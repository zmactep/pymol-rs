# Plan: Python как динамический плагин PyMOL-RS

## Context

Python-поддержка сейчас реализована в `crates/pymol-python/` — PyO3 wheel, общающийся с GUI через IPC. Проблемы:
- Жёсткая связанность: крейт зависит от всех остальных крейтов
- IPC добавляет latency и ограничивает доступ к данным (нет прямого доступа к атомам/цепям/остаткам)
- Нет встроенного Python в GUI — нельзя писать Python-выражения в командной строке
- `runpy` работает через двойной IPC roundtrip

**Цель:** бесшовный опыт как в оригинальном PyMOL — скриптовый движок в командной строке, `run script.py`, прямой доступ к молекулярным данным, cmd.* API без IPC.

**Принципы:**
- **Ядро не знает о Python.** Никаких заглушек, fallback'ов или special-case'ов для Python в `pymol-cmd` / `pymol-gui`. Всё делает плагин через существующие механизмы (команды, dynamic commands, MessageHandler).
- **Единый Python API.** Один и тот же скрипт работает и через `run script.py` (в GUI), и через `python script.py` (standalone). Пакет `pymol_rs` автоматически определяет режим.

---

## Архитектура: два компонента

### 1. `plugins/python/` — Встроенный Python-интерпретатор (новый плагин)

Динамический плагин (cdylib). Встраивает Python через PyO3 `auto-initialize`.

**При загрузке регистрирует через существующий plugin API:**
- `python <expr>` — команда для выполнения Python-выражений. Алиас: `/`
- Регистрирует себя как обработчик `.py` файлов для builtin `run` (через generic file handler hook — см. ниже)

Команда `python` — обычная `impl Command`, регистрируется через `PluginRegistrar::register_command()`.
Алиас `/` задаётся через `fn aliases(&self) -> &[&str] { &["/"] }`.

**Реализует `MessageHandler` с `needs_poll() = true`:**
- Доступ к `SharedContext` для чтения молекулярных данных (объекты, атомы, цепи, остатки)
- `ctx.execute_command()` для записи (color, show, hide, load)
- Результаты команд приходят в следующем poll-цикле

**Доступ к данным (во время `poll()`):**
- **Чтение:** `SharedContext.registry` → `MoleculeObject` → `ObjectMolecule` → атомы, цепи, остатки. Данные клонируются в PyO3-обёртки (snapshot, как текущий `PyAtom`).
- **Запись:** через `ctx.execute_command()` — `cmd.color("red", "chain A")`, `cmd.show("cartoon")`, etc.

**Файловая структура:**
```
plugins/python/
  Cargo.toml          # cdylib; deps: pymol-plugin, pyo3 (auto-initialize), pymol-mol, pymol-select
  src/
    lib.rs            # pymol_plugin!, регистрация команд python, run, /
    engine.rs         # Python interpreter lifecycle, with_gil, eval/exec
    commands.rs       # PythonCommand, RunCommand, SlashCommand (impl Command)
    handler.rs        # PythonHandler: impl MessageHandler (poll-based)
    backend.rs        # PluginBackend — бэкенд для pymol_rs пакета внутри embedded Python
```

**Команда `python` (с алиасом `/`):**

```rust
// commands.rs

/// `python <expr>` — evaluate a Python expression/statement
/// Alias: `/` (e.g., `/print(1+1)`)
pub struct PythonCommand { engine: Arc<Mutex<PythonEngine>> }

impl Command for PythonCommand {
    fn name(&self) -> &str { "python" }
    fn aliases(&self) -> &[&str] { &["/"] }
    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
        let expr = args.rest();  // всё после "python" / "/"
        self.engine.lock().eval(expr)?;
        Ok(())
    }
}
```

**Регистрация в плагине:**
```rust
pymol_plugin! {
    name: "python",
    version: "0.1.0",
    description: "Embedded Python interpreter",
    register: |reg| {
        match PythonEngine::try_init() {
            Ok(engine) => {
                let engine = Arc::new(Mutex::new(engine));
                reg.register_command(PythonCommand { engine: engine.clone() });
                reg.set_message_handler(PythonHandler::new(engine));
                // Register .py file handler for the builtin `run` command
                reg.register_file_handler("py", Box::new(move |path| {
                    engine.lock().exec_file(path)
                }));
            }
            Err(e) => log::warn!("Python plugin: {}", e),
        }
    },
}
```

### 2. `python/` — Unified Python-пакет `pymol_rs`

Pip-installable пакет. PyO3 `extension-module`, собирается maturin. **Без IPC.**

**Ключевая идея: один пакет, два бэкенда.**

Скрипт пользователя всегда одинаков:
```python
from pymol_rs import cmd

cmd.load("protein.pdb")
model = cmd.get_model("protein")
for atom in model.atoms():
    print(f"{atom.chain}:{atom.resn}{atom.resi}.{atom.name} ({atom.coord})")
cmd.color("red", "chain A")
```

Пакет автоматически определяет режим запуска:

```python
# pymol_rs/__init__.py
import sys

if hasattr(sys, '_pymolrs_backend'):
    # Embedded mode: running inside PyMOL-RS GUI (plugin set this up)
    _backend = sys._pymolrs_backend
else:
    # Standalone mode: own Session + CommandExecutor
    from ._pymol_rs import _create_standalone_backend
    _backend = _create_standalone_backend()

cmd = Cmd(_backend)
```

**Бэкенд-интерфейс (Python protocol):**
```python
class CmdBackend:
    """Abstract backend — implemented in Rust (PyO3)."""
    def execute(self, command: str, silent: bool = False) -> None: ...
    def get_names(self) -> list[str]: ...
    def get_model(self, name: str) -> ObjectMolecule: ...
    def count_atoms(self, selection: str) -> int: ...
    def get_view(self) -> list[float]: ...
```

**Standalone backend** (в `python/src/backend.rs`):
```rust
#[pyclass]
struct StandaloneBackend {
    session: Session,
    executor: CommandExecutor,
    needs_redraw: bool,
}

#[pymethods]
impl StandaloneBackend {
    fn execute(&mut self, command: &str, silent: bool) -> PyResult<()> {
        let mut adapter = SessionAdapter {
            session: &mut self.session,
            render_context: None,  // headless
            default_size: (1024, 768),
            needs_redraw: &mut self.needs_redraw,
            async_fetch_fn: None,
        };
        self.executor.do_with_options(&mut adapter, command, true, silent)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    fn get_model(&self, name: &str) -> PyResult<PyObjectMolecule> {
        let mol_obj = self.session.registry.get_molecule(name)
            .ok_or_else(|| PyKeyError::new_err(format!("Object '{}' not found", name)))?;
        Ok(PyObjectMolecule::new(mol_obj.molecule().clone()))
    }

    fn get_names(&self) -> Vec<String> {
        self.session.registry.names().map(|s| s.to_string()).collect()
    }

    fn count_atoms(&self, selection: &str) -> usize { /* ... */ }
}
```

**Plugin backend** (в `plugins/python/src/backend.rs`):
```rust
#[pyclass]
struct PluginBackend {
    // Snapshot of SharedContext data, refreshed each poll()
    names: Arc<Mutex<Vec<String>>>,
    molecules: Arc<Mutex<HashMap<String, ObjectMolecule>>>,  // cloned snapshots
    cmd_queue: Arc<Mutex<Vec<CommandExecRequest>>>,
}

#[pymethods]
impl PluginBackend {
    fn execute(&self, command: &str, silent: bool) -> PyResult<()> {
        self.cmd_queue.lock().push(CommandExecRequest { command: command.into(), silent });
        Ok(())
    }

    fn get_model(&self, name: &str) -> PyResult<PyObjectMolecule> {
        let mols = self.molecules.lock();
        let mol = mols.get(name)
            .ok_or_else(|| PyKeyError::new_err(format!("Object '{}' not found", name)))?;
        Ok(PyObjectMolecule::new(mol.clone()))
    }

    fn get_names(&self) -> Vec<String> {
        self.names.lock().clone()
    }
}
```

Плагин при `poll()`:
1. Обновляет shared snapshots из `SharedContext` (names, molecules)
2. Устанавливает `sys._pymolrs_backend = PluginBackend(...)` (один раз)
3. Выполняет Python-код (скрипты, eval) — Python импортирует `pymol_rs` и находит backend
4. После выполнения сливает `cmd_queue` → `ctx.execute_command()`

**`Cmd` класс (Python, одинаковый для обоих режимов):**
```python
# pymol_rs/cmd.py
class Cmd:
    def __init__(self, backend):
        self._backend = backend

    def load(self, filename, object=None, state=0, format=None, quiet=True):
        cmd_str = f"load {filename}"
        if object: cmd_str += f", {object}"
        self._backend.execute(cmd_str, quiet)

    def color(self, color, selection="all"):
        self._backend.execute(f"color {color}, {selection}")

    def show(self, representation, selection="all"):
        self._backend.execute(f"show {representation}, {selection}")

    def hide(self, representation, selection="all"):
        self._backend.execute(f"hide {representation}, {selection}")

    def get_model(self, name):
        return self._backend.get_model(name)

    def get_names(self):
        return self._backend.get_names()

    def count_atoms(self, selection="all"):
        return self._backend.count_atoms(selection)

    # ... etc — тонкая обёртка над backend
```

**Молекулярные типы (одинаковые в обоих режимах):**
```python
model = cmd.get_model("protein")

# Объект
model.name                    # "protein"
model.atom_count              # 1542
model.get_chains()            # ["A", "B"]
model.get_residue_names()     # ["ALA", "GLY", ...]

# Атомы
for atom in model.atoms():
    atom.name                 # "CA"
    atom.chain                # "A"
    atom.resn                 # "ALA"
    atom.resi                 # "42"
    atom.coord                # (1.0, 2.0, 3.0)
    atom.b_factor             # 15.2
    atom.element.symbol       # "C"

# Селекции
from pymol_rs import selecting
sel = selecting.select(model, "name CA and chain A")
sel.count                     # 128
sel.indices()                 # [0, 4, 11, ...]

# Раскраска, представления
cmd.color("red", "chain A")
cmd.show("cartoon", "all")
cmd.hide("lines", "all")
```

**Файловая структура:**
```
python/
  pyproject.toml
  Cargo.toml                # cdylib; deps: pyo3 (extension-module), pymol-scene, pymol-cmd, pymol-mol, pymol-io, pymol-select, pymol-color, pymol-settings
  src/
    lib.rs                  # PyO3 module init
    backend.rs              # StandaloneBackend (Session + CommandExecutor)
    convert.rs              # Type conversions (перенос)
    error.rs                # Error mapping (перенос)
    mol/                    # PyObjectMolecule, PyAtom, PyBond, ... (перенос)
    selecting/              # Selection bindings (перенос)
    color/                  # Color bindings (перенос)
    io/                     # File I/O bindings (перенос)
    settings/               # Settings bindings (перенос)
  pymol_rs/
    __init__.py             # Auto-detect backend, create cmd
    cmd.py                  # Cmd class — unified API wrapper
    _cli.py                 # CLI entry point (упрощённый)
```

---

## Изменения в существующем коде

Через существующий plugin API:
- `python` (с алиасом `/`) — обычная команда через `PluginRegistrar::register_command()`
- Доступ к данным — через `SharedContext` в `MessageHandler::poll()`
- Выполнение команд — через `PollContext::execute_command()`

**Минимальные generic изменения в ядре (не Python-specific):**

### A. File handler registry для `run` команды

Builtin `run` сейчас обрабатывает только `.pml`. Добавляем generic механизм: плагины могут регистрировать обработчики по расширению файла.

**Файл:** `crates/pymol-plugin/src/registrar.rs`
```rust
pub struct PluginRegistrar {
    // ... existing fields ...
    pub(crate) file_handlers: Vec<(String, Box<dyn Fn(&str) -> Result<(), String> + Send + Sync>)>,
}

impl PluginRegistrar {
    /// Register a file handler for a specific extension.
    /// Used by the builtin `run` command to dispatch non-.pml files.
    pub fn register_file_handler(
        &mut self,
        extension: &str,
        handler: Box<dyn Fn(&str) -> Result<(), String> + Send + Sync>,
    ) {
        self.file_handlers.push((extension.to_string(), handler));
    }
}
```

**Файл:** `crates/pymol-cmd/src/commands/control.rs` — `RunCommand::execute`:
```rust
fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
    let filename = args.get_str(0).ok_or(CmdError::MissingArgument("filename".into()))?;
    let path = expand_path(filename);
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("pml");

    match ext {
        "pml" => {
            let mut engine = ScriptEngine::new();
            engine.run_pml(ctx.viewer, &path)
        }
        other => {
            // Delegate to registered file handler
            if let Some(handler) = ctx.file_handlers.get(other) {
                handler(path.to_str().unwrap())
                    .map_err(|e| CmdError::Generic(e))
            } else {
                Err(CmdError::Generic(format!("No handler for .{} files", other)))
            }
        }
    }
}
```

`CommandContext` получает `file_handlers: &HashMap<String, ...>` — передаётся при создании контекста. Это полностью generic: ядро не знает ни о Python, ни о Lua. Плагин регистрирует `"py"` → handler, другой плагин может зарегистрировать `"lua"` → handler.

### B. Парсер: `/` как алиас команды

Нужно проверить, что парсер корректно обрабатывает `/` как имя команды (через `aliases()`). Если нет — минимальная правка в парсере для спецсимволов (generic, не Python-specific).

---

## Поток данных

### `python print(1+1)` в GUI:
```
"python print(1+1)" → executor → registry.get("python") → PythonCommand →
→ engine.eval("print(1+1)") → with_gil { exec } → stdout captured →
→ ctx.print_info("2")
```

### `/import sys; print(sys.version)` в GUI:
```
"/import sys; ..." → executor → registry.get("/") (alias → "python") →
→ PythonCommand → engine.eval("import sys; ...") → output
```

### `run script.py` в GUI:
```
"run script.py" → executor → builtin RunCommand → ext = "py" →
→ file_handlers["py"](path) → PythonEngine::exec_file →
→ with_gil { exec(script) } → script: `from pymol_rs import cmd` →
→ sys._pymolrs_backend → PluginBackend →
→ cmd.load("x.pdb") → cmd_queue → poll() → host executes
→ cmd.get_model("x") → shared snapshot → clone
```

### `run setup.pml` в GUI (без изменений):
```
"run setup.pml" → executor → builtin RunCommand → ext = "pml" →
→ ScriptEngine::run_pml → execute line by line (как раньше)
```

### Тот же Python-скрипт standalone (`python script.py`):
```
`from pymol_rs import cmd` →
→ sys._pymolrs_backend не найден → StandaloneBackend(Session) →
→ cmd.load("x.pdb") → executor.do_() на owned Session → direct
→ cmd.get_model("x") → session.registry.get_molecule → clone
```

---

## Фазы реализации

### Phase 1: Python-пакет без IPC
1. Перенести `crates/pymol-python/` → `python/`
2. Создать `StandaloneBackend` в `python/src/backend.rs` (Session + CommandExecutor)
3. Создать `pymol_rs/cmd.py` — unified Cmd class
4. Обновить `pymol_rs/__init__.py` — auto-detect backend
5. Добавить `get_model()` → `PyObjectMolecule` через `session.registry`
6. Перенести без изменений: `mol/`, `selecting/`, `color/`, `io/`, `settings/`, `convert.rs`, `error.rs`
7. Удалить: `connection.rs`, `ipc/`, `scripting.rs`
8. Обновить workspace `Cargo.toml`
9. Проверить: `cd python && maturin develop && python -c "from pymol_rs import cmd; ..."`

### Phase 2: Generic file handler в ядре
1. Добавить `register_file_handler()` в `PluginRegistrar` (`crates/pymol-plugin/src/registrar.rs`)
2. Расширить `RunCommand` — dispatch по расширению файла (`crates/pymol-cmd/src/commands/control.rs`)
3. Добавить `file_handlers` в `CommandContext` (`crates/pymol-cmd/src/command.rs`)
4. Wire file handlers из PluginManager в CommandContext (`crates/pymol-gui/src/plugin_manager.rs`)
5. Проверить парсер: `/` как алиас команды (generic правка если нужно)
6. `cargo test` — всё проходит, `run setup.pml` работает как раньше

### Phase 3: Python-плагин
1. `plugins/python/Cargo.toml` — cdylib, PyO3 `auto-initialize`, pymol-plugin, pymol-mol, pymol-select
2. `engine.rs` — PythonEngine (init, eval, exec_file)
3. `commands.rs` — PythonCommand (с алиасом `/`)
4. `backend.rs` — PluginBackend (SharedContext snapshot + cmd_queue)
5. `handler.rs` — PythonHandler: MessageHandler (poll: update snapshots, drain cmd_queue)
6. `lib.rs` — `pymol_plugin!` с register_command(PythonCommand) + register_file_handler("py", ...) + set_message_handler
7. Добавить в workspace, проверить: `make run` → `python print("hello")` → `run test.py`

### Phase 4: Интеграция
1. Проверить единый скрипт: один `test.py` работает через `run test.py` и `python test.py`
2. Проверить доступ к данным: `cmd.get_model()`, итерация атомов, цепей, остатков
3. Проверить раскраску: `cmd.color("red", "chain A")` в обоих режимах

### Phase 5: Очистка
1. Удалить `crates/pymol-python/`
2. Обновить `CLAUDE.md`, `Makefile`

---

## Критичные файлы

| Файл | Роль |
|------|------|
| `crates/pymol-plugin/src/registrar.rs` | Добавить `register_file_handler()` (generic) |
| `crates/pymol-cmd/src/commands/control.rs` | Расширить `RunCommand` — dispatch по расширению |
| `crates/pymol-cmd/src/command.rs` | Добавить `file_handlers` в `CommandContext` |
| `crates/pymol-scene/src/session.rs` | Основа StandaloneBackend |
| `crates/pymol-scene/src/session_adapter.rs` | Мост Session → ViewerLike |
| `crates/pymol-gui/src/plugin_manager.rs` | Wire file handlers из плагинов |
| `crates/pymol-python/src/mol/` | Перенести в python/src/mol/ |
| `crates/pymol-python/src/scripting.rs` | Логика exec_file → перенести в plugins/python |
| `plugins/ipc/src/handler.rs` | Образец MessageHandler с poll() |

## Риски

| Риск | Митигация |
|------|-----------|
| `/` не парсится как имя команды | Проверить парсер; если нужно — минимальная generic правка для спецсимволов |
| GIL блокирует GUI при долгих скриптах |  background thread + `allow_threads` |
| `get_model` клонирует ObjectMolecule | Приемлемо. Для >100K атомов — lazy API позже |
| cmd в плагине асинхронен (next frame) | Для большинства скриптов незаметно |
| libpython не найден | Плагин не загружается, всё остальное работает |
| Два бэкенда с разным поведением | Unified Cmd class + integration tests |

## Верификация

1. `cargo build` — workspace без ошибок
2. `cargo test` — все тесты
3. **Standalone:** `cd python && maturin develop && python -c "from pymol_rs import cmd; cmd.load('protein.pdb'); m = cmd.get_model('protein'); print(m.atom_count)"`
4. **Атомы:** `python -c "...; [print(a.name, a.chain, a.coord) for a in m.atoms()]"`
5. **Раскраска:** `python -c "...; cmd.color('red', 'chain A')"`
6. **Плагин:** `make run` → `python print("hello")` → "hello"
7. **Run .py:** `make run` → `run test.py` → Python-скрипт выполняется с `cmd`
8. **Run .pml:** `make run` → `run setup.pml` → PML работает как раньше (регрессия)
9. **Slash:** `make run` → `/import sys; print(sys.version)` → версия
10. **Единый скрипт:** `test.py` работает и через `run test.py` в GUI, и через `python test.py` standalone
