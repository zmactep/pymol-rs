var f = Object.defineProperty;
var v = (r, e, t) => e in r ? f(r, e, { enumerable: !0, configurable: !0, writable: !0, value: t }) : r[e] = t;
var o = (r, e, t) => v(r, typeof e != "symbol" ? e + "" : e, t);
class w {
  constructor(e) {
    o(this, "container");
    o(this, "pool", []);
    o(this, "activeCount", 0);
    this.container = document.createElement("div"), this.container.style.position = "absolute", this.container.style.top = "0", this.container.style.left = "0", this.container.style.width = "100%", this.container.style.height = "100%", this.container.style.pointerEvents = "none", this.container.style.overflow = "hidden", e.appendChild(this.container);
  }
  update(e, t) {
    for (; this.pool.length < e.length; ) {
      const s = document.createElement("span");
      s.style.position = "absolute", s.style.whiteSpace = "nowrap", s.style.fontFamily = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif", s.style.textShadow = "0 0 3px rgba(0,0,0,0.8), 0 0 6px rgba(0,0,0,0.5)", s.style.display = "none", this.container.appendChild(s), this.pool.push(s);
    }
    for (let s = 0; s < e.length; s++) {
      const i = e[s], n = this.pool[s], a = i.x / t, c = i.y / t;
      n.style.left = `${a}px`, n.style.top = `${c}px`, n.textContent = i.text, i.kind === "atom" ? (n.style.color = "white", n.style.fontSize = "14px", n.style.transform = "translateY(-100%)") : (n.style.color = "#ffff00", n.style.fontSize = "13px", n.style.transform = "translate(-50%, -50%)"), n.style.display = "";
    }
    for (let s = e.length; s < this.activeCount; s++)
      this.pool[s].style.display = "none";
    this.activeCount = e.length;
  }
  destroy() {
    this.container.remove();
  }
}
const y = 25;
class b {
  constructor(e) {
    o(this, "wasm", null);
    o(this, "canvas");
    o(this, "animFrameId", 0);
    o(this, "resizeObserver", null);
    o(this, "labelOverlay");
    o(this, "_deferred", !1);
    o(this, "_revealDuration", 150);
    o(this, "dpr", 1);
    o(this, "clickStart", null);
    /** Called with the raw WASM pick result (PickHitInfo | null) on each click. */
    o(this, "onPick", null);
    this.container = e, getComputedStyle(e).position === "static" && (e.style.position = "relative"), this.canvas = document.createElement("canvas"), this.canvas.id = "pymol-rs-canvas-" + Math.random().toString(36).slice(2, 8), this.canvas.style.width = "100%", this.canvas.style.height = "100%", this.canvas.style.display = "block", this.canvas.tabIndex = 0, e.appendChild(this.canvas), this.labelOverlay = new w(e);
  }
  async init() {
    const e = await import("./pymol_web-BtKWEboG.js");
    await e.default(), this.syncCanvasSize();
    const t = await e.WebViewer.create(this.canvas.id);
    return this.wasm = t, this.bindEvents(), this.startLoop(), t;
  }
  get wasmViewer() {
    return this.wasm;
  }
  get isDeferred() {
    return this._deferred;
  }
  setDeferred(e, t = 150) {
    this._deferred = e, this._revealDuration = t, e && (this.container.style.opacity = "0");
  }
  reveal() {
    return this._deferred ? (this._deferred = !1, new Promise((e) => {
      this.container.style.transition = `opacity ${this._revealDuration}ms ease-in`, this.container.offsetHeight, this.container.style.opacity = "1";
      let t = !1;
      const s = () => {
        t || (t = !0, this.container.removeEventListener("transitionend", s), this.container.style.transition = "", e());
      };
      this.container.addEventListener("transitionend", s), setTimeout(s, this._revealDuration + 50);
    })) : Promise.resolve();
  }
  destroy() {
    var e;
    this.animFrameId && (cancelAnimationFrame(this.animFrameId), this.animFrameId = 0), (e = this.resizeObserver) == null || e.disconnect(), this.labelOverlay.destroy(), this.canvas.remove(), this.wasm = null;
  }
  // ---------------------------------------------------------------------------
  // Render loop
  // ---------------------------------------------------------------------------
  startLoop() {
    let e = performance.now();
    const t = (s) => {
      if (this.animFrameId = requestAnimationFrame(t), !this.wasm) return;
      const i = Math.min((s - e) / 1e3, 0.1);
      if (e = s, this.wasm.process_input(), this.wasm.update_animations(i), this.wasm.needs_redraw()) {
        this.wasm.render_frame();
        const n = this.wasm.get_labels();
        this.labelOverlay.update(
          n ?? [],
          window.devicePixelRatio || 1
        );
      }
    };
    this.animFrameId = requestAnimationFrame(t);
  }
  // ---------------------------------------------------------------------------
  // Event binding
  // ---------------------------------------------------------------------------
  bindEvents() {
    const e = this.canvas;
    this.dpr = window.devicePixelRatio || 1, e.addEventListener("mousedown", (t) => {
      var s;
      t.preventDefault(), e.focus(), (s = this.wasm) == null || s.on_mouse_down(t.offsetX, t.offsetY, t.button, u(t)), t.button === 0 && (this.clickStart = { x: t.offsetX, y: t.offsetY });
    }), e.addEventListener("mousemove", (t) => {
      var s, i;
      (s = this.wasm) == null || s.on_mouse_move(t.offsetX, t.offsetY, u(t)), (i = this.wasm) == null || i.process_hover(t.offsetX * this.dpr, t.offsetY * this.dpr);
    }), e.addEventListener("mouseup", (t) => {
      var s, i, n;
      if ((s = this.wasm) == null || s.on_mouse_up(t.offsetX, t.offsetY, t.button), t.button === 0 && this.clickStart !== null) {
        const a = t.offsetX - this.clickStart.x, c = t.offsetY - this.clickStart.y;
        if (a * a + c * c < y) {
          const l = ((i = this.wasm) == null ? void 0 : i.pick_at_screen(
            t.offsetX * this.dpr,
            t.offsetY * this.dpr
          )) ?? null;
          (n = this.onPick) == null || n.call(this, l);
        }
        this.clickStart = null;
      }
    }), e.addEventListener("mouseleave", () => {
      var t;
      this.clickStart = null, (t = this.wasm) == null || t.process_hover(-1, -1);
    }), e.addEventListener("wheel", (t) => {
      var s;
      t.preventDefault(), (s = this.wasm) == null || s.on_wheel(t.deltaY, u(t));
    }, { passive: !1 }), e.addEventListener("contextmenu", (t) => t.preventDefault()), this.resizeObserver = new ResizeObserver(() => {
      this.dpr = window.devicePixelRatio || 1, this.syncCanvasSize(), this.wasm && this.wasm.resize(this.canvas.width, this.canvas.height);
    }), this.resizeObserver.observe(this.container);
  }
  syncCanvasSize() {
    const e = window.devicePixelRatio || 1, t = this.canvas.getBoundingClientRect();
    this.canvas.width = Math.round(t.width * e), this.canvas.height = Math.round(t.height * e);
  }
}
function u(r) {
  let e = 0;
  return r.shiftKey && (e |= 1), r.ctrlKey && (e |= 2), r.altKey && (e |= 4), r.metaKey && (e |= 8), e;
}
class g {
  constructor() {
    o(this, "listeners", /* @__PURE__ */ new Map());
  }
  on(e, t) {
    this.listeners.has(e) || this.listeners.set(e, /* @__PURE__ */ new Set()), this.listeners.get(e).add(t);
  }
  off(e, t) {
    var s;
    (s = this.listeners.get(e)) == null || s.delete(t);
  }
  emit(e, t) {
    for (const s of this.listeners.get(e) ?? [])
      try {
        s(t);
      } catch (i) {
        console.error(`Event handler error (${e}):`, i);
      }
  }
}
class x {
  constructor(e, t) {
    o(this, "container");
    o(this, "viewer");
    o(this, "output");
    o(this, "input");
    o(this, "history", []);
    o(this, "historyIdx", -1);
    this.container = e, this.viewer = t, e.innerHTML = `
      <div class="repl-header">Command Line</div>
      <div class="repl-output"></div>
      <div class="repl-input-row">
        <span class="repl-prompt">PyMOL&gt;</span>
        <input class="repl-input" type="text" placeholder="Type a command..." spellcheck="false" autocomplete="off" />
      </div>
    `, this.output = e.querySelector(".repl-output"), this.input = e.querySelector(".repl-input"), this.input.addEventListener("keydown", (s) => this.onKey(s));
  }
  async onKey(e) {
    if (e.key === "Enter") {
      const t = this.input.value.trim();
      if (!t) return;
      this.history.push(t), this.historyIdx = this.history.length, this.input.value = "", this.appendLine(`PyMOL> ${t}`, "cmd");
      const s = await this.viewer.executeAsync(t);
      for (const i of s.messages)
        this.appendLine(i.text, i.level);
    } else e.key === "ArrowUp" ? (e.preventDefault(), this.historyIdx > 0 && (this.historyIdx--, this.input.value = this.history[this.historyIdx])) : e.key === "ArrowDown" && (e.preventDefault(), this.historyIdx < this.history.length - 1 ? (this.historyIdx++, this.input.value = this.history[this.historyIdx]) : (this.historyIdx = this.history.length, this.input.value = ""));
  }
  appendLine(e, t) {
    const s = document.createElement("div");
    s.className = `repl-line repl-${t}`, s.textContent = e, this.output.appendChild(s), this.output.scrollTop = this.output.scrollHeight;
  }
  update() {
  }
  destroy() {
    this.container.innerHTML = "";
  }
}
const _ = ["lines", "sticks", "cartoon", "spheres", "surface"];
class L {
  constructor(e, t) {
    o(this, "container");
    o(this, "viewer");
    o(this, "list");
    this.container = e, this.viewer = t, e.innerHTML = `
      <div class="panel-header">Objects</div>
      <div class="object-list"></div>
    `, this.list = e.querySelector(".object-list");
  }
  update() {
    const e = this.viewer.getObjectNames();
    this.list.innerHTML = "";
    for (const s of e) {
      const i = this.viewer.getObjectInfo(s);
      if (!i) continue;
      const n = document.createElement("div");
      n.className = "object-row";
      const a = document.createElement("button");
      a.className = `obj-vis ${i.enabled ? "enabled" : "disabled"}`, a.textContent = i.enabled ? "V" : "-", a.title = i.enabled ? "Hide" : "Show", a.addEventListener("click", () => {
        this.viewer.execute(i.enabled ? `disable ${s}` : `enable ${s}`);
      });
      const c = document.createElement("span");
      if (c.className = "obj-name", c.textContent = i.object_type === "map" ? `${s} (map)` : `${s} (${i.atom_count})`, n.appendChild(a), n.appendChild(c), i.object_type !== "map") {
        const l = document.createElement("span");
        l.className = "obj-reps";
        for (const d of _) {
          const h = document.createElement("button");
          h.className = "rep-btn", h.textContent = d.charAt(0).toUpperCase(), h.title = d, h.addEventListener("click", () => {
            this.viewer.execute(`show ${d}, ${s}`);
          }), l.appendChild(h);
        }
        n.appendChild(l);
      }
      this.list.appendChild(n);
    }
    const t = this.viewer.getSelectionList();
    if (t.length > 0) {
      const s = document.createElement("hr");
      s.className = "object-list-separator", this.list.appendChild(s);
      for (const i of t) {
        const n = document.createElement("div");
        n.className = "object-row selection-row";
        const a = document.createElement("button");
        a.className = `obj-vis selection-vis ${i.visible ? "enabled" : "disabled"}`, a.textContent = i.visible ? "V" : "-", a.title = i.visible ? "Hide indicators" : "Show indicators", a.addEventListener("click", () => {
          this.viewer.execute(`toggle ${i.name}`);
        });
        const c = document.createElement("span");
        c.className = "obj-name selection-name", c.textContent = `(${i.name})`, c.title = i.expression;
        const l = document.createElement("button");
        l.className = "rep-btn selection-delete", l.textContent = "X", l.title = "Delete selection", l.addEventListener("click", () => {
          this.viewer.execute(`deselect ${i.name}`);
        }), n.appendChild(a), n.appendChild(c), n.appendChild(l), this.list.appendChild(n);
      }
    }
  }
  destroy() {
    this.container.innerHTML = "";
  }
}
class C {
  constructor(e, t) {
    o(this, "container");
    o(this, "viewer");
    o(this, "content");
    this.container = e, this.viewer = t, e.innerHTML = `
      <div class="panel-header">Sequence</div>
      <div class="sequence-content"></div>
    `, this.content = e.querySelector(".sequence-content");
  }
  update() {
    const e = this.viewer.getSequenceData();
    this.content.innerHTML = "";
    for (const t of e) {
      const s = document.createElement("div");
      s.className = "seq-chain";
      const i = document.createElement("div");
      i.className = "seq-chain-header", i.textContent = `${t.object_name} / ${t.chain_id || "?"}`, s.appendChild(i);
      const n = document.createElement("div");
      n.className = "seq-residues";
      for (const a of t.residues) {
        const c = document.createElement("span");
        if (c.className = "seq-residue", c.textContent = a.one_letter, c.title = `${a.resn} ${a.resv}`, c.addEventListener("click", () => {
          const l = `${t.object_name} and chain ${t.chain_id} and resi ${a.resv}`;
          this.viewer.execute(`select sele, ${l}`);
        }), n.appendChild(c), a.resv % 10 === 0) {
          const l = document.createElement("span");
          l.className = "seq-spacer", l.textContent = " ", n.appendChild(l);
        }
      }
      s.appendChild(n), this.content.appendChild(s);
    }
  }
  destroy() {
    this.container.innerHTML = "";
  }
}
class E {
  constructor(e, t) {
    o(this, "container");
    o(this, "viewer");
    o(this, "frameLabel", null);
    o(this, "slider", null);
    o(this, "playBtn", null);
    this.container = e, this.viewer = t, e.innerHTML = `
      <div class="panel-header">Movie</div>
      <div class="movie-controls">
        <div class="movie-buttons">
          <button class="movie-btn" data-action="beginning" title="Beginning">|&lt;</button>
          <button class="movie-btn" data-action="backward" title="Step back">&lt;</button>
          <button class="movie-btn movie-play" data-action="play" title="Play/Pause">&#9654;</button>
          <button class="movie-btn" data-action="forward" title="Step forward">&gt;</button>
          <button class="movie-btn" data-action="end" title="End">&gt;|</button>
        </div>
        <div class="movie-timeline">
          <input type="range" class="movie-slider" min="1" max="1" value="1" />
          <span class="movie-frame-label">1 / 1</span>
        </div>
      </div>
    `, this.frameLabel = e.querySelector(".movie-frame-label"), this.slider = e.querySelector(".movie-slider"), this.playBtn = e.querySelector(".movie-play"), e.querySelectorAll(".movie-btn").forEach((s) => {
      s.addEventListener("click", () => {
        switch (s.dataset.action) {
          case "play":
            this.viewer.execute("mplay");
            break;
          case "beginning":
            this.viewer.execute("frame 1");
            break;
          case "end": {
            const n = this.viewer.getMovieState();
            this.viewer.execute(`frame ${n.frame_count}`);
            break;
          }
          case "forward":
            this.viewer.execute("forward");
            break;
          case "backward":
            this.viewer.execute("backward");
            break;
        }
      });
    }), this.slider.addEventListener("input", () => {
      const s = parseInt(this.slider.value, 10);
      this.viewer.execute(`frame ${s}`);
    });
  }
  update() {
    const e = this.viewer.getMovieState();
    this.slider && (this.slider.max = String(Math.max(1, e.frame_count)), this.slider.value = String(e.current_frame + 1)), this.frameLabel && (this.frameLabel.textContent = `${e.current_frame + 1} / ${Math.max(1, e.frame_count)}`), this.playBtn && (this.playBtn.textContent = e.is_playing ? "⏸" : "▶");
  }
  destroy() {
    this.container.innerHTML = "";
  }
}
class S {
  constructor(e, t = {}) {
    o(this, "core");
    o(this, "events", new g());
    o(this, "panels", /* @__PURE__ */ new Map());
    o(this, "options");
    this.options = t, this.core = new b(e), t.defer && this.core.setDeferred(!0, t.revealDuration ?? 150);
  }
  async init() {
    const e = await this.core.init();
    this.options.picking && e.set_picking_enabled(!0), this.core.onPick = (s) => {
      const i = s ?? {
        object_name: null,
        atom_index: null,
        chain: null,
        residue: null,
        expression: null
      };
      this.events.emit("atom-picked", i), this.refreshPanels();
    };
    const t = [];
    if (this.options.layout)
      t.push(...this.options.layout);
    else if (this.options.panels)
      for (const s of this.options.panels)
        t.push({ name: s, slot: "right" });
    for (const s of t) {
      let i;
      if (this.options.slots && (i = this.options.slots[s.slot]), i || (i = document.getElementById("sidebar")), !i) continue;
      const n = document.createElement("div");
      n.className = `pymol-panel pymol-panel-${s.name}`, s.collapsed && n.classList.add("collapsed"), i.appendChild(n);
      let a;
      switch (s.name) {
        case "repl":
          a = new x(n, this);
          break;
        case "objects":
          a = new L(n, this);
          break;
        case "sequence":
          a = new C(n, this);
          break;
        case "movie":
          a = new E(n, this);
          break;
      }
      this.panels.set(s.name, a);
    }
    this.events.emit("ready", {});
  }
  // ---------------------------------------------------------------------------
  // Deferred display
  // ---------------------------------------------------------------------------
  get isDeferred() {
    return this.core.isDeferred;
  }
  async show() {
    await this.core.reveal();
  }
  // ---------------------------------------------------------------------------
  // Picking
  // ---------------------------------------------------------------------------
  /**
   * Enable or disable cursor-based atom picking at runtime.
   *
   * When enabled, left-click picks atoms, updates the `sele` selection, and
   * fires `atom-picked` events. Can also be set at construction time via
   * `ViewerOptions.picking`.
   */
  setPicking(e) {
    var t;
    (t = this.core.wasmViewer) == null || t.set_picking_enabled(e);
  }
  // ---------------------------------------------------------------------------
  // Commands
  // ---------------------------------------------------------------------------
  execute(e) {
    const t = p(e);
    if (t)
      return this.executeAsync(e), { messages: [{ level: "info", text: t.kind === "fetch" ? ` Fetching ${t.code}...` : ` Loading ${t.url}...` }] };
    const s = this.core.wasmViewer;
    if (!s) return { messages: [] };
    const i = s.execute(e);
    return this.events.emit("command-output", i.messages[0] ?? { level: "info", text: "" }), this.refreshPanels(), i;
  }
  async executeAsync(e) {
    const t = p(e);
    if (!t) return this.execute(e);
    try {
      if (t.kind === "fetch") {
        const s = O(t.code, t.format);
        await this.loadUrl(s, { name: t.name, format: t.format });
        const i = { level: "info", text: ` Fetched ${t.code} as "${t.name}"` };
        return this.events.emit("command-output", i), { messages: [i] };
      } else {
        await this.loadUrl(t.url, { name: t.name, format: t.format });
        const s = { level: "info", text: ` Loaded "${t.name}" from URL` };
        return this.events.emit("command-output", s), { messages: [s] };
      }
    } catch (s) {
      const i = { level: "error", text: ` ${s}` };
      return this.events.emit("command-output", i), { messages: [i] };
    }
  }
  loadData(e, t, s) {
    const i = this.core.wasmViewer;
    i && (i.load_data(e, t, s), this.refreshPanels());
  }
  async loadUrl(e, t) {
    var d;
    const s = await fetch(e);
    if (!s.ok) throw new Error(`Fetch failed: ${s.status} ${s.statusText}`);
    const i = new Uint8Array(await s.arrayBuffer());
    let a = new URL(e, location.href).pathname.split("/").pop() ?? "structure";
    a.toLowerCase().endsWith(".gz") && (a = a.slice(0, -3));
    const c = (t == null ? void 0 : t.name) ?? a.replace(/\.[^.]+$/, ""), l = (t == null ? void 0 : t.format) ?? ((d = a.split(".").pop()) == null ? void 0 : d.toLowerCase()) ?? "pdb";
    this.loadData(i, c, l);
  }
  // ---------------------------------------------------------------------------
  // Queries
  // ---------------------------------------------------------------------------
  getObjectNames() {
    const e = this.core.wasmViewer;
    return e ? e.get_object_names() : [];
  }
  getObjectInfo(e) {
    const t = this.core.wasmViewer;
    return t ? t.get_object_info(e) : null;
  }
  getSequenceData() {
    const e = this.core.wasmViewer;
    return e ? e.get_sequence_data() : [];
  }
  getMovieState() {
    const e = this.core.wasmViewer;
    return e ? e.get_movie_state() : { frame_count: 0, current_frame: 0, is_playing: !1 };
  }
  getSelectionList() {
    const e = this.core.wasmViewer;
    return e ? e.get_selection_list() : [];
  }
  countAtoms(e = "all") {
    const t = this.core.wasmViewer;
    return t ? t.count_atoms(e) : 0;
  }
  // ---------------------------------------------------------------------------
  // Events
  // ---------------------------------------------------------------------------
  on(e, t) {
    this.events.on(e, t);
  }
  off(e, t) {
    this.events.off(e, t);
  }
  // ---------------------------------------------------------------------------
  // Panel management
  // ---------------------------------------------------------------------------
  refreshPanels() {
    for (const e of this.panels.values())
      e.update();
    this.events.emit("objects-changed", { names: this.getObjectNames() });
  }
  destroy() {
    for (const e of this.panels.values())
      e.destroy();
    this.panels.clear(), this.core.destroy();
  }
}
function k(r) {
  const e = [], t = {};
  for (const s of r.split(",")) {
    const i = s.trim();
    if (!i) continue;
    const n = i.indexOf("=");
    n !== -1 ? t[i.slice(0, n).trim().toLowerCase()] = i.slice(n + 1).trim() : e.push(i);
  }
  return { positional: e, named: t };
}
function p(r) {
  const e = r.match(/^\s*(\w+)\s+(.*)/s);
  if (!e) return null;
  const t = e[1].toLowerCase(), { positional: s, named: i } = k(e[2]);
  if (t === "fetch" && s.length >= 1) {
    const n = s[0], a = s[1] ?? i.name ?? n, c = s[2] ?? i.type ?? "bcif", l = $(c);
    return { kind: "fetch", code: n, name: a, format: l };
  }
  if (t === "load" && s.length >= 1) {
    const n = s[0];
    if (!/^https?:\/\//i.test(n)) return null;
    const a = s[1] ?? i.object ?? i.name ?? P(n), c = i.format ?? j(n);
    return { kind: "load", url: n, name: a, format: c };
  }
  return null;
}
function $(r) {
  switch (r.toLowerCase()) {
    case "pdb":
      return "pdb";
    case "cif":
    case "mmcif":
      return "cif";
    case "bcif":
    case "binarycif":
      return "bcif";
    default:
      return "bcif";
  }
}
const M = "https://models.rcsb.org", D = "https://files.rcsb.org/download";
function O(r, e) {
  const t = r.toLowerCase();
  return e === "bcif" ? `${M}/${t}.bcif.gz` : `${D}/${t}.${e}.gz`;
}
function P(r) {
  let e = new URL(r, location.href).pathname.split("/").pop() ?? "structure";
  return e.toLowerCase().endsWith(".gz") && (e = e.slice(0, -3)), e.replace(/\.[^.]+$/, "");
}
function j(r) {
  var t;
  let e = new URL(r, location.href).pathname.split("/").pop() ?? "";
  return e.toLowerCase().endsWith(".gz") && (e = e.slice(0, -3)), ((t = e.split(".").pop()) == null ? void 0 : t.toLowerCase()) ?? "pdb";
}
function A(r = "pymol-rs-viewer") {
  customElements.get(r) || customElements.define(
    r,
    class extends HTMLElement {
      constructor() {
        super(...arguments);
        o(this, "viewer", null);
      }
      async connectedCallback() {
        const t = this.attachShadow({ mode: "open" }), s = document.createElement("div");
        s.style.width = "100%", s.style.height = "100%", s.style.position = "relative", t.appendChild(s);
        const i = this.getAttribute("panels"), n = i ? i.split(",").map((h) => h.trim()) : [], a = this.getAttribute("src"), c = this.getAttribute("command"), l = this.getAttribute("defer"), d = l !== null ? l !== "false" : !!(a || c);
        if (this.viewer = new S(s, { panels: n, defer: d }), await this.viewer.init(), a) {
          const h = a.split(/\s+/).filter(Boolean);
          for (const m of h)
            await this.viewer.loadUrl(m);
        }
        if (c)
          for (const h of c.split(";")) {
            const m = h.trim();
            m && await this.viewer.executeAsync(m);
          }
        d && await this.viewer.show();
      }
      disconnectedCallback() {
        var t;
        (t = this.viewer) == null || t.destroy(), this.viewer = null;
      }
    }
  );
}
export {
  S as PyMolRSViewer,
  A as registerElement
};
