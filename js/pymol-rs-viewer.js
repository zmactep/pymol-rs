var f = Object.defineProperty;
var v = (o, e, t) => e in o ? f(o, e, { enumerable: !0, configurable: !0, writable: !0, value: t }) : o[e] = t;
var a = (o, e, t) => v(o, typeof e != "symbol" ? e + "" : e, t);
class w {
  constructor(e) {
    a(this, "container");
    a(this, "pool", []);
    a(this, "activeCount", 0);
    this.container = document.createElement("div"), this.container.style.position = "absolute", this.container.style.top = "0", this.container.style.left = "0", this.container.style.width = "100%", this.container.style.height = "100%", this.container.style.pointerEvents = "none", this.container.style.overflow = "hidden", e.appendChild(this.container);
  }
  update(e, t) {
    for (; this.pool.length < e.length; ) {
      const s = document.createElement("span");
      s.style.position = "absolute", s.style.whiteSpace = "nowrap", s.style.fontFamily = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif", s.style.textShadow = "0 0 3px rgba(0,0,0,0.8), 0 0 6px rgba(0,0,0,0.5)", s.style.display = "none", this.container.appendChild(s), this.pool.push(s);
    }
    for (let s = 0; s < e.length; s++) {
      const i = e[s], n = this.pool[s], r = i.x / t, c = i.y / t;
      n.style.left = `${r}px`, n.style.top = `${c}px`, n.textContent = i.text, i.kind === "atom" ? (n.style.color = "white", n.style.fontSize = "14px", n.style.transform = "translateY(-100%)") : (n.style.color = "#ffff00", n.style.fontSize = "13px", n.style.transform = "translate(-50%, -50%)"), n.style.display = "";
    }
    for (let s = e.length; s < this.activeCount; s++)
      this.pool[s].style.display = "none";
    this.activeCount = e.length;
  }
  destroy() {
    this.container.remove();
  }
}
class y {
  constructor(e) {
    a(this, "wasm", null);
    a(this, "canvas");
    a(this, "animFrameId", 0);
    a(this, "resizeObserver", null);
    a(this, "labelOverlay");
    a(this, "_deferred", !1);
    a(this, "_revealDuration", 150);
    this.container = e, getComputedStyle(e).position === "static" && (e.style.position = "relative"), this.canvas = document.createElement("canvas"), this.canvas.id = "pymol-rs-canvas-" + Math.random().toString(36).slice(2, 8), this.canvas.style.width = "100%", this.canvas.style.height = "100%", this.canvas.style.display = "block", this.canvas.tabIndex = 0, e.appendChild(this.canvas), this.labelOverlay = new w(e);
  }
  async init() {
    const e = await import("./pymol_web-mVc0P39q.js");
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
    e.addEventListener("mousedown", (t) => {
      var s;
      t.preventDefault(), e.focus(), (s = this.wasm) == null || s.on_mouse_down(t.offsetX, t.offsetY, t.button, u(t));
    }), e.addEventListener("mousemove", (t) => {
      var s;
      (s = this.wasm) == null || s.on_mouse_move(t.offsetX, t.offsetY, u(t));
    }), e.addEventListener("mouseup", (t) => {
      var s;
      (s = this.wasm) == null || s.on_mouse_up(t.offsetX, t.offsetY, t.button);
    }), e.addEventListener("wheel", (t) => {
      var s;
      t.preventDefault(), (s = this.wasm) == null || s.on_wheel(t.deltaY, u(t));
    }, { passive: !1 }), e.addEventListener("contextmenu", (t) => t.preventDefault()), this.resizeObserver = new ResizeObserver(() => {
      this.syncCanvasSize(), this.wasm && this.wasm.resize(this.canvas.width, this.canvas.height);
    }), this.resizeObserver.observe(this.container);
  }
  syncCanvasSize() {
    const e = window.devicePixelRatio || 1, t = this.canvas.getBoundingClientRect();
    this.canvas.width = Math.round(t.width * e), this.canvas.height = Math.round(t.height * e);
  }
}
function u(o) {
  let e = 0;
  return o.shiftKey && (e |= 1), o.ctrlKey && (e |= 2), o.altKey && (e |= 4), o.metaKey && (e |= 8), e;
}
class b {
  constructor() {
    a(this, "listeners", /* @__PURE__ */ new Map());
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
class g {
  constructor(e, t) {
    a(this, "container");
    a(this, "viewer");
    a(this, "output");
    a(this, "input");
    a(this, "history", []);
    a(this, "historyIdx", -1);
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
const x = ["lines", "sticks", "cartoon", "spheres", "surface"];
class L {
  constructor(e, t) {
    a(this, "container");
    a(this, "viewer");
    a(this, "list");
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
      const r = document.createElement("button");
      r.className = `obj-vis ${i.enabled ? "enabled" : "disabled"}`, r.textContent = i.enabled ? "V" : "-", r.title = i.enabled ? "Hide" : "Show", r.addEventListener("click", () => {
        this.viewer.execute(i.enabled ? `disable ${s}` : `enable ${s}`);
      });
      const c = document.createElement("span");
      if (c.className = "obj-name", c.textContent = i.object_type === "map" ? `${s} (map)` : `${s} (${i.atom_count})`, n.appendChild(r), n.appendChild(c), i.object_type !== "map") {
        const l = document.createElement("span");
        l.className = "obj-reps";
        for (const d of x) {
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
        const r = document.createElement("button");
        r.className = `obj-vis selection-vis ${i.visible ? "enabled" : "disabled"}`, r.textContent = i.visible ? "V" : "-", r.title = i.visible ? "Hide indicators" : "Show indicators", r.addEventListener("click", () => {
          this.viewer.execute(`toggle ${i.name}`);
        });
        const c = document.createElement("span");
        c.className = "obj-name selection-name", c.textContent = `(${i.name})`, c.title = i.expression;
        const l = document.createElement("button");
        l.className = "rep-btn selection-delete", l.textContent = "X", l.title = "Delete selection", l.addEventListener("click", () => {
          this.viewer.execute(`deselect ${i.name}`);
        }), n.appendChild(r), n.appendChild(c), n.appendChild(l), this.list.appendChild(n);
      }
    }
  }
  destroy() {
    this.container.innerHTML = "";
  }
}
class C {
  constructor(e, t) {
    a(this, "container");
    a(this, "viewer");
    a(this, "content");
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
      for (const r of t.residues) {
        const c = document.createElement("span");
        if (c.className = "seq-residue", c.textContent = r.one_letter, c.title = `${r.resn} ${r.resv}`, c.addEventListener("click", () => {
          const l = `${t.object_name} and chain ${t.chain_id} and resi ${r.resv}`;
          this.viewer.execute(`select sele, ${l}`);
        }), n.appendChild(c), r.resv % 10 === 0) {
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
class _ {
  constructor(e, t) {
    a(this, "container");
    a(this, "viewer");
    a(this, "frameLabel", null);
    a(this, "slider", null);
    a(this, "playBtn", null);
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
class E {
  constructor(e, t = {}) {
    a(this, "core");
    a(this, "events", new b());
    a(this, "panels", /* @__PURE__ */ new Map());
    a(this, "options");
    this.options = t, this.core = new y(e), t.defer && this.core.setDeferred(!0, t.revealDuration ?? 150);
  }
  async init() {
    await this.core.init();
    const e = [];
    if (this.options.layout)
      e.push(...this.options.layout);
    else if (this.options.panels)
      for (const t of this.options.panels)
        e.push({ name: t, slot: "right" });
    for (const t of e) {
      let s;
      if (this.options.slots && (s = this.options.slots[t.slot]), s || (s = document.getElementById("sidebar")), !s) continue;
      const i = document.createElement("div");
      i.className = `pymol-panel pymol-panel-${t.name}`, t.collapsed && i.classList.add("collapsed"), s.appendChild(i);
      let n;
      switch (t.name) {
        case "repl":
          n = new g(i, this);
          break;
        case "objects":
          n = new L(i, this);
          break;
        case "sequence":
          n = new C(i, this);
          break;
        case "movie":
          n = new _(i, this);
          break;
      }
      this.panels.set(t.name, n);
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
        const s = D(t.code, t.format);
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
    let r = new URL(e, location.href).pathname.split("/").pop() ?? "structure";
    r.toLowerCase().endsWith(".gz") && (r = r.slice(0, -3));
    const c = (t == null ? void 0 : t.name) ?? r.replace(/\.[^.]+$/, ""), l = (t == null ? void 0 : t.format) ?? ((d = r.split(".").pop()) == null ? void 0 : d.toLowerCase()) ?? "pdb";
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
function $(o) {
  const e = [], t = {};
  for (const s of o.split(",")) {
    const i = s.trim();
    if (!i) continue;
    const n = i.indexOf("=");
    n !== -1 ? t[i.slice(0, n).trim().toLowerCase()] = i.slice(n + 1).trim() : e.push(i);
  }
  return { positional: e, named: t };
}
function p(o) {
  const e = o.match(/^\s*(\w+)\s+(.*)/s);
  if (!e) return null;
  const t = e[1].toLowerCase(), { positional: s, named: i } = $(e[2]);
  if (t === "fetch" && s.length >= 1) {
    const n = s[0], r = s[1] ?? i.name ?? n, c = s[2] ?? i.type ?? "bcif", l = S(c);
    return { kind: "fetch", code: n, name: r, format: l };
  }
  if (t === "load" && s.length >= 1) {
    const n = s[0];
    if (!/^https?:\/\//i.test(n)) return null;
    const r = s[1] ?? i.object ?? i.name ?? O(n), c = i.format ?? j(n);
    return { kind: "load", url: n, name: r, format: c };
  }
  return null;
}
function S(o) {
  switch (o.toLowerCase()) {
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
const M = "https://models.rcsb.org", k = "https://files.rcsb.org/download";
function D(o, e) {
  const t = o.toLowerCase();
  return e === "bcif" ? `${M}/${t}.bcif.gz` : `${k}/${t}.${e}.gz`;
}
function O(o) {
  let e = new URL(o, location.href).pathname.split("/").pop() ?? "structure";
  return e.toLowerCase().endsWith(".gz") && (e = e.slice(0, -3)), e.replace(/\.[^.]+$/, "");
}
function j(o) {
  var t;
  let e = new URL(o, location.href).pathname.split("/").pop() ?? "";
  return e.toLowerCase().endsWith(".gz") && (e = e.slice(0, -3)), ((t = e.split(".").pop()) == null ? void 0 : t.toLowerCase()) ?? "pdb";
}
function A(o = "pymol-rs-viewer") {
  customElements.get(o) || customElements.define(
    o,
    class extends HTMLElement {
      constructor() {
        super(...arguments);
        a(this, "viewer", null);
      }
      async connectedCallback() {
        const t = this.attachShadow({ mode: "open" }), s = document.createElement("div");
        s.style.width = "100%", s.style.height = "100%", s.style.position = "relative", t.appendChild(s);
        const i = this.getAttribute("panels"), n = i ? i.split(",").map((h) => h.trim()) : [], r = this.getAttribute("src"), c = this.getAttribute("command"), l = this.getAttribute("defer"), d = l !== null ? l !== "false" : !!(r || c);
        if (this.viewer = new E(s, { panels: n, defer: d }), await this.viewer.init(), r) {
          const h = r.split(/\s+/).filter(Boolean);
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
  E as PyMolRSViewer,
  A as registerElement
};
