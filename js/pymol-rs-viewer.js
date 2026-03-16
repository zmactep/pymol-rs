var m = Object.defineProperty;
var u = (r, e, t) => e in r ? m(r, e, { enumerable: !0, configurable: !0, writable: !0, value: t }) : r[e] = t;
var n = (r, e, t) => u(r, typeof e != "symbol" ? e + "" : e, t);
class p {
  constructor(e) {
    n(this, "wasm", null);
    n(this, "canvas");
    n(this, "animFrameId", 0);
    n(this, "resizeObserver", null);
    this.container = e, this.canvas = document.createElement("canvas"), this.canvas.id = "pymol-rs-canvas-" + Math.random().toString(36).slice(2, 8), this.canvas.style.width = "100%", this.canvas.style.height = "100%", this.canvas.style.display = "block", this.canvas.tabIndex = 0, e.appendChild(this.canvas);
  }
  async init() {
    const e = await import("./pymol_web-DZbtLDqN.js");
    await e.default(), this.syncCanvasSize();
    const t = await e.WebViewer.create(this.canvas.id);
    return this.wasm = t, this.bindEvents(), this.startLoop(), t;
  }
  get wasmViewer() {
    return this.wasm;
  }
  destroy() {
    var e;
    this.animFrameId && (cancelAnimationFrame(this.animFrameId), this.animFrameId = 0), (e = this.resizeObserver) == null || e.disconnect(), this.canvas.remove(), this.wasm = null;
  }
  // ---------------------------------------------------------------------------
  // Render loop
  // ---------------------------------------------------------------------------
  startLoop() {
    const e = () => {
      this.animFrameId = requestAnimationFrame(e), this.wasm && (this.wasm.process_input(), this.wasm.needs_redraw() && this.wasm.render_frame());
    };
    this.animFrameId = requestAnimationFrame(e);
  }
  // ---------------------------------------------------------------------------
  // Event binding
  // ---------------------------------------------------------------------------
  bindEvents() {
    const e = this.canvas;
    e.addEventListener("mousedown", (t) => {
      var s;
      t.preventDefault(), e.focus(), (s = this.wasm) == null || s.on_mouse_down(t.offsetX, t.offsetY, t.button, d(t));
    }), e.addEventListener("mousemove", (t) => {
      var s;
      (s = this.wasm) == null || s.on_mouse_move(t.offsetX, t.offsetY, d(t));
    }), e.addEventListener("mouseup", (t) => {
      var s;
      (s = this.wasm) == null || s.on_mouse_up(t.offsetX, t.offsetY, t.button);
    }), e.addEventListener("wheel", (t) => {
      var s;
      t.preventDefault(), (s = this.wasm) == null || s.on_wheel(t.deltaY, d(t));
    }, { passive: !1 }), e.addEventListener("contextmenu", (t) => t.preventDefault()), this.resizeObserver = new ResizeObserver(() => {
      this.syncCanvasSize(), this.wasm && this.wasm.resize(this.canvas.width, this.canvas.height);
    }), this.resizeObserver.observe(this.container);
  }
  syncCanvasSize() {
    const e = window.devicePixelRatio || 1, t = this.canvas.getBoundingClientRect();
    this.canvas.width = Math.round(t.width * e), this.canvas.height = Math.round(t.height * e);
  }
}
function d(r) {
  let e = 0;
  return r.shiftKey && (e |= 1), r.ctrlKey && (e |= 2), r.altKey && (e |= 4), r.metaKey && (e |= 8), e;
}
class v {
  constructor() {
    n(this, "listeners", /* @__PURE__ */ new Map());
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
class f {
  constructor(e, t) {
    n(this, "container");
    n(this, "viewer");
    n(this, "output");
    n(this, "input");
    n(this, "history", []);
    n(this, "historyIdx", -1);
    this.container = e, this.viewer = t, e.innerHTML = `
      <div class="repl-header">Command Line</div>
      <div class="repl-output"></div>
      <div class="repl-input-row">
        <span class="repl-prompt">PyMOL&gt;</span>
        <input class="repl-input" type="text" placeholder="Type a command..." spellcheck="false" autocomplete="off" />
      </div>
    `, this.output = e.querySelector(".repl-output"), this.input = e.querySelector(".repl-input"), this.input.addEventListener("keydown", (s) => this.onKey(s));
  }
  onKey(e) {
    if (e.key === "Enter") {
      const t = this.input.value.trim();
      if (!t) return;
      this.history.push(t), this.historyIdx = this.history.length, this.appendLine(`PyMOL> ${t}`, "cmd");
      const s = this.viewer.execute(t);
      for (const i of s.messages)
        this.appendLine(i.text, i.level);
      this.input.value = "";
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
const w = ["lines", "sticks", "cartoon", "spheres", "surface"];
class b {
  constructor(e, t) {
    n(this, "container");
    n(this, "viewer");
    n(this, "list");
    this.container = e, this.viewer = t, e.innerHTML = `
      <div class="panel-header">Objects</div>
      <div class="object-list"></div>
    `, this.list = e.querySelector(".object-list");
  }
  update() {
    const e = this.viewer.getObjectNames();
    this.list.innerHTML = "";
    for (const t of e) {
      const s = this.viewer.getObjectInfo(t);
      if (!s) continue;
      const i = document.createElement("div");
      i.className = "object-row";
      const a = document.createElement("button");
      a.className = `obj-vis ${s.enabled ? "enabled" : "disabled"}`, a.textContent = s.enabled ? "V" : "-", a.title = s.enabled ? "Hide" : "Show", a.addEventListener("click", () => {
        this.viewer.execute(s.enabled ? `disable ${t}` : `enable ${t}`);
      });
      const o = document.createElement("span");
      o.className = "obj-name", o.textContent = `${t} (${s.atom_count})`;
      const c = document.createElement("span");
      c.className = "obj-reps";
      for (const l of w) {
        const h = document.createElement("button");
        h.className = "rep-btn", h.textContent = l.charAt(0).toUpperCase(), h.title = l, h.addEventListener("click", () => {
          this.viewer.execute(`show ${l}, ${t}`);
        }), c.appendChild(h);
      }
      i.appendChild(a), i.appendChild(o), i.appendChild(c), this.list.appendChild(i);
    }
  }
  destroy() {
    this.container.innerHTML = "";
  }
}
class y {
  constructor(e, t) {
    n(this, "container");
    n(this, "viewer");
    n(this, "content");
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
      const a = document.createElement("div");
      a.className = "seq-residues";
      for (const o of t.residues) {
        const c = document.createElement("span");
        if (c.className = "seq-residue", c.textContent = o.one_letter, c.title = `${o.resn} ${o.resv}`, c.addEventListener("click", () => {
          const l = `${t.object_name} and chain ${t.chain_id} and resi ${o.resv}`;
          this.viewer.execute(`select sele, ${l}`);
        }), a.appendChild(c), o.resv % 10 === 0) {
          const l = document.createElement("span");
          l.className = "seq-spacer", l.textContent = " ", a.appendChild(l);
        }
      }
      s.appendChild(a), this.content.appendChild(s);
    }
  }
  destroy() {
    this.container.innerHTML = "";
  }
}
class g {
  constructor(e, t) {
    n(this, "container");
    n(this, "viewer");
    n(this, "frameLabel", null);
    n(this, "slider", null);
    n(this, "playBtn", null);
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
            const a = this.viewer.getMovieState();
            this.viewer.execute(`frame ${a.frame_count}`);
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
class x {
  constructor(e, t = {}) {
    n(this, "core");
    n(this, "events", new v());
    n(this, "panels", /* @__PURE__ */ new Map());
    n(this, "options");
    this.options = t, this.core = new p(e);
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
      let a;
      switch (t.name) {
        case "repl":
          a = new f(i, this);
          break;
        case "objects":
          a = new b(i, this);
          break;
        case "sequence":
          a = new y(i, this);
          break;
        case "movie":
          a = new g(i, this);
          break;
      }
      this.panels.set(t.name, a);
    }
    this.events.emit("ready", {});
  }
  // ---------------------------------------------------------------------------
  // Commands
  // ---------------------------------------------------------------------------
  execute(e) {
    const t = this.core.wasmViewer;
    if (!t) return { messages: [] };
    const s = t.execute(e);
    return this.events.emit("command-output", s.messages[0] ?? { level: "info", text: "" }), this.refreshPanels(), s;
  }
  loadData(e, t, s) {
    const i = this.core.wasmViewer;
    i && (i.load_data(e, t, s), this.refreshPanels());
  }
  async loadUrl(e, t) {
    var h;
    const s = await fetch(e);
    if (!s.ok) throw new Error(`Fetch failed: ${s.status} ${s.statusText}`);
    const i = new Uint8Array(await s.arrayBuffer()), o = new URL(e, location.href).pathname.split("/").pop() ?? "structure", c = (t == null ? void 0 : t.name) ?? o.replace(/\.[^.]+$/, ""), l = (t == null ? void 0 : t.format) ?? ((h = o.split(".").pop()) == null ? void 0 : h.toLowerCase()) ?? "pdb";
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
function E(r = "pymol-rs-viewer") {
  customElements.get(r) || customElements.define(
    r,
    class extends HTMLElement {
      constructor() {
        super(...arguments);
        n(this, "viewer", null);
      }
      async connectedCallback() {
        const t = this.attachShadow({ mode: "open" }), s = document.createElement("div");
        s.style.width = "100%", s.style.height = "100%", s.style.position = "relative", t.appendChild(s);
        const i = this.getAttribute("panels"), a = i ? i.split(",").map((l) => l.trim()) : [];
        this.viewer = new x(s, { panels: a }), await this.viewer.init();
        const o = this.getAttribute("src");
        o && await this.viewer.loadUrl(o);
        const c = this.getAttribute("command");
        c && this.viewer.execute(c);
      }
      disconnectedCallback() {
        var t;
        (t = this.viewer) == null || t.destroy(), this.viewer = null;
      }
    }
  );
}
export {
  x as PyMolRSViewer,
  E as registerElement
};
