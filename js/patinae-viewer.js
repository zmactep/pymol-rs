var x = Object.defineProperty;
var k = (l, e, t) => e in l ? x(l, e, { enumerable: !0, configurable: !0, writable: !0, value: t }) : l[e] = t;
var a = (l, e, t) => k(l, typeof e != "symbol" ? e + "" : e, t);
class E {
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
const C = 25, L = 240;
class S {
  constructor(e) {
    a(this, "wasm", null);
    a(this, "canvas");
    a(this, "animFrameId", 0);
    a(this, "resizeObserver", null);
    a(this, "labelOverlay");
    a(this, "_deferred", !1);
    a(this, "_revealDuration", 150);
    a(this, "dpr", 1);
    a(this, "clickStart", null);
    a(this, "fatalError", null);
    a(this, "fatalLabel", null);
    a(this, "queuedHover", null);
    a(this, "frameTimes", []);
    a(this, "lastFrameMs", 0);
    a(this, "frameCount", 0);
    a(this, "skipNextFrameTiming", !1);
    /** Called with the raw WASM pick result (PickHitInfo | null) on each click. */
    a(this, "onPick", null);
    /** Called for renderer warnings produced outside direct command execution. */
    a(this, "onOutput", null);
    this.container = e, getComputedStyle(e).position === "static" && (e.style.position = "relative"), this.canvas = document.createElement("canvas"), this.canvas.id = "patinae-canvas-" + Math.random().toString(36).slice(2, 8), this.canvas.style.width = "100%", this.canvas.style.height = "100%", this.canvas.style.display = "block", this.canvas.tabIndex = 0, e.appendChild(this.canvas), this.labelOverlay = new E(e);
  }
  /**
   * Initialise the viewer.
   *
   * @param options.picking — whether to allocate hit-test picking readback
   *   resources (Rg32Uint half-res target + readbacks/reprojection).
   *   Read-only viewers should leave it `false` to save VRAM. Defaults to
   *   `false`. The choice is baked in at WASM construction; changing
   *   later requires re-creating the viewer.
   * @param options.selectionOverlay — whether to draw visible selection /
   *   hover overlays. Defaults to `options.picking`.
   * @param options.memoryProfile — optional renderer memory profile override:
   *   "performance", "balanced", "lite", "manual:<MiB>", or "auto".
   */
  async init(e = {}) {
    const t = await import("./patinae_web-DQxpYDSv.js");
    await t.default(), this.syncCanvasSize();
    const s = e.picking ?? !1, i = e.selectionOverlay ?? s, n = e.memoryProfile ?? "auto", r = t.WebViewer, c = n === "auto" ? await r.create(
      this.canvas.id,
      s,
      i
    ) : await r.createWithMemoryProfile(
      this.canvas.id,
      s,
      i,
      n
    );
    return this.wasm = c, this.fatalError = null, this.fatalLabel = null, this.syncCanvasSize(), this.requireWasm(
      "resize",
      (o) => o.resize(this.canvas.width, this.canvas.height)
    ), this.bindEvents(), this.startLoop(), c;
  }
  get isFailed() {
    return this.fatalError !== null;
  }
  get lastError() {
    return this.fatalError;
  }
  callWasm(e, t) {
    if (this.fatalError !== null) return;
    const s = this.wasm;
    if (s)
      try {
        return t(s);
      } catch (i) {
        this.recordFatal(e, i);
        return;
      }
  }
  requireWasm(e, t) {
    if (this.fatalError !== null)
      throw this.stoppedError(e);
    const s = this.wasm;
    if (!s)
      throw new Error(
        `Patinae viewer is not initialized; cannot call ${e}`
      );
    try {
      return t(s);
    } catch (i) {
      throw this.recordFatal(e, i), i;
    }
  }
  queryWasm(e, t, s) {
    return !this.wasm && this.fatalError === null ? t : this.requireWasm(e, s);
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
  getPerformanceSnapshot() {
    const e = this.frameTimes.length > 0 ? this.frameTimes.reduce((s, i) => s + i, 0) / this.frameTimes.length : 0, t = this.callWasm(
      "get_performance_snapshot",
      (s) => s.get_performance_snapshot()
    ) ?? null;
    return {
      frame_count: this.frameCount,
      avg_frame_ms: e,
      median_frame_ms: y(this.frameTimes, 0.5),
      p95_frame_ms: y(this.frameTimes, 0.95),
      last_frame_ms: this.lastFrameMs,
      wasm: t
    };
  }
  resetPerformanceStats() {
    this.frameTimes = [], this.lastFrameMs = 0, this.frameCount = 0, this.skipNextFrameTiming = !0, this.callWasm(
      "reset_performance_stats",
      (e) => e.reset_performance_stats()
    );
  }
  recordFatal(e, t) {
    this.fatalError === null && (this.fatalError = t, this.fatalLabel = e, this.clickStart = null, this.animFrameId && (cancelAnimationFrame(this.animFrameId), this.animFrameId = 0), console.error(
      `[Patinae] WASM call failed in ${e}; viewer stopped.`,
      t
    ));
  }
  stoppedError(e) {
    const t = this.fatalLabel ?? "unknown";
    return new Error(
      `Patinae viewer stopped after WASM call ${t} failed; cannot call ${e}`,
      { cause: this.fatalError }
    );
  }
  // ---------------------------------------------------------------------------
  // Render loop
  // ---------------------------------------------------------------------------
  startLoop() {
    let e = performance.now();
    const t = (s) => {
      var c;
      if (this.animFrameId = requestAnimationFrame(t), !this.wasm || this.isFailed) return;
      const i = Math.min((s - e) / 1e3, 0.1);
      if (this.lastFrameMs = s - e, this.skipNextFrameTiming ? this.skipNextFrameTiming = !1 : this.lastFrameMs > 0 && (this.frameTimes.push(this.lastFrameMs), this.frameTimes.length > L && this.frameTimes.shift(), this.frameCount++), e = s, this.callWasm("process_input", (o) => o.process_input()), this.isFailed || (this.callWasm("update_animations", (o) => o.update_animations(i)), this.isFailed)) return;
      const n = this.queuedHover;
      if (n && (this.queuedHover = null, this.callWasm(
        "process_hover",
        (o) => o.process_hover(n.x, n.y)
      ), this.isFailed) || (this.callWasm("poll_pending_picks", (o) => o.poll_pending_picks()), this.isFailed)) return;
      if (this.callWasm("needs_redraw", (o) => o.needs_redraw())) {
        if (this.callWasm("render_frame", (h) => h.render_frame()), this.isFailed || (this.drainRendererWarnings(), this.isFailed)) return;
        const o = this.callWasm(
          "get_labels",
          (h) => h.get_labels()
        );
        if (this.isFailed) return;
        this.labelOverlay.update(o ?? [], window.devicePixelRatio || 1);
      }
      const r = this.callWasm(
        "take_completed_pick",
        (o) => o.take_completed_pick()
      );
      this.isFailed || r !== void 0 && ((c = this.onPick) == null || c.call(this, r ?? null));
    };
    this.animFrameId = requestAnimationFrame(t);
  }
  drainRendererWarnings() {
    var t;
    const e = this.callWasm("takeWarnings", (s) => {
      const i = s.takeWarnings;
      return i ? i.call(s) : null;
    });
    if (e)
      for (const s of e.messages ?? [])
        (t = this.onOutput) == null || t.call(this, s);
  }
  // ---------------------------------------------------------------------------
  // Event binding
  // ---------------------------------------------------------------------------
  bindEvents() {
    const e = this.canvas;
    this.dpr = window.devicePixelRatio || 1, e.addEventListener("mousedown", (t) => {
      t.preventDefault(), e.focus(), this.callWasm(
        "on_mouse_down",
        (s) => s.on_mouse_down(t.offsetX, t.offsetY, t.button, f(t))
      ), !this.isFailed && t.button === 0 && (this.clickStart = { x: t.offsetX, y: t.offsetY });
    }), e.addEventListener("mousemove", (t) => {
      this.callWasm(
        "on_mouse_move",
        (s) => s.on_mouse_move(t.offsetX, t.offsetY, f(t))
      ), !this.isFailed && (this.queuedHover = {
        x: t.offsetX * this.dpr,
        y: t.offsetY * this.dpr
      });
    }), e.addEventListener("mouseup", (t) => {
      if (this.callWasm(
        "on_mouse_up",
        (s) => s.on_mouse_up(t.offsetX, t.offsetY, t.button)
      ), !this.isFailed && t.button === 0 && this.clickStart !== null) {
        const s = t.offsetX - this.clickStart.x, i = t.offsetY - this.clickStart.y;
        if (s * s + i * i < C && (this.callWasm(
          "pick_at_screen",
          (n) => n.pick_at_screen(t.offsetX * this.dpr, t.offsetY * this.dpr)
        ), this.isFailed))
          return;
        this.clickStart = null;
      }
    }), e.addEventListener("mouseleave", () => {
      this.clickStart = null, this.queuedHover = { x: -1, y: -1 };
    }), e.addEventListener(
      "wheel",
      (t) => {
        t.preventDefault(), this.callWasm(
          "on_wheel",
          (s) => s.on_wheel(t.deltaY, f(t))
        );
      },
      { passive: !1 }
    ), e.addEventListener("contextmenu", (t) => t.preventDefault()), this.resizeObserver = new ResizeObserver(() => {
      this.dpr = window.devicePixelRatio || 1, this.syncCanvasSize(), this.callWasm(
        "resize",
        (t) => t.resize(this.canvas.width, this.canvas.height)
      );
    }), this.resizeObserver.observe(this.container);
  }
  syncCanvasSize() {
    const e = window.devicePixelRatio || 1, t = this.canvas.getBoundingClientRect();
    this.canvas.width = Math.round(t.width * e), this.canvas.height = Math.round(t.height * e);
  }
}
function f(l) {
  let e = 0;
  return l.shiftKey && (e |= 1), l.ctrlKey && (e |= 2), l.altKey && (e |= 4), l.metaKey && (e |= 8), e;
}
function y(l, e) {
  if (l.length === 0) return 0;
  const t = [...l].sort((i, n) => i - n), s = Math.round((t.length - 1) * Math.min(Math.max(e, 0), 1));
  return t[s] ?? 0;
}
class M {
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
class w {
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
        <span class="repl-prompt">Patinae&gt;</span>
        <input class="repl-input" type="text" placeholder="Type a command..." spellcheck="false" autocomplete="off" />
      </div>
    `, this.output = e.querySelector(".repl-output"), this.input = e.querySelector(".repl-input"), this.input.addEventListener("keydown", (s) => this.onKey(s));
  }
  async onKey(e) {
    if (e.key === "Enter") {
      const t = this.input.value.trim();
      if (!t) return;
      this.history.push(t), this.historyIdx = this.history.length, this.input.value = "", this.appendLine(`Patinae> ${t}`, "cmd");
      const s = await this.viewer.executeAsync(t);
      for (const i of s.messages)
        this.appendOutputMessage(i);
    } else e.key === "ArrowUp" ? (e.preventDefault(), this.historyIdx > 0 && (this.historyIdx--, this.input.value = this.history[this.historyIdx])) : e.key === "ArrowDown" && (e.preventDefault(), this.historyIdx < this.history.length - 1 ? (this.historyIdx++, this.input.value = this.history[this.historyIdx]) : (this.historyIdx = this.history.length, this.input.value = ""));
  }
  appendOutputMessage(e) {
    if (e.level === "clear") {
      this.clearOutput();
      return;
    }
    this.appendLine(e.text, e.level);
  }
  appendLine(e, t) {
    const s = document.createElement("div");
    s.className = `repl-line repl-${t}`, s.textContent = e, this.output.appendChild(s), this.output.scrollTop = this.output.scrollHeight;
  }
  clearOutput() {
    this.output.replaceChildren();
  }
  update() {
  }
  destroy() {
    this.container.innerHTML = "";
  }
}
const F = ["lines", "sticks", "cartoon", "spheres", "surface"];
class $ {
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
        const o = document.createElement("span");
        o.className = "obj-reps";
        for (const h of F) {
          const d = document.createElement("button");
          d.className = "rep-btn", d.textContent = h.charAt(0).toUpperCase(), d.title = h, d.addEventListener("click", () => {
            this.viewer.execute(`show ${h}, ${s}`);
          }), o.appendChild(d);
        }
        n.appendChild(o);
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
        const o = document.createElement("button");
        o.className = "rep-btn selection-delete", o.textContent = "X", o.title = "Delete selection", o.addEventListener("click", () => {
          this.viewer.execute(`deselect ${i.name}`);
        }), n.appendChild(r), n.appendChild(c), n.appendChild(o), this.list.appendChild(n);
      }
    }
  }
  destroy() {
    this.container.innerHTML = "";
  }
}
class P {
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
          const o = `${t.object_name} and chain ${t.chain_id} and resi ${r.resv}`;
          this.viewer.execute(`select sele, ${o}`);
        }), n.appendChild(c), r.resv % 10 === 0) {
          const o = document.createElement("span");
          o.className = "seq-spacer", o.textContent = " ", n.appendChild(o);
        }
      }
      s.appendChild(n), this.content.appendChild(s);
    }
  }
  destroy() {
    this.container.innerHTML = "";
  }
}
class W {
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
class O {
  constructor(e, t = {}) {
    a(this, "core");
    a(this, "events", new M());
    a(this, "panels", /* @__PURE__ */ new Map());
    a(this, "options");
    this.options = t, this.core = new S(e), t.defer && this.core.setDeferred(!0, t.revealDuration ?? 150);
  }
  async init() {
    await this.core.init({
      picking: this.options.picking ?? !1,
      selectionOverlay: this.options.selectionOverlay,
      memoryProfile: this.options.memoryProfile
    }), this.options.picking && this.core.requireWasm("set_picking_enabled", (t) => t.set_picking_enabled(!0)), this.core.onOutput = (t) => this.emitRendererOutput(t), this.core.onPick = (t) => {
      const s = t ?? {
        object_name: null,
        atom_index: null,
        chain: null,
        residue: null,
        expression: null
      };
      this.events.emit("atom-picked", s), this.refreshPanels();
    };
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
      i.className = `patinae-panel patinae-panel-${t.name}`, t.collapsed && i.classList.add("collapsed"), s.appendChild(i);
      let n;
      switch (t.name) {
        case "repl":
          n = new w(i, this);
          break;
        case "objects":
          n = new $(i, this);
          break;
        case "sequence":
          n = new P(i, this);
          break;
        case "movie":
          n = new W(i, this);
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
    this.core.callWasm("set_picking_enabled", (t) => t.set_picking_enabled(e));
  }
  /** Enable or disable the visible selection / hover overlay at runtime. */
  setSelectionOverlay(e) {
    this.core.callWasm(
      "set_selection_overlay_enabled",
      (t) => t.set_selection_overlay_enabled(e)
    );
  }
  // ---------------------------------------------------------------------------
  // Commands
  // ---------------------------------------------------------------------------
  execute(e) {
    const t = b(e);
    if (t)
      return this.executeAsync(e), { messages: [{ level: "info", text: t.kind === "fetch" ? ` Fetching ${t.code}...` : ` Loading ${t.url}...` }] };
    const s = this.core.requireWasm("execute", (i) => i.execute(e));
    return this.emitCommandMessages(s.messages), this.refreshPanels(), s;
  }
  async executeAsync(e) {
    const t = b(e);
    if (!t) return this.execute(e);
    try {
      if (t.kind === "fetch") {
        const s = I(t.code, t.format);
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
    this.core.requireWasm("load_data", (i) => i.load_data(e, t, s)), this.refreshPanels();
  }
  async loadUrl(e, t) {
    var h;
    const s = await fetch(e);
    if (!s.ok) throw new Error(`Fetch failed: ${s.status} ${s.statusText}`);
    const i = new Uint8Array(await s.arrayBuffer());
    let r = new URL(e, location.href).pathname.split("/").pop() ?? "structure";
    r.toLowerCase().endsWith(".gz") && (r = r.slice(0, -3));
    const c = (t == null ? void 0 : t.name) ?? r.replace(/\.[^.]+$/, ""), o = (t == null ? void 0 : t.format) ?? ((h = r.split(".").pop()) == null ? void 0 : h.toLowerCase()) ?? "pdb";
    this.loadData(i, c, o);
  }
  // ---------------------------------------------------------------------------
  // Queries
  // ---------------------------------------------------------------------------
  getObjectNames() {
    return this.core.queryWasm(
      "get_object_names",
      [],
      (e) => e.get_object_names()
    );
  }
  getObjectInfo(e) {
    return this.core.queryWasm(
      "get_object_info",
      null,
      (t) => t.get_object_info(e)
    );
  }
  getSequenceData() {
    return this.core.queryWasm(
      "get_sequence_data",
      [],
      (e) => e.get_sequence_data()
    );
  }
  getMovieState() {
    return this.core.queryWasm(
      "get_movie_state",
      { frame_count: 0, current_frame: 0, is_playing: !1, rock_enabled: !1 },
      (e) => e.get_movie_state()
    );
  }
  getSelectionList() {
    return this.core.queryWasm(
      "get_selection_list",
      [],
      (e) => e.get_selection_list()
    );
  }
  getPerformanceSnapshot() {
    return this.core.getPerformanceSnapshot();
  }
  resetPerformanceStats() {
    this.core.resetPerformanceStats();
  }
  countAtoms(e = "all") {
    return this.core.requireWasm("count_atoms", (t) => t.count_atoms(e));
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
  emitCommandMessages(e) {
    for (const t of e)
      this.events.emit("command-output", t);
  }
  emitRendererOutput(e) {
    this.events.emit("command-output", e);
    const t = this.panels.get("repl");
    t instanceof w && t.appendOutputMessage(e);
  }
  destroy() {
    for (const e of this.panels.values())
      e.destroy();
    this.panels.clear(), this.core.destroy();
  }
}
function q(l) {
  const e = [], t = {};
  for (const s of l.split(",")) {
    const i = s.trim();
    if (!i) continue;
    const n = i.indexOf("=");
    n !== -1 ? t[i.slice(0, n).trim().toLowerCase()] = i.slice(n + 1).trim() : e.push(i);
  }
  return { positional: e, named: t };
}
function b(l) {
  const e = l.match(/^\s*(\w+)\s+(.*)/s);
  if (!e) return null;
  const t = e[1].toLowerCase(), { positional: s, named: i } = q(e[2]);
  if (t === "fetch" && s.length >= 1) {
    const n = s[0], r = s[1] ?? i.name ?? n, c = s[2] ?? i.type ?? "bcif", o = T(c);
    return { kind: "fetch", code: n, name: r, format: o };
  }
  if (t === "load" && s.length >= 1) {
    const n = s[0];
    if (!/^https?:\/\//i.test(n)) return null;
    const r = s[1] ?? i.object ?? i.name ?? j(n), c = i.format ?? N(n);
    return { kind: "load", url: n, name: r, format: c };
  }
  return null;
}
function T(l) {
  switch (l.toLowerCase()) {
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
const A = "https://models.rcsb.org", D = "https://files.rcsb.org/download";
function I(l, e) {
  const t = l.toLowerCase();
  return e === "bcif" ? `${A}/${t}.bcif.gz` : `${D}/${t}.${e}.gz`;
}
function j(l) {
  let e = new URL(l, location.href).pathname.split("/").pop() ?? "structure";
  return e.toLowerCase().endsWith(".gz") && (e = e.slice(0, -3)), e.replace(/\.[^.]+$/, "");
}
function N(l) {
  var t;
  let e = new URL(l, location.href).pathname.split("/").pop() ?? "";
  return e.toLowerCase().endsWith(".gz") && (e = e.slice(0, -3)), ((t = e.split(".").pop()) == null ? void 0 : t.toLowerCase()) ?? "pdb";
}
function H(l = "patinae-viewer") {
  customElements.get(l) || customElements.define(
    l,
    class extends HTMLElement {
      constructor() {
        super(...arguments);
        a(this, "viewer", null);
      }
      async connectedCallback() {
        var v;
        const t = this.attachShadow({ mode: "open" }), s = document.createElement("div");
        s.style.width = "100%", s.style.height = "100%", s.style.position = "relative", t.appendChild(s);
        const i = this.getAttribute("panels"), n = i ? i.split(",").map((u) => u.trim()) : [], r = this.getAttribute("src"), c = this.getAttribute("command"), o = this.getAttribute("defer"), h = o !== null ? o !== "false" : !!(r || c), d = this.getAttribute("selection-overlay"), g = d === null ? void 0 : d !== "false", m = (v = this.getAttribute("memory-profile")) == null ? void 0 : v.trim(), _ = m === "performance" || m === "balanced" || m === "lite" || /^manual:[1-9]\d*$/.test(m ?? "") || m === "auto" ? m : void 0;
        if (this.viewer = new O(s, {
          panels: n,
          defer: h,
          selectionOverlay: g,
          memoryProfile: _
        }), await this.viewer.init(), r) {
          const u = r.split(/\s+/).filter(Boolean);
          for (const p of u)
            await this.viewer.loadUrl(p);
        }
        if (c)
          for (const u of c.split(";")) {
            const p = u.trim();
            p && await this.viewer.executeAsync(p);
          }
        h && await this.viewer.show();
      }
      disconnectedCallback() {
        var t;
        (t = this.viewer) == null || t.destroy(), this.viewer = null;
      }
    }
  );
}
export {
  O as PatinaeViewer,
  H as registerElement
};
