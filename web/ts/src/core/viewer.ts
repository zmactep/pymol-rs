/**
 * ViewerCore — canvas lifecycle, WASM initialization, rAF loop, and event wiring.
 */

import type { WebViewer } from "../../../pkg/pymol_web.js";
import { MOD_SHIFT, MOD_CTRL, MOD_ALT, MOD_META } from "./types.js";

export class ViewerCore {
  private wasm: WebViewer | null = null;
  private canvas: HTMLCanvasElement;
  private animFrameId = 0;
  private resizeObserver: ResizeObserver | null = null;

  constructor(private container: HTMLElement) {
    this.canvas = document.createElement("canvas");
    this.canvas.id = "pymol-rs-canvas-" + Math.random().toString(36).slice(2, 8);
    this.canvas.style.width = "100%";
    this.canvas.style.height = "100%";
    this.canvas.style.display = "block";
    this.canvas.tabIndex = 0; // focusable for keyboard events
    container.appendChild(this.canvas);
  }

  async init(): Promise<WebViewer> {
    // Dynamically import the WASM module
    const wasmModule = await import("../../../pkg/pymol_web.js");
    await wasmModule.default(); // init wasm

    // Match canvas pixel size to layout size
    this.syncCanvasSize();

    const viewer = await wasmModule.WebViewer.create(this.canvas.id);
    this.wasm = viewer;

    this.bindEvents();
    this.startLoop();

    return viewer;
  }

  get wasmViewer(): WebViewer | null {
    return this.wasm;
  }

  destroy(): void {
    if (this.animFrameId) {
      cancelAnimationFrame(this.animFrameId);
      this.animFrameId = 0;
    }
    this.resizeObserver?.disconnect();
    this.canvas.remove();
    this.wasm = null;
  }

  // ---------------------------------------------------------------------------
  // Render loop
  // ---------------------------------------------------------------------------

  private startLoop(): void {
    const loop = () => {
      this.animFrameId = requestAnimationFrame(loop);
      if (!this.wasm) return;

      this.wasm.process_input();

      if (this.wasm.needs_redraw()) {
        this.wasm.render_frame();
      }
    };
    this.animFrameId = requestAnimationFrame(loop);
  }

  // ---------------------------------------------------------------------------
  // Event binding
  // ---------------------------------------------------------------------------

  private bindEvents(): void {
    const c = this.canvas;

    c.addEventListener("mousedown", (e) => {
      e.preventDefault();
      c.focus();
      this.wasm?.on_mouse_down(e.offsetX, e.offsetY, e.button, modBits(e));
    });

    c.addEventListener("mousemove", (e) => {
      this.wasm?.on_mouse_move(e.offsetX, e.offsetY, modBits(e));
    });

    c.addEventListener("mouseup", (e) => {
      this.wasm?.on_mouse_up(e.offsetX, e.offsetY, e.button);
    });

    c.addEventListener("wheel", (e) => {
      e.preventDefault();
      this.wasm?.on_wheel(e.deltaY, modBits(e));
    }, { passive: false });

    c.addEventListener("contextmenu", (e) => e.preventDefault());

    // Resize
    this.resizeObserver = new ResizeObserver(() => {
      this.syncCanvasSize();
      if (this.wasm) {
        this.wasm.resize(this.canvas.width, this.canvas.height);
      }
    });
    this.resizeObserver.observe(this.container);
  }

  private syncCanvasSize(): void {
    const dpr = window.devicePixelRatio || 1;
    const rect = this.canvas.getBoundingClientRect();
    this.canvas.width = Math.round(rect.width * dpr);
    this.canvas.height = Math.round(rect.height * dpr);
  }
}

/** Convert MouseEvent/WheelEvent modifiers to bitmask. */
function modBits(e: MouseEvent | WheelEvent): number {
  let bits = 0;
  if (e.shiftKey) bits |= MOD_SHIFT;
  if (e.ctrlKey) bits |= MOD_CTRL;
  if (e.altKey) bits |= MOD_ALT;
  if (e.metaKey) bits |= MOD_META;
  return bits;
}
