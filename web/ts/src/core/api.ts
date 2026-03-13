/**
 * Public API facade — wraps ViewerCore + panels into a clean interface.
 */

import { ViewerCore } from "./viewer.js";
import { ViewerEvents } from "./events.js";
import type {
  CommandOutput,
  ObjectInfo,
  SequenceChain,
  MovieState,
  ViewerOptions,
  PanelName,
  PanelPlacement,
} from "./types.js";
import type { ViewerEventType, ViewerEventMap } from "./events.js";
import { ReplPanel } from "../panels/repl.js";
import { ObjectListPanel } from "../panels/object-list.js";
import { SequencePanel } from "../panels/sequence.js";
import { MoviePanel } from "../panels/movie.js";

type PanelInstance = ReplPanel | ObjectListPanel | SequencePanel | MoviePanel;

export class PyMolRSViewer {
  private core: ViewerCore;
  private events = new ViewerEvents();
  private panels = new Map<PanelName, PanelInstance>();
  private options: ViewerOptions;

  constructor(container: HTMLElement, options: ViewerOptions = {}) {
    this.options = options;
    this.core = new ViewerCore(container);
  }

  async init(): Promise<void> {
    const wasm = await this.core.init();

    // Build the list of panels to mount
    const placements: PanelPlacement[] = [];

    if (this.options.layout) {
      placements.push(...this.options.layout);
    } else if (this.options.panels) {
      // Legacy mode — all panels go into a sidebar
      for (const name of this.options.panels) {
        placements.push({ name, slot: "right" });
      }
    }

    for (const placement of placements) {
      // Find the target container for this panel
      let target: HTMLElement | null | undefined;
      if (this.options.slots) {
        target = this.options.slots[placement.slot];
      }
      if (!target) {
        target = document.getElementById("sidebar");
      }
      if (!target) continue;

      const panelEl = document.createElement("div");
      panelEl.className = `pymol-panel pymol-panel-${placement.name}`;
      if (placement.collapsed) {
        panelEl.classList.add("collapsed");
      }
      target.appendChild(panelEl);

      let panel: PanelInstance;
      switch (placement.name) {
        case "repl":
          panel = new ReplPanel(panelEl, this);
          break;
        case "objects":
          panel = new ObjectListPanel(panelEl, this);
          break;
        case "sequence":
          panel = new SequencePanel(panelEl, this);
          break;
        case "movie":
          panel = new MoviePanel(panelEl, this);
          break;
      }
      this.panels.set(placement.name, panel);
    }

    this.events.emit("ready", {});
  }

  // ---------------------------------------------------------------------------
  // Commands
  // ---------------------------------------------------------------------------

  execute(command: string): CommandOutput {
    const wasm = this.core.wasmViewer;
    if (!wasm) return { messages: [] };
    const result = wasm.execute(command) as CommandOutput;
    this.events.emit("command-output", result.messages[0] ?? { level: "info", text: "" });
    this.refreshPanels();
    return result;
  }

  loadData(data: Uint8Array, name: string, format: string): void {
    const wasm = this.core.wasmViewer;
    if (!wasm) return;
    wasm.load_data(data, name, format);
    this.refreshPanels();
  }

  async loadUrl(url: string, options?: { name?: string; format?: string }): Promise<void> {
    const resp = await fetch(url);
    if (!resp.ok) throw new Error(`Fetch failed: ${resp.status} ${resp.statusText}`);
    const data = new Uint8Array(await resp.arrayBuffer());

    const urlPath = new URL(url, location.href).pathname;
    const fileName = urlPath.split("/").pop() ?? "structure";
    const name = options?.name ?? fileName.replace(/\.[^.]+$/, "");
    const format = options?.format ?? fileName.split(".").pop()?.toLowerCase() ?? "pdb";

    this.loadData(data, name, format);
  }

  // ---------------------------------------------------------------------------
  // Queries
  // ---------------------------------------------------------------------------

  getObjectNames(): string[] {
    const wasm = this.core.wasmViewer;
    if (!wasm) return [];
    return wasm.get_object_names() as string[];
  }

  getObjectInfo(name: string): ObjectInfo | null {
    const wasm = this.core.wasmViewer;
    if (!wasm) return null;
    return wasm.get_object_info(name) as ObjectInfo | null;
  }

  getSequenceData(): SequenceChain[] {
    const wasm = this.core.wasmViewer;
    if (!wasm) return [];
    return wasm.get_sequence_data() as SequenceChain[];
  }

  getMovieState(): MovieState {
    const wasm = this.core.wasmViewer;
    if (!wasm) return { frame_count: 0, current_frame: 0, is_playing: false };
    return wasm.get_movie_state() as MovieState;
  }

  countAtoms(selection = "all"): number {
    const wasm = this.core.wasmViewer;
    if (!wasm) return 0;
    return wasm.count_atoms(selection);
  }

  // ---------------------------------------------------------------------------
  // Events
  // ---------------------------------------------------------------------------

  on<K extends ViewerEventType>(
    event: K,
    callback: (data: ViewerEventMap[K]) => void
  ): void {
    this.events.on(event, callback);
  }

  off<K extends ViewerEventType>(
    event: K,
    callback: (data: ViewerEventMap[K]) => void
  ): void {
    this.events.off(event, callback);
  }

  // ---------------------------------------------------------------------------
  // Panel management
  // ---------------------------------------------------------------------------

  private refreshPanels(): void {
    for (const panel of this.panels.values()) {
      panel.update();
    }
    this.events.emit("objects-changed", { names: this.getObjectNames() });
  }

  destroy(): void {
    for (const panel of this.panels.values()) {
      panel.destroy();
    }
    this.panels.clear();
    this.core.destroy();
  }
}
