/**
 * Public API facade — wraps ViewerCore + panels into a clean interface.
 */

import { ViewerCore } from "./viewer.js";
import { ViewerEvents } from "./events.js";
import type {
  CommandOutput,
  ObjectInfo,
  PickHitInfo,
  SelectionInfo,
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
    if (options.defer) {
      this.core.setDeferred(true, options.revealDuration ?? 150);
    }
  }

  async init(): Promise<void> {
    const wasm = await this.core.init();

    // Enable picking if requested.
    if (this.options.picking) {
      wasm.set_picking_enabled(true);
    }

    // Forward pick results as typed events.
    this.core.onPick = (hit) => {
      const h = (hit as PickHitInfo | null) ?? {
        object_name: null,
        atom_index: null,
        chain: null,
        residue: null,
        expression: null,
      };
      this.events.emit("atom-picked", h);
      this.refreshPanels();
    };

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
  // Deferred display
  // ---------------------------------------------------------------------------

  get isDeferred(): boolean {
    return this.core.isDeferred;
  }

  async show(): Promise<void> {
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
  setPicking(enabled: boolean): void {
    this.core.wasmViewer?.set_picking_enabled(enabled);
  }

  // ---------------------------------------------------------------------------
  // Commands
  // ---------------------------------------------------------------------------

  execute(command: string): CommandOutput {
    const async_cmd = parseAsyncCommand(command);
    if (async_cmd) {
      // Fire-and-forget; callers who need to await should use executeAsync()
      this.executeAsync(command);
      const text = async_cmd.kind === "fetch"
        ? ` Fetching ${async_cmd.code}...`
        : ` Loading ${async_cmd.url}...`;
      return { messages: [{ level: "info", text }] };
    }
    const wasm = this.core.wasmViewer;
    if (!wasm) return { messages: [] };
    const result = wasm.execute(command) as CommandOutput;
    this.events.emit("command-output", result.messages[0] ?? { level: "info", text: "" });
    this.refreshPanels();
    return result;
  }

  async executeAsync(command: string): Promise<CommandOutput> {
    const cmd = parseAsyncCommand(command);
    if (!cmd) return this.execute(command);

    try {
      if (cmd.kind === "fetch") {
        const url = buildRcsbUrl(cmd.code, cmd.format);
        await this.loadUrl(url, { name: cmd.name, format: cmd.format });
        const msg = { level: "info" as const, text: ` Fetched ${cmd.code} as "${cmd.name}"` };
        this.events.emit("command-output", msg);
        return { messages: [msg] };
      } else {
        await this.loadUrl(cmd.url, { name: cmd.name, format: cmd.format });
        const msg = { level: "info" as const, text: ` Loaded "${cmd.name}" from URL` };
        this.events.emit("command-output", msg);
        return { messages: [msg] };
      }
    } catch (e) {
      const msg = { level: "error" as const, text: ` ${e}` };
      this.events.emit("command-output", msg);
      return { messages: [msg] };
    }
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
    let fileName = urlPath.split("/").pop() ?? "structure";
    // Strip .gz suffix — gzip is handled transparently by the WASM layer
    if (fileName.toLowerCase().endsWith(".gz")) {
      fileName = fileName.slice(0, -3);
    }
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

  getSelectionList(): SelectionInfo[] {
    const wasm = this.core.wasmViewer;
    if (!wasm) return [];
    return wasm.get_selection_list() as SelectionInfo[];
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

// ---------------------------------------------------------------------------
// Command interception helpers
// ---------------------------------------------------------------------------

type AsyncCommand =
  | { kind: "fetch"; code: string; name: string; format: string }
  | { kind: "load"; url: string; name: string; format: string };

/** Parse a PyMOL command string into positional and named args. */
function parsePymolArgs(argsStr: string): { positional: string[]; named: Record<string, string> } {
  const positional: string[] = [];
  const named: Record<string, string> = {};
  for (const part of argsStr.split(",")) {
    const trimmed = part.trim();
    if (!trimmed) continue;
    const eq = trimmed.indexOf("=");
    if (eq !== -1) {
      named[trimmed.slice(0, eq).trim().toLowerCase()] = trimmed.slice(eq + 1).trim();
    } else {
      positional.push(trimmed);
    }
  }
  return { positional, named };
}

/** Detect fetch/load commands that should be handled via JS fetch API. */
function parseAsyncCommand(command: string): AsyncCommand | null {
  const match = command.match(/^\s*(\w+)\s+(.*)/s);
  if (!match) return null;
  const verb = match[1].toLowerCase();
  const { positional, named } = parsePymolArgs(match[2]);

  if (verb === "fetch" && positional.length >= 1) {
    const code = positional[0];
    const name = positional[1] ?? named["name"] ?? code;
    const typeStr = positional[2] ?? named["type"] ?? "bcif";
    const format = normalizeFormat(typeStr);
    return { kind: "fetch", code, name, format };
  }

  if (verb === "load" && positional.length >= 1) {
    const filename = positional[0];
    if (!/^https?:\/\//i.test(filename)) return null; // local path — let WASM handle it
    const name = positional[1] ?? named["object"] ?? named["name"] ?? urlToName(filename);
    const format = named["format"] ?? urlToFormat(filename);
    return { kind: "load", url: filename, name, format };
  }

  return null;
}

function normalizeFormat(s: string): string {
  switch (s.toLowerCase()) {
    case "pdb": return "pdb";
    case "cif": case "mmcif": return "cif";
    case "bcif": case "binarycif": return "bcif";
    default: return "bcif";
  }
}

const RCSB_MODELS = "https://models.rcsb.org";
const RCSB_FILES = "https://files.rcsb.org/download";

function buildRcsbUrl(pdbId: string, format: string): string {
  const id = pdbId.toLowerCase();
  if (format === "bcif") return `${RCSB_MODELS}/${id}.bcif.gz`;
  return `${RCSB_FILES}/${id}.${format}.gz`;
}

function urlToName(url: string): string {
  let fileName = new URL(url, location.href).pathname.split("/").pop() ?? "structure";
  if (fileName.toLowerCase().endsWith(".gz")) fileName = fileName.slice(0, -3);
  return fileName.replace(/\.[^.]+$/, "");
}

function urlToFormat(url: string): string {
  let fileName = new URL(url, location.href).pathname.split("/").pop() ?? "";
  if (fileName.toLowerCase().endsWith(".gz")) fileName = fileName.slice(0, -3);
  return fileName.split(".").pop()?.toLowerCase() ?? "pdb";
}
