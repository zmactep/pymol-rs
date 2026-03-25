/** Types matching the WASM bridge serialization format. */

export interface CommandOutput {
  messages: OutputMessage[];
}

export interface OutputMessage {
  level: "info" | "warning" | "error";
  text: string;
}

export interface ObjectInfo {
  name: string;
  object_type: "molecule" | "map";
  atom_count: number;
  enabled: boolean;
}

export interface SelectionInfo {
  name: string;
  expression: string;
  visible: boolean;
}

export interface SequenceChain {
  object_name: string;
  chain_id: string;
  residues: SequenceResidue[];
}

export interface SequenceResidue {
  resn: string;
  resv: number;
  one_letter: string;
}

export interface MovieState {
  frame_count: number;
  current_frame: number;
  is_playing: boolean;
}

export interface LabelInfo {
  x: number;
  y: number;
  text: string;
  kind: "atom" | "measurement";
}

export type PanelSlot = "top" | "right" | "bottom";

export interface PanelPlacement {
  name: PanelName;
  slot: PanelSlot;
  collapsed?: boolean;
}

/** Information about a picked atom returned by `pick_at_screen`. */
export interface PickHitInfo {
  /** Name of the picked object. */
  object_name: string;
  /** Zero-based atom index within the object, or `null` for non-atom hits. */
  atom_index: number | null;
  /** Chain identifier of the hit atom, or `null`. */
  chain: string | null;
  /** Residue sequence number of the hit atom, or `null`. */
  residue: number | null;
  /** PyMOL selection expression (depends on `mouse_selection_mode`), or `null`. */
  expression: string | null;
}

export interface ViewerOptions {
  /** Simple list — all panels go into a "sidebar" element (legacy). */
  panels?: PanelName[];
  /** Desktop-style layout — panels placed in named slots around the viewport. */
  layout?: PanelPlacement[];
  /** Container elements for each slot (required when using `layout`). */
  slots?: Partial<Record<PanelSlot, HTMLElement>>;
  /** When true, the viewer starts hidden and is revealed only when `show()` is called. */
  defer?: boolean;
  /** Duration of the fade-in transition in milliseconds (default: 150). */
  revealDuration?: number;
  /**
   * When true, left-click picks atoms and updates the `sele` named selection.
   * Fires `atom-picked` events with hit details. Default: false.
   */
  picking?: boolean;
}

export type PanelName = "repl" | "objects" | "sequence" | "movie";

/** Modifier bitmask constants (must match event.rs). */
export const MOD_SHIFT = 1;
export const MOD_CTRL = 2;
export const MOD_ALT = 4;
export const MOD_META = 8;
