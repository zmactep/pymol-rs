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

export interface ViewerOptions {
  /** Simple list — all panels go into a "sidebar" element (legacy). */
  panels?: PanelName[];
  /** Desktop-style layout — panels placed in named slots around the viewport. */
  layout?: PanelPlacement[];
  /** Container elements for each slot (required when using `layout`). */
  slots?: Partial<Record<PanelSlot, HTMLElement>>;
}

export type PanelName = "repl" | "objects" | "sequence" | "movie";

/** Modifier bitmask constants (must match event.rs). */
export const MOD_SHIFT = 1;
export const MOD_CTRL = 2;
export const MOD_ALT = 4;
export const MOD_META = 8;
