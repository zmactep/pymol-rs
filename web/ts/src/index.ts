/**
 * @pymol-rs/viewer — embeddable WebGPU molecular visualization.
 *
 * Usage:
 *   import { PyMolRSViewer } from "@pymol-rs/viewer";
 *   const viewer = new PyMolRSViewer(document.getElementById("container")!);
 *   await viewer.init();
 *   viewer.execute("load https://files.rcsb.org/download/1CRN.pdb");
 */

export { PyMolRSViewer } from "./core/api.js";
export type {
  CommandOutput,
  ObjectInfo,
  SequenceChain,
  SequenceResidue,
  MovieState,
  ViewerOptions,
  PanelName,
  PanelSlot,
  PanelPlacement,
} from "./core/types.js";
export type { ViewerEventType, ViewerEventMap } from "./core/events.js";

// Custom element registration
import { PyMolRSViewer } from "./core/api.js";

/**
 * Register `<pymol-rs-viewer>` as a custom HTML element.
 *
 * Attributes:
 *   src      — URL to load on mount
 *   panels   — comma-separated panel names (repl, objects, sequence, movie)
 *   command  — PyMOL command to run after loading
 */
export function registerElement(tagName = "pymol-rs-viewer"): void {
  if (customElements.get(tagName)) return;

  customElements.define(
    tagName,
    class extends HTMLElement {
      private viewer: PyMolRSViewer | null = null;

      async connectedCallback() {
        const shadow = this.attachShadow({ mode: "open" });

        const wrapper = document.createElement("div");
        wrapper.style.width = "100%";
        wrapper.style.height = "100%";
        wrapper.style.position = "relative";
        shadow.appendChild(wrapper);

        const panelAttr = this.getAttribute("panels");
        const panels = panelAttr
          ? (panelAttr.split(",").map((s) => s.trim()) as Array<"repl" | "objects" | "sequence" | "movie">)
          : [];

        this.viewer = new PyMolRSViewer(wrapper, { panels });
        await this.viewer.init();

        const src = this.getAttribute("src");
        if (src) {
          await this.viewer.loadUrl(src);
        }

        const cmd = this.getAttribute("command");
        if (cmd) {
          this.viewer.execute(cmd);
        }
      }

      disconnectedCallback() {
        this.viewer?.destroy();
        this.viewer = null;
      }
    }
  );
}
