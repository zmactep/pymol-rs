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
 *   src      — URL(s) to load on mount (whitespace-separated for multiple)
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

        const src = this.getAttribute("src");
        const cmd = this.getAttribute("command");
        const deferAttr = this.getAttribute("defer");
        const shouldDefer = deferAttr !== null
          ? deferAttr !== "false"
          : !!(src || cmd);

        this.viewer = new PyMolRSViewer(wrapper, { panels, defer: shouldDefer });
        await this.viewer.init();

        if (src) {
          const urls = src.split(/\s+/).filter(Boolean);
          for (const url of urls) {
            await this.viewer.loadUrl(url);
          }
        }

        if (cmd) {
          // Split on semicolons to support multiple commands;
          // use executeAsync so fetch/load complete before reveal.
          for (const c of cmd.split(";")) {
            const trimmed = c.trim();
            if (trimmed) await this.viewer.executeAsync(trimmed);
          }
        }

        if (shouldDefer) {
          await this.viewer.show();
        }
      }

      disconnectedCallback() {
        this.viewer?.destroy();
        this.viewer = null;
      }
    }
  );
}
