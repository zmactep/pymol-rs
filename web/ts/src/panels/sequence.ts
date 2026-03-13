/**
 * Sequence Viewer panel — displays residue sequences per chain,
 * with click-to-select support.
 */

import type { PyMolRSViewer } from "../core/api.js";

export class SequencePanel {
  private container: HTMLElement;
  private viewer: PyMolRSViewer;
  private content: HTMLElement;

  constructor(container: HTMLElement, viewer: PyMolRSViewer) {
    this.container = container;
    this.viewer = viewer;

    container.innerHTML = `
      <div class="panel-header">Sequence</div>
      <div class="sequence-content"></div>
    `;

    this.content = container.querySelector(".sequence-content")!;
  }

  update(): void {
    const chains = this.viewer.getSequenceData();
    this.content.innerHTML = "";

    for (const chain of chains) {
      const section = document.createElement("div");
      section.className = "seq-chain";

      const header = document.createElement("div");
      header.className = "seq-chain-header";
      header.textContent = `${chain.object_name} / ${chain.chain_id || "?"}`;
      section.appendChild(header);

      const residues = document.createElement("div");
      residues.className = "seq-residues";

      for (const res of chain.residues) {
        const span = document.createElement("span");
        span.className = "seq-residue";
        span.textContent = res.one_letter;
        span.title = `${res.resn} ${res.resv}`;
        span.addEventListener("click", () => {
          const sel = `${chain.object_name} and chain ${chain.chain_id} and resi ${res.resv}`;
          this.viewer.execute(`select sele, ${sel}`);
        });
        residues.appendChild(span);

        // Add space every 10 residues for readability
        if (res.resv % 10 === 0) {
          const spacer = document.createElement("span");
          spacer.className = "seq-spacer";
          spacer.textContent = " ";
          residues.appendChild(spacer);
        }
      }

      section.appendChild(residues);
      this.content.appendChild(section);
    }
  }

  destroy(): void {
    this.container.innerHTML = "";
  }
}
