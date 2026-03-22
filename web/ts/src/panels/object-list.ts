/**
 * Object List panel — shows loaded molecules with visibility toggles
 * and representation buttons.
 */

import type { PyMolRSViewer } from "../core/api.js";

const REPS = ["lines", "sticks", "cartoon", "spheres", "surface"] as const;

export class ObjectListPanel {
  private container: HTMLElement;
  private viewer: PyMolRSViewer;
  private list: HTMLElement;

  constructor(container: HTMLElement, viewer: PyMolRSViewer) {
    this.container = container;
    this.viewer = viewer;

    container.innerHTML = `
      <div class="panel-header">Objects</div>
      <div class="object-list"></div>
    `;

    this.list = container.querySelector(".object-list")!;
  }

  update(): void {
    const names = this.viewer.getObjectNames();
    this.list.innerHTML = "";

    for (const name of names) {
      const info = this.viewer.getObjectInfo(name);
      if (!info) continue;

      const row = document.createElement("div");
      row.className = "object-row";

      // Visibility toggle
      const vis = document.createElement("button");
      vis.className = `obj-vis ${info.enabled ? "enabled" : "disabled"}`;
      vis.textContent = info.enabled ? "V" : "-";
      vis.title = info.enabled ? "Hide" : "Show";
      vis.addEventListener("click", () => {
        this.viewer.execute(info.enabled ? `disable ${name}` : `enable ${name}`);
      });

      // Name
      const label = document.createElement("span");
      label.className = "obj-name";
      label.textContent = info.object_type === "map"
        ? `${name} (map)`
        : `${name} (${info.atom_count})`;

      row.appendChild(vis);
      row.appendChild(label);

      // Representation buttons (molecule-only)
      if (info.object_type !== "map") {
        const reps = document.createElement("span");
        reps.className = "obj-reps";
        for (const rep of REPS) {
          const btn = document.createElement("button");
          btn.className = "rep-btn";
          btn.textContent = rep.charAt(0).toUpperCase();
          btn.title = rep;
          btn.addEventListener("click", () => {
            this.viewer.execute(`show ${rep}, ${name}`);
          });
          reps.appendChild(btn);
        }
        row.appendChild(reps);
      }
      this.list.appendChild(row);
    }

    // Named selections
    const selections = this.viewer.getSelectionList();
    if (selections.length > 0) {
      const sep = document.createElement("hr");
      sep.className = "object-list-separator";
      this.list.appendChild(sep);

      for (const sel of selections) {
        const row = document.createElement("div");
        row.className = "object-row selection-row";

        // Visibility toggle
        const vis = document.createElement("button");
        vis.className = `obj-vis selection-vis ${sel.visible ? "enabled" : "disabled"}`;
        vis.textContent = sel.visible ? "V" : "-";
        vis.title = sel.visible ? "Hide indicators" : "Show indicators";
        vis.addEventListener("click", () => {
          this.viewer.execute(`toggle ${sel.name}`);
        });

        // Name in parentheses (PyMOL convention)
        const label = document.createElement("span");
        label.className = "obj-name selection-name";
        label.textContent = `(${sel.name})`;
        label.title = sel.expression;

        // Delete button
        const del = document.createElement("button");
        del.className = "rep-btn selection-delete";
        del.textContent = "X";
        del.title = "Delete selection";
        del.addEventListener("click", () => {
          this.viewer.execute(`deselect ${sel.name}`);
        });

        row.appendChild(vis);
        row.appendChild(label);
        row.appendChild(del);
        this.list.appendChild(row);
      }
    }
  }

  destroy(): void {
    this.container.innerHTML = "";
  }
}
