/**
 * REPL panel — command input with history and scrollable output log.
 */

import type { PyMolRSViewer } from "../core/api.js";

export class ReplPanel {
  private container: HTMLElement;
  private viewer: PyMolRSViewer;
  private output: HTMLElement;
  private input: HTMLInputElement;
  private history: string[] = [];
  private historyIdx = -1;

  constructor(container: HTMLElement, viewer: PyMolRSViewer) {
    this.container = container;
    this.viewer = viewer;

    container.innerHTML = `
      <div class="repl-header">Command Line</div>
      <div class="repl-output"></div>
      <div class="repl-input-row">
        <span class="repl-prompt">PyMOL&gt;</span>
        <input class="repl-input" type="text" placeholder="Type a command..." spellcheck="false" autocomplete="off" />
      </div>
    `;

    this.output = container.querySelector(".repl-output")!;
    this.input = container.querySelector(".repl-input")!;

    this.input.addEventListener("keydown", (e) => this.onKey(e));
  }

  private onKey(e: KeyboardEvent): void {
    if (e.key === "Enter") {
      const cmd = this.input.value.trim();
      if (!cmd) return;

      this.history.push(cmd);
      this.historyIdx = this.history.length;
      this.appendLine(`PyMOL> ${cmd}`, "cmd");

      const result = this.viewer.execute(cmd);
      for (const msg of result.messages) {
        this.appendLine(msg.text, msg.level);
      }

      this.input.value = "";
    } else if (e.key === "ArrowUp") {
      e.preventDefault();
      if (this.historyIdx > 0) {
        this.historyIdx--;
        this.input.value = this.history[this.historyIdx];
      }
    } else if (e.key === "ArrowDown") {
      e.preventDefault();
      if (this.historyIdx < this.history.length - 1) {
        this.historyIdx++;
        this.input.value = this.history[this.historyIdx];
      } else {
        this.historyIdx = this.history.length;
        this.input.value = "";
      }
    }
  }

  private appendLine(text: string, level: string): void {
    const line = document.createElement("div");
    line.className = `repl-line repl-${level}`;
    line.textContent = text;
    this.output.appendChild(line);
    this.output.scrollTop = this.output.scrollHeight;
  }

  update(): void {
    // REPL doesn't need periodic updates
  }

  destroy(): void {
    this.container.innerHTML = "";
  }
}
