/**
 * Movie panel — animation timeline and playback controls.
 */

import type { PyMolRSViewer } from "../core/api.js";

export class MoviePanel {
  private container: HTMLElement;
  private viewer: PyMolRSViewer;
  private frameLabel: HTMLElement | null = null;
  private slider: HTMLInputElement | null = null;
  private playBtn: HTMLButtonElement | null = null;

  constructor(container: HTMLElement, viewer: PyMolRSViewer) {
    this.container = container;
    this.viewer = viewer;

    container.innerHTML = `
      <div class="panel-header">Movie</div>
      <div class="movie-controls">
        <div class="movie-buttons">
          <button class="movie-btn" data-action="beginning" title="Beginning">|&lt;</button>
          <button class="movie-btn" data-action="backward" title="Step back">&lt;</button>
          <button class="movie-btn movie-play" data-action="play" title="Play/Pause">&#9654;</button>
          <button class="movie-btn" data-action="forward" title="Step forward">&gt;</button>
          <button class="movie-btn" data-action="end" title="End">&gt;|</button>
        </div>
        <div class="movie-timeline">
          <input type="range" class="movie-slider" min="1" max="1" value="1" />
          <span class="movie-frame-label">1 / 1</span>
        </div>
      </div>
    `;

    this.frameLabel = container.querySelector(".movie-frame-label")!;
    this.slider = container.querySelector(".movie-slider")!;
    this.playBtn = container.querySelector(".movie-play")!;

    container.querySelectorAll(".movie-btn").forEach((btn) => {
      btn.addEventListener("click", () => {
        const action = (btn as HTMLElement).dataset.action;
        switch (action) {
          case "play":
            this.viewer.execute("mplay");
            break;
          case "beginning":
            this.viewer.execute("frame 1");
            break;
          case "end": {
            const state = this.viewer.getMovieState();
            this.viewer.execute(`frame ${state.frame_count}`);
            break;
          }
          case "forward":
            this.viewer.execute("forward");
            break;
          case "backward":
            this.viewer.execute("backward");
            break;
        }
      });
    });

    this.slider.addEventListener("input", () => {
      const frame = parseInt(this.slider!.value, 10);
      this.viewer.execute(`frame ${frame}`);
    });
  }

  update(): void {
    const state = this.viewer.getMovieState();
    if (this.slider) {
      this.slider.max = String(Math.max(1, state.frame_count));
      this.slider.value = String(state.current_frame + 1);
    }
    if (this.frameLabel) {
      this.frameLabel.textContent = `${state.current_frame + 1} / ${Math.max(1, state.frame_count)}`;
    }
    if (this.playBtn) {
      this.playBtn.textContent = state.is_playing ? "\u23F8" : "\u25B6";
    }
  }

  destroy(): void {
    this.container.innerHTML = "";
  }
}
