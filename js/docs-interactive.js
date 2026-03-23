// ── Utilities ───────────────────────────────────────────────────────────────

function flashElement(el, className = "flash", ms = 400) {
  el.classList.add(className);
  setTimeout(() => el.classList.remove(className), ms);
}

function execCmdString(viewer, cmdStr) {
  for (const cmd of cmdStr.split(";").map(s => s.trim()).filter(Boolean))
    viewer.execute(cmd);
}

// ── SPA page navigation ─────────────────────────────────────────────────────

/**
 * Initializes SPA-style documentation with hash-routed sections.
 *
 * Each section is a `.doc-page` div with `id` and optional data attributes:
 *   data-src="url"              — structure file(s) to load (whitespace-separated)
 *   data-command="cmd1; cmd2"   — semicolon-separated PyMOL commands
 *
 * Inline commands use `data-cmd` (semicolon-separated) with optional
 * `data-cmd-alt` for toggle behaviour.
 */
export function initDocPages(viewer) {
  // ── Inline command clicks ─────────────────────────────────────────────
  document.addEventListener("click", (e) => {
    const el = e.target.closest("[data-cmd]");
    if (!el) return;
    e.preventDefault();
    const alt = el.dataset.cmdAlt;
    const toggled = el.classList.contains("toggled");
    const cmdStr = (alt && toggled) ? alt : el.dataset.cmd;
    if (alt) el.classList.toggle("toggled");

    // Radio-group: un-toggle siblings in the same group
    const group = el.dataset.cmdGroup;
    if (group && !toggled) {
      document.querySelectorAll(`[data-cmd-group="${group}"]`).forEach(sib => {
        if (sib !== el) sib.classList.remove("toggled");
      });
    }

    execCmdString(viewer, cmdStr);
    flashElement(el);
  });

  // ── Page management ───────────────────────────────────────────────────
  const pages = Array.from(document.querySelectorAll(".doc-page"));
  if (pages.length === 0) return;

  const sidebar = document.querySelector(".docs-sidebar");
  let currentIdx = 0;
  let setupSeq = 0;

  function showPage(idx, pushState = true) {
    if (idx < 0 || idx >= pages.length) return;
    currentIdx = idx;
    const page = pages[idx];

    // Show only the active page
    pages.forEach(p => p.style.display = "none");
    page.style.display = "block";

    // Sidebar highlight
    if (sidebar) {
      sidebar.querySelectorAll("a").forEach((a, i) => {
        a.classList.toggle("active", i === idx);
      });
    }

    // URL hash
    if (pushState) history.pushState(null, "", "#" + page.id);

    // Viewer setup for this section
    setupViewer(page);

    // Prev / Next buttons
    page.querySelectorAll("[data-nav]").forEach(btn => {
      const dir = btn.dataset.nav;
      if (dir === "prev") {
        btn.style.visibility = idx > 0 ? "visible" : "hidden";
        if (idx > 0) btn.textContent = "\u2190 " + (pages[idx - 1].dataset.title || "Previous");
      } else if (dir === "next") {
        btn.style.visibility = idx < pages.length - 1 ? "visible" : "hidden";
        if (idx < pages.length - 1) btn.textContent = (pages[idx + 1].dataset.title || "Next") + " \u2192";
      }
    });

    // Page counter
    const counter = page.querySelector(".page-counter");
    if (counter) counter.textContent = `${idx + 1} of ${pages.length}`;

    // Scroll to top
    window.scrollTo({ top: 0 });
  }

  async function setupViewer(page) {
    const src = page.dataset.src;
    const command = page.dataset.command;
    if (!src && !command) return;

    const seq = ++setupSeq;
    viewer.core.setDeferred(true);
    viewer.execute("delete all");

    if (src) {
      for (const url of src.split(/\s+/).filter(Boolean))
        await viewer.loadUrl(url);
    }

    if (seq !== setupSeq) return;
    if (command) execCmdString(viewer, command);

    await viewer.show();
  }

  // ── Navigation wiring ─────────────────────────────────────────────────
  document.addEventListener("click", (e) => {
    const btn = e.target.closest("[data-nav]");
    if (!btn) return;
    e.preventDefault();
    if (btn.dataset.nav === "prev") showPage(currentIdx - 1);
    else if (btn.dataset.nav === "next") showPage(currentIdx + 1);
  });

  if (sidebar) {
    sidebar.querySelectorAll("a").forEach((a, i) => {
      a.addEventListener("click", (e) => {
        e.preventDefault();
        showPage(i);
      });
    });
  }

  // ── Initial navigation ────────────────────────────────────────────────
  function navigateToHash() {
    const hash = location.hash.slice(1);
    const idx = pages.findIndex(p => p.id === hash);
    showPage(idx >= 0 ? idx : 0, false);
  }

  window.addEventListener("popstate", navigateToHash);
  navigateToHash();
}
