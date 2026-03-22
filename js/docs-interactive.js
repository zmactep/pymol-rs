// ── Inline PDB fragments ────────────────────────────────────────────────────
const FRAGMENTS = {
  "alanine-dipeptide": `ATOM      1  N   ALA A   1       1.458   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.534   1.398   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       4.140   0.338  -0.170  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.521   2.166  -1.236  1.00  0.00           C
ATOM      6  H   ALA A   1       0.458   0.000   0.000  1.00  0.00           H
ATOM      7  HA  ALA A   1       1.631   1.940   0.893  1.00  0.00           H
ATOM      8  HB1 ALA A   1       0.431   2.180  -1.236  1.00  0.00           H
ATOM      9  HB2 ALA A   1       1.893   3.191  -1.236  1.00  0.00           H
ATOM     10  HB3 ALA A   1       1.893   1.646  -2.124  1.00  0.00           H
ATOM     11  N   ALA A   2       4.126   2.574   0.170  1.00  0.00           N
ATOM     12  CA  ALA A   2       5.590   2.668   0.170  1.00  0.00           C
ATOM     13  C   ALA A   2       6.141   4.090   0.170  1.00  0.00           C
ATOM     14  O   ALA A   2       5.390   5.062   0.170  1.00  0.00           O
ATOM     15  CB  ALA A   2       6.120   1.918   1.392  1.00  0.00           C
ATOM     16  OXT ALA A   2       7.370   4.174   0.170  1.00  0.00           O
ATOM     17  H   ALA A   2       3.570   3.400   0.308  1.00  0.00           H
ATOM     18  HA  ALA A   2       5.970   2.178  -0.735  1.00  0.00           H
ATOM     19  HB1 ALA A   2       5.740   0.896   1.392  1.00  0.00           H
ATOM     20  HB2 ALA A   2       7.210   1.918   1.392  1.00  0.00           H
ATOM     21  HB3 ALA A   2       5.740   2.401   2.294  1.00  0.00           H
END`,

  "adenosine": `ATOM      1  P   A   B   1      -0.817   4.551   0.735  1.00  0.00           P
ATOM      2  OP1 A   B   1      -2.093   5.301   0.580  1.00  0.00           O
ATOM      3  OP2 A   B   1       0.346   5.247   0.100  1.00  0.00           O
ATOM      4  O5' A   B   1      -0.467   4.318   2.273  1.00  0.00           O
ATOM      5  C5' A   B   1      -1.469   3.734   3.113  1.00  0.00           C
ATOM      6  C4' A   B   1      -0.883   3.432   4.472  1.00  0.00           C
ATOM      7  O4' A   B   1       0.272   2.582   4.310  1.00  0.00           O
ATOM      8  C3' A   B   1      -0.358   4.633   5.241  1.00  0.00           C
ATOM      9  O3' A   B   1      -1.371   5.283   5.982  1.00  0.00           O
ATOM     10  C2' A   B   1       0.662   4.016   6.195  1.00  0.00           C
ATOM     11  O2' A   B   1       0.058   3.317   7.270  1.00  0.00           O
ATOM     12  C1' A   B   1       1.290   2.925   5.326  1.00  0.00           C
ATOM     13  N9  A   B   1       2.517   3.388   4.640  1.00  0.00           N
ATOM     14  C8  A   B   1       2.620   4.210   3.539  1.00  0.00           C
ATOM     15  N7  A   B   1       3.878   4.434   3.175  1.00  0.00           N
ATOM     16  C5  A   B   1       4.644   3.734   4.085  1.00  0.00           C
ATOM     17  C6  A   B   1       6.048   3.600   4.168  1.00  0.00           C
ATOM     18  N6  A   B   1       6.872   4.152   3.273  1.00  0.00           N
ATOM     19  N1  A   B   1       6.533   2.876   5.196  1.00  0.00           N
ATOM     20  C2  A   B   1       5.680   2.326   6.078  1.00  0.00           C
ATOM     21  N3  A   B   1       4.360   2.371   6.086  1.00  0.00           N
ATOM     22  C4  A   B   1       3.896   3.097   5.058  1.00  0.00           C
END`
};

// ── Utilities ───────────────────────────────────────────────────────────────

function stringToUint8(str) {
  return new TextEncoder().encode(str);
}

function flashElement(el, className = "flash", ms = 400) {
  el.classList.add(className);
  setTimeout(() => el.classList.remove(className), ms);
}

/** Split a semicolon-delimited command string and execute each via viewer. */
async function execCmdString(viewer, cmdStr) {
  const cmds = cmdStr.split(";").map(s => s.trim()).filter(Boolean);
  for (const cmd of cmds) {
    const result = viewer.execute(cmd);
    if (result && result._loadPromise) await result._loadPromise;
  }
}

// ── Load intercept ──────────────────────────────────────────────────────────

const LOAD_RE = /^load\s+(.+)/i;
const LOAD_FRAGMENT_RE = /^load-fragment\s+(.+)/i;

function parseLoadArgs(argsStr) {
  const parts = argsStr.split(",").map(s => s.trim());
  const filename = parts[0];
  const name = parts[1] || filename.replace(/^.*\//, "").replace(/\.[^.]+$/, "");
  const format = parts[3] || filename.replace(/^.*\./, "").toLowerCase();
  return { filename, name, format };
}

function resolveAssetPath(filename, basePath) {
  if (filename.startsWith("http://") || filename.startsWith("https://")) return filename;
  if (!filename.includes("/")) return basePath + filename;
  return filename;
}

export function loadFragment(viewer, name) {
  const pdb = FRAGMENTS[name];
  if (!pdb) { console.warn(`Fragment "${name}" not found`); return; }
  viewer.loadData(stringToUint8(pdb), name, "pdb");
}

/**
 * Patches viewer.execute() to intercept `load` and `load-fragment` commands,
 * converting them to web-compatible loadUrl / loadData calls.
 */
export function patchLoadCommand(viewer, basePath = "../assets/pdb/") {
  const origExecute = viewer.execute.bind(viewer);

  viewer.execute = function(cmd) {
    const fragMatch = cmd.match(LOAD_FRAGMENT_RE);
    if (fragMatch) {
      const name = fragMatch[1].trim();
      loadFragment(viewer, name);
      return { messages: [{ level: "info", text: `Loaded fragment ${name}` }] };
    }
    const match = cmd.match(LOAD_RE);
    if (match) {
      const { filename, name, format } = parseLoadArgs(match[1]);
      const url = resolveAssetPath(filename, basePath);
      const promise = viewer.loadUrl(url, { name, format }).catch(err => {
        console.error(`Failed to load ${url}:`, err);
      });
      return { messages: [{ level: "info", text: `Loading ${name}...` }], _loadPromise: promise };
    }
    return origExecute(cmd);
  };
}

// ── SPA page navigation ─────────────────────────────────────────────────────

/**
 * Initializes SPA-style documentation with hash-routed sections.
 *
 * Each section is a `.doc-page` div with `id` and optional `data-viewer-setup`:
 *   { "pdb": "url", "cmds": [...] }
 *   { "fragment": "name", "cmds": [...] }
 *
 * Inline commands use `data-cmd` (semicolon-separated) with optional
 * `data-cmd-alt` for toggle behaviour.
 */
function forceResync(viewer) {
  if (viewer.core) {
    viewer.core.syncCanvasSize();
    const w = viewer.core.wasmViewer;
    if (w) w.resize(viewer.core.canvas.width, viewer.core.canvas.height);
  }
}

export function initDocPages(viewer) {
  // ── Inline command clicks ─────────────────────────────────────────────
  document.addEventListener("click", async (e) => {
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

    await execCmdString(viewer, cmdStr);
    flashElement(el);
  });

  // ── Page management ───────────────────────────────────────────────────
  const pages = Array.from(document.querySelectorAll(".doc-page"));
  if (pages.length === 0) return;

  const sidebar = document.querySelector(".docs-sidebar");
  let currentIdx = 0;

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

  function setupViewer(page) {
    forceResync(viewer);
    const attr = page.dataset.viewerSetup;
    if (!attr) return;
    try {
      const setup = JSON.parse(attr);

      const allCmds = setup.cmds || [];

      // Wait for WASM to finish computing geometry, then run all commands
      // in original order. View commands (orient/zoom) need correct bounds,
      // so we wait for idle before running anything.
      const runAll = () => {
        const wasm = viewer.core && viewer.core.wasmViewer;
        async function execAll() {
          for (const cmd of allCmds) {
            const result = viewer.execute(cmd);
            if (result && result._loadPromise) await result._loadPromise;
          }
        }
        if (!wasm) { execAll(); return; }
        let idleFrames = 0;
        let sawBusy = false;
        let totalFrames = 0;
        function poll() {
          totalFrames++;
          wasm.process_input();
          if (wasm.needs_redraw()) {
            wasm.render_frame();
            idleFrames = 0;
            sawBusy = true;
          } else {
            if (sawBusy) idleFrames++;
          }
          if (sawBusy ? idleFrames < 3 : totalFrames < 30) {
            requestAnimationFrame(poll);
          } else {
            execAll();
          }
        }
        requestAnimationFrame(poll);
      };

      if (setup.fragment) {
        viewer.execute("delete all");
        loadFragment(viewer, setup.fragment);
        runAll();
      } else if (setup.pdb) {
        viewer.execute("delete all");
        const pdbs = Array.isArray(setup.pdb) ? setup.pdb : [setup.pdb];
        (async () => {
          for (const url of pdbs) {
            await viewer.loadUrl(url, {
              name: url.replace(/^.*\//, "").replace(/\.[^.]+$/, "")
            });
          }
          runAll();
        })();
      } else {
        runAll();
      }
    } catch (e) {
      console.error("Invalid data-viewer-setup JSON:", e);
    }
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

  // Wait for the viewer canvas to reach a stable size before first render
  const canvas = document.querySelector(".docs-viewer-canvas");
  if (canvas) {
    let lastW = 0, lastH = 0, stableFrames = 0;
    function waitForStableSize() {
      const rect = canvas.getBoundingClientRect();
      if (rect.width > 0 && rect.height > 0 && rect.width === lastW && rect.height === lastH) {
        stableFrames++;
      } else {
        stableFrames = 0;
      }
      lastW = rect.width;
      lastH = rect.height;
      if (stableFrames >= 2) {
        forceResync(viewer);
        navigateToHash();
      } else {
        requestAnimationFrame(waitForStableSize);
      }
    }
    requestAnimationFrame(waitForStableSize);
  } else {
    navigateToHash();
  }
}
