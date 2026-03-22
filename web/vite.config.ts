import { defineConfig, Plugin } from "vite";
import { resolve } from "path";
import { copyFileSync } from "fs";

/**
 * Vite plugin that prevents .wasm files from being inlined as base64 data URIs
 * in library mode. Rewrites wasm-bindgen's `new URL('...wasm', import.meta.url)`
 * so Vite doesn't process it, then copies the .wasm file to the output directory.
 */
function wasmExternalPlugin(): Plugin {
  let isBuild = false;

  return {
    name: "wasm-external",
    enforce: "pre",
    config(_cfg, { command }) {
      isBuild = command === "build";
    },
    // Rewrite the new URL() pattern so Vite's asset pipeline doesn't inline it
    transform(code, id) {
      if (!isBuild) return null;
      if (!code.includes("pymol_web_bg.wasm")) return null;

      // Replace: new URL('pymol_web_bg.wasm', import.meta.url)
      // Wrap import.meta.url in an identity function to bypass Vite's static
      // analysis (which would inline the wasm as a data URI in library mode),
      // while preserving the correct module-relative URL at runtime.
      const rewritten = code.replace(
        /new\s+URL\(\s*['"]pymol_web_bg\.wasm['"]\s*,\s*import\.meta\.url\s*\)/g,
        "new URL('pymol_web_bg.wasm', (x => x)(import.meta.url))",
      );

      if (rewritten !== code) {
        return { code: rewritten, map: null };
      }
      return null;
    },
    // Copy the wasm file to dist after build
    closeBundle() {
      if (!isBuild) return;
      try {
        copyFileSync(
          resolve(__dirname, "pkg/pymol_web_bg.wasm"),
          resolve(__dirname, "dist/pymol_web_bg.wasm"),
        );
      } catch {
        // wasm file may not exist during dev
      }
    },
  };
}

export default defineConfig({
  root: ".",
  plugins: [wasmExternalPlugin()],
  build: {
    lib: {
      entry: resolve(__dirname, "ts/src/index.ts"),
      name: "PyMolRS",
      fileName: "pymol-rs-viewer",
      formats: ["es"],
    },
    outDir: "dist",
    emptyOutDir: true,
  },
  server: {
    headers: {
      // Required for SharedArrayBuffer (used by some WASM features)
      "Cross-Origin-Opener-Policy": "same-origin",
      "Cross-Origin-Embedder-Policy": "require-corp",
    },
  },
});
