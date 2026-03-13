import { defineConfig } from "vite";
import { resolve } from "path";

export default defineConfig({
  root: ".",
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
