.PHONY: all build release debug test clean run help \
       python python-release python-dev \
       plugins plugins-install \
       icon app app-full \
       sign sign-full notarize dmg dmg-full \
       bundle-windows \
       web-build web-dev web-clean \
       widget-assets widget-build \
       version

# ── Variables ─────────────────────────────────────────────────────

APP_NAME       := PyMOL-RS
APP_DIR        := target/app/$(APP_NAME).app
BUNDLE_ID      := me.yakovlev.pymol-rs
VERSION        := 0.3.2
ICON_SRC       := images/pymol-rs.png
ICONSET        := target/app/AppIcon.iconset
ICNS           := target/app/AppIcon.icns
BINARY         := target/release/pymol-rs
PLIST_TPL      := macos/Info.plist.in
PYTHON_VERSION := 3.13
ENV_FILE       := .env
PLUGIN_INSTALL_DIR ?= $(HOME)/.pymol-rs/plugins
BUNDLE_DIR     := target/bundle/$(APP_NAME)

ifeq ($(OS),Windows_NT)
MKDIRP = powershell -NoProfile -Command "$$null = New-Item -ItemType Directory -Force -Path"
else
MKDIRP = mkdir -p
endif

# Python installation root (uv-managed preferred — smaller than system Python)
ifeq ($(OS),Windows_NT)
PYTHON_DIST = $(shell powershell -NoProfile -Command "Split-Path (uv python find $(PYTHON_VERSION))")
else
PYTHON_DIST = $(shell uv python find $(PYTHON_VERSION) 2>/dev/null | sed 's|/bin/python[0-9.]*$$||')
endif

# ── Build ─────────────────────────────────────────────────────────

all: release python-release

build: debug

debug:
	cargo build

release:
	cargo build --release

test:
	cargo test

run:
	./target/release/pymol-rs

clean:
	cargo clean
	rm -rf python/target target/wheels target/app target/dmg-stage target/$(APP_NAME).dmg
	rm -rf web/pkg web/dist python/pymol_rs/widget/static

# ── Python ────────────────────────────────────────────────────────

python: python-release

python-release: widget-assets
	cd python && maturin build --release

python-dev: widget-assets
	cd python && maturin develop

# ── Plugins ───────────────────────────────────────────────────────

PLUGIN_CRATES := -p raytracer-plugin -p hello-plugin -p python-plugin -p toolbar-plugin
ifneq ($(OS),Windows_NT)
PLUGIN_CRATES += -p ipc-plugin
endif

plugins:
	cargo build --release $(PLUGIN_CRATES)
	$(MKDIRP) target/release/plugins
ifeq ($(OS),Windows_NT)
	powershell -NoProfile -Command "Copy-Item 'target/release/*_plugin.dll' 'target/release/plugins/' -ErrorAction SilentlyContinue"
	powershell -NoProfile -Command "Copy-Item 'plugins/*/*.deps' 'target/release/plugins/' -ErrorAction SilentlyContinue"
else
	cp target/release/lib*_plugin.dylib target/release/plugins/ 2>/dev/null || \
	cp target/release/lib*_plugin.so    target/release/plugins/ 2>/dev/null || true
endif

plugins-install: plugins
	$(MKDIRP) $(PLUGIN_INSTALL_DIR)
	cp target/release/plugins/* $(PLUGIN_INSTALL_DIR)/

# ── macOS App Bundle ──────────────────────────────────────────────

icon: $(ICON_SRC)
	@echo "── Generating .icns ──"
	@mkdir -p $(ICONSET)
	@for size in 16 32 64 128 256 512 1024; do \
	    sips -z $$size $$size $(ICON_SRC) --out $(ICONSET)/icon_$${size}x$${size}.png >/dev/null 2>&1; \
	done
	@sips -z   32   32 $(ICON_SRC) --out $(ICONSET)/icon_16x16@2x.png    >/dev/null
	@sips -z   64   64 $(ICON_SRC) --out $(ICONSET)/icon_32x32@2x.png    >/dev/null
	@sips -z  256  256 $(ICON_SRC) --out $(ICONSET)/icon_128x128@2x.png  >/dev/null
	@sips -z  512  512 $(ICON_SRC) --out $(ICONSET)/icon_256x256@2x.png  >/dev/null
	@sips -z 1024 1024 $(ICON_SRC) --out $(ICONSET)/icon_512x512@2x.png  >/dev/null
	@rm -f $(ICONSET)/icon_64x64.png $(ICONSET)/icon_1024x1024.png
	iconutil -c icns $(ICONSET) -o $(ICNS)
	@rm -rf $(ICONSET)
	@echo "✓ $(ICNS)"

# Shared: assemble base .app structure (binary + icon + Info.plist)
define assemble-app
	@rm -rf $(APP_DIR)
	@mkdir -p $(APP_DIR)/Contents/MacOS $(APP_DIR)/Contents/Resources
	cp $(BINARY) $(APP_DIR)/Contents/MacOS/pymol-rs
	cp $(ICNS)   $(APP_DIR)/Contents/Resources/AppIcon.icns
	sed -e 's/@APP_NAME@/$(APP_NAME)/g' \
	    -e 's/@BUNDLE_ID@/$(BUNDLE_ID)/g' \
	    -e 's/@VERSION@/$(VERSION)/g' \
	    $(PLIST_TPL) > $(APP_DIR)/Contents/Info.plist
	@rm -f $(ICNS)
endef

# Minimal .app (binary only)
app: release icon
	@echo "── Assembling $(APP_NAME).app ──"
	$(assemble-app)
	@echo "✓ $(APP_DIR)"

# Full .app (binary + plugins + Python + venv)
app-full: release plugins icon python-release
	@echo "── Creating bundled venv (uv + Python $(PYTHON_VERSION)) ──"
	@rm -rf target/app/python-venv
	uv venv --python $(PYTHON_VERSION) target/app/python-venv
	uv pip install --python target/app/python-venv/bin/python3 \
	    $$(ls python/target/wheels/pymol_rs-*.whl | head -1)
	@echo "── Assembling full $(APP_NAME).app ──"
	$(assemble-app)
	@mkdir -p $(APP_DIR)/Contents/PlugIns
	cp target/release/lib*_plugin.dylib $(APP_DIR)/Contents/PlugIns/ 2>/dev/null || true
	cp target/release/lib*_plugin.so    $(APP_DIR)/Contents/PlugIns/ 2>/dev/null || true
	cp -R $(PYTHON_DIST)               $(APP_DIR)/Contents/Resources/python
	cp -R target/app/python-venv        $(APP_DIR)/Contents/Resources/python-venv
	@# Fix venv python symlinks: point to bundled Python, not system Python
	@cd $(APP_DIR)/Contents/Resources/python-venv/bin && \
	  rm -f python python3 python$(PYTHON_VERSION) && \
	  ln -s ../../python/bin/python3 python && \
	  ln -s python python3 && \
	  ln -s python python$(PYTHON_VERSION)
	bash macos/fix-dylib-paths.sh $(APP_DIR)
	@echo "✓ $(APP_DIR) (full bundle)"

# ── Code Signing ──────────────────────────────────────────────────

sign: app
	@echo "── Signing $(APP_NAME).app ──"
	@if [ ! -f $(ENV_FILE) ]; then echo "⚠ Signing skipped: $(ENV_FILE) not found"; exit 0; fi; \
	. ./$(ENV_FILE); \
	if [ -z "$$PYMOL_RS_APPLE_TEAMID" ]; then echo "⚠ Signing skipped: PYMOL_RS_APPLE_TEAMID not set"; exit 0; fi; \
	codesign --force --options runtime \
	    --sign "Developer ID Application: $$PYMOL_RS_APPLE_TEAMID" $(APP_DIR) && \
	codesign --verify --verbose $(APP_DIR) && \
	echo "✓ Signed"

sign-full: app-full
	@echo "── Deep-signing $(APP_NAME).app ──"
	@if [ ! -f $(ENV_FILE) ]; then echo "⚠ Signing skipped: $(ENV_FILE) not found"; exit 0; fi; \
	. ./$(ENV_FILE); \
	if [ -z "$$PYMOL_RS_APPLE_TEAMID" ]; then echo "⚠ Signing skipped: PYMOL_RS_APPLE_TEAMID not set"; exit 0; fi; \
	IDENTITY="Developer ID Application: $$PYMOL_RS_APPLE_TEAMID"; \
	find $(APP_DIR)/Contents/Resources/python \
	     $(APP_DIR)/Contents/Resources/python-venv \
	     $(APP_DIR)/Contents/PlugIns \
	  -type f \( -name '*.dylib' -o -name '*.so' \) 2>/dev/null | \
	  xargs -I{} codesign --force --options runtime --sign "$$IDENTITY" {} && \
	codesign --force --options runtime --sign "$$IDENTITY" \
	  $(APP_DIR)/Contents/Resources/python/bin/python3 && \
	codesign --force --options runtime --sign "$$IDENTITY" $(APP_DIR) && \
	codesign --verify --verbose --deep $(APP_DIR) && \
	echo "✓ Signed (deep)"

# ── Notarize & DMG ───────────────────────────────────────────────

notarize:
	@echo "── Notarizing $(APP_NAME).app ──"
	@if [ ! -f $(ENV_FILE) ]; then echo "⚠ Notarization skipped: $(ENV_FILE) not found"; exit 0; fi; \
	. ./$(ENV_FILE); \
	if [ -z "$$PYMOL_RS_APPLE_EMAIL" ] || [ -z "$$PYMOL_RS_APPLE_TEAMID" ] || [ -z "$$PYMOL_RS_APPLE_APP_PASS" ]; then \
	    echo "⚠ Notarization skipped: missing credentials in $(ENV_FILE)"; exit 0; \
	fi; \
	ditto -c -k --keepParent $(APP_DIR) target/app/$(APP_NAME).zip && \
	xcrun notarytool submit target/app/$(APP_NAME).zip \
	    --apple-id "$$PYMOL_RS_APPLE_EMAIL" \
	    --team-id "$$PYMOL_RS_APPLE_TEAMID" \
	    --password "$$PYMOL_RS_APPLE_APP_PASS" \
	    --wait && \
	xcrun stapler staple $(APP_DIR) && \
	echo "✓ Notarized & stapled" || \
	echo "⚠ Notarization failed — continuing without it"; \
	rm -f target/app/$(APP_NAME).zip

# Shared: create DMG from a clean staging directory (only .app + Applications link)
define create-dmg
	$(MAKE) notarize
	@echo "── Creating DMG ──"
	@rm -rf target/dmg-stage
	@mkdir -p target/dmg-stage
	@cp -R $(APP_DIR) target/dmg-stage/
	@ln -sf /Applications target/dmg-stage/Applications
	hdiutil create -volname "$(APP_NAME)" \
	    -srcfolder target/dmg-stage -ov -format UDZO \
	    target/$(APP_NAME).dmg
	@rm -rf target/dmg-stage
	@echo "✓ target/$(APP_NAME).dmg"
endef

dmg: sign
	$(create-dmg)

dmg-full: sign-full
	$(create-dmg)

# ── Windows Bundle ────────────────────────────────────────────────

bundle-windows: release plugins python-release
	@echo ── Assembling Windows bundle ──
	powershell -NoProfile -Command "if (Test-Path '$(BUNDLE_DIR)') { Remove-Item -Recurse -Force '$(BUNDLE_DIR)' }"
	$(MKDIRP) "$(BUNDLE_DIR)\plugins"
	powershell -NoProfile -Command "Copy-Item 'target/release/pymol-rs.exe' '$(BUNDLE_DIR)/'"
	powershell -NoProfile -Command "Copy-Item 'target/release/*_plugin.dll' '$(BUNDLE_DIR)/plugins/' -ErrorAction SilentlyContinue"
	powershell -NoProfile -Command "Copy-Item 'plugins/*/*.deps' '$(BUNDLE_DIR)/plugins/' -ErrorAction SilentlyContinue"
	powershell -NoProfile -Command "Copy-Item -Recurse '$(PYTHON_DIST)' '$(BUNDLE_DIR)/python'"
	uv venv --python $(PYTHON_VERSION) "$(BUNDLE_DIR)\python-venv"
	powershell -NoProfile -Command "uv pip install --python '$(BUNDLE_DIR)/python-venv/Scripts/python.exe' $$(Get-ChildItem 'python/target/wheels/pymol_rs-*.whl' | Select-Object -First 1).FullName"
	powershell -NoProfile -Command "Copy-Item 'windows/PyMOL-RS.cmd' '$(BUNDLE_DIR)/'"
	@echo ✓ $(BUNDLE_DIR)

# ── Web ──────────────────────────────────────────────────────────

web-build:
	cd web && npm install && npm run build

web-dev:
	cd web && npm install && npm run dev

web-clean:
	rm -rf web/pkg web/dist web/node_modules

# ── Widget ───────────────────────────────────────────────────────

# Copy pre-built web dist into the Python package (requires web-build to have run once)
widget-assets:
ifeq ($(OS),Windows_NT)
	@powershell -NoProfile -Command "if (-not (Test-Path 'web/dist/pymol_web_bg.wasm')) { Write-Host '── Web dist not found, building... ──'; & '$(MAKE)' web-build }"
	$(MKDIRP) python\pymol_rs\widget\static
	powershell -NoProfile -Command "Copy-Item 'web/dist/pymol-rs-viewer.js' 'python/pymol_rs/widget/static/'"
	powershell -NoProfile -Command "Copy-Item 'web/dist/pymol_web_bg.wasm' 'python/pymol_rs/widget/static/'"
	powershell -NoProfile -Command "Get-ChildItem 'web/dist/pymol_web-*.js' | Select-Object -First 1 | Copy-Item -Destination 'python/pymol_rs/widget/static/pymol_web_glue.js'"
else
	@if [ ! -f web/dist/pymol_web_bg.wasm ]; then \
		echo "── Web dist not found, building... ──"; \
		$(MAKE) web-build; \
	fi
	$(MKDIRP) python/pymol_rs/widget/static
	cp web/dist/pymol-rs-viewer.js python/pymol_rs/widget/static/
	cp web/dist/pymol_web_bg.wasm python/pymol_rs/widget/static/
	cp web/dist/pymol_web-*.js python/pymol_rs/widget/static/pymol_web_glue.js
endif

# Full rebuild: web + copy assets
widget-build: web-build widget-assets

# ── Version ──────────────────────────────────────────────────────

version:
	@if [ -z "$(V)" ]; then echo "Usage: make version V=X.Y.Z"; exit 1; fi
	bash scripts/version.sh $(V)

# ── Help ──────────────────────────────────────────────────────────

help:
	@echo "Build:"
	@echo "  all              Build release binaries and Python wheel"
	@echo "  debug / release  Build Rust workspace"
	@echo "  test             Run tests"
	@echo "  run              Launch GUI (release)"
	@echo "  python           Build Python wheel (release)"
	@echo "  python-dev       Install Python package in dev mode"
	@echo "  plugins          Build all plugins (release)"
	@echo "  plugins-install  Build + install plugins to ~/.pymol-rs/plugins"
	@echo "  clean            Clean all build artifacts"
	@echo ""
	@echo "macOS app bundle:"
	@echo "  app              Minimal .app (binary + icon)"
	@echo "  app-full         Full .app (+ plugins + Python + venv)"
	@echo "  sign / sign-full Code-sign the .app"
	@echo "  dmg / dmg-full   Distributable DMG (signed + notarized)"
	@echo ""
	@echo "Windows bundle:"
	@echo "  bundle-windows   Bundle (exe + plugins + Python + venv + launcher)"
	@echo ""
	@echo "Version:"
	@echo "  version V=X.Y.Z Update version across all crates and packages"
	@echo ""
	@echo "Web:"
	@echo "  web-build        Build WASM + TypeScript bundle"
	@echo "  web-dev          Dev server with hot reload"
	@echo "  web-clean        Clean web build artifacts"
	@echo ""
	@echo "Widget:"
	@echo "  widget-build     Build web + copy WASM assets for Jupyter widget"
	@echo ""
	@echo "Signing credentials (.env):"
	@echo "  PYMOL_RS_APPLE_TEAMID    Apple Developer Team ID"
	@echo "  PYMOL_RS_APPLE_EMAIL     Apple ID email"
	@echo "  PYMOL_RS_APPLE_APP_PASS  App-specific password"
