.PHONY: all build release debug python python-release python-dev clean test help \
       app icon sign notarize dmg

# Default target
all: release python-release

# ── Build ──────────────────────────────────────────────────────────

build: debug

debug:
	cargo build

release:
	cargo build --release

python: python-release

python-release:
	cd crates/pymol-python && maturin build --release

python-dev:
	cd crates/pymol-python && maturin develop

test:
	cargo test

clean:
	cargo clean
	rm -rf crates/pymol-python/target
	rm -rf target/wheels
	rm -rf target/app
	rm -rf target/PyMOL-RS.dmg

run:
	./target/release/pymol-rs

# ── macOS App Bundle ───────────────────────────────────────────────

APP_NAME     := PyMOL-RS
APP_DIR      := target/app/$(APP_NAME).app
BUNDLE_ID    := me.yakovlev.pymol-rs
VERSION      := 0.1.0
ICON_SRC     := images/pymol-rs.png
ICONSET      := target/app/AppIcon.iconset
ICNS         := target/app/AppIcon.icns
BINARY       := target/release/pymol-rs
PLIST_TPL    := macos/Info.plist.in

# .env is sourced by sign/notarize/dmg targets at shell level.
# It must define: PYMOL_RS_APPLE_TEAMID, PYMOL_RS_APPLE_EMAIL, PYMOL_RS_APPLE_APP_PASS
ENV_FILE := .env

# Generate .icns from source PNG
icon: $(ICON_SRC)
	@echo "── Generating .icns ──"
	@mkdir -p $(ICONSET)
	sips -z   16   16 $(ICON_SRC) --out $(ICONSET)/icon_16x16.png      >/dev/null
	sips -z   32   32 $(ICON_SRC) --out $(ICONSET)/icon_16x16@2x.png   >/dev/null
	sips -z   32   32 $(ICON_SRC) --out $(ICONSET)/icon_32x32.png      >/dev/null
	sips -z   64   64 $(ICON_SRC) --out $(ICONSET)/icon_32x32@2x.png   >/dev/null
	sips -z  128  128 $(ICON_SRC) --out $(ICONSET)/icon_128x128.png    >/dev/null
	sips -z  256  256 $(ICON_SRC) --out $(ICONSET)/icon_128x128@2x.png >/dev/null
	sips -z  256  256 $(ICON_SRC) --out $(ICONSET)/icon_256x256.png    >/dev/null
	sips -z  512  512 $(ICON_SRC) --out $(ICONSET)/icon_256x256@2x.png >/dev/null
	sips -z  512  512 $(ICON_SRC) --out $(ICONSET)/icon_512x512.png    >/dev/null
	sips -z 1024 1024 $(ICON_SRC) --out $(ICONSET)/icon_512x512@2x.png >/dev/null
	iconutil -c icns $(ICONSET) -o $(ICNS)
	@rm -rf $(ICONSET)
	@echo "✓ $(ICNS)"

# Build .app bundle (release build + icon + Info.plist)
app: release icon
	@echo "── Assembling $(APP_NAME).app ──"
	@mkdir -p $(APP_DIR)/Contents/MacOS
	@mkdir -p $(APP_DIR)/Contents/Resources
	cp $(BINARY) $(APP_DIR)/Contents/MacOS/pymol-rs
	cp $(ICNS)   $(APP_DIR)/Contents/Resources/AppIcon.icns
	sed -e 's/@APP_NAME@/$(APP_NAME)/g' \
	    -e 's/@BUNDLE_ID@/$(BUNDLE_ID)/g' \
	    -e 's/@VERSION@/$(VERSION)/g' \
	    $(PLIST_TPL) > $(APP_DIR)/Contents/Info.plist
	@echo "✓ $(APP_DIR)"

# Code-sign (requires Apple Developer certificate)
sign: app
	@echo "── Signing $(APP_NAME).app ──"
	@test -f $(ENV_FILE) || { echo "ERROR: $(ENV_FILE) not found"; exit 1; }
	. ./$(ENV_FILE) && \
	  test -n "$$PYMOL_RS_APPLE_TEAMID" || { echo "ERROR: PYMOL_RS_APPLE_TEAMID not set in $(ENV_FILE)"; exit 1; } && \
	  codesign --force --options runtime --sign "Developer ID Application: $$PYMOL_RS_APPLE_TEAMID" $(APP_DIR) && \
	  codesign --verify --verbose $(APP_DIR)
	@echo "✓ Signed"

# Notarize with Apple
notarize: sign
	@echo "── Notarizing $(APP_NAME).app ──"
	@test -f $(ENV_FILE) || { echo "ERROR: $(ENV_FILE) not found"; exit 1; }
	. ./$(ENV_FILE) && \
	  test -n "$$PYMOL_RS_APPLE_EMAIL"    || { echo "ERROR: PYMOL_RS_APPLE_EMAIL not set in $(ENV_FILE)"; exit 1; } && \
	  test -n "$$PYMOL_RS_APPLE_TEAMID"   || { echo "ERROR: PYMOL_RS_APPLE_TEAMID not set in $(ENV_FILE)"; exit 1; } && \
	  test -n "$$PYMOL_RS_APPLE_APP_PASS" || { echo "ERROR: PYMOL_RS_APPLE_APP_PASS not set in $(ENV_FILE)"; exit 1; } && \
	  ditto -c -k --keepParent $(APP_DIR) target/app/$(APP_NAME).zip && \
	  xcrun notarytool submit target/app/$(APP_NAME).zip \
	      --apple-id "$$PYMOL_RS_APPLE_EMAIL" \
	      --team-id "$$PYMOL_RS_APPLE_TEAMID" \
	      --password "$$PYMOL_RS_APPLE_APP_PASS" \
	      --wait && \
	  xcrun stapler staple $(APP_DIR) && \
	  rm target/app/$(APP_NAME).zip
	@echo "✓ Notarized & stapled"

# Create distributable DMG
dmg: notarize
	@echo "── Creating DMG ──"
	hdiutil create -volname "$(APP_NAME)" \
	    -srcfolder target/app \
	    -ov -format UDZO \
	    target/$(APP_NAME).dmg
	@echo "✓ target/$(APP_NAME).dmg"

# ── Help ───────────────────────────────────────────────────────────

help:
	@echo "Available targets:"
	@echo "  all            - Build release binaries and Python wheel"
	@echo "  build/debug    - Build Rust workspace (debug)"
	@echo "  release        - Build Rust workspace (release)"
	@echo "  run            - Run release version of GUI interface"
	@echo "  python         - Build Python wheel (release)"
	@echo "  python-dev     - Install Python package in development mode"
	@echo "  test           - Run tests"
	@echo "  clean          - Clean all build artifacts"
	@echo ""
	@echo "  app            - Build macOS .app bundle (release + icon)"
	@echo "  sign           - Code-sign the .app (needs PYMOL_RS_APPLE_TEAMID)"
	@echo "  notarize       - Notarize with Apple (needs _EMAIL, _TEAMID, _APP_PASS)"
	@echo "  dmg            - Create distributable DMG (sign + notarize + package)"
	@echo ""
	@echo "Environment variables for signing:"
	@echo "  PYMOL_RS_APPLE_TEAMID   - Apple Developer Team ID"
	@echo "  PYMOL_RS_APPLE_EMAIL    - Apple ID email"
	@echo "  PYMOL_RS_APPLE_APP_PASS - App-specific password"
	@echo ""
	@echo "  help           - Show this help message"
