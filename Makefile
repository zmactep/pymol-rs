.PHONY: all build release debug python python-release python-dev clean test help

# Default target
all: release python-release

# Build Rust workspace (debug)
build: debug

debug:
	cargo build

# Build Rust workspace (release)
release:
	cargo build --release

# Build Python wheel (release)
python: python-release

python-release:
	cd crates/pymol-python && maturin build --release

# Build and install Python package in development mode
python-dev:
	cd crates/pymol-python && maturin develop

# Run tests
test:
	cargo test

# Clean build artifacts
clean:
	cargo clean
	rm -rf crates/pymol-python/target
	rm -rf target/wheels

# Show help
help:
	@echo "Available targets:"
	@echo "  all            - Build release binaries and Python wheel"
	@echo "  build/debug    - Build Rust workspace (debug)"
	@echo "  release        - Build Rust workspace (release)"
	@echo "  python         - Build Python wheel (release)"
	@echo "  python-dev     - Install Python package in development mode"
	@echo "  test           - Run tests"
	@echo "  clean          - Clean all build artifacts"
	@echo "  help           - Show this help message"
