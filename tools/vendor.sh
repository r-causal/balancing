#!/bin/sh
#
# Produce the vendored source archive shipped with CRAN releases.
#
# Run from the package root whenever src/rust/Cargo.lock changes. The archive
# lets the package build offline: at install time Makevars extracts it and
# copies tools/config.toml into src/.cargo/config.toml to redirect crates-io to
# the vendored directory.
#
# cargo-vendor-filterer is preferred because it strips tests, documentation, and
# platforms the package does not target while keeping valid checksums, which
# keeps the archive within CRAN's size budget. It falls back to plain
# `cargo vendor` when the filterer is not installed.

set -eu

RUST_DIR="src/rust"
VENDOR_DIR="${RUST_DIR}/vendor"
ARCHIVE="${RUST_DIR}/vendor.tar.xz"

cd "$(dirname "$0")/.."

rm -rf "${VENDOR_DIR}"

if command -v cargo-vendor-filterer >/dev/null 2>&1; then
  cargo vendor-filterer \
    --manifest-path "${RUST_DIR}/Cargo.toml" \
    "${VENDOR_DIR}"
else
  echo "cargo-vendor-filterer not found; falling back to cargo vendor" >&2
  cargo vendor \
    --manifest-path "${RUST_DIR}/Cargo.toml" \
    "${VENDOR_DIR}"
fi

tar cJf "${ARCHIVE}" -C "${RUST_DIR}" vendor
rm -rf "${VENDOR_DIR}"

echo "Wrote ${ARCHIVE}"
