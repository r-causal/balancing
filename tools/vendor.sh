#!/bin/sh
#
# Produce the vendored source archive shipped with CRAN releases.
#
# Run from the package root whenever src/rust/Cargo.lock changes. The archive
# lets the package build offline: at install time Makevars extracts it and
# copies tools/config.toml into the config file of a build-local cargo home to
# redirect crates-io to the vendored directory.
#
# cargo-vendor-filterer is required rather than preferred. Plain `cargo vendor`
# has no way to leave out dependency kinds, and for this package that is where
# the size is: criterion, the benchmark harness, is a dev-dependency, and it
# brings plotters, web-sys, windows-sys and their graphs with it. A release build
# compiles none of that. Measured on the graph as of this writing, the unfiltered
# archive is 13.0 MB and the filtered one 7.3 MB, so a fallback to plain
# `cargo vendor` would quietly produce an archive nobody would want to submit.
# The fallback that used to be here did exactly that in CI, where the filterer
# was never installed and so the fallback was the only path ever taken.

set -eu

RUST_DIR="src/rust"
WORKSPACE_MANIFEST="${RUST_DIR}/Cargo.toml"
VENDOR_DIR="${RUST_DIR}/vendor"
ARCHIVE="${RUST_DIR}/vendor.tar.xz"

# The filterer is pointed at the savvy crate rather than at the workspace, and
# that is a workaround rather than a preference. Its dependency-kind filter shells
# out to `cargo tree --prefix none`, which prints a blank line between the trees
# of two workspace members, and its parser treats a line it cannot split into a
# name and a version as fatal: "Invalid output received from cargo tree:" with
# nothing after the colon. This workspace has two members, so the filter fails
# outright when given the workspace manifest. One member's tree carries no blank
# line, and balancing-savvy depends on balancing-core, so its tree is the whole of
# what a release build compiles. The vendored set is still the workspace's; the
# filter only decides which of it is reduced to a stub.
#
# What this gives up is a future workspace member that balancing-savvy does not
# reach. There is none today, and an offline build of one would fail outright
# rather than quietly miss anything.
SAVVY_MANIFEST="${RUST_DIR}/crates/savvy/Cargo.toml"

# Ceiling on the archive, in bytes: roughly a quarter above the size measured when
# this was written. That leaves the graph room to grow and still catches an
# exclusion below quietly ceasing to match, since losing the three osqp ones alone
# puts the archive back over 10 MB. CRAN asks that a source tarball stay near
# 5 MB, so this is a regression guard and not a target; closing the remaining
# distance is its own issue.
ARCHIVE_MAX_BYTES=9000000

cd "$(dirname "$0")/.."

if ! command -v cargo-vendor-filterer >/dev/null 2>&1; then
  echo "cargo-vendor-filterer is not installed, and this script does not produce" >&2
  echo "an archive without it. Install it with:" >&2
  echo "" >&2
  echo "  cargo install cargo-vendor-filterer --locked" >&2
  exit 1
fi

# The archive has to describe the dependency set the committed lockfile describes,
# or an install would compile crate versions no one recorded. `cargo vendor` takes
# --locked for that. cargo-vendor-filterer has no such flag, and passing one is an
# error rather than a no-op, which is how the flag came to be here on a code path
# CI never ran. Asking cargo for the metadata under --locked refuses a stale
# lockfile in the same way, and leaves the vendoring step nothing to re-resolve:
# a lockfile that needs no update does not get one.
cargo metadata --manifest-path "${WORKSPACE_MANIFEST}" --locked --format-version 1 >/dev/null

rm -rf "${VENDOR_DIR}"

# The osqp exclusions are what the size rests on: osqp-sys vendors the whole OSQP
# repository, whose built documentation site is 3 MB of the archive by itself. The
# demo executable's sources are deliberately not excluded. cmake builds it by
# default, since OSQP_BUILD_DEMO_EXE is ON unless the static library is off, so
# dropping osqp/examples fails the compile at add_executable; that was tried. The
# unit tests are guarded both by an option that defaults to off and by a
# top-level-project test, so osqp/tests is unreachable from a build of osqp-sys.
cargo vendor-filterer \
  --manifest-path "${SAVVY_MANIFEST}" \
  --keep-dep-kinds no-dev \
  --exclude-crate-path '*#tests' \
  --exclude-crate-path '*#benches' \
  --exclude-crate-path '*#examples' \
  --exclude-crate-path 'osqp-sys#osqp/site' \
  --exclude-crate-path 'osqp-sys#osqp/docs' \
  --exclude-crate-path 'osqp-sys#osqp/tests' \
  "${VENDOR_DIR}"

tar cJf "${ARCHIVE}" -C "${RUST_DIR}" vendor
rm -rf "${VENDOR_DIR}"

# wc -c rather than stat, whose flag for a file's size differs between GNU and BSD.
ARCHIVE_BYTES="$(wc -c < "${ARCHIVE}" | tr -d ' ')"

if [ "${ARCHIVE_BYTES}" -gt "${ARCHIVE_MAX_BYTES}" ]; then
  echo "${ARCHIVE} is ${ARCHIVE_BYTES} bytes, over the ${ARCHIVE_MAX_BYTES} this" >&2
  echo "script allows. Either the dependency graph grew or one of the exclusions" >&2
  echo "above stopped matching the paths it names, which cargo-vendor-filterer" >&2
  echo "reports as a warning of its own." >&2
  exit 1
fi

echo "Wrote ${ARCHIVE} (${ARCHIVE_BYTES} bytes)"
