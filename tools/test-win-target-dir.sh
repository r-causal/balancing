#!/bin/sh
#
# Exercise the cargo target directory decision in tools/win-target-dir.sh.
#
# Not an R test and not shipped: .Rbuildignore keeps this out of the source
# tarball, while win-target-dir.sh itself has to travel in it. Run it directly,
# from anywhere:
#
#   sh tools/test-win-target-dir.sh
#
# Nothing here is Windows-specific. The decision is string and directory handling
# that behaves the same everywhere, which is the point: a Windows install is the
# only other way to reach this code, and it exercises exactly one of the paths
# below per run.
#
# No set -e and no set -u, deliberately. configure.win runs under neither, and the
# function reads temporary-root variables that are routinely unset, so a harness
# under -u would fail where the real caller does not. Failures are counted instead
# of aborting, so one broken case does not hide the state of the rest.

cd "$(dirname "$0")/.." || exit 1
. ./tools/win-target-dir.sh || exit 1

WORK="$(mktemp -d)" || exit 1
LOG="${WORK}/log"
trap 'rm -rf "${WORK}"' EXIT INT TERM

failures=0
checks=0

# A fresh root per case. Outside the dev profile the relocated tree is named for
# the process id, which is the same for every call in one run of this harness, so
# a root reused across two such cases would have the second decline the tree the
# first created. That is the behaviour being tested elsewhere, and anywhere else
# it would quietly make a case measure something other than its name.
#
# mktemp rather than a counter this function increments, because every call site
# is a command substitution and so runs in a subshell, where an increment is lost
# on return. That version handed every case the same root and looked, from the
# failures, like the code under test refusing to relocate.
case_root() {
  mktemp -d "${WORK}/rXXXXXX"
}

# A directory whose full path is exactly the requested number of characters, so
# that a case about lengths does not rest on how long mktemp's answer happened to
# be on this machine.
root_of_length() {
  rol_pad=""
  rol_want=$(( $1 - ${#WORK} - 1 ))
  while [ "${#rol_pad}" -lt "${rol_want}" ]; do
    rol_pad="${rol_pad}x"
  done
  mkdir -p "${WORK}/${rol_pad}" || exit 1
  printf '%s' "${WORK}/${rol_pad}"
}

# A package directory long enough to force relocation. It is never created: the
# decision measures this string and checksums it, and touches nothing under it.
LONG_PKG="/c/Users/somebody/AppData/Local/Temp/RtmpAbCdEf/R.INSTALL1a2b3c4d5e/balancing"
SHORT_PKG="/c/b"

# Run one case. Every temporary-root variable is set for every case, so that a
# value left behind by an earlier one cannot decide a later one.
run_case() {
  TMPDIR="$1"
  ORIGINAL_TEMP="$2"
  ORIGINAL_TMP="$3"
  TEMP="$4"
  TMP="$5"
  TARGET_DIR=""
  select_target_dir "$6" "$7" >"${LOG}" 2>&1
}

ok() {
  checks=$((checks + 1))
  if [ "$2" = "$3" ]; then
    printf 'ok   %s\n' "$1"
  else
    failures=$((failures + 1))
    printf 'FAIL %s\n       expected: %s\n         actual: %s\n' "$1" "$2" "$3"
  fi
}

# Assert on the shape of the answer rather than on a name holding a process id or
# a checksum, which no test can predict.
ok_under() {
  checks=$((checks + 1))
  case "$3" in
    "$2"/rt-*)
      printf 'ok   %s\n' "$1"
      ;;
    *)
      failures=$((failures + 1))
      printf 'FAIL %s\n       expected: under %s\n         actual: %s\n' "$1" "$2" "$3"
      ;;
  esac
}

ok_logged() {
  checks=$((checks + 1))
  if grep -q -- "$2" "${LOG}"; then
    printf 'ok   %s\n' "$1"
  else
    failures=$((failures + 1))
    printf 'FAIL %s\n       log did not mention: %s\n       log was:\n%s\n' "$1" "$2" "$(sed 's/^/         /' "${LOG}")"
  fi
}

# --- the in-package tree is kept whenever it fits --------------------------

run_case "${WORK}" "" "" "" "" "${SHORT_PKG}" release
ok "a short package directory keeps the in-package tree" '$(CURDIR)/rust/target' "${TARGET_DIR}"

# --- relocation ------------------------------------------------------------

root="$(case_root)"
run_case "${root}" "" "" "" "" "${LONG_PKG}" release
ok_under "a long package directory relocates into TMPDIR" "${root}" "${TARGET_DIR}"
ok "the relocated directory was created" "yes" "$([ -d "${TARGET_DIR}" ] && echo yes || echo no)"

# A root is only worth taking when it is shorter than what it replaces. This one
# is longer than the package directory, so the loop declines it.
run_case "$(root_of_length $(( ${#LONG_PKG} + 40 )))" "" "" "" "" "${LONG_PKG}" release
ok "a root no shorter than the package tree is declined" '$(CURDIR)/rust/target' "${TARGET_DIR}"

# --- the candidate chain ---------------------------------------------------

root="$(case_root)"
run_case "" "" "" "${root}" "$(case_root)" "${LONG_PKG}" release
ok_under "an empty root is skipped for the next one" "${root}" "${TARGET_DIR}"

root="$(case_root)"
run_case "" "${root}" "" "$(case_root)" "" "${LONG_PKG}" release
ok_under "ORIGINAL_TEMP is preferred to the TEMP that replaced it" "${root}" "${TARGET_DIR}"

root="$(case_root)"
run_case "" "" "${root}" "" "$(case_root)" "${LONG_PKG}" release
ok_under "ORIGINAL_TMP is preferred to the TMP that replaced it" "${root}" "${TARGET_DIR}"

root="$(case_root)"
run_case "${WORK}/missing-root" "" "" "${root}" "" "${LONG_PKG}" release
ok_under "a root that does not exist is skipped, not created" "${root}" "${TARGET_DIR}"
ok "and it was not created" "no" "$([ -d "${WORK}/missing-root" ] && echo yes || echo no)"

root="$(case_root)"
run_case "${root}/" "" "" "" "" "${LONG_PKG}" release
ok_under "a trailing slash on a root does not double" "${root}" "${TARGET_DIR}"

# --- the character allowlist -----------------------------------------------

for bad in 'with space' 'with&ersand' 'with$dollar' 'with(paren)' 'with;semi'; do
  mkdir -p "${WORK}/${bad}"
  root="$(case_root)"
  run_case "${WORK}/${bad}" "" "" "${root}" "" "${LONG_PKG}" release
  ok_under "a root holding '${bad}' is refused for the next root" "${root}" "${TARGET_DIR}"
done

# Above ASCII is text, not syntax, and is accepted. An account named José is
# enough to produce such a root, and refusing it would put the build back on the
# in-package tree that the path limit rules out.
mkdir -p "${WORK}/josé"
run_case "${WORK}/josé" "" "" "" "" "${LONG_PKG}" release
ok_under "a root holding a non-ASCII character is accepted" "${WORK}/josé" "${TARGET_DIR}"

# Backslashes are what TEMP and TMP actually arrive with. The separator has to be
# normalised before the allowlist sees it, or every Windows root would be refused
# for holding a backslash.
root="$(case_root)"
run_case "" "" "" "$(printf '%s' "${root}" | tr '/' '\\')" "" "${LONG_PKG}" release
ok_under "a root written with backslashes is normalised, not refused" "${root}" "${TARGET_DIR}"

# --- the two naming schemes ------------------------------------------------

# A development build has to land on the same name twice, or every load_all()
# rebuilds the dependency graph from scratch. One root throughout, since reuse is
# the whole point of the dev name.
dev_root="$(case_root)"
run_case "${dev_root}" "" "" "" "" "${LONG_PKG}" dev
first_dev="${TARGET_DIR}"
run_case "${dev_root}" "" "" "" "" "${LONG_PKG}" dev
ok "a dev build reuses the tree the last dev build of this checkout made" "${first_dev}" "${TARGET_DIR}"
ok_logged "and says that it is reusing it" "reusing the cargo target directory"

# Distinct per checkout, because this project is developed in git worktrees.
run_case "${dev_root}" "" "" "" "" "${LONG_PKG}-other" dev
ok "another checkout gets its own dev tree" "no" "$([ "${TARGET_DIR}" = "${first_dev}" ] && echo yes || echo no)"

# An install must never inherit a tree, so a name already taken is not adopted
# outside the dev profile. Take one, plant a staticlib in it, and ask again with
# the same root offered first and a fresh one behind it.
rel_root="$(case_root)"
next_root="$(case_root)"
run_case "${rel_root}" "" "" "${next_root}" "" "${LONG_PKG}" release
planted="${TARGET_DIR}"
ok_under "a release build relocates" "${rel_root}" "${planted}"
touch "${planted}/libbalancing_savvy.a"
run_case "${rel_root}" "" "" "${next_root}" "" "${LONG_PKG}" release
ok "a release build does not adopt a tree already standing there" "no" "$([ "${TARGET_DIR}" = "${planted}" ] && echo yes || echo no)"
ok_under "it takes the next root instead" "${next_root}" "${TARGET_DIR}"

# A dev tree that cannot be written to is no use either: cargo would fail inside
# the recipe, long after the loop that could have tried another root.
unwritable_root="$(case_root)"
fallback_root="$(case_root)"
run_case "${unwritable_root}" "" "" "" "" "${LONG_PKG}" dev
chmod a-w "${TARGET_DIR}"
run_case "${unwritable_root}" "" "" "${fallback_root}" "" "${LONG_PKG}" dev
ok_under "an unwritable dev tree is skipped for the next root" "${fallback_root}" "${TARGET_DIR}"
ok_logged "and says why" "it is not writable"
chmod u+w "${unwritable_root}"/rt-* 2>/dev/null

# --- reporting -------------------------------------------------------------

root="$(case_root)"
run_case "${root}" "" "" "" "" "${LONG_PKG}" release
ok_logged "the chosen directory is reported with its length" "using cargo target directory"
ok_logged "and the reason it moved" "moved out of the package tree"

# Nothing on offer is short enough, so the shortest available is still taken and
# the log says the limit may be exceeded rather than pretending otherwise. The
# root is sized past the budget but well inside the package tree it replaces.
deep_pkg="${LONG_PKG}$(printf '%0100d' 0)"
run_case "$(root_of_length $(( TARGET_DIR_MAX + 30 )))" "" "" "" "" "${deep_pkg}" release
ok_logged "a relocation that is still too long says so" "still over"

printf '\n%s checks, %s failures\n' "${checks}" "${failures}"
[ "${failures}" -eq 0 ] || exit 1
