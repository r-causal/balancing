#!/bin/sh
#
# Choose where cargo builds on Windows.
#
# Sourced by configure.win, which is the only thing that ships this. It lives in
# its own file so that tools/test-win-target-dir.sh can drive the decision with a
# synthetic environment: the logic below has two naming schemes, a character
# allowlist, a length test, an adoption branch and a fall-through, and a full
# Windows install is a slow and coarse way to exercise any of them. One defect in
# it, a /tmp last resort that no length test could measure, survived two reviews.
#
# POSIX sh with no local variables, so everything here is prefixed wtd_ and reset
# on entry rather than left to whatever a previous call set.
#
# Windows caps a path at 260 characters. osqp-sys, a hard dependency of
# balancing-core because the energy and cfd methods need OSQP, builds OSQP
# through cmake, and cmake writes its compiler probes and object files as much as
# 173 characters below the cargo target directory. So the target directory itself
# has to stay short. Installing from a source tarball unpacks the sources under
# the temporary directory, which leaves the in-package tree near 100 characters
# and puts cmake past the limit; that is why the same sources build from a short
# checkout and fail on r-universe.
#
# Keep the tree inside the package while it fits, because pkgbuild::compile_dll()
# and repeated load_all() reuse the cargo cache and a fresh directory per run
# would rebuild the Rust core from scratch every time. Relocate it into the
# temporary directory otherwise; clean_intermediate removes it for every profile
# but dev.
#
# The budget is the one for aarch64-pc-windows-gnullvm, the deeper of the two
# Windows targets. cmake reaches 146 characters below the target directory before
# the longest name it appends, the 27-character CMakeCXXCompilerABI.cpp.obj, so
# clearing cmake's own CMAKE_OBJECT_PATH_MAX warning threshold of 250 allows 77
# characters here and clearing the roughly 259 the OS turns out to enforce allows
# 85. 76 is under the stricter of the two.
#
# The margin against 85 earns its keep, because the length measured below is only
# an estimate of the length that matters. The limit applies to the path after
# Windows canonicalizes it, and an 8.3 component such as RUNNER~1 grows by three
# characters once cmake resolves it to runneradmin. Rtools' sh also reports a
# drive-letter path in msys form, which happens to be the same length as the
# Windows form make will use, but would not be for a path under the msys root.
TARGET_DIR_MAX=76

# Decide the cargo target directory and set TARGET_DIR to what configure.win
# substitutes into src/Makevars.win: either the literal $(CURDIR)/rust/target for
# the in-package tree, or an absolute path to a relocated one. Everything it
# decides is reported on stdout, since that output is the only record an install
# log carries of which tree a build used.
#
# Arguments: the package directory, and the cargo profile.
select_target_dir() {
  wtd_pkg_dir="$1"
  wtd_profile="$2"

  # make runs in src/, one level below the directory configure.win runs in.
  wtd_in_package="${wtd_pkg_dir}/src/rust/target"
  wtd_relocated=""

  if [ "${#wtd_in_package}" -gt "${TARGET_DIR_MAX}" ]; then
    # What a relocated tree is named depends on the profile, because the two
    # profiles want opposite things from the name.
    #
    # An install needs a name no other install can arrive at. The relocated tree
    # is not scoped to anything else: R exports R_SESSION_TMPDIR on Unix only, so
    # on Windows the roots below name the temporary directory itself and not a
    # per-session directory inside it, and nothing removes what is left there when
    # R exits. A shared name would then let one install inherit another's tree.
    # Two concurrent installs would also race, the first to finish removing the
    # tree the other is building in. The process id settles all of that, at the
    # cost of leaving an interrupted install's tree behind: the name it would be
    # recognised by is gone with the shell that chose it, but nothing can link
    # that tree either and temporary-directory cleaning collects it.
    #
    # A development build needs the opposite. pkgbuild::compile_dll() sets
    # DEBUG=true and installs from the checkout, so a contributor whose checkout
    # is deeper than 60 characters relocates on every load_all(); a name unique
    # per run would rebuild OSQP and the whole dependency graph each time and
    # leave every tree behind, since clean_intermediate keeps the tree for dev.
    # That is the cache loss this whole block exists to avoid. A checksum of the
    # package directory gives a name that is stable from one run to the next and
    # still distinct per checkout, which matters because this project is developed
    # in git worktrees and two checkouts on one machine must not build in the same
    # tree.
    #
    # The staleness a per-install name guards against does not follow the checksum
    # name into the dev profile. A dev tree can only be inherited by a later build
    # of the same checkout, never by an install of another version, and
    # clean_intermediate removes the staticlib for every profile, so the next dev
    # build re-runs cargo and cargo decides for itself what to reuse.
    #
    # The name stays terse, because the whole point of this block is the character
    # budget. The d marks a checksum, which keeps the two schemes from ever naming
    # the same directory: a process id is digits alone.
    wtd_tag="$$"
    if [ "${wtd_profile}" = "dev" ]; then
      # cksum is in POSIX and Rtools carries it, but nothing here depends on that:
      # when it answers with no checksum the per-run name is used instead and the
      # build still runs. That fallback is reported because it costs more than the
      # cache. A dev build then carries the per-run name, and clean_intermediate
      # keeps a dev tree, so every load_all() rebuilds the whole dependency graph
      # and leaves another cargo tree in the temporary directory that nothing later
      # goes back for.
      #
      # The first run of digits in cksum's answer, which keeps the checksum and
      # leaves behind the byte count beside it and any padding around either.
      wtd_digest="$(printf '%s' "${wtd_pkg_dir}" | cksum 2>/dev/null | sed -n 's/^[^0-9]*\([0-9][0-9]*\).*$/\1/p')"
      if [ -n "${wtd_digest}" ]; then
        wtd_tag="d${wtd_digest}"
      else
        echo "cksum reported no checksum of this directory, so the build tree is named for this run alone: no cargo cache will carry over to the next run, and each run will leave its own tree behind in the temporary directory"
      fi
    fi

    # Only the roots the environment names, and no /tmp behind them. /tmp is the one
    # root the length test below cannot measure: under Rtools' sh it is the msys root's
    # own tmp directory, C:/rtools45/tmp, so the string measured here comes out short by
    # the whole msys root, where every other root's msys spelling measures the same as
    # the Windows form the tools end up using. A build tree inside the toolchain
    # installation is also not what a temporary directory is for, since nothing ever
    # collects what is left there. The in-package tree is already the fallback. None of
    # this is about where a value came from: TEMP or TMP can still name /tmp, and for
    # such a root the length reported below is the msys one.
    #
    # ORIGINAL_TEMP and ORIGINAL_TMP come ahead of TEMP and TMP rather than behind
    # them, which is the opposite of how this was first written down. msys2's
    # /etc/profile, which Rtools inherits unchanged, overrides TEMP and TMP to /tmp
    # and exports the values it replaced under those names. So they are set only when
    # TEMP and TMP have already been overridden, and a root behind an overridden TEMP
    # would never be reached: /tmp exists, is writable, and is short, so the loop
    # takes it and stops. Ahead of them, a build started from an Rtools shell puts its
    # tree in a real temporary directory instead of inside the toolchain
    # installation, and a build started anywhere else does not set them and is
    # unaffected.
    for wtd_root in "${TMPDIR}" "${ORIGINAL_TEMP}" "${ORIGINAL_TMP}" "${TEMP}" "${TMP}"; do
      if [ -z "${wtd_root}" ]; then
        continue
      fi

      # TEMP and TMP arrive with backslashes. Cargo, cmake and make all take forward
      # slashes, and make would read a backslash as an escape. Normalise before the
      # test below, so that a separator is not taken for one of the characters that
      # test refuses.
      wtd_normalised="$(printf '%s' "${wtd_root}" | LC_ALL=C tr '\\' '/')"

      # Refuse a root holding a character that make or the shell would read as
      # something other than text, and try the next root instead. Windows allows far
      # more than this in a directory name and a local account is enough to produce it:
      # an account named "Jane Doe" gives a TEMP with a space in it, and one named
      # "A&B" a TEMP with an ampersand.
      #
      # An allowlist rather than a list of the characters known to break something,
      # because the value reaches the generated makefile as TARGET_DIR and make
      # expands it unquoted into PKG_LIBS and into the prerequisites of $(SHLIB),
      # where anything make or the shell reads as a separator or an operator splits
      # it. $ is worse still, since make expands it before any quoting can apply, so
      # rm -Rf "$(TARGET_DIR)" could name a directory that nothing here created. The
      # allowlist covers each of those at once, and covers the next use site that
      # nobody thinks to audit.
      #
      # What the allowlist leaves out is a closed set, which is what makes it
      # auditable rather than a guess: the ASCII whitespace and punctuation, less the
      # six characters a Windows path needs, : / . - _ ~. So the refused characters
      # are space, tab, newline and the other ASCII control characters, delete, and
      # ! " # $ % & ' ( ) * + , ; < = > ? @ [ \ ] ^ ` { | }.
      #
      # Everything above ASCII is accepted, and tr strikes those bytes before the test
      # so that the test judges only what it was written to judge. An account named
      # "José" gives a TEMP with an é in it, and é is not a character make or the shell
      # reads as a separator or an operator, which is the whole criterion here. A
      # refusal is no longer free either: the roots below this one are tried and the
      # in-package tree is the fallback, which builds and warns rather than failing,
      # but that fallback is the path-limit bug this block exists to avoid, so every
      # needless refusal reinstates it. R draws its own line in the same place,
      # shortening a temporary directory with whitespace in it to the 8.3 form and
      # aborting when it cannot, and saying nothing about non-ASCII. LC_ALL=C is there
      # so that the high byte range means bytes under every tr this can meet.
      wtd_ascii="$(printf '%s' "${wtd_normalised}" | LC_ALL=C tr -d '\200-\377')"
      case "${wtd_ascii}" in
        *[!A-Za-z0-9_:/.~-]*)
          echo "skipping temporary directory '${wtd_root}': it holds a character that a make variable cannot carry safely"
          continue
          ;;
      esac

      wtd_candidate="${wtd_normalised%/}/rt-${wtd_tag}"

      # Take the candidate only when moving there helps. The comparison is against
      # the in-package length rather than the budget on purpose: when nothing on
      # offer is short enough, the shortest option is still the best one available.
      if [ "${#wtd_candidate}" -ge "${#wtd_in_package}" ]; then
        continue
      fi

      # mkdir without -p, so a root that does not exist is refused here rather than
      # created, and so an install takes the directory only when it is the one that
      # created it. Process ids are reused, so a tree an interrupted install left
      # behind could carry the name a later install picks, and inheriting a stale
      # cargo tree is what the per-install name exists to prevent.
      if wtd_mkdir_message="$(mkdir "${wtd_candidate}" 2>&1)"; then
        wtd_relocated="${wtd_candidate}"
        break
      fi

      # A dev name is chosen to be the name the last dev build of this checkout
      # chose, so a directory already standing there is that build tree and reusing
      # it is the point of the name. Writability is part of the test, because cargo
      # writes into this directory: adopting one it cannot write to would fail inside
      # the recipe, where the loop that could have tried the next root is long gone.
      # A refusal on writability says as much in its own words, because all mkdir
      # reported for this candidate is that the name exists, which is just as true of
      # the tree that does get adopted.
      if [ "${wtd_profile}" = "dev" ] && [ -d "${wtd_candidate}" ]; then
        if [ -w "${wtd_candidate}" ]; then
          wtd_relocated="${wtd_candidate}"
          echo "reusing the cargo target directory an earlier development build left at '${wtd_candidate}'"
          break
        fi
        echo "skipping cargo target directory '${wtd_candidate}': a development build tree is already there and it is not writable"
        continue
      fi

      # Say what mkdir said rather than discarding it. A root that is gone, a root
      # that refuses writes, and a name taken by something the branch above does not
      # adopt, meaning a name taken outside the dev profile or taken under it by
      # something that is not a directory, all end in the same fallback, and a log that
      # does not say which of them happened leaves no way to tell a leftover directory
      # apart from a misconfigured machine.
      echo "skipping cargo target directory '${wtd_candidate}': ${wtd_mkdir_message:-mkdir gave no reason}"
    done
  fi

  if [ -n "${wtd_relocated}" ]; then
    TARGET_DIR="${wtd_relocated}"
    echo "using cargo target directory: '${TARGET_DIR}' (${#TARGET_DIR} characters)"
    echo "moved out of the package tree, where it would be ${#wtd_in_package} characters, over the ${TARGET_DIR_MAX} that cmake leaves room for"
    if [ "${#TARGET_DIR}" -gt "${TARGET_DIR_MAX}" ]; then
      echo "this is still over ${TARGET_DIR_MAX}, so cmake may run past the Windows path limit"
    fi
  else
    TARGET_DIR='$(CURDIR)/rust/target'
    echo "using cargo target directory: '${wtd_in_package}' (${#wtd_in_package} characters)"
    if [ "${#wtd_in_package}" -gt "${TARGET_DIR_MAX}" ]; then
      echo "no shorter directory was available, so cmake may run past the Windows path limit"
    fi
  fi
}
