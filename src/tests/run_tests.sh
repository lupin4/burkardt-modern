#!/usr/bin/env bash
# burkardt-modern test harness.
#
#   1. DIFFERENTIAL — the same driver is built against the ORIGINAL objects and
#      against the MODERN objects; the two stdouts must be identical. Catches
#      any numeric drift introduced by the modernization.
#   2. ANALYTIC — the driver also checks every value against a closed form and
#      exits non-zero if one is violated. A diff alone cannot catch a bug the
#      two trees SHARE; this layer can.
#
# Both must pass.
set -euo pipefail
cd "$(dirname "$0")/../.."

FC=${FC:-gfortran}

# --- Windows/MSYS2 preflight -------------------------------------------------
# gfortran -cpp needs a writable temp dir. On this platform TMPDIR can resolve
# to C:\WINDOWS\, which is not writable, and the compile dies with
# "Cannot create temporary file ... Permission denied" (exit 127). mkdir.exe
# also lives only in /usr/bin, which is not always on PATH.
case "$(uname -s)" in
  MINGW*|MSYS*|CYGWIN*)
    export PATH="/c/msys64/ucrt64/bin:/c/msys64/usr/bin:$PATH"
    if [ -z "${TMPDIR:-}" ] || [ ! -w "${TMPDIR:-/nonexistent}" ]; then
      export TMPDIR="$(pwd)/build/tmp"
    fi
    ;;
esac
mkdir -p "${TMPDIR:-/tmp}" build/tests

echo "== building both trees =="
make original >/dev/null
make modern   >/dev/null

# Objects the driver needs. core/geometry.f90 owns every routine exercised
# here; add more objects to this list as the driver grows.
OBJ_REL="core/geometry.o"

build_drv () {                      # $1 = tree, $2 = std flag, $3 = output
  local objs=""
  for o in $OBJ_REL; do objs="$objs build/$1/$o"; done
  # shellcheck disable=SC2086
  $FC -O2 "$2" -fall-intrinsics -fno-underscoring -o "$3" src/tests/driver.f90 $objs
}

echo "== linking drivers =="
build_drv original -std=legacy build/tests/drv_orig
build_drv modern   -std=f2018  build/tests/drv_mod

# --- 1. analytic ------------------------------------------------------------
# Run each directly so a non-zero exit (a violated closed form) is caught even
# when the two trees agree with each other.
echo "== analytic check: ORIGINAL =="
orig_rc=0; build/tests/drv_orig > build/tests/out_orig.txt || orig_rc=$?
echo "== analytic check: MODERN =="
mod_rc=0;  build/tests/drv_mod  > build/tests/out_mod.txt  || mod_rc=$?

fail=0
if [ "$orig_rc" -ne 0 ]; then
  echo "ANALYTIC FAIL (original tree):"; grep -A1 EXPECTED build/tests/out_orig.txt || true
  fail=1
fi
if [ "$mod_rc" -ne 0 ]; then
  echo "ANALYTIC FAIL (modern tree):";   grep -A1 EXPECTED build/tests/out_mod.txt || true
  fail=1
fi
[ "$fail" -eq 0 ] && echo "ANALYTIC PASS: both trees satisfy every closed form"

# --- 2. differential --------------------------------------------------------
# Compare with shell builtins only. Depending on which bash Windows resolves
# (MSYS2 vs Git Bash) neither `diff` nor `cmp` is reliably on PATH, and a
# missing tool exits 127 -- which read as "the trees differ" and produced a
# false REGRESSION FAIL on byte-identical output. External tools are now used
# only to RENDER a delta, never to decide pass/fail.
if [ "$(cat build/tests/out_orig.txt)" = "$(cat build/tests/out_mod.txt)" ]; then
  echo "REGRESSION PASS: modern == original ($(grep -c ':' build/tests/out_mod.txt | tr -d ' ') printed values)"
else
  echo "REGRESSION FAIL: modern drifted from original"
  if command -v diff >/dev/null 2>&1; then
    diff -u build/tests/out_orig.txt build/tests/out_mod.txt || true
  else
    echo "  (diff unavailable; original vs modern side by side)"
    paste build/tests/out_orig.txt build/tests/out_mod.txt | awk -F'\t' '$1!=$2 {print "  orig: "$1"\n  mod : "$2}'
  fi
  fail=1
fi

exit "$fail"
