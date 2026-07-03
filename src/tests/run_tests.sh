#!/usr/bin/env bash
# Regression: the MODERN tree must reproduce the ORIGINAL tree's numbers.
# Builds the same driver against both object sets and diffs stdout.
set -euo pipefail
cd "$(dirname "$0")/../.."
FC=${FC:-gfortran}
make original >/dev/null
make modern >/dev/null
mkdir -p build/tests
$FC -O2 -std=legacy -fno-underscoring -o build/tests/drv_orig src/tests/driver.f90 build/original/core/geometry.o
$FC -O2 -std=f2018 -fall-intrinsics -fno-underscoring -o build/tests/drv_mod src/tests/driver.f90 build/modern/core/geometry.o
build/tests/drv_orig > build/tests/out_orig.txt
build/tests/drv_mod  > build/tests/out_mod.txt
if diff -u build/tests/out_orig.txt build/tests/out_mod.txt; then
  echo "REGRESSION PASS: modern == original ($(wc -l < build/tests/out_mod.txt | tr -d ' ') checks)"
else
  echo "REGRESSION FAIL: modern drifted from original"; exit 1
fi
