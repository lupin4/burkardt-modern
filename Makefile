# burkardt-modern — Build & Test
#
# Build original and modern versions, run regression tests
#
# Usage:
#   make all        Build everything
#   make original   Build original sources only
#   make modern     Build modern sources only
#   make test       Regression test: modern == original
#   make bench      Benchmark original vs modern

FC ?= gfortran

# Original: compile as-is (Fortran 90 compat)
FFLAGS_ORIG = -O2 -fPIC -std=legacy -fall-intrinsics -fno-underscoring -cpp

# Modern: Fortran 2018 + parallelism
FFLAGS_MOD = -O2 -march=native -fPIC -std=f2018 -fall-intrinsics \
             -fopenmp -fno-underscoring -cpp

ORIG_SRC = $(shell find src/original -name "*.f90" | sort)
MOD_SRC  = $(shell find src/modern -name "*.f90" 2>/dev/null | sort)

ORIG_OBJ = $(patsubst src/original/%.f90,build/original/%.o,$(ORIG_SRC))
MOD_OBJ  = $(patsubst src/modern/%.f90,build/modern/%.o,$(MOD_SRC))

.PHONY: all original modern test bench clean

all: original modern

original: $(ORIG_OBJ)
	@echo "=== Original: $(words $(ORIG_OBJ)) objects built ==="

modern: $(MOD_OBJ)
	@echo "=== Modern: $(words $(MOD_OBJ)) objects built ==="

build/original/%.o: src/original/%.f90
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS_ORIG) -c $< -o $@

build/modern/%.o: src/modern/%.f90
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS_MOD) -c $< -o $@

test:
	@echo "Running regression tests..."
	@cd src/tests && bash run_tests.sh

bench:
	@echo "Running benchmarks..."
	@cd benchmarks && bash run_bench.sh

clean:
	rm -rf build/
