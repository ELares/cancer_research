# BioFVM external cross-check of the reaction-diffusion solver (#408)

This directory holds the BioFVM driver used to cross-check the Rust
`reaction_diffusion` steady-state supply field against an INDEPENDENT PDE solver
(BioFVM, the diffusion microenvironment of PhysiCell). It is a LOCAL validation
tool, not part of CI: CI reads only the committed result
(`analysis/calibration/rd-biofvm-crosscheck.{json,md}`), produced by
`scripts/validate_rd_vs_biofvm.py`.

`biofvm_rd_check.cpp` solves the same matched scenario as the Rust solver
(steady-state `D grad^2 c - k c = 0`, lambda = sqrt(D/k), Dirichlet vessel sources
clamped to c = 1, no-flux boundaries, uniform decay) and writes the field to
`biofvm_field.csv`.

## Build (one-time, local)

BioFVM is not vendored here (it has its own license). Fetch it and compile the
driver against it. On macOS the toolchain needs two workarounds, both
self-contained (no system installs):

1. **No OpenMP**: Apple clang ships no `libomp`. BioFVM calls no `omp_*` functions
   (only `#pragma omp`, which clang ignores serially), so an empty `omp.h` stub on
   the include path is enough. The driver runs single-threaded (fine at validation
   grid sizes).
2. **libc++ path**: Apple clang does not auto-resolve the SDK's C++ headers here,
   so pass `-isysroot $SDK -isystem $SDK/usr/include/c++/v1` explicitly.

```bash
git clone --depth 1 https://github.com/MathCancer/BioFVM
mkdir -p ompstub && printf '#ifndef _OMP_STUB_H\n#define _OMP_STUB_H\n#endif\n' > ompstub/omp.h
SDK=$(xcrun --show-sdk-path)   # e.g. /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
clang++ -O2 -std=c++11 -w -isysroot "$SDK" -isystem "$SDK/usr/include/c++/v1" \
  -I ompstub -I BioFVM \
  BioFVM/BioFVM_*.cpp BioFVM/pugixml.cpp biofvm_rd_check.cpp -o biofvm_rd_check
```

(On Linux with g++ and OpenMP, just `g++ -O2 -std=c++11 -fopenmp -I BioFVM
BioFVM/BioFVM_*.cpp BioFVM/pugixml.cpp biofvm_rd_check.cpp -o biofvm_rd_check`.)

## Run via the comparison harness

```bash
python3 scripts/validate_rd_vs_biofvm.py --biofvm /path/to/biofvm_rd_check
```

This runs the Rust example (`cargo run --release --example rd_field_dump` in
`simulations/ferroptosis-core`) and the BioFVM driver on matched single- and
two-source scenarios, sweeps BioFVM's time step, and writes the committed result.

## Result

The two independent solvers agree: identical field shape (log-field Pearson
r > 0.99) and effective decay length, with a magnitude offset that is BioFVM's LOD
operator-splitting error and converges to 1 as its time step dt -> 0. See
`analysis/calibration/rd-biofvm-crosscheck.md`.

## A note on the driver

`biofvm_rd_check.cpp` uses canonical BioFVM (`set_density`, `resize_space_uniform`,
assign `diffusion_decay_solver`, `simulate_diffusion_decay`). Dirichlet sources are
re-clamped to 1.0 each step (robust, independent of BioFVM's dirichlet-activation
machinery). Use a small `dt` so the LOD operator-splitting error is small (the
harness sweeps `dt` and shows the convergence).
