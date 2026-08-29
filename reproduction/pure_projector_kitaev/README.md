# Pure-state projector PfQMC: Phase 1 core

This directory contains only the pure Gaussian boundary and equal-time Green
core. Phase 1 deliberately excludes Monte Carlo sweeps, absolute Pfaffian Z2
initialization, MP checkpoints, driven evolution, and a complete walker.

The six tests cover analytic one-site states, identity propagation, random pure
Gaussian constraints, finite-boundary convergence, Kitaev OBC zero-mode/parity
policy, and direct versus periodically thin-QR-stabilized propagation.

Build with the existing Intel oneAPI/MKL and PFAPACK environment:

```bash
export EIGEN3_INCLUDE_DIR=/home/sunxr/software/eigen-3.4.0
export PFQMC_PFAPACK_DIR=/home/sunxr/new-pfqmc-main/inc/pfapack
export PURE_PROJECTOR_BUILD_DIR=/path/to/external/build
bash reproduction/pure_projector_kitaev/build.sh
```

Run the independent groups concurrently:

```bash
export PURE_PROJECTOR_TEST_OUTPUT_DIR=/path/to/external/results
export PURE_PROJECTOR_TEST_JOBS=4
bash reproduction/pure_projector_kitaev/run_tests.sh
```
