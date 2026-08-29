# Pure-state projector PfQMC Phase 1 report

## Scope and provenance

- Work directory: `/home/sunxr/pure-state-projector-phase1`
- Branch: `pure-state-projector-phase1`
- Base: `zero-t-pfqmc@89ebcbc431476dfb71bfecf6b8647a090a81bb3f`
- Tests-first commit: `a950778d1ac45f119feab14f5e1ea983b8260748`
- External build/results: `/home/sunxr/pure-state-projector-phase1-results`

This phase implements only pure Gaussian trial boundaries, pure-projector
equal-time Green construction, thin-QR stabilization, and matrix/Operator
subspace propagation helpers. It does not implement a walker, Monte Carlo
sweeps, absolute Pfaffian Z2 initialization, MP checkpoints, driven evolution,
or performance optimization.

No existing finite-trial-density projector/driven source or result was
modified. No main/develop or parallel issue/proposal-safety/MP worktree was
modified. No scheduler or production physics task was submitted.

## Tests-first record

The complete six-group test harness, build script, parallel runner, and README
were committed before the implementation. The recorded pre-implementation
build exited 4 at the missing `gaussian_trial_state.h` include, as expected.
The two minimal implementation headers were then added.

## GaussianTrialState

`inc/gaussian_trial_state.h` stores the `2N x N` occupied subspace `Phi`,
`G_T = I - 2 Phi Phi^dagger`, dimensions, validation residuals, finite status,
and zero-mode diagnostics.

Construction is supported from a user `Phi` and from a Hermitian,
transpose-skew Majorana Hamiltonian by selecting exactly N negative-energy
eigenvectors. The implementation rejects nonfinite values, invalid dimensions,
failed pure-state invariants, and any eigenvalue within the explicit zero-mode
tolerance. It never selects parity from floating-point noise.

The Kitaev helper calls the existing `SpinlessTvChainUtils::KineticGenerator`
and locates the two edge Majoranas through `majoranaCoord2Idx`. The optional
`epsilon` is added only to a local initial-Hamiltonian copy. It is not stored in
the trial state or used by subsequent propagation. For L=4 and L=6 OBC at
`t=delta=1, mu=0`, epsilon zero is rejected and epsilon `+/-1e-8` selects
opposite Pfaffian parity trial states.

## Pure projector Green

`inc/pure_projector_green.h` implements

`G = I - 2 Phi_R (Phi_L^dagger Phi_R)^(-1) Phi_L^dagger`

through a checked SVD solve, never an explicit matrix inverse. The stable path
first obtains independent rectangular thin QR bases and applies the same
checked overlap solve to `Q_L^dagger Q_R`. Structured results report status,
overlap rank/rcond, solve residual, and involution residual. Singular overlap,
rank-deficient subspaces, low rcond, nonfinite data, or excessive residuals
fail closed.

Matrix and existing `Operator` helpers implement ket left multiplication by B
and bra-subspace left multiplication by B-dagger. The current square UDT code
was not changed.

## Build and parallel test command

```bash
set +u
source /opt/ohpc/pub/apps/intel/oneapi/setvars.sh
set -u
export EIGEN3_INCLUDE_DIR=/home/sunxr/software/eigen-3.4.0
export PFQMC_PFAPACK_DIR=/home/sunxr/new-pfqmc-main/inc/pfapack
export PURE_PROJECTOR_BUILD_DIR=/home/sunxr/pure-state-projector-phase1-results/build
bash reproduction/pure_projector_kitaev/build.sh

export PURE_PROJECTOR_TEST_OUTPUT_DIR=/home/sunxr/pure-state-projector-phase1-results/test_release
export PURE_PROJECTOR_TEST_JOBS=6
bash reproduction/pure_projector_kitaev/run_tests.sh
```

Compilation uses C++17, `PFQMC_SCALE_SAFE_UDT`, Intel oneAPI/MKL, and the
existing PFAPACK libraries. Six independent test processes ran concurrently;
each process used one MKL/OMP thread.

## Results

All six groups pass:

1. Single-site empty/full analytic states: maximum Green error
   `3.1401849173675498e-16`; parities opposite.
2. Identity propagation: maximum direct/thin-QR error
   `6.5851099635466118e-16`; singular and untrusted-rcond cases fail closed.
3. Random legal Gaussian states: maximum invariant error
   `8.8998738271776818e-16`; malformed and nonfinite inputs rejected.
4. Finite-Lambda boundary: monotone convergence through Lambda 24; maximum
   Lambda-24 difference from the pure formula `7.5502493146799674e-11`.
5. Kitaev OBC L=4,6: both unsplit cases report zero-mode ambiguity; all four
   split states pass invariants; epsilon signs select opposite parity. Maximum
   invariant error `6.2803698347350997e-16`.
6. Moderate 24-factor propagation: direct versus periodic thin-QR Green error
   `2.694123373749838e-15`; maximum reported Green residual
   `4.4141059618745704e-15`.

Machine-readable results are in
`reproduction/pure_projector_kitaev/phase1_test_results.csv`; detailed JSON and
stderr files remain in the external results directory.

## Exact file changes

Added:

- `inc/gaussian_trial_state.h`
- `inc/pure_projector_green.h`
- `reproduction/pure_projector_kitaev/README.md`
- `reproduction/pure_projector_kitaev/build.sh`
- `reproduction/pure_projector_kitaev/run_tests.sh`
- `reproduction/pure_projector_kitaev/phase1_core_test.cpp`
- `reproduction/pure_projector_kitaev/phase1_test_results.csv`
- `PURE_PROJECTOR_PHASE1_REPORT.md`

Modified existing production/core files: none.

Phase 1 stops here. The next possible phase is absolute Pfaffian Z2
initialization and an open-boundary static walker; it has not been started.
