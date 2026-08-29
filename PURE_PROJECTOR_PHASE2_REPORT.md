# Pure-state projector PfQMC Phase 2 report

## Scope and provenance

- Branch: `pure-state-projector-phase2`
- Phase 1 base: `a36704b186e7d6231ee4eb820e0917fc54fbf3da`
- Validated source commit: `4f9c8b7f8b27c8bab3e0a7554d81344d76024c9f`
- Validation executable SHA256: `e41f3e7aced915785d8cc736cc5b8ccd697fc80ae0b326b61cb54c6507d359c7`
- TDD red commit: `a650424858a1997b26e324d4418ef80ef5a764bf`
- Implementation commit: `653bd1b16835dadd31de051bf6806a2d796c121b`
- No main/develop merge, no proposal-safety/MP/issue branch merge, no production scan.
- Existing finite-trial-density walker, RNG, SpinlessVOperator, real-Z2, MP oracle and Phase 1 results were not modified.

## Implementation

`PureProjectorWeightEvaluator` starts at `W0=1`, `G0=G_T`, applies slices in their recorded ket action order, uses the existing Pfaffian product routine to fix each square-root branch, and independently propagates the rectangular occupied subspace. It periodically performs thin QR and reconstructs Green with the checked Phase 1 solve. It never calls the closed-contour raw-sign initializer or an MP checkpoint.

The result is fail-closed and includes log absolute weight, complex phase, integer Z2, Green, overlap rank/rcond/residual, Green residual, Pfaffian status, first failing slice, determinant identity error and structured status. `RealZ2` rejects nonfinite, clearly complex, near-zero and untrusted ratios; zero weight returns `z2_sign=0` and `log_abs_weight=-inf`, never NaN or fabricated `+1`. `GenericComplex` retains the full phase. No explicit matrix inverse was added.

`PureStaticProjectorContour` stores independent ket and bra protocols, an explicit middle measurement cut, flattened action-order labels, trial parity and edge splitting. Bra construction uses the strict adjoint and reversed protocol/factor order. There is no `beta_trial`.

`PureProjectorSlowWalker` copies the full candidate configuration, preserves proposal identity and caller-supplied uniform, rebuilds the candidate weight from scratch, makes the Metropolis decision from the log-weight difference, and installs candidate HS/Green/Z2 only after acceptance. It contains no fast update or incremental sweep.

The validation CSV gate performs `flush -> stream check -> close -> fail check` on all four requested streams before printing a complete JSON record.

## Convention issue found and fixed

The noncommuting L=4 oracle exposed the contragredient Majorana spinor convention: a single-particle action chain `B_n ... B_1` corresponds to the reverse dense-spinor factor order `F_1 ... F_n`. L=2 could not expose this because its short bond layers partially commute. The dense oracle and strict bra protocol now encode the verified reverse order. No parity-sign mismatch was found. OBC `t=delta, mu=0` still requires explicit edge splitting/parity; L=6 ED used `epsilon=+1e-8`, parity `+1`.

## Tests

All 8 requested groups passed:

1. Identity contour: `W=1`, Z2 `+1`, `G=G_T`.
2. Trial-state gauge invariance.
3. L=2 full dense-Fock enumeration for OBC/PBC and hs0/hs1.
4. L=4 dense-Fock random configurations: 100 for each OBC/PBC and hs0/hs1 (400 total).
5. Noncommuting ket/bra contour-order test, including sensitivity to deliberately wrong order.
6. Exact L=4 OBC hs1 enumeration over 256 configurations.
7. Four-seed slow MC comparison and same-seed exact reproducibility.
8. Same-Trotter-contour ED comparison: even-L PBC L=4 and explicitly split OBC L=6.

Numerical maxima:

- dense-Fock weight relative error: `8.2295042082128555e-15`
- dense-Fock Green relative error: `7.6138644475598846e-15`
- determinant identity relative error: `1.2804623963124567e-15`
- MC average-sign absolute deviation: `2.2204460492503131e-16`
- MC `S_pi` absolute deviation: `2.8212657952945697e-3`
- MC `R_CDW` absolute deviation: `6.7573621511785031e-3`
- L=4 ED (`S_pi`, `S_pi_dq`, `R_CDW`) differences: `5.5910e-17`, `6.5190e-17`, `1.1005e-16`
- L=6 ED (`S_pi`, `S_pi_dq`, `R_CDW`) differences: `1.1111e-16`, `5.5528e-17`, `6.6622e-16`
- slow-MC minimum overlap rcond: `0.676736`
- slow-MC maximum overlap solve residual: `3.85474e-15`
- zero/untrusted MC proposals: `0`; Pfaffian status: `success`

## Build and run

```bash
source /opt/ohpc/pub/apps/intel/oneapi/setvars.sh
EIGEN3_INCLUDE_DIR=/home/sunxr/software/eigen-3.4.0 \
PFQMC_PFAPACK_DIR=/home/sunxr/new-pfqmc-main/inc/pfapack \
reproduction/pure_projector_kitaev/build.sh
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
reproduction/pure_projector_kitaev/build/phase2_core_test
PURE_PROJECTOR_PHASE2_THREADS=8 \
reproduction/pure_projector_kitaev/run_validation.sh \
reproduction/pure_projector_kitaev/phase2_results
```

## Added or modified files

- `inc/pure_projector_weight.h` (new)
- `inc/pure_projector_static.h` (new)
- `reproduction/pure_projector_kitaev/phase2_core_test.cpp` (new)
- `reproduction/pure_projector_kitaev/phase2_validation.cpp` (new)
- `reproduction/pure_projector_kitaev/run_validation.sh` (new)
- `reproduction/pure_projector_kitaev/build.sh` (modified)
- `reproduction/pure_projector_kitaev/README.md` (modified)
- `dense_configuration_checks.csv` (generated)
- `exact_enumeration_summary.csv` (generated)
- `slow_mc_vs_exact.csv` (generated)
- `pure_projector_ed_comparison.csv` (generated)
- `PURE_PROJECTOR_PHASE2_REPORT.md` (generated)

## Phase 3 gate

Phase 2 provides a bounded slow correctness oracle suitable for checking a later open-boundary forward/backward sweep and fast update. Entry to Phase 3 is technically supported, but this work intentionally does not implement or start Phase 3.
