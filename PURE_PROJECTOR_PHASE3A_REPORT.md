# Pure-state projector PfQMC Phase 3A report

## Provenance and scope

- Branch: `pure-state-projector-phase3a`
- Phase 2 base: `f280ff470a0f74a7d35cf2f9e61e787ae2ec2b50`
- Validated implementation commit: `bb5ad68f9bbfeb0b1ed433ed069674ec4b2384a4`
- Validation executable SHA256: `0c8020958e7da225ce44b1007af6be14f8d005ab708dee2ece1b79ebac916c65`
- TDD red commit: `61415528f84d9834c13422de9ee30a0103034311`
- No main/develop modification or merge. No finite-trial-density walker, MP checkpoint, condition-aware v3, left recovery, driven protocol, production scan, or L>6 run was added.
- The Phase 2 slow walker and Phase 2 evaluator were not modified.

## Implementation

`PureProjectorStackManager` maintains the two open pure-state boundaries,
`Phi_R=B_< Phi_T` and `Phi_L=B_>^dagger Phi_T`. It stores configurable
thin-QR checkpoints, moves a cut factor by factor in the recorded noncommuting
action order, and uses checked factor solves for reverse movement. At stable
cuts it calls the Phase 1 checked formula
`I-2 Q_R (Q_L^dagger Q_R)^-1 Q_L^dagger`; an untrusted overlap fails closed.

`PureProjectorFastWalker` supports `AuditLockstep` and `FastStrict`. A proposal
snapshot records branch, slice, factor, bond, auxiliary index, old/new HS,
configuration hash, fast ratio, and the exact caller-supplied uniform. The
decision path is mutation-free until acceptance.

For a local factor replacement the implementation forms `Delta B Phi_R`,
extracts its numerical rank with a thin SVD, evaluates the determinant lemma
in the resulting rank-2/rank-4 space, and uses the existing Pfaffian product
formula (including the established `(-1)^N` convention) to fix the square-root
and Z2 branch. The accepted Green is obtained with the corresponding Woodbury
update. No explicit matrix inverse or raw/MP full-contour sign correction is
used.

`AuditLockstep` evaluates the Phase 2 full-rebuild reference on every proposal.
`FastStrict` uses the persistent sweep stack and calls the same reference only
for a nonfinite, complex, near-zero, untrusted, or decision-margin alarm. A
failed reference terminates rather than using an untrusted fast value. Absolute
Z2 comes from the Phase 2 initializer and is thereafter transported only by
accepted trusted ratio signs.

## Tests

All 8 required groups passed:

1. Cut round trip: forward to the right boundary and backward to the left restores Green, both subspaces, cut, and configuration hash.
2. Stable rebuild: block sizes 1, 2, and 4 agree with a Phase 2-style full rebuild.
3. Noncommuting left/right orientation: the correct order agrees and an intentionally swapped order is detected.
4. Single proposals: 500 each for L=2 and L=4 using real Kitaev local HS factors; fast ratio, Woodbury Green, full rebuild, and Z2 agree.
5. Directed decisions: positive/negative ratios, rejection, decision-margin fallback, zero/reference failure, complex/untrusted fallback, same uniform, and mutation ordering.
6. Full trajectory lockstep: 1000 proposals each for L=2 and L=4; accept/reject, HS hash, Z2, Green, and observables agree at every step.
7. L=6 static regression: PBC/hs0 and OBC/hs1, four seeds each; 200 proposals per seed. Fast and AuditLockstep averages and acceptance agree exactly; OBC/hs1 remains physical Z2 `+1` throughout.
8. Stabilization stress: 48 moderate noncommuting factors, block sizes 1, 2, 4, and 8; no alarm.

Final maxima and counts:

- fast/slow ratio relative error: `1.5646998427194181e-15`
- fast-updated/full-rebuild Green relative error: `1.0278848524363269e-14`
- stabilization Green relative error: `5.7994357691308648e-15`
- Z2 mismatches: `0`
- trajectory mismatches: `0`
- first stabilization alarm: none (`-1`)
- OBC/hs1 physical Z2 always `+1`: yes

## Timing

Across the timed proposal checks and L=6 regressions (2600 proposals):

- persistent-stack fast path: `1.194136454 s` total, approximately `4.5928e-4 s/proposal`
- Phase 2 full rebuild reference: `2.716115772 s` total, approximately `1.0447e-3 s/proposal`
- fast/slow ratio: `0.4396485843`

This is a bounded correctness measurement only; no mandatory speedup threshold is imposed.

## Build and run

```bash
source /opt/ohpc/pub/apps/intel/oneapi/setvars.sh
EIGEN3_INCLUDE_DIR=/home/sunxr/software/eigen-3.4.0 \
PFQMC_PFAPACK_DIR=/home/sunxr/new-pfqmc-main/inc/pfapack \
reproduction/pure_projector_kitaev/build_phase3a.sh

OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
reproduction/pure_projector_kitaev/build_phase3a/phase3a_core_test

PURE_PROJECTOR_PHASE3A_THREADS=4 \
reproduction/pure_projector_kitaev/run_phase3a_validation.sh \
reproduction/pure_projector_kitaev/phase3a_results
```

The validation prints `status=complete` only after all four requested CSV
streams successfully flush, pass stream checks, close, and pass close checks.

## Added or modified files

- `inc/pure_projector_stack.h` (new)
- `inc/pure_projector_fast.h` (new)
- `reproduction/pure_projector_kitaev/phase3a_core_test.cpp` (new)
- `reproduction/pure_projector_kitaev/phase3a_validation.cpp` (new)
- `reproduction/pure_projector_kitaev/build_phase3a.sh` (new)
- `reproduction/pure_projector_kitaev/run_phase3a_validation.sh` (new)
- `proposal_ratio_checks.csv` (generated)
- `trajectory_lockstep.csv` (generated)
- `block_stabilization.csv` (generated)
- `fast_vs_slow_mc.csv` (generated)
- `PURE_PROJECTOR_PHASE3A_REPORT.md` (generated)

## Phase 3B gate

Phase 3A supplies a clean, proposal-lockstep correctness baseline for a later
production-grade long-contour stack, measurement stream, and larger-L
benchmark. Entry to Phase 3B is technically supported, but Phase 3B was not
started.
