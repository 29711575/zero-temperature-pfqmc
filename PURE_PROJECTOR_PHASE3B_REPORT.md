# Pure-state projector PfQMC Phase 3B report

## Scope and provenance

- Branch: `pure-state-projector-phase3b`
- Base: `026900c9e79515e090adf3b691bbb9a7fa5094e2`
- Validated code commits on HPC: `2139bc60b17061a5b627f2cf499fdbe30140519e`, followed by the explicit parity-policy correction `1c541430c24441efe75c500776d0cb6c91c591df`
- Production driver SHA256: `1e5ddb031186163bb81f05d4c0a5a869f1f802004438ad4b356573f76e74e555`
- Projector type: `pure_state`
- `condition_aware_ratio=false`, `left_recovery=false`, and mutating raw/MP checkpoints are absent/disabled.
- No driven protocol, finite-trial-density walker integration, larger-L scan, or Phase 4 work was added.

The implementation started with a failing production contract test (the driver did not exist), then added the minimum driver and made the contract/I/O test pass. Subsequent failures found by long-contour tests were handled fail-closed and fixed before the final regression run.

## Implementation

The new static driver builds independent ket and bra HS branches, reverses the full bra protocol including noncommuting bond-factor order, and records the flattened action locations. The active stack follows a deterministic forward/backward proposal schedule. Both boundaries are the rectangular `GaussianTrialState::Phi`; stable center measurements and audit references rebuild independently from the two boundaries and do not move or mutate the live sweep cut.

Accepted proposals transactionally replace only the selected factor at the current cut. Thin-QR stabilization and scheduled complete two-sided rebuilds re-anchor the subspaces. The production reference evaluator tracks QR factors and overlap determinants in log form, while retaining the per-factor Pfaffian phase/Z2 rule. It does not form a raw long product and does not apply MP or raw-sign corrections.

Fast proposals check finite/reality/zero/overlap/decision-margin conditions and the updated Green involution residual before mutation. A periodic read-only audit compares the same candidate with the stable full reference. An untrusted value is never used after an alarm; reference failure terminates the run.

The driver supports the requested physical, trial, MC, stabilization, retained-stream, and provenance parameters. It rejects odd-L PBC and mismatched/ambiguous trial parity. The explicit trial parity convention is L-dependent for the tested gapped and edge-split trials: L=6/18 use `+1`, while L=4/12 use `-1`. The existing PBC energy implementation includes the `(L-1,0)` quadratic and interaction boundary bonds; arbitrary `t` is included by a linear decomposition of the existing kinetic generator.

Center-cut output contains average integer Z2 and bin error, signed numerators and sign denominator, `S_pi`, `S_pi_dq`, `R_CDW`, energy, fermion parity, acceptance, complex-energy diagnostics, overlap/stabilization/trust counters, parameters, source/executable provenance, and final HS/RNG hashes. If the sign denominator is unresolved, reweighted values are JSON `null` with an explicit status.

The final `status=complete` record is emitted only after every requested retained CSV write succeeds and the stream passes flush, check, close, and post-close check. `/dev/full` and an invalid output path return nonzero, write an explicit stderr error, and never emit complete.

## Build and test commands

```bash
source /opt/ohpc/pub/apps/intel/oneapi/setvars.sh
export LC_ALL=C LANG=C
EIGEN3_INCLUDE_DIR=/home/sunxr/software/eigen-3.4.0 \
PFQMC_PFAPACK_DIR=/home/sunxr/new-pfqmc-main/inc/pfapack \
reproduction/pure_projector_kitaev/build_phase3b.sh

OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
reproduction/pure_projector_kitaev/build_phase3b/phase3a_core_test

python3 reproduction/pure_projector_kitaev/phase3b_contract_test.py \
  reproduction/pure_projector_kitaev/build_phase3b/pure_projector_driver

python3 reproduction/pure_projector_kitaev/run_phase3b_regression.py core \
  reproduction/pure_projector_kitaev/build_phase3b/pure_projector_driver \
  reproduction/pure_projector_kitaev/phase3b_regression

python3 reproduction/pure_projector_kitaev/run_phase3b_regression.py ed \
  reproduction/pure_projector_kitaev/build_phase3b/pure_projector_driver \
  reproduction/pure_projector_kitaev/phase3b_regression \
  --ed-executable reproduction/pure_projector_kitaev/build_phase3b/pure_projector_ed
```

## Regression results

All completed automated regression groups passed: Phase 1 `6/6`; Phase 2 core `5/5` and validation `8/8`; Phase 3A core `4/4` and validation `8/8`; Phase 3B contract/I/O contract `4/4`; production fast/slow trajectories `2/2`; block sizes `4/4`; I/O cases `4/4`; and ED/theta rows `8/8`. Together with the L=12 and L=18 sign/stability smokes, this is 55 successful checks.

Final Phase 3A lockstep metrics after the Phase 3B changes were ratio error `1.6007632792608383e-15`, Green error `1.0332568796446633e-14`, Z2 mismatches `0`, trajectory mismatches `0`, and OBC/hs1 always `+1`. The final production L=2/4 audit has maximum ratio error `5.5190823337346546e-15`, maximum Green error `3.672157402078089e-14`, and zero trajectory mismatches. Block sizes 1, 2, 4, and 8 produced identical HS hashes and trajectories.

At theta 12, the L=4 PBC fixed-parity comparison differs from ED by `0.0026739899155541202` in energy, zero in parity, `5.55e-17` in `S_pi`, and `2.22e-16` in `R_CDW`. The L=6 OBC comparison differs by `0.05916841652300242` in energy, zero in parity, `0.0026428922599828764` in `S_pi`, and `0.049232082360095575` in `R_CDW`. The endpoint improves over theta 2 in energy for both sizes, and the L=6 endpoint `S_pi` is close to ED, but the finite-statistics sequence is not monotone (notably L=6 theta 8). This is evidence of approach to the fixed-parity ground state, not a high-statistics convergence claim.

OBC/hs1 physical Z2 was `+1` at every measurement in the L=6 regression and the L=12/L=18 long-contour smokes. Their maximum observed Green rebuild error was `5.5581324113188364e-15`; all six smoke slow references succeeded and none failed. The completed short L=6 PBC/hs0 fast/slow regression had preliminary average sign `1.0` for all four seeds.

The bounded L=18 timing smoke (30 proposals including center measurements and three audits) took 41.40 s and used 635348 KiB maximum RSS. The inherited Phase 3A correctness timing gave fast/slow `0.4010`; neither number is a production throughput claim.

## HPC benchmark status

The corrected 36-task array is queued as PBS `153908[]` in its independent result directory. It uses exactly L=6/12/18, theta=L, V=2/4/6, PBC+hs0 and OBC+hs1, two seeds, burn 1000, measurements 2000, dt=0.1, one CPU, 32 GiB, and 72 h. The earlier `153907[]` was cancelled while still queued after the L=12 explicit-parity convention was caught; it produced no benchmark data. Job/task `153850` and its worktree/results were not touched.

Per user direction, this report does not wait for `153908[]`. `production_benchmark.csv` therefore contains all 36 manifest rows with explicit pending/missing result fields; no queued task is represented as completed. The aggregation script will replace those fields after the array finishes.

## Files

Modified:

- `inc/pure_projector_fast.h`
- `inc/pure_projector_stack.h`

Added:

- `reproduction/pure_projector_kitaev/pure_projector_driver.cpp`
- `reproduction/pure_projector_kitaev/pure_projector_ed.cpp`
- `reproduction/pure_projector_kitaev/build_phase3b.sh`
- `reproduction/pure_projector_kitaev/phase3b_contract_test.py`
- `reproduction/pure_projector_kitaev/run_phase3b_regression.py`
- `reproduction/pure_projector_kitaev/prepare_phase3b_benchmark.py`
- `reproduction/pure_projector_kitaev/phase3b_benchmark.pbs`
- `reproduction/pure_projector_kitaev/aggregate_phase3b_benchmark.py`
- `ed_theta_convergence.csv`
- `production_benchmark.csv`
- `block_size_checks.csv`
- `proposal_audit.csv`
- `io_regression.csv`
- `PURE_PROJECTOR_PHASE3B_REPORT.md`

## Readiness

The static driver, safety gate, I/O gate, long-contour smokes, and inherited regressions are ready for the queued benchmark. Phase 4 should **not** start yet: the requested 36-task benchmark is pending, and the theta convergence series would benefit from higher effective sweep statistics before using it as a driven-protocol baseline.
