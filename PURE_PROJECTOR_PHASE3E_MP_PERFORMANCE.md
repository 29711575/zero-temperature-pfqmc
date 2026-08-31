# Pure-state projector Phase 3E: same-proposal MP performance

## Provenance and scope

- Branch: `pure-state-projector-phase3e-mp-performance`
- Base: `081d79d67cbb51f165e7053737f4bb78d0c1a365`
- Implementation commit: `847da9749c0a95c84fd2cbfa9e621b204ae2edae`
- Production executable SHA256: `df8ff62378d8cf50efad0bad72e872daeb0a7b6a93869b83a45209fbe73ad929`
- C++17 with `PFQMC_SCALE_SAFE_UDT`; one CPU per benchmark process.
- No Phase 4 work, high-statistics scan, trust-threshold change, fallback suppression, condition-aware ratio, left recovery, or mutating raw/MP checkpoint was introduced.
- Commit `222a89f0c0f42db38625895bc0d7653a0090f325` was not merged or cherry-picked.

The optimized path is restricted to the production `RealZ2` MP oracle. `GenericComplex` explicitly remains on the legacy evaluator and has a regression asserting identical result and zero cache use.

## Repetition found in the old fallback

The old evaluator already shared the proposal-independent left/right propagation between the pre and post states *within one precision*. The expensive repetition that remained was:

1. every occurrence of a repeated contour factor was converted from double to MP independently;
2. sparse local HS factors were applied as dense MP matrix products;
3. the same canonical double input was rediscovered independently by the 160- and 320-digit stages;
4. the complete MP thin-QR propagation was repeated at each precision.

The accepted/rejected endpoint recovery still constructs the normal double stack after the decision. It was deliberately not bypassed because it is a fail-closed state-restoration gate, not merely diagnostic work.

## Implementation

For each frozen same-proposal fallback, the code now builds one proposal-scoped canonical input:

- canonical action order is preserved as an ordered vector of operator IDs;
- operators are deduplicated only by exact dimensions and exact double matrix entries (hash collisions are checked by exact comparison);
- each precision converts each unique operator once;
- an operator represented exactly as `I + delta` uses sparse row updates when fewer than half its entries differ from identity;
- 160 and 320 digits remain independent numerical evaluations and must satisfy the unchanged consecutive-precision agreement gate; 640 digits is still used on disagreement;
- the cache is destroyed at return from every fallback, so accepted/rejected mutation cannot leave stale MP state.

The same proposal identity and uniform are untouched. All trust thresholds, endpoint solves, Pfaffian branch alignment, real-axis test, zero test, Z2 transport, acceptance rule, and fail-closed exits are unchanged.

New read-only diagnostics report canonical builds/invalidations, exact operator reuse, sparse/dense applications, 160/320/640 timings, conversion, propagation, thin QR, endpoint solve, and local Pfaffian time.

## Frozen Task 4/8

| Task | Ratio | Decision | Final hash | Fallbacks | Failure | Frozen-reference speedup |
|---|---:|---|---:|---:|---:|---:|
| 4 | 0.70741914065425615 | accept | 949208207496548183 | 12 | 0 | 1.449x |
| 8 | 0.0026573889530685144 | reject | 14005514076608941581 | 69 | 0 | 1.474x |

Both ratios, imaginary parts, endpoint rconds, decisions, and final hashes are identical to the base executable. Each replay continued for 64 proposals past the historical failure.

## Controlled true-sweep benchmark

Parameters were `L=6`, `theta=6`, PBC+hs0, one burn sweep plus two measurement sweeps, 720 HS variables and 2160 proposals. The before/after comparison used the same four-process local batch, one thread per process. These runs are timing tests, not physics data.

| V | Fallbacks before/after | Mean fallback before | Mean fallback after | Fallback speedup | Seconds/sweep before | Seconds/sweep after | Sweep speedup |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 285 / 285 | 1.517343 s | 1.054580 s | 1.439x | 147.3996 | 103.3920 | 1.426x |
| 6 | 756 / 756 | 1.504963 s | 1.052246 s | 1.430x | 384.8481 | 270.6604 | 1.422x |

For both V values, acceptance, signed numerators, observables/status, final HS hash, final RNG hash, Z2, and fallback count are identical. Fallback and slow-reference failures are zero. Peak RSS remained about 31 MiB.

The exact-operator cache hit rate is 98.6472% for both runs: across the 160/320 pair, each fallback requested 1922 operators but converted only 26 unique MP matrices (961 requests and 13 conversions per precision). All contour applications in these runs used the exact sparse-delta path.

## MP stage profile

A second two-process run added the thin-QR timer. It reproduced the same hashes and counts. The main MP fractions were:

| V | MP total | thin QR | propagation | endpoint | local Pfaffian | conversion |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 293.221 s | 68.43% | 22.76% | 3.54% | 4.90% | 0.15% |
| 6 | 761.748 s | 67.64% | 23.64% | 3.50% | 4.85% | 0.15% |

The maximum endpoint Green residual was `8.4404869330236453e-09` (V=6); V=4 was `2.5813797673489134e-14`. No endpoint or MP fallback failed.

## Safety and regression results

- Phase 1: 6/6 PASS.
- Phase 2 core: 5/5 PASS; dense/exact/ED validation: 8/8 PASS.
- Phase 3A core: 4/4 PASS; lockstep validation: 8/8 PASS, Z2 mismatch 0 and trajectory mismatch 0.
- Phase 3C core: 9/9 PASS. L12 mirrored tasks 16/17/20/21 each passed 128/128 initializations.
- Production-semantics core: PASS; Phase 3B driver/I/O contract: 4/4 PASS; PBC structural guard/energy: PASS.
- Task 4/8 frozen replay: PASS with the exact expected ratios and decisions.
- Checkpoint/restart: Task 4 and Task 8 continuous/restarted measurement CSVs are byte-identical. Their SHA256 values are respectively `55f1f2d421bf91e109a97c2d44cbfdf0bc5f5236c3458926cf5612dcd5bb8e12` and `c934271b96f77485d5ae4925c47bac9d0cfc7c3a1ac993ecd48494a59d88cc8d`.
- Ratio reciprocity/detailed balance: 1002/1002 rows PASS; maximum reciprocity error `6.421527630581721e-15`; decision failures 0; malformed reference remains fatal.
- OBC+hs1 L=6 short true-sweep: average Z2 `+1`, fallback failures 0, endpoint residual `1.4934400759935883e-14`.
- RealZ2 cache-vs-legacy and GenericComplex-legacy-isolation tests: PASS.

## Build and test commands

On the HPC/toolchain environment:

```bash
source /opt/intel/oneapi/setvars.sh --force
export EIGEN3_INCLUDE_DIR=/path/to/eigen3
export BOOST_INCLUDE_DIR=/path/to/boost
export PFQMC_PFAPACK_DIR=/path/to/pfapack
export PURE_PROJECTOR_PHASE3B_BUILD_DIR=/tmp/pure_phase3e_build
bash reproduction/pure_projector_kitaev/build_phase3b.sh

/tmp/pure_phase3e_build/phase3e_mp_performance_test
/tmp/pure_phase3e_build/phase3e_safety_regression ratio /tmp/ratio_checks
/tmp/pure_phase3e_build/phase3e_safety_regression restart /tmp/restart_checks
```

## Conclusion

The optimization is worth retaining: it reduces controlled true-sweep wall time by about 29.7% without changing any fallback trigger, trajectory, Z2, observable, or trust outcome. The next dominant cost is now clearly MP thin QR (about 68% of fallback time), followed by MP propagation (about 23%).

A checkpoint-based sweep optimization is therefore justified as the next isolated performance step, but it should cache high-precision proposal-independent subspaces keyed by configuration hash/cut and prove precise invalidation after accepted updates. It must retain the current legacy same-proposal oracle as an audit and must not reuse untrusted double live-stack state as an MP trust input.
