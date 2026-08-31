# Pure-projector Phase 3F MP subspace performance

## Scope and provenance

- Branch: `pure-state-projector-phase3f-mp-subspace-performance`
- Base: Phase 3E `e17180d76af8f434add635b3d4d40e4c547b3c64`
- Implementation commit: `c1722dc4f932ab42f34f21a27d5966d03e0bd6f3`
- Validated production executable SHA256: `a3c421f854941b8631a191c0ffdbedb6e6ea29baafa21e11b521d5b5934194e2`
- C++17 and `PFQMC_SCALE_SAFE_UDT`; one CPU per benchmark process.
- No physics formula, trust threshold, precision-agreement rule, fallback trigger,
  Metropolis decision, RNG consumption, Z2 transport, or trajectory was changed.
  GenericComplex, condition-aware ratio, left recovery, and mutating MP/raw sign
  checkpoints remain unchanged/disabled as specified.

The Phase 3E audit showed that its operator cache avoids repeated double-to-MP
operator conversion, but every trusted fallback still propagates both endpoint
subspaces from the trial boundary and repeats thin QR independently at MP160 and
MP320. Those two stages dominated the remaining fallback cost.

## Minimal implementation

`PureMpSubspaceCache` keeps independent MP160, MP320, and MP640 thin-subspace
checkpoints. A checkpoint key contains the configuration hash, canonical contour
order hash, trial-state hash, cut, side (separate ket/right and bra/left stores),
precision (the concrete cache type), and slice count. It never consumes the
double live stack.

For each same-proposal reference the closest unaffected right-prefix and
left-suffix checkpoints are restored, validated for dimensions, key identity,
and MP orthogonality residual, then propagated to the proposal cut in canonical
operator order. Checkpoints remain on the existing eight-factor QR boundaries.
MP160/320/640 caches do not share numerical subspaces, so consecutive-precision
agreement remains independent.

An accepted proposal rekeys checkpoints whose represented contour segment is
provably unaffected and invalidates every affected prefix/suffix. A rejected
proposal performs no invalidation. A configuration/order/trial mismatch clears
the stale cache before use. If checkpoint validation or the new precision path
is untrusted, that precision is cleared and rerun through the Phase 3E oracle;
an untrusted legacy result still fails closed. Cache state is performance-only
and is intentionally not serialized in a restart checkpoint.

The fail-first Phase 3F unit initially did not compile because the cache API and
profile counters were absent. After the minimal implementation it verifies real
restores, exact accepted invalidation, rejected preservation, and deliberate
stale-cache rejection. During the 512-case audit, a structured FNV collision in
the first cache-key draft was exposed by the stale test; the order-sensitive
canonical mixer was corrected before any final benchmark or evidence was kept.

## Oracle and safety validation

`oracle_equivalence.csv` contains 512 comparisons of the new path against both
the Phase 3E cached oracle and the complete no-cache legacy oracle: 259 ket-side,
253 bra-side, 213 block-boundary, 495 accepted, and 17 rejected proposals. All
512 pass. The maximum relative ratio difference is `6.15e-171`, Z2/status/
precision/decision mismatches are zero, and the endpoint residual delta is zero.
There were 1,798 actual checkpoint restores and one deliberately stale cache was
rejected and safely rebuilt. `cache_invalidation_tests.csv` records each cache
lifecycle event.

Frozen regressions remain exact:

- Task 4: ratio `0.70741914065425615`, uniform `0.52928038756354223`, accepted,
  final state hash `949208207496548183` after 64 further proposals.
- Task 8: ratio `0.0026573889530685144`, uniform `0.56577673777691728`, rejected,
  final state hash `14005514076608941581` after 64 further proposals.
- Task 4/8 fallback failures: zero.
- Checkpoint/restart retained measurement files are byte-identical. Their
  SHA256 values are respectively `55f1f2d421bf91e109a97c2d44cbfdf0bc5f5236c3458926cf5612dcd5bb8e12`
  and `c934271b96f77485d5ae4925c47bac9d0cfc7c3a1ac993ecd48494a59d88cc8d`.
- Ratio reciprocity: 1,002 rows pass, maximum error `6.42e-15`, detailed-balance
  failures zero.
- The Phase 3A OBC+hs1 test reports physical Z2 `+1` at every measurement; the
  Phase 3F L=6 control also has average Z2 `+1`, zero fallback failures, and
  endpoint residual `1.49e-14`.

Phase 1 (6/6), Phase 2 core and validation (5/5 and 8/8), Phase 3A core and
validation (4/4 and 8/8), Phase 3B I/O/production contracts, Phase 3C safety,
Phase 3E operator-cache/GenericComplex, production-semantics, frozen replay,
restart, reciprocity, Phase 3F cache unit, and the 512 oracle comparisons all
pass. Dense-Fock maximum weight/Green errors remain `7.50e-15`/`1.64e-14`;
trajectory and Z2 mismatches are zero.

## Performance

The controlled comparison ran V=4 and V=6 simultaneously, one CPU per process,
with the same Phase 3E seeds and true-sweep workload: L=6, theta=6, PBC+hs0,
one burn sweep plus two measurement sweeps. Detailed counters are in
`performance_before_after.csv`.

| Case | Phase 3E seconds/sweep | Phase 3F seconds/sweep | total speedup | fallback speedup | thin-QR speedup | propagation speedup |
|---|---:|---:|---:|---:|---:|---:|
| V=4 | 100.85 | 15.78 | 6.39x | 8.83x | 47.68x | 73.70x |
| V=6 | 259.29 | 35.63 | 7.28x | 8.74x | 46.51x | 69.08x |

Fallback counts are unchanged at 285 and 756, with zero failures and zero legacy
retries in the normal benchmark. Restore hit rates are 98.60% and 98.08%.
Maximum endpoint residuals remain bit-identical to Phase 3E at `2.58e-14` and
`8.44e-9`. Observables, acceptance, final HS hash, and final RNG hash are exact.

Peak RSS changed from 31,908 to 36,564 KiB for V=4 (+14.6%) and from 31,272 to
36,984 KiB for V=6 (+18.3%), about 5--6 MiB per process. The checkpoint stores
peak at about 2.3 MB by the conservative profile estimator. Given the 6.4--7.3x
true-sweep improvement and bounded memory increase, retaining this optimization
is justified.

## Build and files

The validated GCC build used `-no-pie -O2 -fopenmp -std=c++17
-DPFQMC_SCALE_SAFE_UDT`, the project Eigen/Boost includes, the existing static
PFAPACK archives, and sequential MKL. `build_phase3b.sh` now includes both Phase
3F test executables.

Changed production files:

- `inc/pure_projector_mp.h`
- `inc/pure_projector_fast.h`
- `reproduction/pure_projector_kitaev/pure_projector_driver.cpp`
- `reproduction/pure_projector_kitaev/build_phase3b.sh`

New tests/evidence:

- `reproduction/pure_projector_kitaev/phase3f_mp_subspace_test.cpp`
- `reproduction/pure_projector_kitaev/phase3f_oracle_equivalence.cpp`
- `oracle_equivalence.csv`
- `cache_invalidation_tests.csv`
- `performance_before_after.csv`
- `trajectory_hashes.csv`

No HPC production scan was submitted. Phase 4 was not started.
