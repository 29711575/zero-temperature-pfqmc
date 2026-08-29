# Driven `Dexp` comparator fix report

## Scope and provenance

- Branch: `fix-driven-dexp-comparator`
- Clean committed baseline: `4d79bb7ac2b734bc10667d917b914ce316a39a39`
  (`Freeze MP Z2 trust-gate handoff`)
- Issue scope: driven validation/comparator reconstruction ignored the exponent
  in scale-safe `UDT::D`/`UDT::Dexp` storage.
- HPC result directory:
  `/home/sunxr/new-pfqmc-fix-driven-dexp-comparator/reproduction/driven_kitaev/driven_dexp_comparator_fix_20260829/`
- Formal PBS job: `153851.mgt`, 1 CPU, no GPU, 24 h requested walltime.

No `main` or `develop` branch was checked out or modified. No proposal-safety,
I/O, PBC structural, or MP trust-gate worktree was modified. The isolated HPC
clone and branch were created directly from the committed baseline above.

## Fix

Both copies of `driven_fastupdate_check.cpp` now use the shared validation-only
`driven_comparator_reconstruction.h` helpers. The helpers reconstruct the
scaled solve with existing `UDT` accessors:

- `dLargeInverse()` and `dSmallPart()` build the scaled core used by the
  independent comparator Green reconstruction and `ultra_proxy_at()`.
- `actualD()` supplies materialized full scales for diagnostic dense
  reconstruction, scale spread, and the large-scale term in `logdet_full()`.
- `dGreaterThanOne()` selects the large-scale logdet contribution.

There is no direct `f.D` use left in either driven comparator. The production
fast-update kernel, production driver, `src/`, and `inc/` are unchanged.

## Regression results

| Regression | Result | Key result |
|---|---:|---|
| Synthetic `D=(0.5,0.5)`, `Dexp=(21,-19)` | PASS | Legacy reconstruction relative error `0.99999952316295548`; fixed error `0` |
| Comparator reconstruction vs `UDT::onePlusInv()` | PASS | Relative error `0` |
| Synthetic comparator vs direct `2*(I+B)^-1` | PASS | Relative error `0` |
| Synthetic comparator logdet vs direct dense logdet | PASS | Absolute error `0` |
| Pathological `L=6,Vf=6,dt=0.1,flips=300,stb=10` | PASS | Ratio/proxy false errors removed; details below |
| Driven static comparison | PASS | All reported operator/sign/Green/observable differences `0` |
| Driven smoke | PASS | JSON `status=complete` |
| Production source/trajectory check | PASS | No production source diff; Green metrics identical before/after |

Pathological before/after comparison:

| Metric | Before | After |
|---|---:|---:|
| max ratio complex absolute error | `735.16294625667194` | `5.7633405722517015e-6` |
| max ratio magnitude absolute error | `735.16294625667194` | `6.0568731896637473e-7` |
| max comparator solve residual | `210.71238278723098` | `3.5442400737561756e-14` |
| max Green absolute error | `2.2790998949061056e-6` | `2.2790998949061056e-6` |
| max error before stabilization reset | `3.50928396342152e-6` | `3.50928396342152e-6` |

The repaired ratio error is within 2.53 times the Green error, and the repaired
proxy residual is far below that numerical scale. The identical Green and
pre-reset values demonstrate that this validation-only change did not alter
the sampled trajectory or production fast-update behavior. The corrected
reported `D` spread is `4.8590438324853841e93`; the old mantissa-only value
`1.9978317971527826` was not the physical scale spread.

Machine-readable details are in `regression_summary.csv` and
`before_after_comparator.csv`. Raw before/after traces and JSON summaries are
retained in this directory.

## Executable hashes

- Fixed comparator: `23d9d1d38836b95b6bd76b2d6ac7675f3a41edd0e7362a3f2101b0d1d58c03b1`
- Baseline comparator: `f64cc91c67148aded7c4e9e1b7d3e09e284ed1f92bf34426a2b749f0b96cc2cf`
- Production driven smoke executable: `71e013a1ba050834d6d3f8727c29902abf126f19afc4234eb3246ceacfe4669d`
- Static comparator: `dee61ae3080ec6f89b5405ab4e31ffa413180207d6d2390d0aba3a7a2de3200e`
- Synthetic regression: `d798b87b48e245031f847536135d3e4ebb0679129599ee23d7fb790fe0a03263`

The branch remains unmerged for unified integration.
