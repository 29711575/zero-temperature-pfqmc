# PfQMC integration staging status

Staging repository: `/home/sunxr/PfQMC-integration-staging`

Base: `main@2ae3bb7d81c4c4153b773b81dfe897405a94adab`

## Integrated changes

1. `8887c5f3d7ff612f7920ebd69608a25e5dd33ec8` — production-only optimized
   UDT rank-loss guard, sourced from `5b21b110...` via validated clean port
   `b33dc7d...`. It retains the 45-bit lost-bits fail-closed rule, staged
   detailed checking, final `||U^dagger U-I||/sqrt(n) <= 1e-6` gate, and
   finite/bookkeeping diagnostics. Profiling, MP fallback, Householder, and
   stress instrumentation are absent.
2. `f29af53cef565d04ec3598adeeefc443be203245` — output-layer fix from
   `62a1d6a...`; `average_sign=0` remains a valid completed JSON result and
   sign-reweighted observables are emitted as `null` with explicit unresolved
   status.
3. `834ff596cec611bb78a1d546021812c2271c5370` — minimal optional left-sweep
   Green recovery from `5dcf46b...`. It is disabled by default, consumes no
   additional RNG, does not reconsider acceptance or change transported sign,
   preserves the accepted HS/B mutation, and rebuilds at most once per alarmed
   operation.
4. Durable regression source/scripts are consolidated under
   `reproduction/projector_kitaev/core_regression_tests/`. Binaries and outputs
   are directed to external directories and ignored.

No condition-aware/adaptive MP ratio production path and no discrete Z2 sign
transport change is included.

## Independent build and validation

Build directory: `/home/sunxr/PfQMC-integration-staging-build`

Output directory: `/home/sunxr/PfQMC-integration-staging-results`

| Test | Result | Key evidence |
|---|---|---|
| tiny complete HS enumeration | PASS | 256 configurations; max center Green direct error `7.11e-15`; average phase real 1 |
| V=0 Gaussian | PASS | finite result; raw sign imaginary `3.01e-17` |
| L10 control | PASS | guard trigger 0; max fast/full Green `3.18e-15`; sign mismatch 0 |
| task88 short | PASS | guard trigger 0; max lost bits `8.1832`; margin `36.8168` bits; Green `3.05e-14`; sign mismatch 0 |
| task92 short | PASS | guard trigger 0; max lost bits `9.2313`; margin `35.7687` bits; Green `9.18e-13`; sign mismatch 0 |
| guard vs no-guard trajectory | PASS | L10/task88/task92 raw CSV files are byte-identical |
| 3000 local flips | PASS gate | anomalies 0; direct/sequential and direct/stabilized-full ratios pass; max accepted Green error `1.32e-7`; dense determinant-squared diagnostic remains conditioning-only |
| representative all-boundary right Green | PASS | 180 boundaries, finite 180/180, max relative error `2.18e-10` after fixing test snapshot timing |
| UDT normal stress | PASS | max reconstruction `3.73e-16`; max normalized solve residual `7.50e-17` |
| UDT rank-loss synthetic guard | PASS | n12 +/-20 passes; n24 +/-40, n12 +/-500, n52 +/-1500 and +/-2000 all fail closed; no finite nonunitary U escaped |
| reality symmetry / MP ratio | PASS | symmetry and ratio imaginary residual 0 at 80/160-digit small-L samples; `r^2=Q` residual about `1e-157` |
| zero average sign output | PASS | process completed with valid JSON, `average_sign=0`, explicit `unresolved_zero_average_sign`, all reweighted observables `null` |
| optional left recovery, L12 V5 | PASS for recovery correctness; performance caution | finite, end-of-sweep Green equals stabilized full rebuild; hard trajectory triggered frequently, so feature remains default-off/configurable |

The matched CSV result simultaneously verifies that the UDT guard and the
default-off left recovery leave the normal path unchanged. No normal-QMC guard
trigger was observed.

## Test-suite notes

- The local-update CSV retains `r^2/Q` from an unstabilized dense determinant as
  a diagnostic. It is not the gate under poor conditioning; stabilized
  direct/sequential/full ratio agreement and rebuilt Green are the gate.
- The consolidated boundary driver stops at the requested boundary before the
  independent rebuild. An earlier draft let the sweep mutate later HS fields
  before comparison and was discarded as an invalid test definition.
- Previous core-regression coverage (L=6,10,12; OBC/PBC; V=2,4,6; both sweep
  directions) remains recorded in the source branch. The staging smoke is a
  representative independent rebuild, not a rerun of all historical rows.
- The optional recovery threshold is a diagnostic/runtime configuration, not
  a new physical definition. Frequent triggers in a deliberately hard initial
  L12,V5 trajectory are why it is not enabled by default.

## Pending integration decisions

- PBS 153474 v1 remains running. PBS 153573 optimized v3 completed during the
  final inventory check, but its result was deliberately not analyzed here.
  Both require dedicated analysis before any condition-aware ratio production
  logic is considered.
- Discrete Z2 sign transport remains explicitly out of scope.
- Householder and MP UDT fallback remain out of staging; the optimized
  fail-closed guard is sufficient as the current safety candidate because real
  QMC retained roughly 35.9 bits of margin and triggered zero times.

## Classification

Normal-path staging and the UDT guard are **PASS / SAFE TO PORT**. The optional
left recovery is integrated safely as **default off**; its recovery action is
correct, while threshold/rate tuning remains follow-up before any default-on
proposal. Staging must not be merged to main until the user reviews this status
and the protected adaptive-ratio results are separately analyzed.
