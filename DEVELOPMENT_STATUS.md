# PfQMC development status

## Stable main

The stable code base descends linearly from original main
`2ae3bb7d81c4c4153b773b81dfe897405a94adab` through integration hygiene
`55fee175752a896c12e2acf79d2a252c7db74e2f`.

Validated and enabled by default:

- scale-safe UDT with mantissa plus base-2 exponents;
- exponent-aware QR and power-of-two equilibrated factorization/solve;
- staged UDT effective-rank-loss guard (45 lost bits, fail closed);
- final adjoint-orthogonality gate `||U†U-I||/sqrt(n) <= 1e-6`;
- C++17 and `PFQMC_SCALE_SAFE_UDT` across Make, CMake, and projector builds;
- strict JSON with nullable non-finite/unresolved values and feature provenance;
- successful JSON output when `average_sign == 0`, with reweighted observables
  marked unresolved/null;
- consolidated core regression tests.

Present but disabled by default:

- left-sweep Green consistency/rebuild recovery. This remains opt-in so the
  normal trajectory is identical to the validated legacy sweep path.

Not on stable main:

- condition-aware adaptive multiprecision ratio fallback;
- a new discrete Z2 sign-transport implementation;
- MP UDT fallback or exponent-aware Householder prototypes;
- replay/audit event instrumentation.

## Incomplete work and limits

Optimized condition-aware ratio v3 is archived on
`feature/condition-aware-ratio-v3`. PBS 153573 was numerically clean and much
faster than v2, but PBS 153577 failed to capture the required
k132/operator909/aux4 target row (`target_events=0`). V3 is therefore not
approved for targeted validation or stable integration.

The running PBS 153474 v1 executable has an exact SHA-256 but no trustworthy
source commit tying its binary to the currently visible dirty source tree.
`archive/condition-aware-ratio-v1` is explicitly a best-effort source snapshot,
not a reproducible binary claim.

Average sign near zero is not itself a numerical failure. The current stable
main is a numerical implementation baseline, not a final large-size physics
result.

## Recommended workflow

- `main`: stable, validated integration baseline.
- `develop`: start point for unified new work.
- `feature/condition-aware-ratio-v3`: continue ratio validation only after
  fixing the missing target-event observer evidence.
- `archive/*`: read-only historical diagnosis/prototype references.

New development and tests should be performed only in
`/home/sunxr/new-pfqmc-main` or worktrees created from it.
