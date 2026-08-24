# Driven 153416 numerical-failure audit

No Vc fit was performed. The driven contour, tolerance, sampling algorithm, and core PfQMC sources were not changed.

## Static audit of the 208 exit-code-3 tasks

All 208 failures have `failure_reason=numerical_diagnostic_exceeds_tolerance`. Trigger counts are non-exclusive: SIGN_IMAG=198, OBSERVABLE_IMAG=0, GREEN_SKEW=28. Sign corrections are not part of the exit-code-3 predicate and were zero in these failures; no OTHER trigger was found.

| R | planned | failures | failure rate | sign-imag | observable-imag | Green/skew | max sign/tol | max obs/tol | max Green/tol |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 264 | 181 | 68.6% | 177 | 0 | 19 | 69.39 | 0.667 | 14.87 |
| 1.0 | 120 | 27 | 22.5% | 21 | 0 | 9 | 7.47 | 0.017 | 6.72 |
| 2.0 | 72 | 0 | 0.0% | 0 | 0 | 0 | 0.00 | 0 | 0.00 |

The failures are therefore strongly rate-concentrated: R=0.5 accounts for 181/208, R=1 for 27/208, and R=2 for 0/208. Detailed `(L,V,R)` rates and magnitudes are in `driven_failure_audit_by_LVR.csv`.

## Deterministic same-configuration audit (PBS 153420)

The diagnostic copy adds a read-only callback exactly at the production center capture, before later operators in `rightSweep` can change the HS configuration. Direct Pfaffian sign and full-contour Green are evaluated only at the first threshold event, consume no RNG, and use the same boundary/configuration as the tracked sign and fast Green.

| task | L | V | R | seed | first meas. | trigger | tracked/direct ± | mismatch | |Im tracked| | fast/full Green rel. | max fast obs imag | class |
|---:|---:|---:|---:|---:|---:|---|---|---:|---:|---:|---:|---|
| 38 | 26 | 3.80 | 0.5 | 4226003 | 3863 | SIGN_IMAG | -1/-1 | 0 | 1.024e-07 | 6.992e-14 | 3.811e-17 | harmless complex sign drift |
| 443 | 34 | 3.80 | 0.5 | 4234012 | 2456 | SIGN_IMAG | 1/1 | 0 | 4.582e-06 | 6.537e-13 | 6.909e-15 | harmless complex sign drift |
| 46 | 26 | 3.80 | 0.5 | 4226011 | 4620 | SIGN_IMAG | 1/1 | 0 | 1.029e-06 | 5.943e-14 | 2.379e-17 | harmless complex sign drift |
| 187 | 26 | 4.00 | 0.5 | 4226008 | 1547 | SIGN_IMAG | 1/1 | 0 | 5.245e-07 | 1.868e-14 | 2.330e-17 | harmless complex sign drift |
| 324 | 26 | 4.20 | 0.5 | 4226001 | 17 | SIGN_IMAG | 1/1 | 0 | 1.019e-07 | 2.645e-13 | 3.921e-17 | harmless complex sign drift |
| 577 | 34 | 4.00 | 0.5 | 4234002 | 4047 | SIGN_IMAG | -1/-1 | 0 | 1.961e-06 | 2.333e-13 | 1.874e-16 | harmless complex sign drift |

## Conclusion

- All 6/6 representative reruns reproduced a sign-imaginary threshold event deterministically.
- Genuine discrete ±1 sign mismatches: **0/6**.
- At the same HS configuration and boundary, fast/full Green relative errors were at most 6.537e-13 and observable imaginary parts at most 6.909e-15; no numerical instability was seen at those first events.
- The sampled first events support **harmless numerical complex-sign drift with unchanged ±1 eigen-sign**, not a sign-transport mismatch.
- Caveat: 28 original tasks had a Green/skew maximum above tolerance. In task 324, which belongs to both categories, the first event was sign-imag drift and was fully consistent; because the rerun stopped at the first event, the later original Green/skew maximum was not independently frozen. The static audit flags those 28 for possible narrower Green-event follow-up rather than declaring all of them harmless.
