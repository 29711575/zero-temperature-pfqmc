# long first-divergence audit

Parameters: L=6, V=4, PBC, theta=10, beta_trial=8, dt=0.1, hs_scheme=0, guard OFF, seed=984035; 500 burn plus 5000 measurements. Diagnostics are read-only and never correct `q.sign`.

## Coarse localization

No true event was found. All 11,001 boundary-0 sweep-end checks and all 5,000 center-boundary checks made on the same live configuration have matching tracked/direct Pfaffian ± signs and stabilized full-contour Green matrices.

| location | checks | ± mismatches | Green drifts | max Green rel. error |
|---|---:|---:|---:|---:|
| sweep end, boundary 0 | 11001 | 0 | 0 | 7.942e-09 |
| center, same configuration | 5000 | 0 | 0 | 1.023e-11 |

## Precision stage

There is no “last normal sweep → first abnormal sweep” for this seed/run, so no proposal window was activated. Running proposal diagnostics without an abnormal bracket would not be a valid first-divergence localization.

## Conclusion

No PfQMC operation is first to fail in this deterministic run: accepted local updates, stabilization checkpoints, right/left transitions, sweep ends, and center capture remain internally consistent. The earlier 509/2000 audit mismatch compared a sign captured at the center with `getSignRaw()` evaluated after the remainder of the right sweep had changed the HS configuration. It was therefore a cross-configuration diagnostic artifact, not evidence of a transported-sign mismatch.

No core algorithm was modified.
