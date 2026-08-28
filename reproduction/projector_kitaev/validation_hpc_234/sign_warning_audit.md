# sign warning audit — corrected

## Static driver sign usage

At measurement, `rightSweep(center,&g,&sg)` captures the complex transported sign at the center. The driver immediately converts it to the eigen-sign `s = (sg.real() >= 0 ? +1 : -1)`. It accumulates `ss += s`, `sp += s*S_pi_cfg`, and `sd += s*S_pi_dq_cfg`; final values are `S_pi=sp/ss`, `S_pi_dq=sd/ss`, and `R_cdw=1-S_pi_dq/S_pi`. Neither the complex phase, `abs(sg)`, nor the unthresholded value `real(sg)` enters reweighting.

`SIGN_IMAG_WARN` is a collector label when `max(abs(Im(sg))) > 1e-6`. `SIGN_CORRECTION` is a collector label when the driver reports `sign_corrections > 0`.

Every 20 measurements, before the next right sweep, the driver evaluates `z=getSignRaw()` on the current live contour. If `abs(q.sign-z)>1e-2`, it **does modify** `q.sign`, setting it to `+1` or `-1` according to `z.real()`, and increments `sign_corrections`. Otherwise the call is read-only.

## Existing-data comparison

For five severe warning seeds, recomputing from configuration observables using an explicit `sign>=0 ? +1 : -1` gives the same estimator as the stored numerator columns. Maximum absolute difference is 0.000e+00 and maximum difference is 0.000e+00 sigma. Detailed values are in `collected/sign_usage_reweight_compare.csv`.

Therefore the recorded imaginary component does not change S_pi, S_pi_dq, or R_cdw in these data. The corrected same-configuration audit contains zero genuine transported/direct $\pm1$ eigen-sign flips. These warnings are **numerical complex drift with unchanged $\pm1$ eigen-sign**, not sign mismatches in the reweighting estimator.

The targeted records retain only the thresholded direct sign and its imaginary component. They do not retain the real component of the raw complex sign, so `min |Re(s_raw)|` and the maximum normalized raw-sign distance from the nearest $\pm1$ are unavailable and are not inferred.

## Retraction of the earlier mismatch claim

The earlier “509/2000 real mismatches” conclusion is withdrawn. That audit compared `sg` captured at the center with `getSignRaw()` called after the remainder of the right sweep had updated additional HS fields. Those signs belonged to different HS configurations. Correct same-configuration long-run checks found 0 mismatches in 11,001 sweep-end checks and 0 mismatches in 5,000 center checks.

The optional `diagnostic_stride` block of `projector_bins_driver.cpp` has now been corrected to compare the live transported sign with `getSignRaw()` at the same completed-right-sweep configuration and boundary. The validation_hpc_234 production tasks used diagnostic_stride=0, so this diagnostic-only correction does not alter their observables or warning labels.

No core algorithm was modified.
