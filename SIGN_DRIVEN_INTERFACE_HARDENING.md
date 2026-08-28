# Sign and driven interface hardening

## Scope and provenance

- Repository: `/home/sunxr/new-pfqmc-main`
- Branch: `harden-sign-and-driven-interfaces`
- Base: `develop@1956c8c12b14e28f6904314874da62046429a501`
- Condition-aware ratio v3 was not integrated. The UDT numerical algorithm,
  guard thresholds, Green-recovery policy, and Markov/RNG logic were not changed.

## Implemented changes

### Pfaffian status API

`PfaffianResult` now carries a typed status, phase/value, LAPACK `info`, and the
minimum finite pivot magnitude. `pfafWithStatus()`, `signOfPfafWithStatus()`,
`pfaffianForSignOfProductWithStatus()`, and `PfQMC::getSignRawWithStatus()` fail
closed on invalid dimensions, LAPACK errors, zero pivots, and non-finite pivots.
Compatibility wrappers remain, but throw on failure instead of returning `+1`.
Production sign diagnostics use the status-returning API.

The legacy `DenseOperator::signPf_g0_inv` cache is not used by current sign or
propagation code and can be singular for valid OBC dense operators. It now
stores the status and an unavailable (`NaN`) cached value on failure; it neither
fabricates `+1` nor aborts construction of an otherwise valid contour.

### Read-only periodic sign diagnostics

Static projector, bins/static-guard, driven, and maintained regression drivers
no longer overwrite transported `q.sign`. JSON reports comparison counts and
one of `not_sampled`, `agreement`, `mismatch`, `raw_check_untrusted`, or
`raw_check_unavailable`, including LAPACK info and minimum pivot when available.
The legacy `sign_corrections` counter consequently remains zero.

### Driven output and record-only behavior

- `PFQMC_MULTIPRECISION_RECORD_PROXY` installs only the read-only condition
  recorder; the fallback enable flag remains false and Green is never replaced.
- The sign error is computed from bin-level sign means rather than hard-coded.
- Exactly zero average sign completes normally; sign-reweighted observables are
  JSON `null` with status `unresolved_zero_average_sign`.
- All unresolved/non-finite sentinels use standard JSON `null`.
- Projector provenance/feature fields and raw-sign status schema are shared.
- The guard-off diagnostic is named `min_observed_update_denominator`; the
  guard-specific field is `guard_min_prepared_denominator` and is null when the
  guard is disabled.

### Condition proxy and QR audit

`fullContourCoreConditionAtBoundary()` replaces `T.inverse()` with checked
`FullPivLU` solve/rcond and a normalized solve-residual check. This path is
read-only and does not affect a normal trajectory.

The QR diagonal audit records `abs(Im Rii)/max(abs(Re Rii), tiny)` without
changing the existing `abs(real(Rii))` calculation. It covers QR calls from a
short real QMC trajectory and independent random complex matrices.

## Validation

Validation used independent executables under
`/home/sunxr/pfqmc_hardening_validation_20260826`; no production executable or
PBS input was modified.

### Build and unit checks

- Clean validation build: PBS `153592`, complete.
- Pfaffian regular case: `success`; singular case: `lapack_failure`,
  `info=1`; compatibility wrapper threw as required.
- Old inverse and new checked-solve condition proxies agreed to relative
  `2.1271114574762845e-16` at the normal control.
- QR audit: 2160 QMC and 1200 random-complex diagonal samples. Both had maximum
  `abs(Im Rii)/abs(Re Rii) = 0`; no sample exceeded `1e-12`. Thus the existing
  real-diagonal assumption is observed for the LAPACK QR used here, but its
  production calculation was intentionally left unchanged.

### Unified regression

PBS `153593` finished with exit status 0 (walltime 4m18s). Strict parsing loaded
all normal and forced-zero JSON files. Tiny enumeration, Gaussian exact, L10,
task88/task92 short smoke, 3000 local flips, all-boundary Green, left-recovery
harness, UDT normal/rank-loss stress, reality symmetry, and zero-sign tests all
completed. Representative QC:

- tiny maximum direct-center Green error: `7.21e-15`;
- all-boundary maximum Green error: `1.07e-10`, non-finite count 0;
- task88/task92 raw-sign mismatch 0, sign corrections 0, UDT guard triggers 0;
- L10 raw-sign status `agreement`, Green diagnostic maximum `3.22e-15`;
- forced-zero projector, bins, and static-guard outputs completed with
  `average_sign=0`, null reweighted observables, and explicit unresolved status.

The hardened and baseline L10 measurement CSV files were byte-identical:
`SHA256 8427f42f48b30b5a9af5544ca7cf47dac5d212a06eaa4e88051ccfd3b7274460`.

### Static and driven ED smoke

All cases completed with raw-sign status `agreement`, zero mismatch/correction,
finite observables, and no fallback.

| case | observable | QMC | exact | difference / QMC SEM |
|---|---:|---:|---:|---:|
| static L4 V2 | S_pi | 0.0898211 | 0.0910233 | 1.05 |
| static L4 V2 | S_pi+dq | 0.0574876 | 0.0569285 | 1.38 |
| static L4 V2 | R_cdw | 0.359976 | 0.374572 | 1.22 |
| driven static limit L4 V2 | S_pi | 0.117697 | 0.119438 | 0.74 |
| driven static limit L4 V2 | S_pi+dq | 0.0557012 | 0.0538405 | 1.94 |
| driven static limit L4 V2 | R_cdw | 0.526742 | 0.549217 | 1.30 |
| driven ramp L4, V0=0 to V=2 | S_pi | 0.0914858 | 0.0897087 | 0.80 |
| driven ramp L4, V0=0 to V=2 | S_pi+dq | 0.0584724 | 0.0592672 | 1.53 |
| driven ramp L4, V0=0 to V=2 | R_cdw | 0.360859 | 0.339337 | 1.03 |

The forced-zero driven case completed with `average_sign=0` and all reweighted
observables null. Its normal driven smoke reported a nonzero bin-level sign SEM
(`0.0114197`), confirming that `average_sign_err` is no longer hard-coded.

### Record-only trajectory

The fully disabled and `PFQMC_MULTIPRECISION_RECORD_PROXY=1` runs produced
byte-identical bin records:
`SHA256 a7419058ae8f97f7b76874b57b267ac5d3114657991625cc3014ad7ec233ea42`.
The record-only run collected 6484 condition samples while reporting fallback
disabled and fallback count 0. Acceptance, signs, observables, RNG/HS trajectory,
and all QC values were identical.

## Conclusion

PASS. The interface hardening is isolated from the numerical algorithms and
preserves the tested normal trajectory. Pfaffian failure can no longer silently
become `+1`; periodic sign checks are read-only; driven record-only mode cannot
replace Green; zero-sign JSON is standards-compliant; and the checked condition
proxy matches the previous diagnostic on a regular point. Condition-aware ratio
v3 and discrete Z2 sign transport remain outside this branch.
