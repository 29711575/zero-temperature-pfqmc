# Scale-safe UDT validation

Final status: **PASS**. The original long-contour UDT dynamic-range failure is
fixed by the scale-safe UDT implementation. The original failing task88
parameters completed the full burn plus 3000-measurement run without NaN/Inf,
LAPACK failure, sign mismatch, or Green-function divergence.

## Scope and provenance

This validation concerns the static-projector PfQMC implementation for the
interacting Kitaev chain. The numerical change remains isolated from the main
branch: no checkout, merge, rebase, or cherry-pick was performed as part of
this validation.

The HPC execution tree is `/home/sunxr/PfQMC-scale-safe-udt`. The copied HPC
source tree does not contain Git metadata, so a branch name and commit cannot
be reconstructed from that copy. The separately visible main worktree was
checked read-only and is `main` at
`6fb1b303b425771bf1aa3324c555c7e48c22f6b6`; it was not modified. The D PBS
script called the isolated executable below directly, not an executable from
that main worktree.

| artifact | SHA-256 |
|---|---|
| `bin/regression_driver_old` | `483cff314d4abe661de6073081d340b8414483028a73511eea5512e9ecad4b35` |
| `bin/regression_driver_scale_safe` used by A/B/C/D and job 153446 | `95a973a1422949f36819a3614785bf704af41e62e8d1e1805507ce240f986ec3` |
| `src/pfqmc.cpp` validation snapshot | `5fb99d84438ebab7103e281283a9e07e1112ec9476634110e28e82d93e705d0f` |
| `inc/qr_udt.h` validation snapshot | `5525c0c317ad70200713e8bbbf0a272111b4643ac30cdbd53a34b4f4aecc8b6b` |
| `regression_driver.cpp` validation snapshot | `efbecba4444cb9d97c8c483effb51ad2917073c04f00ac53c1fe63702b4e376c` |
| `stage_D_task88_full.pbs` | `2ca8311384b28b8cb0a9e5cc483028e8e80c47a5cc6e4e9ad7b0011d5993c848` |
| `validation_gate.py` | `782a4c82a18c830fdc99911a2b1f67b257cdffabef4effc685493118b1766430` |

The scale-safe path stores each UDT scale as a mantissa plus base-2 exponent,
uses exponent-aware QR for non-materializable dynamic ranges, applies exact
power-of-two equilibration, and factorizes/solves the scaled two-UDT core.
The old path and executable were retained for comparison.

## A: short-contour old/new validation

**PASS.** Algebraic unit tests, tiny enumeration, Gaussian checks, and the
matched `L=26, V=4, theta=26, seed=750049` control agreed between old and new
paths. In the matched control, sign and acceptance agreed row-for-row;
maximum row differences in `S_pi`, `S_pi_dq`, and `R_cdw` were `1.80e-16`,
`3.64e-17`, and `2.13e-14`. Same-configuration sign mismatches were zero.
The v2 control exercised exponent range `[-781,781]` with minimum scaled-core
rcond `4.41e-7` and was row-for-row identical to v1.

## B: frozen 153407 baseline validation

**PASS.** All seven matched-seed cases passed the fail-closed gate. Observable
differences were at roundoff level, average sign and acceptance were exactly
matched, and every same-configuration sign-mismatch count was zero.

| case | exponent range | minimum scaled-core rcond | max fast/full Green relative error |
|---|---:|---:|---:|
| guard V4 off | `[-315,316]` | `6.459877e-8` | `2.994081e-14` |
| guard V4 on | `[-318,319]` | `6.459877e-8` | `2.994081e-14` |
| guard V6 off | `[-393,394]` | `7.125994e-8` | `3.817185e-13` |
| guard V6 on | `[-399,400]` | `1.532817e-7` | `2.354489e-13` |
| interval 5 | `[-537,538]` | `7.872173e-8` | `1.881224e-13` |
| interval 10 | `[-537,538]` | `7.872173e-8` | `2.763382e-13` |
| interval 20 | `[-537,538]` | `7.872173e-8` | `1.084332e-12` |

Small average sign in baseline cases is a physical sign problem and was not
treated as a numerical failure.

## C: original long-contour failures, 500-measurement smoke

**PASS.** Tasks 88, 89, 90, and 92 completed 500 burn plus 500 measurements.
All exercised exponent-aware QR and crossed the original failure regime with
no computational NaN/Inf, LAPACK abort, sign correction, same-configuration
sign mismatch, or Green-function divergence.

| task | `(L, theta, seed)` | exponent range | minimum scaled-core rcond | max fast/full Green relative error | max observable imaginary part |
|---:|---|---:|---:|---:|---:|
| 88 | `(26,39,750049)` | `[-1147,1149]` | `2.118468e-6` | `3.734158e-14` | `1.191542e-14` |
| 89 | `(26,39,750050)` | `[-1159,1160]` | `9.967195e-7` | `1.20e-13` | `5.65e-14` |
| 90 | `(26,39,750051)` | `[-1166,1167]` | `4.479888e-6` | `7.89e-14` | `4.32e-13` |
| 92 | `(34,51,750069)` | `[-1502,1503]` | `1.258375e-6` | `<1e-8` (gate) | `<1e-8` (gate) |

## D: full original task88 run

**PASS.** PBS job `153446.mgt` finished with `Exit_status=0` in `02:27:21`.
Its status file records completion and identifies the same job. Parameters were
`L=26, V=4, theta=39, beta_trial=8, dt=0.1, seed=750049`, with 500 planned burn
steps followed by all 3000 requested measurements. `measurements.csv` contains
exactly 3000 data rows (`measurement=0..2999`).

| QC item | result |
|---|---:|
| scale-safe dense QR calls | `1,342,338` |
| exponent-aware QR calls | `1,573,321` |
| checked solves | `3,405,102` |
| exponent range | `[-1159,1161]` |
| minimum scaled-core rcond | `1.3550875935527368e-6` |
| same-configuration sign mismatches | `0 / 30` diagnostics |
| sign corrections | `0` |
| maximum sign imaginary part | `3.057089225912622e-7` |
| maximum observable imaginary part | `6.556196420877105e-14` |
| maximum fast/full Green relative error | `5.140891191852719e-13` |
| maximum diagnostic `S_pi` absolute difference | `8.708311849403572e-16` |
| maximum diagnostic `R_cdw` absolute difference | `7.294165271787278e-14` |
| acceptance | `0.60708752991453` |
| average sign | `1.0 +/- 0.0` |
| sign recomputations | `150` |
| negative signs | `0` |

All 3000 physical trajectory rows are finite. All 30 scheduled diagnostic rows
are finite; the `nan` values in diagnostic columns on the other 2970 rows are
the defined “diagnostic not sampled on this row” sentinel, not propagated
numerical values. Likewise, `min_update_denominator=inf` is the defined
guard-off/no-proposal-counter sentinel. Neither appears in a computed physical
observable, Green diagnostic, factorization, or solve. Searches of stderr and
the PBS output found no LAPACK failure, LU abort, or other numerical error.

C task88 and D use the same physical parameters, seed, burn, executable, and
RNG definition; only diagnostic cadence differs (50 versus 100). The first 500
physical trajectory rows (`measurement`, sign, signed and unsigned observables,
`R_cdw`, acceptance) agree exactly. Their normalized trajectory SHA-256 is
`ba62a18102f389895e1eb016a5654ee440eabad88eb977b21c473850dea5bcc9`.

## Overall conclusion

**A/B/C/D PASS.** The original long-contour dynamic-range failure has been
reproduced as the validation target and fixed by the scale-safe UDT path. The
new implementation preserves short-contour and frozen-baseline behavior,
survives all original-failure 500-measurement smoke cases, and passes the
original task88 parameters for a complete 500-burn plus 3000-measurement run.
No new unexplained numerical anomaly was observed. This is a validation of the
isolated development implementation only; it does not merge or promote that
implementation into the main branch.
