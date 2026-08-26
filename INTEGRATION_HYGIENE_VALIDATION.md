# Integration staging hygiene validation

## Scope and verdict

Branch: `integration-staging-hygiene`

Parent: `integration-staging@c3ba9655a4905303b7742a0b4ed770db2a175b13`

Verdict: **PASS**.  This change is limited to build-system consistency,
standards-compliant JSON/output handling, explicit build provenance fields,
one UDT guard diagnostic correction, and tests.  It does not change the
ratio, sign-transport, Green-recovery, UDT thresholds, pivoting, or
factorization arithmetic.

All validation outputs and executables were written below
`/home/sunxr/PfQMC-hygiene-validation-20260826`; no production or PBS
executable was overwritten.

## Build consistency

The three supported build entry points now select C++17 and define
`PFQMC_SCALE_SAFE_UDT`:

| Entry point | Clean build | Feature evidence | Executable SHA-256 |
|---|---|---|---|
| Makefile | PASS | command log has `-std=c++17 -DPFQMC_SCALE_SAFE_UDT`; guard failure string present | `723cb039ee57065a928164b21238cf759092e232ea16baebded1e4f2c0b978f9` |
| projector `build.sh` | PASS | script has the same flags; JSON reports scale-safe/guard enabled; guard failure string present | `ba758c8e676f6f3d44530c889ad2e140b54387fac9b93bdb42347eb2d6889eaf` |
| CMake default | PASS | cache has `PFQMC_SCALE_SAFE_UDT:BOOL=ON`; compile database has the macro and `-std=c++17`; guard failure string present | `9472aba6b8f9fc443cce08aa800ffc436aad8d254cff3710820be9bf96a40e1a` |

CMake no longer assigns a compiler after `project()`.  The compiler is chosen
by the caller (`-DCMAKE_CXX_COMPILER=...` or `CXX`).  CMake uses an installed
`MKL::MKL` package when available and the same Intel `-mkl` mechanism as the
other build paths when this cluster's older oneAPI installation has no
`MKLConfig.cmake`.  PFAPACK is an explicit configurable root for all three
builds.  The projector script also supports an external `OUTPUT`, allowing a
clean build without writing into the source tree.

The historical GTest target remains available through
`-DPFQMC_BUILD_TESTS=ON`; it is opt-in so a default core build does not require
an otherwise unused GTest package.  A completely default configure recorded
`PFQMC_SCALE_SAFE_UDT=ON` and `PFQMC_BUILD_TESTS=OFF` and built successfully.

## JSON and zero-average-sign behavior

`projector_driver`, `projector_bins_driver`, and `static_guard_driver` share
`projector_json.h` for nullable finite-number emission and provenance.

At forced `average_sign == 0`, all three drivers:

* return success and emit `"status":"complete"`;
* preserve `"average_sign":0`;
* emit
  `"sign_reweighted_observables_status":"unresolved_zero_average_sign"`;
* emit `null` for `S_pi`, `S_pi_dq`, `R_cdw`, and their unresolved errors;
* never emit non-standard `nan`, `NaN`, `inf`, or `Infinity` JSON tokens.

The regression projector also uses the helper for an unsampled infinite UDT
guard margin, which is now standard JSON `null` rather than `inf`.

A strict parser using `parse_constant` rejection parsed 16 generated JSON
files.  The dedicated zero-sign schema check passed all three forced-zero
outputs.

The forced-zero hook is applied after measurement collection.  Measurement
CSV hashes are identical for normal and forced builds, and the clean
guard-off/guard-on trajectory is identical:

`63c54609064502fd0085abdffc36cde32c5e663c0e865fd6a13d3819c73d72e8`

This single hash is shared by `bins_raw.csv`, `bins_zero_raw.csv`,
`static_guard_raw.csv`, and `static_guard_zero_raw.csv`.

## Provenance schema

Relevant projector JSON now contains:

* `scale_safe_udt: true`;
* `udt_rank_loss_guard: true`;
* `udt_rank_loss_guard_bits: 45`;
* `udt_orthogonality_gate: 1e-6`;
* `left_recovery_enabled: false` for these default-off runs;
* `condition_aware_ratio_enabled: false` (not integrated on this branch);
* `source_commit` and `executable_sha256`.

The last two values are `null` unless a build system supplies reliable
compile-time values.  They are intentionally not guessed from a mutable
working tree or executable path.

## UDT diagnostic correction

The final-Q orthogonality failure now reports the maximum lost-bits value for
the current factorization.  Previously it used the process-global diagnostic
maximum, which could describe an earlier factorization.  The formal 45-bit
rank-loss limit, 32-bit staged precheck, `1e-6` final orthogonality gate, pivot
order, MGS operations, D/T bookkeeping, and fail-closed decisions are
unchanged.

Synthetic guard regression: five cases passed.  The `n=12, +/-20` control
completed with finite unitary output; `n=24, +/-40`, `n=12, +/-500`, and
`n=52, +/-1500/2000` failed closed as expected.  No finite strongly
nonunitary factor was returned.

## Numerical regression

| Test | Result |
|---|---|
| tiny complete HS enumeration | PASS; max direct Green relative error `7.11e-15` |
| V=0 Gaussian | PASS |
| L10 control | PASS; Green diagnostic max `3.18e-15`, sign mismatch 0, UDT trigger 0 |
| task88 short | PASS; Green diagnostic max `3.05e-14`, sign mismatch 0, UDT trigger 0 |
| task92 short | PASS; Green diagnostic max `9.18e-13`, sign mismatch 0, UDT trigger 0 |
| 3000 local flips | PASS; anomalies 0 |
| all-boundary Green | PASS; nonfinite 0, max error `2.18e-10` |
| left-recovery replay | PASS; nonfinite 0, post-recovery max error 0 |
| reality symmetry | PASS |
| scale-safe UDT normal/stress | PASS |

The guard-off and guard-on measurement trajectories in the forced-zero test
case are byte-identical.  Normal QMC runs report zero UDT rank-loss guard
triggers.

## Protected running jobs

No file used by PBS 153474 or the optimized-v3 follow-up was modified.  Their
executable hashes are checked again after the commit in the final handoff.
