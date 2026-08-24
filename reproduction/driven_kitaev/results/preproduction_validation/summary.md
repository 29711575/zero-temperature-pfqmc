# Preproduction validation

Baseline backup: /home/sunxr/PfQMC-main_preproduction_baseline_20260821.tar.gz

| Test | Key result | Status |
|---|---|---|
| Clean build | driven_driver, static_contour_compare, driven_fastupdate_check, static_contact_test, capture_boundary_check; mpiicpc -std=c++14 | PASS |
| Build wrapper | build.sh with mpiicpc/MKL/Eigen/Boost 1.70 produced driven_driver_build_sh | PASS |
| Static regression V=0 | operator, Green, S_pi, S_pi_dq, R_cdw differences all 0 | PASS |
| Static regression V=4 | operator, Green, S_pi, S_pi_dq, R_cdw differences all 0 | PASS |
| Identity contact | S_pi=1/(4L)=0.041666666666666664 at L=6; absolute difference 0 | PASS |
| Capture defense | invalid boundary and pointer-without-boundary rejected; valid midpoint capture succeeds | PASS |
| L=4 ED spot | Vf=3, R=2; z(S_pi)=-0.579, z(S_pi_dq)=0.641, z(R_cdw)=-0.601 | PASS |
| Bin collector | main JSON vs bin collector max difference <= 5.56e-16 for center/error fields | PASS |
| Observable centers | current vs post-cleanup baseline max center difference <= 1.11e-16 | PASS |
| Fast update L=4 | max Green=4.64e-14, ratio=0, post-reset=0 | PASS |
| Fast update L=6 | max Green=1.87e-13, ratio=0, post-reset=0 | PASS |
| L=8 smoke | finite; sign=0.820; acceptance=0.899; max sign imag=2.42e-11; max observable imag=3.10e-15; max green skew=2.75e-13 | PASS |
| Script checks | build/PBS shell syntax passed; production wrappers set bin path and code version | PASS |

## Changes in scope

- Added rightSweep capture-boundary validation and requested-capture completion checks.
- Added projector contour input/count assertions without changing contour ordering.
- Changed static_contour_compare to use StructureFactorCDW directly, and added an independent Identity contact test.
- Added per-bin production records and delete-one-contiguous-bin jackknife errors for S_pi, S_pi_dq, and R_cdw. Central values remain ratio-of-means.
- Added JSON provenance and numerical diagnostics with a 1e-7 production failure tolerance.
- Added a bin collector and capture/contact validation entry points.
- Marked legacy static_note_projector PBC adaptive_guard=true results as historical/non-production.

The observable center values did not change within floating-point roundoff. The repository is ready for production scans using the updated wrappers; no large production task was submitted here.
