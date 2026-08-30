# Pure-projector production semantics review

## Scope and provenance

- Branch: `pure-state-projector-production-semantics-review`
- Base Phase 3C commit: `c901cb50d88c244cc5c6973e0936627899a849d3`
- Semantics implementation commit:
  `62a00f892878d9c8143dce6be2efbf9dffeacfbc`
- Small-system HPC validation-script commit:
  `d5bd8ec5bd9e6a1615ec368e49e946cbb6ae913e`
- Phase 4 was not started.
- The Phase 3C MP same-proposal fallback algorithm, all numeric trust thresholds,
  condition-aware ratio, left recovery, and mutating raw/MP checkpoints were not
  changed.
- Existing HPC results and jobs were not overwritten.  The four redundant local
  L=12 processes corresponding to completed job 153923 were stopped at the
  user's request.

## Burn and measurement units

The old driver called `step()` once per `burn` unit and once per measurement, so
both were proposal counts.  The production default is now `--run-units sweeps`.
A complete sweep visits every HS variable once and alternates forward/reverse
order between sweeps.  `--measurement-stride N` advances N configured units
between retained measurements.  The old behavior remains available only as
`--run-units proposals --measurement-stride 1`.

The completion JSON now records `hs_variable_count`, `burn_proposals`,
`burn_sweep_equivalent`, `measurement_stride`, `measurement_stride_unit`,
`measurement_stride_proposals`, `measurement_proposals`,
`measurement_sweep_equivalent`, and total `proposal_count`.

For the Phase 3B benchmark geometry (`theta=L`, `dt=0.1`), the number of HS
variables is `20 L^2` for even-L PBC and `20 L(L-1)` for OBC.  Thus the former
`burn=1000` meant only the following sweep equivalents:

| L | PBC HS variables | old PBC burn sweeps | OBC HS variables | old OBC burn sweeps |
|---:|---:|---:|---:|---:|
| 6 | 720 | 1.388889 | 600 | 1.666667 |
| 12 | 2880 | 0.347222 | 2640 | 0.378788 |
| 18 | 6480 | 0.154321 | 6120 | 0.163399 |

With the new default, `burn=1000` explicitly means 1000 full sweeps.  No large
production scan was launched with the corrected interpretation.

## Observable normalization

All pure-projector driver, dense/ED oracle, Phase 2 validation, and Phase 3A
validation paths now use one shared definition,

`S(q) = sum_ij exp(i q (i-j)) <rho_i rho_j> / L^2`.

The former duplicated implementations divided by L.  In a frozen L=4 legacy
trajectory, both `S_pi` and `S_pi_dq` changed from 0.25 to 0.0625, exactly the
required additional factor 1/L.  `R_CDW = 1-S_pi_dq/S_pi` remained exactly zero,
as required.  Energy, Z2, acceptance, and the HS/RNG trajectory hashes were
unchanged.

## Physical fermion parity

For block-ordered coordinates `(gamma_0,...,gamma_{L-1},gamma'_0,...,gamma'_{L-1})`,
the reported physical parity is now

`(-i)^L (-1)^[L(L-1)/2] Pf(G)`.

The code records the internal `Pf(-iG)` sign and the block-reordering sign
separately.  `GaussianTrialState::fermionParity()` and production measurements
return the physical label.  The ED program constructs the selected dense-Fock
sector, measures the parity of the actual dense state, verifies that it matches
the explicit trial policy, and reports that measured value rather than echoing
the policy.  Local ED checks gave parity -1 for both L=4 PBC and L=6 OBC test
sectors.

## Corrected checkerboard models

The shared checkerboard bond-count helper now defines:

- even-L PBC: `{L/2, L/2}` (all L bonds, including `(L-1,0)`);
- OBC: `{L/2, (L-1)/2}` (exactly L-1 bonds).

The production driver, Phase 2 dense/exact/ED validation, and Phase 3A fast/slow
validation all use this same helper.  This replaces the incorrect historical
Phase 2/3A counts `{L/2,(L-1)/2}` for PBC and `{(L+1)/2,L/2}` for OBC.

## Small safety fixes

- Pure-projector Green reconstruction always records transpose-skew and diagonal
  residuals.  Tests can explicitly promote them to a fail-closed structure gate;
  production does not silently add a new trust threshold in this review.  The
  physical-parity estimator independently fails closed when these structures are
  untrusted.
- The rank-0 fast-ratio branch now requires a successful checked Green rebuild
  and propagates its overlap diagnostics.
- After every completed proposal, `currentWeight().green` is synchronized to the
  live stack Green.  If synchronization is unavailable it is explicitly cleared,
  never left as a stale pre-update matrix.
- A Phase 3C unit-test fixture contained a non-Gaussian `a==b` pseudo-rotation;
  the new structural diagnostic found it, and the fixture was corrected to a
  genuine plane rotation.

## Regression status

Local checks completed so far:

- Phase 1 core: 6/6 PASS.
- Phase 2 core: 5/5 PASS.
- Phase 3A core: 4/4 PASS.
- Phase 3C core/safety: 9/9 PASS.
- Production semantics core and driver/I/O contract: PASS.
- Phase 3A validation: 1000 proposal comparisons and 2000 trajectory steps;
  maximum ratio relative error `2.3900310206755053e-15`, maximum Green relative
  error `7.745503499262228e-15`, Z2 mismatches 0, trajectory mismatches 0.
- Legacy L=4 new/old trajectory: HS hash `6214351880218834211` and RNG hash
  `11030203536659588101` identical.
- L=6 OBC+hs1 short control: average Z2 1, bin error 0, no trust alarm/fallback.
- L=4/L=6 standalone dense ED: PASS with actual dense-Fock parity -1.
- Frozen Task 4: ratio `0.70741914065425615`, accept, 320-digit
  convergence, final trajectory hash `949208207496548183`.
- Frozen Task 8: ratio `0.0026573889530685144`, reject, 320-digit
  convergence, final trajectory hash `14005514076608941581`.
- Frozen Task 4/8 MP fallback failures: 0; each continued for 64 proposals.
- Phase 3B driver core: fast/audit HS hashes identical for L=2 and L=4;
  block sizes 1/2/4/8 produced the same trajectory and `S_pi`; normal,
  zero-sign, unwritable-path, and `/dev/full` I/O cases all PASS.
- PBC structural guards: odd-L PBC rejection PASS, even-L `B B^-1`
  residual 0, L=4 PBC four-bond energy error `2.3417e-15`, and OBC
  three-bond energy error `8.3969e-16`.
- Corrected Phase 2 dense/exact/ED validation: 8/8 PASS using 4 local CPU
  threads.  Maximum dense weight error `7.4997e-15`, dense Green error
  `1.6405e-14`, and determinant-identity error `1.4048e-15`.
  L=4 ED maximum observable difference was `5.0379e-16`; L=6 was
  `1.7764e-15`.  The largest slow-MC deviations from exact enumeration were
  `0` for average sign, `0.00122493` for `S_pi`, and `0.0134519` for `R_CDW`,
  within the existing statistical envelope.

HPC validation job 153930 was submitted while the local run proceeded, but it
remained queued with CPU time zero.  After the local 8/8 PASS it was cancelled
as redundant; no other HPC job was changed.

## Build

HPC/reference build flags remain C++17 with `PFQMC_SCALE_SAFE_UDT`, Intel MKL,
Eigen, PFAPACK, and OpenMP for the Phase 2 validation executable.  The production
driver executable SHA256 and final Git commits are recorded in the final handoff
section after the verified build.

Verified local `-O2` production driver SHA256:
`d9fd7a40807c68eb20b168781305f459df9741af2d6df5198b5b92c68a344d63`.
The corrected Phase 2 validation executable SHA256 was
`2e4163eacbdb4c9be54aa112001e74ccee8e3c0961845ad2e4af1c54a9bb3ad1`.
