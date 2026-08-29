# Preintegration without proposal-safety

## Scope and isolation

- Work directory: `/home/sunxr/new-pfqmc-preintegration-no-proposal-safety`
- Branch: `preintegration-no-proposal-safety`
- Base: `4d79bb7ac2b734bc10667d917b914ce316a39a39`
- External build/results directory:
  `/home/sunxr/new-pfqmc-preintegration-no-proposal-safety-results`
- No main/develop branch was modified or merged.
- Proposal-safety/#3 commit `d291f374ec328e251a79d67b75fa47311d95f2f6`
  was not integrated. Its worktree and task 153850 were not modified.
- No 150-seed replay, PBS job, production scan, or new zero-temperature work
  was started.
- `condition_aware_ratio=false` and `left_recovery=false`; the production
  ratio/recovery files have zero diff from the MP trust-gate base.

## Production ancestry and integrated commits

The following production commits are direct ancestors and were not
cherry-picked again:

- `0c2df511dff146eb15c9cf395c1f54863422ea9e`
- `89809695c0dea7cb725a0645a0a48bd2f84f378b`
- `14969585b0aee9b24b7bf221834041f11d1ca0df`

The requested C/D/E commits were cherry-picked in order:

| Scope | Source commit | Integrated commit |
|---|---|---|
| C: retained I/O completion gate | `c65c3348bad6c780897585a9a3628e0ab715073a` | `18c101499caa60b9135eb14bc4fa902c914f95c6` |
| D: PBC structural fixes | `3eec30c8d4dbfb9ba6009c6a7c18479ccc4fe734` | `d5a12d5b774693cacbc4bcda610f175e92385a47` |
| E: driven Dexp comparator | `ad1d419e22e2449b81942a8586b069c39e0f0eba` | `0c5d401d9221b8989cafbfac43858f8069cada1d` |

All C/D/E semantic source and test files are byte-equivalent to their source
commits after integration.

## Conflict record

One add/add conflict occurred while applying D: C and D both added a root
`regression_summary.csv`, with incompatible schemas. No core source conflict
occurred.

The minimal resolution preserved both datasets without combining or rewriting
test results:

- C summary -> `IO_COMPLETION_GATE_REGRESSION_SUMMARY.csv`
- D summary -> `PBC_STRUCTURAL_REGRESSION_SUMMARY.csv`

Only the corresponding file references in the two historical reports were
updated. Core behavior and regression semantics are unchanged.

## Regression results

### Existing MP trust-gate full regression

The existing `mp_z2_trust_gate_fix_20260829/build.sh` and
`run_regressions.sh` completed on the login node. All 12 existing groups pass:

- generic complex;
- exact projector and exact finite enumeration;
- OBC+hs1;
- zero-sign JSON;
- sign interface;
- local flips;
- right-boundary Green;
- reality symmetry;
- integration QMC;
- driven static comparison;
- driven smoke.

Key gates:

- OBC+hs1 `z2_average_sign=1`: PASS.
- OBC+hs1 `shadow_trajectory_match=true`: PASS.
- OBC+hs1 `mp_correction_count=0`: PASS.
- Driven smoke `mp_correction_count=0`: PASS.
- Driven smoke `status=complete`: PASS.

### PBC structural gates

- Shared L=3 PBC construction fails fast: PASS.
- Projector L=3 PBC exits 2 with the odd-L guard: PASS.
- Driven L=3 PBC exits 2 with the odd-L guard: PASS.
- L=3 OBC construction remains valid: PASS.
- Even-L L=6 PBC `max ||B B_inv-I||=0`: PASS.
- Even-L L=6 PBC right-boundary Green: 180/180 finite,
  `max_error=2.1837046465355791e-10`: PASS.
- L=4 PBC four-bond direct-Wick interaction-energy relative error
  `2.341683895620198e-15` (tolerance `1e-13`): PASS.
- L=4 OBC three-bond interaction-energy relative error
  `8.3969449632371878e-16`: PASS.

### Retained I/O completion gate and trajectory

- Normal single and dual retained streams: PASS.
- CSV row count, final newline, and parseability: PASS.
- Bin and time-series `/dev/full`: return 2, no complete output, explicit
  flush failure: PASS.
- Unwritable bin and time-series paths: return 2, no complete output: PASS.
- Equivalent MP scale-safe baseline/integrated trajectory SHA-256:
  `956bc2a3a173f0b91c8f3509ad1bc9b640b93e68974990b9fcb56e0433508fac`
  on both binaries: PASS.
- Equivalent bin CSV SHA-256:
  `19792cbcc66058d135c61a589760e4f4ae90539f27b82e212e9abc8e07e7ffe0`
  on both binaries: PASS.
- Physics JSON fields are identical. The raw C helper reports its provenance
  equality row as FAIL only because the expected `source_commit` changed from
  the frozen baseline to the preintegration HEAD; after excluding
  `source_commit`, runtime, and output path metadata, no field differs.

An initial comparison against C's older non-scale-safe build was discarded as
non-equivalent. The authoritative trajectory comparison above uses binaries
built with the same MP scale-safe configuration.

### Driven Dexp comparator

- Synthetic nonzero `Dexp` reconstruction: PASS; fixed reconstruction,
  onePlusInv comparison, direct comparison, and logdet errors are all zero.
- Pathological L=6 comparator: PASS under E's existing scale-aware criteria.
- `max_ratio_complex_abs_error=6.123669014735936e-06` versus
  `max_green_abs_error=4.238120499768675e-06`.
- `max_udt_solve_residual=5.131389506368776e-14`.
- Baseline/fixed comparator trajectory projection is identical, SHA-256
  `e875961ea6e27db15b701edb2c733d23ae563c08efc7e6c025e69e0d24aa3bd8`.

## Actual production diff

Relative to the frozen base, the only production-path changes are:

- `inc/kitaevChain.h`: odd-L PBC guard and PBC boundary interaction energy;
- `reproduction/driven_kitaev/driven_driver.cpp`: retained CSV completion
  gate.

E changes validation/comparator code only. There are no changes in `src/`,
`inc/pfqmc.h`, `inc/operator.h`, `inc/spinless_tV.h`, `main.cpp`, production
ratio/recovery logic, or zero-temperature paths.

## Files intentionally waiting for #3

The following proposal-safety core/interface files remain outside this branch
and must be handled only after #3 is finalized:

- `inc/operator.h`
- `inc/pfqmc.h`
- `inc/spinless_tV.h`
- `src/pfqmc.cpp`
- `reproduction/projector_kitaev/projector_json.h`
- `reproduction/projector_kitaev/real_z2_raw_checker_fix_20260827/projector_json.h`
- `reproduction/projector_kitaev/real_z2_raw_checker_fix_20260827/projector_mp_z2_oracle.h`
- `reproduction/projector_kitaev/real_z2_raw_checker_fix_20260827/projector_real_z2_driver.cpp`

This branch is therefore a clean C+D+E preintegration baseline ready for the
later, separately validated #3 integration.
