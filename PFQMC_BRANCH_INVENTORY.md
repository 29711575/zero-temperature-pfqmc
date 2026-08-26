# PfQMC branch and repository inventory

Snapshot: 2026-08-26 (Asia/Shanghai). This inventory is read-only with respect
to every pre-existing repository. In particular, the worktrees and artifacts
used by PBS 153474 and 153573 were not checked out, rebuilt, or modified.
PBS 153573 completed during the inventory window; its results were not analyzed.

## Repository inventory

| Path | Checked-out branch / HEAD | Base | Status at inventory | Role |
|---|---|---|---|---|
| `/home/sunxr/PfQMC-main` | `fix-condition-aware-ratio` / `d8eebc643fc32ac57c4a77589a52303094feda9a` | local `main=2ae3bb7d81c4c4153b773b81dfe897405a94adab` | dirty: tracked ratio/oracle changes plus untracked diagnostic/results trees; protected because 153474 is running | canonical repository plus the v1/v2 diagnostic worktree; **not modified** |
| `/home/sunxr/PfQMC-condition-ratio-v3` | `optimize-condition-aware-ratio-v3` / `243e660b8341afe84b7354d87cc69cd877eeb492` | remote `main=2ae3bb7...` | clean at inventory | optimized v3 replay used by PBS 153573; completed, not analyzed, **not modified** |
| `/home/sunxr/PfQMC-core-regression-tests` | `optimize-udt-rank-loss-guard` / `5b21b110113bfe3f3f104baf8157cc4e45cf4a3d` | remote `main=2ae3bb7...` | clean | sanitizer/core tests, UDT guard/fallback/Householder audits |
| `/home/sunxr/PfQMC-integration-prep` | `integration-prep` / `5dcf46b77ab86d2d3729068e8449d396326db034` | remote `main=2ae3bb7...` | clean | cleaned UDT-guard port and validated default-off left-sweep recovery |
| `/home/sunxr/PfQMC-driven-vc4-fts` | `driven-vc4-fts` / `c18fa796eb60e619e7cb35514b8e0c8e578bcefa` | remote `main=2ae3bb7...` | untracked benchmark directories | separate driven-project work, outside this static-projector integration |
| `/home/sunxr/PfQMC-v5-sign-compare` | `compare-v5-sign-definition` / `c0dd27401c74f69f186e315a7967447113978e3e` | local `main=2ae3bb7...` | staged/untracked comparison reports/tools | separate sign-definition study; explicitly left untouched at user request |
| `/home/sunxr/PfQMC-integration-staging` | `integration-staging` / created from `2ae3bb7...` | `main=2ae3bb7...` | isolated integration work only | unified staging created by this work |

Tags visible in the main-line clones include `pre-scale-safe-20260824` and
`scale-safe-udt-validated-20260824`; the latter identifies `2ae3bb7...`.

## Development DAG

```text
pre-scale-safe-20260824
        |
        v
2ae3bb7  main: validated scale-safe UDT (A/B/C/D PASS)
        |
        +--> 557b730  near-zero sign audit
              |
              v
            8490e0c  deterministic sign/Green replay
              |
              v
            77c7adf  ratio-reality / Green-gate audit
              |
              v
            1e7c674  Green-recovery prototype
              |
              v
            d8eebc6  MP ratio oracle / condition-aware base
              |\
              | +--> 243e660  optimized condition-aware ratio v3 (153573 complete; analysis pending)
              |
              +--> b33dc7d --> 5dcf46b  integration-prep guard + left recovery
              |
              +--> 837552a --> 62a1d6a --> ... --> ea9e9f09 --> 5b21b110
                       core tests   zero-sign fix        safe/optimized guard

2ae3bb7 --> integration-staging
              +-- 8887c5f  production-only optimized UDT rank-loss guard
              +-- f29af53  zero-average-sign output fix
              +-- 834ff59  optional, default-off left-sweep Green recovery
```

The graph is provenance-oriented: some diagnostic repositories have no local
`main` branch, so their base is the recorded remote/main commit. No history was
rewritten and no source branch was deleted.

## Branch/commit purpose and disposition

| Branch / commit | Purpose | Validation status | Production relevance / staging action |
|---|---|---|---|
| `main@2ae3bb7` | mantissa+base-2 exponent UDT, exponent-aware QR, equilibrated factorization/solve, scale-safe products | A/B/C/D PASS; full original-failure task88 3000-measurement run passed | staging base; retained unchanged |
| `diagnose-near-zero-sign@557b730` | same-configuration sign/Green event capture | audit complete | diagnostic provenance only; do not port |
| `diagnose-sign-green-replay@8490e0c` | deterministic sign drift and Green replay | complete; second-stage-denominator blind spot rejected; rebuild restored Green | diagnostic provenance only |
| `diagnose-ratio-reality-green-gate@77c7adf` | rank-4/sequential ratio reality and structural gate | complete; rank-4 equals sequential; strong complex ratio follows bad Green | diagnostic provenance only |
| `fix-green-recovery-prototype@1e7c674` | same-configuration full-rebuild recovery | deterministic recovery demonstrated; threshold alone not a final production definition | only the later validated left-sweep, default-off subset is staged |
| `diagnose-mp-ratio-oracle@d8eebc6` | boundary/order and converged MP ratio oracle | complete: boundary/order correct, rank-4=sequential, r^2=Q, relevant weight proven real | oracle retained as a test; no production MP ratio logic staged |
| condition-aware v1 / PBS 153474 | early adaptive MP replay | still running at snapshot | wait; executable/PBS/raw protected |
| condition-aware v2 / PBS 153475 | difficult L12,V5 replay | numerical PASS but too expensive (4487 fallbacks; MP dominated walltime) | analysis provenance only; not staged |
| `optimize-condition-aware-ratio-v3@243e660` / PBS 153573 | decouple Green rebuild from MP, start MP at 160, decision-margin/audit reuse | completed during final inventory check; numerical result not analyzed here | await dedicated analysis; no v3 production logic staged |
| `candidate-udt-rank-loss-guard@ea9e9f09` | 45-bit fail-closed rank/orthogonality guard | SAFE TO PORT; real-QMC triggers 0 | superseded by optimized guard |
| `optimize-udt-rank-loss-guard@5b21b110` | staged checks with final adjoint-orthogonality gate | SAFE TO PORT; synthetic failure cases fail closed; matched QMC unchanged | production-only subset staged as `8887c5f` |
| `62a1d6a` | valid JSON when `average_sign==0` | PASS | staged as `f29af53` |
| `integration-prep@5dcf46b` | guard clean port, left recovery, unified regression preparation | guard and deterministic recovery validation complete | minimal left recovery staged default-off; durable tests consolidated |
| `test-core-regression@837552a` and descendants | sanitizer, boundary Green, local algebra, UDT stress, reality symmetry | core numerical tests pass; sanitizer runtime availability remained follow-up | durable source/scripts retained; raw/PBS output excluded |
| `diagnose-udt-householder-feasibility@6be96b7` | distinguish MGS instability from effective rank loss | mixed regime; real QMC has about 35.9-bit margin to 45-bit guard | no Householder or MP fallback staged |

## Running-job provenance held read-only

| Job | Source | Executable provenance | State at snapshot |
|---|---|---|---|
| 153474 | condition-aware v1 in the protected `/home/sunxr/PfQMC-main` diagnostic worktree | executable `/home/sunxr/PfQMC-main/reproduction/projector_kitaev/regression_stress/offcritical_sign_scan/condition_aware_ratio/bin/condition_aware_ratio_driver`, SHA-256 `0f2a6a67f6813611db863ac17695e6d47b0414fb13f501b0d4f4a3096a3c0431`; `replay.pbs` SHA-256 `b25e1d6189cb306d6bf8a10beea40fae88f8f9eb9b2e8b63d348192b15a61d82`; `manifest.csv` SHA-256 `03e88fab1ca1cc35f4632681f25531f6e53fdeede25637c40a253d4109e23d58` | running on one CPU at final check; no recompile or overwrite performed |
| 153573 | `optimize-condition-aware-ratio-v3@243e660b8341afe84b7354d87cc69cd877eeb492` | `/home/sunxr/PfQMC-condition-ratio-v3-build/bin/condition_aware_ratio_driver_v3`, SHA-256 `f97549ac6874fb1b8825de1f9e42cc99ada8b1cc38d548023cd8bdd9fb0f1b59`; `/home/sunxr/PfQMC-condition-ratio-v3-build/replay_v3.pbs`, SHA-256 `12eef8c046e44e60460ac257c8741d9f024ee9b3d8e9efaeda2b33c3bfddb79e` | status marker says complete; not analyzed; no file in its source/build/output chain was changed |

Running-job/result artifacts are intentionally not copied into staging.

## Archive candidates (do not delete yet)

After 153474 is complete and both 153474/153573 analyses are committed, the
near-zero, replay, ratio-gate, early Green-recovery, UDT MP-fallback, and
Householder-feasibility branches can be marked archival. They remain valuable
for provenance and minimal reproducers. `PfQMC-v5-sign-compare` and
`PfQMC-driven-vc4-fts` are separate research lines, not archive candidates of
this integration operation.
