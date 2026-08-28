# Existing static/projector PfQMC average-sign audit

Scope: read-only audit of `/home/sunxr/PfQMC-main/reproduction/projector_kitaev`, excluding every `driven` path.  The companion CSV is the machine-readable table.  Values were recomputed from each available `measurements.csv` `sign` column, not copied from summaries.  Each 25 consecutive measurements is one delete-one jackknife unit, matching `validation_hpc_234/correlated_trotter_jackknife.py`; all equal-length seed streams in a row were concatenated before the jackknife.  `PBC` is `boundary=0`; `OBC` is `boundary=1`.

## L=6 validation_hpc_234 benchmark

| V | PBC <sign> | OBC <sign> |
| ---: | ---: | ---: |
| 0 | 1 (deterministic) | 1 (deterministic) |
| 2 | 0.912367 +/- 0.001683 | 0.905867 +/- 0.001717 |
| 4 | 0.468867 +/- 0.003645 | 0.608867 +/- 0.003230 |
| 6 | 0.448100 +/- 0.003624 | 0.579400 +/- 0.003419 |

Every entry is the complete 12-seed, 5000-measurement/seed benchmark grid (`theta=10`, `beta_trial=8`, `dt=0.1`, guard off).  V=0 has no updates/sign flips and is exactly one, so its error is zero.

## L>=18 static projector inventory

| L | V | BC | theta | beta_trial | dt | guard | seeds | measurements/seed | <sign> | error | status | source |
| ---: | ---: | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | --- | --- |
| 18 | 3.8 | PBC | 8 | 8 | 0.1 | on | 3 | 5000 | 0.044533 | 0.008011 | production/complete | L10_L18_observable_recheck |
| 18 | 4.0 | PBC | 8 | 8 | 0.1 | on | 3 | 5000 | 0.020267 | 0.007962 | production/complete | L10_L18_observable_recheck |
| 18 | 4.0 | PBC | 10 | 8 | 0.1 | on | 3 | 5000 | 0.004933 | 0.008128 | production/complete | L10_L18_observable_recheck |
| 18 | 4.0 | PBC | 12 | 8 | 0.1 | on | 3 | 5000 | 0.018800 | 0.008122 | production/complete | L10_L18_observable_recheck |
| 18 | 4.2 | PBC | 8 | 8 | 0.1 | on | 3 | 5000 | 0.024000 | 0.008213 | production/complete | L10_L18_observable_recheck |
| 18 | 4.4 | PBC | 8 | 8 | 0.1 | on | 3 | 5000 | 0.017200 | 0.008321 | production/complete | L10_L18_observable_recheck |
| 18 | 4.6 | PBC | 8 | 8 | 0.1 | on | 3 | 5000 | 0.016667 | 0.007989 | production/complete | L10_L18_observable_recheck |
| 18 | 3.8 | PBC | 8 | 8 | 0.1 | on | 1 | 5000 | 0.027200 | 0.013745 | incomplete | supp_guard_on_completion |
| 18 | 4.0 | PBC | 10 | 8 | 0.1 | on | 4 | 5000 | 0.001000 | 0.006914 | production/complete | supp_guard_on_completion |
| 18 | 4.0 | PBC | 12 | 8 | 0.1 | on | 4 | 5000 | -0.001800 | 0.007178 | production/complete | supp_guard_on_completion |
| 18 | 4.2 | PBC | 8 | 8 | 0.1 | on | 3 | 5000 | 0.014267 | 0.007985 | production/complete | supp_guard_on_completion |
| 18 | 4.4 | PBC | 8 | 8 | 0.1 | on | 3 | 5000 | 0.026800 | 0.008069 | production/complete | supp_guard_on_completion |
| 18 | 4.6 | PBC | 8 | 8 | 0.1 | on | 4 | 5000 | 0.012000 | 0.006977 | production/complete | supp_guard_on_completion |
| 18 | 3.8 | PBC | 8 | 8 | 0.1 | on | 1 | 300 | 0.026667 | 0.044947 | diagnostic/stress | short_diagnostics/V3.8 |
| 18 | 4.6 | PBC | 8 | 8 | 0.1 | on | 1 | 300 | -0.046667 | 0.054346 | diagnostic/stress | short_diagnostics/V4.6 |

`aligned_oracle/L18_V3.8` and `aligned_oracle/L18_V4.6` contain snapshots/phases but no QMC `measurements.csv` sign series, so they are listed as incomplete in the CSV rather than assigned an ensemble sign.  `regression_stress/raw/gaussian_exact` contains L=20,40,80 PBC/OBC artifacts (six files per L/BC), but no V and no average-sign observable; these are diagnostic/stress only, not PfQMC physics ensembles.

## Campaign provenance and classification

* `results/L10_L18_observable_recheck`: campaign settings and all 42 manifest rows specify `adaptive_guard=true`, threshold 0.1.  Its L=18 V=3.8, 4.0, 4.2, 4.4, 4.6 theta=8 groups each have the complete planned three seeds; V=4 theta=10 and 12 also each have all three planned seeds.  These are `production/complete` under this audit's manifest-completeness rule.
* `results/supp_guard_on_completion`: its 19-row manifest explicitly sets `adaptive_guard=true`.  V=4 theta=10/12 (4 seeds each), V=4.2 (3), V=4.4 (3), and V=4.6 (4) have every manifest row and are `production/complete`.  V=3.8 has only its sole manifest entry: it is retained but marked `incomplete` because it lacks the multi-seed set used for the interaction scan.
* The two 300-measurement short runs and all oracle/stress artifacts are not formal ensemble results.  In particular, do not cite their low/negative central values as physical average signs.

## Conclusions

At L=18 the sign is extremely small near V=4: the guard-on completed theta=10/12 sets are consistent with zero at about 0.007, and even the theta=8 recheck value is only 0.0203 +/- 0.0080.  Moving to V=3.8, 4.2, 4.4, or 4.6 does **not** show a statistically decisive recovery in the available L=18 datasets: the central values remain roughly 0.01--0.045 with 0.008 jackknife errors.  Thus a claimed recovery away from V=4 is not supported by these data.

The largest trustworthy, noncritical-point dataset is L=18 (PBC, guard on, complete) at V=3.8/4.2/4.4/4.6.  No complete L>=20 static PfQMC sign ensemble exists in the repository.  The L=20/40/80 regression artifacts and the one-seed short/oracle runs are diagnostic/incomplete only; they cannot be used as formal physical average-sign values.  No static L>=18 data were found at V=2, 3, 3.5, 5, or 6 (nor at the specifically requested V values beyond the listed 3.8, 4.0, 4.2, 4.4, 4.6).
