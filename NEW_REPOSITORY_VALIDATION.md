# New unified repository validation

Source under test: `main@45c8af962ce3a23daa956ac73ecad63a7dec3d71`
(stable code `55fee175`; the difference is documentation only).

## Builds

| Entry | Result | Executable SHA-256 |
|---|---|---|
| Makefile clean build | PASS | `9d9b1ca5d89fd3c5845263d6bb1c1079679c77d0155fe2da60ba6724e39a1cbf` |
| CMake clean build, scale-safe default ON | PASS | `64ab571d325dc33373ac179523584bbfa4c22bc5ebdf616418b281275b6d19c0` |
| projector `build.sh` | PASS | `42e29fbc6b9150ed7a60625b988c78154661320bacc7df7b79caa30919e2b741` |
| consolidated core drivers | PASS | 14 drivers built in ignored `build/core-regression/bin` |

The first CMake attempt omitted `-DEIGEN3_INCLUDE_DIR` and stopped at missing
`Eigen/Dense`; the corrected clean configure/build passed. No source change was
required.

## Regression job

PBS `153579.mgt` completed successfully using one CPU. The durable runner
executed:

- tiny complete-HS enumeration: PASS;
- V=0 Gaussian: PASS;
- same-contour ED: PASS;
- L10 control: PASS;
- task88/task92 short smokes: PASS;
- 3000 local-update property trials: PASS (`anomalies=0`);
- all-boundary Green: PASS (`nonfinite=0`, maximum relative error
  `1.0674e-10`);
- optional left-sweep recovery replay: PASS;
- normal scale-safe UDT: PASS;
- reality symmetry: PASS;
- zero-average-sign projector/bins/static-guard: PASS;
- strict JSON: 17 normal files and 3 forced-zero files parsed successfully.

L10/task88/task92 all report `scale_safe_udt=true`,
`udt_rank_loss_guard=true`, 45 guard bits, zero guard triggers,
`left_recovery_enabled=false`, and `condition_aware_ratio_enabled=false`.
Maximum observed lost bits were 8.1794 (task88) and 9.4033 (task92), leaving
more than 35 bits of margin to the 45-bit fail gate.

Synthetic UDT outcomes were exactly as required:

- n12, ±20: normal success, no trigger;
- n24, ±40: fail closed;
- n12, ±500: fail closed;
- n52, ±1500 and ±2000: fail closed;
- no finite strongly nonunitary U was emitted.

Forced zero-average-sign JSON completed with `average_sign=0`, status
`unresolved_zero_average_sign`, and sign-reweighted observables encoded as
standard JSON `null`. Raw bins/static-guard measurement CSVs were unchanged by
the forced-zero output path.

## Stable equivalence

`git diff 55fee175..main` contains only root documentation. Core source and test
content are byte-identical to the validated integration hygiene commit, so the
normal numerical path is unchanged. Build/test artifacts are ignored and are
not committed.

Verdict: **PASS — unified stable main is suitable as the base for future work.**
