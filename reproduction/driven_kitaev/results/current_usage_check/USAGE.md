# Current zero-T driven PfQMC usage

## CLI

```text
driven_driver L V0 Vf drive_rate theta_init beta_trial dt delta mu boundary burn measurements seed [adaptive_guard guard_threshold] [timeseries.csv stride]
```

Example verified in `current_usage_check`:

```bash
/tmp/driven_usage_check 4 0 4 2 1 1 0.1 1 0 0 10 20 88001
```

## Parameter rules

- There are 13 required positional parameters.
- `adaptive_guard` and `guard_threshold` must be supplied together.
- `timeseries.csv` and `stride` must be supplied together.
- `theta_init / dt` must be an integer.
- `(Vf - V0) / (drive_rate * dt)` must be an integer.
- The sign of `drive_rate` must point from `V0` to `Vf`.
- `seed` changes the Monte Carlo chain.

## Formal outputs

Use `result.json` as the formal per-task summary and `bin_records.csv` as the
raw pooled-bin input. The checked JSON includes `S_pi`, `S_pi_err`,
`S_pi_dq`, `S_pi_dq_err`, `R_cdw`, `R_cdw_err`, `average_sign`,
`average_sign_err`, `acceptance`, `runtime_seconds`, `max_sign_imag`,
`max_observable_imag`, `sign_corrections`, `adaptive_guard`,
`multiprecision_fallback`, and rebuild/guard diagnostics. `stdout` is captured
as `result.json`; diagnostics/errors are in `stderr.log`.

`R_cdw` uses the production ratio-of-means estimator.

## Minimal HPC PBS example

```bash
#!/bin/bash
#PBS -N driven_usage
#PBS -q workq
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:10:00
#PBS -j oe
set -euo pipefail
cd /home/sunxr/PfQMC-main
source /opt/ohpc/pub/apps/intel/oneapi/setvars.sh >/dev/null
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
out=reproduction/driven_kitaev/results/current_usage_check/pbs
mkdir -p "$out"
export PFQMC_BIN_RECORDS_PATH="$PWD/$out/bin_records.csv"
/tmp/driven_usage_check 4 0 4 2 1 1 0.1 1 0 0 10 20 88004 >"$out/result.json" 2>"$out/stderr.log"
```

This environment and directory pattern was checked with PBS smoke `153380.mgt`.

后续批量驱动脚本应以本文件定义的 CLI 为准，不再依赖旧版 driver 调用方式。
