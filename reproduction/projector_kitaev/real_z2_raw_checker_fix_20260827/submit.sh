#!/usr/bin/env bash
set -euo pipefail;base=$(cd "$(dirname "$0")"&&pwd);cd "$base";mkdir -p results;j=$(qsub production.pbs);f=$(qsub -W depend=afterok:"$j" finalize.pbs);printf 'production=%s\nfinalize=%s\n' "$j" "$f"|tee submitted_jobs.txt
