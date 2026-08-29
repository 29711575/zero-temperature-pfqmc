#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
here=$(cd "$(dirname "$0")" && pwd)
output=${1:-"$here/phase2_results"}
threads=${PURE_PROJECTOR_PHASE2_THREADS:-8}
binary="$here/build/phase2_validation"
test -x "$binary"
mkdir -p "$output"
source_commit=$(git -C "$root" rev-parse HEAD)
executable_sha=$(sha256sum "$binary" | awk '{print $1}')
OMP_NUM_THREADS="$threads" MKL_NUM_THREADS=1 \
  "$binary" "$output" "$source_commit" "$executable_sha"
