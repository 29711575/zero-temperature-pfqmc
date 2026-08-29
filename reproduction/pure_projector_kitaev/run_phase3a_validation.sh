#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
here=$(cd "$(dirname "$0")" && pwd)
output=${1:-"$here/phase3a_results"}
threads=${PURE_PROJECTOR_PHASE3A_THREADS:-4}
binary="$here/build_phase3a/phase3a_validation"
test -x "$binary"
mkdir -p "$output"
source_commit=$(cd "$root" && git rev-parse HEAD)
executable_sha=$(sha256sum "$binary" | awk '{print $1}')
OMP_NUM_THREADS="$threads" MKL_NUM_THREADS=1 \
  "$binary" "$output" "$source_commit" "$executable_sha" "$threads"
