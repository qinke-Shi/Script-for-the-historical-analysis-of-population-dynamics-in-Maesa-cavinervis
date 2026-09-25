#!/usr/bin/env bash
set -euo pipefail

FSC_BIN="${1:-fsc28}"
MAX_JOBS="${2:-1}"
BOOT_ROOT="${3:-M25_boot}"
NUM_SIMS="${4:-100000}"

if ! command -v "$FSC_BIN" >/dev/null 2>&1; then
  echo "Cannot find fsc binary: $FSC_BIN" >&2
  exit 1
fi

if [[ ! -d "$BOOT_ROOT" ]]; then
  echo "Cannot find bootstrap result directory: $BOOT_ROOT" >&2
  exit 1
fi

if ! [[ "$MAX_JOBS" =~ ^[0-9]+$ ]] || [[ "$MAX_JOBS" -lt 1 ]]; then
  echo "MAX_JOBS must be a positive integer" >&2
  exit 1
fi

mapfile -t REPS < <(find "$BOOT_ROOT" -mindepth 1 -maxdepth 1 -type d -name 'M25_boot_[0-9]*' | sort)
if [[ "${#REPS[@]}" -eq 0 ]]; then
  echo "No replicate directories found under $BOOT_ROOT" >&2
  exit 1
fi

run_one() {
  local rep_dir="$1"
  echo "Running ${rep_dir}"
  (
    cd "$rep_dir"
    "$FSC_BIN" \
      -t M25_boot.tpl \
      -e M25_boot.est \
      -n "$NUM_SIMS" \
      -N "$NUM_SIMS" \
      -L 40 \
      -m \
      -M \
      -0 \
      -B12 \
      -c0 \
      --initValues M25_boot.pv \
      -q
  )
}

if [[ "$MAX_JOBS" -eq 1 ]]; then
  for rep_dir in "${REPS[@]}"; do
    run_one "$rep_dir"
  done
else
  export FSC_BIN NUM_SIMS
  export -f run_one
  printf '%s\n' "${REPS[@]}" | xargs -I{} -P "$MAX_JOBS" bash -lc 'run_one "$1"' _ {}
fi
