#!/usr/bin/env bash
set -euo pipefail

FSC_BIN="${1:-fsc28}"

if ! command -v "$FSC_BIN" >/dev/null 2>&1; then
  echo "Cannot find fsc binary: $FSC_BIN" >&2
  exit 1
fi

"$FSC_BIN" -i M25_boot.par -n 100 -j -m -s0 -x -I -q
