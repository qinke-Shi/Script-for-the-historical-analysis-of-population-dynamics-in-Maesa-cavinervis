#!/usr/bin/env bash
set -euo pipefail

BOOT_ROOT="${1:-M25_boot}"

if [[ ! -d "$BOOT_ROOT" ]]; then
  echo "Cannot find bootstrap result directory: $BOOT_ROOT" >&2
  exit 1
fi

shopt -s nullglob
count=0
for rep_dir in "$BOOT_ROOT"/M25_boot_*; do
    [[ -d "$rep_dir" ]] || continue
    cp M25_boot.tpl M25_boot.est M25_boot.pv "$rep_dir"/
    count=$((count + 1))
done

if [[ "$count" -eq 0 ]]; then
  echo "No bootstrap replicate directories found under $BOOT_ROOT" >&2
  exit 1
fi

echo "Prepared $count replicate directories under $BOOT_ROOT"
