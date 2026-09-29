#!/bin/bash
# Step: paralogy pruning (Yang & Smith 2014 pipeline Step 6).
# RT (root-to-tip, needs outgroups) and MI (maximum inclusion, no outgroups needed)
# are both run; compare yields and keep the better one.
# Needs ANY python 2.7 interpreter (see 03_trim_mask.sh); override with PY2=<path>.
#
# usage: bash 04_prune.sh <treeDIR> <prunedDIR_RT> <prunedDIR_MI> <min_ingroup_taxa> <kit_scripts_parent>
# e.g.:  bash scripts/04_prune.sh 2_trees 4_pruned_RT 4_pruned_MI 8 .
set -eu
TR=${1:?tree dir}; OUT_RT=${2:?RT out dir}; OUT_MI=${3:?MI out dir}; MINT=${4:-8}; KIT=${5:-.}
KIT="$(cd "$KIT" && pwd)"   # absolute path robustness
mkdir -p "$OUT_RT" "$OUT_MI"

# ---- pick a python-2 interpreter (same logic as 03_trim_mask.sh) ----
PY2BIN=""
for c in "${PY2:-}" python2.7 python2; do
    [ -n "$c" ] || continue
    if command -v "$c" >/dev/null 2>&1 \
       && "$c" -c 'import sys; sys.exit(0 if sys.version_info[0]==2 else 1)' 2>/dev/null; then
        PY2BIN="$c"; break
    fi
done
if [ -z "$PY2BIN" ] && command -v python >/dev/null 2>&1 \
   && python -c 'import sys; sys.exit(0 if sys.version_info[0]==2 else 1)' 2>/dev/null; then
    PY2BIN=python
fi
if [ -z "$PY2BIN" ]; then
    source activate ys 2>/dev/null || conda activate ys 2>/dev/null || true
    PY2BIN=python
fi
echo "Using python2: $PY2BIN ($("$PY2BIN" --version 2>&1))"

echo "== RT pruning (outgroups: Vvin, Soly) =="
"$PY2BIN" "$KIT/ys_scripts/prune_paralogs_RT.py" "$TR" .tt.mm "$OUT_RT" "$MINT" "$KIT/scripts/taxon_list.txt"

echo "== MI pruning (fallback, no outgroups needed) =="
"$PY2BIN" "$KIT/ys_scripts/prune_paralogs_MI.py" "$TR" .tt.mm 10 10 "$MINT" "$OUT_MI"

echo "RT ortho trees: $(ls "$OUT_RT"/*.tre 2>/dev/null | wc -l)"
echo "MI ortho trees: $(ls "$OUT_MI"/*.tre 2>/dev/null | wc -l)"
