#!/bin/bash
# Step: trim long tips + mask monophyletic same-taxon tips
# (Yang & Smith 2014, pipeline Step 5; genome-annotation version of masking).
# Needs ANY python 2.7 interpreter: system python2.7/python2, an existing
# conda env's python (2.x), or a new "ys" env (fallback). Override with PY2=<path>.
#
# usage: bash 03_trim_mask.sh <treeDIR> <kit_scripts_parent>
# e.g.:  bash scripts/03_trim_mask.sh 2_trees .
set -eu
TR=${1:?tree dir}; KIT=${2:-.}
KIT="$(cd "$KIT" && pwd)"   # absolute: the script cd's into $TR before calling ys_scripts

# ---- pick a python-2 interpreter ----
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
    # last resort: create/use a dedicated conda env (needs network)
    source activate ys 2>/dev/null || conda activate ys 2>/dev/null || true
    PY2BIN=python
fi
echo "Using python2: $PY2BIN ($("$PY2BIN" --version 2>&1))"

cd "$TR"
# 1) trim tips sticking out (> relative cutoff AND >10x sister) or > absolute cutoff
"$PY2BIN" "$KIT/ys_scripts/trim_tips.py" . .tre 10 10
# 2) collapse monophyletic tips from the same taxon, keep shortest terminal branch
"$PY2BIN" "$KIT/ys_scripts/mask_tips_by_taxonID_genomes.py" .
echo "done. .tt: $(ls *.tt 2>/dev/null | grep -v '\.mm$' | wc -l)  .tt.mm: $(ls *.tt.mm 2>/dev/null | wc -l)"
