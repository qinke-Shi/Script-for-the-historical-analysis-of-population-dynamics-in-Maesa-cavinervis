#!/bin/bash
# Build one gene tree per renamed OG alignment with IQ-TREE2 (WAG+G4, no bootstrap,
# matching the Yang & Smith 2014 homolog-tree step which uses unbootstrapped trees
# for paralogy pruning).
#
# usage: bash 02_build_trees.sh <renamedDIR> <treeDIR> <threads>
# e.g.:  bash scripts/02_build_trees.sh 1_renamed 2_trees 24
set -u
SRC=${1:?renamed dir}; DST=${2:?tree out dir}; NT=${3:-24}
mkdir -p "$DST"

build_one() {
    f="$1"; b=$(basename "$f" .fa)
    if [ -s "$DST/$b.tre" ]; then return; fi          # skip finished
    iqtree2 -s "$f" -m WAG+G4 -pre "$DST/$b" -nt 1 -quiet >/dev/null 2>&1 \
      && mv "$DST/$b.treefile" "$DST/$b.tre" \
      || echo "$b" >> "$DST/fail.log"
}
export -f build_one; export DST

ls "$SRC"/OG*.fa | xargs -P "$NT" -I{} bash -c 'build_one "$@"' _ {}

echo "done. trees: $(ls "$DST"/*.tre 2>/dev/null | wc -l)  failed: $(wc -l < "$DST/fail.log" 2>/dev/null || echo 0)"
