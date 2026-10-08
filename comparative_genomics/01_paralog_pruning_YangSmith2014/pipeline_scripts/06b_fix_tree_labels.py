#!/usr/bin/env python3
"""
Fix tip labels in .tre files: re-insert '@' after the 4-letter taxonID at the
start of each tip label (undoing the '@'->'_' replacement that happened during
tree building). After this, tree labels match the renamed FASTA headers
(taxonID@seqID) exactly, as required by the Yang & Smith pruning scripts.

The substitution only fires at label starts (right after '(' or ','), so
taxonID-like words inside the description part of a label are untouched.
Idempotent: labels that already contain '@' are not modified.

Usage: python3 06b_fix_tree_labels.py <tre_dir> <out_dir>
e.g.:  python3 scripts/06b_fix_tree_labels.py 2_trees 2_trees_at
"""
import os, re, sys

TAXA = "Acch|Aeco|Cjap|Dika|Mcav|Prvu|Rhvi|Sdu|Soly|Sxy|Vadu|Vvin"
PAT = re.compile(r"(?<=[(,])(%s)_" % TAXA)

def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]
    os.makedirs(dst, exist_ok=True)
    n, nlabels = 0, 0
    for fn in sorted(os.listdir(src)):
        if not fn.endswith(".tre"):
            continue
        with open(os.path.join(src, fn)) as f:
            text = f.read()
        text, k = PAT.subn(lambda m: m.group(1) + "@", text)
        with open(os.path.join(dst, fn), "w") as f:
            f.write(text)
        n += 1
        nlabels += k
    print("trees fixed : %d" % n)
    print("labels fixed: %d (expect ~35780)" % nlabels)

if __name__ == "__main__":
    main()
