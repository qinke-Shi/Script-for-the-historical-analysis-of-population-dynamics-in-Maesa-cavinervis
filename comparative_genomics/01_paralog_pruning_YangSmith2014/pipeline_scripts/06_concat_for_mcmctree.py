#!/usr/bin/env python3
"""
Concatenate per-ortholog fasta alignments (short taxonID headers, e.g. ">Mcav",
as written by 05_extract_seqs_py3.py) into a sequential PHYLIP file for
MCMCTree.

- Each ortholog partition keeps its own alignment width (widths differ between
  OGs); taxa missing from a partition are filled with gaps ('-').
- Columns occupied (non-gap) by fewer than --occupancy of the 12 taxa are
  dropped (equivalent to phyutility MIN_COLUMN_OCCUPANCY filtering,
  cf. Yang & Smith 2014 pipeline Step 7).
- Sequence headers are taken from the MCMCTree tree file (long names such as
  "Mcav_longest_transcript_pep") so that seqfile and treefile match; the short
  taxonID (part before the first '_') is used to look up sequences.

Usage:
  python3 06_concat_for_mcmctree.py <ortho_fa_dir> <mcmctree_input.tree> \
          <out.phy> [--occupancy 0.3]

Example:
  python3 scripts/06_concat_for_mcmctree.py 5_ortho_fa_MI mcmctree_input.tree \
          6_mcmctree_MI/sens_supermatrix.phy --occupancy 0.3
"""
import os, sys

def read_fasta(path):
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip().split()[0], []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)

def taxon_order_from_tree(tree_path):
    with open(tree_path) as fh:
        nw = fh.read()
    # tips = labels before ':' that are not preceded by ')' (support values);
    # also skip MCMCTree calibration annotations like 'B(0.95, 1.10)'.
    taxa, token, in_tip, depth = [], "", False, 0
    i = 0
    while i < len(nw):
        ch = nw[i]
        if ch == "'":                      # skip quoted MCMCTree annotations
            j = nw.find("'", i + 1)
            i = j + 1 if j != -1 else len(nw)
        elif ch == "(":
            token, in_tip, depth = "", True, depth + 1
        elif ch == ")":
            if token and in_tip and token not in taxa:
                taxa.append(token)
            token, in_tip, depth = "", False, depth - 1
        elif ch == ",":
            if token and in_tip and token not in taxa and depth > 0:
                taxa.append(token)
            token, in_tip = "", True
        elif ch == ":":
            if token and in_tip and token not in taxa and depth > 0:
                taxa.append(token)
            token, in_tip = "", False
        elif ch in " \t\n;":
            pass
        else:
            token += ch
        i += 1
    if token and in_tip and token not in taxa and depth > 0:
        taxa.append(token)
    return taxa

def main():
    argv = sys.argv[1:]
    occ = 0.3
    if "--occupancy" in argv:
        i = argv.index("--occupancy")
        occ = float(argv[i + 1])
        del argv[i:i + 2]
    if len(argv) != 3:
        sys.exit(__doc__)
    fa_dir, tree_path, out_phy = argv
    taxa = taxon_order_from_tree(tree_path)
    print("Taxon order (%d): %s" % (len(taxa), ", ".join(taxa)))

    # long tree-tip name -> short taxonID (used in ortho fasta headers)
    short = {}
    for t in taxa:
        s = t.split("_")[0]
        if s in short.values():
            sys.exit("ERROR: taxonID collision when shortening tree tip names: %s" % s)
        short[t] = s

    # gather alignments; each partition has its own width
    genes, part_w = [], []
    for fn in sorted(os.listdir(fa_dir)):
        if not fn.endswith(".fa"):
            continue
        d = {}
        for head, seq in read_fasta(os.path.join(fa_dir, fn)):
            d[head] = seq
        if not d:
            continue
        ws = set(len(s) for s in d.values())
        if len(ws) != 1:
            sys.exit("ERROR: unequal alignment widths within %s: %s" % (fn, sorted(ws)))
        genes.append(d)
        part_w.append(ws.pop())
    ngenes = len(genes)
    raw_len = sum(part_w)
    print("Genes: %d   concatenated raw length: %d aa" % (ngenes, raw_len))

    # concatenate with gap filling
    concat = {t: [] for t in taxa}
    for g, w in zip(genes, part_w):
        for t in taxa:
            concat[t].append(g.get(short[t], "-" * w))
    concat = {t: "".join(v) for t, v in concat.items()}

    # column occupancy filter
    ntax = float(len(taxa))
    keep = [i for i in range(raw_len)
            if sum(1 for t in taxa if concat[t][i] != "-") / ntax >= occ]
    if not keep:
        sys.exit("ERROR: no columns pass occupancy %.2f" % occ)
    for t in taxa:
        concat[t] = "".join(concat[t][i] for i in keep)
    final_len = len(keep)

    os.makedirs(os.path.dirname(os.path.abspath(out_phy)), exist_ok=True)
    # format identical to the main analysis (SpeciesTreeAlignment.phy):
    # name on its own line, sequence wrapped at 60 chars on following lines
    with open(out_phy, "w") as out:
        out.write(" %d %d\n" % (len(taxa), final_len))
        for t in taxa:
            seq = concat[t]
            out.write(t + "\n")
            for j in range(0, len(seq), 60):
                out.write(seq[j:j + 60] + "\n")

    print("Final length: %d aa (%d columns dropped, occupancy >= %.2f)"
          % (final_len, raw_len - final_len, occ))
    print("Per-taxon gap fraction:")
    for t in taxa:
        print("  %-28s %.3f" % (t, concat[t].count("-") / float(final_len)))
    print("Wrote %s" % out_phy)

if __name__ == "__main__":
    main()
