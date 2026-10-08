#!/usr/bin/env python3
"""
Python-3 equivalent of Yang & Smith's write_ortholog_fasta_files.py
(no Biopython dependency).

Reads the ortholog trees produced by prune_paralogs_RT/MI (*.tre), takes the
tip labels ("taxonID@seqID"), and writes one fasta per ortholog containing the
corresponding sequences with SHORT taxonID headers (e.g. ">Mcav").

Usage:
  python3 05_extract_seqs_py3.py <all.fa> <ortho_tree_dir> <out_dir> <MIN_TAXA>
"""
import os, sys

def read_fasta(path):
    """Yield (header, seq) pairs; seq has no whitespace."""
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

def tip_labels(newick_line):
    """Return tip labels from a newick string.
    Tip labels start right after '(' or ',' and end at ':', ',', ')';
    labels appearing after ')' are internal-node support values, skipped."""
    labels = []
    token = ""
    in_tip = False
    for ch in newick_line:
        if ch == "(":
            token, in_tip = "", True
        elif ch == ",":
            if token and in_tip:
                labels.append(token)
            token, in_tip = "", True
        elif ch == ")":
            if token and in_tip:
                labels.append(token)
            token, in_tip = "", False
        elif ch == ":":
            if token and in_tip:
                labels.append(token)
            token, in_tip = "", False
        elif ch in " \t\n;":
            continue
        else:
            token += ch
    if token and in_tip:
        labels.append(token)
    return labels

def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    all_fa, tree_dir, out_dir, min_taxa = sys.argv[1], sys.argv[2].rstrip("/"), sys.argv[3], int(sys.argv[4])
    os.makedirs(out_dir, exist_ok=True)

    # seqDICT[taxonID][fullID] = seq
    seqDICT = {}
    for head, seq in read_fasta(all_fa):
        taxonID = head.split("@")[0]
        seqDICT.setdefault(taxonID, {})[head] = seq
    print("Loaded %d sequences from %d taxa" % (sum(len(v) for v in seqDICT.values()), len(seqDICT)))

    n_out, n_skip = 0, 0
    for fn in sorted(os.listdir(tree_dir)):
        if not fn.endswith(".tre"):
            continue
        with open(os.path.join(tree_dir, fn)) as fh:
            line = fh.readline()
        labels = tip_labels(line)
        if len(labels) < min_taxa:
            n_skip += 1
            continue
        out = os.path.join(out_dir, fn[:-4] + ".fa")
        with open(out, "w") as fh:
            for lab in labels:
                taxonID = lab.split("@")[0]
                if taxonID not in seqDICT or lab not in seqDICT[taxonID]:
                    sys.exit("ERROR: tip '%s' in %s not found in all.fa "
                             "(check renaming step / trailing whitespace)" % (lab, fn))
                fh.write(">%s\n%s\n" % (taxonID, seqDICT[taxonID][lab]))
        n_out += 1
    print("Written %d ortholog fasta files (%d skipped < %d taxa)" % (n_out, n_skip, min_taxa))

if __name__ == "__main__":
    main()
