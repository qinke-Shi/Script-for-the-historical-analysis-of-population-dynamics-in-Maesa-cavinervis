#!/usr/bin/env python3
"""
Rename OrthoFinder per-OG alignments to Yang & Smith (2014) format:
  header "Species_stem|geneID"  ->  ">taxonID@geneID"
and merge the renamed OG files into one all.fa (input for extract_seqs).

Usage:
  python3 01_rename_and_merge.py <OG_aln_dir> <out_renamed_dir> <out_all.fa> [og_list.txt]

og_list.txt (optional) limits processing to candidate orthogroups (one
"OGxxxxxxx" per line), e.g. OGs with >= 8 of 12 species having exactly one copy.

Example (run from Results_Jan30/):
  python3 scripts/01_rename_and_merge.py MultipleSequenceAlignments 1_renamed 1_renamed/all.fa og_list.txt
"""
import os, re, sys

def clean_id(s, maxlen=60):
    s = re.sub(r"[^A-Za-z0-9_.@-]", "_", s.strip())
    return s[:maxlen] if maxlen else s

def parse_header(h):
    """Return (taxonID, seqID) from an OrthoFinder alignment header."""
    h = h.strip()
    if "|" in h:
        stem, sid = h.split("|", 1)
    else:
        stem, sid = h, h
    taxonID = stem.split("_")[0][:8]          # e.g. Mcav_longest_transcript_pep -> Mcav
    return taxonID, clean_id(sid)

def main():
    argv = sys.argv[1:]
    if len(argv) not in (3, 4):
        sys.exit(__doc__)
    og_dir, out_dir, all_fa = argv[0], argv[1], argv[2]
    og_dir = og_dir.rstrip("/")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(all_fa)) or ".", exist_ok=True)

    wanted = None
    if len(argv) == 4:
        with open(argv[3]) as fh:
            wanted = set(l.strip() for l in fh if l.strip())
        print("OG list: %d candidates" % len(wanted))

    taxons, nfiles, nseqs, nlines = set(), 0, 0, 0
    files = sorted(f for f in os.listdir(og_dir) if f.startswith("OG") and f.endswith(".fa"))
    if wanted is not None:
        have = set(f[:-3] for f in files)
        missing = wanted - have
        if missing:
            print("WARNING: %d listed OGs not found in %s (e.g. %s)"
                  % (len(missing), og_dir, ", ".join(sorted(missing)[:5])))
        files = [f for f in files if f[:-3] in wanted]

    with open(all_fa, "w") as allout:
        for fn in files:
            og = fn[:-3]
            src = os.path.join(og_dir, fn)
            dst = os.path.join(out_dir, og + ".fa")
            with open(src) as fin, open(dst, "w") as fout:
                for line in fin:
                    # lstrip: some source headers carry a leading space
                    if line.lstrip().startswith(">"):
                        taxonID, seqID = parse_header(line.lstrip()[1:])
                        taxons.add(taxonID)
                        fout.write(">%s@%s\n" % (taxonID, seqID))
                        nseqs += 1
                    else:
                        s = line.strip()
                        if s:                          # drop blank/junk lines
                            fout.write(s + "\n")       # always newline-terminated
            nfiles += 1
            # merge into all.fa (each gene belongs to exactly one OG, no duplicates)
            with open(dst) as fin:
                for line in fin:
                    nlines += 1
                    allout.write(line if line.endswith("\n") else line + "\n")

    print("OG files renamed : %d" % nfiles)
    print("Sequences        : %d" % nseqs)
    print("Merged all.fa    : %d lines" % nlines)
    print("Taxon IDs found  : %d -> %s" % (len(taxons), ", ".join(sorted(taxons))))

if __name__ == "__main__":
    main()
