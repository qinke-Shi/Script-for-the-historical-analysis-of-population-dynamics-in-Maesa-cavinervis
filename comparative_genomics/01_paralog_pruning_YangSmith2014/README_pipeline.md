# MCMCTree sensitivity analysis — Yang & Smith (2014) tree-based paralogy pruning

Goal: rebuild the MCMCTree input from **paralog-pruned 1-to-1-derived orthologs**
following the reviewer-recommended pipeline (Yang & Smith 2014,
`phylogenomic_dataset_construction`), then re-run MCMCTree with the **identical**
clock/calibration settings as the main analysis, and compare divergence times.

Main analysis (already done): OrthoFinder 3.0.1b1 `-M msa` SpeciesTreeAlignment
(310,647 aa, STAG tree-based orthology) + MCMCTree (PAML 4.9j, Poisson+Γ4,
independent-rates clock, 5 TimeTree calibrations).
Sensitivity analysis: 2,612 low-copy OGs -> gene trees -> trim/mask tips ->
RT & MI paralogy pruning -> ortholog alignments -> concatenated supermatrix
(column occupancy >= 0.3) -> MCMCTree (same ctl, only seqfile changed).

## 0. Upload & layout

From local Windows machine:

    scp E:\Maesa\orthofinder\tre\mcmc\sensitivity_kit.zip xubo@planthunter:/data_forrest/xubo/shiqinke/Maesa_pop/compare_gene/1_or_input/OrthoFinder/Results_Jan30/

On server:

    cd /data_forrest/xubo/shiqinke/Maesa_pop/compare_gene/1_or_input/OrthoFinder/Results_Jan30/
    unzip sensitivity_kit.zip        # creates ./ys_scripts ./scripts mcmc_sens.ctl mcmctree_input.tree

    # working dirs
    mkdir -p 1_renamed 2_trees 4_pruned_RT 4_pruned_MI 5_ortho_fa_RT 5_ortho_fa_MI 6_mcmctree

Check tools (before make_env):

    which iqtree2 iqtree raxml mafft    # note what exists

## 1. Environment (user-level, no sudo)

    bash scripts/make_env.sh
    # - writes ~/.condarc (Tsinghua mirror; fixes the ConnectionResetError)
    # - creates conda env "ys" (python 2.7, only stdlib needed)
    # - installs iqtree only if no iqtree2/iqtree binary exists

Test:  `conda activate ys && python --version`  (must print Python 2.7.x)

## 2. Candidate OG list, rename alignments -> taxonID@seqID and merge

All commands run from inside `Results_Jan30/` (the kit was unzipped here).

    # candidate OGs: >= 8 of 12 species with exactly one copy (expect 2612)
    awk -F'\t' 'NR>1{n=0; for(i=2;i<=13;i++) if($i==1) n++; if(n>=8) print $1}' \
        Orthogroups/Orthogroups.GeneCount.tsv > og_list.txt
    wc -l og_list.txt

    python3 scripts/01_rename_and_merge.py MultipleSequenceAlignments 1_renamed 1_renamed/all.fa og_list.txt

Verify output line `Taxon IDs found : 12 -> Acch, Aeco, ...` (12 species).
If the count is not 12, stop and check headers: `head -2 MultipleSequenceAlignments/OG0000000.fa`.

## 3. Gene trees (IQ-TREE2, WAG+G4, no bootstrap — as in Y&S homolog step)

    bash scripts/02_build_trees.sh 1_renamed 2_trees 24

~2,612 small trees (only candidate OGs); check `2_trees/fail.log` (empty = fine).
If iqtree2 is missing after step 1, install: `conda install -n base -c bioconda iqtree -y`.

## 4. Trim tips + mask same-taxon monophyletic tips  (Y&S Step 5, python2)

    conda activate ys
    bash scripts/03_trim_mask.sh 2_trees .

## 5. Paralogy pruning  (Y&S Step 6: RT primary, MI fallback)

    bash scripts/04_prune.sh 2_trees 4_pruned_RT 4_pruned_MI 8 .

    # RT = root-to-tip pruning, rooted by outgroups Vvin + Soly
    # MI = maximum inclusion, no outgroups required (used if RT yield is low)

Compare yields; choose the mode with more ortholog trees as the primary.

## 6. Extract ortholog alignments (short taxonID headers)

    python3 scripts/05_extract_seqs_py3.py 1_renamed/all.fa 4_pruned_RT 5_ortho_fa_RT 8
    python3 scripts/05_extract_seqs_py3.py 1_renamed/all.fa 4_pruned_MI 5_ortho_fa_MI 8

## 7. Concatenate -> MCMCTree phylip (gap-fill missing taxa, occupancy >= 0.3)

    python3 scripts/06_concat_for_mcmctree.py 5_ortho_fa_RT mcmctree_input.tree \
            6_mcmctree/sens_supermatrix.phy --occupancy 0.3

(check per-taxon gap fractions in the printed summary; if RT yield was poor,
re-run with 5_ortho_fa_MI.)

## 8. Re-run MCMCTree — identical settings to main analysis

    cd 6_mcmctree
    cp ../mcmc_sens.ctl .
    mcmctree mcmc_sens.ctl          # ctl is identical to the main mcmctree.ctl
                                    # except seqfile/outfile names; usedata=3,
                                    # so generate out.BV and re-run exactly as
                                    # you did for the main analysis (mcmc6)

## 9. Compare & report

- `sens_out.txt` -> FigTree.tre vs main analysis tree: check node ages for the
  5 calibrated nodes and all other nodes (correlation / scatter).
- Copy back for me: `6_mcmctree/sens_supermatrix.phy` header + stats printed in
  step 7, `sens_out.txt`, and the main-vs-sensitivity age comparison.
- These go into the response letter as: "a sensitivity analysis in which
  paralogs were pruned following Yang & Smith (2014) ... produced nearly
  identical divergence-time estimates (r = ...)".

## Notes

- All Yang & Smith python scripts here are the ORIGINAL unmodified py2 scripts
  from `bitbucket.org/yangya/phylogenomic_dataset_construction`; run them in
  the `ys` env. Helper scripts (01/05/06, shell wrappers) are python3/bash.
- Trees use `-m WAG+G4` without bootstrap: pruning uses topology only, and this
  matches the Y&S homolog-tree procedure (RAxML WAG, no bootstrap needed).
- RT taxon list (`scripts/taxon_list.txt`): IN = Acch Aeco Cjap Dika Mcav Prvu
  Rhvi Sdu Sxy Vadu; OUT = Vvin Soly (Vitales/Solanales as outgroups to Ericales).
