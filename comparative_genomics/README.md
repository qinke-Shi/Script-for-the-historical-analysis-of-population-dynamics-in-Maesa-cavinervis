# Comparative genomics: orthology inference, paralog pruning and molecular dating

This directory archives the scripts and configuration files for the phylogenomic
part of the study: orthology inference across the 12 species, the Yang & Smith
(2014) tree-based paralog-pruning sensitivity analysis, and the MCMCTree
divergence dating (main analysis + sensitivity analysis), as described in
Methods 2.3 of the manuscript.

## 0. Orthology inference (`00_orthology_inference/`)

| File | Purpose |
|---|---|
| `standardize_headers.py` | Unify the protein FASTA headers of the 12 species to `>species\|proteinID` |
| `extract_longest_pep.py` | Keep only the longest isoform per gene, producing the `*_longest_transcript_pep.fa` inputs |
| `run_ortho.sh` | OrthoFinder **v3.0.1b1** with `-S diamond -M msa` (`-t 24 -a 20`) over the 12 longest-transcript protein sets (MAFFT is called internally by `-M msa`) |
| `Orthogroups_for_concatenated_alignment_705.txt` | The **705 orthogroups** (≥ 8/12 species with exactly one copy) auto-selected by OrthoFinder; their `SpeciesTreeAlignment` (310,647 aa) is the input of the main MCMCTree analysis |

## 1. Paralog-pruning sensitivity analysis (`01_paralog_pruning_YangSmith2014/`)

Sensitivity analysis following the tree-based paralogy-pruning pipeline of
**Yang & Smith (2014)**; see [`README_pipeline.md`](01_paralog_pruning_YangSmith2014/README_pipeline.md)
for the step-by-step server workflow.

- `og_list.txt` — the **2,612 candidate low-copy orthogroups** used as input
  (orthogroups with ≥ 8 of 12 species having exactly one copy).
- `taxon_list.txt` — ingroup/outgroup (`IN`/`OUT`) designation of the 12 taxa
  (10 ingroups incl. *M. cavinervis*, 2 outgroups: *V. vinifera*, *S. lycopersicum*).
- `pipeline_scripts/` — numbered workflow wrappers:

  | Script | Step |
  |---|---|
  | `01_rename_and_merge.py` | OrthoFinder per-OG alignments → `taxonID@geneID` headers, merged `all.fa` |
  | `02_build_trees.sh` | One IQ-TREE2 gene tree per OG (WAG+G4, unbootstrapped) |
  | `03_trim_mask.sh` | Trim long tips + collapse monophyletic same-taxon tips (python 2.7) |
  | `04_prune.sh` | **MI** (maximum-inclusion) and **RT** (rooted-tree) paralogy pruning |
  | `05_extract_seqs_py3.py` | Extract pruned ortholog sequences |
  | `06_concat_for_mcmctree.py` | Concatenate to PHYLIP (column occupancy ≥ 0.3), headers matched to the MCMCTree tree file |
  | `06b_fix_tree_labels.py` | Harmonise tree labels with the sequence file |
  | `make_env.sh` | One-time server environment setup (python-2.7 env + IQ-TREE2) |

- `ys_scripts/` — the original Yang & Smith (2014)
  `phylogenomic_dataset_construction` scripts invoked by the wrappers
  (`prune_paralogs_MI.py`, `prune_paralogs_RT.py`, `trim_tips.py`,
  `mask_tips_by_taxonID_genomes.py`, `tree_utils.py`, `newick3.py`, `phylo3.py`).

The pruned MI and RT matrices were re-dated with MCMCTree using settings
**identical** to the main analysis; all 11 node ages recovered 95% credibility
intervals fully overlapping those of the main analysis (node-by-node comparison
provided as a supplementary table in the manuscript), indicating that paralog
retention does not affect the divergence-time conclusions.

## 2. Molecular dating (`02_molecular_dating/`)

| File | Purpose |
|---|---|
| `mcmctree.ctl` | Main analysis control file (PAML **4.9j** `mcmctree`: **Poisson + Γ4** rate variation — with `seqtype = 2` (amino acids) `model = 0` means Poisson; the `0:JC69` comment in the ctl applies to nucleotide data only — independent-rates clock `clock = 2`, approximate-likelihood final run `usedata = 2` reading the precomputed `in.BV`, birth–death prior `BDparas = 1 1 0.1`) |
| `mcmctree_input.tree` | 12-taxon guide tree with the five TimeTree-based calibration constraints `B(...)` used in all MCMCTree runs |
| `run_mcmctree.sh` | Job wrapper (adapt the absolute paths before reuse) |
| `mcmc_sens.ctl` | First-pass template (`usedata = 3`): estimates the branch-rates file `out.BV`, which is renamed to `in.BV` for the final runs |
| `mcmc_sens_MI.ctl` / `mcmc_sens_RT.ctl` | Final MI / RT sensitivity runs (`usedata = 2`, reading the estimated `in.BV`); identical to the main analysis except for the pruned alignment |

This setup corresponds to the manuscript description (Methods 2.3): MCMCTree under
the approximate likelihood method, an independent-rates clock, and the
**Poisson+G4** model.

## Requirements

- OrthoFinder v3.0.1b1 (DIAMOND, `-M msa`; MAFFT called internally)
- IQ-TREE2 (gene trees for the pruning pipeline)
- PAML 4.9j (`mcmctree`)
- Python 3 for the wrappers; **Python 2.7** for the original Yang & Smith
  scripts (`PY2=<path>` override, see `03_trim_mask.sh` / `04_prune.sh`)

## Notes

- Scripts are archived byte-identical to the versions executed on the analysis
  server; absolute paths inside the job wrappers (`run_ortho.sh`,
  `run_mcmctree.sh`, `README_pipeline.md`) must be adapted before reuse.
- References: Yang & Smith (2014), *Systematic Biology*; Emms & Kelly (2017),
  *Molecular Biology and Evolution* (STAG/OrthoFinder).
