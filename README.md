# Script for the historical analysis of population dynamics in *Maesa cavinervis*

This repository archives the fastsimcoal2 model definitions and analysis pipeline used to
infer the historical demography (bottlenecks, divergence times and gene flow) of the five
genetic groups of *Maesa cavinervis*, as described in the manuscript.

## Overview of the analysis

- **Data**: 677,266 intergenic SNPs (93 individuals assigned to five genetic groups; two
  NRD samples excluded, 91 individuals used in the final SFS). The dataset was built from
  the LD-pruned Dataset 2 (963,287 SNPs) by re-annotating with SnpEff and retaining
  `intergenic_region` sites (677,266 SNPs).
- **fastsimcoal2 settings**: SNPs treated as 150 bp partially linked DNA loci,
  mutation rate mu = 2.43e-8 per generation (5-year generation time), 100,000 simulations
  and 40 EC-M cycles per run (`-n 100000 -L 40 -M`).
- **Model comparison strategy**: all 27 models were each run once for an initial screen;
  the five best models (M25, M16, M15, M27, M14) were then each re-run with 100 independent
  optimizations to confirm the ranking. **M25** (bottlenecks in Groups 5, 6 and 1, with
  gene flow G2<->G3 and G1<->G3 centred on Group 3) was the best-supported model
  (Delta AIC >= 7,409 against all alternatives); the best run reached
  lnL = -7,164,399.384.
- **Uncertainty**: parameter uncertainty was quantified by a parametric bootstrap:
  100 joint SFS were simulated under the M25 maximum-likelihood point estimates, and all
  21 parameters were re-estimated from each replicate. Bootstrap re-estimates were
  initialised from starting values drawn at random within the prior ranges. The 2.5%/50%/97.5%
  percentiles across the 100 replicates are reported as the 95% confidence intervals
  (manuscript Table S17).

## Repository layout

```
.
├── models/                  # 27 fastsimcoal2 model definitions (M01-M27)
│   ├── M01.tpl ... M27.tpl  # template files: historical model (demography/gene flow)
│   └── M01.est ... M27.est  # parameter files: prior definitions for all parameters
├── batch_100runs/
│   └── run_new_models_linux.sh   # SGE (qsub) script: initial screen of all 27 models,
│                                 # then 100 independent runs for the 5 best models
│                                 # (12 threads per model, 60 threads total)
└── bootstrap_ci/            # M25 parametric-bootstrap CI pipeline
    ├── M25_boot.par         # parameter file anchored to the ML point estimates
                             # (narrow priors, ±5-15% around the ML values)
    ├── 01_generate_bootstrap_100.sh   # simulate 100 joint SFS replicates
    │                        #   fsc28 -i M25_boot.par -n 100 -j -m -s0 -x -I -q
    ├── 02_prepare_reestimate_inputs.sh # assemble per-replicate input directories
    ├── 03_run_reestimate_100.sh       # re-estimate all 21 parameters from each of the
    │                        #   100 simulated SFS (xargs -P32 parallel fsc28 jobs)
    ├── 04_summarize_ci.py             # summarize 100 .bestlhoods files into
    │                        #   bootstrap_ci_summary.tsv (median, 2.5%, 97.5% per parameter)
    └── m25ci2.sh            # SGE wrapper that runs steps 03-04 with randomized
                             # initial values (final, reported version of the pipeline)
```

## Model definitions (M01-M27)

Demes follow the five genetic groups G1-G5; G5, G6, G3, G2, G1 in the template files
correspond to Groups 5/6/3/2/1 with sample sizes 22/50/22/34/30. Migration rates are
forward rates per generation; bottleneck parameters `RES_*` are size-reduction factors,
`TBOT_E` the bottleneck duration and `INC*` post-bottleneck growth/recovery sizes.

| Model | Demographic scenario |
|---|---|
| M01 | Baseline: no migration, no bottleneck (reference for model comparison) |
| M02 | Bottleneck in G5+G6, no migration |
| M03 | No bottleneck, bidirectional migration G2<->G3 |
| M04 | No bottleneck, G1<->G5 migration (collapsed to prior lower bound; not supported) |
| M05 | No bottleneck, G5<->G6 migration (marginal convergence) |
| M06 | No bottleneck, G1<->G6 migration |
| M07 | Bottleneck G5+G6, G1<->G5 migration (collapsed to lower bound) |
| M08 | Bottleneck G5+G6, G5<->G6 migration |
| M09 | Bottleneck G5+G6, G1<->G6 migration |
| M10 | No bottleneck, G5<->G3 migration |
| M11 | Bottleneck G5+G6, G5<->G3 migration |
| M12 | Bottlenecks G5+G6+G1, G1<->G5 migration (marginal convergence) |
| M13 | Bottlenecks G5+G6+G1, G5<->G6 migration |
| M14 | Bottlenecks G5+G6+G1, G1<->G6 migration |
| M15 | No bottleneck, G2<->G3 and G1<->G3 migration |
| M16 | Bottleneck G5+G6, G2<->G3 and G1<->G3 migration |
| M17 | No bottleneck, G1<->G3 migration only |
| M18 | Bottleneck G5+G6, bidirectional G2<->G3 migration |
| M19 | Bottlenecks G5+G6+G1, bidirectional G2<->G3 migration |
| M20 | No bottleneck, G2<->G3 and G5<->G6 migration |
| M21 | Bottleneck G5+G6, G2<->G3 and G5<->G6 migration |
| M22 | Bottlenecks G5+G6+G1, no migration (tests the necessity of gene flow) |
| M23 | Bottlenecks, unidirectional G2->G3 migration only (tests directionality) |
| M24 | Bottlenecks, unidirectional G3->G2 migration only (tests directionality) |
| M25 | **Bottlenecks G5+G6+G1, G2<->G3 and G1<->G3 migration (best-supported model)** |
| M26 | Bottlenecks G5+G6+G1, G2<->G3 and G5<->G6 migration (tests G5<->G6 gene flow) |
| M27 | Bottlenecks G5+G6+G1, G2<->G3 and G1<->G5 migration (tests G1<->G5 gene flow) |

## Key result (M25)

- Best 100-run log-likelihood: **-7,164,399.384**; Delta AIC of the four runner-up models:
  M16 7,409 / M15 15,953 / M27 25,290 / M14 29,180.
- Root split of the five groups: ~581.8 kya (95% CI: 553.0-627.5 kya).
- Prolonged bottlenecks in Groups 5, 6 and 1 ending at ~67.7 kya (95% CI: 65.0-71.3 kya).
- Historical gene flow centred on Group 3 (G2<->G3 and G1<->G3).

## Requirements

- [fastsimcoal2](http://cmpg.unibe.ch/software/fastsimcoal2/) (`fsc28`, Linux 64-bit)
- Python 3 (step 04; standard library only)
- Sun Grid Engine (qsub) for the batch and bootstrap scripts

## Usage

```bash
# 1) Initial screen of all 27 models + 100 independent runs for the 5 best models
qsub run_new_models_linux.sh

# 2) M25 parametric bootstrap CI
bash 01_generate_bootstrap_100.sh /path/to/fsc28   # simulate 100 SFS
bash 02_prepare_reestimate_inputs.sh               # assemble replicate inputs
bash 03_run_reestimate_100.sh /path/to/fsc28 32 M25_boot 100000   # re-estimation
python3 04_summarize_ci.py M25_boot bootstrap_ci_summary.tsv      # CI summary
```

**Note on initial values**: in the final pipeline (`m25ci2.sh`), each bootstrap
re-estimate starts from values drawn at random within the prior ranges. Fixing the
initial values at the ML point estimates instead causes several low migration-rate
parameters to remain at their starting values (flat composite-likelihood surface),
which artificially collapses their bootstrap intervals.
