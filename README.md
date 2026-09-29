# Script for the historical analysis of population dynamics in *Maesa cavinervis*

This repository archives the analysis scripts of the manuscript: (i) the fastsimcoal2
model definitions and analysis pipeline used to infer the historical demography
(bottlenecks, divergence times and gene flow) of the five genetic groups of
*Maesa cavinervis*, and (ii) the comparative genomics pipeline (orthology inference,
Yang & Smith (2014) tree-based paralog pruning, and MCMCTree molecular dating).

**Version note (2026-09-29).** The model definitions previously hosted here were based on an
outdated dataset and an incorrect bottleneck parameterization. This revision replaces all 27
models with the corrected definitions: sample sizes 18/50/22/34/28 and a three-stage
(bootstrap-recovery) bottleneck structure. Every file below is byte-identical to the
configuration that produced the reported screening results and is currently in use for the
100-run refinement.

## Overview of the analysis

- **Data**: intergenic SNPs of the LD-pruned Dataset 2 (`dataset2_LD_intergenic_exNRD`),
  93 individuals assigned to five genetic groups; the two NRD samples are excluded and
  91 individuals enter the final SFS. Observed SFS produced with easySFS
  (projection 18/50/22/34/28 for demes G5/G6/G3/G2/G1).
- **fastsimcoal2 settings**: SNPs treated as 150 bp partially linked DNA loci,
  mutation rate mu = 2.43e-8 per generation (5-year generation time), 100,000 simulations
  and 40 EC-M cycles per estimation (`-n 100000 -L 40 -M`).
- **Model comparison strategy**: all 27 models were each run once for an initial screen
  (`run_screening_27.sh`); the five best models are then re-run with 100 independent
  optimizations each (`run_top5_100run_botr.sh` / `run_100run_fac6new.sh`) to obtain the
  best estimate per model and an empirical run-to-run distribution for uncertainty.
- **Screening result (single-run lnL, AIC = 2k - 2 lnL)**: M25 is the best-supported model,
  followed by M27 (Delta AIC 2,488), M26 (4,064), M24 (5,054) and M13 (6,599); the full
  27-model comparison is reported in manuscript Table S18. All models share the same
  observed-SFS fingerprint (MaxObsLhood = -10,459,502.188), confirming an identical data
  basis across the comparison.

## Bottleneck parameterization (three-stage, botr/recr)

Bottlenecked demes use a proportional three-stage history instead of an instantaneous size
switch, so that the bottleneck Ne and the pre-bottleneck ancestral Ne are explicitly
estimated and exported:

```
modern Ne (NX)
   │  at TBOT_E: size × botr_X            (botr_X ∈ logunif 1e-3 .. 9e-1)
   ▼
bottleneck Ne  PNBOT_X = NX × botr_X      ← exported to .pv/.bestlhoods
   │  duration DURBOT (unif 10 .. 3000)
   ▼
at TBOT_E + DURBOT: size × recr_X         (recr_X ∈ logunif 1.01 .. 50)
   ▼
pre-bottleneck ancestral Ne  PAN_X = PNBOT_X × recr_X   ← exported
```

**Important fastsimcoal2 limitation.** fsc28 parses complex-parameter expressions with
**two operands only**; a three-term expression such as `TIME1 = TBOT_E + DURBOT +
INC1_POSTBOT` is *silently truncated* to `TBOT_E + DURBOT` (the third term is dropped
without any warning). All `.est` files here therefore use the two-step chain

```
TENDBOT_Gx = TBOT_E + DURBOT        (hidden intermediate parameter)
TIME1      = TENDBOT_G5 + INC1_POSTBOT
```

which was verified arithmetically on two independent fsc28 builds (event times in the
run-time `.par` files satisfy TIME1 - TENDBOT = INC1_POSTBOT exactly). Non-bottleneck
models use the equivalent binary form `TIME1 = INC1 + 0`.

## Repository layout

```
.
├── models/                        # 27 fastsimcoal2 model definitions (M01-M27)
│   ├── M01.tpl ... M27.tpl        # template files: historical model (demography/gene flow)
│   └── M01.est ... M27.est        # parameter files: priors + complex-parameter chains
├── batch_100runs/
│   ├── run_screening_27.sh        # SGE (qsub): initial screen, 27 models x 1 run
│   │                              #   (27 x -c 2 = 54 threads)
│   ├── run_top5_100run_botr.sh    # SGE (qsub): 100 independent runs per top-5 model
│   │                              #   (dynamic queue, 30 workers x -c 2 -B 2)
│   └── run_100run_fac6new.sh      # screen/nohup variant of the 100-run batch for
│                                  #   hosts without SGE (16 workers x -c 2 -B 2),
│                                  #   fsc28 v2.7.0.9 command syntax
└── comparative_genomics/          # orthology inference (OrthoFinder v3.0.1b1),
    │                              #   Yang & Smith (2014) paralog pruning, MCMCTree dating
    ├── 00_orthology_inference/    #   header standardisation, longest-isoform extraction,
    │                              #   OrthoFinder run, 705-OG list
    ├── 01_paralog_pruning_YangSmith2014/  # sensitivity pipeline wrappers + original
    │                              #   Y&S scripts + 2,612-OG list
    └── 02_molecular_dating/       #   main + sensitivity MCMCTree control files,
                                   #   calibrated guide tree
```

## Model definitions (M01-M27)

Demes follow the five genetic groups; the deme order in the template files is
G5, G6, G3, G2, G1 = Groups 5/6/3/2/1 with sample sizes 18/50/22/34/28.
Migration rates are forward rates per generation. Bottlenecked models share one bottleneck
episode: start `TBOT_E` (unif 1000-20000), duration `DURBOT`, with per-deme severity
`botr_*`/recovery `recr_*`; post-bottleneck split intervals are `INC1_POSTBOT`, `INC2`,
`INC3`, `INC4` (TIME1-TIME4 chains).

| Model | Demographic scenario |
|---|---|
| M01 | Baseline: no migration, no bottleneck (reference for model comparison) |
| M02 | Bottleneck in G5+G6, no migration |
| M03 | No bottleneck, bidirectional migration G2<->G3 |
| M04 | No bottleneck, G1<->G5 migration (rate collapsed to the lower prior bound; not supported) |
| M05 | No bottleneck, G5<->G6 migration |
| M06 | No bottleneck, G1<->G6 migration |
| M07 | Bottleneck G5+G6, G1<->G5 migration (one rate collapsed to the lower prior bound) |
| M08 | Bottleneck G5+G6, G5<->G6 migration |
| M09 | Bottleneck G5+G6, G1<->G6 migration |
| M10 | No bottleneck, G5<->G3 migration |
| M11 | Bottleneck G5+G6, G5<->G3 migration |
| M12 | Bottlenecks G5+G6+G1, G1<->G5 migration |
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

k (number of searched parameters) ranges from 9 (M01) to 21 (three-bottleneck models with
two migration pairs: M25/M26/M27).

## The 100-run protocol

For each of the five best models (M25, M27, M26, M24, M13), 100 independent `-M`
estimations are executed with distinct random seeds (`seed = model_number * 10000 +
run_number`). Per-run working directories (`run1` ... `run100`) are retained as audit
evidence; each run's final likelihood line is appended to `<model>_allruns.txt`. After all
runs complete, the run with the highest MaxEstLhood provides `best_<model>.bestlhoods`,
`best_<model>.pv` and `best_<model>_maxL.par`, and an AIC summary across the five models
is written. The observed-SFS fingerprint (MaxObsLhood) is checked for every model to
guarantee an identical data basis.

## Requirements

- [fastsimcoal2](http://cmpg.unibe.ch/software/fastsimcoal2/) (Linux 64-bit). Two builds
  were used and produce identical results for these files:
  - the classic `fsc28` build (command line `... -s 0 -M ...` accepted, seed via `-s`);
  - fsc28 **v2.7.0.9**, which **requires `-E 1` whenever `-e` is given** and takes the
    random seed via `-r` (there, `-s` means "output DNA as SNPs" instead):
    `fsc28 -t M.tpl -e M.est -M -E 1 -m -0 -C 10 -n 100000 -L 40 -c 2 -B 2 -r SEED -x -q`
- Python 3 (standard library only) for summary scripts
- Sun Grid Engine (qsub), or GNU screen/nohup on hosts without SGE

## Usage

```bash
# 1) Initial screen of all 27 models (SGE host)
qsub run_screening_27.sh

# 2) 100 independent runs per top-5 model
qsub run_top5_100run_botr.sh          # SGE host (planthunter-style queue)
# or, on a host without SGE (GNU screen):
screen -dmS fsc_100run bash run_100run_fac6new.sh
```

**Note on observed SFS**: fsc28 looks for the `_MAFpop*.obs` / `_jointMAFpop*.obs` files
in the directory where it is invoked; the batch scripts therefore link the 5 marginal and
10 joint observed SFS into every model (or run) directory before execution.

**Note on parameter ranges**: estimation by the Brent algorithm may push parameters beyond
their prior bounds (observed for G6N and DURBOT). Out-of-bound maximum-likelihood values
require widening the corresponding priors before the 100-run refinement
(M25.est v4: G6N upper 400k -> 600k, DURBOT upper 2000 -> 3000).

**Parametric-bootstrap CI**: the bootstrap pipeline of the previous repository version
belonged to the outdated model structure and has been removed. It will be re-derived from
the M25 maximum-likelihood estimates once the 100-run refinement completes.
