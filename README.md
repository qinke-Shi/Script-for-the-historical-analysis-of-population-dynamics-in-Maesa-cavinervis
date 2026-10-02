# Script for the historical analysis of population dynamics in *Maesa cavinervis*

This repository archives the fastsimcoal2 model definitions and analysis pipeline used to
infer the historical demography (bottlenecks, divergence times and gene flow) of the five
genetic groups of *Maesa cavinervis*, as described in the manuscript.

The repository hosts the final 45 model definitions: 27 models (M01-M27) with sample sizes
18/50/22/34/28 and a three-stage (bootstrap-recovery) bottleneck structure, plus an expanded
family **M28-M45** that repeats the 18 bottleneck-model structures with **independently
estimated bottleneck onset (TBOT_Gx) and duration (DURBOT_Gx) for each group**. The `.tpl`
and `.est` files are the final configurations in use for the reported results (descriptive
comments only; prior ranges already incorporate the bounds identified during model tuning).

## Overview of the analysis

- **Data**: intergenic SNPs of the LD-pruned Dataset 2 (`dataset2_LD_intergenic_exNRD`),
  93 individuals assigned to five genetic groups; the two NRD samples are excluded and
  91 individuals enter the final SFS. Observed SFS produced with easySFS
  (projection 18/50/22/34/28 for demes G5/G6/G3/G2/G1).
- **fastsimcoal2 settings**: SNPs treated as 150 bp partially linked DNA loci,
  mutation rate mu = 2.43e-8 per generation (5-year generation time), 100,000 simulations
  and 40 EC-M cycles per estimation (`-n 100000 -L 40 -M`).
- **Model comparison strategy**: all 45 models were each run once for an initial screen
  (`run_screening_27.sh` for M01-M27; `build_indbot_models.py` + the same command for
  M28-M45). The five best models (M25, M27, M45, M43, M26) are then re-run with 100
  independent optimizations each (`100run_top5v2.sh`, serial by rank, 24 workers x 2
  threads) to obtain the best estimate per model and an empirical run-to-run distribution
  for uncertainty.
- **Screening result (single-run lnL, AIC = 2k - 2 lnL)**: M25 is the best-supported model,
  followed by M27 (Delta AIC 2,488), M45 (3,027), M43 (3,633) and M26 (4,064); the full
  45-model comparison is reported in manuscript Table S18. All models share the same
  observed-SFS fingerprint (MaxObsLhood = -10,459,502.188), confirming an identical data
  basis across the comparison.

## Bottleneck parameterization (three-stage, botr/recr)

Bottlenecked demes use a proportional three-stage history instead of an instantaneous size
switch, so that the bottleneck Ne and the pre-bottleneck ancestral Ne are explicitly
estimated and exported:

```
modern Ne (NX)
   │  at TBOT_E: size × botr_X            (botr_X ∈ logunif 1e-3 .. 0.9-0.99, model-specific)
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
├── models/                        # 45 fastsimcoal2 model definitions (M01-M45)
│   ├── M01.tpl ... M45.tpl        # template files: historical model (demography/gene flow)
│   └── M01.est ... M45.est        # parameter files: priors + complex-parameter chains
└── batch_100runs/
    ├── run_screening_27.sh        # SGE (qsub): initial screen, 27 models x 1 run
    │                              #   (27 x -c 2 = 54 threads)
    ├── build_indbot_models.py     # builds M28-M45 (independent bottleneck times) from
    │                              #   M02-M27 by line-level surgery, with validation
    ├── run_top5_100run_botr.sh    # SGE (qsub): 100 independent runs per top-5 model
    │                              #   (dynamic queue, 30 workers x -c 2 -B 2)
    ├── 100run_top5v2.sh           # screen/nohup variant: 100-run batch for the top five
    │                              #   (serial by rank, 24 workers x -c 2 -B 2 = 48 threads),
    │                              #   fsc28 v2.7.0.9 command syntax
    ├── run_100run_fac6new.sh      # earlier screen/nohup 100-run variant (resume-aware)
    ├── m25_boot_run.sh            # M25 parametric-bootstrap pipeline: builds the simulation
    │                              #   par from the ML point estimates (real 150-bp locus
    │                              #   structure), simulates 100 replicate SFS sets,
    │                              #   re-estimates all parameters from each replicate with
    │                              #   random initial values and summarizes the 2.5/50/97.5
    │                              #   percentiles per parameter
    └── plot_six_model_schematics_exNRD.py  # six-panel demographic schematics (top5 + M01)
```

## Model definitions (M01-M45)

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
| M28 | As M02, with independently estimated bottleneck onset and duration for each group |
| M29 | As M07, with independently estimated bottleneck onset and duration for each group |
| M30 | As M08, with independently estimated bottleneck onset and duration for each group |
| M31 | As M09, with independently estimated bottleneck onset and duration for each group |
| M32 | As M11, with independently estimated bottleneck onset and duration for each group |
| M33 | As M16, with independently estimated bottleneck onset and duration for each group |
| M34 | As M18, with independently estimated bottleneck onset and duration for each group |
| M35 | As M21, with independently estimated bottleneck onset and duration for each group |
| M36 | As M12, with independently estimated bottleneck onset and duration for each group |
| M37 | As M13, with independently estimated bottleneck onset and duration for each group |
| M38 | As M14, with independently estimated bottleneck onset and duration for each group |
| M39 | As M19, with independently estimated bottleneck onset and duration for each group |
| M40 | As M22, with independently estimated bottleneck onset and duration for each group |
| M41 | As M23, with independently estimated bottleneck onset and duration for each group |
| M42 | As M24, with independently estimated bottleneck onset and duration for each group |
| M43 | As M25, with independently estimated bottleneck onset and duration for each group |
| M44 | As M26, with independently estimated bottleneck onset and duration for each group |
| M45 | As M27, with independently estimated bottleneck onset and duration for each group |

k (number of searched parameters) ranges from 9 (M01) to 25 (M43/M44/M45).

## The 100-run protocol

For each of the five best models (M25, M27, M45, M43, M26), 100 independent `-M`
estimations are executed with distinct random seeds (`seed = model_number * 10000 +
run_number`). Prior ranges are checked against the screening maximum-likelihood values
before each batch and widened where the search hits a bound; the files hosted here are
the final widened versions. Per-run
working directories (`run1` ... `run100`) are retained as audit
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
their prior bounds (observed for G6N, G1N and DURBOT). Out-of-bound maximum-likelihood
values require widening the corresponding priors before the 100-run refinement; the
`.est` files hosted here are the final versions after these checks.

**Parametric-bootstrap CI**: the pipeline (`m25_boot_run.sh`) simulates 100 joint SFS
replicates under the M25 maximum-likelihood point estimates (using the real locus
structure: 677,266 intergenic SNPs as 150-bp DNA blocks, mu = 2.43e-8), then re-estimates
all parameters from every replicate with random initial values and wide priors, and
reports the 2.5/50/97.5 percentiles per parameter. Random initial values are essential:
fixed ML starting values collapse the CIs on flat composite-likelihood ridges.
