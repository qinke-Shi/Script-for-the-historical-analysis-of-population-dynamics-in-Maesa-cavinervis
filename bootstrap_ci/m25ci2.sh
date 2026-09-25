#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N m25ci2
#$ -l vf=32G,p=32
#$ -pe smp 32
FSC=/home/xubo/Software/fsc28_linux64/fsc28
cd /data_forrest/xubo/shiqinke/Maesa_pop/anly/fas/model_new/M25/bootstrap_ci100_randinit
bash 03_run_reestimate_100.sh "$FSC" 32 M25_boot 100000 > reestimate.log 2>&1
python3 04_summarize_ci.py M25_boot bootstrap_ci_summary.tsv > summarize.log 2>&1
echo DONE_M25CI2
