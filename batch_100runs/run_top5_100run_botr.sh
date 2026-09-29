#!/bin/bash
#$ -N fsc6_top5_botr
#$ -l vf=100G,p=60
#$ -pe smp 60
#$ -j y
#$ -cwd
#$ -o top5_100run_botr.out

# =====================================================================
# Top-5 models 100-run batch, botr/recr bottleneck structure (v2)
# Models: previous top5 (user decision 2026-09-28, screening still running):
#   M25 > M12 > M22 > M13 > M14
# tpl/est source: model6_new (new botr/recr configs, sample 18/50/22/34/28)
# Scheduling: dynamic task queue, 30 workers x (-c 2 -B 2) = 60 threads.
# M25 tasks are placed FIRST in the queue so M25 finishes first.
# Each run gets a distinct seed (-r): model# * 10000 + run#.
# Submit with:  qsub -hold_jid fsc6_scr27 run_top5_100run_botr.sh
#   (auto-starts when the screening job finishes -> no thread oversubscription)
# =====================================================================

SRC_DIR="/data_forrest/xubo/shiqinke/Maesa_pop/anly/fas/model6_new"
WORK_DIR="/data_forrest/xubo/shiqinke/Maesa_pop/anly/fas/model6_top5_100run_botr"
OBS_DIR="/data_forrest/xubo/shiqinke/Maesa_pop/anly/easy_sfs/internal/easySFS_6pop_intergenic/easySFS_5pop_correct/fastsimcoal2"
FSC28="/home/xubo/Software/fsc28_linux64/fsc28"

NUM_RUNS=100
MODELS=("M25" "M12" "M22" "M13" "M14")   # M25 first = highest priority
WORKERS=30                                # 30 runs x 2 threads = 60 threads

mkdir -p "$WORK_DIR" || exit 1
cd "$WORK_DIR" || exit 1
echo "Top5 100-run (botr/recr) started on $(hostname) at $(date)"

# ---- fetch new-structure tpl/est from model6_new ----------------------
for model in "${MODELS[@]}"; do
    mkdir -p "$WORK_DIR/$model"
    cp "$SRC_DIR/$model/${model}.tpl" "$SRC_DIR/$model/${model}.est" "$WORK_DIR/$model/" || {
        echo "MISSING new-structure input files for $model in $SRC_DIR/$model"; exit 1; }
done

# ---- obs symlinks at model-dir level (run dirs link to ../) -----------
# marginals: deme0..4 = groups 5,6,3,2,1 (proj 18,50,22,34,28) -- NOTE:
# 5pop_correct marginals are named {5,6,3,2,1}_MAFpop0.obs (NOT dataset2_*)
for model in "${MODELS[@]}"; do
    cd "$WORK_DIR/$model"
    rm -f *.obs
    ln -sf $OBS_DIR/5_MAFpop0.obs ${model}_MAFpop0.obs
    ln -sf $OBS_DIR/6_MAFpop0.obs ${model}_MAFpop1.obs
    ln -sf $OBS_DIR/3_MAFpop0.obs ${model}_MAFpop2.obs
    ln -sf $OBS_DIR/2_MAFpop0.obs ${model}_MAFpop3.obs
    ln -sf $OBS_DIR/1_MAFpop0.obs ${model}_MAFpop4.obs
    for i in {1..4}; do
        for j in $(seq 0 $((i-1))); do
            ln -sf $OBS_DIR/dataset2_LD_intergenic_exNRD_jointMAFpop${i}_${j}.obs ${model}_jointMAFpop${i}_${j}.obs
        done
    done
done
cd "$WORK_DIR"

# ---- per-run worker ----------------------------------------------------
cat > "$WORK_DIR/run_one.sh" <<'EOS'
#!/bin/bash
M=$1; I=$2
WORK_DIR="/data_forrest/xubo/shiqinke/Maesa_pop/anly/fas/model6_top5_100run_botr"
FSC28="/home/xubo/Software/fsc28_linux64/fsc28"
D="$WORK_DIR/$M/run$I"
mkdir -p "$D"
cd "$D" || exit 1
cp "$WORK_DIR/$M/${M}.tpl" "$WORK_DIR/$M/${M}.est" .
ln -sf ../*.obs .
SEED=$(( 10#${M#M} * 10000 + I ))
$FSC28 -t ${M}.tpl -e ${M}.est -m -0 -C 10 -n 100000 -L 40 -s 0 -M -c 2 -B 2 -r $SEED -q > run.log 2>&1
if [ -f "${M}/${M}.bestlhoods" ]; then
    awk -v r="run$I" 'NR==2 {print r, $0}' "${M}/${M}.bestlhoods" >> "$WORK_DIR/${M}_allruns.txt"
    echo "$(date '+%F %T') $M run$I OK" >> "$WORK_DIR/progress.txt"
else
    echo "$(date '+%F %T') $M run$I FAILED" >> "$WORK_DIR/progress.txt"
fi
EOS
chmod +x "$WORK_DIR/run_one.sh"

# ---- task list: M25's 100 runs first, then M12, M22, M13, M14 --------
: > "$WORK_DIR/tasks.txt"
for model in "${MODELS[@]}"; do
    for run in $(seq 1 $NUM_RUNS); do
        echo "$model $run" >> "$WORK_DIR/tasks.txt"
    done
done

# ---- launch the queue: 30 concurrent runs x 2 threads = 60 threads ---
xargs -P $WORKERS -n 2 "$WORK_DIR/run_one.sh" < "$WORK_DIR/tasks.txt"

# ---- summarize best runs ----------------------------------------------
echo "All 500 runs finished at $(date); extracting best runs..."
for model in "${MODELS[@]}"; do
    cd "$WORK_DIR/$model" || continue
    BEST=""; MAXL=-999999999999
    while read -r rid rest; do
        [ -z "$rid" ] && continue
        L=$(echo "$rest" | awk '{print $(NF-1)}')
        GT=$(awk -v c="$L" -v m="$MAXL" 'BEGIN {print (c > m) ? 1 : 0}')
        if [ "$GT" = "1" ]; then MAXL=$L; BEST=$rid; fi
    done < "${model}_allruns.txt"
    if [ -n "$BEST" ]; then
        cp "$BEST/${model}/${model}.bestlhoods" "best_${model}.bestlhoods"
        cp "$BEST/${model}/${model}.pv" "best_${model}.pv"
        cp "$BEST/${model}/${model}_maxL.par" "best_${model}_maxL.par"
        echo "Best run for $model: $BEST (MaxEstLhood = $MAXL)"
    else
        echo "ERROR: no valid runs for $model"
    fi
done

# ---- AIC summary (k = number of search params (unif/logunif) in est) --
cd "$WORK_DIR"
{
    echo -e "model\tk\tMaxEstLhood\tMaxObsLhood\tAIC\tdAIC"
    for model in "${MODELS[@]}"; do
        B="$model/best_${model}.bestlhoods"
        [ -f "$B" ] || continue
        L=$(awk 'NR==2 {print $(NF-1)}' "$B")
        O=$(awk 'NR==2 {print $NF}' "$B")
        k=$(grep -cE '^[01][[:space:]].*unif' "${model}.est")
        AIC=$(awk -v k="$k" -v l="$L" 'BEGIN {print 2*k - 2*l}')
        echo -e "$model\t$k\t$L\t$O\t$AIC"
    done | sort -t$'\t' -k5,5g | awk -F'\t' 'NR==1 {best=$5} {printf "%s\t%s\t%s\t%s\t%s\t%.1f\n", $1, $2, $3, $4, $5, $5-best}'
} > AIC_summary.txt

touch ALL_DONE
echo "ALL_DONE at $(date)"
cat AIC_summary.txt
