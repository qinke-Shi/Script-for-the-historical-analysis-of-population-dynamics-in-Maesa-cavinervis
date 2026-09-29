#!/bin/bash
#$ -N fsc6_scr27
#$ -l vf=100G,p=60
#$ -pe smp 60
#$ -cwd
#$ -j y
#$ -o fsc6_scr27.out

# 2026-09-27 initial screening rerun: 27 models x 1 run each, corrected sample sizes 18/50/22/34/28
# Parallel scheme: 27 models simultaneously, 2 threads each (-c 2) = 54 threads on a 60-slot node
WORK_DIR="/data_forrest/xubo/shiqinke/Maesa_pop/anly/fas/model6_new"
OBS_DIR="/data_forrest/xubo/shiqinke/Maesa_pop/anly/easy_sfs/internal/easySFS_6pop_intergenic/easySFS_5pop_correct/fastsimcoal2"
FSC28="/home/xubo/Software/fsc28_linux64/fsc28"
MODELS=(M01 M02 M03 M04 M05 M06 M07 M08 M09 M10 M11 M12 M13 M14 M15 M16 M17 M18 M19 M20 M21 M22 M23 M24 M25 M26 M27)

cd "$WORK_DIR" || exit 1
echo "Node $(hostname), screening start: $(date)"

link_obs () {
    local model=$1
    # marginals: deme0..4 = groups 5,6,3,2,1 (proj 18,50,22,34,28)
    ln -sf "$OBS_DIR/5_MAFpop0.obs" "${model}_MAFpop0.obs"
    ln -sf "$OBS_DIR/6_MAFpop0.obs" "${model}_MAFpop1.obs"
    ln -sf "$OBS_DIR/3_MAFpop0.obs" "${model}_MAFpop2.obs"
    ln -sf "$OBS_DIR/2_MAFpop0.obs" "${model}_MAFpop3.obs"
    ln -sf "$OBS_DIR/1_MAFpop0.obs" "${model}_MAFpop4.obs"
    # joint 2D SFS
    local pair
    for pair in 1_0 2_0 2_1 3_0 3_1 3_2 4_0 4_1 4_2 4_3; do
        ln -sf "$OBS_DIR/dataset2_LD_intergenic_exNRD_jointMAFpop${pair}.obs" "${model}_jointMAFpop${pair}.obs"
    done
}

run_model () {
    local model=$1
    cd "$WORK_DIR/$model" || exit 1
    link_obs "$model"
    echo "[$model] start $(date)"
    "$FSC28" -t ${model}.tpl -e ${model}.est -m -0 -C 10 -n 100000 -L 40 -s 0 -M -B 12 -c 2 -q
    local bf="${model}/${model}.bestlhoods"
    if [ -f "$bf" ]; then
        echo "$model $(sed -n 2p "$bf")" >> "$WORK_DIR/screen_progress.txt"
        echo "[$model] done $(date): $(sed -n 2p "$bf")"
    else
        echo "[$model] FAILED: no bestlhoods $(date)"
    fi
}

rm -f "$WORK_DIR/screen_progress.txt"
for m in "${MODELS[@]}"; do
    run_model "$m" > "$WORK_DIR/${m}_screen.log" 2>&1 &
done
wait
touch "$WORK_DIR/SCREEN_DONE"
echo "All 27 screening runs finished: $(date)"

# ---- AIC summary (k = number of search params (unif/logunif) in est) ----
RAW="$WORK_DIR/.aic_raw.txt"
: > "$RAW"
for m in "${MODELS[@]}"; do
    f="$WORK_DIR/$m/$m/$m.bestlhoods"
    k=$(grep -cE '^[01][[:space:]].*unif' "$WORK_DIR/$m/$m.est" || true)
    if [ -f "$f" ]; then
        obs=$(awk 'NR==1{for(i=1;i<=NF;i++)if($i=="MaxObsLhood")o=i}NR==2{print $o}' "$f")
        estl=$(awk 'NR==1{for(i=1;i<=NF;i++)if($i=="MaxEstLhood")e=i}NR==2{print $e}' "$f")
        aic=$(awk -v k="$k" -v l="$estl" 'BEGIN{printf "%.1f", 2*k-2*l}')
        printf '%s\t%s\t%s\t%s\t%s\n' "$m" "$k" "$obs" "$estl" "$aic" >> "$RAW"
    else
        printf '%s\t%s\tFAILED\tFAILED\t0\n' "$m" "$k" >> "$RAW"
    fi
done
{
    printf 'model\tk\tMaxObsLhood\tMaxEstLhood\tAIC\tdAIC\n'
    sort -t$'\t' -k5,5g "$RAW" | awk -F'\t' -v OFS='\t' 'NR==1{min=$5+0}{printf "%s\t%s\t%s\t%s\t%s\t%.1f\n",$1,$2,$3,$4,$5,$5-min}'
} > "$WORK_DIR/AIC_summary.txt"
rm -f "$RAW"
echo "=== AIC summary (k = number of search params in est) ==="
cat "$WORK_DIR/AIC_summary.txt"
