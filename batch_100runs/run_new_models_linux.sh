#!/bin/bash
#$ -N fsc2_100runs
#$ -l vf=100G,p=60
#$ -pe smp 60
#$ -j y
#$ -o fsc2_100run.out

# Configuration for Linux Server
WORK_DIR="/data_forrest/xubo/shiqinke/Maesa_pop/anly/fas/model_new"
OBS_DIR="/data_forrest/xubo/shiqinke/Maesa_pop/anly/easy_sfs/internal/easySFS_6pop_intergenic/easySFS_5pop_correct/fastsimcoal2"
FSC28="/home/xubo/Software/fsc28_linux64/fsc28"

# 5 best models run IN PARALLEL, 12 threads each (5 x 12 = 60 threads total)
NUM_RUNS=100
MODELS=("M25" "M16" "M15" "M27" "M14")
THREADS_PER_MODEL=12

# NOTE: per-run directories are KEPT as audit evidence (never deleted).
#       Previous results are renamed to *_old instead of removed.

cd "$WORK_DIR" || exit 1

echo "Starting fastsimcoal2 100-run batch: 5 models in parallel, 12 threads each (60 total)"
echo "Running on node: $(hostname)"

run_model () {
    local model=$1

    echo "========================================="
    echo "Running $model"
    echo "========================================="

    # Enter model directory
    cd "$WORK_DIR/$model" || exit 1

    # Archive previous results as *_old (keep as evidence, never delete)
    local old f
    for old in run[0-9]*; do
        [ -d "$old" ] || continue
        case "$old" in *_old*) continue ;; esac
        if [ -d "${old}_old" ]; then
            mv "$old" "${old}_old2"
        else
            mv "$old" "${old}_old"
        fi
    done
    for f in best_${model}.* "${model}_allruns.txt"; do
        [ -e "$f" ] && mv "$f" "${f}_old"
    done

    # Clean up any old observation files and link new ones
    rm -f *.obs

    # 5 populations mean Pop0 to Pop4.
    # Link 1D SFS (MAFpopX.obs)
    for i in {0..4}; do
        ln -sf $OBS_DIR/dataset2_LD_intergenic_exNRD_MAFpop${i}.obs ${model}_MAFpop${i}.obs
    done

    # Link 2D SFS (jointMAFpopX_Y.obs) where X > Y
    for i in {1..4}; do
        for j in $(seq 0 $((i-1))); do
            ln -sf $OBS_DIR/dataset2_LD_intergenic_exNRD_jointMAFpop${i}_${j}.obs ${model}_jointMAFpop${i}_${j}.obs
        done
    done

    # Per-run log (kept in model dir): one line per run, for empirical CI / audit
    local allruns="${model}_allruns.txt"

    # Run multiple independent runs to find global maximum likelihood
    local BEST_RUN=""
    local MAX_LHOOD=-999999999999.0

    for run in $(seq 1 $NUM_RUNS); do
        local run_dir="run${run}"
        mkdir -p $run_dir

        # Copy necessary files into run directory
        cp ${model}.tpl ${model}.est $run_dir/

        cd $run_dir || exit 1

        # Use soft links for SFS observation files to save space and inodes
        ln -sf ../*.obs .

        echo "  -> Starting $model run $run / $NUM_RUNS"
        # Run fastsimcoal2:
        # -c 12 : 12 threads per model (5 models in parallel = 60 threads total)
        # -B 12 : 12 batches for Brent maximization
        # -n 100000 : 100,000 simulations
        # -L 40 : 40 loops (ECM cycles)
        # -s 0 : new random seed for every run
        $FSC28 -t ${model}.tpl -e ${model}.est -m -0 -C 10 -n 100000 -L 40 -s 0 -M -B 12 -c $THREADS_PER_MODEL -q

        # Check likelihood
        local lhood_file="${model}/${model}.bestlhoods"
        if [ -f "$lhood_file" ]; then
            # Log this run (run id + full bestlhoods line) for later empirical CI
            awk -v r="run${run}" 'NR==2 {print r, $0}' "$lhood_file" >> "../$allruns"

            # Extract MaxEstLhood (second line, second to last column)
            local curr_lhood=$(awk 'NR==2 {print $(NF-1)}' "$lhood_file")

            # Compare using awk
            local is_greater=$(awk -v c="$curr_lhood" -v m="$MAX_LHOOD" 'BEGIN {print (c > m) ? 1 : 0}')
            if [ "$is_greater" -eq 1 ]; then
                MAX_LHOOD=$curr_lhood
                BEST_RUN="run${run}"
            fi
        fi

        cd ..
    done

    echo "Best run for $model was $BEST_RUN with MaxEstLhood = $MAX_LHOOD"

    # Extract the best run results to the main model directory
    if [ -n "$BEST_RUN" ]; then
        cp $BEST_RUN/${model}/${model}.bestlhoods best_${model}.bestlhoods
        cp $BEST_RUN/${model}/${model}.pv best_${model}.pv
        cp $BEST_RUN/${model}/${model}_maxL.par best_${model}_maxL.par

        echo "Results for $model:"
        cat best_${model}.bestlhoods
    else
        echo "Error: No bestlhoods file generated for $model"
    fi

    cd $WORK_DIR
}

# Launch all 5 models in parallel, each in its own subshell with its own log
for model in "${MODELS[@]}"; do
    run_model "$model" > "${model}_batch.log" 2>&1 &
done

# Wait for all models to finish
wait

echo ""
echo "All models completed!"
