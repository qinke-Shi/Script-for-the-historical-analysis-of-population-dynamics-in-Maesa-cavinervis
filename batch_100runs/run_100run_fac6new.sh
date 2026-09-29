#!/bin/bash
# run_100run_fac6new.sh (GX4-151): 100-run batch in screen
# models SERIAL: M25 -> M27 -> M26 -> M24 -> M13 (M25 first = priority)
# within each model: 16 workers x (-c 2 -B 2) = 32 threads
# fsc28 v2.7.0.9: -E 1 required with -e; seed via -r (NO -s flag)
BASE=/data/hs012/shiqinke/snp1/fac_6new
FSC=/data/hs012/.conda/envs/fastsimocal2/bin/fsc28
NUM_RUNS=100
WORKERS=16
MODELS=("M25" "M27" "M26" "M24" "M13")

cd "$BASE" || exit 1
echo "100-run started $(date '+%F %T') on $(hostname): models ${MODELS[*]}, $WORKERS workers x 2 threads"

cat > "$BASE/run_one_100.sh" <<'EOS'
#!/bin/bash
M=$1
I=$2
BASE=/data/hs012/shiqinke/snp1/fac_6new
FSC=/data/hs012/.conda/envs/fastsimocal2/bin/fsc28
D="$BASE/$M/run$I"
mkdir -p "$D"
cd "$D" || exit 1
cp "$BASE/$M/${M}.tpl" "$BASE/$M/${M}.est" .
ln -sf ../*.obs .
SEED=$(( 10#${M#M} * 10000 + I ))
"$FSC" -t "${M}.tpl" -e "${M}.est" -M -E 1 -m -0 -C 10 -n 100000 -L 40 -c 2 -B 2 -r "$SEED" -x -q > run.log 2>&1
if [ -f "${M}/${M}.bestlhoods" ]; then
  awk -v r="run$I" 'NR==2 {print r, $0}' "${M}/${M}.bestlhoods" >> "$BASE/${M}_allruns.txt"
  echo "$(date '+%F %T') $M run$I OK" >> "$BASE/progress.txt"
else
  echo "$(date '+%F %T') $M run$I FAILED" >> "$BASE/progress.txt"
fi
EOS
chmod +x "$BASE/run_one_100.sh"

: > "$BASE/progress.txt"
for model in "${MODELS[@]}"; do
  echo "$(date '+%F %T') === $model: $NUM_RUNS runs start ===" >> "$BASE/progress.txt"
  : > "$BASE/${model}_allruns.txt"
  seq 1 "$NUM_RUNS" | xargs -P "$WORKERS" -n 1 -I{} "$BASE/run_one_100.sh" "$model" {}
  echo "$(date '+%F %T') === $model: $NUM_RUNS runs done ===" >> "$BASE/progress.txt"
done

echo "All 500 runs finished at $(date '+%F %T'); extracting best runs..."
cd "$BASE" || exit 1
for model in "${MODELS[@]}"; do
  cd "$BASE/$model" || continue
  BEST=""
  MAXL=-999999999999
  while read -r rid rest; do
    [ -z "$rid" ] && continue
    L=$(echo "$rest" | awk '{print $(NF-1)}')
    GT=$(awk -v c="$L" -v m="$MAXL" 'BEGIN {print (c > m) ? 1 : 0}')
    if [ "$GT" = "1" ]; then
      MAXL=$L
      BEST=$rid
    fi
  done < "$BASE/${model}_allruns.txt"
  if [ -n "$BEST" ]; then
    cp "$BEST/${model}/${model}.bestlhoods" "best_${model}.bestlhoods"
    cp "$BEST/${model}/${model}.pv" "best_${model}.pv"
    cp "$BEST/${model}/${model}_maxL.par" "best_${model}_maxL.par"
    echo "Best run for $model: $BEST (MaxEstLhood = $MAXL)" | tee -a "$BASE/progress.txt"
  else
    echo "ERROR: no valid runs for $model" | tee -a "$BASE/progress.txt"
  fi
done

cd "$BASE" || exit 1
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
echo "ALL_DONE at $(date '+%F %T')"
cat AIC_summary.txt
