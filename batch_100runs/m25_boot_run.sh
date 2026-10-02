#!/bin/bash
# m25_boot_run.sh: M25 parametric bootstrap anchored to the pre-extracted best_M25.* files.
# Phases: build sim par -> 10 shards x 10 reps -> assemble 100 reps -> 16 workers re-estimate -> CI summary
BASE=/data/hs012/shiqinke/snp1/fac_6new
FSC=/data/hs012/.conda/envs/fastsimocal2/bin/fsc28
WORKERS=16
BOOT=$BASE/M25_boot

cd "$BASE" || exit 1

# ---------- build simulation par (ML values + real locus structure) ----------
mkdir -p "$BOOT"
cd "$BOOT" || exit 1
awk 'NR==1 { print "// M25 parametric bootstrap par: 677266 intergenic SNPs as 150-bp DNA blocks (mu 2.43e-8 per gen)"; next }
     /^1 0$/ && !d1 { print "677266 0"; d1=1; next }
     /^FREQ 1 0 2.43e-8 OUTEXP$/ && !d2 { print "DNA 150 0 2.43e-8 OUTEXP"; d2=1; next }
     { print }' "$BASE/M25/best_M25_maxL.par" > M25_boot.par
grep -q '^677266 0$' M25_boot.par || { echo "$(date '+%F %T') BOOTSTRAP ABORT: par loci line build failed" >> "$BASE/progress.txt"; exit 1; }
grep -q '^DNA 150 0 2.43e-8 OUTEXP$' M25_boot.par || { echo "$(date '+%F %T') BOOTSTRAP ABORT: par datatype line build failed" >> "$BASE/progress.txt"; exit 1; }

# ---------- simulate 100 replicates (10 parallel shards x 10 reps) ----------
rm -rf "$BOOT"/simshard*
for s in 1 2 3 4 5 6 7 8 9 10; do
  mkdir -p "$BOOT/simshard$s"
  cp "$BOOT/M25_boot.par" "$BOOT/simshard$s/"
  (
    cd "$BOOT/simshard$s" || exit 1
    "$FSC" -i M25_boot.par -n 10 -j -m -s0 -x -I -q -r $((25010 + s)) > sim.log 2>&1
    echo "$(date '+%F %T') sim shard $s rc=$?" >> "$BOOT/sim_shards.log"
  ) &
done
wait
nrep=$(find "$BOOT"/simshard* -type d -name 'M25_boot_[0-9]*' 2>/dev/null | wc -l)
echo "$(date '+%F %T') simulation done: $nrep/100 replicate SFS sets" >> "$BASE/progress.txt"

# ---------- assemble re-estimation directories ----------
cat > "$BASE/run_one_boot.sh" <<'EOI'
#!/bin/bash
d=$1
BASE=/data/hs012/shiqinke/snp1/fac_6new
FSC=/data/hs012/.conda/envs/fastsimocal2/bin/fsc28
cd "$d" || exit 1
n=$(basename "$d")
i=${n#rep}
"$FSC" -t M25_boot.tpl -e M25_boot.est -M -E 1 -m -0 -C 10 -n 100000 -L 40 -c 2 -B 2 -r $((26000 + i)) -x -q > run.log 2>&1
if [ -f M25_boot/M25_boot.bestlhoods ]; then
  awk -v r="$n" 'NR==2 {print r, $0}' M25_boot/M25_boot.bestlhoods >> "$BASE/M25_boot_allruns.txt"
  echo "$(date '+%F %T') BOOT $n OK" >> "$BASE/progress.txt"
else
  echo "$(date '+%F %T') BOOT $n FAILED" >> "$BASE/progress.txt"
fi
EOI
chmod +x "$BASE/run_one_boot.sh"

mkdir -p "$BOOT/reps"
r=0
for s in 1 2 3 4 5 6 7 8 9 10; do
  for repdir in $(find "$BOOT/simshard$s" -type d -name 'M25_boot_[0-9]*' 2>/dev/null); do
    r=$((r + 1))
    d="$BOOT/reps/rep$r"
    mkdir -p "$d"
    cp "$BASE/M25/M25.tpl" "$d/M25_boot.tpl"
    cp "$BASE/M25/M25.est" "$d/M25_boot.est"
    nj=0
    for f in "$repdir"/M25_boot_jointMAFpop*.obs; do
      [ -e "$f" ] || continue
      ln -sf "$f" "$d/$(basename "$f")"
      nj=$((nj + 1))
    done
    if [ "$nj" -lt 10 ]; then
      for f in "$repdir"/M25_boot_jointMAFpop*.txt; do
        [ -e "$f" ] || continue
        b=$(basename "$f" .txt)
        cp "$f" "$d/$b.obs"
        nj=$((nj + 1))
      done
    fi
    for f in "$repdir"/M25_boot_MAFpop*.obs; do
      [ -e "$f" ] && ln -sf "$f" "$d/$(basename "$f")"
    done
    if [ "$nj" -lt 10 ]; then
      echo "$(date '+%F %T') WARNING: $d has only $nj joint obs (incomplete replicate)" >> "$BASE/progress.txt"
      rm -rf "$d"
      r=$((r - 1))
    fi
  done
done
echo "$(date '+%F %T') assembled $r/100 re-estimation directories" >> "$BASE/progress.txt"
if [ "$r" -lt 95 ]; then
  echo "$(date '+%F %T') BOOTSTRAP ABORT: too few complete replicates assembled ($r/100)" >> "$BASE/progress.txt"
  exit 1
fi

# ---------- re-estimation (16 workers x 2 threads, random initial values) ----------
: > "$BASE/M25_boot_allruns.txt"
echo "$(date '+%F %T') bootstrap re-estimation start ($r replicates, $WORKERS workers x 2 threads)" >> "$BASE/progress.txt"
printf '%s\n' "$BOOT"/reps/rep* | xargs -P "$WORKERS" -n 1 "$BASE/run_one_boot.sh"

# ---------- CI summary ----------
cat > "$BOOT/summarize_ci.py" <<'EOP'
import glob, math
from pathlib import Path
BASE = Path('/data/hs012/shiqinke/snp1/fac_6new')
best_files = sorted(glob.glob(str(BASE / 'M25_boot/reps/rep*/M25_boot/M25_boot.bestlhoods')))
rows = []
header = None
for f in best_files:
    lines = [l.strip() for l in open(f, encoding='utf-8') if l.strip()]
    if len(lines) < 2:
        continue
    h = lines[0].split()
    v = lines[1].split()
    if header is None:
        header = h
    row = dict(zip(h, v))
    row['rep'] = Path(f).parts[-3]
    rows.append(row)
def pct(vals, p):
    vals = sorted(vals)
    pos = (len(vals) - 1) * p
    lo = int(math.floor(pos)); hi = int(math.ceil(pos))
    return vals[lo] if lo == hi else vals[lo] + (vals[hi] - vals[lo]) * (pos - lo)
bl = (BASE / 'M25/best_M25.bestlhoods').read_text(encoding='utf-8').splitlines()
h = bl[0].split(); v = bl[1].split()
point = dict(zip(h, v))
with open(BASE / 'M25_boot_ci_summary.tsv', 'w', encoding='utf-8') as w:
    w.write('parameter\tn\tpoint\tpoint_in_CI\tmean\tmedian\tci2.5\tci97.5\n')
    for name in header:
        if name in ('MaxEstLhood', 'MaxObsLhood'):
            continue
        vals = [float(r[name]) for r in rows]
        lo = pct(vals, 0.025); hi = pct(vals, 0.975)
        pt = float(point[name])
        inside = 'yes' if lo <= pt <= hi else 'no'
        w.write('%s\t%d\t%.6g\t%s\t%.6g\t%.6g\t%.6g\t%.6g\n' % (name, len(vals), pt, inside, sum(vals) / len(vals), pct(vals, 0.5), lo, hi))
print('replicates used:', len(rows))
EOP
python3 "$BOOT/summarize_ci.py" > "$BOOT/summarize.log" 2>&1

touch "$BASE/BOOTSTRAP_DONE"
echo "$(date '+%F %T') BOOTSTRAP_DONE; summary:" >> "$BASE/progress.txt"
cat "$BASE/M25_boot_ci_summary.tsv" >> "$BASE/progress.txt" 2>/dev/null
echo "=== M25 BOOTSTRAP PIPELINE COMPLETE ==="
