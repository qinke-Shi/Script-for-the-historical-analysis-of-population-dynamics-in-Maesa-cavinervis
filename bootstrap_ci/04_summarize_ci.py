from __future__ import annotations

from pathlib import Path
import csv
import math
import sys

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('M25_boot')
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path('bootstrap_ci_summary.tsv')

best_files = sorted(ROOT.glob('M25_boot_[0-9]*/M25_boot/M25_boot.bestlhoods'))
if not best_files:
    raise SystemExit(f'No M25_boot.bestlhoods files found under {ROOT}')

rows = []
header = None
for best in best_files:
    lines = [line.strip() for line in best.read_text(encoding='utf-8').splitlines() if line.strip()]
    if len(lines) < 2:
        continue
    current_header = lines[0].split()
    values = lines[1].split()
    if header is None:
        header = current_header
    if current_header != header:
        raise SystemExit(f'Header mismatch in {best}')
    row = {name: float(value) for name, value in zip(header, values)}
    row['replicate'] = best.parent.name
    rows.append(row)

if not rows:
    raise SystemExit('No usable bootstrap rows were found')


def percentile(values: list[float], pct: float) -> float:
    vals = sorted(values)
    if len(vals) == 1:
        return vals[0]
    pos = (len(vals) - 1) * pct
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return vals[lo]
    frac = pos - lo
    return vals[lo] + (vals[hi] - vals[lo]) * frac


excluded = {'MaxEstLhood', 'MaxObsLhood'}
param_names = [name for name in header if name not in excluded]

with OUT.open('w', encoding='utf-8', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['parameter', 'n', 'mean', 'median', 'ci_2.5', 'ci_97.5'])
    for name in param_names:
        values = [row[name] for row in rows]
        mean = sum(values) / len(values)
        writer.writerow([
            name,
            len(values),
            f'{mean:.10g}',
            f'{percentile(values, 0.5):.10g}',
            f'{percentile(values, 0.025):.10g}',
            f'{percentile(values, 0.975):.10g}',
        ])

print(f'Wrote {OUT} from {len(rows)} bootstrap replicates')
