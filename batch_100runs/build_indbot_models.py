# Build the independent-bottleneck-time model family M28-M45 from the shared-episode
# family M01-M27 (line-level surgery, LF output).
# Transformation: the shared TBOT_E / DURBOT become per-deme TBOT_Gx / DURBOT_Gx; the
# three-stage botr/recr structure, PNBOT/PAN exports, and the TIME1 two-step chain
# (TENDBOT_G5 = TBOT_G5+DURBOT_G5 ; TIME1 = TENDBOT_G5+INC1_POSTBOT) are preserved.
# Mapping (source -> target), keeping the family order:
#   M02->M28 M07->M29 M08->M30 M09->M31 M11->M32 M16->M33 M18->M34 M21->M35
#   M12->M36 M13->M37 M14->M38 M19->M39 M22->M40 M23->M41 M24->M42 M25->M43
#   M26->M44 M27->M45
import re
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent

MAP = [
    ("M02", "M28"), ("M07", "M29"), ("M08", "M30"), ("M09", "M31"),
    ("M11", "M32"), ("M16", "M33"), ("M18", "M34"), ("M21", "M35"),
    ("M12", "M36"), ("M13", "M37"), ("M14", "M38"), ("M19", "M39"),
    ("M22", "M40"), ("M23", "M41"), ("M24", "M42"), ("M25", "M43"),
    ("M26", "M44"), ("M27", "M45"),
]
BOT_DEMES = {  # bottlenecked demes per source model (deme index -> group)
    "M02": ["5", "6"], "M07": ["5", "6"], "M08": ["5", "6"], "M09": ["5", "6"],
    "M11": ["5", "6"], "M16": ["5", "6"], "M18": ["5", "6"], "M21": ["5", "6"],
    "M12": ["5", "6", "1"], "M13": ["5", "6", "1"], "M14": ["5", "6", "1"],
    "M19": ["5", "6", "1"], "M22": ["5", "6", "1"], "M23": ["5", "6", "1"],
    "M24": ["5", "6", "1"], "M25": ["5", "6", "1"], "M26": ["5", "6", "1"],
    "M27": ["5", "6", "1"],
}

RE_TPL_TBOT = re.compile(r"^TBOT_E (\d+) (\d+) 0 botr_G(\d) 0 0$")
RE_EST_TBOT = re.compile(r"^1 TBOT_E unif (\S+) (\S+) output$")
RE_EST_DUR = re.compile(r"^1 DURBOT unif (\S+) (\S+) output$")
RE_EST_TEND = re.compile(r"^1 TENDBOT_G(\d) = TBOT_E\+DURBOT hide$")


def build_one(src: str, dst: str) -> dict:
    demes = BOT_DEMES[src]
    src_dir = BASE / src
    dst_dir = BASE / dst
    dst_dir.mkdir(exist_ok=True)

    # ---- tpl ----
    tpl_lines = (src_dir / f"{src}.tpl").read_text(encoding="utf-8").splitlines()
    out = []
    n_tbot = 0
    for i, line in enumerate(tpl_lines):
        m = RE_TPL_TBOT.match(line)
        if m:
            g = m.group(3)
            out.append(f"TBOT_G{g} {m.group(1)} {m.group(2)} 0 botr_G{g} 0 0")
            n_tbot += 1
        elif i == 0:
            out.append(line.replace(src, dst, 1))
        else:
            out.append(line)
    assert n_tbot == len(demes), f"{dst}: tpl TBOT lines {n_tbot} != {len(demes)}"
    (dst_dir / f"{dst}.tpl").write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")

    # ---- est (per-deme lines inherit the source model's prior bounds) ----
    est_lines = (src_dir / f"{src}.est").read_text(encoding="utf-8").splitlines()
    out = []
    header_done = False
    tbot_bounds = dur_bounds = None
    for line in est_lines:
        m = RE_EST_TBOT.match(line)
        if m:
            tbot_bounds = (m.group(1), m.group(2))
            for g in demes:
                out.append(f"1 TBOT_G{g} unif {tbot_bounds[0]} {tbot_bounds[1]} output")
            continue
        m = RE_EST_DUR.match(line)
        if m:
            dur_bounds = (m.group(1), m.group(2))
            for g in demes:
                out.append(f"1 DURBOT_G{g} unif {dur_bounds[0]} {dur_bounds[1]} output")
            continue
        m = RE_EST_TEND.match(line)
        if m:
            g = m.group(1)
            out.append(f"1 TENDBOT_G{g} = TBOT_G{g}+DURBOT_G{g} hide")
            continue
        if not header_done and line.startswith("//"):
            out.append(f"// {dst}: independent-bottleneck-time variant of {src} "
                       "(three-stage botr/recr; PNBOT=NX*botr, PAN=PNBOT*recr; "
                       "per-deme onset TBOT_Gx and duration DURBOT_Gx; TIME1 anchored to TENDBOT_G5)")
            header_done = True
            continue
        out.append(line)
    (dst_dir / f"{dst}.est").write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")

    # ---- validate (param lines only; comments may mention old names) ----
    est_text = "\n".join(l for l in (dst_dir / f"{dst}.est").read_text(encoding="utf-8").splitlines()
                         if not l.strip().startswith("//"))
    tpl_text = (dst_dir / f"{dst}.tpl").read_text(encoding="utf-8")
    errs = []
    if "TBOT_E" in est_text or "TBOT_E" in tpl_text:
        errs.append("TBOT_E residue")
    if re.search(r"\bDURBOT\b(?!_G)", est_text):
        errs.append("standalone DURBOT residue")
    for g in demes:
        if f"1 TBOT_G{g} unif {tbot_bounds[0]} {tbot_bounds[1]} output" not in est_text:
            errs.append(f"missing TBOT_G{g}")
        if f"1 DURBOT_G{g} unif {dur_bounds[0]} {dur_bounds[1]} output" not in est_text:
            errs.append(f"missing DURBOT_G{g}")
        if f"1 TENDBOT_G{g} = TBOT_G{g}+DURBOT_G{g} hide" not in est_text:
            errs.append(f"missing TENDBOT_G{g} chain")
        if f"1 PNBOT_G{g} =" not in est_text or f"1 PAN_G{g} =" not in est_text:
            errs.append(f"missing PNBOT/PAN_G{g}")
    if "1 TIME1 = TENDBOT_G5+INC1_POSTBOT output" not in est_text:
        errs.append("TIME1 two-step chain missing")
    for g in demes:
        if f"TBOT_G{g} " not in tpl_text:
            errs.append(f"tpl missing TBOT_G{g} event")
    ss = re.search(r"Samples sizes and samples age\n(\d+)\n(\d+)\n(\d+)\n(\d+)\n(\d+)\n", tpl_text)
    if not ss or (ss.group(1), ss.group(2), ss.group(3), ss.group(4), ss.group(5)) != ("18", "50", "22", "34", "28"):
        errs.append("sample sizes wrong")
    k = sum(1 for line in est_text.splitlines() if re.match(r"^[01]\s+\w+\s+(unif|logunif)\b", line))
    return {"src": src, "dst": dst, "k": k, "tbot": tbot_bounds, "dur": dur_bounds, "errs": errs}


def main():
    only = sys.argv[1] if len(sys.argv) > 1 else None
    results = []
    for src, dst in MAP:
        if only and dst != only:
            continue
        r = build_one(src, dst)
        results.append(r)
        flag = "OK " if not r["errs"] else "ERR"
        print(f"[{flag}] {r['src']}->{r['dst']} k={r['k']} tbot={r['tbot']} dur={r['dur']} {'; '.join(r['errs'])}")
    if any(r["errs"] for r in results):
        sys.exit(1)
    print(f"built {len(results)} models, all validations passed")


if __name__ == "__main__":
    main()
