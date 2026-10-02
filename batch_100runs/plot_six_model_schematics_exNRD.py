# Schematics of the new top-5 models (M25, M27, M45, M43, M26) + isolation baseline M01,
# styled after E:\Maesa\snp\fsc_k6\model\plot_six_model_schematics.py.
# Results source: 45-model screening (M01-M27 shared bottleneck episode; M28-M45 per-deme
# independent bottleneck onset/duration), 2026-10-01. M25 uses its 100-run best (v4).
# All six share the standard tree: (G6,G5)->G3, (G2,G1), root at TIME4.
# Column order follows the species tree (sequences.phy.treefile): G6, G5, G3, G2, G1.
# Black branches = split order; gray dashed lines = event times (splits + bottleneck onsets);
# red arrows = forward gene flow (each arrow drawn up to the event where its flow window ends).
# Time axis is LINEAR per panel (no sqrt compression).

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "Times New Roman"

GENERATION_TIME = 5

# deme column order per species tree: G6, G5 | G3 | G2, G1
POP_X = {"G6": 0.8, "G5": 1.8, "G3": 2.9, "G2": 4.6, "G1": 5.6}
AXIS_X = 6.35

BASE = Path(__file__).resolve().parent
MODEL_DIR = BASE.parent

MODEL_FILES = {
    "M25": MODEL_DIR / "M25" / "plot" / "best_M25.bestlhoods",
    "M27": MODEL_DIR / "M27" / "M27" / "M27.bestlhoods",
    "M45": MODEL_DIR / "M45" / "M45" / "M45.bestlhoods",
    "M43": MODEL_DIR / "M43" / "M43" / "M43.bestlhoods",
    "M26": MODEL_DIR / "M26" / "M26" / "M26.bestlhoods",
    "M01": MODEL_DIR / "M01" / "M01" / "M01.bestlhoods",
}

# forward-time gene flow (source deme -> target deme) with the event that ends each flow
# window (per the tpl migration-matrix switches / deme mergers)
MODEL_FLOWS = {
    "M25": [("G2", "G3", "TIME3"), ("G3", "G2", "TIME3"),          # G2<->G3
            ("G1", "G3", "TIME3"), ("G3", "G1", "TIME3")],         # G1<->G3
    "M27": [("G2", "G3", "TIME3"), ("G3", "G2", "TIME3"),          # G2<->G3
            ("G1", "G5", "TIME1"), ("G5", "G1", "TIME1")],         # G1<->G5 (ends at G5 merger)
    "M45": [("G2", "G3", "TIME3"), ("G3", "G2", "TIME3"),          # G2<->G3
            ("G1", "G5", "TIME1"), ("G5", "G1", "TIME1")],         # G1<->G5
    "M43": [("G2", "G3", "TIME3"), ("G3", "G2", "TIME3"),          # G2<->G3
            ("G1", "G3", "TIME3"), ("G3", "G1", "TIME3")],         # G1<->G3
    "M26": [("G2", "G3", "TIME3"), ("G3", "G2", "TIME3"),          # G2<->G3
            ("G5", "G6", "TIME1"), ("G6", "G5", "TIME1")],         # G5<->G6 (ends at G6 merger)
    "M01": [],                                                     # isolation baseline: no bottleneck, no gene flow
}


def read_bestlhoods(path: Path) -> dict[str, float]:
    with path.open("r", encoding="utf-8") as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    header, values = rows[0], rows[1]
    return {k.strip(): float(v) for k, v in zip(header, values)}


def model_events(params: dict[str, float]) -> list[tuple[str, float]]:
    events = [(name, params[name]) for name in ("TIME1", "TIME2", "TIME3", "TIME4")]
    for name, value in params.items():
        if name.startswith("TBOT"):   # TBOT_E (shared) or TBOT_Gx (independent onset per deme)
            events.append((name, value))
    return sorted(events, key=lambda item: item[1])


def scaled_y(value: float, max_time: float) -> float:
    return max(value, 0.0) / max_time   # linear time axis (no compression)


def adjusted_label_positions(events, max_time: float) -> list[float]:
    raw = [scaled_y(t, max_time) for _, t in events]
    min_sep = 0.055
    placed = [raw[0]]
    for current in raw[1:]:
        placed.append(max(current, placed[-1] + min_sep))
    overflow = placed[-1] - 1.0
    if overflow > 0:
        placed = [value - overflow for value in placed]
    placed[0] = max(placed[0], 0.04)
    for idx in range(1, len(placed)):
        placed[idx] = max(placed[idx], placed[idx - 1] + min_sep)
    if placed[-1] > 1.02:
        compression = (1.02 - placed[0]) / max(placed[-1] - placed[0], 1e-9)
        placed = [placed[0] + (value - placed[0]) * compression for value in placed]
    return placed


def tree_nodes(params: dict[str, float]) -> dict[str, tuple[float, float]]:
    max_time = params["TIME4"]
    y = {name: scaled_y(params[name], max_time) for name in ("TIME1", "TIME2", "TIME3", "TIME4")}
    x_a56 = (POP_X["G6"] + POP_X["G5"]) / 2
    x_a563 = (x_a56 + POP_X["G3"]) / 2
    x_a21 = (POP_X["G2"] + POP_X["G1"]) / 2
    x_root = (x_a563 + x_a21) / 2
    return {
        **{name: (POP_X[name], 0.0) for name in POP_X},
        "A56": (x_a56, y["TIME1"]),
        "A563": (x_a563, y["TIME2"]),
        "A21": (x_a21, y["TIME3"]),
        "ROOT": (x_root, y["TIME4"]),
    }


def draw_tree(ax, nodes) -> None:
    segments = [
        ("G6", "A56"), ("G5", "A56"), ("A56", "A563"), ("G3", "A563"),
        ("G2", "A21"), ("G1", "A21"), ("A21", "ROOT"), ("A563", "ROOT"),
    ]
    for start, end in segments:
        (x0, y0), (x1, y1) = nodes[start], nodes[end]
        ax.plot([x0, x1], [y0, y1], color="black", lw=2.2, solid_capstyle="round", zorder=3)
    root_x, root_y = nodes["ROOT"]
    ax.plot([root_x, root_x + 0.34], [root_y, min(root_y + 0.08, 1.08)],
            color="black", lw=2.2, solid_capstyle="round", zorder=3)


def draw_event_guides(ax, params: dict[str, float]):
    events = model_events(params)
    max_time = params["TIME4"]
    label_positions = adjusted_label_positions(events, max_time)
    for (label, time_value), label_y in zip(events, label_positions):
        event_y = scaled_y(time_value, max_time)
        ax.hlines(event_y, xmin=0.35, xmax=AXIS_X - 0.08, colors="#9A9A9A", lw=0.9,
                  linestyles="dashed", zorder=1)
        ax.plot([AXIS_X - 0.08, AXIS_X], [event_y, label_y], color="#9A9A9A", lw=0.7, zorder=2)
        ax.text(AXIS_X + 0.12, label_y, label, va="center", ha="left", fontsize=7.6, color="#2F2F2F")
    ax.vlines(AXIS_X, ymin=0.0, ymax=1.06, colors="#8F8F8F", lw=0.9, zorder=1)


def add_arrow(ax, x_from: float, x_to: float, y: float) -> None:
    ax.annotate("", xy=(x_to, y), xytext=(x_from, y),
                arrowprops=dict(arrowstyle="-|>", lw=1.3, color="#CC0000", mutation_scale=11),
                zorder=4)


def draw_migrations(ax, model: str, params: dict[str, float], nodes) -> None:
    flows = MODEL_FLOWS[model]
    if not flows:
        return
    max_time = params["TIME4"]
    levels = {4: [0.30, 0.48, 0.66, 0.84], 2: [0.35, 0.65], 1: [0.50]}[len(flows)]
    for (src, dst, stop_event), frac in zip(flows, levels):
        y_stop = scaled_y(params[stop_event], max_time)
        add_arrow(ax, POP_X[src], POP_X[dst], max(frac * y_stop, 0.04))


def draw_panel(ax, model: str, params: dict[str, float]):
    nodes = tree_nodes(params)
    draw_event_guides(ax, params)
    draw_tree(ax, nodes)
    draw_migrations(ax, model, params, nodes)
    for pop_name, x_value in POP_X.items():
        ax.text(x_value, -0.065, pop_name, ha="center", va="top", fontsize=8.5)
    ax.set_title(model, fontsize=12.5, pad=4)
    ax.set_xlim(0.25, 7.15)
    ax.set_ylim(-0.1, 1.12)
    ax.axis("off")
    return model_events(params)


def main() -> None:
    fig = plt.figure(figsize=(15.6, 9.4))
    grid = fig.add_gridspec(2, 6, wspace=0.28, hspace=0.10)
    axes = [
        fig.add_subplot(grid[0, 0:2]),
        fig.add_subplot(grid[0, 2:4]),
        fig.add_subplot(grid[0, 4:6]),
        fig.add_subplot(grid[1, 0:2]),
        fig.add_subplot(grid[1, 2:4]),
        fig.add_subplot(grid[1, 4:6]),
    ]

    event_table = []
    for ax, (model, path) in zip(axes, MODEL_FILES.items()):
        params = read_bestlhoods(path)
        for event_name, generations in draw_panel(ax, model, params):
            event_table.append({
                "model": model,
                "event": event_name,
                "time_generations": round(generations, 3),
                "time_years": round(generations * GENERATION_TIME, 3),
            })

    fig.text(0.5, 0.028,
             "Black branches show split order; gray dashed lines show event order; red arrows show forward gene flow.",
             ha="center", fontsize=9.5)
    fig.subplots_adjust(left=0.03, right=0.98, top=0.96, bottom=0.08)

    prefix = BASE / "six_model_schematics_top5_M01"
    fig.savefig(prefix.with_suffix(".pdf"))
    fig.savefig(prefix.with_suffix(".svg"))
    fig.savefig(prefix.with_suffix(".png"), dpi=300)
    plt.close(fig)

    tsv_path = prefix.with_name(prefix.name + "_event_times.tsv")
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["model", "event", "time_generations", "time_years"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(event_table)

    for suffix in (".pdf", ".svg", ".png"):
        print(f"Wrote {prefix.with_suffix(suffix)}")
    print(f"Wrote {tsv_path}")


if __name__ == "__main__":
    main()
