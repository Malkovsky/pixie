"""
Plot RmMTree benchmark results (optionally compared with sdsl-lite).

This script reads a JSON produced by `rmm_benchmarks` and, optionally, a JSON
with benchmarks for sdsl-lite. For each operation, it draws a scatter plot
of individual points and a trend line (optionally median-smoothed) for each
implementation on the same figure.

Plots are saved as PNG files and can also be shown interactively.

Examples:
  python3 plot_rmm.py rmm_bench.json --save-dir plots --logx --smooth 3
  python3 plot_rmm.py rmm_bench.json --show
  python3 plot_rmm.py rmm_bench.json --sdsl-json rmm_bench_sdsl.json --save-dir plots --logx --smooth 3
"""

import argparse
import json
import os
import tempfile
from collections import defaultdict
from statistics import median

os.environ.setdefault(
    "MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "pixie-matplotlib")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OPS_ORDER = [
    "rank1",
    "rank0",
    "select1",
    "select0",
    "rank10",
    "select10",
    "excess",
    "fwdsearch",
    "bwdsearch",
    "range_min_query_pos",
    "range_min_query_val",
    "range_max_query_pos",
    "range_max_query_val",
    "mincount",
    "minselect",
    "close",
    "open",
    "enclose",
]


def clean_name(name: str) -> str:
    for suffix in ("_mean", "_median", "_stddev", "_cv"):
        name = name.removesuffix(suffix)
    return name


def load_bench_json(path: str) -> list[dict]:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if isinstance(data, dict) and "benchmarks" in data:
        data = data["benchmarks"]
    rows = []
    for bench in data:
        if not {"name", "cpu_time", "N"}.issubset(bench):
            continue
        if bench.get("run_type") == "aggregate" and bench.get("aggregate_name") != "mean":
            continue
        try:
            n = int(float(bench["N"]))
            cpu_time = float(bench["cpu_time"])
        except (TypeError, ValueError):
            continue
        if n <= 0 or cpu_time <= 0:
            continue
        rows.append({"name": clean_name(str(bench["name"])), "N": n, "cpu_time": cpu_time})
    if not rows:
        raise ValueError(f"{path!r}: no usable benchmark rows found")
    return rows


def ops_in(rows: list[dict]) -> set[str]:
    return {row["name"] for row in rows}


def op_points(rows: list[dict], op: str) -> list[tuple[int, float]]:
    grouped = defaultdict(list)
    for row in rows:
        if row["name"] == op:
            grouped[row["N"]].append(row["cpu_time"])
    return sorted((n, median(values)) for n, values in grouped.items())


def smooth(points: list[tuple[int, float]], window: int) -> list[float]:
    if window <= 1:
        return [time for _, time in points]
    values = [time for _, time in points]
    radius = window // 2
    smoothed = []
    for index in range(len(values)):
        lo = max(0, index - radius)
        hi = min(len(values), index + radius + 1)
        smoothed.append(median(values[lo:hi]))
    return smoothed


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Read JSON with RmMTree benchmark results (and optionally "
            "sdsl-lite results) and plot time per operation versus "
            "sequence size N for each operation."
        ),
        epilog=(
            "Examples:\n"
            "  python3 plot_rmm.py rmm_bench.json --save-dir plots --logx --smooth 3\n"
            "  python3 plot_rmm.py rmm_bench.json --show\n"
            "  python3 plot_rmm.py rmm_bench.json --sdsl-json rmm_bench_sdsl.json --save-dir plots --logx --smooth 3"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "json",
        metavar="JSON",
        help="Path to the JSON file with RmMTree results (output of rmm_benchmarks).",
    )
    ap.add_argument(
        "--sdsl-json",
        metavar="SDSL_JSON",
        help=(
            "Path to JSON file with sdsl-lite benchmark results. "
            "If provided, both implementations will be plotted on the same graphs."
        ),
    )
    ap.add_argument(
        "--save-dir",
        default="plots",
        metavar="DIR",
        help=(
            "Directory to save PNG plots. Will be created if it doesn't exist. "
            "Default: %(default)s"
        ),
    )
    ap.add_argument(
        "--show",
        action="store_true",
        help=(
            "Show plot windows after saving. "
            "By default, plots are only written to disk."
        ),
    )
    ap.add_argument(
        "--logx",
        action="store_true",
        help=(
            "Use a logarithmic X axis (base 2). Handy when N grows in powers of two."
        ),
    )
    ap.add_argument(
        "--smooth",
        type=int,
        default=0,
        metavar="W",
        help=(
            "Median smoothing window size for the trend line. "
            "0 or 1 means no smoothing. Default: %(default)s"
        ),
    )
    args = ap.parse_args()

    os.makedirs(args.save_dir, exist_ok=True)

    df_rmm = load_bench_json(args.json)
    df_sdsl = None
    if args.sdsl_json is not None:
        df_sdsl = load_bench_json(args.sdsl_json)

    rmm_ops = ops_in(df_rmm)
    if df_sdsl is not None:
        sdsl_ops = ops_in(df_sdsl)
        common_ops = rmm_ops & sdsl_ops
        ops_to_plot = [op for op in OPS_ORDER if op in common_ops]
        extra = sorted(common_ops - set(OPS_ORDER))
        ops_to_plot.extend(extra)
    else:
        ops_to_plot = OPS_ORDER

    for op in ops_to_plot:
        d_rmm = op_points(df_rmm, op)
        d_sdsl = []
        if df_sdsl is not None:
            d_sdsl = op_points(df_sdsl, op)
        if not d_rmm or (df_sdsl is not None and not d_sdsl):
            continue
        plt.figure()

        if d_rmm:
            x_rmm = [n for n, _ in d_rmm]
            raw_rmm = [time for _, time in d_rmm]
            y_rmm = smooth(d_rmm, args.smooth)
            plt.scatter(
                x_rmm,
                raw_rmm,
                s=18,
                alpha=0.7,
                linewidths=0,
            )
            plt.plot(
                x_rmm,
                y_rmm,
                marker="o",
                markersize=3,
                linewidth=0.85,
                label="Pixie RmMTree",
            )

        if d_sdsl:
            x_sdsl = [n for n, _ in d_sdsl]
            raw_sdsl = [time for _, time in d_sdsl]
            y_sdsl = smooth(d_sdsl, args.smooth)
            plt.scatter(
                x_sdsl,
                raw_sdsl,
                s=18,
                alpha=0.7,
                linewidths=0.85,
                marker="x",
            )
            plt.plot(
                x_sdsl,
                y_sdsl,
                marker="x",
                markersize=3,
                linewidth=0.85,
                linestyle="--",
                label="sdsl-lite",
            )

        if args.logx:
            plt.xscale("log", base=2)
        plt.xlabel("Sequence size N (bits)")
        plt.ylabel("Time per operation, ns")
        if d_sdsl:
            plt.title(f"RmMTree comparison: {op}")
        else:
            plt.title(f"RmMTree: {op}")
        plt.grid(True, which="both", linestyle="--", alpha=0.4)
        plt.legend(loc="best")
        out = os.path.join(args.save_dir, f"{op}.png")
        plt.savefig(out, bbox_inches="tight", dpi=180)
        if not args.show:
            plt.close()
        print(f"[saved] {out}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
