"""
Plot RmMTree benchmark results.

This script reads a JSON produced by `bench_rmm` and,
for each operation, draws a scatter plot of individual points
and a trend line (optionally median-smoothed).
Plots are saved as PNG files and can also be shown interactively.

Examples:
  python3 plot_rmm.py rmm_bench.json --save-dir plots --logx --smooth 3
  python3 plot_rmm.py rmm_bench.json --show
"""

import argparse
import json
import os
import pandas as pd
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


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Read a JSON with RmMTree benchmark results and plot time per operation "
            "versus sequence size N for each operation."
        ),
        epilog=(
            "Examples:\n"
            "  python3 plot_rmm.py rmm_bench.json --save-dir plots --logx --smooth 3\n"
            "  python3 plot_rmm.py rmm_bench.json --show"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "json",
        metavar="JSON",
        help="Path to the JSON file with results (output of bench_rmm).",
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
    with open(args.json, "r", encoding="utf-8") as f:
        data = json.load(f)
    if isinstance(data, dict) and "benchmarks" in data:
        data = data["benchmarks"]
    df = pd.DataFrame(data)
    df = df.dropna(subset=["cpu_time", "N"])

    for op in OPS_ORDER:
        d = df[df["name"] == op].copy()
        if d.empty:
            continue
        d = d.groupby("N", as_index=False)["cpu_time"].median().sort_values("N")

        yplot = d["cpu_time"]
        if args.smooth and args.smooth > 1:
            d["ns_smooth"] = (
                d["cpu_time"]
                .rolling(window=args.smooth, center=True, min_periods=1)
                .median()
            )
            yplot = d["ns_smooth"]

        plt.figure()
        plt.scatter(d["N"], d["cpu_time"], s=8, alpha=0.3, linewidths=0)
        plt.plot(d["N"], yplot, linewidth=1.5)

        if args.logx:
            plt.xscale("log", base=2)
        plt.xlabel("Sequence size N (bits)")
        plt.ylabel("Time per operation, ns")
        plt.title(f"RmMTree: {op}")
        plt.grid(True, which="both", linestyle="--", alpha=0.4)
        out = os.path.join(args.save_dir, f"{op}.png")
        plt.savefig(out, bbox_inches="tight", dpi=160)
        if not args.show:
            plt.close()
        print(f"[saved] {out}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
