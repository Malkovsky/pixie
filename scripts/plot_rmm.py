"""
Plot RmMTree benchmark results (optionally compared with sdsl-lite).

This script reads a JSON produced by `bench_rmm` and, optionally, a JSON
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


def load_bench_json(path: str) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if isinstance(data, dict) and "benchmarks" in data:
        data = data["benchmarks"]
    df = pd.DataFrame(data)
    required = {"name", "cpu_time", "N"}
    if not required.issubset(df.columns):
        missing = ", ".join(sorted(required - set(df.columns)))
        raise ValueError(
            f"{path!r}: missing required columns in benchmark JSON: {missing}"
        )
    df = df.dropna(subset=["cpu_time", "N"])
    return df


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
        help="Path to the JSON file with RmMTree results (output of bench_rmm).",
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

    rmm_ops = set(df_rmm["name"].dropna().astype(str).unique())
    if df_sdsl is not None:
        sdsl_ops = set(df_sdsl["name"].dropna().astype(str).unique())
        common_ops = rmm_ops & sdsl_ops
        ops_to_plot = [op for op in OPS_ORDER if op in common_ops]
        extra = sorted(common_ops - set(OPS_ORDER))
        ops_to_plot.extend(extra)
    else:
        ops_to_plot = OPS_ORDER

    for op in ops_to_plot:
        d_rmm = df_rmm[df_rmm["name"] == op].copy()
        d_sdsl = None
        if df_sdsl is not None:
            d_sdsl = df_sdsl[df_sdsl["name"] == op].copy()
        if d_rmm.empty or (df_sdsl is not None and (d_sdsl is None or d_sdsl.empty)):
            continue
        plt.figure()

        if not d_rmm.empty:
            d_rmm = (
                d_rmm.groupby("N", as_index=False)["cpu_time"].median().sort_values("N")
            )
            y_rmm = d_rmm["cpu_time"]
            if args.smooth and args.smooth > 1:
                d_rmm["ns_smooth"] = (
                    d_rmm["cpu_time"]
                    .rolling(window=args.smooth, center=True, min_periods=1)
                    .median()
                )
                y_rmm = d_rmm["ns_smooth"]
            plt.scatter(
                d_rmm["N"],
                d_rmm["cpu_time"],
                s=8,
                alpha=0.3,
                linewidths=0,
            )
            plt.plot(
                d_rmm["N"],
                y_rmm,
                linewidth=1.5,
                label="RmMTree",
            )

        if d_sdsl is not None and not d_sdsl.empty:
            d_sdsl = (
                d_sdsl.groupby("N", as_index=False)["cpu_time"]
                .median()
                .sort_values("N")
            )
            y_sdsl = d_sdsl["cpu_time"]
            if args.smooth and args.smooth > 1:
                d_sdsl["ns_smooth"] = (
                    d_sdsl["cpu_time"]
                    .rolling(window=args.smooth, center=True, min_periods=1)
                    .median()
                )
                y_sdsl = d_sdsl["ns_smooth"]
            plt.scatter(
                d_sdsl["N"],
                d_sdsl["cpu_time"],
                s=8,
                alpha=0.3,
                linewidths=0.9,
                marker="x",
            )
            plt.plot(
                d_sdsl["N"],
                y_sdsl,
                linewidth=1.5,
                linestyle="--",
                label="sdsl-lite",
            )

        if args.logx:
            plt.xscale("log", base=2)
        plt.xlabel("Sequence size N (bits)")
        plt.ylabel("Time per operation, ns")
        if d_sdsl is not None and not d_sdsl.empty:
            plt.title(f"RmMTree vs sdsl-lite: {op}")
        else:
            plt.title(f"RmMTree: {op}")
        plt.grid(True, which="both", linestyle="--", alpha=0.4)
        plt.legend(loc="best")
        out = os.path.join(args.save_dir, f"{op}.png")
        plt.savefig(out, bbox_inches="tight", dpi=160)
        if not args.show:
            plt.close()
        print(f"[saved] {out}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
