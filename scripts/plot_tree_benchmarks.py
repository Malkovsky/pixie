#!/usr/bin/env python3
"""Four canonical tree comparisons: fixed succinct configurations, no bands."""
import argparse
import csv
import json
import os
from pathlib import Path
import statistics

os.environ.setdefault("MPLCONFIGDIR", "/tmp/pixie-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_style import apply_style

SERIES = [
    ("adjlist_preorder", "Vectors · preorder", "#777777", "o", "-"),
    ("adjlist", "Vectors · BFS", "#777777", "o", "--"),
    ("pointer_preorder", "Pointers · preorder", "#222222", "s", "-"),
    ("pointer_bfs", "Pointers · BFS", "#222222", "s", "--"),
    ("bp_best", "BP", "#0072B2", "^", "-"),
    ("dfuds_best", "DFUDS", "#D55E00", "D", "-"),
    ("louds_best", "LOUDS", "#009E73", "v", "-"),
]
OPS = [("dfs", "DFS"), ("dfs_children", "DFS + child inspection"),
       ("random_path", "Root-to-leaf path"),
       ("random_path_children", "Path + child inspection")]
SCALE = {"ns": 1., "us": 1e3, "ms": 1e6, "s": 1e9}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("-o", "--output-dir", type=Path, required=True)
    parser.add_argument("--overlay", type=Path, action="append", default=[],
                        help="Additional succinct results; one budget per file")
    args = parser.parse_args()
    groups = {}
    budgets = {}
    limits = set()
    for row in json.loads(args.input.read_text())["benchmarks"]:
        if row.get("error_occurred"):
            raise ValueError("Benchmark input contains errors")
        if row.get("run_type") != "iteration":
            continue
        name = row.get("run_name", row["name"]).split("/")[0]
        nodes = int(row["size_nodes"])
        groups.setdefault((name, nodes), []).append(
            row["cpu_time"] * SCALE[row["time_unit"]])
        if "_best_" in name:
            used = row["cache_and_scratch_bytes"]
            limit = row["cache_limit_bytes"]
            if used > limit:
                raise ValueError(f"Cache budget exceeded: {name}/{nodes}")
            limits.add(limit)
            budgets[name, nodes] = used
    sizes = sorted({n for _, n in groups})
    if len(limits) != 1:
        raise ValueError("Expected a single cache limit")
    budget_kib = next(iter(limits)) / 1024
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    apply_style()
    records = []
    for op, title in OPS:
        fig, ax = plt.subplots(figsize=(10, 6))
        for prefix, label, color, marker, dash in SERIES:
            name = prefix + "_" + op
            values = []
            for n in sizes:
                samples = groups[name, n]
                median = statistics.median(samples)
                records.append([name, n, median, len(samples),
                                budgets.get((name, n), "")])
                values.append(median / n if op.startswith("dfs") else median)
            ax.plot(sizes, values, label=label, color=color,
                    marker=marker, linestyle=dash)
        for overlay_index, source in enumerate(args.overlay):
            extra = {}
            extra_used = {}
            for row in json.loads(source.read_text())["benchmarks"]:
                if row.get("error_occurred"):
                    raise ValueError("Overlay contains benchmark errors")
                if row.get("run_type") != "iteration":
                    continue
                name = row.get("run_name", row["name"]).split("/")[0]
                if "_best_" not in name:
                    continue
                if row["cache_and_scratch_bytes"] > row["cache_limit_bytes"]:
                    raise ValueError("Overlay exceeds cache budget")
                extra.setdefault((name, int(row["size_nodes"]),
                                  row["cache_limit_bytes"]), []).append(
                    row["cpu_time"] * SCALE[row["time_unit"]])
                extra_used[name, int(row["size_nodes"]), row["cache_limit_bytes"]] = row["cache_and_scratch_bytes"]
            for prefix, label, color, marker, _ in SERIES[4:]:
                points = sorted((n, limit, samples)
                                for (name, n, limit), samples in extra.items()
                                if name == prefix + "_" + op)
                if len({limit for _, limit, _ in points}) != 1:
                    raise ValueError("Expected one overlay budget per series")
                xs, ys = [], []
                for n, limit, samples in points:
                    median = statistics.median(samples)
                    xs.append(n)
                    ys.append(median / n if op.startswith("dfs") else median)
                    records.append([prefix + "_" + op + f"_{limit:g}B", n,
                                    median, len(samples),
                                    extra_used[prefix + "_" + op, n, limit]])
                ax.plot(xs, ys, color=color, marker=marker,
                        linestyle=(":", "-.")[overlay_index % 2],
                        markerfacecolor="white",
                        label=f"{label} · {points[0][1]/1024:g} KiB updated")
        ax.set_xscale("log", base=2)
        ax.set_xticks(sizes, [rf"$2^{{{n.bit_length()-1}}}$"
                     if n & (n - 1) == 0 else "10M" for n in sizes])
        ax.set_ylim(bottom=0)
        ax.set_xlabel("Tree size (nodes)")
        ax.set_ylabel("CPU time (ns/node)" if op.startswith("dfs")
                      else "CPU time (ns/path)")
        ax.set_title(title, loc="left", fontweight="bold")
        fig.suptitle((f"Base sweep: {budget_kib:g} KiB · updated budgets shown in legend"
                      if args.overlay else
                      f"Fixed succinct variants · {budget_kib:g} KiB cache + traversal buffers"),
                     fontsize=12)
        ax.legend(loc="upper left", bbox_to_anchor=(0, -.17), ncol=4)
        fig.text(.5, .015, "Median CPU time · construction excluded · lower is better",
                 ha="center", fontsize=9)
        fig.tight_layout(rect=(0, .08, 1, .97))
        for ext in ("png", "svg", "pdf"):
            fig.savefig(output / f"{op}.{ext}", dpi=180)
        plt.close(fig)
    with (output / "plot_data.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["benchmark", "nodes", "median_cpu_ns", "repetitions",
                         "cache_and_scratch_bytes"])
        writer.writerows(records)
    print(f"Saved four plots for {len(sizes)} sizes to {output}")


if __name__ == "__main__":
    main()
