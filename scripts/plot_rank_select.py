#!/usr/bin/env python3
"""Plot rank/select timings and auxiliary space from Google Benchmark JSON.

Produces PNG/SVG figures with one panel per input family and a normalized CSV.
Stored mean/median aggregates take precedence over iterations for each point;
otherwise the selected statistic is calculated from the iterations. Multiple
files are reduced separately, then their point values receive equal weight.
Requires Matplotlib. Example:

  python3 scripts/plot_rank_select.py results.json -o build/rank-select-plots
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import statistics
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

OPERATIONS = ("rank1", "rank0", "select1", "select0", "build")
FAMILIES = (
    "12p5", "50", "87p5", "sparse", "sparse0", "runs", "mixed", "mixed0", "fnbp",
)
FAMILY_LABELS = {
    "12p5": "12.5% ones",
    "50": "50% ones",
    "87p5": "87.5% ones",
    "sparse": "Sparse ones",
    "sparse0": "Sparse zeros",
    "runs": "Runs",
    "mixed": "Mixed",
    "mixed0": "Mixed complement",
    "fnbp": "FNBP",
}
BACKEND_LABELS = {
    "rank_select": "Baseline",
    "rank_select_compressed": "Compressed",
    "rank_select_pasta": "Pasta",
    "rank_select_3star": "3-star",
    "rank_select_succinct_fast": "FastBitVec",
    "rank_select_orzo": "Orzo",
    "rank_select_sdsl": "SDSL (rank-v + MCL)",
    "rank_select_spider": "Spider",
}
OP_PATTERN = r"rank[01]|select[01]|build_both"
REGULAR_NAME = re.compile(rf"^(rank_select(?:_\w+)?)_({OP_PATTERN})_(\w+)$")
FNBP_NAME = re.compile(rf"^(rank_select(?:_\w+)?)_fnbp_({OP_PATTERN})$")
SIZE = re.compile(r"/N:(\d+)(?:/|_|$)")
AGGREGATE_SUFFIX = re.compile(r"_(mean|median|stddev|cv)$")
TIME_TO_NS = {"ns": 1.0, "us": 1_000.0, "ms": 1_000_000.0, "s": 1_000_000_000.0}


def parse_name(row):
    """Return (operation, family, backend, size), or None for unrelated rows."""
    name = row.get("run_name") or row.get("name", "")
    base = name.split("/", 1)[0]
    if "/" not in name:
        base = AGGREGATE_SUFFIX.sub("", base)
    match = FNBP_NAME.fullmatch(base)
    if match:
        backend, operation = match.groups()
        family = "fnbp"
    else:
        match = REGULAR_NAME.fullmatch(base)
        if not match:
            return None
        backend, operation, family = match.groups()
    size_match = SIZE.search(name)
    size = int(size_match[1]) if size_match else int(row.get("N", 0))
    if size <= 0:
        raise ValueError("missing or invalid input size N")
    operation = "build" if operation == "build_both" else operation
    return operation, family, backend, size


def auxiliary_percent(row):
    """Return index bytes / input bytes in percent, excluding source copies."""
    if "aux_bytes" in row and "input_bytes" in row:
        auxiliary, input_bytes = float(row["aux_bytes"]), float(row["input_bytes"])
        if not math.isfinite(input_bytes) or input_bytes <= 0:
            raise ValueError("input_bytes must be positive and finite")
        value = 100 * auxiliary / input_bytes
    elif "aux_bits_per_input_bit" in row:
        value = 100 * float(row["aux_bits_per_input_bit"])
    else:
        return None
    if not math.isfinite(value) or value < 0:
        raise ValueError("auxiliary space must be nonnegative and finite")
    return value


def collect_file(document, metric="wall", aggregate="median", quantity="time"):
    """Return timing points in ns or auxiliary points in percent, plus warnings."""
    if isinstance(document, dict):
        if document.get("context", {}).get("dry_run"):
            raise ValueError("dry-run JSON contains smoke timings, not measurements")
        rows = document.get("benchmarks")
    else:
        rows = document
    if not isinstance(rows, list):
        raise ValueError("expected a Google Benchmark object or a benchmark-row list")
    selected, iterations = defaultdict(list), defaultdict(list)
    warnings = []
    for row in rows:
        if not isinstance(row, dict):
            warnings.append("non-object benchmark row")
            continue
        name = row.get("name", "")
        if row.get("error_occurred"):
            if name.startswith("rank_select"):
                warnings.append(f"{name}: {row.get('error_message', 'benchmark failed')}")
            continue
        suffix = AGGREGATE_SUFFIX.search(name)
        kind = row.get("aggregate_name") or (suffix[1] if suffix else None)
        if kind and kind != aggregate:
            continue
        if row.get("run_type") == "aggregate" and kind != aggregate:
            continue
        try:
            key = parse_name(row)
            if key is None:
                continue
            if quantity == "auxiliary":
                value = auxiliary_percent(row)
                if value is None:
                    continue
            else:
                unit = row.get("time_unit", "ns")
                if unit not in TIME_TO_NS:
                    raise ValueError(f"unsupported time unit {unit!r}")
                value = float(row["real_time" if metric == "wall" else "cpu_time"])
                value *= TIME_TO_NS[unit]
                if not math.isfinite(value) or value <= 0:
                    raise ValueError("timing must be positive and finite")
        except (KeyError, TypeError, ValueError, OverflowError) as error:
            warnings.append(f"{name}: {error}")
            continue
        (selected if kind else iterations)[key].append(value)
    reduce = statistics.median if aggregate == "median" else statistics.mean
    points = {
        key: reduce(selected.get(key) or iterations[key])
        for key in selected.keys() | iterations.keys()
    }
    if quantity == "auxiliary":
        # Use build counters when available; query-only captures can also supply them.
        priority = ("build", "rank1", "rank0", "select1", "select0")
        space = {}
        for key, value in sorted(points.items(), key=lambda item: priority.index(item[0][0])):
            space.setdefault(("auxiliary", *key[1:]), value)
        points = space
    return points, warnings


def load_points(paths, metric, aggregate, include_auxiliary=False):
    combined = defaultdict(list)
    reduce = statistics.median if aggregate == "median" else statistics.mean
    for path in paths:
        document = json.loads(path.read_text(encoding="utf-8"))
        points, warnings = collect_file(document, metric, aggregate)
        if include_auxiliary:
            space, space_warnings = collect_file(document, aggregate=aggregate,
                                               quantity="auxiliary")
            points.update(space)
            warnings.extend(space_warnings)
        for warning in dict.fromkeys(warnings):
            print(f"warning: {path}: {warning}", file=sys.stderr)
        for key, value in points.items():
            combined[key].append(value)
    return {key: reduce(values) for key, values in combined.items()}


def size_label(size):
    if size & (size - 1) == 0:
        return rf"$2^{{{size.bit_length() - 1}}}$"
    return f"{size:,}"


def write_graphs(points, output, metric, aggregate, formats, yscale):
    # Set a writable font-cache location and a headless backend before pyplot.
    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "pixie-matplotlib")
    )
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullLocator, StrMethodFormatter
    import plot_style as style

    style.apply_style()
    families = sorted(
        {key[1] for key in points},
        key=lambda f: (FAMILIES.index(f) if f in FAMILIES else len(FAMILIES), f),
    )
    backend_order = tuple(BACKEND_LABELS)
    backends = sorted(
        {key[2] for key in points},
        key=lambda b: (
            backend_order.index(b) if b in backend_order else len(backend_order), b
        ),
    )
    styles = {
        backend: style.line_style(
            "baseline" if backend == "rank_select" else backend.removeprefix("rank_select_"),
            index,
        )
        for index, backend in enumerate(backends)
    }
    columns = min(style.PANEL_COLUMNS, len(families))
    rows = math.ceil(len(families) / columns)
    for operation in (*OPERATIONS, "auxiliary"):
        fig, axes = plt.subplots(
            rows, columns, squeeze=False, sharey=True,
            figsize=(columns * style.PANEL_SIZE[0], rows * style.PANEL_SIZE[1]),
        )
        handles = {}
        unit = "ms" if operation == "build" else "ns"
        divisor = 1_000_000.0 if operation == "build" else 1.0
        space = operation == "auxiliary"
        axis_scale = "linear" if space and any(
            value == 0 for key, value in points.items() if key[0] == operation
        ) else yscale
        for ax, family in zip(axes.flat, families):
            family_sizes = set()
            for backend in backends:
                values = sorted(
                    (key[3], value / divisor)
                    for key, value in points.items()
                    if key[:3] == (operation, family, backend)
                )
                if not values:
                    continue
                sizes, timings = zip(*values)
                family_sizes.update(sizes)
                line, = ax.plot(
                    sizes, timings,
                    label=BACKEND_LABELS.get(backend, backend), **styles[backend],
                )
                handles[backend] = line
            ax.set_title(FAMILY_LABELS.get(family, family))
            ax.set_xlabel("Input size (bits)")
            action = "build" if operation == "build" else "query"
            ylabel = f"{metric.title()} time per {action}"
            if operation != "build" or axis_scale != "log":
                ylabel += f" ({unit})"
            ax.set_ylabel("Auxiliary index / input (%)" if space else ylabel)
            if family_sizes:
                ax.set_xscale("log", base=2)
                if axis_scale == "log":
                    if operation == "build":
                        style.build_time_y_axis(ax)
                    else:
                        style.log_y_axis(ax, percent=space)
                else:
                    ax.set_yscale(axis_scale)
                    if space:
                        ax.yaxis.set_major_formatter(StrMethodFormatter("{x:g}%"))
                ticks = sorted(family_sizes)
                ax.set_xticks(ticks, [size_label(n) for n in ticks])
                ax.tick_params(axis="y", labelleft=True)
                ax.xaxis.set_minor_locator(NullLocator())
                if operation != "build" or axis_scale != "log":
                    ax.yaxis.set_minor_locator(NullLocator())
            else:
                ax.text(
                    0.5, 0.5, "No measurements", ha="center", va="center",
                    transform=ax.transAxes,
                )
                ax.set_xticks([])
                ax.tick_params(axis="y", left=False, labelleft=False)
        for ax in list(axes.flat)[len(families):]:
            ax.set_visible(False)
        title = "Construction" if operation == "build" else operation.capitalize()
        heading = "Auxiliary index space" if space else f"{title} · {aggregate} {metric} time"
        fig.suptitle(f"Rank/select — {heading}", y=0.995)
        legend_columns = min(style.LEGEND_COLUMNS, columns * 2)
        if handles:
            ordered = [handles[b] for b in backends if b in handles]
            fig.legend(
                ordered, [line.get_label() for line in ordered],
                loc="lower center", ncol=min(len(ordered), legend_columns),
            )
        else:
            print(f"warning: no measurements for {operation}", file=sys.stderr)
        legend_rows = math.ceil(len(handles) / legend_columns) if handles else 0
        bottom = (0.35 + legend_rows * 0.25) / (rows * style.PANEL_SIZE[1])
        fig.tight_layout(rect=(0, bottom, 1, 0.95))
        for fmt in formats:
            fig.savefig(output / f"{operation}.{fmt}", dpi=style.DPI)
        plt.close(fig)


def write_csv(points, output, metric, aggregate):
    with (output / "data.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "operation", "family", "backend", "input_bits",
            "metric", "aggregate", "value", "unit",
        ])
        for (operation, family, backend, size), value in sorted(points.items()):
            build = operation == "build"
            space = operation == "auxiliary"
            writer.writerow([
                operation, family, backend, size, "space" if space else metric, aggregate,
                value / (1_000_000 if build else 1), "%" if space else "ms" if build else "ns",
            ])


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("json", nargs="+", type=Path, help="Google Benchmark JSON file(s)")
    parser.add_argument(
        "-o", "--output-dir", type=Path, default=Path("build/rank-select-plots")
    )
    parser.add_argument("--metric", choices=("wall", "cpu"), default="wall")
    parser.add_argument("--aggregate", choices=("median", "mean"), default="median")
    parser.add_argument(
        "--format", nargs="+", choices=("png", "svg", "pdf"), default=["png", "svg"]
    )
    parser.add_argument("--yscale", choices=("log", "linear"), default="log")
    parser.add_argument(
        "--family", action="append",
        help="Only plot this input family; repeat to select several",
    )
    args = parser.parse_args(argv)
    try:
        points = load_points(args.json, args.metric, args.aggregate, include_auxiliary=True)
        if args.family:
            points = {
                key: value for key, value in points.items() if key[1] in args.family
            }
        if not points:
            parser.error("no valid rank/select timing rows match the inputs and filters")
        args.output_dir.mkdir(parents=True, exist_ok=True)
        write_graphs(
            points, args.output_dir, args.metric, args.aggregate,
            args.format, args.yscale,
        )
        write_csv(points, args.output_dir, args.metric, args.aggregate)
    except (OSError, ValueError, ImportError) as error:
        parser.exit(1, f"error: {error}\n")
    print(
        f"Wrote six figures ({', '.join(args.format)}) and data.csv to "
        f"{args.output_dir}; {len(points)} plotted points."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
