#!/usr/bin/env python3
"""Plot RMQ hybrid variant timings for fully random query rows.

The RMQ query benchmark uses the second positional argument as `max_width`.
When `max_width == N`, query widths are sampled uniformly from `[1, N]` and
left endpoints are sampled uniformly among valid positions for that width.
This script filters a Google Benchmark JSON file to those fully random rows and
plots the hybrid variants against input size.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import statistics
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Sequence, Tuple

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "pixie-matplotlib")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_METHODS = ("rmq_hybrid_btree", "rmq_cartesian_hybrid_btree")
METHOD_LABELS = {
    "rmq_hybrid_btree": "HybridBTree",
    "rmq_cartesian_hybrid_btree": "CartesianHybridBTree",
}
AGGREGATE_SUFFIXES = ("_mean", "_median", "_stddev", "_cv")
TIME_UNIT_TO_NS = {
    "ns": 1.0,
    "us": 1_000.0,
    "ms": 1_000_000.0,
    "s": 1_000_000_000.0,
}
TABLE_SEPARATOR_RE = re.compile(r"^-+$")


class Point(NamedTuple):
    size: int
    time_ns: float


def _clean_name(name: str) -> str:
    for suffix in AGGREGATE_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def _base_name(bench: dict) -> str:
    return _clean_name(bench.get("run_name") or bench.get("name", "")).split(
        "/", 1
    )[0]


def _positional_args_from_name(name: str) -> List[int]:
    args: List[int] = []
    parts = _clean_name(name).split("/")
    for part in parts[1:]:
        if ":" in part:
            continue
        try:
            args.append(int(part))
        except ValueError:
            continue
    return args


def _extract_size_width(bench: dict) -> Tuple[int, int]:
    params = bench.get("params", {})
    size = params.get("N", params.get("n", params.get("size")))
    width = params.get("max_width", params.get("width"))
    if size is not None and width is not None:
        return int(size), int(width)

    if "N" in bench and "max_width" in bench:
        return int(bench["N"]), int(bench["max_width"])

    positional_args = bench.get("args", [])
    if len(positional_args) >= 2:
        return int(positional_args[0]), int(positional_args[1])

    name_args = _positional_args_from_name(
        bench.get("run_name") or bench.get("name", "")
    )
    if len(name_args) >= 2:
        return name_args[0], name_args[1]
    return 0, 0


def _time_ns(bench: dict) -> float:
    time = float(bench.get("cpu_time", bench.get("real_time", 0.0)))
    return time * TIME_UNIT_TO_NS.get(bench.get("time_unit", "ns"), 1.0)


def _is_selected_aggregate(bench: dict, aggregate: str) -> bool:
    name = bench.get("name", "")
    if name.endswith(("_stddev", "_cv")):
        return False
    if bench.get("run_type") != "aggregate":
        return False
    return bench.get("aggregate_name", "") == aggregate or name.endswith(
        f"_{aggregate}"
    )


def _is_iteration_row(bench: dict) -> bool:
    name = bench.get("name", "")
    if any(name.endswith(suffix) for suffix in AGGREGATE_SUFFIXES):
        return False
    return bench.get("run_type") in (None, "iteration")


def _collect_points(
    benchmarks: Iterable[dict], methods: Sequence[str], aggregate: str
) -> Dict[str, List[Point]]:
    aggregate_rows: List[dict] = []
    iteration_rows: List[dict] = []
    method_set = set(methods)
    for bench in benchmarks:
        if _base_name(bench) not in method_set:
            continue
        if _is_selected_aggregate(bench, aggregate):
            aggregate_rows.append(bench)
        elif _is_iteration_row(bench):
            iteration_rows.append(bench)

    rows = aggregate_rows if aggregate_rows else iteration_rows
    grouped: Dict[str, Dict[int, List[float]]] = {}
    for bench in rows:
        size, width = _extract_size_width(bench)
        time_ns = _time_ns(bench)
        if size <= 0 or width != size or time_ns <= 0:
            continue
        method = _base_name(bench)
        grouped.setdefault(method, {}).setdefault(size, []).append(time_ns)

    series: Dict[str, List[Point]] = {}
    for method, by_size in grouped.items():
        label = METHOD_LABELS.get(method, method)
        series[label] = [
            Point(size, statistics.median(times))
            for size, times in sorted(by_size.items())
        ]
    return series


def _parse_int_cell(cell: str) -> int:
    value = cell.strip()
    if value.startswith("2^"):
        return 1 << int(value[2:])
    return int(value)


def _parse_time_cell(cell: str) -> float:
    value = cell.strip()
    if value == "-":
        return 0.0
    return float(value)


def _collect_points_from_rmq_header(
    path: Path, methods: Sequence[str]
) -> Dict[str, List[Point]]:
    selected = set(methods)
    by_label: Dict[str, List[Point]] = {
        "HybridBTree": [],
        "CartesianHybridBTree": [],
    }
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line.startswith("* |"):
            continue
        row = line[1:].strip()
        cells = [cell.strip() for cell in row.strip("|").split("|")]
        if len(cells) < 7:
            continue
        if cells[0] == "N" or TABLE_SEPARATOR_RE.match(cells[0]):
            continue
        try:
            size = _parse_int_cell(cells[0])
            max_width = _parse_int_cell(cells[1])
            hybrid_time = _parse_time_cell(cells[4])
            cartesian_hybrid_time = _parse_time_cell(cells[6])
        except ValueError:
            continue
        if max_width != size:
            continue
        if "rmq_hybrid_btree" in selected and hybrid_time > 0.0:
            by_label["HybridBTree"].append(Point(size, hybrid_time))
        if (
            "rmq_cartesian_hybrid_btree" in selected
            and cartesian_hybrid_time > 0.0
        ):
            by_label["CartesianHybridBTree"].append(
                Point(size, cartesian_hybrid_time)
            )

    return {
        label: sorted(points)
        for label, points in by_label.items()
        if points
    }


def _plot(series: Dict[str, List[Point]], output: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for label, points in sorted(series.items()):
        ax.plot(
            [point.size for point in points],
            [point.time_ns for point in points],
            marker="o",
            markersize=3,
            linewidth=0.8,
            label=label,
        )

    ax.set_xscale("log", base=2)
    ax.set_xlabel("Input size N")
    ax.set_ylabel("CPU time per query (ns)")
    ax.set_title(title)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.35)
    ax.legend(fontsize=8)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def _print_summary(series: Dict[str, List[Point]]) -> None:
    for label, points in sorted(series.items()):
        formatted = ", ".join(
            f"2^{point.size.bit_length() - 1}={point.time_ns:.2f}ns"
            if point.size > 0 and point.size & (point.size - 1) == 0
            else f"{point.size}={point.time_ns:.2f}ns"
            for point in points
        )
        print(f"{label}: {formatted}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Plot rmq_hybrid_btree and rmq_cartesian_hybrid_btree rows where "
            "max_width == N, i.e. fully random RMQ query widths."
        )
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Google Benchmark JSON file, or include/pixie/rmq.h with --from-rmq-header",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("graphs/rmq_hybrid_random_queries.png"),
        help="Output image path",
    )
    parser.add_argument(
        "--title",
        default="RMQ hybrid variants on fully random queries",
        help="Plot title",
    )
    parser.add_argument(
        "--aggregate",
        default="mean",
        choices=["mean", "median"],
        help="Aggregate row to use when repetitions are present",
    )
    parser.add_argument(
        "--method",
        action="append",
        choices=DEFAULT_METHODS,
        help="Benchmark method to include; defaults to both hybrid variants",
    )
    parser.add_argument(
        "--from-rmq-header",
        action="store_true",
        help="Parse the query benchmark table embedded at the top of include/pixie/rmq.h",
    )
    parser.add_argument(
        "--min-size",
        type=int,
        default=0,
        help="Drop points with N smaller than this value",
    )
    args = parser.parse_args()

    methods = tuple(args.method) if args.method else DEFAULT_METHODS
    if args.from_rmq_header:
        series = _collect_points_from_rmq_header(args.input, methods)
    else:
        data = json.loads(args.input.read_text(encoding="utf-8"))
        benchmarks = data.get("benchmarks", data if isinstance(data, list) else [])
        series = _collect_points(benchmarks, methods, args.aggregate)
    if args.min_size > 0:
        series = {
            label: [point for point in points if point.size >= args.min_size]
            for label, points in series.items()
        }
        series = {label: points for label, points in series.items() if points}
    if not series:
        raise SystemExit(
            "No fully random hybrid RMQ rows found. Expected benchmark names "
            "like rmq_hybrid_btree/N/N and rmq_cartesian_hybrid_btree/N/N."
        )
    _plot(series, args.output, args.title)
    _print_summary(series)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
