#!/usr/bin/env python3
"""Plot Google Benchmark timing against evaluated bitvector/tree sizes."""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "pixie-matplotlib")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


Point = Tuple[int, float]


def _clean_name(name: str) -> str:
    for suffix in ("_mean", "_median", "_stddev", "_cv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def _is_timing_row(bench: dict) -> bool:
    name = bench.get("name", "")
    if any(name.endswith(suffix) for suffix in ("_stddev", "_cv")):
        return False
    run_type = bench.get("run_type")
    if run_type == "aggregate":
        return name.endswith("_mean")
    return run_type in (None, "iteration")


def _guess_size_key(benchmarks: Iterable[dict]) -> str:
    candidates = ("n", "N", "tree_size", "size")
    for bench in benchmarks:
        params = bench.get("params", {})
        for key in candidates:
            if key in params:
                return key
        for key in bench.get("ArgNames", []):
            if key in candidates:
                return key
        for key in candidates:
            if key in bench:
                return key
    return "n"


def _extract_size_from_name(name: str, size_key: str) -> int:
    marker = f"{size_key}:"
    if marker not in name:
        return 0
    try:
        return int(name.split(marker, 1)[1].split("/", 1)[0])
    except ValueError:
        return 0


def _extract_size(bench: dict, size_key: str) -> int:
    params = bench.get("params", {})
    if size_key in params:
        return int(params[size_key])
    if size_key in bench:
        return int(bench[size_key])
    args = bench.get("args", [])
    if args:
        return int(args[0])
    return _extract_size_from_name(bench.get("name", ""), size_key)


def _load_series(path: Path, size_key: str | None) -> Dict[str, List[Point]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    benchmarks = data.get("benchmarks", data if isinstance(data, list) else [])
    key = size_key or _guess_size_key(benchmarks)
    series: Dict[str, List[Point]] = {}
    for bench in benchmarks:
        if not _is_timing_row(bench):
            continue
        name = _clean_name(bench.get("name", "unknown")).split("/", 1)[0]
        size = _extract_size(bench, key)
        cpu_time = float(bench.get("cpu_time", bench.get("real_time", 0.0)))
        if size <= 0 or cpu_time <= 0:
            continue
        series.setdefault(name, []).append((size, cpu_time))

    collapsed: Dict[str, List[Point]] = {}
    for name, points in series.items():
        by_size: Dict[int, List[float]] = {}
        for size, time in points:
            by_size.setdefault(size, []).append(time)
        collapsed[name] = sorted(
            (size, sorted(times)[len(times) // 2]) for size, times in by_size.items()
        )
    return collapsed


def _plot(series_by_label: Dict[str, List[Point]], output: Path, title: str) -> None:
    plt.figure(figsize=(10, 6))
    for label, points in sorted(series_by_label.items()):
        xs = [size for size, _ in points]
        ys = [time for _, time in points]
        plt.plot(
            xs,
            ys,
            marker="o",
            markersize=3,
            linewidth=0.85,
            label=label,
        )
    plt.xscale("log", base=2)
    plt.xlabel("Evaluated size")
    plt.ylabel("CPU time (ns)")
    plt.title(title)
    plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.35)
    plt.legend(fontsize=7)
    plt.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output, dpi=180)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Plot bitvector/tree size benchmarks with x-log scale, thin lines, "
            "and markers at evaluated sizes."
        )
    )
    parser.add_argument("json", type=Path, nargs="+", help="Benchmark JSON file(s)")
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("graphs/size_benchmarks.png")
    )
    parser.add_argument("--title", default="Benchmark time vs size")
    parser.add_argument("--size-key", help="Override size key, e.g. n, N, tree_size")
    args = parser.parse_args()

    combined: Dict[str, List[Point]] = {}
    multi_file = len(args.json) > 1
    for path in args.json:
        for name, points in _load_series(path, args.size_key).items():
            label = f"{path.stem}:{name}" if multi_file else name
            combined[label] = points
    if not combined:
        raise SystemExit("No valid benchmark series found.")
    _plot(combined, args.output, args.title)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
