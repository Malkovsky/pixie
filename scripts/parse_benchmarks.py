#!/usr/bin/env python3
"""Parse Google Benchmark JSON and plot CPU time vs container size."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt


def _guess_size_key(benchmarks: Iterable[dict]) -> str:
    candidate_keys = ["tree_size", "size", "n"]
    for key in candidate_keys:
        for bench in benchmarks:
            if "params" in bench and key in bench["params"]:
                return key
            if "ArgNames" in bench and key in bench["ArgNames"]:
                return key
    return "size"


def _extract_size_from_name(name: str, size_key: str) -> int:
    token = f"{size_key}:"
    if token not in name:
        return 0
    try:
        chunk = name.split(token, 1)[1]
        value = chunk.split("/", 1)[0]
        return int(value)
    except (ValueError, IndexError):
        return 0


def _extract_size(bench: dict, size_key: str) -> int:
    if "params" in bench and size_key in bench["params"]:
        return int(bench["params"][size_key])
    if "args" in bench and bench["args"]:
        return int(bench["args"][0])
    name = bench.get("name", "")
    if name:
        return _extract_size_from_name(name, size_key)
    return 0


def _collect_points(
    benchmarks: Iterable[dict], size_key: str
) -> Dict[str, List[Tuple[int, float]]]:
    series: Dict[str, List[Tuple[int, float]]] = {}
    for bench in benchmarks:
        name = bench.get("name", "unknown")
        base_name = name.split("/", 1)[0] if name else "unknown"
        size = _extract_size(bench, size_key)
        cpu_time = float(bench.get("cpu_time", 0.0))
        if cpu_time <= 0:
            cpu_time = float(bench.get("real_time", 0.0))
        if size <= 0 or cpu_time <= 0:
            continue
        series.setdefault(base_name, []).append((size, cpu_time))
    for name in list(series.keys()):
        series[name] = sorted(series[name], key=lambda item: item[0])
    return series


def _plot_series(
    series: Dict[str, List[Tuple[int, float]]],
    output_path: Path,
    title: str,
    size_key: str,
    fmt: str,
) -> None:
    plt.figure(figsize=(10, 6))
    for name, points in series.items():
        sizes = [p[0] for p in points]
        times = [p[1] for p in points]
        plt.plot(sizes, times, marker="^", linewidth=1.0, label=name)

    plt.title(title)
    plt.xlabel(size_key.replace("_", " "))
    plt.ylabel("CPU time (ns)")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, format=fmt, dpi=160)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot CPU time vs container size from Google Benchmark JSON."
    )
    parser.add_argument("input", type=Path, help="Path to benchmark JSON file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("graphs/cpu_time_vs_size.png"),
        help="Output image path",
    )
    parser.add_argument(
        "--title",
        default="CPU time vs container size",
        help="Plot title",
    )
    parser.add_argument(
        "--size-key",
        default="",
        help="Parameter name for size (defaults to auto-detect)",
    )
    parser.add_argument(
        "--format",
        default="png",
        choices=["png", "svg"],
        help="Output image format",
    )
    args = parser.parse_args()

    data = json.loads(args.input.read_text(encoding="utf-8"))
    benchmarks = data.get("benchmarks", [])
    if not benchmarks:
        raise SystemExit("No benchmark data found in JSON.")

    size_key = args.size_key or _guess_size_key(benchmarks)
    series = _collect_points(benchmarks, size_key)
    if not series:
        raise SystemExit("No valid benchmark points found.")

    _plot_series(series, args.output, args.title, size_key, args.format)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
