#!/usr/bin/env python3
"""Render excess-position benchmark JSON as a compact Markdown table."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, Tuple


def _target_from_benchmark(bench: dict) -> str:
    if "params" in bench and "X" in bench["params"]:
        return str(bench["params"]["X"])
    name = bench.get("name", "")
    for suffix in ("_mean", "_median", "_stddev", "_cv"):
        name = name.removesuffix(suffix)
    marker = "/X:"
    if marker in name:
        return name.split(marker, 1)[1].split("/", 1)[0]
    args = bench.get("args", [])
    if args:
        return str(args[0])
    return ""


def _method_from_benchmark(bench: dict) -> str:
    name = bench.get("name", "unknown")
    name = name.split("/", 1)[0]
    return name.removeprefix("BM_ExcessPositions512").lstrip("_") or "Current"


def _is_data_row(bench: dict) -> bool:
    name = bench.get("name", "")
    if not name.startswith("BM_ExcessPositions512"):
        return False
    if any(name.endswith(suffix) for suffix in ("_median", "_stddev", "_cv")):
        return False
    run_type = bench.get("run_type")
    return run_type in (None, "iteration", "aggregate")


def _collect(benchmarks: Iterable[dict]) -> Dict[Tuple[str, str], Tuple[float, float]]:
    rows: Dict[Tuple[str, str], Tuple[float, float]] = {}
    for bench in benchmarks:
        if not _is_data_row(bench):
            continue
        name = bench.get("name", "")
        if bench.get("run_type") == "aggregate" and not name.endswith("_mean"):
            continue
        method = _method_from_benchmark(bench).removesuffix("_mean")
        target = _target_from_benchmark(bench)
        if not target:
            continue
        real_time = float(bench.get("real_time", 0.0))
        cpu_time = float(bench.get("cpu_time", 0.0))
        rows[(method, target)] = (real_time, cpu_time)
    return rows


def _sort_target(value: str) -> int:
    try:
        return int(value)
    except ValueError:
        return 0


def render_markdown(rows: Dict[Tuple[str, str], Tuple[float, float]]) -> str:
    methods = sorted({method for method, _ in rows})
    targets = sorted({target for _, target in rows}, key=_sort_target)
    if not methods or not targets:
        raise SystemExit("No excess benchmark rows found.")

    lines = ["| Method | " + " | ".join(f"X={x}" for x in targets) + " |"]
    lines.append("|---|" + "|".join("---:" for _ in targets) + "|")
    for method in methods:
        cells = []
        for target in targets:
            value = rows.get((method, target))
            cells.append("" if value is None else f"{value[1]:.2f} ns")
        lines.append(f"| {method} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a Markdown CPU-time table from excess benchmark JSON."
    )
    parser.add_argument("json", type=Path, help="Google Benchmark JSON output")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Optional output Markdown file. Defaults to stdout.",
    )
    args = parser.parse_args()

    data = json.loads(args.json.read_text(encoding="utf-8"))
    table = render_markdown(_collect(data.get("benchmarks", [])))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(table + "\n", encoding="utf-8")
    else:
        print(table)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
