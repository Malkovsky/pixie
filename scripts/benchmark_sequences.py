#!/usr/bin/env python3
"""Run every registered sequence benchmark serially, with resumable JSON output.

Results belong outside the repository. A completed group is reused only when
its executable hash, case list, and run settings match the current invocation.
"""

import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys


BINARIES = (
    "permutable_sequence_benchmarks",
    "permutation_benchmarks",
    "bit_sequence_benchmarks",
    "sorted_sequence_benchmarks",
)


def benchmark_literal(name):
    # Google Benchmark's regex dialect rejects Python's escaped spaces.
    return ''.join('\\' + char if char in r'\.^$|?*+()[]{}' else char
                   for char in name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cpu", type=int, default=0)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--min-time", default="0.05s")
    parser.add_argument("--timeout", type=float, default=600)
    parser.add_argument("--binary", action="append", choices=BINARIES)
    parser.add_argument("--filter", help="Python regex selecting registered case names")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--list-only", action="store_true")
    parser.add_argument("--lock", type=Path,
                        default=Path("/tmp/pixie-sequence-benchmarks.lock"))
    args = parser.parse_args()
    if args.repetitions < 1 or args.timeout <= 0:
        parser.error("repetitions and timeout must be positive")

    groups = []
    for name in args.binary or BINARIES:
        binary = (args.build_dir / name).resolve()
        cases = subprocess.check_output(
            [str(binary), "--benchmark_list_tests"], text=True).splitlines()
        cases = [case for case in cases if case.strip()]
        if not cases:
            raise RuntimeError(f"No registered cases: {binary}")
        if args.filter:
            cases = [case for case in cases if re.search(args.filter, case)]
        digest = hashlib.sha256(binary.read_bytes()).hexdigest()
        partitions = {}
        for case in cases:
            operation = case.split("/", 1)[0]
            size = re.search(r"/N:(\d+)(?:/|$)", case)
            key = (operation, size.group(1) if size else None)
            partitions.setdefault(key, []).append(case)
        for (operation, size), selected in partitions.items():
            pattern = "^(" + "|".join(benchmark_literal(case) for case in selected) + ")$"
            groups.append({"binary": str(binary), "sha256": digest,
                           "operation": operation, "N": size,
                           "cases": selected, "filter": pattern,
                           "cpu": args.cpu, "repetitions": args.repetitions,
                           "min_time": args.min_time,
                           "environment": {key: value for key, value in os.environ.items()
                                           if key.startswith("PIXIE_")}})
        print(f"{name}: {len(cases)} cases, {len(partitions)} groups", flush=True)
    if not groups:
        raise RuntimeError("No registered cases matched the selection")
    if args.list_only:
        return 0

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.lock.parent.mkdir(parents=True, exist_ok=True)
    errors = []
    completed_cases = 0
    with args.lock.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        for index, group in enumerate(groups):
            key = hashlib.sha256(json.dumps(group, sort_keys=True).encode()).hexdigest()[:20]
            output = args.output_dir / f"{key}.json"
            metadata = args.output_dir / f"{key}.meta.json"
            log = args.output_dir / f"{key}.log"
            cached = (args.resume and output.exists() and metadata.exists()
                      and json.loads(metadata.read_text()) == group)
            if not cached:
                command = ["taskset", "-c", str(args.cpu), group["binary"],
                           "--benchmark_filter=" + group["filter"],
                           "--benchmark_repetitions=" + str(args.repetitions),
                           "--benchmark_min_time=" + args.min_time,
                           "--benchmark_enable_random_interleaving=true",
                           "--benchmark_report_aggregates_only=true",
                           "--benchmark_out=" + str(output.resolve()),
                           "--benchmark_out_format=json"]
                print(f"[{index + 1}/{len(groups)}] {Path(group['binary']).name} "
                      f"{group['operation']} N={group['N']} ({len(group['cases'])} cases)",
                      flush=True)
                with log.open("w") as stream:
                    subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                   check=True, timeout=args.timeout)
            if not output.exists() or output.stat().st_size == 0:
                raise RuntimeError(f"Benchmark produced no JSON; inspect {log}")
            data = json.loads(output.read_text())
            actual = {row.get("run_name", row["name"]) for row in data["benchmarks"]}
            if actual != set(group["cases"]):
                raise RuntimeError(f"Incomplete or unexpected cases in {output}")
            metadata.write_text(json.dumps(group, indent=2) + "\n")
            for row in data["benchmarks"]:
                if row.get("error_occurred"):
                    errors.append({"name": row["name"], "error": row.get("error_message"),
                                   "file": str(output)})
            completed_cases += len(group["cases"])
            print(f"  completed: {completed_cases} cases; errors/skips: {len(errors)}",
                  flush=True)
    summary = {"groups": len(groups), "cases": completed_cases, "errors": errors}
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    return 1 if errors else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (subprocess.SubprocessError, RuntimeError, OSError, ValueError) as error:
        print(f"Benchmark collection stopped: {error}", file=sys.stderr)
        sys.exit(2)
