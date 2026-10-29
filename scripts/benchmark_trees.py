#!/usr/bin/env python3
"""Collect resumable tree benchmarks with raw repetitions and CSV exports."""

import argparse
import csv
import fcntl
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import time
import zipfile


def atomic_json(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n")
    temporary.replace(path)


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def groups(scope):
    result = [] if scope == "succinct" else [
        ("bp", "adjlist_preorder", False),
        ("bp", "pointer_preorder", False),
        ("louds", "adjlist", False),
        ("louds", "pointer_bfs", False),
    ]
    if scope in ("all", "succinct"):
        result += [
            ("bp", "bp_rmmbtree", True),
            ("dfuds", "dfuds_rmmbtree", True),
            ("louds", "louds_compressed", True),
        ]
    return result


def operations(prefix, succinct, mode="optimized"):
    if succinct and mode == "best":
        return [prefix.split("_")[0] + "_best_" + op for op in
                ("dfs", "dfs_children", "random_path", "random_path_children")]
    if succinct and mode == "unpack-cache":
        family = prefix.split("_")[0]
        return operations(prefix, True, "comparison") + [
            family + "_" + policy + "_" + op
            for policy in ("selective_unpack", "chunk_unpack")
            for op in ("dfs", "dfs_children", "random_path", "random_path_children")]
    if succinct and mode == "comparison":
        return operations(prefix, True, "cache") + operations(prefix, True, "optimized")
    if succinct and mode == "cache":
        family = prefix.split("_")[0]
        return [prefix + "_view_dfs"] + [prefix + "_" + op for op in
            ("dfs_children", "random_path", "random_path_children")] + [
            family + "_clock_" + op for op in
            ("dfs", "dfs_children", "random_path", "random_path_children")]
    if succinct:
        return [
            prefix + "_unpacked_" + op
            for op in ("dfs", "dfs_children", "random_path",
                       "random_path_children")
        ] + [
            prefix + "_indexed_" + op
            for op in ("random_path", "random_path_children")
        ]
    return [
        prefix + "_" + op
        for op in ("dfs", "dfs_children", "random_path",
                   "random_path_children")
    ]


def validate(path, cases, repetitions, cache_bytes=512 * 1024):
    data = json.loads(path.read_text())
    rows = data["benchmarks"]
    if any(row.get("error_occurred") for row in rows):
        raise RuntimeError(f"Benchmark error/skip in {path}")
    actual = {row.get("run_name", row["name"]) for row in rows}
    if actual != set(cases):
        raise RuntimeError(f"Incomplete or unexpected cases in {path}")
    for case in cases:
        samples = [r for r in rows
                   if r.get("run_type") == "iteration"
                   and r.get("run_name", r["name"]) == case]
        if len(samples) != repetitions:
            raise RuntimeError(f"Missing individual repetitions for {case}")
        for sample in samples:
            if sample["iterations"] < 1:
                raise RuntimeError(f"Empty measurement for {case}")
            if "_best_" in case:
                if (sample.get("cache_limit_bytes") != cache_bytes or
                        sample.get("cache_and_scratch_bytes", float("inf")) >
                        sample["cache_limit_bytes"]):
                    raise RuntimeError(f"Invalid combined cache budget for {case}")
    return rows


def export(output, manifest, completed, active=None, error=None):
    rows = []
    for key in completed:
        group = next(g for g in manifest["groups"] if g["key"] == key)
        for row in json.loads((output / f"{key}.json").read_text())["benchmarks"]:
            rows.append(dict(row, family=group["family"],
                             representation=group["representation"],
                             size_exponent=group["exponent"],
                             size_nodes=group["nodes"],
                             source_json=f"{key}.json"))
    atomic_json(output / "combined.json", {"benchmarks": rows})
    columns = sorted({key for row in rows for key in row})
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=columns)
    writer.writeheader()
    writer.writerows(rows)
    temporary = output / "results.csv.tmp"
    temporary.write_text(stream.getvalue())
    temporary.replace(output / "results.csv")
    atomic_json(output / "status.json", {
        "state": "failed" if error else
                 "complete" if len(completed) == len(manifest["groups"]) else
                 "running",
        "completed_groups": len(completed),
        "total_groups": len(manifest["groups"]),
        "completed_cases": sum(len(g["cases"]) for g in manifest["groups"]
                               if g["key"] in completed),
        "total_cases": sum(len(g["cases"]) for g in manifest["groups"]),
        "active_group": active,
        "error": error,
        "updated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=Path("build/benchmarks"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--scope", choices=("all", "references", "succinct"), default="all")
    parser.add_argument("--exponents", type=int, nargs="+",
                        default=list(range(10, 27, 2)))
    parser.add_argument("--nodes", type=int, nargs="+",
                        help="Explicit node counts, replacing --exponents")
    parser.add_argument("--first-size", type=int,
                        help="Measure this selected size before the remaining grid")
    parser.add_argument("--succinct-mode", choices=("best", "optimized", "cache", "comparison", "unpack-cache"),
                        default="best")
    parser.add_argument("--random-interleaving", action="store_true")
    parser.add_argument("--cpu", type=int, default=0)
    parser.add_argument("--cache-bytes", type=int, default=512 * 1024)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--min-time", type=float, default=1.0)
    parser.add_argument("--warmup", type=float, default=0.2)
    parser.add_argument("--timeout", type=float, default=3600)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--list-only", action="store_true")
    args = parser.parse_args()
    if (args.cache_bytes < 0 or args.repetitions < 1 or args.min_time <= 0 or args.warmup < 0
            or args.timeout <= 0 or args.cpu < 0
            or not all(10 <= e <= 26 for e in args.exponents)
            or len(set(args.exponents)) != len(args.exponents)):
        parser.error("Invalid duration, count, CPU, or exponent range")
    sizes = sorted(args.nodes if args.nodes is not None else
                   [1 << e for e in args.exponents])
    if not sizes or len(sizes) != len(set(sizes)) or not all(
            1 << 10 <= n <= 1 << 26 for n in sizes):
        parser.error("Node counts must be unique and in [2^10, 2^26]")
    if args.first_size is not None:
        if args.first_size not in sizes:
            parser.error("--first-size must be a selected node count")
        sizes = [args.first_size] + [n for n in sizes if n != args.first_size]
    build = args.build_dir.resolve()
    output = args.output_dir.resolve()
    config = {
        "cache_bytes": args.cache_bytes,
        "cpu": args.cpu, "repetitions": args.repetitions,
        "min_seconds": args.min_time, "warmup_seconds": args.warmup,
        "tree_seed": 42, "query_seed": 43, "query_pool_bytes": 65536,
        "scope": args.scope, "exponents": sorted(args.exponents),
        "sizes": sizes, "succinct_mode": args.succinct_mode,
        "random_interleaving": args.random_interleaving,
    }
    binaries = {}
    registrations = {}
    for family, _, _ in groups(args.scope):
        if family in binaries:
            continue
        binary = build / f"{family}_tree_benchmarks"
        binaries[family] = {
            "sha256": digest(binary),
            "flags": (build / "CMakeFiles" / f"{family}_tree_benchmarks.dir"
                      / "flags.make").read_text(),
        }
        registrations[family] = subprocess.check_output(
            [str(binary), "--benchmark_list_tests=true"], text=True).splitlines()
    selected = []
    for nodes in sizes:
        exponent = nodes.bit_length() - 1 if nodes & (nodes - 1) == 0 else None
        for family, prefix, succinct in groups(args.scope):
            names = operations(prefix, succinct, args.succinct_mode)
            pattern = "^(" + "|".join(names) + f")/tree_size:{nodes}(/|$)"
            cases = [c for c in registrations[family] if re.search(pattern, c)]
            if len(cases) != len(names):
                raise RuntimeError(f"Missing cases for {prefix}, {nodes} nodes")
            # Registered MinTime would silently override --benchmark_min_time.
            if any("/min_time:" in c for c in cases):
                raise RuntimeError("Rebuild: registered MinTime overrides run settings")
            selected.append({
                "key": (f"n{exponent:02d}" if exponent is not None else
                        f"nodes{nodes}") + f"-{prefix}", "family": family,
                "representation": prefix, "exponent": exponent,
                "nodes": nodes, "cases": cases, "filter": pattern,
            })
    identity = {"config": config, "binaries": binaries, "groups": selected}
    print(f"{len(selected)} groups, {sum(len(g['cases']) for g in selected)} "
          f"cases, {args.repetitions} repetitions", flush=True)
    if args.list_only:
        return
    output.mkdir(parents=True, exist_ok=True)
    with Path("/tmp/pixie-tree-benchmarks.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        manifest_path = output / "manifest.json"
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text())
            if not args.resume or any(manifest[k] != v
                                      for k, v in identity.items()):
                raise RuntimeError("Existing run differs; use a fresh output directory")
        else:
            manifest = dict(identity)
            manifest["created_utc"] = time.strftime(
                "%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            manifest["host"] = subprocess.check_output(["lscpu"], text=True)
            manifest["kernel"] = subprocess.check_output(["uname", "-a"], text=True)
            manifest["compiler"] = subprocess.check_output(
                ["c++", "--version"], text=True)
            manifest["revision"] = subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True).strip()
            manifest["source_hashes"] = {
                str(p): digest(p)
                for base in ("include", "src/benchmarks")
                for p in Path(base).rglob("*")
                if p.suffix in (".h", ".cpp")
            }
            (output / "bin").mkdir(exist_ok=True)
            for family in binaries:
                shutil.copy2(build / f"{family}_tree_benchmarks",
                             output / "bin" / f"{family}_tree_benchmarks")
            shutil.copy2(__file__, output / "collector.py")
            with zipfile.ZipFile(output / "sources.zip", "w",
                                 zipfile.ZIP_DEFLATED) as archive:
                for path in manifest["source_hashes"]:
                    archive.write(path)
            (output / "README.md").write_text(
                "# Tree benchmark sweep\n\n"
                "See `status.json` for completion and the active group. "
                "`manifest.json` records case selection, settings, source hashes, "
                "compiler flags, CPU, and executable hashes.\n\n"
                "Each implementation/size group has raw Google Benchmark JSON "
                "(individual repetitions and aggregates), a console log, exact "
                "command, process resource measurements, and a completion marker. "
                "`combined.json` and `results.csv` contain only completed, "
                "validated groups and are refreshed after every group.\n\n"
                "Time columns retain their `time_unit`; `cpu_time` and "
                "`real_time` are separate. DFS timings are per full traversal; "
                "path timings are per root-to-leaf walk. Construction and bounded "
                "index construction are excluded from query timing. Resource "
                "logs include construction and are not structure-memory counts. "
                "Do not mix aggregate rows with individual repetitions when "
                "computing statistics (`run_type`, `aggregate_name`).\n\n"
                "Trees use seed 42; path pools use seed 43 and shifted offsets. "
                "Inputs are fixed across repetitions for cross-implementation "
                "comparison. Unlike rank/select's default repeated input "
                "generation, this sweep does not vary tree seeds. Adaptive "
                "iteration counts may differ; query counters record work.\n\n"
                "Succinct variants use compressed rank/select and RmMBTree for "
                "BP/DFUDS. DFS uses unpacking. Paths include both unpacked and "
                "bounded-index variants; index bytes are explicit counters. "
                "References include preorder and BFS layouts. No binary RmM "
                "or baseline LOUDS backend comparisons are included.\n\n"
                "Run the collector again with identical arguments plus "
                "`--resume` to reuse completed groups. It rejects changed "
                "executables/settings and validates result hashes. `bin/`, "
                "`sources.zip`, and `collector.py` preserve this run's inputs. "
                "Keep this user-requested dataset out of commits.\n")
            if args.succinct_mode in ("cache", "comparison", "unpack-cache"):
                (output / "README.md").write_text(
                    "# Cache versus references benchmark collection\n\n"
                    "All four workloads compare local-view memoization and "
                    "the experimental CLOCK policy. References use vectors "
                    "and three pointers in both preorder and BFS layouts. "
                    "Comparison mode also retains separately named unpacked "
                    "and fixed-root-index controls.\n\n"
                    "status.json records completion; manifest.json preserves "
                    "settings, commands, source hashes, and binary identity. "
                    "Raw repetitions and aggregates are retained separately "
                    "in each result; use only iteration rows for statistics. "
                    "Process resource logs include setup, which is excluded "
                    "from query timing. The dataset is user-requested and "
                    "should remain outside commits.\n")
            if args.succinct_mode == "best":
                (output / "README.md").write_text(
                    "# Fixed tree configurations: four workloads\n\n"
                    "One BP, DFUDS and LOUDS configuration at every size, with "
                    f"a {args.cache_bytes}-byte combined cache and traversal-buffer cap. "
                    "Constexpr tables are excluded. References, when selected, include vectors "
                    "and three pointers in BFS and preorder layouts.\n\n"
                    "See ../README.md for configuration and budget details. "
                    "status.json records progress; manifest.json records "
                    "settings, source and binary hashes. Individual repetitions "
                    "are preserved in combined.json and results.csv; setup is "
                    "excluded from query timing.\n")
            atomic_json(manifest_path, manifest)
        for family, details in binaries.items():
            if digest(output / "bin" / f"{family}_tree_benchmarks") != details["sha256"]:
                raise RuntimeError("Preserved executable hash mismatch")
        completed = []
        active = None
        try:
            for group in selected:
                key = group["key"]
                result = output / f"{key}.json"
                marker = output / f"{key}.complete.json"
                active = key
                if args.resume and marker.exists():
                    validate(result, group["cases"], args.repetitions, args.cache_bytes)
                    if json.loads(marker.read_text())["sha256"] != digest(result):
                        raise RuntimeError(f"Result hash mismatch: {result}")
                    completed.append(key)
                    continue
                export(output, manifest, completed, active=active)
                command = [
                    "/usr/bin/time", "-v", "-o", str(output / f"{key}.resources.txt"),
                    "taskset", "-c", str(args.cpu),
                    str(output / "bin" / f"{group['family']}_tree_benchmarks"),
                    "--benchmark_filter=" + group["filter"],
                    f"--benchmark_repetitions={args.repetitions}",
                    f"--benchmark_min_time={args.min_time}s",
                    f"--benchmark_min_warmup_time={args.warmup}",
                    "--benchmark_enable_random_interleaving=" +
                    str(args.random_interleaving).lower(),
                    "--benchmark_report_aggregates_only=false",
                    "--benchmark_display_aggregates_only=true",
                    "--benchmark_out_format=json", f"--benchmark_out={result}",
                ]
                atomic_json(output / f"{key}.command.json", command)
                print(f"[{len(completed) + 1}/{len(selected)}] {key}", flush=True)
                with (output / f"{key}.log").open("w") as log:
                    with subprocess.Popen(command, stdout=log, stderr=log,
                                          env=dict(os.environ, PIXIE_TREE_CACHE_BYTES=str(args.cache_bytes)),
                                          start_new_session=True) as process:
                        try:
                            code = process.wait(timeout=args.timeout)
                        except BaseException:
                            os.killpg(process.pid, signal.SIGKILL)
                            process.wait()
                            raise
                        if code:
                            raise RuntimeError(f"{key}: exit {code}; inspect its log")
                validate(result, group["cases"], args.repetitions, args.cache_bytes)
                atomic_json(marker, {"sha256": digest(result)})
                completed.append(key)
                export(output, manifest, completed)
                print(f"  saved {len(group['cases'])} cases", flush=True)
            export(output, manifest, completed)
        except BaseException as error:
            export(output, manifest, completed, active, str(error) or type(error).__name__)
            raise


if __name__ == "__main__":
    main()
