#!/usr/bin/env python3
"""Fresh-process, pinned construction measurements; save each result immediately."""
import argparse
import csv
import hashlib
import json
import math
import random
import statistics
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--output-dir', type=Path, required=True)
parser.add_argument('--exponents', type=int, nargs='+', default=[20,22,24,26])
parser.add_argument('--resume', action='store_true')
succinct = ('bp', 'dfuds', 'louds', 'forest_cover')
trees = (*succinct, 'vectors_bfs', 'vectors_preorder', 'pointers_bfs', 'pointers_preorder', 'bfs_csr')
parser.add_argument('--trees', choices=trees, nargs='+', default=trees)
parser.add_argument('--succinct-budgets', type=int, choices=(0, 32768, 524288),
                    nargs='+', default=(32768, 524288), help='0 omits all optional navigation state')
parser.add_argument('--construction-mode', choices=['consume', 'copy'], default='consume',
                    help='Transfer prepared aligned encoding or copy prepared words')
args = parser.parse_args()
if any(e < 10 or e > 26 for e in args.exponents) or len(set(args.exponents)) != len(args.exponents):
    parser.error('Use distinct exponents between 10 and 26')
root = args.output_dir.resolve()
root.mkdir(parents=True, exist_ok=True)
binary = args.binary.resolve()
previous = json.loads((root / 'manifest.json').read_text()) if args.resume else None
rows = json.loads((root / 'results.json').read_text()) if args.resume else []
digest = hashlib.sha256(binary.read_bytes()).hexdigest()
if previous and (previous['binary_sha256'] != digest or
                 previous.get('construction_mode') != args.construction_mode):
    raise ValueError('Resume requires the same measured binary')
if not args.resume and (root / 'results.json').exists():
    raise ValueError('Output exists; use --resume or a fresh directory')
done = {(r['tree'],r['nodes'],r['budget'],r['repetition']) for r in rows}
jobs = []
for exponent in args.exponents:
    for name in dict.fromkeys(args.trees):
        for budget in (dict.fromkeys(args.succinct_budgets) if name in succinct else (0,)):
            for repetition in range(3):
                job = (name, 2**exponent, budget, repetition)
                if job not in done:
                    jobs.append(job)
random.Random(43).shuffle(jobs)
manifest_jobs = list(dict.fromkeys(tuple(j) for j in
                                 (previous['jobs'] if previous else []) + jobs))
(root / "manifest.json").write_text(json.dumps({
    "jobs": manifest_jobs, "cpu": 0, "seed": 42, "compiler_flags": "-O3 -DNDEBUG -march=native -std=c++20",
    "binary_sha256": digest, "construction_mode": args.construction_mode,
    "clock": "steady_clock elapsed", "repetitions": 3,
}, indent=2))
for i, (name, n, budget, repetition) in enumerate(jobs):
    result = subprocess.run(["taskset", "-c", "0", str(binary), name, str(n), str(budget), args.construction_mode],
                            check=True, capture_output=True, text=True)
    row = json.loads(result.stdout)
    mode = {'forest_cover': 'ordered_source', 'bfs_csr': 'adjacency'}.get(name, args.construction_mode)
    if (row['tree'], row['nodes'], row['budget'], row['construction_mode']) != (name, n, budget, mode):
        raise ValueError(f'Probe returned a different case: {row}')
    row["repetition"] = repetition
    rows.append(row)
    (root / "results.json").write_text(json.dumps(rows, indent=2))
    with (root / "results.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[{i+1}/{len(jobs)}] {name} n={n} budget={budget} total={row['total_s']:.3f}s", flush=True)

# The historical figures use the median of each process's combined time,
# not the sum of separately aggregated tree and cache build times.
groups = {}
for row in rows:
    key = row['tree'], row['budget'], row['nodes']
    groups.setdefault(key, []).append(row)
if sorted(tuple(j) for j in manifest_jobs) != sorted(
        (r['tree'], r['nodes'], r['budget'], r['repetition']) for r in rows):
    raise ValueError('Incomplete or duplicated measurement grid')
summary = []
for (tree, budget, n), group in sorted(groups.items()):
    if sorted(r['repetition'] for r in group) != [0, 1, 2]:
        raise ValueError(f'Incomplete repetitions: {tree}/{budget}/{n}')
    storage = {(r['structure_bytes'] + r['source_bytes'],
                r['persistent_cache_bytes'] + r['scratch_bytes'], r['working_bytes']) for r in group}
    if len(storage) != 1:
        raise ValueError(f'Nonrepeatable storage: {tree}/{budget}/{n}')
    base, charge, total = storage.pop()
    if base <= 0 or total != base + charge or charge > budget:
        raise ValueError(f'Invalid byte accounting: {tree}/{budget}/{n}')
    times = [(r['construct_s'] + r['cache_s']) * 1000 for r in group]
    if not all(math.isfinite(t) and t > 0 for t in times):
        raise ValueError(f'Invalid prepared timing: {tree}/{budget}/{n}')
    ms = statistics.median(times)
    summary.append((tree, budget, n, len(group), ms, base, charge, total, total / n))
with (root / 'summary.csv').open('w', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(('tree', 'budget_bytes', 'nodes', 'repetitions', 'prepared_ms', 'tree_bytes',
                     'charged_bytes', 'working_bytes', 'bytes_per_node'))
    writer.writerows(summary)
