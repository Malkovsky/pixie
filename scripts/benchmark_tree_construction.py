#!/usr/bin/env python3
"""Fresh-process, pinned construction measurements; save each result immediately."""
import argparse
import csv
import hashlib
import json
import random
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--output-dir', type=Path, required=True)
parser.add_argument('--exponents', type=int, nargs='+', default=[20,22,24,26])
parser.add_argument('--resume', action='store_true')
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
    for name in ("bp", "dfuds", "louds", "vectors_bfs", "vectors_preorder",
                 "pointers_bfs", "pointers_preorder"):
        for budget in ((32768, 524288) if name in ("bp", "dfuds", "louds") else (0,)):
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
    row["repetition"] = repetition
    rows.append(row)
    (root / "results.json").write_text(json.dumps(rows, indent=2))
    with (root / "results.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[{i+1}/{len(jobs)}] {name} n={n} budget={budget} total={row['total_s']:.3f}s", flush=True)
