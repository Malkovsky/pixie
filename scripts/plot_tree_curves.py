#!/usr/bin/env python3
"""Render the three offline tree benchmark charts from finished timing JSONs."""
import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import re
from statistics import mean

PREVIEW = Path(__file__).resolve().parents[1] / 'src/docs/images/benchmarks/trees/preview'
FAMILIES = {'random_recursive', 'binary_catalan', 'ordered_plane_catalan'}
CHARTS = {'inspected_dfs': ('Inspected DFS', 'inspected_dfs'),
          'inspected_root_path': ('Inspected root paths', 'inspected_root_path'),
          'random': ('Random endpoint paths', 'tree_path')}
NAME = re.compile(r'(tree_shape|tree_path)/(bp|dfuds|louds|forest_cover|adj_list|pointer)/'
                  r'(inspected_dfs|inspected_root_path|random)/(\w+)/N:(\d+)/seed:42(?:/budget:(\d+))?')


def collect(paths, random):
    groups = defaultdict(dict)
    for path in paths:
        for row in json.loads(path.read_text())['benchmarks']:
            if row.get('error_occurred') or row.get('skipped'):
                raise ValueError(f'Benchmark error in {path}: {row}')
            if row.get('aggregate_name') != 'median':
                continue
            m = NAME.fullmatch(row['run_name'])
            if not m:
                raise ValueError(f'Unexpected timing in {path}: {row["run_name"]}')
            suite, impl, op, family, n, b = m.groups()
            n, b = int(n), int(b or 0)
            if ((op == 'random') != random or (suite == 'tree_path') != random or
                    family not in FAMILIES or b not in (0, 32768, 524288) or
                    row['budget_limit_bytes'] != b):
                raise ValueError(f'Unexpected timing: {row["run_name"]}')
            t = row['cpu_time'] * {'ns': 1, 'us': 1000, 'ms': 10**6, 's': 10**9}[row['time_unit']]
            if op == 'inspected_dfs':
                if row['visited_nodes_per_walk'] != n:
                    raise ValueError(f'Incomplete DFS: {row["run_name"]}')
                t /= n
            elif op == 'inspected_root_path':
                t /= row['paths_per_iteration']
            # Random endpoint timing is already per query, not per pair pool.
            length = None if op == 'inspected_dfs' else row['nodes_per_path']
            g = groups[op, impl, b, n]
            if family in g or not (math.isfinite(t) and t > 0):
                raise ValueError(f'Duplicate or invalid median in {path}: {row["run_name"]}')
            g[family] = t, length
    charts = defaultdict(list)
    for (op, impl, b, n), g in sorted(groups.items()):
        if g.keys() != FAMILIES:
            if not random:
                raise ValueError(f'Incomplete inspected group: {op}/{impl}/{b}/{n}')
            print(f'Omitted {op}/{impl}/{b}/{n}: missing {sorted(FAMILIES - g.keys())}')
            continue
        point = dict(impl=impl, b=b, n=n, t=mean(t for t, _ in g.values()))
        if op != 'inspected_dfs':
            point['length'] = mean(length for _, length in g.values())
        charts[op].append(point)
    return charts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspected-dir', type=Path, required=True)
    parser.add_argument('--random-initial', type=Path, required=True)
    parser.add_argument('--random-shards', type=Path, help='Directory of finished extension shards')
    args = parser.parse_args()
    charts = collect(sorted(args.inspected_dir.glob('*-timing-*.json')), False)
    paths = [args.random_initial]
    if args.random_shards:
        # A finished size summary contains its component shards. Choose one
        # representation per size from a single directory snapshot.
        sizes = defaultdict(list)
        for path in sorted(args.random_shards.iterdir()):
            match = re.fullmatch(r'(p\d+)(?:-\w+-\w+)?\.json', path.name)
            if match:
                sizes[match[1]].append(path)
        for size, shards in sizes.items():
            summary = args.random_shards / f'{size}.json'
            paths.extend([summary] if summary in shards else shards)
    charts.update(collect(paths, True))
    for op, (title, stem) in CHARTS.items():
        data = charts[op]
        if not data:
            raise ValueError(f'No complete medians for {op}')
        counts = Counter(d['n'] for d in data)
        if op != 'random' and counts != {2**e: 12 for e in range(10, 27, 2)}:
            raise ValueError(f'Incomplete inspected timing matrix: {counts}')
        static = 'tree_path_random_budgets' if op == 'random' else stem
        html = f'''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title} | Pixie</title><link rel="stylesheet" href="tree-curves.css">
<main><h1>{title}</h1><p>Equal-weight mean of the random recursive, binary Catalan and ordered plane Catalan tree-family medians. Only sizes with all three medians have dots. Color identifies the implementation. Solid lines are baselines or references, dashed lines use 32 KiB, and dotted lines use 512 KiB. ForestCover's baseline has one decoder slot.</p>
<details id="picker"><summary id="trigger">Visible curves</summary><div id="panel">
<input id="search" type="search" placeholder="Search curves&hellip;" aria-label="Search curves">
<fieldset id="controls"><legend>Visible curves</legend></fieldset>
<button id="all">Show all</button><button id="none">Clear selection</button><button id="done">Done</button></div></details>
<div id="chart"></div><noscript><a href="{static}.svg">View static preview</a></noscript></main>
<script id="points" type="application/json" data-operation="{op}">{json.dumps(data, allow_nan=False, separators=(',', ':'))}</script>
<script src="d3-7.9.0.min.js"></script><script src="plot-0.6.17.min.js"></script><script src="tree-curves.js"></script></html>
'''
        (PREVIEW / f'{stem}_observable.html').write_text(html)
        print(f'{stem}: {len(data)} points; curves by node count: {dict(sorted(counts.items()))}')


if __name__ == '__main__':
    main()
