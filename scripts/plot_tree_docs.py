#!/usr/bin/env python3
"""Export tree figures from local measurement archives without rerunning benchmarks."""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = 'build/tree-refresh-8e72c47'
REVISION = '8e72c4733846eb2d7ae488420f6375c72f1e7735'
SIZES = [2**e for e in range(10, 27, 2)]
PANELS = [
    ('dfs', 'DFS', 'ns/node'),
    ('dfs_children', 'DFS with child inspection', 'ns/node'),
    ('random_path', 'Root-to-leaf path', 'ns/path'),
    ('random_path_children', 'Root-to-leaf path with child inspection', 'ns/path'),
    ('build_prepared', 'Construction from prepared input', 'ms/build'),
    ('memory', 'Memory including reserved traversal scratch', 'bits/node'),
]
# (Google Benchmark prefix, probe tree, figure legend, renderer family, curve, bytes)
SERIES = [
    ('adjlist_preorder', 'vectors_preorder', 'Vectors', 'adj_list', 'preorder', 0),
    ('adjlist', 'vectors_bfs', 'Vectors', 'adj_list', 'BFS', 0),
    ('pointer_preorder', 'pointers_preorder', 'Pointers', 'pointer', 'preorder', 0),
    ('pointer_bfs', 'pointers_bfs', 'Pointers', 'pointer', 'BFS', 0),
] + [(tree + '_best', tree, tree.upper(), tree, f'{b//1024} KiB', b)
     for tree in ('bp', 'dfuds', 'louds') for b in (32768, 524288)]


def check_summary(path, digest):
    # The CSVs have no date/revision column. Pin each to its published vintage.
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError(f'Unexpected measurement summary: {path}')


def prepared_summary(path, digest, cases):
    check_summary(path, digest)
    points = {}
    with path.open(newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            tree, b, n = row['tree'], int(row['budget_bytes']), int(row['nodes'])
            base, charge, total = (int(row[k]) for k in ('tree_bytes', 'charged_bytes', 'working_bytes'))
            ms = float(row['prepared_ms'])
            key = tree, b, n
            if ((tree, b) not in cases or n not in SIZES or key in points or int(row['repetitions']) != 3
                    or not math.isfinite(ms) or ms <= 0 or base <= 0 or not 0 <= charge <= b
                    or total != base + charge or float(row['bytes_per_node']) != total / n):
                raise ValueError(f'Invalid prepared-input row: {row}')
            points[key] = {'build_prepared': ms, 'memory': total / n}
    if points.keys() != {(tree, b, n) for tree, b in cases for n in SIZES}:
        raise ValueError(f'Expected nine measured sizes per curve in {path}')
    return points


def export(archive, forest_summary, uncached_summary, bfs_csr_summary):
    samples = defaultdict(dict)

    def add(panel, tree, budget, n, rep, value):
        key = panel, tree, budget, n
        if rep in samples[key]:
            raise ValueError(f'Duplicate repetition {rep} for {key}')
        samples[key][rep] = value

    for kib in (32, 512):
        run = archive / str(kib)
        if json.loads((run / 'manifest.json').read_text())['revision'] != REVISION:
            raise ValueError(f'Wrong revision in {run}')
        # The 512 KiB run measured the references too. They have no cache budget.
        cases = {f'{prefix}_{op}': (op, tree, b) for prefix, tree, _, _, _, b in SERIES
                 if b == kib*1024 or (kib == 512 and b == 0) for op, _, _ in PANELS[:4]}
        for row in json.loads((run / 'combined.json').read_text())['benchmarks']:
            if row.get('error_occurred'):
                raise ValueError(f'Benchmark error: {row["name"]}')
            if row['run_type'] != 'iteration':
                continue
            op, tree, b = cases[row['run_name'].split('/')[0]]
            if b and (row['cache_limit_bytes'] != b or row['cache_and_scratch_bytes'] > b):
                raise ValueError(f'Invalid traversal cache budget: {row["name"]}')
            n = int(row['size_nodes'])
            ns = row['cpu_time'] * {'ns': 1, 'us': 1000, 'ms': 10**6, 's': 10**9}[row['time_unit']]
            add(op, tree, b, n, row['repetition_index'], ns/n if op.startswith('dfs') else ns)

    for row in json.loads((archive / 'construction/results.json').read_text()):
        tree, b, n, rep = (row[k] for k in ('tree', 'budget', 'nodes', 'repetition'))
        if b and row['persistent_cache_bytes'] + row['scratch_bytes'] > b:
            raise ValueError(f'Invalid construction cache budget: {tree}/{n}')
        add('build_prepared', tree, b, n, rep, (row['construct_s'] + row['cache_s']) * 1000)
        add('memory', tree, b, n, rep, row['working_bytes'] / n)

    # Check each reconstructed median against the exact table saved with the SVGs.
    plotted = {}
    with (archive / 'graphs/plot_data.csv').open(newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            op = row['operation']
            if op == 'build':
                if row['panel_row'] == '1':  # Full setup is archived, not published.
                    continue
                if row['panel_row'] != '0':
                    raise ValueError(f'Unknown build panel: {row["panel_row"]}')
                op = 'build_prepared'
            key = op, row['series'], int(row['budget_bytes']), int(row['nodes'])
            if key in plotted:
                raise ValueError(f'Duplicate plotted point: {key}')
            plotted[key] = float(row['median']), int(row['repetitions'])

    def js(obj):
        return json.dumps(obj, separators=(',', ':'), ensure_ascii=False, allow_nan=False)

    budgeted = [('forest_cover', 32768), ('forest_cover', 524288)]
    uncached = [(tree, 0) for tree in ('bp', 'dfuds', 'louds', 'forest_cover')]
    prepared = prepared_summary(forest_summary,
        'a839ad389e2d90ca5d2c3fe0aa30154a7c7e822c9f08e097bc8731bd3f3bff48', budgeted)
    prepared.update(prepared_summary(uncached_summary,
        'e3796c932bd527359340d47d9e6c7b1c79536ad938d794a2853aad3faf19d954', uncached))

    check_summary(bfs_csr_summary, '556ba66906e2f3b47a2889534d49a5674055eab9b6297217c362a37b42f18e54')
    with bfs_csr_summary.open(newline='', encoding='utf-8') as f:
        csr = list(csv.DictReader(f))
    if [int(row['nodes']) for row in csr] != SIZES:
        raise ValueError('Expected nine ordered BFS CSR measurements')

    metadata = [dict(impl=impl, family=f'{family} (2026-10-30, 8e72c47)', curve=curve, b=b)
                for _, _, family, impl, curve, b in SERIES]
    metadata.append(dict(impl='bfs_csr', family='BFS CSR (2026-11-20, 6ba63bb working tree)',
                         curve='Packed', b=0))
    prepared_cases = budgeted + uncached
    for tree, b in prepared_cases:
        name = 'ForestCover' if tree == 'forest_cover' else tree.upper()
        curve = f'{b//1024} KiB' if b else 'Uncached'
        if not b and tree != 'forest_cover':
            # The picker groups by implementation, under the older run's heading.
            curve += ' (2026-11-20, 6ba63bb)'
        metadata.append(dict(impl=tree, family=f'{name} (2026-11-20, 6ba63bb)',
                             curve=curve, b=b))
    lines = [f'// Generated by scripts/plot_tree_docs.py. Historical archive: {ARCHIVE} at {REVISION}.',
              '// 32/combined.json, 512/combined.json (including references), construction/results.json.',
              '// BFS CSR: bfs-csr-publication/classic-prepared.csv, 2026-11-20 working tree based on 6ba63bb.',
              '// Classic traversal: median CPU time (5 repetitions); prepared build: median wall time (3); storage: working bytes, no cache.',
             '// ForestCover: 2026-11-20 prepared-input probe on 6ba63bb, summary.csv (3 repetitions).',
             '// Uncached: 2026-11-20 probe on 6ba63bb; consume prepared encoding for BP/DFUDS/LOUDS, ordered source for ForestCover (3 repetitions).',
             '// Random recursive tree, seed 42. Each row is one curve in ascending node count.',
             '(() => {', 'const charts = window.pixieTreeCharts ??= {};',
             f'const sizes = {js(SIZES)};',
             'const series = ' + js(metadata) + ';',
             'function chart(key, title, unit, yLabel, yScale, values) {',
             '  charts[key] = {title, unit, yLabel, yScale, points:values.flatMap((row, i) =>',
             '    row.map((t, j) => ({...series[i], n:sizes[j], t})))};', '}']
    for panel, title, unit in PANELS:
        traversal = unit.startswith('ns/')
        reps = 5 if traversal else 3
        metric = 'Median CPU time' if traversal else 'Median wall time' if unit == 'ms/build' else 'Memory'
        key = 'history-' + panel.replace('_', '-')
        args = [key, title, unit, f'{metric} ({unit})', 'linear' if traversal else 'log']
        lines.append('chart(' + ', '.join(map(js, args)) + ', [')
        def append_curve(values):
            # Verify archive medians in their original units; convert only for display.
            if panel == 'memory':
                values = [8 * bytes_per_node for bytes_per_node in values]
            lines.append('  ' + js(values) + ',')

        for _, tree, label, _, curve, b in SERIES:
            values = []
            for n in SIZES:
                iterations = samples.pop((panel, tree, b, n))
                value = median(iterations.values())
                key = panel, f'{label} · {curve}', b, n
                if (set(iterations) != set(range(reps)) or plotted.pop(key) != (value, reps)):
                    raise ValueError(f'Repetitions or median disagree for {key}')
                values.append(value)
            append_curve(values)
        column = {'build_prepared': 'prepared_ms', 'memory': 'bytes_per_node'}.get(panel, panel)
        append_curve([float(row[column]) for row in csr])
        if panel in ('build_prepared', 'memory'):
            for tree, b in prepared_cases:
                append_curve([prepared[tree, b, n][panel] for n in SIZES])
        lines.append(']);')
    if samples or plotted:
        raise ValueError('Unexpected extra measurements')
    return '\n'.join(lines + ['})();', ''])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=ROOT / ARCHIVE)
    parser.add_argument('--forest-summary', type=Path, required=True,
                        help='ForestCover summary.csv from the 2026-11-20 prepared-input run')
    parser.add_argument('--uncached-summary', type=Path, required=True,
                        help='summary.csv from the four-tree uncached prepared-input run')
    parser.add_argument('--bfs-csr-summary', type=Path, required=True,
                        help='classic-prepared.csv from the 2026-11-20 BFS CSR publication run')
    parser.add_argument('--output', type=Path, default=ROOT / 'src/docs/benchmarks/data/tree-history.js')
    args = parser.parse_args()
    data = export(args.input_dir, args.forest_summary, args.uncached_summary, args.bfs_csr_summary)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(data, encoding='utf-8')
