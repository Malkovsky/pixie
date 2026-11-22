#!/usr/bin/env python3
"""Export the matched tree chart measurements from finished Google Benchmark JSON."""
import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
from statistics import mean

from report_tree_shapes import FAMILIES, NS

OUTPUT = Path(__file__).resolve().parents[1] / 'src/docs/benchmarks/data/tree-matched.js'
SIZES = tuple(2**i for i in range(10, 27, 2))
BUDGETS = (0, 32768, 524288)
SUCCINCT = ('bp', 'dfuds', 'louds', 'forest_cover')
REFERENCES = ('adj_list', 'pointer', 'bfs_csr')
ENDPOINT_NS = 'CPU ns per query as reported (ns/path); pairs_per_batch is not a divisor'
NORMALIZATION = {'inspected_dfs': 'CPU ns / visited_nodes_per_walk (ns/node)',
                 'inspected_root_path': 'CPU ns / paths_per_iteration (ns/path)',
                  'tree_path': ENDPOINT_NS}
NAME = re.compile(r'(tree_shape|tree_path)/(\w+)/(inspected_dfs|inspected_root_path|random|sibling)/'
                  r'(\w+)/N:(\d+)/seed:42(?:/budget:(\d+))?')


def export(inspected, random, references, csr):
    groups = defaultdict(dict)
    workloads = {}
    sources = []
    for suite, role, paths in [('tree_shape', 'inspected', inspected), ('tree_shape', 'references', references or ()),
                               ('tree_path', 'random', random), ('tree_path', 'csr', [csr])]:
        for path in paths:
            raw = path.read_bytes()
            doc = json.loads(raw)
            context = doc['context']
            nsizes = set()
            for row in doc['benchmarks']:
                if any(row.get(k) for k in ('error_occurred', 'skipped', 'error_message', 'skip_message')):
                    raise ValueError(f'Failed benchmark in {path}: {row}')
                if row.get('aggregate_name') != 'median':
                    continue
                name = row['run_name']
                m = NAME.fullmatch(name)
                if not m or m[1] != suite or (m[3] in ('random', 'sibling')) != (suite == 'tree_path'):
                    raise ValueError(f'Unexpected workload in {path}: {name}')
                _, impl, op, family, n, b = m.groups()
                if ((role == 'references' and impl not in REFERENCES) or
                        (role == 'inspected' and impl not in SUCCINCT) or
                        (role == 'csr' and impl != 'bfs_csr') or
                        (role == 'random' and (op != 'random' or impl == 'bfs_csr'))):
                    raise ValueError(f'Unexpected input for {role}: {name}')
                if role == 'csr' and op == 'sibling':
                    continue  # The CSR file also contains sibling queries, which are not displayed.
                n, b = int(n), int(b or 0)
                nsizes.add(n)
                if op == 'random':
                    op = 'tree_path'
                usage = row['budget_usage_bytes']
                if (row['budget_limit_bytes'] != b or not math.isfinite(usage) or usage < 0 or
                        usage != int(usage) or (b and usage > b)):
                    raise ValueError(f'Invalid budget: {name}')
                usage = int(usage)
                if (n not in SIZES or family not in FAMILIES or row['aggregate_unit'] != 'time' or
                        row['run_type'] != 'aggregate' or row['threads'] != 1):
                    raise ValueError(f'Invalid size, family or timing row: {name}')
                if suite == 'tree_shape' and row['nodes'] != n:
                    raise ValueError(f'Incorrect tree size: {name}')
                ns = row['cpu_time'] * NS[row['time_unit']]
                divisor = 1
                if op == 'inspected_dfs':
                    divisor = row['visited_nodes_per_walk']
                    if divisor != n or row['inspected_children_per_walk'] != n - 1:
                        raise ValueError(f'Incomplete walk: {name}')
                elif op == 'inspected_root_path':
                    divisor = row['paths_per_iteration']
                    if divisor != 4096:
                        raise ValueError(f'Unexpected root-path batch: {name}')
                elif row['pairs_per_batch'] != 256:
                    raise ValueError(f'Unexpected endpoint pool: {name}')
                sample = [ns / divisor, int(divisor), usage]
                if op != 'inspected_dfs':
                    length = row['nodes_per_path']
                    if not 1 <= length <= n:
                        raise ValueError(f'Invalid path length: {name}')
                    sample.append(length)
                    counts = (length, row['inspected_children_per_path'] if op == 'inspected_root_path'
                              else row['max_nodes_per_path'])
                else:
                    counts = (row['return_events_per_walk'],)
                if workloads.setdefault((op, family, n), counts) != counts:
                    raise ValueError(f'Different traversal workload: {name}')
                g = groups[op, impl, b, n]
                if family in g or row['repetitions'] != 5 or not math.isfinite(ns) or ns <= 0:
                    raise ValueError(f'Duplicate or invalid family median: {name}')
                g[family] = sample
            sources.append(dict(file=path.name, suite=suite, role=role, sizes=sorted(nsizes), date=context['date'],
                                executable=Path(context['executable']).name, library=context['library_version'],
                                sha256=hashlib.sha256(raw).hexdigest()))

    charts = {}
    for op in NORMALIZATION:
        keys = [(impl, b, n) for impl in SUCCINCT for b in BUDGETS for n in SIZES]
        if references is not None or op == 'tree_path':
            keys.extend((impl, 0, n) for impl in REFERENCES for n in SIZES)
        charts[op] = []
        for impl, b, n in keys:
            group = groups.pop((op, impl, b, n), {})
            if group.keys() != set(FAMILIES):
                raise ValueError(f'Missing family medians for {op}/{impl}/{b}/{n}: {set(FAMILIES) - group.keys()}')
            samples = [group[f] for f in FAMILIES]
            point = dict(impl=impl, b=b, n=n, t=mean(s[0] for s in samples), families=samples)
            if op != 'inspected_dfs':
                point['length'] = mean(s[3] for s in samples)
            charts[op].append(point)
    if groups:
        raise ValueError(f'Unexpected benchmark groups: {list(groups)}')
    metadata = dict(seed=42, repetitions=5, pairs_per_batch=256, families=FAMILIES,
                    family_columns=['normalized_cpu_ns', 'divisor', 'budget_usage_bytes', 'nodes_per_path (paths only)'],
                    normalization=NORMALIZATION, aggregation='Equal-weight arithmetic mean of three family medians.',
                    budgets='b is the traversal budget in bytes; 0 denotes baseline, including one ForestCover slot.',
                    sources=sources)
    return dict(metadata=metadata, charts=charts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspected-dir', type=Path, required=True)
    inspected = parser.add_mutually_exclusive_group(required=True)
    inspected.add_argument('--inspected-refs', type=Path, nargs='+', help='Finished reference timing shards')
    inspected.add_argument('--succinct-inspected-only', action='store_true', help='Use the succinct-only inspected matrix')
    parser.add_argument('--csr-endpoints', type=Path, required=True)
    parser.add_argument('--random-initial', type=Path, required=True)
    parser.add_argument('--random-shards', type=Path, help='Directory of finished extension shards')
    args = parser.parse_args()
    paths = [args.random_initial]
    if args.random_shards:
        sizes = defaultdict(list)
        for path in sorted(args.random_shards.iterdir()):
            match = re.fullmatch(r'(p\d+)(?:-\w+-\w+)?\.json', path.name)
            if match:
                sizes[match[1]].append(path)
        for size, shards in sizes.items():
            summary = args.random_shards / f'{size}.json'
            paths.extend([summary] if summary in shards else shards)
    result = export(sorted(args.inspected_dir.glob('*-timing-*.json')), paths, args.inspected_refs,
                     args.csr_endpoints)
    # One point per line keeps diffs reviewable without pretty-printing each number.
    def compact(value):
        return json.dumps(value, allow_nan=False, separators=(',', ':'))
    charts = ',\n'.join(f'"{op}":[\n' + ',\n'.join(map(compact, points)) + '\n]'
                        for op, points in result['charts'].items())
    script = 'globalThis.PIXIE_TREE_MATCHED = {\n"metadata":' + compact(result['metadata'])
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(script + ',\n"charts":{\n' + charts + '\n}};\n', encoding='utf-8')
    for op, points in result['charts'].items():
        print(f'{op}: {len(points)} points at {len({p["n"] for p in points})} sizes')


if __name__ == '__main__':
    main()
