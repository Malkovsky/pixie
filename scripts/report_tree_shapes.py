#!/usr/bin/env python3
"""Report matched tree_shape CPU latencies from timed Google Benchmark JSON."""
import argparse
import json
import math
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

FAMILIES = ('random_recursive', 'binary_catalan', 'ordered_plane_catalan')
CASE = re.compile(r'tree_shape/(\w+)/(inspected_dfs|inspected_root_path)/(\w+)'
                  r'/N:(\d+)/seed:(\d+)/budget:(\d+)')
NS = {'ns': 1, 'us': 1000, 'ms': 1000000, 's': 1000000000}


def collect(paths):
    groups = defaultdict(dict)
    for path in paths:
        samples = defaultdict(list)
        aggregates = {}
        for row in json.loads(path.read_text())['benchmarks']:
            name = row.get('run_name', row['name'])
            match = CASE.fullmatch(name)
            if not match or match[3] not in FAMILIES:
                continue
            impl, work, family, n, seed, budget = match.groups()
            key = ((int(budget), int(n), int(seed), impl, work), family)
            if any(row.get(k) for k in ('skipped', 'error_occurred', 'skip_message', 'error_message')):
                raise ValueError(f'{path}: failed benchmark {name}')
            aggregate = row.get('run_type') == 'aggregate'
            if aggregate and row['aggregate_name'] != 'median':
                continue
            batch = row['paths_per_iteration'] if work == 'inspected_root_path' else 1
            if not math.isfinite(batch) or batch < 1 or batch != int(batch):
                raise ValueError(f'{path}: invalid paths_per_iteration for {name}')
            ns = row['cpu_time'] * NS[row['time_unit']] / batch
            if not math.isfinite(ns) or ns < 0:
                raise ValueError(f'{path}: invalid CPU latency for {name}')
            if aggregate:
                if key in aggregates:
                    raise ValueError(f'{path}: duplicate median for {name}')
                aggregates[key] = ns
            else:
                samples[key].append(ns)

        # Raw rows are normalized before taking the median. Median-only JSON
        # is equivalent because Matched fixes the batch size for each case.
        aggregates.update((key, statistics.median(values)) for key, values in samples.items())
        for (group, family), value in aggregates.items():
            if family in groups[group]:
                raise ValueError(f'Overlapping input datasets: {group} {family}')
            groups[group][family] = value
    return groups


def means(groups, strict=False):
    values = {}
    for group, medians in sorted(groups.items()):
        missing = set(FAMILIES) - medians.keys()
        if missing:
            message = f'{group}: missing {", ".join(sorted(missing))}'
            if strict:
                raise ValueError(message)
            print(message, file=sys.stderr)
            continue
        values[group] = statistics.fmean(medians[f] for f in FAMILIES)
    if not values:
        raise ValueError('No complete tree_shape triplets')
    return values


def render(values):
    lines = ['| Budget (bytes) | Nodes | Seed | Implementation | Workload | Mean CPU latency (ns) |',
             '| ---: | ---: | ---: | --- | --- | ---: |']
    for group, mean in values.items():
        lines.append('| ' + ' | '.join(map(str, (*group, f'{mean:.6g}'))) + ' |')
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', type=Path, nargs='+', help='non-overlapping timed JSON shards')
    parser.add_argument('--strict', action='store_true', help='reject incomplete family groups')
    args = parser.parse_args()
    try:
        print(render(means(collect(args.inputs), strict=args.strict)))
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f'{parser.prog}: {exc}\n')


if __name__ == '__main__':
    main()
