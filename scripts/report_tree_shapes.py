#!/usr/bin/env python3
"""Report matched tree_shape CPU latencies from timed Google Benchmark JSON."""
import argparse
import json
import math
import re
import statistics
import sys
import tempfile
from collections import defaultdict
from itertools import product
from pathlib import Path

FAMILIES = ('random_recursive', 'binary_catalan', 'ordered_plane_catalan')
CASE = re.compile(r'tree_shape/(\w+)/(inspected_dfs|inspected_root_path)/(\w+)'
                  r'/N:(\d+)/seed:(\d+)/budget:(\d+)')
NS = {'ns': 1, 'us': 1000, 'ms': 1000000, 's': 1000000000}
BEGIN = '<!-- tree_shape results begin -->'
END = '<!-- tree_shape results end -->'
SIZES = tuple(2**i for i in range(10, 27, 2))
BUDGETS = (0, 32768, 524288)
IMPLEMENTATIONS = {'bp': 'BP', 'dfuds': 'DFUDS', 'louds': 'LOUDS', 'forest_cover': 'Forest cover'}
WORKLOADS = {'inspected_dfs': ('Whole inspected DFS', 'ms', 1e6),
             'inspected_root_path': ('Inspected root-to-leaf path', 'ns', 1)}


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


def means(groups, publish=False):
    values = {}
    for group, medians in sorted(groups.items()):
        missing = set(FAMILIES) - medians.keys()
        if missing:
            message = f'{group}: missing {", ".join(sorted(missing))}'
            if publish:
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


def updated_page(page, section):
    if page.count(BEGIN) != 1 or page.count(END) != 1:
        raise ValueError('Expected one tree_shape results section')
    before, rest = page.split(BEGIN)
    old, after = rest.split(END)

    def keys(text):
        return {line.rsplit('|', 2)[0] for line in text.splitlines() if line.startswith('| ')}

    if keys(old) - keys(section):
        raise ValueError('Input would remove published rows; supply the full timed dataset')
    return before + BEGIN + '\n\n' + section + '\n\n' + END + after


def publish(path, values):
    expected = set(product(BUDGETS, SIZES, [42], IMPLEMENTATIONS, WORKLOADS))
    if values.keys() != expected:
        raise ValueError(f'Publication matrix mismatch: {len(expected - values.keys())} missing, '
                         f'{len(values.keys() - expected)} unexpected groups')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style as style

    files = {work: f'tree_shape_{work}.svg' for work in WORKLOADS}
    images = [f'\\image html {file} "{WORKLOADS[work][0]} CPU latency by traversal budget"'
              for work, file in files.items()]
    page = updated_page(path.read_text(), '\n\n'.join([*images, render(values)]))
    destination = path.parent.parent / 'images/benchmarks/trees'
    destination.mkdir(parents=True, exist_ok=True)
    style.apply_style()
    matplotlib.rcParams['svg.hashsalt'] = 'pixie-tree-shapes'
    with tempfile.TemporaryDirectory(dir=destination) as tmp:
        for work, (label, unit, divisor) in WORKLOADS.items():
            fig, axes = plt.subplots(1, 3, figsize=(13, 4.4), sharey=True)
            for ax, budget in zip(axes, BUDGETS):
                for i, (impl, title) in enumerate(IMPLEMENTATIONS.items()):
                    ys = [values[budget, n, 42, impl, work] / divisor for n in SIZES]
                    ax.plot(SIZES, ys, label=title, **style.line_style(impl, i))
                ax.set_title(f'{budget // 1024} KiB budget')
                ax.set_xlabel('Tree size (nodes)')
                ax.set_xscale('log', base=2)
                ax.set_xticks(SIZES[::2], [rf'$2^{{{n.bit_length()-1}}}$' for n in SIZES[::2]])
            # DFS covers several orders of magnitude; path times are easier
            # to compare on a linear axis. All budget panels share one scale.
            if work == 'inspected_dfs' and all(v > 0 for k, v in values.items() if k[-1] == work):
                axes[0].set_yscale('log')
            else:
                axes[0].set_ylim(bottom=0)
            axes[0].set_ylabel(f'Mean CPU latency ({unit} per {label.lower()})')
            fig.suptitle(f'{label}: equal-weight mean of three family medians')
            fig.legend(*axes[0].get_legend_handles_labels(), loc='lower center', ncol=4)
            fig.tight_layout(rect=(0, .09, 1, .98))
            fig.savefig(Path(tmp) / files[work], metadata={'Date': None})
            plt.close(fig)
        for file in files.values():
            (Path(tmp) / file).replace(destination / file)
    path.write_text(page)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', type=Path, nargs='+', help='non-overlapping timed JSON shards')
    parser.add_argument('--page', type=Path, help='replace the marked section in the benchmark page')
    args = parser.parse_args()
    try:
        values = means(collect(args.inputs), publish=bool(args.page))
        if args.page:
            publish(args.page, values)
        else:
            print(render(values))
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f'{parser.prog}: {exc}\n')


if __name__ == '__main__':
    main()
