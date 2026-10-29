#!/usr/bin/env python3
"""Validate and report the cache/local/kernel/reference tree comparison sweep."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import zipfile

os.environ.setdefault('MPLCONFIGDIR', '/tmp/pixie-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_style import apply_style

OPS = ['dfs', 'dfs_children', 'random_path', 'random_path_children']
TITLES = ['DFS', 'DFS + child inspection', 'Root-to-leaf path', 'Path + child inspection']
REFS = ['adjlist_preorder', 'pointer_preorder', 'adjlist', 'pointer_bfs']
FAMILIES = ['bp', 'dfuds', 'louds']
BASE = {'bp': 'bp_rmmbtree', 'dfuds': 'dfuds_rmmbtree', 'louds': 'louds_compressed'}
LABEL = {'adjlist_preorder': 'Vectors · preorder', 'adjlist': 'Vectors · BFS',
         'pointer_preorder': 'Pointers · preorder', 'pointer_bfs': 'Pointers · BFS'}
COLORS = {'bp': '#0072B2', 'dfuds': '#D55E00', 'louds': '#009E73',
          'adjlist_preorder': '#777777', 'adjlist': '#777777',
          'pointer_preorder': '#222222', 'pointer_bfs': '#222222'}
SCALE = {'ns': 1., 'us': 1e3, 'ms': 1e6, 's': 1e9}


def name(series, op):
    if series in REFS:
        return series + '_' + op
    family, variant = series.split(':')
    if variant == 'clock':
        return family + '_clock_' + op
    if variant == 'local':
        return BASE[family] + '_' + ('view_dfs' if op == 'dfs' else op)
    return BASE[family] + '_' + variant + '_' + op


def label(series):
    if series in REFS:
        return LABEL[series]
    family, variant = series.split(':')
    return family.upper() + ' · ' + {'clock': 'CLOCK ≤512 KiB', 'local': 'local',
        'unpacked': 'unpacked', 'indexed': 'fixed index ≈1 MiB'}[variant]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    root = args.directory
    manifest = json.loads((root / 'manifest.json').read_text())
    status = json.loads((root / 'status.json').read_text())
    if status['state'] != 'complete':
        raise ValueError('Full report requires a completed sweep')
    groups = {}; raw = []; repetitions = manifest['config']['repetitions']
    for group in manifest['groups']:
        result = root / (group['key'] + '.json')
        marker = json.loads((root / (group['key'] + '.complete.json')).read_text())
        assert hashlib.sha256(result.read_bytes()).hexdigest() == marker['sha256']
        data = json.loads(result.read_text())
        assert not any(r.get('error_occurred') for r in data['benchmarks'])
        rows = [r for r in data['benchmarks'] if r.get('run_type') == 'iteration']
        assert {r.get('run_name', r['name']) for r in rows} == set(group['cases'])
        assert len(rows) == repetitions * len(group['cases'])
        for r in rows:
            assert r['iterations'] > 0 and r['nodes'] == group['nodes']
            assert math.isfinite(r['cpu_time']) and r['cpu_time'] > 0
            if 'visited_nodes_per_iteration' in r and '_dfs_children/' in r['name']:
                assert r['visited_nodes_per_iteration'] == group['nodes']
                assert r['inspected_children_per_iteration'] == group['nodes'] - 1
            if '_clock_' in r['name']:
                assert r['cache_bytes'] <= 512 * 1024
            key = (r.get('run_name', r['name']).split('/')[0], group['nodes'])
            groups.setdefault(key, []).append(r['cpu_time'] * SCALE[r['time_unit']])
            raw.append(r)
    assert all(len(values) == repetitions for values in groups.values())
    points = {key: statistics.median(values) for key, values in groups.items()}
    sizes = sorted(manifest['config']['sizes'])
    assert len(points) == status['total_cases'] == 580
    assert len(raw) == 580 * repetitions
    with (root / 'medians.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['benchmark', 'nodes', 'median_cpu_ns', 'min_cpu_ns', 'max_cpu_ns', 'repetitions'])
        for (benchmark, n), value in sorted(points.items()):
            values = groups[(benchmark, n)]
            writer.writerow([benchmark, n, value, min(values), max(values), len(values)])

    figures = root / 'figures'; figures.mkdir(exist_ok=True)
    apply_style()
    plt.rcParams.update({'font.size': 10})
    def plot(stem, title, series, ops=OPS):
        fig, axes = plt.subplots(1 if len(ops) == 2 else 2, 2, figsize=(12, 5 if len(ops) == 2 else 8))
        for ax, op in zip(axes.flat, ops):
            for item in series:
                ns = [n for n in sizes if (name(item, op), n) in points]
                values = [points[(name(item, op), n)] / (n if op.startswith('dfs') else 1) for n in ns]
                family = item.split(':')[0]
                style = '--' if item in ('adjlist', 'pointer_bfs') or item.endswith(':local') else '-'
                ax.plot(ns, values, label=label(item), color=COLORS[family], linestyle=style,
                        marker='s' if item.startswith('pointer') else 'o', markersize=3, linewidth=1.6)
            ax.set_xscale('log', base=2); ax.set_yscale('log')
            ax.set_xticks([2**e for e in range(10, 27, 2)])
            ax.set_title(TITLES[OPS.index(op)])
            ax.set_xlabel('Nodes'); ax.set_ylabel('ns / node' if op.startswith('dfs') else 'ns / complete path')
            ax.grid(True, which='major', alpha=.2)
        handles, labels = axes.flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', ncol=3, fontsize=9)
        fig.suptitle(title)
        fig.tight_layout(rect=(0, .13 if len(ops)==2 else .10, 1, .95))
        for extension in ['png', 'svg', 'pdf']:
            fig.savefig(figures / f'{stem}.{extension}', dpi=180)
        plt.close(fig)
    plot('clock_vs_references', '512 KiB CLOCK policy versus reference trees', REFS + [f+':clock' for f in FAMILIES])
    plot('local_vs_clock', 'Local memoization versus selective CLOCK admission', [f+':'+v for f in FAMILIES for v in ['local','clock']])
    plot('unpacked_vs_references', 'Existing unpacked kernels versus reference trees', REFS + [f+':unpacked' for f in FAMILIES])
    plot('indexed_paths_vs_references', 'Existing fixed-index paths versus reference trees (≈1 MiB index)', REFS + [f+':indexed' for f in FAMILIES], OPS[2:])

    n = 10000000
    def timing_table(series):
        lines = ['| Implementation | DFS (ms) | DFS + children (ms) | Path (ns) | Path + children (ns) |',
                 '| --- | ---: | ---: | ---: | ---: |']
        for item in series:
            cells = []
            for op in OPS:
                value = points.get((name(item, op), n))
                cells.append('—' if value is None else f'{value / (1e6 if op.startswith("dfs") else 1):.2f}')
            lines.append('| ' + label(item) + ' | ' + ' | '.join(cells) + ' |')
        return '\n'.join(lines)

    relative = []
    for family in FAMILIES:
        matched = REFS[2:] if family == 'louds' else REFS[:2]
        for size in sizes:
            for op in OPS:
                value = points[(name(family+':clock', op), size)]
                for ref in matched:
                    relative.append([family, size, op, ref, points[(name(ref, op), size)] / value])
    with (root/'matched_reference_speedups.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['tree', 'nodes', 'workload', 'reference', 'reference_time_divided_by_clock_time'])
        writer.writerows(relative)
    best = ['| Workload | Fastest reference | Reference time | Fastest measured succinct strategy | Succinct time | Reference / succinct |',
            '| --- | --- | ---: | --- | ---: | ---: |']
    all_succinct = [f+':'+v for f in FAMILIES for v in ['clock','local','unpacked','indexed']]
    for op in OPS:
        ref = min(REFS, key=lambda r: points[(name(r, op), n)])
        candidates = [s for s in all_succinct if (name(s, op), n) in points]
        candidate = min(candidates, key=lambda s: points[(name(s, op), n)])
        r = points[(name(ref,op), n)]; c = points[(name(candidate,op), n)]
        scale = 1e6 if op.startswith('dfs') else 1; unit = 'ms' if scale==1e6 else 'ns'
        best.append(f'| {TITLES[OPS.index(op)]} | {label(ref)} | {r/scale:.2f} {unit} | {label(candidate)} | {c/scale:.2f} {unit} | {r/c:.2f}x |')
    clock_wins = [op for op in OPS if
        min(points[(name(f+':clock', op), n)] for f in FAMILIES) <
        min(points[(name(ref, op), n)] for ref in REFS)]
    conclusion = ('At 10 million nodes, the new 512 KiB CLOCK policy remains slower '
        'than the fastest reference on all four workloads. Existing unpacked '
        'kernels and fixed-index paths must be assessed separately; the winner '
        'table below shows where they beat the references.' if not clock_wins else
        'At 10 million nodes, a CLOCK variant beats the fastest reference for: ' +
        ', '.join(clock_wins) + '.')
    recheck_notes = []
    if manifest.get('rechecks'):
        recheck_notes = ['', '## Confirmation runs', '',
            'The accepted tables and exports use the confirmation runs. '
            'Original records remain under `rechecks/original/`. Query CPU '
            'timings are distinct from whole-process wall time.', '']
        recheck_notes += ['- ' + r['group'] + ': ' + r['reason']
                          for r in manifest['rechecks']]
    lines = ['# Succinct trees versus reference implementations',
        '', conclusion,
        '', 'Fresh comparison of the selective-admission CLOCK policy (512 KiB cap), local-cache controls, existing unpacking kernels and fixed-root-index paths, and both references in preorder and BFS layouts.',
        '', '## Results at 10 million nodes', '', timing_table(REFS+[f+':clock' for f in FAMILIES]),
        '', '## Existing optimized strategies and local controls', '', timing_table([f+':'+v for f in FAMILIES for v in ['local','unpacked','indexed']]),
        '', '## Fastest measured strategies at 10 million nodes', '', '\n'.join(best),
        '', 'This table chooses a winner separately for each workload; it does not imply one implementation/layout wins everywhere. Fixed-index paths use a different, approximately 1 MiB auxiliary index, not the new 512 KiB CLOCK cache.',
        '', '## Scaling graphs', '', 'Lines are medians of five repetitions. There are no shaded bands. DFS plots use ns/node; path plots use ns/complete path. The tables above use ms/complete DFS traversal.',
        '', '![CLOCK versus references](figures/clock_vs_references.png)',
        '', '![Local versus CLOCK](figures/local_vs_clock.png)',
        '', '![Unpacked versus references](figures/unpacked_vs_references.png)',
        '', '![Fixed-index paths versus references](figures/indexed_paths_vs_references.png)',
        '', '## Method and interpretation', '',
        '- Release (`-O3 -DNDEBUG -march=native`), CPU 0, five randomly interleaved repetitions per group, requested minimum 0.2 s and warmup 0.1 s. Long DFS iterations exceed that minimum.',
        '- Ten sizes: 2^10, 2^12, …, 2^26, plus 10,000,000. 70 isolated process groups, 580 cases, 2,900 accepted individual repetitions. No skipped or failed cases.',
        '- Tree seed 42; shifted 64 KiB query pool, seed 43. Construction is excluded. Adaptive iteration counts can differ; path counters retain the average visited/inspected nodes.',
        '- BP and DFUDS use RmMBTree with compressed rank/select and preorder numbering. LOUDS uses compressed rank/select and BFS numbering. Match BP/DFUDS to preorder references, and LOUDS to BFS references when comparing layouts.',
        '- Vector adjacency and contiguous three-pointer records are separate references. Both layouts are measured; retaining all large reference fixtures in one process is avoided.',
        '- CLOCK uses operation-based admission: sequential leaf/parent/sibling/first-child navigation stays local. Degree/indexed-child queries activate shared caching in a view; enumeration then uses cached data when available. Reset clears activation.',
        '- CLOCK rows use ordinary view/path traversal, not unpacking. Unpacked and fixed-index rows remain separately labeled. The latter prebuild a probability-prioritized root-connected index; it is not an eviction cache.',
        '- Cache-byte counters exclude underlying tree storage and allocator bookkeeping. Resource logs measure whole-process setup/query peak memory; they are not bytes-per-node measurements of the data structure. Pointer references also expose their owned bytes explicitly.',
        '- Timing is CPU-time median on this host. Small differences should not be treated as universal rankings; raw repetition ranges are retained.',
        *recheck_notes,
        '', '## Saved data', '',
        '- `manifest.json`: compiler, CPU, source and executable hashes, configuration and selected cases.',
        '- `status.json`: completed collection status.',
        '- Per-group JSON/log/command/resource/completion files retain all raw data.',
        '- `combined.json`, `results.csv`: all validated raw repetitions and aggregates; filter `run_type == iteration` for statistics.',
        '- `medians.csv`: 580 medians and repetition ranges, normalized to nanoseconds.',
        '- `matched_reference_speedups.csv`: CLOCK speedups against each layout-matched reference across all sizes; >1 means CLOCK is faster.',
        '- `bin/`, `sources.zip`, `collector.py`: measured executables and inputs. Figures are also saved as SVG and PDF.',
        '', 'Collection command:', '', '```sh',
        'python3 scripts/benchmark_trees.py --succinct-mode comparison --nodes 1024 4096 16384 65536 262144 1048576 4194304 10000000 16777216 67108864 --first-size 10000000 --random-interleaving --repetitions 5 --min-time 0.2 --warmup 0.1 --output-dir profiling/cache-vs-references-2026-10-24',
        '```', '', 'Append `--resume` to recover completed groups. The collector verifies identical settings/executables and result hashes.', '']
    (root/'REPORT.md').write_text('\n'.join(lines))
    audit = {'rechecked_groups': [r['group'] for r in manifest.get('rechecks', [])], 'groups': len(manifest['groups']), 'cases': len(points), 'raw_repetitions': len(raw),
             'sizes': sizes, 'max_clock_bytes': max(r.get('cache_bytes',0) for r in raw),
             'validation': 'All completed file hashes, expected cases/repetitions, node counts, child-inspection counts, and CLOCK byte bounds verified.'}
    (root/'report_validation.json').write_text(json.dumps(audit,indent=2)+'\n')
    with zipfile.ZipFile(root/'report-and-graphs.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
        for item in ['REPORT.md','medians.csv','matched_reference_speedups.csv','report_validation.json']:
            archive.write(root/item,item)
        for path in figures.iterdir():
            archive.write(path,str(path.relative_to(root)))
    print(timing_table(REFS+[f+':clock' for f in FAMILIES])); print('\n'.join(best))


if __name__ == '__main__':
    main()
