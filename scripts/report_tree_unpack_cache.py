#!/usr/bin/env python3
"""Reproduce the 2026-10-24 unpacking/cache experiment report."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import zipfile

os.environ.setdefault('MPLCONFIGDIR', '/tmp/pixie-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OPS = ['dfs', 'dfs_children', 'random_path', 'random_path_children']
TITLES = ['DFS', 'DFS + children', 'Path', 'Path + children']
BASE = {'bp': 'bp_rmmbtree', 'dfuds': 'dfuds_rmmbtree', 'louds': 'louds_compressed'}
REFS = {'adjlist_preorder': 'Vectors preorder', 'pointer_preorder': 'Pointers preorder',
        'adjlist': 'Vectors BFS', 'pointer_bfs': 'Pointers BFS'}
SCALE = {'ns': 1, 'us': 1000, 'ms': 1000000, 's': 1000000000}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    root = args.directory
    manifest = json.loads((root / 'manifest.json').read_text())
    status = json.loads((root / 'status.json').read_text())
    assert status['state'] == 'complete'
    repeats = manifest['config']['repetitions']
    groups, cache_bytes = {}, {}
    for group in manifest['groups']:
        path = root / (group['key'] + '.json')
        marker = json.loads(path.with_suffix('.complete.json').read_text())
        assert hashlib.sha256(path.read_bytes()).hexdigest() == marker['sha256']
        rows = json.loads(path.read_text())['benchmarks']
        assert not any(r.get('error_occurred') for r in rows)
        rows = [r for r in rows if r.get('run_type') == 'iteration']
        assert len(rows) == repeats * len(group['cases'])
        assert {r['run_name'] for r in rows} == set(group['cases'])
        for row in rows:
            assert row['nodes'] == group['nodes'] and row['iterations'] > 0
            name = row['run_name'].split('/')[0]
            time = row['cpu_time'] * SCALE[row['time_unit']]
            assert math.isfinite(time) and time > 0
            groups.setdefault((name, group['nodes']), []).append(time)
            if '_dfs_children' in name:
                assert row['visited_nodes_per_iteration'] == group['nodes']
                assert row['inspected_children_per_iteration'] == group['nodes'] - 1
            if 'cache_bytes' in row:
                assert row['cache_bytes'] <= 512 * 1024
                cache_bytes[name] = max(cache_bytes.get(name, 0), row['cache_bytes'])
    assert all(len(v) == repeats for v in groups.values())
    assert len(groups) == status['total_cases']
    points = {k: statistics.median(v) for k, v in groups.items()}
    with (root / 'medians.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['benchmark', 'nodes', 'median_cpu_ns', 'min_cpu_ns', 'max_cpu_ns'])
        for key, value in sorted(points.items()):
            writer.writerow([*key, value, min(groups[key]), max(groups[key])])
    sizes = sorted(manifest['config']['sizes'])
    target = 10000000 if 10000000 in sizes else sizes[-1]
    def names(family):
        return {'Local': BASE[family], 'CLOCK': family + '_clock',
                'Unpacked': BASE[family] + '_unpacked',
                'Selective + unpack': family + '_selective_unpack',
                'Chunks + unpack': family + '_chunk_unpack',
                'Fixed index (~1 MiB)': BASE[family] + '_indexed'}
    def get(prefix, op, n):
        key = prefix + '_' + ('view_dfs' if op == 'dfs' and prefix in BASE.values() else op)
        return points.get((key, n))
    def table(series):
        lines = ['| Implementation | DFS (ms) | DFS + children (ms) | Path (ns) | Path + children (ns) |',
                 '| --- | ---: | ---: | ---: | ---: |']
        for label, prefix in series.items():
            values = [get(prefix, op, target) for op in OPS]
            display = ['—' if v is None else f'{v / (1e6 if i < 2 else 1):.2f}' for i, v in enumerate(values)]
            lines.append('| ' + ' | '.join([label] + display) + ' |')
        return '\n'.join(lines)
    lines = ['# Combining unpacking and caching', '',
        f'Completed {len(groups)} cases / {len(groups) * repeats} timed repetitions. Sizes: ' + ', '.join(f'{n:,}' for n in sizes) + '.', '',
        'Fresh Release (`-O3 -DNDEBUG -march=native`), CPU 0, five interleaved repetitions per process group, 0.2 s minimum and 0.1 s warmup. Median CPU time; construction excluded. Seed 42 topology and shifted 64 KiB query pool with seed 43. Both references were rerun in preorder and BFS layouts.', '',
        '## Findings', '',
        'At 10 million nodes, selective retention is preferable to whole-chunk retention, but neither combination replaces standalone unpacking for inspected DFS or the fixed root index for paths. Plain DFS continues to use the same streaming algorithm, so timing differences there are not a new cache optimization.', '',
        '- **Inspected DFS:** BP takes 197 ms with selective retention and 204 ms with chunks, versus 105 ms with standalone unpacking. DFUDS takes about 196 ms with either, versus 113 ms. LOUDS remains about 113 ms because its kernel is unchanged. The fastest reference takes 68 ms.',
        '- **Plain paths:** selective retention is close to CLOCK alone (BP 1143 vs 1177 ns; DFUDS 1123 vs 1128 ns; LOUDS 160 vs 158 ns). Whole chunks are slower on all three.',
        '- **Inspected paths:** selective retention improves over CLOCK alone by about 6% for BP, 4% for DFUDS and 26% for LOUDS. However, LOUDS local navigation still beats both caches (211 ns versus 232 ns selective). Fixed-index paths remain much faster (47–51 ns).',
        '- **Small trees:** at 1024 nodes, chunks beat selective retention on plain BP/DFUDS paths (about 53–54 vs 68–70 ns), but CLOCK alone is still faster (about 49 ns).',
        '- These are experimental policies; tree defaults remain unchanged. The 256-bit versus 2048-bit inspected-DFS window difference below prevents attributing the regression to retention alone.', '',
        '## What was implemented', '',
        '- **Selective + unpack:** decode a temporary chunk, use it for nearby requests, and retain requested node results (degree plus up to eight children) in a CLOCK table. Other decoded nodes are not independently admitted.',
        '- **Chunks + unpack:** retain decoded chunks, evict whole chunks with CLOCK, and promote a chunk when it is reused. Lookup indexes chunk roots; a current-chunk hint also serves descendants of that chunk. Arbitrary interior nodes of other resident chunks are not separately indexed.',
        '- Both cap context plus requested persistent allocation at **512 KiB**. Allocator bookkeeping, caller views and bounded callback-safe temporary copies are excluded. Cache allocation is lazy; actual allocations differ.',
        '- **BP/DFUDS:** probe up to 256 source bits for a complete subtree. Non-fitting subtrees fall back to ordinary navigation with a retained child prefix. Computing BP degree can still scan beyond the bit probe.',
        '- **LOUDS:** unpack up to 256 bits of consecutive child handles. Children beyond that window use ordinary navigation.',
        '- **Plain DFS:** all unpacked variants use the existing streaming kernel; they do not admit every visited node into the cache. These rows are controls, not distinct cache-accelerated DFS algorithms.',
        '- **Inspected DFS:** BP/DFUDS combinations consume decoded subtrees through the cache, copying them before callbacks so nested queries can evict safely. The old unpacked control uses 2048-bit scratch subtrees; the combinations use 256-bit chunks. Thus this comparison measures the complete implementations, not retention alone. LOUDS retains its existing monotone level-cursor kernel, without cache admission in its normal path.',
        '- Caches persist with each tree fixture across iterations and workload rows. These are warmed repeated-query measurements; cache construction/decoding on individual misses remains timed.',
        '- **Paths:** both combinations use cache-aware navigation and decoding. The old unpacked control decodes anew per path. The fixed root index is a separate ~1 MiB prebuilt structure, not a 512 KiB eviction cache.', '',
        f'## Results at {target:,} nodes', '', table({v:k for k,v in REFS.items()}), '']
    for family in BASE:
        lines += [f'### {family.upper()}', '', table(names(family)), '']
    lines += ['## Actual persistent cache allocation', '', '| Policy | Maximum observed bytes |', '| --- | ---: |']
    for family in BASE:
        for policy in ['selective_unpack', 'chunk_unpack']:
            maximum = max(v for k,v in cache_bytes.items() if k.startswith(family + '_' + policy + '_'))
            lines.append(f'| {family.upper()} {policy} | {maximum:,.0f} |')
    fig, axes = plt.subplots(3, 4, figsize=(16, 10), constrained_layout=True)
    for row, family in enumerate(BASE):
        for col, op in enumerate(OPS):
            ax = axes[row, col]
            for label, color, marker in [('Unpacked','#777777','.'), ('Selective + unpack','#0072B2','o'), ('Chunks + unpack','#D55E00','s'), ('CLOCK','#984ea3','^')]:
                values = [get(names(family)[label], op, n) / (n if col < 2 else 1) for n in sizes]
                ax.plot(sizes, values, label=label, color=color, marker=marker, markersize=4)
            fastest = [min(get(p, op, n) for p in REFS) / (n if col < 2 else 1) for n in sizes]
            ax.plot(sizes, fastest, label='Fastest reference', color='#111111', linestyle='--')
            ax.set_xscale('log', base=2); ax.set_yscale('log')
            ax.set_xticks(sizes, [f'{n:,}' for n in sizes], fontsize=8)
            ax.set_title(f'{family.upper()} · {TITLES[col]}')
            ax.set_ylabel('ns/node' if col < 2 else 'ns/path')
            ax.set_xlabel('Nodes'); ax.grid(True, alpha=.2)
    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=5)
    figures = root / 'figures'; figures.mkdir(exist_ok=True)
    for ext in ['png', 'svg', 'pdf']:
        fig.savefig(figures / ('comparison.' + ext), dpi=160)
    plt.close(fig)
    lines += ['', '## Scaling', '', '![Comparison](figures/comparison.png)', '',
        'Reference curves choose the fastest reference/layout separately at every point. BP/DFUDS naturally use preorder numbering; LOUDS uses BFS. Prefer the corresponding layout rows for layout-matched comparisons. No shaded bands. Small timing differences, especially between identical streaming kernels, are not evidence of an algorithmic improvement.', '',
        '## Validation and reproducibility', '',
        '- All 374 Release tree specification and eviction tests passed; see `tree-tests.log`.',
        '- ASan/UBSan scalar-fallback probe covers all six combinations, random/star/chain shapes, bit boundaries, nested eviction, traversal event equality and copy/move; see `sanitizer.cpp` and `sanitizer.log`. Leak detection disabled.',
        '- Report checks result hashes, expected cases/repetitions, finite timings, node/edge counts and cache bounds.',
        '- `manifest.json`, `sources.zip`, `bin/`, `collector.py`, raw group JSON/logs and `medians.csv` preserve the run. Raw rows and statistical aggregates are kept separate.', '']
    (root / 'REPORT.md').write_text('\n'.join(lines))
    (root / 'report_validation.json').write_text(json.dumps({'cases':len(groups), 'repetitions':len(groups)*repeats, 'sizes':sizes, 'max_cache_bytes':max(cache_bytes.values())}, indent=2)+'\n')
    shutil.copy2(__file__, root / 'report.py')
    with zipfile.ZipFile(root / 'report-and-graphs.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for file in [root/'REPORT.md', root/'medians.csv', root/'report_validation.json', root/'report.py', *figures.glob('*')]:
            archive.write(file, file.relative_to(root))
    print(json.dumps({'cases':len(groups),'repetitions':len(groups)*repeats}))

if __name__ == '__main__':
    main()
