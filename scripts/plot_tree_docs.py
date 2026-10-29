#!/usr/bin/env python3
"""Six documentation plots for current tree traversal, construction and space."""
import argparse
import csv
import json
import os
import statistics
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/pixie-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
import plot_style as style

SERIES = [
    ("adjlist_preorder", "vectors_preorder", "Vectors · preorder", "#777777", "o", "-"),
    ("adjlist", "vectors_bfs", "Vectors · BFS", "#777777", "o", "--"),
    ("pointer_preorder", "pointers_preorder", "Pointers · preorder", "#222222", "s", "-"),
    ("pointer_bfs", "pointers_bfs", "Pointers · BFS", "#222222", "s", "--"),
    ("bp_best", "bp", "BP", "#0072B2", "^", "-"),
    ("dfuds_best", "dfuds", "DFUDS", "#D55E00", "D", "-"),
    ("louds_best", "louds", "LOUDS", "#009E73", "v", "-"),
]
OPS = [("dfs", "DFS"), ("dfs_children", "DFS + child inspection"),
       ("random_path", "Root-to-leaf path"),
       ("random_path_children", "Path + child inspection")]
SCALE = {"ns":1, "us":1000, "ms":1000000, "s":1000000000}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--cache32', type=Path, nargs='+', required=True)
    parser.add_argument('--cache512', type=Path, nargs='+', required=True)
    parser.add_argument('--construction', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    query = {}
    inputs = [(0,[args.reference]),(32768,args.cache32),(524288,args.cache512)]
    for budget, paths in inputs:
        seen = set()
        for path in paths:
            current = set()
            for row in json.loads(path.read_text())["benchmarks"]:
                if row.get('error_occurred'):
                    raise ValueError(f'Benchmark error in {path}')
                if row.get('run_type') != 'iteration':
                    continue
                name = row.get('run_name',row['name']).split('/')[0]
                if ('_best_' in name) != bool(budget):
                    continue
                n = int(row['size_nodes'])
                if budget and (row['cache_limit_bytes'] != budget or row['cache_and_scratch_bytes'] > budget):
                    raise ValueError('Invalid runtime cache budget')
                key = budget,name,n
                if key in seen:
                    raise ValueError(f'Overlapping input datasets: {key}')
                current.add(key)
                query.setdefault(key,[]).append(row['cpu_time']*SCALE[row['time_unit']])
            seen.update(current)
    construction = {}
    for row in json.loads(args.construction.read_text()):
        if row['budget'] and row['persistent_cache_bytes'] + row['scratch_bytes'] > row['budget']:
            raise ValueError('Construction cache budget exceeded')
        construction.setdefault((row['budget'],row['tree'],row['nodes']),[]).append(row)
    style.apply_style()
    records = []

    def draw(operation, title):
        build = operation == 'build'
        memory = operation == 'memory'
        fig, axes = plt.subplots(1,2 if build else 1,squeeze=False,
                                 figsize=(13 if build else 10,6.8))
        handles = {}
        sizes = [2**e for e in range(10,27,2)]
        series = [(item,0) for item in SERIES[:4]]
        for prefix,name,label,color,marker,_ in SERIES[4:]:
            for budget in (32768,524288):
                series.append(((prefix,name,f'{label} · {budget//1024} KiB',
                                color,marker,':' if budget==32768 else '-'),budget))
        for row_index, ax in enumerate(axes.flat):
            for (prefix, name, label, color, marker, dash), effective in series:
                ys = []
                available = sorted(n for b,key,n in (construction if build or memory else query)
                                   if b == effective and key == (name if build or memory else prefix+'_'+operation))
                if available != sizes:
                    raise ValueError(f'Expected all even powers 10 through 26 for {label}')
                for n in available:
                    if build or memory:
                        samples = construction[effective,name,n]
                        if len(samples) != 3:
                            raise ValueError('Expected three construction repetitions')
                        if memory:
                            values = [r['working_bytes']/n for r in samples]
                        elif row_index == 0:
                            values = [(r['construct_s']+r['cache_s'])*1000 for r in samples]
                        else:
                            values = [r['total_s']*1000 for r in samples]
                    else:
                        values = query[effective,prefix+'_'+operation,n]
                        if len(values) != 5:
                            raise ValueError('Expected five traversal repetitions')
                        if operation.startswith('dfs'):
                            values = [v/n for v in values]
                    median = statistics.median(values)
                    ys.append(median)
                    records.append([operation,row_index,effective,label,n,median,len(values)])
                line, = ax.plot(available,ys,label=label,color=color,marker=marker,linestyle=dash,
                                markerfacecolor='white' if effective == 32768 else color)
                handles.setdefault(label,line)
            if build:
                ax.set_title('Prepared input' if row_index==0 else 'Full setup')
            ax.set_xlabel('Tree size (nodes)')
            ax.set_xscale('log',base=2)
            ax.set_xticks(sizes,[rf'$2^{{{n.bit_length()-1}}}$' for n in sizes])
            ax.xaxis.set_minor_locator(NullLocator())
            ax.tick_params(axis='y',labelleft=True)
            if build:
                ax.set_ylabel('Wall time per build')
                style.build_time_y_axis(ax)
            elif memory:
                ax.set_ylabel('Bytes/node · including reserved traversal scratch')
                style.log_y_axis(ax)
                ax.yaxis.set_minor_locator(NullLocator())
            else:
                ax.set_ylabel('CPU time (ns/node)' if operation.startswith('dfs') else 'CPU time (ns/path)')
                ax.set_ylim(bottom=0)
        fig.suptitle('Trees — '+title,y=.99)
        fig.legend(list(handles.values()),list(handles),loc='lower center',ncol=4,bbox_to_anchor=(.5,.015))
        fig.tight_layout(rect=(0,.20,1,.96))
        for ext in ('svg','png','pdf'):
            fig.savefig(args.output_dir/f'tree_{operation}.{ext}',dpi=style.DPI)
        plt.close(fig)

    for op,title in OPS:
        draw(op,title+' · median CPU time')
    draw('build','Construction · median wall time')
    draw('memory','Memory usage')
    with (args.output_dir/'plot_data.csv').open('w',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['operation','panel_row','budget_bytes','series','nodes','median','repetitions'])
        writer.writerows(records)
    print('Saved six tree documentation plots')

if __name__ == '__main__':
    main()
