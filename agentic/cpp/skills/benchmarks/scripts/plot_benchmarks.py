#!/usr/bin/env python3
"""Generic Google Benchmark JSON plotter.

Plots timing and (optionally) hardware counter data from Google Benchmark
--benchmark_format=json output. Correctly handles repetition output by
averaging raw iterations and skipping aggregate entries.

Usage:
    python3 plot_benchmarks.py results.json [output_prefix]
    python3 plot_benchmarks.py results.json report --min-n 16 --max-n 1048576
"""
import json
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np

def load_benchmark_json(path):
    with open(path, 'r') as f:
        return json.load(f)

def extract_series(data, metric='time', min_n=None, max_n=None):
    """Extract series from benchmark JSON.

    metric='time'     -> real_time (ns)
    metric='counter'  -> first available hardware counter (e.g. CACHE-MISSES)
    """
    # Determine counter key if requested
    counter_key = None
    if metric == 'counter':
        # Find first counter key from the first iteration entry
        for bench in data.get('benchmarks', []):
            if bench.get('run_type', 'iteration') != 'iteration':
                continue
            for k in bench.keys():
                if k not in ('name', 'run_name', 'run_type', 'repetitions',
                             'repetition_index', 'threads', 'iterations',
                             'real_time', 'cpu_time', 'time_unit',
                             'items_per_second', 'aggregate_name',
                             'aggregate_unit', 'family_index',
                             'per_family_instance_index'):
                    counter_key = k
                    break
            if counter_key:
                break
        if not counter_key:
            return {}

    raw = {}
    for bench in data.get('benchmarks', []):
        if bench.get('run_type', 'iteration') != 'iteration':
            continue

        name = bench['name']
        parts = name.split('/')
        if len(parts) < 2:
            continue

        bench_name = parts[0]
        try:
            n = int(parts[1])
        except ValueError:
            continue

        if min_n is not None and n < min_n:
            continue
        if max_n is not None and n > max_n:
            continue

        if metric == 'time':
            val = bench.get('real_time', bench.get('cpu_time', 0))
        else:
            val = bench.get(counter_key)

        if val is None or val == 0:
            continue

        key = (bench_name, n)
        raw.setdefault(key, []).append(val)

    series = {}
    for (bench_name, n), vals in raw.items():
        series.setdefault(bench_name, []).append((n, sum(vals) / len(vals)))

    for name in series:
        series[name].sort(key=lambda x: x[0])

    return series

def plot_series(series, output_prefix, ylabel, title_suffix):
    if not series:
        print(f"No data to plot for {title_suffix}")
        return

    fig, ax = plt.subplots(figsize=(12, 8))
    colors = plt.cm.tab10(np.linspace(0, 1, len(series)))

    for idx, (name, points) in enumerate(sorted(series.items())):
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        ax.plot(xs, ys, marker='o', markersize=3, label=name, color=colors[idx])

    ax.set_xscale('log')
    ax.set_xlabel('Benchmark parameter n')
    ax.set_ylabel(ylabel)
    ax.set_title(f'Benchmark Results - {title_suffix}')
    ax.legend(loc='upper left', fontsize='small')
    ax.grid(True, which='both', ls='--', alpha=0.5)
    fig.tight_layout()
    fig.savefig(f'{output_prefix}.png', dpi=150)
    fig.savefig(f'{output_prefix}.svg')
    print(f'Saved {output_prefix}.png and {output_prefix}.svg')
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(
        description='Plot Google Benchmark JSON results.')
    parser.add_argument('json_path', help='Path to benchmark JSON file')
    parser.add_argument('output_prefix', nargs='?', default='report',
                        help='Output file prefix (default: report)')
    parser.add_argument('--min-n', type=int, default=None,
                        help='Minimum parameter value to include (inclusive)')
    parser.add_argument('--max-n', type=int, default=None,
                        help='Maximum parameter value to include (inclusive)')
    args = parser.parse_args()

    data = load_benchmark_json(args.json_path)

    # Timing plot
    time_series = extract_series(data, metric='time',
                                 min_n=args.min_n, max_n=args.max_n)
    if time_series:
        plot_series(time_series, args.output_prefix,
                    'Time per iteration (ns)', 'Time')

    # Counter plot
    counter_series = extract_series(data, metric='counter',
                                    min_n=args.min_n, max_n=args.max_n)
    if counter_series:
        # Determine counter name for labels
        counter_name = 'Counter'
        for bench in data.get('benchmarks', []):
            if bench.get('run_type', 'iteration') != 'iteration':
                continue
            for k in bench.keys():
                if k not in ('name', 'run_name', 'run_type', 'repetitions',
                             'repetition_index', 'threads', 'iterations',
                             'real_time', 'cpu_time', 'time_unit',
                             'items_per_second', 'aggregate_name',
                             'aggregate_unit', 'family_index',
                             'per_family_instance_index'):
                    counter_name = k
                    break
            break
        plot_series(counter_series, f'{args.output_prefix}_{counter_name.lower().replace("-", "_")}',
                    f'{counter_name} per iteration', counter_name)

if __name__ == '__main__':
    main()
