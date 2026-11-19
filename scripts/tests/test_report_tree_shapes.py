import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import report_tree_shapes as report


def row(shape, cpu=1, *, work='inspected_dfs', impl='bp', n=1024, seed=42,
        budget=0, unit='ns', **fields):
    name = f'tree_shape/{impl}/{work}/{shape}/N:{n}/seed:{seed}/budget:{budget}'
    return dict(name=name, run_name=name, cpu_time=cpu, time_unit=unit,
                real_time=999, iterations=100, **fields)


class TreeReportTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.path = Path(temp.name) / 'results.json'

    def collect(self, *rows):
        self.path.write_text(json.dumps({'benchmarks': rows}))
        return report.collect([self.path])

    def test_equal_weight_whole_dfs_and_aggregate_only(self):
        a, b, c = report.FAMILIES
        median = dict(run_type='aggregate', aggregate_name='median')
        groups = self.collect(row(a, 1), row(a, 10), row(a, 10000), row(a, 999, **median),
                              row(b, 0.1, unit='us'), row(c, .001, unit='ms', **median),
                              row('balanced', 999), row('catalan', 999))
        self.assertIn('| inspected_dfs | 370 |', report.render(report.means(groups)))

    def test_per_path_before_median(self):
        a, b, c = report.FAMILIES
        def path(shape, time, batch, **kwargs):
            return row(shape, time, work='inspected_root_path', paths_per_iteration=batch, **kwargs)
        groups = self.collect(path(a, 2, 2), path(a, 400, 100), path(a, 120, 10),
                              path(b, 2e-6, 100, unit='s'), path(c, 0.8, 200, unit='us'))
        self.assertIn('| inspected_root_path | 9.33333 |', report.render(report.means(groups)))

    def test_no_mixing_across_group_dimensions(self):
        a, b, c = report.FAMILIES
        for change in [dict(n=4096), dict(budget=32768), dict(impl='dfuds'), dict(seed=7),
                       dict(work='inspected_root_path', paths_per_iteration=4)]:
            with self.subTest(change=change), self.assertRaisesRegex(ValueError, 'missing'):
                report.means(self.collect(row(a), row(b), row(c, **change), row('catalan')),
                             publish=True)

    def test_reject_invalid_data(self):
        for fields in [dict(cpu=float('nan')), dict(error_occurred=True),
                       dict(work='inspected_root_path', paths_per_iteration=0),
                       dict(work='inspected_root_path'), dict(unit='minutes')]:
            with self.subTest(fields=fields), self.assertRaises((ValueError, KeyError)):
                self.collect(row(report.FAMILIES[0], **fields))

    def test_page_preserves_outside_and_rejects_loss_of_rows(self):
        rows = [row(f) for f in report.FAMILIES]
        table = report.render(report.means(self.collect(*rows, *(row(f, budget=32768) for f in report.FAMILIES))))
        page = report.updated_page('before\n' + report.BEGIN + '\npending\n' + report.END + '\nafter\n', table)
        expected = 'before\n' + report.BEGIN + '\n\n' + table + '\n\n' + report.END + '\nafter\n'
        with self.assertRaisesRegex(ValueError, 'remove published rows'):
            report.updated_page(page, report.render(report.means(self.collect(*rows))))
        self.assertEqual(page, expected)

    def test_plot_publication_from_split_workload_shards(self):
        root = self.path.parent
        page = root / 'benchmarks/trees.md'
        page.parent.mkdir()
        pending = f'{report.BEGIN}\npending\n{report.END}'
        page.write_text(pending)
        files = []
        for n in report.SIZES:
            for target, work in [('dfs', 'inspected_dfs'), ('path', 'inspected_root_path')]:
                rows = [row(f, cpu=n * (i + 1), budget=b, impl=impl, work=work, n=n,
                            paths_per_iteration=4) for b in report.BUDGETS
                        for impl in report.IMPLEMENTATIONS for i, f in enumerate(report.FAMILIES)]
                file = root / f'tree_{target}_benchmarks-timing-{n}.json'
                file.write_text(json.dumps({'benchmarks': rows}))
                files.append(file)
        cmd = [sys.executable, '-B', report.__file__, '--page', str(page)]
        result = subprocess.run(cmd + files[:-1], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('12 missing', result.stderr)
        self.assertEqual(page.read_text(), pending)
        self.assertFalse((root / 'images').exists())
        subprocess.run(cmd + files, check=True)
        text = page.read_text()
        self.assertIn('| 0 | 1024 | 42 | bp | inspected_dfs | 2048 |', text)
        self.assertIn('| 0 | 1024 | 42 | bp | inspected_root_path | 512 |', text)
        for work, (label, unit, _) in report.WORKLOADS.items():
            file = f'tree_shape_{work}.svg'
            self.assertIn('\\image html ' + file, text)
            svg = (root / 'images/benchmarks/trees' / file).read_text()
            self.assertIn(f'{unit} per {label.lower()}', svg)
            for title in ['0 KiB budget', '32 KiB budget', '512 KiB budget', *report.IMPLEMENTATIONS.values()]:
                self.assertIn(f'>{title}<', svg)


if __name__ == '__main__':
    unittest.main()
