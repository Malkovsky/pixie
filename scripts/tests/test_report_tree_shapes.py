import json
from statistics import mean
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import report_tree_shapes as report
import plot_tree_curves as matched


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
                             strict=True)

    def test_reject_invalid_data(self):
        for fields in [dict(cpu=float('nan')), dict(error_occurred=True),
                       dict(work='inspected_root_path', paths_per_iteration=0),
                       dict(work='inspected_root_path'), dict(unit='minutes')]:
            with self.subTest(fields=fields), self.assertRaises((ValueError, KeyError)):
                self.collect(row(report.FAMILIES[0], **fields))

    def test_mixed_csr_endpoints(self):
        rows = [dict(run_name=f'tree_path/bfs_csr/{op}/{family}/N:1024/seed:42',
                     run_type='aggregate', aggregate_name='median', aggregate_unit='time',
                     budget_usage_bytes=0, budget_limit_bytes=0, threads=1, repetitions=5,
                     pairs_per_batch=256, nodes_per_path=length, max_nodes_per_path=length,
                     cpu_time=time, time_unit='ns')
                for family in matched.FAMILIES for op, length, time in [('random', 7, 10), ('sibling', 3, 99)]]
        self.path.write_text(json.dumps(dict(context=dict(date='test', executable='bench', library_version='test'),
                                             benchmarks=rows)))
        with patch.multiple(matched, SUCCINCT=(), REFERENCES=('bfs_csr',), SIZES=(1024,)):
            data = matched.export([], [], None, self.path)
        point, = data['charts']['tree_path']
        self.assertEqual(point['families'], [[10, 1, 0, 7]] * 3)
        self.assertNotIn('tree_path_sibling', data['charts'])

    def test_saved_matched_matrix(self):
        script = matched.OUTPUT.read_text()
        data = json.loads(script.removeprefix('globalThis.PIXIE_TREE_MATCHED = ').rstrip(';\n'))
        self.assertEqual(set(data['charts']), {'inspected_dfs', 'inspected_root_path', 'tree_path'})
        self.assertEqual(set(data['metadata']['normalization']), set(data['charts']))
        roles = {s['role'] for s in data['metadata']['sources']}
        self.assertEqual(roles, {'csr', 'references', 'inspected', 'random'})
        for op, points in data['charts'].items():
            expected = {(i, b, n) for i in matched.SUCCINCT for b in matched.BUDGETS for n in matched.SIZES}
            expected.update((i, 0, n) for i in matched.REFERENCES for n in matched.SIZES)
            self.assertEqual(len(points), len(expected))
            self.assertEqual({(p['impl'], p['b'], p['n']) for p in points}, expected)
            lengths = {}
            for p in points:
                samples = p['families']
                self.assertEqual(len(samples), 3)
                self.assertEqual(p['t'], mean(s[0] for s in samples))
                if op != 'inspected_dfs':
                    self.assertEqual(p['length'], mean(s[3] for s in samples))
                    per_family = [s[3] for s in samples]
                    self.assertEqual(lengths.setdefault(p['n'], per_family), per_family)
                divisor = {'inspected_root_path': 4096, 'inspected_dfs': p['n']}.get(op, 1)
                self.assertTrue(all(s[1] == divisor for s in samples))


if __name__ == '__main__':
    unittest.main()
