"""Numerical and schema checks for the rank/select plotting script."""

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plot_rank_select as plots


def row(name, time, unit="ns", aggregate=None):
    result = {"name": name, "real_time": time, "cpu_time": time * 2, "time_unit": unit}
    if aggregate:
        result["run_name"] = name
        result["name"] += "_" + aggregate
        result["run_type"] = "aggregate"
        result["aggregate_name"] = aggregate
    return result


class RankSelectPlotTest(unittest.TestCase):
    def test_names_include_all_operations_fnbp_and_backend_underscores(self):
        for operation in plots.OPERATIONS:
            raw = "build_both" if operation == "build" else operation
            for family in plots.FAMILIES:
                suffix = f"fnbp_{raw}" if family == "fnbp" else f"{raw}_{family}"
                sample = row(f"rank_select_succinct_fast_{suffix}/N:17179869184/min_time:1.000", 1)
                self.assertEqual(plots.parse_name(sample),
                                 (operation, family, "rank_select_succinct_fast", 1 << 34))

    def test_aggregates_preferred_per_point_without_counting_iterations_twice(self):
        baseline = "rank_select_rank1_50/N:1024"
        compressed = "rank_select_compressed_rank1_50/N:1024"
        rows = [row(baseline, 100), row(baseline, 200), row(baseline, 7, aggregate="median"),
                row(baseline, 999, aggregate="mean"), row(baseline, 5, aggregate="stddev"),
                row(baseline, 0.3, aggregate="cv"), row(compressed, 8), row(compressed, 12)]
        points, warnings = plots.collect_file({"benchmarks": rows})
        self.assertEqual(points[("rank1", "50", "rank_select", 1024)], 7)
        self.assertEqual(points[("rank1", "50", "rank_select_compressed", 1024)], 10)
        self.assertFalse(warnings)

    def test_time_units_metric_and_build_csv(self):
        rows = [row("rank_select_rank0_50/N:1024", 0.025, "us"),
                row("rank_select_select0_50/N:1024", 0.000035, "ms"),
                row("rank_select_fnbp_build_both/N:1024", 0.003, "s")]
        points, _ = plots.collect_file(rows, metric="cpu")
        self.assertAlmostEqual(points[("rank0", "50", "rank_select", 1024)], 50)
        self.assertAlmostEqual(points[("select0", "50", "rank_select", 1024)], 70)
        self.assertEqual(points[("build", "fnbp", "rank_select", 1024)], 6_000_000)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            plots.write_csv(points, output, "cpu", "median")
            with (output / "data.csv").open() as stream:
                build = next(r for r in csv.DictReader(stream) if r["operation"] == "build")
            self.assertEqual(build["unit"], "ms")
            self.assertEqual(float(build["value"]), 6)

    def test_auxiliary_percent_excludes_source_copy_and_prefers_build(self):
        query = row("rank_select_orzo_rank1_50/N:1024", 10)
        query.update(aux_bytes=10, input_bytes=128, source_copy_bytes=128)
        build = row("rank_select_orzo_build_both_50/N:1024", 2, "ms")
        build.update(aux_bytes=8, input_bytes=128, backend_owned_bytes=136,
                     source_copy_bytes=128)
        points, warnings = plots.collect_file([query, build], quantity="auxiliary")
        self.assertEqual(points, {("auxiliary", "50", "rank_select_orzo", 1024): 6.25})
        self.assertFalse(warnings)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            plots.write_csv(points, output, "wall", "median")
            with (output / "data.csv").open() as stream:
                sample = next(csv.DictReader(stream))
            self.assertEqual(sample["unit"], "%")
            self.assertEqual(sample["metric"], "space")
            self.assertEqual(float(sample["value"]), 6.25)

    def test_auxiliary_query_fallback_and_aggregate_selection(self):
        name = "rank_select_rank1_50/N:1024"
        a, b = row(name, 1), row(name, 2, aggregate="median")
        a["aux_bits_per_input_bit"] = 0.04
        b["aux_bits_per_input_bit"] = 0.03
        deviation = row(name, 2, aggregate="stddev")
        deviation.update(aux_bytes=0, input_bytes=0)
        points, warnings = plots.collect_file([a, b, deviation], quantity="auxiliary")
        self.assertEqual(points, {("auxiliary", "50", "rank_select", 1024): 3})
        self.assertFalse(warnings)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.json"
            path.write_text(json.dumps({"benchmarks": [a, b]}))
            loaded = plots.load_points([path], "wall", "median", include_auxiliary=True)
            self.assertEqual(loaded[("auxiliary", "50", "rank_select", 1024)], 3)
            self.assertEqual(loaded[("rank1", "50", "rank_select", 1024)], 2)

    def test_auxiliary_missing_invalid_and_zero_counters(self):
        name = "rank_select_build_both_50/N:1024"
        missing, invalid, negative, zero = [row(name, 1) for _ in range(4)]
        invalid.update(aux_bytes=8, input_bytes=0)
        negative["aux_bits_per_input_bit"] = -0.01
        zero.update(aux_bytes=0, input_bytes=128)
        points, warnings = plots.collect_file([missing, invalid, negative, zero],
                                             quantity="auxiliary")
        self.assertEqual(points, {("auxiliary", "50", "rank_select", 1024): 0})
        self.assertEqual(len(warnings), 2)

    def test_legacy_aggregate_names_and_counter_size(self):
        name = "rank_select_orzo_select1_sparse/N:4096/min_time:1.000_median"
        sample = row(name, 11)
        self.assertEqual(plots.collect_file([sample])[0],
                         {("select1", "sparse", "rank_select_orzo", 4096): 11})
        sample = row("rank_select_rank1_50_median", 12)
        sample["N"] = 8192
        self.assertEqual(plots.collect_file([sample])[0],
                         {("rank1", "50", "rank_select", 8192): 12})

    def test_independent_files_receive_equal_weight(self):
        name = "rank_select_rank1_50/N:1024"
        with tempfile.TemporaryDirectory() as directory:
            a, b = Path(directory) / "a.json", Path(directory) / "b.json"
            a.write_text(json.dumps({"benchmarks": [row(name, 1) for _ in range(5)]}))
            b.write_text(json.dumps({"benchmarks": [row(name, 9, aggregate="median")]}))
            self.assertEqual(plots.load_points([a, b], "wall", "median"),
                             {("rank1", "50", "rank_select", 1024): 5})

    def test_failed_invalid_and_unrelated_rows_do_not_create_points(self):
        failed = {"name": "rank_select_orzo_select1_50/N:1024", "error_occurred": True,
                  "error_message": "input has no one bits"}
        rows = [failed, row("rank_select_rank1_50/N:1024", float("nan")),
                row("rank_select_rank1_50/N:1024", -1),
                row("rank_select_rank1_50/N:1024", 10, "minutes"),
                row("unrelated/N:1024", 10), row("rank_select_rank0_50/N:1024", 6)]
        points, warnings = plots.collect_file(rows)
        self.assertEqual(points, {("rank0", "50", "rank_select", 1024): 6})
        self.assertEqual(len(warnings), 4)
        self.assertFalse(any(key[0] == "select0" for key in points))

    def test_mean_option_uses_mean_aggregate_and_raw_mean(self):
        name = "rank_select_select1_50/N:1024"
        points, _ = plots.collect_file([row(name, 1), row(name, 10), row(name, 100)], aggregate="mean")
        self.assertEqual(next(iter(points.values())), 37)
        points, _ = plots.collect_file([row(name, 8, aggregate="mean"), row(name, 99)], aggregate="mean")
        self.assertEqual(next(iter(points.values())), 8)

    def test_dry_runs_and_malformed_documents_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "dry-run"):
            plots.collect_file({"context": {"dry_run": True}, "benchmarks": []})
        with self.assertRaisesRegex(ValueError, "expected"):
            plots.collect_file({"benchmarks": "not a list"})


if __name__ == "__main__":
    unittest.main()
