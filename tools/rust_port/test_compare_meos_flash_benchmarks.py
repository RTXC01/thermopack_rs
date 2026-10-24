"""Negative tests for the committed, offline MEOS timing evidence audit."""
import copy
import json
import unittest

from benchmark_meos_flash import ROOT
from compare_meos_flash_benchmarks import audit


class StoredTimingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = json.loads((ROOT/'docs/rust/meos-flash-baseline.json').read_text(encoding='utf-8'))

    def test_complete_recorded_evidence_passes(self):
        self.assertEqual(audit(self.report), dict(rows=1484,groups=36,measured_samples=5,faster_median_groups=30))

    def test_changed_median_or_sample_is_rejected(self):
        for field in ('median_ns','timed_samples_ns'):
            bad = copy.deepcopy(self.report)
            key = next(iter(bad[field]))
            if field == 'median_ns':bad[field][key] *= 1.1
            else:bad[field][key][0] *= 1.1
            with self.assertRaises(AssertionError):audit(bad)

    def test_wrong_operation_count_cannot_pass_matching_relative_speed(self):
        bad = copy.deepcopy(self.report)
        for label in ('source','native'):
            group = next(iter(bad['raw_checksums_by_sample'][1][label]))
            values = bad['raw_checksums_by_sample'][1][label][group]
            values[1:] = [v*.999 for v in values[1:]]
        with self.assertRaises(AssertionError):audit(bad)

    def test_changed_manifest_or_omitted_group_is_rejected(self):
        bad = copy.deepcopy(self.report);bad['numeric_input_rows'][0][4] *= 1.01
        with self.assertRaises(AssertionError):audit(bad)
        bad = copy.deepcopy(self.report);bad['preflight']['group_counts'].pop(next(iter(bad['preflight']['group_counts'])))
        with self.assertRaises(AssertionError):audit(bad)

    def test_missing_kernel_checks_and_invalid_provenance_are_rejected(self):
        bad = copy.deepcopy(self.report);bad['preflight']['independent_kernel_checks'].pop()
        with self.assertRaises(AssertionError):audit(bad)
        bad = copy.deepcopy(self.report);bad['runtime_inputs_sha256'][next(iter(bad['runtime_inputs_sha256']))] = 'invalid'
        with self.assertRaises(AssertionError):audit(bad)


if __name__ == '__main__':
    unittest.main()
