"""Regression gates for qualified MEOS flash benchmark inputs and checksums."""
import contextlib
import csv
import json
import io
import math
import unittest

from validate_meos_flash import key as flash_case_key

from benchmark_meos_flash import (
    endpoint_coordinates, endpoint_request, expected_checksum, parse_args, parse_rows,
    valid_result, verify_checksum, PROFILES, ROOT,
)


class BenchmarkAuditTests(unittest.TestCase):
    def test_endpoint_classification_uses_original_request_not_equal_volume_slots(self):
        self.assertFalse(endpoint_request(dict(kind='PH')))
        self.assertFalse(endpoint_request(dict(kind='PS', quality=.5)))
        self.assertFalse(endpoint_request(dict(kind='TP', quality=0.)))
        self.assertTrue(endpoint_request(dict(kind='PH', quality=0.)))
        self.assertTrue(endpoint_request(dict(kind='PS', kernel_quality=1.)))

    def test_committed_cohort_distinguishes_endpoint_seeds_from_ordinary_single_phase(self):
        meta = json.loads((ROOT/'tests/fixtures/meos_flash.metadata.json').read_text(encoding='utf-8'))
        requests = {flash_case_key(e['request'],e['result']):e['request'] for e in meta['records'] if e['result'].get('values')}
        endpoints = ordinary = 0
        with (ROOT/'tests/fixtures/meos_flash.csv').open(encoding='utf-8') as handle:
            for row in csv.DictReader(handle,delimiter=';'):
                case = row['case'].split('|')
                if case[1] not in PROFILES or not case[5].startswith(('FLASH:PH:', 'FLASH:PS:')):
                    continue
                endpoint = endpoint_request(requests[row['case']])
                endpoints += endpoint
                ordinary += not endpoint and float(row['liquid_or_overall_volume']) == float(row['vapor_or_overall_volume'])
        self.assertEqual(endpoints,220)
        self.assertEqual(ordinary,192)

    def test_endpoint_projection_uses_bulk_volume_at_tiny_quality(self):
        source = [374.172579, 4053799.082591912, 1., 0., .0001823236455205708, 0., 39192.317466, 157.873038]
        native = [374.172579, 4053799.082591912, 0., 1.3224346108521192e-12,
                  .00018232364550715846, .00022043244508877658, 39192.317466, 157.873038]
        a, b = endpoint_coordinates(source), endpoint_coordinates(native)
        self.assertEqual(a[2:4], b[2:4])
        self.assertLess(abs(a[4]-b[4]), 2e-13)
        self.assertEqual(b[5], 0.)

    def test_vapor_endpoint_uses_weighted_bulk_instead_of_liquid_slot(self):
        q = [300., 1e5, 0., 1.-3e-8, .001, .002, 1000., 10.]
        result = endpoint_coordinates(q)
        self.assertAlmostEqual(result[4], (1.-q[3])*q[4]+q[3]*q[5])
        self.assertGreater(result[4], .0019)
        q[3] = .5
        with self.assertRaises(AssertionError):
            endpoint_coordinates(q)

    def test_unphysical_raw_phase_data_is_rejected_before_canonicalization(self):
        q = [300., 1e5, 0., .5, .001, .002, 1000., 10.]
        self.assertTrue(valid_result(q))
        for index, value in ((2, 3.), (3, -.1), (3, float('nan')), (4, -.001), (5, .0005)):
            bad = list(q)
            bad[index] = value
            self.assertFalse(valid_result(bad))
        self.assertFalse(valid_result(q[:-1]))

    def test_uneven_cycles_check_actual_visit_counts(self):
        checks = {'a':[11.]*8, 'b':[23.]*8, 'c':[37.]*8}
        expected = expected_checksum(checks, list(checks), 10)
        self.assertEqual(expected, [4.*11.+3.*23.+3.*37.]*8)
        verify_checksum(expected, expected)
        with self.assertRaises(AssertionError):
            verify_checksum(expected_checksum(checks, list(checks), 9), expected)

    def test_checksum_checks_all_fields_and_summation_cancellation(self):
        wanted = [0.]*8
        # Small cancellation residuals can arise when summing large terms.
        verify_checksum([1e-8]*8, wanted, [1e7]*8)
        for field in range(8):
            bad = list(wanted)
            bad[field] = 1.
            with self.assertRaises(AssertionError):
                verify_checksum(bad, wanted, [1e7]*8)

    def test_output_parser_rejects_duplicate_nonfinite_or_missing_fields(self):
        good = 'check_0_0_0,'+','.join(['1']*8)
        self.assertEqual(len(parse_rows(good, 'check_', 8)), 1)
        for bad in (good+'\n'+good, good.replace(',1', ',nan', 1), good.rsplit(',', 1)[0]):
            with self.assertRaises(AssertionError):
                parse_rows(bad, 'check_', 8)

    def test_repeats_and_modes_are_validated_before_building(self):
        for args in (['--compiler','unused','--repeats','0'],
                     ['--compiler','unused','--repeats','-1'],
                     ['--compiler','unused','--preflight-only','--measure']):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                parse_args(args)
        self.assertEqual(parse_args(['--compiler','unused','--repeats','3']).repeats, 3)


if __name__ == '__main__':
    unittest.main()
