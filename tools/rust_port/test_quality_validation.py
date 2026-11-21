"""Unit tests for tools/rust_port/validate_quality.py::known_baseline_difference."""
import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_quality import known_baseline_difference


def base_request(eos='PR', beta=0.4, **overrides):
    q = {
        'eos': eos,
        'components': ['C1', 'H2O'],
        'feed': [0.8, 0.2],
        'configuration': 'nrtl',
        'fixed': 'T',
        'value': 350.0,
        'beta': beta,
    }
    q.update(overrides)
    return q


class QualityValidationTests(unittest.TestCase):
    def test_accepts_all_three_documented_nrtl_pairings_within_bounds(self):
        # PR / beta=0.4: fields 6, 8 with abs diff < 1.2e-7
        q_pr_04 = base_request('PR', 0.4)
        self.assertTrue(known_baseline_difference(q_pr_04, [[6, 0.5, 0.5 + 1.0e-7]]))
        self.assertTrue(known_baseline_difference(q_pr_04, [[8, 0.2, 0.2 - 1.1e-7]]))
        self.assertTrue(
            known_baseline_difference(
                q_pr_04,
                [[6, 0.5, 0.5 + 1.0e-7], [8, 0.2, 0.2 - 1.0e-7]],
            )
        )

        # PR / beta=0.1: field 4 relative < 5e-7, fields 7, 8 absolute < 5e-7
        q_pr_01 = base_request('PR', 0.1)
        self.assertTrue(known_baseline_difference(q_pr_01, [[4, 1e-4, 1e-4 * (1.0 + 4.9e-7)]]))
        self.assertTrue(known_baseline_difference(q_pr_01, [[7, 0.3, 0.3 + 4.9e-7]]))
        self.assertTrue(known_baseline_difference(q_pr_01, [[8, 0.7, 0.7 - 4.9e-7]]))
        self.assertTrue(
            known_baseline_difference(
                q_pr_01,
                [
                    [4, 1e-4, 1e-4 * (1.0 + 4.0e-7)],
                    [7, 0.3, 0.3 + 4.0e-7],
                    [8, 0.7, 0.7 - 4.0e-7],
                ],
            )
        )

        # SRK / beta=0.4: fields 3, 4 relative < 5e-7, fields 6, 7, 8 absolute < 5e-7
        q_srk_04 = base_request('SRK', 0.4)
        self.assertTrue(known_baseline_difference(q_srk_04, [[3, 2e-5, 2e-5 * (1.0 + 4.9e-7)]]))
        self.assertTrue(known_baseline_difference(q_srk_04, [[4, 3e-4, 3e-4 * (1.0 - 4.9e-7)]]))
        self.assertTrue(known_baseline_difference(q_srk_04, [[6, 0.8, 0.8 + 4.9e-7]]))
        self.assertTrue(known_baseline_difference(q_srk_04, [[7, 0.2, 0.2 - 4.9e-7]]))
        self.assertTrue(known_baseline_difference(q_srk_04, [[8, 0.9, 0.9 + 4.9e-7]]))
        self.assertTrue(
            known_baseline_difference(
                q_srk_04,
                [
                    [3, 2e-5, 2e-5 * (1.0 + 3.0e-7)],
                    [4, 3e-4, 3e-4 * (1.0 - 3.0e-7)],
                    [6, 0.8, 0.8 + 3.0e-7],
                    [7, 0.2, 0.2 - 3.0e-7],
                    [8, 0.9, 0.9 + 3.0e-7],
                ],
            )
        )

    def test_rejects_srk_0_1_fourth_cartesian_pairing(self):
        q = base_request('SRK', 0.1)
        for idx in (3, 4, 6, 7, 8):
            self.assertFalse(known_baseline_difference(q, [[idx, 1.0, 1.0 + 1e-12]]))

    def test_rejects_alternate_feeds_for_all_three_pairings(self):
        for eos, beta, sample in (
            ('PR', 0.4, [[6, 0.5, 0.5 + 1e-8]]),
            ('PR', 0.1, [[7, 0.5, 0.5 + 1e-8]]),
            ('SRK', 0.4, [[6, 0.5, 0.5 + 1e-8]]),
        ):
            for feed in ([0.5, 0.5], [0.2, 0.8], [0.8, 0.2, 0.0]):
                q = base_request(eos, beta, feed=feed)
                self.assertFalse(known_baseline_difference(q, sample))

    def test_rejects_altered_request_fields(self):
        sample = [[6, 0.5, 0.5 + 1e-8]]
        for overrides in (
            {'eos': 'RK'},
            {'eos': 'VDW'},
            {'beta': 0.9},
            {'beta': 0.5},
            {'components': ['CO2', 'C1']},
            {'components': ['H2O', 'C1']},
            {'configuration': 'hv1'},
            {'configuration': 'hv2'},
            {'configuration': 'classic'},
            {'fixed': 'P'},
            {'fixed': 'K'},
            {'value': 351.0},
        ):
            q = base_request(**overrides)
            self.assertFalse(known_baseline_difference(q, sample))
        omitted = base_request('PR', 0.4)
        del omitted['configuration']
        self.assertFalse(known_baseline_difference(omitted, sample))

    def test_rejects_unexpected_field_indices_per_pairing(self):
        cases = (
            (base_request('PR', 0.4), {6, 8}),
            (base_request('PR', 0.1), {4, 7, 8}),
            (base_request('SRK', 0.4), {3, 4, 6, 7, 8}),
        )
        for q, allowed in cases:
            for idx in range(10):
                if idx not in allowed:
                    self.assertFalse(known_baseline_difference(q, [[idx, 1.0, 1.0 + 1e-12]]))

    def test_rejects_out_of_bound_values_per_pairing(self):
        # PR / beta=0.4: >= 1.2e-7
        q_pr_04 = base_request('PR', 0.4)
        for idx in (6, 8):
            self.assertFalse(known_baseline_difference(q_pr_04, [[idx, 0.0, 1.2e-7]]))
            self.assertFalse(known_baseline_difference(q_pr_04, [[idx, 0.5, 0.5 + 1.3e-7]]))

        # PR / beta=0.1: field 4 >= 5e-7 * max(|a|,|b|), fields 7, 8 >= 5e-7
        q_pr_01 = base_request('PR', 0.1)
        self.assertFalse(known_baseline_difference(q_pr_01, [[4, 1.0, 1.0 + 5e-7 / (1.0 - 5e-7) + 1e-12]]))
        self.assertFalse(known_baseline_difference(q_pr_01, [[4, 1e-4, 1e-4 * (1.0 + 5.1e-7)]]))
        for idx in (7, 8):
            self.assertFalse(known_baseline_difference(q_pr_01, [[idx, 0.0, 5e-7]]))
            self.assertFalse(known_baseline_difference(q_pr_01, [[idx, 0.3, 0.3 + 5.1e-7]]))

        # SRK / beta=0.4: fields 3, 4 >= 5e-7 * max(|a|,|b|), fields 6, 7, 8 >= 5e-7
        q_srk_04 = base_request('SRK', 0.4)
        for idx in (3, 4):
            self.assertFalse(known_baseline_difference(q_srk_04, [[idx, 1e-4, 1e-4 * (1.0 + 5.1e-7)]]))
        for idx in (6, 7, 8):
            self.assertFalse(known_baseline_difference(q_srk_04, [[idx, 0.0, 5e-7]]))
            self.assertFalse(known_baseline_difference(q_srk_04, [[idx, 0.5, 0.5 + 5.1e-7]]))

    def test_rejects_malformed_mixed_and_nonfinite_inputs(self):
        q_pr_04 = base_request('PR', 0.4)
        q_pr_01 = base_request('PR', 0.1)
        q_srk_04 = base_request('SRK', 0.4)

        for q in (q_pr_04, q_pr_01, q_srk_04):
            self.assertFalse(known_baseline_difference(q, []))
            self.assertFalse(known_baseline_difference(q, ['ERROR: solver failed']))

        # Mixed valid + unexpected field or out-of-bound value
        self.assertFalse(
            known_baseline_difference(q_pr_04, [[6, 0.5, 0.5 + 1e-8], [7, 0.5, 0.5 + 1e-8]])
        )
        self.assertFalse(
            known_baseline_difference(q_pr_04, [[6, 0.5, 0.5 + 1e-8], [8, 0.5, 0.5 + 2e-7]])
        )
        self.assertFalse(
            known_baseline_difference(q_pr_01, [[7, 0.5, 0.5 + 1e-8], [3, 1.0, 1.0 + 1e-12]])
        )
        self.assertFalse(
            known_baseline_difference(q_srk_04, [[6, 0.5, 0.5 + 1e-8], [4, 1e-4, 1e-4 * (1.0 + 1e-6)]])
        )

        # Non-finite values (NaN / Inf)
        for bad_val in (math.nan, math.inf, -math.inf):
            self.assertFalse(known_baseline_difference(q_pr_04, [[6, 0.5, bad_val]]))
            self.assertFalse(known_baseline_difference(q_pr_04, [[6, bad_val, 0.5]]))
            self.assertFalse(known_baseline_difference(q_pr_01, [[4, 1.0, bad_val]]))
            self.assertFalse(known_baseline_difference(q_pr_01, [[4, bad_val, bad_val]]))
            self.assertFalse(known_baseline_difference(q_srk_04, [[3, 1.0, bad_val]]))
            self.assertFalse(known_baseline_difference(q_srk_04, [[6, 0.5, bad_val]]))


if __name__ == '__main__':
    unittest.main()
