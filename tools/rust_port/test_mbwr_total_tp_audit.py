"""Negative tests for fail-closed MBWR TP evidence checks."""
import copy
import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_mbwr_total_tp as audit


def valid_tvp_case(kind='TVP'):
    q = dict(T=300., V=.001, n=1., P=100., phase=2, kind=kind)
    values = [0.] * 45
    values[2] = -1e-3
    values[10:14] = [100., 0., -1., 0.]
    residual = [0.] * 10
    residual[2] = 1.
    base = audit.map_tvp(q, values, residual)
    direct = [base[0]]
    for start in (11, 16, 41):
        direct += [base[start]] + base[start + 2:start + 5]
    return q, {
        'CLASSICAL': [0.] * 18,
        'VALUES': values,
        'RESIDUAL': residual,
        'IDEAL': [0., 0., 18., 10., 1.] + [0.] * 10,
        'TVP': direct,
        'RAW_SOUND': [0.],
    }


class EvidenceAuditNegativeTests(unittest.TestCase):
    def singular_gap_fixture(self):
        q = dict(kind='SOURCE', record='mbwrcx21', T=299.97, P=4e6, n=1., phase=2,
                 binding='AR', profile=0, parameters=[300., 5e6, 10000., .1], V=.001)
        raw = {'DENSITY': [4197.516024649706, 2., -855919.5261020884, 1.6763806343078615e-11]}
        expected = [0.] * 22
        expected[0], expected[1] = q['P'], q['V']
        for start in (3, 13, 18):
            expected[start] = 1.
        actual = copy.deepcopy(expected)
        for start in (3, 13, 18):
            actual[start:start+4] = [0.] * 4
        return q, raw, expected, actual

    def test_singular_source_gap_is_classified_for_exact_case(self):
        q, raw, expected, actual = self.singular_gap_fixture()
        self.assertTrue(audit.singular_source_gap(q, raw, expected, actual))

    def test_source_conditioning_does_not_depend_on_record_label(self):
        q, raw, expected, actual = self.singular_gap_fixture()
        q['record'] = 'different-name'
        self.assertTrue(audit.singular_source_gap(q, raw, expected, actual))

    def test_singular_source_gap_allows_coincidentally_finite_native_derivative_flags(self):
        q, raw, expected, actual = self.singular_gap_fixture()
        for start in (3, 13, 18):
            actual[start] = 1.
        self.assertTrue(audit.singular_source_gap(q, raw, expected, actual))

    def test_singular_source_gap_requires_zero_blocks_when_native_flag_is_absent(self):
        q, raw, expected, actual = self.singular_gap_fixture()
        actual[4] = 2.
        self.assertFalse(audit.singular_source_gap(q, raw, expected, actual))

    def test_source_conditioning_rejects_resolved_pressure_slope(self):
        q, raw, _, _ = self.singular_gap_fixture()
        raw['DENSITY'][3] = 1.
        self.assertFalse(audit.unresolved_source_derivative_case(q, raw))

    def test_source_conditioning_rejects_invalid_density(self):
        q, raw, _, _ = self.singular_gap_fixture()
        raw['DENSITY'][0] = 0.
        self.assertFalse(audit.unresolved_source_derivative_case(q, raw))

    def test_singular_source_gap_rejects_closed_source_pressure(self):
        q, raw, expected, actual = self.singular_gap_fixture()
        raw['DENSITY'][2] = q['P']
        self.assertFalse(audit.singular_source_gap(q, raw, expected, actual))

    def test_source_conditioning_rejects_closed_source_pressure(self):
        q, raw, _, _ = self.singular_gap_fixture()
        raw['DENSITY'][2] = q['P']
        self.assertFalse(audit.unresolved_source_derivative_case(q, raw))

    def test_singular_source_gap_rejects_field_drift_outside_declared_blocks(self):
        q, raw, expected, actual = self.singular_gap_fixture()
        actual[0] = 1e6
        self.assertFalse(audit.singular_source_gap(q, raw, expected, actual))

    def test_omitted_case_is_rejected_by_cohort_length(self):
        with self.assertRaisesRegex(AssertionError, '9752'):
            audit.audit({'records': [], 'cases': 9751})

    def test_omitted_unqualified_claim_is_rejected(self):
        computed = [{'case': str(i)} for i in range(93)]
        with self.assertRaisesRegex(AssertionError, 'list differs'):
            audit._verify_unqualified_entries(computed[:-1], computed)

    def test_tampered_expected_cannot_be_hidden_by_matching_native_pair(self):
        q, raw = valid_tvp_case()
        expected = audit._expected(q, raw)
        expected2 = copy.deepcopy(expected)
        expected2['BASELINE'][0] += 1.
        native = copy.deepcopy(expected2)
        self.assertEqual(native, expected2)  # the tampered native/expected pair agrees exactly
        independently_derived = audit._expected(q, raw)
        with self.assertRaisesRegex(AssertionError, 'expected values are not derived'):
            audit._check(audit._same_expected(expected2, independently_derived),
                         'stored expected values are not derived from original evidence')
        self.assertEqual(native, expected2)

    def test_same_expected_ulp_scope_is_strictly_bounded(self):
        q, raw = valid_tvp_case()
        expected = audit._expected(q, raw)
        # Ensure flags are present (1.0) on a synthetic case with non-zero lnphi/Cp/sound.
        want = copy.deepcopy(expected)
        want['PHYSICAL'][41] = -0.25
        want['PHYSICAL'][42] = 1.0
        want['PHYSICAL'][46] = 40.0
        want['PHYSICAL'][47] = 1.0
        want['PHYSICAL'][48] = 50.0
        want['PHYSICAL'][49] = 1.0
        want['PHYSICAL'][50] = 200.0

        # Presence / derivative flag 1-ULP mutation must be rejected.
        for flag_idx in (42, 47, 49):
            bad = copy.deepcopy(want)
            bad['PHYSICAL'][flag_idx] = 1.0 + math.ulp(1.0)
            self.assertFalse(audit._same_expected(bad, want))

        # Unaffected slots 1-ULP mutation must be rejected.
        for view in ('BASELINE', 'PHYSICAL'):
            bad = copy.deepcopy(want)
            v = bad[view][0]
            bad[view][0] = v + math.ulp(abs(v) if v != 0.0 else 1.0)
            self.assertFalse(audit._same_expected(bad, want))

        # Absent Cp/sound placeholder mutation (0.0 -> 1e-300) when flag is 0.0 must be rejected.
        for flag_idx, val_idx in ((47, 48), (49, 50)):
            absent_want = copy.deepcopy(want)
            absent_want['PHYSICAL'][flag_idx] = 0.0
            absent_want['PHYSICAL'][val_idx] = 0.0
            bad = copy.deepcopy(absent_want)
            bad['PHYSICAL'][val_idx] = 1e-300
            self.assertFalse(audit._same_expected(bad, absent_want))

        # lnphi: 1-ULP passes, large drift (1e-12) and sign flip are rejected.
        ok = copy.deepcopy(want)
        ok['PHYSICAL'][41] = want['PHYSICAL'][41] + math.ulp(abs(want['PHYSICAL'][41]))
        self.assertTrue(audit._same_expected(ok, want))
        for bad_val in (want['PHYSICAL'][41] + 1e-12, -want['PHYSICAL'][41]):
            bad = copy.deepcopy(want)
            bad['PHYSICAL'][41] = bad_val
            self.assertFalse(audit._same_expected(bad, want))

        # Present Cp/sound: 2 ULP passes, 8 ULP is rejected.
        for val_idx in (48, 50):
            ok = copy.deepcopy(want)
            ok['PHYSICAL'][val_idx] = want['PHYSICAL'][val_idx] + 2.0 * math.ulp(abs(want['PHYSICAL'][val_idx]))
            self.assertTrue(audit._same_expected(ok, want))
            bad = copy.deepcopy(want)
            bad['PHYSICAL'][val_idx] = want['PHYSICAL'][val_idx] + 8.0 * math.ulp(abs(want['PHYSICAL'][val_idx]))
            self.assertFalse(audit._same_expected(bad, want))

    def test_same_expected_rejects_large_overflow_stored_mutations_and_overflowing_scales(self):
        q, raw = valid_tvp_case()
        want = copy.deepcopy(audit._expected(q, raw))
        want['PHYSICAL'][41] = -0.25
        want['PHYSICAL'][42] = 1.0
        want['PHYSICAL'][46] = 40.0
        want['PHYSICAL'][47] = 1.0
        want['PHYSICAL'][48] = 50.0
        want['PHYSICAL'][49] = 1.0
        want['PHYSICAL'][50] = 61.66642636762747
        scales = (0.0, 2.388105336637237)
        for val_idx in (41, 48, 50):
            bad = copy.deepcopy(want)
            bad['PHYSICAL'][val_idx] = 1e308
            self.assertFalse(audit._same_expected(bad, want, scales))
        # Overflowing sound_factor must fail closed even on a 1-ULP difference.
        one_ulp = copy.deepcopy(want)
        one_ulp['PHYSICAL'][50] = want['PHYSICAL'][50] + math.ulp(want['PHYSICAL'][50])
        self.assertFalse(audit._same_expected(one_ulp, want, (0.0, 1e308)))

    def test_truncated_original_vectors_fail_closed(self):
        with self.assertRaisesRegex(AssertionError, '45 values'):
            audit._raw_vector({'VALUES': [0.] * 44}, 'VALUES')

    def test_nonfinite_expected_is_rejected_by_native_comparison(self):
        q, _ = valid_tvp_case()
        with self.assertRaisesRegex(AssertionError, 'invalid expected'):
            audit._compare(q, {'BASELINE': [math.nan] * 51}, {'BASELINE': [math.nan] * 51})

    def test_type10_substitution_changes_physical_gold(self):
        q, raw = valid_tvp_case()
        correct = audit._expected(q, raw)['PHYSICAL']
        wrong_raw = copy.deepcopy(raw)
        wrong_raw['IDEAL'][3] = 0.
        wrong = audit._expected(q, wrong_raw)['PHYSICAL']
        self.assertNotEqual(correct, wrong)

    def test_density_phase_flag_drift_is_rejected(self):
        q = dict(T=300., V=.001, n=1., P=100., phase=2, kind='SOURCE')
        expected = {'DENSITY': [1., 2., 100., -1.]}
        actual = {'DENSITY': [1., 1., 100., -1.]}
        with self.assertRaisesRegex(AssertionError, 'phase flag drift'):
            audit._compare(q, expected, actual)

    def test_bounded_false_original_root_closure_is_rejected(self):
        q, raw = valid_tvp_case('BOUNDED')
        raw['VALUES'][10] = q['P'] + 1.
        base = audit.map_tvp(q, raw['VALUES'], raw['RESIDUAL'])
        raw['TVP'] = [base[0]]
        for start in (11, 16, 41):
            raw['TVP'] += [base[start]] + base[start + 2:start + 5]
        with self.assertRaisesRegex(AssertionError, 'does not close'):
            audit._expected(q, raw)

    def test_censored_collection_is_rejected(self):
        with self.assertRaisesRegex(AssertionError, 'censored'):
            audit._stdout_vectors({'censored': True, 'returncode': None, 'stdout': ''})


class EvidenceCorpusNegativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report, _ = audit.load_and_audit()

    def test_matching_corrupted_native_and_gold_fail_original_derivation(self):
        report = copy.deepcopy(self.report)
        report['records'][0]['expected']['BASELINE'][0] += 100.
        report['records'][0]['native']['BASELINE'][0] += 100.
        with self.assertRaisesRegex(AssertionError, 'not derived from original'):
            audit.audit(report)

    def test_corrupted_original_and_matching_vectors_fail_raw_process_binding(self):
        report = copy.deepcopy(self.report)
        row = report['records'][0]
        row['original']['CLASSICAL'][2] += 100.
        row['expected']['CLASSICAL'][2] += 100.
        row['native']['CLASSICAL'][2] += 100.
        with self.assertRaisesRegex(AssertionError, 'differs from raw stdout'):
            audit.audit(report)

    def test_conditioning_domain_cannot_hide_changed_entropy_value(self):
        report = copy.deepcopy(self.report)
        case = report['unqualified_source_derivatives'][0]['case']
        row = next(row for row in report['records'] if audit.key(row['request']) == case)
        row['native']['CLASSICAL'][7] += 100.
        with self.assertRaisesRegex(AssertionError, 'CLASSICAL mismatch'):
            audit.audit(report)

    def _clone_with_record(self, rec_idx):
        report = dict(self.report, records=list(self.report['records']))
        report['records'][rec_idx] = copy.deepcopy(self.report['records'][rec_idx])
        return report

    def test_sub_ulp_floor_classical_sign_flip_is_rejected(self):
        report = self._clone_with_record(0)
        row = report['records'][0]
        val = row['expected']['CLASSICAL'][16]
        self.assertNotEqual(val, 0.0)
        row['expected']['CLASSICAL'][16] = -val
        want = audit._expected(row['request'], row['original'])
        self.assertFalse(audit._same_expected(row['expected'], want))
        with self.assertRaisesRegex(AssertionError, 'not derived from original'):
            audit.audit(report)

    def test_corpus_physical_flag_and_value_ulp_mutations_are_rejected(self):
        tvp_index = next(
            i for i, r in enumerate(self.report['records'])
            if r['request']['kind'] == 'TVP'
            and r['expected']['PHYSICAL'][42] == 1.0
            and r['expected']['PHYSICAL'][47] == 1.0
            and r['expected']['PHYSICAL'][49] == 1.0
            and abs(r['expected']['PHYSICAL'][41]) > 1e-6
        )
        row = self.report['records'][tvp_index]
        phys = row['expected']['PHYSICAL']
        residual_fn, sound_factor = audit._physical_scales(row['request'], row['original'])
        # Flag mutation by 1 ULP is rejected (including dlnphi/dT flag 42).
        for flag_idx in (42, 47, 49):
            report = self._clone_with_record(tvp_index)
            report['records'][tvp_index]['expected']['PHYSICAL'][flag_idx] = 1.0 + math.ulp(1.0)
            with self.assertRaisesRegex(AssertionError, 'not derived from original'):
                audit.audit(report)
        # Unaffected slot 1-ULP mutation is rejected.
        for view in ('BASELINE', 'PHYSICAL'):
            report = self._clone_with_record(tvp_index)
            v = report['records'][tvp_index]['expected'][view][0]
            report['records'][tvp_index]['expected'][view][0] = v + math.ulp(abs(v) if v != 0.0 else 1.0)
            with self.assertRaisesRegex(AssertionError, 'not derived from original'):
                audit.audit(report)
        # lnphi: 1-ULP passes _same_expected, while large drift (1e-12) and sign flip fail audit.
        want = audit._expected(row['request'], row['original'])
        ok_expected = copy.deepcopy(row['expected'])
        ok_expected['PHYSICAL'][41] = want['PHYSICAL'][41] + math.ulp(max(abs(want['PHYSICAL'][41]), abs(residual_fn)))
        self.assertTrue(audit._same_expected(ok_expected, want, (residual_fn, sound_factor)))
        for bad_lnphi in (phys[41] + 1e-12, -phys[41]):
            report = self._clone_with_record(tvp_index)
            report['records'][tvp_index]['expected']['PHYSICAL'][41] = bad_lnphi
            with self.assertRaisesRegex(AssertionError, 'not derived from original'):
                audit.audit(report)
        # Present Cp/sound 8-ULP mutation (beyond operand-scaled 4 ULP) is rejected.
        cp_scale = max(abs(phys[48]), abs(phys[46]))
        sound_scale = abs(phys[50]) * sound_factor
        for val_idx, scale in ((48, cp_scale), (50, sound_scale)):
            report = self._clone_with_record(tvp_index)
            report['records'][tvp_index]['expected']['PHYSICAL'][val_idx] = phys[val_idx] + 8.0 * math.ulp(scale)
            with self.assertRaisesRegex(AssertionError, 'not derived from original'):
                audit.audit(report)

    def test_case_487_sound_speed_overflow_mutation_is_rejected(self):
        rec_idx = 486
        row = self.report['records'][rec_idx]
        phys = row['expected']['PHYSICAL']
        self.assertEqual(phys[49], 1.0)
        self.assertLess(abs(phys[50] - 61.66642636762747), 1e-9)
        _, sound_factor = audit._physical_scales(row['request'], row['original'])
        self.assertGreater(sound_factor, 2.0)
        for val_idx in (50, 48, 41):
            report = self._clone_with_record(rec_idx)
            report['records'][rec_idx]['expected']['PHYSICAL'][val_idx] = 1e308
            with self.assertRaisesRegex(AssertionError, 'not derived from original'):
                audit.audit(report)


if __name__ == '__main__':
    unittest.main()
