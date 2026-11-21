"""Regression gates for UV/SV partial-source qualification."""
import copy
import json
import unittest

from reference import ROOT
from validate_meos_isochoric import (
    allowed_outcomes, bounded_endpoint_gap, coordinate_bad, model_data,
    outcome_key, source_closed, source_outcome, valid_native_point,
)


class IsochoricAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = json.loads((ROOT/'tests/fixtures/multiparameter_isochoric.metadata.json').read_text(encoding='utf-8'))
        cls.native = cls.report['independent_native_kernel_checks'][0]['native']
        cls.source = next(e for e in cls.report['records'] if 'values' in e['result'])
        cls.allowed = json.loads((ROOT/'tools/rust_port/meos_isochoric_known_gaps.json').read_text(encoding='utf-8'))['outcomes']

    def test_every_native_record_is_checked(self):
        self.assertEqual(len(self.report['independent_native_kernel_checks']), self.report['attempts'])
        self.assertTrue(all(valid_native_point(e['native']) for e in self.report['independent_native_kernel_checks']))

    def test_native_raw_fractions_cannot_hide_behind_canonical_coordinates(self):
        self.assertTrue(valid_native_point(self.native))
        for index, value in ((12, -.1), (13, float('nan')), (2, 2.), (8, 3.)):
            bad = list(self.native)
            bad[index] = value
            self.assertFalse(valid_native_point(bad))
        self.assertFalse(valid_native_point(self.native[:-1]))

    def test_phase_gap_requires_endpoint_and_closed_coordinates(self):
        source = [300., 1e5, 1., 0., 1., .001, 1000., 10., 0.]
        native = [300., 1e5, 0., 0., 1., .001, 1000., 10., 1.]
        self.assertTrue(bounded_endpoint_gap(source, native, False))
        for index in (0, 1, 5, 6, 7):
            bad = list(native)
            bad[index] *= 1.01
            self.assertFalse(bounded_endpoint_gap(source, bad, False))
        for beta in (.01, .5, .99, float('nan')):
            bad = list(source)
            bad[3:5] = [beta, 1.-beta]
            self.assertFalse(bounded_endpoint_gap(bad, native, False))

    def test_unphysical_fraction_and_truncated_coordinate_arrays_fail(self):
        source = [300., 1e5, 1., -1e-10, 1.+1e-10, .001, 1000., 10., 0.]
        native = [300., 1e5, 0., 0., 1., .001, 1000., 10., 1.]
        self.assertFalse(bounded_endpoint_gap(source, native, False))
        with self.assertRaises(AssertionError):
            coordinate_bad(source, native[:-1], False)

    def test_outcome_ledger_rejects_new_cases_and_changed_failure_types(self):
        self.assertTrue(allowed_outcomes(self.report['source_outcomes'], self.allowed))
        key, value = next(iter(self.allowed.items()))
        self.assertTrue(allowed_outcomes({}, self.allowed))  # source improvements are allowed
        self.assertFalse(allowed_outcomes({key:'censored_timeout'}, self.allowed))
        self.assertFalse(allowed_outcomes({key+'|new':value}, self.allowed))

    def test_outcome_identity_preserves_solver_controls(self):
        q = self.source['request']
        for field in ('minimum_temperature', 'maximum_temperature', 'minimum_pressure', 'maximum_pressure'):
            changed = dict(q)
            changed[field] *= 1.01
            self.assertNotEqual(outcome_key(q), outcome_key(changed))
        changed = dict(q, kind='SV' if q['kind'] == 'UV' else 'UV')
        self.assertNotEqual(outcome_key(q), outcome_key(changed))

    def test_source_closure_rejects_missing_pressure_data_or_specification_drift(self):
        e = copy.deepcopy(self.source)
        self.assertTrue(source_closed(e['request'], e['result'], model_data()))
        e['result']['pressure_errors'] = []
        self.assertFalse(source_closed(e['request'], e['result'], model_data()))
        e = copy.deepcopy(self.source)
        e['result']['values'][6 if e['request']['kind'] == 'UV' else 7] += 1000.
        self.assertFalse(source_closed(e['request'], e['result'], model_data()))

    def test_source_raw_fractions_are_validated_even_when_single_phase(self):
        e = copy.deepcopy(self.source)
        e['result']['raw_result'][3:5] = [-123., 1.]
        self.assertFalse(source_closed(e['request'], e['result'], model_data()))

    def test_source_outcome_ignores_symbolized_backtrace_frames(self):
        sat_with_sv_trace = {
            'reference_error': {
                'output': [
                    ' saturation::sat did not converge',
                    '#5  0x00007f00 in __sv_solver_MOD_twophasesvsinglecomp',
                    '\tat /src/sv_solver.f90:1267',
                ]
            }
        }
        self.assertEqual(source_outcome(sat_with_sv_trace), 'flash:saturation_nonconvergence')
        sv_with_sat_trace = {
            'reference_error': {
                'output': [
                    ' sv_solver::twoPhaseSVsingleComp: No convergence',
                    '\tat /src/saturation.f90:99',
                ]
            }
        }
        self.assertEqual(source_outcome(sv_with_sat_trace), 'flash:sv_nonconvergence')

    def test_outcome_ledger_supports_closed_multi_platform_sets(self):
        declared = {
            'k1': ['flash:uv_nonconvergence', 'endpoint_phase_difference'],
            'k2': ['endpoint_phase_difference', 'flash:sv_nonconvergence'],
        }
        self.assertTrue(allowed_outcomes({'k1': 'flash:uv_nonconvergence'}, declared))
        self.assertTrue(allowed_outcomes({'k1': 'endpoint_phase_difference'}, declared))
        self.assertTrue(allowed_outcomes({'k2': 'endpoint_phase_difference'}, declared))
        self.assertTrue(allowed_outcomes({'k2': 'flash:sv_nonconvergence'}, declared))
        self.assertFalse(allowed_outcomes({'k1': 'censored_timeout'}, declared))
        self.assertFalse(allowed_outcomes({'k1': 'flash:density_iteration'}, declared))
        self.assertFalse(allowed_outcomes({'k2': 'censored_timeout'}, declared))
        self.assertFalse(allowed_outcomes({'k2': 'flash:density_iteration'}, declared))

    def test_source_outcome_trace_only_output_is_unclassified(self):
        trace_only = {
            'reference_error': {
                'output': [
                    '#5  0x1 in __sv_solver_MOD_twophasesvsinglecomp',
                    '\tat /src/sv_solver.f90:1267',
                ]
            }
        }
        self.assertEqual(source_outcome(trace_only), 'flash:unclassified_termination')
        self.assertEqual(source_outcome({'reference_error': {'output': []}}), 'flash:unclassified_termination')
        stage_trace_only = {
            'reference_error': {
                'output': [
                    'THERMOPACK_STAGE=returned_phase_properties',
                    '#0  0x1 in __sv_solver_MOD_twophasesvsinglecomp',
                    '\tat /src/sv_solver.f90:1267',
                ]
            }
        }
        self.assertEqual(source_outcome(stage_trace_only), 'phase_properties:unclassified_termination')
        recorded = self.report['source_outcomes']
        checked = 0
        for e in self.report['records']:
            key = outcome_key(e['request'])
            if key in recorded and recorded[key] != 'endpoint_phase_difference':
                self.assertEqual(source_outcome(e['result']), recorded[key])
                checked += 1
        self.assertGreater(checked, 0)


if __name__ == '__main__':
    unittest.main()
