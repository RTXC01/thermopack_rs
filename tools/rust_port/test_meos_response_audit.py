"""Unit, stability and tolerance regression gates for MEOS response auditing."""
import math
import unittest

from qualify_multiparameter_response import (
    FLOORS, RELATIVE, close, fit_mass, key, metadata,
    physical_reference, requests, source_request,
)


class ResponseAuditTests(unittest.TestCase):
    def test_schedule_keeps_all_models_and_quality_endpoints(self):
        qs = list(requests())
        self.assertEqual(len(qs), 854)
        self.assertEqual(len({key(q) for q in qs}), 854)
        self.assertEqual(len({q['model'] for q in qs}), 14)
        splits = [q for q in qs if q['quality'] is not None]
        self.assertEqual(len(splits), 630)
        self.assertEqual({q['quality'] for q in splits}, {0., .1, .5, .9, 1.})

    def test_si_pressure_input_is_preserved(self):
        q = dict(model='LJ_MEOS:AR', t=150., p=5355486.028939234,
                 quality=.5, phase=0)
        call, scales = source_request(q, metadata())
        self.assertEqual(call['pressure'], q['p'])
        self.assertEqual(scales['p'], 1.)
        self.assertIn(repr(q['p']), key(q))

    def test_reduced_wrapper_scaling_is_consistent_with_neon_fit(self):
        data = metadata()
        q = dict(model='LJTS_REDUCED:LJTS', t=.8, p=.02, quality=.5, phase=0)
        call, s = source_request(q, data)
        self.assertEqual((call['eos'], call['component']), ('LJTS_MEOS', 'NE'))
        neon = data['LJTS_MEOS:NE']
        self.assertAlmostEqual(s['p'] * s['v'], neon['R'] * s['t'])
        self.assertEqual(call['pressure'] / s['p'], q['p'])
        self.assertAlmostEqual(s['sound'] ** 2,
                               neon['R'] * s['t'] / fit_mass(dict(model='LJTS_MEOS:NE'), data))
        self.assertEqual(fit_mass(q, data), 1.)

    def test_small_volume_derivative_drift_cannot_hide_in_pressure_floor(self):
        for index, value in ((2, 1e-8), (3, -1e-10)):
            self.assertFalse(close(value, 2. * value, FLOORS[index], RELATIVE[index]))
        self.assertFalse(close(float('nan'), 1., 1.))
        self.assertFalse(close(float('inf'), float('inf'), 1.))

    def test_physical_sound_rejects_unstable_and_nonfinite_kernels(self):
        state = [0.] * 24
        state[12:15] = [1e5, 100., -1e9]
        state[20:22] = [30., 31.]
        good = physical_reference(300., [.001], [state], None, .04)
        self.assertIsNotNone(good)
        self.assertTrue(math.isfinite(good[4]))
        heavier = physical_reference(300., [.001], [state], None, .16)
        self.assertAlmostEqual(heavier[4] / good[4], .5)
        for index, bad in ((14, 0.), (14, 1.), (20, -1.), (13, float('nan'))):
            altered = list(state)
            altered[index] = bad
            self.assertIsNone(physical_reference(300., [.001], [altered], None, .04))


if __name__ == '__main__':
    unittest.main()

