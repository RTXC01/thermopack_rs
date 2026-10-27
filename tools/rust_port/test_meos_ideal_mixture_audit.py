"""Negative gates for internal-kernel qualification, independent of Fortran."""
import copy
import unittest
from qualify_meos_ideal_mixture import fields, source_closed, differences, requests


class MixtureAuditTests(unittest.TestCase):
    def setUp(self):
        self.q = dict(ids=[0,1],x=[.3,.7],T=450.,P=1e5,phase=2)
        self.result = dict(returncode=0,values=[1.]*len(fields(2)),checks=[1e5,-1e6,50.,2.]*2)
        self.result['values'][-2:] = [.04,.04]

    def test_complete_schedule_and_vector_dimensions(self):
        schedule = list(requests())
        self.assertEqual(len(schedule),734)
        self.assertEqual(len({str(q) for q in schedule}),734)
        self.assertTrue(source_closed(self.q,self.result))
        for n in (1,2,3,5):
            self.assertEqual(len(fields(n)),13+8*n+n*n)

    def test_unclosed_or_unstable_source_cannot_become_fixture(self):
        for slot,bad in ((0,2e5),(1,1.),(2,-1.),(3,9.)):
            result = copy.deepcopy(self.result);result['checks'][slot] = bad
            self.assertFalse(source_closed(self.q,result))
        for bad in (float('nan'),float('inf')):
            result = copy.deepcopy(self.result);result['values'][0] = bad
            self.assertFalse(source_closed(self.q,result))
        for field in ('checks','values'):
            result = copy.deepcopy(self.result);result[field].pop()
            self.assertFalse(source_closed(self.q,result))
        self.assertFalse(source_closed(self.q,dict(censored=True)))
        self.assertFalse(source_closed(self.q,dict(self.result,returncode=1)))

    def test_coordinate_drift_nonfinite_and_truncation_are_rejected(self):
        expected = self.result['values']
        self.assertEqual(differences(expected,expected,2),[])
        for bad in (expected[:-1], [float('nan')]*len(expected), [x*1.01 for x in expected]):
            self.assertTrue(differences(expected,bad,2))


if __name__=='__main__':
    unittest.main()
