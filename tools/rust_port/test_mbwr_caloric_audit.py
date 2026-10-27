"""Mutated evidence must fail even when native and expected vectors agree."""
import copy
import unittest
from audit_mbwr_caloric import audit,load_and_audit
from qualify_mbwr_caloric import differences,floors


class CaloricAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.report,cls.result=load_and_audit()

    def test_complete_original_evidence_and_dimensions(self):
        self.assertEqual(self.result['cases'],5994)
        self.assertEqual(self.result['native_views'],11988)
        self.assertEqual(len(floors(self.report['records'][0]['request'])),45)

    def test_truncation_nonfinite_and_changed_property_cannot_pass(self):
        e=self.report['records'][0];q=e['request'];a=e['source_baseline']
        for bad in (a[:-1],[float('nan')]*45,[x*1.01 for x in a]):
            self.assertTrue(differences(q,a,bad))
        for bad in (a[:-1],[float('nan')]*45,[float('inf')]*45):
            self.assertTrue(differences(q,bad,a))

    def test_changed_original_vector_even_with_matching_native_is_rejected(self):
        bad=copy.deepcopy(self.report)
        bad['records'][0]['source_baseline'][14]+=100
        bad['records'][0]['native'][0][14]+=100
        with self.assertRaises(AssertionError):audit(bad)

    def test_type10_correction_cannot_be_replaced_by_native_gold(self):
        bad=copy.deepcopy(self.report)
        e=next(e for e in bad['records'] if e['source_ideal'][3]==10.)
        e['physical_from_original_jet']=e['native'][1]=e['source_baseline'][:]
        with self.assertRaises(AssertionError):audit(bad)

    def test_missing_request_or_wrong_failure_claim_is_rejected(self):
        bad=copy.deepcopy(self.report);bad['records'].pop()
        with self.assertRaises(AssertionError):audit(bad)
        bad=copy.deepcopy(self.report)
        e=next(e for e in bad['public_constructor_calls'] if 'PUBLIC_RETURNED' not in e['stdout'])
        e['expected_initialization_status']=0
        with self.assertRaises(AssertionError):audit(bad)


if __name__=='__main__': unittest.main()
