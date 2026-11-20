"""Mutated evidence must fail even when native and expected vectors agree."""
import copy
import math
import unittest
from audit_mbwr_caloric import _same_physical,audit,load_and_audit
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

    def _clone_with_record(self,rec_idx):
        bad=dict(self.report,records=list(self.report['records']))
        bad['records'][rec_idx]=copy.deepcopy(self.report['records'][rec_idx])
        return bad

    def test_presence_flag_single_ulp_mutation_is_rejected(self):
        rec_idx=next(
            i for i,e in enumerate(self.report['records'])
            if e['physical_from_original_jet'][39]==1.0
            and e['physical_from_original_jet'][41]==1.0
            and e['physical_from_original_jet'][43]==1.0
        )
        for flag_idx in (39,41,43):
            bad=self._clone_with_record(rec_idx)
            bad['records'][rec_idx]['physical_from_original_jet'][flag_idx]=1.0+math.ulp(1.0)
            with self.assertRaises(AssertionError):audit(bad)

    def test_unaffected_field_single_ulp_mutation_is_rejected(self):
        for idx in (10,38):
            bad=self._clone_with_record(0)
            v=bad['records'][0]['physical_from_original_jet'][idx]
            bad['records'][0]['physical_from_original_jet'][idx]=v+math.ulp(abs(v) if v!=0.0 else 1.0)
            with self.assertRaises(AssertionError):audit(bad)

    def test_present_cp_and_sound_exceeding_4_ulp_is_rejected(self):
        rec_idx=next(
            i for i,e in enumerate(self.report['records'])
            if e['physical_from_original_jet'][39]==1.0 and e['physical_from_original_jet'][43]==1.0
        )
        rec=self.report['records'][rec_idx]
        p=rec['physical_from_original_jet']
        q=rec['request']
        cp_scale=max(abs(p[40]),abs(p[38]))
        denom=abs(p[12]-q['T']*(p[11]*p[11])/p[38])
        factor=max(1.0,abs(p[12])/denom) if denom>0.0 else 1.0
        for idx,scale in ((40,cp_scale),(44,abs(p[44])*factor)):
            bad=self._clone_with_record(rec_idx)
            bad['records'][rec_idx]['physical_from_original_jet'][idx]=p[idx]+8.0*math.ulp(scale)
            with self.assertRaises(AssertionError):audit(bad)

    def test_lnphi_drift_and_sign_flip_bounds(self):
        rec_idx=next(
            i for i,e in enumerate(self.report['records'])
            if e['physical_from_original_jet'][41]==1.0 and abs(e['physical_from_original_jet'][42])>1e-6
        )
        rec=self.report['records'][rec_idx]
        p=rec['physical_from_original_jet']
        residual_fn=rec['source_baseline'][3]-rec['source_ideal'][8]
        # 1-ULP lnphi drift passes.
        ok=p[:]
        ok[42]=p[42]+math.ulp(max(abs(p[42]),abs(residual_fn)))
        self.assertTrue(_same_physical(p,ok,rec['request'],residual_fn))
        # Large lnphi drift (1e-12) is rejected.
        bad=self._clone_with_record(rec_idx)
        bad['records'][rec_idx]['physical_from_original_jet'][42]=p[42]+1e-12
        with self.assertRaises(AssertionError):audit(bad)
        # Sign flip on lnphi is rejected.
        bad=self._clone_with_record(rec_idx)
        bad['records'][rec_idx]['physical_from_original_jet'][42]=-p[42]
        with self.assertRaises(AssertionError):audit(bad)


if __name__=='__main__': unittest.main()
