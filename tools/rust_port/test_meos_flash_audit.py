"""Regression gates for partial-surface MEOS flash qualification."""
import copy
import json
import unittest
from reference import ROOT
from validate_meos_flash import bounded_endpoint_gap,gap_key,specified_pressure,valid_native_point


class AuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        report=json.loads((ROOT/'tests/fixtures/meos_flash.metadata.json').read_text())
        cls.gap=report['native_differences'][0]

    def mutate_native(self,index,value):
        case=copy.deepcopy(self.gap)
        values=list(map(float,case['native'].split(',')));values[index]=value
        case['native']=','.join(map(repr,values));return case

    def test_declared_trace_endpoint_difference_is_bounded(self):
        self.assertTrue(bounded_endpoint_gap(self.gap))

    def test_temperature_or_caloric_drift_is_rejected(self):
        values=list(map(float,self.gap['native'].split(',')))
        for index in (0,1,6,7):
            with self.subTest(index=index):
                self.assertFalse(bounded_endpoint_gap(self.mutate_native(index,values[index]*1.01)))

    def test_bulk_volume_drift_is_rejected(self):
        value=float(self.gap['native'].split(',')[4])
        self.assertFalse(bounded_endpoint_gap(self.mutate_native(4,2.*value)))

    def test_material_fraction_drift_is_rejected(self):
        self.assertFalse(bounded_endpoint_gap(self.mutate_native(8,1.)))

    def test_solver_error_is_never_a_permitted_gap(self):
        case=copy.deepcopy(self.gap);case['native']='ERROR failed'
        self.assertFalse(bounded_endpoint_gap(case))

    def test_case_identity_preserves_grid_dimensions(self):
        q=self.gap['request'];original=gap_key(q)
        for field in ('pressure','minimum_temperature','maximum_temperature','temperature'):
            with self.subTest(field=field):
                changed=dict(q);changed[field]*=1.01;self.assertNotEqual(gap_key(changed),original)
        changed=dict(q,target=q['target']+1e-8);self.assertEqual(gap_key(changed),original)
        self.assertNotEqual(gap_key(dict(q,guess_factor=1.3)),original)

    def test_ordinary_pressure_specification_uses_input(self):
        self.assertEqual(specified_pressure(dict(pressure=100.),dict(values=[300.,200.])),100.)

    def test_quality_pressure_uses_recorded_seed(self):
        self.assertEqual(specified_pressure(dict(pressure=100.,quality=.5),dict(seed=dict(pressure=20.),result_pressure=30.)),20.)

    def test_invalid_raw_fraction_cannot_hide_behind_canonical_endpoint(self):
        point=[300.,1e5,0.,0.,1.,1.,100.,10.,1.,0.,1.,2.]
        self.assertTrue(valid_native_point(point))
        for beta in (-1e-12,1.+1e-12,float('nan')):
            altered=list(point);altered[9]=beta
            self.assertFalse(valid_native_point(altered))

    def test_phase_fraction_inconsistency_is_rejected(self):
        point=[300.,1e5,0.,0.,1.,1.,100.,10.,1.,.5,1.,2.]
        self.assertFalse(valid_native_point(point))
        point[2]=1.;point[3]=.5
        self.assertTrue(valid_native_point(point))


if __name__=='__main__':unittest.main()
