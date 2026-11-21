"""Adversarial checks for the offline MBWR coexistence evidence audit."""
import copy
import gzip
import json
import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_mbwr_coexistence as audit
import generate_data
import inventory


def raw_process(state=(300.,1e6,.001,.02,1e5), check=None):
    if check is None:
        check=[1e6,1e6,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.]
    fmt=lambda x: f'{x:.17e}'
    return {'returncode':0,'stdout':'INIT_RETURNED\nPOINT_STATUS=0 TP='+ ' '.join(map(fmt,state[:2]))+'\nSTATE='+' '.join(map(fmt,state))+'\nCHECK='+' '.join(map(fmt,check))+'\n','stderr':''}


class ParseNegativeTests(unittest.TestCase):
    def test_censored_source_process_rejected(self):
        with self.assertRaisesRegex(AssertionError,'censored'):
            audit._source_binding({'censored':True,'returncode':None,'stdout':''})

    def test_failed_source_process_is_retained_as_initialization_failure(self):
        source={'returncode':1,'stdout':'','stderr':'original init failure'}
        parsed=audit._source_binding(source)
        self.assertEqual(audit.source_category(source,parsed),'initialization_failed')

    def test_truncated_state_and_nonfinite_state_fail_source_classification(self):
        p=raw_process(); p['stdout']=p['stdout'].replace('STATE=','STATE=1.0 ',1)
        with self.assertRaisesRegex(AssertionError,'5 values'):
            audit._source_binding(p)
        p=raw_process(); p['stdout']=p['stdout'].replace('STATE=3.00000000000000000e+02','STATE=nan',1)
        parsed=audit._source_binding(p)
        self.assertEqual(audit.source_category(p,parsed),'nonfinite_state')

    def test_native_state_cannot_be_changed_without_changing_stdout(self):
        row={'state':[300.,1e6,.001,.02,1e5],'roots':[10.,50.],
             'stdout':'3.00000000000000000e+02,1.00000000000000000e+06,1.00000000000000002e-03,2.00000000000000004e-02,1.00000000000000000e+05;1.00000000000000000e+01,5.00000000000000000e+01'}
        audit._native_observation(row,'test')
        row['state'][0]=301.
        with self.assertRaisesRegex(AssertionError,'differs from stdout'):
            audit._native_observation(row,'test')

    def test_original_root_query_must_match_recomputed_request_and_stdout(self):
        q={'record':'mbwrcx1','critical':[300.,5e6,10000.]}
        state=[250.,1e6,.001,.02,1e5]
        query='|'.join(map(str,['mbwrcx1','R',state[0],4*q['critical'][2],state[1],0,0,0]))
        proof={'request':query,'stdout':'ROOTS;1.00000000000000000e+01,5.00000000000000000e+01','roots':[10.,50.]}
        self.assertEqual(audit._native_root_query(proof,q,state),[10.,50.])
        proof['roots']=[10.,51.]
        with self.assertRaisesRegex(AssertionError,'raw root-query stdout'):
            audit._native_root_query(proof,q,state)

    def test_cohort_inputs_rejects_modified_compiled_components_under_historical_hash(self):
        current = dict(audit.source_input_hashes())
        recorded = dict(current)
        recorded['tools/rust_port/generate_data.py'] = audit.HISTORICAL_GENERATE_DATA_SHA256
        current['tools/rust_port/generate_data.py'] = 'f' * 64
        self.assertFalse(audit._cohort_inputs_match(recorded, current, compiled_sha='0' * 64))
        self.assertTrue(
            audit._cohort_inputs_match(
                recorded,
                current,
                compiled_sha=audit.HISTORICAL_COMPILED_COMPONENTS_SHA256,
            )
        )

    def test_generate_data_and_inventory_sorting_is_os_independent(self):
        fluid_files = inventory.inventory()['fluid_files']
        self.assertLess(fluid_files.index('n-Butane.json'), fluid_files.index('Neon.json'))
        self.assertLess(fluid_files.index('pseudo.json'), fluid_files.index('R11.json'))
        generated = generate_data.generate()
        self.assertLess(generated.index('// fluids/n-Butane.json'), generated.index('// fluids/Neon.json'))
        self.assertLess(generated.index('// fluids/pseudo.json'), generated.index('// fluids/R11.json'))

class FullReportNegativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        candidate=Path(__file__).resolve().parents[2]/'.reference/mbwr_coexistence_report.json'
        fixture=Path(__file__).resolve().parents[2]/'tests/fixtures/mbwr_coexistence.audit.json.gz'
        if candidate.exists(): cls.report=json.loads(candidate.read_text(encoding='utf-8'))
        elif fixture.exists(): cls.report=json.loads(gzip.decompress(fixture.read_bytes()))
        else: raise unittest.SkipTest('full coexistence evidence is not yet available')

    def mutated(self, fn):
        report=copy.deepcopy(self.report); fn(report); return report

    def test_baseline_full_report_audits(self):
        audit.audit(self.report)

    def test_changed_native_state_rejected_even_if_stored_classification_is_unchanged(self):
        r=self.mutated(lambda x:x['records'].__getitem__(390)['native']['state'].__setitem__(0,999.))
        with self.assertRaises(AssertionError): audit.audit(r)

    def test_changed_original_state_rejected_against_raw_stdout(self):
        def change(r):
            row=r['records'][390]; row['original']['STATE'][0]+=1.
        with self.assertRaisesRegex(AssertionError,'raw stdout'):
            audit.audit(self.mutated(change))

    def test_paired_source_state_and_original_edit_is_rejected_by_kernel_trace(self):
        def change(r):
            row=r['records'][390]
            row['original']['STATE'][2]*=1.01
            source=row['source']; old=row['original']['STATE'][2]
            # Rewrite only the STATE vector in the raw process, preserving the
            # paired decoded value so raw-to-decoded equality alone cannot pass.
            fields=source['stdout'].splitlines()
            for i,line in enumerate(fields):
                if line.startswith('STATE='):
                    values=line[len('STATE='):].split(); values[2]=f'{old:.17e}'
                    fields[i]='STATE='+' '.join(values)
            source['stdout']='\n'.join(fields)+'\n'
        with self.assertRaisesRegex(AssertionError,'kernel|original source|kernel STATE|classification'):
            audit.audit(self.mutated(change))

    def test_paired_source_check_and_raw_stdout_edit_is_rejected_by_kernel_trace(self):
        def change(r):
            row=r['records'][390]; row['original']['CHECK'][0]+=1000.
            fields=row['source']['stdout'].splitlines()
            for i,line in enumerate(fields):
                if line.startswith('CHECK='):
                    values=line[len('CHECK='):].split(); values[0]=f'{float(values[0])+1000.:.17e}'
                    fields[i]='CHECK='+' '.join(values)
            row['source']['stdout']='\n'.join(fields)+'\n'
        with self.assertRaisesRegex(AssertionError,'kernel CHECK|pressure closure|source classification'):
            audit.audit(self.mutated(change))

    def test_paired_internal_metadata_and_stdout_edit_is_rejected_against_request(self):
        def change(r):
            row=next(x for x in r['records'] if x.get('source') and x['request']['mode']>=3 and x.get('original',{}).get('META'))
            row['original']['META'][4]+=1.
            fields=row['source']['stdout'].splitlines()
            for i,line in enumerate(fields):
                if line.startswith('META='):
                    values=line[len('META='):].split(); values[4]=f'{float(values[4])+1.:.17e}'
                    fields[i]='META='+' '.join(values)
            row['source']['stdout']='\n'.join(fields)+'\n'
        with self.assertRaisesRegex(AssertionError,'META parameter'):
            audit.audit(self.mutated(change))

    def test_paired_native_state_and_stdout_edit_is_rejected_by_kernel_trace(self):
        def change(r):
            row=r['records'][390]
            row['native']['state'][2]*=1.01
            body,roots=row['native']['stdout'].split(';')
            values=body.split(','); values[2]=f'{row["native"]["state"][2]:.17e}'
            row['native']['stdout']=','.join(values)+';'+roots
        with self.assertRaisesRegex(AssertionError,'native kernel request|original kernel|native differs'):
            audit.audit(self.mutated(change))

    def test_changed_root_query_stdout_rejected(self):
        def change(r):
            row=next(x for x in r['records'] if x.get('root_query'))
            row['root_query']['stdout']='ROOTS;1.2345'
        with self.assertRaisesRegex(AssertionError,'stdout'):
            audit.audit(self.mutated(change))

    def test_omitting_a_root_from_raw_query_and_decoded_list_is_rejected(self):
        def change(r):
            row=next(x for x in r['records'] if x.get('root_query') and isinstance(x['root_query'].get('roots'),list) and len(x['root_query']['roots'])>1)
            proof=row['root_query']; proof['roots'].pop()
            proof['stdout']='ROOTS;'+','.join(f'{x:.17e}' for x in proof['roots'])
        with self.assertRaises(AssertionError):
            audit.audit(self.mutated(change))

    def test_native_stdout_is_bound_to_the_decoded_state(self):
        def change(r):
            row=r['records'][0]
            if 'state' in row['native']:
                row['native']['stdout']='0,0,0,0,0;'
            else:
                row['native']['stdout']='ERROR: tampered'
        with self.assertRaises(AssertionError):
            audit.audit(self.mutated(change))

    def test_missing_kernel_check_rejected(self):
        def change(r):
            row=next(x for x in r['records'] if x.get('reference_kernel')); row['reference_kernel']='absent-kernel'
        with self.assertRaisesRegex(AssertionError,'missing original kernel'):
            audit.audit(self.mutated(change))

    def test_missing_authoritative_case_rejected(self):
        def change(r): r['records'].pop()
        with self.assertRaisesRegex(AssertionError,'case count'):
            audit.audit(self.mutated(change))

    def test_unmatched_native_error_is_not_parity_success(self):
        def change(r):
            row=r['records'][390]; row['native']={'error':'ERROR: test','stdout':'ERROR: test'}; row['native_kernel']=None
        with self.assertRaises(AssertionError): audit.audit(self.mutated(change))

    def test_source_classification_is_recomputed(self):
        def change(r):
            row=next(x for x in r['records'] if x.get('source') and x['classification']!='original_solver_match')
            row['classification']='original_solver_match'
        with self.assertRaises(AssertionError): audit.audit(self.mutated(change))

    def test_physical_classification_is_recomputed(self):
        def change(r): r['records'][0]['classification']='native_error_unqualified'
        with self.assertRaisesRegex(AssertionError,'physical classification'):
            audit.audit(self.mutated(change))

    def test_csv_gold_is_rejected_for_nonmatching_original_outcome(self):
        rows=[]
        for e in self.report['records']:
            native=e['native']; observed=','.join(map(repr,native['state'])) if 'state' in native else native['error']
            original=e.get('original',{}).get('STATE') if e['classification']=='original_solver_match' else ''
            gold=original if isinstance(original,str) else ','.join(map(repr,original)) if original else ''
            rows.append(dict(case=e['request']['id'],request=audit.native_key(e['request']),original_solver=gold,
                             native_observed=observed,classification=e['classification']))
        bad=next(i for i,e in enumerate(self.report['records']) if e.get('source') and e['classification']!='original_solver_match')
        rows[bad]['original_solver']='1,2,3,4,5'
        with self.assertRaisesRegex(AssertionError,'source gold'):
            audit._check_fixture_rows(rows,self.report)

    def test_cohort_inputs_rejects_tampered_generate_data_hash(self):
        def change(r):
            r['signature']['cohort_inputs_sha256']['tools/rust_port/generate_data.py'] = '0' * 64
        with self.assertRaisesRegex(AssertionError, 'cohort input hashes'):
            audit.audit(self.mutated(change))

    def test_cohort_inputs_rejects_tampered_other_cohort_hash(self):
        def change(r):
            r['signature']['cohort_inputs_sha256']['tools/rust_port/generate_mbwr_data.py'] = '0' * 64
        with self.assertRaisesRegex(AssertionError, 'cohort input hashes'):
            audit.audit(self.mutated(change))


if __name__=='__main__': unittest.main()
