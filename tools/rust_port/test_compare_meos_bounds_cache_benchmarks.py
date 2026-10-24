"""Synthetic regression tests for the historical MEOS benchmark auditor."""
import copy,hashlib,json,tempfile,unittest
from pathlib import Path

from reference import ROOT
from compare_meos_bounds_cache_benchmarks import compare,digest,normalized

SCRATCH=(ROOT/'.reference').resolve()
BEFORE_REPORT=ROOT/'docs/rust/meos-branches-before-cache-benchmark.json'

class CompareMeosBoundsCacheBenchmarksTests(unittest.TestCase):
    def setUp(self):
        SCRATCH.mkdir(parents=True,exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(dir=SCRATCH)
        self.directory=Path(self.temp.name).resolve()
        self.assertTrue(self.directory.is_relative_to(SCRATCH))
        self.addCleanup(self.temp.cleanup)
        self.before=json.loads(BEFORE_REPORT.read_text(encoding='utf-8'))
        self.after=copy.deepcopy(self.before)
        self.after['native_sources_sha256']['rust/src/multiparameter_intervals.rs']='synthetic-after-hash'
        self.before_path=self._path('before.json')
        self.after_path=self._path('after.json')
        self.output_path=self._path('comparison.json')
        self._write(self.before_path,self.before)
        self._write(self.after_path,self.after)

    def _path(self,name):
        path=(self.directory/name).resolve()
        self.assertTrue(path.is_relative_to(SCRATCH),f'path escaped scratch: {path}')
        return path

    def _write(self,path,value):
        self.assertTrue(path.resolve().is_relative_to(SCRATCH))
        path.write_bytes((json.dumps(value,indent=2)+'\n').replace('\n','\r\n').encode('utf-8'))

    def _compare(self):
        return compare(self.before_path,self.after_path,self.output_path)

    def test_valid_pair_and_lf_normalized_full_report_hashes(self):
        result=self._compare()
        expected_before=hashlib.sha256(self.before_path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
        expected_after=hashlib.sha256(self.after_path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
        self.assertEqual(result['before_report_sha256'],expected_before)
        self.assertEqual(result['after_report_sha256'],expected_after)
        self.assertNotEqual(expected_before,expected_after)
        self.assertEqual(result['normalized_inputs_sha256'],digest(normalized(self.before)))
        self.assertEqual(result['normalized_inputs_sha256'],digest(normalized(self.after)))
        self.assertTrue(self.output_path.is_file())

    def test_nontrivial_cyclic_root_count_sums(self):
        for report in (self.before,self.after):
            report['numeric_cases']['1'][0]['expected'][0]=2.0
            case_count=len(report['numeric_cases']['1'])
            for mode,iterations in (('fresh',100),('prepared',1000)):
                occurrences=(iterations+case_count-1)//case_count
                report['checksums'][f'{mode}_1'][1]+=occurrences
        self._write(self.before_path,self.before);self._write(self.after_path,self.after)
        result=self._compare()
        checks=result['profiles'][0]['checksum_validation']
        self.assertEqual(checks['fresh']['root_count_sum'],113)
        self.assertEqual(checks['prepared']['root_count_sum'],1125)

    def test_iteration_counts_are_not_accepted_as_root_count_sums(self):
        for report in (self.before,self.after):
            report['numeric_cases']['1'][0]['expected'][0]=2.0
            report['checksums']['fresh_1'][1]=100
            report['checksums']['prepared_1'][1]=1000
        self._write(self.before_path,self.before);self._write(self.after_path,self.after)
        with self.assertRaises(AssertionError):self._compare()

    def test_negative_sample_is_rejected(self):
        self.after['samples_ns']['setup_1'][0]=-1
        self._write(self.after_path,self.after)
        with self.assertRaises(AssertionError):self._compare()

    def test_numeric_case_drift_is_rejected(self):
        self.after['numeric_cases']['1'][0]['expected'][1]+=1
        self._write(self.after_path,self.after)
        with self.assertRaises(AssertionError):self._compare()

    def test_checksum_corruption_is_rejected(self):
        self.after['checksums']['fresh_1'][2]+=1
        self._write(self.after_path,self.after)
        with self.assertRaises(AssertionError):self._compare()

    def test_unpermitted_native_source_change_is_rejected(self):
        self.after['native_sources_sha256']['rust/src/lib.rs']='synthetic-extra-change'
        self._write(self.after_path,self.after)
        with self.assertRaises(AssertionError):self._compare()

    def test_incorrect_stored_median_is_rejected(self):
        self.after['median_ns']['setup_1']+=1
        self._write(self.after_path,self.after)
        with self.assertRaises(AssertionError):self._compare()

if __name__=='__main__':
    unittest.main()
