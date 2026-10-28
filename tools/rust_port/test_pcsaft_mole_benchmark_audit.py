import copy
import base64
import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import audit_pcsaft_mole_benchmark as audit


class BenchmarkAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = audit.ROOT / 'docs/rust/pcsaft-mole-performance.json'
        if not path.exists():
            path = audit.ROOT / '.reference/pcsaft-mole-bench/measurements.json'
        cls.record = json.loads(path.read_text(encoding='utf-8'))

    def setUp(self):
        self.data = copy.deepcopy(self.record)

    def rejected(self):
        with self.assertRaises((ValueError, KeyError, TypeError)):
            audit.audit(self.data)

    def test_complete(self):
        result = audit.audit(self.data)
        self.assertEqual((result['profiles'], result['timed_calls'], result['groups']), (32, 1120, 160))

    def test_dropped_sample_with_edited_count(self):
        row = self.data['cases'][0]
        row['samples'].pop(0)
        row['summaries_seconds_per_call']['prepared_baseline']['count'] -= 1
        self.rejected()

    def test_vector_and_checksum_changed_together(self):
        x = self.data['cases'][0]['samples'][0]
        x['result']['vector'][0] += 1
        x['result']['checksum'] += self.data['settings']['iterations']
        x['stdout'] = 'RESULT=' + json.dumps(x['result']) + '\n'
        self.rejected()

    def test_parsed_result_differs_from_stdout(self):
        self.data['cases'][0]['samples'][0]['result']['seconds'] *= 2
        self.rejected()

    def test_duplicate_or_reordered_cases_and_samples(self):
        self.data['cases'][1]['case'] = copy.deepcopy(self.data['cases'][0]['case'])
        self.rejected()
        self.data = copy.deepcopy(self.record)
        x = self.data['cases'][0]['samples']; x[0], x[1] = x[1], x[0]
        self.rejected()

    def test_source_flag_and_tolerance(self):
        self.data['manifest']['source_metadata']['flags'].remove('-frecursive')
        self.rejected()
        self.data = copy.deepcopy(self.record)
        self.data['cases'][0]['preflight']['candidate']['within_tolerance'] = False
        self.rejected()

    def test_hash_mismatch(self):
        key = next(iter(self.data['hashes_after']))
        self.data['hashes_after'][key] = '0' * 64
        self.rejected()

    def test_summary_and_ratio(self):
        self.data['cases'][0]['summaries_seconds_per_call']['prepared_candidate']['median'] *= 2
        self.rejected()
        self.data = copy.deepcopy(self.record)
        self.data['cases'][0]['median_ratios']['baseline_over_candidate_fresh'] *= 2
        self.rejected()

    def test_hidden_timeout(self):
        x = self.data['cases'][0]['samples'][0]
        x['failure'] = 'timeout (process killed and reaped)'
        self.rejected()

    def test_pure_adapter_profile(self):
        self.data['cases'][0]['case']['moles'] = [0.1]
        self.rejected()

    def snapshot_copy(self, directory):
        source = audit.ROOT / 'docs/rust'
        record = directory / 'pcsaft-mole-performance.json'
        for name in (record.name, 'pcsaft-mole-performance.sources.json.gz',
                     'pcsaft-mole-performance.sources.metadata.json'):
            (directory / name).write_bytes((source / name).read_bytes())
        return record

    def test_snapshot_archive_tampering(self):
        with tempfile.TemporaryDirectory() as tmp:
            record = self.snapshot_copy(Path(tmp))
            archive = record.with_name(record.stem + '.sources.json.gz')
            metadata_path = record.with_name(record.stem + '.sources.metadata.json')
            audit.audit_snapshot(record, self.data)
            archive.write_bytes(archive.read_bytes() + b'corrupt')
            with self.assertRaises(ValueError):
                audit.audit_snapshot(record, self.data)
            files = json.loads(gzip.decompress((audit.ROOT / 'docs/rust' / archive.name).read_bytes()))
            name = 'pcsaft_mole_bench.rs'
            modified = base64.b64decode(files[name]) + b'changed'
            files[name] = base64.b64encode(modified).decode('ascii')
            payload = gzip.compress(json.dumps(files).encode())
            archive.write_bytes(payload)
            metadata = json.loads(metadata_path.read_text(encoding='utf-8'))
            metadata['archive_sha256'] = hashlib.sha256(payload).hexdigest()
            metadata['files'][name] = hashlib.sha256(modified).hexdigest()
            metadata_path.write_text(json.dumps(metadata), encoding='utf-8')
            with self.assertRaises(ValueError):
                audit.audit_snapshot(record, self.data)

    def test_snapshot_record_binding(self):
        with tempfile.TemporaryDirectory() as tmp:
            record = self.snapshot_copy(Path(tmp))
            record.write_bytes(record.read_bytes() + b' ')
            with self.assertRaises(ValueError):
                audit.audit_snapshot(record, self.data)


if __name__ == '__main__':
    unittest.main()
