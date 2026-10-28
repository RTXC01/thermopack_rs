import copy
import json
import unittest
from pathlib import Path

import capture_pcsaft_mole_derivatives as audit


class AuditTampering(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = audit.ROOT / "tests/fixtures"
        cls.metadata = json.loads((directory / "pcsaft_mole_derivatives.metadata.json").read_text(encoding="utf-8"))
        cls.archive = [json.loads(line) for line in (directory / "pcsaft_mole_derivatives.oracle.jsonl").read_text(encoding="utf-8").splitlines()]
        cls.fixture = (directory / "pcsaft_mole_derivatives.csv").read_text(encoding="utf-8")
        cls.historical = json.loads((directory / "pcsaft_state.metadata.json").read_text(encoding="utf-8"))
        cls.hashes = {key: audit.digest(audit.ROOT / "tools/rust_port" / name) for key, name in [
            ("probe_source_sha256", "pcsaft_mole_probe.f90"), ("worker_source_sha256", "pcsaft_mole_oracle.py")]}

    def setUp(self):
        self.meta = copy.deepcopy(self.metadata)
        self.rows = copy.deepcopy(self.archive)
        self.csv = self.fixture
        self.history = copy.deepcopy(self.historical)
        self.source_hashes = dict(self.hashes)

    def check(self):
        return audit.audit_data(self.meta, self.rows, self.csv, self.history, self.source_hashes)

    def rejected(self):
        with self.assertRaises(RuntimeError):
            self.check()

    def test_frozen_capture(self):
        self.assertEqual(self.check(), (346, 344, 2))

    def test_reclassified_success_with_adjusted_csv_and_counts(self):
        self.rows[0]["outcome"] = "failure"
        lines = self.csv.splitlines(keepends=True)
        del lines[1]
        self.csv = "".join(lines)
        self.meta["cases"] -= 1
        self.meta["failure_count"] += 1
        self.rejected()

    def test_altered_fixture(self):
        self.csv = self.csv.replace("-1.90699972601948327e-02", "0.00000000000000000e+00", 1)
        self.assertNotEqual(self.csv, self.fixture)
        self.rejected()

    def test_missing_or_reordered_request(self):
        self.rows.pop()
        self.meta["attempted_cases"] -= 1
        self.rejected()
        self.rows = copy.deepcopy(self.archive)
        self.rows[0], self.rows[1] = self.rows[1], self.rows[0]
        self.rejected()

    def test_malformed_values(self):
        for value in ("[NaN, 1]", "[true, 1]", "[1]", "oops"):
            self.rows = copy.deepcopy(self.archive)
            self.rows[0]["stdout"] = "RESULT=" + value + "\n"
            self.rejected()

    def test_false_failure_count(self):
        self.meta["failure_count"] += 1
        self.rejected()

    def test_censored_process(self):
        self.rows[0].update(returncode=None, outcome="timeout")
        self.rejected()

    def test_metadata_and_signature(self):
        self.meta["probe_sha256"] = "invalid"
        self.rejected()
        self.meta = copy.deepcopy(self.metadata)
        self.meta["library_sha256"] = "0" * 64
        self.meta["verified_library_sha256"] = "0" * 64
        self.rejected()
        self.meta = copy.deepcopy(self.metadata)
        self.source_hashes["worker_source_sha256"] = "0" * 64
        self.rejected()


if __name__ == "__main__":
    unittest.main()
