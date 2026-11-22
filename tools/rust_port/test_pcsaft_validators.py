"""Unit tests verifying pinned coverage guards in the PC-SAFT validator."""
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import validate_pcsaft


class FakeOracle:
    def __init__(self, _directory=None):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def call(self, request, timeout=15.0):
        raise AssertionError("oracle should not be called for a truncated grid")


class PcSaftValidatorCoverageTests(unittest.TestCase):
    def test_cases_match_pinned_keys(self):
        produced = frozenset(
            (kind, param_ref, comps, tuple(moles), float(t), float(p), phase)
            for kind, param_ref, comps, moles, t, p, phase in validate_pcsaft.cases()
        )
        self.assertEqual(len(produced), validate_pcsaft.EXPECTED_QUALIFIED_CASE_COUNT)
        self.assertEqual(produced, validate_pcsaft.EXPECTED_CASE_KEYS)

    def test_rdf_cases_match_pinned_keys(self):
        produced = frozenset(
            (comps, tuple(moles), float(t), float(v))
            for comps, moles, t, v in validate_pcsaft.rdf_cases()
        )
        self.assertEqual(len(produced), validate_pcsaft.EXPECTED_RDF_CASE_COUNT)
        self.assertEqual(produced, validate_pcsaft.EXPECTED_RDF_CASE_KEYS)

    def test_truncated_grid_fails_before_oracle_and_fixture_write(self):
        truncated = list(validate_pcsaft.cases())[:-1]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixtures = root / "tests/fixtures"
            fixtures.mkdir(parents=True)
            sentinel = "UNCHANGED_SENTINEL\n"
            csv_path = fixtures / "pcsaft_full.csv"
            meta_path = fixtures / "pcsaft_full.metadata.json"
            csv_path.write_text(sentinel, encoding="utf-8")
            meta_path.write_text(sentinel, encoding="utf-8")
            with (
                mock.patch.object(validate_pcsaft, "ROOT", root),
                mock.patch.object(validate_pcsaft, "cases", return_value=iter(truncated)),
                mock.patch.object(validate_pcsaft, "DeadlineReferenceProcess", FakeOracle),
                mock.patch.object(validate_pcsaft.subprocess, "run") as run_mock,
            ):
                with self.assertRaises(SystemExit) as ctx:
                    validate_pcsaft.validate(None)
                self.assertEqual(ctx.exception.code, 1)
                run_mock.assert_not_called()
            self.assertEqual(csv_path.read_text(encoding="utf-8"), sentinel)
            self.assertEqual(meta_path.read_text(encoding="utf-8"), sentinel)


if __name__ == "__main__":
    unittest.main()
