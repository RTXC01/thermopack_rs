"""Unit tests verifying pinned coverage guards in Lee-Kesler and ExtCSP validators."""
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import validate_extcsp
import validate_lee_kesler


class FakeOracle:
    def __init__(self, _directory=None):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def call(self, request, timeout=10.0):
        return dict(values=[1.0], volume=1.0e-3)


class CspAndLeeKeslerValidatorCoverageTests(unittest.TestCase):
    def test_lee_kesler_cases_match_pinned_keys(self):
        produced = frozenset(
            (comps, ref, tuple(z), float(t), float(p), phase)
            for comps, ref, z, t, p, phase in validate_lee_kesler.cases()
        )
        self.assertEqual(len(produced), validate_lee_kesler.EXPECTED_QUALIFIED_CASE_COUNT)
        self.assertEqual(produced, validate_lee_kesler.EXPECTED_CASE_KEYS)

    def test_extcsp_cases_match_pinned_keys(self):
        produced = frozenset(
            (sh_eos, sh_mix, sh_alpha, ref_eos, ref_comp, param_ref, comps, tuple(moles), float(t), float(p), phase)
            for sh_eos, sh_mix, sh_alpha, ref_eos, ref_comp, param_ref, comps, moles, t, p, phase in validate_extcsp.cases()
        )
        self.assertEqual(len(produced), validate_extcsp.EXPECTED_QUALIFIED_CASE_COUNT)
        self.assertEqual(produced, validate_extcsp.EXPECTED_CASE_KEYS)

    def test_lee_kesler_truncated_grid_fails_before_fixture_write(self):
        truncated = list(validate_lee_kesler.cases())[:-1]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixtures = root / "tests/fixtures"
            fixtures.mkdir(parents=True)
            sentinel = "UNCHANGED_SENTINEL\n"
            csv_path = fixtures / "lee_kesler.csv"
            meta_path = fixtures / "lee_kesler.metadata.json"
            csv_path.write_text(sentinel, encoding="utf-8")
            meta_path.write_text(sentinel, encoding="utf-8")
            with (
                mock.patch.object(validate_lee_kesler, "ROOT", root),
                mock.patch.object(validate_lee_kesler, "cases", return_value=iter(truncated)),
                mock.patch.object(validate_lee_kesler, "DeadlineReferenceProcess", FakeOracle),
                mock.patch.object(validate_lee_kesler.subprocess, "run") as run_mock,
            ):
                with self.assertRaises(SystemExit) as ctx:
                    validate_lee_kesler.validate(None)
                self.assertEqual(ctx.exception.code, 1)
                run_mock.assert_not_called()
            self.assertEqual(csv_path.read_text(encoding="utf-8"), sentinel)
            self.assertEqual(meta_path.read_text(encoding="utf-8"), sentinel)

    def test_extcsp_truncated_grid_fails_before_fixture_write(self):
        truncated = list(validate_extcsp.cases())[:-1]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixtures = root / "tests/fixtures"
            fixtures.mkdir(parents=True)
            sentinel = "UNCHANGED_SENTINEL\n"
            csv_path = fixtures / "extcsp.csv"
            meta_path = fixtures / "extcsp.metadata.json"
            csv_path.write_text(sentinel, encoding="utf-8")
            meta_path.write_text(sentinel, encoding="utf-8")
            with (
                mock.patch.object(validate_extcsp, "ROOT", root),
                mock.patch.object(validate_extcsp, "cases", return_value=iter(truncated)),
                mock.patch.object(validate_extcsp, "DeadlineReferenceProcess", FakeOracle),
                mock.patch.object(validate_extcsp.subprocess, "run") as run_mock,
            ):
                with self.assertRaises(SystemExit) as ctx:
                    validate_extcsp.validate(None)
                self.assertEqual(ctx.exception.code, 1)
                run_mock.assert_not_called()
            self.assertEqual(csv_path.read_text(encoding="utf-8"), sentinel)
            self.assertEqual(meta_path.read_text(encoding="utf-8"), sentinel)


if __name__ == "__main__":
    unittest.main()
