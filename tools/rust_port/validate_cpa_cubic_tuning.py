"""Validate CPA second derivatives/mixing, UMR/VTPR/HV-WS cubic mixing,
apparent_composition, tuning, and complexmodelinit Rust parity against live Fortran."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests/fixtures/cpa_cubic_tuning.csv"
PROBE = Path(__file__).with_name("cpa_cubic_tuning_probe.f90")


def _fmt_vec(vals):
    return " ".join(repr(float(v)) for v in vals)


def cases():
    """Yield (case_id, mode, expected_count, stdin_text) tuples."""

    def cpa_tvn(cid, comps, eos, mix, ref, variant, t, v, moles, extra=None):
        nc = len(moles)
        lines = [
            "CPA_TVN",
            variant,
            comps,
            eos,
            mix,
            ref,
            f"{nc} 0",
            f"{t!r} {v!r}",
            _fmt_vec(moles),
        ]
        if extra is not None:
            lines.append(_fmt_vec(extra))
        return (cid, "CPA_TVN", None, "\n".join(lines) + "\n")

    def cpa_tp(cid, comps, eos, mix, ref, variant, phase, t, p, moles, extra=None):
        nc = len(moles)
        lines = [
            "CPA_TP",
            variant,
            comps,
            eos,
            mix,
            ref,
            f"{nc} {phase}",
            f"{t!r} {p!r}",
            _fmt_vec(moles),
        ]
        if extra is not None:
            lines.append(_fmt_vec(extra))
        return (cid, "CPA_TP", 17 + 8 * nc + nc * nc, "\n".join(lines) + "\n")

    def cubic_tvn(cid, comps, eos, mix, ref, variant, t, v, moles, extra=None):
        nc = len(moles)
        lines = [
            "CUBIC_TVN",
            variant,
            comps,
            eos,
            mix,
            ref,
            f"{nc} 0",
            f"{t!r} {v!r}",
            _fmt_vec(moles),
        ]
        if extra is not None:
            lines.append(_fmt_vec(extra))
        return (cid, "CUBIC_TVN", 6 + 3 * nc + nc * nc, "\n".join(lines) + "\n")

    def cubic_tp(cid, comps, eos, mix, ref, variant, phase, t, p, moles, extra=None):
        nc = len(moles)
        lines = [
            "CUBIC_TP",
            variant,
            comps,
            eos,
            mix,
            ref,
            f"{nc} {phase}",
            f"{t!r} {p!r}",
            _fmt_vec(moles),
        ]
        if extra is not None:
            lines.append(_fmt_vec(extra))
        return (cid, "CUBIC_TP", 17 + 8 * nc + nc * nc, "\n".join(lines) + "\n")

    # ---- 1. CPA TVN 10-block Helmholtz + X_k ----
    yield cpa_tvn("cpa_tvn|H2O|SRK|VDW|DEFAULT|default|320|2.2e-5|1.0", "H2O", "SRK", "Classic", "DEFAULT", "default", 320.0, 2.2e-5, [1.0])
    yield cpa_tvn("cpa_tvn|H2O|SRK|VDW|DEFAULT|default|450|2.5e-3|1.4", "H2O", "SRK", "Classic", "DEFAULT", "default", 450.0, 2.5e-3, [1.4])
    yield cpa_tvn("cpa_tvn|H2O|SRK|VDW|DEFAULT|orig_rdf|350|2.5e-5|1.1", "H2O", "SRK", "Classic", "DEFAULT", "orig_rdf", 350.0, 2.5e-5, [1.1])
    yield cpa_tvn("cpa_tvn|H2O|SRK|VDW|SINTEF|default|360|3.0e-5|1.0", "H2O", "SRK", "Classic", "SINTEF", "default", 360.0, 3.0e-5, [1.0])
    yield cpa_tvn("cpa_tvn|MEOH|SRK|VDW|DEFAULT|default|330|4.5e-5|1.2", "MEOH", "SRK", "Classic", "DEFAULT", "default", 330.0, 4.5e-5, [1.2])
    yield cpa_tvn("cpa_tvn|ETOH|SRK|VDW|DEFAULT|default|360|7.5e-5|1.0", "ETOH", "SRK", "Classic", "DEFAULT", "default", 360.0, 7.5e-5, [1.0])
    yield cpa_tvn("cpa_tvn|NH3|SRK|VDW|3B-SINTEF|default|310|3.5e-5|1.0", "NH3", "SRK", "Classic", "3B-SINTEF", "default", 310.0, 3.5e-5, [1.0])
    yield cpa_tvn("cpa_tvn|NH3|SRK|VDW|4B-SINTEF|default|310|3.5e-5|1.0", "NH3", "SRK", "Classic", "4B-SINTEF", "default", 310.0, 3.5e-5, [1.0])
    yield cpa_tvn("cpa_tvn|H2O,CO2|SRK|VDW|DEFAULT|default|350|1.2e-4|0.7,0.5", "H2O,CO2", "SRK", "Classic", "DEFAULT", "default", 350.0, 1.2e-4, [0.7, 0.5])
    yield cpa_tvn("cpa_tvn|CO2,H2O|SRK|VDW|DEFAULT|default|350|1.2e-4|0.5,0.7", "CO2,H2O", "SRK", "Classic", "DEFAULT", "default", 350.0, 1.2e-4, [0.5, 0.7])
    yield cpa_tvn("cpa_tvn|H2O,ETOH|SRK|VDW|DEFAULT|default|360|6.0e-5|0.4,0.6", "H2O,ETOH", "SRK", "Classic", "DEFAULT", "default", 360.0, 6.0e-5, [0.4, 0.6])
    yield cpa_tvn("cpa_tvn|H2O,MEOH|SRK|VDW|DEFAULT|default|340|4.0e-5|0.55,0.65", "H2O,MEOH", "SRK", "Classic", "DEFAULT", "default", 340.0, 4.0e-5, [0.55, 0.65])
    yield cpa_tvn("cpa_tvn|H2O,MEOH|SRK|VDW|DEFAULT|elliot|340|4.0e-5|0.55,0.65", "H2O,MEOH", "SRK", "Classic", "DEFAULT", "elliot", 340.0, 4.0e-5, [0.55, 0.65])
    yield cpa_tvn("cpa_tvn|H2O,C1|SRK|VDW|DEFAULT|default|380|2.0e-4|0.6,0.4", "H2O,C1", "SRK", "Classic", "DEFAULT", "default", 380.0, 2.0e-4, [0.6, 0.4])
    yield cpa_tvn("cpa_tvn|H2O,C3|SRK|VDW|DEFAULT|default|370|1.5e-4|0.5,0.5", "H2O,C3", "SRK", "Classic", "DEFAULT", "default", 370.0, 1.5e-4, [0.5, 0.5])
    yield cpa_tvn("cpa_tvn|H2O,MEOH,C1|SRK|VDW|DEFAULT|default|390|2.5e-4|0.8,0.5,0.3", "H2O,MEOH,C1", "SRK", "Classic", "DEFAULT", "default", 390.0, 2.5e-4, [0.8, 0.5, 0.3])
    # CPA mixing rules: HV1, HV2, WSCPA
    yield cpa_tvn(
        "cpa_tvn|H2O,MEOH|SRK|HV1|DEFAULT|hv1_tuned|350|4.2e-5|0.5,0.5",
        "H2O,MEOH", "SRK", "Classic", "DEFAULT", "hv1_tuned", 350.0, 4.2e-5, [0.5, 0.5],
        extra=[0.25, 0.25, 180.0, -90.0, 0.12, -0.08, 0.0, 0.0],
    )
    yield cpa_tvn(
        "cpa_tvn|H2O,MEOH|SRK|HV2|DEFAULT|hv2_tuned|350|4.2e-5|0.5,0.5",
        "H2O,MEOH", "SRK", "Classic", "DEFAULT", "hv2_tuned", 350.0, 4.2e-5, [0.5, 0.5],
        extra=[0.25, 0.30, 180.0, -90.0, 0.12, -0.08, 4.0e-4, -3.0e-4],
    )
    yield cpa_tvn(
        "cpa_tvn|H2O,CO2|SRK|WS|DEFAULT|ws_tuned|380|1.8e-4|0.6,0.4",
        "H2O,CO2", "SRK", "WS", "DEFAULT", "ws_tuned", 380.0, 1.8e-4, [0.6, 0.4],
        extra=[0.3, 0.3, 0.08, 0.08, 0.22, -0.14],
    )

    # ---- 2. CPA TP full states + (T, P, n) derivatives ----
    yield cpa_tp("cpa_tp|H2O|SRK|VDW|DEFAULT|default|1|320|1.0e6|1.0", "H2O", "SRK", "Classic", "DEFAULT", "default", 1, 320.0, 1.0e6, [1.0])
    yield cpa_tp("cpa_tp|H2O|SRK|VDW|DEFAULT|default|2|450|2.0e5|1.0", "H2O", "SRK", "Classic", "DEFAULT", "default", 2, 450.0, 2.0e5, [1.0])
    yield cpa_tp("cpa_tp|MEOH|SRK|VDW|DEFAULT|default|1|310|2.0e6|1.0", "MEOH", "SRK", "Classic", "DEFAULT", "default", 1, 310.0, 2.0e6, [1.0])
    yield cpa_tp("cpa_tp|NH3|SRK|VDW|3B-SINTEF|default|2|360|1.0e6|1.0", "NH3", "SRK", "Classic", "3B-SINTEF", "default", 2, 360.0, 1.0e6, [1.0])
    yield cpa_tp("cpa_tp|H2O,CO2|SRK|VDW|DEFAULT|default|2|420|4.0e6|0.4,0.6", "H2O,CO2", "SRK", "Classic", "DEFAULT", "default", 2, 420.0, 4.0e6, [0.4, 0.6])
    yield cpa_tp("cpa_tp|H2O,ETOH|SRK|VDW|DEFAULT|default|1|330|1.5e6|0.55,0.45", "H2O,ETOH", "SRK", "Classic", "DEFAULT", "default", 1, 330.0, 1.5e6, [0.55, 0.45])
    yield cpa_tp("cpa_tp|H2O,ETOH|SRK|VDW|DEFAULT|orig_rdf|1|330|1.5e6|0.55,0.45", "H2O,ETOH", "SRK", "Classic", "DEFAULT", "orig_rdf", 1, 330.0, 1.5e6, [0.55, 0.45])
    yield cpa_tp("cpa_tp|H2O,C1|SRK|VDW|DEFAULT|default|2|400|3.0e6|0.3,0.7", "H2O,C1", "SRK", "Classic", "DEFAULT", "default", 2, 400.0, 3.0e6, [0.3, 0.7])
    yield cpa_tp("cpa_tp|H2O,MEOH,C1|SRK|VDW|DEFAULT|default|2|430|2.0e6|0.4,0.3,0.3", "H2O,MEOH,C1", "SRK", "Classic", "DEFAULT", "default", 2, 430.0, 2.0e6, [0.4, 0.3, 0.3])
    yield cpa_tp(
        "cpa_tp|H2O,CO2|SRK|WS|DEFAULT|ws_tuned|2|420|3.0e6|0.45,0.55",
        "H2O,CO2", "SRK", "WS", "DEFAULT", "ws_tuned", 2, 420.0, 3.0e6, [0.45, 0.55],
        extra=[0.3, 0.3, 0.08, 0.08, 0.22, -0.14],
    )

    # ---- 3. Cubic UMR, VTPR, init_VTPR, HVWongSandler TVN & TP ----
    yield cubic_tvn("cubic_tvn|C1,C3|PR|UMR|DEFAULT|default|310|2.0e-4|0.4,0.6", "C1,C3", "PR", "UMR", "DEFAULT", "default", 310.0, 2.0e-4, [0.4, 0.6])
    yield cubic_tvn("cubic_tvn|C1,C3,NC4|PR|UMR|DEFAULT|default|340|3.2e-4|0.5,0.4,0.3", "C1,C3,NC4", "PR", "UMR", "DEFAULT", "default", 340.0, 3.2e-4, [0.5, 0.4, 0.3])
    yield cubic_tvn("cubic_tvn|C1,C3|PR|VTPR|DEFAULT|default|310|2.0e-4|0.4,0.6", "C1,C3", "PR", "VTPR", "DEFAULT", "default", 310.0, 2.0e-4, [0.4, 0.6])
    yield cubic_tvn("cubic_tvn|C1,C3,NC4|SRK|VTPR|DEFAULT|default|340|3.2e-4|0.5,0.4,0.3", "C1,C3,NC4", "SRK", "VTPR", "DEFAULT", "default", 340.0, 3.2e-4, [0.5, 0.4, 0.3])
    yield cubic_tvn("cubic_tvn|C1,C3|PR|VTPR|DEFAULT|init_vtpr|310|2.0e-4|0.4,0.6", "C1,C3", "PR", "VTPR", "DEFAULT", "init_vtpr", 310.0, 2.0e-4, [0.4, 0.6])
    yield cubic_tvn("cubic_tvn|CO2,C1|PR|HVWONGSANDLER|DEFAULT|default|320|2.5e-4|0.45,0.55", "CO2,C1", "PR", "HVWongSandler", "DEFAULT", "default", 320.0, 2.5e-4, [0.45, 0.55])
    yield cubic_tvn(
        "cubic_tvn|CO2,C1|PR|HVWONGSANDLER|DEFAULT|hv_ws_tuned|320|2.5e-4|0.45,0.55",
        "CO2,C1", "PR", "HVWongSandler", "DEFAULT", "hv_ws_tuned", 320.0, 2.5e-4, [0.45, 0.55],
        extra=[0.28, 0.32, 220.0, -110.0, 0.15, -0.10, 3.0e-4, -2.0e-4, 0.06, 0.02],
    )
    yield cubic_tvn(
        "cubic_tvn|CO2,C1,C3|SRK|HVWONGSANDLER|DEFAULT|hv_ws_tuned|350|3.5e-4|0.4,0.5,0.3",
        "CO2,C1,C3", "SRK", "HVWongSandler", "DEFAULT", "hv_ws_tuned", 350.0, 3.5e-4, [0.4, 0.5, 0.3],
        extra=[0.20, 0.25, 160.0, -80.0, 0.10, -0.05, 2.0e-4, -1.5e-4, 0.04, 0.015],
    )

    yield cubic_tp("cubic_tp|C1,C3|PR|UMR|DEFAULT|default|2|330|3.0e6|0.4,0.6", "C1,C3", "PR", "UMR", "DEFAULT", "default", 2, 330.0, 3.0e6, [0.4, 0.6])
    yield cubic_tp("cubic_tp|C1,C3|PR|VTPR|DEFAULT|default|1|260|5.0e6|0.4,0.6", "C1,C3", "PR", "VTPR", "DEFAULT", "default", 1, 260.0, 5.0e6, [0.4, 0.6])
    yield cubic_tp("cubic_tp|C1,C3|PR|VTPR|DEFAULT|init_vtpr|2|330|3.0e6|0.4,0.6", "C1,C3", "PR", "VTPR", "DEFAULT", "init_vtpr", 2, 330.0, 3.0e6, [0.4, 0.6])
    yield cubic_tp(
        "cubic_tp|CO2,C1|PR|HVWONGSANDLER|DEFAULT|hv_ws_tuned|2|330|4.0e6|0.45,0.55",
        "CO2,C1", "PR", "HVWongSandler", "DEFAULT", "hv_ws_tuned", 2, 330.0, 4.0e6, [0.45, 0.55],
        extra=[0.28, 0.32, 220.0, -110.0, 0.15, -0.10, 3.0e-4, -2.0e-4, 0.06, 0.02],
    )

    # ---- 4. fidel_alpha ----
    for tau12, tau21 in [(0.5, 0.3), (0.0, 0.0), (-0.4, 0.9), (1.2, -0.2)]:
        cid = f"fidel|{tau12}|{tau21}"
        yield (cid, "FIDEL", 1, f"FIDEL\ndefault\n{tau12!r} {tau21!r}\n")

    # ---- 5. lowcase_bij & calcBmixCubic ----
    for cid, comps, variant, b_exp, moles in [
        ("bmix|C1,C3|init_vtpr|1.3333333333333333|0.4,0.6", "C1,C3", "init_vtpr", 4.0 / 3.0, [0.4, 0.6]),
        ("bmix|C1,C3,NC4|umr_exp|2.0|0.3,0.4,0.3", "C1,C3,NC4", "custom", 2.0, [0.3, 0.4, 0.3]),
    ]:
        nc = len(moles)
        stdin = "\n".join([
            "BMIX", variant, comps, "PR", "Classic", "DEFAULT",
            f"{nc} 0", f"300.0 {b_exp!r}", _fmt_vec(moles),
        ]) + "\n"
        yield (cid, "BMIX", 1 + nc * nc, stdin)

    # ---- 6. apparent_composition.f90 ----
    # Case A: H2O + NaCl (nc=2, nce=3, ncsym=1)
    app_a = [
        "APPARENT", "h2o_nacl", "2 3 1", "2",
        "2 2 1.0", "2 3 1.0",
        "310.0 1.0e5",
        "0.92 0.08",
        "-0.4 -1.2 -1.5",
        "0.002 -0.005 0.004",
        "-120.0 -310.0 -290.0",
        "0.8 0.2 -0.1", "0.2 1.5 0.4", "-0.1 0.4 1.7",
        "-0.25 -0.60 -0.70",
        "1.2e-3 -2.3e-3 -1.8e-3",
        "3.1e-8 -4.2e-8 5.0e-8",
        "0.4 -0.1 0.05", "-0.1 0.9 0.25", "0.05 0.25 1.1",
    ]
    yield ("apparent|h2o_nacl", "APPARENT", 3 + 7 * 2 + 2 * 4 + 1, "\n".join(app_a) + "\n")

    # Case B: H2O + NaCl + CaCl2 (nc=3, nce=4, ncsym=1, common ion Cl-)
    app_b = [
        "APPARENT", "h2o_nacl_cacl2", "3 4 1", "4",
        "2 2 1.0", "2 3 1.0", "3 4 1.0", "3 3 2.0",
        "320.0 2.5e5",
        "0.85 0.10 0.05",
        "-0.3 -1.1 -1.4 -2.2",
        "0.001 -0.004 0.003 -0.006",
        "-100.0 -280.0 -250.0 -410.0",
        "0.7 0.1 -0.05 0.02",
        "0.1 1.4 0.3 0.15",
        "-0.05 0.3 1.6 0.25",
        "0.02 0.15 0.25 2.1",
        "-0.20 -0.55 -0.65 -1.10",
        "1.0e-3 -2.0e-3 -1.5e-3 -3.2e-3",
        "2.5e-8 -3.5e-8 4.0e-8 -6.0e-8",
        "0.35 -0.08 0.04 -0.02",
        "-0.08 0.85 0.20 0.10",
        "0.04 0.20 1.05 0.18",
        "-0.02 0.10 0.18 1.40",
    ]
    yield ("apparent|h2o_nacl_cacl2", "APPARENT", 4 + 7 * 3 + 2 * 9 + 1, "\n".join(app_b) + "\n")

    # Case C: Identity (nc=2, nce=2, ncsym=2)
    app_c = [
        "APPARENT", "identity", "2 2 2", "0",
        "300.0 1.0e5",
        "0.4 0.6",
        "-0.5 -0.8",
        "0.001 -0.002",
        "-90.0 -110.0",
        "0.6 0.1", "0.1 0.9",
        "-0.15 -0.22",
        "1.1e-3 -0.9e-3",
        "2.0e-8 -1.5e-8",
        "0.3 -0.05", "-0.05 0.5",
    ]
    yield ("apparent|identity", "APPARENT", 2 + 7 * 2 + 2 * 4 + 1, "\n".join(app_c) + "\n")

    # ---- 7. tuning.f90 ----
    def tuning_case(cid, comps, eos, mix, variant, t, v, moles, extra_lines, n_get):
        nc = len(moles)
        lines = [
            "CUBIC_TVN", variant, comps, eos, mix, "DEFAULT",
            f"{nc} 0", f"{t!r} {v!r}", _fmt_vec(moles),
            *extra_lines,
        ]
        return (cid, "TUNING", n_get + 6 + 3 * nc + nc * nc, "\n".join(lines) + "\n")

    yield tuning_case(
        "tuning|hv2|CO2,C1|PR", "CO2,C1", "PR", "HV2", "tune_hv", 320.0, 2.2e-4, [0.4, 0.6],
        [_fmt_vec([0.22, 0.27, 310.0, -140.0, 0.20, -0.15, 5.0e-4, -4.0e-4])], 8,
    )
    yield tuning_case(
        "tuning|ws|CO2,C1|PR", "CO2,C1", "PR", "WS", "tune_ws", 320.0, 2.2e-4, [0.4, 0.6],
        [_fmt_vec([0.30, 0.35, 0.11, 0.11, 0.25, -0.18])], 6,
    )
    yield tuning_case(
        "tuning|twu_mc_beta|CO2,H2|PR", "CO2,H2", "PR", "Classic", "tune_twu_mc_beta", 180.0, 1.8e-4, [0.5, 0.5],
        [_fmt_vec([0.21, 0.88, 1.95, 0.72, 0.18, -0.09, 3.0696, 12.682])], 6,
    )
    yield tuning_case(
        "tuning|alpha_name|CO2,C1|PR", "CO2,C1", "PR", "Classic", "tune_alpha_name", 280.0, 2.4e-4, [0.45, 0.55],
        [_fmt_vec([0.42, 0.61])], 2,
    )
    yield tuning_case(
        "tuning|kij_lij|CO2,C1|PR", "CO2,C1", "PR", "Classic", "tune_kij_lij", 300.0, 2.0e-4, [0.4, 0.6],
        [_fmt_vec([0.12, 0.05, 0.03])], 4,
    )
    yield tuning_case(
        "tuning|vshift_quad|CO2,C1|PR", "CO2,C1", "PR", "Classic", "tune_vshift", 310.0, 2.1e-4, [0.4, 0.6],
        ["3", _fmt_vec([1.2e-6, -2.5e-9, 1.5e-12, 0.0, 0.0, 0.0])], 7,
    )
    yield tuning_case(
        "tuning|unifac_get|C1,C3|PR", "C1,C3", "PR", "UNIFAC", "tune_unifac_get", 310.0, 2.0e-4, [0.4, 0.6],
        ["1 57"], 3,
    )


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def _extract_single_row(stdout: str, tag: str) -> list[str]:
    rows = [s.split()[1:] for s in stdout.splitlines() if s.startswith(f"{tag} ")]
    if len(rows) != 1:
        raise ValueError(f"Expected exactly 1 {tag} row, got {len(rows)}: {stdout!r}")
    return rows[0]


def verify_defect_init_umr_result(returncode: int, stdout: str, stderr: str) -> None:
    combined = stdout + stderr
    if returncode == 0 or "unknown alpha corr name" not in combined:
        raise RuntimeError(
            f"Expected init_UMR 'unknown alpha corr name' abort, got rc={returncode}, out={combined!r}"
        )


def verify_defect_cpa_hv_init_result(returncode: int, stdout: str, stderr: str) -> None:
    combined = stdout + stderr
    if returncode == 0 or "Wrong GE mixing rule for CPA" not in combined:
        raise RuntimeError(
            f"Expected CPA HV init 'Wrong GE mixing rule for CPA' abort, got rc={returncode}, out={combined!r}"
        )


def verify_defect_cpa_setkijandji_result(returncode: int, stdout: str, stderr: str) -> None:
    combined = stdout + stderr
    if returncode == 0 or "Not able to set binary kij for CPA eos." not in combined:
        raise RuntimeError(
            f"Expected thermopack_setkijandji CPA abort, got rc={returncode}, out={combined!r}"
        )


def verify_defect_cpa_setlijandji_result(returncode: int, stdout: str, stderr: str) -> None:
    combined = stdout + stderr
    if returncode == 0 or "Not able to set binary lij for CPA eos." not in combined:
        raise RuntimeError(
            f"Expected thermopack_setlijandji CPA abort, got rc={returncode}, out={combined!r}"
        )


def verify_defect_cpa_setkij_single_output(stdout: str) -> None:
    tokens = _extract_single_row(stdout, "DEFECT_CPA_SETKIJ")
    if len(tokens) != 2:
        raise ValueError(f"Expected 2 DEFECT_CPA_SETKIJ values, got {tokens!r}")
    vals = [float(x) for x in tokens]
    if not all(math.isfinite(v) for v in vals):
        raise ValueError(f"Non-finite DEFECT_CPA_SETKIJ output: {vals!r}")
    if abs(vals[0] - 0.075) >= 1e-12 or abs(vals[1] - 0.075) < 1e-6:
        raise RuntimeError(f"Unexpected asymmetric CPA setkij values: {vals!r}")


def verify_defect_vtpr_bmix_output(stdout: str) -> None:
    tokens = _extract_single_row(stdout, "DEFECT_VTPR_BMIX")
    if len(tokens) != 3:
        raise ValueError(f"Expected 3 DEFECT_VTPR_BMIX values, got {tokens!r}")
    vals = [float(x) for x in tokens]
    if not all(math.isfinite(v) for v in vals):
        raise ValueError(f"Non-finite DEFECT_VTPR_BMIX output: {vals!r}")
    # calcBmixCubic uses lowcase_bij (s=4/3), whereas cbCalcBmix ignores lowcase_bij -> distinct b
    if abs(vals[0] - vals[1]) < 1e-7:
        raise RuntimeError(
            f"Expected calcBmixCubic ({vals[0]}) to differ from cbCalcBmix ({vals[1]}) in init_VTPR"
        )


def verify_defect_apparent_nan_output(stdout: str) -> None:
    tokens = _extract_single_row(stdout, "DEFECT_APPARENT_NAN")
    if len(tokens) != 3 or tokens != ["F", "T", "F"]:
        raise RuntimeError(f"Expected DEFECT_APPARENT_NAN = [F, T, F], got {tokens!r}")


def verify_defect_scheme_1ea_output(stdout: str) -> None:
    tokens = _extract_single_row(stdout, "DEFECT_1EA")
    if len(tokens) != 4:
        raise ValueError(f"Expected 4 DEFECT_1EA integers, got {tokens!r}")
    vals = [int(x) for x in tokens]
    # numAssocComps == 2 (H2O 4C + CO2 1ea), yet H2O uses sites 1..2 and numAssocSites remains 2
    # (CO2 1ea is omitted from assocIndices_bookkeeping, leaving comp_vs_sites(2,:) unassigned).
    if vals != [2, 2, 1, 2]:
        raise RuntimeError(f"Unexpected DEFECT_1EA bookkeeping flags: {vals!r}")


def verify_upstream_defects(exe: Path, env: dict) -> None:
    """Independently verify confirmed upstream Fortran defects live."""
    def run_defect(name: str):
        return subprocess.run(
            [exe], input=f"DEFECT\n{name}\n", text=True, env=env, capture_output=True, timeout=30
        )

    p = run_defect("init_umr")
    verify_defect_init_umr_result(p.returncode, p.stdout, p.stderr)

    p = run_defect("cpa_hv1_init")
    verify_defect_cpa_hv_init_result(p.returncode, p.stdout, p.stderr)

    p = run_defect("cpa_hv2_init")
    verify_defect_cpa_hv_init_result(p.returncode, p.stdout, p.stderr)

    p = run_defect("cpa_setkijandji")
    verify_defect_cpa_setkijandji_result(p.returncode, p.stdout, p.stderr)

    p = run_defect("cpa_setlijandji")
    verify_defect_cpa_setlijandji_result(p.returncode, p.stdout, p.stderr)

    p = run_defect("cpa_setkij_single")
    if p.returncode != 0:
        raise RuntimeError(("cpa_setkij_single", p.stdout, p.stderr))
    verify_defect_cpa_setkij_single_output(p.stdout)

    p = run_defect("vtpr_bmix_defect")
    if p.returncode != 0:
        raise RuntimeError(("vtpr_bmix_defect", p.stdout, p.stderr))
    verify_defect_vtpr_bmix_output(p.stdout)

    p = run_defect("apparent_nan")
    if p.returncode != 0:
        raise RuntimeError(("apparent_nan", p.stdout, p.stderr))
    verify_defect_apparent_nan_output(p.stdout)

    p = run_defect("scheme_1ea")
    if p.returncode != 0:
        raise RuntimeError(("scheme_1ea", p.stdout, p.stderr))
    verify_defect_scheme_1ea_output(p.stdout)

    print("Confirmed all 7 live upstream Fortran defects (init_umr, cpa_hv1/2_init, cpa_setkij/lijandji, vtpr_bmix, apparent_nan, scheme_1ea)")


def self_test_defect_checks() -> None:
    """Verify defect parsers reject malformed, NaN, and corrected outputs even under python -O."""

    def must_reject(fn, *args):
        try:
            fn(*args)
        except (ValueError, RuntimeError):
            return
        raise RuntimeError(f"Defect check {fn.__name__} unexpectedly accepted {args!r}")

    must_reject(verify_defect_init_umr_result, 0, "unknown alpha corr name", "")
    must_reject(verify_defect_init_umr_result, 1, "other error", "")
    verify_defect_init_umr_result(1, "", "nameToIdx::unknown alpha corr name")

    must_reject(verify_defect_cpa_hv_init_result, 0, "Wrong GE mixing rule for CPA", "")
    must_reject(verify_defect_cpa_hv_init_result, 1, "other error", "")
    verify_defect_cpa_hv_init_result(1, "", "Wrong GE mixing rule for CPA")

    must_reject(verify_defect_cpa_setkijandji_result, 0, "Not able to set binary kij for CPA eos.", "")
    must_reject(verify_defect_cpa_setkijandji_result, 1, "other", "")
    verify_defect_cpa_setkijandji_result(1, "", "Not able to set binary kij for CPA eos.")

    must_reject(verify_defect_cpa_setlijandji_result, 0, "Not able to set binary lij for CPA eos.", "")
    must_reject(verify_defect_cpa_setlijandji_result, 1, "other", "")
    verify_defect_cpa_setlijandji_result(1, "", "Not able to set binary lij for CPA eos.")

    must_reject(verify_defect_cpa_setkij_single_output, "DEFECT_CPA_SETKIJ 0.075 0.075\n")
    must_reject(verify_defect_cpa_setkij_single_output, "DEFECT_CPA_SETKIJ NaN 0.046\n")
    verify_defect_cpa_setkij_single_output("DEFECT_CPA_SETKIJ 0.075 0.04626056\n")

    must_reject(verify_defect_vtpr_bmix_output, "DEFECT_VTPR_BMIX 4.5e-5 4.5e-5 -0.2\n")
    must_reject(verify_defect_vtpr_bmix_output, "DEFECT_VTPR_BMIX NaN 4.5e-5 -0.2\n")
    verify_defect_vtpr_bmix_output("DEFECT_VTPR_BMIX 4.50e-5 4.55e-5 -0.2\n")

    must_reject(verify_defect_apparent_nan_output, "DEFECT_APPARENT_NAN F F F\n")
    must_reject(verify_defect_apparent_nan_output, "DEFECT_APPARENT_NAN F T\n")
    verify_defect_apparent_nan_output("DEFECT_APPARENT_NAN F T F\n")

    must_reject(verify_defect_scheme_1ea_output, "DEFECT_1EA 2 3 1 2\n")
    must_reject(verify_defect_scheme_1ea_output, "DEFECT_1EA 1 2 1 2\n")
    verify_defect_scheme_1ea_output("DEFECT_1EA 2 2 1 2\n")

    print("Self-test passed for all 8 defect validators")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--reference-dir", default=None)
    parser.add_argument("--write-fixture", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test_defect_checks()
        return

    ref = Path(args.reference_dir).resolve() if args.reference_dir else (ROOT / ".reference")
    compiler = Path(shutil.which(args.compiler) or args.compiler).resolve()
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join([str(ref), str(compiler.parent), env["PATH"]])
    lib = ref / ("thermopack.dll" if os.name == "nt" else "libthermopack.so")
    flags = [] if os.name == "nt" else ["-Wl,-rpath," + str(ref)]
    all_cases = list(cases())

    with tempfile.TemporaryDirectory(dir=ref) as tmp:
        exe = Path(tmp) / ("cpa_cubic_tuning_probe.exe" if os.name == "nt" else "cpa_cubic_tuning_probe")
        subprocess.run(
            [compiler, "-fdefault-real-8", "-fdefault-double-8", "-I", ref, PROBE, lib, *flags, "-o", exe],
            check=True,
            env=env,
        )
        verify_upstream_defects(exe, env)
        live = Path(tmp) / "live.csv"
        with live.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter=";", lineterminator="\n")
            writer.writerow(["case_id", "category", "tag", "values"])
            for cid, category, expected_count, stdin in all_cases:
                proc = subprocess.run([exe], input=stdin, text=True, env=env, capture_output=True, timeout=30)
                if proc.returncode:
                    raise RuntimeError((cid, proc.stdout, proc.stderr))
                matched = [
                    s.split()
                    for s in proc.stdout.splitlines()
                    if s.startswith(("TVN ", "TP ", "FIDEL ", "BMIX ", "APPARENT ", "TUNING "))
                ]
                if len(matched) != 1:
                    raise RuntimeError((cid, proc.stdout, proc.stderr))
                tag = matched[0][0]
                values = matched[0][1:]
                if expected_count is not None and len(values) != expected_count:
                    raise ValueError((cid, tag, len(values), expected_count, values))
                if not all(math.isfinite(float(v)) for v in values):
                    raise ValueError((cid, tag, "non-finite", values))
                writer.writerow([cid, category, tag, ",".join(values)])

        env["CPA_CUBIC_TUNING_FIXTURE"] = str(live)
        subprocess.run(
            ["cargo", "test", "--test", "cpa_cubic_tuning", "exact_source_cpa_cubic_tuning"],
            cwd=ROOT,
            env=env,
            check=True,
        )
        if args.write_fixture:
            FIXTURE.write_bytes(live.read_bytes())
            metadata = {
                "cases": len(all_cases),
                "fixture_sha256": sha(FIXTURE),
                "probe_sha256": sha(PROBE),
                "build": json.loads((ref / "metadata.json").read_text(encoding="utf-8")),
            }
            FIXTURE.with_suffix(".metadata.json").write_text(
                json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
            )
    print(f"{len(all_cases)} live CPA / cubic-mixing / apparent / tuning cases passed")


if __name__ == "__main__":
    main()
