"""Exact-source parity qualification for PC-SAFT, sPC-SAFT, PCP-SAFT, and sPCP-SAFT
(including association, multipolar DD/QQ/DQ, and simplified RDF kernels).

Note: Each PCSAFT_FULL row evaluates the 10-block TV residual Helmholtz kernel
(`query_pcsaft_full_fres`) on raw extensive `moles`, and the TP property/derivative
columns on normalized unit-mole `x = moles / sum(moles)` matching ThermoPack's
TP wrapper convention. Direct raw-amount non-unit-mole TP scaling is verified in
`tests/pc_saft_assoc_polar.rs`.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from deadline_reference import DeadlineReferenceProcess
from reference import ROOT

R_GAS = 8.31446261815324

# Independently pinned qualification grid count and exact case identities so a
# silent truncation of cases() or an unexpected reference drop cannot pass.
EXPECTED_QUALIFIED_CASE_COUNT = 35
EXPECTED_CASE_KEYS = frozenset(
    [
        # Non-associating pure and multicomponent (PC-SAFT and sPC-SAFT)
        ("PC-SAFT", "Default", "C3", (1.0,), 300.0, 5.0e5, "L"),
        ("PC-SAFT", "Default", "C3", (1.0,), 300.0, 5.0e5, "V"),
        ("PC-SAFT", "Default", "CO2,C1", (0.3, 0.7), 280.0, 3.0e6, "V"),
        ("sPC-SAFT", "Default", "CO2,C1", (0.3, 0.7), 280.0, 3.0e6, "V"),
        ("PC-SAFT", "Default", "C1,C2,C3", (0.5, 0.3, 0.2), 240.0, 6.0e6, "L"),
        ("sPC-SAFT", "Default", "C1,C2,C3", (1.2, 0.5, 0.3), 320.0, 4.0e6, "V"),
        # Associating pure fluids (4C water, 2B methanol/ethanol, 4C ammonia)
        ("PC-SAFT", "Default", "H2O", (1.0,), 350.0, 1.0e6, "L"),
        ("PC-SAFT", "Default", "H2O", (1.0,), 500.0, 5.0e5, "V"),
        ("sPC-SAFT", "Default", "H2O", (1.0,), 350.0, 1.0e6, "L"),
        ("PC-SAFT", "Default", "MEOH", (1.0,), 320.0, 1.0e6, "L"),
        ("PC-SAFT", "Default", "MEOH", (1.0,), 420.0, 5.0e5, "V"),
        ("PC-SAFT", "Default", "ETOH", (1.0,), 330.0, 1.0e6, "L"),
        ("PC-SAFT", "Default", "NH3", (1.0,), 280.0, 2.0e6, "L"),
        # Associating mixtures (self + cross association, Default + ariComb/geoComb)
        ("PC-SAFT", "Default", "H2O,MEOH", (0.6, 0.4), 340.0, 2.0e6, "L"),
        ("PC-SAFT", "Default", "H2O,MEOH", (0.6, 0.4), 460.0, 5.0e5, "V"),
        ("sPC-SAFT", "Default", "H2O,MEOH", (0.6, 0.4), 340.0, 2.0e6, "L"),
        ("PC-SAFT", "Default", "H2O,ETOH", (0.5, 0.5), 350.0, 1.5e6, "L"),
        ("PC-SAFT", "Default", "NH3,H2O", (0.4, 0.6), 310.0, 2.0e6, "L"),
        ("sPC-SAFT", "Default", "NH3,H2O", (0.4, 0.6), 310.0, 2.0e6, "L"),
        # Associating + inert mixtures (including non-unit total moles)
        ("PC-SAFT", "Default", "H2O,C1", (0.7, 0.3), 400.0, 8.0e6, "L"),
        ("PC-SAFT", "Default", "H2O,CO2", (0.8, 0.2), 360.0, 5.0e6, "L"),
        ("sPC-SAFT", "Default", "H2O,MEOH,C1", (0.8, 0.6, 0.6), 360.0, 6.0e6, "L"),
        # Polar PCP-SAFT & sPCP-SAFT (Q-Q, D-D, D-Q, and polar + associating)
        ("PCP-SAFT", "Default", "CO2", (1.0,), 260.0, 5.0e6, "L"),
        ("PCP-SAFT", "Default", "CO2", (1.0,), 330.0, 4.0e6, "V"),
        ("sPCP-SAFT", "Default", "CO2", (1.0,), 260.0, 5.0e6, "L"),
        ("PCP-SAFT", "Default", "ACETONE", (1.0,), 320.0, 1.0e6, "L"),
        ("PCP-SAFT", "Default", "ACETONE", (1.0,), 420.0, 5.0e5, "V"),
        ("PCP-SAFT", "Default", "H2S", (1.0,), 270.0, 3.0e6, "L"),
        ("PCP-SAFT", "Default", "CO2,ACETONE", (0.4, 0.6), 320.0, 3.0e6, "L"),
        ("sPCP-SAFT", "Default", "CO2,ACETONE", (0.4, 0.6), 320.0, 3.0e6, "L"),
        # Note: H2S is placed first for PCP-SAFT because upstream Fortran
        # `src/saft_association.f90:136,147` omits `eosPCP_SAFT` and reuses
        # `g(1, 1)` for all pairs; placing H2S first makes `g(1, 1)` match the
        # physical `g(H2S, H2S)`. sPCP-SAFT uses pair-independent `g(eta)`.
        ("PCP-SAFT", "Default", "H2S,CO2", (0.5, 0.5), 290.0, 4.0e6, "L"),
        ("sPCP-SAFT", "Default", "CO2,H2S", (0.5, 0.5), 290.0, 4.0e6, "L"),
        ("PCP-SAFT", "Default", "CO2,C1,ACETONE", (0.3, 0.3, 0.4), 330.0, 5.0e6, "L"),
        # Polar + associating (pure PCP-SAFT and mixture sPCP-SAFT)
        ("PCP-SAFT", "Default", "H2O", (1.0,), 360.0, 2.0e6, "L"),
        ("sPCP-SAFT", "Default", "NH3,H2O", (0.4, 0.6), 320.0, 3.0e6, "L"),
    ]
)

EXPECTED_RDF_CASE_COUNT = 3
EXPECTED_RDF_CASE_KEYS = frozenset(
    [
        ("CO2,C1", (0.3, 0.7), 280.0, 2.0e-3),
        ("H2O,MEOH", (0.6, 0.4), 340.0, 3.5e-5),
        ("C1,C2,C3", (0.5, 0.3, 0.2), 260.0, 8.0e-5),
    ]
)


def case_key(request):
    return (
        request["kind"],
        request["parameter_reference"],
        request["components"],
        tuple(request["moles"]),
        float(request["temperature"]),
        float(request["pressure"]),
        request["phase"],
    )


def cases():
    grid = [
        ("PC-SAFT", "Default", "C3", [1.0], [(300.0, 5.0e5, "L"), (300.0, 5.0e5, "V")]),
        ("PC-SAFT", "Default", "CO2,C1", [0.3, 0.7], [(280.0, 3.0e6, "V")]),
        ("sPC-SAFT", "Default", "CO2,C1", [0.3, 0.7], [(280.0, 3.0e6, "V")]),
        ("PC-SAFT", "Default", "C1,C2,C3", [0.5, 0.3, 0.2], [(240.0, 6.0e6, "L")]),
        ("sPC-SAFT", "Default", "C1,C2,C3", [1.2, 0.5, 0.3], [(320.0, 4.0e6, "V")]),
        ("PC-SAFT", "Default", "H2O", [1.0], [(350.0, 1.0e6, "L"), (500.0, 5.0e5, "V")]),
        ("sPC-SAFT", "Default", "H2O", [1.0], [(350.0, 1.0e6, "L")]),
        ("PC-SAFT", "Default", "MEOH", [1.0], [(320.0, 1.0e6, "L"), (420.0, 5.0e5, "V")]),
        ("PC-SAFT", "Default", "ETOH", [1.0], [(330.0, 1.0e6, "L")]),
        ("PC-SAFT", "Default", "NH3", [1.0], [(280.0, 2.0e6, "L")]),
        ("PC-SAFT", "Default", "H2O,MEOH", [0.6, 0.4], [(340.0, 2.0e6, "L"), (460.0, 5.0e5, "V")]),
        ("sPC-SAFT", "Default", "H2O,MEOH", [0.6, 0.4], [(340.0, 2.0e6, "L")]),
        ("PC-SAFT", "Default", "H2O,ETOH", [0.5, 0.5], [(350.0, 1.5e6, "L")]),
        ("PC-SAFT", "Default", "NH3,H2O", [0.4, 0.6], [(310.0, 2.0e6, "L")]),
        ("sPC-SAFT", "Default", "NH3,H2O", [0.4, 0.6], [(310.0, 2.0e6, "L")]),
        ("PC-SAFT", "Default", "H2O,C1", [0.7, 0.3], [(400.0, 8.0e6, "L")]),
        ("PC-SAFT", "Default", "H2O,CO2", [0.8, 0.2], [(360.0, 5.0e6, "L")]),
        ("sPC-SAFT", "Default", "H2O,MEOH,C1", [0.8, 0.6, 0.6], [(360.0, 6.0e6, "L")]),
        ("PCP-SAFT", "Default", "CO2", [1.0], [(260.0, 5.0e6, "L"), (330.0, 4.0e6, "V")]),
        ("sPCP-SAFT", "Default", "CO2", [1.0], [(260.0, 5.0e6, "L")]),
        ("PCP-SAFT", "Default", "ACETONE", [1.0], [(320.0, 1.0e6, "L"), (420.0, 5.0e5, "V")]),
        ("PCP-SAFT", "Default", "H2S", [1.0], [(270.0, 3.0e6, "L")]),
        ("PCP-SAFT", "Default", "CO2,ACETONE", [0.4, 0.6], [(320.0, 3.0e6, "L")]),
        ("sPCP-SAFT", "Default", "CO2,ACETONE", [0.4, 0.6], [(320.0, 3.0e6, "L")]),
        ("PCP-SAFT", "Default", "H2S,CO2", [0.5, 0.5], [(290.0, 4.0e6, "L")]),
        ("sPCP-SAFT", "Default", "CO2,H2S", [0.5, 0.5], [(290.0, 4.0e6, "L")]),
        ("PCP-SAFT", "Default", "CO2,C1,ACETONE", [0.3, 0.3, 0.4], [(330.0, 5.0e6, "L")]),
        ("PCP-SAFT", "Default", "H2O", [1.0], [(360.0, 2.0e6, "L")]),
        ("sPCP-SAFT", "Default", "NH3,H2O", [0.4, 0.6], [(320.0, 3.0e6, "L")]),
    ]
    for kind, param_ref, comps, moles, states in grid:
        for t, p, phase in states:
            yield kind, param_ref, comps, moles, t, p, phase


def rdf_cases():
    yield "CO2,C1", [0.3, 0.7], 280.0, 2.0e-3
    yield "H2O,MEOH", [0.6, 0.4], 340.0, 3.5e-5
    yield "C1,C2,C3", [0.5, 0.3, 0.2], 260.0, 8.0e-5


def verify_case_grid_before_run(case_list, rdf_list):
    keys = frozenset(
        (kind, param_ref, comps, tuple(moles), float(t), float(p), phase)
        for kind, param_ref, comps, moles, t, p, phase in case_list
    )
    missing = sorted(EXPECTED_CASE_KEYS - keys)
    extra = sorted(keys - EXPECTED_CASE_KEYS)
    rdf_keys = frozenset(
        (comps, tuple(moles), float(t), float(v))
        for comps, moles, t, v in rdf_list
    )
    rdf_missing = sorted(EXPECTED_RDF_CASE_KEYS - rdf_keys)
    rdf_extra = sorted(rdf_keys - EXPECTED_RDF_CASE_KEYS)
    if (
        len(case_list) != EXPECTED_QUALIFIED_CASE_COUNT
        or missing
        or extra
        or len(rdf_list) != EXPECTED_RDF_CASE_COUNT
        or rdf_missing
        or rdf_extra
    ):
        raise ValueError(
            json.dumps(
                dict(
                    expected=EXPECTED_QUALIFIED_CASE_COUNT,
                    planned=len(case_list),
                    missing_cases=len(missing),
                    extra_cases=len(extra),
                    expected_rdf=EXPECTED_RDF_CASE_COUNT,
                    planned_rdf=len(rdf_list),
                    missing_rdf=len(rdf_missing),
                    extra_rdf=len(rdf_extra),
                ),
                indent=2,
            )
        )


def verify_expected_coverage(good, failed, censored, unqualified):
    actual_keys = frozenset(case_key(r["request"]) for r in good)
    missing = sorted(EXPECTED_CASE_KEYS - actual_keys)
    extra = sorted(actual_keys - EXPECTED_CASE_KEYS)
    if (
        failed
        or censored
        or unqualified
        or len(good) != EXPECTED_QUALIFIED_CASE_COUNT
        or missing
        or extra
    ):
        raise ValueError(
            json.dumps(
                dict(
                    expected=EXPECTED_QUALIFIED_CASE_COUNT,
                    qualified=len(good),
                    failed=len(failed),
                    censored=len(censored),
                    unqualified=len(unqualified),
                    missing_cases=len(missing),
                    extra_cases=len(extra),
                ),
                indent=2,
            )
        )


def physical_scales(nc, nsites, sumn, t, p, v_tv):
    """Per-column unit-aware physical scales so tiny derivatives are strictly checked."""
    v_tp = v_tv / sumn
    rt = R_GAS * t
    return [
        # F, F_T, F_V, F_n
        sumn,
        sumn / t,
        sumn / v_tv,
        *([1.0] * nc),
        # F_TT, F_VV, F_TV, F_Tn, F_Vn, F_nn
        sumn / (t * t),
        sumn / (v_tv * v_tv),
        sumn / (t * v_tv),
        *([1.0 / t] * nc),
        *([1.0 / v_tv] * nc),
        *([1.0 / sumn] * (nc * nc)),
        # X_k
        *([1.0] * nsites),
        # Z, dZ/dT, dZ/dP, dZ/dn
        1.0,
        1.0 / t,
        1.0 / p,
        *([1.0] * nc),
        # v_tp, dv/dT, dv/dP, dV/dn
        v_tp,
        v_tp / t,
        v_tp / p,
        *([v_tp] * nc),
        # h_res, dh_res/dT, dh_res/dP, dH_res/dn
        rt,
        R_GAS,
        v_tp,
        *([rt] * nc),
        # h, dh/dT, dh/dP, dH/dn
        rt,
        R_GAS,
        v_tp,
        *([rt] * nc),
        # s_res, ds_res/dT, ds_res/dP, dS_res/dn
        R_GAS,
        R_GAS / t,
        v_tp / t,
        *([R_GAS] * nc),
        # s, ds/dT, ds/dP, dS/dn
        R_GAS,
        R_GAS / t,
        v_tp / t,
        *([R_GAS] * nc),
        # ln_phi, dlnphi/dT, dlnphi/dP, dlnphi/dn
        *([1.0] * nc),
        *([1.0 / t] * nc),
        *([1.0 / p] * nc),
        *([1.0] * (nc * nc)),
        # sos_tp
        100.0,
    ]


def rdf_physical_scales(nc, sumn, t, v):
    return [
        1.0,
        1.0 / t,
        1.0 / v,
        *([1.0 / sumn] * nc),
        1.0 / (t * t),
        1.0 / (t * v),
        *([1.0 / (t * sumn)] * nc),
        1.0 / (v * v),
        *([1.0 / (v * sumn)] * nc),
        *([1.0 / (sumn * sumn)] * (nc * nc)),
    ]


def validate(directory):
    case_list = list(cases())
    rdf_list = list(rdf_cases())
    try:
        verify_case_grid_before_run(case_list, rdf_list)
    except ValueError as err:
        print(err)
        raise SystemExit(1) from err

    good = []
    failed = []
    censored = []
    unqualified = []
    rdf_good = []
    with DeadlineReferenceProcess(directory) as oracle:
        for kind, param_ref, comps, moles, t, p, phase in case_list:
            simplified = kind in ("sPC-SAFT", "sPCP-SAFT")
            polar = kind in ("PCP-SAFT", "sPCP-SAFT")
            request = dict(
                mode="PCSAFT_FULL",
                kind=kind,
                simplified=simplified,
                polar=polar,
                parameter_reference=param_ref,
                components=comps,
                moles=moles,
                temperature=t,
                pressure=p,
                phase=phase,
            )
            result = oracle.call(request, timeout=15.0)
            v_tv = result.get("volume", 0.0)
            label = (
                f"PCSAFT_FULL|{kind}|{param_ref}|{comps}|{t}|{p}|{v_tv:.17e}|"
                f"{','.join(map(repr, moles))}|{phase}"
            )
            record = dict(case=label, request=request, result=result)
            if "reference_timeout" in result:
                censored.append(record)
                continue
            if "values" not in result:
                failed.append(record)
                continue
            if not all(map(math.isfinite, result["values"])):
                unqualified.append(record)
                continue
            good.append(record)

        for comps, moles, t, v in rdf_list:
            request = dict(
                mode="SPCSAFT_RDF",
                components=comps,
                moles=moles,
                temperature=t,
                volume=v,
            )
            result = oracle.call(request, timeout=10.0)
            if "values" not in result or not all(map(math.isfinite, result["values"])):
                print(json.dumps(dict(rdf_error=request, result=result), indent=2))
                raise SystemExit(1)
            label = f"SPCSAFT_RDF|{comps}|{t}|{v:.17e}|{','.join(map(repr, moles))}"
            rdf_good.append(dict(case=label, request=request, result=result))

    try:
        verify_expected_coverage(good, failed, censored, unqualified)
    except ValueError as err:
        print(err)
        raise SystemExit(1) from err

    subprocess.run(
        ["cargo", "build", "--release", "--example", "parity", "-j", "2"],
        cwd=ROOT,
        check=True,
    )
    exe = ROOT / "target/release/examples" / ("parity.exe" if os.name == "nt" else "parity")
    all_inputs = [r["case"] for r in good] + [r["case"] for r in rdf_good]
    output = subprocess.run(
        [str(exe)],
        input="\n".join(all_inputs) + "\n",
        capture_output=True,
        text=True,
        check=True,
    ).stdout.splitlines()
    assert len(output) == len(all_inputs)

    full_outputs = output[: len(good)]
    rdf_outputs = output[len(good) :]

    rows = []
    differences = []
    worst_rel = 0.0
    rtol = 1e-8
    atol_scale = 1e-11
    for r, out in zip(good, full_outputs):
        expected = r["result"]["values"]
        rows.append([r["case"], str(r["result"]["nsites"]), *[f"{v:.17e}" for v in expected]])
        if out.startswith("ERROR"):
            differences.append(dict(case=r["case"], rust=out))
            continue
        native = list(map(float, out.split(",")))
        assert len(native) == len(expected), f"{r['case']}: {len(native)} vs {len(expected)}"
        req = r["request"]
        nc = len(req["moles"])
        nsites = int(r["result"]["nsites"])
        sumn = sum(req["moles"])
        v_tv = r["result"]["volume"]
        scales = physical_scales(nc, nsites, sumn, req["temperature"], req["pressure"], v_tv)
        assert len(scales) == len(expected)
        bad = []
        for idx, (a, b, phys) in enumerate(zip(expected, native, scales)):
            atol_i = atol_scale * phys
            denom = max(abs(a), abs(b), 1e-4 * phys)
            err = abs(a - b)
            worst_rel = max(worst_rel, err / denom)
            if not math.isfinite(b) or err > atol_i + rtol * max(abs(a), abs(b)):
                bad.append([idx, a, b, phys])
        if bad:
            differences.append(dict(case=r["case"], quantities=bad))

    rdf_worst_rel = 0.0
    for r, out in zip(rdf_good, rdf_outputs):
        expected = r["result"]["values"]
        native = list(map(float, out.split(",")))
        assert len(native) == len(expected)
        req = r["request"]
        nc = len(req["moles"])
        sumn = sum(req["moles"])
        scales = rdf_physical_scales(nc, sumn, req["temperature"], req["volume"])
        bad = []
        for idx, (a, b, phys) in enumerate(zip(expected, native, scales)):
            atol_i = 1e-12 * phys
            denom = max(abs(a), abs(b), 1e-4 * phys)
            err = abs(a - b)
            rdf_worst_rel = max(rdf_worst_rel, err / denom)
            if not math.isfinite(b) or err > atol_i + 1e-11 * max(abs(a), abs(b)):
                bad.append([idx, a, b, phys])
        if bad:
            differences.append(dict(case=r["case"], quantities=bad))

    if differences:
        print(json.dumps(dict(differences=len(differences), sample=differences[:8]), indent=2))
        raise SystemExit(1)

    path = ROOT / "tests/fixtures/pcsaft_full.csv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "nsites", "values"])
        writer.writerows(rows)

    meta = json.loads((ROOT / ".reference/metadata.json").read_text(encoding="utf-8"))
    meta.update(
        cases=len(rows),
        rdf_cases=len(rdf_good),
        generator="tools/rust_port/validate_pcsaft.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        relative_tolerance=rtol,
        absolute_tolerance_scale=atol_scale,
        tolerance_policy="unit-aware per-quantity physical scales (atol_i = absolute_tolerance_scale * phys_scale_i)",
        reference_errors=failed,
        censored_calls=censored,
        unqualified_returns=unqualified,
        worst_scaled_error=worst_rel,
        worst_rdf_scaled_error=rdf_worst_rel,
    )
    path.with_suffix(".metadata.json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            dict(
                qualified=len(rows),
                rdf_qualified=len(rdf_good),
                failed=len(failed),
                censored=len(censored),
                unqualified=len(unqualified),
                differences=0,
                worst_scaled_error=worst_rel,
                worst_rdf_scaled_error=rdf_worst_rel,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler-directory")
    validate(parser.parse_args().compiler_directory)
