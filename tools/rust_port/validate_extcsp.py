"""Exact-source parity qualification for the Extended CSP equation of state."""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from deadline_reference import DeadlineReferenceProcess
from reference import ROOT


def cases():
    # (shape_eos, shape_mixing, shape_alpha, ref_eos, ref_comp, param_ref, comps, moles, states)
    grid = [
        (
            "SRK",
            "vdW",
            "Classic",
            "NIST_MEOS",
            "C3",
            "DEFAULT",
            "C3",
            [1.0],
            [(280.0, 2.0e6, "L"), (350.0, 1.0e6, "V"), (420.0, 6.0e6, "V")],
        ),
        (
            "SRK",
            "vdW",
            "Classic",
            "NIST_MEOS",
            "C3",
            "DEFAULT",
            "C1,C3",
            [0.6, 0.4],
            [(240.0, 5.0e6, "L"), (320.0, 2.0e6, "V"), (380.0, 8.0e6, "V")],
        ),
        (
            "PR",
            "vdW",
            "Classic",
            "NIST_MEOS",
            "C3",
            "DEFAULT",
            "C1,C3",
            [0.7, 0.3],
            [(240.0, 5.0e6, "L"), (320.0, 2.0e6, "V"), (380.0, 8.0e6, "V")],
        ),
        (
            "SRK",
            "vdW",
            "Classic",
            "NIST_MEOS",
            "C3",
            "DEFAULT",
            "CO2,C1",
            [0.45, 0.55],
            [(240.0, 6.0e6, "L"), (310.0, 3.0e6, "V")],
        ),
        (
            "PR",
            "vdW",
            "Classic",
            "NIST_MEOS",
            "C3",
            "DEFAULT",
            "C1,C2,C3",
            [0.5, 0.3, 0.2],
            [(230.0, 5.0e6, "L"), (320.0, 2.0e6, "V"), (380.0, 8.0e6, "V")],
        ),
        # Non-unit total moles (tests extensive TV scaling directly against csp_calcFres)
        (
            "SRK",
            "vdW",
            "Classic",
            "NIST_MEOS",
            "C3",
            "DEFAULT",
            "C1,C2,C3",
            [1.0, 0.6, 0.4],
            [(300.0, 3.0e6, "V"), (240.0, 6.0e6, "L")],
        ),
        # Activity-based mixing rules through shape EOS
        (
            "SRK",
            "HV1",
            "Classic",
            "NIST_MEOS",
            "C3",
            "DEFAULT",
            "C1,H2O",
            [0.4, 0.6],
            [(350.0, 1.0e6, "V"), (450.0, 5.0e6, "V")],
        ),
        (
            "PR",
            "HV2",
            "Classic",
            "NIST_MEOS",
            "C3",
            "DEFAULT",
            "C1,H2O",
            [0.4, 0.6],
            [(350.0, 1.0e6, "V"), (450.0, 5.0e6, "V")],
        ),
        (
            "PR",
            "NRTL",
            "Classic",
            "NIST_MEOS",
            "C3",
            "Dicko2012",
            "H2S,C3",
            [0.4, 0.6],
            [(350.0, 2.0e6, "V"), (400.0, 5.0e6, "V")],
        ),
        (
            "SRK",
            "UNIFAC",
            "Classic",
            "NIST_MEOS",
            "C3",
            "DEFAULT",
            "C1,C3",
            [0.3, 0.7],
            [(300.0, 1.0e6, "V"), (380.0, 4.0e6, "V")],
        ),
    ]
    for sh_eos, sh_mix, sh_alpha, ref_eos, ref_comp, param_ref, comps, moles, states in grid:
        for t, p, phase in states:
            yield sh_eos, sh_mix, sh_alpha, ref_eos, ref_comp, param_ref, comps, moles, t, p, phase


def validate(directory):
    good = []
    failed = []
    censored = []
    unqualified = []
    with DeadlineReferenceProcess(directory) as oracle:
        for (
            sh_eos,
            sh_mix,
            sh_alpha,
            ref_eos,
            ref_comp,
            param_ref,
            comps,
            moles,
            t,
            p,
            phase,
        ) in cases():
            request = dict(
                mode="EXTCSP",
                shape_eos=sh_eos,
                shape_mixing=sh_mix,
                shape_alpha=sh_alpha,
                ref_eos=ref_eos,
                ref_comp=ref_comp,
                parameter_reference=param_ref,
                components=comps,
                moles=moles,
                temperature=t,
                pressure=p,
                phase=phase,
            )
            result = oracle.call(request, timeout=10.0)
            v_tv = result.get("volume", 0.0)
            label = (
                f"EXTCSP|{sh_eos}|{sh_mix}|{sh_alpha}|{ref_eos}|{ref_comp}|{param_ref}|"
                f"{comps}|{t}|{p}|{v_tv:.17e}|{','.join(map(repr, moles))}|{phase}"
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

    subprocess.run(
        ["cargo", "build", "--release", "--example", "parity", "-j", "2"],
        cwd=ROOT,
        check=True,
    )
    exe = ROOT / "target/release/examples" / ("parity.exe" if os.name == "nt" else "parity")
    output = subprocess.run(
        [str(exe)],
        input="\n".join(r["case"] for r in good) + "\n",
        capture_output=True,
        text=True,
        check=True,
    ).stdout.splitlines()
    assert len(output) == len(good)

    rows = []
    differences = []
    worst_rel = 0.0
    rtol = 1e-7
    atol = 1e-9
    for r, out in zip(good, output):
        expected = r["result"]["values"]
        rows.append([r["case"], *[f"{v:.17e}" for v in expected]])
        if out.startswith("ERROR"):
            differences.append(dict(case=r["case"], rust=out))
            continue
        native = list(map(float, out.split(",")))
        assert len(native) == len(expected)
        bad = []
        for idx, (a, b) in enumerate(zip(expected, native)):
            scale = max(abs(a), abs(b), 1.0)
            err = abs(a - b)
            worst_rel = max(worst_rel, err / scale)
            if not math.isfinite(b) or err > atol + rtol * max(abs(a), abs(b)):
                bad.append([idx, a, b])
        if bad:
            differences.append(dict(case=r["case"], quantities=bad))

    path = ROOT / "tests/fixtures/extcsp.csv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "values"])
        writer.writerows(rows)

    meta = json.loads((ROOT / ".reference/metadata.json").read_text(encoding="utf-8"))
    meta.update(
        cases=len(rows),
        generator="tools/rust_port/validate_extcsp.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        relative_tolerance=rtol,
        absolute_tolerance=atol,
        reference_errors=failed,
        censored_calls=censored,
        unqualified_returns=unqualified,
        worst_scaled_error=worst_rel,
    )
    path.with_suffix(".metadata.json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8"
    )
    (ROOT / ".reference/extcsp-differences.json").write_text(
        json.dumps(differences, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            dict(
                qualified=len(rows),
                failed=len(failed),
                censored=len(censored),
                unqualified=len(unqualified),
                differences=len(differences),
                worst_scaled_error=worst_rel,
            ),
            indent=2,
        )
    )
    if differences:
        print(json.dumps(differences[:8], indent=2))
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler-directory")
    validate(parser.parse_args().compiler_directory)
