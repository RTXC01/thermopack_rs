"""Qualify the UNIFAC mixing rule inside the cubic equation of state.

The original model selects the UNIFAC mixing rule through the cubic
equation of state, so the end-to-end check compares the pressure, its
derivatives, the residual and total caloric properties and the fugacity
coefficients of the exact-source library against the native model.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
from pathlib import Path

from reference import ROOT, ReferenceProcess

MIXTURES = [
    (["C2", "NC4"], [0.5, 0.5]),
    (["C2", "NC4"], [0.9, 0.1]),
    (["NC4", "NC5"], [0.5, 0.5]),
    (["C1", "NC4"], [0.5, 0.5]),
    (["C1", "C3"], [0.3, 0.7]),
    (["H2O", "NC4"], [0.5, 0.5]),
    (["H2O", "C3"], [0.7, 0.3]),
    (["CO2", "NC4"], [0.4, 0.6]),
    (["N2", "CO2"], [0.5, 0.5]),
    (["BENZENE", "NC6"], [0.5, 0.5]),
    (["TOLU", "NC7"], [0.6, 0.4]),
    (["ACETONE", "NC5"], [0.5, 0.5]),
    (["ACETONE", "C1"], [0.5, 0.5]),
    (["H2S", "NC4"], [0.25, 0.75]),
    (["C3", "NC4", "C1"], [0.4, 0.35, 0.25]),
    (["H2O", "CO2", "C3"], [0.5, 0.3, 0.2]),
]
STATES = [
    (250.0, 1e5, "L"),
    (250.0, 1e6, "L"),
    (298.15, 5e5, "L"),
    (298.15, 2e6, "L"),
    (298.15, 1e7, "V"),
    (350.0, 1e6, "L"),
    (350.0, 5e6, "V"),
    (400.0, 2e6, "L"),
    (400.0, 2e7, "V"),
]
# The oracle reports the heat capacity at the ninth position; the native
# parity protocol reports the baseline enthalpy temperature derivative there.
# Both are skipped in the comparison.
LABELS = ["v", "p", "dpdt", "dpdv", "h_res", "s_res", "h", "s", "cp", "ln_phi"]
ABSOLUTE = [1e-14, 1e-2, 1e-6, 1e-2, 1e-7, 1e-9, 1e-7, 1e-9, 1e-7, 1e-9]
RELATIVE = 2e-9


def cases():
    for eos in ("PR", "SRK"):
        for components, composition in MIXTURES:
            for temperature, pressure, phase in STATES:
                yield eos, components, composition, temperature, pressure, phase


def native(lines):
    exe = ROOT / "target/release/examples" / ("parity.exe" if os.name == "nt" else "parity")
    result = subprocess.run(
        [str(exe)], input="\n".join(lines) + "\n", text=True, capture_output=True, check=True
    )
    rows = []
    for line in result.stdout.splitlines():
        rows.append(line if line.startswith("ERROR") else [float(value) for value in line.split(",")])
    assert len(rows) == len(lines), (len(rows), len(lines))
    return rows


def main(directory):
    request_list = list(cases())
    lines = []
    for eos, components, composition, temperature, pressure, phase in request_list:
        lines.append(
            "|".join(
                [
                    eos,
                    ",".join(components),
                    f"{temperature:.17e}",
                    f"{pressure:.17e}",
                    ",".join(f"{value:.17e}" for value in composition),
                    f"{phase}@unifac",
                ]
            )
        )
    native_values = native(lines)

    reference_values = []
    with ReferenceProcess(directory) as oracle:
        for eos, components, composition, temperature, pressure, phase in request_list:
            result = oracle.call(
                {
                    "mode": "V_PROP",
                    "eos": eos,
                    "components": ",".join(components),
                    "temperature": temperature,
                    "pressure": pressure,
                    "feed": list(composition),
                    "phase": phase,
                    "mixing": "UNIFAC",
                }
            )
            # The oracle reports the fugacity derivatives after the fugacity
            # coefficients; the comparison uses the leading quantities.
            values = result.get("values")
            reference_values.append(values[: 9 + len(composition)] if values else None)

    failures = []
    rows = []
    status_rows = []
    worst = {}
    for index, (eos, components, composition, temperature, pressure, phase) in enumerate(request_list):
        case = "|".join(
            [
                eos,
                ",".join(components),
                f"{temperature:.17e}",
                f"{pressure:.17e}",
                ",".join(f"{value:.17e}" for value in composition),
                phase,
                "UNIFAC",
            ]
        )
        expected = reference_values[index]
        actual = native_values[index]
        if expected is None:
            # The original terminated; the native side must report an error.
            if not isinstance(actual, str):
                failures.append(dict(case=case, reason="native value for terminated source",
                                     native=actual))
                status_rows.append([case, 1])
            else:
                status_rows.append([case, 1])
            continue
        if isinstance(actual, str):
            failures.append(dict(case=case, reason="native error", native=actual,
                                 reference=expected))
            status_rows.append([case, 2])
            continue
        bad = []
        for position, (a, b) in enumerate(zip(expected, actual)):
            if position == 8:
                continue
            label = LABELS[position] if position < len(LABELS) else "ln_phi"
            if not math.isfinite(a) or not math.isfinite(b):
                bad.append(label)
                continue
            error = abs(a - b) / max(abs(a), abs(b), 1.0)
            worst[label] = max(worst.get(label, 0.0), error)
            atol = ABSOLUTE[position] if position < len(ABSOLUTE) else ABSOLUTE[-1]
            if abs(a - b) > atol + RELATIVE * max(abs(a), abs(b)):
                bad.append(label)
        if bad:
            failures.append(dict(case=case, quantities=bad, reference=expected, native=actual))
            status_rows.append([case, 3])
            continue
        status_rows.append([case, 0])
        rows.append([case, *expected])

    path = ROOT / "tests/fixtures/unifac_cubic.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "v", "p", "dpdt", "dpdv", "h_res", "s_res", "h", "s", "cp",
                         "ln_phi..."])
        writer.writerows(rows)
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        generator="tools/rust_port/validate_unifac_cubic.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=LABELS,
        worst_scaled_error=worst,        tolerances={label: {"absolute": ABSOLUTE[index] if index < len(ABSOLUTE) else ABSOLUTE[-1],
                            "relative": RELATIVE} for index, label in enumerate(LABELS)},
        convention=(
            "cubic equation of state with the UNIFAC mixing rule; the source reports the heat "
            "capacity at the ninth position and the native parity protocol the baseline enthalpy "
            "temperature derivative, so that position is excluded from the comparison"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    status_path = ROOT / "tests/fixtures/unifac_cubic_status.csv"
    with status_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "source_status"])
        writer.writerows(status_rows)
    status_metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    status_metadata.update(
        cases=len(status_rows),
        generator="tools/rust_port/validate_unifac_cubic.py",
        fixture_sha256=hashlib.sha256(status_path.read_bytes()).hexdigest(),
        quantities=["source_status"],
        convention=(
            "0 compared and matching, 1 the original terminated, 2 the native side reported an "
            "error, 3 a numerical mismatch"
        ),
    )
    status_path.with_suffix(".metadata.json").write_text(json.dumps(status_metadata, indent=2) + "\n")
    (ROOT / ".reference/unifac_cubic_failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    print(json.dumps({"cases": len(request_list), "compared": len(rows),
                      "source_terminations": sum(1 for row in status_rows if row[1] == 1),
                      "native_errors": sum(1 for row in status_rows if row[1] == 2),
                      "mismatches": sum(1 for row in status_rows if row[1] == 3),
                      "failures": len(failures), "worst": worst}, indent=2))
    if failures:
        print(json.dumps(failures[:6], indent=2))
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler-directory")
    main(parser.parse_args().compiler_directory)
