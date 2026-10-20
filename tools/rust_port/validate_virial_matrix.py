"""Compare the native matrix of second virial coefficients against the
exact-source oracle.

The original `secondVirialCoeffMatrix` fills the matrix of the pure and the
cross second virial coefficients of the current model. A reference probe calls
the original routine and the native model builds the same matrix.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the equation of state, the components and the temperature.
CASES = [
    ("PR", ["C1"], 150.0),
    ("PR", ["C1"], 300.0),
    ("PR", ["C1", "C3"], 300.0),
    ("PR", ["C1", "C3"], 150.0),
    ("PR", ["CO2", "C1", "N2"], 300.0),
    ("SRK", ["C1", "C3"], 250.0),
]
# The matrix carries the truncation error of the finite volume estimate of the
# original, and the cross coefficients inherit that error.
RTOL = 1.0e-3
ATOL = 1.0e-8


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
    lines = []
    for eos, components, temperature in CASES:
        lines.append(
            "|".join(
                [
                    "VIRIALMATRIX",
                    eos,
                    ",".join(components),
                    f"{temperature:.17e}",
                ]
            )
        )
    native_values = native(lines)

    with ReferenceProcess(directory) as oracle:
        reference_values = []
        for eos, components, temperature in CASES:
            reference_values.append(
                oracle.call(
                    {
                        "mode": "VIRIAL_MATRIX",
                        "eos": eos,
                        "components": components,
                        "temperature": temperature,
                    }
                )
            )

    failures = []
    rows = []
    worst = 0.0
    for index, case in enumerate(CASES):
        eos, components, temperature = case
        label = f"{eos}|{','.join(components)}|{temperature:.17e}"
        expected = reference_values[index]
        actual = native_values[index]
        if "values" not in expected:
            failures.append(dict(case=label, reason="source did not report values", detail=expected))
            continue
        if isinstance(actual, str):
            failures.append(dict(case=label, reason="native error", native=actual))
            continue
        if len(expected["values"]) != len(actual):
            failures.append(
                dict(case=label, reason="matrix size", source=len(expected["values"]), native=len(actual))
            )
            continue
        bad = []
        for position, (source, native_value) in enumerate(zip(expected["values"], actual)):
            if abs(source - native_value) > ATOL + RTOL * abs(source):
                bad.append(position)
            else:
                worst = max(worst, abs(source - native_value) / max(abs(source), 1.0e-30))
        if bad:
            failures.append(
                dict(
                    case=label,
                    entries=bad,
                    source=expected["values"],
                    native=actual,
                )
            )
            continue
        rows.append([label, *expected["values"]])

    path = ROOT / "tests/fixtures/virial_matrix.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "matrix"])
        for row in rows:
            writer.writerow(
                [row[0], ",".join(f"{value:.17e}" for value in row[1:])]
            )
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_virial_matrix.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original secondVirialCoeffMatrix; the source estimates the pure and the cross "
            "coefficients from the residual Helmholtz energy at a molar volume of 1e7 cubic metre "
            "per mole, and the native model uses the same estimate, so the matrix carries the "
            "truncation error of that estimate"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/virial_matrix_failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    print(
        json.dumps(
            {
                "attempts": len(CASES),
                "compared": len(rows),
                "failures": len(failures),
                "worst_relative_difference": worst,
            },
            indent=2,
        )
    )
    if failures:
        print(json.dumps(failures[:4], indent=2))
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler-directory")
    main(parser.parse_args().compiler_directory)
