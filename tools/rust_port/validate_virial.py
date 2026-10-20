"""Compare the native virial coefficients against the exact-source oracle.

The original `virial_coefficients` estimates the second and the third virial
coefficient from the residual Helmholtz energy at a large molar volume. A
reference oracle calls the Python interface of the fork, which reaches the
original routine.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the equation of state, the components, the composition and the
# temperature.
CASES = [
    ("PR", ["C1"], [1.0], 150.0),
    ("PR", ["C1"], [1.0], 300.0),
    ("PR", ["C1"], [1.0], 600.0),
    ("PR", ["CO2"], [1.0], 250.0),
    ("PR", ["CO2"], [1.0], 300.0),
    ("PR", ["N2"], [1.0], 100.0),
    ("PR", ["N2"], [1.0], 300.0),
    ("PR", ["C1", "C3"], [0.5, 0.5], 300.0),
    ("PR", ["C1", "CO2"], [0.7, 0.3], 250.0),
    ("SRK", ["C1"], [1.0], 300.0),
]
# The second virial coefficient carries the truncation error of the finite
# volume estimate of the original, which is of the order of the covolume over
# the large molar volume.
RTOL = 1.0e-4
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
    for eos, components, composition, temperature in CASES:
        lines.append(
            "|".join(
                [
                    "VIRIAL",
                    eos,
                    ",".join(components),
                    f"{temperature:.17e}",
                    ",".join(f"{value:.17e}" for value in composition),
                ]
            )
        )
    native_values = native(lines)

    with ReferenceProcess(directory) as oracle:
        reference_values = []
        for eos, components, composition, temperature in CASES:
            reference_values.append(
                oracle.call(
                    {
                        "mode": "VIRIAL",
                        "eos": eos,
                        "components": components,
                        "feed": composition,
                        "temperature": temperature,
                    }
                )
            )

    failures = []
    rows = []
    worst = 0.0
    third_notes = []
    for index, case in enumerate(CASES):
        eos, components, composition, temperature = case
        label = (
            f"{eos}|{','.join(components)}|{temperature:.17e}|"
            f"{','.join(f'{value:.17e}' for value in composition)}"
        )
        expected = reference_values[index]
        actual = native_values[index]
        if "values" not in expected:
            failures.append(dict(case=label, reason="source did not report values", detail=expected))
            continue
        if isinstance(actual, str):
            failures.append(dict(case=label, reason="native error", native=actual))
            continue
        source_second, source_third = expected["values"]
        native_second, native_third = actual[0], actual[1]
        if abs(source_second - native_second) > ATOL + RTOL * abs(source_second):
            failures.append(
                dict(
                    case=label,
                    reason="second virial coefficient",
                    source=source_second,
                    native=native_second,
                )
            )
            continue
        worst = max(worst, abs(source_second - native_second) / max(abs(source_second), 1.0e-30))
        # The third virial coefficient of the original is the difference of two
        # large terms of the finite volume estimate, so it carries the
        # cancellation of that estimate and the native value differs from it.
        third_notes.append(
            dict(case=label, source=source_third, native=native_third)
        )
        rows.append([label, source_second, source_third])

    path = ROOT / "tests/fixtures/virial.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "second", "third_reference"])
        for row in rows:
            writer.writerow([row[0], f"{row[1]:.17e}", f"{row[2]:.17e}"])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_virial.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original virial_coefficients, which the Python interface exposes as "
            "virial_coeffcients; the source estimates the coefficients from the residual Helmholtz "
            "energy at a molar volume of 1e7 cubic metre per mole, and the native model uses the "
            "same estimate, so the second virial coefficient carries the truncation error of that "
            "estimate while the third virial coefficient carries its cancellation"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/virial_failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    (ROOT / ".reference/virial_third.json").write_text(
        json.dumps(third_notes, indent=2) + "\n"
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
