"""Compare the native critical points against the exact-source oracle.

The original `calcCriticalTV` locates the critical point of a fluid in the
temperature and the molar volume with the stability matrix of the mixture.
A reference oracle calls the Python interface of the fork, which reaches the
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

# The cases: the equation of state, the component and the composition. The
# critical point of the source needs a single component.
CASES = [
    ("PR", ["C1"], [1.0]),
    ("PR", ["CO2"], [1.0]),
    ("PR", ["C3"], [1.0]),
    ("PR", ["H2O"], [1.0]),
    ("SRK", ["N2"], [1.0]),
    ("SRK", ["C1"], [1.0]),
    ("PR", ["C1", "C3"], [1.0, 0.0]),
]
QUANTITIES = ["temperature", "volume", "pressure"]
# The source solves the critical point with a tolerance of 1e-7.
RTOL = 1.0e-6
ATOL = [1.0e-6, 1.0e-12, 1.0e-3]


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
    for eos, components, composition in CASES:
        lines.append(
            "|".join(
                [
                    "CRITICAL",
                    eos,
                    ",".join(components),
                    ",".join(f"{value:.17e}" for value in composition),
                ]
            )
        )
    native_values = native(lines)

    with ReferenceProcess(directory) as oracle:
        reference_values = []
        for eos, components, composition in CASES:
            reference_values.append(
                oracle.call(
                    {
                        "mode": "CRITICAL",
                        "eos": eos,
                        "components": components,
                        "feed": composition,
                    }
                )
            )

    failures = []
    rows = []
    worst = 0.0
    for index, case in enumerate(CASES):
        eos, components, composition = case
        label = (
            f"{eos}|{','.join(components)}|"
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
        bad = []
        for position, (source, native_value) in enumerate(zip(expected["values"], actual)):
            if abs(source - native_value) > ATOL[position] + RTOL * abs(source):
                bad.append(QUANTITIES[position])
            else:
                worst = max(worst, abs(source - native_value) / max(abs(source), 1.0e-30))
        if bad:
            failures.append(
                dict(
                    case=label,
                    quantities=bad,
                    source=expected["values"],
                    native=actual,
                )
            )
            continue
        rows.append([label, *expected["values"]])

    path = ROOT / "tests/fixtures/critical.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *QUANTITIES])
        for row in rows:
            writer.writerow([row[0], *[f"{value:.17e}" for value in row[1:]]])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_critical.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerances=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original calcCriticalTV, which the Python interface exposes as critical; the "
            "source reports the temperature, the molar volume and the pressure of the critical "
            "point, and the native model solves the two volume derivatives of the pressure, which "
            "are the conditions of the stability matrix for a single component"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/critical_failures.json").write_text(
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
