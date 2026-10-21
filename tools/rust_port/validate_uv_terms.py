"""Compare the native UV theory terms against the exact-source oracle.

The original `uv_a1_u` and `uv_a1_b2` evaluate the first and the second
order terms of the UV theory of the splined Lennard-Jones model with their
partial derivatives. A reference probe calls the original routines.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the segment diameter, the well depth divided by the Boltzmann
# constant, the temperature and the reduced density.
CASES = [
    (3.42e-10, 124.0, 100.0, 0.5),
    (3.42e-10, 124.0, 200.0, 1.0),
    (3.50e-10, 95.0, 300.0, 0.2),
    (3.42e-10, 124.0, 150.0, 1.5),
    (3.50e-10, 95.0, 250.0, 0.8),
]
QUANTITIES = ["a1", "a1_r", "a1_T", "a1_rr", "a1_TT", "a1_rT", "b2", "b2_T", "b2_TT"]
RTOL = 1.0e-11
ATOL = 1.0e-14


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
    with ReferenceProcess(directory) as oracle:
        reference_values = []
        for case in CASES:
            result = oracle.call(
                {
                    "mode": "UV_TERMS",
                    "sigma": case[0],
                    "eps_divk": case[1],
                    "temperature": case[2],
                    "rho": case[3],
                }
            )
            if "values" not in result:
                raise SystemExit(f"the source did not report the terms: {result}")
            reference_values.append(result)

    lines = [
        "UVTERMS|{0:.17e}|{1:.17e}|{2:.17e}|{3:.17e}".format(*case) for case in CASES
    ]
    native_values = native(lines)

    failures = []
    rows = []
    worst = 0.0
    for index, case in enumerate(CASES):
        label = "|".join(f"{value:.17e}" for value in case)
        source = reference_values[index]["values"]
        actual = native_values[index]
        if isinstance(actual, str):
            failures.append(dict(case=label, reason="native error", native=actual))
            continue
        bad = []
        for position, (source_value, native_value) in enumerate(zip(source, actual)):
            if abs(source_value - native_value) > ATOL + RTOL * abs(source_value):
                bad.append((QUANTITIES[position], source_value, native_value))
            else:
                worst = max(
                    worst, abs(source_value - native_value) / max(abs(source_value), 1.0e-30)
                )
        if bad:
            failures.append(dict(case=label, reason="terms", differences=bad[:6]))
            continue
        rows.append([label, *source])

    path = ROOT / "tests/fixtures/uv_terms.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *QUANTITIES])
        for row in rows:
            writer.writerow([row[0], *[f"{value:.17e}" for value in row[1:]]])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_uv_terms.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original uv_a1_u and uv_a1_b2, which evaluate the first and the second order "
            "terms of the UV theory of the splined Lennard-Jones model with their partial "
            "derivatives; the probe uses the hard-sphere diameter of the Barker-Henderson reference"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/uv_terms_failures.json").write_text(
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
