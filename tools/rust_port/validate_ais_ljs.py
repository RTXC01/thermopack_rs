"""Compare the native perturbation terms of the Barker-Henderson theory against
the exact-source oracle.

The original `calc_ais_LJs_ex` and `calc_chi_LJs_ex` evaluate the terms of
the spline of the excess energy. A reference probe calls the original
routines.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the index of the term, the reduced distance of the hard-sphere
# centres and the reduced density.
CASES = [
    (1, 1.05, 0.5),
    (2, 1.10, 1.0),
    (3, 1.02, 0.2),
    (4, 1.20, 1.5),
    (1, 1.30, 0.8),
]
QUANTITIES = ["term", "chi"]
RTOL = 1.0e-13
ATOL = 1.0e-15


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
        for index, reduced_distance, reduced_density in CASES:
            result = oracle.call(
                {
                    "mode": "AIS_LJS",
                    "index": index,
                    "reduced_distance": reduced_distance,
                    "reduced_density": reduced_density,
                }
            )
            if "values" not in result:
                raise SystemExit(f"the source did not report the terms: {result}")
            reference_values.append(result)

    lines = [
        "AISLJS|{0}|{1:.17e}|{2:.17e}".format(*case) for case in CASES
    ]
    native_values = native(lines)

    failures = []
    rows = []
    worst = 0.0
    for index, case in enumerate(CASES):
        label = f"{case[0]}|{case[1]:.17e}|{case[2]:.17e}"
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

    path = ROOT / "tests/fixtures/ais_ljs.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *QUANTITIES])
        for row in rows:
            writer.writerow([row[0], *[f"{value:.17e}" for value in row[1:]]])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_ais_ljs.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original calc_ais_LJs_ex and calc_chi_LJs_ex, which evaluate the terms of the "
            "spline of the excess energy of the Barker-Henderson theory"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/ais_ljs_failures.json").write_text(
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
