"""Compare the native White Bear cavity terms against the exact-source oracle.

The original `minimum_r_cavity_int_LJ`, `cavity_distribution_function` and
`lj_boltzmann_factor` evaluate the cavity terms of the White Bear hard-sphere
reference. A reference probe calls the original routines.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the temperature, the well depth divided by the Boltzmann constant,
# the packing fraction and the reduced separation with its derivatives.
CASES = [
    (100.0, 124.0, 0.30, 1.05, 1.0e-4, 2.0e-7),
    (150.0, 124.0, 0.50, 1.10, -2.0e-4, 3.0e-7),
    (200.0, 124.0, 0.52, 1.12, 5.0e-5, -1.0e-7),
    (300.0, 124.0, 0.60, 1.02, 1.0e-4, 1.0e-7),
    (120.0, 95.0, 0.45, 1.08, 0.0, 0.0),
]
QUANTITIES = [
    "r_min",
    "r_min_t",
    "r_min_tt",
    "y_hs",
    "y_hs_e",
    "y_hs_ee",
    "y_hs_t",
    "y_hs_tt",
    "y_hs_et",
    "e_lj",
    "e_lj_t",
    "e_lj_tt",
    "diameter",
    "diameter_t",
    "diameter_tt",
]
RTOL = 1.0e-11
ATOL = 1.0e-13


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
                    "mode": "WHITE_BEAR",
                    "temperature": case[0],
                    "eps_divk": case[1],
                    "eta": case[2],
                    "r_star": case[3],
                    "r_star_t": case[4],
                    "r_star_tt": case[5],
                }
            )
            if "values" not in result:
                raise SystemExit(f"the source did not report the terms: {result}")
            reference_values.append(result)

    lines = [
        "WHITEBEAR|{0:.17e}|{1:.17e}|{2:.17e}|{3:.17e}|{4:.17e}|{5:.17e}".format(*case)
        for case in CASES
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

    path = ROOT / "tests/fixtures/white_bear.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *QUANTITIES])
        for row in rows:
            writer.writerow([row[0], *[f"{value:.17e}" for value in row[1:]]])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_white_bear.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original minimum_r_cavity_int_LJ, cavity_distribution_function and "
            "lj_boltzmann_factor; the radius of the minimum of the cavity integral, the cavity "
            "correlation function of the hard spheres with its derivatives and the Boltzmann factor "
            "of the Lennard-Jones potential with its derivatives"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/white_bear_failures.json").write_text(
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
