"""Compare the native Lennard-Jones splined terms against the exact-source
oracle.

The original `calc_ai_reduced_LJs_ex` evaluates the reduced perturbation
terms of the Barker-Henderson perturbation theory and the original
`ljs_bh_get_bh_diameter_div_sigma` the reduced hard-sphere diameter. A
reference oracle calls the Python interface of the fork, which reaches both
original routines.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the reduced temperature and the reduced density.
CASES = [
    (1.0, 0.5),
    (0.8, 1.0),
    (1.2, 0.2),
    (0.5, 2.0),
    (0.6, 0.4),
    (0.9, 1.5),
    (1.5, 0.8),
    (0.4, 0.1),
]
QUANTITIES = ["a1", "a2", "a3", "diameter_div_sigma"]
RTOL = 1.0e-10
ATOL = 1.0e-12


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
    lines = [f"LJS|{t_star:.17e}|{rho_star:.17e}" for t_star, rho_star in CASES]
    native_values = native(lines)

    with ReferenceProcess(directory) as oracle:
        reference_values = []
        for t_star, rho_star in CASES:
            reference_values.append(
                oracle.call(
                    {
                        "mode": "LJS",
                        "reduced_temperature": t_star,
                        "reduced_density": rho_star,
                    }
                )
            )

    failures = []
    rows = []
    worst = 0.0
    for index, case in enumerate(CASES):
        t_star, rho_star = case
        label = f"{t_star:.17e}|{rho_star:.17e}"
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
            if abs(source - native_value) > ATOL + RTOL * abs(source):
                bad.append(QUANTITIES[position])
            else:
                worst = max(worst, abs(source - native_value) / max(abs(source), 1.0e-30))
        if bad:
            failures.append(
                dict(case=label, quantities=bad, source=expected["values"], native=actual)
            )
            continue
        rows.append([label, *expected["values"]])

    path = ROOT / "tests/fixtures/ljs.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *QUANTITIES])
        for row in rows:
            writer.writerow([row[0], *[f"{value:.17e}" for value in row[1:]]])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_ljs.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original calc_ai_reduced_LJs_ex and ljs_bh_get_bh_diameter_div_sigma, which the "
            "Python interface exposes as get_pert_a and get_bh_diameter_div_sigma; the reduced "
            "terms are divided by the well depth, its square and its cube, and the diameter is the "
            "Barker-Henderson hard-sphere diameter divided by the segment diameter"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/ljs_failures.json").write_text(json.dumps(failures, indent=2) + "\n")
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
