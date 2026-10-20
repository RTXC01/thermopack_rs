"""Compare the native spline polynomials against the exact-source oracle.

The original `calc_ai_LJs_ex` evaluates one polynomial of the spline of the
excess energy of the splined Lennard-Jones model with its partial
derivatives. A reference probe calls the original routine.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the reduced temperature, the reduced density and the index of the
# polynomial. The first polynomial is the first order term, the second the
# second order term, the third the third order term and the fourth the
# correction of the second order term.
CASES = [
    (1.0, 0.5, 0),
    (1.0, 0.5, 1),
    (1.0, 0.5, 2),
    (1.0, 0.5, 3),
    (0.8, 1.0, 0),
    (0.8, 1.0, 2),
    (1.2, 0.2, 1),
    (0.5, 2.0, 3),
]
QUANTITIES = ["a", "a_r", "a_x", "a_rr", "a_xx", "a_rx"]
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
    lines = [
        f"LJSSPLINE|{t_star:.17e}|{rho_star:.17e}|{index}"
        for t_star, rho_star, index in CASES
    ]
    native_values = native(lines)

    with ReferenceProcess(directory) as oracle:
        reference_values = []
        for t_star, rho_star, index in CASES:
            reference_values.append(
                oracle.call(
                    {
                        "mode": "LJS_SPLINE",
                        "reduced_temperature": t_star,
                        "reduced_density": rho_star,
                        "index": index,
                    }
                )
            )

    failures = []
    rows = []
    worst = 0.0
    for index, case in enumerate(CASES):
        t_star, rho_star, polynomial = case
        label = f"{t_star:.17e}|{rho_star:.17e}|{polynomial}"
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

    path = ROOT / "tests/fixtures/ljs_spline.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *QUANTITIES])
        for row in rows:
            writer.writerow([row[0], *[f"{value:.17e}" for value in row[1:]]])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_ljs_spline.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original calc_ai_LJs_ex, which evaluates one polynomial of the spline of the "
            "excess energy with the prefactor of one; the derivatives are the first and the second "
            "derivative with respect to the reduced temperature, the first and the second "
            "derivative with respect to the reduced density and the mixed derivative"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/ljs_spline_failures.json").write_text(
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
