"""Compare the native u-fraction of the UV theory against the exact-source
oracle.

The original `uv_phi` evaluates the u-fraction of the UV theory with its
partial derivatives. A reference probe calls the original routine.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the well depth divided by the Boltzmann constant, the temperature
# and the reduced density.
CASES = [
    (124.0, 100.0, 0.5),
    (124.0, 200.0, 1.0),
    (95.0, 300.0, 0.2),
    (124.0, 150.0, 1.5),
    (95.0, 250.0, 0.8),
]
QUANTITIES = ["u", "u_r", "u_T", "u_rr", "u_TT", "u_rT"]
RTOL = 1.0e-12
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
                    "mode": "UV_PHI",
                    "eps_divk": case[0],
                    "temperature": case[1],
                    "rho": case[2],
                }
            )
            if "values" not in result:
                raise SystemExit(f"the source did not report the fraction: {result}")
            reference_values.append(result)

    lines = [
        f"UVPHI|{eps_divk:.17e}|{temperature:.17e}|{rho:.17e}"
        for eps_divk, temperature, rho in CASES
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
            failures.append(dict(case=label, reason="fraction", differences=bad[:6]))
            continue
        rows.append([label, *source])

    path = ROOT / "tests/fixtures/uv_phi.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *QUANTITIES])
        for row in rows:
            writer.writerow([row[0], *[f"{value:.17e}" for value in row[1:]]])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_uv_phi.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original uv_phi, which evaluates the u-fraction of the UV theory of the splined "
            "Lennard-Jones model with its partial derivatives"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/uv_phi_failures.json").write_text(
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
