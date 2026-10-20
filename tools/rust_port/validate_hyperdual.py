"""Compare the native hyperdual numbers against the exact-source oracle.

The original `hyperdual_mod` provides a hyperdual number with up to three
differentiation directions and the elementary functions of the calculus. A
reference probe evaluates the original functions and the native model evaluates
the same expressions.
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

FUNCTIONS = [
    "exp",
    "log",
    "sqrt",
    "sin",
    "cos",
    "tan",
    "sinh",
    "cosh",
    "tanh",
    "asin",
    "acos",
    "atan",
]
COMPONENTS = ["f0", "f1", "f2", "f3", "f12", "f13", "f23", "f123"]
# The cases cover the first, second and third differentiation order.
CASES = [
    (1.5, 1.0, 2.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 2),
    (2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1),
    (0.75, 0.5, -1.5, 2.0, 1.0, -0.5, 0.25, 0.125, 3),
    (-0.4, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3),
    (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0),
    (0.25, 2.0, -3.0, 1.5, 0.75, -0.25, 0.5, 0.25, 2),
]
RTOL = 1e-12
ATOL = 1e-14


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
    for case in CASES:
        lines.append("HYPERDUAL|" + ",".join(f"{value:.17e}" for value in case))
    native_values = native(lines)

    with ReferenceProcess(directory) as oracle:
        reference = []
        for case in CASES:
            result = oracle.call(
                {
                    "mode": "HYPERDUAL",
                    "inputs": list(case),
                }
            )
            reference.append(result.get("values"))

    failures = []
    rows = []
    worst = 0.0
    for index, case in enumerate(CASES):
        expected = reference[index]
        actual = native_values[index]
        if expected is None:
            failures.append(dict(case=str(case), reason="source did not report values"))
            continue
        if isinstance(actual, str):
            failures.append(dict(case=str(case), reason="native error", native=actual))
            continue
        bad = []
        for position, (a, b) in enumerate(zip(expected, actual)):
            error = abs(a - b)
            if math.isnan(a) and math.isnan(b):
                # Both sides reproduce the same undefined value.
                continue
            worst = max(worst, error / max(abs(a), 1.0))
            if not math.isfinite(b) or error > ATOL + RTOL * abs(a):
                bad.append(position)
        if bad:
            failures.append(
                dict(case=str(case), quantities=bad, reference=expected, native=actual)
            )
            continue
        rows.append([",".join(f"{value:.17e}" for value in case), *expected])

    path = ROOT / "tests/fixtures/hyperdual.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["inputs", "values"])
        for row in rows:
            writer.writerow([row[0], ",".join(f"{value:.17e}" for value in row[1:])])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_hyperdual.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        functions=FUNCTIONS,
        components=COMPONENTS,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original hyperdual elementary functions; every function contributes its eight "
            "components in the order f0, f1, f2, f3, f12, f13, f23, f123"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/hyperdual_failures.json").write_text(json.dumps(failures, indent=2) + "\n")
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
