"""Compare the native lnf-T solver against the exact-source oracle.

The original `mut_solver` converts a chemical potential specification into a
log fugacity specification at the initial guess and then solves the log
fugacity equation with a damped Newton iteration. A reference probe builds the
source model, evaluates the target at the true state and solves the problem
from a different guess, and the native model does the same.

The chemical potential variant of the original derives its target from the
initial guess, so its solve is a fixed point and the original returns a
degenerate state. The defect ledger records that, and the qualification uses
the log fugacity variant, where the target is an input.
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

MIXTURES = [
    ["C3", "N-H2"],
    ["N-H2", "O-H2"],
    ["C1", "C3"],
    ["CO2", "C3"],
]
# The target states stay in the vapor branch at every temperature, so the
# native model sees a positive pressure, and the guess is a different state.
TARGETS = {
    200.0: [[1000.0, 300.0], [700.0, 200.0]],
    250.0: [[2000.0, 500.0], [1500.0, 400.0]],
    298.15: [[4000.0, 1000.0], [3000.0, 800.0]],
    350.0: [[2000.0, 500.0], [1500.0, 400.0]],
}
TEMPERATURES = [200.0, 250.0, 298.15, 350.0]
# The relative tolerance covers the different Newton iterations of the two
# implementations; the absolute tolerance covers the smallest densities.
RTOL = 1e-6
ATOL = 1e-6


def cases():
    for components in MIXTURES:
        for temperature in TEMPERATURES:
            for target in TARGETS[temperature]:
                guess = [value * 0.9 for value in target]
                yield components, target, guess, temperature


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
    request_list = list(cases())
    with ReferenceProcess(directory) as oracle:
        targets = []
        for components, target, guess, temperature in request_list:
            result = oracle.call(
                {
                    "mode": "MU_T",
                    "kind": "TARGET",
                    "target_kind": 1,
                    "components": components,
                    "densities": target,
                    "temperature": temperature,
                }
            )
            targets.append(result.get("values"))
        reference = []
        for (components, target, guess, temperature), specification in zip(request_list, targets):
            result = oracle.call(
                {
                    "mode": "MU_T",
                    "kind": "LNF",
                    "components": components,
                    "densities": guess,
                    "temperature": temperature,
                    "target": specification,
                }
            )
            reference.append(result)

    lines = []
    for (components, target, guess, temperature), specification in zip(request_list, targets):
        lines.append(
            "|".join(
                [
                    "MUT",
                    ",".join(components),
                    f"{temperature:.17e}",
                    "0",
                    ",".join(f"{value:.17e}" for value in guess),
                    "LNF:" + ",".join(f"{value:.17e}" for value in specification),
                ]
            )
        )
    native_values = native(lines)

    failures = []
    rows = []
    worst = 0.0
    for index, (components, target, guess, temperature) in enumerate(request_list):
        case = "|".join(
            [
                "LNF",
                ",".join(components),
                f"{temperature:.17e}",
                ",".join(f"{value:.17e}" for value in target),
                ",".join(f"{value:.17e}" for value in guess),
            ]
        )
        expected = reference[index]
        actual = native_values[index]
        if "values" not in expected:
            failures.append(dict(case=case, reason="source did not report values"))
            continue
        if isinstance(actual, str):
            failures.append(dict(case=case, reason="native error", native=actual))
            continue
        if expected.get("status", 0) != 0:
            failures.append(
                dict(case=case, reason="source status", status=expected.get("status"))
            )
            continue
        bad = []
        for position, (a, b) in enumerate(zip(expected["values"], actual)):
            error = abs(a - b)
            worst = max(worst, error / max(abs(a), 1.0))
            if error > ATOL + RTOL * abs(a):
                bad.append(position)
        if bad:
            failures.append(
                dict(case=case, quantities=bad, reference=expected["values"], native=actual)
            )
            continue
        rows.append(
            [
                case,
                ",".join(f"{value:.17e}" for value in targets[index]),
                *expected["values"],
            ]
        )

    path = ROOT / "tests/fixtures/mut_solve.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "specification", "densities"])
        for row in rows:
            writer.writerow(
                [
                    row[0],
                    row[1],
                    ",".join(f"{value:.17e}" for value in row[2:]),
                ]
            )
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(request_list),
        generator="tools/rust_port/validate_mut_solve.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original lnf-T solver with the source model as the active equation of state; "
            "the target is the source log fugacity specification at the true state and the "
            "solve starts from a different guess"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/mut_solve_failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    print(
        json.dumps(
            {
                "attempts": len(request_list),
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
