"""Compare the native quadrature rules against the exact-source oracle.

The source offers the Gauss and the Gauss-Kronrod rules in seven orders and
integrates a function together with its differential. A reference probe calls
the original routines.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The quadratures of the source, with the identifier of the routine and the
# name of the native rule.
QUADRATURES = {
    "gauss10": 1,
    "kronrod21": 2,
    "gauss15": 3,
    "kronrod31": 4,
    "gauss20": 5,
    "kronrod41": 6,
    "gauss25": 7,
    "kronrod51": 8,
    "gauss30": 9,
    "kronrod61": 10,
    "gauss40": 11,
    "kronrod81": 12,
    "gauss50": 13,
    "kronrod101": 14,
}
# The integration limits of the test integrand of the probe.
LIMITS = (0.0, 1.0)
RTOL = 1.0e-13
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
        reference_values = {}
        for name, identifier in QUADRATURES.items():
            result = oracle.call(
                {
                    "mode": "QUADRATURE",
                    "quadrature": identifier,
                    "lower": LIMITS[0],
                    "upper": LIMITS[1],
                }
            )
            if "values" not in result:
                raise SystemExit(f"the source did not report the rule: {result}")
            reference_values[name] = result

    lines = [f"QUAD|{index}" for index in range(len(QUADRATURES))]
    native_values = native(lines)

    failures = []
    rows = []
    worst = 0.0
    for index, (name, identifier) in enumerate(QUADRATURES.items()):
        expected = reference_values[name]
        rule = expected["values"]
        count = len(rule) // 2
        positions = rule[:count]
        weights = rule[count:]
        actual = native_values[index]
        if isinstance(actual, str):
            failures.append(dict(case=name, reason="native error", native=actual))
            continue
        if len(actual) != 2 * count + 2:
            failures.append(
                dict(case=name, reason="point count", source=count, native=len(actual) // 2)
            )
            continue
        bad = []
        # The source leaves the last node and the last weight of the Gauss
        # rules with an odd number of points unset, and the value depends on the
        # stack of the call, so the comparison stops before them.
        odd_gauss = name in ("gauss15", "gauss25")
        limit = count - 1 if odd_gauss else count
        for position in range(limit):
            for source, value, quantity in (
                (positions[position], actual[position], "position"),
                (weights[position], actual[count + position], "weight"),
            ):
                if abs(source - value) > ATOL + RTOL * abs(source):
                    bad.append((position, quantity, source, value))
                else:
                    worst = max(worst, abs(source - value) / max(abs(source), 1.0e-30))
        # The integral of the test integrand of the probe and the native.
        integral = expected["integral"]
        differential = expected["differential"]
        native_integral = actual[2 * count]
        native_differential = actual[2 * count + 1]
        if not odd_gauss:
            for source, value, quantity in (
                (integral, native_integral, "integral"),
                (differential, native_differential, "differential"),
            ):
                if abs(source - value) > ATOL + RTOL * abs(source):
                    bad.append((count, quantity, source, value))
                else:
                    worst = max(worst, abs(source - value) / max(abs(source), 1.0e-30))
        if bad:
            failures.append(dict(case=name, reason="rule", differences=bad[:6]))
            continue
        rows.append([name, *positions, *weights, integral, differential])

    path = ROOT / "tests/fixtures/quadrature.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "positions", "weights", "integral", "differential"])
        for row in rows:
            half = (len(row) - 3) // 2
            writer.writerow(
                [
                    row[0],
                    ",".join(f"{value:.17e}" for value in row[1 : 1 + half]),
                    ",".join(f"{value:.17e}" for value in row[1 + half : 1 + 2 * half]),
                    f"{row[-2]:.17e}",
                    f"{row[-1]:.17e}",
                ]
            )
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(QUADRATURES),
        generator="tools/rust_port/validate_quadrature.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original get_quadrature_positions and get_quadrature_weights; the nodes and the "
            "weights of every rule of the source on the unit interval, and the integral of the "
            "square of the integration variable with its differential over the unit interval"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/quadrature_failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    print(
        json.dumps(
            {
                "attempts": len(QUADRATURES),
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
