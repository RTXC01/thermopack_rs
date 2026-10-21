"""Compare the native reduced interaction potential against the exact-source
oracle.

The original `ljs_potential_reduced` evaluates the truncated and shifted
Lennard-Jones potential of the splined model and the reduced force. A
reference probe calls the original routine.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the separations divided by the segment diameter.
CASES = [
    [0.9, 1.0, 1.1, 1.2],
    [1.05, 1.1, 1.2, 1.3, 1.4],
    [1.15, 1.25, 1.3, 1.35, 1.4, 1.45],
    [1.2, 1.3, 1.4, 1.5, 1.6, 1.7],
]
QUANTITIES = ["potential", "force"]
RTOL = 1.0e-12
ATOL = 1.0e-14


def run_native(lines):
    exe = ROOT / "target/release/examples" / ("parity.exe" if os.name == "nt" else "parity")
    result = subprocess.run(
        [str(exe)], input="\n".join(lines) + "\n", text=True, capture_output=True, check=True
    )
    rows = []
    for line in result.stdout.splitlines():
        rows.append(line if line.startswith("ERROR") else line.split(";"))
    assert len(rows) == len(lines), (len(rows), len(lines))
    return rows


def main(directory):
    lines = [
        "LJSPOT|" + ",".join(f"{value:.17e}" for value in separations)
        for separations in CASES
    ]
    native_values = run_native(lines)

    with ReferenceProcess(directory) as oracle:
        reference_values = []
        for separations in CASES:
            reference_values.append(
                oracle.call({"mode": "LJS_POTENTIAL", "separations": separations})
            )

    failures = []
    rows = []
    worst = 0.0
    for index, separations in enumerate(CASES):
        label = ",".join(f"{value:.17e}" for value in separations)
        expected = reference_values[index]
        actual = native_values[index]
        if "values" not in expected:
            failures.append(dict(case=label, reason="source did not report values", detail=expected))
            continue
        if isinstance(actual, str):
            failures.append(dict(case=label, reason="native error", native=actual))
            continue
        bad = []
        for position in range(len(separations)):
            source_potential = expected["values"][2 * position]
            source_force = expected["values"][2 * position + 1]
            native_potential, native_force = [
                float(value) for value in actual[position].split(",")
            ]
            for source, native, quantity in (
                (source_potential, native_potential, "potential"),
                (source_force, native_force, "force"),
            ):
                if abs(source - native) > ATOL + RTOL * abs(source):
                    bad.append((position, quantity, source, native))
                else:
                    worst = max(worst, abs(source - native) / max(abs(source), 1.0e-30))
        if bad:
            failures.append(dict(case=label, reason="potential", differences=bad[:6]))
            continue
        rows.append([label, *expected["values"]])

    path = ROOT / "tests/fixtures/ljs_potential.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "values"])
        for row in rows:
            writer.writerow([row[0], ",".join(f"{value:.17e}" for value in row[1:])])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_ljs_potential.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original ljs_potential_reduced; the source reports the potential divided by the "
            "well depth and the reduced force for every separation divided by the segment diameter"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/ljs_potential_failures.json").write_text(
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
