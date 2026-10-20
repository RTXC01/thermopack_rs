"""Compare the native Joule-Thomson inversion against the exact-source oracle.

The original maps the Joule-Thomson inversion curve and solves the inversion
condition for the temperature at every molar volume of the curve. A reference
oracle calls the Python interface of the fork, which reaches the original
routine, and the native model solves the same condition at the same volumes.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The fluids and the number of points of the inversion curve. The original
# fails to map the curve of methane and of carbon dioxide, and the nitrogen
# curve is the reference.
CASES = [
    ("PR", ["N2"], 40),
    ("SRK", ["N2"], 40),
]
# The relative step of the initial guess of the native solve, so that the
# native solve is a solution and not a repetition of the source value.
GUESS_STEP = 1.02
QUANTITIES = ["temperature", "pressure"]
# The source solves the condition with a tolerance of 1e-10.
ATOL = [1.0e-6, 1.0e-4]
RTOL = 1.0e-7


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
        curves = []
        for eos, components, points in CASES:
            result = oracle.call(
                {
                    "mode": "JTI",
                    "eos": eos,
                    "components": components,
                    "feed": [1.0] * len(components),
                    "nmax": points,
                }
            )
            if "values" not in result:
                raise SystemExit(f"the source did not map the curve: {result}")
            curves.append(result["values"])

    # The first point of the curve is the zero pressure solution of the
    # original, which comes from the virial coefficients and not from the
    # inversion condition, so the comparison starts at the second point.
    lines = []
    labels = []
    for (eos, components, points), curve in zip(CASES, curves):
        temperature, pressure, volume = curve
        for index in range(1, len(temperature)):
            guess = temperature[index] * GUESS_STEP
            lines.append(
                "|".join(
                    [
                        "JTI",
                        eos,
                        ",".join(components),
                        f"{volume[index]:.17e}",
                        f"{guess:.17e}",
                    ]
                )
            )
            labels.append((eos, components, index))
    native_values = native(lines)

    failures = []
    rows = []
    worst = 0.0
    for index, ((eos, components, point), actual) in enumerate(zip(labels, native_values)):
        curve = curves[[case[0] for case in CASES].index(eos)]
        source_temperature = curve[0][point]
        source_pressure = curve[1][point]
        label = f"{eos}|{','.join(components)}|{point}"
        if isinstance(actual, str):
            failures.append(dict(case=label, reason="native error", native=actual))
            continue
        native_temperature, native_pressure = actual[0], actual[2]
        bad = []
        for position, (source, native_value) in enumerate(
            zip((source_temperature, source_pressure), (native_temperature, native_pressure))
        ):
            if abs(source - native_value) > ATOL[position] + RTOL * abs(source):
                bad.append(QUANTITIES[position])
            else:
                worst = max(worst, abs(source - native_value) / max(abs(source), 1.0))
        if bad:
            failures.append(
                dict(
                    case=label,
                    quantities=bad,
                    source=[source_temperature, source_pressure],
                    native=[native_temperature, native_pressure],
                )
            )
            continue
        rows.append([label, source_temperature, source_pressure])

    path = ROOT / "tests/fixtures/joule_thompson.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "temperature", "pressure"])
        for row in rows:
            writer.writerow([row[0], f"{row[1]:.17e}", f"{row[2]:.17e}"])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(lines),
        generator="tools/rust_port/validate_joule_thompson.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerances=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original map_jt_inversion; the source solves the inversion condition for the "
            "temperature at every molar volume of the curve, and the native model solves the same "
            "condition at the same volumes with a perturbed initial guess"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/joule_thompson_failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    print(
        json.dumps(
            {
                "attempts": len(lines),
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
