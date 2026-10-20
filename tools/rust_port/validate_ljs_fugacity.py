"""Compare the native splined Lennard-Jones fugacity against the exact-source
oracle.

The original reports the logarithm of the fugacity of a state through the
Python interface of the fork. A reference oracle asks for it, and the native
model evaluates the residual chemical potential of the single component from
the reduced residual Helmholtz energy.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the temperature in kelvin, the pressure in pascal and the phase.
CASES = [
    (100.0, 1.0e5, 1),
    (150.0, 1.0e5, 1),
    (150.0, 1.0e6, 1),
    (200.0, 1.0e6, 2),
    (200.0, 5.0e6, 1),
    (300.0, 1.0e7, 2),
    (120.0, 5.0e5, 1),
    (250.0, 2.0e7, 1),
]
QUANTITIES = ["log_fugacity"]
RTOL = 1.0e-6
ATOL = 1.0e-9


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
        volumes = []
        for temperature, pressure, phase in CASES:
            state = oracle.call(
                {
                    "mode": "LJS_STATE",
                    "temperature": temperature,
                    "pressure": pressure,
                    "feed": [1.0],
                    "phase": phase,
                }
            )
            volumes.append(state["values"][0] if "values" in state else None)
        reference_values = []
        for (temperature, _, _), volume in zip(CASES, volumes):
            if volume is None:
                reference_values.append({})
                continue
            result = oracle.call(
                {
                    "mode": "LJS_HELMHOLTZ",
                    "temperature": temperature,
                    "volume": volume,
                    "moles": 1.0,
                }
            )
            reference_values.append(result)

    lines = []
    for (temperature, pressure, phase), volume in zip(CASES, volumes):
        if volume is None:
            continue
        lines.append(f"LJSFUGACITY|{temperature:.17e}|{volume:.17e}|1.0")
    native_values = native(lines)

    failures = []
    rows = []
    worst = 0.0
    index = 0
    for (temperature, pressure, phase), expected, volume in zip(
        CASES, reference_values, volumes
    ):
        label = f"{temperature:.17e}|{pressure:.17e}|{phase}"
        if volume is None or "values" not in expected:
            failures.append(dict(case=label, reason="source did not report values", detail=expected))
            continue
        source = expected["values"][2]
        actual = native_values[index]
        index += 1
        if isinstance(actual, str):
            failures.append(dict(case=label, reason="native error", native=actual))
            continue
        native_log_fugacity = actual[0]
        if abs(source - native_log_fugacity) > ATOL + RTOL * abs(source):
            failures.append(
                dict(case=label, reason="log fugacity", source=source, native=native_log_fugacity)
            )
            continue
        worst = max(worst, abs(source - native_log_fugacity) / max(abs(source), 1.0))
        rows.append([label, source, actual[0], volume])

    path = ROOT / "tests/fixtures/ljs_fugacity.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "log_fugacity", "log_fugacity_coefficient", "volume"])
        for row in rows:
            writer.writerow(
                [row[0], f"{row[1]:.17e}", f"{row[2]:.17e}", f"{row[3]:.17e}"]
            )
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_ljs_fugacity.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original composition derivative of the reduced residual Helmholtz energy of the "
            "splined Lennard-Jones model; the source reports the derivative of the reduced energy "
            "with respect to the number of moles, and the native model reports the same "
            "derivative, which is the reduced residual chemical potential of the single component. "
            "The Python interface of the fork reports the logarithm of the fugacity coefficient "
            "of the model as the reduced residual energy itself, which is a different convention"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/ljs_fugacity_failures.json").write_text(
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
