"""Compare the native splined Lennard-Jones caloric state against the
exact-source oracle.

The original reports the enthalpy and the entropy of a state through the
Python interface of the fork. A reference oracle asks for them, and the native
model evaluates the residual part of the caloric properties with the
temperature derivative of the reduced residual Helmholtz energy and the ideal
part with the constant heat capacity of the original.
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
QUANTITIES = ["enthalpy", "entropy"]
# The enthalpy carries the reference energy of the argon record and the
# entropy the reference entropy, so the absolute tolerances are the ones of
# the source reference.
RTOL = 1.0e-8
ATOL = [1.0e-6, 1.0e-9]


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
            if volume_missing(volumes[-1]):
                reference_values.append({})
                continue
            result = oracle.call(
                {
                    "mode": "LJS_CALORIC",
                    "temperature": temperature,
                    "pressure": pressure,
                    "feed": [1.0],
                    "phase": phase,
                }
            )
            reference_values.append(result)

    lines = []
    for (temperature, _, _), volume in zip(CASES, volumes):
        if volume_missing(volume):
            continue
        lines.append(f"LJSCALORIC|{temperature:.17e}|{volume:.17e}|1.0")
    native_values = native(lines)

    failures = []
    rows = []
    worst = 0.0
    index = 0
    for (temperature, pressure, phase), expected, volume in zip(
        CASES, reference_values, volumes
    ):
        label = f"{temperature:.17e}|{pressure:.17e}|{phase}"
        if volume_missing(volume) or "values" not in expected:
            failures.append(dict(case=label, reason="source did not report values", detail=expected))
            continue
        source_h, source_s = expected["values"]
        actual = native_values[index]
        index += 1
        if isinstance(actual, str):
            failures.append(dict(case=label, reason="native error", native=actual))
            continue
        native_h, native_s = actual
        bad = []
        if abs(source_h - native_h) > ATOL[0] + RTOL * abs(source_h):
            bad.append("enthalpy")
        if abs(source_s - native_s) > ATOL[1] + RTOL * abs(source_s):
            bad.append("entropy")
        if not bad:
            worst = max(worst, abs(source_h - native_h) / max(abs(source_h), 1.0))
            worst = max(worst, abs(source_s - native_s) / max(abs(source_s), 1.0))
        if bad:
            failures.append(
                dict(
                    case=label,
                    quantities=bad,
                    source=[source_h, source_s],
                    native=[native_h, native_s],
                )
            )
            continue
        rows.append([label, source_h, source_s, volume])

    path = ROOT / "tests/fixtures/ljs_caloric.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *QUANTITIES, "volume"])
        for row in rows:
            writer.writerow(
                [row[0], f"{row[1]:.17e}", f"{row[2]:.17e}", f"{row[3]:.17e}"]
            )
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_ljs_caloric.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerances=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original state of the splined Lennard-Jones model; the source reports the molar "
            "enthalpy and the molar entropy of the state, and the native model reports the same "
            "quantities with the residual part from the temperature derivative of the reduced "
            "residual Helmholtz energy and the ideal part from the constant heat capacity of the "
            "original"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/ljs_caloric_failures.json").write_text(
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


def volume_missing(volume):
    return volume is None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler-directory")
    main(parser.parse_args().compiler_directory)
