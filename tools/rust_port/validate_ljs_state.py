"""Compare the native splined Lennard-Jones state against the exact-source
oracle.

The original state of the splined Lennard-Jones model follows from its
residual Helmholtz energy. A reference oracle asks the Python interface of the
fork for the molar volume and the compressibility factor of a state, and the
native model evaluates the pressure and the compressibility factor of the same
state.
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
QUANTITIES = ["volume", "compressibility"]
# The native model evaluates the pressure of the state, and the source the
# volume of the same state, so the compressibility factors are compared.
RTOL = 1.0e-6
ATOL = [1.0e-12, 1.0e-9]


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
        for temperature, pressure, phase in CASES:
            reference_values.append(
                oracle.call(
                    {
                        "mode": "LJS_STATE",
                        "temperature": temperature,
                        "pressure": pressure,
                        "feed": [1.0],
                        "phase": phase,
                    }
                )
            )

    lines = []
    for (temperature, pressure, phase), expected in zip(CASES, reference_values):
        if "values" not in expected:
            continue
        volume = expected["values"][0]
        lines.append(f"LJSSTATE|{temperature:.17e}|{volume:.17e}|1.0")
    native_values = native(lines)

    failures = []
    rows = []
    worst = 0.0
    index = 0
    for (temperature, pressure, phase), expected in zip(CASES, reference_values):
        label = f"{temperature:.17e}|{pressure:.17e}|{phase}"
        if "values" not in expected:
            failures.append(dict(case=label, reason="source did not report values", detail=expected))
            continue
        source_volume, source_compressibility = expected["values"]
        actual = native_values[index]
        index += 1
        if isinstance(actual, str):
            failures.append(dict(case=label, reason="native error", native=actual))
            continue
        native_pressure, native_compressibility = actual
        bad = []
        if abs(native_pressure - pressure) > 1.0e-6 * abs(pressure):
            bad.append("pressure")
        if abs(source_compressibility - native_compressibility) > ATOL[1] + RTOL * abs(
            source_compressibility
        ):
            bad.append("compressibility")
        if not bad:
            worst = max(
                worst,
                abs(source_compressibility - native_compressibility)
                    / max(abs(source_compressibility), 1.0e-30),
            )
        if bad:
            failures.append(
                dict(
                    case=label,
                    quantities=bad,
                    source=[source_volume, source_compressibility],
                    native=[native_pressure, native_compressibility],
                )
            )
            continue
        rows.append([label, source_volume, source_compressibility, native_pressure])

    path = ROOT / "tests/fixtures/ljs_state.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "volume", "compressibility", "native_pressure"])
        for row in rows:
            writer.writerow(
                [
                    row[0],
                    f"{row[1]:.17e}",
                    f"{row[2]:.17e}",
                    f"{row[3]:.17e}",
                ]
            )
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_ljs_state.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerances=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original state of the splined Lennard-Jones model; the source reports the molar "
            "volume and the compressibility factor of the state at the requested temperature and "
            "pressure, and the native model evaluates the pressure and the compressibility factor "
            "of the same state"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/ljs_state_failures.json").write_text(
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
