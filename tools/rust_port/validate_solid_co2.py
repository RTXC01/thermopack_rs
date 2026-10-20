"""Compare the native solid carbon dioxide model against the exact-source
Fortran oracle.

The original module `src/co2_gibbs.f90` keeps the reference levels of the
model in module variables. A reference probe shifts those levels with the
source `sco2init` routine and then evaluates the model, and the native model
performs the same shift.
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

# The reference levels of the liquid phase at the triple point. The values
# are the ones of the source reference state of carbon dioxide.
LIQUID_ENTROPY = 205.147591034776
LIQUID_GIBBS = -18000.0
TRIPLE_TEMPERATURE = 216.592
TRIPLE_PRESSURE = 517950.0
# States below and above the triple point temperature.
STATES = [
    (150.0, 1.0e5),
    (180.0, 2.0e5),
    (200.0, 3.0e5),
    (216.592, 5.1795e5),
    (230.0, 1.0e6),
    (250.0, 5.0e6),
    (270.0, 2.0e7),
]
QUANTITIES = [
    "gibbs",
    "dg_dp",
    "d2g_dp2",
    "dg_dt",
    "d2g_dt2",
    "d2g_dtdp",
    "specific_volume",
    "entropy",
    "enthalpy",
    "internal_energy",
    "helmholtz",
    "heat_capacity",
    "speed_of_sound",
]
# The speed of sound is a square root of a small difference and therefore
# needs the loosest tolerance.
ATOL = [1e-6, 1e-16, 1e-26, 1e-9, 1e-14, 1e-17, 1e-16, 1e-9, 1e-6, 1e-6, 1e-6, 1e-9, 1e-6]
RTOL = 1e-11


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
    for temperature, pressure in STATES:
        lines.append(
            "|".join(
                [
                    "SOLIDCO2",
                    f"{LIQUID_ENTROPY:.17e}",
                    f"{LIQUID_GIBBS:.17e}",
                    f"{TRIPLE_TEMPERATURE:.17e}",
                    f"{TRIPLE_PRESSURE:.17e}",
                    f"{temperature:.17e}",
                    f"{pressure:.17e}",
                ]
            )
        )
    native_values = native(lines)

    # The source module keeps its reference levels, so the shift is applied
    # once and the states are then evaluated with the shifted model.
    with ReferenceProcess(directory) as oracle:
        oracle.call(
            {
                "mode": "SOLID_CO2_INIT",
                "liquid_entropy": LIQUID_ENTROPY,
                "liquid_gibbs": LIQUID_GIBBS,
                "temperature": TRIPLE_TEMPERATURE,
                "pressure": TRIPLE_PRESSURE,
            }
        )
        reference_values = []
        for temperature, pressure in STATES:
            result = oracle.call(
                {
                    "mode": "SOLID_CO2",
                    "temperature": temperature,
                    "pressure": pressure,
                }
            )
            reference_values.append(result.get("values"))

    failures = []
    rows = []
    for index, (temperature, pressure) in enumerate(STATES):
        case = "|".join(
            [
                "SOLIDCO2",
                f"{LIQUID_ENTROPY:.17e}",
                f"{LIQUID_GIBBS:.17e}",
                f"{TRIPLE_TEMPERATURE:.17e}",
                f"{TRIPLE_PRESSURE:.17e}",
                f"{temperature:.17e}",
                f"{pressure:.17e}",
            ]
        )
        expected = reference_values[index]
        actual = native_values[index]
        if expected is None or isinstance(actual, str):
            failures.append(dict(case=case, reason="missing values", reference=expected, native=actual))
            continue
        bad = []
        for position, (a, b) in enumerate(zip(expected, actual)):
            atol = ATOL[position] if position < len(ATOL) else ATOL[-1]
            if not math.isfinite(b) or abs(a - b) > atol + RTOL * max(abs(a), abs(b)):
                bad.append(QUANTITIES[position])
        if bad:
            failures.append(dict(case=case, quantities=bad, reference=expected, native=actual))
            continue
        rows.append([case, *expected])

    path = ROOT / "tests/fixtures/solid_co2.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *QUANTITIES])
        writer.writerows(rows)
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        generator="tools/rust_port/validate_solid_co2.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerances=ATOL,
        convention=(
            "the original solid carbon dioxide Gibbs energy model with the source reference "
            "shift of sco2init; the source reports J per mole and the native model uses the "
            "same molar convention"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/solid_co2_failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    print(json.dumps({"cases": len(STATES), "compared": len(rows), "failures": len(failures)}, indent=2))
    if failures:
        print(json.dumps(failures[:4], indent=2))
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler-directory")
    main(parser.parse_args().compiler_directory)
