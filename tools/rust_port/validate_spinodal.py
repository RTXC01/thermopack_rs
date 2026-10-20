"""Compare the native stability limit points against the exact-source oracle.

The original `initial_stab_limit_point` locates the point of the meta
stability limit at a given pressure, where the first volume derivative of the
pressure vanishes. A reference oracle calls the Python interface of the fork,
which reaches the original routine.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the equation of state, the component, the pressure and the phase
# of the limit, the liquid branch or the vapor branch.
CASES = [
    ("PR", ["C1"], 1.0e5, 1),
    ("PR", ["C1"], 1.0e6, 1),
    ("PR", ["C1"], 2.0e6, 1),
    ("PR", ["C1"], 1.0e5, 2),
    ("PR", ["C1"], 1.0e6, 2),
    ("PR", ["C1"], 3.0e6, 2),
    ("PR", ["CO2"], 5.0e5, 1),
    ("PR", ["CO2"], 5.0e6, 1),
    ("PR", ["CO2"], 6.0e6, 2),
    ("PR", ["CO2"], 1.0e6, 2),
    ("SRK", ["C1"], 1.0e6, 1),
    ("SRK", ["C1"], 2.0e6, 2),
]
QUANTITIES = ["temperature", "volume"]
# The source solves the limit with a tolerance of 1e-8 on the volume.
ATOL = [1.0e-6, 1.0e-10]
RTOL = 1.0e-8


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
    for eos, components, pressure, phase in CASES:
        lines.append(
            "|".join(
                [
                    "SPINODAL",
                    eos,
                    ",".join(components),
                    f"{pressure:.17e}",
                    str(phase),
                ]
            )
        )
    native_values = native(lines)

    with ReferenceProcess(directory) as oracle:
        reference_values = []
        for eos, components, pressure, phase in CASES:
            reference_values.append(
                oracle.call(
                    {
                        "mode": "SPINODAL",
                        "eos": eos,
                        "components": components,
                        "feed": [1.0] * len(components),
                        "pressure": pressure,
                        "phase": phase,
                    }
                )
            )

    failures = []
    rows = []
    worst = 0.0
    for index, case in enumerate(CASES):
        eos, components, pressure, phase = case
        label = f"{eos}|{','.join(components)}|{pressure:.17e}|{phase}"
        expected = reference_values[index]
        actual = native_values[index]
        if "values" not in expected:
            failures.append(dict(case=label, reason="source did not report values", detail=expected))
            continue
        if isinstance(actual, str):
            failures.append(dict(case=label, reason="native error", native=actual))
            continue
        source_temperature, source_volume = expected["values"][0], expected["values"][1]
        native_temperature, native_volume = actual[0], actual[1]
        bad = []
        for position, (source, native_value) in enumerate(
            zip((source_temperature, source_volume), (native_temperature, native_volume))
        ):
            if abs(source - native_value) > ATOL[position] + RTOL * abs(source):
                bad.append(QUANTITIES[position])
            else:
                worst = max(worst, abs(source - native_value) / max(abs(source), 1.0e-30))
        if bad:
            failures.append(
                dict(
                    case=label,
                    quantities=bad,
                    source=[source_temperature, source_volume],
                    native=[native_temperature, native_volume],
                )
            )
            continue
        rows.append([label, source_temperature, source_volume])

    path = ROOT / "tests/fixtures/spinodal.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "temperature", "volume"])
        for row in rows:
            writer.writerow([row[0], f"{row[1]:.17e}", f"{row[2]:.17e}"])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_spinodal.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerances=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original initial_stab_limit_point, which the Python interface exposes as "
            "spinodal_point; the point lies on the meta stability limit, where the first volume "
            "derivative of the pressure vanishes, and the source reports the temperature and the "
            "molar volume of the point"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/spinodal_failures.json").write_text(
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
