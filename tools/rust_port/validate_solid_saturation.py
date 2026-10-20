"""Compare the native solid saturation line against the exact-source oracle.

The original `solidPointOnEnvelope` locates the state where a solid phase and
the fluid phase of the equation of state share the fugacity. A reference probe
installs the component in the source globals, initializes the source dry ice
model and calls the original routine. The native model solves the same
equilibrium with a Newton iteration in the logarithm of the pressure or the
temperature.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The specification prescribes the temperature or the pressure, and the other
# variable is solved. The third value is the initial guess of the solved
# variable, close enough to the line for both solvers.
CASES = [
    ("T", 216.592, 517950.0),
    ("T", 216.600, 550000.0),
    ("T", 216.610, 570000.0),
    ("P", 517950.0, 216.692),
    ("P", 600000.0, 216.692),
    ("P", 700000.0, 216.625),
]
# The pressure is solved in a temperature specification and the temperature is
# solved in a pressure specification.
PRESSURE_RTOL = 1e-6
TEMPERATURE_ATOL = 1e-8
TEMPERATURE_RTOL = 1e-9
# The prescribed variable is returned unchanged by both solvers.
PRESCRIBED_RTOL = 1e-12


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
    for kind, value, guess in CASES:
        if kind == "T":
            temperature, pressure, specification = value, guess, f"T:{value:.17e}"
        else:
            temperature, pressure, specification = guess, value, f"P:{value:.17e}"
        lines.append(
            "|".join(
                [
                    "SOLIDSAT",
                    "CO2",
                    f"{temperature:.17e}",
                    f"{pressure:.17e}",
                    "1",
                    specification,
                ]
            )
        )
    native_values = native(lines)

    with ReferenceProcess(directory) as oracle:
        reference_values = []
        for kind, value, guess in CASES:
            if kind == "T":
                temperature, pressure, ln_spec = value, guess, math.log(value)
            else:
                temperature, pressure, ln_spec = guess, value, math.log(value)
            reference_values.append(
                oracle.call(
                    {
                        "mode": "SOLID_SATURATION",
                        "component": "CO2",
                        "temperature": temperature,
                        "pressure": pressure,
                        "phase": 1,
                        "specification": 1 if kind == "T" else 2,
                        "ln_spec": ln_spec,
                    }
                )
            )

    failures = []
    rows = []
    worst = 0.0
    for index, (kind, value, guess) in enumerate(CASES):
        case = f"{kind}:{value:.17e}"
        expected = reference_values[index]
        actual = native_values[index]
        if "values" not in expected:
            failures.append(dict(case=case, reason="source did not report values", detail=expected))
            continue
        if expected.get("status", 0) != 0:
            failures.append(dict(case=case, reason="source status", status=expected.get("status")))
            continue
        if isinstance(actual, str):
            failures.append(dict(case=case, reason="native error", native=actual))
            continue
        source_t, source_p = expected["values"][0], expected["values"][1]
        native_t, native_p = actual[0], actual[1]
        if kind == "T":
            error = abs(source_p - native_p) / abs(source_p)
            prescribed = abs(source_t - native_t) / abs(source_t)
            tolerance = PRESSURE_RTOL
        else:
            error = abs(source_t - native_t) / abs(source_t)
            prescribed = abs(source_p - native_p) / abs(source_p)
            tolerance = TEMPERATURE_ATOL / abs(source_t) + TEMPERATURE_RTOL
        if prescribed > PRESCRIBED_RTOL:
            failures.append(
                dict(case=case, reason="prescribed variable moved", source=[source_t, source_p],
                     native=[native_t, native_p])
            )
            continue
        if error > tolerance:
            failures.append(
                dict(case=case, reason="solved variable", source=[source_t, source_p],
                     native=[native_t, native_p], relative_difference=error)
            )
            continue
        worst = max(worst, error)
        rows.append([case, kind, value, guess, source_t, source_p])

    path = ROOT / "tests/fixtures/solid_saturation.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "specification", "specified", "guess", "temperature", "pressure"])
        for row in rows:
            writer.writerow(
                [
                    row[0],
                    row[1],
                    f"{row[2]:.17e}",
                    f"{row[3]:.17e}",
                    f"{row[4]:.17e}",
                    f"{row[5]:.17e}",
                ]
            )
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_solid_saturation.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        pressure_relative_tolerance=PRESSURE_RTOL,
        temperature_absolute_tolerance=TEMPERATURE_ATOL,
        temperature_relative_tolerance=TEMPERATURE_RTOL,
        worst_relative_difference=worst,
        convention=(
            "the original solidPointOnEnvelope with the source dry ice model; the "
            "specification prescribes the temperature or the pressure and the other variable is "
            "solved"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/solid_saturation_failures.json").write_text(
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
