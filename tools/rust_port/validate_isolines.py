"""Compare the native iso-lines against the exact-source oracle.

The original `isotherm` and `isobar` map a curve of the temperature-pressure
plane with the properties along it. A reference oracle calls the Python
interface of the fork, which reaches the original routines.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the line, the equation of state, the components, the composition,
# the prescribed variable, the limits of the other variable and the number of
# points. The isenthalp and the isentrope add the temperature limits.
CASES = [
    ("isotherm", "PR", ["C1"], [1.0], 150.0, 1.0e5, 1.0e7, 10, None, None),
    ("isotherm", "PR", ["C1"], [1.0], 250.0, 1.0e5, 1.0e7, 10, None, None),
    ("isotherm", "PR", ["CO2"], [1.0], 250.0, 5.0e5, 2.0e7, 12, None, None),
    ("isotherm", "SRK", ["CO2"], [1.0], 280.0, 1.0e6, 1.5e7, 8, None, None),
    ("isobar", "PR", ["C1"], [1.0], 1.0e6, 120.0, 190.0, 8, None, None),
    ("isobar", "PR", ["CO2"], [1.0], 5.0e6, 240.0, 300.0, 8, None, None),
    ("isotherm", "PR", ["C1", "C3"], [0.5, 0.5], 250.0, 1.0e5, 1.0e7, 8, None, None),
    ("isotherm", "PR", ["C1", "CO2"], [0.7, 0.3], 220.0, 1.0e5, 2.0e7, 8, None, None),
    ("isobar", "PR", ["C1", "C3"], [0.5, 0.5], 2.0e6, 180.0, 280.0, 8, None, None),
    ("isenthalp", "PR", ["C1"], [1.0], -78000.0, 1.0e5, 1.0e7, 8, 150.0, 250.0),
    ("isenthalp", "PR", ["CO2"], [1.0], -395800.0, 5.0e5, 1.5e7, 8, 220.0, 300.0),
    ("isentrope", "PR", ["C1"], [1.0], 165.0, 1.0e5, 5.0e5, 8, 150.0, 250.0),
    ("isentrope", "PR", ["C1", "C3"], [0.5, 0.5], 160.0, 1.0e5, 1.0e7, 8, 180.0, 280.0),
]
# The quantities that are compared: the temperature, the pressure, the molar
# volume, the molar entropy and the molar enthalpy of every point.
QUANTITIES = ["temperature", "pressure", "volume", "entropy", "enthalpy"]
RTOL = 1e-6
ATOL = 1e-9


def native(lines):
    exe = ROOT / "target/release/examples" / ("parity.exe" if os.name == "nt" else "parity")
    result = subprocess.run(
        [str(exe)], input="\n".join(lines) + "\n", text=True, capture_output=True, check=True
    )
    rows = []
    for line in result.stdout.splitlines():
        rows.append(line if line.startswith("ERROR") else line.split(";"))
    assert len(rows) == len(lines), (len(rows), len(lines))
    return rows


def close(source, native):
    if math.isnan(source) or math.isnan(native):
        return math.isnan(source) and math.isnan(native)
    return abs(source - native) <= ATOL + RTOL * abs(source)


def main(directory):
    lines = []
    for case in CASES:
        line, eos, components, composition, value, minimum, maximum, points, tmin, tmax = case
        fields = [
            "ISOLINE",
            eos,
            ",".join(components),
            ",".join(f"{item:.17e}" for item in composition),
            f"{value:.17e}",
            f"{minimum:.17e}",
            f"{maximum:.17e}",
        ]
        if tmin is None:
            fields.append(str(points))
        else:
            fields.extend([f"{tmin:.17e}", f"{tmax:.17e}", str(points)])
        fields.append(line)
        lines.append("|".join(fields))
    native_values = native(lines)

    with ReferenceProcess(directory) as oracle:
        reference_values = []
        for case in CASES:
            line, eos, components, composition, value, minimum, maximum, points, tmin, tmax = case
            request = {
                "mode": "ISOLINE",
                "line": line,
                "eos": eos,
                "components": components,
                "feed": composition,
                "points": points,
                "minimum": minimum,
                "maximum": maximum,
            }
            if line == "isobar":
                request["pressure"] = value
            elif line == "isenthalp":
                request["enthalpy"] = value
                request["tmin"] = tmin
                request["tmax"] = tmax
            elif line == "isentrope":
                request["entropy"] = value
                request["tmin"] = tmin
                request["tmax"] = tmax
            else:
                request["temperature"] = value
            reference_values.append(oracle.call(request))

    failures = []
    rows = []
    worst = 0.0
    for index, case in enumerate(CASES):
        line, eos, components, composition, value, minimum, maximum, points, tmin, tmax = case
        label = (
            f"{line}|{eos}|{','.join(components)}|"
            f"{','.join(f'{item:.17e}' for item in composition)}|"
            f"{value:.17e}|{minimum:.17e}|{maximum:.17e}|"
            f"{'none' if tmin is None else f'{tmin:.17e}'}|"
            f"{'none' if tmax is None else f'{tmax:.17e}'}|{points}"
        )
        expected = reference_values[index]
        actual = native_values[index]
        if "values" not in expected:
            failures.append(dict(case=label, reason="source did not report values", detail=expected))
            continue
        if isinstance(actual, str):
            failures.append(dict(case=label, reason="native error", native=actual))
            continue
        source_arrays = expected["values"]
        if len(source_arrays[0]) != len(actual):
            failures.append(
                dict(case=label, reason="point count", source=len(source_arrays[0]), native=len(actual))
            )
            continue
        bad = []
        for point in range(len(actual)):
            native_point = [float(value) for value in actual[point].split(",")]
            # The oracle reports the temperature, the pressure, the volume, the
            # entropy and the enthalpy of every point, and the native point
            # holds the same five quantities.
            source_point = [float(values[point]) for values in source_arrays]
            for position, (source, native_value) in enumerate(zip(source_point, native_point)):
                if not close(source, native_value):
                    bad.append((point, QUANTITIES[position], source, native_value))
                elif not math.isnan(source):
                    worst = max(worst, abs(source - native_value) / max(abs(source), 1.0))
        if bad:
            failures.append(dict(case=label, reason="property", differences=bad[:6]))
            continue
        rows.append([label, source_arrays])

    path = ROOT / "tests/fixtures/isolines.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(
            [
                "case",
                "point",
                "temperature",
                "pressure",
                "volume",
                "entropy",
                "enthalpy",
            ]
        )
        for label, arrays in rows:
            for point in range(len(arrays[0])):
                writer.writerow(
                    [
                        label,
                        str(point),
                        f"{float(arrays[0][point]):.17e}",
                        f"{float(arrays[1][point]):.17e}",
                        f"{float(arrays[2][point]):.17e}",
                        f"{float(arrays[3][point]):.17e}",
                        f"{float(arrays[4][point]):.17e}",
                    ]
                )
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=sum(len(arrays[0]) for _, arrays in rows),
        curves=len(rows),
        attempted_curves=len(CASES),
        generator="tools/rust_port/validate_isolines.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original isotherm and isobar of src/isolines.f90; the fixture holds the "
            "prescribed variable of every point, which is the pressure of an isotherm and the "
            "temperature of an isobar, and the source reports the molar volume, the molar entropy "
            "and the molar enthalpy of the same points"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/isolines_failures.json").write_text(
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
