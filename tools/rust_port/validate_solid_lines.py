"""Compare the native melting and sublimation line correlations against the
exact-source Fortran oracle.

The original routines live in `src/solid_saturation.f90` and read their
parameters from `src/solid_correlation_datadb.f90`. A reference probe installs
the requested component in the source globals, calls the source routine and
restores the globals afterwards.
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

MELTING_COMPONENTS = [
    "NH3", "AR", "CO2", "CO", "CYCLOHEX", "C2", "HE", "H2", "IC4", "IC5", "KR",
    "C1", "MEOH", "NE", "N2", "O2", "P-H2", "C3", "PRLN", "F6S", "XE", "NC4",
    "NC5",
]
SUBLIMATION_COMPONENTS = [
    "NH3", "AR", "CO2", "CO", "C2", "H2", "H2S", "KR", "C1", "NE", "N2", "N2O",
    "O2", "P-H2", "F6S", "XE",
]
# The melting lines are sampled up to a temperature above the triple point and
# the sublimation lines below it.
TEMPERATURE_LIMITS = {
    "MELT": [200.0, 260.0, 300.0, 400.0],
    "SUBL": [50.0, 100.0, 150.0, 200.0],
}
POINTS = 7
# The pressure tolerance covers the absolute and the relative contribution.
# Some melting curves subtract two nearly equal numbers, for example the
# krypton line at its triple point, which amplifies the rounding of the
# exponentiation by a factor of about three thousand. The thin link time
# optimization can also contract a multiplication and an addition into a fused
# operation in one binary and not in another, so the tolerance allows for
# that rounding difference.
RTOL = 1e-10
ATOL = 1e-4


def cases():
    for kind, components in (("MELT", MELTING_COMPONENTS), ("SUBL", SUBLIMATION_COMPONENTS)):
        for component in components:
            for limit in TEMPERATURE_LIMITS[kind]:
                yield kind, component, limit


def native(lines):
    if not lines:
        return []
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
    request_list = list(cases())
    lines = []
    for kind, component, limit in request_list:
        lines.append(
            "|".join(
                [
                    "SOLID",
                    component,
                    f"{limit:.17e}",
                    str(POINTS),
                    "0",
                    kind,
                ]
            )
        )
    native_values = native(lines)

    reference_values = []
    with ReferenceProcess(directory) as oracle:
        for kind, component, limit in request_list:
            result = oracle.call(
                {
                    "mode": "SOLID_LINE",
                    "kind": 0 if kind == "MELT" else 1,
                    "component": component,
                    "temperature": limit,
                    "points": POINTS,
                }
            )
            reference_values.append(result.get("values"))

    failures = []
    rows = []
    status_rows = []
    worst = {}
    for index, (kind, component, limit) in enumerate(request_list):
        case = "|".join([kind, component, f"{limit:.17e}", str(POINTS)])
        expected = reference_values[index]
        actual = native_values[index]
        if expected is None:
            if not isinstance(actual, str):
                failures.append(dict(case=case, reason="native value for terminated source"))
            status_rows.append([case, 1])
            continue
        if isinstance(actual, str):
            failures.append(dict(case=case, reason="native error", native=actual))
            status_rows.append([case, 2])
            continue
        if len(expected) != len(actual):
            failures.append(dict(case=case, reason="length", native=len(actual), reference=len(expected)))
            status_rows.append([case, 3])
            continue
        bad = []
        for position, (a, b) in enumerate(zip(expected, actual)):
            if math.isnan(a) and math.isnan(b):
                # Both sides reproduce the same undefined correlation value.
                continue
            if not math.isfinite(b):
                bad.append(position)
                continue
            error = abs(a - b)
            if error > ATOL + RTOL * max(abs(a), abs(b)):
                bad.append(position)
            if error > worst.get("absolute", 0.0):
                worst["absolute"] = error
                worst["case"] = f"{case} quantity {position}"
        if bad:
            failures.append(dict(case=case, quantities=bad, reference=expected, native=actual))
            status_rows.append([case, 3])
            continue
        status_rows.append([case, 0])
        rows.append([case, *expected])

    path = ROOT / "tests/fixtures/solid_lines.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "temperatures", "pressures"])
        for row in rows:
            half = len(row) - 1
            writer.writerow(
                [
                    row[0],
                    ",".join(f"{value:.17e}" for value in row[1 : 1 + half // 2]),
                    ",".join(f"{value:.17e}" for value in row[1 + half // 2 :]),
                ]
            )
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        generator="tools/rust_port/validate_solid_lines.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=["temperature", "pressure"],
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_difference=worst,
        convention=(
            "the original melting_pressure_correlation and sublimation_pressure_correlation "
            "without scaling to the equation of state; the temperatures and pressures are the "
            "source grid and the source correlation values in Pa"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    status_path = ROOT / "tests/fixtures/solid_lines_status.csv"
    with status_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", "source_status"])
        writer.writerows(status_rows)
    status_metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    status_metadata.update(
        cases=len(status_rows),
        generator="tools/rust_port/validate_solid_lines.py",
        fixture_sha256=hashlib.sha256(status_path.read_bytes()).hexdigest(),
        quantities=["source_status"],
        convention=(
            "0 compared and matching, 1 the original terminated, 2 the native side reported an "
            "error, 3 a numerical mismatch"
        ),
    )
    status_path.with_suffix(".metadata.json").write_text(json.dumps(status_metadata, indent=2) + "\n")
    (ROOT / ".reference/solid_lines_failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    print(
        json.dumps(
            {
                "cases": len(request_list),
                "compared": len(rows),
                "source_terminations": sum(1 for row in status_rows if row[1] == 1),
                "native_errors": sum(1 for row in status_rows if row[1] == 2),
                "mismatches": sum(1 for row in status_rows if row[1] == 3),
                "failures": len(failures),
                "worst_difference": worst,
            },
            indent=2,
        )
    )
    if failures:
        print(json.dumps(failures[:8], indent=2))
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler-directory")
    main(parser.parse_args().compiler_directory)
