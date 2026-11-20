"""Compare the native UNIFAC model against the exact-source Fortran oracle.

The original `Ge_UNIFAC_GH_SG` routine is exercised through a reference
probe that installs the requested component list in the source globals,
calls `init_unifac` and the source excess Gibbs routine, and restores the
globals afterwards. The native side evaluates the same composition with
`thermopack_rs::Unifac`.
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

COMPONENTS = [
    "C2", "C3", "IC4", "NC4", "NC5", "IC5", "NC6", "NC7", "NC8", "NC9",
    "NC10", "NC11", "BENZENE", "TOLUENE", "H2O", "N2", "CO2", "H2S", "C1",
    "ACETONE",
]
VARIANTS = ["UNIFAC", "UMR", "VTPR"]
TEMPERATURES = [200.0, 298.15, 450.0]
# The source reports the excess Gibbs energy per kilomole: its gas constant is
# 1000*R. The native model works per mole, so the native values are scaled
# to the source convention before the comparison.
SOURCE_UNIT_SCALE = 1000.0
# Absolute tolerances per reported quantity: excess Gibbs energy, first and
# second temperature derivative, composition derivatives and the mixed and
# second composition derivatives.
ATOL = [1e-6, 1e-9, 1e-9, 1e-6, 1e-9, 1e-6]
RTOL = 1e-10


def cases():
    """Build the systematic case set."""
    result = []
    for variant in VARIANTS:
        for component in COMPONENTS:
            for temperature in TEMPERATURES:
                result.append(([component], [1.0], temperature, variant))
    for i, first in enumerate(COMPONENTS):
        for second in COMPONENTS[i + 1:]:
            for composition in ([0.5, 0.5], [0.9, 0.1], [0.25, 0.75]):
                result.append(([first, second], list(composition), 298.15, "UNIFAC"))
            for composition in ([0.5, 0.5], [0.8, 0.2]):
                result.append(([first, second], list(composition), 298.15, "UMR"))
                result.append(([first, second], list(composition), 298.15, "VTPR"))
    ternary = [
        (["C3", "NC4", "C1"], [0.4, 0.35, 0.25]),
        (["H2O", "CO2", "C3"], [0.5, 0.3, 0.2]),
        (["C2", "BENZENE", "TOLUENE"], [0.5, 0.3, 0.2]),
        (["NC5", "IC5", "H2S"], [0.6, 0.3, 0.1]),
        (["N2", "CO2", "C1"], [0.2, 0.6, 0.2]),
    ]
    for components, composition in ternary:
        for variant in VARIANTS:
            for temperature in (250.0, 320.0):
                result.append((list(components), list(composition), temperature, variant))
    return result


def key(components, composition, temperature, variant):
    return "|".join(
        [
            "UNIFAC",
            ",".join(components),
            f"{temperature:.17e}",
            ",".join(f"{value:.17e}" for value in composition),
            variant,
        ]
    )


def native(lines):
    """Run the native parity executable for the given request lines."""
    if not lines:
        return []
    exe = ROOT / "target/release/examples" / ("parity.exe" if os.name == "nt" else "parity")
    output = subprocess.run(
        [str(exe)], input="\n".join(lines) + "\n", text=True, capture_output=True, check=True
    ).stdout.splitlines()
    assert len(output) == len(lines), (len(output), len(lines))
    return [line if line.startswith("ERROR") else list(map(float, line.split(","))) for line in output]


def build_fixture(name, rows, quantities, **extra):
    path = ROOT / "tests/fixtures" / name
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *quantities])
        writer.writerows(rows)
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        generator="tools/rust_port/validate_unifac.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=quantities,
        relative_tolerance=RTOL,
        absolute_tolerances=ATOL,
        **extra,
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def main(directory):
    request_list = cases()
    lines = []
    for components, composition, temperature, variant in request_list:
        lines.append(
            "|".join(
                [
                    "UNIFAC",
                    ",".join(components),
                    f"{temperature:.17e}",
                    "0",
                    ",".join(f"{value:.17e}" for value in composition),
                    variant,
                ]
            )
        )
    native_values = native(lines)

    reference_values = []
    reference_errors = []
    with ReferenceProcess(directory) as oracle:
        for index, (components, composition, temperature, variant) in enumerate(request_list):
            request = {
                "mode": "UNIFAC",
                "components": components,
                "composition": composition,
                "temperature": temperature,
                "variant": VARIANTS.index(variant),
            }
            result = oracle.call(request)
            if "values" not in result:
                reference_errors.append((index, result))
                reference_values.append(None)
            else:
                reference_values.append(result["values"])

    failures = []
    rows = []
    status_rows = []
    for index, (components, composition, temperature, variant) in enumerate(request_list):
        expected = reference_values[index]
        actual = native_values[index]
        case = key(components, composition, temperature, variant)
        if expected is None:
            # The original aborted; the native side must not report a value.
            if not isinstance(actual, str):
                failures.append(dict(case=case, reason="native value for aborted source", native=actual))
                status_rows.append([case, 1])
            else:
                status_rows.append([case, 1])
            continue
        if isinstance(actual, str):
            failures.append(dict(case=case, reason="native error", native=actual,
                                 reference=expected))
            status_rows.append([case, 2])
            continue
        if len(actual) != len(expected):
            failures.append(dict(case=case, reason="length", native=len(actual), reference=len(expected)))
            status_rows.append([case, 2])
            continue
        nc = len(components)
        atols = [ATOL[0], ATOL[1], ATOL[2]] + [ATOL[3]] * nc + [ATOL[4]] * nc + [ATOL[5]] * (nc * nc)
        assert len(atols) == len(expected)
        bad = []
        for position, (a, b, atol) in enumerate(zip(expected, actual, atols)):
            scaled = b * SOURCE_UNIT_SCALE
            if not math.isfinite(scaled) or abs(a - scaled) > atol + RTOL * max(abs(a), abs(scaled)):
                bad.append(position)
        if bad:
            failures.append(dict(case=case, quantities=bad, reference=expected, native=actual))
            status_rows.append([case, 3])
            continue
        status_rows.append([case, 0])
        rows.append([case, *expected])

    quantities = ["gibbs", "dg_dt", "d2g_dt2", "dg_dn", "d2g_dn_dt", "d2g_dn_dn"]
    build_fixture(
        "unifac.csv",
        rows,
        quantities,
        convention=(
            "original Ge_UNIFAC_GH_SG with the source component list, group database and "
            "combinatorial terms; mole fractions as the source uses them; the source reports "
            "J per kilomole, the native model reports J per mole and is scaled by 1000 "
            "for the comparison"
        ),
        source_unit_scale=SOURCE_UNIT_SCALE,
        source_status_counts={
            "compared": len(rows),
            "source_terminations": sum(1 for row in status_rows if row[1] == 1),
            "native_errors": sum(1 for row in status_rows if row[1] == 2),
            "mismatches": sum(1 for row in status_rows if row[1] == 3),
        },
    )
    build_fixture(
        "unifac_status.csv",
        status_rows,
        ["source_status"],
        convention=(
            "0 compared and matching, 1 the original terminated, 2 the native side reported "
            "an error, 3 a numerical mismatch"
        ),
    )
    (ROOT / ".reference/unifac_failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    print(
        json.dumps(
            {
                "cases": len(request_list),
                "compared": len(rows),
                "source_terminations": sum(1 for row in status_rows if row[1] == 1),
                "native_errors": sum(1 for row in status_rows if row[1] == 2),
                "mismatches": sum(1 for row in status_rows if row[1] == 3),
                "failures": len(failures),
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
