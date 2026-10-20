"""Compare the native pseudo component construction against the exact-source
oracle.

The original `init_cubic_pseudo` replaces the critical temperature, the
critical pressure, the acentric factor and the molar mass of the components
named `PSEUDO` and then selects the cubic equation of state. A reference
oracle builds the same model through the Python interface of the fork and
evaluates the properties of the mixture.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The pseudo component definitions. The lists are indexed by the position of
# the component, as in the source.
METHANE = (190.564, 4599200.0, 0.01142, 0.01604)
CARBON_DIOXIDE = (304.1282, 7336500.0, 0.22394, 0.04401)
DEFINITIONS = [
    # A single pseudo component.
    (["PSEUDO"], [350.0], [4.0e6], [0.3], [0.05]),
    # Two pseudo components with different critical points.
    (
        ["PSEUDO", "PSEUDO"],
        [350.0, 420.0],
        [4.0e6, 5.0e6],
        [0.3, 0.45],
        [0.05, 0.07],
    ),
    # A pseudo component with methane.
    (
        ["C1", "PSEUDO"],
        [METHANE[0], 350.0],
        [METHANE[1], 4.0e6],
        [METHANE[2], 0.3],
        [METHANE[3], 0.05],
    ),
    # A pseudo component with carbon dioxide.
    (
        ["PSEUDO", "CO2"],
        [350.0, CARBON_DIOXIDE[0]],
        [4.0e6, CARBON_DIOXIDE[1]],
        [0.3, CARBON_DIOXIDE[2]],
        [0.05, CARBON_DIOXIDE[3]],
    ),
]
# The states of every definition. The composition, the temperature, the
# pressure and the phase of the fluid.
STATES = [
    ([1.0], 250.0, 2.0e6, 1),
    ([1.0], 300.0, 1.0e6, 2),
    ([1.0], 280.0, 5.0e5, 2),
]
MIXTURE_STATES = [
    ([0.5, 0.5], 250.0, 2.0e6, 1),
    ([0.5, 0.5], 300.0, 5.0e6, 1),
    ([0.3, 0.7], 280.0, 1.0e7, 1),
    ([0.7, 0.3], 260.0, 3.0e6, 2),
]
EQUATIONS = ["PR", "SRK"]
QUANTITIES = ["volume", "enthalpy", "entropy", "compressibility"]
# The volume is a small difference of the state, the enthalpy and the entropy
# carry the reference levels of the source and the compressibility is a ratio
# of the two.
ATOL = [1e-12, 1e-6, 1e-9, 1e-12]
RTOL = 1e-10


def cases():
    for index, definition in enumerate(DEFINITIONS):
        components, temperatures, pressures, acentric, masses = definition
        for equation in EQUATIONS:
            if len(components) == 1:
                states = STATES
            else:
                states = MIXTURE_STATES
            for composition, temperature, pressure, phase in states:
                yield {
                    "index": index,
                    "equation": equation,
                    "components": components,
                    "critical_temperatures": temperatures,
                    "critical_pressures": pressures,
                    "acentric_factors": acentric,
                    "molar_masses": masses,
                    "composition": composition,
                    "temperature": temperature,
                    "pressure": pressure,
                    "phase": phase,
                }


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
    all_cases = list(cases())
    lines = []
    for case in all_cases:
        lines.append(
            "|".join(
                [
                    "PSEUDOC",
                    case["equation"],
                    ",".join(case["components"]),
                    ",".join(f"{value:.17e}" for value in case["critical_temperatures"]),
                    ",".join(f"{value:.17e}" for value in case["critical_pressures"]),
                    ",".join(f"{value:.17e}" for value in case["acentric_factors"]),
                    ",".join(f"{value:.17e}" for value in case["molar_masses"]),
                    f"{case['temperature']:.17e}",
                    f"{case['pressure']:.17e}",
                    ",".join(f"{value:.17e}" for value in case["composition"]),
                    str(case["phase"]),
                ]
            )
        )
    native_values = native(lines)

    with ReferenceProcess(directory) as oracle:
        reference_values = []
        for case in all_cases:
            reference_values.append(
                oracle.call(
                    {
                        "mode": "PSEUDO_COMPONENT",
                        "eos": case["equation"],
                        "components": case["components"],
                        "critical_temperatures": case["critical_temperatures"],
                        "critical_pressures": case["critical_pressures"],
                        "acentric_factors": case["acentric_factors"],
                        "molar_masses": case["molar_masses"],
                        "temperature": case["temperature"],
                        "pressure": case["pressure"],
                        "feed": case["composition"],
                        "phase": case["phase"],
                    }
                )
            )

    failures = []
    rows = []
    worst = 0.0
    for index, case in enumerate(all_cases):
        label = (
            f"{case['equation']}|{','.join(case['components'])}|"
            f"{case['temperature']:.17e}|{case['pressure']:.17e}|"
            f"{','.join(f'{value:.17e}' for value in case['composition'])}|{case['phase']}"
        )
        expected = reference_values[index]
        actual = native_values[index]
        if "values" not in expected:
            failures.append(dict(case=label, reason="source did not report values", detail=expected))
            continue
        if isinstance(actual, str):
            failures.append(dict(case=label, reason="native error", native=actual))
            continue
        bad = []
        for position, (source, native_value) in enumerate(zip(expected["values"], actual)):
            if not math.isfinite(native_value):
                bad.append(QUANTITIES[position])
                continue
            tolerance = ATOL[position] + RTOL * abs(source)
            error = abs(source - native_value) / max(abs(source), 1.0)
            if abs(source - native_value) > tolerance:
                bad.append(QUANTITIES[position])
            else:
                worst = max(worst, error)
        if bad:
            failures.append(
                dict(case=label, quantities=bad, reference=expected["values"], native=actual)
            )
            continue
        rows.append([label, *expected["values"]])

    path = ROOT / "tests/fixtures/pseudo_components.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(
            [
                "case",
                "equation",
                "components",
                "critical_temperatures",
                "critical_pressures",
                "acentric_factors",
                "molar_masses",
                "temperature",
                "pressure",
                "composition",
                "phase",
                *QUANTITIES,
            ]
        )
        for case, row in zip(all_cases, rows):
            writer.writerow(
                [
                    row[0],
                    case["equation"],
                    ",".join(case["components"]),
                    ",".join(f"{value:.17e}" for value in case["critical_temperatures"]),
                    ",".join(f"{value:.17e}" for value in case["critical_pressures"]),
                    ",".join(f"{value:.17e}" for value in case["acentric_factors"]),
                    ",".join(f"{value:.17e}" for value in case["molar_masses"]),
                    f"{case['temperature']:.17e}",
                    f"{case['pressure']:.17e}",
                    ",".join(f"{value:.17e}" for value in case["composition"]),
                    str(case["phase"]),
                    *[f"{value:.17e}" for value in row[1:]],
                ]
            )
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(all_cases),
        generator="tools/rust_port/validate_pseudo_components.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerances=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original init_cubic_pseudo with the components named PSEUDO; the source reports "
            "the molar volume in cubic metre per mole, the enthalpy and the entropy in joule per "
            "mole and joule per mole per kelvin and the compressibility factor, and the native "
            "model uses the same conventions"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/pseudo_components_failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    print(
        json.dumps(
            {
                "attempts": len(all_cases),
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
