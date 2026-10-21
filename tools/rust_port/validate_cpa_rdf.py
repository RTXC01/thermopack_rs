"""Compare the native hard-sphere contact radial distribution function of the
cubic association model against the exact-source oracle.

The original `saft_rdf::g_rdf_cpa` evaluates the contact function of the
cubic association model. A reference probe calls the original
`master_saft_rdf` for the active model.
"""
import argparse
import csv
import hashlib
import json
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the component, the temperature, the volume and the number of
# moles of every component.
CASES = [
    ("H2O", 300.0, 1.0e-3, [1.0]),
    ("H2O", 350.0, 2.0e-3, [2.0]),
    ("MEOH", 280.0, 5.0e-4, [0.5]),
    ("H2O", 298.15, 1.0e-2, [3.0]),
    ("ETOH", 320.0, 1.5e-3, [0.75]),
]
QUANTITIES = [
    "g",
    "g_t",
    "g_v",
    "g_n",
    "g_tt",
    "g_tv",
    "g_tn",
    "g_vv",
    "g_vn",
    "g_nn",
]
RTOL = 1.0e-12
ATOL = 1.0e-14


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
        for component, temperature, volume, moles in CASES:
            result = oracle.call(
                {
                    "mode": "CPA_RDF",
                    "component": component,
                    "temperature": temperature,
                    "volume": volume,
                    "moles": moles,
                }
            )
            if "values" not in result:
                raise SystemExit(f"the source did not report the function: {result}")
            reference_values.append(result)

    lines = [
        "CPARDF|{0}|{1:.17e}|{2:.17e}|{3}".format(
            case[0], case[1], case[2], ",".join(f"{value:.17e}" for value in case[3])
        )
        for case in CASES
    ]
    native_values = native(lines)

    failures = []
    rows = []
    worst = 0.0
    for index, case in enumerate(CASES):
        label = "{0}|{1:.17e}|{2:.17e}|{3}".format(
            case[0], case[1], case[2], ",".join(f"{value:.17e}" for value in case[3])
        )
        source = reference_values[index]["values"]
        actual = native_values[index]
        if isinstance(actual, str):
            failures.append(dict(case=label, reason="native error", native=actual))
            continue
        bad = []
        if len(source) != len(actual):
            failures.append(
                dict(
                    case=label,
                    reason="number of values",
                    source=len(source),
                    native=len(actual),
                )
            )
            continue
        for position, (source_value, native_value) in enumerate(zip(source, actual)):
            if abs(source_value - native_value) > ATOL + RTOL * abs(source_value):
                bad.append((QUANTITIES[position], source_value, native_value))
            else:
                worst = max(
                    worst, abs(source_value - native_value) / max(abs(source_value), 1.0e-30)
                )
        if bad:
            failures.append(dict(case=label, reason="values", differences=bad[:6]))
            continue
        rows.append([label, *source])

    path = ROOT / "tests/fixtures/cpa_rdf.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *QUANTITIES])
        for row in rows:
            writer.writerow([row[0], *[f"{value:.17e}" for value in row[1:]]])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_cpa_rdf.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original saft_rdf::g_rdf_cpa, which evaluates the hard-sphere contact "
            "radial distribution function of the cubic association model with its volume "
            "and mole number derivatives"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/cpa_rdf_failures.json").write_text(json.dumps(failures, indent=2) + "\n")
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
