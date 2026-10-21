"""Compare the native chain rule of the packing fraction against the
exact-source oracle.

The original `convert_zeta_x_to_TVn` turns the derivatives of a term with
respect to the packing fraction and the reduced distance into the derivatives
with respect to the temperature, the volume and the number of moles. A
reference probe calls the original routine for a fixed test case.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

QUANTITIES = [
    "temperature",
    "volume",
    "moles_1",
    "moles_2",
    "temperature_second",
    "volume_second",
    "temperature_volume",
    "temperature_moles_1",
    "temperature_moles_2",
    "volume_moles_1",
    "volume_moles_2",
    "moles_second_11",
    "moles_second_12",
    "moles_second_21",
    "moles_second_22",
]
RTOL = 1.0e-13
ATOL = 1.0e-15


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
        result = oracle.call({"mode": "CHAIN_RULE"})
    if "values" not in result:
        raise SystemExit(f"the source did not report the derivatives: {result}")
    source = result["values"]
    actual = native(["CHAINRULE"])[0]
    failures = []
    worst = 0.0
    for position, (source_value, native_value) in enumerate(zip(source, actual)):
        if abs(source_value - native_value) > ATOL + RTOL * abs(source_value):
            failures.append(
                dict(
                    case=QUANTITIES[position],
                    source=source_value,
                    native=native_value,
                )
            )
        else:
            worst = max(worst, abs(source_value - native_value) / max(abs(source_value), 1.0e-30))
    rows = [["chain_rule", *[f"{value:.17e}" for value in source]]]

    path = ROOT / "tests/fixtures/chain_rule.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *QUANTITIES])
        for row in rows:
            writer.writerow([row[0], *row[1:]])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=1,
        generator="tools/rust_port/validate_chain_rule.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original convert_zeta_x_to_TVn for a fixed test case with two components; the "
            "native model follows the chain rule of the source for the first and the second "
            "derivatives"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/chain_rule_failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    print(
        json.dumps(
            {
                "attempts": 1,
                "compared": len(rows) if not failures else 0,
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
