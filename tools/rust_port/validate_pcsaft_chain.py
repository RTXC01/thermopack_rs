"""Compare the native reduced Helmholtz energy of the chain term of the
PC-SAFT model against the exact-source oracle.

The original `pc_saft_nonassoc::F_chain_PC_SAFT_TVn` evaluates the reduced
Helmholtz energy of the chain term of the PC-SAFT model and its derivatives
with respect to the temperature, the volume and the mole numbers. A reference
probe calls the published routine of the original.
"""
import argparse
import csv
import hashlib
import json
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the components, the temperature, the volume and the number of
# moles of every component.
CASES = [
    (["CO2"], 300.0, 1.0e-3, [0.1]),
    (["CO2", "C1"], 300.0, 1.0e-3, [0.05, 0.05]),
    (["CO2", "C1"], 320.0, 5.0e-4, [0.1, 0.2]),
    (["N2", "CO2"], 280.0, 2.0e-3, [0.1, 0.1]),
    (["C1", "C2"], 350.0, 1.0e-3, [0.05, 0.05]),
    (["H2O", "MEOH"], 300.0, 1.0e-3, [0.1, 0.1]),
]
QUANTITIES = [
    "F",
    "F_t",
    "F_v",
    "F_n",
    "F_tt",
    "F_tv",
    "F_tn",
    "F_vv",
    "F_vn",
    "F_nn",
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


def label(case):
    return "{0}|{1:.17e}|{2:.17e}|{3}".format(
        ",".join(case[0]),
        case[1],
        case[2],
        ",".join(f"{value:.17e}" for value in case[3]),
    )


def main(directory):
    with ReferenceProcess(directory) as oracle:
        reference_values = []
        for components, temperature, volume, moles in CASES:
            result = oracle.call(
                {
                    "mode": "PCSAFT_CHAIN",
                    "components": ",".join(components),
                    "temperature": temperature,
                    "volume": volume,
                    "moles": moles,
                }
            )
            if "values" not in result:
                raise SystemExit(f"the source did not report the energy: {result}")
            reference_values.append(result)

    lines = [
        "PCSAFTCHAIN|{0}|{1:.17e}|{2:.17e}|{3}".format(
            ",".join(case[0]),
            case[1],
            case[2],
            ",".join(f"{value:.17e}" for value in case[3]),
        )
        for case in CASES
    ]
    native_values = native(lines)

    failures = []
    rows = []
    worst = 0.0
    for index, case in enumerate(CASES):
        name = label(case)
        source = reference_values[index]["values"]
        actual = native_values[index]
        if isinstance(actual, str):
            failures.append(dict(case=name, reason="native error", native=actual))
            continue
        bad = []
        if len(source) != len(actual):
            failures.append(
                dict(
                    case=name,
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
            failures.append(dict(case=name, reason="values", differences=bad[:6]))
            continue
        rows.append([name, *source])

    path = ROOT / "tests/fixtures/pcsaft_chain.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *QUANTITIES])
        for row in rows:
            writer.writerow([row[0], *[f"{value:.17e}" for value in row[1:]]])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_pcsaft_chain.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original pc_saft_nonassoc::F_chain_PC_SAFT_TVn, which evaluates the reduced "
            "Helmholtz energy of the chain term of the PC-SAFT model with its temperature, "
            "volume and mole number derivatives"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/pcsaft_chain_failures.json").write_text(
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
