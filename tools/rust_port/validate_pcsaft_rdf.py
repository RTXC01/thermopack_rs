"""Compare the native hard-sphere contact radial distribution function of the
PC-SAFT model against the exact-source oracle.

The original `pc_saft_nonassoc::g_pc_saft_TVn` applies the chain rule of the
hard-sphere model of Boublik to the packing fractions of the mixture. A
reference probe calls the published routines `calc_dhs`,
`calc_bmcsl_zeta_and_derivatives` and `calc_bmcsl_gij` of the original.
"""
import argparse
import csv
import hashlib
import json
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the components, the temperature, the volume, the number of moles
# of every component and the components of the pair. The indices of the pair
# are zero based, as in the native model, and the probe of the source uses the
# one based indices of Fortran.
CASES = [
    (["CO2"], 300.0, 1.0e-3, [0.1], 0, 0),
    (["CO2", "C1"], 300.0, 1.0e-3, [0.05, 0.05], 0, 0),
    (["CO2", "C1"], 320.0, 5.0e-4, [0.1, 0.2], 0, 1),
    (["N2", "CO2"], 280.0, 2.0e-3, [0.1, 0.1], 1, 1),
    (["C1", "C2"], 350.0, 1.0e-3, [0.05, 0.05], 1, 0),
    (["H2O", "MEOH"], 300.0, 1.0e-3, [0.1, 0.1], 0, 0),
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


def label(case):
    return "{0}|{1:.17e}|{2:.17e}|{3}|{4}|{5}".format(
        ",".join(case[0]),
        case[1],
        case[2],
        ",".join(f"{value:.17e}" for value in case[3]),
        case[4],
        case[5],
    )


def main(directory):
    with ReferenceProcess(directory) as oracle:
        reference_values = []
        for components, temperature, volume, moles, first, second in CASES:
            result = oracle.call(
                {
                    "mode": "PCSAFT_RDF",
                    "components": ",".join(components),
                    "temperature": temperature,
                    "volume": volume,
                    "moles": moles,
                    "first": first + 1,
                    "second": second + 1,
                }
            )
            if "values" not in result:
                raise SystemExit(f"the source did not report the function: {result}")
            reference_values.append(result)

    lines = [
        "PCSAFTRDF|{0}|{1:.17e}|{2:.17e}|{3}|{4}|{5}".format(
            ",".join(case[0]),
            case[1],
            case[2],
            ",".join(f"{value:.17e}" for value in case[3]),
            case[4],
            case[5],
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

    path = ROOT / "tests/fixtures/pcsaft_rdf.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *QUANTITIES])
        for row in rows:
            writer.writerow([row[0], *[f"{value:.17e}" for value in row[1:]]])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_pcsaft_rdf.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original pc_saft_nonassoc::g_pc_saft_TVn, which applies the chain rule "
            "of the hard-sphere model of Boublik to the packing fractions of the mixture "
            "and delivers the pair correlation function with its temperature, volume and "
            "mole number derivatives"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/pcsaft_rdf_failures.json").write_text(json.dumps(failures, indent=2) + "\n")
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
