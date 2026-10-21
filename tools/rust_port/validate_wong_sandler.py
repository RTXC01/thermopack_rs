"""Compare the native mixture coefficients of the Wong-Sandler mixing rule
against the exact-source oracle.

The original `WongSandlerMix` replaces the binary attraction and the binary
covolume of the cubic equation of state with its own coefficients. A reference
probe calls the original routine for the active model.
"""
import argparse
import csv
import hashlib
import json
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the sub-equation of state, the components, the temperature, the
# composition and the binary parameters of the rule.
CASES = [
    ("PR", "CO2,NC4", 300.0, [0.5, 0.5], 0.1, 0.0, 0.0),
    ("PR", "CO2,NC4", 350.0, [0.3, 0.7], 0.05, 0.001, 0.3),
    ("SRK", "CO2,NC4", 280.0, [0.8, 0.2], 0.12, 0.002, 0.5),
    ("PR", "CO2,NC4", 400.0, [0.2, 0.8], 0.0, 0.0, 0.0),
    ("SRK", "CO2,C3", 320.0, [0.6, 0.4], 0.08, 0.0005, 0.25),
    ("VDW", "CO2,NC4", 300.0, [0.4, 0.6], 0.03, 0.0, 0.0),
    ("PR", "CO2,NC4", 320.0, [0.6, 0.4], 0.0, 0.002, 0.0),
    ("SRK", "CO2,C3", 300.0, [0.5, 0.5], 0.1, 0.0, 0.4),
]
QUANTITIES = [
    "a",
    "a_t",
    "a_tt",
    "a_i_1",
    "a_i_2",
    "a_it_1",
    "a_it_2",
    "b",
    "b_t",
    "b_tt",
    "b_i_1",
    "b_i_2",
    "b_it_1",
    "b_it_2",
]
RTOL = 1.0e-10
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
        for eos, components, temperature, composition, kij, tauij, alpha in CASES:
            result = oracle.call(
                {
                    "mode": "WONGSANDLER",
                    "components": components,
                    "eos": eos,
                    "temperature": temperature,
                    "composition": composition,
                    "kij": kij,
                    "tauij": tauij,
                    "alpha": alpha,
                }
            )
            if "values" not in result:
                raise SystemExit(f"the source did not report the coefficients: {result}")
            reference_values.append(result)

    lines = []
    for eos, components, temperature, composition, kij, tauij, alpha in CASES:
        lines.append(
            "WSMIX|{0}|{1}|{2:.17e}|{3}|{4:.17e}|{5:.17e}|{6:.17e}".format(
                eos,
                components,
                temperature,
                ",".join(f"{value:.17e}" for value in composition),
                kij,
                tauij,
                alpha,
            )
        )
    native_values = native(lines)

    failures = []
    rows = []
    worst = 0.0
    for index, case in enumerate(CASES):
        label = "{0}|{1}|{2:.17e}|{3}|{4:.17e}|{5:.17e}|{6:.17e}".format(
            case[0],
            case[1],
            case[2],
            ",".join(f"{value:.17e}" for value in case[3]),
            case[4],
            case[5],
            case[6],
        )
        source = list(reference_values[index]["values"])
        # The cubic equation of state of the original keeps the attraction in
        # pascal litre squared per mole squared and the covolume in litres per
        # mole, the native library uses the SI units.
        source[:7] = [value * 1.0e-6 for value in source[:7]]
        source[7:] = [value * 1.0e-3 for value in source[7:]]
        actual = native_values[index]
        if isinstance(actual, str):
            failures.append(dict(case=label, reason="native error", native=actual))
            continue
        if len(source) != len(actual):
            failures.append(
                dict(case=label, reason="number of values", source=len(source), native=len(actual))
            )
            continue
        bad = []
        for position, (source_value, native_value) in enumerate(zip(source, actual)):
            if abs(source_value - native_value) > ATOL + RTOL * abs(source_value):
                bad.append((QUANTITIES[position], source_value, native_value))
            else:
                worst = max(
                    worst, abs(source_value - native_value) / max(abs(source_value), 1.0e-30)
                )
        if bad:
            failures.append(dict(case=label, reason="values", differences=bad))
            continue
        rows.append([label, *source])

    path = ROOT / "tests/fixtures/wong_sandler.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *QUANTITIES])
        for row in rows:
            writer.writerow([row[0], *[f"{value:.17e}" for value in row[1:]]])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_wong_sandler.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original wong_sandler::WongSandlerMix, which replaces the binary attraction "
            "and the binary covolume of the cubic equation of state with its own coefficients"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/wong_sandler_failures.json").write_text(
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
