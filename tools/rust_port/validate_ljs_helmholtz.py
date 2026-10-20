"""Compare the native splined Lennard-Jones Helmholtz energy against the
exact-source oracle.

The original `calcFresLJs_bh` evaluates the reduced residual Helmholtz energy
of the splined Lennard-Jones model. A reference probe calls the original
routine with the active model of the fork.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess

# The cases: the temperature in kelvin, the volume in cubic metre and the
# number of moles.
CASES = [
    (100.0, 1.0e-3, 1.0),
    (150.0, 1.0e-3, 1.0),
    (200.0, 1.0e-3, 1.0),
    (300.0, 1.0e-3, 1.0),
    (150.0, 5.0e-4, 1.0),
    (200.0, 5.0e-4, 1.0),
    (250.0, 2.0e-4, 1.0),
    (120.0, 1.0e-3, 0.5),
]
QUANTITIES = ["reduced_helmholtz", "reduced_helmholtz_t"]
RTOL = 1.0e-9
ATOL = 1.0e-12


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
    lines = [
        f"LJSHELMHOLTZ|{temperature:.17e}|{volume:.17e}|{moles:.17e}"
        for temperature, volume, moles in CASES
    ]
    native_values = native(lines)

    with ReferenceProcess(directory) as oracle:
        reference_values = []
        for temperature, volume, moles in CASES:
            reference_values.append(
                oracle.call(
                    {
                        "mode": "LJS_HELMHOLTZ",
                        "temperature": temperature,
                        "volume": volume,
                        "moles": moles,
                    }
                )
            )

    failures = []
    rows = []
    worst = 0.0
    for index, case in enumerate(CASES):
        temperature, volume, moles = case
        label = f"{temperature:.17e}|{volume:.17e}|{moles:.17e}"
        expected = reference_values[index]
        actual = native_values[index]
        if "values" not in expected:
            failures.append(dict(case=label, reason="source did not report values", detail=expected))
            continue
        if isinstance(actual, str):
            failures.append(dict(case=label, reason="native error", native=actual))
            continue
        source = expected["values"][0]
        source_t = expected["values"][1]
        native_value = actual[0]
        native_t = actual[1]
        if abs(source - native_value) > ATOL + RTOL * abs(source):
            failures.append(
                dict(case=label, reason="reduced Helmholtz energy", source=source, native=native_value)
            )
            continue
        if abs(source_t - native_t) > ATOL + RTOL * abs(source_t):
            failures.append(
                dict(case=label, reason="reduced Helmholtz energy derivative", source=source_t, native=native_t)
            )
            continue
        worst = max(worst, abs(source - native_value) / max(abs(source), 1.0e-30))
        worst = max(worst, abs(source_t - native_t) / max(abs(source_t), 1.0e-30))
        rows.append([label, source, source_t])

    path = ROOT / "tests/fixtures/ljs_helmholtz.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *QUANTITIES])
        for row in rows:
            writer.writerow([row[0], f"{row[1]:.17e}", f"{row[2]:.17e}"])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(rows),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_ljs_helmholtz.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,
        relative_tolerance=RTOL,
        absolute_tolerance=ATOL,
        worst_relative_difference=worst,
        convention=(
            "the original calcFresLJs_bh; the reduced residual Helmholtz energy of the whole system, "
            "the hard-sphere reference and the three perturbation terms of the original"
        ),
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/ljs_helmholtz_failures.json").write_text(
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
