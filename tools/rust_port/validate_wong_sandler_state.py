"""Compare the native state of a mixture with the Wong-Sandler mixing rule
against the exact-source oracle.

The original `wong_sandler::WongSandlerMix` replaces the binary attraction and
the binary covolume of the cubic equation of state with its own coefficients.
A reference probe evaluates the state of the mixture with the rule at a molar
volume, the native model evaluates the same state.
"""

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reference import ReferenceProcess, ROOT  # noqa: E402

# eos, components, temperature [K], molar volume [m3/mol], composition,
# kij, tauij, alpha
CASES = [
    ("PR", "CO2,NC4", 350.0, 2.0e-4, [0.3, 0.7], 0.05, 0.001, 0.3),
    ("PR", "CO2,NC4", 300.0, 5.0e-3, [0.3, 0.7], 0.0, 0.0, 0.0),
    ("SRK", "CO2,NC4", 350.0, 2.0e-3, [0.8, 0.2], 0.12, 0.002, 0.5),
    ("SRK", "CO2,C3", 320.0, 1.0e-2, [0.6, 0.4], 0.08, 0.0005, 0.25),
    ("PR", "CO2,C3", 300.0, 1.5e-4, [0.5, 0.5], 0.0, 0.0, 0.0),
    ("VDW", "CO2,NC4", 350.0, 1.0e-3, [0.4, 0.6], 0.03, 0.0, 0.0),
]

ABSOLUTE_TOLERANCE = 1.0e-8
RELATIVE_TOLERANCE = 1.0e-9


def native_lines():
    for eos, components, temperature, volume, composition, kij, tauij, alpha in CASES:
        yield (
            "WSSTATE|{0}|{1}|{2:.17e}|{3:.17e}|{4}|{5:.17e}|{6:.17e}|{7:.17e}".format(
                eos,
                components,
                temperature,
                volume,
                ",".join(f"{value:.17e}" for value in composition),
                kij,
                tauij,
                alpha,
            )
        )


def source_values(oracle):
    values = []
    for eos, components, temperature, volume, composition, kij, tauij, alpha in CASES:
        result = oracle.call(
            {
                "mode": "WONGSANDLERSTATE",
                "components": components,
                "eos": eos,
                "temperature": temperature,
                "volume": volume,
                "composition": composition,
                "kij": kij,
                "tauij": tauij,
                "alpha": alpha,
            }
        )
        values.append(list(result["values"]))
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler-directory", help="directory of the runtime DLLs")
    parser.add_argument("--report", help="write the report to this file")
    arguments = parser.parse_args()

    with ReferenceProcess(arguments.compiler_directory) as oracle:
        sources = source_values(oracle)

    subprocess.run(
        ["cargo", "build", "--release", "--example", "parity"],
        cwd=ROOT,
        check=True,
    )
    executable = ROOT / "target/release/examples" / (
        "parity.exe" if os.name == "nt" else "parity"
    )
    process = subprocess.run(
        [str(executable)],
        input="\n".join(native_lines()) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    lines = process.stdout.splitlines()
    assert len(lines) == len(CASES)

    failures = []
    worst = 0.0
    worst_case = None
    compared = 0
    for case, source, line in zip(CASES, sources, lines):
        native = [float(value) for value in line.split(",")]
        assert len(native) == len(source), (len(native), len(source))
        for index, (expected, actual) in enumerate(zip(source, native)):
            difference = abs(expected - actual)
            limit = ABSOLUTE_TOLERANCE + RELATIVE_TOLERANCE * abs(expected)
            relative = difference / abs(expected) if expected != 0.0 else difference
            if relative > worst:
                worst = relative
                worst_case = (case, index, expected, actual)
            compared += 1
            if difference > limit:
                failures.append(
                    {
                        "case": case[0] + "|" + case[1],
                        "index": index,
                        "source": expected,
                        "native": actual,
                        "relative": relative,
                    }
                )

    report = {
        "attempts": len(CASES),
        "compared": compared,
        "failures": len(failures),
        "worst_relative_difference": worst,
        "worst_case": worst_case,
        "absolute_tolerance": ABSOLUTE_TOLERANCE,
        "relative_tolerance": RELATIVE_TOLERANCE,
        "failures_detail": failures,
        "convention": (
            "the state of the mixture with the Wong-Sandler mixing rule at a molar "
            "volume, as the original wong_sandler::WongSandlerMix replaces the binary "
            "attraction and the binary covolume of the cubic equation of state with its "
            "own coefficients"
        ),
    }

    path = ROOT / "tests/fixtures/wong_sandler_state.csv"
    quantities = ["p", "dpdv", "dpdt", "d2pdv2", "h", "s"]
    for index in range(2):
        quantities.append(f"lnphi_{index + 1}")
    for index in range(2):
        quantities.append(f"dlnfugdt_{index + 1}")
    for index in range(2):
        quantities.append(f"dlnfugdp_{index + 1}")
    for first in range(2):
        for second in range(2):
            quantities.append(f"dlnfugdn_{first + 1}_{second + 1}")
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator="\n")
        writer.writerow(["case", *quantities])
        for case, source in zip(CASES, sources):
            label = "WSSTATE|{0}|{1}|{2:.17e}|{3:.17e}|{4}|{5:.17e}|{6:.17e}|{7:.17e}".format(
                case[0],
                case[1],
                case[2],
                case[3],
                ",".join(f"{value:.17e}" for value in case[4]),
                case[5],
                case[6],
                case[7],
            )
            writer.writerow([label, *[f"{value:.17e}" for value in source]])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update(
        cases=len(CASES),
        attempted_cases=len(CASES),
        generator="tools/rust_port/validate_wong_sandler_state.py",
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=quantities,
        relative_tolerance=RELATIVE_TOLERANCE,
        absolute_tolerance=ABSOLUTE_TOLERANCE,
        worst_relative_difference=worst,
    )
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    if arguments.report:
        Path(arguments.report).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
