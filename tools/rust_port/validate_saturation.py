"""Pure saturation differential fixtures from the exact-source reference."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
from reference import ROOT, load


def validate(compiler_directory):
    cubic = load(compiler_directory)
    inputs, expected = [], []
    for eos in ["PR", "SRK", "RK", "VDW"]:
        for name in ["C1", "CO2", "C3"]:
            model = cubic(name, "SRK" if eos == "RK" else eos,
                          alpha="RK" if eos == "RK" else "Classic", volume_shift=False)
            tc = model.get_critical_parameters(1)[0]
            for tr in [0.65, 0.8, 0.95, 0.999]:
                t = tc * tr
                p = model.bubble_pressure(t, [1.0])[0]
                vl = model.specific_volume(t, p, [1.0], model.LIQPH)[0]
                vv = model.specific_volume(t, p, [1.0], model.VAPPH)[0]
                hl = model.enthalpy(t, p, [1.0], model.LIQPH, residual=True)[0]
                hv = model.enthalpy(t, p, [1.0], model.VAPPH, residual=True)[0]
                inputs.append(f"{eos}|{name}|{t}|SAT|1|L")
                expected.append([p, vl, vv, (hv-hl)/(t*(vv-vl))])
    subprocess.run(["cargo", "build", "--release", "--example", "parity"], cwd=ROOT, check=True)
    executable = ROOT / "target/release/examples" / ("parity.exe" if os.name == "nt" else "parity")
    result = subprocess.run([str(executable)], input="\n".join(inputs) + "\n", text=True,
                            capture_output=True, check=True)
    actual = [list(map(float, row.split(','))) for row in result.stdout.splitlines()]
    assert len(actual) == len(expected)
    absolute = [1e-3, 1e-13, 1e-13, 1e-3]
    worst = [0.0] * 4
    for case, ref, native in zip(inputs, expected, actual):
        assert len(native) == len(ref)
        for i, (a, b) in enumerate(zip(ref, native)):
            worst[i] = max(worst[i], abs(a-b)/max(abs(a),abs(b)))
            assert abs(a-b) <= absolute[i] + 2e-7 * max(abs(a), abs(b)), (case, i, a, b)
    output = ROOT / "tests/fixtures/saturation.csv"
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=';', lineterminator='\n')
        writer.writerow(["case", "pressure", "liquid_volume", "vapor_volume", "dpdt"])
        writer.writerows([[case, *row] for case, row in zip(inputs, expected)])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update({"case_count": len(inputs), "worst_relative_error": worst,
                     "tolerance": {"absolute": absolute, "relative": 2e-7},
                     "fixture_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
                     "generator": "tools/rust_port/validate_saturation.py"})
    output.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"cases": len(inputs), "worst_relative_error": worst}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler-directory")
    validate(parser.parse_args().compiler_directory)
