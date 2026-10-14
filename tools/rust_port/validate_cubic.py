"""Generate exact-fork fixtures and compare the Rust parity example.

Run after build_reference.py. Committed fixtures run with cargo test alone.
"""
import argparse
import csv
import hashlib
import json
import os
import subprocess
from reference import ROOT, load


def cases():
    mixtures = [("C1", [1.0]), ("CO2", [1.0]), ("N2", [1.0]),
                ("C3", [1.0]), ("NC4", [1.0]), ("H2O", [1.0]),
                ("C1,CO2", [0.7, 0.3]), ("C1,C2,C3", [0.5, 0.3, 0.2]),
                ("CO2,N2", [0.8, 0.2]), ("CO2,H2O", [0.9, 0.1]),
                ("C1,CO2,N2", [0.0, 0.8, 0.2])]
    for eos in ["PR", "SRK", "RK", "VDW"]:
        for comps, x in mixtures:
            for t in [200.0, 300.0, 500.0]:
                for p in [1e3, 1e5, 1e6, 5e6, 2e7]:
                    for phase in ["L", "V"]:
                        yield eos, comps, t, p, x, phase


def validate(compiler_directory):
    cubic = load(compiler_directory)
    inputs, expected = [], []
    model, previous = None, None
    for eos, comps, t, p, x, phase in cases():
        if previous != (eos, comps):
            model = cubic(comps, "SRK" if eos == "RK" else eos,
                          alpha="RK" if eos == "RK" else "Classic", volume_shift=False)
            previous = eos, comps
        flag = model.LIQPH if phase == "L" else model.VAPPH
        v = model.specific_volume(t, p, x, flag)[0]
        pressure, pt, pv = model.pressure_tv(t, v, x, dpdt=True, dpdv=True)
        hr = model.enthalpy(t, p, x, flag, residual=True)[0]
        sr = model.entropy(t, p, x, flag, residual=True)[0]
        phi = model.thermo(t, p, x, flag)[0]
        h, cp = model.enthalpy(t, p, x, flag, dhdt=True)
        s = model.entropy(t, p, x, flag)[0]
        inputs.append("|".join([eos, comps, str(t), str(p), ",".join(map(str, x)), phase]))
        expected.append([v, pressure, pt, pv, hr, sr, h, s, cp, *map(float, phi)])
    subprocess.run(["cargo", "build", "--release", "--example", "parity"], cwd=ROOT, check=True)
    executable = ROOT / "target/release/examples" / ("parity.exe" if os.name == "nt" else "parity")
    result = subprocess.run([str(executable)], input="\n".join(inputs) + "\n", text=True,
                            capture_output=True, check=True)
    actual = [[float(v) for v in row.split(',')] for row in result.stdout.splitlines()]
    assert len(actual) == len(expected)
    failures, worst = [], {}
    labels = ["v", "p", "dpdt", "dpdv", "h_res", "s_res", "h", "s", "dhdt_baseline"]
    absolute = [1e-14, 1e-2, 1e-6, 1e-2, 1e-7, 1e-9, 1e-7, 1e-9, 1e-7]
    for case, ref, native in zip(inputs, expected, actual):
        assert len(ref) == len(native)
        for i, (a, b) in enumerate(zip(ref, native)):
            label = labels[i] if i < len(labels) else "ln_phi"
            error = abs(a-b) / max(abs(a), abs(b), 1.0)
            worst[label] = max(error, worst.get(label, 0.0))
            atol = absolute[i] if i < len(absolute) else 1e-9
            rtol = 2e-9 if i == 0 else 2e-8
            tolerance = atol + rtol * max(abs(a), abs(b))
            if not abs(a-b) <= tolerance:
                failures.append({"case": case, "quantity": label, "reference": a,
                                 "rust": b, "scaled_error": error})
    output = ROOT / "tests/fixtures"
    output.mkdir(exist_ok=True)
    with (output / "cubic.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=';')
        writer.writerow(["case", *labels, "ln_phi..."])
        writer.writerows([[case, *values] for case, values in zip(inputs, expected)])
    metadata = json.loads((ROOT / ".reference/metadata.json").read_text())
    metadata.update({"case_count": len(inputs), "worst_scaled_error": worst,
                     "tolerances": {label: {"absolute": absolute[i], "relative": 2e-9 if i == 0 else 2e-8}
                                    for i, label in enumerate(labels)} | {"ln_phi": {"absolute": 1e-9, "relative": 2e-8}},
                     "fixture_sha256": hashlib.sha256((output / "cubic.csv").read_bytes()).hexdigest(),
                     "generator": "tools/rust_port/validate_cubic.py"})
    (output / "cubic.metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / ".reference/cubic_failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    print(json.dumps({"cases": len(inputs), "failures": len(failures), "worst": worst}, indent=2))
    if failures:
        print(json.dumps(failures[:12], indent=2))
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler-directory")
    validate(parser.parse_args().compiler_directory)
