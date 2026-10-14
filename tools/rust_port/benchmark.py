"""Compare native Fortran and Rust loops, including checksums and medians."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
from reference import ROOT


def benchmark(compiler, repeats):
    compiler = Path(compiler).resolve()
    env = dict(os.environ)
    env["PATH"] = str(compiler.parent) + os.pathsep + env["PATH"]
    out = ROOT / ".reference"
    executable = out / ("benchmark.exe" if os.name == "nt" else "benchmark")
    flags = ["-O3", "-fdefault-real-8", "-fdefault-double-8", "-frecursive", "-fopenmp"]
    subprocess.run([str(compiler), *flags, "-I", str(out),
                    str(ROOT / "tools/rust_port/benchmark_reference.f90"),
                    "@" + str(out / "objects.rsp"), "-o", str(executable)],
                   env=env, cwd=out, check=True)
    subprocess.run(["cargo", "bench", "--bench", "cubic", "--no-run"], cwd=ROOT, check=True)
    def parse(output):
        return {row.split(',')[0]: list(map(float, row.split(',')[1:]))
                for row in output.splitlines() if ',' in row}
    timings = {}
    for i in range(repeats + 1):
        reference = parse(subprocess.check_output([str(executable)], env=env, text=True))
        native = parse(subprocess.check_output(["cargo", "bench", "--bench", "cubic"], cwd=ROOT, text=True))
        for name, measurements in native.items():
            ref = reference[name.split('_')[0]]
            for a, b in zip(measurements[1:], ref[1:]):
                assert abs(a-b) <= 2e-7 * max(abs(a), abs(b), 1.0), (name, a, b)
            if i > 0:
                timings.setdefault("rust_" + name, []).append(measurements[0])
        if i > 0:
            for name, measurements in reference.items():
                timings.setdefault("fortran_" + name, []).append(measurements[0])
        print(f"Benchmark repetition {i}/{repeats} (0 is warm-up)", flush=True)
    report = {"platform": platform.platform(), "processor": platform.processor(),
              "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
              "reference": json.loads((out / "metadata.json").read_text()),
              "native_sources_sha256": {path.relative_to(ROOT).as_posix():
                  hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
                  for path in sorted((ROOT / "rust").rglob("*.rs"))},
              "cargo_profile": {"optimization": "release", "lto": "thin", "codegen_units": 1},
              "iterations_per_sample": 1000000, "samples_ns_per_state": timings,
              "median_ns_per_state": {name: statistics.median(values) for name, values in timings.items()},
              "workload": "PR CO2/C1/N2 0.7/0.2/0.1, 300 K; pressure+T/V derivatives, vapor volume"}
    (ROOT / "docs/rust/benchmark.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["median_ns_per_state"], indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler", required=True)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    benchmark(args.compiler, args.repeats)
