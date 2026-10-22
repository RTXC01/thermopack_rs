"""Native-executable PS and complete metastable-isentrope comparisons."""
import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess

from reference import ROOT
from benchmark_isochoric import fingerprints
from benchmark_solid_water import parse

ATOL = [2e-5, 1e-10, 2., 2e-5]
RTOL = 2e-6


def benchmark(compiler, repeats, output, preflight_only):
    if repeats < 1:
        raise ValueError('Positive sample count required')
    compiler = Path(compiler).resolve()
    env = dict(os.environ)
    env['PATH'] = str(compiler.parent)+os.pathsep+env['PATH']
    directory = ROOT/'.reference'
    source = ROOT/'tools/rust_port/benchmark_metastable_reference.f90'
    exe = directory/('benchmark_metastable.exe' if os.name == 'nt' else 'benchmark_metastable')
    flags = ['-O3', '-fdefault-real-8', '-fdefault-double-8', '-frecursive', '-fopenmp']
    native_paths = sorted((ROOT/'rust').rglob('*.rs'))
    harness_paths = [Path(__file__).resolve(), source, ROOT/'tools/rust_port/benchmark_isochoric.py', ROOT/'tools/rust_port/benchmark_solid_water.py']
    cargo_paths = [ROOT/'Cargo.toml', ROOT/'Cargo.lock']
    rust, harness, cargo = map(fingerprints, [native_paths, harness_paths, cargo_paths])
    subprocess.run([str(compiler), *flags, '-I', str(directory), str(source), '@'+str(directory/'objects.rsp'), '-o', str(exe)], cwd=directory, env=env, check=True)
    command = ['cargo', 'bench', '--bench', 'metastable', '-j', '2']
    subprocess.run([*command, '--no-run'], cwd=ROOT, check=True)
    ref = parse(subprocess.check_output([str(exe), 'check'], env=env, text=True), ('check_',))
    native = parse(subprocess.check_output([*command, '--', '--check'], cwd=ROOT, text=True), ('check_',))
    assert ref.keys() == native.keys() and len(ref) == 170
    maxima = [0.]*4
    for key, expected in ref.items():
        assert len(expected) == len(native[key]) == 4
        for j, (a, b) in enumerate(zip(expected, native[key])):
            assert math.isfinite(a) and math.isfinite(b) and abs(a-b) <= ATOL[j]+RTOL*max(abs(a), abs(b)), (key, j, a, b)
            maxima[j] = max(maxima[j], abs(a-b))
    samples, checksums = {}, {}
    if not preflight_only:
        for sample in range(repeats+1):
            ref = parse(subprocess.check_output([str(exe)], env=env, text=True), ('source_',))
            native = parse(subprocess.check_output(command, cwd=ROOT, text=True), ('native_',))
            assert len(ref) == len(native) == 4
            for profile in range(1, 5):
                sk, nk = f'source_{profile}', f'native_{profile}'
                expected, actual = ref[sk], native[nk]
                assert len(expected) == len(actual) == 5 and min(expected[0], actual[0]) > 0 and all(map(math.isfinite, expected+actual))
                states = 10000 if profile <= 2 else 240
                for j, (a, b) in enumerate(zip(expected[1:], actual[1:])):
                    assert abs(a-b) <= states*ATOL[j]+RTOL*max(abs(a), abs(b)), (profile, j, a, b)
                if sample:
                    samples.setdefault(sk, []).append(expected[0])
                    samples.setdefault(nk, []).append(actual[0])
            checksums = dict(source=ref, native=native)
            print(f'Metastable benchmark {sample}/{repeats} (0 is warm-up)', flush=True)
    assert rust == fingerprints(sorted((ROOT/'rust').rglob('*.rs'))), 'Native sources changed during measurement'
    assert harness == fingerprints(harness_paths), 'Benchmark inputs changed during measurement'
    assert cargo == fingerprints(cargo_paths), 'Cargo inputs changed during measurement'
    if preflight_only:
        print(json.dumps(dict(source_states=170, native_comparisons=170, maximum_absolute_differences=maxima), indent=2))
        return
    report = dict(measured_at=datetime.now(timezone.utc).isoformat(), platform=platform.platform(), processor=platform.processor(),
        git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        rustc=subprocess.check_output(['rustc', '--version'], text=True).strip(), reference=json.loads((directory/'metadata.json').read_text()),
        compiler=subprocess.check_output([str(compiler), '--version'], text=True).splitlines()[0], compiler_flags=flags,
        native_sources_sha256=rust, benchmark_sources_sha256=harness, cargo_inputs_sha256=cargo,
        source_preflight_states=170, native_preflight_comparisons=170, maximum_preflight_absolute_differences=maxima,
        iterations_per_sample=dict(local_ps=10000, complete_curve=20), samples_ns=samples,
        median_ns={k: statistics.median(v) for k, v in samples.items()}, checksums=checksums,
        cargo_profile=dict(optimization='release', lto='thin', codegen_units=1),
        workload='Profiles 1/2: local homogeneous PS, PR methane/ethane and SRK carbon dioxide/propane, equimolar, 25 T/V-derived specifications, T=250..274 or 350..374 K, V=0.002..0.00248 m3/mol, seeds 1.02T/0.98V. Profiles 3/4: complete 12-point baseline metastable isentropes, PR methane from 180..182 K/6 MPa and equimolar methane/ethane from 250..252 K/10 MPa, 0.1 MPa minimum. Five source temperatures per curve profile. Model construction and target P/S preparation excluded; complete curve spinodal search, equilibrium PS initialization, continuation, returned states and entropy evaluation included. Preflight compares every returned T/V/reported P/actual S; timed samples check all four sums separately. Native local solve returns actual pressure/entropy internally. Fortran post-evaluates these returned quantities. Same numerical outputs; differing nonlinear iteration algorithms. This is scoped to these workloads and does not establish broad-model speedups or full failure-trajectory fidelity.')
    (ROOT/output).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['median_ns'], indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--compiler', required=True)
    p.add_argument('--repeats', type=int, default=5)
    p.add_argument('--output', default='docs/rust/metastable-benchmark.json')
    p.add_argument('--preflight-only', action='store_true')
    a = p.parse_args()
    benchmark(a.compiler, a.repeats, a.output, a.preflight_only)
