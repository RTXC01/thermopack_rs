"""Native-loop PeTS pressure-root timings with all path/state preflights."""
import argparse
import json
import math
import os
import platform
import statistics
import subprocess
from datetime import datetime,timezone
from pathlib import Path
from reference import ROOT
from benchmark_isochoric import fingerprints

WORKLOADS=['liquid','vapor','supercritical']
MODES=['fresh','prepared','baseline']


def benchmark(compiler,repeats,output,preflight_only):
    if repeats<1:raise ValueError('Positive repetition count required')
    compiler=Path(compiler).resolve();env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    out=ROOT/'.reference';source=ROOT/'tools/rust_port/benchmark_pets_roots_reference.f90'
    exe=out/('benchmark_pets_roots.exe' if os.name=='nt' else 'benchmark_pets_roots')
    flags=['-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp']
    native_paths=sorted((ROOT/'rust').rglob('*.rs'));benchmark_paths=[Path(__file__),source,ROOT/'tools/rust_port/benchmark_isochoric.py']
    cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock']
    native_hashes=fingerprints(native_paths);benchmark_hashes=fingerprints(benchmark_paths);cargo_hashes=fingerprints(cargo_paths)
    subprocess.run([str(compiler),*flags,'-I',str(out),str(source),'@'+str(out/'objects.rsp'),'-o',str(exe)],cwd=out,env=env,check=True)
    subprocess.run(['cargo','bench','--bench','pets_roots','--no-run'],cwd=ROOT,check=True)
    native=['cargo','bench','--bench','pets_roots']
    def parse(text,prefix):
        return {s.split(',')[0]:list(map(float,s.split(',')[1:])) for s in text.splitlines() if s.startswith(prefix)}
    reference=parse(subprocess.check_output([str(exe),'check'],env=env,text=True),'check_')
    rust=parse(subprocess.check_output([*native,'--','--check'],cwd=ROOT,text=True),'check_')
    assert len(reference)==300 and len(rust)==900
    maxima=[0.]*4
    for workload in WORKLOADS:
        for i in range(100):
            expected=reference[f'check_{workload}_{i}'];pressure=5e5+i*5e3
            assert len(expected)==4 and all(map(math.isfinite,expected)) and abs(expected[1]-pressure)<1e8*2.220446049250313e-16*pressure and expected[3]<0.
            for mode in MODES:
                actual=rust[f'check_{workload}_{mode}_{i}']
                assert len(actual)==4 and all(map(math.isfinite,actual)) and abs(actual[1]-pressure)<1e8*2.220446049250313e-16*pressure and actual[3]<0.
                for j,(a,b) in enumerate(zip(expected,actual)):
                    assert abs(a-b)<=[1e-10,.1,.01,10.][j]+2e-7*max(abs(a),abs(b)),(workload,mode,i,j,a,b)
                    maxima[j]=max(maxima[j],abs(a-b))
    samples={};checksums={}
    if not preflight_only:
        for i in range(repeats+1):
            ref=parse(subprocess.check_output([str(exe)],env=env,text=True),tuple(w+',' for w in WORKLOADS))
            actual=parse(subprocess.check_output(native,cwd=ROOT,text=True),tuple(w+'_'+m+',' for w in WORKLOADS for m in MODES))
            assert ref.keys()==set(WORKLOADS) and actual.keys()=={w+'_'+m for w in WORKLOADS for m in MODES}
            for workload in WORKLOADS:
                r=ref[workload];assert len(r)==2 and r[0]>0. and all(map(math.isfinite,r))
                if i:samples.setdefault('fortran_'+workload,[]).append(r[0])
                for mode in MODES:
                    name=workload+'_'+mode;v=actual[name]
                    assert len(v)==2 and v[0]>0. and all(map(math.isfinite,v))
                    assert abs(r[1]-v[1])<=1e-12+5e-8*max(abs(r[1]),abs(v[1])),(name,r,v)
                    if i:samples.setdefault('rust_'+name,[]).append(v[0])
            if i:checksums=dict(fortran={k:v[1] for k,v in ref.items()},rust={k:v[1] for k,v in actual.items()})
            print(f'PeTS roots benchmark {i}/{repeats} (0 is warm-up)',flush=True)
    assert native_hashes==fingerprints(sorted((ROOT/'rust').rglob('*.rs'))),'Native sources changed during measurement'
    assert benchmark_hashes==fingerprints(benchmark_paths),'Benchmark sources changed during measurement'
    assert cargo_hashes==fingerprints(cargo_paths),'Cargo inputs changed during measurement'
    if preflight_only:
        print(json.dumps(dict(reference_cases=300,native_path_comparisons=900,maximum_absolute_differences=maxima),indent=2));return
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),processor=platform.processor(),
        rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),reference=json.loads((out/'metadata.json').read_text()),
        reference_benchmark_compiler=subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],reference_benchmark_flags=flags,
        native_sources_sha256=native_hashes,benchmark_sources_sha256=benchmark_hashes,cargo_inputs_sha256=cargo_hashes,
        cargo_profile=dict(optimization='release',lto='thin',codegen_units=1),reference_preflight_cases=300,native_path_preflight_comparisons=900,
        maximum_preflight_absolute_differences=maxima,iterations_per_sample=20000,samples_ns_per_state=samples,
        median_ns_per_state={k:statistics.median(v) for k,v in samples.items()},checksums=checksums,
        workload='Default PeTS; liquid/vapor volume at 120 K and supercritical vapor seed at 300 K; 100 pressures 0.5..0.995 MPa; 20000 native-loop calls per row, model construction excluded. Fortran original seeded nonassociation solver; Rust fresh strict branches, reused isotherm strict branches, and source-compatible last-iterate path. Prepared-isotherm discovery excluded only from prepared timer.')
    destination=ROOT/output;destination.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['median_ns_per_state'],indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler',required=True);parser.add_argument('--repeats',type=int,default=5)
    parser.add_argument('--output',default='docs/rust/benchmark-pets-roots.json');parser.add_argument('--preflight-only',action='store_true')
    args=parser.parse_args();benchmark(args.compiler,args.repeats,args.output,args.preflight_only)
