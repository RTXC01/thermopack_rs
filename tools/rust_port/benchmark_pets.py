"""Native-loop full PeTS Helmholtz derivative benchmark with per-state preflight."""
import argparse,json,math,os,platform,statistics,subprocess
from pathlib import Path
from datetime import datetime,timezone
from reference import ROOT
from benchmark_isochoric import fingerprints


def benchmark(compiler,repeats):
    if repeats<1:raise ValueError('Positive repetition count required')
    compiler=Path(compiler).resolve();env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    out=ROOT/'.reference';source=ROOT/'tools/rust_port/benchmark_pets_reference.f90'
    exe=out/('benchmark_pets.exe' if os.name=='nt' else 'benchmark_pets')
    flags=['-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp']
    native_paths=sorted((ROOT/'rust').rglob('*.rs'));benchmark_paths=[Path(__file__),source,ROOT/'tools/rust_port/benchmark_isochoric.py']
    cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock']
    native_hashes=fingerprints(native_paths);benchmark_hashes=fingerprints(benchmark_paths);cargo_hashes=fingerprints(cargo_paths)
    subprocess.run([str(compiler),*flags,'-I',str(out),str(source),'@'+str(out/'objects.rsp'),'-o',str(exe)],cwd=out,env=env,check=True)
    subprocess.run(['cargo','bench','--bench','pets','--no-run'],cwd=ROOT,check=True)
    rust_command=['cargo','bench','--bench','pets']
    def parse(text,prefix):
        return {line.split(',')[0]:list(map(float,line.split(',')[1:])) for line in text.splitlines() if line.startswith(prefix)}
    reference=parse(subprocess.check_output([str(exe),'check'],env=env,text=True),'check_')
    rust=parse(subprocess.check_output([*rust_command,'--','--check'],cwd=ROOT,text=True),'check_')
    assert len(reference)==1000 and reference.keys()==rust.keys()
    maxima=[0.]*10
    for case,values in reference.items():
        assert len(values)==len(rust[case])==10
        for i,(a,b) in enumerate(zip(values,rust[case])):
            assert math.isfinite(a) and math.isfinite(b) and abs(a-b)<=1e-9+5e-10*max(abs(a),abs(b)),(case,i,a,b)
            maxima[i]=max(maxima[i],abs(a-b))
    samples={};checksums={}
    for i in range(repeats+1):
        ref=parse(subprocess.check_output([str(exe)],env=env,text=True),'helmholtz,')['helmholtz']
        native=parse(subprocess.check_output(rust_command,cwd=ROOT,text=True),('prepared,','with_temperature_preparation,'))
        assert len(ref)==11 and ref[0]>0. and all(map(math.isfinite,ref))
        for name,values in native.items():
            assert len(values)==11 and values[0]>0. and all(map(math.isfinite,values))
            for a,b in zip(ref[1:],values[1:]):assert abs(a-b)<=2e-9*max(abs(a),abs(b),1.),(name,a,b)
            if i:samples.setdefault('rust_'+name,[]).append(values[0])
        assert native.keys()=={'prepared','with_temperature_preparation'}
        if i:
            samples.setdefault('fortran',[]).append(ref[0]);checksums=dict(fortran=ref[1:],rust={k:v[1:] for k,v in native.items()})
        print(f'PeTS benchmark {i}/{repeats} (0 is warm-up)',flush=True)
    assert native_hashes==fingerprints(sorted((ROOT/'rust').rglob('*.rs'))),'Native sources changed during measurement'
    assert benchmark_hashes==fingerprints(benchmark_paths),'Benchmark sources changed during measurement'
    assert cargo_hashes==fingerprints(cargo_paths),'Cargo inputs changed during measurement'
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),processor=platform.processor(),
                rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),reference=json.loads((out/'metadata.json').read_text()),
                reference_benchmark_compiler=subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],reference_benchmark_flags=flags,
                native_sources_sha256=native_hashes,benchmark_sources_sha256=benchmark_hashes,cargo_inputs_sha256=cargo_hashes,
                cargo_profile=dict(optimization='release',lto='thin',codegen_units=1),preflight_cases=1000,iterations_per_sample=1000000,
                maximum_preflight_absolute_differences=maxima,samples_ns_per_state=samples,median_ns_per_state={k:statistics.median(v) for k,v in samples.items()},checksums=checksums,
                workload='Default PeTS (sigma 3.42e-10 m, epsilon/k 136 K); 300 K, 1 mol; 1000 volumes from 5e-5 through 1.049e-3 m3; full F and all nine first/second T,V,n derivatives; native Fortran kernel call, Rust with/without reusable temperature preparation; model construction excluded')
    (ROOT/'docs/rust/benchmark-pets.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['median_ns_per_state'],indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler',required=True);parser.add_argument('--repeats',type=int,default=5)
    args=parser.parse_args();benchmark(args.compiler,args.repeats)
