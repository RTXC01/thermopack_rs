"""Native-loop UV/SV benchmark with identical source-generated specifications.

Every state is compared before timing; checksums are verified on every sample.
Source SV receives TP phase arrays as required by its public low-level API.
Native SV receives the same T/P guess but initializes its own phase amounts.
"""
import argparse,hashlib,json,math,os,platform,statistics,subprocess
from datetime import datetime,timezone
from pathlib import Path
from reference import ROOT


def fingerprints(paths):
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths}


def benchmark(compiler,repeats,tcpr,output):
    if repeats<1:raise ValueError('At least one measured repetition required')
    compiler=Path(compiler).resolve();env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    out=ROOT/'.reference';exe=out/('benchmark_isochoric.exe' if os.name=='nt' else 'benchmark_isochoric')
    flags=['-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp']
    source=ROOT/'tools/rust_port/benchmark_isochoric_reference.f90'
    native_paths=sorted((ROOT/'rust').rglob('*.rs'));benchmark_paths=[Path(__file__),source]
    cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock']
    native_hashes=fingerprints(native_paths);benchmark_hashes=fingerprints(benchmark_paths);cargo_hashes=fingerprints(cargo_paths)
    subprocess.run([str(compiler),*flags,'-I',str(out),str(source),'@'+str(out/'objects.rsp'),'-o',str(exe)],env=env,cwd=out,check=True)
    ref='tcPR' if tcpr else 'DEFAULT'
    text=subprocess.check_output([str(exe),'inputs','unused',ref],env=env,text=True)
    data='\n'.join(line.removeprefix('input,') for line in text.splitlines() if line.startswith('input,'))+'\n'
    initialization_messages=[line for line in text.splitlines() if not line.startswith('input,')]
    inputs=out/'isochoric-benchmark-input.csv';inputs.write_bytes(data.encode())
    rows=[list(map(float,line.split(','))) for line in data.splitlines()];assert len(rows)==100 and all(len(row)==14 for row in rows)
    subprocess.run(['cargo','bench','--bench','isochoric','--no-run'],cwd=ROOT,check=True)
    native=['cargo','bench','--bench','isochoric','--','--inputs',str(inputs)]+(['--tcpr'] if tcpr else [])
    def run_reference(command):return subprocess.check_output([str(exe),command,str(inputs),ref],env=env,text=True)
    def run_rust(check=False):return subprocess.check_output(native+(['--check'] if check else []),cwd=ROOT,text=True)
    def parse(text,prefix):
        return {line.split(',')[0]:list(map(float,line.split(',')[1:])) for line in text.splitlines() if line.startswith(prefix)}
    expected=parse(run_reference('check'),'check_');actual=parse(run_rust(True),'check_')
    assert expected.keys()==actual.keys() and len(expected)==200
    maxima=[0.]*11
    for case,values in expected.items():
        assert len(values)==len(actual[case])==11
        for i,(a,b) in enumerate(zip(values,actual[case])):
            tolerance=0. if i==2 else ([2e-5,.5,0.,2e-7,2e-7]+[2e-7]*6)[i]+2e-7*max(abs(a),abs(b))
            assert math.isfinite(a) and math.isfinite(b) and abs(a-b)<=tolerance,(case,i,a,b,tolerance)
            maxima[i]=max(maxima[i],abs(a-b))
    samples={};checksums={}
    for i in range(repeats+1):
        reference=parse(run_reference('time'),('uv,','sv,'));rust=parse(run_rust(),('uv,','sv,'))
        assert reference.keys()==rust.keys()=={'uv','sv'}
        for name in ['uv','sv']:
            a,b=reference[name],rust[name]
            assert len(a)==len(b)==4 and a[0]>0. and b[0]>0. and all(map(math.isfinite,[*a,*b]))
            for j in range(1,4):assert abs(a[j]-b[j])<=2e-7*max(abs(a[j]),abs(b[j]),1.),(name,j,a,b)
            if i:
                samples.setdefault('fortran_'+name,[]).append(a[0]);samples.setdefault('rust_'+name,[]).append(b[0])
        if i:checksums=dict(fortran={k:v[1:] for k,v in reference.items()},rust={k:v[1:] for k,v in rust.items()})
        print(f'Isochoric benchmark {i}/{repeats} (0 is warm-up)',flush=True)
    assert native_hashes==fingerprints(sorted((ROOT/'rust').rglob('*.rs'))),'Native sources changed during measurements'
    assert benchmark_hashes==fingerprints(benchmark_paths),'Benchmark sources changed during measurements'
    assert cargo_hashes==fingerprints(cargo_paths),'Cargo inputs changed during measurements'
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),parameter_reference=ref,platform=platform.platform(),processor=platform.processor(),
                rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),reference=json.loads((out/'metadata.json').read_text()),
                reference_benchmark_flags=flags,reference_benchmark_compiler=subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],
                native_sources_sha256=native_hashes,benchmark_sources_sha256=benchmark_hashes,cargo_inputs_sha256=cargo_hashes,
                cargo_profile=dict(optimization='release',lto='thin',codegen_units=1),iterations_per_sample=2000,preflight_cases=200,
                reference_initialization_messages=initialization_messages,
                maximum_preflight_absolute_differences=maxima,input_sha256=hashlib.sha256(inputs.read_bytes()).hexdigest(),inputs=rows,
                samples_ns_per_state=samples,median_ns_per_state={k:statistics.median(v) for k,v in samples.items()},checksums=checksums,
                workload='PR C1/C2/C3 0.5/0.3/0.2; 100 source TP states from 180..259.2 K at 1 MPa; identical U/S/V inputs; UV wrapper defaults (zero T/P/phase arrays); SV T/P guess 0.9*T, 1.2*P; source SV gets TP phase arrays, native initializes its own; specification generation and independent closure/parity preflight outside timers')
    destination=ROOT/output if output else ROOT/'docs/rust'/('benchmark-isochoric-tcpr.json' if tcpr else 'benchmark-isochoric.json')
    destination.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['median_ns_per_state'],indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler',required=True);parser.add_argument('--repeats',type=int,default=5)
    parser.add_argument('--tcpr',action='store_true');parser.add_argument('--output')
    args=parser.parse_args();benchmark(args.compiler,args.repeats,args.tcpr,args.output)
