"""Native MBWR density timings, with independent per-state parity preflights."""
import argparse,json,math,os,platform,statistics,subprocess
from pathlib import Path
from datetime import datetime,timezone
from reference import ROOT
from benchmark_isochoric import fingerprints

NAMES=[f'{family}_{workload}_{mode}' for family in ['19','32'] for workload in ['liquid','vapor','supercritical'] for mode in ['fresh','prepared']]
def parse(output,prefixes):
    rows={}
    for line in output.splitlines():
        if line.startswith(prefixes):
            key,*values=line.split(',');assert key not in rows,key
            rows[key]=list(map(float,values))
    return rows

def benchmark(compiler,repeats,output,preflight_only):
    if repeats<1:raise ValueError('Positive repetitions required')
    compiler=Path(compiler).resolve();env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    out=ROOT/'.reference';source=ROOT/'tools/rust_port/benchmark_mbwr_density_reference.f90'
    exe=out/('benchmark_mbwr_density.exe' if os.name=='nt' else 'benchmark_mbwr_density')
    flags=['-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp']
    native=fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    harness_paths=[Path(__file__),source,ROOT/'tools/rust_port/benchmark_isochoric.py']
    harness=fingerprints(harness_paths);cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock'];cargo=fingerprints(cargo_paths)
    subprocess.run([str(compiler),*flags,'-I',str(out),str(source),'@'+str(out/'objects.rsp'),'-o',str(exe)],cwd=out,env=env,check=True)
    subprocess.run(['cargo','bench','--bench','mbwr_density','--no-run'],cwd=ROOT,check=True)
    command=['cargo','bench','--bench','mbwr_density']
    ref=parse(subprocess.check_output([str(exe),'check'],env=env,text=True),('check_',))
    actual=parse(subprocess.check_output([*command,'--','--check'],cwd=ROOT,text=True),('check_',))
    assert ref.keys()==actual.keys() and len(ref)==1200
    maximum=[0.]*4
    for name,a in ref.items():
        b=actual[name];assert len(a)==len(b)==4
        index=int(name.rsplit('_',1)[1]);p=5e5+index*5e3
        assert abs(a[2]-p)<=1e-5+2e-9*p and abs(b[2]-p)<=1e-5+2e-9*p
        assert a[0]>0 and a[3]>=0 and b[0]>0 and b[3]>=0
        for i,(v,w) in enumerate(zip(a,b)):
            assert math.isfinite(w) and abs(v-w)<=[1e-12,0.,1e-5,1e-8][i]+(0. if i==1 else 2e-9)*max(abs(v),abs(w)),(name,i,v,w)
            maximum[i]=max(maximum[i],abs(v-w))
    samples={};checksums={}
    if not preflight_only:
        for repeat in range(repeats+1):
            ref=parse(subprocess.check_output([str(exe)],env=env,text=True),tuple(n+',' for n in NAMES))
            actual=parse(subprocess.check_output(command,cwd=ROOT,text=True),tuple(n+',' for n in NAMES))
            assert ref.keys()==actual.keys()==set(NAMES)
            for name in NAMES:
                r,a=ref[name],actual[name];assert len(r)==len(a)==2
                assert r[0]>0 and a[0]>0 and all(map(math.isfinite,[*r,*a]))
                assert abs(r[1]-a[1])<=1e-8+2e-9*max(abs(r[1]),abs(a[1])),(name,r,a)
                if repeat:
                    samples.setdefault('fortran_'+name,[]).append(r[0]);samples.setdefault('rust_'+name,[]).append(a[0])
            checksums=dict(fortran={k:v[1] for k,v in ref.items()},rust={k:v[1] for k,v in actual.items()})
            print(f'MBWR density benchmark {repeat}/{repeats} (0 is warm-up)',flush=True)
    assert native==fingerprints(sorted((ROOT/'rust').rglob('*.rs'))),'Native sources changed during measurement'
    assert harness==fingerprints(harness_paths),'Benchmark sources changed during measurement'
    assert cargo==fingerprints(cargo_paths),'Cargo inputs changed during measurement'
    if preflight_only:
        print(json.dumps(dict(source_cases=1200,native_cases=1200,maximum_absolute_differences=maximum),indent=2));return
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),processor=platform.processor(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),reference=json.loads((out/'metadata.json').read_text()),
        compiler=subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],compiler_flags=flags,
        native_sources_sha256=native,benchmark_sources_sha256=harness,cargo_inputs_sha256=cargo,
        preflight_cases_per_implementation=1200,maximum_preflight_absolute_differences=maximum,
        iterations_per_sample=20000,samples_ns_per_state=samples,median_ns_per_state={k:statistics.median(v) for k,v in samples.items()},
        checksums=checksums,cargo_profile=dict(optimization='release',lto='thin',codegen_units=1),
        workload='Default methane MBWR19/32; liquid/vapor at 150 K and vapor at 300 K; 100 pressures 0.5..0.995 MPa; 20000 queries per row. Source model/critical/SRK initialization excluded in both implementations. Fresh includes temperature preparation per call (and native explicit-parameter validation/SRK preparation); prepared reuses it. Both return raw density, phase, pressure and pressure slope; checksum includes all four. No FFI or Python loop is timed.')
    (ROOT/output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['median_ns_per_state'],indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler',required=True);p.add_argument('--repeats',type=int,default=5)
    p.add_argument('--output',default='docs/rust/benchmark-mbwr-density.json');p.add_argument('--preflight-only',action='store_true')
    a=p.parse_args();benchmark(a.compiler,a.repeats,a.output,a.preflight_only)
