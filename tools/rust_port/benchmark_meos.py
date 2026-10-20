"""Native pure MEOS state loops, all 24 fields qualified before timing."""
import argparse,json,math,os,platform,statistics,subprocess
from pathlib import Path
from datetime import datetime,timezone
from reference import ROOT
from benchmark_isochoric import fingerprints
from benchmark_solid_water import parse
from validate_meos import ATOL,RTOL

def benchmark(compiler,repeats,output,preflight_only):
    if repeats<1:raise ValueError('Positive sample count required')
    compiler=Path(compiler).resolve();env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    directory=ROOT/'.reference';source=ROOT/'tools/rust_port/benchmark_meos_reference.f90'
    exe=directory/('benchmark_meos.exe' if os.name=='nt' else 'benchmark_meos')
    flags=['-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp']
    rust=fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    harness_paths=[Path(__file__),source,ROOT/'tools/rust_port/validate_meos.py',ROOT/'tools/rust_port/benchmark_isochoric.py',ROOT/'tools/rust_port/benchmark_solid_water.py']
    harness=fingerprints(harness_paths);cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock'];cargo=fingerprints(cargo_paths)
    subprocess.run([str(compiler),*flags,'-I',str(directory),str(source),'@'+str(directory/'objects.rsp'),'-o',str(exe)],cwd=directory,env=env,check=True)
    command=['cargo','bench','--bench','meos']
    subprocess.run([*command,'--no-run'],cwd=ROOT,check=True)
    source_check=parse(subprocess.check_output([str(exe),'check'],env=env,text=True),('check_',))
    native_check=parse(subprocess.check_output([*command,'--','--check'],cwd=ROOT,text=True),('check_',))
    assert len(source_check)==600 and len(native_check)==1200
    maxima=[0.]*24
    for model in range(1,7):
        for i in range(100):
            expected=source_check[f'check_{model}_{i}'];assert len(expected)==24
            for mode in ['fresh','prepared']:
                actual=native_check[f'check_{mode}_{model}_{i}'];assert len(actual)==24
                for j,(a,b) in enumerate(zip(expected,actual)):
                    assert math.isfinite(a) and math.isfinite(b) and abs(a-b)<=ATOL[j]+RTOL*max(abs(a),abs(b)),(mode,model,i,j,a,b)
                    maxima[j]=max(maxima[j],abs(a-b))
    samples={};checksums={}
    if not preflight_only:
        for sample in range(repeats+1):
            ref=parse(subprocess.check_output([str(exe)],env=env,text=True),('source_',))
            native=parse(subprocess.check_output(command,cwd=ROOT,text=True),('fresh_','prepared_'))
            assert len(ref)==6 and len(native)==12
            for model in range(1,7):
                key=f'source_{model}';expected=ref[key];assert expected[0]>0 and all(map(math.isfinite,expected))
                if sample:samples.setdefault(key,[]).append(expected[0])
                for mode in ['fresh','prepared']:
                    key=f'{mode}_{model}';actual=native[key];assert actual[0]>0 and all(map(math.isfinite,actual))
                    assert abs(actual[1]-expected[1])<=1e-5+1e-9*max(abs(actual[1]),abs(expected[1])),(key,actual,expected)
                    if sample:samples.setdefault(key,[]).append(actual[0])
            checksums=dict(source=ref,native=native)
            print(f'MEOS benchmark {sample}/{repeats} (0 is warm-up)',flush=True)
    assert rust==fingerprints(sorted((ROOT/'rust').rglob('*.rs'))),'Native sources changed during measurement'
    assert harness==fingerprints(harness_paths),'Benchmark inputs changed during measurement'
    assert cargo==fingerprints(cargo_paths),'Cargo inputs changed during measurement'
    if preflight_only:
        print(json.dumps(dict(source_states=600,native_comparisons=1200,maximum_absolute_differences=maxima),indent=2));return
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),processor=platform.processor(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),reference=json.loads((directory/'metadata.json').read_text()),
        compiler=subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],compiler_flags=flags,
        native_sources_sha256=rust,benchmark_sources_sha256=harness,cargo_inputs_sha256=cargo,
        source_preflight_states=600,native_preflight_comparisons=1200,maximum_preflight_absolute_differences=maxima,
        iterations_per_sample=100000,samples_ns_per_state=samples,median_ns_per_state={k:statistics.median(v) for k,v in samples.items()},checksums=checksums,
        cargo_profile=dict(optimization='release',lto='thin',codegen_units=1),
        workload='Six pure MEOS profiles: C3, R134A, N-H2, LJ argon, LJTS argon and raw reduced LJTS. T=1.4 Tc; 100 densities delta=0.05..1.535 repeated over 100000 states/sample. All 24 outputs checked before timing and included in checksums; source calls its original ideal/residual jets, pressure, TV alpha, Cv, Cp and sound methods, retaining its internal NIST temperature caches. Rust shares one residual jet and a cached ideal jet. Fresh prepares temperature per state; prepared reuses it outside timer. Model construction, Python and FFI are outside timers. Source-scaled sound is selected for comparison. This excludes phase solvers and critical/dilute-limit workloads.')
    (ROOT/output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['median_ns_per_state'],indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler',required=True);p.add_argument('--repeats',type=int,default=5)
    p.add_argument('--output',default='docs/rust/benchmark-meos.json');p.add_argument('--preflight-only',action='store_true')
    a=p.parse_args();benchmark(a.compiler,a.repeats,a.output,a.preflight_only)
