"""Native-loop ice Ih state benchmark, all twelve outputs checked first."""
import argparse,json,math,os,platform,statistics,subprocess
from pathlib import Path
from datetime import datetime,timezone
from reference import ROOT
from benchmark_isochoric import fingerprints
from validate_solid_water import ATOL,RTOL

def parse(text,prefix):
    result={}
    for line in text.splitlines():
        if line.startswith(prefix):
            name,*values=line.split(',');assert name not in result
            result[name]=list(map(float,values))
    return result
def benchmark(compiler,repeats,output,preflight_only):
    if repeats<1:raise ValueError('Positive sample count required')
    compiler=Path(compiler).resolve();env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    directory=ROOT/'.reference';source=ROOT/'tools/rust_port/benchmark_solid_water_reference.f90'
    exe=directory/('benchmark_solid_water.exe' if os.name=='nt' else 'benchmark_solid_water')
    flags=['-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp']
    rust=fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    harness_paths=[Path(__file__),source,ROOT/'tools/rust_port/validate_solid_water.py',ROOT/'tools/rust_port/benchmark_isochoric.py']
    harness=fingerprints(harness_paths);cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock'];cargo=fingerprints(cargo_paths)
    subprocess.run([str(compiler),*flags,'-I',str(directory),str(source),'@'+str(directory/'objects.rsp'),'-o',str(exe)],cwd=directory,env=env,check=True)
    subprocess.run(['cargo','bench','--bench','solid_water','--no-run'],cwd=ROOT,check=True)
    command=['cargo','bench','--bench','solid_water']
    ref=parse(subprocess.check_output([str(exe),'check'],env=env,text=True),('check_',))
    native=parse(subprocess.check_output([*command,'--','--check'],cwd=ROOT,text=True),('check_',))
    assert len(ref)==100 and len(native)==200
    maximum=[0.]*12
    for i in range(100):
        expected=ref[f'check_{i}'];assert len(expected)==12
        for mode in ['fresh','prepared']:
            actual=native[f'check_{mode}_{i}'];assert len(actual)==12
            for j,(a,b) in enumerate(zip(expected,actual)):
                assert math.isfinite(b) and abs(a-b)<=ATOL[j]+RTOL*max(abs(a),abs(b)),(i,mode,j,a,b)
                maximum[j]=max(maximum[j],abs(a-b))
    samples={};checksums={}
    if not preflight_only:
        for sample in range(repeats+1):
            ref=parse(subprocess.check_output([str(exe)],env=env,text=True),('source,',))['source']
            native=parse(subprocess.check_output(command,cwd=ROOT,text=True),('fresh,','prepared,'))
            assert ref[0]>0 and all(map(math.isfinite,ref))
            for mode in ['fresh','prepared']:
                values=native[mode];assert values[0]>0 and all(map(math.isfinite,values))
                assert abs(values[1]-ref[1])<=1e-6+1e-10*max(abs(values[1]),abs(ref[1])),(mode,values,ref)
                if sample:samples.setdefault('rust_'+mode,[]).append(values[0])
            if sample:samples.setdefault('fortran_source',[]).append(ref[0])
            checksums=dict(source=ref[1],native={k:v[1] for k,v in native.items()})
            print(f'Ice state benchmark {sample}/{repeats} (0 is warm-up)',flush=True)
    assert rust==fingerprints(sorted((ROOT/'rust').rglob('*.rs'))),'Native sources changed during measurement'
    assert harness==fingerprints(harness_paths),'Benchmark inputs changed during measurement'
    assert cargo==fingerprints(cargo_paths),'Cargo inputs changed during measurement'
    if preflight_only:
        print(json.dumps(dict(source_states=100,native_comparisons=200,maximum_absolute_differences=maximum),indent=2));return
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),processor=platform.processor(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),reference=json.loads((directory/'metadata.json').read_text()),
        compiler=subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],compiler_flags=flags,
        native_sources_sha256=rust,benchmark_sources_sha256=harness,cargo_inputs_sha256=cargo,
        preflight_source_states=100,native_preflight_comparisons=200,maximum_preflight_absolute_differences=maximum,
        iterations_per_sample=500000,samples_ns_per_state=samples,median_ns_per_state={k:statistics.median(v) for k,v in samples.items()},
        checksums=checksums,cargo_profile=dict(optimization='release',lto='thin',codegen_units=1),
        workload='Ice Ih at 250 K, 100 pressures 0.1..99.1 MPa; 500000 calls per sample, source/reference initialization excluded. Source evaluates its twelve separate public functions; Rust shares one analytic Gibbs jet, with fresh temperature preparation or reused isotherm. All twelve outputs including the source-scaled Cp are checked and included in checksums. No Python or FFI loop is timed; prepared temperature work excluded only from prepared timer.')
    (ROOT/output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['median_ns_per_state'],indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler',required=True);p.add_argument('--repeats',type=int,default=5)
    p.add_argument('--output',default='docs/rust/benchmark-solid-water.json');p.add_argument('--preflight-only',action='store_true')
    a=p.parse_args();benchmark(a.compiler,a.repeats,a.output,a.preflight_only)
