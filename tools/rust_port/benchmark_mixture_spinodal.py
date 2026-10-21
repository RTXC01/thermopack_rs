"""Native executable comparison of the same default local spinodal solves."""
import argparse,json,math,os,platform,statistics,subprocess
from pathlib import Path
from datetime import datetime,timezone
from reference import ROOT
from benchmark_isochoric import fingerprints
from benchmark_solid_water import parse
from validate_mixture_spinodal import ATOL,RTOL

def benchmark(compiler,repeats,output,preflight_only):
    if repeats<1:raise ValueError('Positive sample count required')
    compiler=Path(compiler).resolve();env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    directory=ROOT/'.reference';source=ROOT/'tools/rust_port/benchmark_mixture_spinodal_reference.f90'
    exe=directory/('benchmark_mixture_spinodal.exe' if os.name=='nt' else 'benchmark_mixture_spinodal')
    flags=['-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp']
    rust=fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    harness_paths=[Path(__file__),source,ROOT/'tools/rust_port/validate_mixture_spinodal.py',ROOT/'tools/rust_port/benchmark_isochoric.py',ROOT/'tools/rust_port/benchmark_solid_water.py']
    harness=fingerprints(harness_paths);cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock'];cargo=fingerprints(cargo_paths)
    subprocess.run([str(compiler),*flags,'-I',str(directory),str(source),'@'+str(directory/'objects.rsp'),'-o',str(exe)],cwd=directory,env=env,check=True)
    command=['cargo','bench','--bench','mixture_spinodal','-j','2']
    subprocess.run([*command,'--no-run'],cwd=ROOT,check=True)
    ref=parse(subprocess.check_output([str(exe),'check'],env=env,text=True),('check_',))
    native=parse(subprocess.check_output([*command,'--','--check'],cwd=ROOT,text=True),('check_',))
    assert len(ref)==len(native) and len(ref)>100
    preflight_count=len(ref)
    maxima=[0.]*3
    for key,expected in ref.items():
        actual=native[key];assert len(expected)==len(actual)==3
        for j,(a,b) in enumerate(zip(expected,actual)):
            assert math.isfinite(a) and math.isfinite(b) and abs(a-b)<=ATOL[j]+RTOL*max(abs(a),abs(b)),(key,j,a,b)
            maxima[j]=max(maxima[j],abs(a-b))
    samples={};checksums={}
    if not preflight_only:
        for sample in range(repeats+1):
            ref=parse(subprocess.check_output([str(exe)],env=env,text=True),('source_',))
            native=parse(subprocess.check_output(command,cwd=ROOT,text=True),('native_',))
            assert len(ref)==len(native)==3
            for profile in range(1,4):
                sk=f'source_{profile}';nk=f'native_{profile}';expected=ref[sk];actual=native[nk]
                assert len(expected)==len(actual)==5 and expected[1]==actual[1] and expected[0]>0 and actual[0]>0
                assert all(map(math.isfinite,expected+actual))
                for j,(a,b) in enumerate(zip(expected[2:],actual[2:])):
                    assert abs(a-b)<=10*ATOL[j]+RTOL*max(abs(a),abs(b)),(nk,j,a,b)
                if sample:samples.setdefault(sk,[]).append(expected[0]);samples.setdefault(nk,[]).append(actual[0])
            checksums=dict(source=ref,native=native)
            print(f'Mixture spinodal benchmark {sample}/{repeats} (0 is warm-up)',flush=True)
    assert rust==fingerprints(sorted((ROOT/'rust').rglob('*.rs'))),'Native sources changed during measurement'
    assert harness==fingerprints(harness_paths),'Benchmark inputs changed during measurement'
    assert cargo==fingerprints(cargo_paths),'Cargo inputs changed during measurement'
    if preflight_only:
        print(json.dumps(dict(source_states=preflight_count,native_comparisons=preflight_count,maximum_absolute_differences=maxima),indent=2));return
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),processor=platform.processor(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),reference=json.loads((directory/'metadata.json').read_text()),
        compiler=subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],compiler_flags=flags,
        native_sources_sha256=rust,benchmark_sources_sha256=harness,cargo_inputs_sha256=cargo,
        source_preflight_states=preflight_count,native_preflight_comparisons=preflight_count,maximum_preflight_absolute_differences=maxima,
        iterations_per_sample=10,samples_ns_per_curve=samples,median_ns_per_curve={k:statistics.median(v) for k,v in samples.items()},checksums=checksums,
        cargo_profile=dict(optimization='release',lto='thin',codegen_units=1),
        workload='Complete mixture spinodal curves at fixed composition: PR/SRK methane+ethane 0.5/0.5 and PR methane+ethane+propane 0.2/0.3/0.5. Initial pressure 0.1 MPa, ln(volume) step 0.04, minimum vapor temperature 0.6 times the composition-weighted critical temperatures. Ten complete curves per sample, five samples after warm-up. Model construction, Python and FFI excluded. Both native executables include initialization, continuation and endpoint refinement. Every T/V/P coordinate and point count is checked before timing and carried in separate timed checksums. Rust uses the explicit baseline API, including the original stale endpoint pressure convention. Physical pressure reporting, other mixtures/models/parameters and source pure grids are outside this measurement.')
    (ROOT/output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['median_ns_per_curve'],indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler',required=True);p.add_argument('--repeats',type=int,default=5)
    p.add_argument('--output',default='docs/rust/benchmark-mixture-spinodal.json');p.add_argument('--preflight-only',action='store_true')
    a=p.parse_args();benchmark(a.compiler,a.repeats,a.output,a.preflight_only)
