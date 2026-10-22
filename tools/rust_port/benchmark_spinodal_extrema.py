"""Native-executable fresh and prepared pure spinodal extrema measurements."""
import argparse,json,math,os,platform,statistics,subprocess
from pathlib import Path
from datetime import datetime,timezone
from reference import ROOT
from benchmark_isochoric import fingerprints
from benchmark_solid_water import parse
from validate_spinodal_extrema import ATOL,RTOL

def benchmark(compiler,repeats,output,preflight_only):
    if repeats<1:raise ValueError('Positive sample count required')
    compiler=Path(compiler).resolve();env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    directory=ROOT/'.reference';source=ROOT/'tools/rust_port/benchmark_spinodal_extrema_reference.f90'
    exe=directory/('benchmark_spinodal_extrema.exe' if os.name=='nt' else 'benchmark_spinodal_extrema')
    flags=['-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp']
    rust=fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    harness_paths=[Path(__file__).resolve(),source,ROOT/'tools/rust_port/validate_spinodal_extrema.py',ROOT/'tools/rust_port/benchmark_isochoric.py',ROOT/'tools/rust_port/benchmark_solid_water.py']
    harness=fingerprints(harness_paths);cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock'];cargo=fingerprints(cargo_paths)
    subprocess.run([str(compiler),*flags,'-I',str(directory),str(source),'@'+str(directory/'objects.rsp'),'-o',str(exe)],cwd=directory,env=env,check=True)
    command=['cargo','bench','--bench','spinodal_extrema','-j','2']
    subprocess.run([*command,'--no-run'],cwd=ROOT,check=True)
    ref=parse(subprocess.check_output([str(exe),'check'],env=env,text=True),('check_',))
    native=parse(subprocess.check_output([*command,'--','--check'],cwd=ROOT,text=True),('check_',))
    assert len(ref)==6 and len(native)==12
    maxima=[0.]*6
    for profile in [1,2]:
        for mode in ['fresh','prepared']:
            for i in range(3):
                expected=ref[f'check_{profile}_{i}'];actual=native[f'check_{profile}_{mode}_{i}']
                assert len(expected)==len(actual)==6
                for j,(a,b) in enumerate(zip(expected,actual)):
                    assert math.isfinite(a) and math.isfinite(b) and abs(a-b)<=ATOL[j%3]+RTOL*max(abs(a),abs(b)),(profile,mode,i,j,a,b)
                    maxima[j]=max(maxima[j],abs(a-b))
    samples={};checksums={}
    if not preflight_only:
        for sample in range(repeats+1):
            ref=parse(subprocess.check_output([str(exe)],env=env,text=True),('source_',))
            native=parse(subprocess.check_output(command,cwd=ROOT,text=True),('native_','setup_'))
            assert len(ref)==2 and len(native)==6
            for profile in [1,2]:
                sk=f'source_{profile}';expected=ref[sk]
                assert len(expected)==7 and expected[0]>0 and all(map(math.isfinite,expected))
                if sample:samples.setdefault(sk,[]).append(expected[0])
                for mode in ['fresh','prepared']:
                    nk=f'native_{profile}_{mode}';actual=native[nk]
                    assert len(actual)==7 and actual[0]>0 and all(map(math.isfinite,actual))
                    for j,(a,b) in enumerate(zip(expected[1:],actual[1:])):
                        assert abs(a-b)<=100*ATOL[j%3]+RTOL*max(abs(a),abs(b)),(nk,j,a,b)
                    if sample:samples.setdefault(nk,[]).append(actual[0])
                setup=f'setup_{profile}'
                assert len(native[setup])==1 and native[setup][0]>0
                if sample:samples.setdefault(setup,[]).append(native[setup][0])
            checksums=dict(source=ref,native=native)
            print(f'Spinodal extrema benchmark {sample}/{repeats} (0 is warm-up)',flush=True)
    assert rust==fingerprints(sorted((ROOT/'rust').rglob('*.rs'))),'Native sources changed during measurement'
    assert harness==fingerprints(harness_paths),'Benchmark inputs changed during measurement'
    assert cargo==fingerprints(cargo_paths),'Cargo inputs changed during measurement'
    if preflight_only:
        print(json.dumps(dict(source_extrema_pairs=6,native_comparisons=12,maximum_absolute_differences=maxima),indent=2));return
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),processor=platform.processor(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),reference=json.loads((directory/'metadata.json').read_text()),
        compiler=subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],compiler_flags=flags,
        native_sources_sha256=rust,benchmark_sources_sha256=harness,cargo_inputs_sha256=cargo,
        source_preflight_extrema_pairs=6,native_preflight_comparisons=12,maximum_preflight_absolute_differences=maxima,
        iterations_per_sample=100,samples_ns=samples,median_ns={k:statistics.median(v) for k,v in samples.items()},checksums=checksums,
        cargo_profile=dict(optimization='release',lto='thin',codegen_units=1),
        workload='Pure spinodal property extrema for PR methane and SRK nitrogen at initial pressure 1 MPa; cycle entropy, enthalpy and internal energy, 100 queries/sample, five samples after one warm-up. Original routine rebuilds two 100-point branches per query. Fresh Rust includes the complete grid and source candidate/refinement/reporting conventions. Prepared Rust caches that grid and all six property arrays; one-time setup is excluded from query timing and recorded separately. Model construction and Python/FFI excluded. Preflight compares every Tmin/Vmin/Pmin/Tmax/Vmax/Pmax and all six sums are checked in every timed sample. These three properties have unambiguous ordering in this workload; pressure-label ties are not benchmarked. Prepared ratios are amortized workflow measurements, not equal one-shot work or broad-model claims. Baseline mixed-coordinate minimum-pressure convention and uncorrected trial volumes are retained.')
    (ROOT/output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['median_ns'],indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler',required=True);p.add_argument('--repeats',type=int,default=5)
    p.add_argument('--output',default='docs/rust/spinodal-extrema-benchmark.json');p.add_argument('--preflight-only',action='store_true')
    a=p.parse_args();benchmark(a.compiler,a.repeats,a.output,a.preflight_only)
