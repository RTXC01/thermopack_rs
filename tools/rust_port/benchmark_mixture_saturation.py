"""Native mixture bubble/dew loops, all returned state fields checked before timing."""
import argparse,json,math,os,platform,statistics,subprocess
from pathlib import Path
from datetime import datetime,timezone
from reference import ROOT
from benchmark_isochoric import fingerprints
from benchmark_solid_water import parse
from validate_mixture_saturation import ATOL,RTOL
ATOL=ATOL+[1e-9]*4

def benchmark(compiler,repeats,output,preflight_only):
    if repeats<1:raise ValueError('Positive sample count required')
    compiler=Path(compiler).resolve();env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    directory=ROOT/'.reference';source=ROOT/'tools/rust_port/benchmark_mixture_saturation_reference.f90'
    exe=directory/('benchmark_mixture_saturation.exe' if os.name=='nt' else 'benchmark_mixture_saturation')
    flags=['-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp']
    rust=fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    harness_paths=[Path(__file__),source,ROOT/'tools/rust_port/validate_mixture_saturation.py',ROOT/'tools/rust_port/benchmark_isochoric.py',ROOT/'tools/rust_port/benchmark_solid_water.py']
    harness=fingerprints(harness_paths);cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock'];cargo=fingerprints(cargo_paths)
    subprocess.run([str(compiler),*flags,'-I',str(directory),str(source),'@'+str(directory/'objects.rsp'),'-o',str(exe)],cwd=directory,env=env,check=True)
    command=['cargo','bench','--bench','mixture_saturation']
    subprocess.run([*command,'--no-run'],cwd=ROOT,check=True)
    source_check=parse(subprocess.check_output([str(exe),'check'],env=env,text=True),('check_',))
    native_check=parse(subprocess.check_output([*command,'--','--check'],cwd=ROOT,text=True),('check_',))
    assert len(source_check)==800 and len(native_check)==800
    maxima=[0.]*8
    for profile in [1,2]:
        for mode in ['BP','DP','BT','DT']:
            for i in range(100):
                key=f'check_{profile}_{mode}_{i}';expected=source_check[key];actual=native_check[key]
                assert len(expected)==len(actual)==8
                for j,(a,b) in enumerate(zip(expected,actual)):
                    assert math.isfinite(a) and math.isfinite(b) and abs(a-b)<=ATOL[j]+RTOL*max(abs(a),abs(b)),(key,j,a,b)
                    maxima[j]=max(maxima[j],abs(a-b))
    samples={};checksums={}
    if not preflight_only:
        for sample in range(repeats+1):
            ref=parse(subprocess.check_output([str(exe)],env=env,text=True),('source_',))
            native=parse(subprocess.check_output(command,cwd=ROOT,text=True),('native_',))
            assert len(ref)==len(native)==8
            for profile in [1,2]:
                for mode in ['BP','DP','BT','DT']:
                    sk=f'source_{profile}_{mode}';nk=f'native_{profile}_{mode}'
                    expected=ref[sk];actual=native[nk]
                    assert len(expected)==len(actual)==9
                    assert expected[0]>0 and actual[0]>0 and all(map(math.isfinite,expected+actual))
                    for j,(a,b) in enumerate(zip(expected[1:],actual[1:])):
                        assert abs(a-b)<=1000*ATOL[j]+RTOL*max(abs(a),abs(b)),(nk,j,a,b)
                    if sample:
                        samples.setdefault(sk,[]).append(expected[0]);samples.setdefault(nk,[]).append(actual[0])
            checksums=dict(source=ref,native=native)
            print(f'Mixture saturation benchmark {sample}/{repeats} (0 is warm-up)',flush=True)
    assert rust==fingerprints(sorted((ROOT/'rust').rglob('*.rs'))),'Native sources changed during measurement'
    assert harness==fingerprints(harness_paths),'Benchmark inputs changed during measurement'
    assert cargo==fingerprints(cargo_paths),'Cargo inputs changed during measurement'
    if preflight_only:
        print(json.dumps(dict(source_states=800,native_comparisons=800,maximum_absolute_differences=maxima),indent=2));return
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),processor=platform.processor(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),reference=json.loads((directory/'metadata.json').read_text()),
        compiler=subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],compiler_flags=flags,
        native_sources_sha256=rust,benchmark_sources_sha256=harness,cargo_inputs_sha256=cargo,
        source_preflight_states=800,native_preflight_comparisons=800,maximum_preflight_absolute_differences=maxima,
        iterations_per_sample=1000,samples_ns_per_state=samples,median_ns_per_state={k:statistics.median(v) for k,v in samples.items()},checksums=checksums,
        cargo_profile=dict(optimization='release',lto='thin',codegen_units=1),
        workload='PR/SRK methane+ethane at fixed phase composition 0.5/0.5. All four bubble/dew pressure/temperature methods; 100 T points 160..199.6 K or P points 0.1..3.07 MPa, 1000 solves/sample, five samples after warm-up. Both original safe wrappers and native default initialization solve each point independently. Model construction, Python and FFI are excluded. T/P, both volumes and all four phase fractions are preflighted and included in separate timed checksums. Rust also computes its branch tangent; the source does not return that tangent. This excludes translated/activity models, critical/retrograde fallback workloads and the full library.')
    (ROOT/output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['median_ns_per_state'],indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler',required=True);p.add_argument('--repeats',type=int,default=5)
    p.add_argument('--output',default='docs/rust/benchmark-mixture-saturation.json');p.add_argument('--preflight-only',action='store_true')
    a=p.parse_args();benchmark(a.compiler,a.repeats,a.output,a.preflight_only)
