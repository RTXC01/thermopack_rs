"""Native executable comparison of the same default local critical solves."""
import argparse,json,math,os,platform,statistics,subprocess
from pathlib import Path
from datetime import datetime,timezone
from reference import ROOT
from benchmark_isochoric import fingerprints
from benchmark_solid_water import parse
from validate_mixture_critical import ATOL,RTOL

def benchmark(compiler,repeats,output,preflight_only):
    if repeats<1:raise ValueError('Positive sample count required')
    compiler=Path(compiler).resolve();env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    directory=ROOT/'.reference';source=ROOT/'tools/rust_port/benchmark_mixture_critical_reference.f90'
    exe=directory/('benchmark_mixture_critical.exe' if os.name=='nt' else 'benchmark_mixture_critical')
    flags=['-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp']
    rust=fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    harness_paths=[Path(__file__),source,ROOT/'tools/rust_port/validate_mixture_critical.py',ROOT/'tools/rust_port/benchmark_isochoric.py',ROOT/'tools/rust_port/benchmark_solid_water.py']
    harness=fingerprints(harness_paths);cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock'];cargo=fingerprints(cargo_paths)
    subprocess.run([str(compiler),*flags,'-I',str(directory),str(source),'@'+str(directory/'objects.rsp'),'-o',str(exe)],cwd=directory,env=env,check=True)
    command=['cargo','bench','--bench','mixture_critical','-j','2']
    subprocess.run([*command,'--no-run'],cwd=ROOT,check=True)
    ref=parse(subprocess.check_output([str(exe),'check'],env=env,text=True),('check_',))
    native=parse(subprocess.check_output([*command,'--','--check'],cwd=ROOT,text=True),('check_',))
    assert len(ref)==len(native)==75
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
                assert len(expected)==len(actual)==4 and expected[0]>0 and actual[0]>0
                assert all(map(math.isfinite,expected+actual))
                for j,(a,b) in enumerate(zip(expected[1:],actual[1:])):
                    assert abs(a-b)<=200*ATOL[j]+RTOL*max(abs(a),abs(b)),(nk,j,a,b)
                if sample:samples.setdefault(sk,[]).append(expected[0]);samples.setdefault(nk,[]).append(actual[0])
            checksums=dict(source=ref,native=native)
            print(f'Mixture critical benchmark {sample}/{repeats} (0 is warm-up)',flush=True)
    assert rust==fingerprints(sorted((ROOT/'rust').rglob('*.rs'))),'Native sources changed during measurement'
    assert harness==fingerprints(harness_paths),'Benchmark inputs changed during measurement'
    assert cargo==fingerprints(cargo_paths),'Cargo inputs changed during measurement'
    if preflight_only:
        print(json.dumps(dict(source_states=75,native_comparisons=75,maximum_absolute_differences=maxima),indent=2));return
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),processor=platform.processor(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),reference=json.loads((directory/'metadata.json').read_text()),
        compiler=subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],compiler_flags=flags,
        native_sources_sha256=rust,benchmark_sources_sha256=harness,cargo_inputs_sha256=cargo,
        source_preflight_states=75,native_preflight_comparisons=75,maximum_preflight_absolute_differences=maxima,
        iterations_per_sample=200,samples_ns_per_state=samples,median_ns_per_state={k:statistics.median(v) for k,v in samples.items()},checksums=checksums,
        cargo_profile=dict(optimization='release',lto='thin',codegen_units=1),
        workload='Default local mixture critical solves, PR and SRK methane/ethane with x(C1)=0.1..0.9 and PR methane/ethane/propane with x(C1)=0.1..0.7, x(C2)=0.2. Each profile uses 25 compositions, 200 solves/sample, model construction outside timing, five samples after warm-up. Both executables include default seed calculation, all iterations, and T/V/P output in separate checksums. Python/FFI are excluded. Native classic mixing uses analytic Hessians and a scalar third-order hyperdual contraction; the source uses its original finite-difference critical equations. Excludes translations, Wong-Sandler, explicit seed/critical-locus/global-stability work and the rest of the library.')
    (ROOT/output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['median_ns_per_state'],indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler',required=True);p.add_argument('--repeats',type=int,default=5)
    p.add_argument('--output',default='docs/rust/benchmark-mixture-critical.json');p.add_argument('--preflight-only',action='store_true')
    a=p.parse_args();benchmark(a.compiler,a.repeats,a.output,a.preflight_only)
