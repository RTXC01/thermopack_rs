"""Complete selected-component pure curves, including both endpoint volumes."""
import argparse,csv,hashlib,json,math,os,platform,statistics,subprocess
from datetime import datetime,timezone
from pathlib import Path
from reference import ROOT
from benchmark_isochoric import fingerprints
from benchmark_solid_water import parse

def benchmark(compiler,repeats,output,preflight_only):
    if repeats<1:raise ValueError('Positive sample count required')
    compiler=Path(compiler).resolve();env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    directory=ROOT/'.reference';source=ROOT/'tools/rust_port/benchmark_pure_curve_reference.f90'
    exe=directory/('benchmark_pure_curve.exe' if os.name=='nt' else 'benchmark_pure_curve')
    flags=['-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp']
    rust=fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    harness_paths=[Path(__file__).resolve(),source,ROOT/'tools/rust_port/benchmark_isochoric.py',ROOT/'tools/rust_port/benchmark_solid_water.py',ROOT/'tests/fixtures/pure_curve.csv',ROOT/'tests/fixtures/pure_curve.metadata.json']
    harness=fingerprints(harness_paths);cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock'];cargo=fingerprints(cargo_paths)
    fixture=ROOT/'tests/fixtures/pure_curve.csv';meta=json.loads(fixture.with_suffix('.metadata.json').read_text())
    assert hashlib.sha256(fixture.read_bytes()).hexdigest()==meta['fixture_sha256']
    assert meta['library_sha256']==json.loads((directory/'metadata.json').read_text())['library_sha256']
    signatures={1:('PR','C1','classic','P','0'),2:('SRK','CO2','classic','T','1'),3:('PR','H2O','translated','P','0'),4:('PR','H2','quantum','P','1')}
    endpoints={}
    for row in csv.DictReader(fixture.open(),delimiter=';'):
        f=row['case'].split('|');o=f[4].split(':')
        for profile,signature in signatures.items():
            if (f[1],f[2],f[3],o[1],o[5])==signature and o[3]=='20000.0' and o[4]=='100' and float(row['critical_volume_radius'])>0.:
                endpoints[profile]=row
    assert len(endpoints)==4
    subprocess.run([str(compiler),*flags,'-I',str(directory),str(source),'@'+str(directory/'objects.rsp'),'-o',str(exe)],cwd=directory,env=env,check=True)
    command=['cargo','bench','--bench','pure_curve','-j','2'];subprocess.run([*command,'--no-run'],cwd=ROOT,check=True)
    ref=parse(subprocess.check_output([str(exe),'check'],env=env,text=True),('check_',))
    native=parse(subprocess.check_output([*command,'--','--check'],cwd=ROOT,text=True),('check_',))
    assert ref.keys()==native.keys() and len(ref)==400
    maxima=[0.]*4;conditioned=0
    for key,expected in ref.items():
        profile,index=map(int,key.split('_')[1:]);actual=native[key];assert len(expected)==len(actual)==4
        for j,(a,b) in enumerate(zip(expected,actual)):
            atol=[2e-6,.1,1e-12,1e-12][j];bound=atol+2e-7*max(abs(a),abs(b))
            if index==99 and j>=2:
                if abs(a-b)>bound:conditioned+=1
                bound+=float(endpoints[profile]['critical_volume_radius'])
                e=endpoints[profile];poly=[float(e[f'critical_polynomial_{c}']) for c in ['a','b','c']]
                for v in [a,b]:
                    z=(v+float(e['critical_shift']))*expected[1]/(8.31446261815324*expected[0])
                    assert abs(math.fsum([z**3,poly[0]*z*z,poly[1]*z,poly[2]]))<=float(e['critical_backward_bound'])
            assert math.isfinite(a) and math.isfinite(b) and abs(a-b)<=bound,(key,j,a,b,bound)
            maxima[j]=max(maxima[j],abs(a-b))
    samples={};checksums={}
    if not preflight_only:
        for sample in range(repeats+1):
            ref=parse(subprocess.check_output([str(exe)],env=env,text=True),('source_',))
            native=parse(subprocess.check_output(command,cwd=ROOT,text=True),('native_',));assert len(ref)==len(native)==4
            for profile in range(1,5):
                sk=f'source_{profile}';nk=f'native_{profile}';expected=ref[sk];actual=native[nk]
                assert len(expected)==len(actual)==6 and min(expected[0],actual[0])>0. and all(map(math.isfinite,expected+actual))
                assert expected[1]==actual[1]==2000.
                for j,(a,b) in enumerate(zip(expected[2:],actual[2:])):
                    bound=2000.*[2e-6,.1,1e-12,1e-12][j]+2e-7*max(abs(a),abs(b))
                    if j>=2:bound+=20.*float(endpoints[profile]['critical_volume_radius'])
                    assert abs(a-b)<=bound,(profile,j,a,b,bound)
                if sample:samples.setdefault(sk,[]).append(expected[0]);samples.setdefault(nk,[]).append(actual[0])
            checksums=dict(source=ref,native=native);print(f'Pure curve benchmark {sample}/{repeats} (0 is warm-up)',flush=True)
    assert rust==fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    assert harness==fingerprints(harness_paths) and cargo==fingerprints(cargo_paths)
    if preflight_only:print(json.dumps(dict(source_states=400,native_states=400,conditioned_endpoint_volumes=conditioned,maximum_absolute_differences=maxima),indent=2));return
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),processor=platform.processor(),git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),reference=json.loads((directory/'metadata.json').read_text()),compiler=subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],compiler_flags=flags,
        native_sources_sha256=rust,benchmark_sources_sha256=harness,cargo_inputs_sha256=cargo,
        source_preflight_states=400,native_preflight_states=400,conditioned_endpoint_volumes=conditioned,endpoint_conditioning=endpoints,maximum_preflight_absolute_differences=maxima,
        iterations_per_sample=20,samples_ns=samples,median_ns={k:statistics.median(v) for k,v in samples.items()},checksums=checksums,cargo_profile=dict(optimization='release',lto='thin',codegen_units=1),
        workload='Complete 100-point pure saturation grids: PR methane, pressure start 0.05 Pc, linear; SRK CO2, temperature start 0.65 Tc, logarithmic; translated PR water, pressure start 0.05 Pc, linear; QuantumCubic H2, pressure start 0.05 Pc, logarithmic. Full source/native grids and both phase volumes compared before timing. Final volume comparisons use source-derived cubic backward-error bounds; all interior outputs and all T/P retain standard tolerances. Twenty complete curves/sample, five samples after one warm-up. Includes cold initial saturation, warm analytic continuation, database-critical endpoint and both returned phase volumes; excludes construction and Python/FFI. Native additionally calculates Clapeyron slopes and explicit errors. Timed point counts and separate T/P/volume sums are checked, with the same endpoint conditioning allowance. Scoped workload measurement, not whole-library speed or exact raw critical-iteration parity.')
    (ROOT/output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['median_ns'],indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler',required=True);p.add_argument('--repeats',type=int,default=5);p.add_argument('--output',default='docs/rust/pure-curve-benchmark.json');p.add_argument('--preflight-only',action='store_true')
    a=p.parse_args();benchmark(a.compiler,a.repeats,a.output,a.preflight_only)
