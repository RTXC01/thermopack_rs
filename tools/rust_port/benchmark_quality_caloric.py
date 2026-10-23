"""Source, compatible Rust and physical Rust H/S corrections from identical seeds."""
import argparse,csv,hashlib,json,math,os,platform,statistics,subprocess
from datetime import datetime,timezone
from pathlib import Path
from reference import ROOT
from benchmark_isochoric import fingerprints
from benchmark_solid_water import parse

def benchmark(compiler,repeats,preflight_only):
    if repeats<1:raise ValueError('Positive repeat count required')
    compiler=Path(compiler).resolve();env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    fixture=ROOT/'tests/fixtures/quality_caloric.csv';meta=json.loads(fixture.with_suffix('.metadata.json').read_text())
    assert hashlib.sha256(fixture.read_bytes()).hexdigest()==meta['fixture_sha256']
    reference=json.loads((ROOT/'.reference/metadata.json').read_text());assert reference['library_sha256']==meta['library_sha256']
    profiles=[('PR','C1,C2','classic'),('SRK','C1,C2','classic'),('PR','C1,C3','translated'),('PR','C1,C3','unifac')]
    inputs=[];counts={}
    rows=list(csv.DictReader(fixture.open(),delimiter=';'))
    for profile,signature in enumerate(profiles,1):
        selected=[]
        for row in rows:
            f=row['case'].split('|');o=f[6].split(':');beta=float(f[5])
            if tuple(f[1:4])!=signature or not 0.<beta<1.:continue
            g=list(map(float,f[7].split(',')));s=6 if o[0]=='H' else 5;spec=float(o[1])
            selected.append([profile,g[0],g[1],beta,s,spec,*g[2:]])
        assert selected;counts[profile]=len(selected);inputs.extend(selected)
    path=ROOT/'.reference/quality-caloric-benchmark-input.txt'
    path.write_text(str(len(inputs))+'\n'+'\n'.join(' '.join(map(repr,row)) for row in inputs)+'\n')
    source=ROOT/'tools/rust_port/benchmark_quality_reference.f90';exe=ROOT/'.reference'/('benchmark_quality_caloric.exe' if os.name=='nt' else 'benchmark_quality_caloric')
    native=fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    harness_paths=[Path(__file__).resolve(),source,ROOT/'tools/rust_port/benchmark_isochoric.py',ROOT/'tools/rust_port/benchmark_solid_water.py',fixture,fixture.with_suffix('.metadata.json')]
    harness=fingerprints(harness_paths);cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock'];cargo=fingerprints(cargo_paths)
    flags=['-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp']
    subprocess.run([str(compiler),*flags,'-I',str(ROOT/'.reference'),str(source),'@'+str(ROOT/'.reference/objects.rsp'),'-o',str(exe)],cwd=ROOT/'.reference',env=env,check=True)
    command=['cargo','bench','--bench','quality_caloric','-j','2'];subprocess.run([*command,'--no-run'],cwd=ROOT,check=True)
    ref=parse(subprocess.check_output([str(exe),'check',str(path),'caloric'],env=env,text=True),('check_',))
    rust=parse(subprocess.check_output([*command,'--','--check'],cwd=ROOT,text=True),('check_',));assert len(rust)==2*len(ref)==2*len(inputs)
    maxima={mode:[0.]*12 for mode in ['baseline','physical']}
    def compare(a,b,mode,scale=1.):
        assert len(a)==len(b)==12
        for j,(x,y) in enumerate(zip(a,b)):
            atol=2e-5 if j==0 else 2. if j==1 else 0. if j==2 else 1e-10 if j<5 else 1e-3 if j==9 else 1e-5 if j==10 else 1e-8 if j==11 else 2e-8
            rtol=2e-6 if mode=='physical' or j==9 else 2e-7
            assert math.isfinite(x) and math.isfinite(y) and abs(x-y)<=scale*atol+rtol*max(abs(x),abs(y)),(j,x,y)
    for key,a in ref.items():
        for mode in maxima:
            b=rust[key.replace('check_',f'check_{mode}_')];compare(a,b,mode)
            maxima[mode]=[max(m,abs(x-y)) for m,x,y in zip(maxima[mode],a,b)]
    samples={};checksums={}
    if not preflight_only:
        for sample in range(repeats+1):
            ref=parse(subprocess.check_output([str(exe),'time',str(path),'caloric'],env=env,text=True),('source_',))
            rust=parse(subprocess.check_output(command,cwd=ROOT,text=True),('baseline_','physical_'));assert len(ref)==4 and len(rust)==8
            for profile in range(1,5):
                a=ref[f'source_{profile}'];assert a[0]>0.
                if sample:samples.setdefault(f'source_{profile}',[]).append(a[0])
                for mode in maxima:
                    b=rust[f'{mode}_{profile}'];assert b[0]>0.;compare(a[1:],b[1:],mode,1000.)
                    if sample:samples.setdefault(f'{mode}_{profile}',[]).append(b[0])
            checksums=dict(source=ref,native=rust);print('Quality benchmark',sample,'/',repeats,flush=True)
    assert native==fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    assert harness==fingerprints(harness_paths) and cargo==fingerprints(cargo_paths)
    if preflight_only:print(json.dumps(dict(states=len(inputs),points_per_profile=counts,maximum_absolute_differences=maxima),indent=2));return
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),processor=platform.processor(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),
        reference=reference,compiler=subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],compiler_flags=flags,
        native_sources_sha256=native,benchmark_sources_sha256=harness,cargo_inputs_sha256=cargo,
        numeric_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),numeric_input_rows=inputs,
        source_preflight_states=len(inputs),native_preflight_states=len(inputs),points_per_profile=counts,maximum_preflight_absolute_differences=maxima,
        iterations_per_sample=1000,samples_ns=samples,median_ns={k:statistics.median(v) for k,v in samples.items()},checksums=checksums,
        workload='Fixed-quality H/S corrections from identical perturbed source coexistence seeds. Four binary workloads: PR methane/ethane, SRK methane/ethane, translated PR methane/propane, PR UNIFAC methane/propane. Every finite interior-quality caloric fixture state in those profiles is compared before timing, including shifted targets. Each sample cycles 1000 corrections, computes both phase volumes/compositions, fixed-quality dP/dT, total H and S, and compares all twelve separate output sums. Source sat_newton/newton_extrapolate plus H/S calls versus baseline_quality_state_with_specification and quality_state_with_specification plus quality_properties. Physical Rust uses the corrected analytic Jacobian; source-compatible Rust preserves the source Jacobian. Excludes model construction, seed generation, parsing and Python/FFI. Five measured samples after one warm-up. No whole-library or cold-initialization speed claim.')
    (ROOT/'docs/rust/quality-caloric-benchmark.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['median_ns'],indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler',required=True);p.add_argument('--repeats',type=int,default=5);p.add_argument('--preflight-only',action='store_true');a=p.parse_args();benchmark(a.compiler,a.repeats,a.preflight_only)
