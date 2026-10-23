"""Full pure MEOS coexistence calls, using only explicitly qualified source states."""
import argparse,csv,hashlib,json,math,os,platform,statistics,subprocess
from datetime import datetime,timezone
from pathlib import Path
from reference import ROOT
from benchmark_isochoric import fingerprints
from benchmark_solid_water import parse

PROFILES=['NIST_MEOS:C3','NIST_MEOS:R134A','LJ_MEOS:KR','LJ_MEOS:AR','LJTS_MEOS:KR']
FLOORS=[2e-7,1e-2,1e-12,1e-12,1e-4]

def benchmark(compiler,repeats,preflight_only,output):
    if repeats<1:raise ValueError('Positive repeat count required')
    compiler=Path(compiler).resolve();env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    fixture=ROOT/'tests/fixtures/meos_saturation.csv';meta=json.loads(fixture.with_suffix('.metadata.json').read_text())
    assert hashlib.sha256(fixture.read_bytes()).hexdigest()==meta['fixture_sha256']
    reference=json.loads((ROOT/'.reference/metadata.json').read_text());assert reference['library_sha256']==meta['library_sha256']
    rows=list(csv.DictReader(fixture.open(),delimiter=';'));inputs=[];expected={};counts={}
    for profile,config in enumerate(PROFILES,1):
        for mode,tag in enumerate(['SAT_P','SAT_T'],1):
            selected=[r for r in rows if r['case'].split('|')[1]==config and r['case'].split('|')[5]==tag]
            assert selected;counts[f'{profile}_{tag}']=len(selected)
            for i,row in enumerate(selected):
                inputs.append([profile,mode,float(row['case'].split('|')[2])])
                expected[f'check_{profile}_{tag}_{i}']=[float(row[k]) for k in ['temperature','pressure','liquid_volume','vapor_volume','dpdt']]
    path=ROOT/'.reference/meos-saturation-benchmark-input.txt'
    path.write_text(str(len(inputs))+'\n'+'\n'.join(' '.join(map(str,row)) for row in inputs)+'\n')
    source=ROOT/'tools/rust_port/benchmark_meos_saturation_reference.f90';exe=ROOT/'.reference'/('benchmark_meos_saturation.exe' if os.name=='nt' else 'benchmark_meos_saturation')
    native=fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    harness_paths=[Path(__file__).resolve(),source,ROOT/'tools/rust_port/benchmark_isochoric.py',ROOT/'tools/rust_port/benchmark_solid_water.py',fixture,fixture.with_suffix('.metadata.json')]
    harness=fingerprints(harness_paths);cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock'];cargo=fingerprints(cargo_paths)
    flags=['-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp']
    subprocess.run([str(compiler),*flags,'-I',str(ROOT/'.reference'),str(source),'@'+str(ROOT/'.reference/objects.rsp'),'-o',str(exe)],cwd=ROOT/'.reference',env=env,check=True)
    command=['cargo','bench','--bench','meos_saturation','-j','2'];subprocess.run([*command,'--no-run'],cwd=ROOT,check=True)
    ref=parse(subprocess.check_output([str(exe),'check',str(path)],env=env,text=True),('check_',))
    rust=parse(subprocess.check_output([*command,'--','--check'],cwd=ROOT,text=True),('check_',))
    assert ref.keys()==rust.keys()==expected.keys()
    maximum=[0.]*5
    def compare(a,b,scale=1.):
        assert len(a)==len(b)==5
        for j,(x,y) in enumerate(zip(a,b)):
            assert math.isfinite(x) and math.isfinite(y) and abs(x-y)<=scale*FLOORS[j]+2e-6*max(abs(x),abs(y)),(j,x,y)
    for key,a in ref.items():
        compare(expected[key],a);compare(a,rust[key])
        maximum=[max(m,abs(x-y)) for m,x,y in zip(maximum,a,rust[key])]
    samples={};checksums={}
    if not preflight_only:
        for sample in range(repeats+1):
            ref=parse(subprocess.check_output([str(exe),'time',str(path)],env=env,text=True),('source_',))
            rust=parse(subprocess.check_output(command,cwd=ROOT,text=True),('native_',))
            assert len(ref)==len(rust)==10
            for suffix in counts:
                a=ref['source_'+suffix];b=rust['native_'+suffix];assert min(a[0],b[0])>0.
                compare(a[1:],b[1:],1000.)
                if sample:
                    samples.setdefault('source_'+suffix,[]).append(a[0]);samples.setdefault('native_'+suffix,[]).append(b[0])
            checksums=dict(source=ref,native=rust);print('MEOS saturation benchmark',sample,'/',repeats,flush=True)
    assert native==fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    assert harness==fingerprints(harness_paths) and cargo==fingerprints(cargo_paths)
    if preflight_only:print(json.dumps(dict(states=len(inputs),points_per_profile=counts,maximum_absolute_differences=maximum),indent=2));return
    medians={k:statistics.median(v) for k,v in samples.items()}
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),processor=platform.processor(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),
        reference=reference,compiler=subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],compiler_flags=flags,
        native_sources_sha256=native,benchmark_sources_sha256=harness,cargo_inputs_sha256=cargo,
        numeric_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),numeric_input_rows=inputs,
        profiles=PROFILES,source_preflight_states=len(inputs),native_preflight_states=len(inputs),points_per_profile=counts,maximum_preflight_absolute_differences=maximum,
        iterations_per_sample=1000,samples_ns=samples,median_ns=medians,checksums=checksums,
        workload='Five pure MEOS models, pressure at fixed T and cold temperature inversion at fixed P timed separately. Every qualified original-source saturation fixture state for these profiles is preflighted; failed and coincident-phase source results are excluded explicitly and retained in the qualification ledger. Each sample cycles 1000 calls and compares five separate output sums: T, P, both phase volumes and Clapeyron slope. Source safe_bubP/safe_bubT plus specificVolume and enthalpy calls versus native saturation_pressure/saturation_temperature. Includes temperature preparation and native fallback when needed; excludes model construction, parsing, Python and FFI. No warm previous-state guess between calls. Five measured samples after one warm-up. Reduced LJTS, other models and source-failed recovery states are outside this timing workload; no whole-library speed claim.')
    Path(output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(median_ns=medians,speedups={s:medians['source_'+s]/medians['native_'+s] for s in counts}),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler',required=True);p.add_argument('--repeats',type=int,default=5);p.add_argument('--preflight-only',action='store_true');p.add_argument('--output',default=ROOT/'docs/rust/meos-saturation-benchmark.json');a=p.parse_args();benchmark(a.compiler,a.repeats,a.preflight_only,a.output)
