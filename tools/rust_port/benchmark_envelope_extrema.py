"""Native-only extrema timing; every output is checked against source-backed roots.

No speedup versus the defective legacy optional extrema path is claimed.
"""
import argparse,csv,hashlib,json,math,platform,statistics,subprocess
from datetime import datetime,timezone
from pathlib import Path
from reference import ROOT
from benchmark_isochoric import fingerprints
from benchmark_solid_water import parse

def benchmark(repeats,preflight_only):
    if repeats<1:raise ValueError('Positive repeat count required')
    fixture=ROOT/'tests/fixtures/envelope_extrema.csv';meta=json.loads(fixture.with_suffix('.metadata.json').read_text())
    assert hashlib.sha256(fixture.read_bytes()).hexdigest()==meta['fixture_sha256']
    profiles=[('PR','C1,C2','classic'),('SRK','C1,C2','classic'),('PR','C1,C3','translated'),('PR','C1,C3','unifac')]
    expected={};inputs={}
    for profile,signature in enumerate(profiles,1):
        selected=[]
        for row in csv.DictReader(fixture.open(),delimiter=';'):
            f=row['case'].split('|')
            if tuple(f[1:4])!=signature or f[8]!='100000.0,15000000.0,80.0,0.0,1.0,0.25,1000':continue
            key=f'check_{profile}_{len(selected)}';selected.append(row['case']);expected[key]=list(map(float,row['source_backed_reference_state'].split(',')))+[0.]
        assert len(selected)==2;inputs[profile]=selected
    native=fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    harness_paths=[Path(__file__).resolve(),ROOT/'tools/rust_port/benchmark_isochoric.py',ROOT/'tools/rust_port/benchmark_solid_water.py',fixture,fixture.with_suffix('.metadata.json')]
    harness=fingerprints(harness_paths);cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock'];cargo=fingerprints(cargo_paths)
    command=['cargo','bench','--bench','envelope_extrema','-j','2'];subprocess.run([*command,'--no-run'],cwd=ROOT,check=True)
    actual=parse(subprocess.check_output([*command,'--','--check'],cwd=ROOT,text=True),('check_',));assert actual.keys()==expected.keys()
    def compare(a,b,scale=1.):
        assert len(a)==len(b)==10
        for i,(x,y) in enumerate(zip(a,b)):
            atol=2e-5 if i==0 else 2. if i==1 else 0. if i==2 else 1e-10 if i<5 else 1e-9 if i==9 else 2e-8
            assert math.isfinite(y) and abs(x-y)<=atol*scale+2e-6*max(abs(x),abs(y)),(i,x,y)
    maxima=[0.]*10
    for key,a in expected.items():
        b=actual[key];compare(a,b);maxima=[max(m,abs(x-y)) for m,x,y in zip(maxima,a,b)]
    samples={};checksums={}
    if not preflight_only:
        for sample in range(repeats+1):
            timed=parse(subprocess.check_output(command,cwd=ROOT,text=True),('native_',));assert len(timed)==4
            for profile in range(1,5):
                name=f'native_{profile}';values=timed[name];assert values[0]>0.
                source_sum=[sum(expected[f'check_{profile}_{i%2}'][j] for i in range(1000)) for j in range(10)]
                compare(source_sum,values[1:],1000.)
                if sample:samples.setdefault(name,[]).append(values[0])
            checksums=timed;print('Native extrema benchmark',sample,'/',repeats,flush=True)
    assert native==fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    assert harness==fingerprints(harness_paths) and cargo==fingerprints(cargo_paths)
    if preflight_only:print(json.dumps(dict(states=len(expected),maximum_absolute_differences=maxima),indent=2));return
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),processor=platform.processor(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),
        reference={k:meta[k] for k in ['baseline','compiler','library_sha256','oracle_sha256']},
        native_sources_sha256=native,benchmark_sources_sha256=harness,cargo_inputs_sha256=cargo,numeric_input_cases=inputs,
        source_backed_reference_preflight_states=len(expected),maximum_preflight_absolute_differences=maxima,
        compared_outputs=['temperature_K','pressure_Pa','source_phase_fraction','liquid_volume_m3_per_mol','vapor_volume_m3_per_mol','liquid_x1','liquid_x2','vapor_y1','vapor_y2','dimensionless_unscaled_stationarity'],
        iterations_per_sample=1000,samples_ns=samples,median_ns={k:statistics.median(v) for k,v in samples.items()},checksums=checksums,
        workload='Native-only warm local P/T extrema refinement, two source-backed guesses per binary profile. Ten outputs, including both phase volumes/compositions and stationarity, are checked before timing and as separate sums in every sample. The full returned state includes dP/dT, which is computed but is singular at a temperature maximum and is not used as a parity/checksum quantity. Five measured samples after one warm-up. Excludes construction, parsing, complete tracing and source-reference generation. No speedup comparison against the defective original optional criconden path.')
    (ROOT/'docs/rust/envelope-extrema-benchmark.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['median_ns'],indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repeats',type=int,default=5);p.add_argument('--preflight-only',action='store_true');a=p.parse_args();benchmark(a.repeats,a.preflight_only)
