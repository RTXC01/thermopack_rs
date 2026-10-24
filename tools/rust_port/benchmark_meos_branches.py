"""Source-backed native-only complete bounded root sets and minimum-Gibbs selection."""
import argparse,csv,hashlib,json,math,platform,statistics,subprocess
from datetime import datetime,timezone
from pathlib import Path
from reference import ROOT
from benchmark_isochoric import fingerprints
from benchmark_solid_water import parse

PROFILES=[('NIST_MEOS:C3',.95),('NIST_MEOS:R134A',1.05),('NIST_MEOS:N-H2',.95),('LJ_MEOS:AR',.95),('LJTS_MEOS:KR',.95)]

def benchmark(repeats,preflight_only):
    if repeats<1:raise ValueError('Positive sample count required')
    fixture=ROOT/'tests/fixtures/meos_pressure_branches.csv';meta=json.loads(fixture.with_suffix('.metadata.json').read_text());assert hashlib.sha256(fixture.read_bytes()).hexdigest()==meta['fixture_sha256']
    metadata=list(csv.DictReader((ROOT/'tests/fixtures/meos_metadata.csv').open(),delimiter=';'))
    # The kernel metadata CSV keeps source constructor Tc in its first numeric field.
    tc={r['case'].split('|')[1]:float(r['Tc']) for r in metadata}
    rows=list(csv.DictReader(fixture.open(),delimiter=';'));expected={};inputs={};counts={}
    for profile,(config,tr) in enumerate(PROFILES,1):
        selected=[r for r in rows if r['case'].split('|')[1]==config and abs(float(r['case'].split('|')[2])/tc[config]-tr)<1e-12];assert selected
        inputs[profile]=[];counts[profile]=len(selected)
        for i,r in enumerate(selected):
            roots=list(map(float,r['stable_densities'].split(',')));values=[float(len(roots)),sum(roots),roots[0],roots[-1],1./float(r['minimum_gibbs_volume'])]
            expected[f'check_{profile}_{i}']=values;inputs[profile].append(dict(case=r['case'],expected=values))
    native=fingerprints(sorted((ROOT/'rust').rglob('*.rs')));cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock'];cargo=fingerprints(cargo_paths)
    harness_paths=[Path(__file__).resolve(),ROOT/'tools/rust_port/benchmark_isochoric.py',ROOT/'tools/rust_port/benchmark_solid_water.py',fixture,fixture.with_suffix('.metadata.json'),ROOT/'tests/fixtures/meos_metadata.csv']
    harness=fingerprints(harness_paths)
    command=['cargo','bench','--bench','meos_branches','-j','2'];subprocess.run([*command,'--no-run'],cwd=ROOT,check=True)
    actual=parse(subprocess.check_output([*command,'--','--check'],cwd=ROOT,text=True),('check_',));assert actual.keys()==expected.keys()
    def compare(a,b,scale=1.):
        assert len(a)==len(b)==5
        for j,(x,y) in enumerate(zip(a,b)):
            tolerance=0. if j==0 else 1e-7*scale+3e-7*max(abs(x),abs(y))
            assert math.isfinite(y) and abs(x-y)<=tolerance,(j,x,y,tolerance)
    maxima=[0.]*5
    for key,a in expected.items():compare(a,actual[key]);maxima=[max(m,abs(x-y)) for m,x,y in zip(maxima,a,actual[key])]
    samples={};checksums={}
    if not preflight_only:
        for sample in range(repeats+1):
            values=parse(subprocess.check_output(command,cwd=ROOT,text=True),('setup_','fresh_','prepared_'));assert len(values)==15
            for profile in inputs:
                setup=values[f'setup_{profile}'][0];assert setup>0.
                if sample:samples.setdefault(f'setup_{profile}',[]).append(setup)
                for mode,iterations in [('fresh',100),('prepared',1000)]:
                    key=f'{mode}_{profile}';v=values[key];assert v[0]>0.
                    wanted=[sum(inputs[profile][i%len(inputs[profile])]['expected'][j] for i in range(iterations)) for j in range(5)]
                    compare(wanted,v[1:],iterations)
                    if sample:samples.setdefault(key,[]).append(v[0])
            checksums=values;print('MEOS branch benchmark',sample,'/',repeats,flush=True)
    assert native==fingerprints(sorted((ROOT/'rust').rglob('*.rs'))) and cargo==fingerprints(cargo_paths) and harness==fingerprints(harness_paths)
    if preflight_only:print(json.dumps(dict(cases=len(expected),points_per_profile=counts,maximum_absolute_differences=maxima),indent=2));return
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),processor=platform.processor(),git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),source_library_sha256=meta['library_sha256'],source_probe=meta['probe'],
        native_sources_sha256=native,cargo_inputs_sha256=cargo,benchmark_sources_sha256=harness,profiles=PROFILES,numeric_cases=inputs,preflight_cases=len(expected),maximum_absolute_differences=maxima,
        iterations_per_sample=dict(fresh=100,prepared=1000,setup=10),samples_ns=samples,median_ns={k:statistics.median(v) for k,v in samples.items()},checksums=checksums,
        workload='Complete stable bounded root set plus minimum-Gibbs selection at fixed T/P, with caller-owned density buffers. Five profiles; every resolved source-backed fixture pressure at each chosen isotherm is preflighted. Outputs are root count, sum of all root densities, extreme densities and minimum-Gibbs density; five separate sums checked per sample. Fresh mode includes interval discovery and destruction; prepared mode reuses branches; setup is timed separately. Model construction and input parsing are excluded. Native-only measurements: no original bulk complete-root API is available here, and legacy phase-switch iteration has different selection semantics. No source speedup claim. Default five measured samples after a warm-up.')
    (ROOT/'docs/rust/meos-branches-benchmark.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['median_ns'],indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repeats',type=int,default=5);p.add_argument('--preflight-only',action='store_true');a=p.parse_args();benchmark(a.repeats,a.preflight_only)
