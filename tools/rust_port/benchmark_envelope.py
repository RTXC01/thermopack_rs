"""Whole-curve native-executable timing with identical adaptive sample grids."""
import argparse,json,math,os,platform,statistics,subprocess
from pathlib import Path
from datetime import datetime,timezone
from reference import ROOT
from benchmark_isochoric import fingerprints
from benchmark_solid_water import parse

def benchmark(compiler,repeats,output,preflight_only,quality=False):
    if repeats<1:raise ValueError('Positive sample count required')
    compiler=Path(compiler).resolve();env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    directory=ROOT/'.reference';source=ROOT/'tools/rust_port/benchmark_envelope_reference.f90'
    exe=directory/('benchmark_envelope.exe' if os.name=='nt' else 'benchmark_envelope')
    flags=['-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp']
    rust=fingerprints(sorted((ROOT/'rust').rglob('*.rs')))
    harness_paths=[Path(__file__).resolve(),source,ROOT/'tools/rust_port/benchmark_isochoric.py',ROOT/'tools/rust_port/benchmark_solid_water.py']
    harness=fingerprints(harness_paths);cargo_paths=[ROOT/'Cargo.toml',ROOT/'Cargo.lock'];cargo=fingerprints(cargo_paths)
    subprocess.run([str(compiler),*flags,'-I',str(directory),str(source),'@'+str(directory/'objects.rsp'),'-o',str(exe)],cwd=directory,env=env,check=True)
    command=['cargo','bench','--bench','quality_envelope' if quality else 'envelope','-j','2']
    subprocess.run([*command,'--no-run'],cwd=ROOT,check=True)
    ref_command=lambda check:[str(exe),'check' if check else 'time','quality' if quality else '']
    ref=parse(subprocess.check_output(ref_command(True),env=env,text=True),('check_',))
    native=parse(subprocess.check_output([*command,'--','--check'],cwd=ROOT,text=True),('check_',))
    assert ref.keys()==native.keys() and (quality or len(ref)==224)
    total=len(ref)
    maxima=[0.]*5;counts={}
    for key,expected in ref.items():
        actual=native[key];assert len(expected)==len(actual)
        profile=key.split('_')[1];counts[profile]=counts.get(profile,0)+1
        for j,(a,b) in enumerate(zip(expected,actual)):
            atol=2e-5 if j==0 else 2. if j==1 else 0. if j==2 else 1e-10 if j<5 else 2e-8
            assert math.isfinite(a) and math.isfinite(b) and abs(a-b)<=atol+2e-6*max(abs(a),abs(b)),(key,j,a,b)
            if j<5:maxima[j]=max(maxima[j],abs(a-b))
    samples={};checksums={}
    if not preflight_only:
        for sample in range(repeats+1):
            ref=parse(subprocess.check_output(ref_command(False),env=env,text=True),('source_',))
            native=parse(subprocess.check_output(command,cwd=ROOT,text=True),('native_',))
            assert len(ref)==len(native)==4
            for profile in range(1,5):
                sk=f'source_{profile}';nk=f'native_{profile}';expected=ref[sk];actual=native[nk]
                assert len(expected)==len(actual) and min(expected[0],actual[0])>0 and all(map(math.isfinite,expected+actual))
                states=10*counts[str(profile)]
                for j,(a,b) in enumerate(zip(expected[1:],actual[1:])):
                    if j in [0,5]:assert a==b,(profile,j,a,b);continue
                    atol=2e-5 if j==1 else 2. if j==2 else 1e-10 if j<5 else 2e-8
                    assert abs(a-b)<=states*atol+2e-6*max(abs(a),abs(b)),(profile,j,a,b)
                if sample:samples.setdefault(sk,[]).append(expected[0]);samples.setdefault(nk,[]).append(actual[0])
            checksums=dict(source=ref,native=native)
            print(f'Envelope benchmark {sample}/{repeats} (0 is warm-up)',flush=True)
    assert rust==fingerprints(sorted((ROOT/'rust').rglob('*.rs'))),'Native sources changed during measurement'
    assert harness==fingerprints(harness_paths),'Benchmark inputs changed during measurement'
    assert cargo==fingerprints(cargo_paths),'Cargo inputs changed during measurement'
    if preflight_only:
        print(json.dumps(dict(source_states=total,native_states=total,points_per_curve=counts,maximum_absolute_differences=maxima),indent=2));return
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),processor=platform.processor(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),reference=json.loads((directory/'metadata.json').read_text()),
        compiler=subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],compiler_flags=flags,
        native_sources_sha256=rust,benchmark_sources_sha256=harness,cargo_inputs_sha256=cargo,
        source_preflight_states=total,native_preflight_states=total,points_per_curve=counts,maximum_preflight_absolute_differences=maxima,
        iterations_per_sample=10,samples_ns=samples,median_ns={k:statistics.median(v) for k,v in samples.items()},checksums=checksums,
        cargo_profile=dict(optimization='release',lto='thin',codegen_units=1),
        workload='Complete mixture VLE envelopes from dew pressure 0.1 MPa back to bubble pressure 0.1 MPa; maximum pressure 15 MPa, minimum temperature 80 K, maximum sensitive log-variable step 0.25. Profiles: PR and SRK equimolar methane/ethane; PR methane/ethane/propane [0.2,0.3,0.5]; PR-UNIFAC equimolar methane/propane. Source and baseline native adaptive grids have identical counts (48/48/59/69) and every preflight T/P/beta/liquid-volume/vapor-volume/liquid-composition/vapor-composition is compared. Ten complete curves/sample, five samples after one warm-up. Timed checks compare point count, phase-flag sum, separate T/P/volume sums and each liquid/vapor composition sum. Model construction and Python/FFI excluded; cold saturation initialization, full continuation/phase formulation switch, all returned phase coordinates and both phase volumes included. Fortran post-evaluates both volumes at its own returned states; native states already carry them. Native additionally computes the saturation slope. Tangents reuse the converged analytic Jacobian and are cached by exact coordinates, specification and formulation. Scoped workloads, not whole-library speed or proof of other adaptive-grid/failure trajectories.')
    if quality:
        report['workload']='Complete fixed-quality loci, initial molar vapor fraction 0.25, through source phase formulation complement to 0.75. Profiles: PR and SRK equimolar methane/ethane, PR methane/ethane/propane [0.2,0.3,0.5], PR UNIFAC equimolar methane/propane. Start/end pressure 0.1 MPa, cap 15 MPa, minimum temperature 80 K, maximum sensitive log-coordinate step 0.25. Entire source/native sampled grids have identical counts/order and every T/P/beta/phase-volume/composition is compared before timing. Ten complete curves/sample, five samples after one warm-up; timed counts, beta sum and separate T/P/volume/composition sums all checked. Includes cold initialization, complete continuation and both reported phase volumes; native also reports fixed-quality slopes/residuals/feed. Excludes model construction and Python/FFI. Scoped whole-curve performance, not complete raw failure fidelity or whole-library speed.'
    output=output or ('docs/rust/quality-envelope-benchmark.json' if quality else 'docs/rust/envelope-benchmark.json')
    (ROOT/output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['median_ns'],indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler',required=True);p.add_argument('--repeats',type=int,default=5)
    p.add_argument('--output');p.add_argument('--preflight-only',action='store_true');p.add_argument('--quality',action='store_true')
    a=p.parse_args();benchmark(a.compiler,a.repeats,a.output,a.preflight_only,a.quality)
