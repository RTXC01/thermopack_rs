"""Native-loop sound benchmark; source-compatible outputs required before timing."""
import argparse,hashlib,json,math,os,platform,statistics,subprocess
from datetime import datetime,timezone
from pathlib import Path
from reference import ROOT


def benchmark(compiler,repeats,tcpr):
    if repeats<1:raise ValueError('At least one measured repetition required')
    compiler=Path(compiler).resolve();env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    out=ROOT/'.reference';exe=out/('benchmark_sound.exe' if os.name=='nt' else 'benchmark_sound')
    flags=['-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp']
    source=ROOT/'tools/rust_port/benchmark_sound_reference.f90'
    subprocess.run([str(compiler),*flags,'-I',str(out),str(source),'@'+str(out/'objects.rsp'),'-o',str(exe)],env=env,cwd=out,check=True)
    subprocess.run(['cargo','bench','--bench','sound','--no-run'],cwd=ROOT,check=True)
    def parse(output,names):
        rows={row.split(',')[0]:list(map(float,row.split(',')[1:])) for row in output.splitlines() if row.split(',')[0] in names}
        assert rows.keys()==names,(rows.keys(),names,output)
        assert all(len(v)==2 and v[0]>0. and all(map(math.isfinite,v)) for v in rows.values())
        return rows
    timings={};checksums={};native=['cargo','bench','--bench','sound']+(['--','--tcpr'] if tcpr else [])
    for i in range(repeats+1):
        reference=parse(subprocess.check_output([str(exe)]+(['tcPR'] if tcpr else []),env=env,text=True),{'tv','equilibrium'})
        rust=parse(subprocess.check_output(native,cwd=ROOT,text=True),{'tv_prepared','tv_with_mixing','equilibrium'})
        for name,values in rust.items():
            ref=reference['tv' if name.startswith('tv_') else name]
            assert abs(ref[1]-values[1])<=2e-7*max(abs(ref[1]),abs(values[1]),1.),(name,ref,values)
            if i:timings.setdefault('rust_'+name,[]).append(values[0])
        if i:
            for name,values in reference.items():timings.setdefault('fortran_'+name,[]).append(values[0])
            checksums={'fortran':{k:v[1] for k,v in reference.items()},'rust':{k:v[1] for k,v in rust.items()}}
        print(f'Sound benchmark {i}/{repeats} (0 is warm-up)',flush=True)
    report=dict(measured_at=datetime.now(timezone.utc).isoformat(),parameter_reference='tcPR' if tcpr else 'DEFAULT',
                platform=platform.platform(),processor=platform.processor(),rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),
                reference=json.loads((out/'metadata.json').read_text()),reference_benchmark_flags=flags,
                native_sources_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in sorted((ROOT/'rust').rglob('*.rs'))},
                benchmark_sources_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in [Path(__file__),source]},
                cargo_profile=dict(optimization='release',lto='thin',codegen_units=1),iterations_per_sample=dict(tv=500000,equilibrium=10000),
                samples_ns_per_state=timings,median_ns_per_state={k:statistics.median(v) for k,v in timings.items()},checksums=checksums,
                workload='PR; TV: CO2/C1/N2 0.7/0.2/0.1 at 300 K, 1000 molar volumes 1e-4..1.099e-3; equilibrium: C1/C2/C3 0.5/0.3/0.2 at 1 MPa, 100 precomputed TP splits from 180..259.2 K; flash initialization outside timer')
    (ROOT/'docs/rust'/('benchmark-sound-tcpr.json' if tcpr else 'benchmark-sound.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['median_ns_per_state'],indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler',required=True);parser.add_argument('--repeats',type=int,default=5);parser.add_argument('--tcpr',action='store_true')
    args=parser.parse_args();benchmark(args.compiler,args.repeats,args.tcpr)
