"""Compare all stored MBWR fits against original coefficient/Helmholtz routines."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe
from generate_mbwr_data import records

# Density alpha (six), extensive F (ten), TV P (four), density P (six).
ABSOLUTE=[2e-8,2e-9,1e-10,2e-10,1e-11,1e-13,
          1e-7,1e-8,1e-7,1e-7,1e-9,1e-8,1e-8,1e-4,1e-7,1e-7,
          .001,.001,.1,.001,.001,.001,1e-5,1e-5,1e-7,1e-9]
RELATIVE=2e-8
QUANTITIES=['alpha','alpha_T','alpha_rho','alpha_TT','alpha_Trho','alpha_rhorho',
    'F','F_T','F_V','F_n','F_TT','F_TV','F_Tn','F_VV','F_Vn','F_nn',
    'P_TV','P_T_at_V','P_V','P_n','P_density','P_T_at_rho','P_rho','P_TT','P_Trho','P_rhorho']


def cases():
    for record in records():
        base=dict(record_index=int(record['source_name'].removeprefix('mbwrcx')),source_name=record['source_name'])
        yield dict(base,mode='MBWR_PARAMETERS',temperature=0.,volume=0.,moles=0.)
        yield dict(base,mode='MBWR_INTEGRATED',temperature=0.,volume=0.,moles=0.)
        for t in [.5,1.,5.,10.,20.,50.,100.,200.,300.,500.,1000.,2000.]:
            for rho in [1e-5,.01,1.,100.,1000.,5000.,10000.,20000.]:
                yield dict(base,mode='MBWR_BASELINE_CORE',temperature=t,volume=1./rho,moles=1.)
        for t in [100.,500.]:
            for rho in [100.,10000.]:
                for n in [.3,3.5]:
                    yield dict(base,mode='MBWR_BASELINE_CORE',temperature=t,volume=n/rho,moles=n)


def encode(request):
    return '|'.join(['MBWR',request['source_name'],str(request['temperature']),str(request['volume']),str(request['moles']),request['mode']])


def validate(directory):
    probe=verify_probe();results=[];errors=[]
    with ReferenceProcess(directory) as oracle:
        for request in cases():
            result=oracle.call(request)
            if 'calculation_error' in result or 'values' not in result:
                errors.append(dict(request=request,**result));continue
            results.append((request,encode(request),result['values']))
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    lines=subprocess.run([str(exe)],input='\n'.join(c for _,c,_ in results)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(lines)==len(results);failures=[];maximum_absolute=[0.]*26;maximum_relative=[0.]*26
    for (request,case,expected),line in zip(results,lines):
        if line.startswith('ERROR'):
            failures.append(dict(case=case,error=line,reference=expected));continue
        actual=list(map(float,line.split(',')));assert len(actual)==len(expected)
        if request['mode'] in ('MBWR_PARAMETERS','MBWR_INTEGRATED'):bad=[i for i,(a,b) in enumerate(zip(expected,actual)) if a!=b]
        else:
            bad=[i for i,(a,b) in enumerate(zip(expected,actual)) if not math.isfinite(b) or abs(a-b)>ABSOLUTE[i]+RELATIVE*max(abs(a),abs(b))]
            for i,(a,b) in enumerate(zip(expected,actual)):
                maximum_absolute[i]=max(maximum_absolute[i],abs(a-b))
                maximum_relative[i]=max(maximum_relative[i],abs(a-b)/max(abs(a),abs(b),1e-300))
        if bad:failures.append(dict(case=case,quantities=bad,reference=expected,rust=actual))
    stable_records=[(r,c.replace('MBWR_BASELINE_CORE','MBWR_CORE'),v) for r,c,v in results if r['mode']=='MBWR_BASELINE_CORE']
    stable_lines=subprocess.run([str(exe)],input='\n'.join(c for _,c,_ in stable_records)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(stable_lines)==len(stable_records);stable_differences=[]
    for (_,case,expected),line in zip(stable_records,stable_lines):
        if line.startswith('ERROR'):
            stable_differences.append(dict(case=case,error=line));continue
        actual=list(map(float,line.split(',')))
        bad=[i for i,(a,b) in enumerate(zip(expected,actual)) if not math.isfinite(b) or abs(a-b)>ABSOLUTE[i]+RELATIVE*max(abs(a),abs(b))]
        if bad:stable_differences.append(dict(case=case,quantities=bad,reference=expected,stable_rust=actual))
    (ROOT/'tests/fixtures/mbwr_stable_differences.json').write_text(json.dumps(dict(
        generator='tools/rust_port/validate_mbwr.py',attempted_cases=len(stable_records),matched_cases=len(stable_records)-len(stable_differences),
        differences=len(stable_differences),failures=stable_differences,
        convention='stable analytic integral versus original integrated-term summation; differences retained without changing parity tolerances; use baseline API for original rounding behavior'),indent=2)+'\n')
    for mode,name in [('MBWR_PARAMETERS','mbwr_parameters'),('MBWR_INTEGRATED','mbwr_integrated'),('MBWR_BASELINE_CORE','mbwr_core')]:
        group=[(r,c,v) for r,c,v in results if r['mode']==mode];path=ROOT/'tests/fixtures'/(name+'.csv')
        with path.open('w',newline='') as handle:
            writer=csv.writer(handle,delimiter=';',lineterminator='\n');writer.writerow(['case','values'])
            writer.writerows([[c,*v] for _,c,v in group])
        metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
        metadata.update(cases=len(group),attempted_cases=len(group)+sum(e['request']['mode']==mode for e in errors),
            reference_errors=[e for e in errors if e['request']['mode']==mode],probe=probe,generator='tools/rust_port/validate_mbwr.py',
            fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),relative_tolerance=RELATIVE if mode=='MBWR_BASELINE_CORE' else 0.,
            absolute_tolerances=ABSOLUTE if mode=='MBWR_BASELINE_CORE' else [0.],
            quantities=QUANTITIES if mode=='MBWR_BASELINE_CORE' else ('kind, raw coefficients, gamma' if mode=='MBWR_PARAMETERS' else 'original integrated alpha coefficients in source order'),
            maximum_absolute_differences=maximum_absolute if mode=='MBWR_BASELINE_CORE' else None,
            maximum_relative_differences=maximum_relative if mode=='MBWR_BASELINE_CORE' else None,
            convention='all original records in source order; minimal source coefficient initialization excludes unused critical and ideal metadata; original computeZCoeff/computeHelmCoeff/alphar_derivatives/Fres/press/MBWR_pressure routines unmodified',
            native_formulation='explicit baseline Helmholtz summation follows original integrated coefficient expressions and order; direct pressure analytic; separate stable integral API; SI public density, mol/L backend; no scratch allocation')
        path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/mbwr_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(finite_cases=len(results),reference_errors=len(errors),failures=len(failures),stable_summation_differences=len(stable_differences)),indent=2))
    if failures:print(json.dumps(failures[:5],indent=2));raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory');validate(parser.parse_args().compiler_directory)
