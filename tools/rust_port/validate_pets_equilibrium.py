"""Qualify PeTS roots, finite raw source iterates, critical and coexistence."""
import argparse,csv,hashlib,json,math,os,subprocess,sys
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe

MODES=['PETS_VOLUME_BASELINE','PETS_VOLUME','PETS_CRITICAL','PETS_SAT','PETS_SAT_T']
ABSOLUTE={'PETS_VOLUME_BASELINE':[1e-9,1e-3,1e-3,1.], 'PETS_VOLUME':[1e-10,1e-10,1e-14],
          'PETS_CRITICAL':[1e-5,.1,1e-11,1e-7], 'PETS_SAT':[1e-5,.01,1e-10,1e-9,.01],
          'PETS_SAT_T':[1e-5,.01,1e-10,1e-9,.01]}


def cases():
    for sigma,eps in [(3.42e-10,136.),(2.9e-10,100.),(4.1e-10,240.)]:
        params=[sigma,eps]
        yield dict(mode='PETS_CRITICAL',parameters=params,temperature=0.,pressure=0.)
        for tr in [.3,.5,.8,1.,1.09,1.5,3.]:
            for p in [.001,1.,1e3,1e5,1e6,5e6,1e8]:
                for phase in ['L','V']:
                    for mode in ['PETS_VOLUME_BASELINE','PETS_VOLUME']:
                        yield dict(mode=mode,parameters=params,temperature=tr*eps,pressure=p,phase=phase)
        # The real source setter does not refresh fallback critical data.
        # Explicitly invoke its public refresh for custom saturation models.
        for tr in [.35,.5,.7,.9,.99,.999]:
            for mode in ['PETS_SAT','PETS_SAT_T']:
                yield dict(mode=mode,parameters=params,temperature=148.11064915309302/136.*tr*eps,
                           pressure=0.,redefine_critical=True)


def encode(request,result):
    mode=request['mode'];p=result['values'][1] if mode=='PETS_SAT_T' else request['pressure']
    return '|'.join(['PETS',','.join(map(str,request['parameters'])),str(request['temperature']),str(p),'1',
                     mode+('@'+request['phase'] if 'phase' in request else '')])


def validate(directory):
    probe=verify_probe();records=[];errors=[]
    with ReferenceProcess(directory) as oracle:
        for request in cases():
            r=dict(request)
            if r['mode']=='PETS_SAT_T':
                initial=oracle.call({**r,'mode':'PETS_SAT'})
                if 'values' not in initial or 'calculation_error' in initial:
                    errors.append(dict(request=r,**initial));continue
                r['pressure']=initial['values'][1]
            result=oracle.call(r)
            if 'calculation_error' in result or 'values' not in result or not all(map(math.isfinite,result['values'])):
                errors.append(dict(request=r,**result));continue
            records.append((r,encode(r,result),result['values']))
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    lines=subprocess.run([str(exe)],input='\n'.join(c for _,c,_ in records)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(lines)==len(records);failures=[]
    for (r,case,expected),line in zip(records,lines):
        if line.startswith('ERROR'):
            failures.append(dict(case=case,error=line,reference=expected));continue
        actual=list(map(float,line.split(',')));assert len(actual)==len(expected)
        for i,(a,b) in enumerate(zip(expected,actual)):
            if not abs(a-b)<=ABSOLUTE[r['mode']][i]+2e-7*max(abs(a),abs(b)):
                failures.append(dict(case=case,quantity=i,reference=a,rust=b))
    names=['pets_volume_baseline','pets_volume','pets_critical','pets_saturation','pets_saturation_temperature']
    for mode,name in zip(MODES,names):
        group=[(case,values) for r,case,values in records if r['mode']==mode];path=ROOT/'tests/fixtures'/(name+'.csv')
        with path.open('w',newline='') as handle:
            writer=csv.writer(handle,delimiter=';',lineterminator='\n');writer.writerow(['case','values'])
            writer.writerows([[c,*v] for c,v in group])
        metadata=json.loads((ROOT/'.reference/metadata.json').read_text());reference_errors=[e for e in errors if e['request']['mode']==mode]
        metadata.update(cases=len(group),attempted_cases=len(group)+len(reference_errors),reference_errors=reference_errors,probe=probe,
                        absolute_tolerances=ABSOLUTE[mode],relative_tolerance=2e-7,
                        source_relative_pressure_closure_tolerance=1e8*sys.float_info.epsilon if mode=='PETS_VOLUME' else None,
                        generator='tools/rust_port/validate_pets_equilibrium.py',fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        convention='raw source iterates retained in baseline; physical-root group requires source pressure closure and negative pressure-volume derivative; explicit critical refresh before saturation')
        path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/pets_equilibrium_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(finite_cases=len(records),reference_errors=len(errors),failures=len(failures)),indent=2))
    if failures:print(json.dumps(failures[:10],indent=2));raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
