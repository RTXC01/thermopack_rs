"""PH/PS temperature recovery and pure latent heat against the exact fork."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
from reference import ROOT,ReferenceProcess


def cases():
    fluids={json.loads(path.read_text())['ident']:json.loads(path.read_text())
            for path in (ROOT/'fluids').glob('*.json')}
    for eos in ['PR','SRK','RK','VDW']:
        for mode in ['PH','PS']:
            for comps,z in [('C1',[1.]),('CO2',[1.]),('C3',[1.]),
                            ('C1,C2',[0.5,0.5]),('C1,C2,C3',[0.5,0.3,0.2]),
                            ('CO2,N2',[0.8,0.2]),('CO2,H2O',[0.9,0.1]),
                            ('C1,CO2,N2',[0.,0.8,0.2])]:
                for t in [150.,220.,300.,450.]:
                    for p in [1e5,1e6,8e6]:
                        yield {'mode':mode,'eos':eos,'components':comps,'feed':z,
                               'temperature':t,'pressure':p}
            for name in ['C1','CO2','C3']:
                for tr in [0.7,0.9]:
                    t=tr*fluids[name]['critical']['temperature']
                    for quality in [0.1,0.5,0.9]:
                        yield {'mode':mode,'eos':eos,'components':name,'feed':[1.],
                               'temperature':t,'quality':quality}


def validate(compiler_directory):
    inputs,expected,swappable=[],[],[]
    reference_errors=[]
    with ReferenceProcess(compiler_directory) as oracle:
        for request in cases():
            result=oracle.call(request)
            if 'reference_error' in result or 'calculation_error' in result:
                reference_errors.append({'request':request,**result})
                continue
            case='|'.join([request['eos'],request['components'],str(result['specification']),
                           str(result['pressure']),','.join(map(str,request['feed'])),request['mode']])
            inputs.append(case)
            expected.append(result['values'])
            swappable.append(result['phase_permutation_allowed'])
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    executable=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    result=subprocess.run([str(executable)],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True)
    actual=result.stdout.splitlines()
    assert len(actual)==len(expected)
    failures=[]
    for case,ref,line,can_swap in zip(inputs,expected,actual,swappable):
        if line.startswith('ERROR'):
            failures.append({'case':case,'error':line})
            continue
        native=list(map(float,line.split(',')))
        assert len(ref)==len(native)
        for i,(a,b) in enumerate(zip(ref,native)):
            tolerance=0.0 if i==1 else 2e-6+2e-7*max(abs(a),abs(b))
            if not abs(a-b)<=tolerance:
                failures.append({'case':case,'quantity':i,'reference':a,'rust':b})
    output=ROOT/'tests/fixtures/caloric_flash.csv'
    with output.open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter=';')
        writer.writerow(['case','phase_permutation_allowed','T,phase,betaV,betaL,x[n],y[n]...'])
        writer.writerows([[case,int(can_swap),*row] for case,row,can_swap in zip(inputs,expected,swappable)])
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update({'cases':len(inputs),'reference_errors':reference_errors,
                     'attempted_cases':len(inputs)+len(reference_errors),
                     'source_closure_scaled_tolerance':1e-8,
                     'absolute_tolerance':2e-6,'relative_tolerance':2e-7,'phase_flag_tolerance':0,
                     'fixture_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
                     'generator':'tools/rust_port/validate_caloric.py'})
    output.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/caloric_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps({'cases':len(inputs),'reference_errors':len(reference_errors),'failures':len(failures)},indent=2))
    if failures:
        print(json.dumps(failures[:16],indent=2))
        raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
