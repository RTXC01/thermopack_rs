"""Qualify equilibrium TP caloric/volume derivatives against source helpers."""
import argparse,csv,hashlib,json,os,subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe
from validate_translation import encode


def cases():
    for eos in ['PR','SRK','RK','VDW']:
        for ref in ['Default'] if eos=='VDW' else ['Default','tcPR']:
            for comps,z in [('C1',[1.]),('CO2',[1.]),('H2',[1.]),('HE',[1.]),('C1,C2,C3',[.5,.3,.2]),
                            ('C1,CO2,N2',[0.,.8,.2]),('CO2,H2O',[.9,.1])]:
                for t in [180.,300.,500.]:
                    for p in [1e3,1e6,8e6]:
                        yield dict(mode='V_RESP',eos=eos,components=comps,feed=z,temperature=t,pressure=p,parameter_reference=ref)


def validate(directory):
    probe=verify_probe();inputs=[];expected=[];errors=[]
    with ReferenceProcess(directory) as oracle:
        for request in cases():
            result=oracle.call(request)
            if 'values' not in result or 'calculation_error' in result:
                errors.append({'request':request,**result});continue
            inputs.append(encode(request));expected.append(result['values'])
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    actual=subprocess.run([str(exe)],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(actual)==len(inputs);failures=[];absolute=[0.,1e-6,1e-10,1e-14]
    for case,reference,line in zip(inputs,expected,actual):
        if line.startswith('ERROR'):
            failures.append(dict(case=case,error=line));continue
        native=list(map(float,line.split(',')));assert len(native)==4
        for i,(a,b) in enumerate(zip(reference,native)):
            tolerance=absolute[i]+(0. if i==0 else 2e-7*max(abs(a),abs(b)))
            if not abs(a-b)<=tolerance:failures.append(dict(case=case,quantity=i,reference=a,rust=b))
    path=ROOT/'tests/fixtures/equilibrium_response.csv'
    with path.open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter=';',lineterminator='\n');writer.writerow(['case','phase,Cp,dVdT,dVdP'])
        writer.writerows([[case,*values] for case,values in zip(inputs,expected)])
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update(cases=len(inputs),attempted_cases=len(inputs)+len(errors),reference_errors=errors,probe=probe,
                    absolute_tolerances=absolute,relative_tolerance=2e-7,phase_flag_tolerance=0,
                    convention='baseline reported Cp; source state_functions helpers after exact-source TP flash',
                    fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),generator='tools/rust_port/validate_response.py')
    path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/response_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(cases=len(inputs),reference_errors=len(errors),failures=len(failures)),indent=2))
    if failures:print(json.dumps(failures[:6],indent=2));raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
