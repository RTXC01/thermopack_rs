"""Compare the source's LV RR substitution step, including flags and stale fields."""
import argparse
import csv
import hashlib
import json
import os
import random
import subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe
from validate_translation import encode as encode_model


def encode(request):
    return encode_model(request)+'@rr:'+','.join(map(str,request['k']))+':'+str(request.get('beta',.5))+':'+str(int(request.get('single_step',False)))


def cases():
    randomizer=random.Random(106251)
    feeds=[('C1,C2,C3',[.5,.3,.2]),('C1,CO2,N2',[.5,.3,.2]),('CO2,H2O',[.9,.1]),
           ('HE,H2,NE',[0.,.5,.5]),('H2,NE',[.6,.4]),('C1,C2',[.2,.8])]
    for eos in ['PR','SRK','RK','VDW']:
        for ref in ['Default'] if eos=='VDW' else ['Default','tcPR','QuantumCubic']:
            for comps,z in feeds:
                for t in [3.,20.,80.] if ref=='QuantumCubic' else [180.,300.,500.]:
                    for p in [1e3,1e6]:
                        n=len(z)
                        for k in [[4.]+[.2]*(n-1),[1.]*n,[.5]*n,[2.]*n,
                                  [10.**randomizer.uniform(-4.,4.) for _ in z]]:
                            yield dict(mode='V_RR',eos=eos,components=comps,feed=z,temperature=t,pressure=p,
                                       parameter_reference=ref,k=k)
                        # At beta=.5 these K values give a zero RR residual for
                        # an equimolar binary; the update stays in the extended
                        # domain for the remaining normalized feeds as well.
                        yield dict(mode='V_RR',eos=eos,components=comps,feed=z,temperature=t,pressure=p,
                                   parameter_reference=ref,k=[2.]+[.5]*(n-1),single_step=True)
    # Near-trivial factors amplify cancellation in extended fractions.
    for eos in ['PR','SRK']:
        for delta in [1e-12,2.**-40]:
            yield dict(mode='V_RR',eos=eos,components='C1,C2',feed=[.5,.5],temperature=500.,pressure=1e3,
                       k=[1.+delta,1.-delta])


def validate(directory):
    probe=verify_probe()
    inputs=[];expected=[];errors=[];attempted=0
    with ReferenceProcess(directory) as oracle:
        for request in cases():
            attempted+=1
            result=oracle.call(request)
            if 'values' not in result or 'calculation_error' in result:
                errors.append({'request':request,**result});continue
            inputs.append(encode(request));expected.append(result['values'])
        # Converged input factors exercise the final RR call and convergence
        # flag. Derive them from the source, retaining finite source failures.
        for eos in ['PR','SRK','RK','VDW']:
            for t in [180.,250.,300.]:
                request=dict(mode='V_TP',eos=eos,components='C1,C2,C3',feed=[.5,.3,.2],temperature=t,pressure=1e6)
                result=oracle.call(request)
                if 'values' not in result or result['values'][0]!=0:continue
                values=result['values'];x=values[3:6];y=values[6:9]
                request.update(mode='V_RR',k=[yi/xi for xi,yi in zip(x,y)],beta=values[1])
                attempted+=1
                result=oracle.call(request)
                if 'values' not in result or 'calculation_error' in result:
                    errors.append({'request':request,**result});continue
                inputs.append(encode(request));expected.append(result['values'])
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    native=subprocess.run([str(exe)],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(native)==len(inputs)
    failures=[]
    for case,ref,line in zip(inputs,expected,native):
        if line.startswith('ERROR'):
            failures.append({'case':case,'error':line});continue
        actual=list(map(float,line.split(',')))
        if len(actual)!=len(ref):
            failures.append({'case':case,'reference':ref,'rust':actual});continue
        for i,(a,b) in enumerate(zip(ref,actual)):
            tolerance=0. if i<2 else 1e-9+2e-8*max(abs(a),abs(b))
            if not abs(a-b)<=tolerance:
                failures.append({'case':case,'quantity':i,'reference':a,'rust':b})
    path=ROOT/'tests/fixtures/rr_iteration.csv'
    with path.open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter=';',lineterminator='\n');writer.writerow(['case','converged,has_solution,betaV,betaL,g,x[n],y[n],K[n],phiL[n],phiV[n]...'])
        writer.writerows([[case,*values] for case,values in zip(inputs,expected)])
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update(cases=len(inputs),attempted_cases=attempted,reference_errors=errors,probe=probe,
                    absolute_tolerance=1e-9,relative_tolerance=2e-8,flag_tolerance=0,
                    false_return_scope='Only convergence and RR status flags qualified when RR has no solution.',
                    fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),generator='tools/rust_port/validate_rr_iteration.py')
    path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/rr_iteration_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(cases=len(inputs),reference_errors=len(errors),failures=len(failures)),indent=2))
    if failures:
        print(json.dumps(failures[:8],indent=2));raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
