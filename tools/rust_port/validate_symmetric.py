"""Qualify native upper Bunch-Kaufman factorization/solve against pinned LAPACK."""
import argparse
import csv
import hashlib
import json
import os
import random
import subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe


def encode(request):
    return '|'.join(['SYMMETRIC',','.join(map(str,request['matrix'])),','.join(map(str,request['rhs'])),'0','0','SYMMETRIC'])


def cases():
    randomizer=random.Random(97055)
    for n in [1,2,3,8,20,32]:
        for _ in range(12):
            m=[[randomizer.uniform(-1.,1.) for _ in range(n)] for _ in range(n)]
            positive=[[sum(m[k][i]*m[k][j] for k in range(n))+(1. if i==j else 0.) for j in range(n)] for i in range(n)]
            indefinite=[[m[min(i,j)][max(i,j)] for j in range(n)] for i in range(n)]
            scales=[10.**randomizer.uniform(-8.,8.) for _ in range(n)]
            scaled=[[positive[i][j]*scales[i]*scales[j] for j in range(n)] for i in range(n)]
            for matrix in [positive,indefinite,scaled]:
                yield dict(mode='SYMMETRIC',matrix=[v for row in matrix for v in row],rhs=[randomizer.uniform(-1.,1.) for _ in range(n)])
    for matrix in [[0.],[2.],[0.,1.,1.,0.],[1.,1.,1.,1.],[0.,0.,0.,1.],
                   [1e-200,0.,0.,1e200],[1.,2.,2.,-1.],
                   [0.,0.,0.,0.,0.,1.,0.,1.,0.],[1.,2.,3.,2.,1.,4.,3.,4.,-1.],
                   [1.,2.,3.,-999.,-1.,4.,-777.,-888.,0.]]:
        n=int(len(matrix)**.5)
        yield dict(mode='SYMMETRIC',matrix=matrix,rhs=[float(i+1) for i in range(n)])


def validate(directory):
    probe=verify_probe();inputs=[];expected=[];errors=[]
    with ReferenceProcess(directory) as oracle:
        for request in cases():
            result=oracle.call(request)
            if 'values' not in result or 'calculation_error' in result:
                errors.append({'request':request,**result});continue
            inputs.append(encode(request));expected.append(result['values'])
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    executable=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    process=subprocess.run([str(executable)],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True)
    lines=process.stdout.splitlines();assert len(lines)==len(inputs)
    failures=[]
    for case,reference,line in zip(inputs,expected,lines):
        if line.startswith('ERROR'):
            failures.append(dict(case=case,error=line));continue
        actual=list(map(float,line.split(',')));assert len(reference)==len(actual)
        n=len(case.split('|')[2].split(','))
        for i,(a,b) in enumerate(zip(reference,actual)):
            tolerance=0. if i<=n else 1e-12+2e-12*max(abs(a),abs(b))
            if not abs(a-b)<=tolerance:failures.append(dict(case=case,quantity=i,reference=a,rust=b))
    path=ROOT/'tests/fixtures/symmetric.csv'
    with path.open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter=';',lineterminator='\n');writer.writerow(['case','success,pivots[n],factor[n*n],solution[n]...'])
        writer.writerows([[case,*values] for case,values in zip(inputs,expected)])
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update(cases=len(inputs),attempted_cases=len(inputs)+len(errors),reference_errors=errors,probe=probe,
                    absolute_tolerance=1e-12,relative_tolerance=2e-12,status_pivot_tolerance=0,
                    factorization='pinned LAPACK DSYTF2 upper, DSYTRS2; unblocked at every size',
                    fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),generator='tools/rust_port/validate_symmetric.py')
    path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/symmetric_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(cases=len(inputs),reference_errors=len(errors),failures=len(failures)),indent=2))
    if failures:print(json.dumps(failures[:4],indent=2));raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
