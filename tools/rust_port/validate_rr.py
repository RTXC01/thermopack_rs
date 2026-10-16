"""Qualify the fork's extended and one-step RR algorithm, without an EOS."""
import argparse,csv,hashlib,json,os,random,subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe


def cases():
    randomizer=random.Random(86021)
    for n in [2,3,8,20]:
        for _ in range(80):
            z=[randomizer.random() for _ in range(n)]
            total=sum(z);z=[v/total for v in z]
            k=[10.**randomizer.uniform(-12.,12.) for _ in range(n)]
            for step in [False,True]:
                yield dict(mode='RR',feed=z,k=k,single_step=step,beta=randomizer.uniform(.05,.95))
    feeds=[[.5,.5],[.2,.8],[.01,.99],[0.,.4,.6]]
    for z in feeds:
        for low,high in [(.1,2.),(.9,1.1),(1e-100,1e100),(0.,2.),(1.,1.),(.1,.9),(2.,3.),(1.-1e-15,1.+1e-15)]:
            k=[high]+[low]*(len(z)-1)
            for step in [False,True]:
                yield dict(mode='RR',feed=z,k=k,single_step=step,beta=.5)


def encode(request):
    return '|'.join(['RR',','.join(map(str,request['feed'])),','.join(map(str,request['k'])),
                     str(request['beta']),str(int(request['single_step'])),'RR'])


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
    process=subprocess.run([str(exe)],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True)
    lines=process.stdout.splitlines();assert len(lines)==len(inputs)
    failures=[]
    for case,ref,line in zip(inputs,expected,lines):
        if line.startswith('ERROR'):
            failures.append({'case':case,'error':line});continue
        native=list(map(float,line.split(',')))
        if len(native)!=len(ref):
            failures.append({'case':case,'reference':ref,'rust':native});continue
        for i,(a,b) in enumerate(zip(ref,native)):
            tol=0. if i==0 else 1e-12+2e-10*max(abs(a),abs(b))
            if not abs(a-b)<=tol:failures.append({'case':case,'quantity':i,'reference':a,'rust':b})
    path=ROOT/'tests/fixtures/rachford_rice.csv'
    with path.open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter=';',lineterminator='\n');writer.writerow(['case','success,betaV,betaL,x[n],y[n]...'])
        writer.writerows([[case,*values] for case,values in zip(inputs,expected)])
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update(cases=len(inputs),attempted_cases=len(inputs)+len(errors),reference_errors=errors,probe=probe,
                    absolute_tolerance=1e-12,relative_tolerance=2e-10,success_flag_tolerance=0,
                    false_return_scope='Only the success flag is qualified; source output fields are undefined on false.',
                    fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),generator='tools/rust_port/validate_rr.py')
    path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/rr_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(cases=len(inputs),reference_errors=len(errors),failures=len(failures)),indent=2))
    if failures:print(json.dumps(failures[:8],indent=2));raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
