"""Qualify all three FMT functionals, source output branches and grid layout."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe

RTOL=2e-10
ATOL=1e-10

def grid():
    for eta in [0.,1e-16,1e-12,1e-8,1e-4,.01,.1,.4,.8,.99,.99999,1.-1e-9,1.-1e-12,1.,1.1,-.1]:
        for scale in [.01,1.,100.]:
            yield [.1*scale,.2*scale,.3*scale,eta,.02*scale,-.04*scale]
        # Uniform hard spheres of radius one, no vector contribution.
        yield [3.*eta/(4.*math.pi),3.*eta/(4.*math.pi),3.*eta,eta,0.,0.]

def validate(directory):
    probe=verify_probe();rows=[];errors=[];attempts=[]
    inputs=list(grid())
    with ReferenceProcess(directory) as oracle:
        for model in [1,2,3]:
            for flags in range(1,8):
                result=oracle.call(dict(mode='FMT',model=model,flags=flags,grid=inputs))
                assert 'rows' in result,result
                assert len(result['rows'])==len(inputs)
                for n,values in zip(inputs,result['rows']):
                    case=f'FMT|{model}|'+','.join(map(str,n))+f'|{flags}|0|STATE'
                    attempts.append(case)
                    if not all(math.isfinite(float(v)) for v in values):
                        errors.append(dict(case=case,values=values));continue
                    rows.append([case,*values])
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    lines=subprocess.run([str(exe)],input='\n'.join(attempts)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(lines)==len(attempts)
    actual=dict(zip(attempts,lines));failures=[];maximum=[0.]*43;scaled=[0.]*43
    for row in rows:
        line=actual[row[0]]
        if line.startswith('ERROR'):
            failures.append(dict(case=row[0],error=line));continue
        values=list(map(float,line.split(',')));assert len(values)==43
        bad=[]
        for i,(a,b) in enumerate(zip(row[1:],values)):
            delta=abs(a-b);maximum[i]=max(maximum[i],delta)
            ratio=delta/(ATOL+RTOL*max(abs(a),abs(b)))
            scaled[i]=max(scaled[i],ratio)
            if not math.isfinite(b) or ratio>1.:bad.append(i)
        if bad:failures.append(dict(case=row[0],quantities=bad,reference=row[1:],rust=values))
    for row in errors:
        if not actual[row['case']].startswith('ERROR'):
            failures.append(dict(case=row['case'],error='native accepted non-finite source state',rust=actual[row['case']]))
    path=ROOT/'tests/fixtures/fmt.csv'
    fields=['Phi',*[f'd{i}' for i in range(6)],*[f'd{i}d{j}' for i in range(6) for j in range(6)]]
    with path.open('w',newline='') as f:
        w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(['case',*fields]);w.writerows(rows)
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update(cases=len(rows),attempted_cases=len(attempts),reference_errors=errors,probe=probe,
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),generator='tools/rust_port/validate_fmt.py',
        quantities=fields,absolute_tolerance=ATOL,relative_tolerance=RTOL,
        maximum_absolute_differences=maximum,maximum_tolerance_fractions=scaled,
        convention='source FMT functional, six weighted densities, full/value/gradient/Hessian optional outputs; source max and hyperdual reciprocal derivative floors retained by baseline API; singular WB/WBII n3=0 returns recorded as native errors')
    path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/fmt_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(cases=len(rows),reference_errors=len(errors),failures=len(failures),worst_tolerance_fraction=max(scaled)),indent=2))
    if failures:print(json.dumps(failures[:2],indent=2));raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');validate(p.parse_args().compiler_directory)
