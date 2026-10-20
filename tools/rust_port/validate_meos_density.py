"""Original pure MEOS density iteration; explicit source failures and deadlines."""
import argparse,csv,hashlib,json,math,os,subprocess
from collections import Counter
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess
from validate_alpha import verify_probe
from validate_meos import MODELS

ATOL=[1e-8,0.,0.,1e-4,1e-5]
RTOL=2e-8
# The density of a Newton step inherits the pressure uncertainty of the
# residual kernel, which the kernel qualification bounds at about 1e-14
# relative. Near the critical point dP/drho is small, so that uncertainty is
# amplified into the density. The density tolerance therefore adds the
# conditioning term EPS_P*|p|/|dP/drho|.
EPS_P=1e-14

def cases(meta):
    tc,pc,rc=meta[:3]
    if pc<=0.:pc=.3*rc*meta[10]*tc
    for tr in [.6,.8,.95,1.,1.05,1.5,3.]:
        for pr in [.0001,.001,.01,.1,.3,1.,3.,10.]:
            for phase in [1,2]:yield tc*tr,pc*pr,phase,8.
    # Defined raw flags above Tc, unsupported flags below Tc, and selected
    # below-triple extrapolations. Deadlines remain censored observations.
    for phase in [-7,3,4]:yield tc*1.5,pc*.1,phase,8.
    for phase in [3,4]:yield tc*.8,pc*.1,phase,8.
    for phase in [1,2]:
        for pr in [.001,1.]:yield tc*.3,pc*pr,phase,2.

def validate(directory):
    probe=verify_probe();finite=[];failures=[];censored=[];attempts=[];source_unqualified=[]
    with DeadlineReferenceProcess(directory) as oracle:
        for e,c in MODELS:
            base=dict(mode='MEOS',eos=e,component=c,temperature=300.,volume=1.)
            init=oracle.call(base);assert 'metadata' in init,init
            for t,p,phase,deadline in cases(init['metadata']):
                case=f'MEOS|{e}:{c}|{t!r}|{p!r}|{phase}|DENSITY';attempts.append(case)
                result=oracle.call(dict(mode='MEOS_DENSITY',eos=e,component=c,temperature=t,pressure=p,phase=phase),timeout=deadline)
                if 'reference_timeout' in result:censored.append(dict(case=case,result=result));continue
                if 'reference_error' in result:failures.append(dict(case=case,result=result));continue
                assert 'values' in result,result
                if any(v is None for v in result['values']):source_unqualified.append(dict(case=case,result=result));continue
                finite.append([case,*result['values']])
            print(f'Density oracle {e}:{c}: {len(finite)} finite, {len(failures)} failed, {len(censored)} censored',flush=True)
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    lines=subprocess.run([str(exe)],input='\n'.join(attempts)+'\n',capture_output=True,text=True,check=True).stdout.splitlines()
    assert len(lines)==len(attempts);native=dict(zip(attempts,lines))
    differences=[];maximum=[0.]*5;counts=Counter();diagnostics=[]
    for row in finite:
        line=native[row[0]]
        if line.startswith('ERROR'):
            differences.append(dict(case=row[0],reference=row[1:],error=line));continue
        actual=list(map(float,line.split(',')));assert len(actual)==12
        bad=[]
        for i,(a,b) in enumerate(zip(row[1:],actual[:5])):
            maximum[i]=max(maximum[i],abs(a-b))
            tolerance=ATOL[i]+RTOL*max(abs(a),abs(b))
            if i==0:
                conditioning=EPS_P*abs(row[4])/abs(row[5]) if row[5] else float('inf')
                tolerance+=conditioning
            if not math.isfinite(b) or abs(a-b)>tolerance:bad.append(i)
        if bad:differences.append(dict(case=row[0],quantities=bad,reference=row[1:],rust=actual))
        counts[f'phase_switches_{int(actual[6])}']+=1
        counts[f'source_error_{int(actual[2])}']+=1
        if actual[7]!=1.:counts['unclosed']+=1
        if actual[8]!=1.:counts['unstable']+=1
        if actual[11]==1.:counts['relaxed_convergence']+=1
        phase=int(row[0].split('|')[4])
        strict_expected=actual[7]==1. and actual[8]==1. and actual[2]==0. and phase in [1,2,3]
        if strict_expected:
            assert math.isfinite(actual[10]) and abs(actual[0]*actual[10]-1.)<1e-12,(row[0],actual)
            counts['strict_volumes']+=1
        else:assert not math.isfinite(actual[10]),(row[0],actual)
        diagnostics.append(dict(case=row[0],iterations=int(actual[5]),phase_switches=int(actual[6]),relative_pressure_error=actual[9]))
    for row in failures+source_unqualified:
        row['rust']=native[row['case']]
        if not row['rust'].startswith('ERROR'):counts['finite_native_for_failed_source']+=1
    for row in censored:row['rust']=native[row['case']]
    path=ROOT/'tests/fixtures/meos_density.csv'
    with path.open('w',newline='') as f:
        w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(['case','rho','phase','ierr','pressure','pressure_rho']);w.writerows(finite)
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update(cases=len(finite),attempted_cases=len(attempts),reference_errors=failures,censored_reference_calls=censored,unqualified_reference_returns=source_unqualified,probe=probe,
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),generator='tools/rust_port/validate_meos_density.py',absolute_tolerances=ATOL,relative_tolerance=RTOL,maximum_absolute_differences=maximum,native_diagnostics=dict(counts),
        convention='original densitySolver with continueOnError=false, ierr present; raw density/phase/status and final original mp_pressure; execution deadlines are censored, intentionally terminated/reaped and not retried; native guard is 10000 steps and 1024 liquid seed doublings; strict volumes require final pressure closure, dP/drho>0, source ierr=0')
    path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/meos_density_differences.json').write_text(json.dumps(differences,indent=2)+'\n')
    (ROOT/'.reference/meos_density_diagnostics.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
    print(json.dumps(dict(attempts=len(attempts),finite=len(finite),reference_errors=len(failures),censored=len(censored),unqualified=len(source_unqualified),differences=len(differences),native_diagnostics=dict(counts)),indent=2))
    if differences:print(json.dumps(differences[:3],indent=2));raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');validate(p.parse_args().compiler_directory)
