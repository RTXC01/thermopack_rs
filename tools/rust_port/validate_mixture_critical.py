"""Qualify local cubic mixture critical points against this exact fork.

Failed, censored, and noncritical source returns remain explicit diagnostics.
The source's own TV mole Hessian independently qualifies returned states.
"""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess

RTOL=1e-6
ATOL=[1e-6,1e-12,1e-2]

def configurations():
    binary=['C1,C2','C1,C3','CO2,C3','N2,C1','NC4,NC5']
    for eos in ['PR','SRK','RK','VDW']:
        for names in binary:
            for x in [.01,.1,.3,.5,.7,.9,.99]:yield eos,names,[x,1-x],'classic'
        for names,feeds in [
            ('C1,C2,C3',[[.2,.3,.5],[.7,.2,.1],[.01,.49,.5],[0,.3,.7],[0,1,0]]),
            ('C1,C2,C3,NC4,NC5',[[.1,.2,.3,.2,.2],[.6,.1,.1,.1,.1]])]:
            for feed in feeds:yield eos,names,feed,'classic'
    for eos in ['PR','SRK']:
        for names in binary[:3]:
            for x in [.1,.5,.9]:yield eos,names,[x,1-x],'translated'
        for names in ['CO2,NC4','CO2,C3']:
            for x in [.1,.5,.9]:yield eos,names,[x,1-x],'ws'

def case_line(eos,names,z,config):
    return f'CRITICAL|{eos}|{names}|'+','.join(map(str,z))+f'|{config}'

def validate(directory):
    cases=list(configurations());good=[];failed=[];unqualified=[];censored=[]
    with DeadlineReferenceProcess(directory) as oracle:
        for eos,names,z,config in cases:
            line=case_line(eos,names,z,config)
            result=oracle.call(dict(mode='CRITICAL',eos=eos,components=names.split(','),feed=z,
                volume_shift=config=='translated',mixing='WongSandler' if config=='ws' else 'vdW',
                stability_diagnostics=True),timeout=8.)
            record=dict(case=line,result=result)
            if 'reference_timeout' in result:censored.append(record)
            elif 'values' not in result:failed.append(record)
            elif (not all(math.isfinite(x) and x>0 for x in result['values'])
                  or abs(result['eigenvalue'])>2e-7 or abs(result['third_derivative'])>2e-6):
                unqualified.append(record)
            else:good.append(record)
    subprocess.run(['cargo','build','--release','--example','parity','-j','2'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    lines=[case_line(*c) for c in cases]
    output=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',capture_output=True,text=True,check=True).stdout.splitlines()
    assert len(output)==len(lines),(len(output),len(lines));native=dict(zip(lines,output))
    differences=[];maxima=[0.]*3;rows=[]
    for record in good:
        expected=record['result']['values'];actual=native[record['case']]
        rows.append([record['case'],*expected,','.join(map(str,record['result']['hessian'])),
                     record['result']['eigenvalue'],record['result']['third_derivative']])
        if actual.startswith('ERROR'):differences.append(dict(record,rust=actual));continue
        actual=list(map(float,actual.split(',')));assert len(actual)==3
        bad=[]
        for i,(a,b) in enumerate(zip(expected,actual)):
            maxima[i]=max(maxima[i],abs(a-b)/max(abs(a),abs(b),ATOL[i]))
            if not math.isfinite(b) or abs(a-b)>ATOL[i]+RTOL*max(abs(a),abs(b)):bad.append(i)
        if bad:differences.append(dict(record,rust=actual,quantities=bad))
    for record in failed+unqualified+censored:record['rust']=native[record['case']]
    path=ROOT/'tests/fixtures/mixture_critical.csv'
    with path.open('w',newline='') as f:
        w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(['case','T','V','P','active_hessian','eigenvalue','third_derivative']);w.writerows(rows)
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update(cases=len(rows),attempted_cases=len(cases),reference_errors=failed,
        unqualified_reference_returns=unqualified,censored_reference_calls=censored,
        generator='tools/rust_port/validate_mixture_critical.py',fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        relative_tolerance=RTOL,absolute_tolerances=ATOL,worst_relative_differences=maxima,
        convention='calcCriticalTV at normalized fixed composition; source TV Hessian checked for a small minimum eigenvalue and third directional derivative; active Hessian also retained independently of the native solved point; raw failed/unqualified/censored returns retained without parity claims')
    path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/mixture_critical_differences.json').write_text(json.dumps(differences,indent=2)+'\n')
    print(json.dumps(dict(attempts=len(cases),qualified=len(good),failed=len(failed),unqualified=len(unqualified),
                         censored=len(censored),differences=len(differences),worst_relative_differences=maxima),indent=2))
    if differences:print(json.dumps(differences[:3],indent=2));raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');validate(p.parse_args().compiler_directory)
