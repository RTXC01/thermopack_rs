"""Qualify fixed-pressure spinodals and complete mixture volume traces."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess
RTOL=2e-6
ATOL=[2e-5,1e-10,2.0]

def configurations():
    for eos in ['PR','SRK','RK','VDW']:
        for names,z in [('C1,C2',[.2,.8]),('C1,C2',[.5,.5]),('C1,C3',[.5,.5]),('CO2,C3',[.5,.5]),
                        ('C1,C2,C3',[.2,.3,.5]),('C1,C2,C3',[0,.3,.7])]:
            for p in [1e6,3e6]:
                for phase in ['L','V']:yield eos,names,z,'classic',phase,p,.04
        for names,z in [('C1,C2',[.5,.5]),('C1,C3',[.5,.5]),('CO2,C3',[.5,.5]),('C1,C2,C3',[.2,.3,.5])]:
            yield eos,names,z,'classic','CURVE',1e5,.04
    for eos in ['PR','SRK']:
        for names,z in [('C1,C2',[.5,.5]),('CO2,C3',[.5,.5]),('C1,C2,C3',[.2,.3,.5])]:
            for phase in ['L','V']:yield eos,names,z,'translated',phase,1e6,.04
            yield eos,names,z,'translated','CURVE',1e5,.04

def line(eos,names,z,config,method,p,step):
    mode='BASE_'+method
    return f'MIXSPIN|{eos}|{names}|'+','.join(map(str,z))+f'|{config}|{mode}|{p}|{step}'

def validate(directory):
    cases=list(configurations());good=[];failed=[];unqualified=[];censored=[];warnings=[]
    with DeadlineReferenceProcess(directory) as oracle:
        for case in cases:
            eos,names,z,config,method,p,step=case;label=line(*case)
            request=dict(mode='SPINODAL_CURVE' if method=='CURVE' else 'SPINODAL',eos=eos,components=names.split(','),
                feed=z,volume_shift=config=='translated',pressure=p,phase=1 if method=='L' else 2,
                log_volume_step=step,stability_diagnostics=True)
            result=oracle.call(request,timeout=8.)
            record=dict(case=label,result=result)
            if 'reference_timeout' in result:censored.append(record)
            elif 'values' not in result:failed.append(record)
            elif any(not all(map(math.isfinite,r)) or r[0]<=0 or r[1]<=0 or abs(r[4])>1e-6 for r in result['values']):unqualified.append(record)
            else:
                good.append(record)
                for i,r in enumerate(result['values']):
                    if abs(r[2]-r[3])>1.0+1e-7*max(abs(r[2]),abs(r[3])):warnings.append(dict(case=label,index=i,values=r))
            print(f'{len(good)}/{len(failed)}/{len(censored)} {label}',flush=True)
    subprocess.run(['cargo','build','--release','--example','parity','-j','2'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    lines=[line(*c) for c in cases]
    output=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',capture_output=True,text=True,check=True).stdout.splitlines()
    assert len(output)==len(lines);native=dict(zip(lines,output));differences=[];rows=[];maxima=[0.]*3
    for record in good:
        label=record['case'];source=record['result']['values'];actual=native[label]
        for i,r in enumerate(source):rows.append([label,i,*r])
        if actual.startswith('ERROR'):differences.append(dict(case=label,rust=actual));continue
        flat=list(map(float,actual.split(',')));actual=[flat[i:i+3] for i in range(0,len(flat),3)]
        if len(actual)!=len(source):differences.append(dict(case=label,source_points=len(source),rust_points=len(actual)));continue
        bad=[]
        for i,(s,r) in enumerate(zip(source,actual)):
            for j,(a,b) in enumerate(zip([s[0],s[1],s[2]],r)):
                maxima[j]=max(maxima[j],abs(a-b)/max(abs(a),abs(b),ATOL[j]))
                if not math.isfinite(b) or abs(a-b)>ATOL[j]+RTOL*max(abs(a),abs(b)):bad.append([i,j,a,b])
        if bad:differences.append(dict(case=label,quantities=bad))
    for r in failed+unqualified+censored:r['rust']=native[r['case']][:500]
    physical_cases=[r for r in good if '|BASE_L|' in r['case'] or '|BASE_V|' in r['case']]
    physical_lines=[r['case'].replace('|BASE_L|','|L|').replace('|BASE_V|','|V|') for r in physical_cases]
    physical_output=subprocess.run([str(exe)],input='\n'.join(physical_lines)+'\n',capture_output=True,text=True,check=True).stdout.splitlines()
    assert len(physical_output)==len(physical_cases)
    branch_differences=[]
    for record,answer in zip(physical_cases,physical_output):
        assert not answer.startswith('ERROR'),(record['case'],answer)
        result=list(map(float,answer.split(',')));s=record['result']['values'][0]
        if any(abs(a-b)>ATOL[j]+RTOL*max(abs(a),abs(b)) for j,(a,b) in enumerate(zip(s[:3],result))):
            branch_differences.append(dict(case=record['case'],source=s[:3],physical_native=result))
    path=ROOT/'tests/fixtures/mixture_spinodal.csv'
    with path.open('w',newline='') as f:
        w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(['case','index','T','V','raw_P','evaluated_P','eigenvalue']);w.writerows(rows)
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update(cases=len(rows),attempted_calls=len(cases),qualified_calls=len(good),reference_errors=failed,
        unqualified_reference_returns=unqualified,censored_reference_calls=censored,reference_pressure_warnings=warnings,
        physical_branch_selection_differences=branch_differences,
        generator='tools/rust_port/validate_mixture_spinodal.py',fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        relative_tolerance=RTOL,absolute_tolerances=ATOL,worst_relative_differences=maxima,
        convention='original initial_stab_limit_point and map_stability_limit at normalized fixed composition; complete mixture volume grids, all source states retained; minimum Helmholtz eigenvalue independently checked; raw reported pressure retained separately from evaluated EOS pressure, including stale temperature endpoint pressures')
    path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/mixture_spinodal_differences.json').write_text(json.dumps(differences,indent=2)+'\n')
    print(json.dumps(dict(attempts=len(cases),qualified_calls=len(good),qualified_states=len(rows),failed=len(failed),unqualified=len(unqualified),
                         censored=len(censored),pressure_warnings=len(warnings),differences=len(differences),physical_selection_differences=len(branch_differences),worst_relative_differences=maxima),indent=2))
    if differences:print(json.dumps(differences[:8],indent=2));raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');validate(p.parse_args().compiler_directory)
