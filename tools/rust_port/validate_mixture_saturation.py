"""Mixture bubble/dew pressure/temperature against the exact-source wrappers."""
import argparse,csv,hashlib,json,math,os,subprocess
from collections import Counter
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess

MIXTURES=[
    ('C1,C2',[130.,180.,220.],[[.2,.8],[.5,.5],[.8,.2]]),
    ('C1,C3',[160.,220.,280.],[[.2,.8],[.5,.5],[.8,.2]]),
    ('CO2,C3',[220.,270.,310.],[[.2,.8],[.5,.5],[.8,.2]]),
    ('C2,C3',[210.,270.,320.],[[.2,.8],[.5,.5]]),
    ('C1,C2,C3',[170.,220.,270.],[[.2,.3,.5],[.7,.2,.1],[.5,.5,0.],[0.,1.,0.]]),
]
RTOL=2e-7
ATOL=[1e-6,1e-2,1e-12,1e-12]

def configurations():
    for eos in ['PR','SRK','RK','VDW']:
        for item in MIXTURES:yield eos,False,False,item
    for eos in ['PR','SRK']:
        for item in MIXTURES[:3]:yield eos,True,False,item
        for item in [MIXTURES[1],('NC4,NC5',[220.,270.,320.],[[.3,.7],[.5,.5]])]:yield eos,False,True,item

def case_line(eos,shift,unifac,comps,spec,feed,method):
    config=f'V_MIX_{method}@CLASSIC@DEFAULT@1' if shift else f'MIX_{method}'+('@unifac' if unifac else '')
    return f'{eos}|{comps}|{spec!r}|0|'+','.join(map(str,feed))+'|'+config

def validate(directory):
    good=[];failed=[];censored=[];unqualified=[];root_warnings=[];all_cases=[]
    with DeadlineReferenceProcess(directory) as oracle:
        def calculate(eos,shift,unifac,comps,spec,feed,method):
            case=case_line(eos,shift,unifac,comps,spec,feed,method);all_cases.append(case)
            request=dict(mode='MIXSAT',eos=eos,components=comps,volume_shift=shift,calculation=method,specification=spec,feed=feed)
            if unifac:request['mixing']='UNIFAC'
            result=oracle.call(request,timeout=8.)
            record=dict(case=case,result=result)
            if 'reference_timeout' in result:censored.append(record);return None
            if 'values' not in result:failed.append(record);return None
            values=result['values']
            reasons=[]
            if not all(map(math.isfinite,values)) or values[0]<=0 or values[1]<=0:reasons.append('non-finite/nonpositive state')
            if result['fugacity_residual']>1e-6:reasons.append('fugacity does not close')
            if result['normalization_error']>1e-7:reasons.append('phase amounts not normalized')
            if result['phase_distance']<=1e-7:reasons.append('coincident or inverted phase volumes')
            n=len(feed)
            gap=max(abs(a-b) for a,b in zip(values[4:4+n],values[4+n:]))
            if result['phase_distance']<1e-4 and gap<1e-4:reasons.append('coalescing phases: boundary coordinates are below log-fugacity numerical resolution')
            if reasons:unqualified.append(dict(record,reasons=reasons));return None
            if result['pressure_error']>1e-7:root_warnings.append(record)
            good.append(record);return result
        for eos,shift,unifac,(comps,temperatures,feeds) in configurations():
            for t in temperatures:
                for feed in feeds:
                    for method,inverse in [('BP','BT'),('DP','DT')]:
                        result=calculate(eos,shift,unifac,comps,t,feed,method)
                        if result is not None:calculate(eos,shift,unifac,comps,result['values'][1],feed,inverse)
            print(f'{eos} shift={shift} UNIFAC={unifac} {comps}: {len(good)} closing, {len(failed)} failed, {len(unqualified)} unqualified, {len(censored)} censored',flush=True)
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    output=subprocess.run([str(exe)],input='\n'.join(all_cases)+'\n',capture_output=True,text=True,check=True).stdout.splitlines()
    assert len(output)==len(all_cases);native=dict(zip(all_cases,output))
    differences=[];maxima=[0.]*5;diagnostics=Counter();rows=[]
    for record in good:
        case=record['case'];reference=record['result']['values'];line=native[case]
        n=(len(reference)-4)//2
        rows.append([case,*reference[:4],','.join(map(str,reference[4:4+n])),','.join(map(str,reference[4+n:]))])
        if line.startswith('ERROR'):
            differences.append(dict(record,rust=line));continue
        actual=list(map(float,line.split(',')));assert len(actual)==7+2*n
        bad=[]
        for i,(a,b) in enumerate(zip(reference,actual[:4]+actual[7:])):
            atol=ATOL[i] if i<4 else 1e-9
            delta=abs(a-b);maxima[min(i,4)]=max(maxima[min(i,4)],delta/max(abs(a),abs(b),atol))
            if not math.isfinite(b) or delta>atol+RTOL*max(abs(a),abs(b)):bad.append(i)
        if bad:differences.append(dict(record,rust=actual,quantities=bad))
        assert actual[5]<3e-10 and actual[2]<actual[3],(case,actual)
        diagnostics[case.split('|')[-1].split('@')[0]]+=1
    for record in failed+unqualified+censored:
        record['rust']=native[record['case']]
        if not record['rust'].startswith('ERROR'):diagnostics['native_solution_for_unqualified_source']+=1
    path=ROOT/'tests/fixtures/mixture_saturation.csv'
    with path.open('w',newline='') as f:
        w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(['case','T','P','v_liquid','v_vapor','liquid','vapor']);w.writerows(rows)
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update(cases=len(rows),attempted_cases=len(all_cases),reference_errors=failed,unqualified_reference_returns=unqualified,censored_reference_calls=censored,reference_pressure_root_warnings=root_warnings,
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),generator='tools/rust_port/validate_mixture_saturation.py',relative_tolerance=RTOL,absolute_tolerances=dict(T=ATOL[0],P=ATOL[1],volume=ATOL[2],composition=1e-9),worst_relative_differences=maxima,native_diagnostics=dict(diagnostics),
        convention='exact-source safe bubble/dew wrappers, fixed-phase composition, phase fugacity/normalization and separation checked before boundary qualification; source volumes from the separate specific_volume routine are compared as returned and pressure-root residual warnings retained separately, not used to discard valid boundary coordinates; four classic cubics, translated PR/SRK and UNIFAC; raw failed/degenerate source returns retained separately; default native analytic Newton follows physical equality rather than raw source trajectories')
    path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/mixture_saturation_differences.json').write_text(json.dumps(differences,indent=2)+'\n')
    print(json.dumps(dict(attempts=len(all_cases),qualified_boundaries=len(good),failed=len(failed),unqualified=len(unqualified),censored=len(censored),source_pressure_root_warnings=len(root_warnings),differences=len(differences),worst_relative_differences=maxima,diagnostics=dict(diagnostics)),indent=2))
    if differences:print(json.dumps(differences[:3],indent=2));raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');validate(p.parse_args().compiler_directory)
