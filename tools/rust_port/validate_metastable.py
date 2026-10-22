"""Local PS solves and source-order metastable isentropes of this exact fork."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess
RTOL=2e-6
ATOL=[2e-5,1e-10,2.0]

def local_cases():
    mixtures=[('C1',[1.]),('CO2',[1.]),('C3',[1.]),('N2',[1.]),('C1,C2',[.5,.5]),('C1,C3',[.5,.5]),
              ('C1,C2,C3',[.2,.3,.5]),('C1,C2,C3',[.5,.5,0.])]
    for eos in ['PR','SRK','RK','VDW']:
        for names,z in mixtures:
            for t in [150.,250.,400.]:
                for v in [8e-5,2e-4,2e-3]:yield dict(eos=eos,components=names.split(','),feed=z,configuration='classic',temperature=t,volume=v)
    for eos in ['PR','SRK']:
        for names,z in [('CO2',[1.]),('C1,C2',[.5,.5]),('C1,C2,C3',[.2,.3,.5])]:
            for t in [150.,250.,400.]:
                for v in [8e-5,2e-4,2e-3]:yield dict(eos=eos,components=names.split(','),feed=z,configuration='translated',temperature=t,volume=v)
        for config,names,z in [('unifac','C1,C3',[.5,.5]),('hv1','C1,H2O',[.5,.5])]:
            for t in [250.,400.]:
                for v in [2e-4,2e-3]:yield dict(eos=eos,components=names.split(','),feed=z,configuration=config,temperature=t,volume=v)
    for names,z in [('H2',[1.]),('H2,NE',[.5,.5])]:
        for t in [20.,50.,120.]:
            for v in [2e-4,2e-3]:yield dict(eos='PR',components=names.split(','),feed=z,configuration='quantum',temperature=t,volume=v,temperature_guess=max(110.,1.02*t))

def curve_cases():
    fluids=[('C1',[180.,220.,260.],6e6),('CO2',[270.,310.,350.],1e7),('C3',[330.,370.,410.],8e6),('N2',[110.,130.,150.],6e6)]
    for eos in ['PR','SRK','RK','VDW']:
        for name,temperatures,p in fluids:
            for t in temperatures:
                for pmin in [1e5,1e6]:yield dict(eos=eos,components=[name],feed=[1.],configuration='classic',temperature=t,volume=0.,initial_pressure=p,minimum_pressure=pmin,points=12)
        for names,z in [('C1,C2',[.5,.5]),('C1,C3',[.5,.5]),('CO2,C3',[.5,.5]),('C1,C2,C3',[.2,.3,.5])]:
            for t in [250.,350.]:yield dict(eos=eos,components=names.split(','),feed=z,configuration='classic',temperature=t,volume=0.,initial_pressure=1e7,minimum_pressure=1e5,points=12)
    for eos in ['PR','SRK']:
        for names,z in [('CO2',[1.]),('C1,C2',[.5,.5]),('C1,C2,C3',[.2,.3,.5])]:
            for t in [250.,350.]:yield dict(eos=eos,components=names.split(','),feed=z,configuration='translated',temperature=t,volume=0.,initial_pressure=1e7,minimum_pressure=1e5,points=12)

def native_line(request,result):
    p=result['pressure'];s=result['entropy'];t=request['temperature'];v=request['volume']
    if request['mode']=='META_PS':method=f"PS:{p!r}:{s!r}:{request.get('temperature_guess',1.02*t)!r}:{request.get('volume_guess',.98*v)!r}"
    else:method=f"CURVE:{p!r}:{s!r}:{request['minimum_pressure']!r}:{request['points']}"
    return f"META|{request['eos']}|{','.join(request['components'])}|{request['configuration']}|{method}|"+','.join(map(str,request['feed']))

def validate(directory):
    good=[];failed=[];unqualified=[];censored=[];warnings=[];attempts=[]
    with DeadlineReferenceProcess(directory) as oracle:
        for mode,cases in [('META_PS',local_cases()),('META_CURVE',curve_cases())]:
            for case in cases:
                request=dict(case,mode=mode);result=oracle.call(request,timeout=8.);record=dict(request=request,result=result);attempts.append(record)
                if 'reference_timeout' in result:censored.append(record);continue
                if 'values' not in result or result.get('status',0):
                    for message in result.get('reference_error',{}).get('output',[]):
                        if message.startswith('THERMOPACK_META_SPEC='):result.update(json.loads(message.split('=',1)[1]))
                    failed.append(record);continue
                values=result['values'];scale_p=max(abs(result['pressure']),1e5);scale_s=max(abs(result['entropy']),1.)
                if mode=='META_PS':
                    reasons=[]
                    if not all(map(math.isfinite,values)) or min(values)<=0:reasons.append('invalid state')
                    if abs(result.get('pressure_error',math.inf))>1e-7*scale_p or abs(result.get('entropy_error',math.inf))>1e-7*scale_s:reasons.append('PS does not close')
                else:
                    reasons=[]
                    if len(values)!=case['points'] or any(not all(map(math.isfinite,row)) or row[0]<=0 or row[1]<=0 for row in values):reasons.append('invalid curve')
                    for i,row in enumerate(values):
                        if abs(row[3])>1e-7*max(abs(row[2]),1e5) or abs(row[4])>1e-7*scale_s:
                            if i in [0,len(values)-1]:warnings.append(dict(request=request,index=i,values=row,entropy=result['entropy']))
                            else:reasons.append('interior PS does not close')
                if reasons:record['reasons']=reasons;unqualified.append(record);continue
                record['case']=native_line(request,result);good.append(record)
            print(mode,'good',len(good),'failed',len(failed),'unqualified',len(unqualified),'censored',len(censored),flush=True)
    subprocess.run(['cargo','build','--release','--example','parity','-j','2'],cwd=ROOT,check=True)
    eligible=[r for r in attempts if 'pressure' in r['result'] and 'entropy' in r['result']]
    lines=[native_line(r['request'],r['result']) for r in eligible]
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    output=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',capture_output=True,text=True,check=True).stdout.splitlines()
    assert len(output)==len(lines);native=dict(zip(lines,output));differences=[];local_rows=[];curve_rows=[];maxima={'local':[0.,0.],'curve':[0.]*3}
    for record in good:
        label=record['case'];result=record['result'];actual=native[label];local=record['request']['mode']=='META_PS'
        if local:local_rows.append([label,*result['values'],result['pressure'],result['entropy']])
        else:
            for i,row in enumerate(result['values']):curve_rows.append([label,i,*row,result['entropy']])
        if actual.startswith('ERROR'):differences.append(dict(case=label,rust=actual));continue
        actual=list(map(float,actual.split(',')));expected=[result['values']] if local else result['values'];width=7 if local else 3
        assert len(actual)==width*len(expected),(label,len(actual),len(expected))
        bad=[]
        for i,row in enumerate(expected):
            for j,(a,b) in enumerate(zip(row[:2] if local else row[:3],actual[i*width:(i+1)*width])):
                maxima['local' if local else 'curve'][j]=max(maxima['local' if local else 'curve'][j],abs(a-b)/max(abs(a),abs(b),ATOL[j]))
                if not math.isfinite(b) or abs(a-b)>ATOL[j]+RTOL*max(abs(a),abs(b)):bad.append([i,j,a,b])
        if bad:differences.append(dict(case=label,quantities=bad))
    for record in failed+unqualified+censored:
        if 'pressure' in record['result'] and 'entropy' in record['result']:record['native']=native[native_line(record['request'],record['result'])][:800]
    meta=json.loads((ROOT/'.reference/metadata.json').read_text())
    for name,rows,header in [('homogeneous_ps',local_rows,['case','T','V','P','S']),('metastable_isentrope',curve_rows,['case','index','T','V','P','pressure_error','entropy_error','S'])]:
        path=ROOT/f'tests/fixtures/{name}.csv'
        with path.open('w',newline='') as f:w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(header);w.writerows(rows)
        data=dict(meta,cases=len(rows),generator='tools/rust_port/validate_metastable.py',fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            relative_tolerance=RTOL,absolute_tolerances=ATOL,reference_errors=failed,unqualified_returns=unqualified,censored_calls=censored,
            source_endpoint_warnings=warnings,worst_relative_differences=maxima['local' if name=='homogeneous_ps' else 'curve'])
        path.with_suffix('.metadata.json').write_text(json.dumps(data,indent=2)+'\n')
    (ROOT/'.reference/metastable_differences.json').write_text(json.dumps(differences,indent=2)+'\n')
    print(json.dumps(dict(attempts=len(attempts),qualified_local=len(local_rows),qualified_curve_states=len(curve_rows),failed=len(failed),unqualified=len(unqualified),
        censored=len(censored),endpoint_warnings=len(warnings),differences=len(differences),worst_relative_differences=maxima),indent=2))
    if differences:print(json.dumps(differences[:6],indent=2));raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');validate(p.parse_args().compiler_directory)
