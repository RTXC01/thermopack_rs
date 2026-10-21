"""Exact-source pure grids and property intersection coordinates/statuses."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess
RTOL=2e-6
ATOL=[2e-5,1e-11,2.0]

def configurations():
    for eos in ['PR','SRK','RK','VDW']:
        for name in ['C1','CO2','C3','N2']:yield eos,name,'classic'
    for eos in ['PR','SRK']:
        for name in ['C1','CO2']:yield eos,name,'translated'

def label(eos,name,config,mode,p):return f'SPINPROP|{eos}|{name}|{config}|{mode}|{p}'

def validate(directory):
    grids=[];roots=[];failed=[];censored=[];p=1e6
    with DeadlineReferenceProcess(directory) as oracle:
        for eos,name,config in configurations():
            base=dict(eos=eos,components=[name],feed=[1.0],volume_shift=config=='translated',pressure=p)
            result=oracle.call(dict(base,mode='PURE_SPINODAL_GRID'),timeout=8.)
            record=dict(case=label(eos,name,config,'GRID',p),result=result)
            if 'reference_timeout' in result:censored.append(record);continue
            if 'values' not in result or result.get('status',0):failed.append(record);continue
            assert len(result['values'])==200
            assert all(all(map(math.isfinite,r)) and r[0]>0 and r[1]>0 for r in result['values'])
            record['case']=label(eos,name,config,f"GRID:{result['vapor_temperature']}",p)
            grids.append(record)
            if (eos,name,config) not in [('PR','C1','classic'),('PR','CO2','classic'),('SRK','C1','classic'),('SRK','N2','classic'),
                                        ('RK','C1','classic'),('VDW','C3','classic'),('PR','CO2','translated'),('SRK','C1','translated')]:continue
            for prop in range(1,7):
                for phase in [0,1,2]:
                    for index in [20,120]:
                        r=oracle.call(dict(base,mode='PURE_SPINODAL_PROP',property=prop,phase=phase,target_point=result['values'][index][:2]),timeout=8.)
                        rec=dict(configuration=[eos,name,config],property=prop,phase=phase,index=index,result=r)
                        if 'reference_timeout' in r:censored.append(rec);continue
                        if 'values' not in r or r.get('status',0):
                            for message in r.get('reference_error',{}).get('output',[]):
                                if message.startswith('THERMOPACK_TARGET='):rec['target']=float(message.split('=',1)[1])
                            failed.append(rec);continue
                        rec['case']=label(eos,name,config,f"ROOT:{prop}:{phase}:{r['target']!r}",p);roots.append(rec)
            print(eos,name,config,len(grids),len(roots),len(failed),len(censored),flush=True)
    subprocess.run(['cargo','build','--release','--example','parity','-j','2'],cwd=ROOT,check=True)
    records=grids+roots;exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    output=subprocess.run([str(exe)],input='\n'.join(r['case'] for r in records)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(output)==len(records);differences=[];grid_rows=[];root_rows=[];warnings=[];maxima={'grid':[0.]*3,'property':[0.]*3};conditioned=[]
    for rec,out in zip(records,output):
        label_=rec['case'];r=rec['result'];is_grid='|GRID:' in label_
        expected=r['values'] if is_grid else [[int(r['found']),*r['values']]]
        if is_grid:
            for i,s in enumerate(expected):grid_rows.append([label_,i,*s])
        else:root_rows.append([label_,int(r['found']),*r['values'],r.get('eigenvalue',0),r.get('property_value',0),r['target']])
        if out.startswith('ERROR'):differences.append(dict(case=label_,rust=out));continue
        a=list(map(float,out.split(',')));width=3 if is_grid else 4
        actual=[a[i:i+width] for i in range(0,len(a),width)];assert len(actual)==len(expected)
        bad=[]
        for i,(s,n) in enumerate(zip(expected,actual)):
            if not is_grid:
                if n[0]!=s[0]:bad.append(['found',s[0],n[0]])
                if not s[0]:continue
                s=s[1:];n=n[1:]
            for j,(v,w) in enumerate(zip(s[:3],n)):
                group=maxima['grid' if is_grid else 'property']
                group[j]=max(group[j],abs(v-w)/max(abs(v),abs(w),ATOL[j]))
                tolerance=ATOL[j]+RTOL*max(abs(v),abs(w))
                if is_grid and i in [99,100] and j==1:
                    # A monic cubic multiple root amplifies coefficient rounding
                    # as eps^(1/3). Only the stored critical entries use this
                    # explicit polynomial-conditioning radius in volume units.
                    radius=8.31446261815324*s[0]/s[2]*(64*2.220446049250313e-16)**(1/3)
                    if abs(v-w)>tolerance:conditioned.append(dict(case=label_,index=i,source=v,native=w,roundoff_radius=radius))
                    tolerance+=radius
                if not math.isfinite(w) or abs(v-w)>tolerance:bad.append([i,j,v,w])
        if bad:differences.append(dict(case=label_,quantities=bad))
        if not is_grid and r['found'] and abs(r['eigenvalue'])>1e-6:warnings.append(rec)
    diagnostics=[]
    for r in failed:
        target=r.get('target',r['result'].get('target'))
        if target is not None:
            eos,name,config=r['configuration'];r['case']=label(eos,name,config,f"ROOT:{r['property']}:{r['phase']}:{target!r}",p)
            diagnostics.append(r)
    if diagnostics:
        answers=subprocess.run([str(exe)],input='\n'.join(r['case'] for r in diagnostics)+'\n',capture_output=True,text=True,check=True).stdout.splitlines()
        assert len(answers)==len(diagnostics)
        for r,a in zip(diagnostics,answers):r['native']=a
    meta=json.loads((ROOT/'.reference/metadata.json').read_text())
    for name,rows,header in [
        ('pure_spinodal_grid',grid_rows,['case','index','T','V','raw_P','evaluated_P','eigenvalue']),
        ('spinodal_property',root_rows,['case','found','T','V','P','eigenvalue','property_value','target'])]:
        path=ROOT/f'tests/fixtures/{name}.csv'
        with path.open('w',newline='') as f:w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(header);w.writerows(rows)
        data=dict(meta,cases=len(rows),generator='tools/rust_port/validate_spinodal_properties.py',fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            relative_tolerance=RTOL,absolute_tolerances=ATOL,reference_errors=failed,censored_calls=censored,
            source_reported_coordinate_warnings=warnings if name=='spinodal_property' else [],worst_relative_differences=maxima['grid' if name=='pure_spinodal_grid' else 'property'])
        if name=='pure_spinodal_grid':data.update(critical_volume_roundoff='(R*T/P)*(64*machine_epsilon)^(1/3), added only to critical entries 99/100',roundoff_conditioned_states=conditioned)
        path.with_suffix('.metadata.json').write_text(json.dumps(data,indent=2)+'\n')
    (ROOT/'.reference/spinodal_property_differences.json').write_text(json.dumps(differences,indent=2)+'\n')
    print(json.dumps(dict(grid_calls=len(grids),grid_states=len(grid_rows),property_calls=len(roots),failed=len(failed),censored=len(censored),
        coordinate_warnings=len(warnings),critical_roundoff_states=len(conditioned),differences=len(differences),worst_relative_differences=maxima),indent=2))
    if differences:print(json.dumps(differences[:6],indent=2));raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');validate(p.parse_args().compiler_directory)
