"""H/S constraints and the exact-source caloric quality Jacobian."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess

def cases():
    mixes=[(['C1','C2'],[.5,.5]),(['C1','C3'],[.5,.5]),(['CO2','C1'],[.8,.2]),(['C1','C2','C3'],[.2,.3,.5]),(['C1','C2','C3'],[.5,.5,0.])]
    for eos in ['PR','SRK','RK','VDW']:
        for names,z in mixes:
            for beta in [0.,.25,.5,.75,1.]:
                for prop in ['H','S']:
                    yield dict(eos=eos,components=names,feed=z,beta=beta,property=prop,base_pressure=1e6)
    for eos in ['PR','SRK']:
        for config in ['translated','unifac']:
            for names,z in [mixes[1],mixes[3]]:
                for beta in [.25,.5,.75]:
                    for prop in ['H','S']:
                        yield dict(eos=eos,components=names,feed=z,beta=beta,property=prop,base_pressure=1e6,configuration=config)
        for config in ['hv1','hv2','nrtl']:
            for beta in [.25,.75]:
                for prop in ['H','S']:
                    yield dict(eos=eos,components=['C1','H2O'],feed=[.8,.2],beta=beta,property=prop,base_pressure=1e6,configuration=config)
    for names in [['H2','NE'],['H2','D2']]:
        for beta in [.25,.5,.75]:
            for prop in ['H','S']:
                yield dict(eos='PR',components=names,feed=[.5,.5],beta=beta,property=prop,base_pressure=1e6,configuration='quantum')
    for prop,offset in [('H',500.),('H',-500.),('S',1.),('S',-1.)]:
        for beta in [.25,.75]:
            yield dict(eos='PR',components=['C1','C2'],feed=[.5,.5],beta=beta,property=prop,base_pressure=1e6,offset=offset)

def label(q,r):
    g=r['guess'];guess=[g['temperature'],g['pressure'],*g['ratios']]
    return '|'.join(['QCAL',q['eos'],','.join(q['components']),q.get('configuration','classic'),','.join(map(repr,q['feed'])),repr(q['beta']),
        f"{q['property']}:{r['target']!r}:0",','.join(map(repr,guess))])

def validate(directory,reuse_source):
    requests=[dict(q,mode='QUALITY_CALORIC') for q in cases()];identity=json.loads((ROOT/'.reference/metadata.json').read_text())
    worker_hash=hashlib.sha256((ROOT/'tools/rust_port/oracle_worker.py').read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    cache=ROOT/'.reference/quality_caloric_source.json'
    if reuse_source:
        saved=json.loads(cache.read_text());assert saved['reference']==identity and saved['oracle_sha256']==worker_hash and saved['requests']==requests
        records=saved['records']
    else:
        records=[]
        with DeadlineReferenceProcess(directory) as oracle:
            for i,q in enumerate(requests):
                records.append(dict(request=q,result=oracle.call(q,timeout=8.)))
                if (i+1)%20==0:print(i+1,'source attempts',flush=True)
        cache.write_text(json.dumps(dict(reference=identity,oracle_sha256=worker_hash,requests=requests,records=records))+'\n')
    good=[];failed=[];censored=[];unqualified=[];warnings=[]
    jacobian_rows=[]
    for e in records:
        q=e['request'];r=e['result'];nc=len(q['feed'])
        if r.get('initial_caloric_jacobian_row') and all(map(math.isfinite,r['initial_caloric_jacobian_row']+r['initial_function'])):
            jacobian_rows.append([label(q,r),','.join(map(repr,r['initial_caloric_jacobian_row'])),','.join(map(repr,r['initial_function'])),','.join(map(repr,r['initial_caloric_finite_difference_row']))])
        if 'reference_timeout' in r:censored.append(e);continue
        if r.get('status')!=0 or not r.get('values'):failed.append(e);continue
        if len(r['values'])!=8+2*nc or not all(map(math.isfinite,r['values'])):
            unqualified.append(e);continue
        if r['fugacity_residual']>2e-7 or r['material_balance_residual']>2e-7 or abs(r['caloric_residual'])>2e-8:
            unqualified.append(e);continue
        good.append(e)
        if max(map(abs,r['pressure_errors']))>1e-3+1e-7*r['values'][1]:warnings.append(e)
    subprocess.run(['cargo','build','--release','--example','parity','-j','2'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    reachable=[e for e in records if 'guess' in e['result']]
    output=subprocess.run([str(exe)],input='\n'.join(label(e['request'],e['result']) for e in reachable)+'\n',capture_output=True,text=True,check=True).stdout.splitlines()
    assert len(output)==len(reachable)
    native_by_label={label(e['request'],e['result']):line for e,line in zip(reachable,output)}
    native=[native_by_label[label(e['request'],e['result'])] for e in good]
    failure_outcomes=[]
    for e in failed:
        if 'guess' not in e['result']:continue
        key=label(e['request'],e['result']);baseline,physical=native_by_label[key].split('|')
        failure_outcomes.append(dict(case=key,source_status=e['result']['status'],baseline_rejected=baseline.startswith('ERROR'),
            physical_solved=not physical.startswith('ERROR'),baseline=baseline,physical=physical))
    by_label={label(e['request'],e['result']):e for e in reachable}
    with DeadlineReferenceProcess(directory) as oracle:
        for result in failure_outcomes:
            if not result['physical_solved']:continue
            e=by_label[result['case']];q=e['request'];r=e['result'];n=len(q['feed']);v=list(map(float,result['physical'].split(',')))
            point=dict(temperature=v[0],pressure=v[1],liquid_volume=v[3],vapor_volume=v[4],liquid=v[5:5+n],vapor=v[5+n:5+2*n])
            checked=oracle.call(dict(q,check_state=point,target=r['target']),timeout=8.)
            result['independent_source_check']=checked
            assert checked.get('fugacity_residual',math.inf)<2e-7 and checked.get('material_balance_residual',math.inf)<2e-7,('independent phase closure',result)
            assert abs(checked['caloric_residual'])<2e-8 and max(map(abs,checked['pressure_errors']))<1e-3+1e-7*v[1],('independent caloric/pressure closure',result)
    rows=[];differences=[];physical_differences=[];maxima=[0.]*4
    for e,line in zip(good,native):
        q=e['request'];r=e['result'];key=label(q,r);baseline,physical=line.split('|');expected=r['values'];nc=len(q['feed'])
        def compare(encoded,physical):
            if encoded.startswith('ERROR'):return [encoded]
            actual=list(map(float,encoded.split(',')));assert len(actual)==len(expected)
            bad=[]
            for j,(a,b) in enumerate(zip(expected,actual)):
                atol=2e-5 if j==0 else 2. if j==1 else 0. if j==2 else 1e-10 if j<5 else 1e-3 if j==5+2*nc else 1e-5 if j==6+2*nc else 1e-8 if j==7+2*nc else 2e-8
                rtol=2e-6 if physical or j==5+2*nc else 2e-7
                if not math.isfinite(b) or abs(a-b)>atol+rtol*max(abs(a),abs(b)):bad.append([j,a,b])
                if j in [0,1,3,4]:maxima[[0,1,3,4].index(j)]=max(maxima[[0,1,3,4].index(j)],abs(a-b)/max(abs(a),abs(b)))
            return bad
        a=compare(baseline,False);b=compare(physical,True)
        assert not physical.startswith('ERROR'),('physical caloric solve',key,physical)
        if not physical.startswith('ERROR'):
            actual=list(map(float,physical.split(',')));v=actual[-2] if q['property']=='H' else actual[-1]
            assert abs(v-r['target'])<=1.1e-9*max(abs(r['target']),1.),('physical target closure',key,v,r['target'])
        if a:differences.append(dict(case=key,differences=a))
        if b:physical_differences.append(dict(case=key,differences=b))
        rows.append([key,','.join(map(repr,expected)),','.join(map(repr,r['initial_caloric_jacobian_row'])),','.join(map(repr,r['initial_function'])),r['fugacity_residual'],r['material_balance_residual'],r['caloric_residual']])
    physical_difference_source_checks=[]
    with DeadlineReferenceProcess(directory) as oracle:
        for entry in physical_differences:
            e=by_label[entry['case']];q=e['request'];r=e['result'];n=len(q['feed'])
            v=list(map(float,native_by_label[entry['case']].split('|')[1].split(',')))
            point=dict(temperature=v[0],pressure=v[1],liquid_volume=v[3],vapor_volume=v[4],liquid=v[5:5+n],vapor=v[5+n:5+2*n])
            checked=oracle.call(dict(q,check_state=point,target=r['target']),timeout=8.)
            physical_difference_source_checks.append(dict(case=entry['case'],native_state=point,source_check=checked))
            assert checked.get('fugacity_residual',math.inf)<2e-7 and checked.get('material_balance_residual',math.inf)<2e-7,('physical difference phase closure',entry,checked)
            assert abs(checked['caloric_residual'])<2e-8 and max(map(abs,checked['pressure_errors']))<1e-3+1e-7*v[1],('physical difference caloric/pressure closure',entry,checked)
    fixture=ROOT/'tests/fixtures/quality_caloric.csv'
    with fixture.open('w',newline='') as f:
        w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(['case','source_state','initial_caloric_jacobian_row','initial_function','fugacity_residual','material_balance_residual','caloric_residual']);w.writerows(rows)
    metadata=dict(identity,cases=len(rows),attempts=len(requests),fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),generator='tools/rust_port/validate_quality_caloric.py',
        reference_errors=failed,censored_calls=censored,unqualified_returns=unqualified,source_pressure_warnings=warnings,
        baseline_differences=differences,physical_differences=physical_differences,worst_state_relative_differences=maxima,
        failed_source_case_native_outcomes=failure_outcomes,
        physical_difference_source_checks=physical_difference_source_checks,
        baseline_relative_tolerance=2e-7,physical_relative_tolerance=2e-6,slope_relative_tolerance=2e-6,
        qualification='Finite closing source sat_newton H/S returns from independently generated fixed-P quality states, perturbed T/P/K guesses and selected shifted H/S targets. Complete failed-call/iteration identity and global stability remain separate.')
    fixture.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    jac_fixture=ROOT/'tests/fixtures/quality_caloric_jacobian.csv'
    with jac_fixture.open('w',newline='') as f:
        w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(['case','initial_caloric_jacobian_row','initial_function','initial_caloric_finite_difference_row']);w.writerows(jacobian_rows)
    jac_meta=dict(identity,cases=len(jacobian_rows),fixture_sha256=hashlib.sha256(jac_fixture.read_bytes()).hexdigest(),generator='tools/rust_port/validate_quality_caloric.py',
        earlier_initialization_failures=[e for e in failed if 'guess' not in e['result']],
        qualification='Source last-row H/S derivatives and full residual at every finite initial guess, including calls whose caloric solve later fails. These are function/derivative comparisons, not successful solve claims.')
    jac_fixture.with_suffix('.metadata.json').write_text(json.dumps(jac_meta,indent=2)+'\n')
    (ROOT/'.reference/quality_caloric_differences.json').write_text(json.dumps(dict(baseline=differences,physical=physical_differences),indent=2)+'\n')
    print(json.dumps(dict(attempts=len(requests),finite_closing=len(good),failures=len(failed),censored=len(censored),unqualified=len(unqualified),warnings=len(warnings),baseline_differences=len(differences),physical_differences=len(physical_differences),worst=maxima),indent=2))
    if differences:print(json.dumps(differences[:5],indent=2))
    if physical_differences:print(json.dumps(physical_differences[:5],indent=2))
    expected_physical={
        ('RK','C1,C3','classic','0.75','S'),
        ('PR','C1,C2,C3','unifac','0.75','S'),
        ('SRK','C1,H2O','hv1','0.25','H'),
    }
    unexpected=[]
    for entry in physical_differences:
        f=entry['case'].split('|');signature=(f[1],f[2],f[3],f[5],f[6].split(':')[0])
        indexes={int(d[0]) for d in entry['differences']}
        if signature not in expected_physical:unexpected.append(entry);continue
        if signature[0]=='RK':
            if not indexes.issubset({4,9}) or any(abs(a-b)>5e-6*max(abs(a),abs(b)) for _,a,b in entry['differences']):unexpected.append(entry)
        elif signature[0]=='PR':
            if indexes!={4} or any(abs(a-b)>3e-6*max(abs(a),abs(b)) for _,a,b in entry['differences']):unexpected.append(entry)
        else:
            if not indexes.issubset({0,1,3,4,9,11}):unexpected.append(entry)
            v=list(map(float,native_by_label[entry['case']].split('|')[1].split(',')))
            if abs(v[0]-135.9567589314)>2e-5 or abs(v[1]-1e6)>2.:unexpected.append(entry)
    if differences or unexpected:raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');p.add_argument('--reuse-source',action='store_true');a=p.parse_args();validate(a.compiler_directory,a.reuse_source)
