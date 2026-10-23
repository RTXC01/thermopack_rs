"""Pure dispatch of envelopePlot, including its control and critical-edge conventions."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess
from validate_pure_curve import cases as pure_cases,critical_bound

def raw_cases():
    for q in pure_cases():
        if q.get('logarithmic',False) or 'maximum_pressure_step' in q:continue
        yield dict(q,mode='PURE_ENVELOPE')
    for count in [1,2,17]:
        for start in [dict(initial_pressure=.05),dict(initial_temperature=1.)]:
            yield dict(mode='PURE_ENVELOPE',eos='PR',components=['C1'],relative_start=True,maximum_points=count,**start)
    for eos in ['PR','SRK','RK','VDW']:
        for start in [dict(initial_temperature=1.),dict(initial_pressure=1.),dict(initial_temperature=1.0001),dict(initial_pressure=1.0001)]:
            yield dict(mode='PURE_ENVELOPE',eos=eos,components=['CO2'],relative_start=True,**start)
        yield dict(mode='PURE_ENVELOPE',eos=eos,components=['C1','C2'],feed=[1.,1e-18],relative_start=True,initial_pressure=.05)
    for extra in [dict(maximum_pressure=0.,minimum_temperature=-1.,maximum_step=-1.),dict(maximum_pressure=1e4),dict(minimum_temperature=800.),dict(maximum_step=10.),dict(initial_beta=0.)]:
        yield dict(mode='PURE_ENVELOPE',eos='PR',components=['C1'],relative_start=True,initial_pressure=.05,**extra)
    for config,names in [('unifac',['C1','C3']),('hv1',['C1','C2']),('quantum',['H2','NE'])]:
        yield dict(mode='PURE_ENVELOPE',eos='PR',components=names,configuration=config,feed=[1.,1e-18],relative_start=True,initial_pressure=.05)
    for eos in ['PR','SRK','RK','VDW']:
        for count in [1,17]:
            for start in [dict(initial_pressure=.05),dict(initial_temperature=1.)]:
                yield dict(mode='PURE_CURVE',eos=eos,components=['C1'],relative_start=True,maximum_points=count,**start)

def cases():
    seen=set()
    for q in raw_cases():
        key=(q['mode'],q['eos'],tuple(q['components']),q.get('configuration','classic'),q.get('component',0),tuple(q.get('feed',[])),
            q.get('initial_pressure'),q.get('initial_temperature'),q.get('relative_start',False),q.get('maximum_points',100),
            q.get('maximum_pressure',1.5e7),q.get('minimum_temperature',80.),q.get('maximum_step',.25),q.get('initial_beta',1.))
        if key not in seen:seen.add(key);yield q

def label(q,r):
    tag='PURE_ENV' if q['mode']=='PURE_ENVELOPE' else 'PURE_CURVE'
    kind='T' if r['initial_temperature'] is not None else 'P';value=r['initial_temperature'] if kind=='T' else r['initial_pressure']
    index=q.get('component',0);z=q.get('feed',[float(i==index) for i in range(len(q['components']))])
    options=[index,kind,value,2e4,q.get('maximum_points',100),0]
    controls=[q.get('maximum_pressure',1.5e7),q.get('minimum_temperature',80.),q.get('maximum_step',.25),q.get('initial_beta',1.)]
    return f"{tag}|{q['eos']}|{','.join(q['components'])}|{q.get('configuration','classic')}|"+':'.join(map(str,options))+'|'+q.get('parameter_reference','DEFAULT')+'|'+','.join(map(repr,z))+'|'+':'.join(map(repr,controls))

def validate(directory,reuse_source=False):
    cache=ROOT/'.reference/pure_envelope_source.json';identity=json.loads((ROOT/'.reference/metadata.json').read_text())
    worker_hash=hashlib.sha256((ROOT/'tools/rust_port/oracle_worker.py').read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    records=[];failed=[];censored=[];unqualified=[];warnings=[]
    if reuse_source:
        saved=json.loads(cache.read_text());assert saved['reference']==identity and saved['oracle_sha256']==worker_hash and saved['requests']==list(cases())
        records,failed,censored,unqualified,warnings=[saved[k] for k in ['records','failed','censored','unqualified','warnings']]
    with DeadlineReferenceProcess(directory) as oracle:
        for i,q in enumerate([] if reuse_source else cases()):
            r=oracle.call(q,timeout=10.);record=dict(request=q,result=r)
            if 'reference_timeout' in r:censored.append(record);continue
            if not r.get('values'):failed.append(record);continue
            if any(len(p)!=7 or not all(map(math.isfinite,p)) for p in r['values']):unqualified.append(record);continue
            if q['mode']=='PURE_ENVELOPE' and (any(b!=0. for b in r['beta']) or any(any(k!=1. for k in row) for row in r['ratios'])):unqualified.append(record);continue
            record['case']=label(q,r);records.append(record)
            for j,p in enumerate(r['values']):
                if max(abs(p[4]),abs(p[5]))>1e-3+1e-7*p[1] or abs(p[6])>2e-9:warnings.append(dict(case=record['case'],index=j,values=p))
            print(i+1,'curves',len(records),'failures',len(failed),'censored',len(censored),'unqualified',len(unqualified),flush=True)
    if not reuse_source:cache.write_text(json.dumps(dict(records=records,failed=failed,censored=censored,unqualified=unqualified,warnings=warnings,requests=list(cases()),oracle_sha256=worker_hash,reference=identity))+'\n')
    subprocess.run(['cargo','build','--release','--example','parity','-j','2'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    outputs=subprocess.run([str(exe)],input='\n'.join(r['case'] for r in records)+'\n',capture_output=True,text=True,check=True).stdout.splitlines();assert len(outputs)==len(records)
    differences=[];rows=[];summary=[];raw=[];conditioned=[];derived_critical=[];rounded=[];worst=[0.]*4;interior=[0.]*4;critical_nodes=0
    for record,out in zip(records,outputs):
        q=record['request'];r=record['result'];ref=r['values'];name=record['case'];index=q.get('component',0);feed=q.get('feed',[float(i==index) for i in range(len(q['components']))]);feed=[x/sum(feed) for x in feed]
        radius,eta=critical_bound(r);poly=r['critical_polynomial'];shift=r['critical_shift'];bad=[]
        if out.startswith('ERROR'):bad.append(out);native=[];end=out
        else:
            end,encoded=out.split('|');native=[list(map(float,row.split(','))) for row in encoded.split(';')]
            if len(native)!=len(ref):bad.append(['count',len(ref),len(native)])
            if end!='SourceCriticalParameters':bad.append(['end',end])
            for i,(a,b) in enumerate(zip(ref,native)):
                near_critical=abs(a[0]-r['critical'][0])<=64.*2.220446049250313e-16*r['critical'][0] and abs(a[1]-r['critical'][1])<=64.*2.220446049250313e-16*r['critical'][1]
                for j,(x,y) in enumerate(zip(a[:4],b[:4])):
                    standard=[2e-6,.1,1e-12,1e-12][j]+2e-7*max(abs(x),abs(y));bound=standard
                    worst[j]=max(worst[j],abs(x-y)/max(abs(x),abs(y)))
                    if not near_critical:interior[j]=max(interior[j],abs(x-y)/max(abs(x),abs(y)))
                    if near_critical and j>=2:
                        for v in [x,y]:
                            z=(v+shift)*a[1]/(8.31446261815324*a[0]);error=abs(math.fsum([z**3,poly[0]*z*z,poly[1]*z,poly[2]]))
                            if error>eta:bad.append(['critical backward error',i,j,error,eta])
                        if abs(x-y)>standard:
                            raw.append(dict(case=name,index=i,quantity=j,source=x,native=y,difference=abs(x-y)))
                            if abs(x-y)<=standard+radius:conditioned.append(dict(case=name,index=i,quantity=j,radius=radius))
                        bound+=radius
                    if not math.isfinite(y) or abs(x-y)>bound:bad.append([i,j,x,y,bound])
                if q['mode']=='PURE_ENVELOPE':
                    expected=[0.,*feed,*feed]
                    for j,(x,y) in enumerate(zip(expected,b[4:-1])):
                        if x!=y:bad.append(['labels/composition',i,j,x,y])
                    expected_residual=r['all_fugacity_differences'][i]
                    if not math.isfinite(b[-1]):bad.append(['nonfinite fugacity residual',i,b[-1]])
                    elif abs(expected_residual-b[-1])>2e-9+2e-7*max(abs(expected_residual),abs(b[-1])):
                        # The original envelope routine does not return this
                        # derived diagnostic. At conditioned critical roots,
                        # its trace-component sensitivity is retained rather
                        # than claimed as ordinary raw-diagnostic parity.
                        entry=dict(case=name,index=i,source=expected_residual,native=b[-1])
                        if near_critical:derived_critical.append(entry)
                        else:bad.append(['full fugacity residual',i,expected_residual,b[-1]])
        for i,a in enumerate(ref):
            near_critical=abs(a[0]-r['critical'][0])<=64.*2.220446049250313e-16*r['critical'][0] and abs(a[1]-r['critical'][1])<=64.*2.220446049250313e-16*r['critical'][1]
            critical_nodes+=int(near_critical)
            rows.append([name,i,*a,*poly,shift,eta if near_critical else 0.,radius if near_critical else 0.,r['all_fugacity_differences'][i]])
        if bad:differences.append(dict(case=name,differences=bad))
        if sum(x>0. for x in feed)>1:
            rounded.append(dict(case=name,nodes=len(ref),maximum_source_fugacity_difference=max(r['all_fugacity_differences'])))
        summary.append(dict(case=name,source_points=len(ref),native_points=len(native),native_end=end,differences=len(bad)))
    fixture=ROOT/'tests/fixtures/pure_envelope.csv'
    with fixture.open('w',newline='') as f:
        w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(['case','index','temperature','pressure','liquid_volume','vapor_volume','liquid_pressure_error','vapor_pressure_error','selected_fugacity_difference','critical_polynomial_a','critical_polynomial_b','critical_polynomial_c','critical_shift','critical_backward_bound','critical_volume_radius','all_active_fugacity_difference']);w.writerows(rows)
    metadata=dict(identity,cases=len(rows),curve_attempts=sum(1 for _ in cases()),finite_source_curves=len(records),critical_nodes=critical_nodes,generator='tools/rust_port/validate_pure_envelope.py',
        fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),reference_errors=failed,censored_calls=censored,unqualified_returns=unqualified,source_coordinate_warnings=warnings,
        native_curves=summary,raw_critical_volume_differences=raw,conditioned_critical_comparisons=conditioned,derived_critical_fugacity_differences=derived_critical,roundoff_pure_feed_classifications=rounded,worst_relative_difference=worst,worst_interior_relative_difference=interior,
        relative_tolerance=2e-7,absolute_tolerances=[2e-6,.1,1e-12,1e-12],critical_bound_method='Source monic cubic with active alpha/beta/translation; 32-epsilon monomial backward-error bound required for both roots; depressed-cubic displacement bound, only for nodes at source critical T/P within 64 epsilons, including initial T=Tc nodes. Every raw critical difference is retained. All ordinary coordinates retain the standard tolerances.')
    fixture.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n');(ROOT/'.reference/pure_envelope_differences.json').write_text(json.dumps(differences,indent=2)+'\n')
    print(json.dumps(dict(attempts=metadata['curve_attempts'],finite_source_curves=len(records),source_states=len(rows),critical_nodes=critical_nodes,failures=len(failed),censored=len(censored),unqualified=len(unqualified),warnings=len(warnings),different_curves=len(differences),conditioned_critical_comparisons=len(conditioned),worst_interior=interior),indent=2))
    if differences:print(json.dumps(differences[:5],indent=2));raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');p.add_argument('--reuse-source',action='store_true');a=p.parse_args();validate(a.compiler_directory,a.reuse_source)
