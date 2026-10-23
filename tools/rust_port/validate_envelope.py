"""Whole mixture envelopes and projection of independent native curves."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess
RTOL=2e-6

def cases():
    mixtures=[('C1,C2',[.5,.5]),('C1,C3',[.5,.5]),('CO2,C1',[.8,.2]),('CO2,C3',[.5,.5]),
        ('C1,C2,C3',[.2,.3,.5]),('C1,C2,C3',[.5,.5,0.]),('C1,C2,C3',[.9,.09,.01])]
    for eos in ['PR','SRK','RK','VDW']:
        for names,z in mixtures:
            for p in [1e5,1e6]:yield dict(eos=eos,components=names.split(','),feed=z,initial_pressure=p)
    for eos in ['PR','SRK']:
        for names,z in mixtures[:3]+mixtures[4:5]:
            for p in [1e5,1e6]:yield dict(eos=eos,components=names.split(','),feed=z,initial_pressure=p,volume_shift=True)
        for names,z in [mixtures[1],mixtures[4]]:
            for p in [1e5,1e6]:yield dict(eos=eos,components=names.split(','),feed=z,initial_pressure=p,mixing='UNIFAC')
        for extra in [dict(maximum_pressure=2e6),dict(minimum_temperature=150.),dict(maximum_points=3),dict(initial_temperature=180.),dict(initial_temperature=180.,initial_beta=0.)]:
            yield dict(eos=eos,components=['C1','C2'],feed=[.5,.5],initial_pressure=1e5,**extra)
    for eos in ['PR','SRK']:
        yield dict(eos=eos,components=['N2','C1'],feed=[.95,.05],initial_pressure=1e5)
        yield dict(eos=eos,components=['C1','C2'],feed=[.999,.001],initial_pressure=1e5)
    for names in [['H2','NE'],['H2','D2']]:
        for p in [1e5,1e6]:yield dict(eos='PR',components=names,feed=[.5,.5],initial_pressure=p,minimum_temperature=2.,quantum=True)

def label(q):
    config='quantum' if q.get('quantum') else 'translated' if q.get('volume_shift') else 'unifac' if q.get('mixing')=='UNIFAC' else 'classic'
    opts=[q['initial_pressure'],q.get('maximum_pressure',1.5e7),q.get('minimum_temperature',80.),q.get('maximum_step',.25),q.get('maximum_points',1000),q.get('initial_temperature',0.),q.get('initial_beta',1.)]
    return f"ENVELOPE|{q['eos']}|{','.join(q['components'])}|{config}|"+','.join(map(repr,q['feed']))+'|'+':'.join(map(repr,opts))

def ending(q,points):
    if len(points)==q.get('maximum_points',1000):return 'PointLimit'
    t,p=points[-1][:2]
    if abs(t-q.get('minimum_temperature',80.))<1e-4:return 'MinimumTemperature'
    if abs(p-q.get('maximum_pressure',1.5e7))<2.:return 'MaximumPressure'
    if len(points)>1 and abs(p-points[0][1])<2.:return 'InitialPressure'
    return 'NoConvergence'

def validate(directory):
    records=[];failed=[];censored=[];unqualified=[];warnings=[]
    with DeadlineReferenceProcess(directory) as oracle:
        for i,q in enumerate(cases()):
            q=dict(q,mode='ENVELOPE');r=oracle.call(q,timeout=12.);record=dict(case=label(q),request=q,result=r)
            if 'reference_timeout' in r:censored.append(record);continue
            if not r.get('values'):failed.append(record);continue
            points=r['values'];nc=len(q['feed']);width=8+2*nc
            if any(len(p)!=width or not all(map(math.isfinite,p)) or min(p[:2]+p[3:5])<=0 for p in points):unqualified.append(record);continue
            record['source_end']=ending(q,points);records.append(record)
            for j,p in enumerate(points):
                reasons=[]
                if max(abs(p[-3]),abs(p[-2]))>2.+1e-7*p[1]:reasons.append('source volume pressure mismatch')
                if p[-1]>1e-7:reasons.append('source VLE residual')
                if reasons:warnings.append(dict(case=record['case'],index=j,reasons=reasons,values=p))
            print(i+1,'curves',len(records),'failed',len(failed),'censored',len(censored),flush=True)
    subprocess.run(['cargo','build','--release','--example','parity','-j','2'],cwd=ROOT,check=True)
    inputs=[]
    for r in records:
        encoded=';'.join(','.join(map(repr,[p[2],p[0],p[1],*b['log_k'],b['best_coordinate']])) for p,b in zip(r['result']['values'],r['result']['bounds']))
        inputs.append(r['case']+'@'+encoded)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    output=subprocess.run([str(exe)],input='\n'.join(inputs)+'\n',capture_output=True,text=True,check=True).stdout.splitlines()
    assert len(output)==len(records)
    rows=[];differences=[];grid_differences=[];whole=[];native_summary=[];curve_rows=[]
    for r,out in zip(records,output):
        source=r['result']['values'];nc=len(r['request']['feed']);label_=r['case'];width=5+2*nc
        if out.startswith('ERROR'):differences.append(dict(case=label_,rust=out));continue
        end,encoded,projected=out.split('|');native=[list(map(float,p.split(','))) for p in encoded.split(';')];projections=projected.split(';')
        assert len(projections)==len(source)
        native_summary.append(dict(case=label_,end=end,source_end=r['source_end'],native_points=len(native),source_points=len(source)))
        if r['source_end']!='NoConvergence' and end!=r['source_end']:differences.append(dict(case=label_,end=[r['source_end'],end]))
        def compare(a,b,index):
            bad=[]
            for j,(x,y) in enumerate(zip(a[:width],b[:width])):
                absolute=2e-5 if j==0 else 2. if j==1 else 0. if j==2 else 1e-10 if j<5 else 2e-8
                tolerance=absolute+RTOL*max(abs(x),abs(y))
                # Raw grid comparison does not use conditioning radii. Source
                # Jacobian backward errors are retained for diagnostics only.
                if not math.isfinite(y) or abs(x-y)>tolerance:bad.append([index,j,x,y])
            return bad
        grid_bad=[]
        if len(native)==len(source):
            for i,(a,b) in enumerate(zip(source,native)):grid_bad.extend(compare(a,b,i))
        else:grid_bad=[['counts',len(source),len(native)]]
        if grid_bad:grid_differences.append(dict(case=label_,differences=grid_bad[:12]))
        else:whole.append(label_)
        curve_rows.append([label_,len(source),r['source_end'],int(not grid_bad)])
        for i,(a,b) in enumerate(zip(source,projections)):
            bound=r['result']['bounds'][i]
            rows.append([label_,i,*a,*bound['correction'],bound['condition'],bound['best_coordinate'],*bound['log_k']])
            if b.startswith('ERROR'):differences.append(dict(case=label_,index=i,rust=b));continue
            b=list(map(float,b.split(',')))
            bad=compare(a,b,i)
            if bad:differences.append(dict(case=label_,projection=bad,bound=bound))
    fixture=ROOT/'tests/fixtures/envelope.csv'
    # Mixtures have different widths. Encode source state and backward-error
    # vectors in two comma-separated columns within the semicolon CSV.
    with fixture.open('w',newline='') as f:
        w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(['case','index','source_state_and_backward_error'])
        for row in rows:w.writerow([row[0],row[1],','.join(map(repr,row[2:]))])
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text());metadata.update(cases=len(rows),curve_attempts=sum(1 for _ in cases()),generator='tools/rust_port/validate_envelope.py',fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),
        relative_tolerance=RTOL,reference_errors=failed,censored_calls=censored,unqualified_returns=unqualified,source_coordinate_warnings=warnings,
        native_curves=native_summary,matching_adaptive_grids=whole,adaptive_grid_differences=grid_differences,
        source_termination_method='Inferred from returned pressure/temperature limits and point count; NoConvergence denotes an unclassified partial prefix because envelopePlot does not return its control returnCause.',
        source_partial_curves=[r for r in native_summary if r['source_end']=='NoConvergence'])
    fixture.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    curve_fixture=ROOT/'tests/fixtures/envelope_curves.csv'
    with curve_fixture.open('w',newline='') as f:
        w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(['case','source_count','source_end','matching_grid']);w.writerows(curve_rows)
    curve_metadata=dict(metadata,cases=len(curve_rows),fixture_sha256=hashlib.sha256(curve_fixture.read_bytes()).hexdigest())
    curve_fixture.with_suffix('.metadata.json').write_text(json.dumps(curve_metadata,indent=2)+'\n')
    (ROOT/'.reference/envelope_differences.json').write_text(json.dumps(differences,indent=2)+'\n')
    print(json.dumps(dict(curve_attempts=metadata['curve_attempts'],finite_source_curves=len(records),source_states=len(rows),failed=len(failed),censored=len(censored),unqualified=len(unqualified),warnings=len(warnings),matching_grids=len(whole),grid_differences=len(grid_differences),differences=len(differences)),indent=2))
    if differences:print(json.dumps(differences[:8],indent=2));raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');validate(p.parse_args().compiler_directory)
