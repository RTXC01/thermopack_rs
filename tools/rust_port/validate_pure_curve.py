"""Selected-component pure saturation grids from the exact-source library."""
import argparse,csv,hashlib,json,math,os,subprocess,sys
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess

def cases():
    starts=[dict(initial_pressure=.05),dict(initial_temperature=.65)]
    for eos in ['PR','SRK','RK','VDW']:
        for name in ['C1','CO2','C3','N2','AR','H2O']:
            for start in starts:
                for logarithmic in [False,True]:
                    yield dict(eos=eos,components=[name],relative_start=True,logarithmic=logarithmic,**start)
    for eos in ['PR','SRK']:
        for name in ['C1','CO2','H2O']:
            for start in starts:
                for logarithmic in [False,True]:
                    yield dict(eos=eos,components=[name],configuration='translated',relative_start=True,logarithmic=logarithmic,**start)
        for config in ['classic','unifac','hv1','hv2','nrtl']:
            for component in range(3):
                for start in starts:
                    yield dict(eos=eos,components=['C1','C2','C3'],component=component,configuration=config,relative_start=True,**start)
    for name in ['H2','D2','HE','NE']:
        for start in starts:
            for logarithmic in [False,True]:
                yield dict(eos='PR',components=[name],configuration='quantum',relative_start=True,logarithmic=logarithmic,**start)
    for count in [2,3,17,500]:
        yield dict(eos='PR',components=['C1'],initial_pressure=.05,relative_start=True,maximum_points=count)
    for step in [1e6,1e7,2.1e6]:
        yield dict(eos='PR',components=['C1'],initial_pressure=.05,relative_start=True,maximum_pressure_step=step)

def label(q,r):
    kind='T' if r['initial_temperature'] is not None else 'P'
    value=r['initial_temperature'] if kind=='T' else r['initial_pressure']
    options=[q.get('component',0),kind,value,q.get('maximum_pressure_step',2e4),q.get('maximum_points',100),int(q.get('logarithmic',False))]
    return f"PURE_CURVE|{q['eos']}|{','.join(q['components'])}|{q.get('configuration','classic')}|"+':'.join(map(str,options))+'|'+q.get('parameter_reference','DEFAULT')

def critical_bound(result):
    a,b,c=result['critical_polynomial'];center=-a/3.
    p=b-a*a/3.;q=2.*a*a*a/27.-a*b/3.+c
    scale=abs(center)**3+abs(a)*center**2+abs(b*center)+abs(c)
    eta=32.*sys.float_info.epsilon*scale
    low=0.;high=1.
    for _ in range(100):
        r=(low+high)*.5
        if r**3-abs(p)*r-abs(q)-eta>0.:high=r
        else:low=r
    t,pressure,_=result['critical']
    return 2.*high*8.31446261815324*t/pressure,eta

def validate(directory,reuse_source=False):
    cache=ROOT/'.reference/pure_curve_source.json'
    oracle_hash=hashlib.sha256((ROOT/'tools/rust_port/oracle_worker.py').read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    identity=json.loads((ROOT/'.reference/metadata.json').read_text())
    records=[];failed=[];censored=[];nonfinite=[];warnings=[]
    if reuse_source:
        saved=json.loads(cache.read_text());assert saved['oracle_sha256']==oracle_hash and saved['reference']==identity
        assert saved['requests']==list(cases())
        records,failed,censored,nonfinite,warnings=[saved[k] for k in ['records','failed','censored','nonfinite','warnings']]
    with DeadlineReferenceProcess(directory) as oracle:
        for index,q in enumerate([] if reuse_source else cases()):
            # Every request constructs a fresh model; only the worker is reused.
            r=oracle.call(dict(q,mode='PURE_CURVE'),timeout=15.)
            record=dict(request=q,result=r)
            if 'reference_timeout' in r:censored.append(record);continue
            if not r.get('values'):failed.append(record);continue
            if any(len(p)!=7 or not all(map(math.isfinite,p)) for p in r['values']):nonfinite.append(record);continue
            record['case']=label(q,r);records.append(record)
            for j,p in enumerate(r['values']):
                if max(abs(p[4]),abs(p[5]))>1e-3+1e-7*p[1] or abs(p[6])>2e-9 or p[2]>=p[3] and j<len(r['values'])-1:
                    warnings.append(dict(case=record['case'],index=j,values=p))
            print(index+1,'curves',len(records),'failed',len(failed),'nonfinite',len(nonfinite),'censored',len(censored),flush=True)
    if not reuse_source:cache.write_text(json.dumps(dict(records=records,failed=failed,censored=censored,nonfinite=nonfinite,warnings=warnings,requests=list(cases()),oracle_sha256=oracle_hash,reference=identity))+'\n')
    subprocess.run(['cargo','build','--release','--example','parity','-j','2'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    outputs=subprocess.run([str(exe)],input='\n'.join(r['case'] for r in records)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(outputs)==len(records)
    differences=[];rows=[];summary=[];worst=[0.]*4;interior_worst=[0.]*4;conditioned=[];raw_endpoint=[]
    for record,out in zip(records,outputs):
        result=record['result'];ref=result['values'];name=record['case'];bad=[]
        radius,eta=critical_bound(result);poly=result['critical_polynomial'];shift=result['critical_shift']
        if out.startswith('ERROR'):bad.append(out);native=[];end=out
        else:
            end,encoded=out.split('|');native=[list(map(float,p.split(','))) for p in encoded.split(';')]
            if len(native)!=len(ref):bad.append(['count',len(ref),len(native)])
            if end!='SourceCriticalParameters':bad.append(['end',end])
            for i,(a,b) in enumerate(zip(ref,native)):
                for j,(x,y) in enumerate(zip(a,b)):
                    if j>=4:break
                    atol=[2e-6,.1,1e-12,1e-12][j];scale=max(abs(x),abs(y))
                    worst[j]=max(worst[j],abs(x-y)/scale)
                    if i<len(ref)-1:interior_worst[j]=max(interior_worst[j],abs(x-y)/scale)
                    error=abs(x-y);standard=atol+2e-7*scale
                    if i==len(ref)-1 and j in [2,3]:
                        for v in [x,y]:
                            z=(v+shift)*a[1]/(8.31446261815324*a[0]);backward=abs(math.fsum([z**3,poly[0]*z*z,poly[1]*z,poly[2]]))
                            if backward>eta:bad.append(['critical backward error',i,j,backward,eta])
                        if error>standard:
                            raw_endpoint.append(dict(case=name,quantity=j,source=x,native=y,difference=error))
                            if error<=standard+radius:conditioned.append(dict(case=name,quantity=j,difference=error,radius=radius))
                        standard+=radius
                    if not math.isfinite(y) or error>standard:bad.append([i,j,x,y])
        if bad:differences.append(dict(case=name,differences=bad))
        summary.append(dict(case=name,source_points=len(ref),native_points=len(native),native_end=end,differences=len(bad)))
        rows.extend([name,i,*p,*(poly+[shift,eta,radius] if i==len(ref)-1 else [0.]*6)] for i,p in enumerate(ref))
    fixture=ROOT/'tests/fixtures/pure_curve.csv'
    with fixture.open('w',newline='') as f:
        w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(['case','index','temperature','pressure','liquid_volume','vapor_volume','liquid_pressure_error','vapor_pressure_error','fugacity_difference','critical_polynomial_a','critical_polynomial_b','critical_polynomial_c','critical_shift','critical_backward_bound','critical_volume_radius']);w.writerows(rows)
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text());metadata.update(cases=len(rows),curve_attempts=sum(1 for _ in cases()),finite_source_curves=len(records),
        generator='tools/rust_port/validate_pure_curve.py',fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),reference_errors=failed,censored_calls=censored,
        nonfinite_returns=nonfinite,source_coordinate_warnings=warnings,native_curves=summary,worst_relative_difference=worst,relative_tolerance=2e-7,absolute_tolerances=[2e-6,.1,1e-12,1e-12],
        worst_interior_relative_difference=interior_worst,
        critical_conditioned_comparisons=conditioned,raw_critical_volume_differences=raw_endpoint,
        critical_bound_method='Only final volumes: source monic cubic z^3+a*z^2+b*z+c, center=-a/3, depressed p/q. eta=32*epsilon*sum(abs(monomials at center)). Both source and native polynomial backward errors must be <= eta. r is the largest positive root of r^3-abs(p)*r-abs(q)-eta=0. The two roots differ by <=2*r, converted to m3/mol with R*T/P; all interior outputs and all T/P keep standard tolerances.')
    fixture.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/pure_curve_differences.json').write_text(json.dumps(differences,indent=2)+'\n')
    print(json.dumps(dict(attempts=metadata['curve_attempts'],source_curves=len(records),source_points=len(rows),failed=len(failed),censored=len(censored),nonfinite=len(nonfinite),warnings=len(warnings),different_curves=len(differences),conditioned_critical_volumes=len(conditioned),worst_relative_difference=worst),indent=2))
    if differences:print(json.dumps(differences[:5],indent=2));raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');p.add_argument('--reuse-source',action='store_true');args=p.parse_args();validate(args.compiler_directory,args.reuse_source)
