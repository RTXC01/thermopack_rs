"""Independent source-backed scalar extrema and raw legacy extrema observations.

The scalar oracle calls the unmodified source EOS/sat_newton with explicit phase
and status arguments. It is NOT the fork's defective optional criconden solver.
"""
import argparse,csv,hashlib,json,math,os,subprocess
from pathlib import Path
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess
from validate_envelope import cases as envelope_cases

def cases():
    yield from envelope_cases()
    for eos in ['PR','SRK']:
        for mixing in ['HV1','HV2','NRTL']:
            for pressure in [1e5,1e6]:
                yield dict(eos=eos,components=['C1','H2O'],feed=[.8,.2],initial_pressure=pressure,mixing=mixing)

def label(q,p):
    config='quantum' if q.get('quantum') else 'translated' if q.get('volume_shift') else {'UNIFAC':'unifac','HV1':'hv1','HV2':'hv2','NRTL':'nrtl'}.get(q.get('mixing'),'classic')
    options=[q['initial_pressure'],q.get('maximum_pressure',1.5e7),q.get('minimum_temperature',2. if q.get('quantum') else 80.),q.get('initial_temperature',0.),q.get('initial_beta',1.),q.get('maximum_step',.25),q.get('maximum_points',1000)]
    return '|'.join(['EXTREMUM',q['eos'],','.join(q['components']),config,','.join(map(repr,q['feed'])),p['kind'],repr(p['point'][2]),','.join(map(repr,p['guess'])),','.join(map(repr,options))])

def validate(directory,reuse_source,allow_known_gaps):
    identity=json.loads((ROOT/'.reference/metadata.json').read_text());requests=list(cases())
    worker_hash=hashlib.sha256((ROOT/'tools/rust_port/oracle_worker.py').read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    cache=ROOT/'.reference/envelope_extrema_source.json'
    if reuse_source:
        saved=json.loads(cache.read_text());assert saved['requests']==requests and saved['oracle_sha256']==worker_hash and saved['reference']==identity
        records,legacy=saved['records'],saved['legacy']
    else:
        records=[];legacy=[]
        with DeadlineReferenceProcess(directory) as oracle:
            for i,q in enumerate(requests):
                records.append(dict(request=q,result=oracle.call(dict(q,mode='ENVELOPE_EXTREMA_REFERENCE'),timeout=12.)))
                if (i+1)%10==0:print(i+1,'independent source curve/reference attempts',flush=True)
        # Legacy criconden may access invalid memory: every observation uses a
        # separate worker, so it cannot poison later reference property calls.
        for q in requests[:28:2]:
            with DeadlineReferenceProcess(directory) as oracle:
                legacy.append(dict(request=q,result=oracle.call(dict(q,mode='ENVELOPE_CRICONDEN'),timeout=8.)))
        cache.write_text(json.dumps(dict(reference=identity,oracle_sha256=worker_hash,requests=requests,records=records,legacy=legacy))+'\n')
    good=[];failed=[];unqualified=[];warnings=[]
    for record in records:
        q=record['request'];r=record['result']
        if not r.get('values'):failed.append(record);continue
        for p in r.get('reference_extrema',[]):
            entry=dict(request=q,result=p)
            if not all(map(math.isfinite,p['point'])) or p['fugacity_residual']>2e-7 or abs(p['scaled_stationarity'])>1e-8 or not p['maximum_verified']:
                unqualified.append(entry);continue
            if max(map(abs,p['pressure_errors']))>2.+1e-7*p['point'][1]:warnings.append(entry)
            good.append(entry)
    subprocess.run(['cargo','build','--release','--example','parity','-j','2'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    rows=subprocess.run([str(exe)],input='\n'.join(label(e['request'],e['result']) for e in good)+'\n',capture_output=True,text=True,check=True).stdout.splitlines()
    assert len(rows)==len(good)
    differences=[];auto_differences=[];fixtures=[];local_matches=0;auto_matches=0;flipped=0
    def compare(expected,encoded,allow_flip):
        if encoded.startswith(('ERROR','NONE')):return [encoded],False
        actual=list(map(float,encoded.split(',')));n=(len(expected)-5)//2;assert len(actual)==len(expected)
        flip=actual[2]!=expected[2]
        if flip:
            if not allow_flip:return ['unexpected phase permutation'],True
            actual=[actual[0],actual[1],1.-actual[2],actual[4],actual[3],*actual[5+n:],*actual[5:5+n]]
        bad=[]
        for i,(a,b) in enumerate(zip(expected,actual)):
            atol=2e-5 if i==0 else 2. if i==1 else 0. if i==2 else 1e-10 if i<5 else 2e-8
            if not math.isfinite(b) or abs(a-b)>atol+2e-6*max(abs(a),abs(b)):bad.append([i,a,b])
        return bad,flip
    for entry,row in zip(good,rows):
        q=entry['request'];p=entry['result'];case=label(q,p);local,auto=row.split('|')
        bad,_=compare(p['point'],local,False);autobad,flip=compare(p['point'],auto,True)
        if bad:differences.append(dict(case=case,differences=bad,native=local))
        else:local_matches+=1
        if autobad:auto_differences.append(dict(case=case,differences=autobad,native=auto))
        else:auto_matches+=1;flipped+=flip
        fixtures.append([case,','.join(map(repr,p['point'])),int(not bad),int(not autobad),int(flip)])
    fixture=ROOT/'tests/fixtures/envelope_extrema.csv'
    with fixture.open('w',newline='') as handle:
        w=csv.writer(handle,delimiter=';',lineterminator='\n');w.writerow(['case','source_backed_reference_state','local_match','automatic_match','automatic_phase_permutation']);w.writerows(fixtures)
    legacy_classification=[]
    for e in legacy:
        q=e['request'];r=e['result'];summary=dict(request=q,returned=bool(r.get('values')),extrema=[])
        for p in r.get('extrema',[]):
            reported=p['reported'];point=p.get('point');reasons=[]
            if min(reported)<=0:reasons.append('unreported coordinate')
            if point:
                n=len(q['feed']);contrast=sum((x-y)**2 for x,y in zip(point[5:5+n],point[5+n:5+2*n]))
                if contrast<1e-12:reasons.append('coincident phases')
                elif abs(p['source_property_stationarity'])/contrast>1e-7:reasons.append('nonstationary source property derivatives')
                if p['physical_fugacity_residual']>2e-7 or p['material_balance_residual']>2e-7:reasons.append('nonclosing coexistence')
                if max(map(abs,point[-3:-1]))>2.+1e-7*point[1]:reasons.append('nonclosing phase pressure')
                if any(abs(x-y)>1e-8*max(abs(x),abs(y)) for x,y in zip(reported,point[:2])):reasons.append('reported coordinate absent from stored curve')
                target=next((v['result']['point'] for v in good if v['request']==q and v['result']['kind']==p['kind']),None)
                matches=target is not None and all(abs(x-y)<=a+2e-6*max(abs(x),abs(y)) for x,y,a in zip(reported,target[:2],[2e-5,2.]))
            else:matches=False
            summary['extrema'].append(dict(kind=p['kind'],reported=reported,reasons=reasons,matches_independent_reference=matches))
        legacy_classification.append(summary)
    metadata=dict(identity,cases=len(fixtures),fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),generator='tools/rust_port/validate_envelope_extrema.py',
        oracle_sha256=worker_hash,curve_attempts=len(requests),reference_failures=failed,
        scalar_reference_errors=[dict(request=e['request'],errors=e['result']['reference_extrema_errors']) for e in records if e['result'].get('reference_extrema_errors')],
        unqualified_reference_returns=unqualified,pressure_warnings=warnings,local_differences=differences,automatic_differences=auto_differences,
        automatic_phase_permutations=flipped,raw_legacy_optional_extrema_observations=legacy,
        legacy_optional_extrema_classification=legacy_classification,
        qualification='Independent scalar bisection of weighted fugacity derivatives using the exact source EOS and sat_newton, with explicit Yphase/ierr arguments. This is not a claim of matching the defective legacy optional envelope criconden path. All raw legacy observations and inner scalar-reference errors remain separate.')
    fixture.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(dict(curves=len(requests),reference_roots=len(good),failed_curves=len(failed),unqualified=len(unqualified),warnings=len(warnings),local_matches=local_matches,automatic_matches=auto_matches,phase_permutations=flipped),indent=2))
    print('Local differences',json.dumps(differences[:8],indent=2));print('Automatic differences',json.dumps(auto_differences[:8],indent=2))
    (ROOT/'.reference/envelope_extrema_differences.json').write_text(json.dumps(dict(local=differences,automatic=auto_differences),indent=2)+'\n')
    if differences or auto_differences:
        known=not auto_differences and len(differences)==1
        if known:
            f=differences[0]['case'].split('|')
            known=(f[1:7]==['PR','C1,C2','classic','0.5,0.5','P','0.0'] and f[8]=='100000.0,15000000.0,80.0,180.0,1.0,0.25,1000'
                and differences[0]['native']=='ERROR:solver did not converge: unresolved envelope extremum phase contrast')
        if not allow_known_gaps or not known:raise SystemExit(1)
        print('One known warm-guess error remains; automatic detection matches this extremum. Not counted as strict local parity.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');p.add_argument('--reuse-source',action='store_true');p.add_argument('--allow-known-gaps',action='store_true');a=p.parse_args();validate(a.compiler_directory,a.reuse_source,a.allow_known_gaps)
