"""Original saturation outcomes and independently checked MBWR coexistence states."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess

from reference import ROOT
from mbwr_coexistence_cases import models, reference_requests, physical_requests, native_key, source_input_hashes
from mbwr_coexistence_oracle import build_reference, collect_cached, encode_request, verify_reference
from qualify_mbwr_caloric import native_inputs, atomic_binary
from validate_meos_flash import atomic_csv, atomic_text, digest

RELATIVE=2e-6
ABSOLUTE=[1e-7,1e-4,1e-11,1e-10,1e-4]


def parse_original(observation):
    text=observation['stdout'];out={}
    for label in ('META','STATE','CHECK'):
        matches=re.findall(r'^'+label+r'=([^\r\n]+)',text,re.M)
        assert len(matches)<=1,(label,matches)
        if matches:out[label]=list(map(float,matches[0].split()))
    point=re.search(r'POINT_STATUS=(-?\d+) TP=([^\r\n]+)',text)
    if point:out['status']=int(point[1]);out['TP']=list(map(float,point[2].split()))
    out['ROOT']=[list(map(float,s.split())) for s in re.findall(r'^ROOT=([^\r\n]+)',text,re.M)]
    return out


def equilibrium_checks(state,check):
    if len(state)!=5 or len(check)!=12 or not all(map(math.isfinite,state+check)):
        return ['missing/nonfinite original state']
    t,p,vl,vv,slope=state
    if min(t,p,vl,vv)<=0 or vv<=vl*(1+1e-6):return ['coincident/invalid phase volumes']
    gap=math.log(vv/vl)
    failures=[]
    for i,bound in ((0,check[10]),(1,check[11])):
        if bound<0 or abs(check[i]-p)>2e-7*p+bound:failures.append('original pressure closure')
    if min(check[2:4])<=0:failures.append('original mechanical stability')
    tolerance=2e-8*gap+64*math.ulp(1.)*max(1.,abs(check[6]),abs(check[7]))
    if abs(check[4]-check[5])>tolerance:failures.append('original fugacity equality')
    if abs(check[6]-check[7])>tolerance:failures.append('original chemical potential equality')
    expected=(-8.31446261815324*t*t*(check[9]-check[8])+p*(vv-vl))/(t*(vv-vl))
    if abs(slope-expected)>1e-8+2e-12*max(abs(slope),abs(expected)):
        failures.append('original Clapeyron identity')
    return failures


def source_category(observation,parsed):
    if observation.get('censored'):return 'censored'
    if 'INIT_RETURNED' not in observation['stdout']:return 'initialization_failed'
    if 'status' not in parsed:return 'calculation_fatal'
    if parsed['status']!=0:return 'nonconverged'
    state=parsed.get('STATE',[])
    if len(state)!=5 or not all(map(math.isfinite,state[:4])):return 'nonfinite_state'
    if min(state[:4])<=0:return 'invalid_phase_state'
    if state[3]<=state[2]*(1+1e-6):return 'coincident_or_swapped'
    if not math.isfinite(state[4]):return 'nonfinite_derived_slope'
    return 'kernel_candidate' if not equilibrium_checks(parsed['STATE'],parsed.get('CHECK',[])) else 'unclosed_or_unstable'


def compare_state(expected,actual):
    if not isinstance(actual,list) or len(actual)!=5 or not all(map(math.isfinite,actual)):
        return ['invalid native state']
    assert len(expected)==5 and all(map(math.isfinite,expected))
    return [(i,a,b) for i,(a,b,floor) in enumerate(zip(actual,expected,ABSOLUTE))
            if abs(a-b)>floor+RELATIVE*max(abs(a),abs(b))]


def run_native(exe,requests,env):
    lines=subprocess.check_output([str(exe)],input='\n'.join(map(native_key,requests))+'\n',text=True,env=env).splitlines()
    assert len(lines)==len(requests)
    out=[]
    for line in lines:
        if line.startswith('ERROR'):out.append(dict(error=line,stdout=line));continue
        a,b=line.split(';');state=list(map(float,a.split(',')));roots=list(map(float,b.split(','))) if b else []
        assert len(state)==5 and all(map(math.isfinite,state+roots))
        assert min(state[:4])>0 and state[3]>state[2]
        out.append(dict(state=state,roots=roots,stdout=line))
    return out


def reference_roots(exe,requests,observations,env):
    queries=[]
    for q,p in zip(requests,observations):
        t,pressure,*_=p['STATE'];maximum=4*q['critical'][2]
        queries.append('|'.join(map(str,[q['record'],'R',t,maximum,pressure,0,0,0])))
    if not queries:return []
    lines=subprocess.check_output([str(exe)],input='\n'.join(queries)+'\n',text=True,env=env).splitlines()
    assert len(lines)==len(queries)
    return [dict(request=query,stdout=line,roots=line if line.startswith('ERROR') else list(map(float,line.split(';')[1].split(','))) if line.split(';')[1] else []) for query,line in zip(queries,lines)]


def kernel_request(model,state,roots):
    q=dict(model,mode=5,refine=0,T=state[0],P=state[1],liquid_volume=state[2],vapor_volume=state[3],roots=roots)
    payload=encode_request(q)
    q['id']='kernel:'+hashlib.sha256(payload.encode()).hexdigest()
    return q


def screen_kernel(q,parsed):
    state=parsed.get('STATE',[]);check=parsed.get('CHECK',[])
    failures=equilibrium_checks(state,check)
    if failures:return failures
    roots=parsed['ROOT']
    if len(roots)!=len(q['roots']):return ['original root count']
    reference=.5*(check[6]+check[7]);t,p,vl,vv,_=state
    for expected,row in zip(q['roots'],roots):
        if len(row)!=6 or not all(map(math.isfinite,row)):return ['invalid original root check']
        rho,pressure,prho,phi,mu,bound=row
        if rho!=expected or abs(pressure-p)>2e-7*p+bound:return ['original root pressure closure']
        if prho < -2e-8*8.31446261815324*t:return ['original root mechanical stability']
        if mu < reference-2e-8*max(1.,abs(reference),abs(mu)):return ['lower original root Gibbs energy']
    return []


def clean(value):
    if isinstance(value,float) and not math.isfinite(value):return str(value)
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
    if isinstance(value,list):return list(map(clean,value))
    return value


def main(args):
    exe,env,signature=build_reference(Path(args.compiler),Path(args.source_build_dir))
    signature['cohort_inputs_sha256']=source_input_hashes()
    requests=reference_requests()
    historic,observations=collect_cached(exe,env,signature,requests,ROOT/'.reference/mbwr_coexistence_source.json',args.reuse_source)
    assert not any(r.get('censored') for r in observations),'censored observations retained; no automatic retry'
    parsed=[parse_original(r) for r in observations]
    categories=[source_category(o,p) for o,p in zip(observations,parsed)]
    print(json.dumps(dict(original_requests=len(requests),outcomes=dict(Counter(categories))),indent=2),flush=True)
    if args.source_only:return
    hashes=native_inputs()
    subprocess.run(['cargo','build','--release','--example','mbwr_coexistence_parity','-j','2'],cwd=ROOT,env=env,check=True)
    native=ROOT/'target/release/examples'/('mbwr_coexistence_parity.exe' if os.name=='nt' else 'mbwr_coexistence_parity')
    native_sha=hashlib.sha256(native.read_bytes()).hexdigest()
    physical=physical_requests();all_requests=physical+requests
    native_values=run_native(native,all_requests,env)
    proposed=[];kernel_map={};native_kernel_keys={}
    for q,result in zip(all_requests,native_values):
        if 'state' in result:
            k=kernel_request(q,result['state'],result['roots'])
            kernel_map.setdefault(k['id'],k);native_kernel_keys[q['id']]=k['id']
    candidates=[i for i,c in enumerate(categories) if c=='kernel_candidate']
    candidate_requests=[requests[i] for i in candidates]
    candidate_parsed=[parsed[i] for i in candidates]
    roots=reference_roots(native,candidate_requests,candidate_parsed,env)
    reference_kernel_keys={};root_queries={}
    for i,proof in zip(candidates,roots):
        root_queries[requests[i]['id']]=proof
        r=proof['roots']
        if isinstance(r,str):continue
        k=kernel_request(requests[i],parsed[i]['STATE'],r)
        kernel_map.setdefault(k['id'],k);reference_kernel_keys[requests[i]['id']]=k['id']
    kernels=list(kernel_map.values())
    kernel_signature,kernel_observations=collect_cached(exe,env,signature,kernels,ROOT/'.reference/mbwr_coexistence_kernels.json',args.reuse_kernels)
    assert not any(r.get('censored') for r in kernel_observations),'censored kernel checks retained; no automatic retry'
    kernel_parsed={q['id']:parse_original(o) for q,o in zip(kernels,kernel_observations)}
    errors=[];records=[];rows=[];counts=Counter()
    for q,result in zip(all_requests,native_values):
        k=native_kernel_keys.get(q['id'])
        if k:
            failures=screen_kernel(kernel_map[k],kernel_parsed[k])
            bad_slope=compare_state(kernel_parsed[k]['STATE'],result['state'])
            if failures or bad_slope:errors.append(dict(case=q['id'],kind='native original-kernel check',failures=failures,differences=bad_slope))
            counts['native_kernel_checked']+=1
        else:counts['native_error_observations']+=1
        if q['context']=='physical':
            records.append(dict(request=q,native=result,kernel=k,classification='kernel_checked_native' if k else 'native_error_unqualified'))
            rows.append([q['id'],native_key(q),'',','.join(map(repr,result['state'])) if k else result['error'],'kernel_checked_native' if k else 'native_error_unqualified'])
    for i,(q,observation,p,category) in enumerate(zip(requests,observations,parsed,categories)):
        result=native_values[len(physical)+i];k=reference_kernel_keys.get(q['id']);classification=category;difference=[]
        if category=='kernel_candidate':
            if k is None:classification='reference_root_query_failed';errors.append(dict(case=q['id'],kind=classification))
            else:
                failures=screen_kernel(kernel_map[k],kernel_parsed[k])
                if failures==['lower original root Gibbs energy']:classification='original_metastable_in_domain'
                elif failures:
                    classification='original_reference_kernel_inconsistent'
                    errors.append(dict(case=q['id'],kind=classification,failures=failures))
                else:
                    maximum=4*q['critical'][2];tc=q['critical'][0]
                    if 1/p['STATE'][2]>maximum or 1/p['STATE'][3]>maximum:classification='original_outside_density_domain'
                    elif q['operation']=='T' and not .5*tc<=p['STATE'][0]<=.9999*tc:classification='original_outside_temperature_domain'
                    else:
                        classification='original_solver_match'
                        difference=compare_state(p['STATE'],result.get('state'))
                        if difference:errors.append(dict(case=q['id'],kind='original solver parity',differences=difference,native=result,original=p['STATE']))
        counts[classification]+=1
        records.append(dict(request=q,source=observation,original=clean(p),classification=classification,native=result,
                            native_kernel=native_kernel_keys.get(q['id']),reference_kernel=k,root_query=root_queries.get(q['id'])))
        reference=','.join(map(repr,p['STATE'])) if classification=='original_solver_match' else ''
        rows.append([q['id'],native_key(q),reference,','.join(map(repr,result['state'])) if 'state' in result else result['error'],classification])
    report=dict(cases=len(all_requests),reference_requests=requests,physical_requests=physical,signature=historic,
        kernel_signature=kernel_signature,kernel_requests=kernels,kernel_processes=kernel_observations,
        native_sources_sha256=hashes,native_executable_sha256=native_sha,generator_sha256=digest(Path(__file__)),
        records=records,errors=errors,counts=dict(counts),relative_tolerance=RELATIVE,absolute_tolerances=ABSOLUTE)
    atomic_text(ROOT/'.reference/mbwr_coexistence_report.json',json.dumps(clean(report),indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(cases=len(all_requests),counts=dict(counts),errors=len(errors)),indent=2),flush=True)
    assert not errors,errors[:6]
    assert hashes==native_inputs() and native_sha==hashlib.sha256(native.read_bytes()).hexdigest()
    verify_reference(signature)
    path=ROOT/'tests/fixtures/mbwr_coexistence.csv'
    atomic_csv(path,['case','request','original_solver','native_observed','classification'],rows)
    report['fixture_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    archive=path.with_suffix('.audit.json.gz')
    atomic_binary(archive,gzip.compress(json.dumps(clean(report),separators=(',',':'),allow_nan=False).encode(),mtime=0))
    metadata={k:v for k,v in report.items() if k not in ('records','kernel_processes','kernel_requests','reference_requests','physical_requests')}
    metadata.update(evidence_file=archive.name,evidence_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),generator='tools/rust_port/qualify_mbwr_coexistence.py')
    atomic_text(path.with_suffix('.metadata.json'),json.dumps(clean(metadata),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--compiler',required=True)
    parser.add_argument('--source-build-dir',required=True)
    parser.add_argument('--reuse-source',action='store_true')
    parser.add_argument('--reuse-kernels',action='store_true')
    parser.add_argument('--source-only',action='store_true')
    main(parser.parse_args())
