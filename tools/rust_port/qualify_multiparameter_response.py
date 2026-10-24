"""Independent exact-source qualification of pure MEOS response and sound.

The public TP wrapper and the direct MEOS TV sound have different mass
conventions. Physical sound is checked from original TV kernels; compatibility
sound is compared to the public TP wrapper. Every error and coordinate
difference fails qualification. Only original-source expected values enter fixtures.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
from collections import Counter
from pathlib import Path
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess
from validate_meos_flash import atomic_csv,atomic_text,digest

FIELDS=['case','has_response','cp','volume_dt','volume_dp','physical_sound','source_tp_sound','pressure']
FLOORS=[0.,1e-8,1e-12,1e-16,1e-7,1e-6,1e-6]
RTOL=2e-6
RELATIVE=[0.,2e-7,2e-7,2e-7,2e-8,RTOL,2e-8]

def finite(v):return isinstance(v,(float,int)) and math.isfinite(v)
def close(a,b,floor,relative=RTOL):return finite(a) and finite(b) and abs(a-b)<=floor+relative*max(abs(a),abs(b))

def fit_mass(q,data):
    result=data[q['model']]['source_molar_mass']
    return result*1e-3 if q['model'].startswith(('LJ_MEOS:','LJTS_MEOS:')) else result

def metadata():
    with (ROOT/'tests/fixtures/meos_metadata.csv').open() as f:
        return {r['case'].split('|')[1]:{k:float(v) for k,v in r.items() if k!='case'} for r in csv.DictReader(f,delimiter=';')}

def requests():
    report=json.loads((ROOT/'tests/fixtures/meos_flash.metadata.json').read_text())
    phases={}
    for e in report['records']:
        q=e['request'];r=e['result']
        if q['kind']=='TP' and 'quality' not in q and r.get('raw_result'):
            phases[(q['eos']+':'+q['component'],q['temperature'],q['pressure'])]=r['raw_result'][2]
    with (ROOT/'tests/fixtures/meos_flash.csv').open() as f:
        for row in csv.DictReader(f,delimiter=';'):
            case=row['case'].split('|')
            if not case[5].startswith('FLASH:TP:'):continue
            model=case[1];t=float(case[2]);p=float(case[3]);phase=phases[model,t,p]
            assert phase in (1,2,4),(row,phase)
            yield dict(model=model,t=t,p=p,quality=None,phase=phase,seed_volumes=[float(row['liquid_or_overall_volume'])])
            if t>metadata()[model]['Tc']:
                for other in (1,2):
                    if other!=phase:yield dict(model=model,t=t,p=p,quality=None,phase=other,seed_volumes=[float(row['liquid_or_overall_volume'])])
    for name in ('meos_saturation.csv','meos_saturation_recovered.csv'):
        with (ROOT/'tests/fixtures'/name).open() as f:
            for row in list(csv.reader(f,delimiter=';'))[1:]:
                case=row[0].split('|')
                if case[5]!='SAT_P':continue
                t,p,vl,vv=map(float,row[1:5])
                for quality in (0.,.1,.5,.9,1.):
                    yield dict(model=case[1],t=t,p=p,quality=quality,phase=0,seed_volumes=[vl,vv])

def key(q):
    mode='HOM:'+{1:'L',2:'V',4:'S'}[q['phase']] if q['quality'] is None else repr(q['quality'])
    return f"MPRESP|{q['model']}|{q['t']!r}|{q['p']!r}|{mode}"

def source_request(q,data):
    eos,component=q['model'].split(':');ts=ps=vs=cs=rs=1.
    if eos=='LJTS_REDUCED':
        m=data['LJTS_MEOS:NE'];ts=m['eps_div_k'];vs=6.02214076e23*m['sigma']**3
        ps=m['R']*ts/vs;cs=math.sqrt(m['R']*ts/(m['source_molar_mass']*1e-3))
        rs=m['R']
        eos,component='LJTS_MEOS','NE'
    return dict(mode='MEOS_EQ',eos=eos,component=component,temperature=q['t']*ts,
                pressure=q['p']*ps,quality=q['quality'],phase=q['phase']),dict(t=ts,p=ps,v=vs,sound=cs,r=rs)

def physical_reference(t,volumes,states,beta,mass):
    if any(not all(finite(s[i]) for i in (12,13,14,16,17,20,21)) or s[14]>=0. or s[20]<=0. for s in states):return None
    if len(states)==1:
        s=states[0];v=volumes[0];pv,pt,cv=s[14],s[13],s[20]
        c2=-v*v*(pv-t*pt*pt/cv)/mass
        return [1.,s[21],-pt/pv,1./pv,math.sqrt(c2)] if c2>0. and s[21]>0. else None
    l,v=states;vl,vv=volumes
    if (vv-vl)/(vv+vl)<1e-6:return None
    slope=(v[17]-l[17])/(t*(vv-vl))
    cv=sum(w*(s[20]-t*(s[13]-slope)**2/s[14]) for w,s in zip((1.-beta,beta),states))
    volume=(1.-beta)*vl+beta*vv
    c2=volume*volume*t*slope*slope/(mass*cv) if cv>0. else -1.
    return [0.,0.,0.,0.,math.sqrt(c2)] if c2>0. and finite(c2) else None

def validate(args):
    data=metadata();qs=list(requests());assert len({key(q) for q in qs})==len(qs)
    inputs=['meos_metadata.csv','meos_flash.csv','meos_flash.metadata.json','meos_saturation.csv','meos_saturation_recovered.csv']
    signature=dict(reference=json.loads((ROOT/'.reference/metadata.json').read_text()),
        probe=json.loads((ROOT/'.reference/probe.metadata.json').read_text()),
        worker_sha256=digest(ROOT/'tools/rust_port/oracle_worker.py'),generator_sha256=digest(Path(__file__)),
        numeric_inputs_sha256={n:digest(ROOT/'tests/fixtures'/n) for n in inputs})
    assert hashlib.sha256((ROOT/'.reference'/('thermopack.dll' if os.name=='nt' else 'libthermopack.so')).read_bytes()).hexdigest()==signature['reference']['library_sha256']
    assert hashlib.sha256((ROOT/'.reference'/('reference_probe.dll' if os.name=='nt' else 'libreference_probe.so')).read_bytes()).hexdigest()==signature['probe']['probe_sha256']
    cache=ROOT/'.reference/meos_response_source.json'
    if args.reuse_source:
        saved=json.loads(cache.read_text());assert saved['signature']==signature,'stale source cache'
        records=saved['records'];assert [r['request'] for r in records]==qs
    else:
        records=[]
        with DeadlineReferenceProcess(args.compiler_directory) as oracle:
            for i,q in enumerate(qs):
                request,scales=source_request(q,data)
                result=oracle.call(request,timeout=4.)
                records.append(dict(request=q,source_request=request,scales=scales,result=result))
                if i%100==0:print(i,'source response/sound calls',flush=True)
        atomic_text(cache,json.dumps(dict(signature=signature,records=records),allow_nan=False)+'\n')
    subprocess.run(['cargo','build','--release','--example','parity','-j','2'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    output=subprocess.run([str(exe)],input='\n'.join(key(q) for q in qs)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(output)==len(qs)
    matches=[];differences=[];source_errors=[];source_unqualified=[];independent=[];native_errors=[]
    kernel_cache={}
    with DeadlineReferenceProcess(args.compiler_directory) as oracle:
        def kernels(q,volumes,t):
            result=[];eos,c=q['model'].split(':')
            for volume in volumes:
                k=(eos,c,t,volume)
                if k not in kernel_cache:
                    checked=oracle.call(dict(mode='MEOS',eos=eos,component=c,temperature=t,volume=volume),timeout=4.)
                    assert 'values' in checked,(q,checked)
                    kernel_cache[k]=checked['values']
                result.append(kernel_cache[k])
            return result
        for e,line in zip(records,output):
            q=e['request'];r=e['result'];s=e['scales'];case=key(q);t=q['t']
            if line.startswith('ERROR'):
                seed=kernels(q,q['seed_volumes'],t)
                mass=fit_mass(q,data)
                expected=physical_reference(t,q['seed_volumes'],seed,q['quality'],mass)
                entry=dict(case=case,native_error=line,source=r,physical_input=expected)
                if expected:native_errors.append(entry)
                else:source_unqualified.append(entry)
                continue
            actual=list(map(float,line.split(',')));assert len(actual)==12 and all(map(finite,actual)),(case,line)
            has,cp,vt,vp,w,baseline,p,mass,vl,vv,beta,phase=actual
            assert has in (0.,1.) and min(w,p,mass,vl,vv)>0. and baseline>=0.,(case,actual)
            assert close(mass,fit_mass(q,data),0.,1e-14),(case,'fit molar mass',mass)
            assert phase in (0.,1.,2.,4.) and (phase==0.)==(q['quality'] is not None),(case,actual)
            if phase==0.:assert 0.<=beta<=1. and beta==q['quality'] and vv>vl,(case,actual)
            else:assert beta in (-1.,0.,1.),(case,actual)
            assert p==q['p'],(case,'pressure input was changed',p,q['p'])
            volumes=[vl,vv] if phase==0. else [vl]
            states=kernels(q,volumes,t)
            expected=physical_reference(t,volumes,states,beta,fit_mass(q,data))
            assert expected is not None,(case,'native result is thermally/mechanically invalid')
            for i,(a,b) in enumerate(zip((has,cp,vt,vp,w),expected)):
                assert close(a,b,[0.,1e-8,1e-12,1e-16,1e-7][i],2e-8),(case,'kernel',i,a,b)
            gas=data[q['model']]['R']
            chemical_gap=None
            if phase==0.:
                chemical_gap=math.log(vv/vl)+states[0][6]-states[1][6]+p*(vl-vv)/(gas*t)
                assert abs(chemical_gap)<2e-8,(case,'native chemical equilibrium',chemical_gap)
            pressure_errors=[st[12]-p for st in states]
            for vol,error in zip(volumes,pressure_errors):
                assert abs(error)<2e-8*p+512.*math.ulp(1.)*gas*t/vol,(case,'native pressure',error)
            independent.append(dict(case=case,kernel_expected=expected,native=actual,pressure_errors=pressure_errors,chemical_potential_gap=chemical_gap))
            if 'reference_timeout' in r or 'reference_error' in r or 'calculation_error' in r:
                source_errors.append(e);continue
            valid=(finite(r.get('sound')) and r['sound']>=0. and all(finite(v) and v<0. for v in r['pressure_dv']))
            if q['quality'] is None:valid=valid and all(map(finite,r['response'])) and r['response'][0]>0.
            else:valid=valid and len(r['volumes'])==2 and (r['volumes'][1]-r['volumes'][0])/(r['volumes'][1]+r['volumes'][0])>=1e-6
            for vol,error in zip(r['volumes'],r['pressure_errors']):
                valid=valid and abs(error)<2e-8*r['pressure']+512.*math.ulp(1.)*data['LJTS_MEOS:NE']['R']*e['source_request']['temperature']/vol
            if not valid:source_unqualified.append(e);continue
            response=[r['response'][0]/s['r'],r['response'][1]*s['t']/s['v'],r['response'][2]*s['p']/s['v']] if q['model'].startswith('LJTS_REDUCED') and has else r['response']
            wanted=[float(has),*(response if has else [0.,0.,0.]),expected[4],r['sound']/s['sound'],r['pressure']/s['p']]
            values=[has,cp,vt,vp,w,baseline,p];bad=[]
            for i,(a,b) in enumerate(zip(values,wanted)):
                floor=FLOORS[i] if i!=6 else 1e-6/s['p']
                if not close(a,b,floor,RELATIVE[i]):bad.append([i,a,b])
            if bad:differences.append(dict(case=case,request=q,native=values,source_expected=wanted,differences=bad,source=r))
            else:matches.append([case,*wanted])
    fixture=ROOT/'tests/fixtures/multiparameter_response_source.csv'
    atomic_csv(fixture,FIELDS,matches)
    meta=dict(signature,cases=len(matches),fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),
        generator='tools/rust_port/qualify_multiparameter_response.py',attempts=len(records),records=records,
        differences=differences,source_errors=source_errors,source_unqualified=source_unqualified,
        independent_native_kernel_checks=independent,native_errors=native_errors,
        relative_tolerances=RELATIVE,absolute_tolerances=FLOORS,
        native_sources_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in sorted((ROOT/'rust').rglob('*.rs'))+[ROOT/'Cargo.toml',ROOT/'Cargo.lock']},
        qualification='Matched original TP response/sound expectations, independent original TV-kernel physical sound. All 854 source calls close; errors, unqualified states and coordinate differences fail the audit.')
    atomic_text(fixture.with_suffix('.metadata.json'),json.dumps(meta,indent=2,allow_nan=False)+'\n')
    atomic_text(ROOT/'.reference/meos_response_report.json',json.dumps(meta,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(attempts=len(records),matched=len(matches),differences=len(differences),source_errors=len(source_errors),source_unqualified=len(source_unqualified),independent=len(independent),native_errors=len(native_errors)),indent=2),flush=True)
    assert not native_errors,native_errors[:3]
    assert not differences, differences[:3]
    assert not source_errors, source_errors[:3]
    assert not source_unqualified, source_unqualified[:3]

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');p.add_argument('--reuse-source',action='store_true');validate(p.parse_args())
