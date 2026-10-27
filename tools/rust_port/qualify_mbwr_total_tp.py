"""Original classical TP dispatch and TVP Hessian qualification for total MBWR."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import subprocess

from reference import ROOT
from audit_mbwr_caloric import load_and_audit
from qualify_mbwr_caloric import R, assemble, atomic_binary, native_inputs
from validate_meos_flash import atomic_csv, atomic_text, digest
from validate_mbwr_density import cases as density_cases

RELATIVE = 2e-7


def key(q):
    return '|'.join([q['record'],repr(q['T']),repr(q['V']),repr(q['n']),q['binding'],str(q['profile']),
                     q['kind'],repr(q['P']),str(q['phase']),','.join(map(repr,q['parameters'])) if q['kind']=='SOURCE' else repr(q.get('maximum_density',60000.))])


def schedule():
    base,_=load_and_audit()
    result=[]
    for i,e in enumerate(base['records']):
        q=dict(e['request']);p=e['source_baseline'][10]
        q.update(kind='TVP' if p>0 else 'INVALID',P=p*(1. if i%2 else .8) if p>0 else 1.,phase=2)
        result.append(q)
    parameters={}
    for q in density_cases():
        parameters.setdefault(q['source_name'],[*q['critical'],q['acentric_factor']])
    for name,a in parameters.items():
        template=next(q for q in result if q['record']==name)
        points=[(tr,pr,phase,1.,0) for tr,pr,phase in itertools.product((.65,.9,.9999,1.0001,1.2,2.),(.01,.8,3.),(1,2))]
        points += [(tr,pr,phase,n,0) for tr,pr,phase,n in itertools.product((.65,1.2),(.01,3.),(1,2),(.3,3.5))]
        if name in ('mbwrcx1','mbwrcx5','mbwrcx10'):
            points += [(1.2,.8,phase,n,profile) for phase,n,profile in itertools.product((1,2),(.3,1.,3.5),(1,2,3))]
        for tr,pr,phase,n,profile in points:
            q=dict(template);q.update(T=a[0]*tr,P=a[1]*pr,V=.001,n=n,phase=phase,profile=profile,kind='SOURCE',parameters=a)
            result.append(q)
    for name in ('mbwrcx1','mbwrcx5','mbwrcx10'):
        a=parameters[name];template=next(q for q in result if q['record']==name)
        for tr,pr,phase,n in itertools.product((.65,.9,.9999,1.0001,1.2,2.),(.01,.8,3.),(1,2,3),(.3,3.5)):
            q=dict(template);q.update(T=a[0]*tr,P=a[1]*pr,V=.001,n=n,phase=phase,profile=0,kind='BOUNDED',maximum_density=60000.)
            result.append(q)
    assert len({key(q) for q in result})==len(result)
    return result


def floors(q,source=False):
    t,v,n,p=q['T'],q['V'],q['n'],q['P'];rt=R*t;e=2e-10
    volume=[v,0.,v/t,v/p,v/n];energy=[n*rt,0.,n*R,n*rt/p,rt]
    entropy=[n*R,0.,n*R/t,n*R/p,R];chemical=[rt,0.,R,rt/p,rt/n]
    intensive=[1.,0.,1./t,1./p,1./n]
    if source: return [e*n*rt/v,e*v]+[e*x for x in energy+entropy+intensive+intensive]
    return [e*n*rt/v]+[e*x for x in volume+energy+energy+entropy+energy+energy+chemical+intensive+intensive+[n*R]]+[0.,e*n*R,0.,2e-8]


def differences(q,expected,actual,source=False):
    length=22 if source else 51
    if isinstance(expected,str):return [] if isinstance(actual,str) and actual.startswith('ERROR') else ['expected error']
    if not isinstance(expected,list) or len(expected)!=length or not all(map(math.isfinite,expected)):return ['invalid expected vector']
    if not isinstance(actual,list) or len(actual)!=length or not all(map(math.isfinite,actual)):return ['invalid native vector']
    q=dict(q,V=expected[1],P=expected[0])
    flags=range(3,22,5) if source else list(range(2,46,5))+[47,49]
    if any(actual[i]!=expected[i] or actual[i] not in (0.,1.) for i in flags):return ['optional flag mismatch']
    return [(i,a,b) for i,(a,b,f) in enumerate(zip(actual,expected,floors(q,source)))
            if abs(a-b)>f+RELATIVE*max(abs(a),abs(b))]


def quantity(values):
    if not math.isfinite(values[0]):return None
    return [values[0],1.,*values[1:]] if all(map(math.isfinite,values[1:])) else [values[0],0.,0.,0.,0.]


def classical(raw):
    assert len(raw)==18
    if not all(math.isfinite(raw[i]) for i in (0,1,2,6,10,14)):return 'ERROR: non-finite original state'
    return raw[:2]+sum((quantity(raw[i:i+4]) for i in (2,6,10,14)),[])


def unresolved_source_derivative_case(q,raw):
    """Source-only conditioning classification; never selected by native error."""
    if q['kind']!='SOURCE':return False
    d=raw.get('DENSITY',[])
    return (len(d)==4 and all(map(math.isfinite,d)) and d[0]>0
            and abs(d[3])<=64*math.ulp(1.)*R*q['T']
            and abs(d[2]-q['P'])>1e-3+1e-6*max(abs(d[2]),q['P']))


def singular_source_gap(q,raw,expected,actual):
    """Validate values/S derivatives but do not qualify extremum H/Z/phi responses.

    The original unclosed zero-slope returns have ill-conditioned implicit
    responses. All raw derivative vectors remain in the evidence. Finite
    triplets or absent triplets can vary with compilation context; none of
    these declared blocks are counted as qualified derivative parity.
    """
    if not unresolved_source_derivative_case(q,raw):return False
    if not isinstance(expected,list) or not isinstance(actual,list) or len(expected)!=22 or len(actual)!=22:return False
    if not all(map(math.isfinite,expected+actual)):return False
    unqualified=set()
    for start in (3,13,18):
        if expected[start] not in (0.,1.) or actual[start] not in (0.,1.):return False
        if actual[start]==0. and actual[start:start+4]!=[0.]*4:return False
        unqualified.update(range(start,start+4))
    margin=floors(dict(q,V=expected[1],P=expected[0]),True)
    return all(abs(a-b)<=floor+RELATIVE*max(abs(a),abs(b)) for i,(a,b,floor) in enumerate(zip(actual,expected,margin)) if i not in unqualified)


def unqualified_entry(q):
    return dict(case=key(q),view='CLASSICAL',unmatched_fields=list(range(3,7))+list(range(13,17))+list(range(18,22)),reason='Unclosed original density return with roundoff-unresolved pressure slope; H/Z/lnphi derivative blocks remain unqualified regardless of coincidental numeric agreement. Original and native vectors are retained. Values and S derivatives must match.')


def map_tvp(q,tv,residual):
    t,v,n=q['T'],q['V'],q['n'];p,pt,pv,pn=tv[10:14];rt=R*t
    regular=abs(pv)>64*math.ulp(1.)*rt*n/v**2
    vt,vp,vn=(-pt/pv,1./pv,-pn/pv) if regular else (0.,0.,0.)
    def mapped(x):return quantity([x[0],x[1]+x[2]*vt,x[2]*vp,x[3]+x[2]*vn]) if regular else [x[0],0.,0.,0.,0.]
    out=[p]+(quantity([v,vt,vp,vn]) if regular else [v,0.,0.,0.,0.])
    for i in (14,18,22,26):out+=mapped(tv[i:i+4])
    out+=quantity([tv[30],-tv[22],v,tv[30]/n])
    out+=quantity([tv[34],-tv[22]/n,v/n,0.])
    z=p*v/(n*rt)
    out+=(quantity([z,z*(vt/v-1./t),z*(1./p+vp/v),0.]) if regular else [z,0.,0.,0.,0.])
    hr=-rt*t*residual[1]-rt*v*residual[2]
    out+=quantity([tv[42],-hr/(n*rt*t),-v*residual[2]/n/p,0.])
    out += [tv[38],tv[39],tv[40],tv[43],tv[44]]
    assert len(out)==51 and all(map(math.isfinite,out))
    return out


def native_run(exe,requests,env):
    lines=subprocess.check_output([str(exe)],input='\n'.join(map(key,requests))+'\n',text=True,env=env).splitlines()
    result=[];cursor=0
    for q in requests:
        kinds=['DENSITY','CLASSICAL'] if q['kind']=='SOURCE' else ['PHYSICAL'] if q['kind']=='BOUNDED' else ['BASELINE','PHYSICAL','CLASSICAL']
        entry={}
        for label in kinds:
            prefix,body=lines[cursor].split(' ',1);cursor+=1;assert prefix==label
            entry[label]=body if body.startswith('ERROR') else list(map(float,body.split(',')))
        result.append(entry)
    assert cursor==len(lines)
    return result


def main(args):
    out=ROOT/'.reference';out.mkdir(exist_ok=True);build=Path(args.source_build_dir).resolve();compiler=Path(args.compiler).resolve()
    env=dict(os.environ);env['PATH']=str(compiler.parent)+os.pathsep+env['PATH']
    identity=json.loads((out/'metadata.json').read_text(encoding='utf-8'))
    assert identity==json.loads((build/'metadata.json').read_text(encoding='utf-8'))
    assert all(hashlib.sha256((ROOT/'src'/(name+'.f90')).read_bytes()).hexdigest()==sha for name,sha in identity['sources'].items())
    dll=out/('thermopack.dll' if os.name=='nt' else 'libthermopack.so');assert hashlib.sha256(dll.read_bytes()).hexdigest()==identity['library_sha256']
    hashes=native_inputs()
    subprocess.run(['cargo','build','--release','--example','mbwr_total_tp_parity','-j','2'],cwd=ROOT,env=env,check=True)
    native=ROOT/'target/release/examples'/('mbwr_total_tp_parity.exe' if os.name=='nt' else 'mbwr_total_tp_parity')
    native_sha=hashlib.sha256(native.read_bytes()).hexdigest();requests=schedule()
    bounded=[q for q in requests if q['kind']=='BOUNDED'];proposals=native_run(native,bounded,env)
    for q,entry in zip(bounded,proposals):
        assert isinstance(entry['PHYSICAL'],list),(key(q),entry)
        q['V']=entry['PHYSICAL'][1]
    manifest=build/'objects.rsp'
    objects=[Path(s.strip().strip('"')) for s in manifest.read_text(encoding='utf-8').splitlines() if s.strip()]
    paths=[p if p.is_absolute() else build/p for p in objects]+list(build.glob('*.mod'))+[manifest]
    object_hashes={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    flags=['-O3','-cpp','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp','-fallow-argument-mismatch','-ffree-line-length-none']
    helper=ROOT/'tools/rust_port/mbwr_total_tp_reference.f90';exe=out/('mbwr_total_tp_reference.exe' if os.name=='nt' else 'mbwr_total_tp_reference')
    subprocess.run([str(compiler),*flags,'-I',str(build),str(helper),'@'+str(manifest),'-o',str(exe)],cwd=out,env=env,check=True)
    signature=dict(reference=identity,object_hashes=object_hashes,helper_sha256=digest(helper),generator_sha256=digest(Path(__file__)),
        support_sha256={name:digest(ROOT/'tools/rust_port'/name) for name in ('qualify_mbwr_caloric.py','validate_mbwr_density.py','audit_mbwr_caloric.py','generate_mbwr_data.py','generate_data.py','reference.py','validate_meos_flash.py')},
        base_evidence_sha256=hashlib.sha256((ROOT/'tests/fixtures/mbwr_caloric.audit.json.gz').read_bytes()).hexdigest(),
        compiler_sha256=hashlib.sha256(compiler.read_bytes()).hexdigest(),compiler_version=subprocess.check_output([str(compiler),'--version'],text=True,env=env).splitlines()[0],
        flags=flags,schedule=requests,executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest())
    cache=out/'mbwr_total_tp_source.json'
    if args.reuse_source:
        saved=json.loads(cache.read_text(encoding='utf-8'));without_exe=lambda x:{k:v for k,v in x.items() if k!='executable_sha256'}
        assert without_exe(saved['signature'])==without_exe(signature),'stale source collection'
        source=saved['source'];signature=saved['signature']
    else:
        text=str(len(requests))+'\n'
        for q in requests:
            mode=1 if q['kind']=='SOURCE' else 2 if q['kind']=='INVALID' else 0
            params=q.get('parameters',[300.,5e6,10000.,.1])
            text+=' '.join(map(str,[q['index'],q['T'],q['V'],q['n'],q['ideal_component'],q['profile'],mode,q['P'],q['phase'],*params]))+'\n'
        try:
            run=subprocess.run([str(exe)],input=text,text=True,capture_output=True,env=env,timeout=60.)
            source=dict(returncode=run.returncode,stdout=run.stdout,stderr=run.stderr)
        except subprocess.TimeoutExpired as e:
            decode=lambda s:s.decode('utf-8',errors='replace') if isinstance(s,bytes) else s or ''
            source=dict(censored=True,timeout=60.,returncode=None,stdout=decode(e.stdout),stderr=decode(e.stderr))
        atomic_text(cache,json.dumps(dict(signature=signature,source=source),allow_nan=False)+'\n')
    assert not source.get('censored'),'censored collection retained; never retried'
    assert source['returncode']==0,source['stderr'][-1000:]
    parsed={}
    for i,label,body in re.findall(r'^CASE=(\d+) (\w+)=([^\r\n]+)',source['stdout'],re.M):
        k=int(i),label;assert k not in parsed;parsed[k]=list(map(float,body.split()))
    grouped={}
    for (i,label),values in parsed.items():grouped.setdefault(i,{})[label]=values
    actual=native_run(native,requests,env);records=[];errors=[];known_mismatches=[];unqualified_source_derivatives=[];rows=[]
    for i,(q,native_values) in enumerate(zip(requests,actual),1):
        expected={};raw=grouped[i]
        expected['CLASSICAL']=classical(raw['CLASSICAL'])
        if unresolved_source_derivative_case(q,raw):unqualified_source_derivatives.append(unqualified_entry(q))
        if q['kind']=='SOURCE':
            d=raw['DENSITY'];expected['DENSITY']=d if all(map(math.isfinite,d)) else 'ERROR: non-finite original density'
        elif q['kind']=='INVALID':expected.update(BASELINE='ERROR: nonpositive EOS pressure',PHYSICAL='ERROR: nonpositive EOS pressure')
        else:
            tv=raw['VALUES'];ideal=raw['IDEAL'];residual=raw['RESIDUAL'];assert len(tv)==45 and len(ideal)==15 and len(residual)==10
            expected['BASELINE']=map_tvp(q,tv,residual)
            f=tv[:10].copy()
            if ideal[3]==10.:f[4]-=q['n']*999.*ideal[4]/(R*q['T']**2)
            physical=assemble(q,f,residual[3],ideal[2]*1e-3)
            expected['PHYSICAL']=map_tvp(q,physical,residual)
            # Separate unchanged TVP entry points verify H/S/lnphi against the jet-derived expectation.
            direct=raw['TVP'];probe=[expected['BASELINE'][0]]
            for start in (11,16,41):probe += [expected['BASELINE'][start]]+expected['BASELINE'][start+2:start+5]
            assert len(direct)==len(probe)==13
            direct_floors=[floors(q)[0]]
            for start in (11,16,41):direct_floors += [floors(q)[start]]+floors(q)[start+2:start+5]
            mismatch=[(j,a,b) for j,(a,b,floor) in enumerate(zip(probe,direct,direct_floors)) if not math.isfinite(b) or abs(a-b)>floor+RELATIVE*max(abs(a),abs(b))]
            if mismatch:errors.append(dict(case=key(q),view='independent original TVP',differences=mismatch))
            if q['kind']=='BOUNDED':
                assert abs(tv[10]-q['P'])<=1e-6+2e-9*q['P'],(key(q),'original EOS pressure closure')
                assert tv[12]<=0.,(key(q),'original mechanical stability')
        for label,value in native_values.items():
            if label=='DENSITY':
                gold=expected[label]
                bad=[] if isinstance(gold,str) and isinstance(value,str) and value.startswith('ERROR') else ['density shape'] if not isinstance(value,list) or not isinstance(gold,list) or len(value)!=4 or len(gold)!=4 else [j for j,(a,b,floor) in enumerate(zip(value,gold,(1e-12,0.,1e-5,1e-8))) if not math.isfinite(a) or abs(a-b)>floor+(0. if j==1 else 2e-9)*max(abs(a),abs(b))]
            else:bad=differences(q,expected[label],value,label=='CLASSICAL')
            if bad and label=='CLASSICAL' and singular_source_gap(q,raw,expected[label],value):
                known_mismatches.append(unqualified_entry(q))
            elif bad:errors.append(dict(case=key(q),view=label,differences=bad))
        clean_raw={label:[str(x) if not math.isfinite(x) else x for x in values] for label,values in raw.items()}
        records.append(dict(request=q,original=clean_raw,expected=expected,native=native_values))
        cell=lambda label:expected.get(label,'') if isinstance(expected.get(label,''),str) else ','.join(map(repr,expected[label]))
        rows.append([key(q),cell('BASELINE'),cell('PHYSICAL'),cell('CLASSICAL'),cell('DENSITY')])
    report=dict(signature=signature,cases=len(requests),errors=errors,known_mismatches=known_mismatches,unqualified_source_derivatives=unqualified_source_derivatives,records=records,source_process=source,
        native_sources_sha256=hashes,native_executable_sha256=native_sha,relative_tolerance=RELATIVE,
        counts={kind:sum(q['kind']==kind for q in requests) for kind in ('TVP','INVALID','SOURCE','BOUNDED')})
    atomic_text(out/'mbwr_total_tp_report.json',json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(cases=len(requests),counts=report['counts'],errors=len(errors),unmatched_source_states=len(known_mismatches),unqualified_source_derivative_states=len(unqualified_source_derivatives)),indent=2),flush=True)
    assert not errors,errors[:5]
    assert len(known_mismatches)<=len(unqualified_source_derivatives)
    assert hashes==native_inputs() and native_sha==hashlib.sha256(native.read_bytes()).hexdigest()
    assert object_hashes=={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    fixture=ROOT/'tests/fixtures/mbwr_total_tp.csv';atomic_csv(fixture,['case','baseline','physical','classical','density'],rows)
    report['fixture_sha256']=hashlib.sha256(fixture.read_bytes()).hexdigest()
    archive=fixture.with_suffix('.audit.json.gz');atomic_binary(archive,gzip.compress(json.dumps(report,separators=(',',':'),allow_nan=False).encode(),mtime=0))
    metadata={k:v for k,v in report.items() if k not in ('records','source_process','signature')}
    metadata['signature']={k:v for k,v in signature.items() if k!='schedule'}
    metadata.update(evidence_file=archive.name,evidence_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),generator='tools/rust_port/qualify_mbwr_total_tp.py')
    atomic_text(fixture.with_suffix('.metadata.json'),json.dumps(metadata,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler',required=True);p.add_argument('--source-build-dir',required=True)
    p.add_argument('--reuse-source',action='store_true');main(p.parse_args())
