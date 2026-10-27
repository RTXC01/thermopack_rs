"""Exact-source total TV MBWR calorics, raw references, and full gradients."""
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
import tempfile

from reference import ROOT
from generate_data import compiled_components
from generate_mbwr_data import records
from validate_meos_flash import atomic_csv, atomic_text, digest

R=8.31446261815324
RELATIVE=2e-7
FIELDS=['F','FT','FV','Fn','FTT','FTV','FTn','FVV','FVn','Fnn',
        'P','PT','PV','Pn']
FIELDS += [q+s for q in ['U','H','S','A','G','mu'] for s in ('','T','V','n')]
FIELDS += ['Cv','Cp_present','Cp','lnphi_present','lnphi','sound_present','sound']
assert len(FIELDS)==45


def requests():
    general,cp=compiled_components(ROOT/'src/compdatadb.f90')
    available={s.upper() for s in general if s in cp}
    for record in records():
        supported=record['identifier'].upper() in available
        for t,rho,n in itertools.product((75.,150.,298.15,300.,500.,900.),(.01,100.,1000.,5000.,15000.),(.3,1.,3.5)):
            yield dict(record=record['source_name'],index=int(record['source_name'][6:]),T=t,V=n/rho,n=n,
                ideal_component=record['identifier'] if supported else 'AR',binding='DEFAULT' if supported else 'AR',profile=0)
    for record in ['mbwrcx1','mbwrcx5','mbwrcx10']:
        component=next(r['identifier'] for r in records() if r['source_name']==record)
        for t,rho,n,profile in itertools.product((298.15,450.),(10.,1000.),(.3,3.5),(1,2,3)):
            yield dict(record=record,index=int(record[6:]),T=t,V=n/rho,n=n,ideal_component=component,binding='DEFAULT',profile=profile)
    for gas in ['N-H2','O-H2','P-H2']:
        for t,rho,n in itertools.product((20.,45.,300.,1000.),(.01,1000.),(.3,1.,3.5)):
            yield dict(record='mbwrcx1',index=1,T=t,V=n/rho,n=n,ideal_component=gas,binding=gas,profile=0)


def key(q):
    return '|'.join([q['record'],repr(q['T']),repr(q['V']),repr(q['n']),q['binding'],str(q['profile'])])


def floors(q):
    t,v,n=q['T'],q['V'],q['n']; r=R; rt=r*t; e=2e-10
    f=[n,n/t,n/v,1.,n/t**2,n/t/v,1./t,n/v**2,1./v,1./n]
    p=[n*rt/v,n*r/v,n*rt/v**2,rt/v]
    energy=[n*rt,n*r,n*rt/v,rt]
    entropy=[n*r,n*r/t,n*r/v,r]
    chemical=[rt,r,rt/v,rt/n]
    return [e*s for s in f+p+energy+energy+entropy+energy+energy+chemical+[n*r]] + [0.,e*n*r,0.,e,0.,2e-8]


def differences(q, expected, actual):
    if len(expected)!=45 or not all(map(math.isfinite,expected)):
        return ['invalid expected vector']
    if len(actual)!=45 or not all(map(math.isfinite,actual)):
        return ['invalid native vector']
    return [(name,a,b) for name,floor,a,b in zip(FIELDS,floors(q),actual,expected)
            if abs(a-b)>floor+RELATIVE*max(abs(a),abs(b))]


def assemble(q,f,residual_fn,mass):
    """Thermodynamic identities applied only to original F/Hessian data."""
    t,v,n=q['T'],q['V'],q['n']; rt=R*t
    a=[rt*f[0],R*(f[0]+t*f[1]),rt*f[2],rt*f[3]]
    p=[-rt*f[2],-R*(f[2]+t*f[5]),-rt*f[7],-rt*f[8]]
    u=[-rt*t*f[1],-R*t*(2*f[1]+t*f[4]),-rt*t*f[5],-rt*t*f[6]]
    s=[-R*(f[0]+t*f[1]),-R*(2*f[1]+t*f[4]),-R*(f[2]+t*f[5]),-R*(f[3]+t*f[6])]
    h=[u[0]+p[0]*v,u[1]+v*p[1],u[2]+p[0]+v*p[2],u[3]+v*p[3]]
    g=[a[0]+p[0]*v,a[1]+v*p[1],a[2]+p[0]+v*p[2],a[3]+v*p[3]]
    mu=[rt*f[3],R*(f[3]+t*f[6]),rt*f[8],rt*f[9]]
    cv=u[1]; cp=None; phi=None; sound=None
    if abs(p[2])>64*math.ulp(1.)*n*rt/v**2:
        cp=cv-t*p[1]**2/p[2]
        if not math.isfinite(cp): cp=None
    if p[0]>0:
        phi=residual_fn-math.log(p[0]*v/(n*rt))
    if cv!=0:
        squared=-v*v/(n*mass)*(p[2]-t*p[1]**2/cv)
        if math.isfinite(squared) and squared>0: sound=math.sqrt(squared)
    result=list(f)+p+u+h+s+a+g+mu+[cv]
    for optional in (cp,phi,sound): result.extend([float(optional is not None),0. if optional is None else optional])
    return result


def native_inputs():
    return {str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/'rust').rglob('*.rs'))+[ROOT/'Cargo.toml',ROOT/'Cargo.lock']}


def public_calls(exe,env):
    with (ROOT/'tests/fixtures/mbwr_initialization_status.csv').open(encoding='utf-8',newline='') as handle:
        status={r['case'].split('|')[1]:int(r['source_status']) for r in csv.DictReader(handle,delimiter=';')}
    result=[]; seen=set()
    general,cp=compiled_components(ROOT/'src/compdatadb.f90')
    available={s.upper() for s in general if s in cp}
    for record in records():
        if record['identifier'].upper() not in available: continue
        model=record['kind'],record['identifier'].upper()
        if model in seen: continue
        seen.add(model)
        index=int(record['source_name'][6:])
        text=f"1\n{index} 298.15 .001 1 {record['identifier']} 0\n"
        try:
            run=subprocess.run([str(exe),'public'],input=text,text=True,capture_output=True,env=env,timeout=4.)
            entry=dict(record=record['source_name'],expected_initialization_status=status[record['source_name']],returncode=run.returncode,stdout=run.stdout,stderr=run.stderr)
        except subprocess.TimeoutExpired as e:
            entry=dict(record=record['source_name'],censored=True,timeout=4.,stdout=str(e.stdout),stderr=str(e.stderr))
        result.append(entry)
    return result


def atomic_binary(path,data):
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix='.'+path.name,delete=False) as handle:
        temporary=Path(handle.name); handle.write(data)
    try: temporary.replace(path)
    finally: temporary.unlink(missing_ok=True)


def main(args):
    out=ROOT/'.reference'; out.mkdir(exist_ok=True)
    build=Path(args.source_build_dir).resolve(); compiler=Path(args.compiler).resolve()
    identity=json.loads((out/'metadata.json').read_text(encoding='utf-8'))
    assert identity==json.loads((build/'metadata.json').read_text(encoding='utf-8'))
    assert all(hashlib.sha256((ROOT/'src'/(name+'.f90')).read_bytes()).hexdigest()==sha for name,sha in identity['sources'].items())
    dll=out/('thermopack.dll' if os.name=='nt' else 'libthermopack.so')
    assert hashlib.sha256(dll.read_bytes()).hexdigest()==identity['library_sha256']
    manifest=build/'objects.rsp'
    objects=[Path(line.strip().strip('"')) for line in manifest.read_text(encoding='utf-8').splitlines() if line.strip()]
    paths=[p if p.is_absolute() else build/p for p in objects]+list(build.glob('*.mod'))+[manifest]
    object_hashes={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    env=dict(os.environ); env['PATH']=str(compiler.parent)+os.pathsep+str(Path.home()/'.cargo/bin')+os.pathsep+env['PATH']
    flags=['-O3','-cpp','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp','-fallow-argument-mismatch','-ffree-line-length-none']
    source=ROOT/'tools/rust_port/mbwr_caloric_reference.f90'
    exe=out/('mbwr_caloric_reference.exe' if os.name=='nt' else 'mbwr_caloric_reference')
    subprocess.run([str(compiler),*flags,'-I',str(build),str(source),'@'+str(manifest),'-o',str(exe)],cwd=out,env=env,check=True)
    schedule=list(requests())
    signature=dict(reference=identity,object_hashes=object_hashes,helper_sha256=digest(source),generator_sha256=digest(Path(__file__)),
        compiler_sha256=hashlib.sha256(compiler.read_bytes()).hexdigest(),compiler_version=subprocess.check_output([str(compiler),'--version'],text=True,env=env).splitlines()[0],
        flags=flags,schedule=schedule,initialization_status_sha256=digest(ROOT/'tests/fixtures/mbwr_initialization_status.csv'),
        executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest())
    cache=out/'mbwr_caloric_source.json'
    if args.reuse_source:
        saved=json.loads(cache.read_text(encoding='utf-8'))
        comparable=lambda s:{k:v for k,v in s.items() if k!='executable_sha256'}
        assert comparable(saved['signature'])==comparable(signature),'stale original collection'
        original=saved['original']; public=saved['public_calls']; signature=saved['signature']
    else:
        text=str(len(schedule))+'\n'+'\n'.join(f"{q['index']} {q['T']} {q['V']} {q['n']} {q['ideal_component']} {q['profile']}" for q in schedule)+'\n'
        # communicate drains output; timeout kills/reaps the batch, never retries.
        try:
            run=subprocess.run([str(exe)],input=text,text=True,capture_output=True,env=env,timeout=45.)
            original=dict(returncode=run.returncode,stdout=run.stdout,stderr=run.stderr)
        except subprocess.TimeoutExpired as e:
            decode=lambda value: value.decode('utf-8',errors='replace') if isinstance(value,bytes) else (value or '')
            original=dict(returncode=None,censored=True,timeout=45.,stdout=decode(e.stdout),stderr=decode(e.stderr))
        public=public_calls(exe,env)
        atomic_text(cache,json.dumps(dict(signature=signature,original=original,public_calls=public),allow_nan=False)+'\n')
    assert not original.get('censored'), 'censored source batch retained in cache; no retry or qualification'
    assert original['returncode']==0,original['stderr'][-1000:]
    parsed={}
    for number,kind,body in re.findall(r'^CASE=(\d+) (VALUES|IDEAL|RAW_SOUND)=([^\r\n]+)',original['stdout'],re.M):
        k=(int(number),kind); assert k not in parsed
        parsed[k]=list(map(float,body.split()))
    assert len(parsed)==3*len(schedule)
    hashes=native_inputs()
    subprocess.run(['cargo','build','--release','--example','mbwr_caloric_parity','-j','2'],cwd=ROOT,env=env,check=True)
    native=ROOT/'target/release/examples'/('mbwr_caloric_parity.exe' if os.name=='nt' else 'mbwr_caloric_parity')
    native_sha=hashlib.sha256(native.read_bytes()).hexdigest()
    lines=subprocess.check_output([str(native)],input='\n'.join(map(key,schedule))+'\n',text=True,env=env).splitlines()
    assert len(lines)==2*len(schedule)
    errors=[]; rows=[]; records_out=[]; maxima=[0.]*45; counts=dict(default=0,explicit=0,reference_edits=0,type10=0,negative_pressure=0)
    for i,q in enumerate(schedule):
        baseline=parsed[i+1,'VALUES']; ideal=parsed[i+1,'IDEAL']; rawsound=parsed[i+1,'RAW_SOUND']
        assert len(baseline)==45 and len(ideal)==15 and len(rawsound)==1
        assert all(map(math.isfinite,baseline+ideal)),(key(q),'non-finite original TV state')
        assert all(baseline[j] in (0.,1.) for j in (39,41,43))
        physical_f=baseline[:10].copy()
        if ideal[3]==10.:
            physical_f[4]-=q['n']*999.*ideal[4]/(R*q['T']**2)
            counts['type10']+=1
        physical=assemble(q,physical_f,baseline[3]-ideal[8],ideal[2]*1e-3)
        # ideal[5:] are the ten original ideal jet fields; Fn is slot 8.
        outputs=[]
        for label,expected,line in [('baseline',baseline,lines[2*i]),('physical',physical,lines[2*i+1])]:
            if line.startswith('ERROR'):
                errors.append(dict(case=key(q),mode=label,native_error=line)); outputs.append(line); continue
            actual=list(map(float,line.split(','))); outputs.append(actual)
            bad=differences(q,expected,actual)
            if bad:errors.append(dict(case=key(q),mode=label,differences=bad))
            else:
                for j,(a,b) in enumerate(zip(actual,expected)): maxima[j]=max(maxima[j],abs(a-b))
        counts['default' if q['binding']=='DEFAULT' else 'explicit']+=1
        counts['reference_edits']+=q['profile']!=0
        counts['negative_pressure']+=baseline[10]<=0
        records_out.append(dict(request=q,source_baseline=baseline,source_ideal=ideal,source_raw_sound=[str(x) if not math.isfinite(x) else x for x in rawsound],
            physical_from_original_jet=physical,native=outputs))
        rows.append([key(q),','.join(map(repr,baseline)),','.join(map(repr,physical))])
    public_matches=0; constructor_failures=0
    for entry in public:
        assert not entry.get('censored'),entry
        assert 'PUBLIC_BEGIN' in entry['stdout'],entry
        if 'PUBLIC_RETURNED' not in entry['stdout']:
            assert entry['expected_initialization_status']==1,entry
            constructor_failures+=1; continue
        assert entry['returncode']==0,entry
        match=re.search(r'^CASE=1 VALUES=([^\r\n]+)',entry['stdout'],re.M);assert match,entry
        values=list(map(float,match[1].split()))
        q=next(q for q in schedule if q['record']==entry['record'] and q['T']==298.15 and q['n']==1. and q['V']==.001 and q['profile']==0)
        expected=parsed[schedule.index(q)+1,'VALUES']
        bad=differences(q,expected,values);assert not bad,(entry['record'],bad)
        entry['values']=values;public_matches+=1
    report=dict(signature=signature,attempts=len(schedule),cases=len(schedule),errors=errors,records=records_out,counts=counts,
        fields=FIELDS,relative_tolerance=RELATIVE,maximum_absolute_differences=maxima,native_sources_sha256=hashes,native_executable_sha256=native_sha,
        source_process=original,public_constructor_calls=public,public_tv_matches=public_matches,public_constructor_failures=constructor_failures,
        qualification='Original TV dispatch with raw ideal constants and all first TVn derivatives. Coefficient-only allocation bypasses unrelated critical initialization; absent general ideal metadata requires an explicit ideal component. Physical expectations use the original total Helmholtz Hessian with the documented type10 Cp correction. No equilibrium/fit-accuracy qualification.')
    atomic_text(out/'mbwr_caloric_report.json',json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(attempts=len(schedule),errors=len(errors),counts=counts,public_tv_matches=public_matches,public_constructor_failures=constructor_failures),indent=2),flush=True)
    assert not errors, errors[:4]
    assert hashes==native_inputs()
    assert native_sha==hashlib.sha256(native.read_bytes()).hexdigest()
    assert object_hashes=={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    fixture=ROOT/'tests/fixtures/mbwr_caloric.csv'
    atomic_csv(fixture,['case','source_baseline','physical_from_original_jet'],rows)
    report.update(generator='tools/rust_port/qualify_mbwr_caloric.py',fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest())
    archive=fixture.with_suffix('.audit.json.gz')
    atomic_binary(archive,gzip.compress(json.dumps(report,separators=(',',':'),allow_nan=False).encode('utf-8'),mtime=0))
    metadata={k:v for k,v in report.items() if k not in ('records','source_process','signature')}
    metadata['signature']={k:v for k,v in signature.items() if k!='schedule'}
    metadata.update(evidence_file=archive.name,evidence_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
        absolute_tolerance='2e-10 times ideal natural scale of each quantity/derivative; flags exact; sound floor 2e-8 m/s. See floors(request) in generator.')
    atomic_text(fixture.with_suffix('.metadata.json'),json.dumps(metadata,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler',required=True);p.add_argument('--source-build-dir',required=True)
    p.add_argument('--reuse-source',action='store_true');main(p.parse_args())
