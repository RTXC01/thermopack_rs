"""Qualify internal residual TP mixture kernels; public construction is broken."""
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import subprocess

from reference import ROOT
from deadline_reference import DeadlineReferenceProcess
from validate_meos_flash import atomic_csv, atomic_text, digest

NAMES = ['C3', 'R134A', 'N-H2', 'O-H2', 'P-H2']
R = 8.31446261815324
RELATIVE = 2e-8


def fields(n):
    result = []
    def block(names, floors):
        result.extend(zip(names, floors))
    block(['Z','ZT','ZP'], [1e-10,1e-12,1e-15])
    block([f'Zx{i}' for i in range(n)], [1e-10]*n)
    block(['H','HT','HP'], [2e-7,1e-8,1e-10])
    block([f'Hx{i}' for i in range(n)], [2e-7]*n)
    block(['S','ST','SP'], [2e-9,1e-10,1e-12])
    block([f'Sx{i}' for i in range(n)], [2e-9]*n)
    block(['G','GT','baseline_GP','physical_GP'], [2e-7,2e-9,2e-7,1e-10])
    block([f'Gx{i}' for i in range(n)], [2e-7]*n)
    block([f'phi{i}' for i in range(n)], [1e-10]*n)
    block([f'lnphiT{i}' for i in range(n)], [1e-12]*n)
    block([f'lnphiP{i}' for i in range(n)], [1e-15]*n)
    block([f'lnphiN{i}_{j}' for j in range(n) for i in range(n)], [0.]*(n*n))
    block([f'V{i}' for i in range(n)], [1e-12]*n)
    assert len(result) == 13+8*n+n*n
    return result


def requests():
    configurations = [([i],[1.]) for i in range(5)]
    configurations += [(list(pair),[x,1.-x]) for pair in itertools.combinations(range(5),2) for x in (0.,.5,1.)]
    configurations += [([2,3,4],x) for x in ([.2,.3,.5],[1.,0.,0.],[0.,0.,1.])]
    configurations += [(list(range(5)),[float(i==j) for i in range(5)]) for j in range(5)]
    configurations += [(list(range(5)),x) for x in ([.2]*5,[.05,.1,.2,.25,.4])]
    for ids,x in configurations:
        for t,p,phase in itertools.product((450.,600.),(1e4,1e5,1e6,1e7),(1,2)):
            yield dict(ids=ids,x=x,T=t,P=p,phase=phase)
    for t,p,phase in itertools.product((45.,60.),(1e5,1e6,1e7),(1,2)):
        yield dict(ids=[2,3,4],x=[.2,.3,.5],T=t,P=p,phase=phase)
    for phase in (1,2):
        yield dict(ids=[0,1],x=[.3,.7],T=180.,P=5e6,phase=phase)


def key(q):
    return '|'.join([repr(q['T']),repr(q['P']),str(q['phase']),
        ','.join(NAMES[i] for i in q['ids']),','.join(map(repr,q['x']))])


def collect(exe, q, env):
    text = f"{len(q['ids'])} {q['T']} {q['P']} {q['phase']}\n"
    text += ' '.join(str(i+1) for i in q['ids'])+'\n'+' '.join(map(repr,q['x']))+'\n'
    try:
        run = subprocess.run([str(exe)], input=text, text=True, capture_output=True, env=env, timeout=4.)
    except subprocess.TimeoutExpired as e:
        return dict(censored=True,timeout=4.,stdout=str(e.stdout),stderr=str(e.stderr))
    output = dict(returncode=run.returncode,stdout=run.stdout,stderr=run.stderr)
    for prefix,name in [('VALUES=','values'),('CHECKS=','checks')]:
        selected = [s for s in run.stdout.splitlines() if s.startswith(prefix)]
        if len(selected)==1:
            output[name] = [float(s) for s in selected[0][len(prefix):].split()]
    return output


def source_closed(q, result):
    n = len(q['ids'])
    if result.get('returncode') != 0 or len(result.get('values',[])) != len(fields(n)) or len(result.get('checks',[])) != 4*n:
        return False
    if not all(map(math.isfinite,result['values']+result['checks'])):
        return False
    for i,v in enumerate(result['values'][-n:]):
        p,pv,cv,phase = result['checks'][4*i:4*i+4]
        if v<=0. or abs(p-q['P']) > 2e-8*q['P']+512.*math.ulp(1.)*R*q['T']/v or pv>=0. or cv<=0. or phase not in (1.,2.,4.):
            return False
    return True


def differences(expected, actual, n):
    if len(actual)!=len(fields(n)) or not all(map(math.isfinite,actual)):
        return ['invalid native vector']
    return [(name,a,b) for (name,floor),a,b in zip(fields(n),actual,expected)
            if abs(a-b)>floor+RELATIVE*max(abs(a),abs(b))]


def native_inputs():
    paths = sorted((ROOT/'rust').rglob('*.rs'))+[ROOT/'Cargo.toml',ROOT/'Cargo.lock']
    return {str(p.relative_to(ROOT)):digest(p) for p in paths}


def main(args):
    out = ROOT/'.reference'; out.mkdir(exist_ok=True)
    compiler = Path(args.compiler).resolve()
    source_build = Path(args.source_build_dir).resolve()
    source = ROOT/'tools/rust_port/meos_ideal_mixture_reference.f90'
    manifest = source_build/'objects.rsp'
    objects = [Path(s.strip().strip('"')) for s in manifest.read_text(encoding='utf-8').splitlines() if s.strip()]
    objects = [p if p.is_absolute() else source_build/p for p in objects]
    modules = list(source_build.glob('*.mod'))
    provenance = {str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in objects+modules+[manifest]}
    identity = json.loads((out/'metadata.json').read_text(encoding='utf-8'))
    assert identity==json.loads((source_build/'metadata.json').read_text(encoding='utf-8'))
    assert identity['library_sha256']==hashlib.sha256((out/('thermopack.dll' if os.name=='nt' else 'libthermopack.so')).read_bytes()).hexdigest()
    probe_identity = json.loads((out/'probe.metadata.json').read_text(encoding='utf-8'))
    probe = out/('reference_probe.dll' if os.name=='nt' else 'libreference_probe.so')
    assert probe_identity['library_sha256']==identity['library_sha256']
    assert probe_identity['probe_sha256']==hashlib.sha256(probe.read_bytes()).hexdigest()
    assert probe_identity['source_sha256']==digest(ROOT/'tools/rust_port/reference_probe.f90')
    assert all(hashlib.sha256((ROOT/'src'/(name+'.f90')).read_bytes()).hexdigest()==sha
               for name,sha in identity['sources'].items()), 'compiled original source differs from checkout'
    env = dict(os.environ)
    env['PATH'] = str(compiler.parent)+os.pathsep+str(Path.home()/'.cargo/bin')+os.pathsep+env['PATH']
    flags = ['-O3','-cpp','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp',
             '-fallow-argument-mismatch','-ffree-line-length-none']
    exe = out/('meos_ideal_mixture_reference.exe' if os.name=='nt' else 'meos_ideal_mixture_reference')
    subprocess.run([str(compiler),*flags,'-I',str(source_build),str(source),'@'+str(manifest),'-o',str(exe)],env=env,cwd=out,check=True)
    schedule = list(requests())
    signature = dict(reference=identity,original_objects_sha256=provenance,source_sha256=digest(source),
        generator_sha256=digest(Path(__file__)),requests=schedule,flags=flags,
        compiler_sha256=hashlib.sha256(compiler.read_bytes()).hexdigest(),
        compiler_version=subprocess.check_output([str(compiler),'--version'],env=env,text=True).splitlines()[0],
        executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest())
    # Relinking can encode a timestamp; source reuse binds all actual source,
    # object and input identities while retaining the historical executable hash.
    cache = out/'meos_ideal_mixture_source.json'
    if args.reuse_source:
        saved = json.loads(cache.read_text(encoding='utf-8'))
        comparable = lambda x: {k:v for k,v in x.items() if k!='executable_sha256'}
        assert comparable(saved['signature'])==comparable(signature), 'stale original-source collection'
        records,constructor = saved['records'],saved['public_constructor']
        signature = saved['signature']
    else:
        run = subprocess.run([str(exe),'public'],env=env,capture_output=True,text=True,timeout=4.)
        constructor = dict(returncode=run.returncode,stdout=run.stdout,stderr=run.stderr)
        records = [dict(request=q,source=collect(exe,q,env)) for q in schedule]
        atomic_text(cache,json.dumps(dict(signature=signature,public_constructor=constructor,records=records),allow_nan=False)+'\n')
    assert 'CONSTRUCTOR_BEGIN' in constructor['stdout'] and 'CONSTRUCTOR_RETURNED' not in constructor['stdout']
    assert 'Wrong input for multiparameter EoS.' in constructor['stdout']+constructor['stderr']
    native_sources = native_inputs()
    subprocess.run(['cargo','build','--release','--example','meos_ideal_mixture_parity','-j','2'],cwd=ROOT,env=env,check=True)
    native = ROOT/'target/release/examples'/('meos_ideal_mixture_parity.exe' if os.name=='nt' else 'meos_ideal_mixture_parity')
    native_executable = hashlib.sha256(native.read_bytes()).hexdigest()
    kernel_paths = [out/'reference_probe.dll' if os.name=='nt' else out/'libreference_probe.so',out/'probe.metadata.json',
                   ROOT/'tools/rust_port/oracle_worker.py',ROOT/'tools/rust_port/reference.py',ROOT/'tools/rust_port/deadline_reference.py']
    kernel_inputs = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in kernel_paths}
    lines = subprocess.check_output([str(native)],input='\n'.join(map(key,schedule))+'\n',text=True,env=env).splitlines()
    assert len(lines)==len(records)
    errors = []; matches = []; independent = []; original_states = {}
    with DeadlineReferenceProcess(compiler.parent) as oracle:
        for entry,line in zip(records,lines):
            q = entry['request']; n = len(q['ids']); case = key(q)
            if not source_closed(q,entry['source']):
                errors.append(dict(case=case,kind='unclosed_source')); continue
            if line.startswith('ERROR'):
                errors.append(dict(case=case,kind='native_error',output=line)); continue
            values = list(map(float,line.split(','))); entry['native']=values
            bad = differences(entry['source']['values'],values,n)
            if bad:
                errors.append(dict(case=case,kind='coordinates',differences=bad)); continue
            lookup = dict(zip((s[0] for s in fields(n)),values))
            for i,(component,v) in enumerate(zip(q['ids'],values[-n:])):
                kernel_key = (NAMES[component],q['T'],v)
                if kernel_key not in original_states:
                    kernel = oracle.call(dict(mode='MEOS',eos='NIST_MEOS',component=NAMES[component],temperature=q['T'],volume=v))
                    assert 'values' in kernel, kernel
                    original_states[kernel_key]=kernel['values']
                st = original_states[kernel_key]; t,p=q['T'],q['P']; rt=R*t
                assert v>0. and st[14]<0. and st[20]>0.
                assert abs(st[12]-p)<=2e-8*p+512.*math.ulp(1.)*rt/v, (case,'original pressure closure',st[12],p)
                expected = dict(Hx=rt*st[7]+p*v-rt,Sx=R*(st[7]-st[6]+math.log(p*v/rt)),
                                Gx=rt*(st[6]+p*v/rt-1.-math.log(p*v/rt)))
                for label,target in expected.items():
                    floor = 2e-9 if label=='Sx' else 2e-7
                    actual = lookup[f'{label}{i}']
                    assert abs(actual-target)<=floor+RELATIVE*max(abs(actual),abs(target)), (case,label,actual,target)
                independent.append(dict(case=case,component=i,original_state=kernel_key))
            matches.append([case,','.join(map(repr,entry['source']['values'])),
                            ','.join(repr(floor) for _,floor in fields(n))])
    report = dict(signature=signature,public_constructor=constructor,attempts=len(records),matched=len(matches),
        errors=errors,records=records,independent_native_kernel_checks=independent,
        original_kernel_states=[dict(component=k[0],T=k[1],V=k[2],values=v) for k,v in original_states.items()],
        fields_by_component_count={str(n):fields(n) for n in sorted({len(q['ids']) for q in schedule})},
        relative_tolerance=RELATIVE,qualification='Internal residual TP kernels only. Public NIST_MEOS_MIX constructor fails in unchanged original; helper constructs the correct dynamic type with supported pure labels. No total caloric or phase-solver parity claim.',
        native_sources_sha256=native_sources,native_executable_sha256=native_executable,
        kernel_oracle_inputs_sha256=kernel_inputs)
    atomic_text(out/'meos_ideal_mixture_report.json',json.dumps(report,indent=2,allow_nan=False)+'\n')
    assert not errors, errors[:5]
    assert provenance=={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in objects+modules+[manifest]}
    assert kernel_inputs=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in kernel_paths}
    assert native_executable==hashlib.sha256(native.read_bytes()).hexdigest()
    fixture = ROOT/'tests/fixtures/meos_ideal_mixture.csv'
    assert native_sources==native_inputs(), 'native source changed during qualification'
    atomic_csv(fixture,['case','values','absolute_tolerances'],matches)
    report.update(cases=len(matches),fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),generator='tools/rust_port/qualify_meos_ideal_mixture.py')
    atomic_text(fixture.with_suffix('.metadata.json'),json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(attempts=len(records),matched=len(matches),independent=len(independent),original_states=len(original_states),errors=len(errors)),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--compiler',required=True); p.add_argument('--source-build-dir',required=True)
    p.add_argument('--reuse-source',action='store_true'); main(p.parse_args())
