"""Qualify full QuantumCubic TP/PH/PS and pure coexistence at cryogenic T."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
from reference import ROOT
from validate_translated_flash import validate


def recheck_stored():
    """Compare frozen source snapshots, preserving every raw source result."""
    inputs=[];metadata={}
    for mode in ['tp','ph','ps','sat']:
        path=ROOT/'tests/fixtures'/('quantum_'+mode+'.csv')
        info=json.loads(path.with_suffix('.metadata.json').read_text())
        assert hashlib.sha256(path.read_bytes()).hexdigest()==info['fixture_sha256']
        metadata[mode]=info
        with path.open(newline='') as handle:
            for row in list(csv.reader(handle,delimiter=';'))[1:]:
                inputs.append((mode,row[0],list(map(float,row[2:]))))
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    executable=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    process=subprocess.run([str(executable)],input='\n'.join(row[1] for row in inputs)+'\n',text=True,capture_output=True,check=True)
    lines=process.stdout.splitlines();assert len(lines)==len(inputs)
    failures=[];counts={mode:dict(matched_cases=0,failed_cases=0) for mode in metadata}
    for (mode,case,expected),line in zip(inputs,lines):
        failure=None
        if line.startswith('ERROR'):
            failure=dict(mode=mode,case=case,error=line,reference=expected)
        else:
            actual=list(map(float,line.split(',')));assert len(actual)==len(expected)
            bad=[];info=metadata[mode]
            for i,(a,b) in enumerate(zip(expected,actual)):
                if mode=='sat':
                    tolerance=info['saturation_absolute_tolerances'][i]+info['relative_tolerance']*max(abs(a),abs(b))
                else:
                    flag=0 if mode=='tp' else 1
                    tolerance=info['phase_flag_tolerance'] if i==flag else info['absolute_tolerance']+info['relative_tolerance']*max(abs(a),abs(b))
                if not abs(a-b)<=tolerance: bad.append(i)
            if bad: failure=dict(mode=mode,case=case,quantities=bad,reference=expected,rust=actual)
        counts[mode]['failed_cases' if failure else 'matched_cases']+=1
        if failure: failures.append(failure)
    sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
             for p in sorted((ROOT/'rust').rglob('*.rs'))}
    compiler=subprocess.check_output(['rustc','--version'],text=True).strip()
    for mode,info in metadata.items():
        info['native_comparison']={**counts[mode],'native_sources_sha256':sources,'rustc':compiler}
        (ROOT/'tests/fixtures'/('quantum_'+mode+'.metadata.json')).write_text(json.dumps(info,indent=2)+'\n')
    report=dict(native_sources_sha256=sources,rustc=compiler,summary=counts,failures=failures,
                source_fixture_sha256={mode:info['fixture_sha256'] for mode,info in metadata.items()},
                generator='tools/rust_port/validate_quantum_flash.py --fixtures-only')
    (ROOT/'tests/fixtures/quantum_native_differences.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'native_comparison':counts,'unresolved_cases':len(failures)},indent=2))
    if failures: raise SystemExit(1)


def cases():
    fluids={json.loads(path.read_text())['ident']:json.loads(path.read_text())
            for path in (ROOT/'fluids').glob('*.json')}
    for eos in ['PR','SRK','RK']:
        for comps,z in [('HE',[1.]),('NE',[1.]),('H2',[1.]),('D2',[1.]),('O-H2',[1.]),('P-H2',[1.]),
                        ('N-H2',[1.]),('E-H2',[1.]),('HE,NE',[.4,.6]),('H2,NE',[.6,.4]),
                        ('HE,C1',[.1,.9]),('HE,H2,NE',[0.,.5,.5]),('C1,CO2,N2',[.5,.3,.2])]:
            config=dict(parameter_reference='QuantumCubic',operations=[['bounds',.5,999.]])
            for t in [3.,8.,20.,40.,80.,300.]:
                for p in [1e3,1e5,1e6]:
                    for mode in ['V_TP','V_PH','V_PS']:
                        yield dict(mode=mode,eos=eos,components=comps,feed=z,temperature=t,pressure=p,**config)
            if len(z)==1 and comps not in ['N-H2','E-H2']:
                for tr in [.7,.9]:
                    t=tr*fluids[comps]['critical']['temperature']
                    yield dict(mode='V_SAT',eos=eos,components=comps,feed=z,temperature=t,pressure=0.,**config)
                    for quality in [.1,.5,.9]:
                        for mode in ['V_PH','V_PS']:
                            yield dict(mode=mode,eos=eos,components=comps,feed=z,temperature=t,quality=quality,**config)
        for t in [150.,220.,300.,450.]:
            for p in [1e5,1e6,8e6]:
                for mode in ['V_TP','V_PH','V_PS']:
                    operations=[['init','Default'],['kij',0,1,.07],['lij',1,2,.04]]
                    operations += [['beta',i,'QUANTUM',[3.+i,12.]] for i in range(3)]
                    yield dict(mode=mode,eos=eos,components='C1,C2,C3',feed=[.5,.3,.2],temperature=t,
                               pressure=p,operations=operations)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory')
    parser.add_argument('--fixtures-only',action='store_true',help='Recheck frozen source results without Fortran')
    args=parser.parse_args()
    if args.fixtures_only:
        recheck_stored()
    else:
        validate(args.compiler_directory,cases(),'quantum','tools/rust_port/validate_quantum_flash.py')
