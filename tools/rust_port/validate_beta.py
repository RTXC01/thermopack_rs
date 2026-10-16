"""Qualify beta factors and full QuantumCubic properties against this fork."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe
from validate_translation import encode


def parameters():
    for name in ['HE','H2','NE']:
        for ref in ['X'*100+'QuantumCubic','X'*89+'QuantumCubic','X'*88+'QuantumCubic','Default']:
            for t in [2.,20.,100.,300.]:
                yield dict(mode='V_BETA',eos='PR',components=name,feed=[1.],temperature=t,parameter_reference=ref)
    names=[json.loads(path.read_text())['ident'] for path in sorted((ROOT/'fluids').glob('*.json'))]
    for eos in ['PR','SRK','RK']:
        for name in names:
            for ref in ['Default','QuantumCubic','QuantumCubic/tcPR/DEFAULT']:
                for t in [2.,20.,100.,300.]:
                    yield dict(mode='V_BETA',eos=eos,components=name,feed=[1.],temperature=t,parameter_reference=ref,
                               quantum_helper=eos=='PR' and ref=='QuantumCubic/tcPR/DEFAULT')
    for eos in ['PR','SRK','RK','VDW']:
        for corr,params in [('CLASSIC',[0.,0.]),('QUANTUM',[0.,0.]),('QUANTUM',[1.25,2.]),('QUANTUM',[10.,100.])]:
            for t in [2.,20.,100.,300.]:
                yield dict(mode='V_BETA',eos=eos,components='C1,HE,H2',feed=[.5,.3,.2],temperature=t,
                           operations=[['beta',i,corr,params] for i in range(3)])


def properties():
    feeds=[('HE',[1.]),('NE',[1.]),('H2',[1.]),('D2',[1.]),('O-H2',[1.]),('P-H2',[1.]),
           ('N-H2',[1.]),('E-H2',[1.]),('HE,NE',[.4,.6]),('H2,NE',[.6,.4]),
           ('HE,C1',[.1,.9]),('HE,H2,NE',[0.,.5,.5]),('C1,CO2,N2',[.5,.3,.2])]
    for eos in ['PR','SRK','RK']:
        for comps,z in feeds:
            for t in [3.,8.,20.,40.,80.,300.]:
                for p in [1e3,1e5,1e6,1e7]:
                    for phase in ['L','V']:
                        for ref in ['QuantumCubic','QuantumCubic/tcPR/DEFAULT']:
                            yield dict(mode='V_PROP',eos=eos,components=comps,feed=z,temperature=t,pressure=p,
                                       phase=phase,parameter_reference=ref,
                                       quantum_helper=eos=='PR' and ref=='QuantumCubic/tcPR/DEFAULT')
    for eos in ['PR','SRK','RK','VDW']:
        for comps,z in [('C1,C2,C3',[.5,.3,.2]),('HE,H2,NE',[.2,.4,.4])]:
            for corr,params in [('CLASSIC',[0.,0.]),('QUANTUM',[0.,0.]),('QUANTUM',[1.25,2.]),('QUANTUM',[10.,100.])]:
                operations=[['init','Default']] if eos!='VDW' else []
                operations += [['set',i,6,[2e-6,1e-8,-1e-11,1e-14,-1e-17,1e-20]] for i in range(3)]
                operations += [['beta',i,corr,[v*(1.+.2*i) for v in params]] for i in range(3)]
                operations += [['kij',0,1,.07],['lij',1,2,.04]]
                for t in [100.,300.,700.]:
                    for p in [1e3,1e5,1e6,1e7]:
                        for phase in ['L','V']:
                            yield dict(mode='V_PROP',eos=eos,components=comps,feed=z,temperature=t,pressure=p,
                                       phase=phase,operations=operations)


def validate(compiler_directory):
    probe=verify_probe()
    inputs,expected,modes,errors=[],[],[],[]
    attempted={'V_BETA':0,'V_PROP':0}
    with ReferenceProcess(compiler_directory) as oracle:
        for request in [*parameters(),*properties()]:
            attempted[request['mode']]+=1
            result=oracle.call(request)
            if 'values' not in result or 'calculation_error' in result:
                errors.append({'request':request,**result});continue
            inputs.append(encode(request));expected.append(result['values']);modes.append(request['mode'])
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    executable=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    process=subprocess.run([str(executable)],input='\n'.join(inputs)+'\n',text=True,capture_output=True)
    if process.returncode:
        (ROOT/'.reference/beta_native_output.txt').write_text(process.stdout+'\n'+process.stderr)
        raise RuntimeError('Native parity process failed; inspect .reference/beta_native_output.txt')
    actual=process.stdout.splitlines();assert len(actual)==len(expected)
    grouped={'V_BETA':[],'V_PROP':[]};failures=[]
    absolute=[1e-14,1e-2,1e-6,1e-2,1e-7,1e-9,1e-7,1e-9,1e-7]
    for case,mode,ref,line in zip(inputs,modes,expected,actual):
        grouped[mode].append((case,ref))
        if line.startswith('ERROR'):
            failures.append({'case':case,'error':line});continue
        native=list(map(float,line.split(',')));assert len(native)==len(ref)
        n=len(case.split('|')[1].split(','))
        for i,(a,b) in enumerate(zip(ref,native)):
            if mode=='V_BETA':
                tolerance=0. if i%9 in [0,1,2,6,7,8] else 1e-14+2e-12*max(abs(a),abs(b))
            else:
                ab=absolute[i] if i<9 else 1e-9 if i<9+2*n else 1e-11 if i<9+3*n else 1e-8
                rel=2e-9 if i==0 else 2e-8 if i<9+n else 2e-7
                tolerance=ab+rel*max(abs(a),abs(b))
            if not abs(a-b)<=tolerance:
                failures.append({'case':case,'quantity':i,'reference':a,'rust':b})
    for mode,name in [('V_BETA','beta_parameters'),('V_PROP','beta_properties')]:
        output=ROOT/'tests/fixtures'/(name+'.csv')
        with output.open('w',newline='') as handle:
            writer=csv.writer(handle,delimiter=';',lineterminator='\n');writer.writerow(['case','values...'])
            writer.writerows([[case,*values] for case,values in grouped[mode]])
        metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
        metadata.update({'cases':len(grouped[mode]),'attempted_cases':attempted[mode],'probe':probe,
                         'reference_errors':[e for e in errors if e['request']['mode']==mode],
                         'property_absolute_tolerances':absolute,'parameter_absolute_tolerance':1e-14,
                         'parameter_relative_tolerance':2e-12,'exact_parameter_indices_per_component':[0,1,2,6,7,8],
                         'volume_relative_tolerance':2e-9,'property_relative_tolerance':2e-8,'fugacity_relative_tolerance':2e-7,
                         'fugacity_derivative_absolute_tolerances':{'dt':1e-9,'dp':1e-11,'dn':1e-8},
                         'fixture_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'generator':'tools/rust_port/validate_beta.py'})
        output.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/beta_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps({'cases':{mode:len(rows) for mode,rows in grouped.items()},'reference_errors':len(errors),
                      'failures':len(failures)},indent=2))
    if failures:
        print(json.dumps(failures[:20],indent=2));raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
