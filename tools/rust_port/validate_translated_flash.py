"""Qualify translated TP/PH/PS flashes and pure saturation against this fork."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
from reference import ROOT,ReferenceProcess
from validate_translation import encode


def cases():
    fluids={json.loads(path.read_text())['ident']:json.loads(path.read_text())
            for path in (ROOT/'fluids').glob('*.json')}
    for eos in ['PR','SRK','RK']:
        for config in [dict(volume_shift=True),dict(parameter_reference='tcPR'),
                       dict(alpha='TWU',parameter_reference='tcRK',volume_shift=True),
                       dict(custom=True)]:
            for comps,z in [('C1',[1.]),('CO2',[1.]),('C3',[1.]),('C1,C2',[.5,.5]),
                            ('C1,C2,C3',[.5,.3,.2]),('CO2,N2',[.8,.2]),
                            ('CO2,H2O',[.9,.1]),('C1,CO2,N2',[0.,.8,.2])]:
                options={k:v for k,v in config.items() if k!='custom'}
                if config.get('custom'):
                    options['operations']=[['init','Default']]+[
                        ['set',i,6,[v*(1.+i*.2) for v in [2e-6,1e-8,-1e-11,1e-14,-1e-17,1e-20]]]
                        for i in range(len(z))]
                for t in [150.,220.,300.,450.]:
                    for p in [1e5,1e6,8e6]:
                        for mode in ['V_TP','V_PH','V_PS']:
                            yield dict(mode=mode,eos=eos,components=comps,feed=z,temperature=t,
                                       pressure=p,**options)
                if len(z)==1:
                    for tr in [.7,.9]:
                        t=tr*fluids[comps]['critical']['temperature']
                        yield dict(mode='V_SAT',eos=eos,components=comps,feed=z,temperature=t,
                                   pressure=0.,**options)
                        for quality in [.1,.5,.9]:
                            for mode in ['V_PH','V_PS']:
                                yield dict(mode=mode,eos=eos,components=comps,feed=z,temperature=t,
                                           quality=quality,**options)


def validate(compiler_directory):
    inputs,expected,swappable,errors=[],[],[],[]
    with ReferenceProcess(compiler_directory) as oracle:
        for request in cases():
            result=oracle.call(request)
            if 'values' not in result or 'calculation_error' in result:
                errors.append({'request':request,**result})
                continue
            data=dict(request)
            if request['mode'] in ('V_PH','V_PS'):
                data.update(temperature=result['specification'],pressure=result['pressure'])
            inputs.append(encode(data));expected.append(result['values'])
            swappable.append(result['phase_permutation_allowed'])
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    executable=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    process=subprocess.run([str(executable)],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True)
    actual=process.stdout.splitlines()
    assert len(actual)==len(expected)
    failures=[]
    grouped={'V_TP':[],'V_PH':[],'V_PS':[],'V_SAT':[]}
    for case,ref,line,can_swap in zip(inputs,expected,actual,swappable):
        mode=case.split('|')[5].split('@')[0]
        grouped[mode].append((case,can_swap,ref))
        if line.startswith('ERROR'):
            failures.append({'case':case,'error':line});continue
        native=list(map(float,line.split(',')))
        assert len(native)==len(ref)
        for i,(a,b) in enumerate(zip(ref,native)):
            if mode=='V_SAT':
                absolute=[1e-6,1e-3,1e-13,1e-11,1e-3][i]
                tolerance=absolute+2e-7*max(abs(a),abs(b))
            else:
                flag_index=0 if mode=='V_TP' else 1
                tolerance=0. if i==flag_index else 2e-6+2e-7*max(abs(a),abs(b))
            if not abs(a-b)<=tolerance:
                failures.append({'case':case,'quantity':i,'reference':a,'rust':b,
                                 'source_phase_permutation_allowed':can_swap})
    for mode,rows in grouped.items():
        output=ROOT/'tests/fixtures'/('translated_'+mode[2:].lower()+'.csv')
        with output.open('w',newline='') as handle:
            writer=csv.writer(handle,delimiter=';',lineterminator='\n')
            writer.writerow(['case','phase_permutation_allowed','values...'])
            writer.writerows([[case,int(swap),*values] for case,swap,values in rows])
        mode_errors=[e for e in errors if e['request']['mode']==mode]
        metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
        metadata.update({'cases':len(rows),'attempted_cases':len(rows)+len(mode_errors),
                         'reference_errors':mode_errors,'source_closure_scaled_tolerance':1e-8,
                         'absolute_tolerance':2e-6,'relative_tolerance':2e-7,'phase_flag_tolerance':0,
                         'saturation_absolute_tolerances':[1e-6,1e-3,1e-13,1e-11,1e-3],
                         'fixture_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
                         'generator':'tools/rust_port/validate_translated_flash.py'})
        output.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/translated_flash_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps({'cases':{mode:len(rows) for mode,rows in grouped.items()},
                      'reference_errors':len(errors),'failures':len(failures)},indent=2))
    if failures:
        print(json.dumps(failures[:20],indent=2));raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
