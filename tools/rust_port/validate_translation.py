"""Qualify volume-translation parameters, polynomial derivatives and full EOS."""
import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe


def parameter_cases():
    names=[json.loads(path.read_text())['ident'] for path in sorted((ROOT/'fluids').glob('*.json'))]
    for eos in ['PR','SRK']:
        for name in names:
            for ref in ['Default','tcPR']:
                for t in [20.,298.15,700.]:
                    yield {'mode':'V_CI','eos':eos,'components':name,'parameter_reference':ref,
                           'volume_shift':True,'temperature':t,'feed':[1.]}
    source=(ROOT/'src/compdatadb.f90').read_text()
    for block in re.findall(r'cidatadb\((.*?)^\s*\)',source,re.M|re.S):
        text=dict(re.findall(r'(\w+)\s*=\s*"([^"]*)"',block))
        numeric=dict(re.findall(r'\b(\w+)\s*=\s*([+-]?[\d.]+(?:[eEdD][+-]?\d+)?)',block))
        coefficients=[float(numeric.get(key,'0').replace('D','e').replace('d','e'))
                      for key in ['ciA','ciB','ciC','ciDD','ciE','ciF']]
        for t in [20.,298.15,1000.]:
            yield {'mode':'V_CI','eos':text['eosid'],'components':text['cid'],'temperature':t,'feed':[1.],
                   'parameter_reference':'tcPR' if text['cid']=='BUTANAL' else 'Default',
                   'operations':[['init','Default'],['set',0,int(numeric['c_type']),coefficients]]}
    for eos in ['PR','SRK','RK']:
        for ref in ['Default','tcPR','tcRK','missing']:
            for kind in [1,2,3,6]:
                operations=[['init',ref],['set',1,kind,[2e-6,1e-8,-1e-11,1e-14,-1e-17,1e-20]]]
                for disable in [False,True]:
                    yield {'mode':'V_CI','eos':eos,'components':'C1,CO2,H2O','temperature':300.,'feed':[.5,.3,.2],
                           'operations':operations+([['disable']] if disable else [])}


def property_cases():
    # init_cubic's len=100 reference cannot match Cp database len=40 labels.
    # Empty/slash references also exercise alpha/ci exact-token fallback.
    for eos in ['PR','SRK','RK']:
        for name in ['CO2','NC16']:
            for ref in ['', ' ', '/', 'Default/']:
                for t in [298.15,500.]:
                    for phase in ['L','V']:
                        yield dict(mode='V_PROP',eos=eos,components=name,feed=[1.],temperature=t,
                                   pressure=1e6,phase=phase,alpha='TWU',parameter_reference=ref,volume_shift=True)
    configurations=[{'parameter_reference':'Default','volume_shift':True},
                    {'parameter_reference':'tcPR'},
                    {'alpha':'TWU','parameter_reference':'tcRK','volume_shift':True}]
    for kind in [1,2,3,6]:
        configurations.append({'operations':[['init','Default'],['set',-1,kind,[2e-6,1e-8,-1e-11,1e-14,-1e-17,1e-20]]]})
    for eos in ['PR','SRK','RK']:
        for comps,z in [('C1',[1.]),('CO2',[1.]),('C3',[1.]),('C1,C2,C3',[.5,.3,.2]),
                        ('CO2,H2O',[.9,.1]),('C1,CO2,N2',[0.,.8,.2])]:
            for config in configurations:
                operations=[]
                for op in config.get('operations',[]):
                    if op[0]=='set' and op[1]==-1:
                        operations.extend([['set',i,op[2],[v*(1.+i*.2) for v in op[3]]] for i in range(len(z))])
                    else:
                        operations.append(op)
                for t in [200.,300.,500.]:
                    for p in [1e3,1e5,5e6,2e7]:
                        for phase in ['L','V']:
                            yield {'mode':'V_PROP','eos':eos,'components':comps,'feed':z,'temperature':t,
                                   'pressure':p,'phase':phase,**config,'operations':operations}


def encode(request):
    name='V_'+request['phase'] if request['mode']=='V_PROP' else request['mode']
    config='@'.join([name,request.get('alpha','RK' if request['eos']=='RK' else 'Classic'),
                     request.get('parameter_reference','Default'),'1' if request.get('volume_shift') else '0'])
    for op in request.get('operations',[]):
        if op[0]=='set':
            config+='@set:'+str(op[1])+':'+str(op[2])+':'+','.join(map(str,op[3]))
        else:
            config+='@'+':'.join(map(str,op))
    return '|'.join([request['eos'],request['components'],str(request['temperature']),str(request.get('pressure',0.)),
                     ','.join(map(str,request['feed'])),config])


def validate(compiler_directory):
    probe=verify_probe()
    inputs,expected,modes,errors=[],[],[],[]
    with ReferenceProcess(compiler_directory) as oracle:
        for request in [*parameter_cases(),*property_cases()]:
            result=oracle.call(request)
            inputs.append(encode(request));modes.append(request['mode'])
            if 'values' not in result or 'calculation_error' in result:
                errors.append({'request':request,**result});expected.append(None)
            else:
                expected.append(result['values'])
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    executable=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    result=subprocess.run([str(executable)],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True)
    actual=result.stdout.splitlines()
    assert len(actual)==len(expected)
    failures=[]
    grouped={'V_CI':[],'V_PROP':[]}
    property_absolute=[1e-14,1e-2,1e-6,1e-2,1e-7,1e-9,1e-7,1e-9,1e-7]
    for case,mode,ref,line in zip(inputs,modes,expected,actual):
        if ref is None:
            continue
        grouped[mode].append((case,ref))
        if line.startswith('ERROR'):
            failures.append({'case':case,'error':line});continue
        native=list(map(float,line.split(',')))
        assert len(native)==len(ref)
        n=len(case.split('|')[1].split(','))
        for i,(a,b) in enumerate(zip(ref,native)):
            if mode=='V_CI':
                tolerance=0. if i%12 in (0,1,11) else 1e-22+2e-12*max(abs(a),abs(b))
            else:
                absolute=property_absolute[i] if i<9 else 1e-9 if i<9+2*n else 1e-11 if i<9+3*n else 1e-8
                relative=2e-9 if i==0 else 2e-8 if i<9+n else 2e-7
                tolerance=absolute+relative*max(abs(a),abs(b))
            if not abs(a-b)<=tolerance:
                failures.append({'case':case,'quantity':i,'reference':a,'rust':b})
    for mode,name in [('V_CI','translation_parameters'),('V_PROP','translation_properties')]:
        output=ROOT/'tests/fixtures'/(name+'.csv')
        with output.open('w',newline='') as handle:
            writer=csv.writer(handle,delimiter=';');writer.writerow(['case','values...'])
            writer.writerows([[case,*values] for case,values in grouped[mode]])
        metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
        metadata.update({'cases':len(grouped[mode]),'attempted_cases':modes.count(mode),'probe':probe,
                         'reference_errors':[e for e in errors if e['request']['mode']==mode],
                         'property_absolute_tolerances':property_absolute,'parameter_relative_tolerance':2e-12,
                         'parameter_absolute_tolerance':1e-22,'exact_parameter_indices_per_component':[0,1,11],
                         'volume_relative_tolerance':2e-9,'property_relative_tolerance':2e-8,
                         'fugacity_relative_tolerance':2e-7,
                         'fugacity_derivative_absolute_tolerances':{'dt':1e-9,'dp':1e-11,'dn':1e-8},
                         'fixture_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
                         'generator':'tools/rust_port/validate_translation.py'})
        output.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/translation_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps({'parameters':len(grouped['V_CI']),'properties':len(grouped['V_PROP']),
                      'reference_errors':len(errors),'failures':len(failures)},indent=2))
    if failures:
        print(json.dumps(failures[:20],indent=2));raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
