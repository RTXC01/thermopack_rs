"""Qualify source alpha initialization, explicit setters and derivatives."""
import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
from reference import ROOT, ReferenceProcess
from generate_data import compiled_components

NAMES = ['CLASSIC','TWU','MC','GERG','CLASSICFIT','UMR','GD','RK','SOAVE','PR','PT','SW','VDW','PR78']


def records():
    source = (ROOT/'src/compdatadb.f90').read_text()
    pattern = r'type\s*\(alphadatadb\),\s*parameter\s*::\s*(twu|mc)(\d+)\s*=\s*&\s*alphadatadb\((.*?)^\s*\)'
    for family, index, block in re.findall(pattern, source, re.M|re.S|re.I):
        fields = dict(re.findall(r'(\w+)\s*=\s*"([^"]*)"',block))
        yield family.lower(),fields


def cases():
    general,_ = compiled_components(ROOT/'src/compdatadb.f90')
    for eos in ['PR','SRK','RK','VDW']:
        for name in ['C1','CO2','H2O','NC12','ETOH']:
            for alpha in NAMES:
                temperatures = [.8] if alpha=='UMR' and name!='H2O' else [.4,.8,.999999,1.,1.000001,1.5,4.]
                for tr in temperatures:
                    yield {'mode':'ALPHA','eos':eos,'components':name,'alpha':alpha,
                           'temperature':tr*general[name]['Tc']}
    for family,record in records():
        if record['eosid'] not in ['PR','SRK','VDW']:
            continue
        correlations = ['TWU','MC'] if family=='twu' else ['UMR'] if record['cid']=='H2O' else []
        for alpha in correlations:
            for tr in [.8,1.,1.5]:
                yield {'mode':'ALPHA','eos':record['eosid'],'components':record['cid'],'alpha':alpha,
                       'alpha_reference':record['ref'],'temperature':tr*general[record['cid']]['Tc']}
    for eos in ['PR','SRK','RK','VDW']:
        for name in ['C1','CO2','H2O']:
            for alpha in ['TWU','MC','UMR','GERG']:
                for ref in ['missing','tcPR/tcRK/DEFAULT','tcRK/tcPR/DEFAULT','DEFAULT/tcPR','tcPR']:
                    yield {'mode':'ALPHA','eos':eos,'components':name,'alpha':alpha,
                           'alpha_reference':ref,'temperature':.8*general[name]['Tc']}
    for eos in ['PR','SRK','RK','VDW']:
        for name in ['C1','H2O','NC12']:
            for alpha in NAMES:
                parameters = [.2,.9,2.1] if alpha=='TWU' else [.8,.2,-.1] if alpha in ['MC','GERG'] else [.45] if alpha=='SW' else [.7]
                for tr in ([.8] if alpha in ['UMR','GERG'] else [.4,.8,.999999,1.,1.000001,1.5,4.]):
                    yield {'mode':'ALPHA','eos':eos,'components':name,'alpha':'Classic',
                           'setters':[[0,alpha,parameters]],'temperature':tr*general[name]['Tc']}
    for eos in ['PR','SRK']:
        for alpha in ['TWU','MC','PR78']:
            for t in [180.,300.,500.]:
                yield {'mode':'ALPHA','eos':eos,'components':'C1,CO2,H2O','alpha':alpha,
                       'alpha_reference':'tcPR/tcRK/DEFAULT','setters':[[1,'MC',[.8,.2,-.1]]],
                       'temperature':t}


def encode(request):
    config = request.get('alpha_reference','Default')
    for index,name,params in request.get('setters',[]):
        config += '@'+str(index)+':'+name+':'+','.join(map(str,params))
    return '|'.join([request['eos'],request['components'],str(request['temperature']),request['alpha'],config,'ALPHA'])


def verify_probe():
    probe = json.loads((ROOT/'.reference/probe.metadata.json').read_text())
    library = ROOT/'.reference'/('thermopack.dll' if os.name=='nt' else 'libthermopack.so')
    helper = ROOT/'.reference'/('reference_probe.dll' if os.name=='nt' else 'libreference_probe.so')
    assert probe['library_sha256']==hashlib.sha256(library.read_bytes()).hexdigest(), 'Rebuild probe'
    assert probe['probe_sha256']==hashlib.sha256(helper.read_bytes()).hexdigest(), 'Rebuild probe'
    assert probe['source_sha256']==hashlib.sha256((ROOT/'tools/rust_port/reference_probe.f90').read_bytes().replace(b'\r\n',b'\n')).hexdigest(), 'Rebuild probe'
    return probe


def validate(compiler_directory):
    probe=verify_probe()
    inputs,expected,errors = [],[],[]
    with ReferenceProcess(compiler_directory) as oracle:
        for request in cases():
            result=oracle.call(request)
            inputs.append(encode(request))
            if 'values' not in result or 'calculation_error' in result:
                expected.append(None)
                errors.append({'request':request,**result})
            else:
                expected.append(result['values'])
    (ROOT/'.reference/alpha_oracle.json').write_text(json.dumps({'inputs':inputs,'expected':expected,'errors':errors},indent=2)+'\n')
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    executable=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    result=subprocess.run([str(executable)],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True)
    actual=result.stdout.splitlines()
    assert len(actual)==len(expected)
    failures=[]
    absolute=[0.,1e-14,1e-14,1e-14,1e-13,1e-13,1e-16]
    relative=2e-12
    numerical=[]
    for case,ref,row in zip(inputs,expected,actual):
        if ref is None:
            if not row.startswith('ERROR'):
                failures.append({'case':case,'error':'Rust returned values for a failed source call','rust':row})
            continue
        numerical.append((case,ref))
        if row.startswith('ERROR'):
            failures.append({'case':case,'error':row})
            continue
        native=list(map(float,row.split(',')))
        assert len(native)==len(ref)
        for i,(a,b) in enumerate(zip(ref,native)):
            tolerance=0. if i%7==0 else absolute[i%7]+relative*max(abs(a),abs(b))
            if not abs(a-b)<=tolerance:
                failures.append({'case':case,'quantity':i,'reference':a,'rust':b})
    output=ROOT/'tests/fixtures/alpha.csv'
    with output.open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter=';')
        writer.writerow(['case','method,params[3],alpha,dt,dtt per component...'])
        writer.writerows([[case,*values] for case,values in numerical])
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update({'attempted_cases':len(inputs),'cases':len(numerical),'reference_errors':errors,
                     'scope':'tpInitAlphaCorr, setSingleAlphaCorr and cbCalcAlphaTerm; excludes eoslibinit translated/quantum constructor overrides',
                     'absolute_tolerance_per_component':absolute,'relative_tolerance':relative,
                     'probe':probe,'fixture_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
                     'generator':'tools/rust_port/validate_alpha.py'})
    output.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/alpha_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps({'attempted_cases':len(inputs),'cases':len(numerical),'reference_errors':len(errors),'failures':len(failures)},indent=2))
    if failures:
        print(json.dumps(failures[:15],indent=2))
        raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
