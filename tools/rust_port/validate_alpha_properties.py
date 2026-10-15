"""Qualify alpha variants through EOS properties and fugacity derivatives."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe


def cases():
    configurations=[('TWU','tcPR/tcRK/DEFAULT',None),('MC','tcPR/tcRK/DEFAULT',None),('PR78','DEFAULT',None)]
    for name in ['CLASSIC','TWU','MC','CLASSICFIT','GD','RK','SOAVE','PR','PT','SW','VDW','PR78']:
        parameters=[.2,.9,2.1] if name=='TWU' else [.8,.2,-.1] if name=='MC' else [.45] if name=='SW' else [.7]
        configurations.append(('CLASSIC','DEFAULT',(name,parameters)))
    for eos in ['PR','SRK','RK','VDW']:
        for comps,z in [('C1',[1.]),('CO2,H2O',[.9,.1]),('C1,C2,C3',[.5,.3,.2])]:
            for alpha,ref,setter in configurations:
                for t in [250.,400.,700.]:
                    for p in [1e3,1e6,2e7]:
                        for phase in ['L','V']:
                            request={'mode':'ALPHA_PROPERTIES','eos':eos,'components':comps,'feed':z,
                                     'alpha':alpha,'alpha_reference':ref,'temperature':t,'pressure':p,'phase':phase}
                            if setter:
                                request['setters']=[[i,*setter] for i in range(len(z))]
                            yield request


def encode(request):
    mode='@'.join(['A_'+request['phase'],request['alpha'],request['alpha_reference']])
    for index,name,parameters in request.get('setters',[]):
        mode+='@'+str(index)+':'+name+':'+','.join(map(str,parameters))
    return '|'.join([request['eos'],request['components'],str(request['temperature']),str(request['pressure']),
                     ','.join(map(str,request['feed'])),mode])


def validate(compiler_directory):
    probe=verify_probe()
    inputs,expected,errors=[],[],[]
    with ReferenceProcess(compiler_directory) as oracle:
        for request in cases():
            result=oracle.call(request)
            if 'values' not in result or 'calculation_error' in result:
                errors.append({'request':request,**result})
            else:
                inputs.append(encode(request))
                expected.append(result['values'])
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    executable=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    result=subprocess.run([str(executable)],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True)
    actual=result.stdout.splitlines()
    assert len(actual)==len(expected)
    failures=[]
    property_absolute=[1e-14,1e-2,1e-6,1e-2,1e-7,1e-9,1e-7,1e-9,1e-7]
    for case,ref,line in zip(inputs,expected,actual):
        native=list(map(float,line.split(',')))
        assert len(native)==len(ref)
        n=len(case.split('|')[1].split(','))
        for i,(a,b) in enumerate(zip(ref,native)):
            absolute=property_absolute[i] if i<9 else 1e-9 if i<9+2*n else 1e-11 if i<9+3*n else 1e-8
            relative=2e-9 if i==0 else 2e-8 if i<9+n else 2e-7
            if not abs(a-b)<=absolute+relative*max(abs(a),abs(b)):
                failures.append({'case':case,'quantity':i,'reference':a,'rust':b})
    output=ROOT/'tests/fixtures/alpha_properties.csv'
    with output.open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter=';',lineterminator='\n')
        writer.writerow(['case','v,p,dpdt,dpdv,hr,sr,h,s,dhdt_baseline,phi[n],phiT[n],phiP[n],phin[n*n]'])
        writer.writerows([[case,*values] for case,values in zip(inputs,expected)])
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update({'attempted_cases':len(inputs)+len(errors),'cases':len(inputs),'reference_errors':errors,
                     'probe':probe,
                     'scope':'conventional cubic EOS with alpha kernel selection; no volume translation or quantum beta',
                     'property_absolute_tolerances':property_absolute,'volume_relative_tolerance':2e-9,
                     'property_relative_tolerance':2e-8,'fugacity_relative_tolerance':2e-7,
                     'fugacity_derivative_absolute_tolerances':{'dt':1e-9,'dp':1e-11,'dn':1e-8},
                     'fixture_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
                     'generator':'tools/rust_port/validate_alpha_properties.py'})
    output.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/alpha_property_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps({'cases':len(inputs),'reference_errors':len(errors),'failures':len(failures)},indent=2))
    if failures:
        print(json.dumps(failures[:15],indent=2))
        raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
