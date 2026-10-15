"""Compare fixed-composition critical and pure saturation inversion helpers."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
from reference import ROOT,ReferenceProcess


def validate(compiler_directory):
    data={d['ident']:d for path in (ROOT/'fluids').glob('*.json')
          for d in [json.loads(path.read_text())]}
    inputs,expected=[],[]
    with ReferenceProcess(compiler_directory) as oracle:
        for eos in ['PR','SRK','RK','VDW']:
            for names,z in [('C1',[1.]),('CO2',[1.]),('C3',[1.]),
                            ('C1,C2',[0.5,0.5]),('CO2,H2O',[0.9,0.1]),
                            ('C1,C2,C3',[0.5,0.3,0.2]),('C1,CO2,N2',[0.,0.8,0.2])]:
                result=oracle.call({'mode':'PSEUDO','eos':eos,'components':names,'feed':z})
                assert 'values' in result,result
                inputs.append('|'.join([eos,names,'0','PSEUDO',','.join(map(str,z)),'V']))
                expected.append(result['values'])
            for name in ['C1','CO2','C3']:
                for tr in [0.65,0.8,0.95]:
                    result=oracle.call({'mode':'SAT_T','eos':eos,'components':name,'feed':[1.],
                                        'temperature':tr*data[name]['critical']['temperature']})
                    assert 'values' in result,result
                    inputs.append('|'.join([eos,name,str(result['pressure']),'SAT_T','1','V']))
                    expected.append(result['values'])
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    executable=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    result=subprocess.run([str(executable)],input='\n'.join(inputs)+'\n',capture_output=True,text=True,check=True)
    actual=[list(map(float,row.split(','))) for row in result.stdout.splitlines()]
    assert len(actual)==len(expected)
    for case,ref,native in zip(inputs,expected,actual):
        assert len(ref)==len(native)
        for i,(a,b) in enumerate(zip(ref,native)):
            atol=1e-12 if i in (2,3) else 1e-6
            assert abs(a-b)<=atol+2e-8*max(abs(a),abs(b)),(case,i,a,b)
    output=ROOT/'tests/fixtures/initialization.csv'
    with output.open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter=';',lineterminator='\n')
        writer.writerow(['case','values...'])
        writer.writerows([[case,*row] for case,row in zip(inputs,expected)])
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update({'cases':len(inputs),'relative_tolerance':2e-8,
                     'absolute_tolerance_by_index':[1e-6,1e-6,1e-12,1e-12,1e-6],
                     'fixture_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
                     'generator':'tools/rust_port/validate_initialization.py'})
    output.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps({'cases':len(inputs),'failures':0}))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
