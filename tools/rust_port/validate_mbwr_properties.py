"""Qualify original MBWR G/H/S/Z/lnphi values and TP response formulas."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe
from validate_mbwr import cases as kernel_cases
from generate_mbwr_data import records

ABSOLUTE=[.001,1e-6,1e-6,1e-5,1e-9,1e-6,1e-8,1e-7,1e-10,1e-8,
          1e-10,1e-9,1e-14,1e-9,1e-8,1e-8,1e-14,1e-8]
RELATIVE=2e-7
QUANTITIES=['P','Gres','Hres','Hres_T','Hres_P','Hres_n','Sres','Sres_T','Sres_P','Sres_n',
    'Z','Z_T','Z_P','Z_n','lnphi','lnphi_T','lnphi_P','lnphi_n']


def cases():
    for case in kernel_cases():
        if case['mode']=='MBWR_BASELINE_CORE':yield dict(case,mode='MBWR_TP_PROPERTIES')
    for record in records():
        for t in [100.,500.]:
            for rho in [100.,10000.]:
                for n in [.3,3.5]:
                    for factor in [.8,1.2]:
                        yield dict(mode='MBWR_TP_PROPERTIES',record_index=int(record['source_name'].removeprefix('mbwrcx')),
                            source_name=record['source_name'],temperature=t,volume=n/rho,moles=n,off_eos=True,pressure_factor=factor)


def encode(request,result):
    return '|'.join(['MBWR',request['source_name'],str(request['temperature']),str(request['volume']),str(request['moles']),
                     'MBWR_TP_PROPERTIES@'+str(result['values'][0])])


def validate(directory):
    probe=verify_probe();records=[];errors=[]
    with ReferenceProcess(directory) as oracle:
        for request in cases():
            result=oracle.call(request)
            if 'calculation_error' in result or 'values' not in result:
                errors.append(dict(request=request,**result));continue
            records.append((request,encode(request,result),result['values']))
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    lines=subprocess.run([str(exe)],input='\n'.join(c for _,c,_ in records)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(lines)==len(records);failures=[];maximum=[0.]*18
    for (_,case,expected),line in zip(records,lines):
        if line.startswith('ERROR'):
            failures.append(dict(case=case,error=line,reference=expected));continue
        actual=list(map(float,line.split(',')));assert len(actual)==len(expected)
        bad=[]
        for i,(a,b) in enumerate(zip(expected,actual)):
            maximum[i]=max(maximum[i],abs(a-b))
            if not math.isfinite(b) or abs(a-b)>ABSOLUTE[i]+RELATIVE*max(abs(a),abs(b)):bad.append(i)
        if bad:failures.append(dict(case=case,quantities=bad,reference=expected,rust=actual))
    path=ROOT/'tests/fixtures/mbwr_tp_properties.csv'
    with path.open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter=';',lineterminator='\n');writer.writerow(['case',*QUANTITIES])
        writer.writerows([[case,*values] for _,case,values in records])
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update(cases=len(records),attempted_cases=len(records)+len(errors),reference_errors=errors,probe=probe,
        generator='tools/rust_port/validate_mbwr_properties.py',fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,absolute_tolerances=ABSOLUTE,relative_tolerance=RELATIVE,maximum_absolute_differences=maximum,
        off_eos_cases=sum(r.get('off_eos',False) for r,_,_ in records),
        source_nonpositive_pressure_cases=sum(float(e['values'][0])<=0 for e in errors if 'values' in e),
        execution_context='oracle and native parity driver run on process main threads; standalone native source-snapshot test uses the same context; host/runtime floating environment is not modified',
        convention='original MBWR state functions at explicit P; includes source entropy pressure-response signs; nonclosing P inputs intentional; primitive units SI extensive, Z and lnphi intensive')
    path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/mbwr_properties_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(finite_cases=len(records),reference_errors=len(errors),failures=len(failures)),indent=2))
    if failures:print(json.dumps(failures[:5],indent=2));raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory');validate(parser.parse_args().compiler_directory)
