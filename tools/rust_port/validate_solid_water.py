"""Qualify the original ice Ih Gibbs functions and owned reference shifts."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe

CONFIGS=[[],[(100.,-10000.)],[(185.,-16000.)],[(-50.,2000.)],[(100.,-10000.),(100.,-10000.)],[(100.,-10000.),(185.,-16000.)]]
TEMPERATURES=[0.,.01,1.,20.,50.,100.,150.,200.,250.,273.152519,273.16,300.,500.]
PRESSURES=[-1e6,0.,611.657,1e5,101325.,1e6,1e7,1e8,2e8,1e9]
QUANTITIES=['G','G_P','G_PP','G_T','G_TT','G_TP','V','S','H','U','A','source_Cp']
ATOL=[1e-6,1e-16,1e-25,1e-9,1e-12,1e-17,1e-16,1e-9,1e-6,1e-6,1e-6,1e-10]
RTOL=1e-11

def validate(directory):
    probe=verify_probe();rows=[];errors=[]
    for shifts in CONFIGS:
        config=':'.join(','.join(map(str,pair)) for pair in shifts) or 'raw'
        # A new process starts from the original module's default references.
        with ReferenceProcess(directory) as oracle:
            for s,g in shifts:
                result=oracle.call(dict(mode='SOLID_H2O_INIT',entropy=s,gibbs=g))
                assert result.get('initialized'),result
            for t in TEMPERATURES:
                for p in PRESSURES:
                    case='|'.join(['SOLIDH2O',config,str(t),str(p),'0','STATE'])
                    result=oracle.call(dict(mode='SOLID_H2O',temperature=t,pressure=p))
                    if 'values' not in result or 'calculation_error' in result:
                        errors.append(dict(case=case,result=result));continue
                    rows.append([case,*result['values']])
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    output=subprocess.run([str(exe)],input='\n'.join(r[0] for r in rows)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(output)==len(rows)
    failures=[];maximum=[0.]*12
    for row,line in zip(rows,output):
        if line.startswith('ERROR'):
            failures.append(dict(case=row[0],error=line,reference=row[1:]));continue
        actual=list(map(float,line.split(',')));assert len(actual)==12
        bad=[]
        for i,(a,b) in enumerate(zip(row[1:],actual)):
            maximum[i]=max(maximum[i],abs(a-b))
            if not math.isfinite(b) or abs(a-b)>ATOL[i]+RTOL*max(abs(a),abs(b)):bad.append(i)
        if bad:failures.append(dict(case=row[0],quantities=bad,reference=row[1:],rust=actual))
    path=ROOT/'tests/fixtures/solid_water.csv'
    with path.open('w',newline='') as f:
        w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(['case',*QUANTITIES]);w.writerows(rows)
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update(cases=len(rows),attempted_cases=len(rows)+len(errors),reference_errors=errors,probe=probe,
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),generator='tools/rust_port/validate_solid_water.py',
        quantities=QUANTITIES,absolute_tolerances=ATOL,relative_tolerance=RTOL,maximum_absolute_differences=maximum,
        convention='original raw Gibbs functions and sho2_init sequences, SI molar units; Cp explicitly preserves the original second MW multiplication; T=0 and extrapolated states included without asserting phase stability')
    path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/solid_water_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(cases=len(rows),reference_errors=len(errors),failures=len(failures),maximum_absolute_differences=maximum),indent=2))
    if failures:print(json.dumps(failures[:6],indent=2));raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');validate(p.parse_args().compiler_directory)
