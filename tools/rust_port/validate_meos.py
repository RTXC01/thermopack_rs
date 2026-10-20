"""Qualify pure NIST/LJ multiparameter kernels and homogeneous TV properties."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe

MODELS=[('NIST_MEOS',c) for c in ['C3','R134A','N-H2','O-H2','P-H2']]+[(e,c) for e in ['LJ_MEOS','LJTS_MEOS'] for c in ['NE','AR','KR','XE']]+[('LJTS_REDUCED','LJTS')]
TR=[.45,.65,.8,.95,.9999,1.,1.0001,1.2,1.8,3.,8.]
DELTA=[1e-8,1e-4,.01,.1,.4,.85,1.,1.15,1.8,2.6,3.2]
AMOUNTS=[.01,1.,3.7]
JETS=['alpha','tau','tautau','delta','deltadelta','deltatau']
EXT=['F','F_T','F_V','F_n','F_TT','F_TV','F_Tn','F_VV','F_Vn','F_nn']
FIELDS=[*[f'ideal_{s}' for s in JETS],*[f'residual_{s}' for s in JETS],'P','P_T','P_V','U','S','H','A','G','Cv','Cp','lnphi','source_sound',*[f'residual_{s}' for s in EXT],*[f'ideal_{s}' for s in EXT],'ancillary_liquid','ancillary_vapor']
ATOL=[1e-11]*12+[1e-5,1e-7,1e-2]+[1e-6,1e-9,1e-6,1e-6,1e-6,1e-9,1e-8,1e-10,1e-8]+[1e-10,1e-10,1e-6,1e-10,1e-10,1e-6,1e-10,1e-1,1e-6,1e-10 ]*2+[1e-10,1e-10]
RTOL=2e-10
META=['Tc','Pc','rhoc','acf','source_molar_mass','Ttriple','Ptriple','rholiq_triple','rhovap_triple','Rfit','R','maxT','maxP','sigma','eps_div_k','upPol','upExp','upExpExp']

def key(e,c,t,v,n,mode='STATE'):
    return f'MEOS|{e}:{c}|{t!r}|{v!r}|{n!r}|{mode}'

def validate(directory):
    probe=verify_probe();rows=[];metadata_rows=[];errors=[];undefined=[0]*46
    with ReferenceProcess(directory) as oracle:
        for e,c in MODELS:
            request=dict(mode='MEOS',eos=e,component=c,temperature=300.,volume=1.,amount=1.)
            initial=oracle.call(request);assert 'metadata' in initial,initial
            m=initial['metadata'];metadata_rows.append([key(e,c,300.,1.,1.,'META'),*m])
            for tr in TR:
                for delta in DELTA:
                    t=m[0]*tr;v=1./(m[2]*delta)
                    for n in AMOUNTS:
                        case=key(e,c,t,v,n)
                        result=oracle.call(dict(request,temperature=t,volume=v,amount=n))
                        if 'values' not in result:
                            errors.append(dict(case=case,result=result));continue
                        values=result['values'];assert len(values)==46
                        for j,vv in enumerate(values):
                            if vv is None:undefined[j]+=1
                        rows.append([case,*values])
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    queries=metadata_rows+rows
    output=subprocess.run([str(exe)],input='\n'.join(r[0] for r in queries)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(output)==len(queries)
    failures=[];maximum=[0.]*46;fractions=[0.]*46
    for index,(row,line) in enumerate(zip(queries,output)):
        if line.startswith('ERROR'):
            failures.append(dict(case=row[0],error=line));continue
        actual=list(map(float,line.split(',')));assert len(actual)==len(row)-1
        bad=[]
        for j,(a,b) in enumerate(zip(row[1:],actual)):
            if a is None:
                if math.isfinite(b):bad.append(j)
                continue
            atol=ATOL[j] if index>=len(metadata_rows) else 1e-10
            delta=abs(a-b);fraction=delta/(atol+RTOL*max(abs(a),abs(b)))
            if index>=len(metadata_rows):
                maximum[j]=max(maximum[j],delta);fractions[j]=max(fractions[j],fraction)
            if not math.isfinite(b) or fraction>1.:bad.append(j)
        if bad:failures.append(dict(case=row[0],quantities=bad,reference=row[1:],rust=actual))
    base=json.loads((ROOT/'.reference/metadata.json').read_text())
    for label,data,fields in [('meos',rows,FIELDS),('meos_metadata',metadata_rows,META)]:
        path=ROOT/f'tests/fixtures/{label}.csv'
        with path.open('w',newline='') as f:
            w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(['case',*fields]);w.writerows(data)
        metadata=dict(base,cases=len(data),attempted_cases=len(rows)+len(errors) if label=='meos' else len(data),reference_errors=errors if label=='meos' else [],probe=probe,
            fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),generator='tools/rust_port/validate_meos.py',quantities=fields,
            absolute_tolerances=ATOL if label=='meos' else [1e-10]*18,relative_tolerance=RTOL,
            undefined_quantities=undefined if label=='meos' else [],maximum_absolute_differences=maximum if label=='meos' else [],maximum_tolerance_fractions=fractions if label=='meos' else [],
            convention='pure default-reference NIST and LJ MEOS kernels, explicit homogeneous TV states including unstable/extrapolated states; empty fields are source non-finite properties; LJ sound uses the source meos mass field, with physical SI sound exposed separately; source metadata max/fit counts only qualified for LJ')
        path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/meos_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(cases=len(rows),metadata=len(metadata_rows),reference_errors=len(errors),undefined_quantities=dict(zip(FIELDS,undefined)),differences=len(failures),worst_tolerance_fraction=max(fractions)),indent=2))
    if failures:print(json.dumps(failures[:2],indent=2));raise SystemExit(1)
    if errors:raise SystemExit('Unexpected reference process failures')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');validate(p.parse_args().compiler_directory)
