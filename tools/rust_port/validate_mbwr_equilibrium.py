"""Qualify MBWR metadata, local critical equations and bounded stable roots.

The independent critical solve uses only the original pressure/derivative
routine, differencing its density column. Failed source iterates are retained.
"""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe
from generate_mbwr_data import records,stored_properties
from generate_data import compiled_components

def key(record,t,x,y,mode):
    return '|'.join(['MBWR',record['source_name'],str(t),str(x),str(y),mode])

def close(a,b,atol,rtol):
    return math.isfinite(b) and abs(a-b)<=atol+rtol*max(abs(a),abs(b))

def original_critical(oracle,record,t,rho):
    """Damped log-space Newton; no calls to the native critical solver."""
    index=int(record['source_name'][6:])
    def evaluate(t,r):
        result=oracle.call(dict(mode='MBWR_CRITICAL_JET',record_index=index,temperature=t,density=r))
        if 'values' not in result or 'calculation_error' in result or 'reference_error' in result:
            raise ArithmeticError(str(result))
        return result['values']
    for iteration in range(80):
        q=evaluate(t,rho)
        norm=math.hypot(q[1]/(8.31446261815324*t),rho*q[2]/(8.31446261815324*t))
        if norm<2e-10 and q[0]>0:
            return [t,q[0],rho,q[0]/(rho*8.31446261815324*t)]
        h=rho*1e-4
        hi=evaluate(t,rho+h);lo=evaluate(t,rho-h)
        a,b,c,d=t*q[3],rho*q[2],t*q[4],rho*(hi[2]-lo[2])/(2*h)
        det=a*d-b*c
        if not math.isfinite(det) or det==0:break
        st=(-q[1]*d+b*q[2])/det;sr=(-a*q[2]+q[1]*c)/det
        length=min(1.,.4/max(abs(st),abs(sr)))
        accepted=False
        for _ in range(24):
            nt=t*math.exp(length*st);nr=rho*math.exp(length*sr)
            nq=evaluate(nt,nr)
            nn=math.hypot(nq[1]/(8.31446261815324*nt),nr*nq[2]/(8.31446261815324*nt))
            if nn<norm:
                t,rho=nt,nr;accepted=True;break
            length*=.5
        if not accepted:break
    raise ArithmeticError(dict(reason='independent original-pressure critical solve did not converge',temperature=t,density=rho,norm=norm))

def run_native(cases):
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    output=subprocess.run([str(exe)],input='\n'.join(cases)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(output)==len(cases),(len(output),len(cases))
    return [line if line.startswith('ERROR') else list(map(float,line.split(','))) if line else [] for line in output]

def fixture(name,rows,quantities,probe,**extra):
    path=ROOT/'tests/fixtures'/name
    with path.open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter=';',lineterminator='\n')
        writer.writerow(['case',*quantities]);writer.writerows(rows)
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update(cases=len(rows),probe=probe,generator='tools/rust_port/validate_mbwr_equilibrium.py',
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),quantities=quantities,**extra)
    path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')

def validate(directory):
    probe=verify_probe();metadata_rows=[];critical_rows=[];density_rows=[]
    critical_errors=[];density_errors=[];failed_iterates=[];selection_differences=[]
    critical_unavailable=[];density_unavailable=[]
    general,_=compiled_components(ROOT/'src/compdatadb.f90')
    tables={(int(kind),identifier):v for kind,identifier,v in stored_properties() if kind!='triple'}
    inputs=[]
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    with ReferenceProcess(directory) as oracle:
        for record in records():
            index=int(record['source_name'][6:]);id=record['identifier']
            result=oracle.call(dict(mode='MBWR_METADATA',record_index=index))
            assert 'calculation_error' not in result and 'values' in result,result
            metadata_rows.append([key(record,0,0,0,'MBWR_METADATA'),*result['values']])
            stored=tables.get((record['kind'],id));component=general.get(id)
            if not component:
                density_unavailable.append(dict(record=record['source_name'],reason='no compiled source acentric factor'))
            if stored:
                initial=[stored[0],stored[1],stored[2]*1000.]
            elif component:
                initial=[component['Tc'],component['Pc'],component['Pc']/(component['Zc']*8.31446261815324*component['Tc'])]
            else:
                critical_unavailable.append(dict(record=record['source_name'],reason='no source table or compiled component initial values'))
                continue
            # Physical critical points are solved from identical explicit seeds.
            try:
                critical=original_critical(oracle,record,initial[0],initial[2])
                critical_rows.append([key(record,initial[0],initial[2],0,'MBWR_CRITICAL'),*critical])
            except ArithmeticError as error:
                critical_errors.append(dict(record=record['source_name'],initial=initial,error=str(error)))
            # Exercise original density iterations with explicit source metadata.
            # Unknown global fluids have no source acentric factor; record them.
            if not component:
                continue
            maximum=max(4.*initial[2],2000./math.sqrt(record['parameters'][-1]))
            for tr in [.65,.9,.9999,1.0001,1.2,2.]:
                for pr in [.01,.1,.8,1.,1.2,3.]:
                    for phase in [1,2]:
                        request=dict(mode='MBWR_DENSITY',record_index=index,temperature=initial[0]*tr,
                            pressure=initial[1]*pr,phase=phase,critical=initial,acentric_factor=component['acf'])
                        result=oracle.call(request)
                        if 'values' not in result or 'calculation_error' in result:
                            density_errors.append(dict(request=request,result=result));continue
                        rho,found,p,dpdrho=result['values']
                        native_case=key(record,request['temperature'],request['pressure'],maximum,'MBWR_ROOTS')
                        if not (0.<rho<=maximum and dpdrho>=0. and close(p,request['pressure'],.01,2e-6)):
                            failed_iterates.append(dict(request=request,values=result['values'],maximum_density=maximum))
                            continue
                        density_rows.append([native_case,rho,phase,found,p,dpdrho])
                        inputs.append(dict(case=native_case,phase=phase,found=found,request=request,reference=result['values']))
    failures=[]
    for rows,atol,rtol in [(metadata_rows,[1e-12]*8,2e-14),(critical_rows,[2e-5,.1,.01,1e-7],2e-7)]:
        for row,actual in zip(rows,run_native([r[0] for r in rows])):
            if isinstance(actual,str) or len(actual)!=len(row)-1 or any(not close(a,b,atol[i],rtol) for i,(a,b) in enumerate(zip(row[1:],actual))):
                failures.append(dict(case=row[0],reference=row[1:],rust=actual))
    for row,info,roots in zip(density_rows,inputs,run_native([r[0] for r in density_rows])):
        if isinstance(roots,str) or not any(close(row[1],r,1e-3,2e-6) for r in roots):
            failures.append(dict(**info,rust=roots));continue
        selected=roots[-1] if info['phase']==1 else roots[0]
        if not close(row[1],selected,1e-3,2e-6):
            selection_differences.append(dict(**info,rust_roots=roots,selected_density=selected))
    fixture('mbwr_metadata.csv',metadata_rows,['critical_present','Tc','Pc','rhoc','Zc','triple_present','Ttriple','Ptriple'],probe,
        absolute_tolerance=1e-12,relative_tolerance=2e-14,convention='original explicit critical/triple tables and lookup presence flags')
    fixture('mbwr_critical.csv',critical_rows,['Tc','Pc','rhoc','Zc'],probe,reference_errors=critical_errors,unavailable_records=critical_unavailable,
        absolute_tolerances=[2e-5,.1,.01,1e-7],relative_tolerance=2e-7,
        convention='independent log-space Newton of original MBWR_pressure first/second density derivatives; no native solver calls; not the coarse source critical grid routine')
    fixture('mbwr_roots.csv',density_rows,['density','phase_in','phase_found','pressure','pressure_rho'],probe,
        reference_errors=density_errors,unavailable_records=density_unavailable,unqualified_returns=failed_iterates,
        phase_selection_differences=selection_differences,absolute_tolerance=1e-3,relative_tolerance=2e-6,
        convention='membership in native bounded stable root set, not source iteration trajectory or requested phase identity; explicit source metadata or compiled component initial guesses')
    (ROOT/'.reference/mbwr_equilibrium_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(metadata=len(metadata_rows),critical=len(critical_rows),critical_errors=len(critical_errors),critical_unavailable=len(critical_unavailable),
        qualified_roots=len(density_rows),density_errors=len(density_errors),unqualified_returns=len(failed_iterates),
        density_unavailable=len(density_unavailable),phase_selection_differences=len(selection_differences),failures=len(failures)),indent=2))
    if failures:
        print(json.dumps(failures[:8],indent=2));raise SystemExit(1)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory');validate(parser.parse_args().compiler_directory)
