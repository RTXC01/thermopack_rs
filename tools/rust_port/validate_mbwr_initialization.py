"""Compare original default MBWR initialization and its coarse critical grid."""
import argparse,csv,hashlib,json,math,os,re,subprocess
from reference import ROOT,ReferenceProcess
from generate_mbwr_data import records,stored_properties
from validate_alpha import verify_probe

INITIAL_ATOL=[1e-8,1e-3,1e-6,1e-10,1e-14,1e-14,1e-14,1e-12]
GRID_ATOL=[1e-8,1e-3,1e-6]
RTOL=2e-10
ADDRESS=re.compile(r'0x[0-9a-fA-F]+')

def normalize_output(output):
    """Drop ASLR-dependent load addresses from recorded source backtraces.

    The frame count and the final diagnostic message are retained; only the
    per-frame instruction addresses are replaced so that regenerating the
    fixture in another process produces identical metadata.
    """
    return [ADDRESS.sub('0xADDR',line) for line in output]

def normalize_result(result):
    """Apply `normalize_output` to every recorded source output list."""
    if isinstance(result,dict):
        return {k:(normalize_output(v) if k=='output' and isinstance(v,list) else normalize_result(v)) for k,v in result.items()}
    if isinstance(result,list):
        return [normalize_result(v) for v in result]
    return result

def encode(record,mode):return '|'.join(['MBWR',record['source_name'],'0','0','0',mode])
def native(cases):
    if not cases:return []
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    lines=subprocess.run([str(exe)],input='\n'.join(cases)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(lines)==len(cases)
    return [s if s.startswith('ERROR') else list(map(float,s.split(','))) for s in lines]
def fixture(name,rows,quantities,probe,**extra):
    path=ROOT/'tests/fixtures'/name
    with path.open('w',newline='') as f:
        w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(['case',*quantities]);w.writerows(rows)
    meta=json.loads((ROOT/'.reference/metadata.json').read_text())
    meta.update(cases=len(rows),probe=probe,generator='tools/rust_port/validate_mbwr_initialization.py',
        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),quantities=quantities,relative_tolerance=RTOL,**extra)
    path.with_suffix('.metadata.json').write_text(json.dumps(meta,indent=2)+'\n')

def validate(directory):
    probe=verify_probe();initial=[];initial_errors=[];undefined=[];grid=[];grid_errors=[]
    tables={(int(k),id):v for k,id,v in stored_properties() if k!='triple'}
    with ReferenceProcess(directory) as oracle:
        for record in records():
            request=dict(mode='MBWR_INITIALIZATION',record_index=int(record['source_name'][6:]))
            result=oracle.call(request);key=encode(record,request['mode'])
            if 'undefined_source_initialization' in result:
                undefined.append(dict(case=key,request=request,record=record['identifier'],result=normalize_result(result)));continue
            if 'values' not in result or 'calculation_error' in result:
                initial_errors.append(dict(case=key,request=request,result=normalize_result(result)));continue
            initial.append([key,*result['values']])
        for record in records():
            table=tables.get((record['kind'],record['identifier']))
            if not table:continue
            offsets=[(0.,0.)]
            if record['source_name'] in ['mbwrcx1','mbwrcx58']:
                offsets +=[(t,r) for t in [-.05,.05] for r in [-50.,50.]]
            for dt,dr in offsets:
                values=[table[0]+dt,table[1],table[2]*1000.+dr]
                request=dict(mode='MBWR_CRITICAL_GRID',record_index=int(record['source_name'][6:]),initial=values)
                key=encode(record,'MBWR_CRITICAL_GRID@'+','.join(map(str,values)))
                result=oracle.call(request)
                if 'values' not in result or 'calculation_error' in result:
                    grid_errors.append(dict(case=key,request=request,result=normalize_result(result)));continue
                grid.append([key,*result['values']])
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    failures=[]
    for rows,atol in [(initial,INITIAL_ATOL),(grid,GRID_ATOL)]:
        for row,actual in zip(rows,native([r[0] for r in rows])):
            if isinstance(actual,str) or len(actual)!=len(row)-1:
                failures.append(dict(case=row[0],reference=row[1:],rust=actual));continue
            bad=[i for i,(a,b) in enumerate(zip(row[1:],actual)) if not math.isfinite(b) or abs(a-b)>atol[i]+RTOL*max(abs(a),abs(b))]
            if bad:failures.append(dict(case=row[0],quantities=bad,reference=row[1:],rust=actual))
    for errors in [initial_errors,undefined,grid_errors]:
        for error,result in zip(errors,native([e['case'] for e in errors])):
            error['native_result']=result
            if not isinstance(result,str):failures.append(dict(reference_error=error,rust=result))
    fixture('mbwr_initialization.csv',initial,['Tc','Pc','rhoc','Zc','acentric_factor','b_srk','m_srk','a0_srk'],probe,
        absolute_tolerances=INITIAL_ATOL,reference_errors=initial_errors,undefined_initializers=undefined,
        attempted_records=65,convention='original initializeMBWRmodel/readDbParameters, including table selection and coarse critical grid; unknown component IDs are diagnosed before original compdb(0) access; SI output for density and SRK b/a0')
    statuses=[[row[0],0] for row in initial]+[[e['case'],1] for e in initial_errors]+[[e['case'],2] for e in undefined]
    statuses.sort(key=lambda r:int(r[0].split('|')[1][6:]))
    fixture('mbwr_initialization_status.csv',statuses,['source_status'],probe,absolute_tolerances=[0.],
        convention='0 finite original initializer, 1 original Fortran termination, 2 undefined original compdb(0) access diagnosed before execution; native initialization must return a value only for status 0')
    fixture('mbwr_critical_grid.csv',grid,['Tc','Pc','rhoc'],probe,absolute_tolerances=GRID_ATOL,reference_errors=grid_errors,
        convention='original MBWR_criticalParameters with explicit SI initial T/P/rho, converted to backend mol/L; 0.001 K and 1 mol/m3 increments, original loop order and stopping rules')
    (ROOT/'.reference/mbwr_initialization_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(default_initializations=len(initial),reference_terminations=len(initial_errors),
        undefined_initializers=len(undefined),grid_results=len(grid),grid_errors=len(grid_errors),failures=len(failures)),indent=2))
    if failures:print(json.dumps(failures[:8],indent=2));raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');validate(p.parse_args().compiler_directory)
