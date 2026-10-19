"""Compare raw MBWR_density iterates and phase flags, including bad returns."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe
from generate_data import compiled_components
from generate_mbwr_data import records,stored_properties

ABSOLUTE=[1e-12,0.,1e-5,1e-8]
RELATIVE=2e-9
QUANTITIES=['density','phase_found','pressure','pressure_rho']

def cases():
    general,_=compiled_components(ROOT/'src/compdatadb.f90')
    tables={(int(kind),id):v for kind,id,v in stored_properties() if kind!='triple'}
    for record in records():
        stored=tables.get((record['kind'],record['identifier']));component=general.get(record['identifier'])
        if stored:critical=[stored[0],stored[1],stored[2]*1000.]
        elif component:critical=[component['Tc'],component['Pc'],component['Pc']/(component['Zc']*8.31446261815324*component['Tc'])]
        else:critical=[300.,5e6,10000.]
        acf=component['acf'] if component else .1
        points=[(tr,critical[1]*pr) for tr in [.65,.9,.9999,1.0001,1.2,2.] for pr in [.01,.1,.8,1.,1.2,3.]]
        points += [(tr,p) for tr in [.7,1.2] for p in [1e-4,.01,99.,100.]]
        points += [(tr,critical[1]*pr) for tr in [.9499,.9501,1.0099,1.0101,1.0499,1.0501] for pr in [.9499,.9501,1.0499,1.0501]]
        for tr,p in points:
            for phase in [1,2]:
                yield dict(mode='MBWR_DENSITY',record_index=int(record['source_name'][6:]),source_name=record['source_name'],
                    temperature=critical[0]*tr,pressure=p,phase=phase,critical=critical,acentric_factor=acf,
                    explicit_synthetic_metadata=component is None)

def encode(request):
    config=','.join(map(str,[*request['critical'],request['acentric_factor']]))
    return '|'.join(['MBWR',request['source_name'],str(request['temperature']),str(request['pressure']),
        'L' if request['phase']==1 else 'V','MBWR_DENSITY@'+config])

def validate(directory):
    probe=verify_probe();rows=[];errors=[];requests=[]
    with ReferenceProcess(directory) as oracle:
        for request in cases():
            result=oracle.call(request)
            if 'values' not in result or 'calculation_error' in result or 'reference_error' in result:
                errors.append(dict(request=request,result=result));continue
            requests.append(request);rows.append([encode(request),*result['values']])
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    output=subprocess.run([str(exe)],input='\n'.join(r[0] for r in rows)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(output)==len(rows)
    failures=[];paths={};maximum=[0.]*4
    for row,request,line in zip(rows,requests,output):
        if line.startswith('ERROR'):
            failures.append(dict(case=row[0],request=request,reference=row[1:],error=line));continue
        actual=list(map(float,line.split(',')));assert len(actual)==5
        paths[int(actual[4])]=paths.get(int(actual[4]),0)+1
        bad=[]
        for i,(a,b) in enumerate(zip(row[1:],actual)):
            maximum[i]=max(maximum[i],abs(a-b))
            if not math.isfinite(b) or abs(a-b)>ABSOLUTE[i]+(0. if i==1 else RELATIVE)*max(abs(a),abs(b)):bad.append(i)
        if bad:failures.append(dict(case=row[0],request=request,quantities=bad,reference=row[1:],rust=actual))
    if errors:
        outputs=subprocess.run([str(exe)],input='\n'.join(encode(e['request']) for e in errors)+'\n',
            text=True,capture_output=True,check=True).stdout.splitlines()
        assert len(outputs)==len(errors)
        for error,line in zip(errors,outputs):
            error['native_result']=line
            if not line.startswith('ERROR'):failures.append(dict(reference_error=error,native_result=line))
    path=ROOT/'tests/fixtures/mbwr_density.csv'
    with path.open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter=';',lineterminator='\n');writer.writerow(['case',*QUANTITIES]);writer.writerows(rows)
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update(cases=len(rows),attempted_cases=len(rows)+len(errors),probe=probe,reference_errors=errors,
        generator='tools/rust_port/validate_mbwr_density.py',fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        quantities=QUANTITIES,absolute_tolerances=ABSOLUTE,relative_tolerance=RELATIVE,phase_flag_tolerance=0,
        maximum_absolute_differences=maximum,native_path_counts=paths,
        synthetic_initialization_cases=sum(r['explicit_synthetic_metadata'] for r in requests),
        source_pressure_not_closed=sum(abs(row[3]-r['pressure'])>1e-3+1e-6*max(abs(row[3]),r['pressure']) for row,r in zip(rows,requests)),
        convention='all finite original raw iterates including unclosed/unstable returns; explicit critical/acentric initialization; 17 legacy fits use synthetic acentric factor 0.1, with synthetic critical metadata where neither source table nor component data exists; source meta_extrem absent')
    path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/mbwr_density_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(finite_cases=len(rows),reference_errors=len(errors),native_paths=paths,failures=len(failures)),indent=2))
    if failures:print(json.dumps(failures[:8],indent=2));raise SystemExit(1)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory');validate(parser.parse_args().compiler_directory)
