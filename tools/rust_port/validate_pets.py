"""Qualify PeTS density/extensive Helmholtz derivatives, TV values and potential."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe


def cases():
    for sigma,eps in [(3.42e-10,136.),(2.9e-10,100.),(4.1e-10,240.),(3.42e-10,1.)]:
        for tr in [.02,.3,.7,1.,1.5,3.,10.]:
            t=tr*eps;d=sigma*(1.-.127112544*math.exp(-3.052785558*eps/t))
            for eta in [1e-10,1e-6,.05,.2,.4,.65,.85]:
                rho=eta/(math.pi/6.*d**3*6.02214076e23)
                for n in [.3,1.,3.5]:
                    yield dict(mode='PETS_CORE',parameters=[sigma,eps],temperature=t,volume=n/rho,moles=n)
        for reduced in [[.5,.9,1.,1.1,2.,2.49999,2.5,3.,4.], [3.,1.,.8], [1.,3.,.8], [.7,1.2,1.1,2.]]:
            yield dict(mode='PETS_POTENTIAL',parameters=[sigma,eps],distances=[r*sigma for r in reduced])


def encode(request):
    params=','.join(map(str,request['parameters']))
    if request['mode']=='PETS_CORE':
        return '|'.join(['PETS',params,*map(str,[request['temperature'],request['volume'],request['moles']]),request['mode']])
    return '|'.join(['PETS',params,','.join(map(str,request['distances'])),'0','0',request['mode']])


def absolute_tolerances(mode):
    if mode=='PETS_POTENTIAL':return None
    return [1e-22,1e-23,1e-24]+[1e-11,1e-11,1e-16,1e-10,1e-15,1e-20]*3+[
        1e-10,1e-10,1e-8,1e-10,1e-9,1e-7,1e-9,1e-3,1e-7,1e-9]*2+[
        1e-3,1e-4,1.,1e-3,1e-5,1e-7,1e-5,1e-5,1e-6,1e-5,1e-6,1e-15,1e-11,1e-9,1e-6]


def validate(directory):
    probe=verify_probe();records=[];errors=[]
    with ReferenceProcess(directory) as oracle:
        for request in cases():
            result=oracle.call(request)
            if 'calculation_error' in result or 'values' not in result:
                errors.append(dict(request=request,**result));continue
            records.append((request,encode(request),result['values']))
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    lines=subprocess.run([str(exe)],input='\n'.join(case for _,case,_ in records)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(lines)==len(records);failures=[]
    for (request,case,values),line in zip(records,lines):
        if line.startswith('ERROR'):
            failures.append(dict(case=case,error=line));continue
        actual=list(map(float,line.split(',')));assert len(actual)==len(values)
        atol=absolute_tolerances(request['mode'])
        for i,(a,b) in enumerate(zip(values,actual)):
            absolute=atol[i] if atol else (1e-10 if i<len(values)//2 else 1e-3)
            if not abs(a-b)<=absolute+5e-10*max(abs(a),abs(b)):
                failures.append(dict(case=case,quantity=i,reference=a,rust=b))
    for mode,name in [('PETS_CORE','pets_core'),('PETS_POTENTIAL','pets_potential')]:
        group=[(case,values) for request,case,values in records if request['mode']==mode]
        path=ROOT/'tests/fixtures'/(name+'.csv')
        with path.open('w',newline='') as handle:
            writer=csv.writer(handle,delimiter=';',lineterminator='\n');writer.writerow(['case','diameter,alpha,alphaHS,alphaDisp,F,tail,P,U,S,A,H,Cv,mu,sound,B,C,lnphi,Cp' if mode=='PETS_CORE' else 'potential,force'])
            writer.writerows([[case,*values] for case,values in group])
        metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
        metadata.update(cases=len(group),attempted_cases=len(group)+sum(e['request']['mode']==mode for e in errors),
                        reference_errors=[e for e in errors if e['request']['mode']==mode],probe=probe,
                        absolute_tolerances=absolute_tolerances(mode),potential_absolute_tolerance=1e-10,force_absolute_tolerance=1e-3,
                        relative_tolerance=5e-10,fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        quantity_order=(['diameter','diameter_T','diameter_TT']+
                                        [family+suffix for family in ['alpha','alpha_hs','alpha_disp'] for suffix in ['', '_T','_rho','_TT','_Trho','_rhorho']]+
                                        [family+suffix for family in ['F','correction'] for suffix in ['', '_T','_V','_n','_TT','_TV','_Tn','_VV','_Vn','_nn']]+
                                        ['P','P_T','P_V','P_n','U','S','A','H','Cv','mu','sound','B_finite_V','C_finite_V','lnphi_positive_pressure','Cp']) if mode=='PETS_CORE' else None,
                        nonpositive_pressure_cases=sum(values[41]<=0. for _,values in group) if mode=='PETS_CORE' else None,
                        zero_sound_cases=sum(values[51]==0. for _,values in group) if mode=='PETS_CORE' else None,
                        generator='tools/rust_port/validate_pets.py',
                        convention='actual Fortran setter; alpha density and extensive F derivatives; source first-cutoff vector termination')
        path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/pets_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(finite_cases=len(records),reference_errors=len(errors),failures=len(failures)),indent=2))
    if failures:print(json.dumps(failures[:8],indent=2));raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
