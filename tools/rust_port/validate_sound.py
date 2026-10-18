"""Qualify equilibrium and TV sound speed against this exact ThermoPack fork."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT,ReferenceProcess
from validate_response import cases as response_cases
from validate_translation import encode as encode_model


def encode(request,result):
    r=dict(request)
    if r['mode']=='V_SOUND_TV':r['pressure']=result['volume']
    if 'quality' in r:r['pressure']=result['pressure']
    text=encode_model(r)
    if 'quality' in r:text+='@quality:'+str(r['quality'])
    if 'roots' in r:text+='@roots:'+':'.join(r['roots'])
    if r['mode']=='V_SOUND_TV':text+='@amount:'+str(r.get('amount',1.))
    return text


def cases():
    for request in response_cases():
        yield {**request,'mode':'V_SOUND'}
        for phase in ['L','V']:
            yield {**request,'mode':'V_SOUND_TV','phase':phase,'amount':3.5}
    for eos in ['PR','SRK','RK','VDW']:
        for ref in ['Default'] if eos=='VDW' else ['Default','tcPR']:
            for comps,t in [('C1',150.),('CO2',250.),('C3',300.),('C1,C2,N2',150.)]:
                z=[1.] if ',' not in comps else [1.,0.,0.]
                for quality in [0.,.01,.5,.99,1.]:
                    yield dict(mode='V_SOUND',eos=eos,components=comps,feed=z,temperature=t,pressure=0.,quality=quality,parameter_reference=ref)
    for eos in ['PR','SRK','RK']:
        for comps,z in [('HE',[1.]),('H2',[1.]),('NE',[1.]),('H2,NE',[.6,.4]),('HE,H2,NE',[0.,.5,.5]),('C1,CO2,N2',[.5,.3,.2])]:
            for t in [3.,20.,80.,300.]:
                for p in [1e3,1e6]:
                    request=dict(eos=eos,components=comps,feed=z,temperature=t,pressure=p,parameter_reference='QuantumCubic')
                    yield dict(mode='V_SOUND',**request)
                    for phase in ['L','V']:yield dict(mode='V_SOUND_TV',phase=phase,amount=.3,**request)
        for comps,t in [('HE',3.),('H2',20.),('NE',30.)]:
            for quality in [0.,.01,.5,.99,1.]:
                yield dict(mode='V_SOUND',eos=eos,components=comps,feed=[1.],temperature=t,pressure=0.,quality=quality,parameter_reference='QuantumCubic')
    for eos in ['PR','SRK','RK','VDW']:
        for t in [180.,250.,350.]:
            for roots in [['L','L'],['V','V'],['V','L']]:
                yield dict(mode='V_SOUND',eos=eos,components='C1,C2,C3',feed=[.5,.3,.2],temperature=t,pressure=1e6,roots=roots)
        for name,t in [('C1',100.),('CO2',180.)]:
            for v in [5e-5,6e-5,1e-4,2e-4,1e-3]:
                yield dict(mode='V_SOUND_TV',eos=eos,components=name,feed=[1.],temperature=t,volume=v,pressure=0.)


    for eos in ['PR','SRK','RK','VDW']:
        operations=[] if eos=='VDW' else [['init','Default']]+[['set',i,6,[2e-6,1e-8,-1e-11,1e-14,-1e-17,1e-20]] for i in range(3)]
        operations += [['beta',i,'QUANTUM',[3.+i,12.]] for i in range(3)]+[['kij',0,1,.07],['lij',1,2,.04]]
        for t in [150.,300.,700.]:
            for p in [1e3,1e6]:
                request=dict(eos=eos,components='C1,C2,C3',feed=[.5,.3,.2],temperature=t,pressure=p,operations=operations)
                yield dict(mode='V_SOUND',**request)
                for phase in ['L','V']:yield dict(mode='V_SOUND_TV',phase=phase,amount=7.3,**request)


def validate(directory):
    inputs=[];expected=[];errors=[];modes=[];nonpositive_pressure=0;nonpositive_positive_sound=0
    with ReferenceProcess(directory) as oracle:
        for request in cases():
            result=oracle.call(request)
            if 'values' not in result or 'calculation_error' in result or not all(map(math.isfinite,result['values'])):
                errors.append({'request':request,**result});continue
            inputs.append(encode(request,result));expected.append(result['values']);modes.append(request['mode'])
            if request['mode']=='V_SOUND_TV' and result['pressure']<=0.:
                nonpositive_pressure+=1
                if result['values'][0]>0.:nonpositive_positive_sound+=1
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    actual=subprocess.run([str(exe)],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(actual)==len(inputs);failures=[];groups={mode:[] for mode in set(modes)}
    for case,reference,line,mode in zip(inputs,expected,actual,modes):
        groups[mode].append((case,reference))
        if line.startswith('ERROR'):
            failures.append(dict(case=case,error=line,reference=reference));continue
        native=list(map(float,line.split(',')));assert len(native)==len(reference)
        for i,(a,b) in enumerate(zip(reference,native)):
            tolerance=0. if mode=='V_SOUND' and i==0 else 1e-5+2e-7*max(abs(a),abs(b))
            if not abs(a-b)<=tolerance:failures.append(dict(case=case,quantity=i,reference=a,rust=b))
    for mode,name in [('V_SOUND','equilibrium_sound'),('V_SOUND_TV','tv_sound')]:
        path=ROOT/'tests/fixtures'/(name+'.csv')
        with path.open('w',newline='') as handle:
            writer=csv.writer(handle,delimiter=';',lineterminator='\n');writer.writerow(['case','phase,sound' if mode=='V_SOUND' else 'sound'])
            writer.writerows([[case,*values] for case,values in groups[mode]])
        metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
        failures_for_mode=[error for error in errors if error['request']['mode']==mode]
        metadata.update(cases=len(groups[mode]),attempted_cases=len(groups[mode])+len(failures_for_mode),reference_errors=failures_for_mode,
                        absolute_tolerance=1e-5,relative_tolerance=2e-7,phase_flag_tolerance=0,
                        convention='source Cp/Cv and stiffness/fallback conventions',
                        tv_total_mole_numbers=[.3,1.,3.5,7.3] if mode=='V_SOUND_TV' else None,
                        nonpositive_pressure_cases=nonpositive_pressure if mode=='V_SOUND_TV' else None,
                        nonpositive_pressure_positive_sound_cases=nonpositive_positive_sound if mode=='V_SOUND_TV' else None,
                        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),generator='tools/rust_port/validate_sound.py')
        path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/sound_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(cases={k:len(v) for k,v in groups.items()},reference_errors=len(errors),failures=len(failures)),indent=2))
    if failures:print(json.dumps(failures[:8],indent=2));raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
