"""Compare native TV caloric properties and UV/SV flashes with the exact fork.

UV exercises both the zero-initialized Python wrapper and explicit guesses.
SV calls the public Fortran dispatcher with TP-initialized phase arrays and
perturbed T/P because this fork comments out its Python SV wrapper.
Every finite output is retained, including failed specification closure.
"""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT,ReferenceProcess
from validate_alpha import verify_probe
from validate_translation import encode as encode_model


def cases():
    for eos in ['PR','SRK','RK','VDW']:
        for ref in ['Default'] if eos=='VDW' else ['Default','tcPR']:
            for comps,z in [('C1',[1.]),('CO2',[1.]),('H2',[1.]),('C1,C2,C3',[.5,.3,.2]),
                            ('C1,CO2,N2',[0.,.8,.2])]:
                for t,p in [(180.,1e6),(300.,1e6),(200.,8e6),(500.,1e3)]:
                    for mode,guessed in [('V_UV',False),('V_UV',True),('V_SV',True)]:
                        yield dict(mode=mode,eos=eos,components=comps,feed=z,temperature=t,pressure=p,
                                   parameter_reference=ref,guessed=guessed)
            for comps,t in [('C1',150.),('CO2',250.),('C1,C2,N2',150.)]:
                z=[1.] if ',' not in comps else [1.,0.,0.]
                for quality in [.1,.5,.9]:
                    for mode in ['V_UV','V_SV']:
                        yield dict(mode=mode,eos=eos,components=comps,feed=z,temperature=t,pressure=0.,
                                   parameter_reference=ref,quality=quality)
            for t in [150.,300.,700.]:
                for v in [8e-5,1e-4,1e-3]:
                    yield dict(mode='V_CALORIC_TV',eos=eos,components='C1,C2,C3',feed=[.5,.3,.2],
                               temperature=t,pressure=0.,volume=v,parameter_reference=ref)
        operations=[] if eos=='VDW' else [['init','Default']]+[['set',i,6,[2e-6,1e-8,-1e-11,1e-14,-1e-17,1e-20]] for i in range(3)]
        operations += [['beta',i,'QUANTUM',[3.+i,12.]] for i in range(3)]+[['kij',0,1,.07],['lij',1,2,.04]]
        for t in [150.,300.,700.]:
            for v in [8e-5,1e-3]:
                yield dict(mode='V_CALORIC_TV',eos=eos,components='C1,C2,C3',feed=[.5,.3,.2],
                           temperature=t,pressure=0.,volume=v,operations=operations)
    for eos in ['PR','SRK','RK']:
        for comp in ['HE','H2','NE']:
            for t in [3.,20.,300.]:
                for v in [1e-3,1e-2]:
                    yield dict(mode='V_CALORIC_TV',eos=eos,components=comp,feed=[1.],temperature=t,
                               pressure=0.,volume=v,parameter_reference='QuantumCubic')


def encode(request,result):
    r=dict(request)
    if request['mode']=='V_CALORIC_TV':r['pressure']=result['volume']
    else:r.update(temperature=result['specification'],pressure=result['volume'])
    line=encode_model(r)
    if result.get('temperature_guess') is not None:
        line+='@guess:'+str(result['temperature_guess'])+':'+str(result['pressure_guess'])
    return line


def tolerance(index,mode,a,b):
    if mode=='V_CALORIC_TV':return [1e-4,1e-6,1e-8][index]+2e-9*max(abs(a),abs(b))
    if index==2:return 0.
    return ([2e-5,.5,0.,2e-7,2e-7]+[2e-7]*50)[index]+2e-7*max(abs(a),abs(b))


def validate(directory):
    probe=verify_probe();records=[];errors=[]
    with ReferenceProcess(directory) as oracle:
        for request in cases():
            result=oracle.call(request)
            if 'calculation_error' in result or 'values' not in result or not all(map(math.isfinite,result['values'])):
                errors.append(dict(request=request,**result));continue
            records.append(dict(case=encode(request,result),reference=result['values'],request=request,**{k:v for k,v in result.items() if k!='values'}))
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    lines=subprocess.run([str(exe)],input='\n'.join(r['case'] for r in records)+'\n',text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(lines)==len(records)
    failures=[]
    for record,line in zip(records,lines):
        mode=record['request']['mode']
        if line.startswith('ERROR'):
            failures.append(dict(case=record['case'],error=line,reference=record['reference']));continue
        actual=list(map(float,line.split(',')));expected=record['reference'];assert len(actual)==len(expected)
        bad=[i for i,(a,b) in enumerate(zip(expected,actual)) if not abs(a-b)<=tolerance(i,mode,a,b)]
        if bad:failures.append(dict(case=record['case'],fields=bad,reference=expected,rust=actual))
        if mode!='V_CALORIC_TV' and (not all(map(math.isfinite,[record['scaled_specification_error'],record['relative_volume_error']]))
                                    or abs(record['scaled_specification_error'])>1e-7 or abs(record['relative_volume_error'])>1e-7):
            failures.append(dict(case=record['case'],error='source returned state fails specification closure',
                                 scaled_specification_error=record['scaled_specification_error'],relative_volume_error=record['relative_volume_error']))
    for modes,name in [(['V_CALORIC_TV'],'caloric_tv'),(['V_UV','V_SV'],'isochoric_flash')]:
        group=[r for r in records if r['request']['mode'] in modes];path=ROOT/'tests/fixtures'/(name+'.csv')
        with path.open('w',newline='') as handle:
            writer=csv.writer(handle,delimiter=';',lineterminator='\n');writer.writerow(['case','P,U,S' if name=='caloric_tv' else 'T,P,two_phase,betaV,betaL,x,y'])
            writer.writerows([[r['case'],*r['reference']] for r in group])
        metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
        metadata.update(cases=len(group),reference_errors=[e for e in errors if e['request']['mode'] in modes],
                        attempted_cases=sum(r['request']['mode'] in modes for r in records)+sum(e['request']['mode'] in modes for e in errors),
                        probe=probe,generator='tools/rust_port/validate_isochoric.py',
                        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        convention='molar TV pressure, internal energy and entropy' if name=='caloric_tv' else 'single-phase output canonicalized to x=z, y=0, betaV=0, betaL=1; two-phase arrays retained',
                        tolerances={'TV':{'absolute':[1e-4,1e-6,1e-8],'relative':2e-9},
                                    'flash':{'absolute':[2e-5,.5,0.,2e-7,2e-7],'composition_absolute':2e-7,'relative':2e-7,'two_phase_relative':0}},
                        raw_phase_codes=sorted({r['raw_result'][2] for r in group if 'raw_result' in r}),
                        maximum_source_scaled_specification_error=max((abs(r.get('scaled_specification_error',0)) for r in group),default=0),
                        maximum_source_relative_volume_error=max((abs(r.get('relative_volume_error',0)) for r in group),default=0),
                        comparisons=[f for f in failures if f['case'] in {r['case'] for r in group}])
        path.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/isochoric_records.json').write_text(json.dumps(records,indent=2)+'\n')
    (ROOT/'.reference/isochoric_failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps(dict(finite_cases=len(records),reference_errors=len(errors),failures=len(failures)),indent=2))
    if failures:print(json.dumps(failures[:6],indent=2));raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
