"""Line-oriented subprocess for a library that can terminate its caller."""
import argparse
import json
import sys
import warnings
import math
import os
from ctypes import c_double,c_int,c_char_p,CDLL,POINTER,byref
from reference import ROOT,load


def overall(model,flash,method):
    t,p,z=flash.T,flash.p,flash.z
    calculate=getattr(model,method)
    if flash.phase==model.TWOPH:
        return (flash.betaL*calculate(t,p,flash.x,model.LIQPH)[0]
                +flash.betaV*calculate(t,p,flash.y,model.VAPPH)[0])
    return calculate(t,p,z,model.LIQPH if flash.phase==model.SINGLEPH else flash.phase)[0]


def phase_permutation_allowed(model,result):
    if result.phase!=model.TWOPH:
        return False
    for composition in [result.x,result.y]:
        vl=model.specific_volume(result.T,result.p,composition,model.LIQPH)[0]
        vv=model.specific_volume(result.T,result.p,composition,model.VAPPH)[0]
        if abs(vl-vv)>1e-10*max(abs(vl),abs(vv)):
            return False
    return True


def load_probe():
    probe=CDLL(str(ROOT/'.reference'/('reference_probe.dll' if os.name=='nt' else 'libreference_probe.so')))
    probe.query_alpha.argtypes=[c_int,c_double,POINTER(c_int),POINTER(c_double),POINTER(c_double)]
    probe.query_alpha.restype=None
    probe.initialize_alpha.argtypes=[c_char_p,c_int,c_char_p,c_int]
    probe.initialize_alpha.restype=None
    probe.query_translation.argtypes=[c_int,c_double,POINTER(c_int),POINTER(c_int),POINTER(c_double),POINTER(c_double),POINTER(c_int)]
    probe.query_translation.restype=None
    probe.query_beta.argtypes=[c_int,c_double,POINTER(c_int),POINTER(c_double),POINTER(c_double),POINTER(c_double)]
    probe.query_beta.restype=None
    probe.query_rr.argtypes=[c_int,POINTER(c_double),POINTER(c_double),c_int,POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_int)]
    probe.query_rr.restype=None
    probe.query_rr_iteration.argtypes=[c_int,c_double,c_double,POINTER(c_double),POINTER(c_double),c_int,
                                      *([POINTER(c_double)]*8),POINTER(c_int)]
    probe.query_rr_iteration.restype=None
    probe.query_symmetric.argtypes=[c_int,POINTER(c_double),POINTER(c_double),POINTER(c_int),POINTER(c_int)]
    probe.query_symmetric.restype=None
    probe.query_response.argtypes=[c_int,c_double,c_double,POINTER(c_double),c_double,c_double,
                                   POINTER(c_double),POINTER(c_double),c_int,POINTER(c_double)]
    probe.query_response.restype=None
    probe.query_sv.argtypes=[c_int,*([POINTER(c_double)]*9),POINTER(c_int),POINTER(c_int)]
    probe.query_sv.restype=None
    probe.query_pets.argtypes=[*([c_double]*5),POINTER(c_double)]
    probe.query_pets.restype=None
    probe.query_pets_potential.argtypes=[c_int,POINTER(c_double),c_double,c_double,POINTER(c_double),POINTER(c_double)]
    probe.query_pets_potential.restype=None
    probe.set_pets_parameters.argtypes=[c_double,c_double]
    probe.set_pets_parameters.restype=None
    probe.query_pets_volume.argtypes=[*([c_double]*4),c_int,POINTER(c_double)]
    probe.query_pets_volume.restype=None
    return probe


def main(compiler_directory):
    cubic=load(compiler_directory)
    model,previous=None,None
    probe=None
    for line in sys.stdin:
        request=json.loads(line)
        if request['mode']=='SYMMETRIC':
            if probe is None: probe=load_probe()
            n=len(request['rhs'])
            a=(c_double*(n*n))(*request['matrix']);b=(c_double*n)(*request['rhs'])
            pivots=(c_int*n)();status=c_int()
            probe.query_symmetric(n,a,b,pivots,byref(status))
            values=[float(status.value),*map(float,pivots),*map(float,a),*map(float,b)]
            if all(map(math.isfinite,values)):
                print('THERMOPACK_RESULT='+json.dumps({'values':values}),flush=True)
            else:
                print('THERMOPACK_RESULT='+json.dumps({'calculation_error':'non-finite symmetric factor/solution',
                      'values':list(map(str,values))}),flush=True)
            continue
        if request['mode']=='RR':
            if probe is None: probe=load_probe()
            n=len(request['feed'])
            z=(c_double*n)(*request['feed']);k=(c_double*n)(*request['k'])
            x=(c_double*n)();y=(c_double*n)()
            beta=c_double(request.get('beta',.5));betal=c_double(.5);status=c_int()
            probe.query_rr(n,z,k,int(request.get('single_step',False)),byref(beta),byref(betal),x,y,byref(status))
            values=[1.,beta.value,betal.value,*map(float,x),*map(float,y)] if status.value else [0.]
            if all(map(math.isfinite,values)):
                print('THERMOPACK_RESULT='+json.dumps({'values':values}),flush=True)
            else:
                print('THERMOPACK_RESULT='+json.dumps({'calculation_error':'non-finite RR output','values':list(map(str,values))}),flush=True)
            continue
        if request['mode'].startswith('PETS_'):
            sigma,eps=request.get('parameters',[3.42e-10,136.])
            key=('PETS',sigma,eps,request.get('redefine_critical',False))
            if previous!=key:
                from thermopack.pets import pets
                model=pets()
                if probe is None:probe=load_probe()
                model.activate();probe.set_pets_parameters(sigma,eps)
                if request.get('redefine_critical',False):
                    model.redefine_critical_parameters(Tc_initials=[1.1*eps],vc_initials=[math.pi/6.*sigma**3*6.02214076e23/.15])
                previous=key
            if probe is None:probe=load_probe()
            model.activate();sigma,eps=request.get('parameters',[3.42e-10,136.])
            if request['mode'] in ('PETS_UV','PETS_SV'):
                from types import SimpleNamespace
                probe.set_pets_parameters(sigma,eps)
                model.set_tmin(request.get('minimum_temperature',1.))
                model.set_tmax(request.get('maximum_temperature',2000.))
                try:
                    t,p=request['temperature'],request['pressure']
                    if 'quality' in request:
                        p=model.bubble_pressure(t,[1.])[0]
                        seed=SimpleNamespace(T=t,p=p,z=[1.],x=[1.],y=[1.],betaV=request['quality'],
                                             betaL=1.-request['quality'],phase=model.TWOPH)
                    else:seed=model.two_phase_tpflash(t,p,[1.])
                    volume=overall(model,seed,'specific_volume')
                    entropy=overall(model,seed,'entropy')
                    energy=overall(model,seed,'enthalpy')-p*volume
                    target=energy if request['mode']=='PETS_UV' else entropy
                    guess_t=t*request.get('guess_factor',1.2)
                    guess_p=p*request.get('pressure_factor',.8)
                    if request['mode']=='PETS_UV':
                        result=model.two_phase_uvflash([1.],energy,volume,
                            temp=guess_t if request.get('guessed',True) else None,
                            press=guess_p if request.get('guessed',True) else None)
                        status=0
                    else:
                        tv,pv,bv,bl=c_double(guess_t),c_double(guess_p),c_double(seed.betaV),c_double(seed.betaL)
                        zv=(c_double*1)(1.);xv=(c_double*1)(*seed.x);yv=(c_double*1)(*seed.y)
                        sv,vv=c_double(entropy),c_double(volume);phase=c_int(seed.phase);code=c_int()
                        probe.query_sv(1,byref(tv),byref(pv),zv,byref(bv),byref(bl),xv,yv,
                                       byref(sv),byref(vv),byref(phase),byref(code))
                        status=code.value
                        result=SimpleNamespace(T=tv.value,p=pv.value,z=[1.],x=list(xv),y=list(yv),
                                               betaV=bv.value,betaL=bl.value,phase=phase.value)
                    raw=[float(result.T),float(result.p),int(result.phase),float(result.betaV),float(result.betaL),
                         *map(float,result.x),*map(float,result.y)]
                    two=result.phase==model.TWOPH
                    actual_v=overall(model,result,'specific_volume')
                    h=overall(model,result,'enthalpy');s=overall(model,result,'entropy');u=h-result.p*actual_v
                    caloric_error=((u if request['mode']=='PETS_UV' else s)-target)/(8.31446261815324*(result.T if request['mode']=='PETS_UV' else 1.))
                    volume_error=(actual_v-volume)/volume
                    phases=[model.LIQPH,model.VAPPH] if two else [model.LIQPH if result.phase==model.SINGLEPH else result.phase]
                    closure=[]
                    for phase_id in phases:
                        v=model.specific_volume(result.T,result.p,[1.],phase_id)[0]
                        a,b=model.pressure_tv(result.T,v,[1.],dpdv=True)
                        closure.append(dict(pressure=float(a),dv=float(b)))
                    canonical_two=two and 1e-8 < result.betaV < 1.-1e-8
                    values=[float(result.T),float(result.p),float(canonical_two),float(result.betaV) if canonical_two else 0.,
                            float(actual_v),float(u),float(s),float(h)]
                    extra=dict(specification=target,volume=volume,raw_result=raw,status=status,
                               scaled_specification_error=caloric_error,relative_volume_error=volume_error,pressure_closure=closure,
                               temperature_guess=guess_t if request.get('guessed',True) or request['mode']=='PETS_SV' else None,
                               pressure_guess=guess_p if request.get('guessed',True) or request['mode']=='PETS_SV' else None)
                    if status!=0 or not all(map(math.isfinite,values+[caloric_error,volume_error])) or abs(caloric_error)>1e-7 or abs(volume_error)>1e-7 or any(
                            not math.isfinite(c['pressure']) or not math.isfinite(c['dv']) or abs(c['pressure']-result.p)>1e8*sys.float_info.epsilon*result.p or c['dv']>=0. for c in closure):
                        print('THERMOPACK_RESULT='+json.dumps(dict(calculation_error='source volume flash fails status, finite/stable pressure or U/S/V closure',
                              values=list(map(str,values)),**extra)),flush=True);continue
                    print('THERMOPACK_RESULT='+json.dumps(dict(values=values,**extra)),flush=True);continue
                except Exception as error:
                    if str(error) not in ('bubble_pressure calculation failed',):raise
                    print('THERMOPACK_RESULT='+json.dumps(dict(calculation_error=str(error))),flush=True);continue
            elif request['mode'] in ('PETS_TP','PETS_PH','PETS_PS'):
                from types import SimpleNamespace
                probe.set_pets_parameters(sigma,eps)
                model.set_tmin(request.get('minimum_temperature',1.))
                model.set_tmax(request.get('maximum_temperature',2000.))
                try:
                    t,p=request['temperature'],request['pressure']
                    if 'quality' in request:
                        p=model.bubble_pressure(t,[1.])[0]
                        seed=SimpleNamespace(T=t,p=p,z=[1.],x=[1.],y=[1.],betaV=request['quality'],
                                             betaL=1.-request['quality'],phase=model.TWOPH)
                    else:seed=model.two_phase_tpflash(t,p,[1.])
                    target=None
                    with warnings.catch_warnings(record=True) as messages:
                        warnings.simplefilter('always')
                        if request['mode']=='PETS_TP':result=seed
                        else:
                            method='enthalpy' if request['mode']=='PETS_PH' else 'entropy'
                            target=overall(model,seed,method)
                            result=getattr(model,'two_phase_phflash' if method=='enthalpy' else 'two_phase_psflash')(
                                p,[1.],target,temp=t*request.get('guess_factor',1.2) if request.get('guessed',True) else None)
                    raw=[float(result.T),float(result.p),int(result.phase),float(result.betaV),float(result.betaL),
                         *map(float,result.x),*map(float,result.y)]
                    two=result.phase==model.TWOPH
                    if two:
                        vl=model.specific_volume(result.T,p,[1.],model.LIQPH)[0]
                        vv=model.specific_volume(result.T,p,[1.],model.VAPPH)[0]
                        phase_volumes=[vl,vv]
                    else:
                        v=model.specific_volume(result.T,p,[1.],model.LIQPH if result.phase==model.SINGLEPH else result.phase)[0]
                        vl,vv=v,v;phase_volumes=[v]
                    closure=[model.pressure_tv(result.T,v,[1.],dpdv=True) for v in phase_volumes]
                    h=overall(model,result,'enthalpy');s=overall(model,result,'entropy')
                    actual=h if request['mode']=='PETS_PH' else s
                    scaled_error=0. if target is None else (actual-target)/(8.31446261815324*(result.T if request['mode']=='PETS_PH' else 1.))
                    canonical_two=two and 1e-8 < result.betaV < 1.-1e-8
                    overall_v=result.betaL*vl+result.betaV*vv if two else vl
                    values=[float(result.T),p,float(canonical_two),float(result.betaV) if canonical_two else 0.,
                            float(vl if canonical_two else overall_v),float(vv if canonical_two else overall_v),float(h),float(s)]
                    extra=dict(specification=target,raw_result=raw,scaled_specification_error=scaled_error,
                               pressure_closure=[dict(pressure=float(a),dv=float(b)) for a,b in closure],
                               warnings=[str(m.message) for m in messages])
                    if not all(map(math.isfinite,values)) or abs(scaled_error)>1e-6 or any(
                            not math.isfinite(a) or not math.isfinite(b) or abs(a-p)>1e8*sys.float_info.epsilon*p or b>=0. for a,b in closure):
                        print('THERMOPACK_RESULT='+json.dumps(dict(calculation_error='source flash fails finite/stable pressure or caloric closure',
                              values=list(map(str,values)),**extra)),flush=True);continue
                    print('THERMOPACK_RESULT='+json.dumps(dict(values=values,**extra)),flush=True);continue
                except Exception as error:
                    if str(error) not in ('PH flash calculation failed','PS flash calculation failed','bubble_pressure calculation failed'):
                        raise
                    print('THERMOPACK_RESULT='+json.dumps(dict(calculation_error=str(error))),flush=True);continue
            elif request['mode'] in ('PETS_VOLUME_BASELINE','PETS_VOLUME','PETS_CRITICAL','PETS_SAT','PETS_SAT_T'):
                probe.set_pets_parameters(sigma,eps)
                try:
                    if request['mode']=='PETS_CRITICAL':
                        guess_t=1.1*eps;guess_v=math.pi/6.*sigma**3*6.02214076e23/.15
                        t,v,p=model.critical([1.],temp=guess_t,v=guess_v,tol=1e-9)
                        values=[t,p,v,p*v/(8.31446261815324*t)]
                    elif request['mode']=='PETS_VOLUME_BASELINE':
                        output=(c_double*4)();phase=model.LIQPH if request['phase']=='L' else model.VAPPH
                        probe.query_pets_volume(request['temperature'],request['pressure'],sigma,eps,phase,output)
                        values=list(map(float,output))
                    elif request['mode']=='PETS_VOLUME':
                        t,p=request['temperature'],request['pressure'];phase=model.LIQPH if request['phase']=='L' else model.VAPPH
                        v,vt,vp=model.specific_volume(t,p,[1.],phase,dvdt=True,dvdp=True)
                        actual,pv=model.pressure_tv(t,v,[1.],dpdv=True)
                        if abs(actual-p)>=1e8*sys.float_info.epsilon*p or pv>=0.:
                            print('THERMOPACK_RESULT='+json.dumps({'calculation_error':'source volume is mechanically unstable' if pv>=0. else 'source volume fails pressure closure',
                                  'volume':float(v),'returned_pressure':float(actual),'pressure_dv':float(pv)}),flush=True);continue
                        values=[float(v),float(vt),float(vp)]
                    else:
                        if request['mode']=='PETS_SAT':
                            t=request['temperature'];p=model.bubble_pressure(t,[1.])[0]
                        else:
                            p=request['pressure'];t=model.bubble_temperature(p,[1.])[0]
                        vl=model.specific_volume(t,p,[1.],model.LIQPH)[0]
                        vv=model.specific_volume(t,p,[1.],model.VAPPH)[0]
                        hl=model.enthalpy(t,p,[1.],model.LIQPH)[0];hv=model.enthalpy(t,p,[1.],model.VAPPH)[0]
                        values=[float(t),float(p),float(vl),float(vv),float((hv-hl)/(t*(vv-vl)))]
                except Exception as error:
                    if str(error) not in ('critical calculation failed','bubble_pressure calculation failed','bubble_temperature calculation failed'):
                        raise
                    print('THERMOPACK_RESULT='+json.dumps({'calculation_error':str(error)}),flush=True);continue
            elif request['mode']=='PETS_CORE':
                output=(c_double*56)()
                probe.query_pets(request['temperature'],request['volume'],request['moles'],sigma,eps,output)
                values=list(map(float,output))
            elif request['mode']=='PETS_POTENTIAL':
                n=len(request['distances']);r=(c_double*n)(*request['distances']);out=(c_double*n)();force=(c_double*n)()
                probe.query_pets_potential(n,r,sigma,eps,out,force)
                values=[*map(float,out),*map(float,force)]
            else:raise ValueError('unknown PeTS mode')
            if not all(map(math.isfinite,values)):
                print('THERMOPACK_RESULT='+json.dumps({'calculation_error':'non-finite PeTS output','values':list(map(str,values))}),flush=True)
            else:print('THERMOPACK_RESULT='+json.dumps({'values':values}),flush=True)
            continue
        eos,comps=request['eos'],request['components']
        alpha=request.get('alpha','RK' if eos=='RK' else 'Classic')
        ref=request.get('parameter_reference',request.get('alpha_reference','Default'))
        shift=request.get('volume_shift',False)
        quantum_helper=request.get('quantum_helper',False)
        key=(eos,comps,alpha,ref,shift,quantum_helper,'alpha_kernel' if request['mode'] in ('ALPHA','ALPHA_PROPERTIES') else 'public')
        if previous!=key:
            # ALPHA qualifies tpInitAlphaCorr separately from eoslibinit's
            # forced tcPR translation / QuantumCubic critical-data overrides.
            alpha_mode=request['mode'] in ('ALPHA','ALPHA_PROPERTIES')
            base_ref=('tcPR' if 'BUTANAL' in comps.split(',') else 'Default') if alpha_mode else ref
            if quantum_helper:
                from thermopack.quantum_cubic import qcubic
                model=qcubic(comps)
            else:
                model=cubic(comps,"SRK" if eos=="RK" else eos,
                            alpha='Classic' if alpha_mode else alpha,parameter_reference=base_ref,volume_shift=shift)
            previous=key
        extra={}
        if request['mode'].startswith('V_'):
            for operation in request.get('operations',[]):
                if operation[0]=='init':
                    model.init_peneloux_volume_translation(operation[1])
                elif operation[0]=='set':
                    model.set_ci(operation[1]+1,*operation[3],ci_type=operation[2])
                elif operation[0]=='disable':
                    model.disable_volume_translation()
                elif operation[0]=='beta':
                    model.set_beta_corr(operation[1]+1,operation[2],operation[3])
                elif operation[0] in ('kij','lij'):
                    getattr(model,'set_'+operation[0])(operation[1]+1,operation[2]+1,operation[3])
                elif operation[0]=='bounds':
                    model.set_tmin(operation[1]);model.set_tmax(operation[2])
                else:
                    raise ValueError('unknown translation operation')
                previous=None
        if request['mode']=='V_CALORIC_TV':
            t,v,z=request['temperature'],request['volume'],request['feed']
            values=[float(model.pressure_tv(t,v,z)[0]),float(model.internal_energy_tv(t,v,z)[0]),
                    float(model.entropy_tv(t,v,z)[0])]
            extra={'volume':v}
            swappable=False
        elif request['mode'] in ('V_UV','V_SV'):
            from types import SimpleNamespace
            t,p,z=request['temperature'],request['pressure'],request['feed']
            if 'quality' in request:
                p=model.bubble_pressure(t,z)[0]
                seed=SimpleNamespace(T=t,p=p,z=z,x=z,y=z,betaV=request['quality'],
                                     betaL=1.-request['quality'],phase=model.TWOPH)
            else:
                seed=model.two_phase_tpflash(t,p,z)
            v=overall(model,seed,'specific_volume')
            s=overall(model,seed,'entropy')
            u=overall(model,seed,'enthalpy')-p*v
            spec=s if request['mode']=='V_SV' else u
            guess_t=t*request.get('guess_factor',.9)
            guess_p=p*request.get('pressure_factor',1.2)
            if request['mode']=='V_UV':
                result=model.two_phase_uvflash(z,u,v,
                    temp=guess_t if request.get('guessed',False) else None,
                    press=guess_p if request.get('guessed',False) else None)
                status=0
            else:
                if probe is None:probe=load_probe()
                model.activate();n=model.nc
                tv,pv,bv,bl=c_double(guess_t),c_double(guess_p),c_double(seed.betaV),c_double(seed.betaL)
                zv=(c_double*n)(*z);xv=(c_double*n)(*seed.x);yv=(c_double*n)(*seed.y)
                sv,vv=c_double(s),c_double(v);phase=c_int(seed.phase);code=c_int()
                probe.query_sv(n,byref(tv),byref(pv),zv,byref(bv),byref(bl),xv,yv,
                               byref(sv),byref(vv),byref(phase),byref(code))
                status=code.value
                result=SimpleNamespace(T=tv.value,p=pv.value,z=z,x=list(xv),y=list(yv),
                                       betaV=bv.value,betaL=bl.value,phase=phase.value)
            raw=[float(result.T),float(result.p),int(result.phase),float(result.betaV),float(result.betaL),
                 *map(float,result.x),*map(float,result.y)]
            actual_v=overall(model,result,'specific_volume')
            actual=overall(model,result,'entropy') if request['mode']=='V_SV' else overall(model,result,'enthalpy')-result.p*actual_v
            two=result.phase==model.TWOPH
            values=[float(result.T),float(result.p),float(two),
                    float(result.betaV) if two else 0.,float(result.betaL) if two else 1.,
                    *map(float,result.x if two else z),*map(float,result.y if two else [0.]*len(z))]
            extra={'specification':spec,'volume':v,'raw_result':raw,'status':status,
                   'scaled_specification_error':(actual-spec)/(8.31446261815324*(result.T if request['mode']=='V_UV' else 1.)),
                   'relative_volume_error':(actual_v-v)/v,
                   'temperature_guess':guess_t if request['mode']=='V_SV' or request.get('guessed') else None,
                   'pressure_guess':guess_p if request['mode']=='V_SV' or request.get('guessed') else None}
            if status!=0 or not all(map(math.isfinite,values)):
                print('THERMOPACK_RESULT='+json.dumps({'calculation_error':'non-finite volume flash or nonzero status',
                      'values':list(map(str,values)),**extra}),flush=True)
                continue
            swappable=phase_permutation_allowed(model,result)
        elif request['mode']=='V_SOUND_TV':
            t,z=request['temperature'],request['feed']
            v=request.get('volume')
            if v is None:
                v=model.specific_volume(t,request['pressure'],z,model.LIQPH if request['phase']=='L' else model.VAPPH)[0]
            amount=request.get('amount',1.)
            values=[float(model.speed_of_sound_tv(t,v*amount,[zi*amount for zi in z]))]
            extra={'volume':v,'pressure':float(model.pressure_tv(t,v,z)[0])}
            swappable=False
        elif request['mode']=='V_SOUND':
            t,p,z=request['temperature'],request['pressure'],request['feed']
            if 'quality' in request:
                from types import SimpleNamespace
                p=model.bubble_pressure(t,z)[0]
                result=SimpleNamespace(x=z,y=z,betaV=request['quality'],betaL=1.-request['quality'],phase=model.TWOPH)
                extra={'pressure':p}
            else:
                result=model.two_phase_tpflash(t,p,z)
            if 'roots' in request:
                model.activate();n=model.nc
                tv,pv,bv,blv=c_double(t),c_double(p),c_double(result.betaV),c_double(result.betaL)
                xv=(c_double*n)(*result.x);yv=(c_double*n)(*result.y);zv=(c_double*n)(*z);phase=c_int(result.phase)
                codes={'L':model.LIQPH,'V':model.VAPPH}
                roots=(c_int*2)(*[codes[root] for root in reversed(request['roots'])])
                routine=model.s_sos_sound_velocity_2ph
                routine.argtypes=[POINTER(c_double)]*7+[POINTER(c_int),POINTER(c_int)]
                routine.restype=c_double
                sound=routine(byref(tv),byref(pv),xv,yv,zv,byref(bv),byref(blv),byref(phase),roots)
            else:
                sound=model.speed_of_sound(t,p,result.x,result.y,z,result.betaV,result.betaL,result.phase)
            values=[float(result.phase),float(sound)]
            swappable=False
        elif request['mode']=='V_RESP':
            if probe is None:probe=load_probe()
            t,p,z=request['temperature'],request['pressure'],request['feed']
            result=model.two_phase_tpflash(t,p,z)
            n=model.nc;feed=(c_double*n)(*z);x=(c_double*n)(*result.x);y=(c_double*n)(*result.y)
            derivatives=(c_double*3)()
            model.activate()
            probe.query_response(n,t,p,feed,result.betaV,result.betaL,x,y,result.phase,derivatives)
            values=[float(result.phase),*map(float,derivatives)]
            if not all(map(math.isfinite,values)):
                print('THERMOPACK_RESULT='+json.dumps({'calculation_error':'non-finite equilibrium response',
                      'values':list(map(str,values))}),flush=True)
                continue
            swappable=False
        elif request['mode']=='V_RR':
            if probe is None: probe=load_probe()
            model.activate()
            n=model.nc
            z=(c_double*n)(*request['feed']);k=(c_double*n)(*request['k'])
            x=(c_double*n)();y=(c_double*n)();kout=(c_double*n)()
            fugl=(c_double*n)();fugv=(c_double*n)()
            beta=c_double(request.get('beta',.5));betal=c_double(.5);g=c_double()
            flags=(c_int*2)()
            probe.query_rr_iteration(n,request['temperature'],request['pressure'],z,k,
                                     int(request.get('single_step',False)),byref(beta),byref(betal),
                                     x,y,kout,fugl,fugv,byref(g),flags)
            values=list(map(float,flags))
            if flags[1]:
                values.extend([beta.value,betal.value,g.value,*map(float,x),*map(float,y),
                               *map(float,kout),*map(float,fugl),*map(float,fugv)])
            if not all(map(math.isfinite,values)):
                print('THERMOPACK_RESULT='+json.dumps({'calculation_error':'non-finite RR iteration output',
                      'values':list(map(str,values))}),flush=True)
                continue
            swappable=False
        elif request['mode']=='V_BETA':
            if probe is None:
                probe=load_probe()
            model.activate()
            values=[]
            for index in range(model.nc):
                method=c_int()
                params=(c_double*2)();derivatives=(c_double*3)();critical=(c_double*3)()
                probe.query_beta(index+1,request['temperature'],byref(method),params,derivatives,critical)
                values.extend([method.value,*map(float,params),*map(float,derivatives),*map(float,critical)])
            swappable=False
        elif request['mode']=='V_CI':
            if probe is None:
                probe=load_probe()
            model.activate()
            values=[]
            for index in range(model.nc):
                enabled,kind,dependent=c_int(),c_int(),c_int()
                parameters=(c_double*6)()
                derivatives=(c_double*3)()
                probe.query_translation(index+1,request['temperature'],byref(enabled),byref(kind),parameters,derivatives,byref(dependent))
                values.extend([enabled.value,kind.value,*map(float,parameters),*map(float,derivatives),dependent.value])
            swappable=False
        elif request['mode']=='V_PROP':
            t,p,z=request['temperature'],request['pressure'],request['feed']
            phase=model.LIQPH if request['phase']=='L' else model.VAPPH
            v=model.specific_volume(t,p,z,phase)[0]
            pressure,pt,pv=model.pressure_tv(t,v,z,dpdt=True,dpdv=True)
            hr=model.enthalpy(t,p,z,phase,residual=True)[0]
            sr=model.entropy(t,p,z,phase,residual=True)[0]
            h,cp=model.enthalpy(t,p,z,phase,dhdt=True)
            s=model.entropy(t,p,z,phase)[0]
            phi,dt,dp,dn=model.thermo(t,p,z,phase,dlnfugdt=True,dlnfugdp=True,dlnfugdn=True)
            values=list(map(float,[v,pressure,pt,pv,hr,sr,h,s,cp,*phi,*dt,*dp,*dn.flatten()]))
            if not all(map(math.isfinite,values)):
                print('THERMOPACK_RESULT='+json.dumps({'calculation_error':'non-finite translated property',
                      'values':list(map(str,values))}),flush=True)
                continue
            swappable=False
        elif request['mode'] in ('ALPHA','ALPHA_PROPERTIES'):
            if probe is None:
                probe=load_probe()
            model.activate()
            probe.initialize_alpha(alpha.encode('ascii'),len(alpha),ref.encode('ascii'),len(ref))
            for index,name,parameters in request.get('setters',[]):
                model.set_alpha_corr(index+1,name,parameters)
                previous=None
            model.activate()
            values=[]
            for index in range(model.nc):
                method=c_int()
                parameters=(c_double*3)()
                derivatives=(c_double*3)()
                probe.query_alpha(index+1,request['temperature'],byref(method),parameters,derivatives)
                values.extend([method.value,*map(float,parameters),*map(float,derivatives)])
            if not all(map(math.isfinite,values)):
                print('THERMOPACK_RESULT='+json.dumps({'calculation_error':'non-finite alpha property',
                      'values':list(map(str,values))}),flush=True)
                continue
            if request['mode']=='ALPHA_PROPERTIES':
                t,p,z=request['temperature'],request['pressure'],request['feed']
                phase=model.LIQPH if request['phase']=='L' else model.VAPPH
                v=model.specific_volume(t,p,z,phase)[0]
                pressure,pt,pv=model.pressure_tv(t,v,z,dpdt=True,dpdv=True)
                hr=model.enthalpy(t,p,z,phase,residual=True)[0]
                sr=model.entropy(t,p,z,phase,residual=True)[0]
                h,cp=model.enthalpy(t,p,z,phase,dhdt=True)
                s=model.entropy(t,p,z,phase)[0]
                phi,dt,dp,dn=model.thermo(t,p,z,phase,dlnfugdt=True,dlnfugdp=True,dlnfugdn=True)
                values=list(map(float,[v,pressure,pt,pv,hr,sr,h,s,cp,*phi,*dt,*dp,*dn.flatten()]))
                if not all(map(math.isfinite,values)):
                    print('THERMOPACK_RESULT='+json.dumps({'calculation_error':'non-finite alpha EOS property',
                          'values':list(map(str,values))}),flush=True)
                    continue
            swappable=False
        elif request['mode'] in ('IDEAL','IDEAL_CONFIG'):
            if request['mode']=='IDEAL_CONFIG':
                # Every configuration case starts from fresh default offsets.
                model=cubic(comps,'PR',volume_shift=False)
                previous=None
                for action in request['actions']:
                    if action[0]=='cp':
                        model.set_ideal_cp(1,action[1],action[2])
                    elif action[0]=='h':
                        model.set_enthalpy_of_formation(1,action[1])
                    elif action[0]=='s':
                        model.set_standard_entropy(1,action[1],
                            reference_pressure='1ATM' if action[2]==101325. else '1BAR')
                    else:
                        raise ValueError('unknown ideal configuration action')
            t,p=request['temperature'],request['pressure']
            h,cp=model.idealenthalpysingle(t,1,dhdt=True)
            s,st,sp=model.idealentropysingle(t,p,1,dsdt=True,dsdp=True)
            values=list(map(float,[h,s,cp,st,sp]))
            if not all(map(math.isfinite,values)):
                print('THERMOPACK_RESULT='+json.dumps({'calculation_error':'non-finite ideal property',
                      'values':list(map(str,values))}),flush=True)
                continue
            if request['mode']=='IDEAL_CONFIG':
                corr,parameters=model.get_ideal_cp(1)
                values.extend([float(model.get_enthalpy_of_formation(1)),
                               float(model.get_standard_entropy(1)),int(corr),*map(float,parameters)])
            swappable=False
        elif request['mode']=='PSEUDO':
            model.activate()
            composition=(c_double*len(request['feed']))(*request['feed'])
            t,p,w,z,v=(c_double() for _ in range(5))
            routine=getattr(model.tp,model.get_export_name('eos','pseudo'))
            routine.argtypes=[POINTER(c_double)]*6
            routine.restype=None
            routine(composition,byref(t),byref(p),byref(w),byref(z),byref(v))
            values=[t.value,p.value,v.value,z.value]
            swappable=False
        elif request['mode'] in ('SAT_T','V_SAT'):
            z=request['feed']
            p=model.bubble_pressure(request['temperature'],z)[0]
            t=model.dew_temperature(p,z)[0]
            vl=model.specific_volume(t,p,z,model.LIQPH)[0]
            vv=model.specific_volume(t,p,z,model.VAPPH)[0]
            hl=model.enthalpy(t,p,z,model.LIQPH)[0]
            hv=model.enthalpy(t,p,z,model.VAPPH)[0]
            values=[float(t),float(p),float(vl),float(vv),float((hv-hl)/(t*(vv-vl)))]
            swappable=False
            extra={'pressure':float(p)}
        elif request['mode'] in ('TP','V_TP','V_TP_AUDIT'):
            result=model.two_phase_tpflash(request['temperature'],request['pressure'],request['feed'])
            values=[int(result.phase),float(result.betaV),float(result.betaL),
                    *map(float,result.x),*map(float,result.y)]
            swappable=phase_permutation_allowed(model,result)
            if request['mode']=='V_TP_AUDIT':
                t,p,z=request['temperature'],request['pressure'],request['feed']
                def gibbs(x,phase):
                    phi=model.thermo(t,p,x,phase,ophase=True)[0]
                    return sum(xi*(math.log(xi)+fi) for xi,fi in zip(x,phi) if xi>0)
                feed_g=gibbs(z,model.MINGIBBSPH)
                if result.phase==model.TWOPH:
                    split_g=result.betaL*gibbs(result.x,model.LIQPH)+result.betaV*gibbs(result.y,model.VAPPH)
                    l=model.thermo(t,p,result.x,model.LIQPH)[0]
                    v=model.thermo(t,p,result.y,model.VAPPH)[0]
                    residual=[math.log(yi/xi)+vi-li for xi,yi,vi,li,zi in zip(result.x,result.y,v,l,z) if zi>0 and xi>0 and yi>0]
                    balance=max(abs(result.betaL*xi+result.betaV*yi-zi) for xi,yi,zi in zip(result.x,result.y,z))
                    extra={'audit':{'feed_gibbs_rt':float(feed_g),'split_gibbs_rt':float(split_g),
                           'gibbs_increase_rt':float(split_g-feed_g),'chemical_potential_residual':residual,
                           'max_chemical_potential_residual':max(map(abs,residual)),'max_material_balance_residual':float(balance)}}
        elif request['mode'] in ('PH','PS','V_PH','V_PS'):
            mode=request['mode'][-2:]
            t,z=request['temperature'],request['feed']
            quantity='enthalpy' if mode=='PH' else 'entropy'
            if 'quality' in request:
                p=model.bubble_pressure(t,z)[0]
                calculate=getattr(model,quantity)
                value=(1-request['quality'])*calculate(t,p,z,model.LIQPH)[0]
                value+=request['quality']*calculate(t,p,z,model.VAPPH)[0]
            else:
                p=request['pressure']
                feed=model.two_phase_tpflash(t,p,z)
                value=overall(model,feed,quantity)
            try:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    calculate=model.two_phase_phflash if mode=='PH' else model.two_phase_psflash
                    result=calculate(p,z,value)
            except Exception as error:
                if str(error) not in ('PH flash calculation failed','PS flash calculation failed'):
                    raise
                print('THERMOPACK_RESULT='+json.dumps({'calculation_error':str(error),
                      'specification':float(value),'pressure':float(p)}),flush=True)
                continue
            if caught:
                if any(str(w.message)!='PH solver not fully converged' for w in caught):
                    raise RuntimeError('Unexpected oracle warning: '+str([str(w.message) for w in caught]))
                print('THERMOPACK_RESULT='+json.dumps({'calculation_error':[str(w.message) for w in caught],
                      'specification':float(value),'pressure':float(p)}),flush=True)
                continue
            actual=overall(model,result,quantity)
            scale=8.31446261815324*(result.T if mode=='PH' else 1.0)
            if not math.isfinite(actual) or abs(actual-value)>1e-8*scale:
                print('THERMOPACK_RESULT='+json.dumps({'calculation_error':'returned state does not close caloric specification',
                      'specification':float(value),'pressure':float(p),'returned_temperature':float(result.T),
                      'returned_property':float(actual),'residual':float(actual-value)}),flush=True)
                continue
            values=[float(result.T),int(result.phase),float(result.betaV),float(result.betaL),
                    *map(float,result.x),*map(float,result.y)]
            swappable=phase_permutation_allowed(model,result)
            extra={'specification':float(value),'pressure':float(p)}
        else:
            raise ValueError('unsupported oracle mode')
        print('THERMOPACK_RESULT='+json.dumps({'values':values,'phase_permutation_allowed':swappable,**extra}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--compiler-directory')
    main(parser.parse_args().compiler_directory)
