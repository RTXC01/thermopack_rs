"""Line-oriented subprocess for a library that can terminate its caller."""
import argparse
import json
import sys
import warnings
import math
from ctypes import c_double,POINTER,byref
from reference import load


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


def main(compiler_directory):
    cubic=load(compiler_directory)
    model,previous=None,None
    for line in sys.stdin:
        request=json.loads(line)
        eos,comps=request['eos'],request['components']
        if previous!=(eos,comps):
            model=cubic(comps,"SRK" if eos=="RK" else eos,
                        alpha="RK" if eos=="RK" else "Classic",volume_shift=False)
            previous=eos,comps
        extra={}
        if request['mode'] in ('IDEAL','IDEAL_CONFIG'):
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
        elif request['mode']=='SAT_T':
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
        elif request['mode']=='TP':
            result=model.two_phase_tpflash(request['temperature'],request['pressure'],request['feed'])
            values=[int(result.phase),float(result.betaV),float(result.betaL),
                    *map(float,result.x),*map(float,result.y)]
            swappable=phase_permutation_allowed(model,result)
        elif request['mode'] in ('PH','PS'):
            mode=request['mode']
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
