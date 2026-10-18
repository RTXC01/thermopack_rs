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
    return probe


def main(compiler_directory):
    cubic=load(compiler_directory)
    model,previous=None,None
    probe=None
    for line in sys.stdin:
        request=json.loads(line)
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
        if request['mode']=='V_RR':
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
