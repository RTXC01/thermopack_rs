"""Line-oriented subprocess for a library that can terminate its caller."""
import argparse
import json
import sys
import warnings
import math
import os
from ctypes import c_double,c_int,c_char,c_char_p,CDLL,POINTER,byref
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
    probe.query_meos.argtypes=[c_double,c_double,c_double,c_int,POINTER(c_double),POINTER(c_double)]
    probe.query_meos.restype=None
    probe.query_meos_density.argtypes=[c_double,c_double,c_int,c_int,POINTER(c_double)]
    probe.query_meos_density.restype=None
    probe.query_fmt.argtypes=[c_int,c_int,POINTER(c_double),c_int,POINTER(c_double)]
    probe.query_fmt.restype=None
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
    probe.query_mbwr_parameters.argtypes=[c_int,POINTER(c_double)]
    probe.query_mbwr_parameters.restype=None
    probe.query_mbwr.argtypes=[c_int,c_double,c_double,c_double,POINTER(c_double)]
    probe.query_mbwr.restype=None
    probe.query_mbwr_integrated.argtypes=[c_int,POINTER(c_double)]
    probe.query_mbwr_integrated.restype=None
    probe.query_mbwr_properties.argtypes=[c_int,c_double,c_double,c_double,c_int,c_double,POINTER(c_double)]
    probe.query_mbwr_properties.restype=None
    probe.query_mbwr_metadata.argtypes=[c_int,POINTER(c_double)]
    probe.query_mbwr_metadata.restype=None
    probe.query_mbwr_critical_jet.argtypes=[c_int,c_double,c_double,POINTER(c_double)]
    probe.query_mbwr_critical_jet.restype=None
    probe.query_mbwr_density.argtypes=[c_int,c_double,c_double,c_int,POINTER(c_double),c_double,POINTER(c_double)]
    probe.query_mbwr_density.restype=None
    probe.query_mbwr_initialization.argtypes=[c_int,POINTER(c_double),POINTER(c_int)]
    probe.query_mbwr_initialization.restype=None
    probe.query_mbwr_critical_grid.argtypes=[c_int,POINTER(c_double),POINTER(c_double)]
    probe.query_mbwr_critical_grid.restype=None
    probe.query_unifac.argtypes=[c_int,POINTER(c_char),c_int,c_double,POINTER(c_double),POINTER(c_double),POINTER(c_int)]
    probe.query_unifac.restype=None
    probe.query_solid_line.argtypes=[c_int,POINTER(c_char),c_double,c_int,POINTER(c_double),POINTER(c_double),POINTER(c_int)]
    probe.query_solid_line.restype=None
    probe.query_solid_co2_init.argtypes=[c_double,c_double,c_double,c_double]
    probe.query_solid_co2_init.restype=None
    probe.query_solid_co2.argtypes=[c_double,c_double,POINTER(c_double),POINTER(c_int)]
    probe.query_solid_co2.restype=None
    probe.query_solid_h2o_init.argtypes=[c_double,c_double]
    probe.query_solid_h2o_init.restype=None
    probe.query_solid_h2o.argtypes=[c_double,c_double,POINTER(c_double)]
    probe.query_solid_h2o.restype=None
    probe.query_mut_solve.argtypes=[c_int,c_int,POINTER(c_char),POINTER(c_double),c_double,POINTER(c_double),POINTER(c_double),POINTER(c_int)]
    probe.query_mut_solve.restype=None
    probe.query_mut_target.argtypes=[c_int,c_double,POINTER(c_char),POINTER(c_double),c_int,POINTER(c_double),POINTER(c_int)]
    probe.query_mut_target.restype=None
    probe.query_hyperdual.argtypes=[POINTER(c_double),POINTER(c_double),POINTER(c_int)]
    probe.query_hyperdual.restype=None
    probe.query_solid_saturation.argtypes=[POINTER(c_char),c_double,c_double,c_int,c_int,c_double,POINTER(c_double),POINTER(c_int)]
    probe.query_solid_saturation.restype=None
    probe.query_cpa_rdf.argtypes=[c_double,c_double,POINTER(c_double),POINTER(c_double),POINTER(c_int)]
    probe.query_cpa_rdf.restype=None
    probe.query_wong_sandler.argtypes=[c_double,POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_int)]
    probe.query_wong_sandler.restype=None
    return probe


def main(compiler_directory):
    cubic=load(compiler_directory)
    model,previous=None,None
    probe=None
    for line in sys.stdin:
        request=json.loads(line)
        if request['mode'] in ['MEOS','MEOS_DENSITY']:
            if probe is None:probe=load_probe()
            reduced=int(request['eos']=='LJTS_REDUCED')
            key=('MEOS',request['eos'],request['component'])
            if not reduced:
                if previous!=key:
                    from thermopack.multiparameter import multiparam
                    model=multiparam(request['component'],request['eos']);previous=key
                model.activate()
            if request['mode']=='MEOS_DENSITY':
                output=(c_double*5)()
                probe.query_meos_density(request['temperature'],request['pressure'],request['phase'],reduced,output)
                print('THERMOPACK_RESULT='+json.dumps(dict(values=[float(v) if math.isfinite(v) else None for v in output])),flush=True)
                continue
            output=(c_double*46)();metadata=(c_double*18)()
            probe.query_meos(request['temperature'],request['volume'],request.get('amount',1.),reduced,output,metadata)
            values=[float(v) if math.isfinite(v) else None for v in output]
            print('THERMOPACK_RESULT='+json.dumps(dict(values=values,metadata=list(map(float,metadata)))),flush=True)
            continue
        if request['mode']=='FMT':
            if probe is None:probe=load_probe()
            grid=request['grid'];size=len(grid)
            inputs=(c_double*(6*size))(*(row[i] for i in range(6) for row in grid))
            output=(c_double*(43*size))()
            probe.query_fmt(request['model'],size,inputs,request['flags'],output)
            rows=[list(map(float,output[i*43:(i+1)*43])) for i in range(size)]
            # Non-finite entries are reported per row, never dropped from the matrix.
            print('THERMOPACK_RESULT='+json.dumps(dict(rows=[v if all(map(math.isfinite,v)) else list(map(str,v)) for v in rows])),flush=True)
            continue
        if request['mode']=='SOLID_H2O_INIT':
            if probe is None:probe=load_probe()
            probe.query_solid_h2o_init(request['entropy'],request['gibbs'])
            print('THERMOPACK_RESULT='+json.dumps(dict(initialized=True)),flush=True);continue
        if request['mode']=='SOLID_H2O':
            if probe is None:probe=load_probe()
            output=(c_double*12)();probe.query_solid_h2o(request['temperature'],request['pressure'],output)
            values=list(map(float,output))
            if all(map(math.isfinite,values)):print('THERMOPACK_RESULT='+json.dumps(dict(values=values)),flush=True)
            else:print('THERMOPACK_RESULT='+json.dumps(dict(calculation_error='non-finite source ice state',values=list(map(str,values)))),flush=True)
            continue
        if request['mode']=='SOLID_CO2_INIT':
            if probe is None:probe=load_probe()
            probe.query_solid_co2_init(request['liquid_entropy'],request['liquid_gibbs'],
                                       request['temperature'],request['pressure'])
            print('THERMOPACK_RESULT='+json.dumps(dict(values=[])),flush=True)
            continue
        if request['mode']=='SOLID_CO2':
            if probe is None:probe=load_probe()
            output=(c_double*13)();count=c_int()
            probe.query_solid_co2(request['temperature'],request['pressure'],output,byref(count))
            values=list(map(float,output))[:count.value]
            print('THERMOPACK_RESULT='+json.dumps(dict(values=values)),flush=True)
            continue
        if request['mode']=='HYPERDUAL':
            if probe is None:probe=load_probe()
            inputs=(c_double*9)(*request['inputs'])
            values=(c_double*96)();count=c_int()
            probe.query_hyperdual(inputs,values,byref(count))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:count.value])),flush=True)
            continue
        if request['mode']=='SOLID_SATURATION':
            if probe is None:probe=load_probe()
            component=request['component']
            from thermopack.cubic import cubic
            key=('SOLID_SATURATION',component,request.get('eos','PR'))
            if previous!=key:
                model=cubic(component,request.get('eos','PR'));previous=key
            model.activate()
            names=bytearray(b' '*20)
            encoded=component.encode()[:20]
            names[:len(encoded)]=encoded
            result=(c_double*2)();status=c_int()
            probe.query_solid_saturation((c_char*20).from_buffer_copy(bytes(names)),
                                        request['temperature'],request['pressure'],
                                        request['phase'],request['specification'],
                                        request['ln_spec'],result,byref(status))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,result)),status=status.value)),flush=True)
            continue
        if request['mode']=='CRITICAL':
            eos,comps=request['eos'],request['components']
            from thermopack.cubic import cubic
            model=cubic(",".join(comps),'SRK' if eos=='RK' else eos,volume_shift=request.get('volume_shift',False),
                        mixing=request.get('mixing','vdW'),alpha='RK' if eos=='RK' else 'Classic')
            previous=None
            model.activate()
            try:
                temperature,volume,pressure=model.critical(request['feed'],
                                        request.get('temperature',0.0),
                                        request.get('volume',0.0),
                                        request.get('tolerance',1.0e-7))
            except Exception as error:
                print('THERMOPACK_RESULT='+json.dumps(dict(calculation_error=str(error))),flush=True)
                continue
            values=[float(temperature),float(volume),float(pressure)]
            diagnostics={}
            if request.get('stability_diagnostics'):
                import numpy as np
                z=np.array(request['feed']);active=np.flatnonzero(z>0)
                weights=np.sqrt(z[active])
                hessian=model.fugacity_tv(temperature,volume,z,dlnphidn=True)[1]
                matrix=hessian[np.ix_(active,active)]*weights[:,None]*weights[None,:]
                eigenvalues,vectors=np.linalg.eigh(matrix)
                u=vectors[:,0]
                if u[np.argmax(abs(u))]<0:u=-u
                direction=np.zeros(len(z));direction[active]=weights*u
                step=min([1e-3]+[.1*z[i]/abs(direction[i]) for i in active if direction[i]!=0])
                def contraction(s):
                    h=model.fugacity_tv(temperature,volume,z+s*direction,dlnphidn=True)[1]
                    return sum(direction[i]*h[i,j]*direction[j] for i in active for j in active)
                third=(contraction(-2*step)-8*contraction(-step)+8*contraction(step)-contraction(2*step))/(12*step)
                diagnostics=dict(hessian=[float(hessian[i,j]) for i in active for j in active],
                                 eigenvalue=float(eigenvalues[0]),third_derivative=float(third))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=values,**diagnostics)),flush=True)
            continue
        if request['mode']=='QUADRATURE':
            if probe is None:probe=load_probe()
            probe.query_quadrature.argtypes=[c_int,POINTER(c_double),POINTER(c_int)]
            probe.query_quadrature.restype=None
            probe.query_quadrature_integral.argtypes=[c_int,c_double,c_double,POINTER(c_double),POINTER(c_int)]
            probe.query_quadrature_integral.restype=None
            quad=request['quadrature']
            values=(c_double*202)();number=c_int()
            probe.query_quadrature(quad,values,byref(number))
            rule=list(map(float,values))[:number.value]
            integrals=(c_double*2)();count=c_int()
            probe.query_quadrature_integral(quad,request['lower'],request['upper'],integrals,byref(count))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=rule,integral=float(integrals[0]),
                  differential=float(integrals[1]))),flush=True)
            continue
        if request['mode']=='AIS_LJS':
            if probe is None:probe=load_probe()
            probe.query_ais_ljs.argtypes=[c_int,c_double,c_double,POINTER(c_double),POINTER(c_int)]
            probe.query_ais_ljs.restype=None
            values=(c_double*2)();number=c_int()
            probe.query_ais_ljs(request['index'],request['reduced_distance'],
                               request['reduced_density'],values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='CPA_RDF':
            if probe is None:probe=load_probe()
            key=('CPA_RDF',request['component'])
            if previous!=key:
                from thermopack.cpa import cpa
                model=cpa(request['component']);previous=key
            model.activate()
            probe.query_cpa_rdf.argtypes=[c_double,c_double,POINTER(c_double),POINTER(c_double),POINTER(c_int)]
            probe.query_cpa_rdf.restype=None
            count=len(request['moles'])
            amounts=(c_double*16)(*request['moles'],*([0.0]*(16-count)))
            values=(c_double*256)();number=c_int()
            probe.query_cpa_rdf(request['temperature'],request['volume'],amounts,values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='PCSAFT_STATE':
            if probe is None:probe=load_probe()
            key=('PCSAFT_STATE',request['components'])
            if previous!=key:
                from thermopack.pcsaft import pcsaft
                model=pcsaft(request['components']);previous=key
            model.activate()
            probe.query_pcsaft_state.argtypes=[c_double,c_double,POINTER(c_double),
                                                POINTER(c_double),POINTER(c_int)]
            probe.query_pcsaft_state.restype=None
            count=len(request['moles'])
            amounts=(c_double*16)(*request['moles'],*([0.0]*(16-count)))
            values=(c_double*256)();number=c_int()
            probe.query_pcsaft_state(request['temperature'],request['volume'],amounts,
                                     values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='PCSAFT_CHAIN':
            if probe is None:probe=load_probe()
            key=('PCSAFT_CHAIN',request['components'])
            if previous!=key:
                from thermopack.pcsaft import pcsaft
                model=pcsaft(request['components']);previous=key
            model.activate()
            probe.query_pcsaft_chain.argtypes=[c_double,c_double,POINTER(c_double),
                                                POINTER(c_double),POINTER(c_int)]
            probe.query_pcsaft_chain.restype=None
            count=len(request['moles'])
            amounts=(c_double*16)(*request['moles'],*([0.0]*(16-count)))
            values=(c_double*256)();number=c_int()
            probe.query_pcsaft_chain(request['temperature'],request['volume'],amounts,
                                     values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='PCSAFT_HS':
            if probe is None:probe=load_probe()
            key=('PCSAFT_HS',request['components'])
            if previous!=key:
                from thermopack.pcsaft import pcsaft
                model=pcsaft(request['components']);previous=key
            model.activate()
            probe.query_pcsaft_hs.argtypes=[c_double,c_double,POINTER(c_double),
                                             POINTER(c_double),POINTER(c_int)]
            probe.query_pcsaft_hs.restype=None
            count=len(request['moles'])
            amounts=(c_double*16)(*request['moles'],*([0.0]*(16-count)))
            values=(c_double*256)();number=c_int()
            probe.query_pcsaft_hs(request['temperature'],request['volume'],amounts,
                                  values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='PCSAFT_RDF':
            if probe is None:probe=load_probe()
            key=('PCSAFT_RDF',request['components'])
            if previous!=key:
                from thermopack.pcsaft import pcsaft
                model=pcsaft(request['components']);previous=key
            model.activate()
            probe.query_pcsaft_rdf.argtypes=[c_double,c_double,POINTER(c_double),c_int,c_int,
                                             POINTER(c_double),POINTER(c_int)]
            probe.query_pcsaft_rdf.restype=None
            count=len(request['moles'])
            amounts=(c_double*16)(*request['moles'],*([0.0]*(16-count)))
            values=(c_double*256)();number=c_int()
            probe.query_pcsaft_rdf(request['temperature'],request['volume'],amounts,
                                   request['first'],request['second'],values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='WONGSANDLER':
            if probe is None:probe=load_probe()
            key=('WONGSANDLER',request['components'],request.get('eos','PR'))
            if previous!=key:
                from thermopack.cubic import cubic
                model=cubic(request['components'],request.get('eos','PR'),mixing='WongSandler')
                previous=key
            model.activate()
            probe.query_wong_sandler.argtypes=[c_double,POINTER(c_double),POINTER(c_double),
                                               POINTER(c_double),POINTER(c_double),POINTER(c_double),
                                               POINTER(c_int)]
            probe.query_wong_sandler.restype=None
            def _padded(values):
                data=list(values)+[0.0]*(16-len(values))
                return (c_double*16)(*data)
            composition=_padded(request['composition'])
            kij=_padded([request['kij']])
            tauij=_padded([request['tauij']])
            alpha=_padded([request['alpha']])
            values=(c_double*256)();number=c_int()
            probe.query_wong_sandler(request['temperature'],composition,kij,tauij,alpha,
                                     values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='WONGSANDLERSTATE':
            if probe is None:probe=load_probe()
            key=('WONGSANDLER',request['components'],request.get('eos','PR'))
            if previous!=key:
                from thermopack.cubic import cubic
                model=cubic(request['components'],request.get('eos','PR'),mixing='WongSandler')
                previous=key
            model.activate()
            probe.query_wong_sandler_state.argtypes=[c_double,c_double,POINTER(c_double),
                                                    POINTER(c_double),POINTER(c_double),
                                                    POINTER(c_double),POINTER(c_double),
                                                    POINTER(c_int)]
            probe.query_wong_sandler_state.restype=None
            def _padded_state(values):
                data=list(values)+[0.0]*(16-len(values))
                return (c_double*16)(*data)
            composition=_padded_state(request['composition'])
            kij=_padded_state([request['kij']])
            tauij=_padded_state([request['tauij']])
            alpha=_padded_state([request['alpha']])
            values=(c_double*256)();number=c_int()
            probe.query_wong_sandler_state(request['temperature'],request['volume'],
                                          composition,kij,tauij,alpha,values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='CHAIN_RULE':
            if probe is None:probe=load_probe()
            probe.query_chain_rule.argtypes=[POINTER(c_double),POINTER(c_int)]
            probe.query_chain_rule.restype=None
            values=(c_double*15)();number=c_int()
            probe.query_chain_rule(values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='UV_VIRIAL_DIFFERENCE':
            if probe is None:probe=load_probe()
            probe.query_uv_virial_difference.argtypes=[c_double,c_double,c_double,c_double,
                                                     POINTER(c_double),POINTER(c_int)]
            probe.query_uv_virial_difference.restype=None
            values=(c_double*12)();number=c_int()
            probe.query_uv_virial_difference(request['sigma'],request['eps_divk'],
                                            request['temperature'],request['rho'],
                                            values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='WCA_TERM':
            if probe is None:probe=load_probe()
            probe.query_wca_term.argtypes=[c_int,c_double,c_double,c_double,
                                           POINTER(c_double),POINTER(c_int)]
            probe.query_wca_term.restype=None
            values=(c_double*6)();number=c_int()
            probe.query_wca_term(request['index'],request['eps_divk'],
                                request['reduced_distance'],request['reduced_density'],
                                values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='UV_VIRIALS':
            if probe is None:probe=load_probe()
            probe.query_uv_virials.argtypes=[c_double,c_double,POINTER(c_double),POINTER(c_int)]
            probe.query_uv_virials.restype=None
            values=(c_double*9)();number=c_int()
            probe.query_uv_virials(request['eps_divk'],request['temperature'],values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='UV_TERMS':
            if probe is None:probe=load_probe()
            probe.query_uv_terms.argtypes=[c_double,c_double,c_double,c_double,POINTER(c_double),POINTER(c_int)]
            probe.query_uv_terms.restype=None
            values=(c_double*9)();number=c_int()
            probe.query_uv_terms(request['sigma'],request['eps_divk'],request['temperature'],
                                request['rho'],values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='UV_PHI':
            if probe is None:probe=load_probe()
            probe.query_uv_phi.argtypes=[c_double,c_double,c_double,POINTER(c_double),POINTER(c_int)]
            probe.query_uv_phi.restype=None
            values=(c_double*6)();number=c_int()
            probe.query_uv_phi(request['eps_divk'],request['temperature'],request['rho'],
                              values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='WHITE_BEAR':
            if probe is None:probe=load_probe()
            probe.query_white_bear.argtypes=[c_double,c_double,c_double,c_double,c_double,c_double,
                                            POINTER(c_double),POINTER(c_int)]
            probe.query_white_bear.restype=None
            values=(c_double*15)();number=c_int()
            probe.query_white_bear(request['temperature'],request['eps_divk'],request['eta'],
                                   request['r_star'],request['r_star_t'],request['r_star_tt'],
                                   values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='LJS_WCA_TERM':
            if probe is None:probe=load_probe()
            probe.query_ljs_wca_term.argtypes=[c_int,c_double,c_double,POINTER(c_double),POINTER(c_int)]
            probe.query_ljs_wca_term.restype=None
            values=(c_double*1)();number=c_int()
            probe.query_ljs_wca_term(request['index'],request['reduced_temperature'],
                                   request['reduced_density'],values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='LJS_POTENTIAL':
            if probe is None:probe=load_probe()
            probe.query_ljs_potential.argtypes=[c_int,POINTER(c_double),POINTER(c_double),POINTER(c_int)]
            probe.query_ljs_potential.restype=None
            radii=request['separations']
            count=len(radii)
            radius_buffer=(c_double*32)()
            for index,radius in enumerate(radii):
                radius_buffer[index]=radius
            values=(c_double*64)();number=c_int()
            probe.query_ljs_potential(count,radius_buffer,values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='LJS_CALORIC':
            from thermopack.ljs_bh import ljs_bh
            model=ljs_bh(request.get('parameter_reference','Default'),2.0)
            temperature=request['temperature']
            pressure=request['pressure']
            z=request['feed']
            phase=model.LIQPH if request['phase']==1 else model.VAPPH
            try:
                enthalpy=model.enthalpy(temperature,pressure,z,phase)[0]
                entropy=model.entropy(temperature,pressure,z,phase)[0]
            except Exception as error:
                print('THERMOPACK_RESULT='+json.dumps(dict(calculation_error=str(error))),flush=True)
                continue
            values=[float(enthalpy),float(entropy)]
            print('THERMOPACK_RESULT='+json.dumps(dict(values=values)),flush=True)
            continue
        if request['mode']=='LJS_FUGACITY':
            from thermopack.ljs_bh import ljs_bh
            model=ljs_bh(request.get('parameter_reference','Default'),2.0)
            temperature=request['temperature']
            pressure=request['pressure']
            z=request['feed']
            phase=model.LIQPH if request['phase']==1 else model.VAPPH
            try:
                logf=model.thermo(temperature,pressure,z,phase)[0]
            except Exception as error:
                print('THERMOPACK_RESULT='+json.dumps(dict(calculation_error=str(error))),flush=True)
                continue
            print('THERMOPACK_RESULT='+json.dumps(dict(values=[float(logf[0])])),flush=True)
            continue
        if request['mode']=='LJS_STATE':
            from thermopack.ljs_bh import ljs_bh
            model=ljs_bh(request.get('parameter_reference','Default'),2.0)
            temperature=request['temperature']
            pressure=request['pressure']
            z=request['feed']
            phase=model.LIQPH if request['phase']==1 else model.VAPPH
            try:
                volume=model.specific_volume(temperature,pressure,z,phase)[0]
                compressibility=model.zfac(temperature,pressure,z,phase)[0]
            except Exception as error:
                print('THERMOPACK_RESULT='+json.dumps(dict(calculation_error=str(error))),flush=True)
                continue
            values=[float(volume),float(compressibility)]
            print('THERMOPACK_RESULT='+json.dumps(dict(values=values)),flush=True)
            continue
        if request['mode']=='LJS_HELMHOLTZ':
            if probe is None:probe=load_probe()
            probe.query_ljs_helmholtz.argtypes=[c_double,c_double,c_double,POINTER(c_double),POINTER(c_int)]
            probe.query_ljs_helmholtz.restype=None
            from thermopack.ljs_bh import ljs_bh
            model=ljs_bh(request.get('parameter_reference','Default'),2.0)
            control=request.get('control')
            if control:
                model.model_control(**control)
            values=(c_double*3)();number=c_int()
            probe.query_ljs_helmholtz(request['temperature'],request['volume'],
                                     request['moles'],values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='LJS_SPLINE':
            if probe is None:probe=load_probe()
            probe.query_ljs_spline.argtypes=[c_double,c_double,c_int,POINTER(c_double),POINTER(c_int)]
            probe.query_ljs_spline.restype=None
            from thermopack.ljs_bh import ljs_bh
            model=ljs_bh(request.get('parameter_reference','Default'),2.0)
            values=(c_double*6)();number=c_int()
            probe.query_ljs_spline(request['reduced_temperature'],
                                   request['reduced_density'],
                                   request['index'],values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='LJS':
            from thermopack.ljs_bh import ljs_bh
            model=ljs_bh(request.get('parameter_reference','Default'),2.0)
            t_star=request['reduced_temperature']
            rho_star=request['reduced_density']
            first,second,third=model.get_pert_a(t_star,rho_star)
            diameter=model.get_bh_diameter_div_sigma(t_star)
            values=[float(first),float(second),float(third),float(diameter)]
            print('THERMOPACK_RESULT='+json.dumps(dict(values=values)),flush=True)
            continue
        if request['mode']=='THIRD_VIRIAL_MATRIX':
            if probe is None:probe=load_probe()
            probe.query_binary_third_virial_matrix.argtypes=[c_double,POINTER(c_double),POINTER(c_int)]
            probe.query_binary_third_virial_matrix.restype=None
            eos,comps=request['eos'],request['components']
            from thermopack.cubic import cubic
            model=cubic(",".join(comps),eos)
            model.activate()
            values=(c_double*4)();number=c_int()
            probe.query_binary_third_virial_matrix(request['temperature'],values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='VIRIAL_MATRIX':
            if probe is None:probe=load_probe()
            probe.query_second_virial_matrix.argtypes=[c_double,POINTER(c_double),POINTER(c_int)]
            probe.query_second_virial_matrix.restype=None
            eos,comps=request['eos'],request['components']
            from thermopack.cubic import cubic
            model=cubic(",".join(comps),eos)
            model.activate()
            values=(c_double*16)();number=c_int()
            probe.query_second_virial_matrix(request['temperature'],values,byref(number))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
            continue
        if request['mode']=='VIRIAL':
            eos,comps=request['eos'],request['components']
            from thermopack.cubic import cubic
            model=cubic(",".join(comps),eos)
            model.activate()
            second,third=model.virial_coeffcients(request['temperature'],request['feed'])
            values=[float(second),float(third)]
            print('THERMOPACK_RESULT='+json.dumps(dict(values=values)),flush=True)
            continue
        if request['mode']=='JTI':
            eos,comps=request['eos'],request['components']
            from thermopack.cubic import cubic
            model=cubic(",".join(comps),eos)
            model.activate()
            z=request['feed']
            try:
                temperature,pressure,volume=model.joule_thompson_inversion(z,request.get('nmax',100))
            except Exception as error:
                print('THERMOPACK_RESULT='+json.dumps(dict(calculation_error=str(error))),flush=True)
                continue
            arrays=[list(map(float,values)) for values in (temperature,pressure,volume)]
            print('THERMOPACK_RESULT='+json.dumps(dict(values=arrays)),flush=True)
            continue
        if request['mode']=='GE_DERIVATIVES':
            from thermopack.cubic import cubic
            model=cubic(','.join(request['components']),request['eos'],mixing=request['mixing'],
                parameter_reference=request.get('parameter_reference','DEFAULT'),volume_shift=request.get('volume_shift',False))
            previous=None;t=request['temperature'];p=request['pressure'];z=request['feed'];phase=model.LIQPH if request['phase']=='L' else model.VAPPH
            phi,dt,dp,dn=model.thermo(t,p,z,phase,dlnfugdt=True,dlnfugdp=True,dlnfugdn=True)
            v=float(model.specific_volume(t,p,z,phase)[0]);error=float(model.pressure_tv(t,v,z)[0])-p
            print('THERMOPACK_RESULT='+json.dumps(dict(values=[*map(float,dt),*map(float,dp),*map(float,dn.flatten())],volume=v,pressure_error=error)),flush=True)
            continue
        if request['mode']=='PURE_CURVE':
            from thermopack.cubic import cubic
            config=request.get('configuration','classic');eos=request['eos'];names=request['components']
            if config=='quantum':
                from thermopack.quantum_cubic import qcubic
                model=qcubic(','.join(names))
            else:
                model=cubic(','.join(names),'SRK' if eos=='RK' else eos,alpha='RK' if eos=='RK' else 'Classic',
                    mixing={'unifac':'UNIFAC','hv1':'HV1','hv2':'HV2','nrtl':'NRTL'}.get(config,'vdW'),
                    volume_shift=config=='translated',parameter_reference=request.get('parameter_reference','DEFAULT'))
            previous=None;component=request.get('component',0);z=[float(i==component) for i in range(len(names))]
            tc,vc,pc=map(float,model.get_critical_parameters(component+1))
            critical=[tc,pc,vc];value=request.get('initial_pressure')
            if value is not None and request.get('relative_start',False):value*=pc
            temp=request.get('initial_temperature')
            if temp is not None and request.get('relative_start',False):temp*=critical[0]
            t,p,vl,vv=model.get_pure_fluid_saturation_curve(value,initial_temperature=temp,i=component+1,
                max_delta_press=request.get('maximum_pressure_step',2e4),nmax=request.get('maximum_points',100),
                log_linear_grid=request.get('logarithmic',False))
            points=[]
            for ti,pi,li,vi in zip(t,p,vl,vv):
                phi_l=model.thermo(ti,pi,z,model.LIQPH)[0][component]
                phi_v=model.thermo(ti,pi,z,model.VAPPH)[0][component]
                pl=model.pressure_tv(ti,li,z)[0];pv=model.pressure_tv(ti,vi,z)[0]
                points.append(list(map(float,[ti,pi,li,vi,pl-pi,pv-pi,phi_l-phi_v])))
            # Critical cubic coefficients in the source's L/mol convention.
            a=float(model.get_energy_constants()[component]);b=float(model.get_covolumes()[component]);kr=1000.*8.31446261815324
            if probe is None:probe=load_probe()
            method=c_int();params_a=(c_double*3)();alpha_values=(c_double*3)()
            probe.query_alpha(component+1,tc,byref(method),params_a,alpha_values)
            params_b=(c_double*2)();beta_values=(c_double*3)();critical_data=(c_double*3)()
            probe.query_beta(component+1,tc,byref(method),params_b,beta_values,critical_data)
            a*=alpha_values[0];b*=beta_values[0]
            aa=a*pc/(tc*kr*tc*kr);bb=b*pc/(tc*kr)
            m1,m2=(-1.+math.sqrt(2.),-1.-math.sqrt(2.)) if eos=='PR' else (0.,-1.) if eos in ('SRK','RK') else (0.,0.)
            pp=bb*(-m1-m2-1.)-1.;qq=bb*((m1+m2)*(bb+1.)+bb*m1*m2)+aa;rr=-(bb**3+bb**2)*(m1*m2)-aa*bb
            shift=0.
            if config in ('translated','quantum'):
                if probe is None:probe=load_probe()
                enabled,kind,dependent=c_int(),c_int(),c_int();params=(c_double*6)();derivatives=(c_double*3)()
                probe.query_translation(component+1,tc,byref(enabled),byref(kind),params,derivatives,byref(dependent));shift=float(derivatives[0])
            print('THERMOPACK_RESULT='+json.dumps(dict(values=points,critical=critical,initial_pressure=value,initial_temperature=temp,
                critical_polynomial=[pp,qq,rr],critical_shift=shift)),flush=True)
            continue
        if request['mode']=='ENVELOPE':
            import numpy as np
            from thermopack.cubic import cubic
            eos=request['eos'];names=request['components'];z=np.array(request['feed']);nc=len(z);nmax=request.get('maximum_points',1000)
            if request.get('quantum'):
                from thermopack.quantum_cubic import qcubic
                model=qcubic(','.join(names))
            else:
                model=cubic(','.join(names),'SRK' if eos=='RK' else eos,alpha='RK' if eos=='RK' else 'Classic',
                    volume_shift=request.get('volume_shift',False),mixing=request.get('mixing','vdW'))
            previous=None;model.activate()
            amounts=(c_double*nc)(*z);tc=c_double(request.get('initial_temperature',0.));pc=c_double(request['initial_pressure'])
            spec=c_int(2 if 'initial_temperature' in request else 1);beta=c_double(request.get('initial_beta',1.))
            maximum=c_double(request.get('maximum_pressure',1.5e7));limit=c_int(nmax);count=c_int()
            temperatures=(c_double*nmax)();pressures=(c_double*nmax)();ratios=(c_double*(nmax*nc))();fractions=(c_double*nmax)()
            ds=c_double(request.get('maximum_step',.25));tmin=c_double(request.get('minimum_temperature',80.))
            calculate=model.s_envelope_plot
            calculate.argtypes=[POINTER(c_double)]*3+[POINTER(c_int)]+[POINTER(c_double)]*2+[POINTER(c_int)]+[POINTER(c_double)]*4+[POINTER(c_int)]+[POINTER(c_double)]*3+[POINTER(c_int)]+[POINTER(c_double)]*2
            calculate.restype=None
            calculate(amounts,byref(tc),byref(pc),byref(spec),byref(beta),byref(maximum),byref(limit),temperatures,pressures,ratios,fractions,byref(count),None,None,byref(ds),None,byref(tmin),None)
            points=[];bounds=[]
            fun=getattr(model.tp,model.get_export_name('saturation','sat_fun_newton'));fun.argtypes=[POINTER(c_double)]*3;fun.restype=None
            jac=getattr(model.tp,model.get_export_name('saturation','sat_diff_newton'));jac.argtypes=[POINTER(c_double)]*3;jac.restype=None
            for i in range(count.value):
                t=temperatures[i];p=pressures[i];b=fractions[i];k=np.array([ratios[j*nmax+i] for j in range(nc)])
                x=z/(1-b+b*k);y=k*x;x/=x.sum();y/=y.sum()
                vl=float(model.specific_volume(t,p,x,model.LIQPH)[0]);vv=float(model.specific_volume(t,p,y,model.VAPPH)[0])
                pl=float(model.pressure_tv(t,vl,x)[0]);pv=float(model.pressure_tv(t,vv,y)[0])
                # Independent source Newton backward error at a nontrivial
                # locked phase-ratio coordinate, including near critical states.
                u=np.r_[np.log(k),math.log(t),math.log(p)];active=np.flatnonzero(z>0);locked=int(active[np.argmax(np.abs(u[active]))])
                param=np.r_[z,b,float(locked+1),u[locked],2.]
                uc=(c_double*(nc+2))(*u);params=(c_double*(nc+4))(*param);f=(c_double*(nc+2))();j=(c_double*((nc+2)**2))()
                fun(f,uc,params);jac(j,uc,params)
                matrix=np.array(j).reshape((nc+2,nc+2),order='F');residual=np.array(f)
                try:
                    correction=np.abs(np.linalg.solve(matrix,residual));condition=float(np.linalg.cond(matrix))
                    rhs=np.zeros(nc+2);rhs[-1]=1.;tangent=np.linalg.solve(matrix,rhs);best=int(np.argmax(np.abs(tangent)))
                except np.linalg.LinAlgError:correction=np.full(nc+2,math.inf);condition=math.inf;best=locked
                points.append([t,p,b,vl,vv,*x,*y,pl-p,pv-p,float(np.max(np.abs(residual)))])
                bounds.append(dict(correction=correction.tolist(),condition=condition,locked=locked,log_k=u[:nc].tolist(),best_coordinate=best))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=points,bounds=bounds,count=count.value)),flush=True)
            continue
        if request['mode'] in ('META_PS','META_CURVE'):
            from thermopack.cubic import cubic
            eos=request['eos'];comps=request['components'];z=request['feed']
            config=request.get('configuration','classic')
            if config=='quantum':
                from thermopack.quantum_cubic import qcubic
                model=qcubic(','.join(comps))
            else:
                model=cubic(','.join(comps),'SRK' if eos=='RK' else eos,
                    alpha='RK' if eos=='RK' else 'Classic',volume_shift=config=='translated',
                    mixing={'unifac':'UNIFAC','hv1':'HV1'}.get(config,'vdW'))
            previous=None;model.activate()
            t=request['temperature'];v=request['volume']
            if request['mode']=='META_CURVE' and 'initial_pressure' in request:
                v=float(model.specific_volume(t,request['initial_pressure'],z,model.LIQPH)[0])
            initial_volume=float(v)
            entropy=float(model.entropy_tv(t,v,z)[0]);pressure=float(model.pressure_tv(t,v,z)[0])
            if request['mode']=='META_CURVE':pressure=request.get('initial_pressure',pressure)
            print('THERMOPACK_META_SPEC='+json.dumps(dict(pressure=pressure,entropy=entropy)),flush=True)
            if request['mode']=='META_PS':
                tg=request.get('temperature_guess',t*1.02);vg=request.get('volume_guess',v*.98)
                pc=c_double(pressure);sc=c_double(entropy);amounts=(c_double*len(z))(*z)
                tc=c_double(tg);vc=c_double(vg);status=c_int()
                calculate=getattr(model.tp,model.get_export_name('spinodal','tv_meta_ps'))
                calculate.argtypes=[POINTER(c_double)]*5+[POINTER(c_int)];calculate.restype=None
                calculate(byref(pc),byref(sc),amounts,byref(tc),byref(vc),byref(status))
                result=dict(values=[tc.value,vc.value],pressure=pressure,entropy=entropy,status=status.value)
                if math.isfinite(tc.value) and math.isfinite(vc.value) and tc.value>0 and vc.value>0:
                    result['pressure_error']=float(model.pressure_tv(tc.value,vc.value,z)[0])-pressure
                    result['entropy_error']=float(model.entropy_tv(tc.value,vc.value,z)[0])-entropy
            else:
                try:
                    temperatures,volumes,pressures=model.map_meta_isentrope(z,pressure,entropy,request['minimum_pressure'],request['points'])
                    values=[]
                    for t,v,p in zip(temperatures,volumes,pressures):
                        pe=float(model.pressure_tv(t,v,z)[0])-float(p)
                        se=float(model.entropy_tv(t,v,z)[0])-entropy
                        values.append(list(map(float,[t,v,p,pe,se])))
                    result=dict(values=values,pressure=pressure,entropy=entropy,initial_volume=initial_volume)
                except Exception as error:result=dict(calculation_error=str(error),pressure=pressure,entropy=entropy)
            print('THERMOPACK_RESULT='+json.dumps(result),flush=True)
            continue
        if request['mode'] in ('PURE_SPINODAL_GRID','PURE_SPINODAL_PROP','PURE_SPINODAL_EXTREMA'):
            from thermopack.cubic import cubic
            import numpy as np
            eos=request['eos'];names=request['components'];z=request['feed'];n=len(z)
            model=cubic(','.join(names),'SRK' if eos=='RK' else eos,
                        alpha='RK' if eos=='RK' else 'Classic',volume_shift=request.get('volume_shift',False))
            previous=None;model.activate()
            p0=c_double(request['pressure']);amounts=(c_double*n)(*z);status=c_int()
            if request['mode']=='PURE_SPINODAL_EXTREMA':
                prop=c_int(request['property']);coordinates=[c_double() for _ in range(6)]
                calculate=getattr(model.tp,model.get_export_name('spinodal','locate_spinodal_prop_min_max_pure_fluid'))
                calculate.argtypes=[POINTER(c_int),*([POINTER(c_double)]*8),POINTER(c_int)];calculate.restype=None
                calculate(byref(prop),amounts,byref(p0),*[byref(v) for v in coordinates],byref(status))
                result=dict(values=[v.value for v in coordinates],status=status.value)
                if not status.value and all(math.isfinite(v.value) for v in coordinates):
                    diagnostics=[]
                    for offset in [0,3]:
                        t,v,p=[c.value for c in coordinates[offset:offset+3]]
                        actual=float(model.pressure_tv(t,v,z)[0])
                        h=model.fugacity_tv(t,v,z,dlnphidn=True)[1]
                        value={2:lambda:v,4:lambda:t}.get(prop.value)
                        if value is None:
                            method={1:'entropy_tv',3:'enthalpy_tv',5:'pressure_tv',6:'internal_energy_tv'}[prop.value]
                            value=lambda:float(getattr(model,method)(t,v,z)[0])
                        diagnostics.append([actual,float(h[0,0]),value()])
                    result['diagnostics']=diagnostics
                    if prop.value==5:
                        a=[v.value for v in coordinates]
                        result['alternative_pressures']=[float(model.pressure_tv(a[0],a[4],z)[0]),float(model.pressure_tv(a[0],a[1],z)[0])]
            elif request['mode']=='PURE_SPINODAL_GRID':
                vt=request.get('vapor_temperature',.6*model.critical_temperature(1))
                ts,vs,ps=model.spinodal(z,initial_pressure=request['pressure'],
                    initial_liquid_temperature=request.get('liquid_temperature'),min_temperature_vapor=vt)
                data=[]
                for t,v,p in zip(ts,vs,ps):
                    h=model.fugacity_tv(t,v,z,dlnphidn=True)[1]
                    active=np.flatnonzero(np.array(z)>0);w=np.sqrt(np.array(z)[active])
                    lam=float(np.linalg.eigvalsh(h[np.ix_(active,active)]*w[:,None]*w[None,:])[0])
                    actual=float(model.pressure_tv(t,v,z)[0]);data.append([t,v,p,actual,lam])
                result=dict(values=data,status=status.value,vapor_temperature=vt)
            else:
                prop=c_int(request['property']);phase=c_int(request['phase'])
                def property_value(t,v):
                    if prop.value==2:return v
                    if prop.value==4:return t
                    method={1:'entropy_tv',3:'enthalpy_tv',5:'pressure_tv',6:'internal_energy_tv'}[prop.value]
                    return float(getattr(model,method)(t,v,z)[0])
                target=c_double(property_value(*request['target_point']) if 'target_point' in request else request['target'])
                # Preserve the exact specification even when the Fortran call
                # terminates before returning its structured result.
                print('THERMOPACK_TARGET='+repr(target.value),flush=True)
                t=c_double();v=c_double();p=c_double();found=c_int()
                calculate=getattr(model.tp,model.get_export_name('spinodal','locate_spinodal_prop_pure_fluid'))
                calculate.argtypes=[POINTER(c_int),POINTER(c_double),POINTER(c_double),POINTER(c_int),*([POINTER(c_double)]*4),POINTER(c_int),POINTER(c_int),POINTER(c_double)]
                calculate.restype=None
                calculate(byref(prop),amounts,byref(target),byref(phase),byref(p0),byref(t),byref(v),byref(p),byref(found),byref(status),None)
                result=dict(values=[t.value,v.value,p.value],found=bool(found.value),status=status.value,target=target.value)
                if found.value:
                    result['eigenvalue']=float(model.fugacity_tv(t.value,v.value,z,dlnphidn=True)[1][0,0])
                    values=[float(model.entropy_tv(t.value,v.value,z)[0]),v.value,
                            float(model.enthalpy_tv(t.value,v.value,z)[0]),t.value,
                            float(model.pressure_tv(t.value,v.value,z)[0]),float(model.internal_energy_tv(t.value,v.value,z)[0])]
                    result['property_value']=values[prop.value-1]
            print('THERMOPACK_RESULT='+json.dumps(result),flush=True)
            continue
        if request['mode'] in ('SPINODAL','SPINODAL_CURVE'):
            eos,comps=request['eos'],request['components']
            from thermopack.cubic import cubic
            model=cubic(",".join(comps),'SRK' if eos=='RK' else eos,
                        alpha='RK' if eos=='RK' else 'Classic',volume_shift=request.get('volume_shift',False))
            previous=None
            model.activate()
            z=request['feed']
            try:
                if request['mode']=='SPINODAL_CURVE':
                    temperature,volume,pressure=model.spinodal(z,initial_pressure=request['pressure'],
                        dlnv=request.get('log_volume_step'),initial_liquid_temperature=request.get('temperature'),
                        min_temperature_vapor=request.get('minimum_temperature'))
                else:
                    phase=model.LIQPH if request['phase']==1 else model.VAPPH
                    t,v=model.spinodal_point(z,request['pressure'],phase,temperature=request.get('temperature'))
                    temperature=[t];volume=[v];pressure=[request['pressure']]
                if request.get('stability_diagnostics'):
                    import numpy as np
                    z=np.array(z);active=np.flatnonzero(z>0);weights=np.sqrt(z[active])
                    data=[]
                    for t,v,p in zip(temperature,volume,pressure):
                        h=model.fugacity_tv(t,v,z,dlnphidn=True)[1]
                        matrix=h[np.ix_(active,active)]*weights[:,None]*weights[None,:]
                        lam=float(np.linalg.eigvalsh(matrix)[0])
                        actual=float(model.pressure_tv(t,v,z)[0])
                        data.append(list(map(float,[t,v,p,actual,lam])))
                    result=dict(values=data)
                else:result=dict(values=list(map(float,[temperature[0],volume[0],pressure[0]])))
            except Exception as error:result=dict(calculation_error=str(error))
            print('THERMOPACK_RESULT='+json.dumps(result),flush=True)
            continue
        if request['mode']=='ISOLINE':
            eos,comps=request['eos'],request['components']
            from thermopack.cubic import cubic
            model=cubic(",".join(comps),"SRK" if eos=="RK" else eos)
            model.activate()
            z=request['feed']
            points=request.get('points',20)
            if request['line']=='isobar':
                t,v,s,h=model.get_isobar(request['pressure'],z,
                                        request['minimum'],request['maximum'],points)
                p=[request['pressure']]*len(t)
            elif request['line']=='isenthalp':
                t,p,v,s=model.get_isenthalp(request['enthalpy'],z,
                                            request['minimum'],request['maximum'],
                                            request['tmin'],request['tmax'],points)
                h=[request['enthalpy']]*len(t)
            elif request['line']=='isentrope':
                t,p,v,h=model.get_isentrope(request['entropy'],z,
                                            request['minimum'],request['maximum'],
                                            request['tmin'],request['tmax'],points)
                s=[request['entropy']]*len(t)
            else:
                p,v,s,h=model.get_isotherm(request['temperature'],z,
                                         request['minimum'],request['maximum'],points)
                t=[request['temperature']]*len(p)
            arrays=[list(map(float,values)) for values in (t,p,v,s,h)]
            # The source writes a not-a-number point into a grid of one point,
            # so the isolines may carry that value.
            print('THERMOPACK_RESULT='+json.dumps(dict(values=arrays)),flush=True)
            continue
        if request['mode']=='PSEUDO_COMPONENT':
            from thermopack.cubic import cubic
            comps=request['components']
            model=cubic(",".join(comps),request.get('eos','PR'))
            model.init_pseudo(",".join(comps),request['critical_temperatures'],
                              request['critical_pressures'],request['acentric_factors'],
                              request.get('molar_masses'))
            model.activate()
            t,p,z=request['temperature'],request['pressure'],request['feed']
            phase=model.LIQPH if request['phase']==1 else model.VAPPH
            v=model.specific_volume(t,p,z,phase)[0]
            h=model.enthalpy(t,p,z,phase)[0]
            s=model.entropy(t,p,z,phase)[0]
            compressibility=model.zfac(t,p,z,phase)[0]
            values=[float(v),float(h),float(s),float(compressibility)]
            if not all(map(math.isfinite,values)):
                print('THERMOPACK_RESULT='+json.dumps(dict(calculation_error='non-finite pseudo property',
                      values=list(map(str,values)))),flush=True)
                continue
            print('THERMOPACK_RESULT='+json.dumps(dict(values=values)),flush=True)
            continue
        if request['mode']=='MU_T':
            if probe is None:probe=load_probe()
            components=request['components']
            count=len(components)
            from thermopack.cubic import cubic
            key=('MU_T',tuple(components),request.get('eos','PR'))
            if previous!=key:
                model=cubic(",".join(components),request.get('eos','PR'));previous=key
            model.activate()
            names=bytearray(b' '*20*count)
            for index,identifier in enumerate(components):
                names[index*20:index*20+20]=identifier.encode()[:20].ljust(20,b' ')
            densities=(c_double*count)(*request['densities'])
            if request['kind']=='TARGET':
                values=(c_double*count)();number=c_int()
                probe.query_mut_target(count,request['temperature'],
                                       (c_char*len(names)).from_buffer_copy(bytes(names)),
                                       densities,int(request.get('target_kind',0)),values,byref(number))
                print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values))[:number.value])),flush=True)
                continue
            target=(c_double*count)(*request['target'])
            values=(c_double*count)();status=c_int()
            probe.query_mut_solve(count,0 if request['kind']=='MU' else 1,
                                  (c_char*len(names)).from_buffer_copy(bytes(names)),
                                  densities,request['temperature'],target,values,byref(status))
            print('THERMOPACK_RESULT='+json.dumps(dict(values=list(map(float,values)),status=status.value)),flush=True)
            continue
        if request['mode']=='SOLID_LINE':
            if probe is None:probe=load_probe()
            identifier=request['component']
            names=bytearray(b' '*20)
            encoded=identifier.encode()[:20]
            names[:len(encoded)]=encoded
            points=int(request['points'])
            temperatures=(c_double*points)();pressures=(c_double*points)();count=c_int()
            probe.query_solid_line(int(request['kind']),(c_char*20).from_buffer_copy(bytes(names)),
                                  request['temperature'],points,temperatures,pressures,byref(count))
            values=list(map(float,temperatures))[:count.value]+list(map(float,pressures))[:count.value]
            print('THERMOPACK_RESULT='+json.dumps(dict(values=values)),flush=True)
            continue
        if request['mode']=='UNIFAC':
            if probe is None:probe=load_probe()
            components=request['components']
            names=bytearray(b' '*512)
            for index,identifier in enumerate(components):
                encoded=identifier.encode()[:20].ljust(20,b' ')
                names[index*20:index*20+20]=encoded
            composition=(c_double*len(components))(*request['composition'])
            output=(c_double*128)();count=c_int()
            probe.query_unifac(len(components),(c_char*512).from_buffer_copy(bytes(names)),
                               int(request.get('variant',0)),request['temperature'],
                               composition,output,byref(count))
            values=list(map(float,output))[:count.value]
            print('THERMOPACK_RESULT='+json.dumps(dict(values=values)),flush=True)
            continue
        if request['mode'].startswith('MBWR_'):
            if probe is None:probe=load_probe()
            record=request['record_index']
            if request['mode']=='MBWR_INITIALIZATION':
                output=(c_double*8)();status=c_int();probe.query_mbwr_initialization(record,output,byref(status))
                if status.value:
                    print('THERMOPACK_RESULT='+json.dumps(dict(undefined_source_initialization='getCompDBindex returned 0; original initializer reads compdb(0)',probe_status=status.value)),flush=True)
                    continue
                values=list(map(float,output))
            elif request['mode']=='MBWR_CRITICAL_GRID':
                initial=(c_double*3)(*request['initial'])
                output=(c_double*3)();probe.query_mbwr_critical_grid(record,initial,output)
                values=list(map(float,output))
            elif request['mode']=='MBWR_METADATA':
                output=(c_double*8)();probe.query_mbwr_metadata(record,output)
                values=list(map(float,output))
            elif request['mode']=='MBWR_CRITICAL_JET':
                output=(c_double*5)();probe.query_mbwr_critical_jet(record,request['temperature'],request['density'],output)
                values=list(map(float,output))
            elif request['mode']=='MBWR_DENSITY':
                critical=(c_double*3)(*request['critical'])
                output=(c_double*4)();probe.query_mbwr_density(record,request['temperature'],request['pressure'],
                    request['phase'],critical,request['acentric_factor'],output)
                values=list(map(float,output))
            elif request['mode']=='MBWR_TP_PROPERTIES':
                output=(c_double*18)();probe.query_mbwr_properties(record,request['temperature'],request['volume'],request['moles'],
                    int(request.get('off_eos',False)),request.get('pressure_factor',1.),output)
                values=list(map(float,output))
            elif request['mode']=='MBWR_INTEGRATED':
                output=(c_double*40)();probe.query_mbwr_integrated(record,output)
                values=list(map(float,output))[:40 if record<=11 else 22]
            elif request['mode']=='MBWR_PARAMETERS':
                output=(c_double*34)();probe.query_mbwr_parameters(record,output)
                values=list(map(float,output))[:int(output[0])+2]
            else:
                output=(c_double*26)();probe.query_mbwr(record,request['temperature'],request['volume'],request['moles'],output)
                values=list(map(float,output))
            if all(map(math.isfinite,values)):
                print('THERMOPACK_RESULT='+json.dumps(dict(values=values)),flush=True)
            else:print('THERMOPACK_RESULT='+json.dumps(dict(calculation_error='non-finite original MBWR kernel output',values=list(map(str,values)))),flush=True)
            continue
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
        mixing=request.get('mixing')
        key=(eos,comps,alpha,ref,shift,quantum_helper,mixing,'alpha_kernel' if request['mode'] in ('ALPHA','ALPHA_PROPERTIES') else 'public')
        if previous!=key:
            # ALPHA qualifies tpInitAlphaCorr separately from eoslibinit's
            # forced tcPR translation / QuantumCubic critical-data overrides.
            alpha_mode=request['mode'] in ('ALPHA','ALPHA_PROPERTIES')
            base_ref=('tcPR' if 'BUTANAL' in comps.split(',') else 'Default') if alpha_mode else ref
            mixing_option={'mixing':mixing} if mixing else {}
            if quantum_helper:
                from thermopack.quantum_cubic import qcubic
                model=qcubic(comps)
            else:
                model=cubic(comps,"SRK" if eos=="RK" else eos,
                            alpha='Classic' if alpha_mode else alpha,parameter_reference=base_ref,volume_shift=shift,
                            **mixing_option)
            previous=key
        extra={}
        if request['mode']=='MIXSAT':
            import numpy as np
            method=request['calculation'];z=request['feed'];known=request['specification']
            try:
                if method=='BP':p,y=model.bubble_pressure(known,z);t=known;x=np.array(z)
                elif method=='DP':p,x=model.dew_pressure(known,z);t=known;y=np.array(z)
                elif method=='BT':t,y=model.bubble_temperature(known,z);p=known;x=np.array(z)
                elif method=='DT':t,x=model.dew_temperature(known,z);p=known;y=np.array(z)
                else:raise ValueError('Unknown saturation method')
                vl=float(model.specific_volume(t,p,x,model.LIQPH)[0]);vv=float(model.specific_volume(t,p,y,model.VAPPH)[0])
                phi_l=model.thermo(t,p,x,model.LIQPH)[0];phi_v=model.thermo(t,p,y,model.VAPPH)[0]
                active=[i for i,v in enumerate(z) if v>0]
                if any(x[i]<=0 or y[i]<=0 for i in active):residual=float('inf')
                else:residual=max(abs(math.log(x[i])+phi_l[i]-math.log(y[i])-phi_v[i]) for i in active)
                pl=float(model.pressure_tv(t,vl,x)[0]);pv=float(model.pressure_tv(t,vv,y)[0])
                values=[float(t),float(p),vl,vv,*map(float,x),*map(float,y)]
                print('THERMOPACK_RESULT='+json.dumps(dict(values=values,fugacity_residual=float(residual),pressure_error=max(abs(pl-p),abs(pv-p))/max(1.,abs(p)),normalization_error=max(abs(sum(x)-1.),abs(sum(y)-1.)),phase_distance=(vv-vl)/max(vv,vl))),flush=True)
            except Exception as exc:
                print('THERMOPACK_RESULT='+json.dumps(dict(calculation_error=str(exc))),flush=True)
            continue
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
