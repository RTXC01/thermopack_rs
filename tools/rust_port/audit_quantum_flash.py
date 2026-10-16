"""Record independently checked defects in returned exact-fork TP states."""
import argparse
import json
from reference import ROOT,ReferenceProcess


def audit(compiler_directory):
    cases=[('PR','H2,NE',[.6,.4],3.,1e3),('PR','HE,H2,NE',[0.,.5,.5],3.,1e3)]
    cases += [(eos,'C1,CO2,N2',[.5,.3,.2],80.,1e5) for eos in ['PR','SRK','RK']]
    results=[]
    with ReferenceProcess(compiler_directory) as oracle:
        for eos,comps,z,t,p in cases:
            request=dict(mode='V_TP_AUDIT',eos=eos,components=comps,feed=z,temperature=t,pressure=p,
                         parameter_reference='QuantumCubic',operations=[['bounds',.5,999.]])
            result=oracle.call(request)
            assert 'audit' in result,result
            check=result['audit']
            if t==3.:
                assert check['gibbs_increase_rt']>0.1 and check['max_chemical_potential_residual']<1e-8
            else:
                assert result['values'][1]<1e-12 and check['max_chemical_potential_residual']>1.
            results.append({'request':request,**result})
    output=dict(reference=json.loads((ROOT/'.reference/metadata.json').read_text()),results=results,
                generator='tools/rust_port/audit_quantum_flash.py',
                gibbs_convention='sum(x_i*(ln(x_i)+ln(phi_i))); linear ideal reference terms cancel at fixed feed T,P')
    (ROOT/'tests/fixtures/quantum_flash_audit.json').write_text(json.dumps(output,indent=2)+'\n')
    print(f'{len(results)} source-returned states audited: stationary higher-Gibbs splits and nonclosing vanishing phases')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory')
    audit(parser.parse_args().compiler_directory)
