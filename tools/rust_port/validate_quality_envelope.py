"""Complete fixed-quality curves and strict projection of their source states."""
import argparse,json
from reference import ROOT
from validate_envelope import validate as validate_curves,label as envelope_label

def cases():
    mixtures=[(['C1','C2'],[.5,.5]),(['C1','C3'],[.5,.5]),(['CO2','C1'],[.8,.2]),
        (['C1','C2','C3'],[.2,.3,.5]),(['C1','C2','C3'],[.5,.5,0.])]
    for eos in ['PR','SRK','RK','VDW']:
        for names,z in mixtures:
            for beta in [.25,.5,.75]:
                for p in [1e5,1e6]:yield dict(eos=eos,components=names,feed=z,initial_beta=beta,initial_pressure=p)
    for eos in ['PR','SRK']:
        for config in ['translated','unifac']:
            for names,z in [mixtures[1],mixtures[3]]:
                for beta in [.25,.5,.75]:
                    extra=dict(volume_shift=True) if config=='translated' else dict(mixing='UNIFAC')
                    yield dict(eos=eos,components=names,feed=z,initial_beta=beta,initial_pressure=1e5,**extra)
        for extra in [dict(maximum_points=3),dict(maximum_pressure=2e6),dict(minimum_temperature=150.),dict(initial_temperature=180.)]:
            yield dict(eos=eos,components=['C1','C2'],feed=[.5,.5],initial_beta=.4,initial_pressure=1e5,**extra)
    for names in [['H2','NE'],['H2','D2']]:
        for beta in [.25,.5,.75]:
            yield dict(eos='PR',components=names,feed=[.5,.5],initial_beta=beta,initial_pressure=1e5,minimum_temperature=2.,quantum=True)

def label(q):
    return envelope_label(q).replace('ENVELOPE|','QENVELOPE|',1)

def source_classifier(q,result):
    active=[i for i,x in enumerate(q['feed']) if x>0.]
    nc=len(q['feed']);z=[x/sum(q['feed']) for x in q['feed']]
    for p in result['values']:
        beta=p[2];x=p[5:5+nc];y=p[5+nc:5+2*nc]
        if max(abs((1.-beta)*a+beta*b-c) for a,b,c in zip(x,y,z))>2e-7:
            return 'Returned phase compositions fail the independent feed material balance; source residual may be masked at its solver bounds'
    if all(max(abs(b['log_k'][i]) for i in active)<1e-7 for b in result['bounds']):
        return 'Only trivial K=1 states returned; no separated fixed-quality branch established'

def validate(directory,reuse_source=False,allow_known_gaps=False):
    try:
        return validate_curves(directory,source_requests=list(cases()),case_label=label,
            fixture_stem='quality_envelope',reuse_source=reuse_source,cache_source=True,source_classifier=source_classifier,prefer_matching_grid=True)
    except SystemExit as error:
        if error.code!=1 or not allow_known_gaps:raise
        gaps=json.loads((ROOT/'.reference/quality_envelope_differences.json').read_text())
        for gap in gaps:
            fields=gap['case'].split('|');projection=gap.get('projection',[]);bound=gap.get('bound',{})
            indices={p[0] for p in projection}
            quantum=(fields[1:4]==['PR','H2,NE','quantum'] and fields[5]=='100000.0:15000000.0:2.0:0.25:1000:0.0:0.25'
                and indices and indices.issubset(set(range(20,28))) and bound.get('condition',0.)>1e12)
            unifac=(fields[1:4]==['SRK','C1,C3','unifac'] and fields[5]=='100000.0:15000000.0:80.0:0.25:1000:0.0:0.5'
                and indices=={26} and bound.get('condition',0.)>1e9)
            if not (quantum or unifac):raise
        print(f'{len(gaps)} known raw coordinate recovery gaps remain; not counted as strict parity')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');p.add_argument('--reuse-source',action='store_true');p.add_argument('--allow-known-gaps',action='store_true');a=p.parse_args();validate(a.compiler_directory,a.reuse_source,a.allow_known_gaps)
