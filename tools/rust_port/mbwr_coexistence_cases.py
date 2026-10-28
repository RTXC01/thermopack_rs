"""Fixed all-fit coexistence cohorts, independent of native solver outcomes."""
import csv
import hashlib
from reference import ROOT
from generate_data import compiled_components
from generate_mbwr_data import records
from validate_mbwr_density import cases as density_cases
from validate_meos_flash import digest


def fixture(name):
    with (ROOT/'tests/fixtures'/name).open(encoding='utf-8',newline='') as handle:
        return list(csv.DictReader(handle,delimiter=';'))


def models():
    general, cp = compiled_components(ROOT/'src/compdatadb.f90')
    available = {name.upper() for name in general if name in cp}
    critical = {row['case'].split('|')[1]:[float(row[k]) for k in ('Tc','Pc','rhoc')]
                for row in fixture('mbwr_critical.csv')}
    initialized = {row['case'].split('|')[1]:[float(row[k]) for k in ('Tc','Pc','rhoc','acentric_factor')]
                   for row in fixture('mbwr_initialization.csv')}
    status = {row['case'].split('|')[1]:int(row['source_status'])
              for row in fixture('mbwr_initialization_status.csv')}
    supplied = {}
    for q in density_cases():
        supplied.setdefault(q['source_name'],[*q['critical'],q['acentric_factor']])
    result=[];seen=set()
    for r in records():
        name=r['source_name'];supported=r['identifier'].upper() in available
        identity=r['kind'],r['identifier'].upper()
        canonical=supported and identity not in seen
        if supported:seen.add(identity)
        parameters=initialized.get(name,supplied[name])
        active=critical.get(name,parameters[:3])
        result.append(dict(record=name,index=int(name[6:]),component=r['identifier'] if supported else 'AR',
            parameters=parameters,critical=active,refine=int(name in critical),public=canonical,
            source_initialization_status=status[name],general_metadata=supported,
            density_metadata='original_initialization' if name in initialized else 'explicit_supplied',
            critical_metadata='independent_original_EOS' if name in critical else 'explicit_supplied'))
    assert len(result)==65 and sum(m['public'] for m in result)==46
    return result


def reference_requests():
    out=[]
    for model in models():
        for context in ('public','internal'):
            if context=='public' and not model['public']:continue
            tc,pc,_=model['critical']
            for tr in (.7,.9,.99):
                q=dict(model,context=context,operation='P',T=tr*tc,P=0.,mode=1 if context=='public' else 3)
                q['id']=f"{context}:{model['record']}:P:{tr}"
                out.append(q)
            # All 49 independently located critical fits participate in inverse
            # calls; unavailable critical metadata is not fabricated as default.
            if model['refine']:
                for pr in (.01,.2,.8):
                    q=dict(model,context=context,operation='T',T=0.,P=pr*pc,mode=2 if context=='public' else 4)
                    q['id']=f"{context}:{model['record']}:T:{pr}"
                    out.append(q)
    assert len(out)==618 and len({q['id'] for q in out})==618
    return out


def physical_requests():
    out=[]
    for model in models():
        tc,_,rho=model['critical']
        for tr in (.65,.8,.95,.99,.9999,1.0001):
            q=dict(model,context='physical',operation='P',T=tr*tc,P=0.)
            q['maximum_density']=(4. if tr<.99 else 2.)*rho
            q['id']=f"physical:{model['record']}:P:{tr}"
            out.append(q)
    assert len(out)==390
    return out


def native_key(q):
    tc,_,rho=q['critical'];maximum=q.get('maximum_density',4.*rho)
    return '|'.join(map(str,[q['record'],q['operation'],q['T'] if q['operation']=='P' else q['P'],
                             maximum,.5*tc,.9999*tc,0,0]))


def source_input_hashes():
    names=('mbwr_critical.csv','mbwr_initialization.csv','mbwr_initialization_status.csv')
    result={f'tests/fixtures/{name}':hashlib.sha256((ROOT/'tests/fixtures'/name).read_bytes()).hexdigest() for name in names}
    for name in ('mbwr_coexistence_cases.py','generate_data.py','generate_mbwr_data.py','validate_mbwr_density.py'):
        result['tools/rust_port/'+name]=digest(ROOT/'tools/rust_port'/name)
    return result
