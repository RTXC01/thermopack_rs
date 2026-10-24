"""Pure MEOS TP/PH/PS qualification with retained raw source outcomes.

Coordinate differences are reported explicitly. --allow-known-gaps permits
reporting a partially qualified surface; it never reclassifies a difference
as agreement. All source attempts remain in the metadata audit trail.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
import io
import tempfile

from reference import ROOT
from deadline_reference import DeadlineReferenceProcess
from pathlib import Path


def atomic_text(path,text):
    # Windows indexers may map a large fixture, preventing truncation. Replace
    # a closed sibling file atomically so readers keep a complete old snapshot.
    with tempfile.NamedTemporaryFile(mode='w',encoding='utf-8',newline='\n',dir=path.parent,
                                     prefix='.'+path.name+'.',suffix='.tmp',delete=False) as handle:
        temporary=Path(handle.name);handle.write(text)
    try:temporary.replace(path)
    finally:temporary.unlink(missing_ok=True)


def atomic_csv(path,header,rows):
    buffer=io.StringIO(newline='')
    writer=csv.writer(buffer,delimiter=';',lineterminator='\n')
    writer.writerow(header);writer.writerows(rows)
    atomic_text(path,buffer.getvalue())


def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def finite(value):
    return isinstance(value,(int,float)) and math.isfinite(value)


def valid_native_point(point):
    if len(point)!=12 or not all(map(finite,point)):return False
    t,p,two,beta,vl,vv,_,_,raw_two,raw_beta,raw_vl,raw_vv=point
    if min(t,p,vl,vv,raw_vl,raw_vv)<=0. or two not in (0.,1.) or raw_two not in (0.,1.):return False
    if raw_two:
        if not 0.<=raw_beta<=1. or (raw_vv-raw_vl)/(raw_vv+raw_vl)<1e-6:return False
    elif raw_beta not in (-1.,0.,1.):return False
    canonical=bool(raw_two) and 1e-8<raw_beta<1.-1e-8
    return bool(two)==canonical and beta==(raw_beta if canonical else 0.)


def gap_key(q):
    with (ROOT/'tests/fixtures/meos_metadata.csv').open() as handle:
        row=next(r for r in csv.DictReader(handle,delimiter=';') if r['case'].split('|')[1]==f"{q['eos']}:{q['component']}")
    tc=float(row['Tc']);pressure_scale=float(row['R'])*tc*float(row['rhoc'])
    return '|'.join(map(str,[q['eos'],q['component'],q['kind'],round(q['temperature']/tc,6),
                             round(q['pressure']/pressure_scale,8),q.get('kernel_quality',q.get('quality')),q.get('guess_factor'),
                             round(q['minimum_temperature']/tc,6),round(q['maximum_temperature']/tc,6)]))


def specified_pressure(q,r):
    # Saturation-quality requests intentionally replace the placeholder P.
    # Ordinary TP and direct caloric specifications must retain the input P.
    return r['seed']['pressure'] if 'quality' in q and r.get('seed') else q['pressure']


def bounded_endpoint_gap(e):
    """Allow only tiny trace-fraction changes with matching bulk state."""
    if e['native'].startswith('ERROR'):return False
    source=e['source'];a=source['values'];b=list(map(float,e['native'].split(',')))
    if b[8] and (b[11]-b[10])/(b[11]+b[10])<1e-6:return False
    for i in (0,1,6,7):
        floor=[2e-7,1e-4,0.,2e-7,1e-12,1e-12,2e-5,2e-7][i]
        if abs(a[i]-b[i])>floor+2e-6*max(abs(a[i]),abs(b[i])):return False
    source_volume=(1.-a[3])*a[4]+a[3]*a[5] if a[2] else a[4]
    native_volume=(1.-b[9])*b[10]+b[9]*b[11] if b[8] else b[4]
    if abs(source_volume-native_volume)>1e-12+2e-6*max(abs(source_volume),abs(native_volume)):return False
    beta=b[9] if b[8] else e['request'].get('kernel_quality')
    return beta is not None and abs(source['raw_result'][3]-beta)<=2e-6


def requests(smoke=False):
    with (ROOT/'tests/fixtures/meos_metadata.csv').open() as handle:
        rows=list(csv.DictReader(handle,delimiter=';'))
    for row in rows:
        eos,component=row['case'].split('|')[1].split(':')
        if smoke and component not in ('C3','N-H2','LJTS'):continue
        tc,pc=float(row['Tc']),float(row['Pc'])
        if pc==0.:
            # Dimensionless critical pressure obtained from the identical NE
            # LJTS fit. This is the stored-parameter pressure, not a refinement.
            ne=next(r for r in rows if 'LJTS_MEOS:NE|' in r['case'])
            vscale=6.02214076e23*float(ne['sigma'])**3
            pc=float(ne['Pc'])*vscale/(float(ne['R'])*float(ne['eps_div_k']))
        base=dict(mode='MEOS_FLASH',eos=eos,component=component,
                  minimum_temperature=.3*tc,maximum_temperature=3.*tc)
        points=[(.55,.03),(.8,.03),(.8,.5),(.999,.5),(1.001,1.1),(1.2,.03),(1.2,1.1),(2.,3.)]
        if smoke:points=[(.8,.03),(1.2,1.1)]
        for tr,pr in points:
            for kind in ('TP','PH','PS'):
                q=dict(base,kind=kind,temperature=tr*tc,pressure=pr*pc)
                yield q
                if kind!='TP':yield dict(q,guess_factor=1.2)
        for tr in ([.8] if smoke else [.5,.8,.95,.999]):
            for quality in ([.5] if smoke else [0.,.1,.5,.9,1.]):
                for kind in ('PH','PS'):
                    q=dict(base,kind=kind,temperature=tr*tc,pressure=pc,quality=quality)
                    yield q
                    if not smoke:yield dict(q,guess_factor=1.2)
            yield dict(base,kind='TP',temperature=tr*tc,pressure=pc,quality=.5)
        if not smoke:
            # Original physical default bounds, including hydrogen/reduced
            # saturation below Tmin. The lever rule bypasses these bounds.
            tscale=1.
            if eos=='LJTS_REDUCED':
                tscale=float(next(r for r in rows if 'LJTS_MEOS:NE|' in r['case'])['eps_div_k'])
            for kind in ('PH','PS'):
                yield dict(base,kind=kind,temperature=.8*tc,pressure=pc,quality=.5,
                           minimum_temperature=80./tscale,maximum_temperature=999./tscale)
            # Independently qualified coexistence coordinates supply source
            # TV-kernel caloric targets. This exercises the actual PH/PS calls
            # even if the source bubble-pressure seed cannot be generated.
            seeds=[]
            for filename in ('meos_saturation.csv','meos_saturation_recovered.csv'):
                with (ROOT/'tests/fixtures'/filename).open() as handle:
                    for fields in list(csv.reader(handle,delimiter=';'))[1:]:
                        parts=fields[0].split('|')
                        if parts[1]!=f'{eos}:{component}' or parts[5]!='SAT_P':continue
                        t,p,vl,vv=map(float,fields[1:5])
                        if any(abs(t/tc-tr)<1e-10 for tr in (.999,.9999)):
                            seeds.append((t,p,vl,vv))
            assert len(seeds)==2,(eos,component,seeds)
            for t,p,vl,vv in seeds:
                for quality in [0.,.5,1.]:
                    for kind in ('PH','PS'):
                        q=dict(base,kind=kind,temperature=t,pressure=p,
                               kernel_seed_volumes=[vl,vv],kernel_quality=quality)
                        yield q
                        yield dict(q,guess_factor=1.2)


def key(request,result):
    q=request
    guess='none' if 'guess_factor' not in q else repr(q['temperature']*q['guess_factor'])
    target=result.get('specification',q.get('target',0.))
    p=specified_pressure(q,result)
    control=f"FLASH:{q['kind']}:{target}:{guess}:{q['minimum_temperature']}:{q['maximum_temperature']}"
    return f"MEOS|{q['eos']}:{q['component']}|{q['temperature']}|{p}|1|{control}"


def validate(args):
    identity=json.loads((ROOT/'.reference/metadata.json').read_text())
    signature=dict(reference=identity,oracle_sha256=digest(ROOT/'tools/rust_port/oracle_worker.py'),
                   generator_sha256=digest(Path(__file__)),smoke=args.smoke,
                   metadata_fixture_sha256=digest(ROOT/'tests/fixtures/meos_metadata.csv'),
                   probe_sha256=hashlib.sha256((ROOT/'.reference'/('reference_probe.dll' if os.name=='nt' else 'libreference_probe.so')).read_bytes()).hexdigest(),
                   source_seed_fixtures={name:digest(ROOT/'tests/fixtures'/name) for name in ('meos_saturation.csv','meos_saturation_recovered.csv')})
    cache=ROOT/'.reference'/('meos_flash_smoke_source.json' if args.smoke else 'meos_flash_source.json')
    if args.reuse_source:
        saved=json.loads(cache.read_text());assert saved['signature']==signature,'stale source cache'
        records=saved['records']
    else:
        records=[];kernel_seeds={}
        with DeadlineReferenceProcess(args.compiler_directory) as oracle:
            previous=None
            for q in requests(args.smoke):
                model=(q['eos'],q['component'])
                if model!=previous:
                    print(*model,len(records),'source attempts',flush=True);previous=model
                if 'kernel_seed_volumes' in q:
                    seed_key=(*model,q['temperature'])
                    if seed_key not in kernel_seeds:
                        points=[]
                        for volume in q['kernel_seed_volumes']:
                            checked=oracle.call(dict(mode='MEOS',eos=q['eos'],component=q['component'],temperature=q['temperature'],volume=volume))
                            assert 'values' in checked,checked
                            points.append(checked['values'])
                        kernel_seeds[seed_key]=points
                    index=17 if q['kind']=='PH' else 16
                    left,right=kernel_seeds[seed_key];beta=q['kernel_quality']
                    q['target']=(1.-beta)*left[index]+beta*right[index]
                    q['source_kernel_caloric_endpoints']=[left[index],right[index]]
                records.append(dict(request=q,result=oracle.call(q,timeout=4.)))
        atomic_text(cache,json.dumps(dict(signature=signature,records=records),allow_nan=False)+'\n')
    qualified=[];unqualified=[];errors=[];censored=[]
    for e in records:
        r=e['result']
        if 'reference_timeout' in r:censored.append(e)
        elif 'calculation_error' in r or 'reference_error' in r or not r.get('values'):errors.append(e)
        elif (not r['pressure_stable'] or not finite(r['scaled_specification_error']) or abs(r['scaled_specification_error'])>2e-7
              or abs(r['values'][1]-specified_pressure(e['request'],r))>8.*math.ulp(1.)*abs(specified_pressure(e['request'],r))
              or r['material_balance_error']>1e-8
              or not all(map(finite,r['raw_result']))
              or (r['raw_result'][2]==0 and not all(0.<=v<=1. for v in r['raw_result'][3:5]))
              or not all(map(finite,r['values']))):unqualified.append(e)
        else:qualified.append(e)
    subprocess.run(['cargo','build','--release','--example','parity','-j','2'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    recoverable=[e for e in errors+unqualified if finite(e['result'].get('specification',e['request'].get('target')))]
    evaluated=qualified+recoverable
    cases=[key(e['request'],e['result']) for e in evaluated]
    output=subprocess.run([str(exe)],input='\n'.join(cases)+'\n',capture_output=True,text=True,check=True).stdout.splitlines()
    assert len(output)==len(cases),(len(output),len(cases))
    matches=[];differences=[];maximum=[0.]*8;native_points=[];degenerate=[];recoveries=[];recovery_errors=[]
    for number,(case,e,encoded) in enumerate(zip(cases,evaluated,output)):
        recovering=number>=len(qualified)
        expected=e['result'].get('values',[]);bad=[]
        if encoded.startswith('ERROR'):bad=[encoded]
        else:
            actual=list(map(float,encoded.split(',')));assert valid_native_point(actual),(case,'invalid native phase data',actual)
            native_points.append((case,e,actual))
            if recovering:
                recoveries.append([case,*actual[:8],'source_error' if 'calculation_error' in e['result'] or 'reference_error' in e['result'] else 'source_nonclosing'])
                continue
            for i,(a,b) in enumerate(zip(expected,actual)):
                floor=[2e-7,1e-4,0.,2e-7,1e-12,1e-12,2e-5,2e-7][i]
                maximum[i]=max(maximum[i],abs(a-b)/max(abs(a),abs(b),1e-300))
                if not math.isfinite(b) or abs(a-b)>floor+2e-6*max(abs(a),abs(b)):bad.append([i,a,b])
        if recovering:
            recovery_errors.append(dict(case=case,native=encoded,source=e['result']));continue
        if e['request']['kind']=='TP' and 'quality' in e['request']:
            # At a pure coexistence T/P, quality and hence H/S/V are not
            # specified. Keep all raw endpoint choices, and separately check
            # the returned native phase states against source kernels.
            degenerate.append(dict(case=case,request=e['request'],source=e['result'],native=encoded,differences=bad))
        elif bad:differences.append(dict(case=case,request=e['request'],source=e['result'],native=encoded,differences=bad))
        else:matches.append([case,*expected])
    # Evaluate returned native phases through independent original TV kernels;
    # this does not rely on the original TP/PH/PS solver reaching the same root.
    independent=[];kernel_cache={}
    with DeadlineReferenceProcess(args.compiler_directory) as oracle:
        for case,e,point in native_points:
            t,p,_,_,_,_,h,s,two,beta,vl,vv=point;q=e['request'];states=[]
            for volume in ([vl,vv] if two else [vl]):
                state_key=(q['eos'],q['component'],t,volume)
                if state_key not in kernel_cache:
                    checked=oracle.call(dict(mode='MEOS',eos=q['eos'],component=q['component'],temperature=t,volume=volume))
                    assert 'values' in checked,(case,checked)
                    kernel_cache[state_key]=checked
                states.append(kernel_cache[state_key])
            r=states[0]['metadata'][10];weights=[1.-beta,beta] if two else [1.]
            source_h=sum(w*st['values'][17] for w,st in zip(weights,states))
            source_s=sum(w*st['values'][16] for w,st in zip(weights,states))
            if two:
                assert (vv-vl)/(vv+vl)>=1e-6,(case,'coincident native phases')
                assert abs(states[0]['values'][22]-states[1]['values'][22])<2e-8,(case,'native fugacity closure')
            caloric_errors=[(h-source_h)/(r*t),(s-source_s)/r]
            assert max(map(abs,caloric_errors))<2e-10,(case,caloric_errors)
            pressure_errors=[]
            for volume,st in zip(([vl,vv] if two else [vl]),states):
                computed=st['values'][12];slope=st['values'][14]
                assert finite(computed) and finite(slope) and slope<0.,(case,computed,slope)
                allowance=2e-8*p+512.*math.ulp(1.)*r*t/volume
                assert abs(computed-p)<allowance,(case,computed,p,allowance)
                pressure_errors.append(computed-p)
            target=e['result'].get('specification',q.get('target'))
            spec_error=None if target is None else ((source_h-target)/(r*t) if q['kind']=='PH' else (source_s-target)/r)
            assert spec_error is None or abs(spec_error)<2e-10,(case,spec_error)
            independent.append(dict(case=case,pressure_errors=pressure_errors,
                                    caloric_kernel_errors=caloric_errors,scaled_specification_error=spec_error))
    report=dict(signature,attempts=len(records),qualified_source_returns=len(qualified),matched=len(matches),
                native_sources_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in sorted((ROOT/'rust').rglob('*.rs'))+[ROOT/'Cargo.toml',ROOT/'Cargo.lock']},
                reference_errors=errors,censored_calls=censored,unqualified_returns=unqualified,
                native_differences=differences,maximum_relative_differences=maximum,
                records=records,independent_native_kernel_checks=independent,relative_tolerance=2e-6,
                qualification_status='partial_surface_with_explicit_source_outcomes_and_endpoint_gaps',
                coexistence_tp_endpoint_choices=degenerate,
                recovered_native_states=len(recoveries),native_recovery_errors=recovery_errors,
                absolute_tolerances=[2e-7,1e-4,0.,2e-7,1e-12,1e-12,2e-5,2e-7],
                qualification='Canonical phase properties of DEFAULT caloric reference states; raw phase flags retained. Reduced LJTS is a scaled physical NE LJTS reference, not a public reduced flash wrapper. Difference records are explicit remaining gaps; matched fixture alone is not full surface parity.')
    atomic_text(ROOT/'.reference/meos_flash_report.json',json.dumps(report,indent=2,allow_nan=False)+'\n')
    if not args.smoke:
        fixture=ROOT/'tests/fixtures/meos_flash.csv'
        atomic_csv(fixture,['case','T','P','interior_two_phase','quality','liquid_or_overall_volume','vapor_or_overall_volume','H','S'],matches)
        report.update(cases=len(matches),fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),generator='tools/rust_port/validate_meos_flash.py')
        atomic_text(fixture.with_suffix('.metadata.json'),json.dumps(report,indent=2,allow_nan=False)+'\n')
        fixture=ROOT/'tests/fixtures/meos_flash_coexistence_tp.csv'
        atomic_csv(fixture,['case','T','P','interior_two_phase','quality','liquid_or_overall_volume','vapor_or_overall_volume','H','S'],
                   ([e['case'],*e['source']['values']] for e in degenerate))
        atomic_text(fixture.with_suffix('.metadata.json'),json.dumps(dict(reference=identity,
            cases=len(degenerate),fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),generator='tools/rust_port/validate_meos_flash.py',
            qualification='Pure coexistence TP states: native phase membership and source-kernel closure, not source endpoint/quality trajectory parity. Full raw endpoint choices retained in meos_flash.metadata.json.'),indent=2)+'\n')
        fixture=ROOT/'tests/fixtures/meos_flash_recovered.csv'
        atomic_csv(fixture,['case','T','P','interior_two_phase','quality','liquid_or_overall_volume','vapor_or_overall_volume','H','S','source_outcome'],recoveries)
        atomic_text(fixture.with_suffix('.metadata.json'),json.dumps(dict(reference=identity,
            cases=len(recoveries),fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),generator='tools/rust_port/validate_meos_flash.py',
            qualification='Native solutions checked against original pressure/caloric TV kernels where public PH/PS fails or returns a nonclosing specification. These are not original solver coordinate parity claims.'),indent=2)+'\n')
    print(json.dumps(dict(attempts=len(records),source_qualified=len(qualified),matched=len(matches),
                         differences=len(differences),source_errors=len(errors),unqualified=len(unqualified),censored=len(censored),maximum=maximum),indent=2),flush=True)
    for e in differences[:6]:print(json.dumps(dict(case=e['case'],differences=e['differences'])),flush=True)
    # A known endpoint fraction difference must never turn into a permitted
    # native solver error. Regenerating fixtures cannot bless new gaps.
    assert not [e for e in differences if e['native'].startswith('ERROR')], 'native errors on qualified source cases'
    assert not recovery_errors, 'native recovery regressions'
    if differences:
        if not args.allow_known_gaps:raise SystemExit(1)
        allowed=json.loads((ROOT/'tools/rust_port/meos_flash_known_gaps.json').read_text())['case_keys']
        unexpected=[e['case'] for e in differences if gap_key(e['request']) not in allowed or not bounded_endpoint_gap(e)]
        assert not unexpected,('new coordinate gaps',unexpected)
    assert not censored,'new censored source calls'
    if errors or unqualified:
        assert args.allow_known_gaps,'source errors/unqualified returns require explicit partial-surface mode'
        allowed=json.loads((ROOT/'tools/rust_port/meos_flash_known_gaps.json').read_text())['source_outcomes']
        unexpected=[]
        for outcome,entries in [('error',errors),('unqualified',unqualified)]:
            for e in entries:
                tag='|'.join([gap_key(e['request']),outcome,e['result'].get('stage','terminal'),e['result'].get('calculation_error','')])
                if tag not in allowed:unexpected.append(tag)
        assert not unexpected,('new source outcomes',unexpected)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--compiler-directory')
    parser.add_argument('--reuse-source',action='store_true');parser.add_argument('--smoke',action='store_true')
    parser.add_argument('--allow-known-gaps',action='store_true');validate(parser.parse_args())
