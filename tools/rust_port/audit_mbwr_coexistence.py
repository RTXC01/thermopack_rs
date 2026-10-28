"""Offline integrity audit for MBWR coexistence evidence."""
import csv
import gzip
import hashlib
import json
import math
import re
from pathlib import Path

from reference import ROOT
from mbwr_coexistence_cases import reference_requests, physical_requests, native_key, source_input_hashes
from qualify_mbwr_coexistence import (RELATIVE, ABSOLUTE, clean, parse_original,
    equilibrium_checks, source_category, compare_state, screen_kernel, kernel_request,
    encode_request)

HASH_RE = re.compile(r'^[a-f0-9]{64}$')
EXPECTED = 1008


def check(ok, message):
    if not ok:
        raise AssertionError(message)


def _finite_vector(value, n, label):
    check(isinstance(value, list) and len(value) == n, f'{label} must contain {n} values')
    out=[]
    for x in value:
        if isinstance(x, str):
            check(x in ('nan','inf','-inf'), f'invalid nonfinite encoding in {label}')
            out.append({'nan':math.nan,'inf':math.inf,'-inf':-math.inf}[x])
        else:
            check(type(x) in (int,float), f'invalid number in {label}')
            out.append(float(x))
    return out


def _same_raw(a,b):
    return (math.isnan(a) and math.isnan(b)) or a == b


def _native_observation(result, label):
    """Bind decoded native values to the unmodified line emitted by the executable."""
    raw=result.get('stdout')
    check(isinstance(raw,str) and '\n' not in raw and '\r' not in raw,f'{label} missing raw native stdout')
    if 'error' in result:
        check(set(result)=={'error','stdout'},f'{label} native error record has unexpected fields')
        check(result['error']==raw and raw.startswith('ERROR'),f'{label} native error differs from stdout')
        return None
    check(set(result)=={'state','roots','stdout'},f'{label} native success record has unexpected fields')
    state=_finite_vector(result.get('state'),5,label+' state')
    roots=result.get('roots')
    check(isinstance(roots,list) and all(type(x) in (int,float) and math.isfinite(x) for x in roots),f'{label} roots invalid')
    try:
        body,root_body=raw.split(';')
        parsed_state=[float(x) for x in body.split(',')]
        parsed_roots=[float(x) for x in root_body.split(',')] if root_body else []
    except (ValueError,TypeError) as exc:
        raise AssertionError(f'{label} raw native stdout malformed') from exc
    check(len(parsed_state)==5 and len(parsed_roots)==len(roots),f'{label} raw native stdout shape mismatch')
    check(all(_same_raw(a,float(b)) for a,b in zip(state,parsed_state)),f'{label} native state differs from stdout')
    check(all(_same_raw(float(a),b) for a,b in zip(roots,parsed_roots)),f'{label} native roots differ from stdout')
    return state


def _native_root_query(proof, request, source_state):
    check(isinstance(proof,dict), 'missing original root-query trace')
    check(set(proof)=={'request','stdout','roots'},'original root-query trace has unexpected fields')
    maximum=4*request['critical'][2]
    query='|'.join(map(str,[request['record'],'R',source_state[0],maximum,source_state[1],0,0,0]))
    check(proof.get('request')==query,'original root-query request differs from source STATE')
    raw=proof.get('stdout')
    check(isinstance(raw,str) and '\n' not in raw and '\r' not in raw,'missing raw original root-query stdout')
    roots=proof.get('roots')
    if raw.startswith('ERROR:'):
        check(roots==raw,'root-query error differs from raw stdout')
        return None
    check(raw.startswith('ROOTS;'),'malformed original root-query stdout')
    body=raw[len('ROOTS;'):]
    parsed=[float(x) for x in body.split(',')] if body else []
    check(isinstance(roots,list) and len(roots)==len(parsed),'original root list shape differs from stdout')
    check(all(type(x) in (int,float) and math.isfinite(x) and _same_raw(float(x),p) for x,p in zip(roots,parsed)),
          'original roots differ from raw root-query stdout')
    return parsed


def _source_binding(source):
    check(not source.get('censored'), 'source collection was censored')
    check(type(source.get('returncode')) is int, 'source process return code is missing')
    out=source.get('stdout','')
    check(isinstance(out,str),'source stdout is missing')
    for label in ('META','STATE','CHECK'):
        check(len(re.findall(r'^'+label+r'=',out,re.M))<=1,f'duplicate source {label} line')
    check(len(re.findall(r'^POINT_STATUS=',out,re.M))<=1,'duplicate source POINT_STATUS line')
    check(out.count('INIT_RETURNED')<=1,'duplicate source initialization marker')
    parsed=parse_original(source)
    for label,n in (('META',10),('STATE',5),('CHECK',12)):
        if label in parsed:
            _finite_vector(clean(parsed[label]),n,'parsed '+label)
    if 'TP' in parsed:
        check(len(parsed['TP'])==2,'original TP point must contain two values')
        if 'STATE' in parsed and all(math.isfinite(x) for x in parsed['TP']+parsed['STATE'][:2]):
            check(parsed['TP']==parsed['STATE'][:2],'original TP point differs from reported STATE')
    return parsed


def _sig_hashes(report):
    sig=report.get('signature',{})
    check(isinstance(sig,dict) and sig, 'missing source signature')
    for field in ('compiler_sha256','helper_sha256','oracle_module_sha256','executable_sha256','reference_library_sha256'):
        check(HASH_RE.fullmatch(str(sig.get(field,''))) is not None, f'invalid signature hash {field}')
    for field in ('source_hashes','object_hashes','dependency_hashes'):
        values=sig.get(field)
        check(isinstance(values,dict) and values and all(HASH_RE.fullmatch(str(x)) for x in values.values()),f'invalid {field}')
    check(isinstance(sig.get('flags'),list) and sig.get('environment',{}).get('OMP_NUM_THREADS')=='1','invalid compiler execution signature')
    cohort=sig.get('cohort_inputs_sha256')
    check(isinstance(cohort,dict) and cohort and all(HASH_RE.fullmatch(str(x)) for x in cohort.values()),'invalid authoritative cohort input hashes')
    for field in ('native_sources_sha256',):
        values=report.get(field,{})
        check(isinstance(values,dict) and values and all(HASH_RE.fullmatch(str(x)) for x in values.values()),f'invalid {field}')
    check(HASH_RE.fullmatch(str(report.get('native_executable_sha256',''))) is not None,'invalid native executable hash')


def _check_source_rows(report, reqs):
    process=report.get('source_processes', report.get('source_process'))
    if isinstance(process,list):
        check(len(process)==len(reqs),'source process count mismatch')
        return [_source_binding(x) for x in process]
    # Current qualifier stores the isolated observations on each source record.
    recs=report.get('records',[])
    source_rows=[r.get('source') for r in recs if isinstance(r,dict) and 'source' in r]
    check(len(source_rows)==len(reqs),'source process count mismatch')
    return [_source_binding(x) for x in source_rows]


def _check_fixture_rows(rows, report):
    check(len(rows)==EXPECTED,'CSV row count mismatch')
    expected_columns=['case','request','original_solver','native_observed','classification']
    if rows:
        check(list(rows[0].keys())==expected_columns,'CSV columns differ from required evidence schema')
    for row,e in zip(rows,report['records']):
        q=e['request']; check(row['case']==q['id'] and row['request']==native_key(q),'CSV request binding mismatch')
        check(row['classification']==e['classification'],'CSV classification mismatch')
        want=e['original']['STATE'] if e['classification']=='original_solver_match' else ''
        check(row['original_solver']==(want if isinstance(want,str) else ','.join(map(repr,want))),'CSV source gold mismatch')
        native=e['native']; got=','.join(map(repr,native['state'])) if 'state' in native else native['error']
        check(row['native_observed']==got,'CSV native observation mismatch')


def audit(report):
    physical=physical_requests(); refs=reference_requests(); allreq=physical+refs
    records=report.get('records',[])
    check(report.get('cases')==EXPECTED and len(records)==EXPECTED,'case count is not 1008')
    check(report.get('physical_requests')==physical and report.get('reference_requests')==refs,'authoritative cohort requests differ')
    check([e.get('request') for e in records]==allreq,'record requests differ from authoritative cohorts')
    check(report.get('relative_tolerance')==RELATIVE and report.get('absolute_tolerances')==ABSOLUTE,'comparison tolerances differ')
    check(report.get('signature',{}).get('cohort_inputs_sha256')==source_input_hashes(),'cohort input hashes differ from current authoritative inputs')
    _sig_hashes(report)
    ks=report.get('kernel_signature')
    check(isinstance(ks,dict) and ks,'missing original kernel signature')
    for field in ('compiler_sha256','helper_sha256','oracle_module_sha256','executable_sha256','reference_library_sha256'):
        check(HASH_RE.fullmatch(str(ks.get(field,''))) is not None,f'invalid kernel signature hash {field}')
    for field in ('source_hashes','object_hashes','dependency_hashes'):
        values=ks.get(field)
        check(isinstance(values,dict) and values and all(HASH_RE.fullmatch(str(x)) for x in values.values()),f'invalid kernel {field}')
    check(isinstance(ks.get('flags'),list) and ks.get('environment',{}).get('OMP_NUM_THREADS')=='1','invalid kernel compiler execution signature')
    check(HASH_RE.fullmatch(str(report.get('generator_sha256',''))) is not None,'invalid qualifier generator hash')
    kernels=report.get('kernel_requests',[]); observations=report.get('kernel_processes',[])
    check(len(kernels)==len(observations),'kernel request/process count mismatch')
    kernel_raw={}
    for q,o in zip(kernels,observations):
        check(q.get('mode')==5 and q.get('id'),'invalid kernel request')
        check(q['id']=='kernel:'+hashlib.sha256(encode_request(q).encode()).hexdigest(),'kernel ID does not bind encoded request')
        check(o.get('returncode')==0,'original kernel process failed')
        parsed=_source_binding(o); kernel_raw[q['id']]=(q,parsed)
    check(len({q['id'] for q in kernels})==len(kernels),'duplicate kernel query')
    counts={}
    used_kernel_ids=set()
    for i,(e,q) in enumerate(zip(records,allreq)):
        check(e.get('request')==q,f'record {i} request changed')
        result=e.get('native',{})
        native_state=_native_observation(result,f'record {i}')
        native_kid=e.get('kernel') if q['context']=='physical' else e.get('native_kernel')
        if native_state is not None:
            check(native_kid in kernel_raw,f'record {i} missing original kernel check for native state')
            nkq,nkp=kernel_raw[native_kid]
            used_kernel_ids.add(native_kid)
            check(nkq['T']==native_state[0] and nkq['P']==native_state[1] and nkq['liquid_volume']==native_state[2] and nkq['vapor_volume']==native_state[3],f'record {i} native kernel request differs from native state')
            check(nkq['roots']==result.get('roots'),f'record {i} native kernel roots differ from native roots')
            check(encode_request(nkq)==encode_request(kernel_request(q,native_state,result['roots'])),f'record {i} kernel request is not the exact native state and roots')
            check(not screen_kernel(nkq,nkp),f'record {i} original kernel fails native state')
            check(not compare_state(nkp['STATE'],native_state),f'record {i} native state differs from original kernel STATE')
            counts['native_kernel_checked']=counts.get('native_kernel_checked',0)+1
        else:
            check(native_kid is None,f'record {i} has kernel check without native state')
            counts['native_error_observations']=counts.get('native_error_observations',0)+1
        if q['context']=='physical':
            classification='kernel_checked_native' if native_state is not None else 'native_error_unqualified'
            check(e.get('classification')==classification,f'record {i} physical classification is wrong')
        else:
            source=e.get('source')
            check(isinstance(source,dict),f'record {i} missing original process')
            parsed=_source_binding(source)
            stored=e.get('original')
            check(stored==clean(parsed),f'record {i} stored source fields differ from raw stdout')
            if q['mode']>=3 and 'META' in parsed:
                # Internal calls explicitly assign the MBWR critical constants
                # from the signed request before META is emitted.
                for j,(actual,expected,floor) in enumerate(zip(parsed['META'][4:8],q['parameters'],(1e-10,1e-7,1e-10,1e-12))):
                    check(abs(actual-expected)<=floor+1e-12*max(abs(actual),abs(expected)),f'record {i} internal META parameter {j} differs from request')
            if parsed.get('status')==0 and parsed.get('STATE') and all(math.isfinite(x) for x in parsed['STATE'][:2]):
                expected_input=q['T'] if q['operation']=='P' else q['P']
                observed_input=parsed['STATE'][0] if q['operation']=='P' else parsed['STATE'][1]
                check(observed_input==expected_input,f'record {i} original solver state differs from requested input')
            category=source_category(source,parsed)
            classification=category
            kid=e.get('reference_kernel')
            if category=='kernel_candidate':
                root_proof=e.get('root_query')
                native_roots=_native_root_query(root_proof,q,parsed['STATE'])
                check(native_roots is not None,f'record {i} original root query failed')
                check(e.get('reference_kernel') is not None,f'record {i} missing original kernel check')
                check(kid in kernel_raw,f'record {i} missing original kernel check')
                kq,kp=kernel_raw[kid]
                used_kernel_ids.add(kid)
                check(kq['record']==q['record'] and kq['T']==parsed['STATE'][0] and kq['P']==parsed['STATE'][1],f'record {i} kernel request state changed')
                check(kq['liquid_volume']==parsed['STATE'][2] and kq['vapor_volume']==parsed['STATE'][3],f'record {i} kernel volumes changed')
                check(kp['STATE'][:4]==parsed['STATE'][:4],f'record {i} kernel STATE differs from original source STATE')
                check(kp['CHECK']==parsed['CHECK'],f'record {i} kernel CHECK differs from original source CHECK')
                check(kq['roots']==native_roots,f'record {i} kernel roots differ from original root-query result')
                check(encode_request(kq)==encode_request(kernel_request(q,parsed['STATE'],native_roots)),f'record {i} original kernel request differs from source STATE and queried roots')
                kernel_failures=screen_kernel(kq,kp)
                if kernel_failures==['lower original root Gibbs energy']:
                    classification='original_metastable_in_domain'
                else:
                    check(not kernel_failures,f'record {i} original kernel checks fail: {kernel_failures}')
                maximum=4*q['critical'][2]
                if classification!='original_metastable_in_domain' and (1/parsed['STATE'][2]>maximum or 1/parsed['STATE'][3]>maximum): classification='original_outside_density_domain'
                elif classification!='original_metastable_in_domain' and q['operation']=='T' and not .5*q['critical'][0]<=parsed['STATE'][0]<=.9999*q['critical'][0]: classification='original_outside_temperature_domain'
                elif classification!='original_metastable_in_domain':
                    classification='original_solver_match'
                    diffs=compare_state(parsed['STATE'],result.get('state'))
                    check(not diffs,f'record {i} native differs from original solver')
            check(e.get('classification')==classification,f'record {i} stored source classification is wrong')
            counts[classification]=counts.get(classification,0)+1
    check(used_kernel_ids==set(kernel_raw),'kernel query set contains unreferenced or missing evidence')
    check(report.get('counts')==counts,'classification counts differ from recomputation')
    check(not report.get('errors'),'qualifier report contains errors')
    return {'cases':len(records),'counts':counts,'kernel_checks':len(kernels)}


def load_and_audit():
    fixture=ROOT/'tests/fixtures/mbwr_coexistence.csv'; metadata_path=fixture.with_suffix('.metadata.json')
    metadata=json.loads(metadata_path.read_text(encoding='utf-8')); archive=fixture.with_suffix('.audit.json.gz')
    check(metadata.get('evidence_file')==archive.name,'metadata archive name mismatch')
    fb=fixture.read_bytes(); ab=archive.read_bytes()
    check(hashlib.sha256(fb).hexdigest()==metadata.get('fixture_sha256'),'fixture hash mismatch')
    check(hashlib.sha256(ab).hexdigest()==metadata.get('evidence_sha256'),'archive hash mismatch')
    report=json.loads(gzip.decompress(ab)); check(report.get('fixture_sha256')==metadata['fixture_sha256'],'report fixture hash mismatch')
    result=audit(report)
    for field in ('cases','counts','errors','native_sources_sha256','native_executable_sha256','relative_tolerance','absolute_tolerances','generator_sha256'):
        check(report.get(field)==metadata.get(field),f'metadata {field} mismatch')
    check(report.get('signature')==metadata.get('signature'),'metadata signature mismatch')
    check(report.get('kernel_signature')==metadata.get('kernel_signature'),'metadata kernel signature mismatch')
    with fixture.open(encoding='utf-8',newline='') as f: rows=list(csv.DictReader(f,delimiter=';'))
    _check_fixture_rows(rows,report)
    result['fixture_sha256']=metadata['fixture_sha256']; return report,result


if __name__=='__main__':
    _,result=load_and_audit(); print(json.dumps(result,sort_keys=True))
