"""Audit compressed MBWR caloric evidence and golden vectors without Fortran."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re

from reference import ROOT
from qualify_mbwr_caloric import FIELDS,RELATIVE,assemble,differences,key


def audit(report):
    records=report['records']; schedule=report['signature']['schedule']
    assert len(records)==len(schedule)==report['attempts']==report['cases']==5994
    assert [e['request'] for e in records]==schedule
    assert len({key(q) for q in schedule})==len(schedule)
    assert report['fields']==FIELDS and report['relative_tolerance']==RELATIVE
    assert not report['errors'] and report['source_process']['returncode']==0
    hashes=report['native_sources_sha256']|report['signature']['object_hashes']
    assert hashes and all(re.fullmatch('[a-f0-9]{64}',s) for s in hashes.values())
    assert re.fullmatch('[a-f0-9]{64}',report['native_executable_sha256'])
    raw={}
    for i,kind,body in re.findall(r'^CASE=(\d+) (VALUES|IDEAL|RAW_SOUND)=([^\r\n]+)',report['source_process']['stdout'],re.M):
        k=(int(i),kind);assert k not in raw;raw[k]=list(map(float,body.split()))
    assert len(raw)==3*len(records)
    counts=dict(default=0,explicit=0,reference_edits=0,type10=0,negative_pressure=0)
    for i,e in enumerate(records):
        q=e['request']; source=e['source_baseline']; ideal=e['source_ideal']
        assert len(source)==45 and len(ideal)==15 and all(map(math.isfinite,source+ideal))
        assert source==raw[i+1,'VALUES'] and ideal==raw[i+1,'IDEAL']
        assert all(source[j] in (0.,1.) for j in (39,41,43))
        f=source[:10].copy()
        if ideal[3]==10.:
            f[4]-=q['n']*999.*ideal[4]/(8.31446261815324*q['T']**2)
            counts['type10']+=1
        expected=assemble(q,f,source[3]-ideal[8],ideal[2]*1e-3)
        assert expected==e['physical_from_original_jet']
        assert len(e['native'])==2
        for target,actual in zip((source,expected),e['native']):
            assert isinstance(actual,list) and not differences(q,target,actual)
            assert all(actual[j] in (0.,1.) for j in (39,41,43))
        counts['default' if q['binding']=='DEFAULT' else 'explicit']+=1
        counts['reference_edits']+=q['profile']!=0
        counts['negative_pressure']+=source[10]<=0
    assert counts==report['counts']
    matched=failed=0
    for e in report['public_constructor_calls']:
        assert not e.get('censored') and 'PUBLIC_BEGIN' in e['stdout']
        if 'PUBLIC_RETURNED' not in e['stdout']:
            assert e['expected_initialization_status']==1;failed+=1
        else:
            assert e['returncode']==0
            match=re.search(r'^CASE=1 VALUES=([^\r\n]+)',e['stdout'],re.M);assert match
            values=list(map(float,match[1].split()));assert values==e['values']
            target=next(r for r in records if r['request']['record']==e['record'] and r['request']['T']==298.15
                        and r['request']['V']==.001 and r['request']['n']==1. and r['request']['profile']==0)
            assert not differences(target['request'],target['source_baseline'],values)
            matched+=1
    assert matched==report['public_tv_matches']==30
    assert failed==report['public_constructor_failures']==16
    return dict(cases=len(records),native_views=2*len(records),public_matches=matched,public_constructor_failures=failed,**counts)


def load_and_audit():
    fixture=ROOT/'tests/fixtures/mbwr_caloric.csv'
    metadata=json.loads(fixture.with_suffix('.metadata.json').read_text(encoding='utf-8'))
    assert hashlib.sha256(fixture.read_bytes()).hexdigest()==metadata['fixture_sha256']
    archive=fixture.parent/metadata['evidence_file']
    assert archive.parent==fixture.parent and archive.name=='mbwr_caloric.audit.json.gz'
    data=archive.read_bytes();assert hashlib.sha256(data).hexdigest()==metadata['evidence_sha256']
    report=json.loads(gzip.decompress(data))
    assert report['fixture_sha256']==metadata['fixture_sha256']
    result=audit(report)
    for field in ('cases','attempts','fields','counts','native_sources_sha256','native_executable_sha256','public_tv_matches','public_constructor_failures'):
        assert report[field]==metadata[field]
    with fixture.open(encoding='utf-8',newline='') as handle:
        rows=list(csv.DictReader(handle,delimiter=';'))
    assert len(rows)==len(report['records'])
    for row,e in zip(rows,report['records']):
        assert row['case']==key(e['request'])
        assert list(map(float,row['source_baseline'].split(',')))==e['source_baseline']
        assert list(map(float,row['physical_from_original_jet'].split(',')))==e['physical_from_original_jet']
    return report,result


if __name__=='__main__':
    _,result=load_and_audit();print(json.dumps(result,sort_keys=True))
