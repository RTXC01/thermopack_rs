"""Audit paired historical native MEOS branch benchmark reports."""
import argparse,hashlib,json,math,statistics
from pathlib import Path
from reference import ROOT

BEFORE=ROOT/'docs/rust/meos-branches-before-cache-benchmark.json'
AFTER=ROOT/'docs/rust/meos-branches-benchmark.json'
OUTPUT=ROOT/'docs/rust/meos-bounds-cache-comparison.json'
UNCHANGED=('source_library_sha256','source_probe','platform','processor','rustc','profiles','numeric_cases','preflight_cases','iterations_per_sample','cargo_inputs_sha256','benchmark_sources_sha256')
CHANGED_NATIVE='rust/src/multiparameter_intervals.rs'

def load(path):
    report=json.loads(Path(path).read_text(encoding='utf-8'))
    assert isinstance(report,dict),f'{path}: expected JSON object'
    return report

def digest(value):
    payload=json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()
    return hashlib.sha256(payload).hexdigest()

def finite_positive(value,label):
    assert isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value) and value>0,label

def validate_report(report,label):
    profiles=report['profiles'];assert len(profiles)==5,f'{label}: expected five profiles'
    assert report['preflight_cases']==40,f'{label}: expected forty preflight cases'
    assert report['iterations_per_sample']['setup']==10,f'{label}: unexpected setup iterations'
    samples=report['samples_ns'];medians=report['median_ns']
    expected={f'{mode}_{i}' for i in range(1,6) for mode in ('setup','fresh','prepared')}
    assert set(samples)==set(medians)==expected,f'{label}: sample/median keys differ from profiles'
    for key,values in samples.items():
        assert isinstance(values,list) and len(values)==5,f'{label}: {key} must have five samples'
        for value in values:finite_positive(value,f'{label}: {key} has invalid sample')
        median=statistics.median(values);stored=medians[key]
        assert isinstance(stored,(int,float)) and math.isfinite(stored) and stored==median,f'{label}: incorrect median for {key}'
    for mode,iters in (('fresh',100),('prepared',1000)):
        assert report['iterations_per_sample'][mode]==iters,f'{label}: unexpected {mode} iterations'
        for i in range(1,6):
            values=report['checksums'][f'{mode}_{i}']
            assert isinstance(values,list) and len(values)==6,f'{label}: invalid {mode}_{i} checksum vector'
            assert all(isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x) for x in values),f'{label}: non-finite checksum'
            numeric_cases=report['numeric_cases'][str(i)]
            expected_count=sum(numeric_cases[n%len(numeric_cases)]['expected'][0] for n in range(iters))
            assert values[1]==expected_count,f'{label}: {mode}_{i} root-count sum mismatch'

def normalized(report):
    """Comparable input identity, omitting run outputs and the permitted source edit."""
    result={key:report[key] for key in UNCHANGED}
    native=report['native_sources_sha256'];assert isinstance(native,dict)
    result['native_sources_sha256']={k:v for k,v in native.items() if k!=CHANGED_NATIVE}
    return result

def compare(before_path=BEFORE,after_path=AFTER,output=OUTPUT):
    before=load(before_path);after=load(after_path)
    validate_report(before,'before');validate_report(after,'after')
    for key in UNCHANGED:
        assert before[key]==after[key],f'benchmark input changed: {key}'
    a=before['native_sources_sha256'];b=after['native_sources_sha256']
    assert set(a)==set(b), 'native source file-key sets differ'
    changed={key for key in a if a[key]!=b[key]}
    assert changed=={CHANGED_NATIVE},f'expected only {CHANGED_NATIVE} to differ, got {sorted(changed)}'
    profiles=[]
    for i,profile in enumerate(before['profiles'],1):
        modes={}
        for mode in ('setup','fresh','prepared'):
            key=f'{mode}_{i}';old=before['samples_ns'][key];new=after['samples_ns'][key]
            modes[mode]=dict(before_median_ns=before['median_ns'][key],after_median_ns=after['median_ns'][key],
                before_range_ns=[min(old),max(old)],after_range_ns=[min(new),max(new)],
                median_speedup=before['median_ns'][key]/after['median_ns'][key],
                all_after_samples_beat_all_before=max(new)<min(old))
        checksum={}
        for mode,iters in (('fresh',100),('prepared',1000)):
            old=before['checksums'][f'{mode}_{i}'][1:];new=after['checksums'][f'{mode}_{i}'][1:]
            assert old[0]==new[0],f'{mode}_{i}: root-count sum changed'
            for field,(x,y) in enumerate(zip(old[1:],new[1:]),1):
                tolerance=1e-9*iters+1e-12*max(abs(x),abs(y))
                assert abs(x-y)<=tolerance,(f'{mode}_{i} checksum field {field}',x,y,tolerance)
            checksum[mode]=dict(root_count_sum=old[0],iterations=iters,maximum_absolute_difference=max(abs(x-y) for x,y in zip(old[1:],new[1:])))
        profiles.append(dict(profile=profile,setup=modes['setup'],fresh=modes['fresh'],prepared=modes['prepared'],checksum_validation=checksum))
    before_norm=normalized(before);after_norm=normalized(after)
    assert before_norm==after_norm,'normalized benchmark input identities differ'
    def report_sha(path):return hashlib.sha256(Path(path).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    def relative(path):
        try:return Path(path).resolve().relative_to(ROOT.resolve()).as_posix()
        except ValueError:return Path(path).resolve().as_posix()
    result=dict(scope='Native-only forty source-backed root/Gibbs states across five isotherms; unchanged numerical inputs. These timing pressures each have one stable root, while their isotherms still require branch discovery. No Fortran or multi-root-query speedup claim.',
        generator=Path(__file__).resolve().relative_to(ROOT.resolve()).as_posix(),
        before_report=relative(before_path),after_report=relative(after_path),
        before_report_sha256=report_sha(before_path),after_report_sha256=report_sha(after_path),
        normalized_inputs_sha256=digest(before_norm),
        native_source_files_changed=sorted(changed),profiles=profiles)
    Path(output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--before',default=BEFORE);parser.add_argument('--after',default=AFTER);parser.add_argument('--output',default=OUTPUT)
    args=parser.parse_args();print(json.dumps(compare(args.before,args.after,args.output),indent=2))
