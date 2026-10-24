"""Audit stored MEOS flash timing evidence without a Fortran runtime."""
import argparse
import hashlib
import json
import math
import re
import statistics
from pathlib import Path

from benchmark_meos_flash import (
    ROOT, endpoint_coordinates, expected_checksum, valid_result, verify_checksum,
)


def audit(report):
    assert report['schema_version'] == 1
    assert re.fullmatch(r'[0-9a-f]{40,64}', report['git_head'])
    for name in ('native_sources_sha256', 'benchmark_sources_sha256',
                 'cargo_inputs_sha256', 'runtime_inputs_sha256'):
        assert report[name]
        assert all(re.fullmatch(r'[0-9a-f]{64}', value) for value in report[name].values())
    runtime = report['runtime_inputs_sha256']
    basename = lambda name: name.replace('\\', '/').rsplit('/', 1)[-1]
    library = {value for name,value in runtime.items() if basename(name) in ('thermopack.dll','libthermopack.so')}
    probe = {value for name,value in runtime.items() if basename(name) in ('reference_probe.dll','libreference_probe.so')}
    assert library == {report['reference_metadata']['library_sha256']}
    assert probe == {report['probe_metadata']['probe_sha256']}
    assert report['probe_metadata']['library_sha256'] in library
    rows = report['numeric_input_rows']; preflight = report['preflight']
    assert rows and preflight['rows'] == len(rows)
    manifest = str(len(rows))+'\n'+'\n'.join(' '.join(format(x,'.17g') for x in row) for row in rows)+'\n'
    assert hashlib.sha256(manifest.encode()).hexdigest() == report['numeric_input_sha256']
    groups = {}; endpoints = ordinary = 0
    outputs = preflight['outputs']
    source = outputs['source']; native = outputs['native']
    for row in rows:
        assert len(row) == 18 and all(map(math.isfinite,row))
        profile, lane, kind = map(int,row[:3])
        assert row[:3] == [profile,lane,kind] and 0 <= profile < len(report['profile_names'])
        assert lane in report['lane_names'].values()
        assert kind == (0 if lane == 0 else 1 if lane in (1,2) else 2 if lane in (3,4) else 5)
        group = f'{profile}_{lane}'; keys = groups.setdefault(group,[])
        key = f'check_{group}_{len(keys)}'; keys.append(key)
        context = preflight['request_contexts'][key]
        endpoint = context['seed_quality'] in (0.,1.) and lane in (1,2,3,4)
        if lane != 5:
            assert context['endpoint_seed'] == endpoint
            endpoints += endpoint
            ordinary += lane in (1,2,3,4) and not endpoint and row[14] == row[15]
        a, b = source[key][:], native[key][:]
        assert valid_result(a,lane == 5) and valid_result(b,lane == 5)
        if endpoint:
            a,b = endpoint_coordinates(a),endpoint_coordinates(b)
        elif lane != 5 and a[2] != 0. and b[2] != 0.:
            a[3] = b[3] = 0.
        floors = [2e-7,2e-5,0.,2e-7,2e-5,2e-7,2e-5,2e-7] if lane == 5 else [2e-7,2e-5,0.,2e-7,1e-12,1e-12,2e-5,2e-7]
        assert all(abs(x-y) <= floor+2e-6*max(abs(x),abs(y)) for x,y,floor in zip(a,b,floors))
    keys = {key for selected in groups.values() for key in selected}
    assert source.keys() == native.keys() == keys
    assert preflight['request_contexts'].keys() == keys
    assert {k.removeprefix('status_') for k in outputs['source_status']} == {k.removeprefix('check_') for k in keys}
    assert all(v in ([-1.],[0.]) for v in outputs['source_status'].values())
    assert preflight['group_counts'] == {key:len(value) for key,value in groups.items()}
    assert endpoints == preflight['audit']['endpoint_seed_requests']
    assert ordinary == preflight['audit']['ordinary_homogeneous_caloric_rows']
    checks = preflight['independent_kernel_checks']
    assert {(e['case'],e['implementation']) for e in checks} == {(key,label) for key in keys for label in ('source','native')}
    assert len(checks) == 2*len(rows)
    assert all(all(math.isfinite(v) and abs(v) <= 2e-8 for v in e['scaled_caloric_errors']) for e in checks)
    for state in preflight['original_kernel_states']:
        assert state['temperature'] > 0. and state['volume'] > 0.
        assert len(state['values']) == 7 and all(map(math.isfinite,state['values']))
        assert state['values'][2] < 0. and state['values'][6] > 0.
    samples = report['timed_samples_ns']; medians = report['median_ns']
    expected_names = {label+'_'+group for label in ('source','native') for group in groups}
    assert samples.keys() == medians.keys() == expected_names
    lengths = {len(values) for values in samples.values()}
    assert len(lengths) == 1 and min(lengths) >= 1
    count = lengths.pop(); raw = report['raw_checksums_by_sample']
    assert len(raw) == count+1 and report['iterations_per_sample'] == 1000
    for i, sample in enumerate(raw):
        assert sample['sample_index'] == i and sample['warmup'] == (i == 0)
        assert sample['process_order'] == (['native','source'] if i % 2 else ['source','native'])
        for label, values_by_case in (('source',source),('native',native)):
            assert sample[label].keys() == groups.keys()
            for group, selected in groups.items():
                values = sample[label][group]
                assert len(values) == 9 and all(map(math.isfinite,values)) and values[0] > 0.
                n = len(selected)
                magnitude = [math.fsum(abs(values_by_case[k][j])*(1000//n+(idx < 1000 % n))
                                      for idx,k in enumerate(selected)) for j in range(8)]
                verify_checksum(values[1:],expected_checksum(values_by_case,selected),magnitude)
                if i:
                    assert samples[label+'_'+group][i-1] == values[0]
    for name,values in samples.items():
        assert all(math.isfinite(v) and v > 0. for v in values)
        assert statistics.median(values) == medians[name]
    return dict(rows=len(rows),groups=len(groups),measured_samples=count,
                faster_median_groups=sum(medians['source_'+g] > medians['native_'+g] for g in groups))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reports',nargs='*',type=Path,default=[ROOT/'docs/rust/meos-flash-baseline.json'])
    for path in parser.parse_args().reports:
        result = audit(json.loads(path.read_text(encoding='utf-8')))
        print(path.name,json.dumps(result,sort_keys=True))
