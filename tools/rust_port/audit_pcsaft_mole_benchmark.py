"""Audit the captured PC-SAFT mixture F_n timing record without executing it."""
import argparse
import base64
import gzip
import hashlib
import json
import math
import re
import statistics
from pathlib import Path

import benchmark_pcsaft_mole_derivatives as bench

ROOT = bench.ROOT
GROUPS = [('prepared', 'baseline'), ('prepared', 'candidate'), ('prepared', 'original'), ('fresh', 'baseline'), ('fresh', 'candidate')]
SNAPSHOT_FILES = {'benchmark_pcsaft_mole_derivatives.py', 'pcsaft_mole_bench.rs',
                  'pcsaft_mole_bench.f90', 'pcsaft_mole_bench_worker.py'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def raw_result(record, count):
    require(isinstance(record, dict) and record.get('returncode') == 0 and 'failure' not in record, 'failed process')
    lines = [x[7:] for x in record.get('stdout', '').splitlines() if x.startswith('RESULT=')]
    require(len(lines) == 1, 'missing or duplicate RESULT')
    result = json.loads(lines[0], parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
    require(result == record.get('result') and isinstance(result, dict), 'raw and parsed result differ')
    require(bench.valid_vector(result.get('vector'), count), 'invalid vector')
    return result


def audit_snapshot(record_path, data):
    """Bind the adjacent source archive to the captured record and manifest."""
    record_path = Path(record_path)
    metadata_path = record_path.with_name(record_path.stem + '.sources.metadata.json')
    metadata = json.loads(metadata_path.read_text(encoding='utf-8'))
    require(metadata['record'] == record_path.name, 'snapshot record name')
    record_bytes = record_path.read_bytes().replace(b'\r\n', b'\n')
    require(hashlib.sha256(record_bytes).hexdigest() == metadata['record_normalized_sha256'], 'snapshot record hash')
    archive_name = record_path.stem + '.sources.json.gz'
    require(metadata['archive'] == archive_name, 'snapshot archive name')
    archive = record_path.with_name(archive_name).read_bytes()
    require(hashlib.sha256(archive).hexdigest() == metadata['archive_sha256'], 'snapshot archive hash')
    require(metadata['encoding'] == 'base64 of exact measured input bytes', 'snapshot encoding')
    files = json.loads(gzip.decompress(archive))
    require(set(files) == SNAPSHOT_FILES and set(metadata['files']) == SNAPSHOT_FILES, 'snapshot files')
    hashes = data['manifest']['hashes']
    for name in SNAPSHOT_FILES:
        contents = base64.b64decode(files[name], validate=True)
        digest = hashlib.sha256(contents).hexdigest()
        matching = [v for k, v in hashes.items() if k.replace('\\', '/').split('/')[-1] == name]
        require(len(matching) == 1 and digest == metadata['files'][name] == matching[0], 'snapshot source hash: ' + name)


def audit(data, require_current_inputs=False):
    require(data.get('schema') == 1, 'record schema')
    m = data['manifest']
    require(m.get('schema') == 1, 'manifest schema')
    serialized = json.dumps(m, indent=2, allow_nan=False) + '\n'
    # The manifest was written in text mode on Windows (CRLF); accept either
    # platform's exact text-mode encoding, never a canonicalized JSON digest.
    hashes = {hashlib.sha256(serialized.replace('\n', nl).encode('utf-8')).hexdigest() for nl in ('\n', '\r\n')}
    require(data['manifest_sha256'] in hashes, 'manifest SHA mismatch')
    h = m['hashes']
    require(isinstance(h, dict) and bool(h) and all(isinstance(k, str) and isinstance(v, str) and re.fullmatch('[0-9a-f]{64}', v) for k, v in h.items()), 'invalid hashes')
    require(data['hashes_before'] == h == data['hashes_after'] and data['hashes_unchanged'] is True, 'input hashes changed')
    meta = m['source_metadata']
    require(meta['baseline'] == 'd68c794' and meta['library_sha256'] == h[m['source_dll']], 'source DLL identity')
    require(meta['flags'] == ['-O3', '-cpp', '-fdefault-real-8', '-fdefault-double-8', '-frecursive', '-fopenmp', '-fallow-argument-mismatch', '-ffree-line-length-none'], 'source flags')
    require(meta['compiler'] in m['versions']['gfortran'] and m['versions']['cargo'].startswith('cargo ') and m['versions']['rustc'].startswith('rustc ') and 'commit-hash:' in m['versions']['rustc'], 'compiler identity')
    require(set(m['git']) == {m['baseline'], m['candidate']}, 'git roots')
    for state in m['git'].values():
        require(re.fullmatch('[0-9a-f]{40}', state['head']) is not None and all(isinstance(state[k], str) for k in ('branch', 'status', 'tracked_rust_diff')), 'git state')
    require(m['git'][m['baseline']]['head'] == m['git'][m['candidate']]['head'] and m['git'][m['baseline']]['head'].startswith('9bd3c68'), 'git revision')
    require(not m['git'][m['baseline']]['tracked_rust_diff'] and bool(m['git'][m['candidate']]['tracked_rust_diff']), 'baseline/candidate state')
    cmds = m['commands']
    require(len(cmds) == 5, 'build command count')
    for i, label in ((0, 'baseline'), (2, 'candidate')):
        require(cmds[i]['cwd'] == m[label] and cmds[i]['argv'][1:] == ['build', '--release', '--lib', '--message-format=json', '-j1'], 'cargo command')
        c = cmds[i+1]['argv']; b = m['builds'][label]
        require(cmds[i+1]['cwd'] == m[label] and all(x in c for x in ('--edition=2024', '-O', '-C', 'lto=thin', '--extern', 'thermopack_rs='+b['rlib'], '-L', '-o', b['exe'])) and any(x.endswith('pcsaft_mole_bench.rs') for x in c), 'Rust build command')
        require(all(b[x] in h for x in ('exe', 'rlib')), 'missing build hash')
    c = cmds[4]['argv']
    require(cmds[4]['cwd'] and all(x in c for x in ('-O3', '-fdefault-real-8', '-fdefault-double-8', '-frecursive', '-shared', m['source_dll'], '-o', m['helper'])) and any(x.endswith('pcsaft_mole_bench.f90') for x in c), 'Fortran command')
    require(m['helper'] in h, 'helper hash')
    for suffix in ('benchmark_pcsaft_mole_derivatives.py', 'pcsaft_mole_bench.rs', 'pcsaft_mole_bench.f90', 'pcsaft_mole_bench_worker.py', 'pc_saft_nonassoc.f90', 'reference.py', 'metadata.json', 'thermopack_var.mod', 'pc_saft_nonassoc.mod'):
        require(any(k.lower().endswith(suffix.lower()) for k in h), 'missing input hash: '+suffix)
    if require_current_inputs:
        bench.validate_hashes(h)
        for root, state in m['git'].items():
            require(bench.git_state(root) == state, 'current git state differs')
    s = data['settings']
    require(s['preflight_only'] is False and type(s['samples']) is int and s['samples'] >= 7 and type(s['iterations']) is int and s['iterations'] > 0 and type(s['warmup']) is int and s['warmup'] > 0, 'invalid measurement settings')
    rows = data['cases']; expected = list(bench.cases())
    require(len(rows) == len(expected) == 32 and [r['case'] for r in rows] == expected, 'case order or scientific scope')
    ratios = {key: [] for key in ('baseline_over_candidate_prepared', 'original_over_candidate_prepared', 'baseline_over_candidate_fresh')}
    for row in rows:
        nc = len(row['case']['moles']); pre = row['preflight']
        require(set(pre) == {'original', 'baseline', 'candidate'}, 'preflight implementations')
        vectors = {name: raw_result(pre[name], nc)['vector'] for name in pre}
        for name, tol in (('baseline', 2e-6), ('candidate', 2e-11)):
            require(pre[name].get('within_tolerance') is True and all(abs(a-b) <= tol*(1+abs(a)) for a,b in zip(vectors['original'], vectors[name])), 'preflight tolerance')
        for name in pre:
            require(bench.sample_error(pre[name], vectors['original'], vectors[name], nc, 2e-6 if name == 'baseline' else 2e-11, 1) is None, 'invalid preflight result')
        order = []
        for sample in range(s['samples']):
            names = ['baseline', 'candidate', 'original']
            if (sample // 3) % 2: names.reverse()
            shift = sample % 3; names = names[shift:] + names[:shift]
            order.extend((sample, mode, name) for mode in ('prepared', 'fresh') for name in names if mode != 'fresh' or name != 'original')
        samples = row['samples']
        require([(x.get('sample'), x.get('mode'), x.get('implementation')) for x in samples] == order, 'missing, duplicate, or out-of-order samples')
        grouped = {g: [] for g in GROUPS}
        for x in samples:
            name = x['implementation']; result = raw_result(x, nc)
            require(x.get('validation_error', 'missing') is None and bench.sample_error(x, vectors['original'], vectors[name], nc, 2e-6 if name == 'baseline' else 2e-11, s['iterations']) is None, 'invalid timed sample')
            grouped[(x['mode'], name)].append(result['seconds'] / s['iterations'])
        require(row.get('timing_valid') is True, 'timing invalid')
        summaries = row['summaries_seconds_per_call']
        require(set(summaries) == {a+'_'+b for a,b in GROUPS}, 'summary groups')
        medians = {}
        for (mode,name), times in grouped.items():
            require(len(times) == s['samples'], 'incomplete group')
            values = dict(count=len(times), median=statistics.median(times), minimum=min(times), maximum=max(times), stdev=statistics.stdev(times))
            actual = summaries[mode+'_'+name]
            require(set(actual) == set(values) and all(type(actual[k]) in (int,float) and math.isfinite(actual[k]) and math.isclose(actual[k], v, rel_tol=1e-12, abs_tol=1e-20) for k,v in values.items()), 'summary mismatch')
            medians[mode+'_'+name] = values['median']
        computed = dict(baseline_over_candidate_prepared=medians['prepared_baseline']/medians['prepared_candidate'], original_over_candidate_prepared=medians['prepared_original']/medians['prepared_candidate'], baseline_over_candidate_fresh=medians['fresh_baseline']/medians['fresh_candidate'])
        require(set(row['median_ratios']) == set(computed), 'ratio keys')
        for k,v in computed.items():
            require(type(row['median_ratios'][k]) in (int,float) and math.isclose(row['median_ratios'][k], v, rel_tol=1e-12), 'ratio mismatch')
            ratios[k].append(v)
    return dict(profiles=len(rows), preflights=3*len(rows), timed_calls=len(rows)*s['samples']*5, groups=len(rows)*5, ratios=ratios)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('record', nargs='?', type=Path, default=ROOT/'docs/rust/pcsaft-mole-performance.json')
    p.add_argument('--require-current-inputs', action='store_true')
    args = p.parse_args()
    data = json.loads(args.record.read_text(encoding='utf-8'))
    audit_snapshot(args.record, data)
    report = audit(data, args.require_current_inputs)
    print(f"Audited {report['profiles']} profiles, {report['preflights']} preflights, {report['timed_calls']} timed calls, {report['groups']} complete groups.")
    for key, values in report['ratios'].items():
        print(f'{key}: min={min(values):.4g}, median={statistics.median(values):.4g}, max={max(values):.4g}')


if __name__ == '__main__':
    main()
