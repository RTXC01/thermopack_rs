"""Offline integrity audit for total-MBWR TP parity evidence."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re

from reference import ROOT
from qualify_mbwr_total_tp import (RELATIVE, classical, differences, floors, key,
                                   map_tvp, schedule, singular_source_gap,
                                   unresolved_source_derivative_case, unqualified_entry)
from qualify_mbwr_caloric import R, assemble

COUNTS = {'TVP': 4959, 'INVALID': 1035, 'SOURCE': 3434, 'BOUNDED': 324}
EXPECTED_CASES = sum(COUNTS.values())
RAW_LENGTHS = {'CLASSICAL': 18, 'VALUES': 45, 'RESIDUAL': 10, 'IDEAL': 15,
               'TVP': 13, 'RAW_SOUND': 1, 'DENSITY': 4}
HASH_RE = re.compile(r'[a-f0-9]{64}')


def _check(condition, message):
    if not condition:
        raise AssertionError(message)


def _raw_float(x):
    """Decode the report's JSON-safe representation of an IEEE value."""
    if isinstance(x, str):
        _check(x in ('nan', 'inf', '-inf'), f'bad nonfinite encoding: {x!r}')
        return {'nan': math.nan, 'inf': math.inf, '-inf': -math.inf}[x]
    _check(type(x) in (int, float), f'raw value is not numeric: {x!r}')
    return float(x)


def _raw_vector(raw, label):
    _check(label in raw, f'missing original {label}')
    values = raw[label]
    _check(isinstance(values, list) and len(values) == RAW_LENGTHS[label],
           f'original {label} must contain {RAW_LENGTHS[label]} values')
    return [_raw_float(x) for x in values]


def _same_encoded(a, b):
    """Require storage to faithfully encode the source stdout values."""
    return (math.isnan(a) and math.isnan(b)) or a == b


def _same_sign(a, b):
    return a == 0.0 or b == 0.0 or (a > 0.0) == (b > 0.0)


def _within_4_ulp(actual, expected, scale):
    if not (math.isfinite(scale) and scale >= 0.0):
        return False
    bound = 4.0 * math.ulp(scale)
    return math.isfinite(bound) and _same_sign(actual, expected) and abs(actual - expected) <= bound


def _physical_scales(q, raw):
    if q['kind'] not in ('TVP', 'BOUNDED'):
        return None
    tv = _raw_vector(raw, 'VALUES')
    ideal = _raw_vector(raw, 'IDEAL')
    residual = _raw_vector(raw, 'RESIDUAL')
    f = tv[:10].copy()
    if ideal[3] == 10.:
        f[4] -= q['n'] * 999. * ideal[4] / (R * q['T'] ** 2)
    p1 = -R * (f[2] + q['T'] * f[5])
    p2 = -R * q['T'] * f[7]
    cv = -R * q['T'] * (2.0 * f[1] + q['T'] * f[4])
    if cv != 0.0:
        denom = abs(p2 - q['T'] * (p1 * p1) / cv)
        sound_factor = max(1.0, abs(p2) / denom) if denom > 0.0 else 1.0
    else:
        sound_factor = 1.0
    return (residual[3], sound_factor)


def _same_expected(stored, want, scales=None):
    """Require stored expected views to match derived values, allowing operand-scaled 4 ULP only on PHYSICAL lnphi and present Cp/sound."""
    if not isinstance(stored, dict) or set(stored) != set(want):
        return False
    residual_fn, sound_factor = scales if scales is not None else (0.0, 1.0)
    for label, expected in want.items():
        actual = stored[label]
        if isinstance(expected, str):
            if actual != expected:
                return False
        else:
            if not isinstance(actual, list) or len(actual) != len(expected):
                return False
            for j, (a, b) in enumerate(zip(actual, expected)):
                if type(a) not in (int, float) or type(b) not in (int, float):
                    return False
                if not (math.isfinite(a) and math.isfinite(b)):
                    return False
                if a == b:
                    continue
                if label == 'PHYSICAL' and j == 41:
                    if not _within_4_ulp(a, b, max(abs(b), abs(residual_fn))):
                        return False
                elif label == 'PHYSICAL' and j == 48 and expected[47] == 1.0 and actual[47] == 1.0:
                    if not _within_4_ulp(a, b, max(abs(b), abs(expected[46]))):
                        return False
                elif label == 'PHYSICAL' and j == 50 and expected[49] == 1.0 and actual[49] == 1.0:
                    if not _within_4_ulp(a, b, abs(b) * sound_factor):
                        return False
                else:
                    return False
    return True


def _stdout_vectors(source):
    _check(not source.get('censored'), 'source collection was censored')
    _check(source.get('returncode') == 0, 'source collection failed')
    vectors = {}
    for m in re.finditer(r'^CASE=(\d+) (CLASSICAL|VALUES|RESIDUAL|IDEAL|TVP|RAW_SOUND|DENSITY)=([^\r\n]+)$', source.get('stdout', ''), re.M):
        i, label, body = int(m[1]), m[2], m[3]
        k = (i, label)
        _check(k not in vectors, f'duplicate raw stdout vector {k}')
        tokens = body.split()
        _check(len(tokens) == RAW_LENGTHS[label], f'raw stdout {k} has wrong length')
        try:
            vectors[k] = [float(x.replace('D', 'E')) for x in tokens]
        except ValueError as exc:
            raise AssertionError(f'invalid raw stdout number in {k}') from exc
    return vectors


def _expected(q, raw):
    c = _raw_vector(raw, 'CLASSICAL')
    out = {'CLASSICAL': classical(c)}
    if q['kind'] == 'SOURCE':
        d = _raw_vector(raw, 'DENSITY')
        out['DENSITY'] = d if all(map(math.isfinite, d)) else 'ERROR: non-finite original density'
    elif q['kind'] == 'INVALID':
        out.update(BASELINE='ERROR: nonpositive EOS pressure', PHYSICAL='ERROR: nonpositive EOS pressure')
    else:
        tv = _raw_vector(raw, 'VALUES')
        ideal = _raw_vector(raw, 'IDEAL')
        residual = _raw_vector(raw, 'RESIDUAL')
        base = map_tvp(q, tv, residual)
        out['BASELINE'] = base
        f = tv[:10].copy()
        if ideal[3] == 10.:
            f[4] -= q['n'] * 999. * ideal[4] / (R * q['T'] ** 2)
        physical = assemble(q, f, residual[3], ideal[2] * 1e-3)
        out['PHYSICAL'] = map_tvp(q, physical, residual)
        direct = _raw_vector(raw, 'TVP')
        probe = [base[0]]
        for start in (11, 16, 41):
            probe += [base[start]] + base[start + 2:start + 5]
        _check(len(probe) == 13, 'internal direct TVP probe shape')
        direct_floors = [floors(q)[0]]
        for start in (11, 16, 41):
            direct_floors += [floors(q)[start]] + floors(q)[start + 2:start + 5]
        for j, (want, got, floor) in enumerate(zip(probe, direct, direct_floors)):
            _check(math.isfinite(got) and abs(want-got) <= floor + RELATIVE*max(abs(want), abs(got)),
                   f'independent original TVP mismatch at {j}')
        if q['kind'] == 'BOUNDED':
            _check(abs(tv[10] - q['P']) <= 1e-6 + 2e-9*q['P'], 'original EOS pressure does not close')
            _check(tv[12] <= 0., 'original EOS state is mechanically unstable')
    return out


def _compare(q, expected, actual):
    _check(set(actual) == set(expected), f'{q["kind"]} output views differ from required views')
    for label, want in expected.items():
        got = actual[label]
        if isinstance(want, str):
            _check(isinstance(got, str) and got.startswith('ERROR'), f'{label}: expected retained source error')
        elif label == 'DENSITY':
            _check(isinstance(got, list) and len(got) == 4, 'DENSITY native vector shape')
            _check(isinstance(want, list) and len(want) == 4, 'DENSITY original vector shape')
            _check(got[1] == want[1], 'DENSITY phase flag drift')
            floors_ = (1e-12, 0., 1e-5, 1e-8)
            for j, (a, b, floor) in enumerate(zip(got, want, floors_)):
                _check(math.isfinite(a) and math.isfinite(b), f'DENSITY nonfinite at {j}')
                _check(abs(a-b) <= floor + (0. if j == 1 else 2e-9)*max(abs(a), abs(b)), f'DENSITY mismatch at {j}')
        else:
            bad = differences(q, want, got, label == 'CLASSICAL')
            _check(not bad, f'{label} mismatch: {bad[:3]}')


def _verify_unqualified_entries(claimed, computed):
    _check(len(computed) == 93, f'expected 93 raw-conditioned source cases, found {len(computed)}')
    _check(claimed == computed, 'unqualified source derivative list differs from original raw evidence')


def audit(report):
    records = report['records']
    requested = schedule()
    _check(report.get('cases') == EXPECTED_CASES and len(records) == EXPECTED_CASES, 'case count is not 9752')
    sig = report['signature']; recorded_schedule = sig['schedule']
    _check(len(recorded_schedule) == EXPECTED_CASES, 'signature schedule length mismatch')
    for i, (authoritative, recorded) in enumerate(zip(requested, recorded_schedule), 1):
        if authoritative['kind'] == 'BOUNDED':
            _check(set(authoritative)-{'V'} == set(recorded)-{'V'} and
                   all(authoritative[k] == recorded[k] for k in authoritative if k != 'V'),
                   f'bounded request {i} differs from authoritative schedule')
            _check(type(recorded.get('V')) in (int, float) and math.isfinite(recorded['V']) and recorded['V'] > 0.,
                   f'bounded proposal {i} has invalid volume')
        else:
            _check(recorded == authoritative, f'request {i} differs from authoritative schedule')
    _check([e['request'] for e in records] == recorded_schedule, 'record requests differ from signed schedule')
    _check(len({key(q) for q in requested}) == EXPECTED_CASES, 'duplicate request in cohort')
    counts = {k: sum(q['kind'] == k for q in requested) for k in COUNTS}
    _check(counts == COUNTS and report.get('counts') == COUNTS, 'cohort kind counts mismatch')
    _check(report.get('relative_tolerance') == RELATIVE, 'relative tolerance mismatch')
    _check(not report.get('errors'), 'qualifier report contains errors')
    source = report['source_process']
    _check(not source.get('censored') and source.get('returncode') == 0, 'source collection failed or censored')
    hashes = report['native_sources_sha256'] | sig['object_hashes'] | sig.get('support_sha256', {})
    _check(bool(hashes) and all(isinstance(v, str) and HASH_RE.fullmatch(v) for v in hashes.values()), 'invalid source/object hash')
    _check(HASH_RE.fullmatch(report['native_executable_sha256']) is not None, 'invalid executable hash')
    raw_stdout = _stdout_vectors(source)
    declared_gaps = []
    unqualified = []
    for i, e in enumerate(records, 1):
        q = e['request']; kind = q['kind']; raw = e['original']
        needed = {'CLASSICAL', 'DENSITY'} if kind == 'SOURCE' else {'CLASSICAL'}
        if kind in ('TVP', 'BOUNDED'):
            needed |= {'VALUES', 'RESIDUAL', 'IDEAL', 'TVP', 'RAW_SOUND'}
        _check(set(raw) == needed, f'case {i}: original raw labels differ from required evidence')
        for label in needed:
            stored = _raw_vector(raw, label)
            parsed = raw_stdout.get((i, label))
            _check(parsed is not None and len(parsed) == len(stored), f'case {i}: missing raw stdout {label}')
            _check(all(_same_encoded(a, b) for a, b in zip(stored, parsed)), f'case {i}: stored {label} differs from raw stdout')
        want = _expected(q, raw)
        scales = _physical_scales(q, raw)
        _check(_same_expected(e.get('expected'), want, scales), f'case {i}: stored expected values are not derived from original evidence')
        if unresolved_source_derivative_case(q, raw):
            unqualified.append(unqualified_entry(q))
        native_views = ({'DENSITY', 'CLASSICAL'} if kind == 'SOURCE' else
                        {'PHYSICAL'} if kind == 'BOUNDED' else
                        {'BASELINE', 'PHYSICAL', 'CLASSICAL'})
        _check(set(e['native']) == native_views, f'case {i}: native output views differ from required views')
        # Actual differences inside the source-only conditioning domain are
        # recorded explicitly; the original gold remains unchanged.
        for label in native_views:
            expected_values = want[label]
            try:
                _compare(q, {label: expected_values}, {label: e['native'][label]})
            except AssertionError:
                if label != 'CLASSICAL' or not singular_source_gap(q, raw, expected_values, e['native'][label]):
                    raise
                declared_gaps.append(unqualified_entry(q))
    _check(len(raw_stdout) == sum(len({'CLASSICAL','DENSITY'} if q['kind']=='SOURCE' else {'CLASSICAL'} |
                                             ({'VALUES','RESIDUAL','IDEAL','TVP','RAW_SOUND'} if q['kind'] in ('TVP','BOUNDED') else set()))
                                     for q in requested), 'unexpected or missing raw stdout rows')
    _verify_unqualified_entries(report.get('unqualified_source_derivatives'), unqualified)
    _check(report.get('known_mismatches') == declared_gaps,
           'declared known mismatches differ from independently classified native differences')
    return {'cases': len(records), 'counts': counts, 'raw_vectors': len(raw_stdout),
            'known_mismatches': len(declared_gaps), 'unqualified_source_derivatives': len(unqualified),
            'qualified_cases': len(records) - len(unqualified)}


def load_and_audit():
    fixture = ROOT/'tests/fixtures/mbwr_total_tp.csv'
    metadata_path = ROOT/'tests/fixtures/mbwr_total_tp.metadata.json'
    metadata = json.loads(metadata_path.read_text(encoding='utf-8'))
    archive = fixture.with_suffix('.audit.json.gz')
    _check(metadata.get('evidence_file') == archive.name, 'metadata evidence archive name mismatch')
    csv_bytes = fixture.read_bytes(); archive_bytes = archive.read_bytes()
    _check(hashlib.sha256(csv_bytes).hexdigest() == metadata.get('fixture_sha256'), 'fixture SHA-256 mismatch')
    _check(hashlib.sha256(archive_bytes).hexdigest() == metadata.get('evidence_sha256'), 'archive SHA-256 mismatch')
    report = json.loads(gzip.decompress(archive_bytes))
    _check(report.get('fixture_sha256') == metadata['fixture_sha256'], 'report fixture SHA-256 mismatch')
    result = audit(report)
    for field in ('cases','counts','errors','known_mismatches','unqualified_source_derivatives','native_sources_sha256','native_executable_sha256','relative_tolerance'):
        _check(report[field] == metadata[field], f'metadata {field} mismatch')
    _check({k:v for k,v in report['signature'].items() if k != 'schedule'} == metadata['signature'], 'metadata signature mismatch')
    with fixture.open(encoding='utf-8', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter=';'))
    _check(len(rows) == EXPECTED_CASES, 'fixture CSV row count mismatch')
    for row, e in zip(rows, report['records']):
        _check(row['case'] == key(e['request']), 'fixture CSV case key mismatch')
        for column, label in (('baseline','BASELINE'),('physical','PHYSICAL'),('classical','CLASSICAL'),('density','DENSITY')):
            value = e['expected'].get(label, '')
            cell = value if isinstance(value, str) else ','.join(map(repr, value))
            _check(row[column] == cell, f'fixture CSV {column} is not bound to evidence')
    result['fixture_sha256'] = metadata['fixture_sha256']
    return report, result


if __name__ == '__main__':
    _, result = load_and_audit()
    print(json.dumps(result, sort_keys=True))
