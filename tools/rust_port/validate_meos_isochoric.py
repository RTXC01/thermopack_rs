"""Pure MEOS UV/SV audit with original-kernel targets and retained failures.

Every positive case runs natively, including source failures. Partial source
qualification needs an exact declared outcome ledger; it never accepts native
errors or physical coordinate drift.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
from pathlib import Path

from reference import ROOT
from deadline_reference import DeadlineReferenceProcess
from validate_meos_flash import atomic_csv, atomic_text, digest

FIELDS = ['case', 'T', 'P', 'two_phase', 'beta_v', 'beta_l', 'V', 'U', 'S', 'raw_phase']


def finite(x):
    return isinstance(x, (int, float)) and math.isfinite(x)


def model_data():
    with (ROOT/'tests/fixtures/meos_metadata.csv').open(encoding='utf-8') as f:
        return {r['case'].split('|')[1]: {k: float(v) for k, v in r.items() if k != 'case'}
                for r in csv.DictReader(f, delimiter=';')}


def seed_schedule(data):
    for model, m in data.items():
        t = 1.2*m['Tc'] if model.startswith('LJTS_REDUCED') else max(1.2*m['Tc'], 90.)
        yield dict(model=model, t=t, p=None, volumes=[3./m['rhoc']], quality=None, tag='homogeneous')
    for name in ('meos_saturation.csv', 'meos_saturation_recovered.csv'):
        with (ROOT/'tests/fixtures'/name).open(encoding='utf-8') as f:
            for row in list(csv.reader(f, delimiter=';'))[1:]:
                case = row[0].split('|')
                if case[5] != 'SAT_P':
                    continue
                t, p, vl, vv = map(float, row[1:5])
                tr = t/data[case[1]]['Tc']
                if not any(abs(tr-x) < 1e-10 for x in (.5, .9, .999)):
                    continue
                for beta in (0., .5, 1.):
                    yield dict(model=case[1], t=t, p=p, volumes=[vl, vv], quality=beta,
                               tag=f'saturation:{tr:.3f}:{beta}')


def kernel(oracle, model, t, v):
    eos, component = model.split(':')
    result = oracle.call(dict(mode='MEOS', eos=eos, component=component,
                              temperature=t, volume=v), timeout=4.)
    assert 'values' in result, (model, t, v, result)
    return result['values']


def equilibrium_kernel(data, model, t, p, volumes, beta, states):
    r = data[model]['R']
    assert len(volumes) == len(states) in (1, 2)
    assert all(all(finite(s[i]) for i in (6, 12, 13, 14, 15, 16, 17, 20))
               and s[14] < 0. and s[20] > 0. for s in states)
    for v, s in zip(volumes, states):
        assert abs(s[12]-p) <= 2e-8*p+512.*math.ulp(1.)*r*t/v, ('pressure', s[12], p)
    weights = [1.] if len(states) == 1 else [1.-beta, beta]
    volume = sum(w*v for w, v in zip(weights, volumes))
    u = sum(w*s[15] for w, s in zip(weights, states))
    entropy = sum(w*s[16] for w, s in zip(weights, states))
    gap = None
    if len(states) == 1:
        cv = states[0][20]
    else:
        vl, vv = volumes
        assert (vv-vl)/(vv+vl) >= 1e-6
        gap = math.log(vv/vl)+states[0][6]-states[1][6]+p*(vl-vv)/(r*t)
        assert abs(gap) < 2e-8, ('chemical potential', gap)
        slope = (states[1][17]-states[0][17])/(t*(vv-vl))
        cv = sum(w*(s[20]-t*(s[13]-slope)**2/s[14]) for w, s in zip(weights, states))
    assert finite(cv) and cv > 0.
    return dict(volume=volume, u=u, s=entropy, cv=cv, chemical_potential_gap=gap)


def prepare(oracle, data):
    result = []
    cache = {}
    for seed in seed_schedule(data):
        model, t = seed['model'], seed['t']
        states = []
        for v in seed['volumes']:
            k = model, t, v
            if k not in cache:
                cache[k] = kernel(oracle, model, t, v)
            states.append(cache[k])
        p = states[0][12] if seed['p'] is None else seed['p']
        seed = dict(seed, p=p, kernels=states,
                    expected=equilibrium_kernel(data, model, t, p, seed['volumes'], seed['quality'], states))
        result.append(seed)
    assert len(result) == 140, len(result)  # fourteen homogeneous + fourteen*3*3 saturation seeds
    assert len({(s['model'], s['tag']) for s in result}) == len(result)
    return result


def requests(seeds, data):
    for seed in seeds:
        eos, component = seed['model'].split(':')
        reduced = eos == 'LJTS_REDUCED'
        tc = data[seed['model']]['Tc']
        for kind in ('UV', 'SV'):
            for guessed in (False, True):
                q = dict(mode='MEOS_ISO', eos=eos, component=component, kind=kind,
                         temperature=seed['t'], volume=seed['expected']['volume'],
                         target=seed['expected']['u' if kind == 'UV' else 's'],
                         minimum_temperature=.3*tc if reduced else (80. if seed['quality'] is None else 2.),
                         maximum_temperature=3.*tc if reduced else 999.,
                         minimum_pressure=1e-8 if reduced else 10., maximum_pressure=70. if reduced else 1e8,
                         temperature_guess=.93*seed['t'] if guessed else None,
                         pressure_guess=1.2*seed['p'] if guessed else None,
                         seed_phase=4 if seed['quality'] is None else 0,
                         seed_beta=0. if seed['quality'] is None else seed['quality'],
                         seed_tag=seed['tag'])
                yield q, seed


def encode(q):
    control = ':'.join(map(str, ['ISO', q['kind'], q['target'],
        'none' if q['temperature_guess'] is None else q['temperature_guess'],
        'none' if q['pressure_guess'] is None else q['pressure_guess'],
        q['minimum_temperature'], q['maximum_temperature'], q['minimum_pressure'], q['maximum_pressure']]))
    return f"MEOS|{q['eos']}:{q['component']}|{q['temperature']!r}|{q['volume']!r}|1|{control}"


def outcome_key(q):
    return '|'.join(map(str, [q['eos'], q['component'], q['kind'], q['seed_tag'],
        q['temperature_guess'] is not None, q['minimum_temperature'], q['maximum_temperature'],
        q['minimum_pressure'], q['maximum_pressure']]))


def source_outcome(r):
    if 'reference_timeout' in r:
        return 'censored_timeout'
    if 'reference_error' in r:
        err = r['reference_error']
        if isinstance(err, dict) and isinstance(err.get('output'), list):
            output = err['output']
            stage = 'phase_properties' if any('THERMOPACK_STAGE=returned_phase_properties' in str(line) for line in output) else 'flash'
            filtered = [str(line) for line in output if not re.match(r'^\s*(?:#\d+\s+|at\s+)', str(line))]
            text = ('\n'.join(filtered) if filtered else str(r)).lower()
        else:
            stage = 'phase_properties' if 'THERMOPACK_STAGE=returned_phase_properties' in str(r) else 'flash'
            text = str(r).lower()
        for name, token in [('above_critical_saturation', 't > tc'), ('density_iteration', 'densitysolver'),
                            ('sv_nonconvergence', 'svsinglecomp'), ('uv_nonconvergence', 'uvsinglecomp'),
                            ('saturation_nonconvergence', 'saturation')]:
            if token in text:
                return stage+':'+name
        return stage+':unclassified_termination'
    if 'calculation_error' in r:
        return r.get('stage', 'unknown')+':calculation_error'
    return 'returned'


def source_closed(q, r, data):
    if source_outcome(r) != 'returned' or not isinstance(r.get('values'), list):
        return False
    a = r['values']; raw = r.get('raw_result', [])
    if len(a) != 8 or len(raw) != 7 or not all(map(finite, a+raw)):
        return False
    t, p, two, bv, bl, v, u, s = a
    if two not in (0., 1.) or min(t, p, v) <= 0. or raw[2] not in (0, 1, 2, 4) or (raw[2] == 0) != bool(two):
        return False
    raw_bv, raw_bl = raw[3:5]
    fractions_closed = 0. <= raw_bv <= 1. and 0. <= raw_bl <= 1. and abs(raw_bv+raw_bl-1.) <= 1e-10
    if raw[2] == 0 and (not fractions_closed or a[3:5] != raw[3:5]):
        return False
    if raw[2] == 1 and (raw_bv, raw_bl) != (0., 1.):
        return False
    if raw[2] == 2 and (raw_bv, raw_bl) != (1., 0.):
        return False
    if raw[2] == 4 and not (fractions_closed or (raw_bv, raw_bl) == (-1., -1.)):
        return False
    if not two and a[3:5] != [0., 1.]:
        return False
    if q['kind'] == 'SV' and r.get('status') != 0:
        return False
    if two and (min(bv, bl) < 0. or max(bv, bl) > 1. or abs(bv+bl-1.) > 1e-10 or raw[5:] != [1., 1.]):
        return False
    gas = data[q['eos']+':'+q['component']]['R']
    residual = (u if q['kind'] == 'UV' else s)-q['target']
    scale = gas*(t if q['kind'] == 'UV' else 1.)
    return (abs(residual) <= 2e-7*scale and abs(v/q['volume']-1.) <= 2e-7
        and q['minimum_temperature'] <= t <= q['maximum_temperature']
        and q['minimum_pressure'] <= p <= q['maximum_pressure']
        and len(r.get('pressure_dv', [])) == (2 if two else 1)
        and all(finite(x) and x < 0. for x in r['pressure_dv'])
        and len(r.get('pressure_errors', [])) == (2 if two else 1)
        and all(finite(e) and abs(e) <= 2e-7*p for e in r['pressure_errors']))


def valid_native_point(a):
    if len(a) != 14 or not all(map(finite, a)):
        return False
    t, p, two, bv, bl, v, u, s, phase, vl, vv, cv, raw_bv, raw_bl = a
    if two not in (0., 1.) or phase not in (0., 1., 2., 4.) or (phase == 0.) != bool(two):
        return False
    if min(t, p, v, vl, vv, cv) <= 0. or not (0. <= bv <= 1. and 0. <= bl <= 1.) or abs(bv+bl-1.) > 1e-12:
        return False
    if two:
        return (raw_bv == bv and raw_bl == bl and vv > vl
                and (vv-vl)/(vv+vl) >= 1e-6)
    expected = {1.:(0., 1.), 2.:(1., 0.), 4.:(-1., -1.)}[phase]
    return (bv, bl) == (0., 1.) and (raw_bv, raw_bl) == expected and vl == vv


def outcome_allowed(value, expected):
    if isinstance(expected, list):
        return value in expected
    return expected == value


def allowed_outcomes(actual, declared):
    return all(outcome_allowed(value, declared.get(key)) for key, value in actual.items())


def coordinate_bad(a, b, reduced):
    floors = [2e-8 if reduced else 2e-6, 1e-10 if reduced else 1e-3,
              0., 3e-6, 3e-6, 1e-8 if reduced else 1e-12,
              1e-8 if reduced else 1e-4, 1e-8 if reduced else 1e-6, 0.]
    assert len(a) == len(b) == 9, 'complete coordinate arrays required'
    return [[i, x, y] for i, (x, y) in enumerate(zip(a, b))
            if not finite(x) or not finite(y) or abs(x-y) > floors[i]+(0. if i in (2, 8) else 3e-7)*max(abs(x), abs(y))]


def bounded_endpoint_gap(wanted, got, reduced):
    if len(wanted) != 9 or len(got) != 9 or not all(map(finite, wanted+got)):
        return False
    for row in (wanted, got):
        if row[2] not in (0., 1.) or row[8] not in (0., 1., 2., 4.) or (row[8] == 0.) != bool(row[2]):
            return False
        if not (0. <= row[3] <= 1. and 0. <= row[4] <= 1.) or abs(row[3]+row[4]-1.) > 1e-10:
            return False
    endpoint = ((wanted[2] == 1. and min(wanted[3], wanted[4]) <= 2e-6)
                or (got[2] == 1. and min(got[3], got[4]) <= 2e-6))
    return (endpoint and wanted[2] != got[2]
            and not any(b[0] in (0, 1, 5, 6, 7) for b in coordinate_bad(wanted, got, reduced)))


def main(args):
    data = model_data()
    inputs = ['meos_metadata.csv', 'meos_saturation.csv', 'meos_saturation_recovered.csv']
    signature = dict(reference=json.loads((ROOT/'.reference/metadata.json').read_text(encoding='utf-8')),
                     probe=json.loads((ROOT/'.reference/probe.metadata.json').read_text(encoding='utf-8')),
                     generator_sha256=digest(Path(__file__)), worker_sha256=digest(ROOT/'tools/rust_port/oracle_worker.py'),
                     inputs_sha256={n:digest(ROOT/'tests/fixtures'/n) for n in inputs})
    assert hashlib.sha256((ROOT/'.reference'/('thermopack.dll' if os.name == 'nt' else 'libthermopack.so')).read_bytes()).hexdigest() == signature['reference']['library_sha256']
    assert hashlib.sha256((ROOT/'.reference'/('reference_probe.dll' if os.name == 'nt' else 'libreference_probe.so')).read_bytes()).hexdigest() == signature['probe']['probe_sha256']
    cache = ROOT/'.reference/meos_isochoric_source_v2.json'
    if args.reuse_source:
        saved = json.loads(cache.read_text(encoding='utf-8'))
        # Audit code may evolve without rerunning fatal source solvers. Verify
        # the unchanged source binary/worker/input identities, then recompute
        # every target with original TV kernels and compare actual call inputs.
        identity = lambda x: {k:v for k,v in x.items() if k != 'generator_sha256'}
        assert identity(saved['signature']) == identity(signature), 'stale source cache'
        with DeadlineReferenceProcess(args.compiler_directory) as oracle:
            seeds = prepare(oracle, data)
        schedule = list(requests(seeds, data))
        assert [e['request'] for e in saved['records']] == [q for q, seed in schedule], 'changed source call inputs'
        records = [dict(old, seed=seed) for old, (q, seed) in zip(saved['records'], schedule)]
        collection_signature = saved['signature']
    else:
        with DeadlineReferenceProcess(args.compiler_directory) as oracle:
            seeds = prepare(oracle, data)
            records = []
            for i, (q, seed) in enumerate(requests(seeds, data)):
                records.append(dict(request=q, seed=seed, result=oracle.call(q, timeout=4.)))
                if i % 50 == 0:
                    print(i, 'source UV/SV calls', flush=True)
        atomic_text(cache, json.dumps(dict(signature=signature, seeds=seeds, records=records), allow_nan=False)+'\n')
        collection_signature = signature
    if args.source_only:
        print('Source records ready:', len(records), 'attempts with independently prepared targets', flush=True)
        return
    subprocess.run(['cargo', 'build', '--release', '--example', 'parity', '-j', '2'], cwd=ROOT, check=True)
    exe = ROOT/'target/release/examples'/('parity.exe' if os.name == 'nt' else 'parity')
    output = subprocess.run([str(exe)], input='\n'.join(encode(e['request']) for e in records)+'\n',
                            text=True, capture_output=True, check=True).stdout.splitlines()
    assert len(output) == len(records)
    atomic_text(ROOT/'.reference/meos_isochoric_native_raw.json', json.dumps(dict(signature=signature, output=output), allow_nan=False)+'\n')
    matches, phase_gaps, recoveries, differences, native_errors, independent, outcomes = [], [], [], [], [], [], {}
    independent_source, unqualified_source_kernels = [], []
    with DeadlineReferenceProcess(args.compiler_directory) as oracle:
        for e, line in zip(records, output):
            q, r = e['request'], e['result']; case = encode(q)
            if line.startswith('ERROR'):
                native_errors.append(dict(case=case, native_error=line, source=r)); continue
            a = list(map(float, line.split(',')))
            assert valid_native_point(a), (case, line)
            t, p, two, bv, bl, v, u, s, phase, vl, vv, cv, raw_bv, raw_bl = a
            volumes = [vl, vv] if two else [vl]
            model = q['eos']+':'+q['component']; gas = data[model]['R']
            checked = equilibrium_kernel(data, model, t, p, volumes, bv, [kernel(oracle, model, t, x) for x in volumes])
            assert abs(v/q['volume']-1.) < 2e-8, (case, 'native specified V closure', v, q['volume'])
            for label, x, y, scale in [('V', v, checked['volume'], q['volume']), ('U', u, checked['u'], gas*t),
                    ('S', s, checked['s'], gas), ('specified caloric target', u if q['kind'] == 'UV' else s, q['target'], gas*(t if q['kind'] == 'UV' else 1.)),
                    ('Cv', cv, checked['cv'], max(abs(cv), gas))]:
                assert abs(x-y) <= 3e-10*scale, (case, label, x, y)
            assert q['minimum_temperature'] <= t <= q['maximum_temperature'] and q['minimum_pressure'] <= p <= q['maximum_pressure'], (case, 'bounds', a)
            independent.append(dict(case=case, native=a, source_kernel_check=checked))
            got = a[:8]+[phase]; reduced = q['eos'] == 'LJTS_REDUCED'
            qualified_source = source_closed(q, r, data)
            kernel_failure = None
            if qualified_source:
                source_t, source_p = r['values'][:2]
                source_states = [kernel(oracle, model, source_t, x) for x in r['phase_volumes']]
                try:
                    source_check = equilibrium_kernel(data, model, source_t, source_p,
                        r['phase_volumes'], r['values'][3], source_states)
                except AssertionError as error:
                    kernel_failure = str(error)
                    unqualified_source_kernels.append(dict(case=case, reason=kernel_failure, kernels=source_states))
                    qualified_source = False
                else:
                    independent_source.append(dict(case=case, source_kernel_check=source_check))
            if not qualified_source:
                outcome = source_outcome(r)
                outcomes[outcome_key(q)] = ('returned_kernel_nonclosing' if kernel_failure is not None else 'returned_nonclosing') if outcome == 'returned' else outcome
                seed = e['seed']; expected = seed['expected']
                wanted = [seed['t'], seed['p'], two, bv, bl, expected['volume'], expected['u'], expected['s'], phase]
                bad = coordinate_bad(wanted, got, reduced)
                if bad:
                    differences.append(dict(case=case, native=got, expected_seed=wanted, differences=bad, source=r))
                else:
                    recoveries.append([case, seed['t'], seed['p'], expected['volume'], expected['u'], expected['s']])
                continue
            wanted = r['values']+[r['raw_result'][2]]
            bad = coordinate_bad(wanted, got, reduced)
            if not bad:
                matches.append([case, *wanted]); continue
            if bounded_endpoint_gap(wanted, got, reduced):
                outcomes[outcome_key(q)] = 'endpoint_phase_difference'
                phase_gaps.append([case, *wanted])
            else:
                differences.append(dict(case=case, source=wanted, native=got, differences=bad))
    meta = dict(signature, source_collection_signature=collection_signature, attempts=len(records), records=records, independent_native_kernel_checks=independent,
                source_outcomes=outcomes, differences=differences, native_errors=native_errors,
                independent_source_kernel_checks=independent_source, unqualified_source_kernels=unqualified_source_kernels,
                component_database_sha256=digest(ROOT/'src/compdatadb.f90'),
                comparison_relative_tolerance=3e-7, native_scaled_caloric_tolerance=3e-10, native_relative_volume_tolerance=2e-8,
                strict_matches=len(matches), endpoint_phase_differences=len(phase_gaps), native_recoveries=len(recoveries),
                native_sources_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in sorted((ROOT/'rust').rglob('*.rs'))+[ROOT/'Cargo.toml', ROOT/'Cargo.lock']},
                qualification='Partial pure UV/SV grid. Source failures and endpoint phase differences remain explicit gaps; every native output is independently closed using original TV kernels.')
    for name, rows, header in [('multiparameter_isochoric', matches, FIELDS),
            ('multiparameter_isochoric_phase_gaps', phase_gaps, FIELDS),
            ('multiparameter_isochoric_recovered', recoveries, ['case', 'seed_T', 'seed_P', 'V', 'U', 'S'])]:
        fixture = ROOT/'tests/fixtures'/f'{name}.csv'
        atomic_csv(fixture, header, rows)
        info = meta if name == 'multiparameter_isochoric' else dict(signature, qualification=meta['qualification'])
        atomic_text(fixture.with_suffix('.metadata.json'), json.dumps(dict(info, cases=len(rows),
            fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(), generator='tools/rust_port/validate_meos_isochoric.py'), indent=2, allow_nan=False)+'\n')
    atomic_text(ROOT/'.reference/meos_isochoric_report.json', json.dumps(meta, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(attempts=len(records), strict_matches=len(matches), phase_gaps=len(phase_gaps),
        recoveries=len(recoveries), native_errors=len(native_errors), differences=len(differences),
        independent=len(independent), independent_source=len(independent_source)), indent=2), flush=True)
    assert not native_errors, native_errors[:3]
    assert not differences, differences[:3]
    if outcomes:
        if not args.allow_known_gaps:
            raise SystemExit('Partial source qualification; inspect complete report and exact declared outcomes')
        allowed = json.loads((ROOT/'tools/rust_port/meos_isochoric_known_gaps.json').read_text(encoding='utf-8'))
        if not allowed_outcomes(outcomes, allowed['outcomes']):
            unexpected = {
                k: {'actual': v, 'declared': allowed['outcomes'].get(k)}
                for k, v in sorted(outcomes.items())
                if not outcome_allowed(v, allowed['outcomes'].get(k))
            }
            print(json.dumps({'unexpected_outcomes_count': len(unexpected), 'unexpected_outcomes': unexpected}, indent=2), flush=True)
            raise AssertionError(f'new or changed source/phase outcomes require review ({len(unexpected)} unexpected)')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--compiler-directory')
    parser.add_argument('--reuse-source', action='store_true')
    parser.add_argument('--source-only', action='store_true')
    parser.add_argument('--allow-known-gaps', action='store_true')
    main(parser.parse_args())
