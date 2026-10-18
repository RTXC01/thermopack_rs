"""Compare pure PeTS flashes against isolated exact-source public interfaces."""
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess

from reference import ROOT, ReferenceProcess
from validate_alpha import verify_probe

ABSOLUTE = [1e-5, .01, 0., 1e-7, 1e-9, 1e-9, 1e-3, 1e-5]
RELATIVE = 2e-7


def cases():
    for sigma, eps in [(3.42e-10, 136.), (2.9e-10, 100.), (4.1e-10, 240.)]:
        common = dict(parameters=[sigma, eps], redefine_critical=True)
        for tr in [.5, .8, 1., 1.09, 1.1, 1.5, 3.]:
            for p in [.001, 1., 1e3, 1e5, 1e6, 5e6, 1e8]:
                yield dict(common, mode='PETS_TP', temperature=tr*eps, pressure=p)
                for mode in ['PETS_PH', 'PETS_PS']:
                    for guessed in [False, True]:
                        yield dict(common, mode=mode, temperature=tr*eps, pressure=p, guessed=guessed)
        for tr in [.5, .7, .9, .99, .999]:
            for beta in [0., .1, .5, .9, 1.]:
                for mode in ['PETS_PH', 'PETS_PS']:
                    for guessed in [False, True]:
                        yield dict(common, mode=mode, temperature=148.11064915309302/136.*eps*tr,
                                   pressure=0., quality=beta, guessed=guessed)


def encode(request, result):
    guess = str(request['temperature']*request.get('guess_factor', 1.2)) if request.get('guessed', True) else 'none'
    return '|'.join(['PETS', ','.join(map(str, request['parameters'])), str(request['temperature']),
                     str(result['values'][1]), '1', request['mode']+'@'+str(result['specification'])+'@'+guess])


def validate(directory):
    probe = verify_probe()
    records, errors = [], []
    with ReferenceProcess(directory) as oracle:
        for request in cases():
            result = oracle.call(request)
            if 'calculation_error' in result or 'values' not in result:
                errors.append(dict(request=request, **result))
            else:
                records.append((request, encode(request, result), result))
    subprocess.run(['cargo', 'build', '--release', '--example', 'parity'], cwd=ROOT, check=True)
    exe = ROOT/'target/release/examples'/('parity.exe' if os.name == 'nt' else 'parity')
    lines = subprocess.run([str(exe)], input='\n'.join(c for _, c, _ in records)+'\n', text=True,
                           capture_output=True, check=True).stdout.splitlines()
    assert len(lines) == len(records)
    failures = []
    for (request, case, result), line in zip(records, lines):
        if line.startswith('ERROR'):
            failures.append(dict(case=case, error=line, reference=result['values']))
            continue
        actual = list(map(float, line.split(',')))
        expected = result['values']
        assert len(actual) == len(expected)
        for i, (a, b) in enumerate(zip(expected, actual)):
            if not math.isfinite(b) or abs(a-b) > ABSOLUTE[i]+RELATIVE*max(abs(a), abs(b)):
                failures.append(dict(case=case, quantity=i, reference=a, rust=b))
    for mode in ['PETS_TP', 'PETS_PH', 'PETS_PS']:
        group = [(request, case, result) for request, case, result in records if request['mode'] == mode]
        path = ROOT/'tests/fixtures'/(mode.lower()+'.csv')
        with path.open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter=';', lineterminator='\n')
            writer.writerow(['case', 'T', 'P', 'two_phase', 'betaV', 'vL_or_v', 'vV_or_v', 'H', 'S'])
            writer.writerows([[case, *result['values']] for _, case, result in group])
        source_errors = [e for e in errors if e['request']['mode'] == mode]
        metadata = json.loads((ROOT/'.reference/metadata.json').read_text())
        metadata.update(cases=len(group), attempted_cases=len(group)+len(source_errors), reference_errors=source_errors,
                        probe=probe, generator='tools/rust_port/validate_pets_flash.py',
                        absolute_tolerances=ABSOLUTE, relative_tolerance=RELATIVE,
                        fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        convention='canonical homogeneous phase amounts/volumes, including quality within 1e-8 of either endpoint; original phase flags and arrays preserved per case; finite/stable pressure and caloric closure required',
                        source_results=[dict(case=case, **{k: v for k, v in result.items() if k != 'values'}) for _, case, result in group])
        path.with_suffix('.metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    (ROOT/'.reference/pets_flash_failures.json').write_text(json.dumps(failures, indent=2)+'\n')
    print(json.dumps(dict(finite_cases=len(records), reference_errors=len(errors), failures=len(failures),
                         by_mode={m: sum(r['mode'] == m for r, _, _ in records) for m in ['PETS_TP', 'PETS_PH', 'PETS_PS']}), indent=2))
    if failures:
        print(json.dumps(failures[:10], indent=2))
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
