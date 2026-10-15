"""Qualify default ideal caloric correlations for every fluid in this fork."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
from reference import ROOT, ReferenceProcess


def cases():
    for path in sorted((ROOT / 'fluids').glob('*.json')):
        name = json.loads(path.read_text())['ident']
        for t in [20., 30., 30.000001, 80., 150., 298.15, 500., 1000.]:
            for p in [1e3, 1e5, 1e7]:
                yield {'mode': 'IDEAL', 'eos': 'PR', 'components': name,
                       'temperature': t, 'pressure': p}


def validate(compiler_directory):
    inputs, expected, reference_errors = [], [], []
    with ReferenceProcess(compiler_directory) as oracle:
        for request in cases():
            result = oracle.call(request)
            if 'values' not in result or 'calculation_error' in result:
                reference_errors.append({'request': request, **result})
                continue
            inputs.append('|'.join(['PR', request['components'], str(request['temperature']),
                                   str(request['pressure']), '1', 'IDEAL']))
            expected.append(result['values'])
    subprocess.run(['cargo', 'build', '--release', '--example', 'parity'], cwd=ROOT, check=True)
    executable = ROOT / 'target/release/examples' / ('parity.exe' if os.name == 'nt' else 'parity')
    result = subprocess.run([str(executable)], input='\n'.join(inputs)+'\n', text=True,
                            capture_output=True, check=True)
    actual = result.stdout.splitlines()
    assert len(actual) == len(expected)
    absolute = [1e-7, 1e-9, 1e-9, 1e-10, 1e-12]
    relative = 2e-8
    failures = []
    for case, ref, row in zip(inputs, expected, actual):
        if row.startswith('ERROR'):
            failures.append({'case': case, 'error': row})
            continue
        native = list(map(float, row.split(',')))
        assert len(native) == len(ref) == len(absolute)
        for i, (a, b) in enumerate(zip(ref, native)):
            if not abs(a-b) <= absolute[i] + relative * max(abs(a), abs(b)):
                failures.append({'case': case, 'quantity': i, 'reference': a, 'rust': b})
    output = ROOT / 'tests/fixtures/ideal.csv'
    with output.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter=';', lineterminator='\n')
        writer.writerow(['case', 'h', 's', 'cp', 'dsdt_baseline', 'dsdp'])
        writer.writerows([[case, *values] for case, values in zip(inputs, expected)])
    metadata = json.loads((ROOT / '.reference/metadata.json').read_text())
    metadata.update({'attempted_cases': len(inputs) + len(reference_errors),
                     'cases': len(inputs), 'reference_errors': reference_errors,
                     'absolute_tolerance': absolute, 'relative_tolerance': relative,
                     'fixture_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
                     'generator': 'tools/rust_port/validate_ideal.py'})
    output.with_suffix('.metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    (ROOT / '.reference/ideal_failures.json').write_text(json.dumps(failures, indent=2)+'\n')
    print(json.dumps({'cases': len(inputs), 'reference_errors': len(reference_errors),
                      'failures': len(failures)}, indent=2))
    if failures:
        print(json.dumps(failures[:12], indent=2))
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
