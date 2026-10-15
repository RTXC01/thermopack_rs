"""Qualify all supported Cp setters and ordered reference-state changes."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
from reference import ROOT, ReferenceProcess

COEFFICIENTS = {
    1: [8., .02, -1e-5, 2e-9],
    2: [3.7, 1.2, .001, 1e-6, -1e-10, 2e-14],
    3: [.2, 1e-3, -1e-7],
    4: [25., .02, -1e-5, 2e-9],
    5: [1.2, .001, -1e-7, 2e-11, 10.],
    6: [1.2, .001, -1e-7, 2e-11, -1e-15],
    7: [30000., 20000., 700., 15000., 300.],
    8: [2.5, .001, -1e-7, 2e-11, -1e-15],
    9: [25., .02, 100., 20.],
    10: [],
    11: [30., 2., .5, -.02, .001],
    12: [2.5, 1.5, 500., .7, 2000.],
}


def actions(corr):
    cp = ['cp', corr, COEFFICIENTS[corr]]
    return [[cp], [cp, ['h', -4321.]], [cp, ['s', 222., 101325.]],
            [cp, ['h', 0.], ['s', 0., 1e5]],
            [['h', -4321.], cp], [['s', 222., 1e5], cp]]


def encode(action):
    if action[0] == 'cp':
        return ':'.join(['cp', str(action[1]), ','.join(map(str, action[2]))])
    return ':'.join(map(str, action))


def cases():
    for corr in range(1, 13):
        names = ['N-H2', 'E-H2', 'P-H2', 'O-H2'] if corr == 10 else ['C1']
        for name in names:
            for sequence in actions(corr):
                for t in [80., 298.15, 500., 1000.]:
                    for p in [1e3, 1e5]:
                        yield {'mode': 'IDEAL_CONFIG', 'eos': 'PR', 'components': name,
                               'actions': sequence, 'temperature': t, 'pressure': p}


def validate(compiler_directory):
    inputs, expected = [], []
    with ReferenceProcess(compiler_directory) as oracle:
        for request in cases():
            result = oracle.call(request)
            assert 'values' in result and 'calculation_error' not in result, (request, result)
            inputs.append('|'.join(['PR', request['components'], str(request['temperature']),
                                    str(request['pressure']), '/'.join(map(encode, request['actions'])),
                                    'IDEAL_CONFIG']))
            expected.append(result['values'])
    subprocess.run(['cargo', 'build', '--release', '--example', 'parity'], cwd=ROOT, check=True)
    executable = ROOT / 'target/release/examples' / ('parity.exe' if os.name == 'nt' else 'parity')
    result = subprocess.run([str(executable)], input='\n'.join(inputs)+'\n', text=True,
                            capture_output=True, check=True)
    actual = result.stdout.splitlines()
    assert len(actual) == len(expected)
    absolute = [1e-7, 1e-9, 1e-9, 1e-10, 1e-12, 1e-10, 1e-10] + [0.] * 22
    relative = 2e-8
    failures = []
    for case, ref, row in zip(inputs, expected, actual):
        native = list(map(float, row.split(',')))
        assert len(native) == len(ref) == len(absolute)
        for i, (a, b) in enumerate(zip(ref, native)):
            tolerance = absolute[i] + relative * max(abs(a), abs(b)) if i < 7 else 0.
            if not abs(a-b) <= tolerance:
                failures.append({'case': case, 'quantity': i, 'reference': a, 'rust': b})
    output = ROOT / 'tests/fixtures/ideal_configuration.csv'
    with output.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter=';', lineterminator='\n')
        writer.writerow(['case', 'h', 's', 'cp', 'dsdt_baseline', 'dsdp', 'formation_enthalpy',
                         'standard_entropy', 'correlation', 'coefficients...'])
        writer.writerows([[case, *values] for case, values in zip(inputs, expected)])
    metadata = json.loads((ROOT / '.reference/metadata.json').read_text())
    metadata.update({'cases': len(inputs), 'reference_errors': [],
                     'absolute_tolerance': absolute, 'relative_tolerance': relative,
                     'exact_quantities_from_index': 7,
                     'fixture_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
                     'generator': 'tools/rust_port/validate_ideal_configuration.py'})
    output.with_suffix('.metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    (ROOT / '.reference/ideal_configuration_failures.json').write_text(json.dumps(failures, indent=2)+'\n')
    print(json.dumps({'cases': len(inputs), 'failures': len(failures)}, indent=2))
    if failures:
        print(json.dumps(failures[:12], indent=2))
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
