"""Check the committed numerical fixture bytes without a Fortran runtime."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def verify():
    checked = 0
    unverified = set((ROOT/'tests/fixtures').glob('*.csv'))
    for metadata in sorted((ROOT/'tests/fixtures').glob('*.metadata.json')):
        data = json.loads(metadata.read_text(encoding='utf-8'))
        if 'fixture_sha256' not in data:
            continue
        fixture = metadata.with_name(metadata.name.removesuffix('.metadata.json')+'.csv')
        content = fixture.read_bytes()
        if b'\r' in content or not content.endswith(b'\n'):
            raise ValueError(f'{fixture.name}: fixtures require LF line endings')
        if hashlib.sha256(content).hexdigest() != data['fixture_sha256']:
            raise ValueError(f'{fixture.name}: SHA-256 differs from metadata')
        with fixture.open(encoding='utf-8', newline='') as handle:
            count = sum(1 for _ in csv.reader(handle, delimiter=';'))-1
        if count != data.get('cases', data.get('case_count')):
            raise ValueError(f'{fixture.name}: row count differs from metadata')
        unverified.discard(fixture)
        checked += 1
    if not checked:
        raise ValueError('No reference fixtures found')
    if unverified:
        raise ValueError('Fixtures missing integrity metadata: '+str(sorted(unverified)))
    print(f'{checked} fixture hashes, row counts and LF line endings verified')


if __name__ == '__main__':
    verify()
