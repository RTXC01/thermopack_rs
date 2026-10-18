"""Compare raw fork-order TP results against all frozen source TP snapshots."""
import csv
import hashlib
import json
import os
import subprocess
from reference import ROOT


def validate():
    cases=[];fixture_hashes={}
    for name in ['tp_flash','translated_tp','quantum_tp']:
        path=ROOT/'tests/fixtures'/(name+'.csv')
        info=json.loads(path.with_suffix('.metadata.json').read_text())
        digest=hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest==info['fixture_sha256']
        fixture_hashes[name]=digest
        for row in list(csv.reader(path.open(),delimiter=';'))[1:]:
            fields=row[0].split('|')
            config=fields[5].split('@');config[0]+='_BASELINE'
            fields[5]='@'.join(config)
            # Flash snapshots record a swap-permitted flag before properties;
            # this comparison always requires the original raw ordering.
            cases.append((name,'|'.join(fields),list(map(float,row[2:]))))
    subprocess.run(['cargo','build','--release','--example','parity'],cwd=ROOT,check=True)
    executable=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    process=subprocess.run([str(executable)],input='\n'.join(row[1] for row in cases)+'\n',text=True,capture_output=True,check=True)
    lines=process.stdout.splitlines();assert len(lines)==len(cases)
    counts={name:dict(matched_cases=0,failed_cases=0) for name in fixture_hashes};failures=[]
    for (name,case,expected),line in zip(cases,lines):
        failure=None
        if line.startswith('ERROR'):
            failure=dict(group=name,case=case,error=line,reference=expected)
        else:
            actual=list(map(float,line.split(',')));assert len(actual)==len(expected)
            bad=[i for i,(a,b) in enumerate(zip(expected,actual)) if not abs(a-b)<=(0. if i==0 else 2e-7+2e-7*max(abs(a),abs(b)))]
            if bad:failure=dict(group=name,case=case,quantities=bad,reference=expected,rust=actual)
        counts[name]['failed_cases' if failure else 'matched_cases']+=1
        if failure:failures.append(failure)
    report=dict(summary=counts,failures=failures,source_fixture_sha256=fixture_hashes,
                absolute_tolerance=2e-7,relative_tolerance=2e-7,phase_flag_tolerance=0,
                phase_permutation_allowed=False,
                native_sources_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                                       for p in sorted((ROOT/'rust').rglob('*.rs'))},
                rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),generator='tools/rust_port/validate_baseline_tp.py')
    (ROOT/'tests/fixtures/baseline_tp_differences.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(summary=counts,unresolved_cases=len(failures)),indent=2))
    if failures:
        print(json.dumps(failures[:8],indent=2));raise SystemExit(1)


if __name__=='__main__':validate()
