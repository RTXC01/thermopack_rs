"""Qualify cubic fugacity T/P/mole-number derivatives against the fork."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
from reference import ROOT, load
from validate_cubic import cases


def validate(compiler_directory):
    cubic=load(compiler_directory)
    inputs,expected=[],[]
    model,previous=None,None
    for eos,comps,t,p,x,phase in cases():
        if previous!=(eos,comps):
            model=cubic(comps,"SRK" if eos=="RK" else eos,
                        alpha="RK" if eos=="RK" else "Classic",volume_shift=False)
            previous=eos,comps
        flag=model.LIQPH if phase=="L" else model.VAPPH
        _,dt,dp,dn=model.thermo(t,p,x,flag,dlnfugdt=True,dlnfugdp=True,dlnfugdn=True)
        inputs.append("|".join([eos,comps,str(t),str(p),",".join(map(str,x)),"D_"+phase]))
        expected.append([*dt,*dp,*dn.flatten()])
    subprocess.run(["cargo","build","--release","--example","parity"],cwd=ROOT,check=True)
    executable=ROOT/"target/release/examples"/("parity.exe" if os.name=="nt" else "parity")
    result=subprocess.run([str(executable)],input="\n".join(inputs)+"\n",text=True,
                          capture_output=True,check=True)
    actual=[list(map(float,row.split(','))) for row in result.stdout.splitlines()]
    assert len(actual)==len(expected)
    failures=[]
    worst={"dt":0.0,"dp":0.0,"dn":0.0}
    for case,ref,native in zip(inputs,expected,actual):
        n=len(case.split('|')[1].split(','))
        assert len(ref)==len(native)
        for i,(a,b) in enumerate(zip(ref,native)):
            label="dt" if i<n else "dp" if i<2*n else "dn"
            absolute=1e-9 if label=="dt" else 1e-11 if label=="dp" else 1e-8
            error=abs(a-b)/max(abs(a),abs(b),1.0)
            worst[label]=max(worst[label],error)
            if not abs(a-b)<=absolute+2e-7*max(abs(a),abs(b)):
                failures.append({"case":case,"quantity":label,"reference":a,"rust":b})
    output=ROOT/"tests/fixtures/fugacity_derivatives.csv"
    with output.open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter=';',lineterminator='\n')
        writer.writerow(["case","dt[n],dp[n],dn[n*n]..."])
        writer.writerows([[case,*row] for case,row in zip(inputs,expected)])
    metadata=json.loads((ROOT/".reference/metadata.json").read_text())
    metadata.update({"cases":len(inputs),"worst_scaled_error":worst,
                     "absolute_tolerances":{"dt":1e-9,"dp":1e-11,"dn":1e-8},"relative_tolerance":2e-7,
                     "fixture_sha256":hashlib.sha256(output.read_bytes()).hexdigest(),
                     "generator":"tools/rust_port/validate_derivatives.py"})
    output.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/".reference/derivative_failures.json").write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps({"cases":len(inputs),"failures":len(failures),"worst":worst},indent=2))
    if failures:
        print(json.dumps(failures[:10],indent=2))
        raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
