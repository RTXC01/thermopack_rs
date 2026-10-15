"""Exact-source TP phase decisions, phase fractions and compositions."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
from reference import ROOT,ReferenceProcess
from validate_cubic import cases as property_cases


def cases():
    seen=set()
    for eos,comps,t,p,z,_ in property_cases():
        key=(eos,comps,t,p,tuple(z))
        if key not in seen:
            seen.add(key)
            yield eos,comps,t,p,z
    for eos in ["PR","SRK","RK","VDW"]:
        for comps,z,temps,pressures in [
            ("C1,C2",[0.5,0.5],[110.,150.,180.,220.,260.,300.,400.],[1e4,1e5,5e5,1e6,2e6,5e6,1e7]),
            ("C1,C2,C3",[0.5,0.3,0.2],[180.,220.,260.,300.],[5e5,1e6,2e6,5e6]),
            ("CO2,N2",[0.8,0.2],[250.,280.,300.,310.,320.],[1e6,4e6,6e6,8e6,1e7]),
            ("C1,C2,NC5,NC10",[0.5,0.2,0.2,0.1],[200.,300.,400.],[1e5,1e6,5e6]),
        ]:
            for t in temps:
                for p in pressures:
                    key=(eos,comps,t,p,tuple(z))
                    if key not in seen:
                        seen.add(key)
                        yield eos,comps,t,p,z
    # The exact difficult TP points from unittests/test_tp_flash.pf.
    for t,p,z in [
        (302.46137724669268,7705980.6837845352,
         [0.98199899999065965,6.0000000025901721e-3,1.0000000003932113e-6,
          5.0000000019329808e-4,4.0000000015823992e-4,9.9999999995230694e-5,
          1.0000000006101501e-2,1.0000000003014362e-3]),
        (302.46161994265759,7706002.2603671746,
         [0.98199899997992279,6.0000000055675665e-3,1.0000000008460226e-6,
          5.0000000041550564e-4,4.0000000034013847e-4,9.9999999989748805e-5,
          1.0000000013115396e-2,1.0000000006479682e-3]),
    ]:
        yield "PR","CO2,N2,O2,AR,CO,H2S,H2,C1",t,p,z


def validate(compiler_directory):
    inputs,expected,swappable=[],[],[]
    reference_failures=[]
    with ReferenceProcess(compiler_directory) as oracle:
        for eos,comps,t,p,z in cases():
            case='|'.join([eos,comps,str(t),str(p),','.join(map(str,z)),"TP"])
            result=oracle.call({'mode':'TP','eos':eos,'components':comps,
                                'temperature':t,'pressure':p,'feed':z})
            inputs.append(case)
            expected.append(result.get('values'))
            swappable.append(result.get('phase_permutation_allowed',False))
            if 'reference_error' in result:
                reference_failures.append({'case':case,**result['reference_error']})
    subprocess.run(["cargo","build","--release","--example","parity"],cwd=ROOT,check=True)
    executable=ROOT/"target/release/examples"/("parity.exe" if os.name=="nt" else "parity")
    result=subprocess.run([str(executable)],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True)
    actual=result.stdout.splitlines()
    assert len(actual)==len(expected)
    failures=[]
    permutations=[]
    max_error=0.0
    for case,ref,line,can_swap in zip(inputs,expected,actual,swappable):
        if ref is None:
            continue
        if line.startswith("ERROR"):
            failures.append({"case":case,"error":line})
            continue
        native=list(map(float,line.split(',')))
        assert len(native)==len(ref)
        def matches(values):
            return all(abs(a-b)<= (0.0 if i==0 else 2e-7+2e-7*max(abs(a),abs(b)))
                       for i,(a,b) in enumerate(zip(ref,values)))
        if can_swap and not matches(native) and native[0]==0:
            n=(len(native)-3)//2
            swapped=[native[0],native[2],native[1],*native[3+n:],*native[3:3+n]]
            if matches(swapped):
                permutations.append(case)
                native=swapped
        for i,(a,b) in enumerate(zip(ref,native)):
            tolerance=0.0 if i==0 else 2e-7+2e-7*max(abs(a),abs(b))
            max_error=max(max_error,abs(a-b))
            if not abs(a-b)<=tolerance:
                failures.append({"case":case,"quantity":i,"reference":a,"rust":b})
    output=ROOT/"tests/fixtures/tp_flash.csv"
    with output.open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter=';',lineterminator='\n')
        writer.writerow(["case","phase_permutation_allowed","phase,betaV,betaL,x[n],y[n]..."])
        writer.writerows([[case,int(can_swap),*row] for case,row,can_swap in zip(inputs,expected,swappable) if row is not None])
    metadata=json.loads((ROOT/".reference/metadata.json").read_text())
    metadata.update({"attempted_cases":len(inputs),"cases":len(inputs)-len(reference_failures),
                     "reference_failures":reference_failures,"max_absolute_error_after_permutation":max_error,
                     "phase_label_permutations":permutations,
                     "absolute_tolerance":2e-7,"relative_tolerance":2e-7,"phase_flag_tolerance":0,
                     "fixture_sha256":hashlib.sha256(output.read_bytes()).hexdigest(),
                     "generator":"tools/rust_port/validate_flash.py"})
    output.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/".reference/flash_failures.json").write_text(json.dumps(failures,indent=2)+'\n')
    print(json.dumps({"attempted_cases":len(inputs),"reference_errors":len(reference_failures),
                      "phase_label_permutations":len(permutations),"failures":len(failures),
                      "max_absolute_error_after_permutation":max_error},indent=2))
    if failures:
        print(json.dumps(failures[:16],indent=2))
        raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--compiler-directory')
    validate(parser.parse_args().compiler_directory)
