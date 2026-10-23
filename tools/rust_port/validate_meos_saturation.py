"""Public pure MEOS coexistence, including a scaled source reference for reduced LJTS."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess
from validate_meos import MODELS

TR=[.45,.5,.6,.75,.9,.95,.99,.999,.9999]

def finite(value):
    return isinstance(value,(float,int)) and math.isfinite(value)

def key(q):
    mode='SAT_P' if 'temperature' in q else 'SAT_T';value=q.get('temperature',q.get('pressure'))
    return f"MEOS|{q['eos']}:{q['component']}|{value!r}|0|1|{mode}"

def validate(directory,reuse_source):
    identity=json.loads((ROOT/'.reference/metadata.json').read_text())
    worker=hashlib.sha256((ROOT/'tools/rust_port/oracle_worker.py').read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    cache=ROOT/'.reference/meos_saturation_source.json'
    if reuse_source:
        saved=json.loads(cache.read_text());assert saved['reference']==identity and saved['oracle_sha256']==worker
        records=saved['records']
    else:
        records=[]
        with DeadlineReferenceProcess(directory) as oracle:
            for eos,component in MODELS:
                initial=oracle.call(dict(mode='MEOS',eos=eos,component=component,temperature=300.,volume=1.));assert 'metadata' in initial,initial
                for tr in TR:
                    q=dict(mode='MEOS_SAT',eos=eos,component=component,temperature=tr*initial['metadata'][0]);r=oracle.call(q,timeout=8.)
                    records.append(dict(request=q,result=r))
                    if r.get('values') and all(map(finite,r['values'])) and r['values'][1]>0:
                        q=dict(mode='MEOS_SAT',eos=eos,component=component,pressure=r['values'][1]);records.append(dict(request=q,result=oracle.call(q,timeout=8.)))
                print(eos,component,len(records),'source attempts',flush=True)
        cache.write_text(json.dumps(dict(reference=identity,oracle_sha256=worker,records=records))+'\n')
    good=[];failed=[];unqualified=[];warnings=[];censored=[]
    for e in records:
        r=e['result'];v=r.get('values')
        if 'reference_timeout' in r:censored.append(e);continue
        if not v:failed.append(e);continue
        if not all(map(finite,v)) or min(v[:4])<=0 or (v[3]-v[2])/(v[3]+v[2])<1e-6 or not r['mechanically_stable'] or abs(r['fugacity_residual'])>2e-7:
            reasons=[]
            if not all(map(finite,v)):reasons.append('nonfinite_output')
            if all(map(finite,v[2:4])) and v[3]+v[2]>0 and (v[3]-v[2])/(v[3]+v[2])<1e-6:reasons.append('unresolved_phase_volume_gap')
            if not r['mechanically_stable']:reasons.append('unstable_endpoint')
            if not finite(r['fugacity_residual']) or abs(r['fugacity_residual'])>2e-7:reasons.append('chemical_potential_closure')
            unqualified.append(dict(e,qualification_reasons=reasons));continue
        if max(map(abs,r['pressure_errors']))>1e-5+1e-7*v[1]:warnings.append(e)
        good.append(e)
    subprocess.run(['cargo','build','--release','--example','parity','-j','2'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    output=subprocess.run([str(exe)],input='\n'.join(key(e['request']) for e in records)+'\n',capture_output=True,text=True,check=True).stdout.splitlines();assert len(output)==len(records)
    native=dict(zip((key(e['request']) for e in records),output))
    independent=[];native_errors=[];recovered=[]
    good_keys={key(e['request']) for e in good};failed_keys={key(e['request']) for e in failed}
    with DeadlineReferenceProcess(directory) as oracle:
        for e in records:
            q=e['request'];case=key(q);encoded=native[case]
            if encoded.startswith('ERROR'):
                native_errors.append(dict(case=case,native_error=encoded,original_result=e['result']));continue
            point=list(map(float,encoded.split(',')))
            checked=oracle.call(dict(mode='MEOS_COEX_CHECK',eos=q['eos'],component=q['component'],point=point),timeout=8.)
            if 'pressure_errors' not in checked:raise AssertionError(('independent source check failed',case,checked))
            assert checked['mechanically_stable'] and abs(checked['chemical_potential_residual'])<2e-9,(case,checked)
            barrier=checked['nonconvex_barrier'];assert barrier and point[2]<barrier['volume']<point[3] and barrier['pressure_dv']>0.,(case,checked)
            ptol=5e-7+2e-8*point[1]+512.*math.ulp(1.)*checked['gas_constant']*point[0]/point[2]
            assert max(map(abs,checked['pressure_errors']))<ptol,(case,checked,ptol)
            assert abs(checked['dpdt']-point[4])<1e-7+2e-7*max(abs(checked['dpdt']),abs(point[4])),(case,checked,point)
            independent.append(dict(case=case,native_point=point,source_check=checked))
            if case not in good_keys:
                recovered.append([case,*point,'source_error' if case in failed_keys else 'source_unqualified',checked['chemical_potential_residual'],*checked['pressure_errors']])
    differences=[];rows=[];maximum=[0.]*5
    for e in good:
        encoded=native[key(e['request'])]
        q=e['request'];r=e['result'];case=key(q);expected=r['values'];bad=[]
        if encoded.startswith('ERROR'):bad=[encoded]
        else:
            actual=list(map(float,encoded.split(',')));assert len(actual)==5
            for i,(a,b) in enumerate(zip(expected,actual)):
                floor=[2e-7,1e-2,1e-12,1e-12,1e-4][i]
                maximum[i]=max(maximum[i],abs(a-b)/max(abs(a),abs(b)))
                if not math.isfinite(b) or abs(a-b)>floor+2e-6*max(abs(a),abs(b)):bad.append([i,a,b])
        if bad:differences.append(dict(case=case,differences=bad))
        rows.append([case,*expected])
    fixture=ROOT/'tests/fixtures/meos_saturation.csv'
    with fixture.open('w',newline='') as handle:
        w=csv.writer(handle,delimiter=';',lineterminator='\n');w.writerow(['case','temperature','pressure','liquid_volume','vapor_volume','dpdt']);w.writerows(rows)
    meta=dict(identity,cases=len(rows),fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),generator='tools/rust_port/validate_meos_saturation.py',
        oracle_sha256=worker,attempts=len(records),reference_errors=failed,censored_calls=censored,unqualified_returns=unqualified,pressure_warnings=warnings,native_differences=differences,
        native_errors=native_errors,independent_native_state_checks=independent,
        relative_tolerance=2e-6,absolute_tolerances=[2e-7,1e-2,1e-12,1e-12,1e-4],maximum_relative_differences=maximum,
        reference_method='Public bubble_pressure/bubble_temperature and source phase volumes, fugacity, pressure and latent heat. Reduced LJTS is scaled from the identical LJTS_MEOS NE fit using source sigma, epsilon and gas constant; no public reduced saturation wrapper is claimed.')
    fixture.with_suffix('.metadata.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')
    recovery_fixture=ROOT/'tests/fixtures/meos_saturation_recovered.csv'
    with recovery_fixture.open('w',newline='') as handle:
        w=csv.writer(handle,delimiter=';',lineterminator='\n');w.writerow(['case','native_temperature','native_pressure','native_liquid_volume','native_vapor_volume','native_dpdt','source_outcome','source_mu_residual','source_liquid_pressure_error','source_vapor_pressure_error']);w.writerows(recovered)
    recovery_fixture.with_suffix('.metadata.json').write_text(json.dumps(dict(identity,cases=len(recovered),fixture_sha256=hashlib.sha256(recovery_fixture.read_bytes()).hexdigest(),generator='tools/rust_port/validate_meos_saturation.py',
        qualification='Native coexistence points independently checked against the exact-source Helmholtz/pressure/caloric kernels, where the original public saturation solver fails or returns an unqualified state. These are not original-solver coordinate parity claims.'),indent=2)+'\n')
    (ROOT/'.reference/meos_saturation_differences.json').write_text(json.dumps(differences,indent=2)+'\n')
    print(json.dumps(dict(attempts=len(records),qualified=len(good),errors=len(failed),censored=len(censored),unqualified=len(unqualified),warnings=len(warnings),differences=len(differences),native_errors=len(native_errors),independently_checked=len(independent),recovered=len(recovered),maximum=maximum),indent=2))
    if differences:print(json.dumps(differences[:12],indent=2));raise SystemExit(1)
    assert not native_errors,native_errors

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');p.add_argument('--reuse-source',action='store_true');a=p.parse_args();validate(a.compiler_directory,a.reuse_source)
