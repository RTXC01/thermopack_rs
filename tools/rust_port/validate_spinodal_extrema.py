"""Pure spinodal extrema selection and raw coordinates from this exact fork."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess

RTOL=2e-6
ATOL=[2e-5,1e-11,2.]

def cases():
    for eos in ['PR','SRK','RK','VDW']:
        for name in ['C1','CO2','C3','N2','H2','H2O','MEOH','ETOH','MEG']:
            for p in [1e5,1e6]:
                for prop in range(1,7):yield eos,name,'classic',p,prop
    for eos in ['PR','SRK']:
        for name in ['C1','CO2','H2O']:
            for p in [1e5,1e6]:
                for prop in range(1,7):yield eos,name,'translated',p,prop

def validate(directory):
    good=[];failed=[];censored=[];unqualified=[];warnings=[]
    with DeadlineReferenceProcess(directory) as oracle:
        for eos,name,config,p,prop in cases():
            label=f'SPINPROP|{eos}|{name}|{config}|EXTREMA:{prop}|{p}'
            request=dict(mode='PURE_SPINODAL_EXTREMA',eos=eos,components=[name],feed=[1.],volume_shift=config=='translated',pressure=p,property=prop)
            result=oracle.call(request,timeout=8.)
            record=dict(case=label,request=request,result=result)
            if 'reference_timeout' in result:censored.append(record);continue
            if 'values' not in result or result.get('status',0):failed.append(record);continue
            values=result['values']
            if len(values)!=6 or not all(map(math.isfinite,values)) or min(values[0],values[1],values[3],values[4])<=0:
                unqualified.append(record);continue
            good.append(record)
            for i,(actual,lam,value) in enumerate(result['diagnostics']):
                state=values[3*i:3*i+3]
                reasons=[]
                if abs(actual-state[2])>2.+1e-7*max(abs(actual),abs(state[2])):reasons.append('reported pressure differs from EOS coordinates')
                if abs(lam)>1e-6:reasons.append('reported coordinates are off spinodal')
                if reasons:warnings.append(dict(case=label,endpoint=i,reasons=reasons,coordinates=state,actual_pressure=actual,eigenvalue=lam,property_value=value))
    subprocess.run(['cargo','build','--release','--example','parity','-j','2'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    output=subprocess.run([str(exe)],input='\n'.join(r['case'] for r in good)+'\n',capture_output=True,text=True,check=True).stdout.splitlines()
    assert len(output)==len(good)
    rows=[];differences=[];maxima=[0.]*6;ambiguous=[];permutations=[]
    for r,out in zip(good,output):
        result=r['result'];original=result['values'];diagnostics=result['diagnostics']
        tie=r['request']['property']==5 and abs(diagnostics[0][2]-diagnostics[1][2])<=2.+1e-7*max(abs(diagnostics[0][2]),abs(diagnostics[1][2]))
        alternative=result.get('alternative_pressures',[0.,0.])
        rows.append([r['case'],*original,*[v for d in diagnostics for v in d],int(tie),*alternative])
        if tie:ambiguous.append(r['case'])
        if out.startswith('ERROR'):differences.append(dict(case=r['case'],rust=out));continue
        native=list(map(float,out.split(',')));assert len(native)==7
        if native[6]!=int(tie):differences.append(dict(case=r['case'],ambiguity=[int(tie),native[6]]))
        expected=original
        if tie:
            swapped=[original[3],original[4],alternative[0],original[0],original[1],alternative[1]]
            def discrepancy(values):return max(abs(a-b)/(ATOL[j%3]+RTOL*max(abs(a),abs(b))) for j,(a,b) in enumerate(zip(values,native)))
            if discrepancy(swapped)<discrepancy(original):expected=swapped;permutations.append(r['case'])
        bad=[]
        for j,(a,b) in enumerate(zip(expected,native)):
            maxima[j]=max(maxima[j],abs(a-b)/max(abs(a),abs(b),ATOL[j%3]))
            if not math.isfinite(b) or abs(a-b)>ATOL[j%3]+RTOL*max(abs(a),abs(b)):bad.append([j,a,b])
        if bad:differences.append(dict(case=r['case'],quantities=bad))
    fixture=ROOT/'tests/fixtures/spinodal_extrema.csv'
    with fixture.open('w',newline='') as f:
        w=csv.writer(f,delimiter=';',lineterminator='\n')
        w.writerow(['case','T_min','V_min','P_min','T_max','V_max','P_max','EOS_P_min','lambda_min','property_min','EOS_P_max','lambda_max','property_max','ordering_ambiguous','alternative_P_min','alternative_P_max'])
        w.writerows(rows)
    metadata=json.loads((ROOT/'.reference/metadata.json').read_text())
    metadata.update(cases=len(rows),attempts=sum(1 for _ in cases()),generator='tools/rust_port/validate_spinodal_extrema.py',fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),
        relative_tolerance=RTOL,absolute_tolerances=ATOL,reference_errors=failed,censored_calls=censored,unqualified_returns=unqualified,
        source_coordinate_warnings=warnings,ambiguous_pressure_orderings=ambiguous,primary_order_differences=permutations,worst_relative_differences=maxima)
    fixture.with_suffix('.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (ROOT/'.reference/extrema_differences.json').write_text(json.dumps(differences,indent=2)+'\n')
    print(json.dumps(dict(attempts=metadata['attempts'],finite_source_returns=len(good),unambiguous_returns=len(good)-len(ambiguous),ambiguous_candidates=len(ambiguous),primary_order_differences=len(permutations),failed=len(failed),censored=len(censored),unqualified=len(unqualified),warnings=len(warnings),differences=len(differences),worst_relative_differences=maxima),indent=2))
    if differences:print(json.dumps(differences[:10],indent=2));raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');validate(p.parse_args().compiler_directory)
