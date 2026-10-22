"""Exact-source TP fugacity derivatives for activity-based cubic mixing."""
import argparse,csv,hashlib,json,math,os,subprocess
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess

def models():
    for eos in ['PR','SRK']:
        for names,z in [('C1,C3',[.3,.7]),('C1,C2,C3',[.2,.3,.5]),('C1,C2,C3',[.5,.5,0.])]:yield eos,names,z,'UNIFAC','DEFAULT',False
        yield eos,'C1,C3',[.3,.7],'UNIFAC','DEFAULT',True
    for rule in ['HV1','HV2']:
        yield 'SRK','C1,H2O',[.4,.6],rule,'DEFAULT',False
        yield 'PR','C1,H2O',[.4,.6],rule,'DEFAULT',False
    yield 'PR','H2O,MEG',[.3,.7],'HV1','Hemmingsen2011',False
    yield 'PR','H2S,C3',[.4,.6],'NRTL','Dicko2012',False

def validate(directory):
    good=[];failed=[];censored=[];unqualified=[];warnings=[]
    with DeadlineReferenceProcess(directory) as oracle:
        for eos,names,z,mix,reference,translated in models():
            for t in [250.,300.,400.]:
                for p in [1e5,1e6,1e7]:
                    for phase in ['L','V']:
                        label=f"GEDER|{eos}|{names}|{'translated' if translated else 'classic'}|"+','.join(map(repr,z))+f'|{t}:{p}:{phase}:{mix}:{reference}'
                        request=dict(mode='GE_DERIVATIVES',eos=eos,components=names.split(','),feed=z,mixing=mix,parameter_reference=reference,volume_shift=translated,temperature=t,pressure=p,phase=phase)
                        result=oracle.call(request,timeout=8.);record=dict(case=label,request=request,result=result)
                        if 'reference_timeout' in result:censored.append(record);continue
                        if 'values' not in result:failed.append(record);continue
                        if not all(map(math.isfinite,result['values'])):unqualified.append(record);continue
                        if abs(result['pressure_error'])>2.+1e-7*p:warnings.append(record)
                        good.append(record)
    subprocess.run(['cargo','build','--release','--example','parity','-j','2'],cwd=ROOT,check=True)
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    output=subprocess.run([str(exe)],input='\n'.join(r['case'] for r in good)+'\n',capture_output=True,text=True,check=True).stdout.splitlines();assert len(output)==len(good)
    rows=[];differences=[];maxima=[0.]*3
    for r,out in zip(good,output):
        expected=r['result']['values'];rows.append([r['case'],*expected]);nc=len(r['request']['feed'])
        if out.startswith('ERROR'):differences.append(dict(case=r['case'],rust=out));continue
        native=list(map(float,out.split(',')));assert len(native)==len(expected)
        bad=[]
        for i,(a,b) in enumerate(zip(expected,native)):
            group=0 if i<nc else 1 if i<2*nc else 2;atol=[1e-9,1e-11,1e-8][group]
            maxima[group]=max(maxima[group],abs(a-b)/max(abs(a),abs(b),1.))
            if not math.isfinite(b) or abs(a-b)>atol+2e-7*max(abs(a),abs(b)):bad.append([i,a,b])
        if bad:differences.append(dict(case=r['case'],quantities=bad))
    path=ROOT/'tests/fixtures/ge_fugacity_derivatives.csv'
    with path.open('w',newline='') as f:w=csv.writer(f,delimiter=';',lineterminator='\n');w.writerow(['case','dt[n],dp[n],dn[n*n]']);w.writerows(rows)
    meta=json.loads((ROOT/'.reference/metadata.json').read_text());meta.update(cases=len(rows),generator='tools/rust_port/validate_ge_derivatives.py',fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        relative_tolerance=2e-7,absolute_tolerances=[1e-9,1e-11,1e-8],reference_errors=failed,censored_calls=censored,unqualified_returns=unqualified,source_pressure_warnings=warnings,worst_scaled_errors=maxima)
    path.with_suffix('.metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    (ROOT/'.reference/ge-derivative-differences.json').write_text(json.dumps(differences,indent=2)+'\n')
    print(json.dumps(dict(qualified=len(rows),failed=len(failed),censored=len(censored),unqualified=len(unqualified),warnings=len(warnings),differences=len(differences),worst_scaled_errors=maxima),indent=2))
    if differences:print(json.dumps(differences[:8],indent=2));raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');validate(p.parse_args().compiler_directory)
