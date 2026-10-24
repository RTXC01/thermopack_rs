"""Analytic MEOS density jets, source-EOS critical points and stable pressure roots.

Reference critical/branch searches use only original mp_pressure outputs, not
native derivative polynomials or native root-search code.
"""
import argparse,csv,hashlib,inspect,json,math,os,subprocess
from collections import defaultdict
from pathlib import Path
from reference import ROOT
from deadline_reference import DeadlineReferenceProcess
from validate_alpha import verify_probe
from validate_meos import MODELS

JET_FIELDS=['pressure','pressure_rho','pressure_rhorho','pressure_rhorhorho','pressure_T','pressure_Trho','pressure_Trhorho']

def case(e,c,t,x,mode):return f'MEOS|{e}:{c}|{t!r}|{x!r}|1|{mode}'

class Source:
    def __init__(self,oracle,e,c,meta):self.oracle,self.e,self.c,self.meta=oracle,e,c,meta
    def pressure(self,t,rhos):
        r=self.oracle.call(dict(mode='MEOS_PRESSURE_GRID',eos=self.e,component=self.c,temperature=t,densities=rhos),timeout=15.)
        if 'rows' not in r or any(v is None for row in r['rows'] for v in row):raise ArithmeticError(str(r))
        return r['rows']
    def jet(self,t,rho):
        h=5e-4*self.meta[2]
        central=rho>=2*h
        offsets=[-2,-1,0,1,2] if central else [0,1,2,3,4]
        first=[1,-8,0,8,-1] if central else [-25,48,-36,16,-3]
        second=[-1,16,-30,16,-1] if central else [35,-104,114,-56,11]
        estimates=[];rounding=[]
        for step in [h,.5*h]:
            rows=self.pressure(t,[rho+i*step for i in offsets]);base=rows[2 if central else 0]
            values=[base[0],base[1],0.,0.,base[2],0.,0.];uncertainty=[0.]*7
            for index,column,weights,power in [(2,1,first,1),(3,1,second,2),(5,2,first,1),(6,2,second,2)]:
                values[index]=sum(w*r[column] for w,r in zip(weights,rows))/(12.*step**power)
                uncertainty[index]=128.*math.ulp(1.)*sum(abs(w*r[column]) for w,r in zip(weights,rows))/(12.*step**power)
            estimates.append(values);rounding.append(uncertainty)
        error=[4.*abs(a-b)+u+v for a,b,u,v in zip(*estimates,*rounding)]
        return estimates[1],error
    def critical(self,t,rho):
        for iteration in range(80):
            q,u=self.jet(t,rho);scale=self.meta[10]*t;norm=math.hypot(q[1]/scale,rho*q[2]/scale)
            if norm<max(5e-10,2.*rho*u[2]/scale):
                if q[0]<=0 or q[3]<=0:raise ArithmeticError(dict(reason='nonordinary stationary pressure point',temperature=t,density=rho,jet=q,uncertainty=u))
                a,b,c,d=t*q[5],rho*q[2],t*q[6],rho*q[3];det=a*d-b*c
                ep=abs(q[1])+2e-11*scale+128.*math.ulp(1.)*scale;ee=abs(q[2])+u[2]+2e-11*scale/rho
                ut=t*(abs(d)*ep+abs(b)*ee)/abs(det);ur=rho*(abs(c)*ep+abs(a)*ee)/abs(det)
                return dict(values=[t,q[0],rho,q[0]/(rho*scale)],jet=q,uncertainty=u,coordinate_uncertainty=[ut,ur],iterations=iteration)
            a,b,c,d=t*q[5],rho*q[2],t*q[6],rho*q[3];det=a*d-b*c
            if not math.isfinite(det) or det==0:break
            st=(-q[1]*d+b*q[2])/det;sr=(-a*q[2]+q[1]*c)/det;length=min(1.,.3/max(abs(st),abs(sr)));accepted=False
            for _ in range(24):
                nt=t*math.exp(length*st);nr=rho*math.exp(length*sr);nq,nu=self.jet(nt,nr)
                nn=math.hypot(nq[1]/(self.meta[10]*nt),nr*nq[2]/(self.meta[10]*nt))
                if nn<norm:t,rho=nt,nr;accepted=True;break
                length*=.5
            if not accepted:break
        raise ArithmeticError(dict(reason='source-pressure critical iteration',temperature=t,density=rho,norm=norm))
    def bisect(self,t,lo,hi,quantity,target=0.):
        a=self.pressure(t,[lo])[0][quantity]-target
        for _ in range(80):
            mid=lo+(hi-lo)*.5
            if mid==lo or mid==hi:break
            b=self.pressure(t,[mid])[0][quantity]-target
            if b==0.:return mid
            if (a<0.)==(b<0.):lo=mid;a=b
            else:hi=mid
        return lo+(hi-lo)*.5
    def branches(self,t,maximum,critical,n):
        rc=self.meta[2]
        grid=[maximum*i/n for i in range(n+1)]+[rc*10.**(-14.+14.*i/128.) for i in range(129)]
        for point in critical:
            rho=point['values'][2];grid.append(rho)
            for k in range(2,11):grid.extend([rho*(1.-10.**-k),rho*(1.+10.**-k)])
        grid=sorted(set(r for r in grid if 0.<=r<=maximum));data=self.pressure(t,grid);spinodals=[]
        for i in range(len(grid)-1):
            a,b=data[i][1],data[i+1][1]
            if (a<0.)!=(b<0.):spinodals.append(self.bisect(t,grid[i],grid[i+1],1))
        unique=[]
        for rho in spinodals:
            if not unique or abs(rho-unique[-1])>1e-9*rc:unique.append(rho)
        bounds=[0.,*unique,maximum];intervals=[]
        for lo,hi in zip(bounds,bounds[1:]):
            if self.pressure(t,[(lo+hi)*.5])[0][1]>0.:
                a,b=self.pressure(t,[lo,hi]);intervals.append([lo,hi,a[0],b[0]])
        return dict(spinodals=unique,intervals=intervals,grid_size=len(grid))
    def roots(self,t,branches,p):
        roots=[]
        for lo,hi,pl,ph in branches['intervals']:
            if pl<=p<=ph:roots.append(self.bisect(t,lo,hi,0,p))
        return roots

def native(cases):
    exe=ROOT/'target/release/examples'/('parity.exe' if os.name=='nt' else 'parity')
    output=subprocess.run([str(exe)],input='\n'.join(cases)+'\n',text=True,capture_output=True,check=True).stdout.splitlines();assert len(output)==len(cases)
    return [line if line.startswith('ERROR') else list(map(float,line.split(','))) for line in output]

def write_fixture(name,fields,rows,identity,**extra):
    path=ROOT/'tests/fixtures'/name
    with path.open('w',newline='') as handle:
        w=csv.writer(handle,delimiter=';',lineterminator='\n');w.writerow(fields);w.writerows(rows)
    path.with_suffix('.metadata.json').write_text(json.dumps(dict(identity,cases=len(rows),fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),generator='tools/rust_port/validate_meos_equilibrium.py',**extra),indent=2,allow_nan=False)+'\n')

def validate(directory,reuse_source):
    identity=json.loads((ROOT/'.reference/metadata.json').read_text());probe=verify_probe()
    worker=hashlib.sha256((ROOT/'tools/rust_port/oracle_worker.py').read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    source_algorithm=hashlib.sha256(inspect.getsource(Source).replace('\r\n','\n').encode()).hexdigest()
    generator=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    cache=ROOT/'.reference/meos_equilibrium_source.json'
    legacy=list(csv.DictReader((ROOT/'tests/fixtures/meos_density.csv').open(),delimiter=';'))
    legacy_hash=hashlib.sha256((ROOT/'tests/fixtures/meos_density.csv').read_bytes()).hexdigest()
    if reuse_source:
        saved=json.loads(cache.read_text());assert saved['reference']==identity and saved['probe']==probe and saved['worker']==worker and saved['source_algorithm_sha256']==source_algorithm and saved['generator_sha256']==generator and saved['legacy_sha256']==legacy_hash;records=saved['records']
    else:
        records=[]
        with DeadlineReferenceProcess(directory) as oracle:
            for e,c in MODELS:
                m=oracle.call(dict(mode='MEOS',eos=e,component=c,temperature=300.,volume=1.))['metadata'];source=Source(oracle,e,c,m);jets=[];critical=[];errors=[]
                for tr in [.45,.7,.95,1.,1.3,2.]:
                    for delta in [0.,.05,.5,.95,1.,1.05,2.,4.]:
                        t=tr*m[0];rho=delta*m[2];q,u=source.jet(t,rho);jets.append(dict(case=case(e,c,t,rho,'PRESSURE_JET'),values=q,uncertainty=u))
                for ratio in [1.,.95,1.05]:
                    k=case(e,c,m[0],ratio*m[2],'CRITICAL')
                    try:critical.append(dict(case=k,**source.critical(m[0],ratio*m[2])))
                    except ArithmeticError as error:errors.append(dict(case=k,error=str(error)))
                primary=next((r for r in critical if float(r['case'].split('|')[3])==m[2]),None)
                automatic=primary or max(critical,key=lambda r:r['values'][0]);automatic=dict(automatic,case=case(e,c,0.,0.,'CRITICAL_AUTO'))
                selected=[r for r in legacy if r['case'].split('|')[1]==f'{e}:{c}'];groups=defaultdict(set)
                for row in selected:
                    f=row['case'].split('|');groups[float(f[2])].add(float(f[3]))
                isotherms=[]
                for t,pressures in sorted(groups.items()):
                    maximum=8.*m[2];coarse=source.branches(t,maximum,critical,2048);fine=source.branches(t,maximum,critical,4096)
                    assert len(coarse['spinodals'])==len(fine['spinodals']),(e,c,t,coarse,fine)
                    assert all(abs(a-b)<2e-7*m[2] for a,b in zip(coarse['spinodals'],fine['spinodals'])),(e,c,t,coarse,fine)
                    roots=[dict(pressure=p,densities=source.roots(t,fine,p)) for p in sorted(pressures)]
                    isotherms.append(dict(temperature=t,maximum_density=maximum,coarse=coarse,reference=fine,roots=roots))
                records.append(dict(eos=e,component=c,metadata=m,jets=jets,critical=critical,critical_errors=errors,automatic=automatic,isotherms=isotherms))
                print(e,c,len(jets),'jets',len(critical),'ordinary critical seeds',len(isotherms),'isotherms',flush=True)
        cache.write_text(json.dumps(dict(reference=identity,probe=probe,worker=worker,source_algorithm_sha256=source_algorithm,generator_sha256=generator,legacy_sha256=legacy_hash,records=records),allow_nan=False)+'\n')
    subprocess.run(['cargo','build','--release','--example','parity','-j','2'],cwd=ROOT,check=True)
    with DeadlineReferenceProcess(directory) as checker:
        compare(records,identity,probe,worker,legacy,checker)

def compare(records,identity,probe,worker,legacy,checker):
    jet_rows=[];critical_rows=[];root_rows=[];spinodal_rows=[];differences=[];unresolved=[];diagnostics=[];legacy_membership=[];selection_differences=[]
    independent_checks=[];roundoff_limited=[]
    for record in records:
        e,c,m=record['eos'],record['component'],record['metadata']
        js=record['jets'];actual=native([r['case'] for r in js])
        for r,q in zip(js,actual):
            if isinstance(q,str):differences.append(dict(case=r['case'],native_error=q));continue
            scales=[m[10]*m[0]*m[2],m[10]*m[0],m[10]*m[0]/m[2],m[10]*m[0]/m[2]**2,m[10]*m[2],m[10],m[10]/m[2]]
            tolerances=[u+2e-8*s+2e-7*abs(v) for u,s,v in zip(r['uncertainty'],scales,r['values'])]
            bad=[i for i,(a,b,tol) in enumerate(zip(r['values'],q,tolerances)) if not math.isfinite(b) or abs(a-b)>tol]
            if bad:differences.append(dict(case=r['case'],quantities=bad,reference=r,native=q,tolerances=tolerances))
            jet_rows.append([r['case'],*r['values'],*tolerances])
        cs=[*record['critical'],record['automatic']];actual=native([r['case'] for r in cs])
        for r,q in zip(cs,actual):
            expected=r['values'];ut,ur=r['coordinate_uncertainty'];floors=[2e-9+8*ut,1e-3+3e-8*abs(expected[1]),1e-6+8*ur,2e-8+8*ur*expected[3]/expected[2]]
            if isinstance(q,str) or any(not math.isfinite(b) or abs(a-b)>tol for a,b,tol in zip(expected,q,floors)):
                differences.append(dict(case=r['case'],reference=r,native=q,tolerances=floors))
            critical_rows.append([r['case'],*expected,*floors])
        if record['critical_errors']:
            rejected=native([r['case'] for r in record['critical_errors']])
            for r,q in zip(record['critical_errors'],rejected):
                if not isinstance(q,str):differences.append(dict(case=r['case'],nonordinary_source_but_native_finite=q))
        source=Source(checker,e,c,m)
        for iso in record['isotherms']:
            t=iso['temperature'];maximum=iso['maximum_density'];pressures=[r['pressure'] for r in iso['roots']]
            k=f'MEOS|{e}:{c}|{t!r}|{maximum!r}|'+','.join(map(repr,pressures))+'|BRANCHES';encoded=native([k])[0]
            if isinstance(encoded,str):differences.append(dict(case=k,native_error=encoded));continue
            count=int(encoded[0]);brackets=[encoded[1+5*i:6+5*i] for i in range(count)];offset=1+5*count
            for rho in iso['reference']['spinodals']:
                matched=[b for b in brackets if b[0]-2e-7*m[2]<=rho<=b[1]+2e-7*m[2]]
                if not matched:differences.append(dict(case=k,missing_spinodal=rho,brackets=brackets))
                spinodal_rows.append([case(e,c,t,maximum,'SPINODALS'),rho])
            for b in brackets:
                if b[4]==1. and sum(b[0]-2e-7*m[2]<=rho<=b[1]+2e-7*m[2] for rho in iso['reference']['spinodals'])!=1:
                    differences.append(dict(case=k,unmatched_isolated_spinodal=b,reference=iso['reference']['spinodals']))
            for row in iso['roots']:
                p=row['pressure'];expected=row['densities'];n=int(encoded[offset]);offset+=1
                root_case=case(e,c,t,maximum,f'ROOTS:{p!r}')
                if n<0:
                    phases=encoded[offset:offset+3];offset+=3
                    boxes=[b for b in brackets if b[2]<=p<=b[3]]
                    limited=[]
                    for rho,values in zip(expected,source.pressure(t,expected)):
                        tolerance=2e-8*p+256.*math.ulp(1.)*m[10]*t*rho
                        if abs(values[0]-p)>tolerance:
                            neighbors=source.pressure(t,[math.nextafter(rho,0.),rho,math.nextafter(rho,math.inf)])
                            limited.append(dict(density=rho,source_pressure=values[0],pressure_error=values[0]-p,strict_tolerance=tolerance,
                                source_pressure_quantum=math.ulp(rho)*abs(values[1]),neighbor_pressures=[r[0] for r in neighbors]))
                    if not boxes and not(n==-2 and limited):differences.append(dict(case=root_case,unexplained_native_error=True,raw_source_candidates=limited))
                    resolved_volumes=[]
                    for index,value in enumerate(phases):
                        if not math.isfinite(value):resolved_volumes.append(None);continue
                        rho=1./value;target=expected[0 if index==0 else -1];q=source.pressure(t,[rho])[0]
                        if index==2 or abs(rho-target)>2e-8*m[2]+2e-7*abs(rho) or q[1]<=0. or abs(q[0]-p)>4e-8*p+512.*math.ulp(1.)*m[10]*t*rho:
                            differences.append(dict(case=root_case,unqualified_partial_phase=index,volume=value,source=q))
                        resolved_volumes.append(1./target)
                    unresolved.append(dict(case=root_case,reference_densities=expected,brackets=boxes,native_reason='density_precision' if n==-2 else 'spinodal_window',unqualified_source_candidates=limited,resolved_extreme_volumes=resolved_volumes));continue
                q=encoded[offset:offset+n];allowances=encoded[offset+n:offset+2*n];volumes=encoded[offset+2*n:offset+2*n+3];offset+=2*n+3
                if n!=len(expected) or any(abs(a-b)>2e-8*m[2]+2e-7*max(abs(a),abs(b)) for a,b in zip(expected,q)):
                    differences.append(dict(case=root_case,reference_densities=expected,native_densities=q))
                scores=[]
                checked=source.pressure(t,q)
                for rho,allowance,pressure_values in zip(q,allowances,checked):
                    raw=checker.call(dict(mode='MEOS',eos=e,component=c,temperature=t,volume=1./rho))
                    assert 'values' in raw,raw
                    alpha=raw['values'][6];score=alpha+math.log(rho/m[2])+p/(rho*m[10]*t);scores.append(score)
                    actual_p,derivative,_=pressure_values;ordinary=4e-12*p+64.*math.ulp(1.)*(m[10]*t*rho+p)
                    entry=dict(case=root_case,density=rho,source_pressure=actual_p,source_pressure_rho=derivative,source_reduced_gibbs=score,pressure_error=actual_p-p,native_pressure_allowance=allowance)
                    independent_checks.append(entry)
                    cap=2e-8*p+256.*math.ulp(1.)*m[10]*t*rho
                    if derivative<=0. or abs(actual_p-p)>cap+min(4.*allowance,cap):differences.append(dict(source_root_check=entry))
                    if abs(actual_p-p)>ordinary:roundoff_limited.append(entry)
                if q:
                    if not all(map(math.isfinite,volumes)) or abs(volumes[0]*q[0]-1.)>1e-9 or abs(volumes[1]*q[-1]-1.)>1e-9:
                        differences.append(dict(case=root_case,incorrect_extreme_selection=volumes,roots=q))
                    else:
                        choice=min(range(n),key=lambda i:abs(q[i]*volumes[2]-1.))
                        if abs(q[choice]*volumes[2]-1.)>1e-9 or scores[choice]>min(scores)+1e-10+512.*math.ulp(1.)*max(map(abs,scores)):
                            differences.append(dict(case=root_case,incorrect_gibbs_selection=volumes[2],roots=q,source_scores=scores))
                source_minimum=1./expected[min(range(n),key=lambda i:scores[i])] if expected else None
                root_rows.append([root_case,','.join(map(repr,expected)),1./expected[0] if expected else None,1./expected[-1] if expected else None,source_minimum])
                for legacy_row in legacy:
                    f=legacy_row['case'].split('|')
                    if f[1]!=f'{e}:{c}' or float(f[2])!=t or float(f[3])!=p:continue
                    rho=float(legacy_row['rho']);dp=float(legacy_row['pressure_rho'])
                    tolerance=2e-8*m[2]+2e-7*abs(rho)+1e-13*abs(p/dp) if dp else 2e-8*m[2]+2e-7*abs(rho)
                    if dp>0. and not any(abs(rho-a)<=tolerance for a in q):differences.append(dict(case=legacy_row['case'],missing_legacy_density=rho,native=q))
                    legacy_membership.append(legacy_row['case'])
                    requested=int(f[4]);selected=volumes[1 if requested==1 else 0] if requested in [1,2] and q else None
                    if selected is not None and abs(1./selected-rho)>tolerance:selection_differences.append(dict(case=legacy_row['case'],legacy_density=rho,selected_density=1./selected))
            assert offset==len(encoded),(offset,len(encoded))
            diagnostics.append(dict(eos=e,component=c,temperature=t,brackets=brackets,source=iso['reference']))
        print(e,c,'native checks complete',flush=True)
    common=dict(identity,probe=probe,oracle_sha256=worker,generator_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest(),legacy_density_fixture_sha256=hashlib.sha256((ROOT/'tests/fixtures/meos_density.csv').read_bytes()).hexdigest())
    write_fixture('meos_pressure_jet.csv',['case',*JET_FIELDS,*['tolerance_'+s for s in JET_FIELDS]],jet_rows,common,
        method='Original mp_pressure P/P_rho/P_T; independent five-point density differences at two step sizes for higher/mixed derivatives. Stored per-field tolerances include step-size disagreement, rounding amplification and explicit small kernel floors.')
    write_fixture('meos_critical.csv',['case','temperature','pressure','density','compressibility','tolerance_temperature','tolerance_pressure','tolerance_density','tolerance_compressibility'],critical_rows,common,
        source_attempts=sum(len(r['critical'])+len(r['critical_errors']) for r in records),source_errors=[e for r in records for e in r['critical_errors']],
        references=[c for r in records for c in [*r['critical'],r['automatic']]],method='Independent damped critical solve through original pressure and first derivatives; density columns are five-point source differences. These are source-EOS-backed equations, not direct original critical-solver trajectory parity.')
    write_fixture('meos_pressure_branches.csv',['case','stable_densities','vapor_volume','liquid_volume','minimum_gibbs_volume'],root_rows,common,
        isotherms=diagnostics,unresolved_pressure_windows=unresolved,legacy_membership_count=len(legacy_membership),legacy_selection_differences=selection_differences,
        independent_root_checks=independent_checks,roundoff_limited_pressure_checks=roundoff_limited,
        method='Independent source-pressure/stiffness bisections on doubled 2048/4096 grids, log-dilute points and independently source-solved critical-density neighborhoods. Every finite existing source density iterate in a resolved pressure query is checked for root-set membership. Native unresolved windows are retained, not silently omitted.')
    write_fixture('meos_spinodals.csv',['case','density'],spinodal_rows,common,method='All source pressure-stiffness sign changes on the independent doubled-grid/critical-focused search; native isolated boxes must match, and unresolved boxes are reported separately.')
    write_fixture('meos_pressure_unresolved.csv',['case','reason','maximum_source_relative_pressure_error','maximum_source_pressure_quantum','resolved_vapor_volume','resolved_liquid_volume','resolved_minimum_gibbs_volume'],
        [[r['case'],r['native_reason'],max([abs(x['pressure_error'])/float(r['case'].split(':')[-1]) for x in r['unqualified_source_candidates']] or [0.]),max([x['source_pressure_quantum'] for x in r['unqualified_source_candidates']] or [0.]),*r['resolved_extreme_volumes']] for r in unresolved],common,
        evidence=unresolved,method='Strict whole-root-set queries reject pressure windows or density-precision-limited branches. Raw independent source candidates, source pressure errors, neighbor pressures and density quanta remain explicit; other source phase iterates are not declared invalid.')
    (ROOT/'.reference/meos_equilibrium_differences.json').write_text(json.dumps(differences,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(jets=len(jet_rows),critical=len(critical_rows),resolved_pressure_states=len(root_rows),unresolved_pressure_states=len(unresolved),spinodals=len(spinodal_rows),legacy_membership=len(legacy_membership),legacy_selection_differences=len(selection_differences),independent_roots=len(independent_checks),roundoff_limited=len(roundoff_limited),differences=len(differences)),indent=2))
    if differences:print(json.dumps(differences[:10],indent=2));raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compiler-directory');p.add_argument('--reuse-source',action='store_true');a=p.parse_args();validate(a.compiler_directory,a.reuse_source)
