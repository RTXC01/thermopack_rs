"""Original-source qualified pure MEOS flash benchmarking, with opt-in timing."""
import argparse
import csv
import hashlib
import json
import math
import os
import platform
import statistics
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from deadline_reference import DeadlineReferenceProcess

ROOT = Path(__file__).resolve().parents[2]
PROFILES = ['NIST_MEOS:C3', 'NIST_MEOS:R134A', 'LJ_MEOS:KR', 'LJ_MEOS:AR',
            'LJTS_MEOS:KR', 'LJTS_REDUCED:LJTS']
LANES = {'TP':0, 'PH_cold':1, 'PH_guessed':2, 'PS_cold':3, 'PS_guessed':4,
         'prepared_scalar_TV_HS':5}
ITERATIONS = 1000


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hashes(paths):
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
            for p in paths}


def parse_rows(text, prefix, fields):
    result = {}
    for line in text.splitlines():
        if not line.startswith(prefix):
            continue
        key, *parts = line.split(',')
        assert key not in result, ('duplicate output key', key)
        values = list(map(float, parts))
        assert len(values) == fields and all(map(math.isfinite, values)), (key, values)
        result[key] = values
    return result


def valid_result(q, property_lane=False):
    if len(q) != 8 or not all(map(math.isfinite, q)) or min(q[:2]) <= 0.:
        return False
    phase, beta, vl, vv = q[2:6]
    if property_lane:
        return phase == 4. and beta == -1.
    if phase == 0.:
        return 0. <= beta <= 1. and 0. < vl < vv
    if phase == 1.:
        return beta == 0. and vl > 0. and vv == 0.
    if phase == 2.:
        return beta == 1. and vl == 0. and vv > 0.
    return phase == 4. and beta in (-1., 0., 1.) and vl > 0. and vv == 0.


def endpoint_coordinates(q):
    assert valid_result(q)
    t, p, phase, beta, vl, vv, h, s = q
    if phase == 0.:
        assert min(beta, 1.-beta) <= 2e-6, ('material split cannot be hidden as an endpoint', q)
        volume = (1.-beta)*vl+beta*vv
    else:
        volume = vv if phase == 2. else vl
    return [t, p, 4., 0., volume, 0., h, s]


def expected_checksum(checks, keys, iterations=ITERATIONS):
    assert keys and iterations > 0
    n = len(keys)
    return [math.fsum(checks[key][j]*(iterations//n+(i < iterations % n))
                      for i, key in enumerate(keys)) for j in range(8)]


def verify_checksum(actual, expected, magnitude=None, iterations=ITERATIONS):
    assert len(actual) == len(expected) == 8
    magnitude = [0.]*8 if magnitude is None else magnitude
    assert len(magnitude) == 8 and iterations > 0
    for j, (got, want) in enumerate(zip(actual, expected)):
        # Bound naive repeated summation by absolute term magnitudes too,
        # so cancellation is not confused with a changed workload.
        tolerance = 1e-11+2e-12*max(abs(got), abs(want))+4.*math.ulp(1.)*iterations*magnitude[j]
        assert math.isfinite(got) and abs(got-want) <= tolerance, (j, got, want, tolerance)


def fingerprints(paths):
    return {str(p.resolve()):digest(p) for p in paths}


def parse_args(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--compiler', required=True)
    ap.add_argument('--source-build-dir')
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument('--preflight-only', action='store_true')
    mode.add_argument('--measure', action='store_true')
    ap.add_argument('--repeats', type=int, default=5)
    ap.add_argument('--output', default=str(ROOT/'docs/rust/meos-flash-benchmark.json'))
    args = ap.parse_args(argv)
    if args.repeats < 1:
        ap.error('--repeats must be positive')
    return args


def main():
    args = parse_args()
    git_head = subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    if args.measure:
        output = Path(args.output)
        if output.exists():
            raise SystemExit('Existing timing reports are immutable; choose a fresh --output path')
        output.parent.mkdir(parents=True, exist_ok=True)
    fixture = ROOT/'tests/fixtures/meos_flash.csv'
    meta = json.loads(fixture.with_suffix('.metadata.json').read_text(encoding='utf-8'))
    assert digest(fixture) == meta['fixture_sha256']
    coexist = ROOT/'tests/fixtures/meos_flash_coexistence_tp.csv'
    recovered = ROOT/'tests/fixtures/meos_flash_recovered.csv'
    metadata_fixture = ROOT/'tests/fixtures/meos_metadata.csv'
    def read_rows(path):
        with path.open(encoding='utf-8', newline='') as f:
            return list(csv.DictReader(f, delimiter=';'))
    metadata = {r['case'].split('|')[1]:{k:float(v) for k,v in r.items() if k != 'case'}
                for r in read_rows(metadata_fixture)}
    ne = metadata['LJTS_MEOS:NE']
    # The original worker explicitly scales its default physical neon model.
    # Fail on parameter drift instead of silently benchmarking different units.
    assert (ne['eps_div_k'], ne['sigma'], ne['R']) == (39.83, 2.8e-10, 8.31446261815324)
    gaps_path = ROOT/'tools/rust_port/meos_flash_known_gaps.json'
    gap_keys = json.loads(gaps_path.read_text(encoding='utf-8'))['case_keys']
    gap_cells = {tuple(key.split('|')[:5]) for key in gap_keys}
    audit = dict(undefined_pure_coexistence_tp=sum(r['case'].split('|')[1] in PROFILES for r in read_rows(coexist)),
                 source_recovery_rows_outside_primary_fixture=sum(r['case'].split('|')[1] in PROFILES for r in read_rows(recovered)),
                 declared_coordinate_gap_cases=len(gap_keys), gap_grid_intersections=0,
                 exact_caloric_endpoint_rows=0, nearcritical_caloric_rows=0)
    rows = []
    for row in read_rows(fixture):
        case = row['case'].split('|'); profile = case[1]; control = case[5].split(':')
        if profile not in PROFILES:
            continue
        mode = control[1]
        assert mode in ('TP', 'PH', 'PS')
        vals = [float(row[k]) for k in ('T','P','quality','liquid_or_overall_volume','vapor_or_overall_volume','H','S')]
        assert all(map(math.isfinite, vals))
        if mode == 'TP' and vals[3] != vals[4]:
            audit['undefined_pure_coexistence_tp'] += 1
            continue
        md = metadata[profile]
        tr = float(case[2])/md['Tc']
        pr = float(case[3])/(md['R']*md['Tc']*md['rhoc'])
        eos, component = profile.split(':')
        if any(c[0:3] == (eos, component, mode) and abs(tr-float(c[3])) <= 2e-5
               and abs(pr-float(c[4])) <= 2e-5 for c in gap_cells):
            audit['gap_grid_intersections'] += 1
        kind = {'TP':0,'PH':1,'PS':2}[mode]
        has_guess = mode != 'TP' and control[3] != 'none'
        guess = float(control[3]) if has_guess else 0.
        lane = {'TP':0,'PH':1,'PS':3}[mode]+int(has_guess)
        tmin, tmax = map(float, control[-2:])
        expected = [vals[0], vals[1], 0., vals[2], vals[3], vals[4], vals[5], vals[6]]
        r = [PROFILES.index(profile), lane, kind,
             float(case[2]) if mode == 'TP' else guess, float(case[3]),
             0. if mode == 'TP' else float(control[2]), guess, int(has_guess), tmin, tmax, *expected]
        rows.append(r)
        if mode != 'TP':
            audit['exact_caloric_endpoint_rows'] += int(vals[3] == vals[4])
            audit['nearcritical_caloric_rows'] += int(abs(tr-1.) < .002)
        prop = r.copy(); prop[1] = prop[2] = 5; prop[3:5] = vals[:2]
        rows.append(prop)
    assert rows
    rows.sort(key=lambda r:(r[0],r[1],r[3],r[4],r[7]))
    out = ROOT/'.reference'; out.mkdir(exist_ok=True)
    source_build = Path(args.source_build_dir).resolve() if args.source_build_dir else out
    manifest = out/'meos-flash-benchmark.txt'
    manifest.write_text(str(len(rows))+'\n'+'\n'.join(' '.join(format(x,'.17g') for x in r) for r in rows)+'\n', encoding='utf-8', newline='\n')
    source = ROOT/'tools/rust_port/benchmark_meos_flash_reference.f90'
    bench = ROOT/'rust/benches/meos_flash.rs'
    sources = [source, bench, Path(__file__), ROOT/'tools/rust_port/reference_probe.f90', fixture,
        fixture.with_suffix('.metadata.json'), metadata_fixture, coexist, coexist.with_suffix('.metadata.json'),
        recovered, ROOT/'tests/fixtures/meos_saturation.csv', ROOT/'tests/fixtures/meos_saturation_recovered.csv', gaps_path,
        ROOT/'tools/rust_port/oracle_worker.py', ROOT/'tools/rust_port/reference.py', ROOT/'tools/rust_port/deadline_reference.py']
    sourcehash = hashes(sources); rustpaths = sorted((ROOT/'rust').rglob('*.rs'))
    rusthash = hashes(rustpaths)
    cargo_candidates = [ROOT/'Cargo.toml', ROOT/'Cargo.lock', ROOT/'.cargo/config', ROOT/'.cargo/config.toml',
                        ROOT/'rust-toolchain', ROOT/'rust-toolchain.toml']
    cargopaths = [p for p in cargo_candidates if p.exists()]; cargohash = hashes(cargopaths)
    identity_path = out/'metadata.json'; probe_identity_path = out/'probe.metadata.json'
    identity = json.loads(identity_path.read_text(encoding='utf-8'))
    probe_identity = json.loads(probe_identity_path.read_text(encoding='utf-8'))
    library = out/('thermopack.dll' if os.name == 'nt' else 'libthermopack.so')
    probe = out/('reference_probe.dll' if os.name == 'nt' else 'libreference_probe.so')
    assert identity['library_sha256'] == meta['reference']['library_sha256'] == digest(library)
    assert probe_identity['library_sha256'] == digest(library)
    assert probe_identity['probe_sha256'] == digest(probe)
    compiler = Path(args.compiler).resolve(); rust_bin = Path.home()/'.cargo/bin'
    cargo = rust_bin/('cargo.exe' if os.name == 'nt' else 'cargo')
    rustc = rust_bin/('rustc.exe' if os.name == 'nt' else 'rustc')
    env = dict(os.environ); env['PATH'] = str(compiler.parent)+os.pathsep+str(rust_bin)+os.pathsep+env['PATH']
    flags = ['-O3','-cpp','-fdefault-real-8','-fdefault-double-8','-frecursive','-fopenmp',
             '-fallow-argument-mismatch','-ffree-line-length-none']
    source_exe = out/('benchmark_meos_flash.exe' if os.name == 'nt' else 'benchmark_meos_flash')
    object_manifest = source_build/'objects.rsp'
    object_paths = []
    for line in object_manifest.read_text(encoding='utf-8').splitlines():
        if line.strip():
            p = Path(line.strip().strip('"'))
            object_paths.append(p.resolve() if p.is_absolute() else (source_build/p).resolve())
    assert len(object_paths) == len(set(object_paths)) and object_paths
    module_paths = sorted(set(source_build.glob('*.mod')) | set(out.glob('*.mod')))
    objects_before = fingerprints(object_paths+module_paths+[object_manifest])
    subprocess.run([str(compiler),*flags,'-I',str(out),'-I',str(source_build),str(source),
                    '@'+str(object_manifest),'-o',str(source_exe)], cwd=out, env=env, check=True)
    build = subprocess.run([str(cargo),'bench','--bench','meos_flash','-j','1','--no-run','--message-format=json'],
                           cwd=ROOT, env=env, text=True, capture_output=True, check=True)
    artifacts = [json.loads(line) for line in build.stdout.splitlines() if line.startswith('{')]
    native_exe = Path(next(r['executable'] for r in artifacts if r.get('reason') == 'compiler-artifact'
                          and r['target']['name'] == 'meos_flash' and r.get('executable')))
    binaries = [source_exe, native_exe, library, probe, identity_path, probe_identity_path, manifest]
    frozen = fingerprints(binaries+object_paths+module_paths+[object_manifest])
    assert all(frozen[k] == v for k,v in objects_before.items())
    ref_text = subprocess.check_output([str(source_exe),'check',str(manifest)], env=env, text=True)
    native_text = subprocess.check_output([str(native_exe),str(manifest),'--check'], env=env, text=True)
    ref = parse_rows(ref_text, 'check_', 8); native = parse_rows(native_text, 'check_', 8)
    statuses = parse_rows(ref_text, 'status_', 1)
    assert ref.keys() == native.keys() and len(ref) == len(rows)
    assert {k.removeprefix('status_') for k in statuses} == {k.removeprefix('check_') for k in ref}
    groups = {}; row_by_key = {}; maxima = [0.]*8; endpoint_pairs = []; warnings = []
    for r in rows:
        group = int(r[0]), int(r[1]); selected = groups.setdefault(group, [])
        key = f'check_{group[0]}_{group[1]}_{len(selected)}'; selected.append(key); row_by_key[key] = r
        x, y = ref[key][:], native[key][:]; prop = group[1] == 5
        assert valid_result(x, prop) and valid_result(y, prop), (key, x, y)
        status = statuses[key.replace('check_', 'status_', 1)][0]
        assert status in (-1., 0.), (key, status)
        if status: warnings.append(dict(case=key, status=status))
        if not prop:
            for point in (x, y):
                assert abs(point[1]-r[4]) <= 512.*math.ulp(1.)*max(abs(r[4]), 1e-30), (key, 'requested pressure changed', point[1], r[4])
            if r[2] in (1,2):
                gas = metadata[PROFILES[group[0]]]['R']
                scale = gas*x[0] if r[2] == 1 else gas
                assert abs((x[6] if r[2] == 1 else x[7])-r[5]) <= 2e-7*scale, (key, 'source caloric specification')
                native_scale = gas*y[0] if r[2] == 1 else gas
                assert abs((y[6] if r[2] == 1 else y[7])-r[5]) <= 3e-10*native_scale, (key, 'native caloric specification')
        endpoint = group[1] in (1,2,3,4) and r[14] == r[15]
        if endpoint:
            endpoint_pairs.append(dict(case=key, source=x[:], native=y[:]))
            x, y = endpoint_coordinates(x), endpoint_coordinates(y)
        elif not prop and x[2] != 0. and y[2] != 0.:
            x[3] = y[3] = 0.  # one-phase quality is undefined
        floors = [2e-7,2e-5,0.,2e-7,2e-5,2e-7,2e-5,2e-7] if prop else [2e-7,2e-5,0.,2e-7,1e-12,1e-12,2e-5,2e-7]
        for j,(a,b) in enumerate(zip(x,y)):
            assert abs(a-b) <= floors[j]+2e-6*max(abs(a),abs(b)), (key,j,a,b)
            maxima[j] = max(maxima[j], abs(a-b))
        if not prop:
            # Canonical fixture columns hold bulk V when the source is homogeneous.
            saved = r[10:18]
            actual = endpoint_coordinates(ref[key]) if endpoint else ref[key]
            if endpoint:
                expected_fields = (0,1,4,6,7)
            elif actual[2] == 0.:
                expected_fields = (0,1,4,5,6,7)
            else:
                expected_fields = (0,1,5,6,7) if actual[2] == 2. else (0,1,4,6,7)
            for j in expected_fields:
                assert abs(actual[j]-saved[j]) <= floors[j]+2e-6*max(abs(actual[j]),abs(saved[j])), (key,'source fixture',j,actual[j],saved[j])
    kernel_cache = {}; kernel_checks = []
    with DeadlineReferenceProcess(str(compiler.parent)) as oracle:
        for (profile, lane), keys in groups.items():
            model = PROFILES[profile]; eos, component = model.split(':'); gas = metadata[model]['R']
            for key in keys:
                r = row_by_key[key]
                for label, q in [('source',ref[key]), ('native',native[key])]:
                    t,p,phase,beta,vl,vv,h,entropy = q
                    if lane == 5:
                        volumes, weights = r[14:16], None
                    elif phase == 0.:
                        volumes, weights = [vl,vv], [1.-beta,beta]
                    else:
                        volumes, weights = [vv if phase == 2. else vl], [1.]
                    states = []
                    for volume in volumes:
                        cache_key = model,t,volume
                        if cache_key not in kernel_cache:
                            result = oracle.call(dict(mode='MEOS',eos=eos,component=component,
                                temperature=t,volume=volume),timeout=4.)
                            assert 'values' in result, (key,label,result)
                            st = result['values']
                            kernel_cache[cache_key] = [st[j] for j in (6,12,14,15,16,17,20)]
                        state = kernel_cache[cache_key]
                        assert all(v is not None and math.isfinite(v) for v in state)
                        assert state[2] < 0. and state[6] > 0., (key,label,'unstable phase',state)
                        states.append(state)
                        if lane != 5:
                            assert abs(state[1]-p) <= 2e-8*p+512.*math.ulp(1.)*gas*t/volume, (key,label,'EOS pressure',state[1],p)
                    if lane == 5:
                        wanted = [states[0][5],states[0][4],states[1][5],states[1][4]]
                        actual = q[4:8]; scales = [gas*t,gas,gas*t,gas]
                    else:
                        wanted = [sum(w*st[5] for w,st in zip(weights,states)),sum(w*st[4] for w,st in zip(weights,states))]
                        actual = [h,entropy]; scales = [gas*t,gas]
                        if phase == 0.:
                            mu = math.log(vv/vl)+states[0][0]-states[1][0]+p*(vl-vv)/(gas*t)
                            assert abs(mu) <= 2e-8, (key,label,'chemical equilibrium',mu)
                    residuals = [(a-b)/scale for a,b,scale in zip(actual,wanted,scales)]
                    assert max(map(abs,residuals)) <= 2e-8, (key,label,'original kernel H/S',residuals)
                    kernel_checks.append(dict(case=key,implementation=label,scaled_caloric_errors=residuals))
    def assert_frozen():
        assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip() == git_head
        assert hashes(sources) == sourcehash and hashes(rustpaths) == rusthash and hashes(cargopaths) == cargohash
        assert sorted((ROOT/'rust').rglob('*.rs')) == rustpaths
        assert [p for p in cargo_candidates if p.exists()] == cargopaths
        assert sorted(set(source_build.glob('*.mod')) | set(out.glob('*.mod'))) == module_paths
        assert fingerprints(binaries+object_paths+module_paths+[object_manifest]) == frozen
    assert_frozen()
    preflight = dict(schema_version=1, rows=len(rows), group_counts={f'{p}_{m}':len(keys) for (p,m),keys in groups.items()},
                     maximum_absolute_differences=maxima, canonicalized_endpoints=endpoint_pairs,
                     source_warnings=warnings, audit=audit, independent_kernel_checks=kernel_checks,
                     original_kernel_states=[dict(model=k[0],temperature=k[1],volume=k[2],values=v) for k,v in kernel_cache.items()],
                     kernel_fields=['residual_alpha','P','P_V','U','S','H','Cv'],
                     outputs=dict(source=ref,native=native,source_status=statuses))
    (out/'meos-flash-benchmark-preflight.json').write_text(json.dumps(preflight,indent=2)+'\n', encoding='utf-8')
    if not args.measure:
        print(json.dumps(dict(preflight_rows=len(rows), timed_eligible_rows=len(rows), audit=audit,
                              source_warnings=len(warnings), independent_kernel_checks=len(kernel_checks), timing='not measured'),indent=2))
        return
    samples = {}; checksums = []
    for sample in range(args.repeats+1):
        runners = [('source', [str(source_exe),'time',str(manifest)], 'source_'),
                   ('native', [str(native_exe),str(manifest)], 'native_')]
        if sample % 2: runners.reverse()
        times = {}
        for label, command, prefix in runners:
            raw = parse_rows(subprocess.check_output(command, cwd=ROOT, env=env, text=True), prefix, 9)
            times[label] = {k.removeprefix(prefix):v for k,v in raw.items()}
        assert times['source'].keys() == times['native'].keys()
        assert set(times['source']) == {f'{p}_{m}' for p,m in groups}
        for group, keys in groups.items():
            key = f'{group[0]}_{group[1]}'
            for label, expected in [('source',ref), ('native',native)]:
                values = times[label][key]
                assert values[0] > 0.
                magnitude = [math.fsum(abs(expected[k][j])*(ITERATIONS//len(keys)+(i < ITERATIONS % len(keys)))
                                      for i,k in enumerate(keys)) for j in range(8)]
                verify_checksum(values[1:], expected_checksum(expected, keys), magnitude)
                if sample: samples.setdefault(label+'_'+key, []).append(values[0])
        checksums.append(dict(sample_index=sample, warmup=sample == 0,
                              process_order=[r[0] for r in runners], **times))
        print(f'MEOS flash benchmark sample {sample}/{args.repeats} (0 warm-up)', flush=True)
    assert_frozen()
    report = dict(schema_version=1, measured_at=datetime.now(timezone.utc).isoformat(), platform=platform.platform(), processor=platform.processor(),
        git_head=git_head,
        rustc=subprocess.check_output([str(rustc),'-Vv'],env=env,text=True).strip(), rustflags=env.get('RUSTFLAGS'),
        compiler=subprocess.check_output([str(compiler),'--version'],env=env,text=True).splitlines()[0], compiler_flags=flags,
        reference_metadata=identity, probe_metadata=probe_identity, native_sources_sha256=rusthash,
        benchmark_sources_sha256=sourcehash, cargo_inputs_sha256=cargohash, runtime_inputs_sha256=frozen,
        source_build_dir=str(source_build), numeric_input_sha256=digest(manifest), numeric_input_rows=rows,
        preflight=preflight, profile_names=PROFILES, lane_names=LANES, iterations_per_sample=ITERATIONS,
        timed_samples_ns=samples, median_ns={k:statistics.median(v) for k,v in samples.items()}, raw_checksums_by_sample=checksums,
        timing_definitions='Parsing, construction and process startup are excluded. Each model is initialized once per worker. TP/PH/PS include solving and equilibrium-property evaluation. Cold/guessed caloric calls are separate lanes with identical bounds and supplied guesses. The scalar TV lane makes four H/S calls at fixture phase volumes, with Rust temperature preparation outside its timer. Source activate_model is included per operation; bounds are set only for PH/PS. Original physical neon LJTS is scaled to reduced units. Qualified caloric endpoints and nearcritical rows are included with explicit bulk-volume canonicalization; every raw pair is retained. Pure TP without specified quality at coexistence and source-failure/known-gap-only fixtures are outside this common successful cohort. There are no per-input mismatch exclusions. Each profile/lane group cycles 1000 fresh operations in fixed manifest order. One warm-up precedes measured samples; process order alternates. Aggregate checksums are checked against independently preflighted per-case outputs with exact visit counts. All source, Rust, Cargo, binary, object, metadata and manifest fingerprints are checked after preflight and after timing.')
    with Path(args.output).open('x',encoding='utf-8',newline='\n') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(dict(median_ns=report['median_ns']),indent=2))


if __name__ == '__main__':
    main()
