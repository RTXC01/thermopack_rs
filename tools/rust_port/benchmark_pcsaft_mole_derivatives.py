"""Build and measure the internal nonassociation PC-SAFT F_n vector.

Run `build --help` for a sequential, manifest-producing build. Run `--help`
for measurement; --manifest is mandatory and all frozen hashes are rechecked.
For mixtures fugacity_coefficients(rho, 1) returns F_n. The pure special case
is excluded because that adapter does not return F_n for one component.
Prepared times reuse model and prepared state; fresh times construct a prepared
state per gradient (but reuse the model). The original times include its normal
internal cache updates and its model is initialized outside the timed loop.
"""
import argparse
import hashlib
import json
import math
import os
import platform
import statistics
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
COMPONENTS = ['CO2', 'C1', 'C2', 'C3', 'NC4', 'NC5', 'NC6', 'NC7']

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def checked(command, cwd=None):
    result = subprocess.run([str(x) for x in command], cwd=cwd, text=True,
                            encoding='utf-8', errors='replace', capture_output=True)
    if result.returncode:
        raise RuntimeError(f"command failed: {command}\n{result.stdout}\n{result.stderr}")
    return result.stdout

def tree_files(root):
    root = Path(root)
    paths = sorted((root / 'rust/src').rglob('*.rs'))
    paths += [root / 'Cargo.toml', root / 'Cargo.lock']
    return [p.resolve() for p in paths if p.exists()]

def git_state(root):
    return dict(head=checked(['git', 'rev-parse', 'HEAD'], root).strip(),
                branch=checked(['git', 'rev-parse', '--abbrev-ref', 'HEAD'], root).strip(),
                tracked_rust_diff=checked(['git', 'diff', '--binary', 'HEAD', '--',
                                           'rust', 'Cargo.toml', 'Cargo.lock'], root),
                status=checked(['git', 'status', '--porcelain=v1', '--untracked-files=all',
                                '--', 'rust', 'Cargo.toml', 'Cargo.lock'], root))

def hashes(paths):
    return {str(Path(path).resolve()): digest(path) for path in paths}

def validate_hashes(expected):
    changed = []
    for path, expected_hash in expected.items():
        if not Path(path).is_file() or digest(path) != expected_hash:
            changed.append(path)
    if changed:
        raise RuntimeError('frozen files missing or changed: ' + ', '.join(changed))

def build_main(argv):
    p = argparse.ArgumentParser(description='Build sequentially and freeze a benchmark manifest')
    p.add_argument('--root', required=True, type=Path, help='baseline Rust checkout')
    p.add_argument('--candidate', type=Path, default=ROOT)
    p.add_argument('--source-root', required=True, type=Path)
    p.add_argument('--baseline-head', required=True, help='expected baseline Git commit')
    p.add_argument('--candidate-head', required=True, help='expected candidate Git commit')
    p.add_argument('--compiler-dir', required=True, type=Path)
    p.add_argument('--cargo', required=True, type=Path)
    p.add_argument('--rustc', required=True, type=Path)
    p.add_argument('--gfortran', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path, help='build manifest JSON')
    a = p.parse_args(argv)
    baseline, candidate, source = (a.root.resolve(), a.candidate.resolve(), a.source_root.resolve())
    cargo, rustc, gfortran = (a.cargo.resolve(), a.rustc.resolve(), a.gfortran.resolve())
    work = candidate / '.reference/pcsaft-mole-bench'
    work.mkdir(parents=True, exist_ok=True)
    dll = source / '.reference/thermopack.dll'
    metadata_path = source / '.reference/metadata.json'
    metadata = json.loads(metadata_path.read_text(encoding='utf-8'))
    if metadata.get('baseline') != 'd68c794' or metadata.get('library_sha256') != digest(dll):
        raise RuntimeError('original reference metadata or DLL hash mismatch')
    modules = [source / '.reference' / x for x in ('thermopack_var.mod', 'pc_saft_nonassoc.mod')]
    scripts = [HERE / x for x in ('benchmark_pcsaft_mole_derivatives.py', 'pcsaft_mole_bench.rs',
               'pcsaft_mole_bench.f90', 'pcsaft_mole_bench_worker.py')]
    source_files = [source / 'src/pc_saft_nonassoc.f90', source / 'tools/rust_port/reference.py', metadata_path]
    source_files += sorted((source / 'addon/pycThermopack/thermopack').rglob('*.py'))
    frozen_inputs = tree_files(baseline) + tree_files(candidate) + [dll, *modules, *scripts, *source_files, cargo, rustc, gfortran]
    before = hashes(frozen_inputs)
    states_before = {str(root): git_state(root) for root in (baseline, candidate)}
    for root, expected in [(baseline, a.baseline_head), (candidate, a.candidate_head)]:
        actual = states_before[str(root)]['head']
        if not actual.startswith(expected) or len(expected) < 7:
            raise RuntimeError(f'Git identity mismatch for {root}: {actual} versus {expected}')
    os.environ['RUSTC'] = str(rustc)
    commands, builds = [], {}
    def execute(cmd, cwd=None):
        commands.append(dict(argv=list(map(str,cmd)), cwd=str(cwd) if cwd else None))
        return checked(cmd, cwd)
    for label, root in [('baseline',baseline), ('candidate',candidate)]:
        cmd = [cargo, 'build', '--release', '--lib', '--message-format=json', '-j1']
        cargo_output = execute(cmd, root)
        matches = []
        for line in cargo_output.splitlines():
            try:
                message = json.loads(line)
            except json.JSONDecodeError:
                continue
            if message.get('reason') == 'compiler-artifact' and message.get('target', {}).get('name') == 'thermopack_rs':
                matches.extend(Path(f).resolve() for f in message.get('filenames', []) if f.endswith('.rlib'))
        if len(matches) != 1:
            raise RuntimeError(f'{label}: expected one cargo rlib, found {matches}')
        rlib = matches[0]
        exe = work / (label + ('.exe' if os.name == 'nt' else ''))
        execute([rustc, '--edition=2024', '-O', '-C', 'lto=thin', '--extern',
                 'thermopack_rs=' + str(rlib), '-L', 'dependency=' + str(rlib.parent if rlib.parent.name == 'deps' else rlib.parent / 'deps'),
                 HERE / 'pcsaft_mole_bench.rs', '-o', exe], root)
        builds[label] = dict(rlib=str(rlib), exe=str(exe))
    helper = work / ('pcsaft_mole_bench.dll' if os.name == 'nt' else 'pcsaft_mole_bench.so')
    execute([gfortran, '-O3', '-fdefault-real-8', '-fdefault-double-8', '-frecursive', '-shared',
             '-I', source / '.reference', HERE / 'pcsaft_mole_bench.f90', dll, '-o', helper], work)
    validate_hashes(before)
    states_after = {str(root): git_state(root) for root in (baseline, candidate)}
    if states_before != states_after:
        raise RuntimeError('Git source state changed during build')
    outputs = [Path(builds[x][key]) for x in builds for key in ('exe','rlib')] + [helper]
    manifest = dict(schema=1, baseline=str(baseline), candidate=str(candidate), source_root=str(source),
                    compiler_dir=str(a.compiler_dir.resolve()), source_dll=str(dll), helper=str(helper),
                    builds=builds, git=states_before, commands=commands, source_metadata=metadata, environment=dict(RUSTC=str(rustc), RUSTFLAGS=os.environ.get('RUSTFLAGS'), CARGO_ENCODED_RUSTFLAGS=os.environ.get('CARGO_ENCODED_RUSTFLAGS'), CARGO_TARGET_DIR=os.environ.get('CARGO_TARGET_DIR')),
                    versions=dict(cargo=checked([cargo,'--version']), rustc=checked([rustc,'--version','--verbose']),
                                  gfortran=checked([gfortran,'--version'])),
                    hashes={**before, **hashes(outputs)})
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('Wrote build manifest:', a.output)

def finite_float(text):
    value = float(text)
    if not math.isfinite(value):
        raise ValueError('nonfinite number')
    return value

def call(command, request, timeout):
    try:
        result = subprocess.run(command, input=json.dumps(request) + '\n', text=True,
                                encoding='utf-8', errors='replace', capture_output=True, timeout=timeout)
        lines = [s[7:] for s in result.stdout.splitlines() if s.startswith('RESULT=')]
        record = dict(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
        if result.returncode == 0 and len(lines) == 1:
            try:
                record['result'] = json.loads(lines[0], parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)), parse_float=lambda x: finite_float(x))
            except (ValueError, TypeError) as e:
                record['failure'] = 'malformed result JSON: ' + str(e)
        else:
            record['failure'] = 'process failed or did not return exactly one result'
        return record
    except subprocess.TimeoutExpired as e:
        return dict(failure='timeout (process killed and reaped)', timeout_seconds=timeout,
                    stdout=str(e.stdout or ''), stderr=str(e.stderr or ''))

def result_vector(record):
    result = record.get('result')
    return result.get('vector') if isinstance(result, dict) else None

def valid_vector(value, count):
    return (count > 1 and isinstance(value, list) and len(value) == count and
            all(type(x) in (int, float) and math.isfinite(x) for x in value))

def sample_error(record, original, preflight, count, tolerance, iterations):
    result = record.get('result')
    if record.get('returncode') != 0 or not isinstance(result, dict):
        return 'failed or missing result'
    vector = result.get('vector')
    if not valid_vector(vector, count) or not valid_vector(original, count) or not valid_vector(preflight, count):
        return 'invalid vector shape or nonfinite component'
    for v, o, p in zip(vector, original, preflight):
        if abs(v-o) > tolerance*(1+abs(o)) or abs(v-p) > tolerance*(1+abs(p)):
            return 'vector differs from preflight or original'
    seconds, checksum = result.get('seconds'), result.get('checksum')
    if type(seconds) not in (float,int) or not math.isfinite(seconds) or seconds <= 0:
        return 'invalid elapsed time'
    if type(checksum) not in (float,int) or not math.isfinite(checksum):
        return 'invalid checksum'
    terms = [(i+1)*v for i,v in enumerate(vector)]
    expected = iterations * math.fsum(terms)
    # Repeated sequential summation can incur O(iterations * nc * epsilon) relative error.
    bound = 16 * sys.float_info.epsilon * iterations * count * iterations * math.fsum(map(abs, terms)) + 1e-12
    if abs(checksum-expected) > bound:
        return 'checksum inconsistent with full vector and iteration count'
    return None

def validate_manifest(manifest):
    required = ('schema','baseline','candidate','source_root','compiler_dir','source_dll','helper',
                'builds','git','commands','versions','hashes','source_metadata')
    if not isinstance(manifest, dict) or any(k not in manifest for k in required) or manifest['schema'] != 1:
        raise ValueError('missing or unsupported build manifest fields')
    cmds = manifest['commands']
    if not isinstance(cmds, list) or len(cmds) != 5:
        raise ValueError('invalid recorded build commands')
    for index in (0,2):
        command = cmds[index].get('argv', [])
        if command[1:] != ['build','--release','--lib','--message-format=json','-j1']:
            raise ValueError('invalid cargo command')
    for index in (1,3):
        command = cmds[index].get('argv', [])
        if not all(flag in command for flag in ('--edition=2024','-O','-C','lto=thin','--extern','-L')):
            raise ValueError('invalid Rust harness command')
    if not all(flag in cmds[4].get('argv', []) for flag in ('-O3','-fdefault-real-8','-fdefault-double-8','-frecursive','-shared')):
        raise ValueError('invalid Fortran helper command')
    if manifest['source_metadata'].get('baseline') != 'd68c794' or manifest['source_metadata'].get('library_sha256') != digest(manifest['source_dll']):
        raise ValueError('original DLL does not match source metadata')
    for label in ('baseline','candidate'):
        for key in ('exe','rlib'):
            if manifest['builds'][label][key] not in manifest['hashes']:
                raise ValueError('build output missing hash')
    for key in ('source_dll','helper'):
        if manifest[key] not in manifest['hashes']:
            raise ValueError('required source or helper hash missing')

def cpu_name():
    if os.name == 'nt':
        try:
            import winreg
            with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, r'HARDWARE\DESCRIPTION\System\CentralProcessor\0') as key:
                return winreg.QueryValueEx(key, 'ProcessorNameString')[0]
        except OSError:
            pass
    return platform.processor() or os.environ.get('PROCESSOR_IDENTIFIER', 'unknown')

def cases():
    for count in (2, 3, 4, 8):
        for simplified in (False, True):
            for density in ('dilute', 'dense'):
                for zero in (False, True):
                    n = [0.1 / count] * count
                    if zero:
                        n[0] = 0.0
                        n[1:] = [0.1 / (count-1)] * (count-1)
                    yield dict(components=','.join(COMPONENTS[:count]), temperature=350.,
                               volume=1e-3 if density == 'dilute' else 1e-4,
                               moles=n, simplified=simplified, density=density, zero_feed=zero)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--preflight-only', action='store_true')
    p.add_argument('--samples', type=int, default=7)
    p.add_argument('--iterations', type=int, default=10000)
    p.add_argument('--warmup', type=int, default=1000)
    p.add_argument('--timeout', type=float, default=120.)
    args = p.parse_args()
    if args.samples < 7 or args.iterations < 1 or args.warmup < 1:
        p.error('require at least seven samples and positive iterations and warmup')
    manifest = json.loads(args.manifest.read_text(encoding='utf-8'))
    validate_manifest(manifest)
    validate_hashes(manifest['hashes'])
    before = manifest['hashes']
    for root, state in manifest['git'].items():
        if git_state(root) != state:
            raise RuntimeError('Git source state changed since build: ' + root)
    args.baseline_bin = Path(manifest['builds']['baseline']['exe'])
    args.candidate_bin = Path(manifest['builds']['candidate']['exe'])
    args.source_root = Path(manifest['source_root'])
    args.compiler_dir = Path(manifest['compiler_dir'])
    args.helper = Path(manifest['helper'])
    out = dict(schema=1, manifest=manifest, manifest_sha256=digest(args.manifest), settings=vars(args).copy(), host=dict(platform=platform.platform(),
               processor=platform.processor(), cpu_name=cpu_name(), python=sys.version), hashes_before=before, cases=[])
    out['settings'] = {k: str(v) if isinstance(v, Path) else v for k,v in out['settings'].items()}
    def run(case, which, mode='prepared', iterations=None):
        req = dict(case, warmup=args.warmup, iterations=iterations or args.iterations)
        if which == 'original':
            cmd = [sys.executable, str(HERE / 'pcsaft_mole_bench_worker.py'), str(args.source_root.resolve()),
                   str(args.compiler_dir.resolve()), str(args.helper.resolve())]
        else:
            exe = args.baseline_bin if which == 'baseline' else args.candidate_bin
            cmd = [str(exe.resolve()), 'simplified' if case['simplified'] else 'standard', case['components'],
                   str(case['temperature']), str(case['volume']), ','.join(map(str,case['moles'])), mode,
                   str(args.warmup), str(req['iterations'])]
        return call(cmd, req, args.timeout)
    for case in cases():
        row = dict(case=case, preflight={}, samples=[])
        for name in ('original', 'baseline', 'candidate'):
            row['preflight'][name] = run(case, name, iterations=1)
        ok = len(case['moles']) > 1
        original = result_vector(row['preflight']['original'])
        ok = ok and valid_vector(original, len(case['moles'])) and row['preflight']['original'].get('returncode') == 0
        for name, tolerance in [('baseline', 2e-6), ('candidate', 2e-11)]:
            vector = result_vector(row['preflight'][name])
            passed = bool(ok and row['preflight'][name].get('returncode') == 0 and valid_vector(vector, len(case['moles'])) and
                all(math.isfinite(a) and math.isfinite(b) and abs(a-b) <= tolerance*(1+abs(a)) for a,b in zip(original,vector)))
            row['preflight'][name]['within_tolerance'] = passed
            ok &= passed
        if ok and not args.preflight_only:
            for sample in range(args.samples):
                order = ['baseline', 'candidate', 'original']
                if (sample // 3) % 2: order.reverse()
                shift = sample % 3
                order = order[shift:] + order[:shift]
                for mode in ('prepared', 'fresh'):
                    for name in order:
                        if name == 'original' and mode == 'fresh':
                            continue
                        row['samples'].append(dict(sample=sample, implementation=name, mode=mode, **run(case,name,mode)))
            summaries = {}
            groups = [('prepared','baseline'), ('prepared','candidate'), ('prepared','original'), ('fresh','baseline'), ('fresh','candidate')]
            valid = True
            for mode, name in groups:
                group = [x for x in row['samples'] if x['mode']==mode and x['implementation']==name]
                if len(group) != args.samples:
                    valid = False
                for sample in group:
                    tolerance = 2e-6 if name == 'baseline' else 2e-11
                    preflight = result_vector(row['preflight'][name])
                    sample['validation_error'] = sample_error(sample, original, preflight, len(case['moles']), tolerance, args.iterations)
                    valid &= sample['validation_error'] is None
                if len(group) == args.samples and all(x.get('validation_error') is None for x in group):
                    times = [x['result']['seconds']/args.iterations for x in group]
                    summaries[mode+'_'+name] = dict(count=len(times), median=statistics.median(times), minimum=min(times), maximum=max(times), stdev=statistics.stdev(times))
            row['timing_valid'] = valid
            row['summaries_seconds_per_call'] = summaries
            row['median_ratios'] = {
                'baseline_over_candidate_prepared': summaries['prepared_baseline']['median'] / summaries['prepared_candidate']['median'],
                'original_over_candidate_prepared': summaries['prepared_original']['median'] / summaries['prepared_candidate']['median'],
                'baseline_over_candidate_fresh': summaries['fresh_baseline']['median'] / summaries['fresh_candidate']['median'],
            } if valid else {}
        out['cases'].append(row)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(out, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    out['hashes_after'] = {path: digest(path) if Path(path).is_file() else None for path in before}
    out['hashes_unchanged'] = (out['hashes_after'] == before and digest(args.manifest) == out['manifest_sha256']
        and all(git_state(root) == state for root, state in manifest['git'].items()))
    args.output.write_text(json.dumps(out, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    if not out['hashes_unchanged'] or any(not r['preflight'][n].get('within_tolerance') for r in out['cases'] for n in ('baseline','candidate')) or (not args.preflight_only and any(not r.get('timing_valid', False) for r in out['cases'])):
        sys.exit(1)

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'build':
        build_main(sys.argv[2:])
    else:
        main()
