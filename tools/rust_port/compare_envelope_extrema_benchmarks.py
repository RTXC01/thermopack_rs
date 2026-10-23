"""Audit the controlled before/after native extrema measurements."""
import argparse
import hashlib
import json
import math
import statistics
from pathlib import Path
from reference import ROOT


def compare(before_path, after_path, output_path):
    before = json.loads(before_path.read_text())
    after = json.loads(after_path.read_text())
    for field in [
        'reference', 'platform', 'processor', 'rustc', 'numeric_input_cases',
        'source_backed_reference_preflight_states', 'iterations_per_sample',
        'compared_outputs', 'benchmark_sources_sha256', 'cargo_inputs_sha256',
    ]:
        assert before[field] == after[field], field
    assert before['native_sources_sha256'].keys() == after['native_sources_sha256'].keys()
    changed = sorted(path for path, value in before['native_sources_sha256'].items()
                     if value != after['native_sources_sha256'][path])
    assert changed == ['rust/src/envelope_extrema.rs'], changed
    profiles = {}
    maximum_differences = [0.] * 10
    iterations = after['iterations_per_sample']
    assert before['samples_ns'].keys() == after['samples_ns'].keys()
    for name in before['samples_ns']:
        old_samples, new_samples = before['samples_ns'][name], after['samples_ns'][name]
        assert len(old_samples) == len(new_samples) >= 5
        assert all(math.isfinite(v) and v > 0. for v in old_samples + new_samples)
        old, new = statistics.median(old_samples), statistics.median(new_samples)
        assert old == before['median_ns'][name] and new == after['median_ns'][name]
        a, b = before['checksums'][name][1:], after['checksums'][name][1:]
        assert len(a) == len(b) == 10
        for j, (x, y) in enumerate(zip(a, b)):
            floor = (2e-5 if j == 0 else 2. if j == 1 else 0. if j == 2 else
                     1e-10 if j < 5 else 1e-9 if j == 9 else 2e-8)
            assert math.isfinite(x) and math.isfinite(y)
            assert abs(x-y) <= floor*iterations + 2e-6*max(abs(x), abs(y)), (name, j, x, y)
            maximum_differences[j] = max(maximum_differences[j], abs(x-y))
        profiles[name] = dict(before_median_ns=old, after_median_ns=new,
                              speedup=old/new, before_range_ns=[min(old_samples), max(old_samples)],
                              after_range_ns=[min(new_samples), max(new_samples)],
                              all_after_samples_faster_than_all_before=max(new_samples)<min(old_samples))
    digest = lambda path: hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
    report = dict(
        generator='tools/rust_port/compare_envelope_extrema_benchmarks.py',
        before_report=before_path.relative_to(ROOT).as_posix(), before_report_sha256=digest(before_path),
        after_report=after_path.relative_to(ROOT).as_posix(), after_report_sha256=digest(after_path),
        changed_native_sources=changed, identical_benchmark_and_cargo_inputs=True,
        identical_source_backed_reference_and_numeric_cases=True,
        compared_outputs=after['compared_outputs'],
        maximum_absolute_checksum_differences=maximum_differences, profiles=profiles,
        scope='Eight warm native P/T extrema refinements across four binary profiles. '
              'Ratios compare the previous Rust fixed-seed method with cached nearby predictors; '
              'they are not speedups against the original defective optional source solver.')
    output_path.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(profiles, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--before', type=Path, default=ROOT/'docs/rust/envelope-extrema-before-cache-benchmark.json')
    parser.add_argument('--after', type=Path, default=ROOT/'docs/rust/envelope-extrema-benchmark.json')
    parser.add_argument('--output', type=Path, default=ROOT/'docs/rust/envelope-extrema-cache-comparison.json')
    args = parser.parse_args()
    compare(args.before.resolve(), args.after.resolve(), args.output.resolve())
