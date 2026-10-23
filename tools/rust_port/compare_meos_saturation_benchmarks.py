"""Audit the nested/coupled inverse builds on identical source-qualified inputs."""
import hashlib,json,math,statistics
from reference import ROOT


def compare():
    paths=[ROOT/'docs/rust/meos-saturation-before-coupled-inverse-benchmark.json',ROOT/'docs/rust/meos-saturation-benchmark.json']
    before,after=[json.loads(p.read_text()) for p in paths]
    for field in ['reference','platform','processor','rustc','compiler','compiler_flags',
                  'numeric_input_rows','numeric_input_sha256','profiles','points_per_profile',
                  'iterations_per_sample','cargo_inputs_sha256','source_preflight_states','native_preflight_states']:
        assert before[field]==after[field],field
    assert before['native_sources_sha256'].keys()==after['native_sources_sha256'].keys()
    changed=[p for p,v in before['native_sources_sha256'].items() if v!=after['native_sources_sha256'][p]]
    assert changed==['rust/src/multiparameter_saturation.rs'],changed
    a,b=before['benchmark_sources_sha256'],after['benchmark_sources_sha256'];assert a.keys()==b.keys()
    changed_evidence=[p for p in a if a[p]!=b[p]]
    assert changed_evidence==['tests/fixtures/meos_saturation.metadata.json'],changed_evidence
    # The numeric source fixture is unchanged. Metadata records new independent
    # source-kernel checks of the optimized native coordinates, not new inputs.
    fixture='tests/fixtures/meos_saturation.csv'
    assert a[fixture]==b[fixture]
    profiles={}
    for suffix in after['points_per_profile']:
        for implementation in ['native','source']:
            name=implementation+'_'+suffix
            old,new=before['samples_ns'][name],after['samples_ns'][name]
            assert len(old)==len(new)==5 and all(math.isfinite(v) and v>0. for v in old+new)
            assert statistics.median(old)==before['median_ns'][name]
            assert statistics.median(new)==after['median_ns'][name]
            key='native' if implementation=='native' else 'source'
            x=before['checksums'][key][name][1:];y=after['checksums'][key][name][1:]
            assert len(x)==len(y)==5
            for j,(u,v) in enumerate(zip(x,y)):
                floor=[2e-7,1e-2,1e-12,1e-12,1e-4][j]*1000.
                assert math.isfinite(u) and math.isfinite(v) and abs(u-v)<floor+2e-6*max(abs(u),abs(v)),(name,j,u,v)
        name='native_'+suffix
        profiles[suffix]=dict(before_native_median_ns=before['median_ns'][name],after_native_median_ns=after['median_ns'][name],
            native_build_speedup=before['median_ns'][name]/after['median_ns'][name],
            source_over_current_native=after['median_ns']['source_'+suffix]/after['median_ns'][name],
            before_native_range_ns=[min(before['samples_ns'][name]),max(before['samples_ns'][name])],
            after_native_range_ns=[min(after['samples_ns'][name]),max(after['samples_ns'][name])])
    digest=lambda p:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    patch=ROOT/'docs/rust/meos-saturation-coupled-inverse.patch'
    report=dict(generator='tools/rust_port/compare_meos_saturation_benchmarks.py',
        before_report=paths[0].relative_to(ROOT).as_posix(),before_report_sha256=digest(paths[0]),
        after_report=paths[1].relative_to(ROOT).as_posix(),after_report_sha256=digest(paths[1]),
        changed_native_sources=changed,changed_qualification_evidence=changed_evidence,
        baseline_reproduction_patch=patch.relative_to(ROOT).as_posix(),baseline_reproduction_patch_sha256=digest(patch),
        baseline_module_sha256=before['native_sources_sha256']['rust/src/multiparameter_saturation.rs'],
        identical_source_fixture_numeric_inputs_benchmark_code_and_cargo=True,profiles=profiles,
        scope='Five profiles and 76 qualified source states. Only the native coexistence module changes between builds. The public saturation-pressure method retains its algorithm; compiler optimization can change its timing too. Three cold inverse profiles improve substantially; LJ full krypton is similar and LJ full argon regresses relative to the nested Rust build. Separate source ratios apply only to the current measured workload. Independent source-kernel checks in metadata are refreshed for the optimized native points; the original numeric source fixture is unchanged.')
    (ROOT/'docs/rust/meos-saturation-inverse-comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(profiles,indent=2))


if __name__=='__main__':compare()
