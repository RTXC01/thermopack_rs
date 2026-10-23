# Cubic phase-fraction curves

`Cubic::quality_envelope` and `baseline_quality_envelope` continue the
[fixed-quality equations](quality.md) through T/P turns and phase interchange.
They return full phase compositions, volumes, fractions, slopes, the original
feed and explicit termination. These are local two-phase loci; global
multiphase stability is not established.

## API and conventions

```rust
use thermopack_rs::{Cubic, CubicKind};
use thermopack_rs::envelope::{EnvelopeOptions, EnvelopeEnd};
# fn main() -> Result<(), thermopack_rs::Error> {
let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"])?;
let curve = model.quality_envelope(&[0.5, 0.5], 0.25, EnvelopeOptions::default())?;
assert_eq!(curve.end, EnvelopeEnd::InitialPressure);
# Ok(())
# }
```

The fraction is molar and the feed must be normalized. Options use the same
K/Pa/log-coordinate conventions as `phase_envelope`. The requested fraction
chooses the nearest bubble/dew initializer; `options.initial_kind` is unused.
Starts fix T or P. Exact pure feeds and beta=0/1 dispatch to the established
pure or bubble/dew envelope algorithms. Baseline pure classification retains
the source rounding-level convention and reports beta=0 there.

The phase interchange complements beta and inverts K, preserving the original
feed. This locus follows the unordered pair of phase amounts: returned beta
can be the requested fraction or its complement. It must not be assumed
constant across a critical crossing. Physical reporting orders the phases by
volume and complements beta when swapping labels. Baseline reporting retains
the fork's labels. Both representations satisfy the lever rule.

Physical tracing closes maximum-pressure and temperature boundaries, validates
its starting state and retains accepted prefixes when correction fails. The
baseline retains the source's deferred minimum-temperature override: a start
below that optional exit temperature can be accepted. A point cap returns an
explicit `PointLimit`, not an assertion that the entire locus is complete.

## Continuation and critical recovery

The predictor uses the most sensitive log K/T/P coordinate, with direction
inferred from the previous accepted point. Analytic phase-composition
derivatives supply the tangent and Newton corrector. Step lengths adapt to
iteration counts, while accepted jumps are limited to 25 K and 1 MPa.

The step attempts follow the source order: full, half, double, then critical
recovery. The recovery locks a nontrivial log K, holds predicted T/P and tries
the opposite phase formulation. Extra short attempts before this recovery
were found to change quantum H2/D2 trajectories: they advanced too close to
the critical point, making the subsequent branch jump exceed 1 MPa. Preserving
the order closes all three qualified H2/D2 fraction curves. Both source and
physical Rust paths have regression checks for those complete returns.

Near the pressure cap the fork can repeat its last state because
`exp(log(Pmax)) < Pmax` by rounding. The baseline preserves those repetitions;
the physical path closes the cap explicitly. The stored source matrix
therefore has 8,019 states but 6,059 distinct thermodynamic coordinate tuples,
including 1,960 repeated cap states.

The baseline corrector also preserves immediate source exits for a trivial
K=1 trial and for trials at the default pressure/temperature bounds. The fork
zeroes its equations there and reports ierr=-2/-1. Continuing a line search
through an out-of-domain trial changed failure behavior and made the UNIFAC
curve about 51 ms locally. Rejecting that trial immediately reduced a local
profile to about 4.3 ms. The default pressure interval is 10 Pa to 100 MPa;
physical correction retains its separate domain policy. Exact-source probes
and regression tests check both immediate exits.

## Qualification and open differences

The source matrix attempts 158 curves. All return finite arrays, but two
PR/SRK CO2/methane beta=0.25 returns are unqualified. They exit at the global
999 K solver bound with phase compositions that fail the independent feed
material balance. The source residual is zeroed at that bound. Their raw
arrays and classification reasons remain in metadata, rather than being
counted as successful VLE curves.

The other 156 source returns store 8,019 states. They cover four classical
cubic geometries, binary/ternary/zero feeds, translated PR/SRK, UNIFAC and
quantum H2/Ne and H2/D2; fractions 0.25/0.5/0.75; pressure and temperature
starts, point caps, pressure limits and temperature exits.

All counts, order, fractions, T/P, phase volumes and compositions match on
**148 complete sampled grids**. The other eight sampled grids or partial
trajectories remain different. Independent native traces recover **8,010
stored source states** directly or by warm single-coordinate correction.
Two H2/Ne source returns are unclassified partial prefixes; they are not
counted as complete curves.

The immediate source-bound exit optimization in this checkpoint's predecessor
already raised the matching-grid count from 133 to 148 in the checked-in
fixtures. The caloric follow-up corrects the stale test assertion and prose;
the freshly regenerated source comparison leaves those fixture bytes unchanged.

Nine coordinate recovery gaps remain outside ordinary tolerance:

* One near-coincident SRK UNIFAC methane/propane node at beta=0.5 has a
  source Jacobian condition estimate about 4.05e9. The pressure difference is
  about 18.83 Pa, just beyond its ordinary bound; volumes also differ.
* Eight nodes of the partial quantum H2/Ne beta=0.25 return have condition
  estimates above 1e12 and nearly coincident phase compositions. Their raw
  T/volume/composition differences are retained.

These nine states are not counted as strict parity or assigned a rigorous
conditioning allowance. Two additional warm projection differences occur on
grids that already pass every direct coordinate comparison. Those secondary
corrector differences remain separate from curve-state recovery.

Ordinary comparisons retain relative 2e-6, with absolute floors 2e-5 K, 2 Pa,
1e-10 m³/mol and 2e-8 for composition. No tolerance is expanded to hide the
remaining gaps. The strict validator exits unsuccessfully while they remain.
`--allow-known-gaps` permits a regression run only for the explicitly recorded
node identities and conditioning ranges, and still saves/reports every gap.
It does not turn them into matching states.

Independent tests close pressure/fugacity/material balance, preserve exact
zero species and feed, check density ordering and complemented fractions,
close physical pressure/temperature bounds and compare all matching source
grids. Existing bubble/dew, pure-grid and fixed-quality point regressions
remain part of the checks.

HV/NRTL/Wong-Sandler curve qualification, caloric constraints, global stability,
other EOS families, critical insertion, exact raw failure/iteration behavior
and the nine raw coordinate gaps remain pending. Full-library parity is not
established by this checkpoint.

## Performance

Local Windows release medians per complete curve:

| Workload | States | Source | Rust | Speedup |
|---|---:|---:|---:|---:|
| PR methane/ethane | 44 | 3.809 ms | 0.841 ms | 4.53x |
| SRK methane/ethane | 45 | 6.639 ms | 0.822 ms | 8.08x |
| PR methane/ethane/propane | 46 | 7.035 ms | 0.978 ms | 7.19x |
| UNIFAC PR methane/propane | 66 | 17.340 ms | 4.169 ms | 4.16x |

`benchmark_envelope.py --quality` compares four complete curves from an initial
fraction of 0.25: PR and SRK methane/ethane, PR methane/ethane/propane and PR
UNIFAC methane/propane. Their source/native grids have identical counts
44/45/46/66. All **201 states** are compared before timing, including both
volumes, fractions and every phase composition.

Each sample computes ten complete curves. Counts, beta sum and separate
T/P/volume/composition sums are checked again. Five measured samples follow
one warm-up. Model construction and Python/FFI are excluded; cold initial
state, complete continuation/phase interchange and both volumes are included.
Native states additionally carry slopes, residuals and the feed. The
[report](quality-envelope-benchmark.json) retains all samples, source identity
and native/harness/Cargo fingerprints. These are scoped workloads, not
whole-library performance or evidence that the remaining raw gaps are closed.

The [before-bounds report](quality-envelope-before-bounds-benchmark.json)
retains the slower implementation's samples and source fingerprints. The
benchmark's `--initial-only` Rust mode separates cold initialization from
whole-curve cost; initialization was below 1 ms in all four local profiles.
The before/after manifests differ only in `rust/src/quality.rs`; harness and
Cargo hashes match, and every timed native output sum is identical. The
UNIFAC median falls from 51.014 ms to 4.169 ms, a 12.24x native improvement.
