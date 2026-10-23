# Mixture envelope pressure and temperature maxima

`Cubic::envelope_extrema` detects interior pressure and temperature peaks in a
sampled bubble/dew curve and refines them to stationary coexistence states.
The largest successfully refined candidates are returned with full phase
compositions, volumes and a stationarity residual. Failed candidates remain
explicit. A partial or metastable curve does not prove global maxima or phase
stability.

## API

```rust
use thermopack_rs::{Cubic, CubicKind};
use thermopack_rs::envelope::EnvelopeOptions;
# fn main() -> Result<(), thermopack_rs::Error> {
let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"])?;
let curve = model.phase_envelope(&[0.5, 0.5], EnvelopeOptions::default())?;
let maxima = model.envelope_extrema(&[0.5, 0.5], &curve)?;
let maximum_pressure = maxima.pressure.unwrap().state;
let maximum_temperature = maxima.temperature.unwrap().state;
assert!(maxima.failures.is_empty());
# Ok(())
# }
```

`refine_envelope_extremum(feed, guess, EnvelopeExtremum::Pressure)` or
`::Temperature` performs a local refinement from a `MixtureSaturation` guess.
It preserves the supplied fixed-phase formulation. A local stationary root
alone does not establish a maximum. Automatic analysis tries the peak and
neighboring guesses, then selects the largest refined candidate of each kind.
These methods do not insert points into the original trace or change existing
envelope sampling. They accept at least two active feed components; pure
critical endpoints are handled by the separate pure-envelope API.

Automatic candidates must also be maxima relative to corrected neighboring
states on both sides, using a relative 5e-5 change of the free coordinate.
An unresolved neighbor or a stationary minimum is retained as a failed
candidate instead of being reported as a maximum.

Missing peaks use `None`, including curves cut short by point or pressure
limits. Endpoints are not relabeled as stationary extrema. The automatic
analysis retains candidate index, kind and error for failed refinements. A
formulation switch may jump across closely spaced turns; in that case the
local search interval permits one neighboring cell beyond the sampled span.

## Numerical method

The fixed feed is the parent phase and the other phase is incipient. With
unnormalized incipient weights w_i, the last equations are

```text
pressure extremum:     F = T * sum_i w_i*(dT ln(phi_inc_i)-dT ln(phi_parent_i))
temperature extremum:  F = P * sum_i w_i*(dP ln(phi_inc_i)-dP ln(phi_parent_i)).
```

These are the weighted fugacity-derivative conditions used by the fork's
`criconden_fun_newton`, made dimensionless. On a separated coexistence state
they locate zero dP/dT or zero dT/dP. Both also vanish as the phases coincide.
The solver divides F by the squared normalized phase-composition distance Q.
This prevents small raw residuals at a coalescing state from being accepted as
a pressure or temperature maximum. The reported `stationarity_residual` is
the unscaled, dimensionless F.

The first path brackets stationarity while correcting coexistence at fixed T
for a pressure extremum or fixed P for a temperature extremum. Analytic
tangents predict the inner state; safeguarded secant/bisection closes the
outer equation. Inner coexistence tolerance is 1e-12, and the scaled outer
residual is below 1e-9. Alternative predictors handle changes in phase-ratio
orientation. The augmented Newton fallback retains analytic coexistence and
material-balance rows and centrally differences only the stationarity row.
It solves to 1e-10 and restricts steps to the warm ratio orientation.

Q below 1e-12 is rejected. A separate arithmetic cancellation indicator,
16*epsilon*scale*sum(abs(w_i)*(abs(derivative_inc_i)+abs(derivative_parent_i)))/Q,
rejects values above 1e-8. This is a resolution guard, not an interval proof
covering every EOS derivative error. General critical regularization and
arbitrarily small phase contrast remain open. Prepared phase data and
analytic T/P/mole derivatives reuse the existing native envelope machinery.

## Defect in the fork's optional extrema path

In `src/saturation_curve.f90`, `fun_criconden` calls

```fortran
iter = sat_newton(Z,K,t,p,beta,s,log(x),ierr)
```

The eighth argument is `Yphase`, not `ierr`; the local integer is uninitialized.
The Jacobian then uses that value as a phase flag. This is undefined source
behavior. Direct observations of the unmodified reference are retained
separately, including access violations and nonclosing or
nonstationary returns. No native timing is presented as a speedup over this
defective path, and its raw failure/trajectory behavior is not qualified.

The qualification oracle therefore uses an **independent source-backed
reference algorithm**: scalar bisection through the exact source EOS and
`sat_newton`, explicitly passing the vapor-phase and status arguments. Source
residual/Jacobian functions polish inner coexistence because the primitive's
1e-8 stopping tolerance is too loose for a coordinate reference near critical
points. This helper is not the original optional extrema solver. It also
checks source phase pressure and compares states on both sides of each root
to verify a maximum. Every helper failure and all direct legacy observations
remain in metadata; the original source/library is unchanged.

The retained local sample has 14 direct legacy calls: eight access violations
and six returned curves containing twelve reported extrema. Ten reported T/P
pairs agree with the independent coordinate reference, but only seven also
pass the independent stationarity/coexistence classification. Five are
nonstationary under the scaled derivative check, including two nonclosing
coexistence states at the default 10 Pa and 100 MPa pressure bounds. These are
observations of this build, not portable guarantees about undefined behavior.

## Qualification

The matrix attempts 110 source curves and obtains 147 independently verified
maxima. All 110 ordinary source traces return; 29 curves have one or more
recorded scalar-reference search errors. No finite reference maxima are
discarded and none have phase-pressure warnings. Coverage includes four
classic cubic geometries, translated PR/SRK, UNIFAC, binary/ternary/zero
feeds, cryogenic quantum cases, activity-model attempts and curve limits.
An attempted model/configuration without a verified maximum is not counted as
a qualified extremum; detailed successful cases and errors remain in metadata.

| Configuration | Verified reference maxima |
|---|---:|
| Classic PR / SRK / RK / VDW | 31 / 31 / 24 / 14 |
| Translated PR / SRK | 15 / 15 |
| UNIFAC PR / SRK | 8 / 8 |
| Quantum H2/Ne | 1 |

The HV1/HV2/NRTL and H2/D2 attempts do not qualify any extrema in this matrix;
their limits, absent peaks or failed reference searches remain explicit.
There are 39 scalar-reference errors across the 29 affected curves.

Automatic native analysis matches all 147 reference maxima. In 97 comparisons
the physical trace's density ordering reverses the source's phase labels;
the comparison explicitly complements the phase flag and swaps both volumes
and complete compositions. These are qualified phase permutations, rather
than exact raw phase-array identity. Tests also preserve zero components,
close phase pressure/fugacity and check neighboring states on both sides of
every automatically recovered maximum. An off-root test checks every row of
the augmented Jacobian for classic, translated, UNIFAC, Huron-Vidal and quantum
states, in both fixed-phase formulations.

Local warm correction strictly matches 146 of 147 supplied guesses in the
local release build. PR methane/ethane, from the source's 180 K dew-start trace,
has one pressure-maximum guess that reaches the phase-contrast resolution
guard. Automatic analysis recovers the same maximum from a neighboring guess.
This local error is retained, not counted as strict local parity. The validator
fails by default; `--allow-known-gaps` permits only this exact named error and
requires every automatic maximum to match. It does not permit coordinate
differences or additional failures.

Coordinate comparisons use relative 2e-6, plus absolute floors 2e-5 K, 2 Pa,
1e-10 m³/mol and 2e-8 for phase mole fractions. The phase flag must match after
the explicitly recorded permutation. Source-backed reference coexistence is
polished below 2e-12 and scaled stationarity below 2e-9. Both source and native
neighbor checks verify local maxima independently of coordinate agreement.

Complete raw optional-solver fidelity, insertion into exactly matching legacy
sample grids, global stability, extrema of critical loci and other EOS families
remain pending. Existing [envelope](envelope.md),
[phase-fraction](quality-envelopes.md) and [pure endpoint](pure-envelope.md)
qualifications remain separate.

## Performance

The [native-only report](envelope-extrema-benchmark.json) times warm refinement
from eight source-backed guesses: pressure and temperature maxima for PR/SRK
methane/ethane, translated PR methane/propane and PR UNIFAC methane/propane.
All ten output quantities are compared before timing, including phase
compositions/volumes and stationarity. Each sample performs 1000 refinements
and checks every separate output sum. Five measured samples follow one warm-up.
Model construction, parsing, complete tracing and reference generation are
excluded. This is an absolute native measurement, with no source speedup claim.

Local Windows release medians per warm refinement:

| Profile | Preflight states | Native median |
|---|---:|---:|
| PR methane/ethane | 2 | 493.53 µs |
| SRK methane/ethane | 2 | 520.77 µs |
| Translated PR methane/propane | 2 | 305.88 µs |
| UNIFAC PR methane/propane | 2 | 1546.74 µs |

The report lists the ten compared fields explicitly. Branch dP/dT is computed
as part of the returned state, but is singular at a temperature maximum and
is not a coordinate-parity or checksum quantity. Inputs, raw timings, separate
checksums and native/harness/Cargo fingerprints are retained and audited after
timing. These measurements establish a reproducible starting point for further
optimization.

The scalar path currently predicts every trial from one fixed corrected seed.
Reusing nearer corrected states and their tangents is an optimization
candidate, subject to the same separation/maximum checks and a fresh benchmark.
