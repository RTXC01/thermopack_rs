# Mixture VLE phase envelopes

The cubic port now follows mixture bubble/dew curves through temperature and
pressure turns and the near-critical formulation switch. It solves the
fugacity and phase-normalization equations in native Rust. The full source
adaptive grid and failure trajectory are qualified only for the subsets below.

## APIs and units

`Cubic::phase_envelope(z, options)` returns a `PhaseEnvelope` containing
`MixtureSaturation` points, a termination reason and an optional failure.
The default start is a dew point at 0.1 MPa, the maximum pressure is 15 MPa,
the maximum logarithmic coordinate step is 0.25 and the point limit is 1000.
The minimum temperature defaults to the model's bound; positive lower user
limits are clamped to that bound. Compositions must be normalized mole
fractions, with at least two active species. Exact zero species remain absent.

`EnvelopeConstraint::Pressure` takes Pa, `Temperature` takes K and
`LogK { component, value }` uses the zero-based index in the original feed
(which must refer to an active species) and dimensionless ln(y_i/x_i).
Automatic tracing starts at T or P. A log K
constraint is available for the warm corrector, not automatic initialization.
Volumes are m³/mol, and `dpdt` is Pa/K along the fixed-parent-composition
boundary. That slope can diverge at a temperature turn.

```rust
use thermopack_rs::{Cubic, CubicKind};
use thermopack_rs::envelope::{EnvelopeOptions, EnvelopeEnd};
# fn main() -> Result<(), thermopack_rs::Error> {
let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"])?;
let curve = model.phase_envelope(&[0.5, 0.5], EnvelopeOptions::default())?;
assert_eq!(curve.end, EnvelopeEnd::InitialPressure);
# Ok(())
# }
```

The physical API uses an absolute 1e-9 equation tolerance, orders returned
phases by molar volume, updates the fixed-phase `kind` accordingly, and closes
requested pressure/temperature boundaries. `InitialPressure`, `MaximumPressure`
and `MinimumTemperature` identify the reached boundary. `PointLimit` is an
explicit bounded prefix. `NoConvergence` retains the accepted prefix with
`failure=Some(error)`. Invalid input or failed initial-point construction
returns an error before a trace exists.

`Cubic::baseline_phase_envelope` retains the fork's 1e-8 equation tolerance,
formulation labels and strict pressure-limit comparison. A label can switch
while both EOS compositions have one volume root, before physical phase
density order changes. Consequently baseline liquid/vapor volume labels can
be reversed near the critical region. The physical API's reporting order is
a separate convention; it does not establish global phase stability.

`Cubic::refine_envelope_state(seed, constraint)` keeps the seed's specified
phase composition fixed and corrects T/P or one log K. It uses the analytic
tangent as predictor, with a stationary-seed fallback. This is a local solve;
use a seed on the intended branch. The corrector rejects trivial coalesced
solutions and reports failure rather than silently returning K=1.

## Equations and continuation

With parent composition z and vapor fraction beta equal to zero or one,

```text
x_i = z_i / (1 - beta + beta K_i)
y_i = K_i x_i
ln K_i + ln phi_V,i - ln phi_L,i = 0
sum(y_i - x_i) = 0
```

Active log K values, log T and log P are the unknowns. An additional equation
fixes the selected continuation coordinate. Analytic T/P/mole fugacity
derivatives form the augmented Jacobian. Normalization derivatives account
for the changing incipient phase. Each predictor follows the most sensitive
logarithmic variable, so a temperature or pressure turn does not require a
cold restart. The safeguarded Newton corrector uses bounded logarithmic steps
and backtracking. Step size adapts to iteration count, with retries for large
T/P jumps or failed correction. Near criticality a phase-ratio coordinate
avoids fixing an almost singular T or P direction.

Initialization first estimates coexistence with both compositions fixed to
the feed, then corrects the actual mixture boundary. This source-guided
pseudo-saturation seed matters in quantum H2/Ne at 1 MPa: a Wilson-first local
Newton solve can select a different dense-fluid solution near 29.90 K, while
the fork's gas/liquid envelope begins at 33.05 K. A dedicated regression checks
the selected branch and the complete pressure loop.

The converged analytic Jacobian also supplies the continuation tangent.
An exact-coordinate cache avoids recalculating it for point reporting and
the next predictor. Its key includes the fixed coordinate and formulation;
phase changes invalidate it. Root-count checks reconstruct the current
composition rather than relying on a previous failed trial. Tests independently
check Jacobians, cache invalidation and normalized coexistence equations.
Models own their parameters; trace workspaces and caches are call-local.

## Qualification against this fork

`tools/rust_port/validate_envelope.py` makes 98 isolated exact-source calls:
classic PR/SRK/RK/VDW, translated PR/SRK, UNIFAC and quantum cubic; binary and
ternary feeds; zero species; a nearly pure feed; 0.1/1 MPa starts; temperature
starts; bubble starts; low-temperature, pressure and point limits. Quantum
cases include H2/Ne and H2/D2 with a 2 K model bound.

All 98 source calls returned finite nonempty arrays. There were no fatal calls,
censored calls, nonfinite source states or pressure/fugacity closure warnings.
Termination is inferred from the arrays, because the wrapper does not return
its internal control cause: 87 return to the starting pressure, two reach a
temperature limit, four reach their point limit, and five are unclassified
partial prefixes. A finite prefix is not claimed to be a completed envelope.

| Comparison | Evidence |
|---|---|
| Source state recovery | All 6185 sampled T/P/phase-flag/volume/composition states match after warm projection from independent native curves |
| Full sampled-grid parity | 71 curves / 5193 stored states match their counts, order, labels and coordinates |
| Remaining grid differences | 27 curves retain count, sampling or partial-trajectory differences in metadata |
| Physical trace checks | Every returned node in the same matrix is checked for fugacity equality, pressure closure, zero species and reporting order |

Projection selects a nearby native state on the same composition branch and
locks one source coordinate. All other variables are solved by Rust; source
coordinates are not substituted as output. Source Jacobian sensitivities choose
a regular comparison coordinate, and source backward errors are retained as
diagnostics only. They are not used to enlarge tolerances. The comparison uses
relative 2e-6, with absolute floors 2e-5 K, 2 Pa, 1e-10 m³/mol and 2e-8 for
composition. Phase flags are checked as well.

The arrays include two pressure-cap stagnation cases that each fill the
1000-point buffer with repeated terminal states. There are 4225 distinct
stored state tuples within the 98 curves. Reporting all 6185 samples preserves
the fork's output; it is not a claim of 6185 independent thermodynamic states.

The fixtures and metadata retain every source state, partial return, source
identity, Jacobian backward-error diagnostic, comparison result and hash.
Tests recover every state from a separately traced curve and additionally
check the physical API. A source partial return is allowed to remain partial;
its accepted physical states must still close the equations.

## Source pressure-cap convention

The source clamps a logarithmic predictor to log(Pmax), then tests
exp(log(Pmax)) >= Pmax. At 2 MPa on this platform, the exponential rounds below
the limit, so the comparison misses the boundary. Repeated identical predictors
can fill all 1000 entries. Baseline tracing preserves this behavior; physical
tracing recognizes the logarithmic boundary and returns `MaximumPressure`.
The source also can stop just above the model Tmin when correction fails;
the physical tracer attempts the exact temperature endpoint.

## Limits of the milestone

This is local VLE continuation, with one phase composition fixed. Global
multiphase stability, liquid/liquid/vapor and solid envelopes, general fixed
quality, pure-fluid envelope routing, cricondenbar/cricondentherm refinement,
critical-point insertion, other EOS families and complete raw iteration/error
fidelity remain pending. Existing local critical-point APIs are separate;
a maximum over sampled nodes is not an exact criconden point.

Full raw sampled-grid parity is specifically **not established for 27 of
these 98 curves**, even though every observed source state is recovered.
Broader alpha/beta/reference combinations and critical degeneracies require
their own qualification. See the [coverage ledger](coverage.md).

## Performance

`benchmark_envelope.py` times complete curves with identical source/native
adaptive grids: PR/SRK methane/ethane, a PR ternary mixture and PR-UNIFAC
methane/propane, with 48/48/59/69 points. All 224 states are compared before
timing, including both phase volumes and compositions. Each timed sample
checks point counts, phase-label sums, T/P/volume sums and every composition
sum separately. Cold initialization and the entire continuation are included;
model construction and Python/FFI are excluded. The source post-evaluates both
volumes; Rust carries them and additionally computes branch slopes.

The 2026-10-23 local medians, per complete curve, are:

| Curve | Fortran ms | Rust ms | Speedup |
|---|---:|---:|---:|
| PR methane/ethane | 5.010 | 1.901 | 2.64× |
| SRK methane/ethane | 5.490 | 1.900 | 2.89× |
| PR methane/ethane/propane | 6.154 | 1.245 | 4.94× |
| PR-UNIFAC methane/propane | 16.987 | 15.681 | 1.08× |

The UNIFAC workload is only slightly faster and remains an optimization
target. The [machine-readable report](envelope-benchmark.json) records compiler settings, five samples after one
warm-up, hashes of native/harness/Cargo inputs and the exact-source identity.
These are scoped workload measurements, not whole-library performance claims.
