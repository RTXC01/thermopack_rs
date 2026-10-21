# Mixture spinodal points and curves

The cubic model now solves local mixture stability limits from the smallest
eigenvalue of `sqrt(z) A_nn/(RT) sqrt(z)` at fixed temperature and total volume.
The existing `spinodal::spinodal_point` accepts normalized mixtures and selects
liquid/vapor mechanical seeds before refining pressure and the zero eigenvalue.
An embedded pure fluid uses its actual composition; the earlier helper built
an invalid all-ones vector for a multi-component model with zero species.

## Native interfaces

* `PreparedState::mixture_stability(V)` returns the minimum eigenvalue and its
  unit direction, omitting the critical third derivative. Absent species have
  zero direction entries. The analytic TV mole Hessian and native symmetric
  eigenpair are shared with the [critical solver](mixture-critical.md).
* `spinodal_point_with_guess` refines T/V for a specified positive pressure
  from an explicit local seed; no global phase classification is imposed.
* `spinodal_temperature_at_volume` and `spinodal_volume_at_temperature` solve
  a zero eigenvalue for either scalar specification, including negative
  homogeneous pressures where the TV state remains defined.
* `trace_spinodal` maps increasing volume from a liquid-side pressure boundary,
  through temperature turning points toward vapor. Options select starting
  pressure, optional initial liquid temperature, minimum vapor temperature,
  positive log-volume step and a maximum point count. Defaults are 0.1 MPa,
  step 0.02, 0.6 times the composition-weighted component critical temperatures,
  and 1000 points. The result explicitly distinguishes reaching the minimum
  temperature, returning to the initial pressure, and an incomplete point limit.

The pressure solve uses log T/V coordinates, analytic pressure derivatives,
centered eigenvalue derivatives, bounded steps and backtracking. Scalar solves
also use logarithmic variables. Residuals must be below 2e-10. Continuation
predicts the temperature slope implicitly, corrects the zero eigenvalue at the
new volume, and refines either endpoint to its specification. Pressure outputs
are evaluated from the returned physical state. Pure curves use the same volume
continuation; the fork's separate pure 200-point temperature grids are pending.

```rust
use thermopack_rs::{Cubic, CubicKind, Phase};
use thermopack_rs::spinodal::spinodal_point;
use thermopack_rs::mixture_spinodal::{trace_spinodal, SpinodalTraceOptions};
# fn main() -> Result<(), thermopack_rs::Error> {
let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"])?;
let point = spinodal_point(&model, 1e6, Phase::Vapor, &[0.5, 0.5])?;
assert!((point.temperature - 154.02227).abs() < 1e-4);
let curve = trace_spinodal(&model, &[0.5, 0.5], SpinodalTraceOptions::default())?;
assert!(curve.points.len() > 100);
# Ok(())
# }
```

## Source compatibility

The fork's `nestedLoopVT` chooses the sign of its volume walk from the pressure
error alone. On the vapor branch the pressure slope has the opposite direction,
so this walk can cross the turning region and return the liquid-side boundary.
All 54 vapor requests in the qualified mixture grid show that selection. For
example PR C1+C2 at equal fractions and 1 MPa returns approximately 228.55 K
and 7.221e-5 m³/mol in the fork; the independently closing vapor-side boundary
is approximately 154.02 K and 6.566e-4 m³/mol.

`baseline_mixture_spinodal_point` preserves the original pseudo-saturation seed,
phase-dependent offset, mechanical extremum and volume-walk sign, followed by
native refinement. It matches both requested source phase flags in this grid.
It bounds work to 1000 volume steps; original iteration trajectories and all
failed-call behavior are not claimed identical.

`baseline_trace_spinodal` also preserves pressure reporting for mixture curves.
The first pressure is the requested starting pressure, including when an
explicit initial temperature changes the actual pressure. At a temperature
endpoint the fork replaces T and V but retains P evaluated before refinement.
The native compatibility curve retains that pressure too; `trace_spinodal`
reports P evaluated at the final T/V. These conventions are explicit because
those raw source pressure entries do not close against the EOS. The baseline
curve rejects the separate unported pure sampling convention.

## Qualification and limits

`validate_mixture_spinodal.py` compares **1310 states from all 130 attempted
source calls**, comprising 108 fixed-pressure points and 22 complete mixture
curves. PR, SRK, RK (source SRK with RK alpha), van der Waals, translated PR/SRK,
binary/ternary mixtures and zero species are covered. Curves use a log-volume
step of 0.04. All point counts and all T/V/raw-P coordinates match at relative
2e-6 with absolute floors 2e-5 K, 1e-10 m³/mol and 2 Pa. Worst observed relative
differences are approximately 1.16e-8, 1.24e-8 and 1.44e-7, respectively.

Every source state is independently checked for a minimum stability eigenvalue
below 1e-6. All raw pressures and actual EOS pressures are retained; 22 endpoint
pressure warnings are recorded. There are no failed, unqualified or censored
source calls in this grid. The 54 physical vapor-selection differences remain
explicit diagnostics, separately from compatibility parity. Native tests also
check local pressure/eigenvalue closure, temperature stability crossings, scalar
inverse specifications, explicit seeds, both endpoint types, partial curves,
pure/zero-species/translated limits and invalid input.

This is local stability, not a global equilibrium or nucleation calculation.
Other activity-rule Hessians, other EOS families, critical degeneracies,
spinodal-property intersections, metastable isentropes and full raw trajectory
fidelity remain pending. Default pressure initialization requires a mechanical
seed below the fixed-composition pseudo-critical pressure; explicit local
seeds can address other roots without promising their global classification.

Reproduce with the exact-source build:

```text
python tools/rust_port/validate_mixture_spinodal.py
cargo test --test mixture_spinodal --test mixture_spinodal_source
python tools/rust_port/benchmark_mixture_spinodal.py --compiler <gfortran>
```

The [native executable benchmark](benchmark-mixture-spinodal.json) measures
complete baseline-compatible curves for PR/SRK C1+C2 and PR C1+C2+C3. Both
executables include initialization, every continuation step and endpoint
refinement, with model construction and Python/FFI excluded. All coordinates
and point counts are checked before timing; each returned quantity has a
separate timed checksum. Source, Cargo and harness hashes are pinned, and the
benchmark rejects edits during measurement. The timings apply to these curves,
not to the entire library or all spinodal modes.

Windows release medians, five samples after warm-up and ten curves per sample:

| Complete curve | Original | Native | Speedup |
|---|---:|---:|---:|
| PR C1+C2 | 2.458 ms | 1.313 ms | 1.87× |
| SRK C1+C2 | 2.585 ms | 1.273 ms | 2.03× |
| PR C1+C2+C3 | 3.311 ms | 1.786 ms | 1.85× |

All 175 preflight states match, with maximum absolute T/V/P differences of
2.10e-6 K, 5.31e-13 m³/mol and 0.112 Pa. Counts and all three output sums are
also checked during timing. The new stability path avoids critical third
derivatives, and mechanical volume bisection stops when floating-point bounds
can no longer move instead of spending the remaining 200 iterations there.
