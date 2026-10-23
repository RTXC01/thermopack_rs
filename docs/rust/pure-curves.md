# Pure-fluid saturation curves

The native cubic implementation now maps a selected component's coexistence
curve from a specified temperature or pressure to its critical endpoint.
`saturation_curve::singleCompSaturation` and
`pure_fluid_saturation_wrapper` are the source specification. Original Fortran
and parameter databases are unchanged.

## API and units

```rust
use thermopack_rs::{Cubic, CubicKind};
use thermopack_rs::pure_saturation_curve::{PureSaturationOptions, PureSaturationEnd};
# fn main() -> Result<(), thermopack_rs::Error> {
let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"])?;
let curve = model.pure_saturation_curve(1, PureSaturationOptions::default())?;
assert_eq!(curve.end, PureSaturationEnd::CriticalPoint);
# Ok(())
# }
```

The component index is zero-based in the existing model. All other species
are set to exact zero without reconstructing the model or replacing its
parameters. The default start is 0.1 MPa, with a requested pressure spacing
of 20 kPa and at most 100 points. `PureSaturationStart::Temperature` uses K;
`Pressure` uses Pa. Returned volumes are m³/mol and slopes are Pa/K.

The point limit redistributes the entire grid and includes both endpoints;
it does not return an early bounded prefix. A failed interior or endpoint
calculation does retain the accepted prefix with `NoConvergence` and a
failure. Invalid inputs or initial-state construction fail before a curve
exists. At least two points and a positive finite start are required. Linear
grids require a positive finite pressure spacing; logarithmic grids ignore it,
including zero or nonfinite values, as in the source.

`Cubic::pure_saturation_curve` ends at the critical point of the active EOS.
Its final phase volumes coincide at the pressure inflection, and the slope
is the fixed-volume pressure derivative there. Interior phase roots remain
distinct. Existing cubic pure saturation methods do not enforce caloric
temperature bounds; this selected-component curve has the same convention.

`Cubic::baseline_pure_saturation_curve` retains the fork's database-critical
T/P endpoint, raw cubic volume iteration and 1e-10 fugacity tolerance. It
returns `SourceCriticalParameters` at that endpoint. A modified model need
not be critical at its database endpoint. Complete raw iteration identity is
not claimed, especially at the nearly triple critical root.

## Grid and solver

For a pressure start, the grid is uniform in P; for a temperature start, it
is uniform in T. With `logarithmic=true`, it is uniform in the logarithm of
that same specified coordinate. As in the source, the number of points is
determined by the pressure interval even for temperature grids:

```text
linear:       min(round((Pc-P0)/requested_pressure_step) + 2, maximum_points)
logarithmic:  min(round((ln(Pc)-ln(P0))/0.02) + 2, maximum_points)
```

Positive half-integers round upward, matching Fortran NINT. The point cap
can enlarge the actual spacing. The logarithmic grid uses the source's
fixed 0.02 log-pressure request and ignores the pressure-spacing value.

Cold initialization uses existing congruent pure saturation. Each interior
point predicts from the previous Clapeyron slope, then corrects one logarithmic
variable with analytic fugacity derivatives, bounded steps and backtracking:

```text
F   = ln(phi_L) - ln(phi_V)
F_T = (H_V^res - H_L^res)/(R*T²)
F_P = (V_L - V_V)/(R*T)
```

The physical solver uses 1e-12 fugacity tolerance. A coarse predictor can
leave the two-root region and produce identical roots with F=0. Both solvers
reject this trivial interior solution and fall back to bracketed saturation.
Coarse 3-point and large-pressure-step source cases cover this regression.
Independent T/P perturbations check the analytic equations, including
translation and selected components in activity-based models.

## Qualification

`validate_pure_curve.py` attempts 203 isolated source calls. It obtains 201
finite complete curves containing **19,352 stored states**. Two VDW nitrogen
pressure starts terminate during source initialization; neither is counted as
successful parity. There are no censored calls, nonfinite returned states or
source pressure/fugacity closure warnings in this matrix.

The matrix includes PR/SRK/RK/VDW, methane/CO2/propane/nitrogen/argon/water,
translated PR/SRK, quantum H2/D2/He/Ne, and each selected component of a
ternary model with vdW/UNIFAC/HV1/HV2/NRTL mixing. It covers T/P starts,
linear/logarithmic grids, point limits 2/3/17/100/500 and coarse spacings.
Every source grid count, all T/P and all interior volumes pass relative
2e-7 with absolute floors 2e-6 K, 0.1 Pa and 1e-12 m³/mol.
Observed worst relative interior differences are 9.23e-12 in T, 8.50e-11
in P, 7.02e-12 in liquid volume and 7.96e-11 in vapor volume.

All **402 final volume comparisons** need an explicit critical-root
conditioning allowance; none is silently reported as passing the ordinary
volume tolerance. Their largest raw relative difference is about 1.68e-5.
The critical endpoint temperatures and pressures keep ordinary tolerances.

### Critical volume bound

The oracle records the source monic compressibility polynomial
`f(z)=z³+a*z²+b*z+c`, including active alpha/beta and volume translation.
With center `z0=-a/3`, the depressed polynomial is `w³+p*w+q`. Set
`eta=32*epsilon*sum(abs(monomials at z0))`. **Both source and native endpoint
volumes must have a polynomial backward error at most eta.**

For any such root, `|w|³ <= |p|*|w| + |q| + eta`. The largest positive root r
of `r³-|p|*r-|q|-eta=0` therefore bounds its displacement from z0. The two
volume roots differ by at most `2*r*R*T/P`. That source-derived radius is
stored per endpoint and recomputed by the fixture tests. It is used only
for final volumes, with separate raw differences retained in metadata.
It never changes the interior or T/P tolerances and does not claim identical
raw critical iterations.

The independent physical API is checked throughout the same 201-case matrix
for distinct interior phases, fugacity equality and coincident endpoint
volumes. Additional tests verify pressure closure, grid spacing, endpoint
first/second volume derivatives, selected-component equivalence and invalid
inputs. Fixture hashes, row counts and LF line endings are audited.

## Performance and remaining scope

`benchmark_pure_curve.py` compares four complete 100-point grids before
timing: PR methane, SRK CO2, translated PR water and quantum H2. All 400
T/P/volume states are checked, with the same explicit final-volume
conditioning. Timed checks compare counts and separate T/P/volume sums.
Cold initialization, complete continuation and both endpoint volumes are
included; model construction and Python/FFI are excluded. Rust additionally
calculates slopes. Reports retain individual samples and source/input hashes.
The native timing calls `baseline_pure_saturation_curve` so both implementations
return the source endpoint convention.

Local 2026-10-23 medians per complete curve, with five measured samples after
one warm-up and twenty curves per sample:

| Curve | Fortran ms | Rust ms | Speedup |
|---|---:|---:|---:|
| PR methane, linear P grid | 6.113 | 0.904 | 6.76× |
| SRK CO2, logarithmic T grid | 5.399 | 0.766 | 7.05× |
| Translated PR water, linear P grid | 6.417 | 0.967 | 6.63× |
| Quantum H2, logarithmic P grid | 6.681 | 1.292 | 5.17× |

These qualify these workloads, not the speed of the entire library. The
[machine-readable report](pure-curve-benchmark.json) includes compiler flags,
per-sample checksums, source identity, endpoint bounds and all native/harness/
Cargo input hashes.

Pure routing through the general envelope API, non-cubic pure curves,
pure-fluid azeotrope end-point searches, broader parameter modifications
and complete raw failure/iteration behavior remain pending. See the
[coverage ledger](coverage.md) and [overall progress estimate](progress.md).
