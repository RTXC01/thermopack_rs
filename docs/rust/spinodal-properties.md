# Pure spinodal grids and property searches

The cubic port now implements the fork's two 100-point pure spinodal grids,
pure property-intersection selection and a reusable search workspace. A separate
physical search finds every bracketed property crossing on a supplied pure or
mixture stability curve.

`SpinodalProperty` covers entropy (J/mol/K), volume (m³/mol), enthalpy (J/mol),
temperature (K), pressure (Pa) and internal energy (J/mol). All compositions are
normalized mole fractions; the model owns the ideal caloric reference levels.

## Grids and physical roots

`baseline_pure_spinodal_curve(model, z, P0, liquid_T, vapor_T)` returns 100
equally spaced liquid temperatures from the starting boundary to stored Tc,
followed by the 100 vapor points in reverse order. Optional temperatures change
the branch starts before sampling. Without a vapor temperature the vapor end
is at P0, as used internally by the source property searches.

`baseline_trace_spinodal` now routes pure compositions through this grid; its
default vapor minimum is 0.6 times the composition-weighted stored critical
temperatures, matching the public Python mapper. Pure sampling ignores the
log-volume increment. An explicit smaller native point limit returns a prefix
with `PointLimit` status. Physical `trace_spinodal` keeps volume continuation.

Entries 99 and 100 retain stored Tc/Pc and the original phase-volume convention.
`PreparedState::baseline_volume` ports `cbCalcZfac` gflag=1 and its bounded
Halley/relaxation iteration. It returns the finite last iterate without the
physical root filters of `volume`; status/trajectory identity and general
parameter qualification are not implied. Near a multiple root, coefficient
rounding changes the last volume more than a usual relative tolerance allows.

The mechanical vapor-volume bracket now starts at the fixed-composition
critical volume, where the subcritical pressure slope is positive. The old
bracket started at the covolume: near Tc it could skip both extrema and converge
to the repulsive pole. The van der Waals propane near-critical grid exposed
this defect. The new grid tests check every noncritical state against the local
stability condition and reject that spurious pole.

## Intersection APIs

`spinodal_intersections(model, z, curve, property, target)` uses a traced curve
to bracket and refine crossings. Every returned state closes both the property
and minimum Helmholtz eigenvalue equations. Temperature is refined at fixed
volume, with interpolation along the bracket for the initial guess. This works
with pure and mixture curves and returns all bracketed crossings. An incomplete
or coarse curve cannot establish the absence of other or tangent crossings.

`baseline_pure_spinodal_intersection` preserves the fork's pure branch selection:
liquid, vapor, or both (`None`), with vapor precedence when both have a crossing.
It retains the last monotonic segment after a property extremum and uses the
source King/Pegasus refinement and its 1e-8 objective tolerance. Grid-target
equality is stabilized within relative 2e-8 to avoid changing the reporting
interval solely through independent grid iteration roundoff.

The source refinement predicts T from the trial V, corrects V at fixed T to
evaluate the property, but returns the **uncorrected trial V**. The compatibility
API preserves this reporting convention. A source success therefore does not
ensure that its returned state is on the spinodal or closes the property.
The source also normalizes its volume objective by 1, giving an absolute
1e-8 m³/mol objective tolerance. The physical search scales volume by its target.

## Reusing a search

`Cubic::prepare_pure_spinodal_search` owns a grid and all six property arrays.
Its borrowed immutable model prevents parameter or reference-state changes
while it is alive. Each `locate` query reuses those data, including source branch
and coordinate-reporting conventions. Construction is explicit; one query does
not silently pay for rebuilding the grid on every later request.

```rust
use thermopack_rs::{Cubic, CubicKind, Phase};
use thermopack_rs::spinodal_properties::SpinodalProperty;
# fn main() -> Result<(), thermopack_rs::Error> {
let model = Cubic::new(CubicKind::PengRobinson, &["C1"])?;
let search = model.prepare_pure_spinodal_search(&[1.0], 1e6, None)?;
let point = search.locate(SpinodalProperty::Temperature, 180.0,
                          Some(Phase::Liquid))?.expect("crossing");
assert!((point.temperature - 180.0).abs() < 1e-5);
# Ok(())
# }
```

## Qualification

`validate_spinodal_properties.py` compares the exact-source library:

* Twenty full grids, **4000 states**: PR/SRK/RK/VDW methane, CO₂, propane and
  nitrogen, plus translated PR/SRK methane and CO₂; P0=1 MPa and vapor end 0.6 Tc.
  Every returned coordinate is retained, including both critical entries.
* **276 successful intersection returns** from 288 attempts, all six properties,
  liquid/vapor/both selection, crossing and no-crossing statuses, and targets
  constructed from source grid states. Twelve van der Waals caloric searches
  abort the source process through invalid-volume evaluation; their exact
  targets, exit logs and native responses are retained. Native diagnostics are
  four errors, four no-crossing returns and four finite returns, without a
  parity claim for calls that never return from the source.
* There are no censored calls. At 145 successful returned coordinates the source
  minimum eigenvalue exceeds 1e-6 in magnitude; these raw reporting warnings are
  recorded without replacing or discarding the source result.

Nominal T/V/P comparison is relative 2e-6 with absolute floors 2e-5 K,
1e-11 m³/mol and 2 Pa. Only the critical entries 99/100 add the volume radius
`(R*T/P)*(64*machine_epsilon)^(1/3)`, reflecting a monic cubic multiple root's
coefficient-rounding sensitivity. Thirty-one comparisons require that radius;
the largest critical-volume relative difference is 1.68e-5. Metadata records
each such comparison, including both values and its radius. Ordinary grid
states and all property outputs keep the nominal tolerance.

Rust tests compare fresh and cached queries independently with every source
return. Physical tests check all six constraints, pressure and stability closure,
both pure/mixture curves, no crossing, invalid input, ordinary baseline volume
roots and the narrow vapor bracket. Qualification uses one-component model
instances and default caloric references; raw embedded-pure/mixed-model grids,
other alpha/beta/reference settings, all failure trajectories and other EOS
families remain unqualified. Mixture source-order property-intersection selection,
property extrema APIs and metastable isentropes remain pending.

Reproduce after building the exact-source reference:

```text
python tools/rust_port/validate_spinodal_properties.py
cargo test --test spinodal_properties
python tools/rust_port/benchmark_spinodal_properties.py --compiler <gfortran>
```

The [benchmark report](benchmark-spinodal-properties.json) compares temperature
intersections on PR methane and SRK nitrogen liquid spinodals: 25 targets from
0.94 to 0.98 Tc, 100 queries per sample, five samples after warm-up. Fresh native
queries include the complete grid construction. Prepared queries exclude the
one-time grid/property preparation, which is timed separately as `setup_*`.
The original public routine remaps both branches for each query. Prepared gains
are amortized workflow gains, not equal work for a single fresh call. Model
construction, Python and FFI are outside every timer. All 50 source states and
100 native preflight results, plus separate T/V/P timed checksums, are compared.
All source, Cargo and harness hashes are pinned; edits during measurement fail.

Windows release medians, five samples after warm-up:

| Liquid temperature query | Source | Native fresh | Native prepared | Preparation |
|---|---:|---:|---:|---:|
| PR C1 | 7.890 ms | 0.3625 ms | 18.16 µs | 0.5456 ms |
| SRK N2 | 2.872 ms | 0.3394 ms | 17.77 µs | 0.4264 ms |

Fresh queries are **8.46–21.77×** faster in these workloads. Prepared queries
are **161.6–434.5×** faster after the separately measured setup; setup amortizes
after two queries versus fresh native searches here. Every preflight output
matches: maximum absolute T/V/P differences are 2.85e-14 K, 2.27e-15 m³/mol and
2.24e-8 Pa. These are workload-specific API measurements, not whole-library or
all-property performance claims.
