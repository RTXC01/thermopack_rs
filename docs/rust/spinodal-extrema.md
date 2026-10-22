# Pure spinodal property extrema

`baseline_pure_spinodal_extrema(model, z, property, P0)` ports
`spinodal::locate_spinodal_prop_min_max_pure_fluid`. The reusable workspace
`Cubic::prepare_pure_spinodal_search` also exposes `extrema(property)`, so
repeated queries reuse the same grid and all six stored property arrays.

The six properties are entropy, volume, enthalpy, temperature, pressure and
internal energy, with the SI molar units documented in the
[spinodal property notes](spinodal-properties.md). Only pure compositions are
accepted. Positive initial pressure is required; reported pressures can be
negative. These APIs preserve the source's candidate selection and reporting
conventions. They do not promise global extrema over the whole spinodal.

## Selection and reporting

For each of the two 100-point source grids, the routine selects the last
discrete turning segment. If it finds one, it refines the zero property
derivative with respect to pressure, using the existing source King/Pegasus
refiner and its absolute 1e-8 objective tolerance. If a branch is monotone,
the **first point of that branch** is its candidate. The endpoint at Tc is
not included as an alternative. Thus even the temperature maximum can be
below the critical temperature.

The two selected candidates are sorted by their evaluated property values.
Refinement retains the original uncorrected trial volume, as in the property
intersection API. The source then evaluates `P_min` using **T_max and V_min**;
`P_max` uses T_max and V_max. The native result reproduces that mixed-coordinate
minimum pressure. Recompute EOS pressure at each returned T/V when actual
state pressure is needed. Physical property-intersection and stability APIs
remain available separately.

`BaselineSpinodalExtrema::ordering_ambiguous` identifies pressure candidates
whose property values are indistinguishable within 2 Pa plus relative 1e-7.
Both monotone branches start at P0; their tiny pressure residuals determine
which one the original strict comparison calls minimum. Independent solver
roundoff can reverse those labels. The native comparison retains its strict
ordering and sets the ambiguity flag. Candidate coordinates are qualified in
either order in this case. Raw bit-for-bit label parity is not established;
callers must not interpret an ambiguous ordering as a resolved inequality.

The returned pressure convention depends on which candidate is named maximum.
For ambiguous cases the fixture retains the source ordering and additionally
evaluates both pressure formulas for the alternative assignment using the
source EOS. This checks the complete candidate/reporting calculation without
silently discarding the 50 observed primary-order differences.

## Reuse and optimization

```rust
use thermopack_rs::{Cubic, CubicKind};
use thermopack_rs::spinodal_properties::SpinodalProperty;
# fn main() -> Result<(), thermopack_rs::Error> {
let model = Cubic::new(CubicKind::PengRobinson, &["C1"])?;
let workspace = model.prepare_pure_spinodal_search(&[1.0], 1e6, None)?;
let entropy = workspace.extrema(SpinodalProperty::Entropy)?;
let energy = workspace.extrema(SpinodalProperty::InternalEnergy)?;
assert!(!entropy.ordering_ambiguous && !energy.ordering_ambiguous);
# Ok(())
# }
```

The borrowed workspace prevents model mutation during reuse. Each selection
scan uses fixed stack arrays, with no temporary heap allocation for its
indices, branch points or property-value arrays. Cached monotone candidates
reuse their property values. Both reported pressures use one prepared
temperature/composition state, avoiding duplicate mixture preparation.
Refinement and state preparation still allocate their normal workspaces;
the whole query is not advertised as allocation-free.

## Exact-source qualification

`tools/rust_port/validate_spinodal_extrema.py` makes 504 isolated attempts:
four classical cubics, nine fluids (C1, CO2, C3, N2, H2, H2O, MEOH, ETOH, MEG),
two starting pressures (0.1 and 1 MPa), all six properties, plus translated
PR/SRK methane, carbon dioxide and water.

| Source outcome | Count | Native evidence |
|---|---:|---|
| Unambiguous finite returns | 379 | All six ordered T/V/P coordinates match |
| Pressure-tied candidate pairs | 83 | Both candidates and both reporting assignments checked; 50 primary label differences |
| Source failures | 42 | 36 nonzero statuses and 6 aborts retained |
| Censored / unqualified finite returns | 0 / 0 | Eight-second per-call deadline |

Tolerances are relative 2e-6 with absolute floors of 2e-5 K, 1e-11 m³/mol and
2 Pa. The largest relative coordinate difference, after explicitly allowing
the pressure-tied candidate permutation, is 9.90e-7. The metadata preserves
all failures, ambiguity assignments, exact library/source identities and
fixture SHA-256. The comparison does not qualify those permutations as full
raw ordered-result parity.

There are 613 coordinate warnings: 462 pressure-coordinate mismatches and
151 off-spinodal reported volumes. These remain separate from coordinate
matching. Tests run fresh and cached queries for every finite return, verify
the pressure formula independently, and confirm caloric stationarity by
differencing the property along corrected spinodal volumes. They also cover
the monotone first-point convention and invalid inputs. Other model families,
broad degeneracy, original solver trajectories and pressure-tie raw ordering
remain pending.

## Performance

The native-executable benchmark cycles entropy, enthalpy and internal energy
for PR methane and SRK nitrogen at P0=1 MPa. It compares all six reported
coordinates of six source extrema pairs and twelve fresh/cached native pairs
before timing. Each timed sample checks six independent sums over 100 queries.
One-time cached setup is measured separately and excluded from cached query
timing. Model construction and Python/FFI are excluded throughout.

`spinodal-extrema-benchmark.json` records five samples after one warm-up,
compiler/settings and hashes of native sources, harness inputs and Cargo
files. Cached speedups are amortized workflow measurements, and apply only to
the stated workloads; pressure-order ties are not part of the timing workload.

Windows x86-64 medians, Rust 1.99 release with thin LTO and one codegen unit
versus gfortran 14.3 at O3:

| Model | Source | Native fresh | Native cached | One-time setup |
|---|---:|---:|---:|---:|
| PR methane | 8.336 ms | 0.4754 ms | 35.16 µs | 0.5663 ms |
| SRK nitrogen | 3.208 ms | 0.4636 ms | 33.40 µs | 0.4298 ms |

Fresh speedups are 6.92–17.54×; cached speedups are 96.07–237.10× after the
separately reported setup. The largest preflight coordinate discrepancies are
1.33e-7 K, 1.81e-12 m³/mol and 0.0378 Pa. Every preflight coordinate and timed
sum passes the documented tolerance. These are scoped API measurements,
not whole-library or equal-work one-shot cached speedups.
