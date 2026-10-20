# Ice Ih Gibbs energy

`SolidWater` ports the twelve public thermodynamic functions and `sho2_init`
reference shift from `src/h2o_gibbs.f90`. The source implements the
Feistel–Wagner ice Ih Gibbs equation with two complex temperature terms and
pressure-dependent polynomial coefficients. The native kernel uses a small
private complex type, principal logarithms and scaled real/complex division;
it adds no runtime dependency or unsafe code.

```rust
use thermopack_rs::SolidWater;
let mut ice = SolidWater::new();
ice.initialize(100.0, -10000.0)?;
let isotherm = ice.prepare(250.0)?;
let state = isotherm.state(1e5)?;
assert!(state.volume > 0.0);
# Ok::<(), thermopack_rs::Error>(())
```

`initialize(liquid_entropy, liquid_gibbs)` takes molar liquid reference values
at the source's fixed 273.152519 K / 101325 Pa melting reference. It first
sets the solid entropy to S_liquid - DeltaS_melt, then aligns Gibbs energy with
the liquid, preserving the source update order. These reference conditions
differ from the equation's 273.160 K / 611.657 Pa reduction constants. Each model
owns its two reference offsets; shifts do not mutate another model or global
EOS. Invalid shifts leave the model unchanged.

`prepare(T)` caches the temperature functions, their first two derivatives and
complex logarithms. Pressure sweeps allocate no scratch storage. `gibbs(P)` on
the prepared state returns G, G_T, G_P, G_TT, G_TP and G_PP. `state(P)` shares
that jet across volume, entropy, enthalpy, internal/Helmholtz energy and Cp.
Scalar convenience methods are also available on the model. Quantities use
K, Pa and molar SI units; G/H/U/A use J/mol and entropy/Cp use J/(mol K).
Zero temperature is supported. Negative/non-finite T and non-finite P or
computed Gibbs jets are errors. Finite negative pressures and extrapolated
temperatures are evaluated without asserting that ice Ih is the stable phase.

## Heat-capacity convention

`sh2o_d2gdT2` already converts the Gibbs derivative to molar units, but the
original `sh2o_heat_capacity` multiplies `-T*G_TT` by the molar mass again.
`heat_capacity` and `SolidWaterState::heat_capacity_p` return the physical
derivative Cp=-T*G_TT. `baseline_heat_capacity` and the state's
`baseline_heat_capacity_p` preserve the source's extra factor 0.01801528.
Reference comparisons explicitly select this baseline quantity. The physical
quantity is checked independently against dH/dT at fixed pressure.

The source water sound-speed routine is commented out; no source-compatible
sound function is claimed here. Solid saturation, fluid-reference initialization
from an EOS, solid mixtures/phase selection and three-phase/hydrate interfaces
remain pending.

## Qualification

All 780 attempted source states are finite and match twelve outputs: G and its
five first/second derivatives, V, S, H, U, A and source-scaled Cp. The matrix
uses thirteen temperatures from 0 to 500 K, ten pressures from -1 to 1000 MPa,
and six reference scenarios (raw, three shifted references, repeated shift and
two different sequential shifts). Each scenario starts in a fresh oracle
process to avoid module-reference leakage. Extrapolated points qualify the
equation's numerical behavior, not phase stability.

The relative tolerance is 1e-11 with quantity-specific absolute floors in
fixture metadata. Independent native tests verify all Gibbs derivatives by
finite differences, Maxwell/caloric identities, physical Cp, reference
alignment/idempotence, model independence, prepared/scalar equivalence, zero T
and invalid input/transactional behavior. Hosted CI has not been observed here.

## Native-loop state measurements

The [recorded benchmark](benchmark-solid-water.json) compares native executable
loops at 250 K and 100 pressures from 0.1 to 99.1 MPa. Both implementations
produce all twelve source outputs, including source-scaled Cp; every output is
checked in 100 source/200 native preflight states and included in timed checksums.
Initialization, Python and FFI are outside the timers. Each sample contains
500000 states; five samples follow a warm-up. Source/harness/Cargo fingerprints
must remain unchanged during measurement.

| Implementation | Median ns/state | Source time / native time |
|---|---:|---:|
| Original twelve Fortran functions | 3691.67 | 1 |
| Rust shared Gibbs jet, fresh temperature preparation | 848.07 | 4.35 |
| Rust shared Gibbs jet, reused temperature preparation | 40.44 | 91.28 |

The source has separate public functions that repeat temperature work. The
native state API shares it across properties; prepared mode additionally moves
the fixed-temperature work outside the timer. This measures the full-state
isotherm workload on this Windows GNU host, not individual scalar functions,
arbitrary temperature trajectories, phase solvers or whole-library speed.

```
python tools/rust_port/validate_solid_water.py --compiler-directory /path/to/compiler/bin
cargo test --test solid_water
python tools/rust_port/benchmark_solid_water.py --compiler /path/to/gfortran
```
