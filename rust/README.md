# Native Rust ThermoPack

Experimental native port; see [coverage](../docs/rust/coverage.md) before use.
The original Fortran/Python implementation is retained as the reference.

`Mbwr` provides MBWR-19/32 pressure and residual Helmholtz kernels with all
65 stored fits. It offers stable integrals and explicit source-order methods;
see [MBWR scope and qualification](../docs/rust/mbwr.md).
Its [residual TP state functions](../docs/rust/mbwr-properties.md) now provide
G/H/S/Z/lnphi and local TP responses, with explicit source-compatible methods.
[Bounded pressure branches and critical points](../docs/rust/mbwr-equilibrium.md)
and [source density iteration](../docs/rust/mbwr-density.md) have separate APIs
and qualifications. [Default initialization](../docs/rust/mbwr-initialization.md)
also ports the source coarse critical search. Full caloric and flash/CSP integration
remain pending.

```
cargo test --all-targets
cargo run --release --example cubic
cargo bench --bench cubic
```

The crate has no runtime dependencies or Fortran calls. `Cubic::new` embeds the
fork's fluid and default PR/SRK interaction database. `prepare` validates and
caches T/composition parameters. Pressure scans, roots, residual properties and
fugacity into caller-owned buffers allocate no memory after preparation.
`Pets::new()` provides the pure-argon PeTS Helmholtz/TV kernel; its
`prepare(T)` caches temperature terms for repeated densities. PeTS accepts
total volume and mole number, and returns extensive caloric quantities. See
[PeTS scope and qualification](../docs/rust/pets.md) for derivatives, sound,
potentials and virials. [PeTS equilibrium](../docs/rust/pets-equilibrium.md)
adds reusable pressure-root branches, critical points, pure saturation and
explicit source last-iterate compatibility. [PeTS flashes](../docs/rust/pets-flash.md)
add pure TP/PH/PS, specified saturation quality and molar equilibrium properties.
[PeTS fixed-volume equilibrium](../docs/rust/pets-isochoric.md) adds analytic
equilibrium Cv and pure UV/SV. Remaining inherited interfaces are pending.

Temperatures: K; pressures: Pa; volumes: m3/mol; energies: J/mol;
compositions: normalized mole fractions. Root selection is explicit. A
minimum-Gibbs root at fixed composition is not a mixture stability test.

```rust
use thermopack_rs::{Cubic, CubicKind, Phase};

let model = Cubic::new(CubicKind::PengRobinson, &["CO2", "C1"])?;
let state = model.prepare(300.0, &[0.8, 0.2])?;
let volume = state.volume(5e6, Phase::Vapor)?;
let pressure = state.pressure(volume)?;
let enthalpy = state.enthalpy(volume)?;
let mut ln_phi = [0.0; 2];
state.ln_fugacity_into(volume, &mut ln_phi)?;
```

Pure cubic saturation: `model.saturation_pressure(temperature)` and
`model.saturation_temperature(pressure)` for a model with one component.

Mixture bubble/dew pressure and temperature solves and warm branch refinement
are available on `Cubic`; see [mixture saturation](../docs/rust/mixture-saturation.md)
for fixed-phase composition, analytic slopes and the qualification limits.
Returns coexistence temperature/pressure, phase volumes and Clapeyron slope.
`model.tp_flash(T,P,z)` performs VLE stability checks and a TP split.
`model.ph_flash(P,z,H)` and `model.ps_flash(P,z,S)` recover temperature and
phase amounts from molar enthalpy or entropy. Use
`model.equilibrium_properties(&flash)` for overall V/H/S/U and equilibrium Cp,
including phase redistribution. Caloric options set temperature bounds and
initial guesses; defaults are 80 to 999 K. `pseudo_critical(z)` is a
fixed-composition estimate, not the true critical point of a mixture.

Mixture [spinodal points and curves](../docs/rust/mixture-spinodal.md) use the
Helmholtz stability eigenvalue, with explicit source-compatible and physical
branch/pressure conventions. Partial curves report their point-limit status.
Pure grids and [property searches](../docs/rust/spinodal-properties.md) include
fresh/source-compatible queries, reusable prepared searches and physical
crossings of all six scalar property constraints on pure or mixture curves.
`model.ideal_gas_mut(i)` configures owned ideal Cp parameters and H/S reference
states for component i (zero based). A Cp setter preserves integration constants;
reference setters recalculate them as in the source. Other models and the global
database remain independent. `components()` exposes active component metadata, while
`ideal_gas(i)` exposes the active caloric parameters.
`Cubic::with_alpha` selects an alpha correlation and parameter reference;
`set_alpha(i, correlation, coefficients)` sets one component explicitly.
Effective parameters and fallback status are available through `alpha(i)`.
See [alpha conventions](../docs/rust/alpha.md), including source initialization
quirks and the separate scope of translated/quantum constructor behavior.
`Cubic::with_options` selects alpha, ci and binary parameter references
and optional volume translation; tcPR references force TWU and translation as
in the source constructor. Component ci parameters and enable/disable controls
are owned by each model. See [volume translation](../docs/rust/volume-translation.md).
`set_beta(i, Beta::new(...))` configures covolume factors and their derivatives.
The `QuantumCubic` reference forces quantum beta/TWU/translation and source
critical-data overrides. Model component records now contain effective critical
data, independently of the global database. See [beta and QuantumCubic](../docs/rust/beta.md)
for qualified properties and pending cryogenic flash compatibility.
`Cubic::quantum(names)` also matches `qcubic`'s combined parameter priorities
and 2 K minimum. Model-local bounds are available through `temperature_bounds`
and `set_temperature_bounds`; `caloric_options()` inherits them for per-call
customization.
`state.fugacity_derivatives_into` provides analytic T/P/n derivatives into
caller-owned buffers. See [phase equilibrium conventions](../docs/rust/phase-equilibrium.md),
including the Single phase fraction sentinel and a known raw ordering mismatch.
The sound-speed and heat-capacity API uses physical derivatives;
`baseline_enthalpy_derivative_tp` explicitly preserves a vdW derivative defect
and the ideal-H2 Cp convention in this fork. `Component::ideal` exposes the
source Cp and a separate `enthalpy_dt` for the derivative of its H values.
See [numerical differences](../docs/rust/baseline-defects.md).

See [validation](../docs/rust/validation.md),
[benchmarks](../docs/rust/performance.md) and the
[migration design](../docs/rust/README.md). Remaining SAFT capabilities, mu-T flashes, phase envelopes,
solids/hydrates and existing language-interface compatibility remain pending.

Regenerate data with `python tools/rust_port/generate_data.py`.

The native upper symmetric-indefinite solver is derived from pinned LAPACK
DSYTF2/DSYTRS2 and reference BLAS algorithms. Its BSD-3-Clause notices and terms
are retained in [LAPACK-LICENSE](LAPACK-LICENSE); ThermoPack-derived code retains
Apache-2.0 attribution. See [symmetric-solver notes](../docs/rust/symmetric-solver.md)
for workspace requirements and qualification scope.

`equilibrium_response(&flash)` returns physical analytic Cp and volume T/P
derivatives, including phase redistribution. `baseline_equilibrium_response`
retains the fork's reported Cp conventions. Pure latent states have no
independent TP response. See [equilibrium-response notes](../docs/rust/equilibrium-response.md)
for derived U/Cv/isentropic identities and qualification limits.

`equilibrium_speed_of_sound(&flash)` includes thermal/mechanical/compositional
phase redistribution, including pure coexistence. `speed_of_sound_tv(T,V,n)`
accepts total volume and mole numbers for a single-phase TV evaluation.
The corresponding `baseline_...` methods preserve source caloric and fallback
conventions; `..._with_phases` supports fluid-root overrides. See
[sound qualification and benchmarks](../docs/rust/sound.md).

`state.caloric_tv(v)` evaluates molar P, U, S and physical Cv, including
homogeneous negative-pressure states. `uv_flash(z,U,V)` and `sv_flash(z,S,V)`
recover equilibrium with pure latent quality and absent feed components.
`volume_flash_options()` supplies model-local temperature bounds and optional
T/P guesses. See [isochoric qualification](../docs/rust/isochoric.md) for the
409 finite source flash comparisons, 141 TV comparisons, canonical single-phase
output conventions and remaining source trajectory/cryogenic limits.
Native-loop [UV/SV measurements](../docs/rust/isochoric.md#native-loop-measurements)
include common source-generated specifications and per-state parity checks;
the qualified default/tcPR workloads show UV speedups of about 27x and SV
speedups of 2.52x/2.34x. These are workload-specific measurements.

The cubic model also carries the curve and the limit families of the source:
`isolines::isotherm`, `isolines::isobar`, `isolines::isenthalp` and
`isolines::isentrope` map the iso-lines, `spinodal::spinodal_point` reports
the point of the meta stability limit at a given pressure,
`joule_thompson::inversion_point` reports the inversion temperature at a
given molar volume, `critical::critical_point` reports a local pure or mixture
critical point, and `Cubic::new_pseudo` builds the equation of state from
user supplied critical data. `state.virial_coefficients` reports the second
and the third virial coefficient. The solid carbon dioxide model carries the
solid-fluid line in `solid_saturation::solid_co2_saturation`, and the
`solid_co2_with_fluid_reference` constructor places the solid on the liquid at
the triple point of the source. See [iso-line notes](../docs/rust/isolines.md),
[stability limit notes](../docs/rust/spinodal.md),
[pseudo component notes](../docs/rust/pseudo-components.md) and
[solid notes](../docs/rust/solid.md).

The splined Lennard-Jones model with the Barker-Henderson perturbation theory
lives in `ljs`: `LjsBh::new` builds the argon record of the source database,
`LjsBh::with_parameters` replaces the segment diameter and the well depth,
`LjsBh::set_control` switches the terms of the model, and the model reports
the reduced perturbation terms, the Barker-Henderson diameter, the spline
polynomials, the reduced residual Helmholtz energy and its derivatives, the
pressure of a state, the reduced residual chemical potential, the caloric state
and the reduced interaction potential. See
[Lennard-Jones splined notes](../docs/rust/ljs.md).
