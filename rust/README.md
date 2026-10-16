# Native Rust ThermoPack

Experimental native port; see [coverage](../docs/rust/coverage.md) before use.
The original Fortran/Python implementation is retained as the reference.

```
cargo test --all-targets
cargo run --release --example cubic
cargo bench --bench cubic
```

The crate has no runtime dependencies or Fortran calls. `Cubic::new` embeds the
fork's fluid and default PR/SRK interaction database. `prepare` validates and
caches T/composition parameters. Pressure scans, roots, residual properties and
fugacity into caller-owned buffers allocate no memory after preparation.

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
Returns coexistence temperature/pressure, phase volumes and Clapeyron slope.
`model.tp_flash(T,P,z)` performs VLE stability checks and a TP split.
`model.ph_flash(P,z,H)` and `model.ps_flash(P,z,S)` recover temperature and
phase amounts from molar enthalpy or entropy. Use
`model.equilibrium_properties(&flash)` for overall V/H/S/U and equilibrium Cp,
including phase redistribution. Caloric options set temperature bounds and
initial guesses; defaults are 80 to 999 K. `pseudo_critical(z)` is a
fixed-composition estimate, not the true critical point of a mixture.
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
[migration design](../docs/rust/README.md). SAFT, UV/SV/mu-T flashes, phase envelopes,
solids/hydrates and existing language-interface compatibility remain pending.

Regenerate data with `python tools/rust_port/generate_data.py`.
