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

Pure cubic saturation: `model.saturation_pressure(temperature)` for a model
with one component. Returns coexistence pressure, phase volumes and Clapeyron
slope. The sound-speed and heat-capacity API uses physical derivatives;
`baseline_enthalpy_derivative_tp` explicitly preserves a vdW derivative defect
in this fork. See [numerical differences](../docs/rust/baseline-defects.md).

See [validation](../docs/rust/validation.md),
[benchmarks](../docs/rust/performance.md) and the
[migration design](../docs/rust/README.md). SAFT, flashes, phase envelopes,
solids/hydrates and existing language-interface compatibility remain pending.

Regenerate data with `python tools/rust_port/generate_data.py`.
