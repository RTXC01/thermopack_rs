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

Regenerate data with `python tools/rust_port/generate_data.py`.
