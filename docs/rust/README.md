# ThermoPack Rust migration

Baseline: `d68c794` (the fork checked out on 2026-10-14). Original
Fortran, Python and C++ sources remain the numerical specification.
The Rust implementation is additive until every capability has passed
the acceptance criteria below. Apache-2.0 attribution is retained.

## Acceptance criteria

* Native Rust implementation of each model, solver and public capability.
* Differential tests against a library compiled from this exact baseline;
  a PyPI wheel alone does not establish parity with this fork.
* Explicit absolute and relative tolerances per quantity, with physical
  identities, derivative checks and challenging phase-boundary cases.
* Benchmarks that compare identical models, parameters, states, compiler
  optimization and outputs. Rust by itself does not guarantee a speedup.
* Reproducible builds, API documentation, examples and migration notes.

## Architecture and order

1. Inventory public interfaces, model families, databases and tests.
2. Embed source-derived fluid data, constants and interaction parameters;
   port numerical primitives and conventional cubic models.
3. Ideal caloric properties, analytic Helmholtz derivatives, fugacity,
   phase selection, saturation and stable flash calculations.
4. General solvers and remaining cubic mixing/alpha/translation variants.
5. CPA, PC-SAFT, SAFT-VR Mie/quantum, multiparameter, MBWR, CSP,
   Lee-Kesler, PETS, Lennard-Jones spline and specialized solid/hydrate models.
6. Envelopes, critical/spinodal curves, DFT, remaining derivatives and
   interfaces; parity and performance qualification of the entire surface.

Model instances own immutable parameters. Call-specific workspaces own
scratch arrays. Shared mutable active-model state is excluded from the native
API, allowing independent models to run concurrently. The numerical core uses
`f64`, SI molar units and explicit errors. Stable algorithms precede SIMD or
parallelization. Allocation-free prepared states and precomputed pair constants
are the first optimization targets. No fast-math or reduced precision by default.

## Traceable development

`main` remains the upstream baseline. Work lives on `codex/rust-port` with
small, ordered commits: inventory/design, numerical/data foundations, model
families, validation and measured optimization. Commit descriptions identify
source routines and validation. The generated inventory records every Python
method and Fortran source file; the coverage ledger distinguishes native,
validated, and pending capabilities. A wrapper calling Fortran is a reference
adapter, not a completed Rust port.

See [coverage.md](coverage.md), [inventory.json](inventory.json),
[numerics.md](numerics.md), [validation.md](validation.md),
[performance.md](performance.md) and [baseline-defects.md](baseline-defects.md).
The [phase equilibrium notes](phase-equilibrium.md) describe stability,
derivatives, TP result conventions and qualification limits.
Regenerate inventory with
`python tools/rust_port/inventory.py`.
