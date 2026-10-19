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
derivatives, TP/PH/PS result conventions and qualification limits.
The [alpha notes](alpha.md) describe temperature correlations, setters,
reference selection and the boundary between alpha and translated/quantum models.
The [volume translation notes](volume-translation.md) cover ci polynomials,
reference-aware constructors, physical volumes and caloric corrections.
The [beta notes](beta.md) cover quantum covolume factors, full QuantumCubic
critical-data overrides and the remaining cryogenic solver differences.
The [extended RR notes](rachford-rice.md) cover source-compatible negative-flash
bounds, independent phase fractions and one-step initialization.
The [symmetric-solver notes](symmetric-solver.md) cover the native upper
Bunch–Kaufman port, caller-owned buffers and pinned LAPACK qualification.
The [equilibrium-response notes](equilibrium-response.md) cover analytic
volume/caloric T/P derivatives, phase redistribution and source Cp conventions.
The [sound notes](sound.md) cover fluid SV acoustics, extensive TV inputs,
physical/source conventions and native-loop benchmark results.
The [isochoric notes](isochoric.md) cover direct TV caloric properties, coupled
UV/SV equations, pure latent quality and qualified source comparison scope.
The [PeTS notes](pets.md) cover the pure-fluid Helmholtz/TV kernel, potential,
LJ correction, source-order virials, temperature preparation and remaining
higher-level interfaces.
The [PeTS equilibrium notes](pets-equilibrium.md) cover source-compatible volume
iteration, strict pressure roots, critical scaling and pure coexistence.
The [PeTS flash notes](pets-flash.md) cover pure TP/PH/PS, specified quality,
equilibrium caloric properties and explicit endpoint comparison conventions.
The [PeTS fixed-volume notes](pets-isochoric.md) cover equilibrium TV states,
analytic Cv with phase redistribution, pure UV/SV and their qualification.
The [PeTS root measurements](pets-root-performance.md) distinguish fresh,
prepared and source-iteration paths and identify branch-discovery overhead.
The [MBWR kernel notes](mbwr.md) cover all stored Bender/MBWR32 fits, pressure,
residual derivatives, source summation fidelity and the stable integral path.
The [MBWR state-function notes](mbwr-properties.md) cover residual G/H/S/Z/lnphi,
physical/source TP responses, entropy signs and raw-rounding test context.
The [MBWR equilibrium notes](mbwr-equilibrium.md) cover bounded pressure branches,
derivative-hierarchy spinodals, stored critical metadata and analytic critical points.
The [PC-SAFT notes](pc-saft.md) cover PC-SAFT and sPC-SAFT Helmholtz kernels,
Lorentz-Berthelot combining rules, bounded pressure-volume branches and pure coexistence.
The [SAFT-VR Mie notes](saft-vrmie.md) cover SAFT-VR Mie and SAFT-VRQ Mie Helmholtz kernels,
Barker-Henderson numerical quadrature, BMCSL reference, exact contact RDF derivation and mixture scaling.
The [CPA notes](cpa.md) cover the Cubic-Plus-Association equation of state, Wertheim multi-site association,
simplified contact RDF, quadratic closed branches and non-associating component mixing conventions.
The [Excess Gibbs notes](excess-gibbs.md) cover Huron-Vidal (HV1, HV2) and NRTL activity model mixing rules,
infinite-pressure excess Helmholtz energy, exact analytical derivatives and binary parameter databases.
The [Multiparameter notes](multiparameter.md) cover reference-quality NIST multiparameter equations of state,
fundamental reduced Helmholtz energy formulations, exact first/second derivative jets, and ancillary density equations.
The [Lee-Kesler notes](lee-kesler.md) cover the Lee-Kesler corresponding states equation of state,
Plöcker critical mixing rules, analytical compositional derivatives, and BWR generalized fluids.
Regenerate inventory with
`python tools/rust_port/inventory.py`.
