# Fundamental measure theory

`fmt::FmtModel` ports Rosenfeld, White Bear and White Bear II from
`src/fmt.f90`. The source model codes 1/2/3 are accepted by `from_code`.
The six weighted densities are ordered n0, n1, n2, n3, nv1, nv2. The last
two are the source's scalar representation of vector weighted densities in
one spatial dimension. General vector/tensor contractions are not added.
The result is reduced excess free-energy density Phi, with no kB*T factor.
The gradient has six entries and the symmetric Hessian has 36 entries.

```rust
use thermopack_rs::fmt::FmtModel;
let model = FmtModel::WhiteBearII;
let state = model.energy_density([0.1, 0.2, 0.3, 0.1, 0., 0.])?;
assert!(state.value > 0.);
# Ok::<(), thermopack_rs::Error>(())
```

The analytic path separates Phi into three coefficient functions of n3
and low-degree polynomials in the other five weighted densities. It computes
one second-order jet for each coefficient and assembles the complete Hessian
explicitly. This replaces 21 repeated two-direction hyperdual evaluations.
Both point and caller-owned grid APIs allocate no scratch storage. A grid
length mismatch is an error before any output write; a failed point leaves
the already-written prefix updated and the remaining suffix untouched.

## Source conventions and compatibility

Both paths retain `max(1-n3,1e-12)`, including the constant derivative branch
on equality. Finite negative weighted densities and n3 beyond close packing
are evaluated as algebraic extrapolations, not certified physical states.
Non-finite inputs/results are errors. White Bear formulas divide by n3 or
n3 squared: n3=0 is singular in the source and is an error here, including
the vacuum. No dilute-series regularization is claimed.

The source's `PowerHyperDualReal` evaluates the original reciprocal value,
but floors the base used in its derivative powers at magnitude 1e-15.
This affects coefficients with tiny denominators in dilute or near-close-
packing states. `energy_density` uses mathematical reciprocal derivatives;
`baseline_energy_density` follows the original 21 hyperdual evaluations,
operand order and derivative floor. Each has a corresponding `_grid` method.
The baseline is a numerical compatibility path, not the optimized derivative
path. The raw White Bear formulas can suffer cancellation at small n3 in
both paths; mathematical derivatives alone do not regularize those formulas.

The Windows Rust and source math runtimes also round negative integer powers
differently. The baseline calculates squares/cubes with compensated products
and corrects their reciprocals using fused multiply-add residuals. This
resolves the observed cancellation-amplified differences without changing
fixture tolerances or calling the Fortran library from production Rust.

## Qualification

The exact-source probe evaluates a 64-point grid for all three models and
all seven nonempty combinations of optional Phi/gradient/Hessian outputs.
All 1288 finite returns match all selected quantities with 2e-10 relative
and 1e-10 absolute tolerances. All 56 non-finite source returns (White Bear
variants at n3=0) are retained in metadata and produce native errors.
The grid includes sixteen n3 values from -0.1 to 1.1, n3 as small as 1e-16,
both sides of the denominator-clamp region, three independent density scales,
nonzero signed vector terms and uniform-sphere bulk states. Source rows are
transposed through a true multi-point Fortran call, checking grid layout.

Native tests independently difference the analytic gradient and Hessian at
nonsingular states, check exact Hessian symmetry, the Rosenfeld bulk excess
energy and the Carnahan-Starling limit of both White Bear variants, and cover
invalid inputs, model codes, grid lengths and the explicit derivative-floor
distinction. These qualify the FMT kernels; spatial convolutions, other DFT
functionals, SAFT density-functional integration and Python/C wrappers remain
pending. Hosted CI has not been observed here.

## Native grid measurements

The [recorded benchmark](benchmark-fmt.json) times grids of 100 points with
n3=0.05..0.842 and varying other weighted densities, including nonzero vector
terms. Each sample evaluates 1000 grids; five samples follow a warm-up.
All 43 fields are compared at 300 source points and 600 native points before
measurement and included in timed checksums. Both implementations write full
output grids into preallocated buffers. Python, FFI and buffer initialization
are outside the timers; kernel execution, output writes and checksum work
are inside them. Source/harness/Cargo fingerprints guard measurement inputs.

| Functional | Fortran ns/point | Analytic Rust ns/point | Baseline Rust ns/point | Source / analytic | Source / baseline |
|---|---:|---:|---:|---:|---:|
| Rosenfeld | 9022.20 | 130.59 | 971.90 | 69.09 | 9.28 |
| White Bear | 10545.08 | 139.84 | 1049.38 | 75.41 | 10.05 |
| White Bear II | 17483.80 | 156.16 | 1420.33 | 111.96 | 12.31 |

The analytic path avoids repeated hyperdual evaluations. The baseline retains
those evaluations but needs only the four active second-order fields, rather
than the source's general type with third-order fields and runtime order
branches. These measurements apply to this full-Hessian grid workload on the
recorded Windows GNU host, not value-only calls, singular states, a complete
DFT solver or the entire thermodynamics library.

```
python tools/rust_port/validate_fmt.py --compiler-directory /path/to/compiler/bin
cargo test --test fmt
python tools/rust_port/benchmark_fmt.py --compiler /path/to/gfortran
```
