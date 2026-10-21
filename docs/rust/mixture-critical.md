# Mixture critical points

`critical::critical_point(model, z)` now solves true local cubic mixture
critical points as well as the existing pure pressure inflections. The input
is a normalized mole-fraction vector. Results contain T (K), molar V (m³/mol)
and P (Pa). This is a local solve, not a complete search of a critical locus.

```rust
use thermopack_rs::{Cubic, CubicKind, critical::critical_point};
# fn main() -> Result<(), thermopack_rs::Error> {
let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"])?;
let point = critical_point(&model, &[0.5, 0.5])?;
assert!((point.temperature - 265.85774).abs() < 1e-4);
let conditions = model.mixture_critical_conditions(
    point.temperature, point.volume, &[0.5, 0.5])?;
assert!(conditions.eigenvalue.abs() < 2e-9);
# Ok(())
# }
```

`mixture_critical::critical_point_with_guess` accepts explicit positive T and V
seeds. Default initialization starts near the pseudo-critical temperature and
four covolumes, with alternate seeds after failed local solves. Explicit seeds
select a local basin; multiple critical points are not enumerated.

## Equations and derivatives

At fixed temperature and **total** volume, let H be the mole Hessian of A/(RT).
For one mole at composition z, the scaled stability matrix is
`B_ij = sqrt(z_i) H_ij sqrt(z_j)`. The solver drives its smallest eigenvalue to
zero and drives the third derivative along `d_i = sqrt(z_i) u_i` to zero,
where u is the corresponding unit eigenvector. The largest absolute entry of
u is chosen positive to keep the third condition's sign consistent.

The fork uses the same two conditions in `critical::calcCriticalTV`, evaluating
the third derivative through symmetric perturbations of chemical potentials
at step 1e-4. Native classic mixing instead differentiates the extensive
Helmholtz energy in a single scalar hyperdual direction of order three.
This avoids a third-order composition tensor and four extra state preparations
for each condition evaluation. It includes temperature-dependent covolumes,
binary attraction/covolume interactions and linear-in-amount volume translation.
Wong-Sandler uses a five-point directional derivative of the analytic Hessian
with a positivity-limited step of at most 1e-3.

The local Newton iteration uses ln(T), ln(V), a centered numerical Jacobian
with step 2e-4, a maximum coordinate step of 0.3, and backtracking. Both residuals
must be below 2e-9. Invalid states and nonconvergence return errors. The active
symmetric matrix uses native Jacobi rotations; absent species are excluded
exactly instead of inserting the source's artificial 1e8 diagonal. The public
direction retains zero entries for those species.

`PreparedState::log_fugacity_tv_into(V, log_f, H)` exposes physical log fugacity
in Pa and the analytic mole Hessian with caller-owned output buffers, for a
one-mole state. It uses fixed-volume expressions directly, so neither positive
pressure nor division by the isothermal volume derivative is required. It
remains regular at a mechanical spinodal. Zero components have log fugacity
minus infinity and an infinite ideal diagonal; their cross derivatives remain
finite. At total amount N, use a normalized state at V/N and divide its Hessian
by N. The source TV interface includes component reference terms depending
only on temperature: its mole Hessian is identical, but its raw values are not
this physical log-fugacity convention.

## Qualification

`validate_mixture_critical.py` attempts 198 cases with the exact-source library:
PR, SRK, RK (the fork's SRK plus RK alpha), van der Waals, translated PR/SRK,
and Wong-Sandler with its default zero binary parameters. Five binary mixtures
span fractions 0.01–0.99; ternary and five-component mixtures include absent
components and an embedded pure fluid.

All **197 successful source returns** match in T, V and P at relative 1e-6 with
absolute floors 1e-6 K, 1e-12 m³/mol and 0.01 Pa. The observed largest relative
differences are 6.22e-8, 3.66e-7 and 3.46e-7, respectively. The source conditions
are independently checked: minimum eigenvalue below 2e-7 and third directional
derivative below 2e-6. There are no censored or unqualified returns. One default
van der Waals five-component solve fails in the source; the native local solve
succeeds, retained as a diagnostic without claiming parity for that call.

Fixtures retain each source active Hessian. Rust compares those at the
**source** T/V, independently of its own solved coordinates, at 1e-9 relative
and absolute, plus eigenvalues and third derivatives. Independent tests check
mole finite differences, physical TP/TV fugacity identity, negative-pressure
states, null-vector residuals, a temperature crossing of local stability,
finite differences of the third condition, component permutation, translation
invariance, explicit perturbed seeds, pure/zero-component limits and invalid
inputs. The eigensolver is also checked on known spectra through 20 variables.

## Limits

These conditions locate local mixture critical states. They do not establish
global multiphase stability, a positive relaxed fourth-order condition, critical
end points or a complete critical locus. A mixture critical point generally has
a negative nonzero fixed-composition pressure-volume derivative; using only a
pressure inflection would give a different result. Other activity mixing
Hessians and other equation-of-state families remain pending. Quantum critical
states and custom beta/interaction configurations have not been
reference-qualified by this grid. Original Newton
trajectories and failed-call semantics are not reproduced by this physical API.

Reproduce with the exact-source build and the tools below:

```text
python tools/rust_port/validate_mixture_critical.py
cargo test --test mixture_critical
python tools/rust_port/benchmark_mixture_critical.py --compiler <gfortran>
```

The benchmark compares complete default critical solves in separate native
executables for PR/SRK C1+C2 and PR C1+C2+C3, 25 compositions per profile,
200 solves per sample and five samples after warm-up. Model construction and
Python/FFI are outside timing. All T/V/P outputs are checked individually before
timing and accumulated in three separate checksums during timing. The report
pins source, Cargo and harness hashes and rejects edits during measurement.
The measured Windows release medians are:

| Workload | Original | Native | Speedup |
|---|---:|---:|---:|
| PR C1+C2 | 190.98 µs | 43.93 µs | 4.35× |
| SRK C1+C2 | 200.85 µs | 43.88 µs | 4.58× |
| PR C1+C2+C3 | 216.37 µs | 61.39 µs | 3.52× |

All 75 preflight states match within the documented tolerances; maximum
absolute T/V/P differences are 2.07e-5 K, 4.88e-11 m³/mol and 1.82 Pa.
See [the measured report](benchmark-mixture-critical.json) for platform,
samples and hashes. These workload measurements do not measure the entire library.
