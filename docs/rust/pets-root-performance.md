# PeTS pressure-root native-loop measurements

The initial [report](benchmark-pets-roots-before.json) measures the production
root implementation from `18cca1b` with the new benchmark harness. Each row
performs 20000 volume queries in native code. Five measured samples follow a
discarded warm-up. Both models use default PeTS parameters. Model construction
is outside the timer; per-call parameter preparation is included in fresh and
source-iteration Rust rows. Isotherm branch discovery is excluded only from
the explicitly prepared row.

| Workload | Original Fortran ns/query | Rust fresh strict | Rust prepared strict | Rust source iteration |
|---|---:|---:|---:|---:|
| Liquid seed, 120 K | 40201.55 | 50387.41 | 1671.22 | 1142.43 |
| Vapor seed, 120 K | 19812.03 | 50385.62 | 1683.19 | 761.43 |
| Vapor seed, 300 K | 16758.18 | 41901.65 | 563.33 | 679.77 |

Each isotherm sweeps 100 pressures from 0.5 to 0.995 MPa. The original Fortran
routine is `pc_saft_nonassoc_volume_solver`. Rust's source-iteration method
reproduces that seeded last-iterate algorithm. The strict paths discover
stable branches and enforce pressure closure; prepared queries reuse that
discovery. These paths have different validation/search responsibilities, so
the table does not imply identical internal work. Every compared state does
have a closing stable source root and the same selected volume.

Before timing, all 300 source states and 900 Rust path/state combinations
compare V, P, dP/dT and dP/dV. Each independently closes pressure at the source
criterion of 1e8 machine epsilons and has negative dP/dV. Every sample also
checks its volume checksum against Fortran. Metadata records maximum per-state
differences, every timing sample, compiler/platform details, original-source
metadata and Rust/harness/Cargo hashes. Mutating those inputs during a run
invalidates the measurement. No fast-math or reduced precision is used.

The initial fresh strict path is slower than the original on this workload.
Its logarithmic/dense density grid, sorting and spinodal refinement are repeated at
every temperature preparation. Reusing an isotherm already removes most of
that cost. The measurement identifies fresh branch discovery as an optimization
target; it is not a claim of a completed library-wide performance goal.

## Density-only root optimization

The [optimized report](benchmark-pets-roots.json) uses the same harness, states,
iteration counts and tolerances. The original sampling points and sign-change
refinements are retained, including low-temperature extra branches. A process-
wide immutable cache stores only the parameter-independent sorted density
grid. Each prepared isotherm still discovers its own temperature-dependent
spinodals. Stable interval endpoint pressures and roundoff scales are cached
inside that isotherm, avoiding repeated boundary evaluations in each query.

Root discovery/refinement now evaluates only the two density derivatives
needed for pressure and its density slope. It avoids constructing and checking
unused thermal/mixed derivatives at every density point. The full public
Helmholtz/TV path retains all derivative and finite-value checks. Pressure,
slope and the cancellation-aware scale are checked in the root kernel.
Temperature preparation also caches the common dispersion multiplier.
An independent comparison against the full kernel covers 528 combinations
of parameters, temperature and density, including zero density, cryogenic
extrapolation and the dense boundary.

| Workload | Fortran ns/query in optimized run | Rust fresh strict | Rust prepared strict | Rust source iteration |
|---|---:|---:|---:|---:|
| Liquid seed, 120 K | 39437.20 | 26954.05 | 1078.20 | 1083.91 |
| Vapor seed, 120 K | 26531.11 | 26884.31 | 1060.40 | 734.89 |
| Vapor seed, 300 K | 21811.79 | 20433.80 | 308.05 | 638.45 |

Fresh Rust queries improve 1.87x, 1.87x and 2.05x over the initial Rust report;
prepared queries improve 1.55x, 1.59x and 1.83x. The fresh vapor row remains
roughly tied with Fortran. Fortran timing also varies between the two runs,
so use the raw samples and within-run comparisons when assessing ratios.
The grid search remains an important cost for changing temperatures; this
stage does not establish uniformly faster cold solves or complete root
isolation beyond the documented numerical sampling.

The existing [Helmholtz benchmark](pets.md) measures the derivative kernel,
not this root algorithm. The numbers here do not qualify saturation, flashes,
different temperatures/pressures, broad degeneracy or other EOS families.

```
python tools/rust_port/benchmark_pets_roots.py --compiler /path/to/gfortran --output docs/rust/benchmark-pets-roots-before.json
python tools/rust_port/benchmark_pets_roots.py --compiler /path/to/gfortran --preflight-only
```
