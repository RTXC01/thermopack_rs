# Homogeneous PS and metastable isentropes

The cubic port solves local pressure/entropy equations and continues these
roots down a pressure grid. It also implements the fork's complete metastable
isentrope mapping convention, including equilibrium initialization and raw
spinodal endpoint reporting. These methods use native Rust throughout.

## Local solver

`Cubic::homogeneous_ps(P, S, z, T_guess, V_guess)` solves two equations at a
fixed composition. Pressure is Pa, entropy J/mol/K, temperature K and molar
volume m³/mol. Compositions are normalized by the prepared-state interface;
the model's ideal caloric references apply.

The solve uses logarithmic T/V variables, an analytic 2×2 Jacobian, bounded
logarithmic steps and backtracking. Maxwell's relations give S_T=Cv/T and
S_V=P_T. Residuals are scaled by max(|P|, 100000 Pa) and max(|S|, 1 J/mol/K),
with a 1e-10 convergence tolerance. Each step uses the physical Cv from the
native caloric evaluator, including the model's temperature-dependent terms.
No finite-difference derivative or equilibrium flash is needed for a local
solve. Both final equation residuals and the iteration count are returned.

This is a local root solve: the initial coordinates determine the branch.
Negative pressure and locally unstable states are supported. A converged
root does not establish global phase stability. The source's seed ±100 K
and volume factor-of-ten bounds, nonlinear trajectories and failure codes
are not reproduced.

## Physical continuation

`Cubic::homogeneous_isentrope(z, S, seed, options)` refines an explicit initial
homogeneous state and continues it toward `minimum_pressure`. The analytic
constant-entropy tangent predicts T/V; local PS correction closes both
equations. Pressure steps subdivide when correction fails or either logarithmic
coordinate changes by more than 0.4. Requested grid points are returned;
accepted intermediate steps remain internal.

By default the initial state must be locally stable. A bracketed sign change
of the minimum TV Helmholtz stability eigenvalue terminates the curve at a
refined spinodal with absolute eigenvalue below 1e-9. This includes mixture
composition instability, which can precede a mechanical dP/dV=0 boundary.
`stop_at_spinodal=false` disables this test and permits unstable continuation.
The supplied minimum pressure may be negative. The method reports an error
when correction or subdivision cannot converge; it does not return an
unqualified partial curve as success.

```rust
use thermopack_rs::{Cubic, CubicKind, Phase};
use thermopack_rs::metastable::IsentropeOptions;
# fn main() -> Result<(), thermopack_rs::Error> {
let model = Cubic::new(CubicKind::PengRobinson, &["C1"])?;
let state = model.prepare(180.0, &[1.0])?;
let v = state.volume(6e6, Phase::Liquid)?;
let s = state.caloric_tv(v)?.entropy;
let seed = model.homogeneous_ps(6e6, s, &[1.0], 180.0, v)?;
let curve = model.homogeneous_isentrope(&[1.0], s, seed,
    IsentropeOptions { points: 20, ..Default::default() })?;
# Ok(())
# }
```

Eigenvalue termination uses the current cubic stability Hessian scope:
classic PR/SRK/RK/VDW, translated PR/SRK and supported Wong-Sandler parameters.
Other activity-model Hessians remain pending, although their local PS equations
can be solved. Pure and mixture tests independently check pressure/entropy
closure and stability on both sides of the termination condition. Coarse
two-point pure traces test pressure subdivision. A local sign-change search
does not prove that an arbitrarily coarse path has no unbracketed stability
crossings or that the selected root is globally stable.

## Fork mapping convention

`Cubic::baseline_metastable_isentrope(z, P0, S, Pmin, count)` follows
`spinodal::map_meta_isentrope`. It first locates the last source-selected
spinodal entropy crossing. Mixture lookup now uses
`baseline_mixture_spinodal_intersection`: a complete source-order volume trace
to the model Tmin, the last strict property sign change and the source
property refiner. Only entropy selection is qualified indirectly here;
standalone mixture searches for the other five properties remain pending.

The initial state comes from the native equilibrium PS flash. If it contains
two phases, the larger phase fraction selects a volume root evaluated at the
**feed** composition. A single phase flag selects the liquid root as in the
source. This first state is reported without local entropy correction.
The endpoint retains the source's uncorrected spinodal trial volume. Interior
points use the source linear pressure grid and analytic tangent, followed by
the native PS corrector. The result carries actual EOS pressure and entropy
errors separately from reported pressure. Source first/last states may fail
closure even when the routine succeeds.

`count=1` returns the initial point with `PointLimit`; the original divides
by zero for this case. Positive Pmin<P0 and at least one point are required.
Full source error/last-iterate behavior, near-critical branch choices and
equilibrium flash failures are outside this qualification.

## Qualification against this fork

`tools/rust_port/validate_metastable.py` isolates source failures and imposes
an eight-second per-call deadline. Committed metadata records every attempt,
failure, endpoint warning, exact source/library identity and fixture hash.

| Comparison | Attempts | Qualified returns | Source failures |
|---|---:|---:|---:|
| Local homogeneous PS | 370 | 298 states | 72 nonzero VDW statuses |
| Complete metastable maps | 140 | 97 curves / 1164 states | 34 aborts, 9 errors |

There were no censored calls or unqualified finite returns. All qualified
coordinates match. Local worst relative T/V differences are 2.55e-11 and
3.79e-11; complete-map T/V/P differences are 8.37e-9, 2.09e-8 and 2.43e-8.
The comparison uses relative 2e-6 with absolute floors 2e-5 K, 1e-10 m³/mol and
2 Pa. The local grid covers classic PR/SRK/RK, translated PR/SRK, UNIFAC,
HV1 and quantum cubic hydrogen/neon, including negative pressures. The
classical VDW attempts all return nonzero source statuses and do not qualify
source-result parity for that solver. Complete maps cover successful classic
PR/SRK/RK/VDW and translated PR/SRK pure, binary and ternary cases.

The 25 first/last coordinate warnings remain explicit. Tests compare all raw
coordinates and separately require native interior equation closure. Independent
tests cover negative-pressure/unstable local roots, zero feed species, analytic
Jacobian identities by finite differences, pure and mixture physical endpoints,
invalid inputs and point counts. This work also exposed and corrected the
previous Rust unrecorded HV-pair attraction/default-interaction defects;
see [baseline conventions](baseline-defects.md).

## Performance

The native-executable benchmark checks all 170 preflight T/V/P/S states and
all four returned-quantity sums during every timed sample. Target preparation
and model construction are excluded from both executables. Curve timing
includes spinodal search, equilibrium initialization, continuation and entropy
evaluation for every returned point. The local solve timing includes returned
pressure and entropy evaluation on both sides. Results and workload details
are stored in `metastable-benchmark.json`; speedups apply only to its workloads.

Five samples after one warm-up, Windows x86-64, Rust 1.99 release with thin
LTO and one codegen unit versus gfortran 14.3 at O3:

| Workload | Source median | Rust median | Ratio |
|---|---:|---:|---:|
| Local PS, PR methane/ethane | 29.36 µs | 5.80 µs | 5.06× |
| Local PS, SRK carbon dioxide/propane | 31.62 µs | 6.80 µs | 4.65× |
| Complete 12-point PR methane map | 5.24 ms | 0.56 ms | 9.34× |
| Complete 12-point PR methane/ethane map | 9.41 ms | 2.77 ms | 3.40× |

The benchmark uses 25 local specifications per mixture and five initial
temperatures per mapping workload. Each sample runs 10000 local solves or
20 complete maps. Native and reference sources, benchmark inputs and Cargo
settings are fingerprinted; the runner rejects changes during measurement.
Numerical equivalence of these outputs does not imply identical nonlinear
iteration trajectories or performance for other models or inputs.
