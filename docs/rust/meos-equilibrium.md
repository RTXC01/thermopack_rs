# Multiparameter pressure branches and local critical points

All fourteen native NIST/Lennard-Jones variants provide analytic density
derivatives, local ordinary critical points and reusable bounded pressure
branches. These build on the existing [state kernels](multiparameter.md) and
[pure coexistence methods](meos-saturation.md). The original phase-switch
density iteration and its strict `PreparedMultiparameter::volume` wrapper
retain their previous behavior and qualification.

```rust
use thermopack_rs::{Multiparameter, Phase};
let model = Multiparameter::new("C3")?;
let critical = model.critical_point()?;
let branches = model.prepare_equilibrium(280.0, 8.0 * model.rc())?;
let mut densities = vec![0.0; branches.maximum_root_count()];
let count = branches.densities_into(1e6, &mut densities)?;
let stable_volume = branches.volume(1e6, Phase::MinimumGibbs)?;
assert!(count > 0 && critical.temperature > 280.0 && stable_volume > 0.0);
# Ok::<(), thermopack_rs::Error>(())
```

Physical models use SI molar units; reduced LJTS keeps its original reduced
units and R=1. The maximum density is an explicit search boundary. It does
not assert the validity of an extrapolated EOS, or enumerate roots beyond it.

## Analytic pressure derivatives

`PreparedMultiparameter::pressure_density_jet(rho)` returns P, P_rho,
P_rhorho, P_rhorhorho, P_T, P_Trho and P_Trhorho, including rho=0. It uses
fixed-size stack storage and no numerical differences or heap allocation.
The established pressure/first-derivative summation is retained; the extra
derivatives use an exponential-polynomial recurrence.

For a residual term `c(T)*delta^d*exp(q(delta,T))`, define

```text
H0 = d*delta^(d+1) + q_delta*delta^(d+2)
H(k+1) = derivative(Hk, delta) + q_delta*Hk
p* = delta + sum(c*exp(q)*H0),   p* = P/(R*T*rhoc)
```

H1, H2 and H3 provide the first three density derivatives of p*. The
exponent is zero, `-delta^l`, or
`-eta*(delta-epsilon)^2-beta*(tau-gamma)^2`. At fixed density, each pressure
term's temperature derivative has weight `1-t` for polynomial/damped terms
or `1-t+2*beta*tau*(tau-gamma)` for Gaussian terms. Units follow by restoring
R, T and rhoc. Polynomial capacities are checked against all compiled terms.

`critical_point_from(T,rho)` solves P_rho=P_rhorho=0 with an analytic
log-temperature/log-density Newton Jacobian, bounded steps and backtracking.
Acceptance requires a normalized residual below `2e-11`, positive pressure
and positive third density derivative. It is a local ordinary fluid critical
candidate, without a global uniqueness or phase-stability claim.

`critical_point()` first uses the stored reducing T and density. If that
guess fails or reaches a nonordinary stationary point, it tries 0.95 and
1.05 times the reducing density and selects the higher-temperature ordinary
result. The full LJ fit needs this policy: independent source-pressure solves
find two nearby ordinary candidates and a nonordinary central candidate.
For LJ argon, the ordinary temperatures are about 154.1628000 and
154.1628405 K, with volumes about 7.65412e-5 and 7.52532e-5 m3/mol. Both local
results are retained in the fixture. Other seeds can select another basin.
This qualification concerns the source EOS equations; the original generic
critical solver's trajectories, statuses and default candidate selection are
separate pending compatibility work.

## Branch discovery and queries

`prepare_equilibrium(T,maximum_density)` and prepared `equilibrium(maximum_density)`
discover pressure-stiffness zeros over the whole specified density interval.
They allocate derivative polynomials, interval bounds, a work list and the
resulting branch data once. Subsequent scalar `volume` queries and
`densities_into` calls allocate no memory. `densities` is the convenience
method that allocates its result vector.

Polynomial coefficients and Horner evaluation use outward interval arithmetic.
The pressure/stiffness bounds are intersected with centered second-order
Taylor bounds; curvature uses a centered first-order bound. These avoid the
large dependency inflation of an uncentered polynomial bound near a flat
critical region. A sign-definite stiffness bound excludes stationary points;
a sign-definite curvature bound establishes monotonic stiffness, allowing
single spinodals to be isolated by interval Newton contraction.

The search uses no fixed density grid. Work is bounded at 100000 cells and
depth 60. Regions that reach the density-resolution limit or have unresolved
stiffness at their center remain explicit boxes. `spinodal_brackets()` exposes
each density/pressure enclosure and whether one simple spinodal is isolated.
Queries that intersect an enclosure's pressure window return an error. The
algorithm does not invent a critical spinodal count or silently drop an
unresolved branch. Synthetic tests cover pairs only about 2e-6 apart and an
unresolved nearly double contact.

The bounds use an explicit four-ulp outward guard for platform `exp`.
Completeness inside the domain is conditional on that transcendental rounding
assumption and the documented arithmetic; this is not a formal proof about
every platform's math library. Source-backed tests independently check the
actual thermodynamic results.

Each stable interval uses bracketed Newton pressure inversion with at most
160 iterations. Ordinary pressure tolerance is `4e-12*P` plus a 64-epsilon
ideal-pressure scale; the density step or bracket must also converge. At a
floating-point stall, neighboring densities are checked and an individual-term
cancellation allowance may be used, capped by `2e-8*P+256*epsilon*R*T*rho`.
`pressure_roundoff(rho)` exposes the uncapped allowance for diagnostics. A
large allowance alone never makes a badly unclosed pressure acceptable.

Roots are sorted by density. Vapor/liquid select the lowest/highest density;
minimum Gibbs compares `alpha_res+ln(delta)+P/(rho*R*T)` across all stable
roots. The common ideal temperature/reference terms cancel. These explicit
selection rules can differ from the original iterator's phase switching.
Full solid stability, mixtures and global phase diagrams remain separate.

## Exact-fork qualification

The source probe adds only batched calls to the original `mp_pressure`.
The source library and its EOS implementations are unchanged. Reference
derivatives use five-point differences at two step sizes, with stored
per-field allowances for step disagreement and rounding amplification.
Independent critical searches use those original pressure derivatives.

| Evidence | Result |
|---|---|
| Pressure jets across all fourteen variants | 672 states, all seven quantities match |
| Three source-EOS critical guesses per variant | 38 ordinary results; four full-LJ central guesses are explicitly nonordinary |
| Ordinary local and default critical comparisons | 52 matches, with coordinate conditioning recorded |
| Independent bounded pressure search | 112 isotherms; doubled 2048/4096 grids, dilute log points and independently solved critical neighborhoods agree |
| Source spinodals | All 212 have native enclosures; isolated boxes match individually |
| Complete pressure-root queries | 804 source-backed sets match; 1110 native roots independently checked |
| Existing finite source density iterates in resolved sets | All 1653 are members of the native sets |
| Source phase-switch versus explicit extreme selection | 67 differences retained in metadata |

Every returned root is checked again by the original pressure/stiffness and
Helmholtz kernels. Source Gibbs scores verify native minimum-Gibbs selection;
source pressure and stiffness verify closure and mechanical stability.
Four R134A extrapolated roots use the bounded cancellation allowance and
remain explicit in metadata. The independent source pressure check is capped
at `4e-8*P+512*epsilon*R*T*rho`.

Eight further R134A complete-set requests at 0.3 or 0.6 reducing temperature
are density-precision errors. A source-bisection candidate in an interior
branch fails meaningful pressure closure: at one 4059.28 Pa request, its
pressure changes by roughly 4.9e14 Pa per representable density step. Accepting
that candidate merely because its density bracket collapsed would give a
misleading result. Raw candidates, actual source pressures, neighboring
pressures and density quanta remain in the unresolved fixture/metadata.

Those eight requests also contain ordinary outer roots. Explicit vapor or
liquid queries can still return a closing extreme branch when it is
resolvable; these results are source-checked and tested separately. The
complete root-set and minimum-Gibbs queries retain their errors. Sixteen
existing finite phase iterates belong to these requests, so they are kept
outside the 1653 complete-set membership count. The original phase iterator
and its 1694-attempt regression fixture remain available and unchanged.

The fixtures separate jets, ordinary critical equations, spinodals, complete
root sets and unresolved requests. All raw outcomes stay visible. This does
not qualify all densities, alternate reference initializations, mixture
behavior, generic critical-solver compatibility or complete multiparameter
flash/phase-diagram behavior.

## Measurement and reproduction

The [native benchmark](meos-branches-benchmark.json) compares fresh discovery
with reused branches for complete root sets plus minimum-Gibbs selection.
Five separate outputs are compared against the source-backed fixture before
timing and accumulated separately in every sample. Preparation and query
costs are reported separately. The original phase iterator has different
selection semantics, so this record makes no source-speedup claim.

| Profile | Discovery ms | Fresh complete query ms | Reused complete query us |
|---|---:|---:|---:|
| NIST C3, 0.95 Tc | 2.769 | 2.827 | 12.50 |
| NIST R134A, 1.05 Tc | 0.736 | 0.739 | 13.51 |
| NIST N-H2, 0.95 Tc | 0.841 | 0.864 | 7.65 |
| LJ argon, 0.95 Tc | 5.457 | 5.487 | 19.20 |
| LJTS krypton, 0.95 Tc | 2.585 | 2.597 | 14.10 |

Each profile preflights all eight matching source-backed fixture pressures.
A complete query requests all roots into caller storage and separately asks
for minimum-Gibbs selection; both API calls are timed. Fresh samples cycle
100 such queries, reused samples 1000, and discovery samples 10 preparations.
Five measured samples follow one warm-up. Model construction and parsing are
excluded. Reuse removes the dominant discovery cost; cold discovery remains
an optimization target. These timings cover this host and workload.

```text
python tools/rust_port/build_probe.py --compiler /path/to/gfortran
python tools/rust_port/validate_meos_equilibrium.py --compiler-directory /path/to/compiler/bin
cargo test --test meos_equilibrium
python tools/rust_port/benchmark_meos_branches.py --preflight-only
python tools/rust_port/benchmark_meos_branches.py
```

Source-cache reuse checks the library, probe, worker, entire generator and
legacy input-fixture hashes. Hosted CI results have not been observed here.
