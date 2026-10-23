# Pure multiparameter coexistence

`Multiparameter::saturation_pressure(T)` returns temperature, pressure, liquid
and vapor molar volumes, and the Clapeyron slope. The same query is available
on `PreparedMultiparameter`. `saturation_temperature(P)` inverts that curve
from a cold `0.8 Tc` guess; `saturation_temperature_with_guess(P,T)` accepts an
explicit initial temperature. All calculations are native Rust.

The supported models are the five NIST fluids C3, R134A, N-H2, O-H2 and P-H2;
full and truncated-shifted Lennard-Jones for Ne, Ar, Kr and Xe; and reduced
LJTS. Physical models use K, Pa and m3/mol. Reduced LJTS retains its original
reduced units and R=1, as described in [the kernel qualification](multiparameter.md).

```rust
use thermopack_rs::Multiparameter;
let propane = Multiparameter::new("C3")?;
let point = propane.saturation_pressure(280.0)?;
let inverse = propane.saturation_temperature(point.pressure)?;
assert!((inverse.temperature - 280.0).abs() < 1e-5);
# Ok::<(), thermopack_rs::Error>(())
```

## Equations and safeguards

The unknowns are `u=[ln(delta_l),ln(delta_v)]`, where `delta=rho/rhoc`.
For residual reduced Helmholtz energy `a`, write `a_d=delta*a_delta`,
`a_dd=delta^2*a_deltadelta` and `a_t=tau*a_tau`, with `tau=Tc/T`.
At each density define `p*=delta*(1+a_d)` and `D=1+2*a_d+a_dd`.
The pressure scale is `R*T*rhoc`. The two residuals and their analytic
Jacobian are

```text
F = [p*_l - p*_v,
     u_l - u_v + a_l - a_v + a_d,l - a_d,v]
J = [delta_l*D_l, -delta_v*D_v;
     D_l,         -D_v]
```

The second equation is chemical-potential equality divided by RT. The ideal
temperature and reference constants cancel exactly at fixed T. The prepared
temperature cache is reused at every density, and the two-density evaluation,
Newton step, backtracking and continuation use fixed-size storage without
heap allocation. No repeated TP phase-volume search is needed inside this
coupled solve.

Ancillary densities seed the solution. Stable trial branches must satisfy
`delta_l>1>delta_v>0`, `delta_l<=16` and `D_l,D_v>0`. Newton log-density steps
are capped at 0.5, with at most 20 backtracks and 100 iterations. The pressure
residual tolerance is `4e-11*abs(p*_v)` plus a 32-epsilon cancellation floor
computed from both phase pressure terms. The chemical residual divided by
the log-density gap must be below `2e-11`. Solving the raw chemical equation,
while checking its normalized residual at convergence, avoids an unhelpful
quotient-Jacobian direction toward the hydrogen vapor spinodal.

Accepted solutions have positive pressure, a log-density gap of at least
`1e-6`, and a sampled point with `D<0` between the stable endpoints. The
nonconvexity check tries the reducing density, arithmetic density midpoint
and a 32-interval density grid, stopping at the first negative stiffness.
It prevents reporting an unresolved coincident pair in a monotone region.
It is a conservative sampled check, not an exhaustive root or global phase
stability analysis. At low temperatures the reducing density itself can have
positive stiffness, so checking only that point incorrectly rejects valid
propane coexistence.

If the initial solve fails, a lower point at `0.9*T` seeds continuation toward
the target. Each step halves the remaining distance to stored Tc, with bounded
temperature backtracking if needed. This recovers the ortho-hydrogen point
at `T/Tc=0.9999` without loosening the convergence tolerances.

Pressure uses the vapor value to avoid dense-phase cancellation. The latent
enthalpy is `R*T*(a_t,v-a_t,l)+P*(Vv-Vl)`; dividing by `T*(Vv-Vl)` gives
`dP/dT`.

Temperature inversion first solves coexistence at the initial guess. It then
corrects log temperature and both log densities together, adding the equation
`ln(p*_v*R*T*rhoc/Ptarget)=0`. The extra Jacobian column for log temperature is
`[-delta_l*a_dt,l+delta_v*a_dt,v,
-(a_t,l-a_t,v+a_dt,l-a_dt,v)]`, where `a_dt=tau*delta*a_taudelta`.
The pressure-specification row is
`[0,delta_v*D_v/p*_v,1-delta_v*a_dt,v/p*_v]`. A stack-allocated 3x3 system
uses the existing equilibrated linear solver. Density steps retain the 0.5
cap, log-temperature steps are capped at 0.15, and every trial requires
stable branches and a temperature below stored Tc. Final acceptance keeps
the same coexistence tolerances and requires `abs(ln(Psat/Ptarget))<2e-10`.
The final point passes the same phase-gap and nonconvexity checks as a
pressure query.

The original nested Clapeyron temperature iteration remains a bounded
fallback if coupled correction fails. This fallback recomputes coexistence
at each trial temperature and backtracks failures toward the last valid
point. The faster primary path avoids fully resolving both densities for
every intermediate temperature. Both analytic Jacobians are independently
differenced off root.

## Exact-fork evidence

The validator makes 232 public source saturation attempts over all 14 model
variants at `T/Tc=0.45,0.5,0.6,0.75,0.9,0.95,0.99,0.999,0.9999`, followed by
inverse requests wherever the source pressure call returns finite outputs.
All attempts remain in the metadata:

| Source outcome | Count | Native qualification |
|---|---:|---|
| Finite, distinct, stable and chemically closed phases | 164 | All five returned fields match |
| Reported saturation error | 24 | 20 pressure and 4 temperature errors; Rust solutions independently checked |
| Unqualified returned state | 44 | All are inverse calls with coincident phase volumes and nonfinite slope; Rust solutions independently checked |
| Censored deadline | 0 | Separate classification retained |

The 164 direct comparisons use relative tolerance `2e-6` and absolute floors
`[2e-7 K, 1e-2 Pa, 1e-12 m3/mol, 1e-12 m3/mol, 1e-4 Pa/K]` for T, P, Vl,
Vv and slope. Reduced LJTS uses the corresponding reduced units. Maximum
relative differences are `2.10e-11,6.02e-10,2.50e-7,3.51e-7,2.95e-7`.
All 164 source states also pass independent pressure checks.

All 232 native results, including the 68 recoveries, are evaluated again by
the unmodified source Helmholtz, pressure and caloric kernels at the native
phase volumes. The source confirms stable endpoints, a nonconvex point between
them, equal chemical potential, pressure closure and the Clapeyron slope.
The largest source chemical residual is `3.07e-12` in RT units. The 68 recovery
coordinates are **not** claims that the original saturation solver returned
the same coordinates. They have a separate fixture and explicit source outcome.

The public reduced LJTS saturation wrapper is unavailable in the fork. Its
reference coordinates are scaled from the identical physical LJTS Ne fit,
using source epsilon/k, sigma, gas constant and Avogadro's constant:
`Tscale=epsilon/k`, `Vscale=N_A*sigma^3`, `Pscale=R*Tscale/Vscale`.
Independent reduced-state checks use the original reduced `init_LJTS` object
directly. No new source EOS or saturation algorithm is substituted.

Tests also check all 14 models at all nine temperatures, cold inversion of
each native pressure, pressure/fugacity/latent-entropy closure, and independent
curve differences for the Clapeyron slope at 0.8 and 0.99 Tc. The coupled
Jacobian is independently differenced off root. Invalid values and the
unresolved critical endpoint return errors.

## Remaining scope and reproduction

The branch restriction uses the model's stored reducing density and
temperature. These methods do not refine a physical critical point or
enumerate all pressure roots. Inversion rejects pressures at or above the
EOS pressure at stored Tc/rhoc. Very narrow or unsupported coexistence pairs
return errors. Below-triple points are fluid extrapolations: solid stability
is outside this API. Full global phase stability, alternate reference
initialization, mixtures and multiparameter TP/PH/PS/UV/SV flashes and phase
diagrams remain pending.

The [benchmark record](meos-saturation-benchmark.json) compares complete
five-output pressure and inverse-temperature queries in native loops, with
all 76 selected source states compared before timing. Its five profiles are
NIST propane, NIST R134A, LJ krypton, LJ argon and LJTS krypton. Source-failed
and coincident-phase states are explicitly outside those timing workloads.
Hydrogen and reduced LJTS performance are not measured by this benchmark.
See [the performance summary](performance.md) for the measured medians.

```text
python tools/rust_port/validate_meos_saturation.py --compiler-directory /path/to/compiler/bin
cargo test --test meos_saturation
python tools/rust_port/benchmark_meos_saturation.py --compiler /path/to/gfortran --preflight-only
python tools/rust_port/benchmark_meos_saturation.py --compiler /path/to/gfortran
```
