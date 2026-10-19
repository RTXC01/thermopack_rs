# MBWR bounded pressure branches and critical points

`Mbwr::prepare_equilibrium(T, maximum_density)` isolates mechanically stable
pressure branches over `[0, maximum_density]` in mol/m3. Its `volume(P, Phase)`
solves for molar volume (m3/mol) at positive pressure. `Phase::Vapor` selects
the lowest stable density root (largest volume), `Phase::Liquid` selects the
highest stable density root (smallest volume), and `Phase::MinimumGibbs`
compares pure residual Gibbs energies across all mechanically stable roots in
the domain.

```rust
use thermopack_rs::{Mbwr, MbwrKind, Phase};
let model = Mbwr::new(MbwrKind::Mbwr32, "C1")?;
let critical = model.critical_point()?;
assert!((critical.temperature - 190.56).abs() < 0.1);
let branches = model.prepare_equilibrium(150.0, 35000.0)?;
let vapor_volume = branches.volume(1e6, Phase::Vapor)?;
let liquid_volume = branches.volume(1e6, Phase::Liquid)?;
assert!((model.pressure_tv(150.0, vapor_volume, 1.0)?.value - 1e6).abs() < 0.01);
assert!((model.pressure_tv(150.0, liquid_volume, 1.0)?.value - 1e6).abs() < 0.01);
# Ok::<(), thermopack_rs::Error>(())
```

## Spinodal discovery via derivative hierarchy

Because MBWR equations of state feature polynomial and Gaussian-exponential terms
`P(T, r) = P_poly(T, r) + exp(-gamma * r^2) * P_exp(T, r)`, finding spinodals
(`dP/drho = 0`) cannot rely on fixed grid sampling without risking missing narrow
spinodal pairs or near-critical extrema.

The implementation applies Rolle's theorem recursively to the Gaussian-polynomial
system:
1. Differentiating `exp(x^2) * A(x) + B(x)` until `B(x)` vanishes.
2. The remaining polynomial roots are found analytically/hierarchically.
3. Earlier derivative roots are recovered in reverse order by bracketed bisections
   between known stationary points.
4. Evaluating `exp(-x^2) * g(x)` avoids overflow while strictly preserving signs,
   zeros, and ordering.

This resolves narrow spinodal pairs and tangencies without a fixed density grid.
The intervals where `dP/drho >= 0` are extracted as mechanically stable branches,
and bracketed Newton iterations find the precise density roots for any given pressure.

## Critical and triple point properties

* `source_critical_point()` provides the original tabulated approximate critical
  properties (Tc, Pc, Vc, Zc) from `src/mbwrdata.f90` (`STORED_PROPERTIES`).
* `source_triple_point()` provides original CoolProp triple-point (T, P) metadata.
* `critical_point()` and `critical_point_from(T_init, rho_init)` solve the true
  critical point equations `dP/drho = 0` and `d^2P/drho^2 = 0` via a damped Newton
  iteration in `(log T, log rho)` using analytic pressure derivatives through third
  density order and second density / first temperature order (`critical_jet`).
  Positive pressure and stable curvature `d^3P/drho^3 > 0` are required.

## Testing and verification

Unit and integration tests in `tests/mbwr_equilibrium.rs` verify:
* Exact retrieval of stored critical and triple points for all Bender19 and MBWR32 fits.
* Convergence of the analytic critical point equations, verifying zero first and second
  density derivatives of pressure to within `1e-6 * RT`.
* Branch isolation, spinodal density ordering, and exact pressure closure for both
  vapor and liquid roots.
* Minimum Gibbs energy selection and supercritical monotonic isotherm behavior.
* Comprehensive rejection of non-physical/non-positive inputs.
