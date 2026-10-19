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
2. The remaining polynomial roots are isolated through their derivative hierarchy.
3. Earlier derivative roots are recovered in reverse order by bracketed bisections
   between known stationary points.
4. Evaluating `exp(-x^2) * g(x)` avoids overflow while strictly preserving signs,
   zeros, and ordering.

This resolves narrow spinodal pairs and tangencies without a fixed density grid.
The intervals where `dP/drho >= 0` are extracted as mechanically stable branches,
and bracketed Newton iterations find the precise density roots for any given pressure.
The hierarchy remains subject to floating-point rounding, particularly at repeated
roots. Absolute Horner bounds reject search ranges that overflow its polynomial
arithmetic. An upper density bound is required explicitly; extrapolated equations
can have extra stable loops, so a globally physical branch is not inferred from
mechanical stability alone. Prepared branches can be reused across pressures.

## Critical and triple point properties

* `source_critical_point()` provides the original tabulated approximate critical
  properties (Tc, Pc, Vc, Zc) from `src/mbwrdata.f90` (`STORED_PROPERTIES`).
* `source_triple_point()` provides original CoolProp triple-point (T, P) metadata.
* `critical_point()` and `critical_point_from(T_init, rho_init)` solve the true
  critical point equations `dP/drho = 0` and `d^2P/drho^2 = 0` via a damped Newton
  iteration in `(log T, log rho)` using analytic pressure derivatives through third
  density order and second density / first temperature order (`critical_jet`).
  Positive pressure and stable curvature `d^3P/drho^3 > 0` are required.
  Default seeds use the source MBWR table where present, otherwise the compiled
  component Tc and Pc/(Zc*R*Tc). Gamma is a fitted exponential coefficient and
  does not generally estimate critical density; using it previously prevented
  R21 convergence. Sixteen legacy records have neither original source table nor
  compiled component seeds; these require explicit initial values. The solve is
  local, and multiple critical points of extrapolated/custom equations are possible.
  It does not reproduce the original coarse `MBWR_criticalParameters` grid search.

## Testing and verification

`validate_mbwr_equilibrium.py` calls original Fortran routines through the
unmodified-source reference library. Its committed CSV fixtures and metadata cover:

* All 65 stored records' critical/triple table presence and values.
* 49 local critical points from independent Newton solves using the original
  `MBWR_pressure` derivatives. The reference density Jacobian column is differenced;
  the native solver uses analytic third-order derivatives. All comparisons pass.
  The 16 records lacking default seeds are recorded separately, not counted as
  successful or failed numerical calls.
* 3456 original `MBWR_density` calls, spanning 48 records with compiled component
  acentric factors, six reduced temperatures (0.65 to 2), six reduced pressures
  (0.01 to 3), and both input phases. Initialization passes source critical tables
  where available and compiled component guesses otherwise; it bypasses the
  unrelated initializer's coarse critical grid search. All calls return finite.
* 3406 of those density outputs close pressure, have nonnegative P_rho, and lie
  inside the explicit native domain. Every one belongs to the native stable root
  set. CSV also preserves requested/reported source phases, returned pressure and
  its density derivative. The 50 unqualified iterates remain in metadata with raw
  requests/results. The other 17 records lack a compiled source acentric factor
  and are listed as unavailable for this initialized-density qualification.

Root-set membership is the qualified comparison. Sixty finite source results
differ from selecting the lowest/highest stable density at the requested phase;
55 also report a source fallback phase change; five retain the requested flag
while returning the opposite endpoint root. All are retained with complete
native root sets. All 50 unqualified returns fail pressure closure; 21 also have
negative pressure slope (these counts overlap). They are not phase-selection
parity successes. Source iteration paths, fallback phase conventions and the
coarse critical search remain pending. Full ideal/caloric, saturation, flash and
CSP integration are also pending; no MBWR speedup is claimed here.

Independent native tests verify critical equations, derivative jets for all 65
fits, near-critical spinodal pairs at T/Tc=1-1e-7 for all 49 critical fixtures,
pressure closure/stability, minimum Gibbs selection, custom ideal coefficients,
search-domain endpoints and overflow rejection. Quantity-specific tolerances are
recorded in fixture metadata; they do not imply bitwise cross-platform identity.

```
python tools/rust_port/validate_mbwr_equilibrium.py --compiler-directory /path/to/compiler/bin
cargo test --test mbwr_equilibrium
```
