# PeTS equilibrium at fixed volume and UV/SV flashes

`Pets::equilibrium_tv(T,V)` returns fluid equilibrium for positive temperature
and molar volume (m3/mol). The result contains a `FlashResult`, molar V/H/S/U
and isobaric Cp in `properties`, and the equilibrium heat capacity at constant
volume in `heat_capacity_v` (J/mol/K). This is distinct from the homogeneous
`caloric_tv(T,total_V,n)` kernel, which can evaluate unstable EOS loops.

Below Tc, the solver computes pure coexistence. If V lies between the saturated
volumes, the lever rule fixes vapor fraction:

```
betaV = (V - VL) / (VV - VL)
```

The equilibrium U and S include the changing phase amounts. Their analytic
fixed-volume derivative follows Maxwell and Clapeyron relations:

```
dVi/dT along saturation = (dPsat/dT - (dP/dT)Vi) / (dP/dV)Ti
Cv_equilibrium = sum_i beta_i * [Cv_i - T*((dP/dT)Vi - dPsat/dT)^2/(dP/dV)Ti]
(dS/dT)V = Cv_equilibrium/T
```

This expression avoids subtracting ideal caloric reference offsets. It includes
phase redistribution and remains meaningful in pure coexistence, where the
isobaric equilibrium heat capacity is undefined. At a saturated endpoint the
reported Cv is the derivative from the coexistence side. Independent centered
differences of equilibrium U and S verify the interior response.

Outside coexistence, a mechanically stable homogeneous state is required.
Subcritical states also verify the minimum-Gibbs pressure branch. Equilibrium
packing fraction is limited to 0.75, matching pressure-root work; the raw TV
kernel retains its larger domain. Pressure and Cv must be positive and finite.
The branch-discovery and low-temperature extrapolation limitations in
[PeTS equilibrium](pets-equilibrium.md) still apply. No solid phase is considered.

`uv_flash(U,V)` and `sv_flash(S,V)` recover temperature, pressure and phase
amounts from molar U (J/mol) or S (J/mol/K) and V. The pure problem is reduced
to a bracketed temperature solve along this equilibrium isochore, using Cv or
Cv/T for safeguarded Newton steps. It avoids solving a coupled mixture split
for a one-species model. Homogeneous outputs use `Single`; coexistence outputs
include endpoints with explicit quality. Full source flag/trajectory identity
is a separate pending requirement.

```rust
use thermopack_rs::Pets;
let model = Pets::new()?;
let state = model.equilibrium_tv(100.0, 1e-3)?;
let recovered = model.uv_flash(state.properties.internal_energy, 1e-3)?;
assert!((recovered.temperature - 100.0).abs() < 1e-6);
# Ok::<(), thermopack_rs::Error>(())
```

The `_with_options` variants use `VolumeFlashOptions`. Defaults are 80–999 K
and final pressure bounds 10 Pa–100 MPa. Temperature guesses are clamped to
the temperature bounds. A failed initial evaluation tries the critical
temperature, then the upper bound. Temperature trials can retreat toward a
valid endpoint. The scalar Newton iteration limit excludes initial bracket
expansion. Invalid specifications, bounds, guesses, non-finite properties or
nonpositive Cv return errors. Final pressure bounds are enforced; intermediate
positive-pressure states can lie outside those bounds while bracketing T.
The pressure guess is validated but is unnecessary for this pure temperature
solve. The shared composition-search `tp` options do not apply.

The owned ideal-gas Cp and H/S references apply to homogeneous and coexistence
properties and flashes. No original-library calls or runtime dependencies are
introduced. Native tests cover parameter scaling, U/S/V recovery, independent
fixed-volume response, caloric setter ownership, bounds and non-finite input.

## Source qualification

The matrix attempts 666 calls across three sigma/epsilon pairs. Homogeneous
seeds cover seven T/epsilon values (0.5–3) and seven pressures (0.001 Pa–50 MPa).
Coexistence seeds cover five T/Tc values (0.5–0.999), five qualities including
both endpoints and the same three UV/SV initialization modes. All 399
qualifying results match: 224 UV and 175 SV. The remaining 267 calls comprise
260 Fortran process terminations and seven returned UV states that fail the
qualification below. These remain recorded failures, never parity successes.

The exact-source validator calls the original public UV dispatcher both with
and without T/P guesses. This fork comments out its Python SV wrapper, so SV
uses the existing probe into `twoPhaseSVflash`, with TP-initialized phase arrays
and perturbed T/P guesses. Custom PeTS parameters explicitly refresh source
critical metadata before flash calls.

Source results must have successful status, finite quantities, negative
homogeneous dP/dV, pressure closure at the source volume tolerance of 1e8
machine epsilons, scaled U/(RT) or S/R closure within 1e-7 and relative V closure
within 1e-7. Failed calls and available diagnostics are retained. Comparisons
cover T/P, interior two-phase status/quality, V/U/S/H. Quality within 1e-8 of
an endpoint is compared as a homogeneous physical state; metadata retains raw
source flags, fractions, compositions and closure errors. This is the same
explicit endpoint convention as the [PeTS PH/PS fixtures](pets-flash.md).

No UV/SV speedup is claimed by this implementation stage. The historical TV
kernel benchmark measures a different operation. Complete inherited PeTS
interfaces, broad degeneracy and the remaining ThermoPack model families
remain pending in the [coverage ledger](coverage.md).

```
python tools/rust_port/validate_pets_isochoric.py --compiler-directory /path/to/compiler/bin
cargo test --test pets_isochoric
```
