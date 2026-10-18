# Fixed-volume caloric properties and equilibrium

`PreparedState::caloric_tv(v)` returns molar pressure, internal energy, entropy
and Cv at T and physical molar volume v. Units are Pa, J/mol and J/(mol K).
The analytic Helmholtz expressions include temperature-dependent covolume beta
and volume translation. This path also evaluates homogeneous negative-pressure
states with EOS volume above the covolume; these are useful initialization
states and do not assert stable equilibrium. Cv differentiates actual caloric
values, including the fork's ideal-H2 values; source-reported defective heat
capacities are described in [baseline-defects.md](baseline-defects.md).

`Cubic::uv_flash(z, u, v)` and `sv_flash(z, s, v)` recover an equilibrium
`FlashResult` from normalized mole fractions and molar U or S and V. Options
are available through `volume_flash_options()` and the `_with_options` methods.
They inherit model temperature bounds, use 10 Pa–100 MPa pressure bounds by
default, and accept optional positive T/P guesses. Search bounds limit every
coupled step. Failed convergence returns an error.

```rust
use thermopack_rs::{Cubic, CubicKind};
let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2", "C3"])?;
let z = [0.5, 0.3, 0.2];
let feed = model.tp_flash(200.0, 1e6, &z)?;
let properties = model.equilibrium_properties(&feed)?;
let recovered = model.uv_flash(&z, properties.internal_energy, properties.volume)?;
assert!((recovered.temperature - 200.0).abs() < 1e-6);
# Ok::<(), thermopack_rs::Error>(())
```

The coupled equations follow `uv_solver`, `sv_solver` and `state_functions`.
Unknowns are component vapor amounts, log temperature and log pressure.
Residuals are phase chemical-potential differences, `(Utarget+P*Vtarget-H)/(RT)`
or `(Starget-S)/R`, and `P*(V-Vtarget)/(RT)`. The analytic symmetric matrix is
solved with the native Bunch–Kaufman implementation. Steps preserve material
balance, retain trace amounts by updating the smaller phase amount, and use
bounded residual-reducing backtracking. Pure latent states retain a quality
unknown; a TP flash cannot determine their quality.

Initialization tries a homogeneous TV caloric solve, ordinary TP states and
finite phase-amount seeds. Fixed-composition saturation estimates are computed
only when an earlier seed has failed. Final U/S and V closure is independently
checked. For mixtures, a strict TP flash verifies the equilibrium branch and
Gibbs energy. Its stopping tolerance is separate from the coupled residual's
default 1e-12 tolerance: the tighter UV/SV state is retained after verification.
The implementation uses physical caloric derivatives, including the corrected
vdW/H2 derivatives, so it does not reproduce every original iteration trajectory.

Qualification contains 409 finite source flash results from 546 attempts:
230 UV and 179 SV. It covers PR/SRK/RK/vdW, default and tcPR construction,
pure methane/CO2/H2, a hydrocarbon ternary, absent feed components, compressed
and dilute homogeneous states, and pure qualities 0.1/0.5/0.9. UV exercises
the Python wrapper's zero-initialized arrays with default and explicit T/P
guesses. The fork comments out its Python SV wrapper, so a probe calls the
public Fortran SV dispatcher with TP-initialized phase arrays and perturbed
T/P. Native SV receives only T/P guesses and computes its own phase seeds.

Single-phase comparisons canonicalize source output to x=z, y=0, betaV=0,
betaL=1. The source has inconsistent single-phase array/flag conventions;
native homogeneous results use TP phase labels and fraction sentinels, and
raw source single-phase output compatibility remains pending. Two-phase arrays,
fractions, T and P are compared directly. Temperature tolerance is
2e-5 K + 2e-7 relative, pressure 0.5 Pa + 2e-7 relative, fractions/compositions
2e-7 absolute + 2e-7 relative; the two-phase flag must agree exactly.
Independent native tests require much tighter caloric/volume closure and
check component balance and chemical equilibrium.

The 137 unsuccessful source calls (131 terminations and six invalid/nonzero
status returns) remain in fixture metadata. None establish parity. Maximum
scaled caloric and relative volume closure errors among retained source
outputs are 2.35e-9 and 1.50e-9. Another 141 direct TV source outputs qualify
P/U/S, including cryogenic QuantumCubic and custom quintic ci, quantum beta
and interactions. Independent finite differences check Cv, S_T=Cv/T and
S_V=P_T with temperature-dependent beta and ci. Broad cryogenic UV/SV,
critical degeneracy, complete source retry/output fidelity and other model
families remain pending. No UV/SV speedup has yet been measured.

Reproduce after building the exact-source library and probe:

```
python tools/rust_port/validate_isochoric.py --compiler-directory /path/to/compiler/bin
cargo test --test isochoric
```
