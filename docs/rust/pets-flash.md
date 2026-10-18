# Pure PeTS TP, PH and PS flashes

`Pets::tp_flash(T,P)` selects the minimum-Gibbs homogeneous pressure root.
At pure coexistence, T and P do not determine vapor quality. Use
`saturated_flash(T,betaV)` to specify the molar vapor fraction in [0,1].
The latter returns `TwoPhase` even at the two saturated endpoints.

`ph_flash(P,H)` and `ps_flash(P,S)` recover temperature and phase amounts
from molar enthalpy (J/mol) or entropy (J/mol/K). Below the critical pressure,
they first determine saturation temperature and the liquid/vapor caloric
limits. Specifications between those limits determine quality directly.
Outside the latent interval, a safeguarded scalar Newton solve uses the
homogeneous heat capacity and minimum-Gibbs pressure root. The scalar
temperature solver is shared with the existing cubic implementation.

```rust
use thermopack_rs::Pets;
let model = Pets::new()?;
let flash = model.saturated_flash(100.0, 0.4)?;
let state = model.equilibrium_properties(&flash)?;
let recovered = model.ph_flash(flash.pressure, state.enthalpy)?;
assert!((recovered.vapor_fraction - 0.4).abs() < 1e-8);
# Ok::<(), thermopack_rs::Error>(())
```

The `_with_options` variants accept `CaloricFlashOptions`, including initial
temperature, bounds, scaled residual tolerance and iteration limit. Defaults
are 80–999 K. The composition-search `tp` member of those shared options does
not apply to this pure-fluid solver. Invalid specifications, bounds, guesses
and off-range targets return errors. These methods make no Fortran calls.

`equilibrium_properties(&flash)` returns molar V/H/S/U and equilibrium Cp.
Pure coexistence has undefined isobaric heat capacity, represented by `None`.
Homogeneous Cp is available only when finite and positive. Phase volumes are
recomputed for the current model; public cached volume fields are not trusted.
Two-phase inputs require pure normalized compositions, finite physical amounts
that sum to one, and pressure within 1e-8 relative of this model's coexistence
curve. Properties use the owned ideal-gas reference and parameter state.

TP uses `Single` with the existing -1 fraction sentinels when only one stable
root is available. With several roots, the selected smallest/largest root is
liquid/vapor. PH/PS label homogeneous subcritical states by their side of the
latent interval; supercritical results use `Single`. Both pure composition
arrays are populated for PH/PS, following the source pure caloric interfaces.
Full reproduction of source phase labels and solver trajectories remains pending.

## Exact-source qualification

`validate_pets_flash.py` attempts 1035 public source calculations using three
sigma/epsilon parameter pairs. Homogeneous inputs cover seven reduced
temperatures (T/epsilon=0.5–3) and seven pressures (0.001 Pa–100 MPa).
Latent inputs cover five T/Tc values (0.5–0.999), five qualities including both
endpoints, and PH/PS with and without an initial temperature. Custom source
parameters explicitly refresh source critical metadata before flash work.

All 923 qualifying results match: 126 TP, 395 PH and 402 PS. The other 112
calls comprise seven PH calculation failures and 105 returned states failing
finite, mechanically stable pressure closure or caloric closure. Metadata
retains all failed requests and available diagnostics. A source pressure root
must meet its own relative criterion, 1e8 machine epsilons, and negative dP/dV.
Caloric closure must be within 1e-6 in H/(RT) or S/R. Those qualification
criteria are separate from the final quantity-specific parity tolerances.

The comparison includes T/P, interior two-phase status/quality, occupied
volume(s), H and S. Quality within 1e-8 of zero or one is compared as a
homogeneous saturated endpoint: the unused phase volume and phase amount do
not define a different physical state. This resolves 76 cases where the
source reports `TwoPhase` exactly at an endpoint and Rust can return a
homogeneous phase due to root/inverse-saturation rounding. The original
source flag, both fractions, both composition arrays, closure diagnostics
and warnings are preserved for every finite result. Raw flag identity is
not claimed by these fixtures.

Native tests independently recover temperatures from H/S for 108 homogeneous
states, close latent quality and U=H-PV, and reject malformed/off-curve states
and invalid options. Caloric setter tests also exercise owned nondefault Cp
and H/S references. Fixtures use 2e-7 relative tolerance with absolute
tolerances documented per quantity. No new flash performance claim is made.
The earlier [TV kernel benchmark](pets.md) does not measure these solvers.
PeTS UV/SV, complete inherited interfaces and broad critical/low-temperature
degeneracy remain pending; the full library migration is still in progress.

```
python tools/rust_port/validate_pets_flash.py --compiler-directory /path/to/compiler/bin
cargo test --test pets_flash
```
