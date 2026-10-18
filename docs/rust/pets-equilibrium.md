# PeTS pressure roots and pure-fluid coexistence

`Pets::volume(T,P,Phase)` solves for molar volume (m3/mol) at positive pressure.
`prepare_equilibrium(T)` returns reusable isotherm branches; its `volume` and
`volumes` methods reuse spinodals across pressure queries. Liquid selects the
smallest stable volume, vapor the largest, and `MinimumGibbs` compares pure
Gibbs energies across the available stable roots. These are homogeneous
pure-fluid root choices; no mixture stability calculation is implied.

```rust
use thermopack_rs::{Pets, Phase};
let model = Pets::new()?;
let isotherm = model.prepare_equilibrium(120.0)?;
let liquid = isotherm.volume(1e6, Phase::Liquid)?;
let vapor = isotherm.volume(1e6, Phase::Vapor)?;
let saturation = model.saturation_pressure(120.0)?;
let inverse = model.saturation_temperature(saturation.pressure)?;
# Ok::<(), thermopack_rs::Error>(())
```

The admissible packing interval for pressure-root work is [0,0.75], following
the source SAFT solver's upper bound. The raw TV kernel remains available up
to packing fraction one. Spinodal discovery combines logarithmic dilute and
linear dense sampling, includes the critical packing fraction to resolve
near-critical loops, and refines derivative sign changes. Root refinement uses
bracketed Newton steps on stable intervals and a pressure roundoff bound that
includes cancellation between hard-sphere and dispersion derivatives. This
is a numerical search over the represented branches, not a proof of complete
root isolation for every unphysical low-temperature extrapolation.

Some low-temperature isotherms have extra dense spinodals from the mathematical
EOS. The saturation bracket uses the first vapor/liquid spinodal pair. Root
queries preserve all stable intervals detected within the admissible packing
range. No solid-phase competition is considered.

`critical_point()` returns true pure-fluid T/P/V/Z. A fixed-size fourth-order
Taylor expansion in packing fraction supplies analytic density derivatives
for the two critical equations. The temperature column of the Newton system
is differenced; density columns are analytic. A thread-safe cache stores the
parameter-independent reduced critical point, never an active model. Sigma
and epsilon/k rescale it for each model, so setters cannot leave stale critical
data. The default critical point is approximately 148.110649153 K,
4.789851142 MPa and 7.791769622e-5 m3/mol. Its density curvature is independently
checked with the existing Helmholtz kernel.

`saturation_pressure(T)` brackets pure coexistence and equates Gibbs energies
using the two phase roots. The chemical-potential difference is evaluated
without subtracting ideal caloric reference constants. Clapeyron slopes use
residual thermal derivatives and imposed P*deltaV, also cancelling ideal
reference offsets analytically. `saturation_temperature(P)` inverts this curve
with bracketed Newton steps. At/above the critical point these APIs return an
error; the exact representable critical T/P pair has an explicit critical-root
limit. Independent tests verify root T/P response, chemical equilibrium,
Clapeyron slopes, inverse recovery, parameter scaling and reference invariance.

## Source iteration compatibility

`baseline_volume_iteration(T,P,Phase)` ports
`pc_saft_nonassoc_volume_solver` with its liquid/vapor initialization, packing
floor 1e-10, ceiling 0.75, Newton/bisection order and stopping conditions.
`PetsVolumeIteration` exposes the finite last volume, iteration count, separate
volume/pressure convergence flags and relative pressure error. `MinimumGibbs`
is not an initialization phase for this compatibility method.

The original routine returns its last iterate even when pressure has not
converged. At 68 K and 0.001 Pa, its vapor seed can finish on a liquid-like
volume instead of the dilute vapor root. The compatibility method preserves
that result and reports pressure nonconvergence. The strict solver returns
the large vapor volume and independently closes pressure. These two behaviors
are intentionally available through separate methods.

Qualification attempts 627 source calls across three sigma/epsilon pairs,
seven temperatures, seven pressures from 0.001 Pa to 100 MPa and both source
phase seeds. All 294 finite raw source iterates match, including nonconverged
returns. A separate root-response comparison requires independently recomputed
pressure closure at the original relative tolerance, `1e8*epsilon_f64`, and
negative dP/dV. All 237 qualifying roots and their T/P derivatives match.
The other 57 source returns remain in metadata and in the raw-iterate fixtures;
they are not silently promoted to converged solutions.

All three source critical points, 18 saturation-pressure states and 18 inverse
states match. Saturation temperatures span 0.35–0.999 of Tc. The reference uses
the actual Fortran parameter setter and explicitly calls its public critical
refresh before custom-parameter saturation: that setter alone leaves source
fallback critical metadata stale. The reference critical calculation uses
explicit T/V starting guesses; its numerical V can be less accurate than the
native analytic-density critical equations.

Fixture metadata records exact-source/probe hashes and quantity-specific
absolute tolerances with 2e-7 relative tolerance. The source and native EOS
kernel remain independent implementations. The original Fortran files are
unmodified. Broad critical degeneracy, complete source trajectory/phase-flag
fidelity, flash solvers and the remaining inherited SAFT interfaces remain
pending. No pressure-root or saturation speedup is claimed yet; the separate
[kernel benchmark](pets.md) measures the TV derivative path only.

```
python tools/rust_port/validate_pets_equilibrium.py --compiler-directory /path/to/compiler/bin
cargo test --test pets_equilibrium
```
