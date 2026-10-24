# Pure multiparameter fixed-volume equilibrium

`Multiparameter::equilibrium_tv(T,V)` evaluates the pure liquid/vapor branch
state, mixed properties and equilibrium Cv. UV and SV invert that fixed-volume
path using analytic Cv (or Cv/T), local bracket expansion, Newton steps and
bisection. No fixed temperature grid is scanned. Homogeneous states require
positive pressure and mechanical/thermal stability; coexistence uses the
qualified native saturation pair and lever-rule volume. This does not prove
global stability of every extrapolated EOS branch.

The two-phase heat capacity includes quality redistribution at fixed volume:

```text
Cv_eq,i = Cv_i - T (P_T,i - dP_sat/dT)^2 / P_V,i
Cv_eq = (1-beta) Cv_eq,L + beta Cv_eq,V
```

Endpoint states retain the coexistence-side response when marked `TwoPhase`.
That differs from the homogeneous-side response at the same limiting point.

## Bounds and guesses

`volume_flash_options()` supplies the source package solution bounds of
80..999 K and 10 Pa..100 MPa for physical models. The LJ fit maxT/maxP fields
do not replace those package bounds. Reduced LJTS uses 0.3..3 Tc and
1e-8..70 in native reduced units. Inputs otherwise use SI and the existing
DEFAULT caloric reference convention.

The component-database triple temperature repairs an out-of-range UV initial
guess; it is **not** a lower limit on the final UV temperature. In the original
`twoPhaseUVsingleComp`, the outer seed adjustment reads component Ttr, but
`singleCompUV_Tv` and `singleCompUV` use package temperature bounds for their
solutions. This matters for krypton targets below database Ttr. SV uses the
package bounds. Caller temperature bounds constrain the complete native solve.

Pressure bounds constrain the accepted solution. Intermediate positive,
mechanically stable fixed-volume trials can exceed Pmax while descending from
a critical-temperature seed to a valid cold dense root. The source homogeneous
UV solver likewise checks pressure after solving its fixed-volume temperature
path. Invalid final pressure is rejected. A pressure guess is validated but
unused in native iteration, since P is evaluated directly from T/V.

## Current qualification

The reviewed audit has **560 source attempts** across fourteen variants:
56 homogeneous controls and 504 latent targets at T/Tc=0.5, 0.9 and 0.999,
qualities 0, 0.5 and 1, UV/SV and supplied/absent guesses. Targets are computed
with original TV kernels at previously qualified saturation volumes, including
independently qualified saturation recoveries. They are retained before any
solver call, so a source termination cannot erase its target.

All **560 native outputs** independently close against original TV kernels
with no native failures or physical coordinate mismatches. Checks include
specified V and U/S, pressure, mechanical and thermal stability, mixed caloric
properties, equilibrium Cv and coexistence chemical potential. Native caloric
closure uses a scaled 3e-10 tolerance in R T or R; relative volume closure is
2e-8. Comparisons to independently closed source coordinates use explicit
per-quantity tolerances recorded in the validator.

Of the source attempts, **240 return qualified states**: **181 strict
coordinate/phase matches** and **59 endpoint phase-flag differences**.
These differences remain compatibility gaps: selecting the homogeneous or
coexistence flag can change downstream sound and response. They are not
classified as equivalent representations. Broad endpoint trajectory fidelity
remains pending.

The remaining **320 source calls terminate**: 305 saturation nonconvergences,
12 SV nonconvergences, two UV nonconvergences and one density iteration limit.
All 240 returned source states also pass independent original TV-kernel
pressure, thermal stability and coexistence chemical-potential checks.
The 320 finite native recoveries are independently checked using the same
original kernels. Full requests, raw source output, status where exposed,
warnings, native results and fingerprints are retained in metadata. The public
UV wrapper exposes no status argument; its status is recorded as null, rather
than inventing a zero status. Reduced wrapper inputs/outputs scale physical
neon LJTS; direct original reduced TV kernels check native outputs.

The audit uses four-second per-call deadlines. A timeout is censored, killed
and reaped and never silently retried. Source reuse verifies binary/worker/input
identities, recomputes all targets with original TV kernels, and requires the
actual call schedule to match. Historical collection and current audit generator
fingerprints are recorded separately.

```text
python tools/rust_port/validate_meos_isochoric.py --compiler-directory /path/to/compiler/bin
python tools/rust_port/validate_meos_isochoric.py --reuse-source --allow-known-gaps --compiler-directory /path/to/compiler/bin
cargo test --test multiparameter_isochoric
```

Default validation fails on this partial source surface. The explicit
`--allow-known-gaps` mode accepts only exact registered case keys and unchanged
outcome classes; new errors or enlarged endpoint differences fail. Improved
source outcomes may remove a gap. Native errors and physical mismatches always
fail. No censored calls occur in the committed grid.

Standalone Rust replay covers all 560 cases without Fortran and verifies native
raw fractions, bounds, V and caloric specification closure. Endpoint coordinates
are replayed, while exact endpoint flag comparison belongs to the live source
audit in its recorded arithmetic environment: another floating-point build may
select a different side. The 181 full matches and 59 flag differences above
are observations of that qualification build, not a promise of identical raw
endpoint flags across builds. Non-endpoint phase and quality values are checked
in fixture replay. Original raw phase flags remain in the committed audit.

Regression tests cover all fourteen variants, independent energy/entropy Cv
differences, pressure-limit initialization, sub-triple UV roots, endpoints and
bounds. Python negative tests reject incomplete outputs, invalid raw fractions,
nonclosing source states, coordinate drift and undeclared outcomes. CI runs
standalone replay and the explicit partial-source gate. Global stability,
complete endpoint trajectory/flag fidelity and broad degeneracy remain pending.
No performance measurement or speedup claim is made for these methods.
