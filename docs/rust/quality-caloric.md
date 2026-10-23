# Caloric fixed-quality constraints

The cubic fixed-quality corrector accepts total molar enthalpy or entropy,
alongside T, P and log K. It solves for temperature, pressure and both phase
compositions at a fixed feed and molar vapor fraction. This is a local
coexistence calculation; it does not establish global phase stability or a
unique solution when the caloric specification intersects several branches.

## API

```rust
use thermopack_rs::{Cubic, CubicKind};
use thermopack_rs::quality::QualitySpecification;
# fn main() -> Result<(), thermopack_rs::Error> {
let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"])?;
let initial = model.quality_temperature(1e6, &[0.5, 0.5], 0.25)?;
let target = model.quality_properties(&initial)?.enthalpy + 500.0;
let next = model.refine_quality_specification(
    &initial, QualitySpecification::Enthalpy(target)
)?;
assert_eq!(next.feed, initial.feed);
assert!((model.quality_properties(&next)?.enthalpy-target).abs() < 1e-4);
# Ok(())
# }
```

`quality_state_with_specification` takes a `QualityGuess` with T, P and y/x
ratios. `refine_quality_specification` reuses a previous state's original feed
and phase amounts. Existing `quality_state_with_guess` and
`refine_quality_state` signatures remain available for T/P/log-K constraints.
`baseline_quality_state_with_specification` preserves the fork's numerical
conventions for qualified source comparisons.

H and U use J/mol, S uses J/(mol K), and V uses m³/mol. `quality_properties`
reports lever-rule H/S/U/V at the supplied state's phase volumes and pressure,
including the model's ideal-gas reference levels. It does not perform another
equilibrium solve. Tests shift both reference levels to make the specified H
or S zero, check exact absent species, and verify H=U+PV.

For a pure feed the physical method brackets temperature locally around the
warm guess and bisects along the existing pure coexistence calculation, below
the EOS critical point. Embedded pure feeds and phase fractions 0, 0.3 and 1
are tested. The source-compatible primitive continues to reject pure
log-K equations. Pure caloric critical endpoints and arbitrary cold H/S
initialization remain outside this qualification.

## Analytic derivatives and source convention

The final equation is `(property-target)/max(abs(target),1)`. The physical
Jacobian differentiates each normalized phase composition away from the
material-balance root. With raw phase sums Sx/Sy, q_i=Xraw_i*Yraw_i/z_i,
phase molar properties h_l/h_v and partial molar properties hbar_l/hbar_v,
its log-K derivative is

```text
beta*(1-beta)*q_i * ((hbar_v_i-h_v)/Sy - (hbar_l_i-h_l)/Sx) / scale.
```

The same expression applies to entropy. T/P derivatives use analytic Cp,
volume and fugacity derivatives. Phase preparation is reused from the
fugacity calculation. The physical Cp follows Cv-T*P_T²/P_V.

The fork's `saturation.f90` Jacobian uses the **sum** of the vapor and liquid
partial molar caloric properties in its log-K row, rather than the derivative
above. The baseline Rust path preserves that sum and the fork's existing VDW
and type-10 hydrogen Cp derivative conventions. The physical path uses the
normalized derivative. This distinction affects convergence and sometimes
branch selection; treating these paths as interchangeable would hide a real
source defect.

The fixture `quality_caloric_jacobian.csv` retains 280 initial source residuals,
analytic caloric rows and independently differenced source-function rows,
including guesses whose source solve subsequently fails. Rust's baseline row
matches the source analytic row; the physical row matches finite differences
of the source function. Additional off-root finite-difference tests cover
classic, translated, UNIFAC, Huron-Vidal, zero-species and quantum cases,
including quality endpoints.

## Qualification and explicit differences

`validate_quality_caloric.py` makes 292 unique source attempts. There are 203
finite closing H/S returns, 12 earlier fixed-P quality initialization failures,
and 77 caloric solver failures. There are no censored calls, unqualified finite
returns or returned pressure warnings. The exact-source library identity and
all failures are retained in fixture metadata.

The baseline matches all 203 finite states and rejects all 77 reachable source
caloric failures. This establishes returned-state and error-category behavior
for this matrix, not identical failed vectors, iteration counts or raw status
integers. The physical solver converges for all 280 reachable guesses. Each
of the 77 additional roots is independently checked through the source's
fugacity, pressure and H/S calculations before being counted as recovered.

The matrix covers PR, SRK, RK and VDW; binary, ternary and zero feeds; beta=0/1
and interior values; translated cubics, UNIFAC, HV1/HV2/NRTL; quantum H2/Ne and
H2/D2; and selected H targets shifted by ±500 J/mol or S by ±1 J/(mol K).

Of the 203 finite source states, 200 physical results also match all coordinate
comparisons. Three differences remain explicit, without widening the ordinary
tolerance:

* RK methane/propane, beta=0.75, S: vapor volume differs by about 2.7e-6
  relative and dP/dT by 2.1e-6 relative.
* PR UNIFAC methane/ethane/propane, beta=0.75, S: vapor volume differs by
  about 2.06e-6 relative. These first two cases reflect the source's looser
  caloric closure and are bounded individually in regression checks.
* SRK HV1 methane/water, beta=0.25, H: source-compatible correction reaches
  approximately 142.112 K and 1.688 MPa, while the physical correction stays
  near the warm branch at 135.957 K and 1 MPa. Both close the same H, phase
  equilibrium and material balance. They have different entropy. All three
  differing physical states pass independent source property/phase checks.

Ordinary source comparisons use relative 2e-7 for the baseline and 2e-6 for
the physical path, with absolute floors 2e-5 K, 2 Pa, 1e-10 m³/mol, 2e-8 for
mole fractions, 1e-5 J/mol for H and 1e-8 J/(mol K) for S. Slopes use relative
2e-6 plus 1e-3 Pa/K. Physical target closure is checked separately at
1.1e-9*max(abs(target),1); the independent source checks use normalized
caloric residual 2e-8, phase/material residual 2e-7 and pressure residual
1e-3 Pa + 1e-7*P. The validator saves raw differences and only permits the
three named, bounded observations. A passing regression check does not count
them as strict coordinate parity.

Complete H/S curve tracing, general coexistence constraints, broad critical
regularization, nondefault source global controls, global stability and other
EOS families remain pending. The [T/P/log-K point](quality.md) and
[phase-fraction curve](quality-envelopes.md) qualifications remain separate.

## Performance

The [native-executable report](quality-caloric-benchmark.json) compares source,
source-compatible Rust and physical Rust from identical perturbed coexistence
seeds. Four binary profiles include every finite interior-quality H/S fixture
state for PR/SRK methane/ethane, translated PR methane/propane and PR UNIFAC
methane/propane. All twelve returned quantities are checked before timing,
including both phase compositions and volumes, dP/dT, total H and total S.

Each sample cycles 1000 corrections and compares all twelve separate output
sums. Five measured samples follow one warm-up. Model construction, input
parsing, seed generation and Python/FFI are excluded. The physical path's
corrected Jacobian changes the work needed to converge; this is part of the
measured algorithm. It is not a whole-library or cold-initialization claim.

Local Windows release medians per correction, including H/S reporting:

| Workload | States | Source | Compatible Rust | Physical Rust | Source/compatible | Source/physical |
|---|---:|---:|---:|---:|---:|---:|
| PR methane/ethane | 13 | 1349.07 µs | 165.41 µs | 31.30 µs | 8.16x | 43.11x |
| SRK methane/ethane | 6 | 1197.66 µs | 143.10 µs | 25.78 µs | 8.37x | 46.46x |
| Translated PR methane/propane | 4 | 986.85 µs | 128.22 µs | 30.84 µs | 7.70x | 32.00x |
| UNIFAC PR methane/propane | 4 | 3637.07 µs | 617.34 µs | 89.07 µs | 5.89x | 40.84x |

The 27 preflight states are the same for all three executables. The report
retains every input, sample, separate checksum, source-library identity and
native/harness/Cargo fingerprint; fingerprints are checked again after timing.
