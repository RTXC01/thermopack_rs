# MBWR-19 and MBWR-32 kernels

`Mbwr` ports the pure-fluid Bender (19 coefficient) and MBWR-32 equations from
this fork. `Mbwr::new(MbwrKind,identifier)` selects the first matching source
record, including compiled fluid aliases. `from_record("mbwrcxN")` selects an
explicit fit; `mbwr_records()` exposes all 65 records in unchanged source order.
There are 54 MBWR19 and 11 MBWR32 records. Three MBWR32 methane fits are stored;
the ordinary methane constructor selects `mbwrcx1`, as the original lookup does.
Legacy identifiers remain available without requiring global ideal-gas data.

`with_parameters(kind,coefficients,gamma)` owns a custom 19/32 coefficient
array; values must be finite and gamma positive. Raw parameters use the source
layout/units, including gamma in (L/mol)^2. `parameters()` returns those values
without conversion. `record()` identifies database-backed models; custom
parameters have no database record. Models and clones own independent arrays.
The generated [data manifest](mbwr-data.json) records original input hashes.
The generator verifies both declaration and lookup-array order.

```rust
use thermopack_rs::{Mbwr, MbwrKind};
let model = Mbwr::new(MbwrKind::Mbwr32, "methane")?;
let prepared = model.prepare(300.0)?;
let pressure = prepared.pressure_tv(1e-3, 1.0)?;
let residual = prepared.residual_helmholtz(1e-3, 1.0)?;
let original = prepared.baseline_residual_helmholtz(1e-3, 1.0)?;
# Ok::<(), thermopack_rs::Error>(())
```

All public state inputs use SI: T in K, density in mol/m3, total volume in m3,
and mole number in mol. The polynomial backend uses r=rho/1000 in mol/L.
`pressure_density` returns P (Pa) and first/second T/rho derivatives, including
the mixed derivative. `pressure_tv` returns P and T/total-V/n derivatives.
`residual_density` returns dimensionless alpha=Ares/(nRT) and all first/second
T/rho derivatives. `residual_helmholtz` returns extensive F=Ares/(RT), in mol,
with all nine first/second T,V,n derivatives. It uses the shared pure-fluid
extensive chain rule. These are homogeneous mathematical EOS evaluations;
negative pressure and unstable loops remain observable. No phase stability
or fitted-temperature validity is asserted by these kernel calls.

## Stable and source-order formulations

Pressure follows the source polynomial plus exp(-gamma*r^2) terms. Temperature
preparation groups pressure and integrated-alpha coefficients and their two
thermal derivatives. Simultaneous Horner evaluation returns polynomial values
and density derivatives. No scratch allocation or Fortran calls are needed.

The stable alpha path integrates each exponential pressure term analytically.
For its m-th term it uses:

```
I_m(r) = integral(0..r) s^(2m-1) exp(-gamma*s^2) ds
alpha_exp = sum_m [E_m(T)/(1000*R*T)] * I_m(r)
```

A positive incomplete-gamma series avoids cancellation for gamma*r^2<4;
above that threshold a recurrence uses `expm1` for the first integral.
Density derivatives come directly from the integrand. Temperature derivatives
use prepared coefficients. At zero density this gives exact analytic limits:
alpha, alpha_T and alpha_TT vanish, pressure vanishes and dP/drho=RT.
`virial_coefficients(T)` returns analytic B=alpha_rho(0), C=alpha_rhorho(0).
The original density derivative routine divides by density at zero; the native
zero-density limit is an explicit extension.

The original implementation expands that integral into 22/40 terms, including
large density-independent polynomial/exponential pairs. Their subtraction can
lose significant digits at extreme dilute/cryogenic extrapolations. The native
`baseline_residual_density`, `baseline_residual_helmholtz` and
`baseline_pressure_tv` methods retain the original integrated coefficients,
term order and derivative expressions. `baseline_helmholtz_terms()` exposes
the read-only expansion for inspection. The generator translates the original
coefficient expressions and preserves their arithmetic tree. Gamma powers
five/six use the original compiler's multiplication chains; a final coefficient
ulp can otherwise be amplified by cancellation. All stored integrated
coefficients compare exactly in the pinned reference qualification.

Source-order thermal coefficients are initialized lazily inside a prepared
state and reused across densities. The default stable pressure/integral path
does not pay for their power evaluations. The cache belongs to that immutable
temperature/model state and contains no active global EOS. Source-order zero
density also uses the analytic native limit. Non-finite results produce errors.

## Qualification and remaining work

The reference probe reads the original 65 records and calls the unmodified
`computeZCoeff`, `computeHelmCoeff`, `alphar_derivatives`, `MBWR_Fres`,
`MBWR_press`, `makeParam` and `MBWR_pressure` routines. Initialization loads only
coefficients/schema/counts needed by these kernels. It skips unrelated critical
iteration and ideal metadata: the tested routines do not use them, and several
legacy records have no global fluid/critical entry. This is kernel qualification,
not qualification of the complete original constructor or its density solver.

All 65 raw parameter rows and all 65 integrated-coefficient rows match exactly.
All 6760 core states match 26 quantities through the source-order methods:
six density-alpha derivatives, ten extensive-F quantities, four TV pressure
quantities and six direct pressure-density quantities. The matrix covers all
fits, twelve temperatures from 0.5 to 2000 K, eight densities from 1e-5 to
20000 mol/m3, and extra n=0.3/3.5 extensive-scaling cases. All 6890 source calls
are finite; no failed calls were removed. Quantity-specific absolute tolerances
and 2e-8 relative tolerance are recorded in fixture metadata.

The default stable integral is also compared over the full matrix. Its 318
source-summation differences are retained in the
[comparison report](../../tests/fixtures/mbwr_stable_differences.json),
without widening tolerances. They distinguish the two numerical formulations;
they are not parity successes for the stable path. Use the `baseline_*` methods
when migrating code that relies on the original rounding behavior. Native
finite-difference tests independently verify all density/thermal derivatives,
Helmholtz pressure identities, extensive Euler/scaling, the integral switch,
zero-density limits, perfect-gas coefficients, ownership and invalid inputs.

This stage does not provide MBWR ideal/full caloric state functions, density
iteration, critical metadata/search, saturation, flashes or CSP reference-fluid
integration. Those remain pending within the full ThermoPack migration. No
MBWR speedup is claimed before a dedicated native-loop benchmark.

```
python tools/rust_port/generate_mbwr_data.py
cargo fmt --all
python tools/rust_port/validate_mbwr.py --compiler-directory /path/to/compiler/bin
cargo test --test mbwr
```
