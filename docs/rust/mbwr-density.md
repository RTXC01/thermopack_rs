# MBWR source density iteration

`Mbwr::baseline_density(T,P,phase,initialization)` ports the default
`MBWR_density` routine from `src/mbwr.f90`, with `meta_extrem` absent. It
returns the source's raw finite last iterate and phase flag, including unclosed
or unstable returns. This compatibility API complements the separately
[qualified bounded stable roots](mbwr-equilibrium.md).

```rust
use thermopack_rs::{Mbwr, MbwrKind, Phase};
let model = Mbwr::new(MbwrKind::Mbwr32, "C1")?;
let initialization = model.source_density_parameters()?;
let prepared = model.prepare_baseline_density(150.0, initialization)?;
let result = prepared.density(1e6, Phase::Vapor)?;
assert!(result.pressure_converged && result.mechanically_stable);
# Ok::<(), thermopack_rs::Error>(())
```

All public inputs and returned densities use SI: K, Pa and mol/m3.
`MbwrDensityParameters::new(Tc,Pc,rhoc,acentric_factor)` supplies the explicit
initialization data. Parameters are validated again during preparation even
if callers construct the public struct directly. `source_density_parameters()`
uses the original explicit MBWR critical table and compiled component acentric
factor. It returns an error when either is unavailable. It does not silently
substitute the independently refined physical critical point.

The source initializer's coarse critical grid search remains separate and
pending. For records without a critical table entry, callers can pass its
result, compiled component estimates, or explicit custom initialization data.
Different initialization data can change source fallback behavior. The
qualification below states which data it passes; it does not qualify the full
original initializer for every legacy record.

## Source arithmetic and control flow

The port retains:

* `MBWR_coef(1)` temperature expressions and `MBWR_pressure` summation order,
  including MBWR19 subtraction/unit conversion, Gaussian terms and original
  integer density powers. The compatibility pressure loop is separate from
  the physical kernel's Horner evaluation.
* `LiquidDensitySRKGuess`, its discriminant branches, special critical region,
  ideal estimate above T/Tc=1.05, and density scaling constants.
* `newton_density`'s phase-existence checks, 15-iteration cap and stale pressure
  deviation stopping test; initial curvature switches and retry in the other
  phase; and the first fallback's `barenewton` 20-iteration condition.
* The GUNDER=1 `MBWR_density_tplib` fallback, its quadratic step, reductions,
  second-start search, loose pressure check, and raw final phase convention.
  Its first pressure evaluation uses the liquid initial density even for a
  vapor-start density. The port retains this original arithmetic.

`MbwrDensityIteration` exposes density, reported phase, returned pressure and
P_rho, relative pressure error, derived `pressure_converged` and
`mechanically_stable` flags, plus a native `MbwrDensityPath` diagnostic.
`pressure_converged` tests the returned pressure's relative error against
1e-6; the source Newton loop does not necessarily enforce that on its final
iterate. Mechanical stability requires positive density and nonnegative P_rho;
neither flag establishes thermodynamic global stability. Non-finite final
outputs return an error. `Phase::MinimumGibbs` is not a source iteration seed;
use the bounded branch API for that selection.

`PreparedMbwrDensity` owns all fixed-temperature coefficients and SRK
initialization constants. It survives dropping the original model. Its density
queries allocate no scratch memory, and multiple prepared states share no
mutable global model. `source_parameters()` exposes read-only gamma/BP/BE
coefficients in the original backend's mol/L units. No speedup is claimed
until a dedicated native-loop comparison is measured.

## Reference qualification

`validate_mbwr_density.py` invokes the original `MBWR_density` and
`MBWR_pressure` through the exact-source probe with explicit initialization.
All 65 stored fits participate in 8840 attempted calls, spanning liquid and
vapor input flags, T/Tc=0.65 to 2, P/Pc=0.01 to 3, very low pressures down to
1e-4 Pa, the P=100 Pa initial-density cutoff and both sides of the SRK
critical/ideal-seed thresholds.

For 48 records, initialization uses the MBWR critical table where present,
otherwise compiled component Tc/Pc/Zc estimates, with the compiled acentric
factor. Seventeen legacy fits lack that compiled factor and explicitly use
0.1; sixteen also lack critical metadata and explicitly use
(Tc,Pc,rhoc)=(300 K,5 MPa,10000 mol/m3). These 2312 calls qualify supplied
initialization behavior, not default physical fluid metadata.

All 8839 finite original results match density, exact reported phase flag,
pressure and P_rho. All five paths execute: 4764 requested-phase Newton,
1881 initial-curvature switches, 1564 other-phase retries, 84 first fallbacks
and 546 legacy fallbacks. These are native diagnostics, not independently
instrumented source path counts. The source's 438 unclosed returns are included
in the passing raw comparison. The single source NaN return is retained in
metadata and checked for a native error.

Absolute tolerances are [1e-12 mol/m3, 0, 1e-5 Pa, 1e-8 Pa/(mol/m3)] with
2e-9 relative tolerance for the numerical fields and zero phase tolerance.
Comparisons do not widen tolerances for the bad iterates. The standalone source
snapshot test runs on process main to match the reference oracle's recorded
floating-point context; normal threaded tests cover domains, custom ideal
coefficients, ownership and prepared/convenience API equivalence.

This closes the raw-output and phase-choice gap identified in the previous
bounded-root qualification for the documented explicit initialization cases.
Inner trajectory identity, `meta_extrem`, full default initialization,
`MBWR_stableDensity`, full ideal/caloric properties, coexistence/flashes and CSP
integration remain pending within the full migration.

```
python tools/rust_port/validate_mbwr_density.py --compiler-directory /path/to/compiler/bin
cargo test --test mbwr_density --test mbwr_source_density
```
