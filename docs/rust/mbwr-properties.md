# MBWR residual state functions and TP responses

`Mbwr::residual_tp(T,total_V,n)` evaluates a homogeneous state with pressure
determined by the EOS. It returns residual G/H (J), S (J/K), compressibility Z
and ln(phi). These residuals compare the real and ideal states at the same
T/P/n. They require positive pressure; negative-pressure raw Helmholtz/pressure
kernels remain separately available. A local TP response does not establish
phase stability or select an equilibrium density root.

H, S, Z and ln(phi) use `MbwrTpQuantity`: a scalar value and an optional
`MbwrTpDerivatives { dt, dp, dn }` block, holding the other T/P/n variables
fixed. H/S are extensive; Z/ln(phi) are intensive. Physical pure n responses
are H/n, S/n and exactly zero for the two intensive quantities. The ideal-gas
Cp/reference database is unnecessary for these residual-only functions.

```rust
use thermopack_rs::{Mbwr, MbwrKind};
let model = Mbwr::new(MbwrKind::Mbwr32, "C1")?;
let state = model.residual_tp(300.0, 1e-3, 1.0)?;
assert!((state.enthalpy.value - 300.0 * state.entropy.value
         - state.gibbs_energy).abs() < 1e-7);
let derivative = state.entropy.derivatives.unwrap();
# Ok::<(), thermopack_rs::Error>(())
```

The implementation shares one residual Helmholtz snapshot across all state
functions. Physical values use the stable integral and `ln_1p` for the dilute
compressibility correction. G=nRT*ln(phi) and H-TS=G are checked independently.
Pressure response gives V_T=-P_T/P_V; physical S_T=H_T/T and S_P=nR/P-V_T.
At a singular or roundoff-unresolved P_V, H/S/Z keep their finite values with
`None` derivatives. Pure fugacity first derivatives can remain defined there.
The physical roundoff threshold is 64 machine epsilons times nRT/V^2.
Non-finite scalar values and invalid input domains return errors.

## Original source conventions

`baseline_residual_tp(T,P,total_V,n)` ports `MBWR_Gres`, `MBWR_Hres`,
`MBWR_Sres`, `MBWR_zfac` and `MBWR_lnphi`, including their pressure expressions
and arithmetic order. It accepts an explicit positive P, as those original
routines do, even when it does not close the EOS. Its derivatives reproduce
the fork's expressions at that input; an off-EOS input is not asserted to be
a closed thermodynamic state. Baseline response blocks become `None` when
non-finite, without suppressing otherwise finite scalar values.

The source entropy routine has pressure-response signs different from the
other routines. In the ideal-coefficient limit, Sres=0 but its T/P derivatives
are -2*n*R/T and 2*V/T. The baseline method preserves these; the physical
method has zero ideal residual responses. Source units in the comments are
occasionally misleading: actual code and extensive-scaling qualification use
J and J/K for total G/H/S, not kJ/mol.

Pure Z_n and ln(phi)_n are mathematically zero, but the original expressions
subtract rounded large intermediates. Their small residues can therefore
depend strongly on the host math runtime and floating environment at extreme
extrapolations. On this Windows host, exp(-0.96) differs by one ulp between
process-main and newly created Rust worker threads. That is amplified in
these near-zero residues. The source snapshot executable has `harness=false`
and runs on process main, matching the reference oracle and serial parity
driver. Normal worker-thread tests verify the physical derivatives and exact
intensivity. Production code does not change the caller's floating environment;
raw cancellation residues are not promised to be bitwise portable across math
runtimes or execution contexts.

## Qualification

The reference adapter invokes the five original routines without modifying
Fortran sources. It shares the same coefficient-only initialization as the
[kernel qualification](mbwr.md); no critical or ideal metadata is needed.
All 5983 finite states match 18 outputs through the baseline API, including
all original derivative expressions and tiny n-response residues in the
recorded process-main context. The 7800 attempted calls span every stored fit,
twelve T values from 0.5 to 2000 K, eight densities from 1e-5 to 20000 mol/m3,
and n=0.3/1/3.5. Another 1040 requests within that total intentionally supply
P=abs(P_EOS)*factor+1000 Pa with factor 0.8/1.2 to qualify explicit P behavior.

All 1817 other calls have nonpositive EOS pressure and non-finite logarithmic
TP residuals. Their requests and raw outputs remain in metadata. They are not
parity successes or silently dropped. Quantity-specific absolute tolerances
and 2e-7 relative tolerance are fixed in metadata; no tolerance was widened to
hide entropy signs or roundoff artifacts.

Physical responses are independently checked by local pressure-root finite
differences for every fit, including all T/P/n derivatives of H/S/Z/ln(phi),
G_T=-S, G_P=V-nRT/P, H-TS=G, G=nRT*ln(phi), extensive scaling and perfect-gas
limits. A custom singular EOS tests retained values/optional responses.
Their local Newton helper is test-only. Bounded physical density roots and
critical points have separate [equilibrium qualification](mbwr-equilibrium.md).
Source density-iteration fidelity, ideal/full caloric properties and flash/CSP
integration remain pending. No new speedup is claimed here.

```
python tools/rust_port/validate_mbwr_properties.py --compiler-directory /path/to/compiler/bin
cargo test --test mbwr_properties --test mbwr_source_properties
```
