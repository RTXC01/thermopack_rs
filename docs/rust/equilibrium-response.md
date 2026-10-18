# Equilibrium temperature/pressure response

`Cubic::equilibrium_response(&flash)` returns analytic equilibrium Cp, dV/dT
at fixed P, and dV/dP at fixed T. Phase amount and composition changes are
included. `EquilibriumResponse` provides derived Cv, the isentropic volume
pressure derivative, and internal-energy T/P derivatives. Pass the temperature
and pressure of the evaluated flash to those identities. These derivatives are
foundations for UV/SV flashes and equilibrium sound speed; those solvers are
not yet implemented by this change.

The response assumes a converged equilibrium split; it does not establish
stability. For a pure two-phase state, T and P are constrained to coexistence,
so independent TP response is undefined and the method returns `None`.
Endpoint splits and singular, non-finite or nonpositive-Cp responses also return
`None`. Non-finite input, invalid composition and invalid fractions return errors.

For a mixture split, let H be the chemical-potential Hessian in vapor-amount
coordinates, dT the vapor-minus-liquid log-fugacity temperature derivative, and
dP its pressure counterpart. Eliminating phase redistribution gives:

* `Cp = Cp_frozen + R*T^2 * dT^T H^-1 dT`.
* `V_T = V_T_frozen - R*T * dP^T H^-1 dT`.
* `V_P = V_P_frozen - R*T * dP^T H^-1 dP`.

Frozen phase volume derivatives follow by implicit differentiation of the EOS.
All corrections include temperature-dependent covolume and volume translation
through the existing analytic phase derivatives. Zero-feed components are
excluded from the response solve. The new routine does not change the existing
PH/PS solver or `equilibrium_properties` implementation.

`baseline_equilibrium_response` instead retains the fork's reported vdW and
type-10 H2 Cp conventions, including finite negative Cp. This explicit baseline
response need not differentiate returned H/S. It compares with the source's
`state_functions::dhdt_twoPhase`, `dvdt_twoPhase` and `dvdp_twoPhase` helpers;
single-phase calls compare with `eos::enthalpy` and `specificVolume` derivatives.
Their near-critical fallback trajectories are not fully reproduced by the
Schur-complement elimination, and full raw parity is not claimed beyond the
qualified cases.

The source grid attempts 441 cases across PR/SRK/RK/vdW; default and tcPR
construction; seven feeds including H2, He, water, zero fractions and a
hydrocarbon ternary; 180/300/500 K and 1 kPa/1 MPa/8 MPa. All 435 finite source
returns match at the fixture's recorded tolerances; six source terminations
remain recorded failures. Independent native central differences in 72 states
check physical H, V and U responses and `S_P = -V_T`, with nonzero kij/lij and
volume translation. Pure latent response and zero-feed handling have separate
tests. Broader quantum/critical/multiphase response qualification remains pending.
