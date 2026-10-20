# Equilibrium temperature/pressure response

`Cubic::equilibrium_response(&flash)` returns analytic equilibrium Cp, dV/dT
at fixed P, and dV/dP at fixed T. Phase amount and composition changes are
included. `EquilibriumResponse` provides derived Cv, the isentropic volume
pressure derivative, and internal-energy T/P derivatives. Pass the temperature
and pressure of the evaluated flash to those identities. These derivatives are
foundations for UV/SV flashes. The separate [sound implementation](sound.md)
uses the full SV matrix to also handle pure coexistence; UV/SV flashes remain
pending.

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

## Virial coefficients

The second and the third virial coefficients follow the definition of the
original, the pressure expansion

$$
\frac{P}{R T} = \rho + B \rho^2 + C \rho^3 + O(\rho^4)
$$

at the zero density. The original estimates them from the residual Helmholtz
energy at a molar volume of 1e7 cubic metre per mole, and the native method
`PreparedState::virial_coefficients` uses the same estimate.

`tools/rust_port/validate_virial.py` compares the native model against the
exact-source library. The oracle calls the Python interface of the fork, which
reaches the original `virial_coefficients`.

* 10 cases: methane, carbon dioxide and nitrogen, two mixtures, in the
  Peng-Robinson and the Soave-Redlich-Kwong equations of state.
* The second virial coefficient matches with a relative 1e-4 and an absolute
  1e-8 cubic metre per mole, and the observed worst relative difference is
  1.8e-4. The difference is the truncation error of the finite volume
  estimate, which the native model and the original share.
* The third virial coefficient of the original is the difference of two large
  terms of that estimate, so it carries its cancellation and it is much larger
  than the analytic value of the equation of state. The native value follows
  the analytic value, and the two differ. The fixture holds the value of the
  original as a reference and the test does not compare it.

The committed fixture is `tests/fixtures/virial.csv`; the test is
`tests/virial.rs`.

## The matrix of the second virial coefficients

The original `secondVirialCoeffMatrix` fills the matrix of the pure and the
cross second virial coefficients of the model. The native method
`Cubic::second_virial_coefficient_matrix` builds the same matrix: the pure
coefficients come from the pure fluids and the cross coefficient of a pair
follows from the virial coefficient of the equal mixture of the pair and the
pure coefficients.

`tools/rust_port/validate_virial_matrix.py` compares the native model against
the exact-source library. A reference probe calls the original routine.

* 6 cases: methane, the methane-propane mixture and the
  carbon dioxide-methane-nitrogen mixture in the Peng-Robinson equation of
  state and the methane-propane mixture in the Soave-Redlich-Kwong equation
  of state.
* The matrix matches with a relative 1e-3 and an absolute 1e-8 cubic metre
  per mole, and the observed worst relative difference is 1.5e-4. The
  difference is the truncation error of the finite volume estimate of the
  original, which the cross coefficients amplify.

The committed fixture is `tests/fixtures/virial_matrix.csv`; the test is
`tests/virial_matrix.rs`.