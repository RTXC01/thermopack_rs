# Lennard-Jones splined model

The splined Lennard-Jones model of the fork approximates the excess energy of
the Lennard-Jones fluid with the Barker-Henderson perturbation theory. The
spline is a set of polynomials in the reduced temperature and the reduced
density, and the model adds the hard-sphere reference, the first, the second
and the third order term of the perturbation theory and the correction of the
second order term.

## Source scope

| Original routine | Rust |
| --- | --- |
| `src/lj_splined.f90:calc_ai_reduced_LJs_ex` | `ljs::LjsBh::reduced_perturbation_terms` |
| `src/lj_splined.f90:calc_ai_LJs_ex` | `ljs::LjsBh::spline_polynomial` |
| `src/lj_splined.f90:ljs_bh_get_bh_diameter_div_sigma` | `ljs::LjsBh::diameter_div_sigma` |
| `src/lj_splined.f90:init_LJs_bh` | `ljs::LjsBh::new`, `ljs::LjsBh::with_parameters` |
| `src/lj_splined.f90:ljs_bh_model_control` | `ljs::LjsBh::set_control`, `ljs::LjsControl` |
| `src/lj_splined.f90:calcFresLJs_bh` | `ljs::LjsBh::residual_helmholtz_energy` |
| `src/lj_splined.f90:calc_LJs_WCA_ai_tr` | `ljs_wca::perturbation_term` |
| `src/lj_splined.f90:calcFresLJs_bh` (volume derivative) | `ljs::LjsBh::pressure` |
| `src/lj_splined.f90:calcFresLJs_bh` (temperature derivative) | `ljs::LjsBh::residual_helmholtz_energy_derivative` |
| `src/lj_splined.f90:calcFresLJs_bh` (mole number derivative) | `ljs::LjsBh::fugacity_coefficient` |
| `src/lj_splined.f90:calcFresLJs_bh` (state functions) | `ljs::LjsBh::state` |
| `src/lj_splined.f90:specific_volume` | `ljs::LjsBh::volume` (native) |
| `src/lj_splined.f90:ljs_potential_reduced` | `ljs::LjsBh::reduced_potential` |
| `src/lj_splined.f90:calc_LJs_bh_zeta` | `ljs::LjsBh::packing_fraction` |

The native model provides the reduced perturbation terms and the reduced
Barker-Henderson hard-sphere diameter, the polynomials of the spline with
their partial derivatives, the parameters and the switches of the model, the
reduced residual Helmholtz energy with its temperature derivative, the pressure
of a state, the reduced residual chemical potential and the caloric state. The
fugacity of a state and the volume solve are not yet ported.

## Method

The reduced perturbation terms follow the source. The first order term is the
first polynomial of the spline with the prefactor of the source, the second
order term is the second polynomial multiplied by the correction polynomial and
the isothermal compressibility factor of the hard-sphere reference, and the
third order term is the third polynomial. The coefficients of the spline are
the four by fifteen array of the source, and the native model keeps them in
the order of that array.

The Barker-Henderson diameter comes from the diameter routine of the
SAFT-VR Mie model of the native library, with the segment diameter, the well
depth, the exponents and the prefactor of the source.

The reduced residual Helmholtz energy follows the source. The hard-sphere
reference is the Carnahan-Starling expression of the packing fraction of the
Barker-Henderson diameter, the first order term is the first polynomial of
the spline with the prefactor of the source, the second order term is the
second polynomial with the correction polynomial and the isothermal
compressibility factor of the hard-sphere reference, and the third order term
is the third polynomial. The source divides the three terms by the first, the
second and the third power of the temperature, and the native model follows
that convention.

The model control of the source can replace the third order term with the
Lafitte form, which follows the dispersion packing fraction and the polynomials
of the van der Waals attractive constant of the source. The native model
implements that form as well.

The database of the source holds the argon record with the segment diameter
3.42e-10 metre and the well depth 124.0 kelvin, and the native model uses the
same values.

## Validation

`tools/rust_port/validate_ljs.py` compares the native model against the
exact-source library. The oracle calls the Python interface of the fork, which
reaches the original routines through `get_pert_a` and
`get_bh_diameter_div_sigma`.

* 8 cases: the reduced temperature and the reduced density over the range of
  the dilute and the dense fluid.
* The compared quantities are the three reduced perturbation terms and the
  reduced diameter. All 8 cases match with a relative 1e-10 and an absolute
  1e-12, and the observed worst relative difference is 1.5e-12.

The committed fixture is `tests/fixtures/ljs.csv`; the test is `tests/ljs.rs`.

`tools/rust_port/validate_ljs_spline.py` compares the polynomials of the
spline and their partial derivatives against the exact-source library. The
reference probe calls the original `calc_ai_LJs_ex` with the reduced
temperature of the model.

* 8 cases: the four polynomials of the spline over the reduced temperature and
  the reduced density.
* The compared quantities are the polynomial, the first and the second
  derivative with respect to the reduced temperature, the first and the second
  derivative with respect to the reduced density and the mixed derivative. All
  8 cases match with a relative 1e-10, and the observed worst relative
  difference is 3.1e-14.

The committed fixture is `tests/fixtures/ljs_spline.csv`; the test is
`tests/ljs.rs`.

`tools/rust_port/validate_ljs_helmholtz.py` compares the reduced residual
Helmholtz energy against the exact-source library. The reference probe calls
the original `calcFresLJs_bh` with the active model of the fork.

* 12 cases: the temperature, the volume and the number of moles of the single
  component over the dilute and the dense fluid, and the same states with the
  Lafitte form of the third order term.
* The compared quantities are the reduced residual Helmholtz energy of the
  whole system and its temperature derivative. All 12 cases match with a
  relative 1e-9, and the observed worst relative difference is 4.4e-14.

The committed fixture is `tests/fixtures/ljs_helmholtz.csv`; the test is
`tests/ljs.rs`.

`tools/rust_port/validate_ljs_state.py` compares the state of the model
against the exact-source library. The reference oracle asks the Python
interface of the fork for the molar volume and the compressibility factor of a
state, and the native model evaluates the pressure and the compressibility
factor of the same state with the analytic volume derivative of the energy.

* 8 cases: the temperature, the pressure and the phase of the single
  component over the dilute and the dense fluid.
* The compared quantities are the pressure of the state and the compressibility
  factor. All 8 cases match with a relative 1e-6, and the observed worst
  relative difference is 1.9e-13.

The committed fixture is `tests/fixtures/ljs_state.csv`; the test is
`tests/ljs.rs`.

`tools/rust_port/validate_ljs_fugacity.py` compares the reduced residual
chemical potential against the exact-source library. The reference probe calls
the original `calcFresLJs_bh` with the derivative with respect to the number of
moles, and the native model evaluates the same derivative, the reduced
residual chemical potential of the single component.

* 8 cases: the temperature, the pressure and the phase of the single
  component over the dilute and the dense fluid.
* The compared quantity is the reduced residual chemical potential. All 8 cases
  match with a relative 1e-6, and the observed worst relative difference is
  1.1e-15.

The Python interface of the fork reports the logarithm of the fugacity
coefficient of this model as the reduced residual energy of the state, which is
a different convention from the reduced residual chemical potential. The
qualification therefore compares the derivative of the original with the same
derivative of the native model.

The committed fixture is `tests/fixtures/ljs_fugacity.csv`; the test is
`tests/ljs.rs`.

`tools/rust_port/validate_ljs_caloric.py` compares the caloric state against
the exact-source library. The reference oracle asks the Python interface of the
fork for the molar enthalpy and the molar entropy of the state.

* 8 cases: the temperature, the pressure and the phase of the single
  component over the dilute and the dense fluid.
* The compared quantities are the molar enthalpy and the molar entropy of the
  state. All 8 cases match with a relative 1e-8, and the observed worst
  relative difference is 4.2e-15.

The native model builds the ideal part of the state from the argon record and
the constant heat capacity of the original, and the residual part from the
temperature derivative of the reduced residual Helmholtz energy.

The committed fixture is `tests/fixtures/ljs_caloric.csv`; the test is
`tests/ljs.rs`.

`tools/rust_port/validate_ljs_potential.py` compares the reduced interaction
potential against the exact-source library. The reference probe calls the
original `ljs_potential_reduced`.

* 4 cases: the separations divided by the segment diameter over the attractive
  branch, the minimum and the smoothing of the cutoff.
* The compared quantities are the potential divided by the well depth and the
  reduced force. All 4 cases match with a relative 1e-12, and the observed
  worst relative difference is 2.8e-15.

The committed fixture is `tests/fixtures/ljs_potential.csv`; the test is
`tests/ljs.rs`.

## The perturbation terms of the White Bear family

The White Bear family replaces the perturbation terms of the
Barker-Henderson theory with the soft repulsion of the White Bear model. The
source evaluates the four terms with polynomials in the reduced distance of the
hard-sphere centres and the reduced density, and the native module
`ljs_wca` evaluates the same terms with their first and second derivatives.

`tools/rust_port/validate_ljs_wca_term.py` compares the native terms against
the exact-source library. The reference probe calls the original
`calc_LJs_WCA_ai_tr`.

* 8 cases: the four terms over the reduced temperature and the reduced
  density.
* The compared quantity is the term. All 8 cases match with a relative 1e-12,
  and the observed worst relative difference is zero.

The committed fixture is `tests/fixtures/ljs_wca_term.csv`; the test is
`tests/ljs_wca.rs`.

`tools/rust_port/validate_wca_term.py` compares the same terms with the
prefactor of the source, the well depth of the model instead of the reduced
form, against the original `calc_LJs_WCA_ai`.

* 5 cases: the index of the term, the well depth, the reduced distance and the
  reduced density.
* The compared quantities are the term and its five partial derivatives. All 5
  cases match with a relative 1e-12, and the observed worst relative difference
  is 1.0e-15.

The committed fixture is `tests/fixtures/wca_term.csv`; the test is
`tests/ljs_wca.rs`.

`tools/rust_port/validate_uv_virial_difference.py` compares the difference of
the second virial coefficient of the UV theory and the first order term of the
theory against the original `uv_delta_b2_overall` and `calc_uv_WCA`.

* 5 cases: the segment diameter, the well depth, the temperature and the
  reduced density.
* The compared quantities are the twelve values of the two routines. All 5 cases
  match with a relative 1e-11, and the observed worst relative difference is
  1.5e-15.

The committed fixture is `tests/fixtures/uv_virial_difference.csv`; the test is
`tests/uv_phi.rs`.

## The u-fraction of the UV theory

The UV theory of the source describes the Lennard-Jones fluid with the
u-fraction of the theory. The native module `ljs` provides the fraction of the
source and the temperature dependent form of it, with their partial
derivatives.

`tools/rust_port/validate_uv_phi.py` compares the native fraction against the
exact-source library. The reference probe calls the original `uv_phi`.

* 5 cases: the well depth, the temperature and the reduced density.
* The compared quantities are the fraction and its five partial derivatives.
  All 5 cases match with a relative 1e-12, and the observed worst relative
  difference is 4.1e-15.

The committed fixture is `tests/fixtures/uv_phi.csv`; the test is
`tests/uv_phi.rs`. The temperature dependent form of the source is a private
routine, so the native model provides it without a qualification of its own.

`tools/rust_port/validate_uv_terms.py` compares the first and the second order
terms of the UV theory against the exact-source library. The reference probe
calls the original `uv_a1_u` and `uv_a1_b2` with the hard-sphere diameter of
the Barker-Henderson reference.

* 5 cases: the segment diameter, the well depth, the temperature and the
  reduced density.
* The compared quantities are the nine values of the two routines. All 5 cases
  match with a relative 1e-11, and the observed worst relative difference is
  2.5e-16.

The committed fixture is `tests/fixtures/uv_terms.csv`; the test is
`tests/uv_phi.rs`.

`tools/rust_port/validate_uv_terms.py` also compares the q-function of the UV
theory, the hard-sphere second virial coefficient and the second virial
coefficient of the Lennard-Jones spline against the original `uv_q_of_T`,
`uv_hs_b2` and `ljs_b2`.

* 4 cases: the well depth and the temperature.
* The compared quantities are the nine values of the three routines. All 4 cases
  match with a relative 1e-11, and the observed worst relative difference is
  3.7e-16.

The committed fixture is `tests/fixtures/uv_virials.csv`; the test is
`tests/uv_phi.rs`.

## The cavity terms of the White Bear reference

The White Bear family adds the cavity of the hard-sphere reference and the soft
repulsion of the model. The source evaluates the hard-sphere diameter of the
reference, the radius of the minimum of the cavity integral, the cavity
correlation function of the hard spheres with its derivatives and the Boltzmann
factor of the Lennard-Jones potential with its derivatives. The native module
`hardsphere_wca` evaluates the same quantities.

`tools/rust_port/validate_white_bear.py` compares the native terms against the
exact-source library. The reference probe calls the original
`calc_dhs_WCA`, `minimum_r_cavity_int_LJ`, `cavity_distribution_function` and
`lj_boltzmann_factor`.

* 5 cases: the temperature, the well depth, the packing fraction and the
  reduced separation, on both sides of the packing fraction of 0.52 where the
  source changes the correlation.
* The compared quantities are the fifteen values of the four routines. All 5
  cases match with a relative 1e-11, and the observed worst relative difference
  is 1.6e-14.

The committed fixture is `tests/fixtures/white_bear.csv`; the test is
`tests/white_bear.rs`.

## Remaining scope

* The cavity integral of the source, and with it the soft repulsion of the
  White Bear model, is not yet ported, and the native model provides the
  perturbation terms and the cavity terms only.

* The volume solve of the model is native. The pressure of the model falls
  with the molar volume over parts of the range, because the source scales the
  perturbation terms with the well depth and the temperature, and the isotherm
  therefore has more than one crossing with a pressure. The original solves
  the volume with a bracketed solver that starts from the volume of the
  saturation state or from the ideal gas volume, and it selects one of the
  crossings. The native model brackets the crossings on a logarithmic grid and
  takes the smallest one for the liquid phase and the largest one for the vapor
  phase, which can be another crossing than the original picks.
* The WCA and the UV variants of the source, the hard-sphere reference of the
  White Bear family and the convolution of the fundamental measure theory are
  not yet ported.
* The WCA and the UV variants of the source, the hard-sphere reference of the
  White Bear family and the convolution of the fundamental measure theory are
  not yet ported.
* The models of the source use the diameter of the Barker-Henderson reference
  and the mixing rules of the SAFT-VR Mie family, and the native model covers
  the single component case only.
