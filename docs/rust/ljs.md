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
| `src/lj_splined.f90:calcFresLJs_bh` | not yet ported |

The native model provides the reduced perturbation terms and the reduced
Barker-Henderson hard-sphere diameter, the polynomials of the spline with
their partial derivatives, and the parameters of the model. The residual
Helmholtz energy of the model, with its temperature, volume and composition
derivatives, is not yet ported.

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

## Remaining scope

* The residual Helmholtz energy of the model and its derivatives are not yet
  ported, and the native model does not yet evaluate a state of the splined
  Lennard-Jones fluid.
* The WCA and the UV variants of the source, the hard-sphere reference of the
  White Bear family and the convolution of the fundamental measure theory are
  not yet ported.
* The models of the source use the diameter of the Barker-Henderson reference
  and the mixing rules of the SAFT-VR Mie family, and the native model covers
  the single component case only.
