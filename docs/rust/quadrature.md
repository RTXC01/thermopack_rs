# Quadratures

The source offers the Gauss and the Gauss-Kronrod rules in seven orders and
integrates a function together with its differential. The native model offers
the same rules with the same nodes and the same weights.

## Source scope

| Original routine | Rust |
| --- | --- |
| `src/quadratures.f90:get_quadrature_positions` | `quadrature::Quadrature::positions` |
| `src/quadratures.f90:get_quadrature_weights` | `quadrature::Quadrature::weights` |
| `src/quadratures.f90:calc_quadrature_error` | `quadrature::quadrature_error` |
| `src/quadratures.f90:integrate_with_differential` | `quadrature::integrate_with_differential` |
| `src/quadratures.f90:integrate_with_differential_const_bound` | `quadrature::integrate_with_differential_const_bound` |

The nodes and the weights of the source are parameter arrays of
`src/quadratures.f90`. The script `tools/rust_port/generate_quadrature_data.py`
parses them and writes `rust/src/quadrature_data.rs`, in the same way that the
other generated databases of the native library are produced.

## Method

The rules are built the way the source builds them. The Gauss rules take every
second node of the Gauss-Kronrod rule of the matching order and mirror the
weights of the first half, and the Gauss-Kronrod rules use the nodes and the
weights of their own arrays. The integration evaluates the function at the
nodes, scales the values with the half width of the interval and accumulates
the differential of the differential with the derivative of the nodes, as the
source does.

The error estimate compares the integral of the rule with the integral of the
Gauss rule of the same order on the even nodes. The Gauss rules have no
estimate and report zero, as in the source.

## The unset last point of the odd Gauss rules

The source mirrors the nodes and the weights of the Gauss rules with a loop
over the first half of the points, which leaves the last point of the rules
with an odd number of points unset. The value depends on the stack of the
call, so the integral of the Gauss rules with fifteen and with twenty-five
points carries an undefined contribution. The native model follows the
construction of the source and leaves the last point at zero, and the
qualification compares the rules with an odd number of points up to that
point.

## Validation

`tools/rust_port/validate_quadrature.py` compares the native rules against the
exact-source library. The reference probe calls the original routines and
integrates the square of the integration variable with its differential over the
unit interval.

* 14 cases: the seven Gauss rules and the seven Gauss-Kronrod rules.
* The compared quantities are the nodes, the weights, the integral and the
  differential. All 14 cases match with a relative 1e-13, and the observed worst
  relative difference is 3.3e-16.

The committed fixture is `tests/fixtures/quadrature.csv`; the test is
`tests/quadrature.rs`.

## Remaining scope

* The adaptive integration of the source is a fixed rule with an error report,
  and the native model offers the same.
* The quadratures of the source are used by the convolution of the fundamental
  measure theory and by the White Bear models, and those parts are pending.
