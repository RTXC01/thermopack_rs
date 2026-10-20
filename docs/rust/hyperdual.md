# Hyperdual numbers

The crate provides hyperdual numbers with up to three differentiation
directions, ported from `src/hyperdual_mod.f90`.

## Model

A hyperdual number carries the value and the partial derivatives with respect
to up to three independent directions:

```text
f = f0 + f1 e1 + f2 e2 + f3 e3 + f12 e1 e2 + f13 e1 e3 + f23 e2 e3 + f123 e1 e2 e3
```

The overall order of the differentiation travels with the number, so a number
of order one never evaluates the higher terms and a number of order three
evaluates all of them. The default order is two, and the third order has to be
activated explicitly.

The module provides the elementary functions of the calculus with their
hyperdual counterparts:

| Function | Derivatives |
|---|---|
| `exp` | value and the first, second and third order derivatives |
| `ln`, `log10` | the inverse value with the chain rule |
| `sqrt` | the real power with the source's guard against a vanishing base |
| `sin`, `cos`, `tan` | the trigonometric functions |
| `sinh`, `cosh`, `tanh` | the hyperbolic functions, built from the exponential |
| `asin`, `acos`, `atan` | the inverse trigonometric functions |
| `atan2` | the arc tangent of a ratio, with the first order derivatives |

The operators follow the source: the addition and the subtraction add the
components, the multiplication applies the product rule, the division
multiplies with the inverse and the integer power repeats the multiplication.
The absolute value keeps the sign of the derivatives.

## Source limitation reproduced

The third order derivative of the arc tangent in the original writes
`deriv**3**v1%f1`, which is a double power: the exponent of the derivative is
three raised to the first derivative of the argument. The native model
preserves that expression, and the qualification confirms that the resulting
component agrees with the exact-source library.

## Validation

`tools/rust_port/validate_hyperdual.py` compares the native functions against
the exact-source library through a reference probe that evaluates the original
module.

* 6 states covering the first, second, third and zero differentiation order,
  with twelve functions and eight components each.
* All 6 states agree with the exact-source library, and the largest observed
  relative difference is zero.
* The committed fixture is `tests/fixtures/hyperdual.csv`; the test is
  `tests/hyperdual.rs`.

## Remaining scope

* The general-purpose quadratures of `src/quadratures.f90` are not ported. The
  SAFT-VR Mie model uses its own allocation-free Gauss-Legendre quadrature.
* The hyperdual module is not yet used by the mixture models, which compute
  their derivatives analytically.
