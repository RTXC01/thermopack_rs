//! Hyperdual numbers with up to three differentiation directions.
//! Ported from `src/hyperdual_mod.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! A hyperdual number carries the value and the partial derivatives with
//! respect to up to three independent directions:
//!
//! ```text
//! f = f0 + f1 e1 + f2 e2 + f3 e3 + f12 e1 e2 + f13 e1 e3 + f23 e2 e3 + f123 e1 e2 e3
//! ```
//!
//! The overall order of the differentiation is kept with the number, so a
//! number of order one never evaluates the higher terms.

use std::ops::{Add, Div, Mul, Neg, Sub};

/// A hyperdual number with up to three differentiation directions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HyperDual {
    /// The value.
    pub f0: f64,
    /// Derivative with respect to the first direction.
    pub f1: f64,
    /// Derivative with respect to the second direction.
    pub f2: f64,
    /// Derivative with respect to the third direction.
    pub f3: f64,
    /// Mixed derivative of the first and second directions.
    pub f12: f64,
    /// Mixed derivative of the first and third directions.
    pub f13: f64,
    /// Mixed derivative of the second and third directions.
    pub f23: f64,
    /// Mixed derivative of all three directions.
    pub f123: f64,
    /// The overall order of the differentiation. The default is two, and three
    /// has to be activated explicitly.
    pub order: u8,
}

impl Default for HyperDual {
    fn default() -> Self {
        Self::real(0.0)
    }
}

impl HyperDual {
    /// A hyperdual number without derivatives.
    pub fn real(value: f64) -> Self {
        Self {
            f0: value,
            f1: 0.0,
            f2: 0.0,
            f3: 0.0,
            f12: 0.0,
            f13: 0.0,
            f23: 0.0,
            f123: 0.0,
            order: 2,
        }
    }

    /// A hyperdual number of the given differentiation order.
    pub fn with_order(value: f64, order: u8) -> Self {
        let mut number = Self::real(value);
        number.order = order;
        number
    }

    /// A hyperdual number with a first order derivative.
    pub fn first(value: f64, derivative: f64) -> Self {
        let mut number = Self::real(value);
        number.f1 = derivative;
        number
    }

    /// A hyperdual number with two first order derivatives.
    pub fn second(value: f64, first: f64, second: f64) -> Self {
        let mut number = Self::real(value);
        number.f1 = first;
        number.f2 = second;
        number
    }

    /// The absolute value, keeping the sign of the derivatives.
    pub fn abs(&self) -> Self {
        let mut result = *self;
        if self.f0 < 0.0 {
            result.f0 = -self.f0;
            if result.order > 0 {
                result.f1 = -self.f1;
            }
            if result.order > 1 {
                result.f2 = -self.f2;
                result.f12 = -self.f12;
            }
            if result.order == 3 {
                result.f3 = -self.f3;
                result.f13 = -self.f13;
                result.f23 = -self.f23;
                result.f123 = -self.f123;
            }
        }
        result
    }

    /// The real power, following the source's guard against a vanishing base.
    pub fn pow_real(&self, exponent: f64) -> Self {
        const TOLERANCE: f64 = 1.0e-15;
        let mut base = self.f0;
        if base.abs() < TOLERANCE {
            base = if base >= 0.0 { TOLERANCE } else { -TOLERANCE };
        }
        let mut result = Self::real(self.f0.powf(exponent));
        result.order = self.order;
        if result.order > 0 {
            let derivative = exponent * base.powf(exponent - 1.0);
            result.f1 = self.f1 * derivative;
            if result.order > 1 {
                let second = exponent * (exponent - 1.0) * base.powf(exponent - 2.0);
                result.f2 = self.f2 * derivative;
                result.f12 = self.f12 * derivative + self.f1 * self.f2 * second;
                if result.order == 3 {
                    result.f3 = self.f3 * derivative;
                    result.f13 = self.f13 * derivative + self.f1 * self.f3 * second;
                    result.f23 = self.f23 * derivative + self.f2 * self.f3 * second;
                    result.f123 = self.f123 * derivative
                        + second * (self.f12 * self.f3 + self.f13 * self.f2 + self.f23 * self.f1)
                        + exponent
                            * (exponent - 1.0)
                            * (exponent - 2.0)
                            * base.powf(exponent - 3.0)
                            * self.f1
                            * self.f2
                            * self.f3;
                }
            }
        }
        result
    }

    /// The integer power.
    pub fn pow_int(&self, exponent: i32) -> Self {
        let mut result = Self::real(1.0);
        for _ in 0..exponent.unsigned_abs() {
            result = result * *self;
        }
        if exponent < 0 {
            result = Self::real(1.0) / result;
        }
        result
    }

    /// The exponential.
    pub fn exp(&self) -> Self {
        let value = self.f0.exp();
        let mut result = Self::real(value);
        result.order = self.order;
        if result.order > 0 {
            result.f1 = value * self.f1;
        }
        if result.order > 1 {
            result.f2 = value * self.f2;
            result.f12 = value * (self.f12 + self.f1 * self.f2);
        }
        if result.order == 3 {
            result.f3 = value * self.f3;
            result.f13 = value * (self.f13 + self.f1 * self.f3);
            result.f23 = value * (self.f23 + self.f2 * self.f3);
            result.f123 = value
                * (self.f123
                    + self.f12 * self.f3
                    + self.f13 * self.f2
                    + self.f23 * self.f1
                    + self.f1 * self.f2 * self.f3);
        }
        result
    }

    /// The natural logarithm.
    pub fn ln(&self) -> Self {
        let mut result = Self::real(self.f0.ln());
        result.order = self.order;
        if result.order > 0 {
            let inverse = 1.0 / self.f0;
            let first = self.f1 * inverse;
            result.f1 = first;
            if result.order > 1 {
                let second = self.f2 * inverse;
                let mixed = self.f12 * inverse;
                result.f2 = second;
                result.f12 = mixed - (first * second);
                if result.order == 3 {
                    let third = self.f3 * inverse;
                    let mixed13 = self.f13 * inverse;
                    let mixed23 = self.f23 * inverse;
                    result.f3 = third;
                    result.f13 = mixed13 - (first * third);
                    result.f23 = mixed23 - (second * third);
                    result.f123 = self.f123 * inverse
                        - (mixed * third)
                        - (mixed13 * second)
                        - (mixed23 * first)
                        + 2.0 * (first * second * third);
                }
            }
        }
        result
    }

    /// The base ten logarithm.
    pub fn log10(&self) -> Self {
        self.ln() / 10f64.ln()
    }

    /// The square root.
    pub fn sqrt(&self) -> Self {
        self.pow_real(0.5)
    }

    /// The sine.
    pub fn sin(&self) -> Self {
        let value = self.f0.sin();
        let mut result = Self::real(value);
        result.order = self.order;
        if result.order > 0 {
            let derivative = self.f0.cos();
            result.f1 = derivative * self.f1;
            if result.order > 1 {
                result.f2 = derivative * self.f2;
                result.f12 = derivative * self.f12 - value * self.f1 * self.f2;
                if result.order == 3 {
                    result.f3 = derivative * self.f3;
                    result.f13 = derivative * self.f13 - value * self.f1 * self.f3;
                    result.f23 = derivative * self.f23 - value * self.f2 * self.f3;
                    result.f123 = derivative * (self.f123 - self.f1 * self.f2 * self.f3)
                        - value * (self.f12 * self.f3 - self.f13 * self.f2 - self.f23 * self.f1);
                }
            }
        }
        result
    }

    /// The cosine.
    pub fn cos(&self) -> Self {
        let value = self.f0.cos();
        let mut result = Self::real(value);
        result.order = self.order;
        if result.order > 0 {
            let derivative = -self.f0.sin();
            result.f1 = derivative * self.f1;
            if result.order > 1 {
                result.f2 = derivative * self.f2;
                result.f12 = derivative * self.f12 - value * self.f1 * self.f2;
                if result.order == 3 {
                    result.f3 = derivative * self.f3;
                    result.f13 = derivative * self.f13 - value * self.f1 * self.f3;
                    result.f23 = derivative * self.f23 - value * self.f2 * self.f3;
                    result.f123 = derivative * (self.f123 - self.f1 * self.f2 * self.f3)
                        - value * (self.f12 * self.f3 + self.f13 * self.f2 + self.f23 * self.f1);
                }
            }
        }
        result
    }

    /// The tangent.
    pub fn tan(&self) -> Self {
        let value = self.f0.tan();
        let mut result = Self::real(value);
        result.order = self.order;
        if result.order > 0 {
            let derivative = value * value + 1.0;
            result.f1 = derivative * self.f1;
            if result.order > 1 {
                result.f2 = derivative * self.f2;
                result.f12 = derivative * self.f12 + self.f1 * self.f2 * 2.0 * value * derivative;
                if result.order == 3 {
                    result.f3 = derivative * self.f3;
                    result.f13 =
                        derivative * self.f13 + self.f1 * self.f3 * 2.0 * value * derivative;
                    result.f23 =
                        derivative * self.f23 + self.f2 * self.f3 * 2.0 * value * derivative;
                    result.f123 = derivative * self.f123
                        + (self.f12 * self.f3 + self.f13 * self.f2 + self.f23 * self.f1)
                            * 2.0
                            * value
                            * derivative
                        + 2.0 * derivative * (3.0 * derivative - 2.0) * self.f1 * self.f2 * self.f3;
                }
            }
        }
        result
    }

    /// The hyperbolic sine.
    pub fn sinh(&self) -> Self {
        let first = self.exp();
        let second = (-*self).exp();
        (first - second) * 0.5
    }

    /// The hyperbolic cosine.
    pub fn cosh(&self) -> Self {
        let first = self.exp();
        let second = (-*self).exp();
        (first + second) * 0.5
    }

    /// The hyperbolic tangent.
    pub fn tanh(&self) -> Self {
        let first = self.exp();
        let second = (-*self).exp();
        (first - second) / (first + second)
    }

    /// The arc sine.
    pub fn asin(&self) -> Self {
        let mut result = Self::real(self.f0.asin());
        result.order = self.order;
        if result.order > 0 {
            let inverse = 1.0 / (1.0 - self.f0 * self.f0);
            let derivative = inverse.sqrt();
            result.f1 = derivative * self.f1;
            if result.order > 1 {
                let second = self.f0 * inverse * derivative;
                result.f2 = derivative * self.f2;
                result.f12 = derivative * self.f12 + second * self.f1 * self.f2;
                if result.order == 3 {
                    result.f3 = derivative * self.f3;
                    result.f13 = derivative * self.f13 + second * self.f1 * self.f3;
                    result.f23 = derivative * self.f23 + second * self.f2 * self.f3;
                    result.f123 = derivative * self.f123
                        + second * (self.f12 * self.f3 + self.f13 * self.f2 + self.f23 * self.f1)
                        + (2.0 * self.f0 * self.f0 + 1.0)
                            * inverse
                            * inverse
                            * derivative
                            * self.f1
                            * self.f2
                            * self.f3;
                }
            }
        }
        result
    }

    /// The arc cosine.
    pub fn acos(&self) -> Self {
        let mut result = Self::real(self.f0.acos());
        result.order = self.order;
        if result.order > 0 {
            let inverse = 1.0 / (1.0 - self.f0 * self.f0);
            let derivative = -inverse.sqrt();
            result.f1 = derivative * self.f1;
            if result.order > 1 {
                let second = self.f0 * inverse * derivative;
                result.f2 = derivative * self.f2;
                result.f12 = derivative * self.f12 + second * self.f1 * self.f2;
                if result.order == 3 {
                    result.f3 = derivative * self.f3;
                    result.f13 = derivative * self.f13 + second * self.f1 * self.f3;
                    result.f23 = derivative * self.f23 + second * self.f2 * self.f3;
                    result.f123 = derivative * self.f123
                        + second * (self.f12 * self.f3 + self.f13 * self.f2 + self.f23 * self.f1)
                        + (2.0 * self.f0 * self.f0 + 1.0)
                            * inverse
                            * inverse
                            * derivative
                            * self.f1
                            * self.f2
                            * self.f3;
                }
            }
        }
        result
    }

    /// The arc tangent.
    pub fn atan(&self) -> Self {
        let value = self.f0.atan();
        let mut result = Self::real(value);
        result.order = self.order;
        if result.order > 0 {
            let inverse = 1.0 / (1.0 + self.f0 * self.f0);
            let derivative = inverse;
            result.f1 = derivative * self.f1;
            if result.order > 1 {
                let second = -2.0 * self.f0 * inverse * derivative;
                result.f2 = derivative * self.f2;
                result.f12 = derivative * self.f12 + second * self.f1 * self.f2;
                if result.order == 3 {
                    result.f3 = derivative * self.f3;
                    result.f13 = derivative * self.f13 + second * self.f1 * self.f3;
                    result.f23 = derivative * self.f23 + second * self.f2 * self.f3;
                    // The original writes `deriv**3**v1%f1` here, which is a
                    // double power: the exponent of the derivative is
                    // 3 raised to the first derivative of the argument. The
                    // original behavior is preserved.
                    result.f123 = derivative * self.f123
                        + second * (self.f12 * self.f3 + self.f13 * self.f2 + self.f23 * self.f1)
                        + (6.0 * self.f0 * self.f0 - 2.0)
                            * derivative.powf(3.0_f64.powf(self.f1))
                            * self.f2
                            * self.f3;
                }
            }
        }
        result
    }

    /// The arc tangent of the ratio of two numbers, using the first order
    /// derivatives only.
    pub fn atan2(&self, other: &HyperDual) -> HyperDual {
        let (a, b) = (self.f0, self.f1);
        let (c, d) = (other.f0, other.f1);
        let mut result = Self::real(a.atan2(c));
        result.f1 = (c * b - a * d) / (a * a + c * c);
        result
    }
}

impl Add for HyperDual {
    type Output = HyperDual;
    fn add(self, other: HyperDual) -> HyperDual {
        let mut result = self;
        result.order = self.order.max(other.order);
        result.f0 = self.f0 + other.f0;
        if result.order > 0 {
            result.f1 = self.f1 + other.f1;
        }
        if result.order > 1 {
            result.f2 = self.f2 + other.f2;
            result.f12 = self.f12 + other.f12;
        }
        if result.order == 3 {
            result.f3 = self.f3 + other.f3;
            result.f13 = self.f13 + other.f13;
            result.f23 = self.f23 + other.f23;
            result.f123 = self.f123 + other.f123;
        }
        result
    }
}

impl Add<f64> for HyperDual {
    type Output = HyperDual;
    fn add(self, other: f64) -> HyperDual {
        let mut result = self;
        result.f0 = self.f0 + other;
        result
    }
}

impl Sub for HyperDual {
    type Output = HyperDual;
    fn sub(self, other: HyperDual) -> HyperDual {
        let mut result = self;
        result.order = self.order.max(other.order);
        result.f0 = self.f0 - other.f0;
        if result.order > 0 {
            result.f1 = self.f1 - other.f1;
        }
        if result.order > 1 {
            result.f2 = self.f2 - other.f2;
            result.f12 = self.f12 - other.f12;
        }
        if result.order == 3 {
            result.f3 = self.f3 - other.f3;
            result.f13 = self.f13 - other.f13;
            result.f23 = self.f23 - other.f23;
            result.f123 = self.f123 - other.f123;
        }
        result
    }
}

impl Neg for HyperDual {
    type Output = HyperDual;
    fn neg(self) -> HyperDual {
        let mut result = self;
        result.f0 = -self.f0;
        if self.order > 0 {
            result.f1 = -self.f1;
        }
        if self.order > 1 {
            result.f2 = -self.f2;
            result.f12 = -self.f12;
        }
        if self.order == 3 {
            result.f3 = -self.f3;
            result.f13 = -self.f13;
            result.f23 = -self.f23;
            result.f123 = -self.f123;
        }
        result
    }
}

impl Mul for HyperDual {
    type Output = HyperDual;
    fn mul(self, other: HyperDual) -> HyperDual {
        let mut result = Self::real(0.0);
        result.f0 = self.f0 * other.f0;
        result.order = self.order.max(other.order);
        if result.order > 0 {
            result.f1 = self.f0 * other.f1 + self.f1 * other.f0;
        }
        if result.order > 1 {
            result.f2 = self.f0 * other.f2 + self.f2 * other.f0;
            result.f12 =
                self.f0 * other.f12 + self.f1 * other.f2 + self.f2 * other.f1 + self.f12 * other.f0;
        }
        if result.order == 3 {
            result.f3 = self.f0 * other.f3 + self.f3 * other.f0;
            result.f13 =
                self.f0 * other.f13 + self.f1 * other.f3 + self.f3 * other.f1 + self.f13 * other.f0;
            result.f23 =
                self.f0 * other.f23 + self.f2 * other.f3 + self.f3 * other.f2 + self.f23 * other.f0;
            result.f123 = self.f0 * other.f123
                + self.f12 * other.f3
                + self.f3 * other.f12
                + self.f123 * other.f0
                + self.f23 * other.f1
                + self.f1 * other.f23
                + self.f13 * other.f2
                + self.f2 * other.f13;
        }
        result
    }
}

impl Mul<f64> for HyperDual {
    type Output = HyperDual;
    fn mul(self, other: f64) -> HyperDual {
        let mut result = self;
        result.f0 = self.f0 * other;
        if self.order > 0 {
            result.f1 = self.f1 * other;
        }
        if self.order > 1 {
            result.f2 = self.f2 * other;
            result.f12 = self.f12 * other;
        }
        if self.order == 3 {
            result.f3 = self.f3 * other;
            result.f13 = self.f13 * other;
            result.f23 = self.f23 * other;
            result.f123 = self.f123 * other;
        }
        result
    }
}

impl Div for HyperDual {
    type Output = HyperDual;
    fn div(self, other: HyperDual) -> HyperDual {
        self * other.pow_real(-1.0)
    }
}

impl Div<f64> for HyperDual {
    type Output = HyperDual;
    fn div(self, other: f64) -> HyperDual {
        self * (1.0 / other)
    }
}

impl Sub<f64> for HyperDual {
    type Output = Self;
    fn sub(self, rhs: f64) -> Self {
        self + (-rhs)
    }
}
impl Add<HyperDual> for f64 {
    type Output = HyperDual;
    fn add(self, rhs: HyperDual) -> HyperDual {
        rhs + self
    }
}
impl Sub<HyperDual> for f64 {
    type Output = HyperDual;
    fn sub(self, rhs: HyperDual) -> HyperDual {
        -rhs + self
    }
}
impl Mul<HyperDual> for f64 {
    type Output = HyperDual;
    fn mul(self, rhs: HyperDual) -> HyperDual {
        rhs * self
    }
}
impl Div<HyperDual> for f64 {
    type Output = HyperDual;
    fn div(self, rhs: HyperDual) -> HyperDual {
        rhs.pow_real(-1.) * self
    }
}
impl std::ops::AddAssign for HyperDual {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}
impl std::ops::SubAssign for HyperDual {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

/// The sum of a slice of hyperdual numbers.
pub fn sum(values: &[HyperDual]) -> HyperDual {
    let mut total = HyperDual::real(0.0);
    for value in values {
        total = total + *value;
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_first_derivative_follows_the_function() {
        let x = HyperDual::first(2.0, 1.0);
        let value = x.exp();
        assert!((value.f0 - 2.0_f64.exp()).abs() < 1e-15);
        assert!((value.f1 - 2.0_f64.exp()).abs() < 1e-15);
    }

    #[test]
    fn the_second_derivative_follows_the_function() {
        let x = HyperDual::second(1.5, 1.0, 2.0);
        let value = x.ln();
        // The chain rule multiplies the derivative of the argument with the
        // derivative of the logarithm.
        assert!((value.f1 - 1.0 / 1.5).abs() < 1e-15);
        assert!((value.f2 - 2.0 / 1.5).abs() < 1e-15);
        // The mixed derivative carries the product rule of the two directions.
        let expected_mixed = -(1.0 / 1.5) * (2.0 / 1.5);
        assert!((value.f12 - expected_mixed).abs() < 1e-14);
    }

    #[test]
    fn the_third_order_derivative_follows_the_function() {
        let mut x = HyperDual::with_order(1.25, 3);
        x.f1 = 1.0;
        x.f2 = 1.0;
        x.f3 = 1.0;
        let value = x.exp();
        let base = 1.25_f64.exp();
        assert!((value.f0 - base).abs() < 1e-15);
        assert!((value.f1 - base).abs() < 1e-15);
        assert!((value.f2 - base).abs() < 1e-15);
        assert!((value.f3 - base).abs() < 1e-15);
        // The mixed derivative of the exponential carries the product rule of
        // the two directions.
        assert!((value.f12 - base).abs() < 1e-14);
        assert!((value.f13 - base).abs() < 1e-14);
        assert!((value.f23 - base).abs() < 1e-14);
        assert!((value.f123 - base).abs() < 1e-13);
    }

    #[test]
    fn the_multiplication_follows_the_product_rule() {
        let first = HyperDual::second(3.0, 1.0, 2.0);
        let second = HyperDual::second(4.0, 5.0, 6.0);
        let product = first * second;
        assert!((product.f0 - 12.0).abs() < 1e-15);
        assert!((product.f1 - (4.0 * 1.0 + 3.0 * 5.0)).abs() < 1e-15);
        assert!((product.f2 - (4.0 * 2.0 + 3.0 * 6.0)).abs() < 1e-15);
        // The mixed derivative follows the product rule of the two directions.
        let expected_mixed = first.f0 * second.f12
            + first.f1 * second.f2
            + first.f2 * second.f1
            + first.f12 * second.f0;
        assert!((product.f12 - expected_mixed).abs() < 1e-14);
        assert!(expected_mixed.abs() > 1.0);
    }

    #[test]
    fn the_order_limits_the_evaluation() {
        let mut x = HyperDual::with_order(1.0, 1);
        x.f1 = 1.0;
        x.f2 = 1.0;
        x.f12 = 1.0;
        let value = x.exp();
        assert_eq!(value.f0, 1.0_f64.exp());
        assert_eq!(value.f1, 1.0_f64.exp());
        assert_eq!(value.f2, 0.0);
        assert_eq!(value.f12, 0.0);
    }

    #[test]
    fn the_inverse_and_the_power_agree() {
        let x = HyperDual::first(4.0, 1.0);
        let inverse = HyperDual::real(1.0) / x;
        assert!((inverse.f0 - 0.25).abs() < 1e-15);
        assert!((inverse.f1 + 0.0625).abs() < 1e-15);
        let squared = x.pow_int(2);
        assert!((squared.f0 - 16.0).abs() < 1e-14);
        assert!((squared.f1 - 8.0).abs() < 1e-14);
    }
}
