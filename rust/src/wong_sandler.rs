//! The Wong-Sandler mixing rule of the cubic equation of state.
//!
//! The original `src/wong_sandler.f90` replaces the binary attraction and the
//! binary covolume of the cubic equation of state with its own coefficients.
//! The rule needs the binary interaction parameters and the parameters of the
//! second order term of the perturbation theory, both of which the original
//! stores as a fraction of the temperature.

use crate::excess_gibbs::ExcessGibbsJet;
use crate::{Error, Result};

/// The degree of the polynomials of a fraction.
const DEGREE: usize = 2;

/// A fraction of a temperature dependent parameter.
///
/// The value is the ratio of two quadratic polynomials in the temperature, as
/// in the original `cubic_eos::Fraction`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fraction {
    /// The coefficients of the numerator polynomial.
    pub numerator: [f64; DEGREE + 1],
    /// The coefficients of the denominator polynomial.
    pub denominator: [f64; DEGREE + 1],
}

impl Default for Fraction {
    /// The zero fraction, which the original uses for a missing parameter.
    fn default() -> Self {
        Self {
            numerator: [0.0; DEGREE + 1],
            denominator: [1.0, 0.0, 0.0],
        }
    }
}

impl Fraction {
    /// A constant fraction.
    pub fn constant(value: f64) -> Self {
        Self {
            numerator: [value, 0.0, 0.0],
            denominator: [1.0, 0.0, 0.0],
        }
    }

    /// The value of the fraction and its first and second derivatives.
    ///
    /// The original `excess_gibbs::GetFraction` returns zero when the
    /// denominator vanishes.
    pub fn evaluate(&self, x: f64) -> (f64, f64, f64) {
        let (numerator, numerator_d, numerator_dd) = polynomial(&self.numerator, x);
        let (denominator, denominator_d, denominator_dd) = polynomial(&self.denominator, x);
        if denominator == 0.0 {
            return (0.0, 0.0, 0.0);
        }
        let value = numerator / denominator;
        let derivative =
            (numerator_d * denominator - numerator * denominator_d) / (denominator * denominator);
        let second = (numerator_dd * denominator * denominator
            - numerator * denominator_dd * denominator
            - 2.0 * numerator_d * denominator * denominator_d
            + 2.0 * numerator * denominator_d * denominator_d)
            / (denominator * denominator * denominator);
        (value, derivative, second)
    }
}

/// The value of a quadratic polynomial and its first and second derivatives.
///
/// The original `excess_gibbs::GetPoly` evaluates the polynomial and the two
/// derivatives.
fn polynomial(coefficients: &[f64; DEGREE + 1], x: f64) -> (f64, f64, f64) {
    if x == 0.0 {
        return (coefficients[0], coefficients[1], 0.5 * coefficients[2]);
    }
    let mut value = 0.0;
    let mut derivative = 0.0;
    let mut second = 0.0;
    for (power, &coefficient) in coefficients.iter().enumerate() {
        let term = x.powi(power as i32);
        let power = power as f64;
        value += coefficient * term;
        derivative += power * coefficient * term / x;
        second += power * (power - 1.0) * coefficient * term / (x * x);
    }
    (value, derivative, second)
}

/// The binary parameters of the Wong-Sandler mixing rule.
///
/// The matrices hold the interaction of the component `j` in the component
/// `i` at the index `i * count + j`, as in the original.
#[derive(Debug, Clone)]
pub struct WongSandlerData {
    /// The number of components.
    pub count: usize,
    /// The binary interaction parameters, a fraction of the temperature.
    pub kij: Vec<Fraction>,
    /// The parameters of the second order term, a fraction of the temperature.
    pub tauij: Vec<Fraction>,
    /// The exponents of the second order term.
    pub alphaij: Vec<f64>,
}

impl WongSandlerData {
    /// New data for the given number of components.
    ///
    /// Every parameter starts as the zero fraction, as in the original.
    pub fn new(count: usize) -> Self {
        Self {
            count,
            kij: vec![Fraction::default(); count * count],
            tauij: vec![Fraction::default(); count * count],
            alphaij: vec![0.0; count * count],
        }
    }

    /// Set the binary interaction parameter of a pair.
    pub fn set_kij(&mut self, i: usize, j: usize, fraction: Fraction) -> Result<()> {
        let index = self.index(i, j)?;
        self.kij[index] = fraction;
        Ok(())
    }

    /// Set the parameter of the second order term of a pair.
    pub fn set_tauij(&mut self, i: usize, j: usize, fraction: Fraction) -> Result<()> {
        let index = self.index(i, j)?;
        self.tauij[index] = fraction;
        Ok(())
    }

    /// Set the exponent of the second order term of a pair.
    pub fn set_alphaij(&mut self, i: usize, j: usize, value: f64) -> Result<()> {
        let index = self.index(i, j)?;
        self.alphaij[index] = value;
        Ok(())
    }

    /// The index of the interaction of the component `j` in the component `i`.
    fn index(&self, i: usize, j: usize) -> Result<usize> {
        if i >= self.count || j >= self.count {
            return Err(Error::InvalidInput("component index out of range"));
        }
        Ok(i * self.count + j)
    }
}

/// The second mole derivatives of the mixture coefficients of the
/// Wong-Sandler rule.
///
/// The prepared state of the cubic equation of state needs them, because the
/// classic mixing rules derive the second mole derivatives of the fugacity
/// from the binary coefficients of the attraction and of the covolume.
#[derive(Debug, Clone)]
pub struct PreparedWongSandler {
    /// The second mole derivatives of the attraction coefficient.
    pub a_ij: Vec<f64>,
    /// The second mole derivatives of the covolume coefficient.
    pub b_ij: Vec<f64>,
}

/// The mixture coefficients of the Wong-Sandler mixing rule.
///
/// The attraction and the covolume are the extensive coefficients of the
/// original, the partial derivatives are with respect to the number of moles
/// of the components.
#[derive(Debug, Clone)]
pub struct WongSandlerCoefficients {
    /// The attraction coefficient.
    pub a: f64,
    /// The temperature derivative of the attraction coefficient.
    pub a_t: f64,
    /// The second temperature derivative of the attraction coefficient.
    pub a_tt: f64,
    /// The partial derivatives of the attraction coefficient.
    pub a_i: Vec<f64>,
    /// The mixed temperature and mole derivatives of the attraction
    /// coefficient.
    pub a_it: Vec<f64>,
    /// The second mole derivatives of the attraction coefficient.
    pub a_ij: Vec<f64>,
    /// The covolume coefficient.
    pub b: f64,
    /// The temperature derivative of the covolume coefficient.
    pub b_t: f64,
    /// The second temperature derivative of the covolume coefficient.
    pub b_tt: f64,
    /// The partial derivatives of the covolume coefficient.
    pub b_i: Vec<f64>,
    /// The mixed temperature and mole derivatives of the covolume coefficient.
    pub b_it: Vec<f64>,
    /// The second mole derivatives of the covolume coefficient.
    pub b_ij: Vec<f64>,
}

/// The mixture coefficients of the Wong-Sandler mixing rule.
///
/// The attraction parameters of the components are the products of the
/// critical attraction and the alpha function, the covolumes are the
/// covolumes of the components. The excess Gibbs energy at the infinite
/// pressure limit and its derivatives come from the excess Gibbs model of the
/// mixture. The limiting constant `limit_c` is the one of the infinite
/// pressure limit of the sub-equation of state.
#[allow(clippy::too_many_arguments)]
pub fn wong_sandler_mix(
    data: &WongSandlerData,
    temperature: f64,
    composition: &[f64],
    a_pure: &[f64],
    a_pure_t: &[f64],
    a_pure_tt: &[f64],
    b_pure: &[f64],
    limit_c: f64,
    excess: &ExcessGibbsJet,
) -> Result<WongSandlerCoefficients> {
    let count = data.count;
    if composition.len() != count
        || a_pure.len() != count
        || a_pure_t.len() != count
        || a_pure_tt.len() != count
        || b_pure.len() != count
    {
        return Err(Error::InvalidInput("component count mismatch"));
    }
    if limit_c == 0.0 {
        return Err(Error::InvalidInput("the limiting constant vanishes"));
    }
    let gas_constant = crate::R;
    let gas_temperature = gas_constant * temperature;
    let mut total = 0.0;
    for &share in composition {
        total += share;
    }
    if total <= 0.0 {
        return Err(Error::InvalidInput("the composition must be positive"));
    }
    // The reduced covolume of every component and its temperature
    // derivatives.
    let mut reduced = vec![0.0; count];
    let mut reduced_t = vec![0.0; count];
    let mut reduced_tt = vec![0.0; count];
    for i in 0..count {
        reduced[i] = b_pure[i] - a_pure[i] / gas_temperature;
        reduced_t[i] = -a_pure_t[i] / gas_temperature + a_pure[i] / (gas_temperature * temperature);
        reduced_tt[i] = -a_pure_tt[i] / gas_temperature
            + 2.0 * a_pure_t[i] / (gas_temperature * temperature)
            - 2.0 * a_pure[i] / (gas_temperature * temperature * temperature);
    }
    // The pair coefficients and their temperature derivatives.
    let mut pair = vec![0.0; count * count];
    let mut pair_t = vec![0.0; count * count];
    let mut pair_tt = vec![0.0; count * count];
    for i in 0..count {
        for j in 0..count {
            let index = i * count + j;
            let (kij, kij_t, kij_tt) = data.kij[index].evaluate(temperature);
            pair[index] = 0.5 * (reduced[i] + reduced[j]) * (1.0 - kij);
            pair_t[index] = 0.5 * (reduced_t[i] + reduced_t[j]) * (1.0 - kij)
                - 0.5 * (reduced[i] + reduced[j]) * kij_t;
            pair_tt[index] = 0.5 * (reduced_tt[i] + reduced_tt[j]) * (1.0 - kij)
                - (reduced_t[i] + reduced_t[j]) * kij_t
                - 0.5 * (reduced[i] + reduced[j]) * kij_tt;
        }
    }
    // The sums of the appendix A8 of the original.
    let mut q = 0.0;
    let mut q_t = 0.0;
    let mut q_tt = 0.0;
    let mut q_i = vec![0.0; count];
    let mut q_it = vec![0.0; count];
    let mut q_ij = vec![0.0; count * count];
    for i in 0..count {
        for j in 0..count {
            let index = i * count + j;
            let share = composition[i] * composition[j];
            q += share * pair[index];
            q_t += share * pair_t[index];
            q_tt += share * pair_tt[index];
            q_it[i] += 2.0 * pair_t[index] * composition[j];
        }
    }
    for k in 0..count {
        for i in 0..count {
            q_i[k] += composition[i] * (pair[i * count + k] + pair[k * count + i]);
            q_ij[i * count + k] = pair[i * count + k] + pair[k * count + i];
        }
    }
    // The sums of the appendix A9 of the original.
    let mut d = 0.0;
    let mut d_t = 0.0;
    let mut d_tt = 0.0;
    let mut d_i = vec![0.0; count];
    let mut d_it = vec![0.0; count];
    let mut d_ij = vec![0.0; count * count];
    for i in 0..count {
        d += composition[i] * a_pure[i] / b_pure[i];
        d_t += composition[i] * a_pure_t[i] / b_pure[i];
        d_tt += composition[i] * a_pure_tt[i] / b_pure[i];
        d_i[i] = (a_pure[i] / b_pure[i] - excess.dg_dn[i] / limit_c) / gas_temperature;
        d_it[i] = (a_pure_t[i] / b_pure[i] - excess.d2g_dn_dt[i] / limit_c - d_i[i] * gas_constant)
            / gas_temperature;
        for j in 0..count {
            d_ij[i * count + j] = -excess.d2g_dn_dn[i * count + j] / (limit_c * gas_temperature);
        }
    }
    d = (d - excess.gibbs / limit_c) / gas_temperature;
    d_t = (d_t - excess.dg_dt / limit_c - d * gas_constant) / gas_temperature;
    d_tt = (d_tt - excess.d2g_dt2 / limit_c - 2.0 * d_t * gas_constant) / gas_temperature;
    // The covolume of the mixture and its derivatives.
    let denominator = total - d;
    let b = q / denominator;
    let b_t = (q_t + d_t * b) / denominator;
    let b_tt = (q_tt + 2.0 * d_t * b_t + d_tt * b) / denominator;
    let mut b_i = vec![0.0; count];
    let mut b_it = vec![0.0; count];
    for i in 0..count {
        b_i[i] = (q_i[i] + b * (d_i[i] - 1.0)) / denominator;
        b_it[i] = (q_it[i] + b_i[i] * d_t + b_t * (d_i[i] - 1.0) + b * d_it[i]) / denominator;
    }
    let mut b_ij = vec![0.0; count * count];
    for i in 0..count {
        for j in 0..count {
            b_ij[i * count + j] = (q_ij[i * count + j]
                + b_i[i] * (d_i[j] - 1.0)
                + b_i[j] * (d_i[i] - 1.0)
                + b * d_ij[i * count + j])
                / denominator;
        }
    }
    // The attraction of the mixture and its derivatives.
    let a = q * d * gas_temperature / denominator;
    let a_t = gas_constant * (b_t * d * temperature + b * d_t * temperature + b * d);
    let a_tt = gas_constant
        * (b_tt * d * temperature
            + 2.0 * b_t * d_t * temperature
            + 2.0 * b_t * d
            + 2.0 * b * d_t
            + b * d_tt * temperature);
    let mut a_i = vec![0.0; count];
    let mut a_it = vec![0.0; count];
    for k in 0..count {
        a_i[k] = gas_temperature * (b_i[k] * d + b * d_i[k]);
        a_it[k] = gas_constant
            * (b_it[k] * d * temperature
                + b_i[k] * d_t * temperature
                + b_i[k] * d
                + b_t * d_i[k] * temperature
                + b * d_it[k] * temperature
                + b * d_i[k]);
    }
    let mut a_ij = vec![0.0; count * count];
    for i in 0..count {
        for j in 0..count {
            a_ij[i * count + j] = gas_temperature
                * (b_ij[i * count + j] * d
                    + b_i[i] * d_i[j]
                    + b_i[j] * d_i[i]
                    + b * d_ij[i * count + j]);
        }
    }
    Ok(WongSandlerCoefficients {
        a,
        a_t,
        a_tt,
        a_i,
        a_it,
        a_ij,
        b,
        b_t,
        b_tt,
        b_i,
        b_it,
        b_ij,
    })
}
