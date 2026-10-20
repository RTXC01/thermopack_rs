//! Lennard-Jones splined model with the Barker-Henderson perturbation theory.
//! Ported from `src/lj_splined.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! The original splines the excess energy of the Lennard-Jones fluid with
//! the Barker-Henderson perturbation theory. The spline of the excess energy
//! is a set of polynomials in the reduced temperature and the reduced density,
//! and the original exposes the reduced perturbation terms and the
//! Barker-Henderson diameter through the Python interface of the fork.
//!
//! The native model provides the same two quantities, the reduced
//! perturbation terms and the Barker-Henderson diameter, and it prepares the
//! parameters of the model.

use crate::error::Result;

/// The polynomial coefficients of the source, in the order of the source
/// array. The source stores them in a four by fifteen array, and the first
/// polynomial of a row holds its first five coefficients, the second
/// polynomial the next five and the third the last five.
const POLYNOMIALS: [f64; 60] = [
    0.039813685,
    -0.106920241,
    1.291256810,
    0.000000000,
    0.469671161,
    -0.295578203,
    -2.763380451,
    21.299989640,
    -0.346120315,
    0.343721691,
    1.416116081,
    -11.284586250,
    -0.344778684,
    0.248954221,
    0.331595545,
    -0.258500011,
    -0.535278660,
    0.358663104,
    -0.248745451,
    0.996821621,
    0.000000000,
    0.000000000,
    0.007804235,
    0.000000000,
    -7.861555167,
    7.084518499,
    2.799119819,
    -540.982471100,
    10.084157010,
    -9.814532406,
    -4.455465501,
    599.182974200,
    0.518650855,
    0.172122707,
    1.855450051,
    -116.677210100,
    0.000000000,
    0.000000000,
    -0.350181778,
    3.804716369,
    0.000000000,
    0.000000000,
    0.000000000,
    0.000000000,
    32.718743790,
    -28.998566930,
    0.000000000,
    3131.407636000,
    -51.573804390,
    47.108198520,
    0.000000000,
    -4126.665342000,
    6.550964581,
    -7.263318791,
    0.000000000,
    1013.877746000,
    0.000000000,
    0.000000000,
    0.000000000,
    -70.079287230,
];

/// The segment diameter of the argon record of the source database in metre.
const ARGON_SIGMA: f64 = 3.42e-10;

/// The well depth of the argon record of the source database in kelvin.
const ARGON_EPS: f64 = 124.0;

/// The repulsive exponent of the model.
const LAMBDA_A: f64 = 6.0;

/// The attractive exponent of the model.
const LAMBDA_R: f64 = 12.0;

/// The prefactor of the potential of the model.
const C_FAC: f64 = 4.0;

/// The Lennard-Jones splined model with the Barker-Henderson perturbation
/// theory.
#[derive(Debug, Clone, Copy)]
pub struct LjsBh {
    sigma: f64,
    eps: f64,
}

impl Default for LjsBh {
    fn default() -> Self {
        Self::new()
    }
}

impl LjsBh {
    /// The model of the source database, the argon record.
    pub fn new() -> Self {
        Self {
            sigma: ARGON_SIGMA,
            eps: ARGON_EPS,
        }
    }

    /// The model with the segment diameter and the well depth of the user.
    pub fn with_parameters(sigma: f64, eps: f64) -> Self {
        Self { sigma, eps }
    }

    /// The segment diameter of the model in metre.
    pub fn sigma(&self) -> f64 {
        self.sigma
    }

    /// The well depth of the model divided by the Boltzmann constant in
    /// kelvin.
    pub fn eps(&self) -> f64 {
        self.eps
    }

    /// The Barker-Henderson hard-sphere diameter of the model and its first
    /// and second temperature derivatives in metre.
    pub fn diameter(&self, temperature: f64) -> Result<(f64, f64, f64)> {
        crate::saft_vrmie::calc_hardsphere_diameter_pure(
            temperature,
            self.sigma,
            self.eps,
            LAMBDA_A,
            LAMBDA_R,
            C_FAC,
            1.0,
            0,
        )
    }

    /// The Barker-Henderson hard-sphere diameter divided by the segment
    /// diameter at a reduced temperature.
    pub fn diameter_div_sigma(&self, reduced_temperature: f64) -> Result<f64> {
        let temperature = reduced_temperature * self.eps;
        let (diameter, _, _) = self.diameter(temperature)?;
        Ok(diameter / self.sigma)
    }

    /// The reduced perturbation terms of the model.
    ///
    /// The reduced temperature and the reduced density are dimensionless, and
    /// the terms are the first, the second and the third order term of the
    /// perturbation theory divided by the well depth, the well depth squared
    /// and the well depth cubed.
    pub fn reduced_perturbation_terms(
        &self,
        reduced_temperature: f64,
        reduced_density: f64,
    ) -> Result<(f64, f64, f64)> {
        let temperature = reduced_temperature * self.eps;
        let (diameter, _, _) = self.diameter(temperature)?;
        let x0 = self.sigma / diameter;
        let tau = 1.0 - 1.0 / x0;
        let first = 2.0 * std::f64::consts::PI * polynomial(0, tau, reduced_density);
        let second = -std::f64::consts::PI * polynomial(1, tau, reduced_density);
        let correction = polynomial(3, tau, reduced_density);
        let eta = std::f64::consts::PI * reduced_density / (6.0 * x0 * x0 * x0);
        let compressibility = khs(eta);
        let second = second * (1.0 + correction) * compressibility;
        let third = polynomial(2, tau, reduced_density);
        Ok((first, second, third))
    }
}

/// One polynomial of the spline of the excess energy.
///
/// The source stores the coefficients in a four by fifteen array, and the
/// polynomial of the term `index` uses the coefficients of the columns zero
/// to four, five to nine and ten to fourteen of that row.
fn polynomial(index: usize, tau: f64, rho: f64) -> f64 {
    // The powers of the reduced density, the fifth power first.
    let mut powers = [0.0; 5];
    powers[4] = rho;
    for position in (0..4).rev() {
        powers[position] = powers[position + 1] * rho;
    }
    // The three polynomials in the reduced temperature.
    let g = [1.0, tau, tau * tau];
    let mut value = 0.0;
    for (position, block) in [0usize, 5, 10].into_iter().enumerate() {
        let mut sum = 0.0;
        for column in 0..5 {
            sum += POLYNOMIALS[(block + column) * 4 + index] * powers[column];
        }
        value += sum * g[position];
    }
    value
}

/// The isothermal compressibility factor of the hard-sphere reference, the
/// Carnahan-Starling expression of the original.
fn khs(eta: f64) -> f64 {
    let denominator =
        1.0 + 4.0 * eta + 4.0 * eta * eta - 4.0 * eta * eta * eta + eta * eta * eta * eta;
    let one_minus = 1.0 - eta;
    one_minus.powi(4) / denominator
}
