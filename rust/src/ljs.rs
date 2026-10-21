//! Lennard-Jones splined model with the Barker-Henderson perturbation theory.
//! Ported from `src/lj_splined.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! The original splines the excess energy of the Lennard-Jones fluid with
//! the Barker-Henderson perturbation theory. The spline of the excess energy
//! is a set of polynomials in the reduced temperature and the reduced density.
//!
//! The native model provides the parameters and the switches of the model,
//! the Barker-Henderson hard-sphere diameter, the polynomials of the spline
//! with their partial derivatives, the reduced perturbation terms, the reduced
//! residual Helmholtz energy with its temperature and mole number
//! derivatives, the pressure of a state, the reduced residual chemical
//! potential, the caloric state and the reduced interaction potential.

use crate::Phase;
use crate::error::{Error, Result};

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

/// The van der Waals attractive constant of the Lafitte form of the third
/// order term.
const ALPHA_LJS: f64 = 0.5357287350120993;

/// The state of the model at a temperature and a volume.
#[derive(Debug, Clone, Copy)]
pub struct LjsState {
    /// Pressure of the state in pascal.
    pub pressure: f64,
    /// Compressibility factor of the state.
    pub compressibility: f64,
    /// Molar reduced residual Helmholtz energy of the state.
    pub residual_helmholtz: f64,
    /// Molar residual internal energy in joule per mole.
    pub u_residual: f64,
    /// Molar residual enthalpy in joule per mole.
    pub h_residual: f64,
    /// Molar residual entropy in joule per mole and kelvin.
    pub s_residual: f64,
    /// Molar enthalpy of the state in joule per mole.
    pub h_total: f64,
    /// Molar entropy of the state in joule per mole and kelvin.
    pub s_total: f64,
}

/// The Lennard-Jones splined model with the Barker-Henderson perturbation
/// theory.
#[derive(Debug, Clone, Copy)]
pub struct LjsBh {
    sigma: f64,
    eps: f64,
    control: LjsControl,
}

/// The terms of the model that the original lets the user switch on and off.
#[derive(Debug, Clone, Copy)]
pub struct LjsControl {
    /// Whether the hard-sphere reference is included.
    pub enable_hard_sphere: bool,
    /// Whether the first order perturbation term is included.
    pub enable_a1: bool,
    /// Whether the second order perturbation term is included.
    pub enable_a2: bool,
    /// Whether the third order perturbation term is included.
    pub enable_a3: bool,
    /// Whether the correction of the second order term is included.
    pub enable_chi_correction: bool,
    /// Whether the Lafitte form of the third order term is used.
    pub use_lafitte_a3: bool,
}

impl Default for LjsControl {
    fn default() -> Self {
        Self {
            enable_hard_sphere: true,
            enable_a1: true,
            enable_a2: true,
            enable_a3: true,
            enable_chi_correction: true,
            use_lafitte_a3: false,
        }
    }
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
            control: LjsControl::default(),
        }
    }

    /// The model with the segment diameter and the well depth of the user.
    pub fn with_parameters(sigma: f64, eps: f64) -> Self {
        Self {
            sigma,
            eps,
            control: LjsControl::default(),
        }
    }

    /// The terms of the model that the user switches.
    pub fn control(&self) -> LjsControl {
        self.control
    }

    /// Switches the terms of the model, as the original does.
    pub fn set_control(&mut self, control: LjsControl) {
        self.control = control;
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

    /// The logarithm of the fugacity coefficient of the model.
    ///
    /// The residual chemical potential of a single component follows from the
    /// reduced residual Helmholtz energy and from its derivative with respect
    /// to the number of moles, because the reduced density and the packing
    /// fraction follow the number of moles at a fixed temperature and volume.
    pub fn fugacity_coefficient(&self, temperature: f64, volume: f64, moles: f64) -> Result<f64> {
        let (energy, _, mole_derivative) = self.energy_derivatives(temperature, volume, moles)?;
        Ok(energy + moles * mole_derivative)
    }

    /// The packing fraction of the model at a temperature and a density.
    ///
    /// The source uses the hypothetical pure fluid packing fraction of the
    /// Barker-Henderson diameter.
    pub fn packing_fraction(&self, temperature: f64, density: f64) -> Result<f64> {
        crate::error::positive(density, "the LJS density")?;
        let (diameter, _, _) = self.diameter(temperature)?;
        Ok(std::f64::consts::PI / 6.0 * crate::N_A * diameter.powi(3) * density)
    }

    /// The reduced interaction potential of the model and the reduced force.
    ///
    /// The separations are divided by the segment diameter, the potential is
    /// divided by the well depth and the force is the negative derivative of
    /// the potential with respect to the separation. The potential of the
    /// source is the truncated and shifted Lennard-Jones potential with a
    /// smoothing of the cutoff, and it is zero beyond the cutoff.
    pub fn reduced_potential(&self, separations: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let xs = (26.0_f64 / 7.0).powf(1.0 / 6.0);
        let xc = 67.0 / 48.0 * xs;
        let a = -24192.0 / 3211.0 / (xs * xs);
        let b = -387072.0 / 61009.0 / xs.powi(3);
        let mut potential = Vec::with_capacity(separations.len());
        let mut force = Vec::with_capacity(separations.len());
        for &separation in separations {
            if separation < xs {
                let inverse = 1.0 / separation.powi(6);
                potential.push(4.0 * (inverse * inverse - inverse));
                // The force is the negative derivative of the potential with
                // respect to the separation.
                force.push((inverse / separation) * (48.0 * inverse - 24.0));
            } else if separation < xc {
                let difference = separation - xc;
                potential
                    .push(a * difference * difference + b * difference * difference * difference);
                force.push(-(2.0 * a * difference + 3.0 * b * difference * difference));
            } else {
                potential.push(0.0);
                force.push(0.0);
            }
        }
        (potential, force)
    }

    /// The molar volume of the model at a temperature and a pressure.
    ///
    /// The phase selects the liquid branch, which has the smaller molar
    /// volume, and the vapor branch, which has the larger molar volume. The
    /// solve brackets the pressure difference on a logarithmic grid of molar
    /// volumes and refines the bracket with a bisection.
    pub fn volume(&self, temperature: f64, pressure: f64, phase: Phase) -> Result<f64> {
        crate::error::positive(temperature, "the LJS temperature")?;
        crate::error::positive(pressure, "the LJS pressure")?;
        let lower: f64 = 1.0e-8;
        let upper: f64 = 1.0e2;
        // The pressure of the state grows with the volume far from the
        // critical region, so the grid brackets the roots of the difference.
        let mut points = Vec::with_capacity(256);
        let mut volumes = Vec::with_capacity(256);
        for index in 0..256 {
            let volume = lower * (upper / lower).powf(index as f64 / 255.0);
            let value = self.pressure(temperature, volume, 1.0)?.0 - pressure;
            points.push(value);
            volumes.push(volume);
        }
        let mut roots = Vec::new();
        for index in 0..points.len() - 1 {
            if points[index] == 0.0 {
                roots.push(volumes[index]);
            } else if points[index].signum() != points[index + 1].signum() {
                let mut low = volumes[index];
                let mut high = volumes[index + 1];
                for _ in 0..200 {
                    let middle = 0.5 * (low + high);
                    let value = self.pressure(temperature, middle, 1.0)?.0 - pressure;
                    if value.signum() == points[index].signum() {
                        low = middle;
                    } else {
                        high = middle;
                    }
                }
                roots.push(0.5 * (low + high));
            }
        }
        match phase {
            Phase::Liquid => roots.first().copied().ok_or(Error::NotConverged(
                "the liquid volume of the splined Lennard-Jones model",
            )),
            _ => roots.last().copied().ok_or(Error::NotConverged(
                "the vapor volume of the splined Lennard-Jones model",
            )),
        }
    }

    /// The state of the model at a temperature and a volume.
    ///
    /// The ideal part of the state follows the original, which builds the
    /// model on the argon record and fixes the heat capacity of the ideal gas
    /// to a constant of 2.5 times the gas constant.
    pub fn state(&self, temperature: f64, volume: f64, moles: f64) -> Result<LjsState> {
        let (energy, energy_t) =
            self.residual_helmholtz_energy_derivative(temperature, volume, moles)?;
        let (pressure, compressibility) = self.pressure(temperature, volume, moles)?;
        let u_residual = -crate::R * temperature * temperature * energy_t;
        let h_residual = u_residual + crate::R * temperature * (compressibility - 1.0);
        let s_residual =
            -crate::R * (temperature * energy_t + energy) + crate::R * compressibility.ln();
        let ideal = self.ideal()?;
        let properties = ideal.properties(temperature, pressure)?;
        Ok(LjsState {
            pressure,
            compressibility,
            residual_helmholtz: energy,
            u_residual,
            h_residual,
            s_residual,
            h_total: properties.enthalpy + h_residual,
            s_total: properties.entropy + s_residual,
        })
    }

    /// The ideal gas part of the model.
    ///
    /// The original builds the model on the argon record and replaces the heat
    /// capacity of the ideal gas with a constant of 2.5 times the gas
    /// constant, and the reference energies stay those of the record.
    fn ideal(&self) -> Result<crate::ideal::IdealGas> {
        let component = crate::component::component("AR")?;
        let mut gas = crate::ideal::IdealGas::new(component)?;
        let mut coefficients = [0.0; 21];
        coefficients[0] = 2.5;
        gas.set_heat_capacity(8, &coefficients)?;
        Ok(gas)
    }

    /// The pressure of the model at a temperature and a volume.
    ///
    /// The temperature is in kelvin, the volume in cubic metre and the number
    /// of moles is the amount of the single component. The method returns the
    /// pressure and the compressibility factor. The volume derivative of the
    /// reduced residual Helmholtz energy follows the chain rule of the source,
    /// the reduced density and the packing fraction depend on the volume, and
    /// the reduced temperature does not.
    pub fn pressure(&self, temperature: f64, volume: f64, moles: f64) -> Result<(f64, f64)> {
        crate::error::positive(temperature, "the LJS temperature")?;
        crate::error::positive(volume, "the LJS volume")?;
        crate::error::positive(moles, "the LJS number of moles")?;
        let density = moles / volume;
        let (diameter, _, _) = self.diameter(temperature)?;
        let tau = 1.0 - diameter / self.sigma;
        let reduced_density = crate::N_A * self.sigma.powi(3) * density;
        let reduced_density_v = -crate::N_A * self.sigma.powi(3) * moles / (volume * volume);
        let eta = std::f64::consts::PI / 6.0 * crate::N_A * diameter.powi(3) * density;
        let eta_v = -eta / volume;
        // The volume derivative of the molar reduced residual energy.
        let mut derivative = 0.0;
        if self.control.enable_hard_sphere {
            derivative += (4.0 - 2.0 * eta) / (1.0 - eta).powi(3) * eta_v;
        }
        if self.control.enable_a1 {
            let (_, value_r, ..) = polynomial_with_derivatives(0, tau, reduced_density, 1.0);
            derivative +=
                2.0 * std::f64::consts::PI * self.eps * value_r * reduced_density_v / temperature;
        }
        if self.control.enable_a2 {
            let (value, value_r, ..) = polynomial_with_derivatives(1, tau, reduced_density, 1.0);
            let (correction, correction_r, ..) =
                polynomial_with_derivatives(3, tau, reduced_density, 1.0);
            let mut second = value;
            let mut second_v = value_r * reduced_density_v;
            if self.control.enable_chi_correction {
                second *= 1.0 + correction;
                second_v = second_v * (1.0 + correction) + value * correction_r * reduced_density_v;
            }
            let compressibility = khs(eta);
            let compressibility_v = khs_eta(eta) * eta_v;
            second_v = second_v * compressibility + second * compressibility_v;
            derivative += -std::f64::consts::PI * self.eps * self.eps * second_v
                / (temperature * temperature);
        }
        if self.control.enable_a3 {
            if self.control.use_lafitte_a3 {
                // The Lafitte form of the third order term follows the
                // dispersion packing fraction and the polynomials of the
                // van der Waals attractive constant of the source.
                let zeta = std::f64::consts::PI / 6.0 * crate::N_A * self.sigma.powi(3) * density;
                let f_alpha = crate::saft_vrmie::alpha_coefficients(ALPHA_LJS);
                let exponent = f_alpha[4] * zeta + f_alpha[5] * zeta * zeta;
                let factor = exponent.exp();
                let slope = -self.eps.powi(3)
                    * (f_alpha[3] * factor
                        + zeta * (f_alpha[4] + 2.0 * f_alpha[5] * zeta) * factor);
                derivative += slope * (-zeta / volume) / temperature.powi(3);
            } else {
                let (_, value_r, ..) = polynomial_with_derivatives(2, tau, reduced_density, 1.0);
                derivative += self.eps.powi(3) * value_r * reduced_density_v / temperature.powi(3);
            }
        }
        // The pressure of the ideal gas and the pressure of the residual
        // energy, which is the negative volume derivative of the residual
        // energy.
        let ideal = crate::R * temperature * moles / volume;
        let residual = -crate::R * temperature * moles * derivative;
        let pressure = ideal + residual;
        let compressibility = pressure / (density * crate::R * temperature);
        Ok((pressure, compressibility))
    }

    /// The molar reduced residual Helmholtz energy of the model and its
    /// temperature derivative.
    ///
    /// The temperature is in kelvin, the volume in cubic metre and the number
    /// of moles is the amount of the single component. The energy is the
    /// reduced residual energy of the original, the hard-sphere reference and
    /// the three perturbation terms. The derivative follows the chain rule of
    /// the original: the reduced temperature depends on the temperature through
    /// the Barker-Henderson diameter, and the reduced density does not.
    pub fn residual_helmholtz_energy_derivative(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<(f64, f64)> {
        let (energy, derivative, _) = self.energy_derivatives(temperature, volume, moles)?;
        Ok((energy, derivative))
    }

    /// The molar reduced residual Helmholtz energy of the model with its
    /// temperature and mole number derivatives.
    fn energy_derivatives(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<(f64, f64, f64)> {
        crate::error::positive(temperature, "the LJS temperature")?;
        crate::error::positive(volume, "the LJS volume")?;
        crate::error::positive(moles, "the LJS number of moles")?;
        let density = moles / volume;
        let (diameter, diameter_t, _) = self.diameter(temperature)?;
        let tau = 1.0 - diameter / self.sigma;
        let tau_t = -diameter_t / self.sigma;
        let reduced_density = crate::N_A * self.sigma.powi(3) * density;
        let reduced_density_n = crate::N_A * self.sigma.powi(3) / volume;
        let eta = std::f64::consts::PI / 6.0 * crate::N_A * diameter.powi(3) * density;
        let eta_t = std::f64::consts::PI / 6.0
            * crate::N_A
            * 3.0
            * diameter
            * diameter
            * diameter_t
            * density;
        let eta_n = std::f64::consts::PI / 6.0 * crate::N_A * diameter.powi(3) / volume;
        let mut energy = 0.0;
        let mut derivative = 0.0;
        let mut mole_derivative = 0.0;
        if self.control.enable_hard_sphere {
            let slope = (4.0 - 2.0 * eta) / (1.0 - eta).powi(3);
            energy += (4.0 * eta - 3.0 * eta * eta) / (1.0 - eta).powi(2);
            derivative += slope * eta_t;
            mole_derivative += slope * eta_n;
        }
        if self.control.enable_a1 {
            let (value, value_r, value_x, ..) =
                polynomial_with_derivatives(0, tau, reduced_density, 1.0);
            let first = 2.0 * std::f64::consts::PI * self.eps * value;
            let first_t = 2.0 * std::f64::consts::PI * self.eps * value_x * tau_t;
            energy += first / temperature;
            derivative += first_t / temperature - first / (temperature * temperature);
            mole_derivative +=
                2.0 * std::f64::consts::PI * self.eps * value_r * reduced_density_n / temperature;
        }
        if self.control.enable_a2 {
            let (value, value_r, value_x, ..) =
                polynomial_with_derivatives(1, tau, reduced_density, 1.0);
            let (correction, correction_r, correction_x, ..) =
                polynomial_with_derivatives(3, tau, reduced_density, 1.0);
            let mut second = -std::f64::consts::PI * self.eps * self.eps * value;
            let mut second_t = -std::f64::consts::PI * self.eps * self.eps * value_x * tau_t;
            let mut second_r = -std::f64::consts::PI * self.eps * self.eps * value_r;
            if self.control.enable_chi_correction {
                let scaled = 1.0 + correction;
                let scaled_t = correction_x * tau_t;
                let scaled_r = correction_r;
                second *= scaled;
                second_t = second_t * scaled
                    + (-std::f64::consts::PI * self.eps * self.eps * value) * scaled_t;
                second_r = second_r * scaled
                    + (-std::f64::consts::PI * self.eps * self.eps * value) * scaled_r;
            }
            let compressibility = khs(eta);
            let compressibility_t = khs_eta(eta) * eta_t;
            let unscaled = second;
            second *= compressibility;
            second_t = second_t * compressibility + unscaled * compressibility_t;
            energy += second / (temperature * temperature);
            derivative +=
                second_t / (temperature * temperature) - 2.0 * second / temperature.powi(3);
            // The compressibility factor follows the packing fraction, and the
            // packing fraction follows the number of moles.
            mole_derivative += second_r * reduced_density_n * compressibility
                / (temperature * temperature)
                + unscaled * khs_eta(eta) * eta_n / (temperature * temperature);
        }
        if self.control.enable_a3 {
            let (third, third_t, third_n) = if self.control.use_lafitte_a3 {
                // The Lafitte form of the third order term follows the
                // dispersion packing fraction and the polynomials of the
                // van der Waals attractive constant of the source.
                let zeta = std::f64::consts::PI / 6.0 * crate::N_A * self.sigma.powi(3) * density;
                let zeta_n = std::f64::consts::PI / 6.0 * crate::N_A * self.sigma.powi(3) / volume;
                let f_alpha = crate::saft_vrmie::alpha_coefficients(ALPHA_LJS);
                let exponent = f_alpha[4] * zeta + f_alpha[5] * zeta * zeta;
                let factor = exponent.exp();
                let value = -self.eps.powi(3) * f_alpha[3] * zeta * factor;
                let slope = -self.eps.powi(3)
                    * (f_alpha[3] * factor
                        + zeta * (f_alpha[4] + 2.0 * f_alpha[5] * zeta) * factor);
                (value, 0.0, slope * zeta_n)
            } else {
                let (value, value_r, value_x, ..) =
                    polynomial_with_derivatives(2, tau, reduced_density, 1.0);
                let third = self.eps.powi(3) * value;
                let third_t = self.eps.powi(3) * value_x * tau_t;
                let third_n = self.eps.powi(3) * value_r * reduced_density_n;
                (third, third_t, third_n)
            };
            energy += third / temperature.powi(3);
            derivative += third_t / temperature.powi(3) - 3.0 * third / temperature.powi(4);
            mole_derivative += third_n / temperature.powi(3);
        }
        Ok((energy, derivative, mole_derivative))
    }

    /// The molar reduced residual Helmholtz energy of the model.
    ///
    /// The temperature is in kelvin, the volume in cubic metre and the number
    /// of moles is the amount of the single component. The energy is the
    /// reduced residual energy of the original, the hard-sphere reference and
    /// the three perturbation terms.
    pub fn residual_helmholtz_energy(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<f64> {
        Ok(self
            .residual_helmholtz_energy_derivative(temperature, volume, moles)?
            .0)
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

    /// One polynomial of the spline of the excess energy and its partial
    /// derivatives.
    ///
    /// The derivatives are the ones of the source: the first and the second
    /// derivative with respect to the reduced temperature, the first and the
    /// second derivative with respect to the reduced density and the mixed
    /// derivative.
    pub fn spline_polynomial(
        &self,
        index: usize,
        reduced_temperature: f64,
        reduced_density: f64,
    ) -> Result<(f64, f64, f64, f64, f64, f64)> {
        let temperature = reduced_temperature * self.eps;
        let (diameter, _, _) = self.diameter(temperature)?;
        let x0 = self.sigma / diameter;
        let tau = 1.0 - 1.0 / x0;
        Ok(polynomial_with_derivatives(
            index,
            tau,
            reduced_density,
            1.0,
        ))
    }
}

/// One polynomial of the spline of the excess energy.
///
/// The source stores the coefficients in a four by fifteen array, and the
/// polynomial of the term `index` uses the coefficients of the columns zero
/// to four, five to nine and ten to fourteen of that row.
fn polynomial(index: usize, tau: f64, rho: f64) -> f64 {
    polynomial_with_derivatives(index, tau, rho, 1.0).0
}

/// One polynomial of the spline of the excess energy with its partial
/// derivatives and the prefactor of the source.
fn polynomial_with_derivatives(
    index: usize,
    tau: f64,
    rho: f64,
    prefactor: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    // The powers of the reduced density and their first two derivatives, the
    // fifth power first.
    let mut powers = [0.0; 5];
    powers[4] = rho;
    for position in (0..4).rev() {
        powers[position] = powers[position + 1] * rho;
    }
    let mut first = [0.0; 5];
    first[4] = 1.0;
    first[3] = 2.0 * rho;
    first[2] = 3.0 * powers[3];
    first[1] = 4.0 * powers[2];
    first[0] = 5.0 * powers[1];
    let mut second = [0.0; 5];
    second[4] = 0.0;
    second[3] = 2.0;
    second[2] = 6.0 * rho;
    second[1] = 12.0 * powers[3];
    second[0] = 20.0 * powers[2];
    // The three polynomials in the reduced temperature and their first two
    // derivatives.
    let g = [
        [1.0, 0.0, 0.0],
        [tau, 1.0, 0.0],
        [tau * tau, 2.0 * tau, 2.0],
    ];
    let mut value = 0.0;
    let mut value_x = 0.0;
    let mut value_r = 0.0;
    let mut value_xx = 0.0;
    let mut value_rr = 0.0;
    let mut value_rx = 0.0;
    for (position, block) in [0usize, 5, 10].into_iter().enumerate() {
        let mut sum = 0.0;
        let mut sum_first = 0.0;
        let mut sum_second = 0.0;
        for column in 0..5 {
            let coefficient = POLYNOMIALS[(block + column) * 4 + index];
            sum += coefficient * powers[column];
            sum_first += coefficient * first[column];
            sum_second += coefficient * second[column];
        }
        value += sum * g[position][0];
        value_x += sum * g[position][1];
        value_r += sum_first * g[position][0];
        value_xx += sum * g[position][2];
        value_rr += sum_second * g[position][0];
        value_rx += sum_first * g[position][1];
    }
    (
        prefactor * value,
        prefactor * value_r,
        prefactor * value_x,
        prefactor * value_rr,
        prefactor * value_xx,
        prefactor * value_rx,
    )
}

/// The isothermal compressibility factor of the hard-sphere reference, the
/// Carnahan-Starling expression of the original.
fn khs(eta: f64) -> f64 {
    let denominator =
        1.0 + 4.0 * eta + 4.0 * eta * eta - 4.0 * eta * eta * eta + eta * eta * eta * eta;
    let one_minus = 1.0 - eta;
    one_minus.powi(4) / denominator
}

/// The derivative of the isothermal compressibility factor of the hard-sphere
/// reference with respect to the packing fraction.
fn khs_eta(eta: f64) -> f64 {
    let eta2 = eta * eta;
    let denominator = 1.0 + 4.0 * eta + 4.0 * eta2 - 4.0 * eta2 * eta + eta2 * eta2;
    4.0 * (eta2 - 5.0 * eta - 2.0) * (1.0 - eta).powi(3) / (denominator * denominator)
}
/// The u-fraction of the UV theory and its partial derivatives.
///
/// The well depth is divided by the Boltzmann constant, the temperature is in
/// kelvin and the density is the reduced density. The method returns the
/// fraction, its derivative with respect to the reduced density, its
/// derivative with respect to the temperature, the second derivative with
/// respect to the reduced density and the mixed derivative, in the order of the
/// source. The fraction of the source does not depend on the temperature, so
/// the temperature derivatives vanish.
pub fn uv_phi(_eps_divk: f64, _temperature: f64, rho: f64) -> (f64, f64, f64, f64, f64, f64) {
    const C: [f64; 5] = [0.11072527, 2.5809656, 1.2333385, 2.3534453, 4.1241290];
    let x = C[0] * rho + C[1] * rho.powf(C[2]) + C[3] * rho.powf(C[4]);
    let x_r = C[0] + C[1] * C[2] * rho.powf(C[2] - 1.0) + C[3] * C[4] * rho.powf(C[4] - 1.0);
    let x_rr = C[1] * C[2] * (C[2] - 1.0) * rho.powf(C[2] - 2.0)
        + C[3] * C[4] * (C[4] - 1.0) * rho.powf(C[4] - 2.0);
    let e2x = (2.0 * x).exp();
    let tanh = (e2x - 1.0) / (e2x + 1.0);
    let denominator = (e2x + 1.0).powi(2);
    let tanh_r = 4.0 * e2x * x_r / denominator;
    let tanh_rr = 2.0 * tanh_r * x_r + 4.0 * e2x * x_rr / denominator
        - 16.0 * e2x * e2x * x_r * x_r / (denominator * (e2x + 1.0));
    (tanh, tanh_r, 0.0, tanh_rr, 0.0, 0.0)
}

/// The temperature dependent u-fraction of the UV theory and its partial
/// derivatives.
///
/// The method follows the source and reports the same quantities as
/// `uv_phi`.
pub fn uv_phi_temp(eps_divk: f64, temperature: f64, rho: f64) -> (f64, f64, f64, f64, f64, f64) {
    const C: [f64; 7] = [
        1.01236, 0.99164, 2.63226, 1.21219, 2.25005, 3.91411, 0.436593,
    ];
    let beta = eps_divk / temperature;
    let beta_t = -eps_divk / (temperature * temperature);
    let beta_tt = 2.0 * eps_divk / temperature.powi(3);
    let sqrt_c7 = C[6].sqrt();
    let denominator_sq = 1.0 + C[6] * beta * beta;
    let denominator = denominator_sq.sqrt();
    let sigma_beta = sqrt_c7 * beta / denominator;
    let sigma_beta_b =
        sqrt_c7 / denominator - sqrt_c7 * C[6] * beta * beta / (denominator * denominator_sq);
    let sigma_beta_bb = -3.0 * sqrt_c7 * C[6] * beta / (denominator * denominator_sq)
        + 3.0 * sqrt_c7 * C[6] * C[6] * beta.powi(3)
            / (denominator * denominator_sq * denominator_sq);
    let prefactor = C[0] + sigma_beta * (C[1] - C[0]);
    let x = C[2] * rho.powf(C[3]) + C[4] * rho.powf(C[5]);
    let x_r = C[2] * C[3] * rho.powf(C[3] - 1.0) + C[4] * C[5] * rho.powf(C[5] - 1.0);
    let x_rr = C[2] * C[3] * (C[3] - 1.0) * rho.powf(C[3] - 2.0)
        + C[4] * C[5] * (C[5] - 1.0) * rho.powf(C[5] - 2.0);
    let e2x = (2.0 * x).exp();
    let tanh = (e2x - 1.0) / (e2x + 1.0);
    let denominator = (e2x + 1.0).powi(2);
    let tanh_r = 4.0 * e2x * x_r / denominator;
    let tanh_rr = 2.0 * tanh_r * x_r + 4.0 * e2x * x_rr / denominator
        - 16.0 * e2x * e2x * x_r * x_r / (denominator * (e2x + 1.0));
    let value = prefactor * tanh;
    let value_b = sigma_beta_b * (C[1] - C[0]) * tanh;
    let value_bb = sigma_beta_bb * (C[1] - C[0]) * tanh;
    let value_r = prefactor * tanh_r;
    let value_rr = prefactor * tanh_rr;
    let value_rb = sigma_beta_b * (C[1] - C[0]) * tanh_r;
    (
        value,
        value_r,
        value_b * beta_t,
        value_rr,
        value_bb * beta_t * beta_t + value_b * beta_tt,
        value_rb * beta_t,
    )
}
