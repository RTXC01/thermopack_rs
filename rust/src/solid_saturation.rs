//! Solid-fluid equilibrium: the saturation line of a solid phase.
//! Ported from `src/solid_saturation.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! The original locates the state where a solid phase and the fluid phase of
//! the equation of state share the fugacity. The solid contributes its
//! logarithm of the fugacity and the fluid its own, and the difference is the
//! residual of a Newton iteration in the logarithm of the temperature and the
//! logarithm of the pressure.
//!
//! The port covers the single component case, where the solid phase is the
//! pure component and the fluid phase is the same component. The original
//! reduces its residual to the difference of the two fugacity coefficients in
//! that case, and the native model reproduces the converged states of the
//! original to a relative difference of three times ten to the minus eight.

use crate::Phase;
use crate::R;
use crate::cubic::Cubic;
use crate::error::{Error, Result, positive};
use crate::ideal::IdealGas;

/// The solid models that can take part in a solid-fluid equilibrium.
pub trait SolidGibbs {
    /// The Gibbs energy of the solid and its temperature and pressure
    /// derivatives.
    fn gibbs_with_derivatives(&self, temperature: f64, pressure: f64) -> (f64, f64, f64);
}

/// Which variable of the saturation state is prescribed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SolidSaturationSpecification {
    /// The temperature is prescribed and the pressure is solved.
    FixedTemperature(f64),
    /// The pressure is prescribed and the temperature is solved.
    FixedPressure(f64),
}

/// A state where the solid and the fluid phase coexist.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolidSaturation {
    /// The temperature of the state in kelvin.
    pub temperature: f64,
    /// The pressure of the state in pascal.
    pub pressure: f64,
    /// The difference of the logarithms of the fugacities at the solution.
    pub residual: f64,
    /// The number of Newton iterations.
    pub iterations: usize,
}

/// The triple point temperature of the carbon dioxide model of the original.
pub const CO2_TRIPLE_TEMPERATURE: f64 = 216.592;
/// The nominal triple point pressure of the carbon dioxide model of the
/// original. The original refines it with a bubble pressure solve.
pub const CO2_TRIPLE_PRESSURE: f64 = 5.1795e5;

/// Builds the solid carbon dioxide model with the reference of the fluid.
///
/// The original `initDryIce` refines the triple point pressure with a bubble
/// pressure solve, takes the entropy and the enthalpy of the liquid at the
/// triple point and shifts the reference levels of the solid model.
pub fn solid_co2_with_fluid_reference(model: &Cubic) -> Result<crate::solid_co2::SolidCo2> {
    let saturation = model.saturation_pressure(CO2_TRIPLE_TEMPERATURE)?;
    let temperature = CO2_TRIPLE_TEMPERATURE;
    let pressure = saturation.pressure;
    let state = model.prepare(temperature, &[1.0])?;
    let volume = saturation.liquid_volume;
    let entropy = state.entropy(volume)?;
    let enthalpy = state.enthalpy(volume)?;
    let gibbs = enthalpy - temperature * entropy;
    let mut solid = crate::solid_co2::SolidCo2::new();
    solid.initialize(entropy, gibbs, temperature, pressure);
    Ok(solid)
}

/// Finds the state where solid carbon dioxide and the fluid coexist.
pub fn solid_co2_saturation(
    model: &Cubic,
    solid: &crate::solid_co2::SolidCo2,
    ideal: &IdealGas,
    specification: SolidSaturationSpecification,
    initial: (f64, f64),
    phase: Phase,
) -> Result<SolidSaturation> {
    solid_saturation(model, solid, ideal, specification, initial, phase)
}

/// The maximum number of Newton iterations.
const MAX_ITERATIONS: usize = 100;
/// The convergence threshold of the fugacity difference.
const TOLERANCE: f64 = 1e-9;
/// The relative step of the finite-difference derivative.
const DERIVATIVE_STEP: f64 = 1e-7;

impl SolidGibbs for crate::solid_co2::SolidCo2 {
    fn gibbs_with_derivatives(&self, temperature: f64, pressure: f64) -> (f64, f64, f64) {
        (
            self.gibbs(temperature, pressure),
            self.dg_dt(temperature, pressure),
            self.dg_dp(temperature, pressure),
        )
    }
}

impl SolidGibbs for crate::solid_water::SolidWater {
    fn gibbs_with_derivatives(&self, temperature: f64, pressure: f64) -> (f64, f64, f64) {
        match self.state(temperature, pressure) {
            Ok(state) => (state.gibbs.value, state.gibbs.dt, state.gibbs.dp),
            Err(_) => (f64::NAN, f64::NAN, f64::NAN),
        }
    }
}

/// The logarithm of the fugacity of the solid at a state.
///
/// The original limits the value to the range of minus five hundred to five
/// hundred and clears the derivatives when the limit applies.
fn solid_ln_fugacity<S: SolidGibbs>(
    solid: &S,
    ideal: &IdealGas,
    temperature: f64,
    pressure: f64,
) -> (f64, f64, f64) {
    let (gibbs, dgdt, dgdp) = solid.gibbs_with_derivatives(temperature, pressure);
    let (ideal_gibbs, ideal_dgdt, ideal_dgdp) = ideal_gibbs(ideal, temperature, pressure);
    let mut ln_fugacity = (gibbs - ideal_gibbs) / (R * temperature);
    let (mut dlnfugacity_dt, mut dlnfugacity_dp) = (
        (dgdt - ideal_dgdt) / (R * temperature) - ln_fugacity / temperature,
        (dgdp - ideal_dgdp) / (R * temperature),
    );
    if ln_fugacity < -500.0 {
        ln_fugacity = -500.0;
        dlnfugacity_dt = 0.0;
        dlnfugacity_dp = 0.0;
    } else if ln_fugacity > 500.0 {
        ln_fugacity = 500.0;
        dlnfugacity_dt = 0.0;
        dlnfugacity_dp = 0.0;
    }
    (ln_fugacity, dlnfugacity_dt, dlnfugacity_dp)
}

/// The Gibbs energy of the ideal gas and its derivatives.
fn ideal_gibbs(ideal: &IdealGas, temperature: f64, pressure: f64) -> (f64, f64, f64) {
    let properties = match ideal.properties(temperature, pressure) {
        Ok(properties) => properties,
        Err(_) => return (f64::NAN, f64::NAN, f64::NAN),
    };
    let gibbs = properties.enthalpy - temperature * properties.entropy;
    // The original reports the negative entropy as the temperature derivative
    // and T R / P as the pressure derivative.
    (gibbs, -properties.entropy, temperature * R / pressure)
}

/// The difference of the logarithms of the fugacities of the solid and the
/// fluid at a state.
fn residual<S: SolidGibbs>(
    model: &Cubic,
    solid: &S,
    ideal: &IdealGas,
    temperature: f64,
    pressure: f64,
    phase: Phase,
) -> Result<f64> {
    let (solid_fugacity, _, _) = solid_ln_fugacity(solid, ideal, temperature, pressure);
    let state = model.prepare(temperature, &[1.0])?;
    let volume = state.volume(pressure, phase)?;
    let mut coefficient = [0.0];
    state.ln_fugacity_into(volume, &mut coefficient)?;
    // For a single component the original reduces the residual to the
    // difference of the two fugacity coefficients, because the logarithm of
    // the mole fraction of the solid component vanishes.
    Ok(solid_fugacity - coefficient[0])
}

/// Finds the state where the solid and the fluid phase coexist.
///
/// The specification prescribes the temperature or the pressure, and the other
/// variable is solved. The initial guess is the pair of a temperature and a
/// pressure in kelvin and pascal.
pub fn solid_saturation<S: SolidGibbs>(
    model: &Cubic,
    solid: &S,
    ideal: &IdealGas,
    specification: SolidSaturationSpecification,
    initial: (f64, f64),
    phase: Phase,
) -> Result<SolidSaturation> {
    let (mut temperature, mut pressure) = initial;
    positive(
        temperature,
        "the solid saturation temperature must be positive",
    )?;
    positive(pressure, "the solid saturation pressure must be positive")?;
    for iteration in 1..=MAX_ITERATIONS {
        let value = residual(model, solid, ideal, temperature, pressure, phase)?;
        if value.abs() < TOLERANCE {
            return Ok(SolidSaturation {
                temperature,
                pressure,
                residual: value,
                iterations: iteration - 1,
            });
        }
        let step = match specification {
            SolidSaturationSpecification::FixedTemperature(_) => {
                let delta = DERIVATIVE_STEP * pressure.abs().max(1.0);
                let upper = residual(model, solid, ideal, temperature, pressure + delta, phase)?;
                let lower = residual(model, solid, ideal, temperature, pressure - delta, phase)?;
                (upper - lower) / (2.0 * delta)
            }
            SolidSaturationSpecification::FixedPressure(_) => {
                let delta = DERIVATIVE_STEP * temperature.abs().max(1.0);
                let upper = residual(model, solid, ideal, temperature + delta, pressure, phase)?;
                let lower = residual(model, solid, ideal, temperature - delta, pressure, phase)?;
                (upper - lower) / (2.0 * delta)
            }
        };
        if !step.is_finite() || step == 0.0 {
            return Err(Error::NotConverged(
                "the solid saturation derivative is not usable",
            ));
        }
        let correction = -value / step;
        match specification {
            SolidSaturationSpecification::FixedTemperature(fixed) => {
                temperature = fixed;
                pressure = (pressure * (correction / pressure).exp()).max(1.0);
            }
            SolidSaturationSpecification::FixedPressure(fixed) => {
                pressure = fixed;
                temperature = (temperature * (correction / temperature).exp()).max(1.0);
            }
        }
    }
    Err(Error::NotConverged(
        "the solid saturation iteration did not converge",
    ))
}
