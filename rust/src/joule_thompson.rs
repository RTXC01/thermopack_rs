//! Joule-Thompson inversion.
//! Ported from `src/joule_thompson_inversion.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! The inversion temperature of a state is the temperature at which the
//! Joule-Thomson coefficient vanishes, the enthalpy of a throttling process
//! does not change. The original solves the condition
//!
//! ```text
//! T (dP/dT)_v + v (dP/dv)_T = 0
//!
//! ```
//!
//! for the temperature at a given molar volume, and it maps the inversion
//! curve with that solution. The native model solves the same condition with
//! a Newton iteration in the logarithm of the temperature, and it uses the
//! same centred difference for the derivative of the condition, as the
//! original does when it extrapolates along the curve.

use crate::cubic::Cubic;
use crate::error::{Error, Result, positive};

/// A point of the Joule-Thomson inversion curve.
#[derive(Debug, Clone, Copy)]
pub struct InversionPoint {
    /// Temperature of the point in kelvin.
    pub temperature: f64,
    /// Molar volume of the point in cubic metre per mole.
    pub volume: f64,
    /// Pressure of the point in pascal.
    pub pressure: f64,
}

/// The relative tolerance of the inversion condition.
const TOLERANCE: f64 = 1.0e-12;
/// The maximum number of the Newton iterations.
const ITERATIONS: usize = 100;

/// The relative step of the centred difference of the condition.
const DIFFERENCE_STEP: f64 = 1.0e-6;

/// The inversion temperature at a given molar volume.
///
/// The volume is in cubic metre per mole, the temperature is the initial
/// guess in kelvin and the composition holds the mole fractions. The
/// condition of the original is solved in the logarithm of the temperature,
/// and the bounds of the flash keep the temperature inside the range of the
/// model.
pub fn inversion_temperature(
    model: &Cubic,
    volume: f64,
    temperature_guess: f64,
    composition: &[f64],
) -> Result<f64> {
    positive(volume, "the inversion volume")?;
    positive(temperature_guess, "the inversion temperature")?;
    if composition.len() != model.components().len() {
        return Err(Error::InvalidInput("composition mismatch"));
    }
    let bounds = model.caloric_options();
    let mut temperature =
        temperature_guess.clamp(bounds.minimum_temperature, bounds.maximum_temperature);
    for _ in 0..ITERATIONS {
        let residual = condition(model, temperature, volume, composition)?;
        let scale = scale(model, temperature, volume, composition)?;
        if residual.abs() <= TOLERANCE * scale {
            return Ok(temperature);
        }
        let above = condition(
            model,
            temperature * (1.0 + DIFFERENCE_STEP),
            volume,
            composition,
        )?;
        let below = condition(
            model,
            temperature * (1.0 - DIFFERENCE_STEP),
            volume,
            composition,
        )?;
        let slope = (above - below) / (2.0 * DIFFERENCE_STEP);
        if slope == 0.0 {
            return Err(Error::NotConverged("the Joule-Thompson inversion"));
        }
        let step = -(residual / slope).exp_m1();
        let mut candidate = temperature * (1.0 + step);
        if !(bounds.minimum_temperature..=bounds.maximum_temperature).contains(&candidate) {
            candidate = 0.5 * (temperature + bounds.minimum_temperature.max(temperature - 10.0));
        }
        if (candidate - temperature).abs() < 1.0e-12 {
            return Ok(temperature);
        }
        temperature = candidate;
    }
    Err(Error::NotConverged("the Joule-Thompson inversion"))
}

/// The inversion point at a given molar volume.
///
/// The point holds the temperature, the molar volume and the pressure of the
/// state.
pub fn inversion_point(
    model: &Cubic,
    volume: f64,
    temperature_guess: f64,
    composition: &[f64],
) -> Result<InversionPoint> {
    let temperature = inversion_temperature(model, volume, temperature_guess, composition)?;
    let state = model.prepare(temperature, composition)?;
    let pressure = state.pressure(volume)?;
    Ok(InversionPoint {
        temperature,
        volume,
        pressure: pressure.value,
    })
}

/// The condition of the inversion at a given temperature and volume.
fn condition(model: &Cubic, temperature: f64, volume: f64, composition: &[f64]) -> Result<f64> {
    let state = model.prepare(temperature, composition)?;
    let pressure = state.pressure(volume)?;
    Ok(temperature * pressure.dt + volume * pressure.dv)
}

/// The scale of the condition, the sum of the magnitudes of its two terms.
fn scale(model: &Cubic, temperature: f64, volume: f64, composition: &[f64]) -> Result<f64> {
    let state = model.prepare(temperature, composition)?;
    let pressure = state.pressure(volume)?;
    Ok((temperature * pressure.dt).abs() + (volume * pressure.dv).abs())
}
