//! Critical points of a pure fluid.
//! Ported from `src/critical.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! The critical point of a pure fluid is the state where the isotherm has an
//! inflection, so the first and the second volume derivative of the pressure
//! vanish. The original solves the two conditions with the stability matrix
//! of the mixture, and for a single component the two conditions of the
//! matrix are the volume derivatives of the pressure.

use crate::R;
use crate::cubic::Cubic;
use crate::error::{Error, Result};

/// A critical point of a pure fluid.
#[derive(Debug, Clone, Copy)]
pub struct CriticalPoint {
    /// Temperature of the critical point in kelvin.
    pub temperature: f64,
    /// Molar volume of the critical point in cubic metre per mole.
    pub volume: f64,
    /// Pressure of the critical point in pascal.
    pub pressure: f64,
}

/// The maximum number of the Newton iterations.
const ITERATIONS: usize = 100;

/// The relative step of the finite difference of the Jacobian.
const STEP: f64 = 1.0e-6;

/// The relative tolerance of the two conditions.
const TOLERANCE: f64 = 1.0e-10;

/// The critical point of a pure fluid.
///
/// The composition holds the mole fractions, which must describe a single
/// component. The initial guess is the pseudo critical state of the model,
/// and the two conditions of the critical point are solved with a Newton
/// iteration in the temperature and the molar volume.
pub fn critical_point(model: &Cubic, composition: &[f64]) -> Result<CriticalPoint> {
    if composition.len() != model.components().len() {
        return Err(Error::InvalidInput("composition mismatch"));
    }
    let largest = composition
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let sum: f64 = composition.iter().sum();
    if composition.len() > 1 && largest != sum {
        return Err(Error::Unsupported(
            "the critical point is only ported for a single component",
        ));
    }
    let guess = model.pseudo_critical(composition)?;
    let mut temperature = guess.temperature;
    let mut volume = guess.volume;
    for _ in 0..ITERATIONS {
        let (first, second) = conditions(model, temperature, volume, composition)?;
        let (scale_first, scale_second) = scales(model, temperature, volume, composition)?;
        if first.abs() <= TOLERANCE * scale_first && second.abs() <= TOLERANCE * scale_second {
            let state = model.prepare(temperature, composition)?;
            let pressure = state.pressure(volume)?;
            return Ok(CriticalPoint {
                temperature,
                volume,
                pressure: pressure.value,
            });
        }
        // The Jacobian of the two conditions with respect to the temperature
        // and the volume.
        let delta_temperature = STEP * temperature.abs();
        let delta_volume = STEP * volume.abs();
        let (first_up, second_up) =
            conditions(model, temperature + delta_temperature, volume, composition)?;
        let (first_down, second_down) =
            conditions(model, temperature - delta_temperature, volume, composition)?;
        let (first_right, second_right) =
            conditions(model, temperature, volume + delta_volume, composition)?;
        let (first_left, second_left) =
            conditions(model, temperature, volume - delta_volume, composition)?;
        let a11 = (first_up - first_down) / (2.0 * delta_temperature);
        let a12 = (first_right - first_left) / (2.0 * delta_volume);
        let a21 = (second_up - second_down) / (2.0 * delta_temperature);
        let a22 = (second_right - second_left) / (2.0 * delta_volume);
        let determinant = a11 * a22 - a12 * a21;
        if determinant == 0.0 || !determinant.is_finite() {
            return Err(Error::NotConverged("the critical point"));
        }
        let step_temperature = (-first * a22 + second * a12) / determinant;
        let step_volume = (-second * a11 + first * a21) / determinant;
        temperature += step_temperature;
        volume += step_volume;
        if temperature <= 0.0 || volume <= 0.0 {
            return Err(Error::NotConverged("the critical point"));
        }
    }
    Err(Error::NotConverged("the critical point"))
}

/// The two conditions of the critical point at a given temperature and volume.
fn conditions(
    model: &Cubic,
    temperature: f64,
    volume: f64,
    composition: &[f64],
) -> Result<(f64, f64)> {
    let state = model.prepare(temperature, composition)?;
    let pressure = state.pressure(volume)?;
    Ok((pressure.dv, state.pressure_second_derivative(volume)?))
}

/// The scales of the two conditions, the magnitudes of their terms.
fn scales(model: &Cubic, temperature: f64, volume: f64, composition: &[f64]) -> Result<(f64, f64)> {
    let state = model.prepare(temperature, composition)?;
    let pressure = state.pressure(volume)?;
    let first = R * temperature / (volume * volume);
    let second = 2.0 * R * temperature / (volume * volume * volume);
    Ok((
        first.abs() + pressure.dv.abs(),
        second.abs() + state.pressure_second_derivative(volume)?.abs(),
    ))
}
