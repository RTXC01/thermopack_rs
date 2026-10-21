//! Spinodal limits of pure fluids and mixtures.
//! Ported from `src/spinodal.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! For pure fluids, the original maps the meta-stability
//! limit where the isotherm has an extremum, so the first volume derivative
//! of the pressure vanishes. The original reports the point of that limit at
//! a given pressure, and the native model reports the same point.
//!
//! The original solves the limit with a nested loop that walks the limit line
//! from an initial point, and the native model solves the two conditions of
//! the point directly, the pressure of the state and the vanishing first
//! volume derivative of the pressure.

use crate::Phase;
use crate::cubic::Cubic;
use crate::error::{Error, Result, positive};

/// A point on a local stability limit at fixed composition.
#[derive(Debug, Clone, Copy)]
pub struct SpinodalPoint {
    /// Temperature of the point in kelvin.
    pub temperature: f64,
    /// Molar volume of the point in cubic metre per mole.
    pub volume: f64,
    /// Pressure of the point in pascal.
    pub pressure: f64,
}

/// The number of bisection steps of the volume and the temperature solves.
const ITERATIONS: usize = 200;

/// The relative tolerance of the pressure solve.
const PRESSURE_TOLERANCE: f64 = 1.0e-10;

/// The local stability limit at a given positive pressure.
///
/// For mixtures this uses the minimum Helmholtz stability eigenvalue, with
/// separate mechanical liquid/vapor seeds. See `mixture_spinodal` for explicit
/// seeds, curve continuation and the original pressure-walk compatibility API.
/// For a pure fluid the phase selects the liquid limit, where the molar volume
/// is smaller than the critical volume and the isotherm has its minimum, or the vapor limit,
/// where it is larger and the isotherm has its maximum. The limit of the
/// liquid branch reaches the critical pressure at the critical temperature
/// and falls towards the zero temperature, and the limit of the vapor branch
/// reaches the critical pressure as well and falls to the zero pressure.
pub fn spinodal_point(
    model: &Cubic,
    pressure: f64,
    phase: Phase,
    composition: &[f64],
) -> Result<SpinodalPoint> {
    positive(pressure, "the spinodal pressure")?;
    if composition.len() != model.components().len()
        || composition.iter().any(|&z| !z.is_finite() || z < 0.0)
        || (composition.iter().sum::<f64>() - 1.0).abs() > 1e-10
    {
        return Err(Error::InvalidInput(
            "normalized spinodal composition required",
        ));
    }
    if phase == Phase::MinimumGibbs {
        return Err(Error::InvalidInput(
            "a liquid or vapor spinodal branch is required",
        ));
    }
    if composition.iter().filter(|&&x| x > 0.0).count() > 1 {
        return crate::mixture_spinodal::mixture_spinodal_point(
            model,
            pressure,
            phase,
            composition,
        );
    }
    mechanical_spinodal_point(model, pressure, phase, composition)
}

pub(crate) fn mechanical_spinodal_point(
    model: &Cubic,
    pressure: f64,
    phase: Phase,
    composition: &[f64],
) -> Result<SpinodalPoint> {
    let critical = model.pseudo_critical(composition)?;
    let liquid = phase == Phase::Liquid;
    // The pressure of the limit is at or below the critical pressure.
    if pressure > critical.pressure {
        return Err(Error::Unsupported(
            "the stability limit of a pure fluid lies below the critical pressure",
        ));
    }
    // The pressure of the limit rises with the temperature on both branches.
    let limit_pressure = |temperature: f64| -> Result<f64> {
        let volume = limit_volume(model, temperature, critical.volume, liquid, composition)?;
        let state = model.prepare(temperature, composition)?;
        Ok(state.pressure(volume)?.value)
    };
    // Find the lower end of the bracket, where the pressure of the limit is
    // below the requested pressure.
    let mut lower = 0.5 * critical.temperature;
    let mut difference = limit_pressure(lower)? - pressure;
    for _ in 0..ITERATIONS {
        if difference < 0.0 {
            break;
        }
        lower *= 0.5;
        difference = limit_pressure(lower)? - pressure;
    }
    let mut upper = critical.temperature;
    let mut temperature = 0.5 * (lower + upper);
    for _ in 0..ITERATIONS {
        let volume = limit_volume(model, temperature, critical.volume, liquid, composition)?;
        let state = model.prepare(temperature, composition)?;
        let current = state.pressure(volume)?;
        let residual = current.value - pressure;
        if residual.abs() <= PRESSURE_TOLERANCE * pressure {
            return Ok(SpinodalPoint {
                temperature,
                volume,
                pressure: current.value,
            });
        }
        // A higher temperature carries a higher pressure on the limit.
        if residual < 0.0 {
            lower = temperature;
        } else {
            upper = temperature;
        }
        temperature = 0.5 * (lower + upper);
    }
    Err(Error::NotConverged("the spinodal point"))
}

/// The molar volume of the stability limit at a given temperature.
///
/// The first volume derivative of the pressure is negative below the liquid
/// limit, positive between the two limits and negative above the vapor
/// limit, so the root of each branch is bracketed by the critical volume.
pub(crate) fn limit_volume(
    model: &Cubic,
    temperature: f64,
    critical_volume: f64,
    liquid: bool,
    composition: &[f64],
) -> Result<f64> {
    let state = model.prepare(temperature, composition)?;
    let mut lower = state.covolume() * (1.0 + 1.0e-10);
    let mut upper = if liquid {
        critical_volume - state.volume_translation().value
    } else {
        critical_volume * 100.0
    };
    if liquid {
        // The first volume derivative of the pressure is negative above the
        // covolume and positive at the critical volume.
        if state.pressure(lower)?.dv >= 0.0 {
            return Err(Error::NotConverged("the liquid limit volume bracket"));
        }
    } else {
        // The first volume derivative of the pressure is positive at the
        // critical volume and negative at a large volume.
        let mut value = state.pressure(upper)?.dv;
        for _ in 0..100 {
            if value <= 0.0 {
                break;
            }
            upper *= 10.0;
            value = state.pressure(upper)?.dv;
        }
        if value > 0.0 {
            return Err(Error::NotConverged("the vapor limit volume bracket"));
        }
    }
    for _ in 0..ITERATIONS {
        let middle = 0.5 * (lower + upper);
        if middle == lower || middle == upper {
            break;
        }
        let derivative = state.pressure(middle)?.dv;
        // The derivative rises with the volume on the liquid branch and falls
        // on the vapor branch.
        if (liquid && derivative < 0.0) || (!liquid && derivative > 0.0) {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    Ok(0.5 * (lower + upper))
}
