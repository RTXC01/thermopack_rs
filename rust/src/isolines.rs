//! Iso-lines in the temperature-pressure plane.
//! Ported from `src/isolines.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! The original maps an isotherm and an isobar of a pure fluid or of a
//! mixture. The pure fluid branch splits the pressure or the temperature grid
//! at the saturation state of the fluid, so that every point of the grid is
//! evaluated in the phase that is stable there, and the mixture branch runs a
//! temperature-pressure flash at every point of the grid and weights the
//! properties of the two phases.
//!
//! The original maps an isenthalp and an isentrope with a pressure-enthalpy
//! and a pressure-entropy flash. Those two lines are not yet ported, because
//! the native model has no pressure-enthalpy and no pressure-entropy flash.

use crate::Phase;
use crate::cubic::Cubic;
use crate::error::{Error, Result};
use crate::flash::FlashPhase;

/// A curve in the temperature-pressure plane with the properties along it.
#[derive(Debug, Clone)]
pub struct Isoline {
    /// Temperature of every point in kelvin.
    pub temperature: Vec<f64>,
    /// Pressure of every point in pascal.
    pub pressure: Vec<f64>,
    /// Molar volume of every point in cubic metre per mole.
    pub volume: Vec<f64>,
    /// Molar entropy of every point in joule per mole and kelvin.
    pub entropy: Vec<f64>,
    /// Molar enthalpy of every point in joule per mole.
    pub enthalpy: Vec<f64>,
}

impl Isoline {
    fn with_points(points: usize) -> Self {
        Self {
            temperature: vec![0.0; points],
            pressure: vec![0.0; points],
            volume: vec![0.0; points],
            entropy: vec![0.0; points],
            enthalpy: vec![0.0; points],
        }
    }
}

/// The original `isSingleComp`: a single component model, or a composition
/// whose largest mole fraction equals the sum of the mole fractions.
fn is_single_component(composition: &[f64]) -> bool {
    if composition.len() == 1 {
        return true;
    }
    let sum: f64 = composition.iter().sum();
    let largest = composition
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    largest == sum
}

/// The index of the component with the largest mole fraction, as `maxloc` of
/// the original returns the first of the largest values.
fn dominant_component(composition: &[f64]) -> usize {
    let mut index = 0;
    let mut largest = f64::NEG_INFINITY;
    for (position, &value) in composition.iter().enumerate() {
        if value > largest {
            largest = value;
            index = position;
        }
    }
    index
}

/// The equidistant points of the original `linspace`.
///
/// The original divides by the number of intervals, so a single point yields
/// a not-a-number value. The native model reproduces that value, because the
/// split of a pure fluid grid at the saturation state can leave a single
/// point in one of the two branches.
fn linspace(minimum: f64, maximum: f64, points: usize) -> Vec<f64> {
    if points == 0 {
        return Vec::new();
    }
    if points == 1 {
        return vec![f64::NAN];
    }
    let step = (maximum - minimum) / (points as f64 - 1.0);
    (0..points)
        .map(|index| minimum + step * index as f64)
        .collect()
}

/// The number of points of the vapor branch of a pure fluid, as the original
/// computes it from the saturation pressure.
fn vapor_points(points: usize, saturation: f64, minimum: f64, maximum: f64) -> usize {
    if maximum <= minimum {
        return points;
    }
    let share = points as f64 * (saturation - minimum) / (maximum - minimum);
    let mut count = share.ceil() as usize;
    if count > (points as f64 * 0.5).round() as usize {
        count = count.saturating_sub(1);
    }
    count.min(points)
}

/// Maps an isotherm of the model.
///
/// The temperature is in kelvin, the pressure limits are in pascal and the
/// composition holds the mole fractions. The returned points carry the
/// temperature, the pressure, the molar volume, the molar entropy and the
/// molar enthalpy.
pub fn isotherm(
    model: &Cubic,
    temperature: f64,
    minimum_pressure: f64,
    maximum_pressure: f64,
    composition: &[f64],
    points: usize,
) -> Result<Isoline> {
    crate::error::positive(temperature, "the isotherm temperature")?;
    crate::error::positive(minimum_pressure, "the minimum isotherm pressure")?;
    crate::error::positive(maximum_pressure, "the maximum isotherm pressure")?;
    if points < 2 {
        return Err(Error::InvalidInput("an isotherm needs at least two points"));
    }
    if composition.len() != model.components().len() {
        return Err(Error::InvalidInput("composition mismatch"));
    }
    if is_single_component(composition) {
        let pure = dominant_component(composition);
        let critical = model.components()[pure].critical_temperature;
        let mut line = Isoline::with_points(points);
        let mut pressures = Vec::with_capacity(points);
        let mut phases = Vec::with_capacity(points);
        if temperature >= critical {
            pressures.extend(linspace(minimum_pressure, maximum_pressure, points));
            phases.extend(vec![Phase::Vapor; points]);
        } else {
            let saturation = model.pseudo_saturation_pressure(temperature, composition)?;
            let bubble = saturation.pressure;
            if bubble <= minimum_pressure {
                pressures.extend(linspace(minimum_pressure, maximum_pressure, points));
                phases.extend(vec![Phase::Liquid; points]);
            } else if bubble >= maximum_pressure {
                pressures.extend(linspace(minimum_pressure, maximum_pressure, points));
                phases.extend(vec![Phase::Vapor; points]);
            } else {
                let vapor = vapor_points(points, bubble, minimum_pressure, maximum_pressure);
                let liquid = points - vapor;
                pressures.extend(linspace(minimum_pressure, bubble, vapor));
                pressures.extend(linspace(bubble, maximum_pressure, liquid));
                phases.extend(vec![Phase::Vapor; vapor]);
                phases.extend(vec![Phase::Liquid; liquid]);
            }
        }
        for (index, &pressure) in pressures.iter().enumerate() {
            line.temperature[index] = temperature;
            line.pressure[index] = pressure;
            if !pressure.is_finite() {
                // The original writes a not-a-number pressure into the grid
                // and the properties of that point follow.
                line.volume[index] = f64::NAN;
                line.enthalpy[index] = f64::NAN;
                line.entropy[index] = f64::NAN;
                continue;
            }
            let state = model.prepare(temperature, composition)?;
            let volume = state.volume(pressure, phases[index])?;
            line.volume[index] = volume;
            line.enthalpy[index] = state.enthalpy(volume)?;
            line.entropy[index] = state.entropy(volume)?;
        }
        Ok(line)
    } else {
        let mut line = Isoline::with_points(points);
        let pressures = linspace(minimum_pressure, maximum_pressure, points);
        for (index, &pressure) in pressures.iter().enumerate() {
            let flash = model.tp_flash(temperature, pressure, composition)?;
            line.temperature[index] = temperature;
            line.pressure[index] = pressure;
            if flash.phase == FlashPhase::TwoPhase {
                let (volume, enthalpy, entropy) =
                    two_phase_properties(model, temperature, pressure, &flash)?;
                line.volume[index] = volume;
                line.enthalpy[index] = enthalpy;
                line.entropy[index] = entropy;
            } else {
                let phase = if flash.phase == FlashPhase::Vapor {
                    Phase::Vapor
                } else {
                    Phase::Liquid
                };
                let state = model.prepare(temperature, composition)?;
                let volume = state.volume(pressure, phase)?;
                line.volume[index] = volume;
                line.enthalpy[index] = state.enthalpy(volume)?;
                line.entropy[index] = state.entropy(volume)?;
            }
        }
        Ok(line)
    }
}

/// Maps an isobar of the model.
///
/// The pressure is in pascal, the temperature limits are in kelvin and the
/// composition holds the mole fractions. The returned points carry the
/// temperature, the pressure, the molar volume, the molar entropy and the
/// molar enthalpy.
pub fn isobar(
    model: &Cubic,
    pressure: f64,
    minimum_temperature: f64,
    maximum_temperature: f64,
    composition: &[f64],
    points: usize,
) -> Result<Isoline> {
    crate::error::positive(pressure, "the isobar pressure")?;
    crate::error::positive(minimum_temperature, "the minimum isobar temperature")?;
    crate::error::positive(maximum_temperature, "the maximum isobar temperature")?;
    if points < 2 {
        return Err(Error::InvalidInput("an isobar needs at least two points"));
    }
    if composition.len() != model.components().len() {
        return Err(Error::InvalidInput("composition mismatch"));
    }
    if is_single_component(composition) {
        let pure = dominant_component(composition);
        let critical = model.components()[pure].critical_pressure;
        let mut line = Isoline::with_points(points);
        let mut temperatures = Vec::with_capacity(points);
        let mut phases = Vec::with_capacity(points);
        if pressure >= critical {
            temperatures.extend(linspace(minimum_temperature, maximum_temperature, points));
            phases.extend(vec![Phase::Liquid; points]);
        } else {
            let saturation = model.pseudo_saturation_temperature(pressure, composition)?;
            let bubble = saturation.temperature;
            if bubble <= minimum_temperature {
                temperatures.extend(linspace(minimum_temperature, maximum_temperature, points));
                phases.extend(vec![Phase::Vapor; points]);
            } else if bubble >= maximum_temperature {
                temperatures.extend(linspace(minimum_temperature, maximum_temperature, points));
                phases.extend(vec![Phase::Liquid; points]);
            } else {
                let liquid = vapor_points(points, bubble, minimum_temperature, maximum_temperature);
                let vapor = points - liquid;
                temperatures.extend(linspace(minimum_temperature, bubble, liquid));
                temperatures.extend(linspace(bubble, maximum_temperature, vapor));
                phases.extend(vec![Phase::Liquid; liquid]);
                phases.extend(vec![Phase::Vapor; vapor]);
            }
        }
        for (index, &temperature) in temperatures.iter().enumerate() {
            line.temperature[index] = temperature;
            line.pressure[index] = pressure;
            if !temperature.is_finite() {
                line.volume[index] = f64::NAN;
                line.enthalpy[index] = f64::NAN;
                line.entropy[index] = f64::NAN;
                continue;
            }
            let state = model.prepare(temperature, composition)?;
            let volume = state.volume(pressure, phases[index])?;
            line.volume[index] = volume;
            line.enthalpy[index] = state.enthalpy(volume)?;
            line.entropy[index] = state.entropy(volume)?;
        }
        Ok(line)
    } else {
        let mut line = Isoline::with_points(points);
        let temperatures = linspace(minimum_temperature, maximum_temperature, points);
        for (index, &temperature) in temperatures.iter().enumerate() {
            let flash = model.tp_flash(temperature, pressure, composition)?;
            line.temperature[index] = temperature;
            line.pressure[index] = pressure;
            if flash.phase == FlashPhase::TwoPhase {
                let (volume, enthalpy, entropy) =
                    two_phase_properties(model, temperature, pressure, &flash)?;
                line.volume[index] = volume;
                line.enthalpy[index] = enthalpy;
                line.entropy[index] = entropy;
            } else {
                let phase = if flash.phase == FlashPhase::Vapor {
                    Phase::Vapor
                } else {
                    Phase::Liquid
                };
                let state = model.prepare(temperature, composition)?;
                let volume = state.volume(pressure, phase)?;
                line.volume[index] = volume;
                line.enthalpy[index] = state.enthalpy(volume)?;
                line.entropy[index] = state.entropy(volume)?;
            }
        }
        Ok(line)
    }
}

/// The specification and the limits of a caloric iso-line.
#[derive(Debug, Clone, Copy)]
pub struct CaloricLine {
    /// The specified molar enthalpy in joule per mole, or the specified molar
    /// entropy in joule per mole and kelvin.
    pub specification: f64,
    /// The lowest pressure of the grid in pascal.
    pub minimum_pressure: f64,
    /// The highest pressure of the grid in pascal.
    pub maximum_pressure: f64,
    /// The lowest temperature of the flash in kelvin.
    pub minimum_temperature: f64,
    /// The highest temperature of the flash in kelvin.
    pub maximum_temperature: f64,
}

/// Maps an isenthalp of the model.
///
/// The specification of the line is the molar enthalpy in joule per mole, the
/// limits bound the pressure grid in pascal and the temperature limits bound
/// the flash in kelvin, and the composition holds the mole fractions. Every
/// point of the pressure grid carries a pressure-enthalpy flash, and the
/// points that the flash cannot solve are dropped, as in the original.
pub fn isenthalp(
    model: &Cubic,
    line: CaloricLine,
    composition: &[f64],
    points: usize,
) -> Result<Isoline> {
    crate::error::positive(line.minimum_pressure, "the minimum isenthalp pressure")?;
    crate::error::positive(line.maximum_pressure, "the maximum isenthalp pressure")?;
    crate::error::positive(
        line.minimum_temperature,
        "the minimum isenthalp temperature",
    )?;
    crate::error::positive(
        line.maximum_temperature,
        "the maximum isenthalp temperature",
    )?;
    if points < 2 {
        return Err(Error::InvalidInput(
            "an isenthalp needs at least two points",
        ));
    }
    if composition.len() != model.components().len() {
        return Err(Error::InvalidInput("composition mismatch"));
    }
    let pressures = linspace(line.minimum_pressure, line.maximum_pressure, points);
    let mut result = Isoline::with_points(0);
    let mut temperature = 0.5 * (line.maximum_temperature + line.minimum_temperature);
    for &pressure in &pressures {
        let options = crate::caloric::CaloricFlashOptions {
            temperature_guess: Some(temperature),
            minimum_temperature: line.minimum_temperature,
            maximum_temperature: line.maximum_temperature,
            ..Default::default()
        };
        let flash =
            match model.ph_flash_with_options(pressure, composition, line.specification, options) {
                Ok(flash) => flash,
                Err(_) => continue,
            };
        temperature = flash.temperature;
        let (volume, _, entropy) = flash_properties(model, &flash, pressure)?;
        result.temperature.push(flash.temperature);
        result.pressure.push(pressure);
        result.volume.push(volume);
        result.entropy.push(entropy);
        result.enthalpy.push(line.specification);
    }
    Ok(result)
}

/// Maps an isentrope of the model.
///
/// The specification of the line is the molar entropy in joule per mole and
/// kelvin, the limits bound the pressure grid in pascal and the temperature
/// limits bound the flash in kelvin, and the composition holds the mole
/// fractions. Every point of the pressure grid carries a pressure-entropy
/// flash, and the points that the flash cannot solve are dropped, as in the
/// original.
pub fn isentrope(
    model: &Cubic,
    line: CaloricLine,
    composition: &[f64],
    points: usize,
) -> Result<Isoline> {
    crate::error::positive(line.minimum_pressure, "the minimum isentrope pressure")?;
    crate::error::positive(line.maximum_pressure, "the maximum isentrope pressure")?;
    crate::error::positive(
        line.minimum_temperature,
        "the minimum isentrope temperature",
    )?;
    crate::error::positive(
        line.maximum_temperature,
        "the maximum isentrope temperature",
    )?;
    if points < 2 {
        return Err(Error::InvalidInput(
            "an isentrope needs at least two points",
        ));
    }
    if composition.len() != model.components().len() {
        return Err(Error::InvalidInput("composition mismatch"));
    }
    let pressures = linspace(line.minimum_pressure, line.maximum_pressure, points);
    let mut result = Isoline::with_points(0);
    let mut temperature = 0.5 * (line.maximum_temperature + line.minimum_temperature);
    for &pressure in &pressures {
        let options = crate::caloric::CaloricFlashOptions {
            temperature_guess: Some(temperature),
            minimum_temperature: line.minimum_temperature,
            maximum_temperature: line.maximum_temperature,
            ..Default::default()
        };
        let flash =
            match model.ps_flash_with_options(pressure, composition, line.specification, options) {
                Ok(flash) => flash,
                Err(_) => continue,
            };
        temperature = flash.temperature;
        let (volume, enthalpy, _) = flash_properties(model, &flash, pressure)?;
        result.temperature.push(flash.temperature);
        result.pressure.push(pressure);
        result.volume.push(volume);
        result.entropy.push(line.specification);
        result.enthalpy.push(enthalpy);
    }
    Ok(result)
}

/// The phase weighted molar volume, enthalpy and entropy of a flash state.
fn flash_properties(
    model: &Cubic,
    flash: &crate::flash::FlashResult,
    pressure: f64,
) -> Result<(f64, f64, f64)> {
    if flash.phase == FlashPhase::TwoPhase {
        two_phase_properties(model, flash.temperature, pressure, flash)
    } else {
        let phase = if flash.phase == FlashPhase::Vapor {
            Phase::Vapor
        } else {
            Phase::Liquid
        };
        let state = model.prepare(flash.temperature, &flash.feed)?;
        let volume = state.volume(pressure, phase)?;
        Ok((volume, state.enthalpy(volume)?, state.entropy(volume)?))
    }
}

/// The phase weighted molar volume, enthalpy and entropy of a two phase state.
fn two_phase_properties(
    model: &Cubic,
    temperature: f64,
    pressure: f64,
    flash: &crate::flash::FlashResult,
) -> Result<(f64, f64, f64)> {
    let vapor_state = model.prepare(temperature, &flash.vapor)?;
    let vapor_volume = vapor_state.volume(pressure, Phase::Vapor)?;
    let vapor_enthalpy = vapor_state.enthalpy(vapor_volume)?;
    let vapor_entropy = vapor_state.entropy(vapor_volume)?;
    let liquid_state = model.prepare(temperature, &flash.liquid)?;
    let liquid_volume = liquid_state.volume(pressure, Phase::Liquid)?;
    let liquid_enthalpy = liquid_state.enthalpy(liquid_volume)?;
    let liquid_entropy = liquid_state.entropy(liquid_volume)?;
    let vapor_fraction = flash.vapor_fraction;
    let liquid_fraction = flash.liquid_fraction;
    Ok((
        vapor_fraction * vapor_volume + liquid_fraction * liquid_volume,
        vapor_fraction * vapor_enthalpy + liquid_fraction * liquid_enthalpy,
        vapor_fraction * vapor_entropy + liquid_fraction * liquid_entropy,
    ))
}
