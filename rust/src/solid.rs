//! Melting and sublimation line correlations of the original fork.
//! Ported from `src/solid_saturation.f90` with the parameter database
//! generated into `crate::solid_data` from `src/solid_correlation_datadb.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! ```
//! use thermopack_rs::solid;
//! let line = solid::melting_line("CO2", 300.0, 5).unwrap();
//! assert_eq!(line.temperature.len(), 5);
//! assert!(line.pressure[0] > 0.0);
//! ```

use crate::error::{Error, Result};
use crate::solid_data::{
    SolidCorrelation, lookup_melting_correlation, lookup_sublimation_correlation,
};

/// A melting or sublimation line sampled on a temperature grid.
#[derive(Debug, Clone)]
pub struct SolidLine {
    /// Temperatures of the line in K.
    pub temperature: Vec<f64>,
    /// Pressures of the line in Pa.
    pub pressure: Vec<f64>,
}

/// Melting line of a component.
///
/// The temperature grid starts at the triple point and reaches
/// `max_temperature`, limited by the correlation's maximum temperature and
/// by 1.1 times the triple point temperature. The source applies the same
/// limits.
pub fn melting_line(component: &str, max_temperature: f64, points: usize) -> Result<SolidLine> {
    let record = lookup_melting_correlation(component).ok_or_else(|| {
        Error::UnknownComponent(format!("no melting line correlation for {component}"))
    })?;
    if points < 2 {
        return Err(Error::InvalidInput("at least two line points required"));
    }
    let triple = record.triple_temperature;
    let limit = max_temperature
        .max(1.1 * triple)
        .min(record.maximum_temperature);
    let mut temperature = Vec::with_capacity(points);
    let mut pressure = Vec::with_capacity(points);
    for index in 0..points {
        let t = triple + (limit - triple) * index as f64 / (points - 1) as f64;
        temperature.push(t);
        pressure.push(melting_pressure(record, t)?);
    }
    Ok(SolidLine {
        temperature,
        pressure,
    })
}

/// Sublimation line of a component.
///
/// The temperature grid starts at `min_temperature`, limited by the
/// correlation's minimum temperature and by 0.9 times the triple point
/// temperature, and reaches the triple point.
pub fn sublimation_line(component: &str, min_temperature: f64, points: usize) -> Result<SolidLine> {
    let record = lookup_sublimation_correlation(component).ok_or_else(|| {
        Error::UnknownComponent(format!("no sublimation line correlation for {component}"))
    })?;
    if points < 2 {
        return Err(Error::InvalidInput("at least two line points required"));
    }
    let triple = record.triple_temperature;
    let start = min_temperature
        .min(0.9 * triple)
        .max(record.minimum_temperature);
    let mut temperature = Vec::with_capacity(points);
    let mut pressure = Vec::with_capacity(points);
    for index in 0..points {
        let t = start + (triple - start) * index as f64 / (points - 1) as f64;
        temperature.push(t);
        pressure.push(sublimation_pressure(record, t)?);
    }
    Ok(SolidLine {
        temperature,
        pressure,
    })
}

/// Melting pressure of a component at a temperature in Pa.
///
/// Follows `melting_pressure` of `src/solid_saturation.f90`.
pub fn melting_pressure(record: &SolidCorrelation, temperature: f64) -> Result<f64> {
    let t_red = temperature / record.reducing_temperature;
    let sum = match &record.correlation[..2] {
        "ML" => melting_sum(record, t_red),
        "MP" => {
            // The source splits the correlation at 22 K.
            if temperature > 22.0 {
                record.coeff[0] + record.coeff[1] * temperature.powf(record.exponents[1])
            } else {
                record.coeff[2] + record.coeff[3] * temperature.powf(record.exponents[3])
            }
        }
        other => {
            return Err(Error::UnsupportedDetail(format!(
                "melting line correlation: {other}"
            )));
        }
    };
    Ok(record.reducing_pressure * apply_exponent(record, sum, t_red)?)
}

/// Sublimation pressure of a component at a temperature in Pa.
///
/// Follows `sublimation_pressure` of `src/solid_saturation.f90`.
pub fn sublimation_pressure(record: &SolidCorrelation, temperature: f64) -> Result<f64> {
    let t_red = temperature / record.reducing_temperature;
    let sum = sublimation_sum(record, t_red);
    Ok(record.reducing_pressure * apply_exponent(record, sum, t_red)?)
}

/// The reduced-temperature sum of a melting line correlation.
fn melting_sum(record: &SolidCorrelation, t_red: f64) -> f64 {
    let mut sum = 0.0;
    let mut j = 0;
    for _ in 0..record.n_coeff_1 {
        sum += record.coeff[j] * t_red.powf(record.exponents[j]);
        j += 1;
    }
    for _ in 0..record.n_coeff_2 {
        sum += record.coeff[j] * (t_red - 1.0).powf(record.exponents[j]);
        j += 1;
    }
    for _ in 0..record.n_coeff_3 {
        sum += record.coeff[j] * t_red.ln().powf(record.exponents[j]);
        j += 1;
    }
    for _ in 0..(record.n_coeff - record.n_coeff_1 - record.n_coeff_2 - record.n_coeff_3) {
        sum += record.coeff[j] * (t_red.powf(record.exponents[j]) - 1.0);
        j += 1;
    }
    sum
}

/// The reduced-temperature sum of a sublimation line correlation.
fn sublimation_sum(record: &SolidCorrelation, t_red: f64) -> f64 {
    let mut sum = 0.0;
    let mut j = 0;
    for _ in 0..record.n_coeff_1 {
        sum += record.coeff[j] * t_red.powf(record.exponents[j]);
        j += 1;
    }
    for _ in 0..record.n_coeff_2 {
        sum += record.coeff[j] * (1.0 - t_red).powf(record.exponents[j]);
        j += 1;
    }
    for _ in 0..record.n_coeff_3 {
        sum += record.coeff[j] * t_red.ln().powf(record.exponents[j]);
        j += 1;
    }
    for _ in 0..(record.n_coeff - record.n_coeff_1 - record.n_coeff_2 - record.n_coeff_3) {
        sum += record.coeff[j] * (t_red.powf(record.exponents[j]) - 1.0);
        j += 1;
    }
    sum
}

/// The final exponentiation of a correlation, selected by its last character.
///
/// The original stops with an error for any other label. The `S2` label of
/// the hydrogen sulfide record has no fourth character and therefore stops
/// there as well.
fn apply_exponent(record: &SolidCorrelation, sum: f64, t_red: f64) -> Result<f64> {
    match record.correlation.as_bytes().get(3) {
        Some(b'1') => Ok(sum),
        Some(b'2') => Ok(sum.exp()),
        Some(b'3') => Ok((sum / t_red).exp()),
        _ => Err(Error::UnsupportedDetail(format!(
            "solid line correlation: {}",
            record.correlation
        ))),
    }
}
