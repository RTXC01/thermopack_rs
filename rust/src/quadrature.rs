//! Quadrature rules and the integration of a function with its differential.
//! Ported from `src/quadratures.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! The source offers the Gauss and the Gauss-Kronrod rules in seven orders
//! and integrates a function together with its differential. The native model
//! offers the same rules with the same nodes and the same weights, and the
//! same integration, which also returns the differential of the differential.

use crate::error::{Error, Result};
use crate::quadrature_data as data;

/// The quadratures of the source.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Quadrature {
    /// The Gauss rule with ten points.
    Gauss10,
    /// The Gauss rule with fifteen points.
    Gauss15,
    /// The Gauss rule with twenty points.
    Gauss20,
    /// The Gauss rule with twenty-five points.
    Gauss25,
    /// The Gauss rule with thirty points.
    Gauss30,
    /// The Gauss rule with forty points.
    Gauss40,
    /// The Gauss rule with fifty points.
    Gauss50,
    /// The Gauss-Kronrod rule with twenty-one points.
    Kronrod21,
    /// The Gauss-Kronrod rule with thirty-one points.
    Kronrod31,
    /// The Gauss-Kronrod rule with forty-one points.
    Kronrod41,
    /// The Gauss-Kronrod rule with fifty-one points.
    Kronrod51,
    /// The Gauss-Kronrod rule with sixty-one points.
    Kronrod61,
    /// The Gauss-Kronrod rule with eighty-one points.
    Kronrod81,
    /// The Gauss-Kronrod rule with one hundred and one points.
    Kronrod101,
}

impl Quadrature {
    /// The number of points of the rule.
    pub fn points(&self) -> usize {
        match self {
            Quadrature::Gauss10 => 10,
            Quadrature::Gauss15 => 15,
            Quadrature::Gauss20 => 20,
            Quadrature::Gauss25 => 25,
            Quadrature::Gauss30 => 30,
            Quadrature::Gauss40 => 40,
            Quadrature::Gauss50 => 50,
            Quadrature::Kronrod21 => 21,
            Quadrature::Kronrod31 => 31,
            Quadrature::Kronrod41 => 41,
            Quadrature::Kronrod51 => 51,
            Quadrature::Kronrod61 => 61,
            Quadrature::Kronrod81 => 81,
            Quadrature::Kronrod101 => 101,
        }
    }

    /// The Gauss rule that the error estimate of the rule uses.
    fn gauss_partner(&self) -> Option<Quadrature> {
        match self {
            Quadrature::Kronrod21 => Some(Quadrature::Gauss10),
            Quadrature::Kronrod31 => Some(Quadrature::Gauss15),
            Quadrature::Kronrod41 => Some(Quadrature::Gauss20),
            Quadrature::Kronrod51 => Some(Quadrature::Gauss25),
            Quadrature::Kronrod61 => Some(Quadrature::Gauss30),
            Quadrature::Kronrod81 => Some(Quadrature::Gauss40),
            Quadrature::Kronrod101 => Some(Quadrature::Gauss50),
            _ => None,
        }
    }

    /// The nodes of the rule on the unit interval.
    pub fn positions(&self) -> Vec<f64> {
        let (half, nodes): (usize, &[f64]) = match self {
            Quadrature::Gauss10 => (10, &data::KRONROD_21_NODES),
            Quadrature::Gauss15 => (15, &data::KRONROD_31_NODES),
            Quadrature::Gauss20 => (20, &data::KRONROD_41_NODES),
            Quadrature::Gauss25 => (25, &data::KRONROD_51_NODES),
            Quadrature::Gauss30 => (30, &data::KRONROD_61_NODES),
            Quadrature::Gauss40 => (40, &data::KRONROD_81_NODES),
            Quadrature::Gauss50 => (50, &data::KRONROD_101_NODES),
            Quadrature::Kronrod21 => (10, &data::KRONROD_21_NODES),
            Quadrature::Kronrod31 => (15, &data::KRONROD_31_NODES),
            Quadrature::Kronrod41 => (20, &data::KRONROD_41_NODES),
            Quadrature::Kronrod51 => (25, &data::KRONROD_51_NODES),
            Quadrature::Kronrod61 => (30, &data::KRONROD_61_NODES),
            Quadrature::Kronrod81 => (40, &data::KRONROD_81_NODES),
            Quadrature::Kronrod101 => (50, &data::KRONROD_101_NODES),
        };
        if self.gauss_partner().is_some() {
            // The Gauss-Kronrod rules carry the nodes of the rule on the
            // whole interval.
            let mut positions = vec![0.0; 2 * half + 1];
            positions[half] = nodes[half];
            for index in 0..half {
                positions[index] = -nodes[index];
                positions[half + 1 + index] = nodes[half - 1 - index];
            }
            positions
        } else {
            // The Gauss rules take every second node of the Kronrod rule.
            let mut positions = vec![0.0; half];
            for index in 0..half / 2 {
                positions[index] = -nodes[2 * index + 1];
                positions[half / 2 + index] = nodes[half - 1 - 2 * index];
            }
            positions
        }
    }

    /// The weights of the rule on the unit interval.
    pub fn weights(&self) -> Vec<f64> {
        let (half, source): (usize, &[f64]) = match self {
            Quadrature::Gauss10 => (10, &data::GAUSS_10_WEIGHTS),
            Quadrature::Gauss15 => (15, &data::GAUSS_15_WEIGHTS),
            Quadrature::Gauss20 => (20, &data::GAUSS_20_WEIGHTS),
            Quadrature::Gauss25 => (25, &data::GAUSS_25_WEIGHTS),
            Quadrature::Gauss30 => (30, &data::GAUSS_30_WEIGHTS),
            Quadrature::Gauss40 => (40, &data::GAUSS_40_WEIGHTS),
            Quadrature::Gauss50 => (50, &data::GAUSS_50_WEIGHTS),
            Quadrature::Kronrod21 => (10, &data::KRONROD_21_WEIGHTS),
            Quadrature::Kronrod31 => (15, &data::KRONROD_31_WEIGHTS),
            Quadrature::Kronrod41 => (20, &data::KRONROD_41_WEIGHTS),
            Quadrature::Kronrod51 => (25, &data::KRONROD_51_WEIGHTS),
            Quadrature::Kronrod61 => (30, &data::KRONROD_61_WEIGHTS),
            Quadrature::Kronrod81 => (40, &data::KRONROD_81_WEIGHTS),
            Quadrature::Kronrod101 => (50, &data::KRONROD_101_WEIGHTS),
        };
        if self.gauss_partner().is_some() {
            // The Gauss-Kronrod rules carry the weights of the rule on the
            // whole interval.
            let mut weights = vec![0.0; 2 * half + 1];
            weights[half] = source[half];
            for index in 0..half {
                weights[index] = source[index];
                weights[half + 1 + index] = source[half - 1 - index];
            }
            weights
        } else {
            // The Gauss rules mirror the weights of the first half.
            let mut weights = vec![0.0; half];
            for index in 0..half / 2 {
                weights[index] = source[index];
                weights[half / 2 + index] = source[half / 2 - 1 - index];
            }
            weights
        }
    }
}

/// The relative error estimate of a quadrature.
///
/// The estimate compares the integral of the rule with the integral of the
/// Gauss rule of the same order on the even nodes, as the source does. The
/// Gauss rules have no estimate and report zero.
pub fn quadrature_error(quadrature: Quadrature, values: &[f64], integral: f64) -> Result<f64> {
    let partner = match quadrature.gauss_partner() {
        Some(partner) => partner,
        None => return Ok(0.0),
    };
    let weights = partner.weights();
    if values.len() < weights.len() {
        return Err(Error::InvalidInput(
            "the quadrature error needs the values of every node",
        ));
    }
    let mut gauss = 0.0;
    for (index, weight) in weights.iter().enumerate() {
        gauss += weight * values[2 * index + 1];
    }
    Ok((gauss - integral) / integral)
}

/// Integrates a function and its differential over the limits.
///
/// The function returns the value of the integrand, the value of its
/// differential and the derivative of the differential with respect to the
/// integration variable. The limits of the differential are the limits of the
/// source, and the method returns the integral and the differential of the
/// integral.
pub fn integrate_with_differential<F>(
    quadrature: Quadrature,
    mut function: F,
    lower: f64,
    upper: f64,
    differential_lower: f64,
    differential_upper: f64,
) -> Result<(f64, f64)>
where
    F: FnMut(f64) -> Result<(f64, f64, f64)>,
{
    let positions = quadrature.positions();
    let weights = quadrature.weights();
    let half = 0.5 * (upper - lower);
    let differential_half = 0.5 * (differential_lower - differential_upper);
    let mut integral = 0.0;
    let mut differential = 0.0;
    for index in 0..positions.len() {
        let point = half * positions[index] + half + lower;
        let differential_point = (positions[index] + 1.0) * differential_half + differential_lower;
        let (value, value_differential, value_differential_r) = function(point)?;
        let prefactor = half * weights[index];
        let prefactor_differential = differential_half * weights[index];
        integral += prefactor * value;
        differential += prefactor * value_differential
            + prefactor_differential * value
            + prefactor * value_differential_r * differential_point;
    }
    Ok((integral, differential))
}

/// Integrates a function over the limits with a constant bound.
pub fn integrate_with_differential_const_bound<F>(
    quadrature: Quadrature,
    mut function: F,
    lower: f64,
    upper: f64,
) -> Result<(f64, f64)>
where
    F: FnMut(f64) -> Result<(f64, f64)>,
{
    let positions = quadrature.positions();
    let weights = quadrature.weights();
    let half = 0.5 * (upper - lower);
    let mut integral = 0.0;
    let mut differential = 0.0;
    for index in 0..positions.len() {
        let point = half * positions[index] + half + lower;
        let (value, value_differential) = function(point)?;
        let prefactor = half * weights[index];
        integral += prefactor * value;
        differential += prefactor * value_differential;
    }
    Ok((integral, differential))
}
