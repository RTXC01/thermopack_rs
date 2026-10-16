//! Temperature-dependent covolume factors from cbbeta.f90.
//!
//! ```
//! use thermopack_rs::{Cubic, CubicKind, CubicOptions, Phase};
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let model = Cubic::with_options(CubicKind::PengRobinson, &["HE", "H2"], CubicOptions {
//!     parameter_reference: "QuantumCubic", ..Default::default()
//! })?;
//! let beta = model.beta(1)?;
//! let tc = model.components()[1].critical_temperature;
//! assert!((beta.evaluate(tc, tc)?.value - 1.0).abs() < 1e-14);
//! let state = model.prepare(40.0, &[0.5, 0.5])?;
//! let volume = state.volume(1e5, Phase::Vapor)?;
//! assert!(volume > 0.0);
//! # Ok(())
//! # }
//! ```
use crate::error::positive;
use crate::{Component, Error, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum BetaCorrelation {
    Classic = 1,
    Quantum = 2,
}

impl BetaCorrelation {
    pub fn from_name(name: &str) -> Result<Self> {
        match name.trim().to_ascii_uppercase().as_str() {
            "CLASSIC" => Ok(Self::Classic),
            "QUANTUM" => Ok(Self::Quantum),
            _ => Err(Error::Unsupported("covolume beta correlation")),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Beta {
    pub correlation: BetaCorrelation,
    pub parameters: [f64; 2],
}
impl Default for Beta {
    fn default() -> Self {
        Self {
            correlation: BetaCorrelation::Classic,
            parameters: [0.; 2],
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct BetaValue {
    /// Dimensionless covolume factor beta.
    pub value: f64,
    /// First temperature derivative [1/K].
    pub dt: f64,
    /// Second temperature derivative [1/K2].
    pub dtt: f64,
}

impl Beta {
    pub fn new(correlation: BetaCorrelation, parameters: [f64; 2]) -> Result<Self> {
        if parameters.iter().any(|p| !p.is_finite()) {
            return Err(Error::InvalidInput("finite beta parameters required"));
        }
        Ok(Self {
            correlation,
            parameters,
        })
    }
    pub(crate) fn quantum(component: &Component) -> Self {
        let parameters = match component.ident.to_ascii_uppercase().as_str() {
            "HE" => [1.4912, 3.2634],
            "NE" => [0.4673, 2.4634],
            "H2" => [3.0696, 12.682],
            "D2" => [1.6501, 7.309],
            "P-H2" => [3.0710, 12.688],
            "O-H2" => [3.0702, 12.684],
            _ => [0., 0.],
        };
        Self {
            correlation: BetaCorrelation::Quantum,
            parameters,
        }
    }
    #[inline]
    pub fn evaluate(&self, temperature: f64, critical_temperature: f64) -> Result<BetaValue> {
        positive(temperature, "temperature")?;
        positive(critical_temperature, "critical temperature")?;
        if self.correlation == BetaCorrelation::Classic {
            return Ok(BetaValue {
                value: 1.,
                dt: 0.,
                dtt: 0.,
            });
        }
        let [xi, adjusted] = self.parameters;
        let inverse = 1. / (temperature + adjusted);
        let delta = xi * inverse;
        let delta_t = -xi * inverse * inverse;
        let delta_tt = 2. * xi * inverse * inverse * inverse;
        let factor = 1. / (1. + xi / (critical_temperature + adjusted)).powi(3);
        let one = 1. + delta;
        let value = factor * one.powi(3);
        let dt = factor * 3. * one * one * delta_t;
        let dtt = factor * (6. * one * delta_t * delta_t + 3. * one * one * delta_tt);
        if !value.is_finite() || !dt.is_finite() || !dtt.is_finite() {
            return Err(Error::InvalidInput("beta outside finite numerical range"));
        }
        Ok(BetaValue { value, dt, dtt })
    }
}

pub(crate) fn quantum_component(component: &mut Component) {
    let data = match component.ident.to_ascii_uppercase().as_str() {
        "HE" => Some((5.1953, 2.276e5, 4.002602)),
        "H2" => Some((33.145, 12.964e5, 2.01588)),
        "NE" => Some((44.4, 2661630.8085948555, 20.179)),
        "D2" => Some((38.34, 16.796e5, 4.0282)),
        "O-H2" => Some((33.22, 13.1065e5, 2.01594)),
        "P-H2" => Some((32.938, 12.858e5, 2.01588)),
        _ => None,
    };
    if let Some((t, p, mass)) = data {
        component.critical_temperature = t;
        component.critical_pressure = p;
        component.molar_mass = mass;
    }
}

#[derive(Debug, Clone)]
pub(crate) struct PreparedBeta {
    pub components: Vec<BetaValue>,
    pub dt: f64,
    pub dtt: f64,
    pub partial_t: Vec<f64>,
}
