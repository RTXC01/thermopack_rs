//! Peneloux volume translation and its temperature derivatives.
//!
//! ```
//! use thermopack_rs::{Cubic, CubicKind, CubicOptions, Phase};
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let model = Cubic::with_options(CubicKind::PengRobinson, &["CO2"], CubicOptions {
//!     parameter_reference: "tcPR",
//!     ..Default::default()
//! })?;
//! assert!(model.volume_translation_enabled());
//! let state = model.prepare(300.0, &[1.0])?;
//! let physical_volume = state.volume(5e6, Phase::Vapor)?;
//! let eos_volume = physical_volume + state.volume_translation().value;
//! assert!(eos_volume > state.eos_covolume());
//! # Ok(())
//! # }
//! ```
use crate::error::positive;
use crate::{Component, CubicKind, Error, R, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum TranslationKind {
    Constant = 1,
    Linear = 2,
    Quadratic = 3,
    Quintic = 6,
}
impl TranslationKind {
    pub fn from_code(code: u8) -> Result<Self> {
        match code {
            1 => Ok(Self::Constant),
            2 => Ok(Self::Linear),
            3 => Ok(Self::Quadratic),
            6 => Ok(Self::Quintic),
            _ => Err(Error::Unsupported("volume translation correlation")),
        }
    }
    fn degree(self) -> usize {
        self as usize - 1
    }
}

#[derive(Debug, Clone, Copy)]
pub struct TranslationValue {
    pub value: f64,
    pub dt: f64,
    pub dtt: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct VolumeTranslation {
    pub kind: TranslationKind,
    pub coefficients: [f64; 6],
    pub database_reference: Option<&'static str>,
    pub rackett_fallback: bool,
}
impl Default for VolumeTranslation {
    fn default() -> Self {
        Self {
            kind: TranslationKind::Constant,
            coefficients: [0.; 6],
            database_reference: None,
            rackett_fallback: false,
        }
    }
}
impl VolumeTranslation {
    pub fn new(kind: TranslationKind, coefficients: [f64; 6]) -> Result<Self> {
        if coefficients.iter().any(|c| !c.is_finite()) {
            return Err(Error::InvalidInput(
                "finite volume translation coefficients required",
            ));
        }
        Ok(Self {
            kind,
            coefficients,
            ..Default::default()
        })
    }
    #[inline]
    pub fn evaluate(&self, t: f64) -> Result<TranslationValue> {
        positive(t, "temperature")?;
        let degree = self.kind.degree();
        let (mut value, mut dt, mut dtt) = (self.coefficients[degree], 0., 0.);
        for &coefficient in self.coefficients[..degree].iter().rev() {
            dtt = dtt * t + 2. * dt;
            dt = dt * t + value;
            value = value * t + coefficient;
        }
        if !value.is_finite() || !dt.is_finite() || !dtt.is_finite() {
            return Err(Error::InvalidInput(
                "volume translation outside numerical range",
            ));
        }
        Ok(TranslationValue { value, dt, dtt })
    }
    pub(crate) fn initialize(
        kind: CubicKind,
        component: &Component,
        reference: &str,
    ) -> Result<Self> {
        let eos = if kind == CubicKind::RedlichKwong {
            "SRK"
        } else {
            kind.name()
        };
        if let Some(record) = crate::component::translation_record(eos, component.ident, reference)
        {
            return Ok(Self {
                kind: TranslationKind::from_code(record.correlation)?,
                coefficients: record.coefficients,
                database_reference: Some(record.reference),
                rackett_fallback: false,
            });
        }
        let mut result = Self {
            rackett_fallback: true,
            ..Self::default()
        };
        if component.rackett_factor > 0. {
            let (scale, intercept) = match kind {
                CubicKind::PengRobinson => (0.50033, 0.25969),
                CubicKind::SoaveRedlichKwong | CubicKind::RedlichKwong => (0.40768, 0.29441),
                _ => {
                    return Err(Error::Unsupported(
                        "Rackett translation requires PR or SRK geometry",
                    ));
                }
            };
            result.coefficients[0] = scale * R * component.critical_temperature
                / component.critical_pressure
                * (intercept - component.rackett_factor);
        }
        Ok(result)
    }
}

#[derive(Debug, Clone)]
pub(crate) struct PreparedTranslation {
    pub c: f64,
    pub ct: f64,
    pub ctt: f64,
    pub components: Vec<TranslationValue>,
}
