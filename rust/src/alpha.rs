//! Cubic temperature-dependent attraction correlations from cbalpha.f90.
//! Initialization follows the fork's parameter-selection and fallback rules.
//!
//! ```
//! use thermopack_rs::{Cubic, CubicKind, AlphaCorrelation, Phase};
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let mut model = Cubic::with_alpha(CubicKind::PengRobinson, &["CO2", "C1"],
//!     AlphaCorrelation::Twu, "tcPR")?;
//! assert_eq!(model.alpha(0)?.correlation, AlphaCorrelation::Twu);
//! model.set_alpha(1, AlphaCorrelation::MathiasCopeman, &[0.8, 0.2, -0.1])?;
//! let state = model.prepare(300.0, &[0.8, 0.2])?;
//! let volume = state.volume(5e6, Phase::Vapor)?;
//! assert!(volume > state.covolume());
//! # Ok(())
//! # }
//! ```
use crate::error::positive;
use crate::{Component, CubicKind, Error, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum AlphaCorrelation {
    Classic = 1,
    Twu = 2,
    MathiasCopeman = 3,
    Gerg = 4,
    ClassicFit = 5,
    Umr = 6,
    GraboskiDaubert = 7,
    RedlichKwong = 8,
    Soave = 9,
    PengRobinson = 10,
    PatelTeja = 11,
    SchmidtWenzel = 12,
    VanDerWaals = 13,
    PengRobinson78 = 14,
}
impl AlphaCorrelation {
    pub fn from_name(name: &str) -> Result<Self> {
        match name.trim().to_ascii_uppercase().as_str() {
            "CLASSIC" => Ok(Self::Classic),
            "TWU" => Ok(Self::Twu),
            "MC" => Ok(Self::MathiasCopeman),
            "GERG" => Ok(Self::Gerg),
            "CLASSICFIT" => Ok(Self::ClassicFit),
            "UMR" => Ok(Self::Umr),
            "GD" => Ok(Self::GraboskiDaubert),
            "RK" => Ok(Self::RedlichKwong),
            "SOAVE" => Ok(Self::Soave),
            "PR" => Ok(Self::PengRobinson),
            "PT" => Ok(Self::PatelTeja),
            "SW" => Ok(Self::SchmidtWenzel),
            "VDW" => Ok(Self::VanDerWaals),
            "PR78" => Ok(Self::PengRobinson78),
            _ => Err(Error::InvalidInput("unknown alpha correlation")),
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            Self::Classic => "CLASSIC",
            Self::Twu => "TWU",
            Self::MathiasCopeman => "MC",
            Self::Gerg => "GERG",
            Self::ClassicFit => "CLASSICFIT",
            Self::Umr => "UMR",
            Self::GraboskiDaubert => "GD",
            Self::RedlichKwong => "RK",
            Self::Soave => "SOAVE",
            Self::PengRobinson => "PR",
            Self::PatelTeja => "PT",
            Self::SchmidtWenzel => "SW",
            Self::VanDerWaals => "VDW",
            Self::PengRobinson78 => "PR78",
        }
    }
    pub fn parameter_count(self) -> usize {
        if matches!(self, Self::Twu | Self::MathiasCopeman | Self::Gerg) {
            3
        } else {
            1
        }
    }

    /// Acentric-factor correlation from get_alpha_param_acentric. This is
    /// separate from initialization: several source init paths leave m=0.
    pub fn acentric_parameter(self, w: f64) -> Result<f64> {
        if !w.is_finite() {
            return Err(Error::InvalidInput("finite acentric factor required"));
        }
        let [a, b, c, d, e] = match self {
            Self::Soave => [0.48, 1.574, -0.176, 0., 0.],
            Self::PengRobinson => [0.37464, 1.54226, -0.26992, 0., 0.],
            Self::GraboskiDaubert => [0.48508, 1.555171, -0.15613, 0., 0.],
            Self::VanDerWaals | Self::RedlichKwong => [1., 0., 0., 0., 0.],
            Self::SchmidtWenzel => [0.465, 1.347, -0.528, 0., 0.],
            Self::PatelTeja => [0.452413, 1.30982, -0.295937, 0., 0.],
            Self::Umr => [0.384401, 1.52276, -0.213808, 0.034616, -0.001976],
            Self::PengRobinson78 if w > 0.491 => [0.379642, 1.48503, -0.164423, 0.016666, 0.],
            Self::PengRobinson78 => [0.37464, 1.54226, -0.26992, 0., 0.],
            _ => {
                return Err(Error::InvalidInput(
                    "alpha has no acentric-factor correlation",
                ));
            }
        };
        Ok(a + b * w + c * w * w + d * w.powi(3) + e * w.powi(4))
    }
}

#[derive(Debug, Clone, Copy)]
pub struct AlphaValue {
    pub value: f64,
    pub dt: f64,
    pub dtt: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Alpha {
    pub correlation: AlphaCorrelation,
    pub parameters: [f64; 3],
    pub requested: AlphaCorrelation,
    pub used_fallback: bool,
    pub database_reference: Option<&'static str>,
}
impl Alpha {
    pub fn new(correlation: AlphaCorrelation, parameters: &[f64]) -> Result<Self> {
        if parameters.len() != correlation.parameter_count()
            || parameters.iter().any(|p| !p.is_finite())
        {
            return Err(Error::InvalidInput(
                "finite alpha parameters matching the correlation required",
            ));
        }
        let mut data = [0.; 3];
        data[..parameters.len()].copy_from_slice(parameters);
        Ok(Self {
            correlation,
            parameters: data,
            requested: correlation,
            used_fallback: false,
            database_reference: None,
        })
    }

    pub(crate) fn initialize(
        kind: CubicKind,
        component: &Component,
        requested: AlphaCorrelation,
        reference: &str,
    ) -> Result<Self> {
        use AlphaCorrelation as C;
        let classic = match kind {
            CubicKind::PengRobinson => C::PengRobinson,
            CubicKind::VanDerWaals => C::VanDerWaals,
            _ => C::Soave,
        };
        let eos = if kind == CubicKind::RedlichKwong {
            "SRK"
        } else {
            kind.name()
        };
        let mut selected = Self {
            correlation: requested,
            parameters: [0.; 3],
            requested,
            used_fallback: false,
            database_reference: None,
        };
        match requested {
            C::Classic => {
                selected.correlation = classic;
                selected.parameters[0] = classic.acentric_parameter(component.acentric_factor)?;
            }
            C::PengRobinson78 => {
                selected.parameters[0] = requested.acentric_parameter(component.acentric_factor)?
            }
            C::Twu | C::MathiasCopeman => {
                // Both paths use the TWU table in this fork, including MC.
                if let Some(record) = select_record(eos, component.ident, false, reference)
                    && record.parameters.iter().any(|p| *p != 0.)
                {
                    selected.parameters = record.parameters;
                    selected.database_reference = Some(record.reference);
                } else {
                    selected.used_fallback = true;
                }
            }
            C::Umr if component.ident == "H2O" => {
                selected.correlation = C::MathiasCopeman;
                if let Some(record) = select_record(eos, component.ident, true, reference)
                    && record.parameters.iter().any(|p| *p != 0.)
                {
                    selected.parameters = record.parameters;
                    selected.database_reference = Some(record.reference);
                } else {
                    selected.used_fallback = true;
                }
            }
            C::Umr => {
                selected.parameters[0] = requested.acentric_parameter(component.acentric_factor)?
            }
            C::Gerg => selected.used_fallback = true,
            _ => {} // The source leaves explicit SOAVE/PR/GD/PT/SW/FIT parameters zero.
        }
        if selected.used_fallback {
            selected.correlation = classic;
            selected.parameters = [
                classic.acentric_parameter(component.acentric_factor)?,
                0.,
                0.,
            ];
        }
        Ok(selected)
    }

    pub fn evaluate(&self, t: f64, tc: f64) -> Result<AlphaValue> {
        positive(t, "temperature")?;
        positive(tc, "critical temperature")?;
        use AlphaCorrelation as C;
        let tr = t / tc;
        let [a, b, c] = self.parameters;
        let (value, dt, dtt) = match self.correlation {
            C::VanDerWaals => (1., 0., 0.),
            C::RedlichKwong => {
                let alpha = 1. / tr.sqrt();
                (alpha, -0.5 * alpha / t, 0.75 * alpha / (t * t))
            }
            C::Twu => {
                let b1 = c * (b - 1.);
                let b2 = c * b;
                let a1 = tr.powf(b1);
                let a2 = tr.powf(b2);
                let alpha = a1 * (a * (1. - a2)).exp();
                (
                    alpha,
                    alpha * (b1 - a * b2 * a2) / t,
                    alpha
                        * (b1 * b1 - b1 - 2. * b1 * a * a2 * b2 - a * a2 * b2 * b2
                            + a * a2 * b2
                            + a * a * a2 * a2 * b2 * b2)
                        / (t * t),
                )
            }
            C::MathiasCopeman => {
                let (b, c) = if tr < 1. { (b, c) } else { (0., 0.) };
                let y = 1. - tr.sqrt();
                let yt = -0.5 / (tr.sqrt() * tc);
                let ytt = 0.25 / (tr.powf(1.5) * tc * tc);
                let f = 1. + a * y + b * y * y + c * y * y * y;
                let ft = a * yt + 2. * b * yt * y + 3. * c * yt * y * y;
                let ftt = a * ytt
                    + 2. * b * (ytt * y + yt * yt)
                    + 3. * c * (ytt * y * y + 2. * yt * yt * y);
                (f * f, 2. * ft * f, 2. * ftt * f + 2. * ft * ft)
            }
            C::SchmidtWenzel => {
                let sqrt = tr.sqrt();
                let arg = 5. * tr.min(1.) - 3. * a - 1.;
                let k = a + arg * arg / 70.;
                let (kt, ktt) = if tr > 1. {
                    (0., 0.)
                } else {
                    (arg / 7., 5. / 7.)
                };
                let f = 1. + k * (1. - sqrt);
                let ft = kt * (1. - sqrt) - 0.5 * k / sqrt;
                let ftt = ktt * (1. - sqrt) + 0.25 * k / (tr * sqrt) - kt / sqrt;
                (
                    f * f,
                    2. * f * ft / tc,
                    (2. * ft * ft + 2. * f * ftt) / (tc * tc),
                )
            }
            C::Gerg => {
                return Err(Error::Unsupported(
                    "explicit GERG alpha is unimplemented in the source",
                ));
            }
            C::Umr => {
                return Err(Error::Unsupported(
                    "non-water UMR alpha is unimplemented in the source",
                ));
            }
            _ => {
                let yy = -2. * a * (a + 1.);
                let sqrt = (t * tc).sqrt();
                (
                    1. + a * a * (tr - 1.) + yy * (tr.sqrt() - 1.),
                    a * a / tc + 0.5 * yy / sqrt,
                    -yy / (4. * t * sqrt),
                )
            }
        };
        if !value.is_finite() || !dt.is_finite() || !dtt.is_finite() {
            return Err(Error::InvalidInput("alpha outside finite numerical range"));
        }
        Ok(AlphaValue { value, dt, dtt })
    }
}

fn select_record(
    eos: &str,
    component: &str,
    mc: bool,
    reference: &str,
) -> Option<&'static crate::component::AlphaRecord> {
    let mut best = None;
    let mut fallback = None;
    let mut rank = usize::MAX;
    for record in crate::component::alpha_records()
        .iter()
        .filter(|r| r.eos == eos && r.component == component && r.mathias_copeman == mc)
    {
        let matching = reference
            .split('/')
            .map(str::trim)
            .filter(|s| !s.eq_ignore_ascii_case("DEFAULT"))
            .position(|s| {
                record
                    .reference
                    .split('/')
                    .any(|r| s.eq_ignore_ascii_case(r.trim()))
            });
        if let Some(position) = matching {
            if position < rank {
                best = Some(record);
                rank = position;
            }
        } else if record.reference.to_ascii_uppercase().contains("DEFAULT") {
            fallback = Some(record);
        }
    }
    best.or(fallback)
}
