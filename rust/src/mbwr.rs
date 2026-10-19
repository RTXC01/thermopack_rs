//! Pure-fluid MBWR-19 (Bender) and MBWR-32 pressure/Helmholtz kernels.
//! All public densities use mol/m3; the original polynomial backend uses mol/L.
//!
//! ```
//! use thermopack_rs::{Mbwr, MbwrKind};
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let model = Mbwr::new(MbwrKind::Mbwr32, "methane")?;
//! let prepared = model.prepare(300.0)?;
//! let pressure = prepared.pressure_tv(1e-3, 1.0)?;
//! assert!(pressure.value > 0.0);
//! let original_order = prepared.baseline_residual_helmholtz(1e-3, 1.0)?;
//! let stable_integral = prepared.residual_helmholtz(1e-3, 1.0)?;
//! assert!((original_order.value - stable_integral.value).abs() < 1e-8);
//! # Ok(())
//! # }
//! ```
use crate::{DensityHelmholtz, Error, PureHelmholtz, PurePressure, R, Result, error::positive};
use std::sync::OnceLock;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MbwrKind {
    Bender19,
    Mbwr32,
}
impl MbwrKind {
    pub fn coefficient_count(self) -> usize {
        match self {
            Self::Bender19 => 19,
            Self::Mbwr32 => 32,
        }
    }
    fn terms(self) -> &'static [Term] {
        match self {
            Self::Bender19 => data::TERMS_19,
            Self::Mbwr32 => data::TERMS_32,
        }
    }
    fn exponential_count(self) -> usize {
        if self == Self::Bender19 { 2 } else { 6 }
    }
}

/// An unchanged source database record. Parameters include gamma as the final
/// value. Source order is significant: three MBWR32 methane fits are stored.
#[derive(Debug, Clone, Copy)]
pub struct MbwrRecord {
    pub kind: MbwrKind,
    pub source_name: &'static str,
    pub identifier: &'static str,
    pub equation_number: i32,
    pub set_number: i32,
    pub low_property: i32,
    pub high_property: i32,
    pub parameters: &'static [f64],
}
#[derive(Clone, Copy)]
struct Term {
    density_power: usize,
    temperature_power: f64,
    exponential: bool,
}
#[derive(Debug, Clone, Copy, Default)]
/// Original integrated alpha term, using T in K and density in mol/L.
pub struct MbwrHelmholtzTerm {
    pub coefficient: f64,
    pub density_power: i32,
    pub temperature_power: f64,
    pub exponential: bool,
}
#[path = "mbwr_data.rs"]
mod data;
#[path = "mbwr_properties.rs"]
mod properties;
pub use properties::{MbwrResidualTp, MbwrTpDerivatives, MbwrTpQuantity};
pub fn mbwr_records() -> &'static [MbwrRecord] {
    data::RECORDS
}

/// Owned coefficients. This stage provides homogeneous residual kernels;
/// critical points, density iteration and full caloric/flash integration follow.
#[derive(Debug, Clone)]
pub struct Mbwr {
    kind: MbwrKind,
    coefficients: [f64; 32],
    gamma: f64,
    record: Option<&'static MbwrRecord>,
    source_terms: [MbwrHelmholtzTerm; 40],
    source_count: usize,
}

/// Pressure (Pa) differentiated in T (K) and molar density (mol/m3).
#[derive(Debug, Clone, Copy)]
pub struct MbwrPressure {
    pub value: f64,
    pub dt: f64,
    pub drho: f64,
    pub dtt: f64,
    pub dtrho: f64,
    pub drhorho: f64,
}

/// Temperature-prepared pressure and integrated Helmholtz coefficients.
/// Density sweeps allocate no scratch memory and evaluate one exponential.
#[derive(Debug, Clone)]
pub struct PreparedMbwr<'a> {
    model: &'a Mbwr,
    temperature: f64,
    pressure_polynomial: [[f64; 3]; 8],
    pressure_exponential: [[f64; 3]; 6],
    alpha_polynomial: [[f64; 3]; 8],
    alpha_exponential: [[f64; 3]; 6],
    source_coefficients: OnceLock<[[f64; 3]; 40]>,
}

impl Mbwr {
    /// First matching database record, as in the source lookup. Compiled fluid
    /// aliases are accepted in addition to the legacy MBWR identifiers.
    pub fn new(kind: MbwrKind, identifier: &str) -> Result<Self> {
        let canonical = crate::component(identifier)
            .ok()
            .map(|c| c.ident)
            .unwrap_or(identifier);
        let record = data::RECORDS
            .iter()
            .find(|r| {
                r.kind == kind
                    && (r.identifier.eq_ignore_ascii_case(identifier)
                        || r.identifier.eq_ignore_ascii_case(canonical))
            })
            .ok_or_else(|| Error::UnknownComponent(identifier.to_owned()))?;
        Self::from_record(record.source_name)
    }
    /// Select any stored fit explicitly, including alternate methane records.
    pub fn from_record(source_name: &str) -> Result<Self> {
        let record = data::RECORDS
            .iter()
            .find(|r| r.source_name.eq_ignore_ascii_case(source_name))
            .ok_or_else(|| Error::UnknownComponent(source_name.to_owned()))?;
        let count = record.kind.coefficient_count();
        let mut model = Self::with_parameters(
            record.kind,
            &record.parameters[..count],
            record.parameters[count],
        )?;
        model.record = Some(record);
        Ok(model)
    }
    pub fn with_parameters(kind: MbwrKind, coefficients: &[f64], gamma: f64) -> Result<Self> {
        positive(gamma, "MBWR exponential gamma")?;
        if coefficients.len() != kind.coefficient_count()
            || coefficients.iter().any(|v| !v.is_finite())
        {
            return Err(Error::InvalidInput(
                "matching finite MBWR coefficient array required",
            ));
        }
        let mut values = [0.; 32];
        values[..coefficients.len()].copy_from_slice(coefficients);
        let (source_terms, source_count) = data::source_helmholtz(kind, coefficients, gamma);
        Ok(Self {
            kind,
            coefficients: values,
            gamma,
            record: None,
            source_terms,
            source_count,
        })
    }
    pub fn kind(&self) -> MbwrKind {
        self.kind
    }
    pub fn record(&self) -> Option<&'static MbwrRecord> {
        self.record
    }
    pub fn parameters(&self) -> (&[f64], f64) {
        (
            &self.coefficients[..self.kind.coefficient_count()],
            self.gamma,
        )
    }
    /// Original integrated alpha terms c*T^j*r^i, optionally multiplied by
    /// exp(-gamma*r^2), with backend density r in mol/L. Read-only source order.
    pub fn baseline_helmholtz_terms(&self) -> &[MbwrHelmholtzTerm] {
        &self.source_terms[..self.source_count]
    }
    pub fn prepare(&self, t: f64) -> Result<PreparedMbwr<'_>> {
        positive(t, "temperature")?;
        let inverse = 1. / t;
        let powers = [
            t,
            t.sqrt(),
            1.,
            inverse,
            inverse * inverse,
            inverse.powi(3),
            inverse.powi(4),
        ];
        let mut polynomial = [[0.; 3]; 8];
        let mut exponential = [[0.; 3]; 6];
        for (index, term) in self.kind.terms().iter().enumerate() {
            let j = term.temperature_power;
            let ti = if j == 1. {
                0
            } else if j == 0.5 {
                1
            } else if j == 0. {
                2
            } else {
                (-j) as usize + 2
            };
            let mut c = self.coefficients[index];
            if self.kind == MbwrKind::Bender19 {
                c *= if (1..=4).contains(&index) {
                    -1000.
                } else {
                    1000.
                };
            }
            let value = c * powers[ti];
            let values = [
                value,
                j * value * inverse,
                j * (j - 1.) * value * inverse * inverse,
            ];
            let output = if term.exponential {
                &mut exponential[(term.density_power - 3) / 2]
            } else {
                &mut polynomial[term.density_power - 2]
            };
            for (a, b) in output.iter_mut().zip(values) {
                *a += b;
            }
        }
        let integrate = |c: [f64; 3]| -> [f64; 3] {
            let value = c[0] / (1000. * R * t);
            [
                value,
                c[1] / (1000. * R * t) - value * inverse,
                c[2] / (1000. * R * t) - 2. * c[1] / (1000. * R * t) * inverse
                    + 2. * value * inverse * inverse,
            ]
        };
        let result = PreparedMbwr {
            model: self,
            temperature: t,
            pressure_polynomial: polynomial,
            pressure_exponential: exponential,
            alpha_polynomial: polynomial.map(integrate),
            alpha_exponential: exponential.map(integrate),
            source_coefficients: OnceLock::new(),
        };
        if result
            .pressure_polynomial
            .iter()
            .chain(&result.pressure_exponential)
            .chain(&result.alpha_polynomial)
            .chain(&result.alpha_exponential)
            .flatten()
            .any(|v| !v.is_finite())
        {
            return Err(Error::InvalidInput(
                "non-finite MBWR temperature coefficients",
            ));
        }
        Ok(result)
    }
    pub fn pressure_density(&self, t: f64, rho: f64) -> Result<MbwrPressure> {
        self.prepare(t)?.pressure_density(rho)
    }
    pub fn residual_density(&self, t: f64, rho: f64) -> Result<DensityHelmholtz> {
        self.prepare(t)?.residual_density(rho)
    }
    /// Original integrated-coefficient summation order, including cancellation
    /// at extreme dilute/cryogenic extrapolations. Zero density uses the limit.
    pub fn baseline_residual_density(&self, t: f64, rho: f64) -> Result<DensityHelmholtz> {
        self.prepare(t)?.baseline_residual_density(rho)
    }
    pub fn baseline_residual_helmholtz(&self, t: f64, v: f64, n: f64) -> Result<PureHelmholtz> {
        self.prepare(t)?.baseline_residual_helmholtz(v, n)
    }
    pub fn baseline_pressure_tv(&self, t: f64, v: f64, n: f64) -> Result<PurePressure> {
        self.prepare(t)?.baseline_pressure_tv(v, n)
    }
    pub fn residual_helmholtz(&self, t: f64, v: f64, n: f64) -> Result<PureHelmholtz> {
        self.prepare(t)?.residual_helmholtz(v, n)
    }
    pub fn pressure_tv(&self, t: f64, v: f64, n: f64) -> Result<PurePressure> {
        self.prepare(t)?.pressure_tv(v, n)
    }
    pub fn virial_coefficients(&self, t: f64) -> Result<(f64, f64)> {
        let a = self.residual_density(t, 0.)?;
        Ok((a.drho, a.drhorho))
    }
}

// Simultaneous polynomial value and first/second derivatives; ascending powers.
fn horner(coefficients: &[[f64; 3]], component: usize, x: f64) -> [f64; 3] {
    let mut a = [0.; 3];
    for coefficient in coefficients.iter().rev() {
        a[2] = a[2] * x + 2. * a[1];
        a[1] = a[1] * x + a[0];
        a[0] = a[0] * x + coefficient[component];
    }
    a
}

// Integral from 0 to r of s^(2m-1)*exp(-gamma*s^2), m=1..count.
// A positive series avoids cancellation in the dilute incomplete-gamma tail.
fn exponential_integrals(r: f64, gamma: f64, count: usize, exp: f64) -> [f64; 6] {
    let r2 = r * r;
    let x = gamma * r2;
    let mut output = [0.; 6];
    if x < 4. {
        let mut power = r2;
        for m in 1..=count {
            let mut term = 1. / m as f64;
            let mut sum = term;
            for k in 1..=80 {
                term *= x / (m + k) as f64;
                sum += term;
                if term <= f64::EPSILON * sum {
                    break;
                }
            }
            output[m - 1] = 0.5 * power * exp * sum;
            power *= r2;
        }
    } else {
        output[0] = -(-x).exp_m1() / (2. * gamma);
        let mut power = r2;
        for m in 2..=count {
            output[m - 1] = ((m - 1) as f64 * output[m - 2] - 0.5 * power * exp) / gamma;
            power *= r2;
        }
    }
    output
}

impl PreparedMbwr<'_> {
    pub fn baseline_residual_density(&self, rho: f64) -> Result<DensityHelmholtz> {
        if !rho.is_finite() || rho < 0. {
            return Err(Error::InvalidInput(
                "nonnegative finite molar density required",
            ));
        }
        if rho == 0. {
            return self.residual_density(rho);
        }
        let source_coefficients = self.source_coefficients.get_or_init(|| {
            let t = self.temperature;
            let mut coefficients = [[0.; 3]; 40];
            for (out, term) in coefficients
                .iter_mut()
                .zip(&self.model.source_terms)
                .take(self.model.source_count)
            {
                let j = term.temperature_power;
                let c = term.coefficient;
                *out = [
                    c * t.powf(j),
                    c * j * t.powf(j - 1.),
                    c * j * (j - 1.) * t.powf(j - 2.),
                ];
            }
            coefficients
        });
        let r = rho * 1e-3;
        let inverse = 1. / r;
        let inverse2 = inverse * inverse;
        let two_gamma = 2. * self.model.gamma;
        let two_gamma_r = two_gamma * r;
        let exp = (-self.model.gamma * r.powi(2)).exp();
        let mut result = [[0.; 3]; 3];
        for (k, row) in result.iter_mut().enumerate() {
            let mut poly = [0.; 3];
            let mut expo = [0.; 3];
            for (term, c) in self
                .model
                .source_terms
                .iter()
                .zip(source_coefficients)
                .take(self.model.source_count)
            {
                let i = term.density_power;
                let ii = i as f64;
                let value = c[k] * r.powi(i);
                if term.exponential {
                    expo[0] += value;
                    expo[1] += value * (ii * inverse - two_gamma_r);
                    expo[2] += value
                        * ((i * (i - 1)) as f64 * inverse2
                            + two_gamma * (-1. - 2. * ii + two_gamma_r * r));
                } else {
                    poly[0] += value;
                    poly[1] += ii * value;
                    poly[2] += (i * (i - 1)) as f64 * value;
                }
            }
            *row = [
                poly[0] + expo[0] * exp,
                (poly[1] * inverse + expo[1] * exp) * 1e-3,
                (poly[2] * inverse2 + expo[2] * exp) * 1e-6,
            ];
        }
        let a = DensityHelmholtz {
            value: result[0][0],
            dt: result[1][0],
            drho: result[0][1],
            dtt: result[2][0],
            dtrho: result[1][1],
            drhorho: result[0][2],
        };
        if a.values().iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput(
                "non-finite baseline MBWR Helmholtz derivatives",
            ));
        }
        Ok(a)
    }
    pub fn baseline_residual_helmholtz(&self, v: f64, n: f64) -> Result<PureHelmholtz> {
        positive(v, "total volume")?;
        positive(n, "mole number")?;
        let a = self.baseline_residual_density(n / v)?;
        // Retain the original converter's arithmetic order. Algebraic
        // rearrangement changes the cancellation artifacts in pure TP n responses.
        let f = PureHelmholtz {
            value: n * a.value,
            dt: n * a.dt,
            dv: -a.drho * n.powi(2) / (v * v),
            dn: a.value + n / v * a.drho,
            dtt: n * a.dtt,
            dtv: -n.powi(2) / (v * v) * a.dtrho,
            dtn: a.dt + n / v * a.dtrho,
            dvv: n.powi(2) / (v * v * v) * (2. * a.drho + n / v * a.drhorho),
            dvn: -(n / (v * v)) * (2. * a.drho + n / v * a.drhorho),
            dnn: 2. / v * a.drho + n / v.powi(2) * a.drhorho,
        };
        if f.values().iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput(
                "non-finite baseline MBWR extensive derivatives",
            ));
        }
        Ok(f)
    }
    pub fn baseline_pressure_tv(&self, v: f64, n: f64) -> Result<PurePressure> {
        let f = self.baseline_residual_helmholtz(v, n)?;
        let rt = R * self.temperature;
        let dv = -rt * (n / (v * v) + f.dvv);
        let p = PurePressure {
            value: rt * (n / v - f.dv),
            dt: R * (n / v - f.dv) - rt * f.dtv,
            dv,
            dn: -dv * v / n,
        };
        if [p.value, p.dt, p.dv, p.dn].iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput("non-finite baseline MBWR pressure"));
        }
        Ok(p)
    }
    pub fn pressure_density(&self, rho: f64) -> Result<MbwrPressure> {
        if !rho.is_finite() || rho < 0. {
            return Err(Error::InvalidInput(
                "nonnegative finite molar density required",
            ));
        }
        let r = rho * 1e-3;
        let r2 = r * r;
        let gamma = self.model.gamma;
        let exp = (-gamma * r2).exp();
        let evaluate = |component| {
            let b = horner(&self.pressure_polynomial, component, r);
            let c = horner(&self.pressure_exponential, component, r2);
            let q = r * r2 * c[0];
            let qr = 3. * r2 * c[0] + 2. * r2 * r2 * c[1];
            let qrr = 6. * r * c[0] + 14. * r * r2 * c[1] + 4. * r * r2 * r2 * c[2];
            [
                r2 * b[0] + exp * q,
                2. * r * b[0] + r2 * b[1] + exp * (qr - 2. * gamma * r * q),
                2. * b[0]
                    + 4. * r * b[1]
                    + r2 * b[2]
                    + exp
                        * (qrr - 4. * gamma * r * qr + (4. * gamma * gamma * r2 - 2. * gamma) * q),
            ]
        };
        let a = evaluate(0);
        let at = evaluate(1);
        let att = evaluate(2);
        let p = MbwrPressure {
            value: R * self.temperature * rho + a[0],
            dt: R * rho + at[0],
            drho: R * self.temperature + a[1] * 1e-3,
            dtt: att[0],
            dtrho: R + at[1] * 1e-3,
            drhorho: a[2] * 1e-6,
        };
        if [p.value, p.dt, p.drho, p.dtt, p.dtrho, p.drhorho]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(Error::InvalidInput("non-finite MBWR pressure derivatives"));
        }
        Ok(p)
    }
    pub fn residual_density(&self, rho: f64) -> Result<DensityHelmholtz> {
        if !rho.is_finite() || rho < 0. {
            return Err(Error::InvalidInput(
                "nonnegative finite molar density required",
            ));
        }
        let r = rho * 1e-3;
        let r2 = r * r;
        let gamma = self.model.gamma;
        let exp = (-gamma * r2).exp();
        let integrals = exponential_integrals(r, gamma, self.model.kind.exponential_count(), exp);
        let mut thermal = [0.; 3];
        let mut ar = 0.;
        let mut arr = 0.;
        let mut atr = 0.;
        for (k, thermal_value) in thermal.iter_mut().enumerate() {
            let mut poly = 0.;
            for (i, c) in self.alpha_polynomial.iter().enumerate().rev() {
                poly = poly * r + c[k] / (i + 1) as f64;
            }
            *thermal_value = r * poly
                + self
                    .alpha_exponential
                    .iter()
                    .zip(integrals)
                    .map(|(c, v)| c[k] * v)
                    .sum::<f64>();
            if k < 2 {
                let b = horner(&self.alpha_polynomial, k, r);
                let c = horner(&self.alpha_exponential, k, r2);
                let first = b[0] + exp * r * c[0];
                if k == 0 {
                    ar = first;
                    arr = b[1] + exp * ((1. - 2. * gamma * r2) * c[0] + 2. * r2 * c[1]);
                } else {
                    atr = first;
                }
            }
        }
        let result = DensityHelmholtz {
            value: thermal[0],
            dt: thermal[1],
            drho: ar * 1e-3,
            dtt: thermal[2],
            dtrho: atr * 1e-3,
            drhorho: arr * 1e-6,
        };
        if result.values().iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput("non-finite MBWR Helmholtz derivatives"));
        }
        Ok(result)
    }
    pub fn residual_helmholtz(&self, v: f64, n: f64) -> Result<PureHelmholtz> {
        positive(v, "total volume")?;
        positive(n, "mole number")?;
        self.residual_density(n / v)?.extensive(v, n)
    }
    pub fn pressure_tv(&self, v: f64, n: f64) -> Result<PurePressure> {
        positive(v, "total volume")?;
        positive(n, "mole number")?;
        let rho = n / v;
        let p = self.pressure_density(rho)?;
        let result = PurePressure {
            value: p.value,
            dt: p.dt,
            dv: -p.drho * rho / v,
            dn: p.drho / v,
        };
        if [result.value, result.dt, result.dv, result.dn]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(Error::InvalidInput("non-finite MBWR TV pressure"));
        }
        Ok(result)
    }
}
