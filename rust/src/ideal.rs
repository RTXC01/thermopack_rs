//! Analytic ideal-gas caloric correlations from ideal.f90.
//!
//! ```
//! use thermopack_rs::{Cubic, CubicKind, EntropyReference};
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let mut model = Cubic::new(CubicKind::PengRobinson, &["C1"])?;
//! let ideal = model.ideal_gas_mut(0)?; // Rust component indices start at zero.
//! ideal.set_heat_capacity(8, &[2.5])?;
//! ideal.set_formation_enthalpy(0.0)?;
//! ideal.set_standard_entropy(150.0, EntropyReference::Bar)?;
//! let values = ideal.properties(298.15, 1e5)?;
//! assert!(values.enthalpy.abs() < 1e-8);
//! assert!((values.entropy - 150.0).abs() < 1e-8);
//! # Ok(())
//! # }
//! ```
use crate::error::positive;
use crate::{Component, Error, R, Result};

#[derive(Debug, Clone, Copy)]
pub struct IdealProperties {
    /// Source Cp correlation [J/(mol K)]. Type 10 in this fork differs from
    /// the derivative of its enthalpy by a factor 1000; see `enthalpy_dt`.
    pub cp: f64,
    /// Actual dH/dT of the returned values. Native caloric solves use this.
    pub enthalpy_dt: f64,
    pub enthalpy: f64,
    pub entropy: f64,
}

/// Reference pressure accepted by the fork's standard-entropy setter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EntropyReference {
    Bar,
    Atmosphere,
}
impl EntropyReference {
    pub fn pressure(self) -> f64 {
        match self {
            Self::Bar => 100000.0,
            Self::Atmosphere => 101325.0,
        }
    }
}

/// Owned ideal-gas parameters with cached integration constants. Changing Cp
/// preserves those constants, as in set_ideal_cp_correlation. Changing a
/// reference value recalculates them, as in set_reference_energy.
#[derive(Debug, Clone)]
pub struct IdealGas {
    data: Component,
    offsets: (f64, f64),
}
impl IdealGas {
    pub fn new(component: &Component) -> Result<Self> {
        Ok(Self {
            data: *component,
            offsets: component.ideal_offsets()?.unwrap_or((0.0, 0.0)),
        })
    }

    /// Build an ideal gas with explicit enthalpy and entropy integration
    /// constants [J/mol] and [J/(mol K)]. This preserves raw source constants
    /// when initialization copies a component without reference normalization.
    pub fn with_integration_constants(
        component: &Component,
        enthalpy_offset: f64,
        entropy_offset: f64,
    ) -> Result<Self> {
        if !enthalpy_offset.is_finite() || !entropy_offset.is_finite() {
            return Err(Error::InvalidInput(
                "finite ideal-gas integration constants required",
            ));
        }
        Ok(Self {
            data: *component,
            offsets: (enthalpy_offset, entropy_offset),
        })
    }

    pub fn properties(&self, temperature: f64, pressure: f64) -> Result<IdealProperties> {
        self.data
            .ideal_with_offsets(temperature, pressure, self.offsets)
    }

    /// Correlation identifier and all 21 coefficients, including zero padding.
    pub fn heat_capacity_correlation(&self) -> (u8, &[f64; 21]) {
        (self.data.cp_type, &self.data.cp)
    }

    /// Set a supported correlation. Short arrays are zero padded. Unlike the
    /// Python wrapper, oversized arrays return an error instead of truncating.
    /// Integration constants remain unchanged until a reference setter is used.
    pub fn set_heat_capacity(&mut self, correlation: u8, coefficients: &[f64]) -> Result<()> {
        if !(1..=12).contains(&correlation) {
            return Err(Error::Unsupported("ideal heat-capacity correlation"));
        }
        if coefficients.len() > 21 || coefficients.iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput(
                "at most 21 finite Cp coefficients required",
            ));
        }
        if correlation == 10 && !["P-H2", "N-H2", "O-H2", "E-H2"].contains(&self.data.ident) {
            return Err(Error::InvalidInput("type-10 Cp requires a hydrogen isomer"));
        }
        self.data.cp_type = correlation;
        self.data.cp.fill(0.0);
        self.data.cp[..coefficients.len()].copy_from_slice(coefficients);
        Ok(())
    }

    pub fn set_standard_entropy(
        &mut self,
        entropy: f64,
        reference: EntropyReference,
    ) -> Result<()> {
        if !entropy.is_finite() {
            return Err(Error::InvalidInput("finite standard entropy required"));
        }
        let mut data = self.data;
        data.reference_entropy = entropy;
        data.reference_pressure = reference.pressure();
        self.update_reference(data)
    }

    /// Preserve get_standard_entropy's reference-pressure conversion, including
    /// its minus sign for 1ATM. Use properties(298.15, 1e5) for actual S at 1bar.
    pub fn standard_entropy(&self) -> f64 {
        self.data.reference_entropy - R * (self.data.reference_pressure / 1e5).ln()
    }

    pub fn set_formation_enthalpy(&mut self, enthalpy: f64) -> Result<()> {
        if !enthalpy.is_finite() {
            return Err(Error::InvalidInput("finite formation enthalpy required"));
        }
        let mut data = self.data;
        data.reference_enthalpy = enthalpy;
        self.update_reference(data)
    }

    pub fn formation_enthalpy(&self) -> f64 {
        self.data.reference_enthalpy
    }

    fn update_reference(&mut self, data: Component) -> Result<()> {
        if let Some(offsets) = data.ideal_offsets()? {
            self.offsets = offsets;
        }
        self.data = data;
        Ok(())
    }
}

impl Component {
    /// Pure ideal gas at T (K), P (Pa), with the fork's reference state.
    pub fn ideal(&self, temperature: f64, pressure: f64) -> Result<IdealProperties> {
        self.ideal_with_offsets(
            temperature,
            pressure,
            self.ideal_offsets()?.unwrap_or((0.0, 0.0)),
        )
    }

    fn ideal_with_offsets(
        &self,
        temperature: f64,
        pressure: f64,
        offsets: (f64, f64),
    ) -> Result<IdealProperties> {
        positive(temperature, "temperature")?;
        positive(pressure, "pressure")?;
        let (cp, h, s) = self.ideal_primitives(temperature)?;
        Ok(IdealProperties {
            cp,
            enthalpy_dt: cp * if self.cp_type == 10 { 1000.0 } else { 1.0 },
            enthalpy: h + offsets.0,
            entropy: s + offsets.1 - R * pressure.ln(),
        })
    }

    fn ideal_offsets(&self) -> Result<Option<(f64, f64)>> {
        let valid_reference = (self.reference_enthalpy != 0.0 || self.reference_entropy != 0.0)
            && self.reference_enthalpy.abs() < 1e15
            && self.reference_entropy.abs() < 1e15;
        if valid_reference {
            let (_, h0, s0) = self.ideal_primitives(298.15)?;
            Ok(Some((
                self.reference_enthalpy - h0,
                self.reference_entropy - s0 + R * self.reference_pressure.ln(),
            )))
        } else {
            Ok(None)
        }
    }

    fn ideal_primitives(&self, t: f64) -> Result<(f64, f64, f64)> {
        let c = &self.cp;
        match self.cp_type {
            1 | 2 | 3 | 4 | 5 | 6 | 8 | 11 => {
                let mut p = [0.0; 5];
                let mut inv = 0.0;
                let (scale, ts) = match self.cp_type {
                    1 => {
                        p[..4].copy_from_slice(&c[..4]);
                        (4.1868, t)
                    }
                    2 => {
                        for i in 0..5 {
                            p[i] = (i + 1) as f64 * c[i + 1];
                        }
                        (self.molar_mass, t)
                    }
                    3 => {
                        p[..3].copy_from_slice(&c[..3]);
                        (4.1868 * self.molar_mass, t * 1.8)
                    }
                    4 => {
                        p[..4].copy_from_slice(&c[..4]);
                        (1.0, t)
                    }
                    5 => {
                        p[..4].copy_from_slice(&c[..4]);
                        inv = c[4];
                        (self.molar_mass, t)
                    }
                    6 => {
                        p.copy_from_slice(&c[..5]);
                        (self.molar_mass, t)
                    }
                    8 => {
                        p.copy_from_slice(&c[..5]);
                        (R, t)
                    }
                    11 => {
                        p[..4].copy_from_slice(&c[..4]);
                        inv = c[4];
                        (1.0, t * 1e-3)
                    }
                    _ => unreachable!(),
                };
                let cp = p.iter().rev().fold(0.0, |acc, ci| acc * ts + ci) + inv / (ts * ts);
                let mut h = -inv / ts;
                let mut s = p[0] * ts.ln() - inv / (2.0 * ts * ts);
                for (i, pi) in p.iter().enumerate() {
                    h += pi * crate::numerics::ipow(ts, (i + 1) as i32) / (i + 1) as f64;
                    if i > 0 {
                        s += pi * crate::numerics::ipow(ts, i as i32) / i as f64;
                    }
                }
                let hscale = match self.cp_type {
                    3 => 1.0 / 1.8,
                    11 => 1e3,
                    _ => 1.0,
                };
                let h_constant = if self.cp_type == 2 { c[0] } else { 0.0 };
                Ok((scale * cp, scale * (h * hscale + h_constant), scale * s))
            }
            7 => {
                let x = c[2].abs() / t;
                let y = c[4].abs() / t;
                // Share two Boltzmann factors across Cp, H and S instead of
                // separately calling sinh/cosh/tanh and logarithmic helpers.
                let qx = (-2.0 * x).exp();
                let qy = (-2.0 * y).exp();
                let dx = if x < 0.5 {
                    -(-2.0 * x).exp_m1()
                } else {
                    1.0 - qx
                };
                let dy = 1.0 + qy;
                let cp = c[0]
                    + 4.0 * c[1] * x * x * qx / (dx * dx)
                    + 4.0 * c[3] * y * y * qy / (dy * dy);
                let h = c[0] * t + c[1] * c[2].abs() * (1.0 + qx) / dx
                    - c[3] * c[4].abs() * (1.0 - qy) / dy;
                // These include the source's raw entropy integration constant;
                // it remains observable when a Cp setter preserves old offsets.
                let log_dx = if x < 0.5 { dx.ln() } else { (-qx).ln_1p() };
                let s = c[0] * t.ln()
                    + c[1] * (2.0 * x * qx / dx - log_dx)
                    + c[3] * (2.0 * y * qy / dy + qy.ln_1p());
                Ok((cp * 1e-3, h * 1e-3, s * 1e-3))
            }
            9 => {
                positive(t + c[3], "Mogensen logarithm domain")?;
                let cp = c[0] + c[1] * t + c[2] / (t + c[3]);
                let h = c[0] * t + c[1] * t * t / 2.0 + c[2] * (t + c[3]).ln();
                let s = c[0] * t.ln()
                    + c[1] * t
                    + if c[3] == 0.0 {
                        -c[2] / t
                    } else {
                        c[2] / c[3] * (t.ln() - (t + c[3]).ln())
                    };
                Ok((cp, h, s))
            }
            12 => {
                let (mut cp, mut h, mut s) = (c[0], c[0] * t, c[0] * t.ln());
                for pair in c[1..].as_chunks::<2>().0 {
                    let (amplitude, theta) = (pair[0], pair[1]);
                    if amplitude.abs() < 1e-12 {
                        break;
                    }
                    let x = theta / t;
                    let e = (-x).exp();
                    let denom = -(-x).exp_m1();
                    cp += amplitude * x * x * e / (denom * denom);
                    h += amplitude * theta * e / denom;
                    s += amplitude * (x * e / denom - denom.ln());
                }
                Ok((cp * R, h * R, s * R))
            }
            10 => crate::hydrogen::primitives(self.ident, t),
            _ => Err(Error::Unsupported(
                "component has no default ideal heat-capacity correlation",
            )),
        }
    }
}
