//! Fundamental-measure free-energy densities from `src/fmt.f90`.
//! Weighted densities are ordered n0, n1, n2, n3, nv1, nv2, with the
//! source's one-dimensional vector contraction. No temperature factor is
//! included: the result is the reduced excess free-energy density Phi.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! ```
//! use thermopack_rs::fmt::FmtModel;
//! let state = FmtModel::WhiteBearII.energy_density([0.1, 0.2, 0.3, 0.1, 0., 0.])?;
//! assert!(state.value > 0.);
//! assert_eq!(state.hessian[1][3], state.hessian[3][1]);
//! # Ok::<(), thermopack_rs::Error>(())
//! ```
use crate::{Error, Result};
use std::f64::consts::PI;

#[path = "fmt_baseline.rs"]
mod baseline;

/// Model identifiers correspond to the original FMT_ROSENFELD/WB/WBII codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum FmtModel {
    Rosenfeld = 1,
    WhiteBear = 2,
    WhiteBearII = 3,
}

#[derive(Debug, Clone, Copy)]
pub struct FmtDensity {
    pub value: f64,
    /// Derivatives with respect to n0, n1, n2, n3, nv1, nv2.
    pub gradient: [f64; 6],
    pub hessian: [[f64; 6]; 6],
}

// A single-variable second-order jet is sufficient: only the coefficient
// functions depend on n3. The remaining weighted-density algebra is explicit.
#[derive(Clone, Copy)]
struct J {
    v: f64,
    d: f64,
    dd: f64,
}
impl J {
    fn constant(v: f64) -> Self {
        Self { v, d: 0., dd: 0. }
    }
    fn scale(self, c: f64) -> Self {
        Self {
            v: self.v * c,
            d: self.d * c,
            dd: self.dd * c,
        }
    }
    fn inverse(self) -> Self {
        let x = self.v;
        Self {
            v: 1. / self.v,
            d: -self.d / (x * x),
            dd: -self.dd / (x * x) + self.d * self.d * 2. / (x * x * x),
        }
    }
    fn log(self) -> Self {
        let d = self.d / self.v;
        Self {
            v: self.v.ln(),
            d,
            dd: self.dd / self.v - d * d,
        }
    }
}
impl std::ops::Add for J {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            v: self.v + b.v,
            d: self.d + b.d,
            dd: self.dd + b.dd,
        }
    }
}
impl std::ops::Sub for J {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self {
            v: self.v - b.v,
            d: self.d - b.d,
            dd: self.dd - b.dd,
        }
    }
}
impl std::ops::Mul for J {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self {
            v: self.v * b.v,
            d: self.v * b.d + self.d * b.v,
            dd: self.v * b.dd + 2. * self.d * b.d + self.dd * b.v,
        }
    }
}

impl FmtModel {
    pub fn from_code(code: i32) -> Result<Self> {
        match code {
            1 => Ok(Self::Rosenfeld),
            2 => Ok(Self::WhiteBear),
            3 => Ok(Self::WhiteBearII),
            _ => Err(Error::InvalidInput("unknown FMT model")),
        }
    }

    /// Evaluate Phi and its first/second derivatives, without allocations.
    /// Keeps the source max(1-n3,1e-12) branch and raw coefficient formulas.
    /// White Bear variants are undefined at n3=0; no dilute series is implied.
    pub fn energy_density(self, weighted: [f64; 6]) -> Result<FmtDensity> {
        self.evaluate(weighted)
    }

    /// Preserve the original hyperdual reciprocal-derivative floor as well.
    /// Follows the original 21 hyperdual evaluations and operand ordering.
    /// Near tiny denominators this differs from mathematical differentiation.
    pub fn baseline_energy_density(self, weighted: [f64; 6]) -> Result<FmtDensity> {
        baseline::evaluate(self, weighted)
    }

    /// Evaluate a grid into caller-owned output; lengths must agree. On a
    /// failed point, the preceding output prefix has already been updated.
    pub fn energy_density_grid(
        self,
        weighted: &[[f64; 6]],
        output: &mut [FmtDensity],
    ) -> Result<()> {
        self.grid(weighted, output, false)
    }

    /// Source-compatible grid evaluation into caller-owned output. Failure
    /// semantics match `energy_density_grid`.
    pub fn baseline_energy_density_grid(
        self,
        weighted: &[[f64; 6]],
        output: &mut [FmtDensity],
    ) -> Result<()> {
        self.grid(weighted, output, true)
    }

    fn grid(self, weighted: &[[f64; 6]], output: &mut [FmtDensity], baseline: bool) -> Result<()> {
        if weighted.len() != output.len() {
            return Err(Error::InvalidInput("FMT input/output grid lengths differ"));
        }
        for (&n, out) in weighted.iter().zip(output) {
            *out = if baseline {
                self.baseline_energy_density(n)?
            } else {
                self.energy_density(n)?
            };
        }
        Ok(())
    }

    fn evaluate(self, n: [f64; 6]) -> Result<FmtDensity> {
        if n.iter().any(|x| !x.is_finite()) {
            return Err(Error::InvalidInput(
                "finite FMT weighted densities required",
            ));
        }
        let eta = J {
            v: n[3],
            d: 1.,
            dd: 0.,
        };
        let one = J::constant(1.);
        // max_dr chooses the constant branch also on exact equality.
        let u = if 1. - n[3] > 1e-12 {
            one - eta
        } else {
            J::constant(1e-12)
        };
        let log = u.log();
        let u2 = u * u;
        let (b, c) = match self {
            Self::Rosenfeld => (u.inverse(), u2.scale(24. * PI).inverse()),
            Self::WhiteBear => (
                u.inverse(),
                (eta + u2 * log) * ((eta * eta).scale(36. * PI) * u2).inverse(),
            ),
            Self::WhiteBearII => {
                let eta2 = eta * eta;
                let phi2 =
                    (J::constant(2.) - eta + (u.scale(2.) * log) * eta.inverse()).scale(1. / 3.);
                let phi3 =
                    (eta.scale(2.) - eta2.scale(3.) + eta.scale(2.) * eta2 + u2.scale(2.) * log)
                        .scale(1. / 3.)
                        * eta2.inverse();
                (
                    (one + phi2) * u.inverse(),
                    (one - phi3) * u2.scale(24. * PI).inverse(),
                )
            }
        };
        let s = n[1] * n[2] - n[4] * n[5];
        let q = n[2] * n[2] * n[2] - 3. * n[2] * n[5] * n[5];
        let q2 = 3. * n[2] * n[2] - 3. * n[5] * n[5];
        let q5 = -6. * n[2] * n[5];
        let gradient = [
            -log.v,
            n[2] * b.v,
            n[1] * b.v + q2 * c.v,
            -n[0] * log.d + s * b.d + q * c.d,
            -n[5] * b.v,
            -n[4] * b.v + q5 * c.v,
        ];
        let mut h = [[0.; 6]; 6];
        h[1][2] = b.v;
        h[4][5] = -b.v;
        h[2][2] = 6. * n[2] * c.v;
        h[2][5] = -6. * n[5] * c.v;
        h[5][5] = -6. * n[2] * c.v;
        let mixed = [
            -log.d,
            n[2] * b.d,
            n[1] * b.d + q2 * c.d,
            -n[0] * log.dd + s * b.dd + q * c.dd,
            -n[5] * b.d,
            -n[4] * b.d + q5 * c.d,
        ];
        for (i, v) in mixed.into_iter().enumerate() {
            h[i][3] = v;
            h[3][i] = v;
        }
        h[2][1] = h[1][2];
        h[5][4] = h[4][5];
        h[5][2] = h[2][5];
        let value = -n[0] * log.v + s * b.v + q * c.v;
        if !value.is_finite()
            || gradient
                .iter()
                .chain(h.iter().flatten())
                .any(|v| !v.is_finite())
        {
            return Err(Error::InvalidInput("non-finite FMT energy or derivatives"));
        }
        Ok(FmtDensity {
            value,
            gradient,
            hessian: h,
        })
    }
}
