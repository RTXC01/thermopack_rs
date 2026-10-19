//! Residual pure-fluid TP state functions from mbwr_additional.f90.
//!
//! ```
//! use thermopack_rs::{Mbwr, MbwrKind};
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let model = Mbwr::new(MbwrKind::Mbwr32, "C1")?;
//! let state = model.residual_tp(300.0, 1e-3, 1.0)?;
//! assert!((state.enthalpy.value - 300.0 * state.entropy.value
//!          - state.gibbs_energy).abs() < 1e-7);
//! let original = model.baseline_residual_tp(300.0, state.pressure, 1e-3, 1.0)?;
//! assert!(original.entropy.derivatives.is_some());
//! # Ok(())
//! # }
//! ```
use super::*;

/// Derivatives at fixed remaining T/P/n variables. H/S are extensive;
/// compressibility and fugacity quantities are intensive.
#[derive(Debug, Clone, Copy)]
pub struct MbwrTpDerivatives {
    pub dt: f64,
    pub dp: f64,
    pub dn: f64,
}
#[derive(Debug, Clone, Copy)]
pub struct MbwrTpQuantity {
    pub value: f64,
    /// None if any requested local TP response is singular/non-finite.
    pub derivatives: Option<MbwrTpDerivatives>,
}
/// Residual G/H (J), S (J/K), compressibility and ln(phi) for a pure state.
/// Ideal reference is at the same T/P/n, not at the same T/V/n.
#[derive(Debug, Clone, Copy)]
pub struct MbwrResidualTp {
    pub pressure: f64,
    pub gibbs_energy: f64,
    pub enthalpy: MbwrTpQuantity,
    pub entropy: MbwrTpQuantity,
    pub compressibility: MbwrTpQuantity,
    pub ln_fugacity_coefficient: MbwrTpQuantity,
}

fn quantity(value: f64, dt: f64, dp: f64, dn: f64) -> MbwrTpQuantity {
    MbwrTpQuantity {
        value,
        derivatives: if [dt, dp, dn].iter().all(|v| v.is_finite()) {
            Some(MbwrTpDerivatives { dt, dp, dn })
        } else {
            None
        },
    }
}

impl Mbwr {
    /// Positive-pressure homogeneous residual properties and physical local
    /// TP derivatives. The EOS determines pressure from total V and n.
    pub fn residual_tp(&self, t: f64, v: f64, n: f64) -> Result<MbwrResidualTp> {
        self.prepare(t)?.residual_tp(v, n)
    }
    /// Original MBWR state functions at explicit P, including the fork's
    /// entropy-response sign convention. P need not close the EOS.
    pub fn baseline_residual_tp(&self, t: f64, p: f64, v: f64, n: f64) -> Result<MbwrResidualTp> {
        self.prepare(t)?.baseline_residual_tp(p, v, n)
    }
}

impl PreparedMbwr<'_> {
    pub fn residual_tp(&self, v: f64, n: f64) -> Result<MbwrResidualTp> {
        let f = self.residual_helmholtz(v, n)?;
        let t = self.temperature;
        let rt = R * t;
        let delta_z = -v * f.dv / n;
        let z = 1. + delta_z;
        let p = n * rt / v * z;
        positive(p, "positive MBWR TP pressure")?;
        let log_z = delta_z.ln_1p();
        let pt = p / t - rt * f.dtv;
        let pv = -rt * (f.dvv + n / (v * v));
        let vt = -pt / pv;
        let h = -R * t * t * f.dt + n * rt * delta_z;
        let s = R * (-f.value - t * f.dt + n * log_z);
        let phi = f.dn - log_z;
        let ht = t * (vt * pt - R * (2. * f.dt + t * f.dtt + n / t));
        let hp = v - t * vt;
        let mut result = MbwrResidualTp {
            pressure: p,
            gibbs_energy: n * rt * phi,
            enthalpy: quantity(h, ht, hp, h / n),
            entropy: quantity(s, ht / t, n * R / p - vt, s / n),
            compressibility: quantity(z, -z * (1. / t - vt / v), z * (1. / p + 1. / (v * pv)), 0.),
            ln_fugacity_coefficient: quantity(phi, -h / (n * R * t * t), v / (n * rt) - 1. / p, 0.),
        };
        if pv.abs() <= 64. * f64::EPSILON * rt * n / (v * v) {
            result.enthalpy.derivatives = None;
            result.entropy.derivatives = None;
            result.compressibility.derivatives = None;
        }
        finite_values(result)
    }

    pub fn baseline_residual_tp(&self, p: f64, v: f64, n: f64) -> Result<MbwrResidualTp> {
        positive(p, "positive MBWR TP pressure")?;
        let f = self.baseline_residual_helmholtz(v, n)?;
        let t = self.temperature;
        let z = p * v / (n * R * t);
        let log_z = z.ln();
        let pt = p / t - R * t * f.dtv;
        let h_pv = -R * t * (f.dvv + n / (v * v));
        let h_pn = R * t * (-f.dvn + 1. / v);
        let h_vt = -pt / h_pv;
        let h_vn = -h_pn / h_pv;
        let h = -R * t * t * f.dt + p * v - n * R * t;
        let hq = quantity(
            h,
            t * (h_vt * pt - R * (2. * f.dt + t * f.dtt + n / t)),
            v - t * h_vt,
            t * (h_vn * pt - R * (f.dtn * t + 1.)),
        );
        // Preserve the signs in MBWR_Sres separately. They differ from the
        // actual EOS response used by MBWR_Hres/zfac/lnphi and residual_tp.
        let s_pv = -R * (t * f.dvv - n * t / (v * v));
        let s_pn = -R * (t * f.dvn + t / v);
        let s_vt = -pt / s_pv;
        let s_vn = -s_pn / s_pv;
        let sq = quantity(
            R * (-f.value - t * f.dt + n * log_z),
            s_vt * pt - R * (2. * f.dt + t * f.dtt + n / t),
            n * R / p - s_vt,
            s_vn * pt - R * (f.dn + t * f.dtn + 1. - log_z),
        );
        let pv = -R * t * f.dvv - n * R * t / (v * v);
        let pn = -R * t * f.dvn + R * t / v;
        let vt = -pt / pv;
        let vn = -pn / pv;
        let zq = quantity(
            z,
            -z * (1. / t - vt / v),
            z * (1. / p + 1. / (v * pv)),
            -z * (1. / n - vn / v),
        );
        let phi = quantity(
            f.dn - log_z,
            f.dtn + (1. - vn * pt / R) / t,
            vn / (R * t) - 1. / p,
            f.dnn + 1. / n + pv * vn.powi(2) / (R * t),
        );
        finite_values(MbwrResidualTp {
            pressure: p,
            gibbs_energy: R * t * f.value + p * v - n * R * t * (1. + log_z),
            enthalpy: hq,
            entropy: sq,
            compressibility: zq,
            ln_fugacity_coefficient: phi,
        })
    }
}

fn finite_values(p: MbwrResidualTp) -> Result<MbwrResidualTp> {
    if [
        p.pressure,
        p.gibbs_energy,
        p.enthalpy.value,
        p.entropy.value,
        p.compressibility.value,
        p.ln_fugacity_coefficient.value,
    ]
    .iter()
    .any(|v| !v.is_finite())
    {
        return Err(Error::InvalidInput("non-finite MBWR residual TP values"));
    }
    Ok(p)
}
