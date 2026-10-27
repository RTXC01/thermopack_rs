//! Total MBWR properties in physical TP coordinates and the distinct source
//! single-component TP convention.
//!
//! Ported from ThermoPack `mbwr_additional.f90`, `single_component.f90`, and
//! `eostv.f90`. Copyright thermotools and Contributors. Apache-2.0.
//!
//! ```
//! use thermopack_rs::{MbwrFluid, MbwrKind, Phase};
//! let fluid = MbwrFluid::new(MbwrKind::Mbwr32, "C1")?;
//! let prepared = fluid.prepare_tp(150., 30000.)?;
//! let q = prepared.caloric_tp(1e6, 1., Phase::Vapor)?;
//! assert!((q.pressure-1e6).abs() < 1e-3);
//! assert!((q.gibbs_energy.derivatives.unwrap().dp-q.volume.value).abs() < 1e-12);
//! # Ok::<(), thermopack_rs::Error>(())
//! ```
use super::*;
use crate::{
    MbwrDensityIteration, MbwrDensityParameters, MbwrEquilibrium, MbwrTpDerivatives,
    MbwrTpQuantity, Phase,
};

/// Total-state TP caloric properties. Extensive quantities use the supplied
/// total mole number; volume is total cubic metres and pressure is Pa.
#[derive(Debug, Clone, Copy)]
pub struct MbwrTpProperties {
    pub pressure: f64,
    pub volume: MbwrTpQuantity,
    pub helmholtz_energy: MbwrTpQuantity,
    pub internal_energy: MbwrTpQuantity,
    pub enthalpy: MbwrTpQuantity,
    pub entropy: MbwrTpQuantity,
    pub gibbs_energy: MbwrTpQuantity,
    pub chemical_potential: MbwrTpQuantity,
    pub compressibility: MbwrTpQuantity,
    pub ln_fugacity_coefficient: MbwrTpQuantity,
    pub heat_capacity_v: f64,
    pub heat_capacity_p: Option<f64>,
    pub speed_of_sound: Option<f64>,
}

/// Outputs of the original source TP routines at an explicitly imposed P.
/// H/S values and derivatives follow `single_component.f90`; residual Z and
/// ln(phi) retain their source TP implementation. `pressure` need not close
/// the EOS at `volume`. Non-finite native derivative triplets are absent.
/// Last-bit differences at an unclosed extremum can therefore differ from
/// huge finite derivative records returned by the fork.
#[derive(Debug, Clone, Copy)]
pub struct MbwrSourceTpProperties {
    pub pressure: f64,
    pub volume: f64,
    pub enthalpy: MbwrTpQuantity,
    pub entropy: MbwrTpQuantity,
    pub compressibility: MbwrTpQuantity,
    pub ln_fugacity_coefficient: MbwrTpQuantity,
}

/// Combined prepared physical caloric and bounded stable-root evaluators.
#[derive(Debug, Clone)]
pub struct PreparedMbwrFluidTp<'a> {
    pub(super) fluid: PreparedMbwrFluid<'a>,
    equilibrium: MbwrEquilibrium<'a>,
}

/// Prepared source-order TP evaluator, including the original density path.
#[derive(Debug, Clone)]
pub struct PreparedMbwrSourceTp<'a> {
    pub(super) fluid: PreparedMbwrFluid<'a>,
    density: crate::PreparedMbwrDensity,
}

#[derive(Debug, Clone)]
pub struct MbwrSourceTpState {
    pub iteration: MbwrDensityIteration,
    pub properties: Result<MbwrSourceTpProperties>,
}

impl MbwrFluid {
    /// Physical TP properties obtained from the total TV equations by the
    /// implicit EOS response. Requires positive finite pressure.
    pub fn caloric_tvp(&self, t: f64, v: f64, n: f64) -> Result<MbwrTpProperties> {
        self.prepare(t)?.caloric_tvp(v, n)
    }
    /// Prepare bounded stable MBWR pressure roots at temperature T.
    pub fn prepare_tp(&self, t: f64, maximum_density: f64) -> Result<PreparedMbwrFluidTp<'_>> {
        let fluid = self.prepare(t)?;
        let equilibrium = fluid.eos.equilibrium(maximum_density)?;
        Ok(PreparedMbwrFluidTp { fluid, equilibrium })
    }
    /// Physical branch-selected TP evaluator. `phase` chooses the low/high
    /// stable root or the minimum residual Gibbs root.
    pub fn caloric_tp(
        &self,
        t: f64,
        p: f64,
        n: f64,
        phase: Phase,
        maximum_density: f64,
    ) -> Result<MbwrTpProperties> {
        self.prepare_tp(t, maximum_density)?.caloric_tp(p, n, phase)
    }
    /// Classical source TP properties at an imposed pressure and total volume.
    pub fn baseline_caloric_tp_at_volume(
        &self,
        t: f64,
        p: f64,
        v: f64,
        n: f64,
    ) -> Result<MbwrSourceTpProperties> {
        self.prepare(t)?.baseline_caloric_tp_at_volume(p, v, n)
    }
    /// Source-order TV Hessian mapped to physical TP derivatives.
    pub fn baseline_caloric_tvp(&self, t: f64, v: f64, n: f64) -> Result<MbwrTpProperties> {
        self.prepare(t)?.baseline_caloric_tvp(v, n)
    }
    /// Prepare the source density path with explicit critical data.
    pub fn prepare_baseline_tp(
        &self,
        t: f64,
        parameters: MbwrDensityParameters,
    ) -> Result<PreparedMbwrSourceTp<'_>> {
        Ok(PreparedMbwrSourceTp {
            fluid: self.prepare(t)?,
            density: self.eos.prepare_baseline_density(t, parameters)?,
        })
    }
    /// Source density and TP properties using original cached density defaults.
    pub fn baseline_caloric_tp(
        &self,
        t: f64,
        p: f64,
        n: f64,
        phase: Phase,
    ) -> Result<MbwrSourceTpState> {
        let params = self.eos.source_density_parameters()?;
        self.baseline_caloric_tp_with_parameters(t, p, n, phase, params)
    }
    /// Source density and properties using explicit initialization data.
    pub fn baseline_caloric_tp_with_parameters(
        &self,
        t: f64,
        p: f64,
        n: f64,
        phase: Phase,
        parameters: MbwrDensityParameters,
    ) -> Result<MbwrSourceTpState> {
        self.prepare_baseline_tp(t, parameters)?
            .caloric_tp(p, n, phase)
    }
}

impl PreparedMbwrFluid<'_> {
    /// Physical TV state mapped to TP partial derivatives. Singular pressure
    /// response preserves values and suppresses derivatives requiring V(T,P,n).
    pub fn caloric_tvp(&self, v: f64, n: f64) -> Result<MbwrTpProperties> {
        let (properties, residual) = self.evaluate_with_residual(v, n, false)?;
        self.map_tv_to_tp(properties, residual, v, n)
    }
    /// Source-order TV Hessian mapped to physical TP derivatives.
    pub fn baseline_caloric_tvp(&self, v: f64, n: f64) -> Result<MbwrTpProperties> {
        let (properties, residual) = self.evaluate_with_residual(v, n, true)?;
        self.map_tv_to_tp(properties, residual, v, n)
    }
    fn map_tv_to_tp(
        &self,
        q: MbwrTvProperties,
        residual: crate::PureHelmholtz,
        v: f64,
        n: f64,
    ) -> Result<MbwrTpProperties> {
        positive(q.pressure.value, "positive MBWR TP pressure")?;
        let t = self.temperature;
        let pv = q.pressure.dv;
        let threshold = 64.0 * f64::EPSILON * R * t * n / (v * v);
        let regular = pv.is_finite() && pv.abs() > threshold;
        let vt = if regular { -q.pressure.dt / pv } else { 0.0 };
        let vp = if regular { 1.0 / pv } else { 0.0 };
        let vn = if regular { -q.pressure.dn / pv } else { 0.0 };
        let conv = |x: MbwrTvQuantity| MbwrTpQuantity {
            value: x.value,
            derivatives: if regular {
                let dt = x.dt + x.dv * vt;
                let dp = x.dv * vp;
                let dn = x.dn + x.dv * vn;
                [dt, dp, dn]
                    .iter()
                    .all(|x| x.is_finite())
                    .then_some(MbwrTpDerivatives { dt, dp, dn })
            } else {
                None
            },
        };
        let volume = MbwrTpQuantity {
            value: v,
            derivatives: regular
                .then_some(MbwrTpDerivatives {
                    dt: vt,
                    dp: vp,
                    dn: vn,
                })
                .filter(finite_tp_derivatives),
        };
        let z = q.pressure.value * v / (n * R * t);
        let zq = MbwrTpQuantity {
            value: z,
            derivatives: regular
                .then_some(MbwrTpDerivatives {
                    dt: z * (vt / v - 1.0 / t),
                    dp: z * (1.0 / q.pressure.value + vp / v),
                    dn: 0.0,
                })
                .filter(finite_tp_derivatives),
        };
        let phi_value = q.ln_fugacity_coefficient.ok_or(Error::InvalidInput(
            "MBWR fugacity requires positive pressure",
        ))?;
        let delta_z = -v * residual.dv / n;
        let residual_enthalpy = -R * t * t * residual.dt - R * t * v * residual.dv;
        let phi = MbwrTpQuantity {
            value: phi_value,
            derivatives: Some(MbwrTpDerivatives {
                dt: -residual_enthalpy / (n * R * t * t),
                dp: delta_z / q.pressure.value,
                dn: 0.0,
            })
            .filter(finite_tp_derivatives),
        };
        // Exact pure-fluid Gibbs identities remain valid at a spinodal even
        // though the volume response itself is singular.
        let gibbs_energy = MbwrTpQuantity {
            value: q.gibbs_energy.value,
            derivatives: Some(MbwrTpDerivatives {
                dt: -q.entropy.value,
                dp: v,
                dn: q.chemical_potential.value,
            })
            .filter(finite_tp_derivatives),
        };
        let chemical_potential = MbwrTpQuantity {
            value: q.chemical_potential.value,
            derivatives: Some(MbwrTpDerivatives {
                dt: -q.entropy.value / n,
                dp: v / n,
                dn: 0.0,
            })
            .filter(finite_tp_derivatives),
        };
        Ok(MbwrTpProperties {
            pressure: q.pressure.value,
            volume,
            helmholtz_energy: conv(q.helmholtz_energy),
            internal_energy: conv(q.internal_energy),
            enthalpy: conv(q.enthalpy),
            entropy: conv(q.entropy),
            gibbs_energy,
            chemical_potential,
            compressibility: zq,
            ln_fugacity_coefficient: phi,
            heat_capacity_v: q.heat_capacity_v,
            heat_capacity_p: q.heat_capacity_p,
            speed_of_sound: q.speed_of_sound,
        })
    }
    /// Source-order classical TP total H/S at imposed pressure.
    pub fn baseline_caloric_tp_at_volume(
        &self,
        p: f64,
        v: f64,
        n: f64,
    ) -> Result<MbwrSourceTpProperties> {
        positive(n, "mole number")?;
        let r = self.eos.baseline_residual_tp(p, v, n)?;
        let ideal_h = n * self.ideal.enthalpy;
        let ideal_s = n * (self.ideal.entropy - R * p.ln());
        let mut h = r.enthalpy;
        h.value += ideal_h;
        if !h.value.is_finite() {
            return Err(Error::InvalidInput("non-finite total source TP enthalpy"));
        }
        h.derivatives = h.derivatives.and_then(|mut d| {
            d.dt += n * self.source_cp;
            finite_tp_derivatives(&d).then_some(d)
        });
        let mut s = r.entropy;
        s.value += ideal_s;
        if !s.value.is_finite() {
            return Err(Error::InvalidInput("non-finite total source TP entropy"));
        }
        s.derivatives = s.derivatives.and_then(|mut d| {
            d.dt += n * self.source_cp / self.temperature;
            finite_tp_derivatives(&d).then_some(d)
        });
        Ok(MbwrSourceTpProperties {
            pressure: p,
            volume: v,
            enthalpy: h,
            entropy: s,
            compressibility: r.compressibility,
            ln_fugacity_coefficient: r.ln_fugacity_coefficient,
        })
    }
}

impl PreparedMbwrFluidTp<'_> {
    pub fn caloric_tvp(&self, v: f64, n: f64) -> Result<MbwrTpProperties> {
        self.fluid.caloric_tvp(v, n)
    }
    pub fn baseline_caloric_tvp(&self, v: f64, n: f64) -> Result<MbwrTpProperties> {
        self.fluid.baseline_caloric_tvp(v, n)
    }
    /// Choose the requested stable root and evaluate at total volume n*Vm.
    pub fn caloric_tp(&self, p: f64, n: f64, phase: Phase) -> Result<MbwrTpProperties> {
        positive(n, "mole number")?;
        let vm = self.equilibrium.volume(p, phase)?;
        self.fluid.caloric_tvp(n * vm, n)
    }
}

impl PreparedMbwrSourceTp<'_> {
    pub fn baseline_caloric_tp_at_volume(
        &self,
        p: f64,
        v: f64,
        n: f64,
    ) -> Result<MbwrSourceTpProperties> {
        self.fluid.baseline_caloric_tp_at_volume(p, v, n)
    }
    /// Preserve raw finite source density diagnostics even if caloric assembly fails.
    pub fn caloric_tp(&self, p: f64, n: f64, phase: Phase) -> Result<MbwrSourceTpState> {
        positive(n, "mole number")?;
        let iteration = self.density.density(p, phase)?;
        let v = n / iteration.density;
        let properties = self.fluid.baseline_caloric_tp_at_volume(p, v, n);
        Ok(MbwrSourceTpState {
            iteration,
            properties,
        })
    }
}

fn finite_tp_derivatives(d: &MbwrTpDerivatives) -> bool {
    [d.dt, d.dp, d.dn].iter().all(|x| x.is_finite())
}
