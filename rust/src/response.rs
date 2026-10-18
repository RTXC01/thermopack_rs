//! Analytic equilibrium response after eliminating phase amount/composition
//! changes with the chemical-potential Hessian.
use super::*;

/// TP response at the evaluated state. Derived identities require that state's
/// temperature/pressure. Baseline Cp retains source defects, so identities using
/// it need not differentiate returned H/S values.
#[derive(Debug, Clone, Copy)]
pub struct EquilibriumResponse {
    /// Equilibrium Cp, J/(mol K). Physical convention: dH/dT at fixed pressure.
    pub heat_capacity: f64,
    /// Equilibrium dV/dT at fixed pressure, m3/(mol K).
    pub volume_dt: f64,
    /// Equilibrium dV/dP at fixed temperature, m3/(mol Pa).
    pub volume_dp: f64,
}
impl EquilibriumResponse {
    pub fn isentropic_volume_dp(self, temperature: f64) -> f64 {
        self.volume_dp + temperature * self.volume_dt * self.volume_dt / self.heat_capacity
    }
    pub fn heat_capacity_v(self, temperature: f64) -> f64 {
        self.heat_capacity + temperature * self.volume_dt * self.volume_dt / self.volume_dp
    }
    pub fn internal_energy_dt(self, pressure: f64) -> f64 {
        self.heat_capacity - pressure * self.volume_dt
    }
    pub fn internal_energy_dp(self, temperature: f64, pressure: f64) -> f64 {
        -temperature * self.volume_dt - pressure * self.volume_dp
    }
}

impl Cubic {
    /// Analytic TP response of a converged physical equilibrium state.
    ///
    /// Includes phase redistribution in Cp, dV/dT and dV/dP. None denotes a
    /// pure two-phase state (T and P are constrained to coexistence), an endpoint
    /// split, or a singular/non-finite/nonpositive-Cp response. This does not test stability.
    /// Heat capacities differentiate actual native caloric values, including
    /// the explicitly documented differences from source Cp conventions.
    pub fn equilibrium_response(&self, flash: &FlashResult) -> Result<Option<EquilibriumResponse>> {
        self.equilibrium_response_impl(flash, false)
    }

    /// TP response retaining the fork's reported vdW/type-10 H2 Cp defects.
    /// May report negative Cp; this need not differentiate returned H/S.
    pub fn baseline_equilibrium_response(
        &self,
        flash: &FlashResult,
    ) -> Result<Option<EquilibriumResponse>> {
        self.equilibrium_response_impl(flash, true)
    }

    fn equilibrium_response_impl(
        &self,
        flash: &FlashResult,
        baseline: bool,
    ) -> Result<Option<EquilibriumResponse>> {
        let t = flash.temperature;
        let p = flash.pressure;
        positive(p, "pressure")?;
        let feed = self.prepare(t, &flash.feed)?;
        if flash.phase != FlashPhase::TwoPhase {
            let volume = feed.volume(p, Phase::MinimumGibbs)?;
            let cp = if baseline {
                feed.baseline_enthalpy_derivative_tp(volume)?
            } else {
                phase_caloric(&feed, volume)?.2
            };
            let pressure = feed.pressure(volume)?;
            return Ok(finite_response(
                cp,
                -pressure.dt / pressure.dv,
                1. / pressure.dv,
                baseline,
            ));
        }
        let bv = flash.vapor_fraction;
        let bl = flash.liquid_fraction;
        if !bv.is_finite() || !bl.is_finite() || bv < 0. || bl < 0. || (bv + bl - 1.).abs() > 1e-10
        {
            return Err(Error::InvalidInput(
                "physical equilibrium phase fractions required",
            ));
        }
        if bv == 0. || bl == 0. {
            return Ok(None);
        }
        let active: Vec<_> = flash
            .feed
            .iter()
            .enumerate()
            .filter(|(_, z)| **z > 0.)
            .map(|(i, _)| i)
            .collect();
        if active.len() == 1 {
            return Ok(None);
        }
        let liquid = self.prepare(t, &flash.liquid)?;
        let vapor = self.prepare(t, &flash.vapor)?;
        let vl = liquid.volume(p, Phase::Liquid)?;
        let vv = vapor.volume(p, Phase::Vapor)?;
        let cpl = if baseline {
            liquid.baseline_enthalpy_derivative_tp(vl)?
        } else {
            phase_caloric(&liquid, vl)?.2
        };
        let cpv = if baseline {
            vapor.baseline_enthalpy_derivative_tp(vv)?
        } else {
            phase_caloric(&vapor, vv)?.2
        };
        let pl = liquid.pressure(vl)?;
        let pv = vapor.pressure(vv)?;
        let n = flash.feed.len();
        let m = active.len();
        let mut dtl = vec![0.; n];
        let mut dtv = vec![0.; n];
        let mut dpl = vec![0.; n];
        let mut dpv = vec![0.; n];
        let mut dnl = vec![0.; n * n];
        let mut dnv = vec![0.; n * n];
        liquid.fugacity_derivatives_into(vl, &mut dtl, &mut dpl, &mut dnl)?;
        vapor.fugacity_derivatives_into(vv, &mut dtv, &mut dpv, &mut dnv)?;
        let delta_t: Vec<_> = active.iter().map(|&i| dtv[i] - dtl[i]).collect();
        let delta_p: Vec<_> = active.iter().map(|&i| dpv[i] - dpl[i]).collect();
        let mut matrix = vec![0.; m * m];
        for (ii, &i) in active.iter().enumerate() {
            for (jj, &j) in active.iter().enumerate() {
                matrix[ii * m + jj] = (dnv[i * n + j] - 1.) / bv + (dnl[i * n + j] - 1.) / bl;
            }
            matrix[ii * m + ii] += 1. / (bv * flash.vapor[i]) + 1. / (bl * flash.liquid[i]);
        }
        let mut t_matrix = matrix.clone();
        let mut response_t = delta_t.clone();
        let mut response_p = delta_p.clone();
        if numerics::solve_linear_in_place(&mut t_matrix, &mut response_t).is_err()
            || numerics::solve_linear_in_place(&mut matrix, &mut response_p).is_err()
        {
            return Ok(None);
        }
        let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
        let cp = bl * cpl + bv * cpv + R * t * t * dot(&delta_t, &response_t);
        let vt = -bl * pl.dt / pl.dv - bv * pv.dt / pv.dv - R * t * dot(&delta_p, &response_t);
        let vp = bl / pl.dv + bv / pv.dv - R * t * dot(&delta_p, &response_p);
        Ok(finite_response(cp, vt, vp, baseline))
    }
}

fn finite_response(
    heat_capacity: f64,
    volume_dt: f64,
    volume_dp: f64,
    baseline: bool,
) -> Option<EquilibriumResponse> {
    if [heat_capacity, volume_dt, volume_dp]
        .iter()
        .all(|v| v.is_finite())
        && (baseline || heat_capacity > 0.)
    {
        Some(EquilibriumResponse {
            heat_capacity,
            volume_dt,
            volume_dp,
        })
    } else {
        None
    }
}
