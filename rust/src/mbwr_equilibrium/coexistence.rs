// Copyright (c) 2022-present thermotools and Contributors.
// SPDX-License-Identifier: Apache-2.0
//
//! Pure MBWR coexistence on analytically isolated mechanical branches.
//! This is a native Rust implementation based on the MBWR pressure and
//! residual Helmholtz equations; see the repository LICENSE for terms. Global
//! stability comparisons apply only to stable roots inside the caller's density
//! search domain; that domain is not a fitted EOS validity claim.
//!
//! ```
//! use thermopack_rs::Mbwr;
//! let model = Mbwr::from_record("mbwrcx1")?;
//! let sat = model.saturation_pressure(150., 45_000.)?;
//! assert!(sat.liquid_volume < sat.vapor_volume);
//! assert!(sat.pressure > 0.);
//! # Ok::<(), thermopack_rs::Error>(())
//! ```
use super::*;
use crate::Saturation;

impl MbwrEquilibrium<'_> {
    /// Solve residual chemical-potential equality from the vapor branch to any
    /// higher stable branch, selecting a globally stable candidate in-domain.
    pub fn saturation_pressure(&self) -> Result<Saturation> {
        self.saturation_pressure_impl(None)
    }

    /// As above, constrained to the caller's explicit pressure search window in Pa.
    pub fn saturation_pressure_bounded(
        &self,
        min_pressure: f64,
        max_pressure: f64,
    ) -> Result<Saturation> {
        positive(min_pressure, "minimum saturation pressure")?;
        positive(max_pressure, "maximum saturation pressure")?;
        if min_pressure >= max_pressure {
            return Err(Error::InvalidInput(
                "saturation pressure bounds must be increasing",
            ));
        }
        self.saturation_pressure_impl(Some((min_pressure, max_pressure)))
    }

    fn saturation_pressure_impl(&self, requested: Option<(f64, f64)>) -> Result<Saturation> {
        if self.intervals.len() < 2 || self.spinodals.len() < 2 {
            return Err(Error::NotConverged("distinct MBWR vapor-liquid branches"));
        }
        let vapor = &self.intervals[0];
        let mut candidates: Vec<(Saturation, f64, f64)> = Vec::new();
        let mut last_error = Error::NotConverged("MBWR coexistence pressure bracket");
        for liquid in self.intervals.iter().skip(1) {
            match self.solve_pair(vapor, liquid, requested) {
                Ok((sat, rv, rl)) => {
                    if self.globally_stable_candidate(sat.pressure, rv, rl)? {
                        candidates.push((sat, rv, rl));
                    }
                }
                Err(e @ (Error::NoPhysicalRoot | Error::NotConverged(_))) => last_error = e,
                Err(e) => return Err(e),
            }
        }
        if candidates.is_empty() {
            return Err(
                if last_error == Error::NotConverged("MBWR coexistence pressure bracket") {
                    Error::NoPhysicalRoot
                } else {
                    last_error
                },
            );
        }
        candidates.sort_by(|a, b| a.0.pressure.total_cmp(&b.0.pressure));
        if candidates.len() > 1 {
            let p0 = candidates[0].0.pressure;
            if candidates
                .iter()
                .skip(1)
                .any(|c| (c.0.pressure - p0).abs() > 2e-10 * p0.max(c.0.pressure))
            {
                return Err(Error::InvalidInput(
                    "multiple globally stable MBWR coexistence candidates",
                ));
            }
            candidates.sort_by(|a, b| b.2.total_cmp(&a.2));
        }
        Ok(candidates.remove(0).0)
    }

    fn solve_pair(
        &self,
        vapor: &Interval,
        liquid: &Interval,
        requested: Option<(f64, f64)>,
    ) -> Result<(Saturation, f64, f64)> {
        let mut lo = vapor.low_pressure.max(liquid.low_pressure);
        let mut hi = vapor.high_pressure.min(liquid.high_pressure);
        if let Some((a, b)) = requested {
            lo = lo.max(a);
            hi = hi.min(b);
        }
        // The dilute branch starts at P=0. Keep the auto search above its f64
        // pressure floor; explicit windows remain available at lower pressures.
        if requested.is_none() {
            lo = lo.max(64. * f64::EPSILON * R * self.prepared.temperature * self.maximum_density);
        }
        positive(lo, "positive MBWR coexistence pressure window")?;
        if hi.partial_cmp(&lo) != Some(std::cmp::Ordering::Greater) {
            return Err(Error::NotConverged(
                "overlapping MBWR branch pressure window",
            ));
        }
        let residual = |p: f64| -> Result<(f64, f64, f64)> {
            let rv = self.root(p, vapor)?;
            let rl = self.root(p, liquid)?;
            if rl <= rv {
                return Err(Error::NotConverged("separated MBWR coexistence phases"));
            }
            let av = self.prepared.residual_density(rv)?;
            let al = self.prepared.residual_density(rl)?;
            let f = (rl / rv).ln() + al.value - av.value + rl * al.drho - rv * av.drho;
            if !f.is_finite() {
                return Err(Error::InvalidInput(
                    "non-finite MBWR coexistence chemical potential",
                ));
            }
            Ok((f, rv, rl))
        };
        let (mut fl, _, _) = residual(lo)?;
        let (fh, _, _) = residual(hi)?;
        let mut a = lo.ln();
        let mut b = hi.ln();
        if fl == 0. {
            let q = residual(lo)?;
            return Ok((self.saturation_at(lo, q)?, q.1, q.2));
        }
        if fh == 0. {
            let q = residual(hi)?;
            return Ok((self.saturation_at(hi, q)?, q.1, q.2));
        }
        if fl.signum() == fh.signum() {
            return Err(Error::NotConverged(
                "MBWR coexistence chemical potential bracket",
            ));
        }
        let mut x = a + (b - a) * 0.5;
        for _ in 0..120 {
            let p = x.exp();
            let q = residual(p)?;
            let (_, rv, rl) = q;
            // Normalize tolerance by the phase separation to reject false
            // convergence as the two outer roots approach one another.
            let phase_gap = ((rl - rv) / rv).ln_1p();
            let tol = 2e-12 * phase_gap.max(1e-8);
            if q.0.abs() <= tol {
                return Ok((self.saturation_at(p, q)?, q.1, q.2));
            }
            if b - a <= 8. * f64::EPSILON * x.abs().max(1.) {
                return Err(Error::NotConverged(
                    "MBWR coexistence pressure roundoff before chemical-potential closure",
                ));
            }
            let same = q.0.signum() == fl.signum();
            if same {
                a = x;
                fl = q.0;
            } else {
                b = x;
            }
            let derivative = p * ((1. / rl) - (1. / rv)) / (R * self.prepared.temperature);
            let newton = x - q.0 / derivative;
            x = if derivative.is_finite()
                && derivative < 0.
                && newton.is_finite()
                && newton > a
                && newton < b
            {
                newton
            } else {
                a + (b - a) * 0.5
            };
        }
        Err(Error::NotConverged("MBWR saturation pressure"))
    }

    fn globally_stable_candidate(&self, p: f64, rv: f64, rl: f64) -> Result<bool> {
        let mu = |rho: f64| -> Result<f64> {
            let a = self.prepared.residual_density(rho)?;
            Ok(rho.ln() + a.value + rho * a.drho)
        };
        let muv = mu(rv)?;
        let mul = mu(rl)?;
        let reference = 0.5 * (muv + mul);
        let gap = ((rl - rv) / rv).ln_1p();
        let tol = 2e-12 * gap.max(1e-8);
        for rho in self.densities(p)? {
            if (rho - rv).abs() <= 64. * f64::EPSILON * rv.max(rho)
                || (rho - rl).abs() <= 64. * f64::EPSILON * rl.max(rho)
            {
                continue;
            }
            if mu(rho)? < reference - tol {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn saturation_at(&self, p: f64, q: (f64, f64, f64)) -> Result<Saturation> {
        let (_, rv, rl) = q;
        if rl <= rv {
            return Err(Error::NotConverged("separated MBWR coexistence phases"));
        }
        let phase_gap = ((rl - rv) / rv).ln_1p();
        if q.0.abs() > 2e-12 * phase_gap.max(1e-8) {
            return Err(Error::NotConverged(
                "MBWR coexistence chemical-potential closure",
            ));
        }
        let pv = self.prepared.pressure_density(rv)?;
        let pl = self.prepared.pressure_density(rl)?;
        let tolerance_v = 2e-10 * p + 64. * f64::EPSILON * R * self.prepared.temperature * rv;
        let tolerance_l = 2e-10 * p + 64. * f64::EPSILON * R * self.prepared.temperature * rl;
        // Permit the full pressure-sum cancellation bound only when the ideal
        // roundoff floor cannot close a density root already at machine limit.
        let closure_v = if (pv.value - p).abs() <= tolerance_v {
            tolerance_v
        } else {
            2e-10 * p + self.prepared.pressure_roundoff_bound(rv).unwrap_or(0.)
        };
        let closure_l = if (pl.value - p).abs() <= tolerance_l {
            tolerance_l
        } else {
            2e-10 * p + self.prepared.pressure_roundoff_bound(rl).unwrap_or(0.)
        };
        if !p.is_finite()
            || p <= 0.
            || (pv.value - p).abs() > closure_v
            || (pl.value - p).abs() > closure_l
            || pv.drho <= 0.
            || pl.drho <= 0.
        {
            return Err(Error::NotConverged(
                "MBWR coexistence pressure closure or mechanical stability",
            ));
        }
        let t = self.prepared.temperature;
        let av = self.prepared.residual_density(rv)?;
        let al = self.prepared.residual_density(rl)?;
        let hv = -R * t * t * av.dt;
        let hl = -R * t * t * al.dt;
        let latent = hv - hl + p * (1. / rv - 1. / rl);
        let dpdt = latent / (t * (1. / rv - 1. / rl));
        if !dpdt.is_finite() {
            return Err(Error::InvalidInput("non-finite MBWR Clapeyron slope"));
        }
        Ok(Saturation {
            temperature: t,
            pressure: p,
            liquid_volume: 1. / rl,
            vapor_volume: 1. / rv,
            dpdt,
        })
    }
}

impl Mbwr {
    /// Pure-fluid coexistence at fixed temperature; max_density bounds branch search [mol/m3].
    pub fn saturation_pressure(&self, temperature: f64, max_density: f64) -> Result<Saturation> {
        self.prepare_equilibrium(temperature, max_density)?
            .saturation_pressure()
    }
    /// Fixed-temperature coexistence constrained to an explicit pressure window in Pa.
    pub fn saturation_pressure_bounded(
        &self,
        temperature: f64,
        min_pressure: f64,
        max_pressure: f64,
        max_density: f64,
    ) -> Result<Saturation> {
        self.prepare_equilibrium(temperature, max_density)?
            .saturation_pressure_bounded(min_pressure, max_pressure)
    }
    /// Fixed-pressure inverse inside the caller's temperature bracket in K.
    /// Both endpoints and the intervening interval must lie on one valid,
    /// continuous, monotone coexistence branch; branch events are not detected.
    pub fn saturation_temperature(
        &self,
        pressure: f64,
        min_temperature: f64,
        max_temperature: f64,
        max_density: f64,
    ) -> Result<Saturation> {
        positive(pressure, "saturation pressure")?;
        positive(min_temperature, "minimum saturation temperature")?;
        positive(max_temperature, "maximum saturation temperature")?;
        if min_temperature >= max_temperature {
            return Err(Error::InvalidInput(
                "saturation temperature bounds must be increasing",
            ));
        }
        let mut lo = min_temperature;
        let mut hi = max_temperature;
        let sl = self.saturation_pressure(lo, max_density)?;
        let sh = self.saturation_pressure(hi, max_density)?;
        let pressure_tolerance = 2e-11 * pressure + 64. * f64::EPSILON * pressure;
        if (pressure - sl.pressure).abs() <= pressure_tolerance {
            return Ok(sl);
        }
        if (pressure - sh.pressure).abs() <= pressure_tolerance {
            return Ok(sh);
        }
        if pressure < sl.pressure.min(sh.pressure) || pressure > sl.pressure.max(sh.pressure) {
            return Err(Error::NotConverged(
                "MBWR inverse saturation temperature bracket",
            ));
        }
        let increasing = sh.pressure > sl.pressure;
        let mut t = lo + (hi - lo) * 0.5;
        for _ in 0..100 {
            let s = self.saturation_pressure(t, max_density)?;
            if (s.pressure - pressure).abs() <= pressure_tolerance {
                return Ok(s);
            }
            if (s.pressure < pressure) == increasing {
                lo = t;
            } else {
                hi = t;
            }
            // Use the analytic Clapeyron derivative only while its Newton step
            // stays inside the caller's bracket; otherwise bisect.
            let newton = t + (pressure - s.pressure) / s.dpdt;
            t = if s.dpdt.is_finite()
                && s.dpdt != 0.
                && newton.is_finite()
                && newton > lo
                && newton < hi
            {
                newton
            } else {
                lo + (hi - lo) * 0.5
            };
            if hi - lo <= 8. * f64::EPSILON * t {
                return Err(Error::NotConverged(
                    "MBWR inverse saturation temperature roundoff before pressure closure",
                ));
            }
        }
        Err(Error::NotConverged("MBWR saturation temperature"))
    }
}

impl crate::MbwrFluid {
    pub fn saturation_pressure(&self, temperature: f64, max_density: f64) -> Result<Saturation> {
        self.mbwr().saturation_pressure(temperature, max_density)
    }
    pub fn saturation_pressure_bounded(
        &self,
        temperature: f64,
        min_pressure: f64,
        max_pressure: f64,
        max_density: f64,
    ) -> Result<Saturation> {
        self.mbwr().saturation_pressure_bounded(
            temperature,
            min_pressure,
            max_pressure,
            max_density,
        )
    }
    pub fn saturation_temperature(
        &self,
        pressure: f64,
        min_temperature: f64,
        max_temperature: f64,
        max_density: f64,
    ) -> Result<Saturation> {
        self.mbwr()
            .saturation_temperature(pressure, min_temperature, max_temperature, max_density)
    }
}
