//! Pure-fluid coexistence in logarithmic liquid/vapor density coordinates.
use super::*;
use crate::Saturation;

struct CoexistenceEvaluation {
    residual: [f64; 2],
    jacobian: [f64; 4],
    liquid: HelmholtzJet,
    vapor: HelmholtzJet,
    liquid_delta: f64,
    vapor_delta: f64,
    vapor_pressure: f64,
    liquid_stiffness: f64,
    vapor_stiffness: f64,
    pressure_roundoff: f64,
}

impl PreparedMultiparameter<'_> {
    // Columns are log liquid density, log vapor density and log temperature.
    // The pressure specification adds one equation to the coexistence system.
    fn isobaric_equations(
        &self,
        e: &CoexistenceEvaluation,
        pressure: f64,
    ) -> Result<([f64; 3], [f64; 9])> {
        positive(e.vapor_pressure, "multiparameter inverse vapor pressure")?;
        let log_pressure =
            (e.vapor_pressure * self.eos.gas_constant() * self.temperature * self.eos.rc()
                / pressure)
                .ln();
        let pt = -e.liquid_delta * e.liquid.d_deltatau + e.vapor_delta * e.vapor.d_deltatau;
        let mt = -(e.liquid.d_tau - e.vapor.d_tau + e.liquid.d_deltatau - e.vapor.d_deltatau);
        let density_response = e.vapor_delta * e.vapor_stiffness / e.vapor_pressure;
        let temperature_response = 1. - e.vapor_delta * e.vapor.d_deltatau / e.vapor_pressure;
        Ok((
            [e.residual[0], e.residual[1], log_pressure],
            [
                e.jacobian[0],
                e.jacobian[1],
                pt,
                e.jacobian[2],
                e.jacobian[3],
                mt,
                0.,
                density_response,
                temperature_response,
            ],
        ))
    }

    fn coexistence_equations(&self, u: [f64; 2]) -> Result<CoexistenceEvaluation> {
        let dl = u[0].exp();
        let dv = u[1].exp();
        if !dl.is_finite() || !dv.is_finite() || dl <= 1. || dv >= 1. || dv <= 0. || dl > 16. {
            return Err(Error::NotConverged(
                "multiparameter coexistence density bounds",
            ));
        }
        let l = self.residual_at_delta(dl)?;
        let v = self.residual_at_delta(dv)?;
        let pl = dl * (1. + l.d_delta);
        let pv = dv * (1. + v.d_delta);
        let kl = 1. + 2. * l.d_delta + l.d_deltadelta;
        let kv = 1. + 2. * v.d_delta + v.d_deltadelta;
        let gap = u[0] - u[1];
        let mu = gap + (l.alpha - v.alpha) + (l.d_delta - v.d_delta);
        let residual = [pl - pv, mu];
        let jacobian = [dl * kl, -dv * kv, kl, -kv];
        if residual.iter().chain(&jacobian).any(|v| !v.is_finite()) {
            return Err(Error::NotConverged(
                "nonfinite multiparameter coexistence equations",
            ));
        }
        Ok(CoexistenceEvaluation {
            residual,
            jacobian,
            liquid: l,
            vapor: v,
            liquid_delta: dl,
            vapor_delta: dv,
            vapor_pressure: pv,
            liquid_stiffness: kl,
            vapor_stiffness: kv,
            pressure_roundoff: 32.
                * f64::EPSILON
                * (dl * (1. + l.d_delta.abs()) + dv * (1. + v.d_delta.abs())),
        })
    }

    /// Pure fluid coexistence at this prepared temperature. The liquid/vapor
    /// density branches straddle the model's reducing density. Ancillary
    /// estimates seed an analytic two-variable Newton solve; no TP root
    /// selection or ideal-gas reference-level subtraction is required.
    pub fn saturation_pressure(&self) -> Result<Saturation> {
        if self.temperature >= self.eos.tc() {
            return Err(Error::InvalidInput(
                "multiparameter coexistence requires T below Tc",
            ));
        }
        let gap = (1. - self.temperature / self.eos.tc()).sqrt();
        let dl = self.sat_delta_estimate(Phase::Liquid)?.max(1. + gap);
        let dv = self
            .sat_delta_estimate(Phase::Vapor)?
            .min((1. - gap).max(1e-12));
        let initial = [dl.ln(), dv.ln()];
        if let Ok(q) = self.saturation_from_log_densities(initial) {
            return Ok(q);
        }
        // Follow the stable branches from a less critical temperature if the
        // ancillary start cannot reach the narrow target coexistence pair.
        let mut t = 0.9 * self.temperature;
        let lower = self.eos.prepare(t)?;
        let width = (1. - t / self.eos.tc()).sqrt();
        let mut point = lower.saturation_from_log_densities([
            lower
                .sat_delta_estimate(Phase::Liquid)?
                .max(1. + width)
                .ln(),
            lower
                .sat_delta_estimate(Phase::Vapor)?
                .min((1. - width).max(1e-12))
                .ln(),
        ])?;
        for _ in 0..80 {
            let mut next = self
                .temperature
                .min(self.eos.tc() - 0.5 * (self.eos.tc() - t));
            let initial = [
                (1. / (self.eos.rc() * point.liquid_volume)).ln(),
                (1. / (self.eos.rc() * point.vapor_volume)).ln(),
            ];
            let mut accepted = None;
            for _ in 0..16 {
                if next <= t {
                    break;
                }
                if let Ok(q) = self
                    .eos
                    .prepare(next)?
                    .saturation_from_log_densities(initial)
                {
                    accepted = Some(q);
                    break;
                }
                next = 0.5 * (t + next);
            }
            point = accepted.ok_or(Error::NotConverged(
                "multiparameter coexistence continuation",
            ))?;
            t = point.temperature;
            if t == self.temperature {
                return Ok(point);
            }
        }
        Err(Error::NotConverged(
            "multiparameter coexistence continuation limit",
        ))
    }

    fn saturation_from_log_densities(&self, mut u: [f64; 2]) -> Result<Saturation> {
        for _ in 0..80 {
            let e = self.coexistence_equations(u)?;
            if e.liquid_stiffness > 0. && e.vapor_stiffness > 0. {
                break;
            }
            if e.liquid_stiffness <= 0. {
                u[0] += 0.05;
            }
            if e.vapor_stiffness <= 0. {
                u[1] -= 0.1;
            }
        }
        for _ in 0..100 {
            let e = self.coexistence_equations(u)?;
            let pressure_tolerance = e.pressure_roundoff + 4e-11 * e.vapor_pressure.abs();
            let chemical_tolerance = 2e-11;
            if e.residual[0].abs() < pressure_tolerance
                && e.residual[1].abs() / (u[0] - u[1]) < chemical_tolerance
            {
                if e.liquid_stiffness <= 0. || e.vapor_stiffness <= 0. || e.vapor_pressure <= 0. {
                    return Err(Error::NotConverged(
                        "unstable multiparameter coexistence branches",
                    ));
                }
                let mut barrier = false;
                for delta in [1., 0.5 * (e.liquid_delta + e.vapor_delta)]
                    .into_iter()
                    .chain((1..32).filter(|i| *i != 16).map(|i| {
                        e.vapor_delta + (e.liquid_delta - e.vapor_delta) * f64::from(i) / 32.
                    }))
                {
                    let middle = self.residual_at_delta(delta)?;
                    if 1. + 2. * middle.d_delta + middle.d_deltadelta < 0. {
                        barrier = true;
                        break;
                    }
                }
                if !barrier || u[0] - u[1] < 1e-6 {
                    return Err(Error::NotConverged(
                        "unresolved multiparameter coexistence barrier",
                    ));
                }
                let r = self.eos.gas_constant();
                let t = self.temperature;
                let rc = self.eos.rc();
                let vl = 1. / (rc * e.liquid_delta);
                let vv = 1. / (rc * e.vapor_delta);
                let p = e.vapor_pressure * r * t * rc;
                let latent = r * t * (e.vapor.d_tau - e.liquid.d_tau) + p * (vv - vl);
                let dpdt = latent / (t * (vv - vl));
                positive(p, "multiparameter saturation pressure")?;
                positive(dpdt, "multiparameter saturation slope")?;
                return Ok(Saturation {
                    temperature: t,
                    pressure: p,
                    liquid_volume: vl,
                    vapor_volume: vv,
                    dpdt,
                });
            }
            let [a, b, c, d] = e.jacobian;
            let det = a * d - b * c;
            if !det.is_finite() || det == 0. {
                return Err(Error::NotConverged(
                    "singular multiparameter coexistence Jacobian",
                ));
            }
            let mut step = [
                (-d * e.residual[0] + b * e.residual[1]) / det,
                (c * e.residual[0] - a * e.residual[1]) / det,
            ];
            let limit = (step[0].abs().max(step[1].abs()) / 0.5).max(1.);
            step.iter_mut().for_each(|v| *v /= limit);
            let norm = e.residual[0].abs().max(e.residual[1].abs());
            let mut fraction = 1.;
            let mut accepted = None;
            for _ in 0..20 {
                let trial = [u[0] + fraction * step[0], u[1] + fraction * step[1]];
                if let Ok(q) = self.coexistence_equations(trial)
                    && q.liquid_stiffness > 0.
                    && q.vapor_stiffness > 0.
                    && q.residual[0].abs().max(q.residual[1].abs()) < norm
                {
                    accepted = Some(trial);
                    break;
                }
                fraction *= 0.5;
            }
            u = accepted.ok_or(Error::NotConverged(
                "multiparameter coexistence backtracking",
            ))?;
        }
        Err(Error::NotConverged("multiparameter coexistence iteration"))
    }
}

impl Multiparameter {
    /// Subcritical pure-fluid coexistence pressure, phase volumes and
    /// Clapeyron slope. This is fluid coexistence, not solid-phase stability.
    /// Stable phase branches must straddle the model's reducing density, with
    /// a resolved nonconvex region between them. Critical endpoints and
    /// unresolved or unsupported branches return an error.
    ///
    /// ```
    /// use thermopack_rs::Multiparameter;
    /// let propane = Multiparameter::new("C3")?;
    /// let point = propane.saturation_pressure(280.)?;
    /// assert!(point.liquid_volume < point.vapor_volume);
    /// let inverse = propane.saturation_temperature(point.pressure)?;
    /// assert!((inverse.temperature - 280.).abs() < 1e-5);
    /// # Ok::<(), thermopack_rs::Error>(())
    /// ```
    pub fn saturation_pressure(&self, temperature: f64) -> Result<Saturation> {
        self.prepare(temperature)?.saturation_pressure()
    }

    pub fn saturation_temperature(&self, pressure: f64) -> Result<Saturation> {
        self.saturation_temperature_with_guess(pressure, 0.8 * self.tc())
    }

    /// Invert the native coexistence curve using an analytic coupled solve for
    /// temperature and both densities. A bounded Clapeyron temperature iteration
    /// remains the fallback if coupled correction cannot resolve a valid pair.
    pub fn saturation_temperature_with_guess(
        &self,
        pressure: f64,
        temperature: f64,
    ) -> Result<Saturation> {
        positive(pressure, "multiparameter coexistence pressure")?;
        positive(temperature, "multiparameter coexistence temperature guess")?;
        let critical_pressure = self.prepare(self.tc())?.pressure(1. / self.rc())?.value;
        if pressure >= critical_pressure || temperature >= self.tc() {
            return Err(Error::InvalidInput(
                "multiparameter inverse coexistence below critical point required",
            ));
        }
        let point = self.saturation_pressure(temperature)?;
        if let Ok(q) = self.isobaric_coexistence(pressure, point) {
            return Ok(q);
        }
        self.nested_saturation_temperature(pressure, point)
    }

    fn isobaric_coexistence(&self, pressure: f64, point: Saturation) -> Result<Saturation> {
        let mut u = [
            (1. / (self.rc() * point.liquid_volume)).ln(),
            (1. / (self.rc() * point.vapor_volume)).ln(),
            point.temperature.ln(),
        ];
        for _ in 0..100 {
            let prepared = self.prepare(u[2].exp())?;
            let e = prepared.coexistence_equations([u[0], u[1]])?;
            let (residual, mut matrix) = prepared.isobaric_equations(&e, pressure)?;
            if residual[2].abs() < 2e-10
                && residual[0].abs() < e.pressure_roundoff + 4e-11 * e.vapor_pressure.abs()
                && residual[1].abs() / (u[0] - u[1]) < 2e-11
            {
                return prepared.saturation_from_log_densities([u[0], u[1]]);
            }
            let mut step = residual.map(|v| -v);
            crate::numerics::solve_linear_in_place(&mut matrix, &mut step)?;
            let limit = (step[0].abs().max(step[1].abs()) / 0.5)
                .max(step[2].abs() / 0.15)
                .max(1.);
            step.iter_mut().for_each(|v| *v /= limit);
            let norm = residual.into_iter().map(f64::abs).fold(0., f64::max);
            let mut fraction = 1.;
            let mut accepted = None;
            for _ in 0..24 {
                let trial = std::array::from_fn(|i| u[i] + fraction * step[i]);
                let t = trial[2].exp();
                if t > self.tc() * 1e-6
                    && t < self.tc()
                    && let Ok(p) = self.prepare(t)
                    && let Ok(e) = p.coexistence_equations([trial[0], trial[1]])
                    && e.liquid_stiffness > 0.
                    && e.vapor_stiffness > 0.
                    && let Ok((r, _)) = p.isobaric_equations(&e, pressure)
                    && r.into_iter().map(f64::abs).fold(0., f64::max) < norm
                {
                    accepted = Some(trial);
                    break;
                }
                fraction *= 0.5;
            }
            u = accepted.ok_or(Error::NotConverged(
                "multiparameter isobaric coexistence backtracking",
            ))?;
        }
        Err(Error::NotConverged(
            "multiparameter isobaric coexistence iteration",
        ))
    }

    fn nested_saturation_temperature(
        &self,
        pressure: f64,
        mut point: Saturation,
    ) -> Result<Saturation> {
        let (mut lo, mut hi) = (self.tc() * 1e-6, self.tc() * (1. - 1e-10));
        for _ in 0..100 {
            let f = (point.pressure / pressure).ln();
            if f.abs() < 2e-10 {
                return Ok(point);
            }
            if f > 0. {
                hi = point.temperature;
            } else {
                lo = point.temperature;
            }
            let newton = point.temperature - f * point.pressure / point.dpdt;
            let mut next = if newton > lo && newton < hi {
                newton
            } else {
                0.5 * (lo + hi)
            };
            let mut accepted = None;
            for _ in 0..20 {
                if let Ok(q) = self.saturation_pressure(next) {
                    accepted = Some(q);
                    break;
                }
                next = 0.5 * (next + point.temperature);
            }
            point = accepted.ok_or(Error::NotConverged(
                "multiparameter inverse coexistence backtracking",
            ))?;
        }
        Err(Error::NotConverged(
            "multiparameter inverse coexistence iteration",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn coupled_coexistence_jacobian_matches_independent_differences() {
        for model in [
            Multiparameter::new("C3").unwrap(),
            Multiparameter::new("R134A").unwrap(),
            Multiparameter::new("P-H2").unwrap(),
            Multiparameter::lennard_jones("NE", LennardJonesKind::Full).unwrap(),
            Multiparameter::reduced_lennard_jones_truncated_shifted(),
        ] {
            let s = model.prepare(0.8 * model.tc()).unwrap();
            let u = [2_f64.ln(), 0.05_f64.ln()];
            let e = s.coexistence_equations(u).unwrap();
            for j in 0..2 {
                let mut up = u;
                let mut um = u;
                let h = 2e-6;
                up[j] += h;
                um[j] -= h;
                let a = s.coexistence_equations(up).unwrap();
                let b = s.coexistence_equations(um).unwrap();
                for i in 0..2 {
                    let finite = (a.residual[i] - b.residual[i]) / (2. * h);
                    let actual = e.jacobian[2 * i + j];
                    assert!(
                        (finite - actual).abs() < 2e-7 + 2e-6 * finite.abs().max(actual.abs()),
                        "{} ({i},{j}): {finite} vs {actual}",
                        model.record.ident
                    );
                }
            }
        }
    }

    #[test]
    fn isobaric_coexistence_jacobian_matches_off_root_differences() {
        for model in [
            Multiparameter::new("C3").unwrap(),
            Multiparameter::new("R134A").unwrap(),
            Multiparameter::new("N-H2").unwrap(),
            Multiparameter::lennard_jones("KR", LennardJonesKind::Full).unwrap(),
            Multiparameter::reduced_lennard_jones_truncated_shifted(),
        ] {
            let u = [2_f64.ln(), 0.05_f64.ln(), (0.8 * model.tc()).ln()];
            let pressure = 0.2 * model.gas_constant() * model.tc() * model.rc();
            let evaluate = |u: [f64; 3]| {
                let p = model.prepare(u[2].exp()).unwrap();
                let e = p.coexistence_equations([u[0], u[1]]).unwrap();
                p.isobaric_equations(&e, pressure).unwrap()
            };
            let (_, jacobian) = evaluate(u);
            for j in 0..3 {
                let (mut up, mut um) = (u, u);
                let h = 2e-6;
                up[j] += h;
                um[j] -= h;
                let (a, _) = evaluate(up);
                let (b, _) = evaluate(um);
                for i in 0..3 {
                    let finite = (a[i] - b[i]) / (2. * h);
                    assert!(
                        (finite - jacobian[3 * i + j]).abs() < 3e-7 + 3e-6 * finite.abs(),
                        "{} ({i},{j}): {finite} vs {}",
                        model.record.ident,
                        jacobian[3 * i + j]
                    );
                }
            }
        }
    }
}
