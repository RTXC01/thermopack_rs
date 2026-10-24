//! Analytic density derivatives and local pure-fluid critical equations.
use super::*;
use crate::PureCriticalPoint;

#[path = "multiparameter_intervals.rs"]
mod intervals;
pub use intervals::{MultiparameterEquilibrium, MultiparameterSpinodalBracket};

const COEFFICIENTS: usize = 40;

/// Pressure and analytic density/temperature derivatives in the model's units.
/// Density is mol/m3 for physical models and the original reduced density for
/// reduced LJTS. Temperature derivatives hold density fixed.
#[derive(Debug, Clone, Copy)]
pub struct MultiparameterPressureJet {
    pub value: f64,
    pub drho: f64,
    pub drhorho: f64,
    pub drhorhorho: f64,
    pub dt: f64,
    pub dtrho: f64,
    pub dtrhorho: f64,
}

#[derive(Debug, Clone, Copy)]
enum DensityExponent {
    Plain,
    Power(i32),
    Gaussian {
        eta: f64,
        epsilon: f64,
        temperature: f64,
    },
}

#[derive(Debug, Clone)]
struct DensityTerm {
    density_power: usize,
    polynomial: [[f64; COEFFICIENTS]; 4],
    lengths: [usize; 4],
    exponent: DensityExponent,
    coefficient: f64,
    temperature_factor: f64,
    family: usize,
}

impl DensityTerm {
    fn new(
        d: i32,
        exponent: DensityExponent,
        coefficient: f64,
        temperature_factor: f64,
        family: usize,
    ) -> Result<Self> {
        if d < 0 {
            return Err(Error::InvalidInput("negative multiparameter density power"));
        }
        let mut q = [0.; COEFFICIENTS];
        let degree = match exponent {
            DensityExponent::Plain => 0,
            DensityExponent::Power(l) => {
                if l <= 0 || l as usize >= COEFFICIENTS {
                    return Err(Error::InvalidInput("multiparameter damping degree"));
                }
                q[l as usize - 1] = -f64::from(l);
                l as usize - 1
            }
            DensityExponent::Gaussian { eta, epsilon, .. } => {
                q[0] = 2. * eta * epsilon;
                q[1] = -2. * eta;
                1
            }
        };
        let d = d as usize;
        if d + 2 + degree + 3 * degree >= COEFFICIENTS {
            return Err(Error::InvalidInput(
                "multiparameter pressure polynomial capacity",
            ));
        }
        let mut polynomial = [[0.; COEFFICIENTS]; 4];
        polynomial[0][d + 1] = d as f64;
        for i in 0..=degree {
            polynomial[0][d + 2 + i] += q[i];
        }
        let mut lengths = [d + 3 + degree; 4];
        for k in 1..4 {
            lengths[k] = (lengths[k - 1] + degree).min(COEFFICIENTS);
            for i in 0..lengths[k - 1] {
                if i > 0 {
                    polynomial[k][i - 1] += i as f64 * polynomial[k - 1][i];
                }
                for (j, &v) in q[..=degree].iter().enumerate() {
                    polynomial[k][i + j] += v * polynomial[k - 1][i];
                }
            }
        }
        if polynomial.iter().flatten().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput(
                "multiparameter pressure coefficient overflow",
            ));
        }
        Ok(Self {
            density_power: d,
            polynomial,
            lengths,
            exponent,
            coefficient,
            temperature_factor,
            family,
        })
    }

    fn evaluate(&self, delta: f64) -> [f64; 4] {
        let exponent = match self.exponent {
            DensityExponent::Plain => 0.,
            DensityExponent::Power(l) => -delta.powi(l),
            DensityExponent::Gaussian {
                eta,
                epsilon,
                temperature,
            } => -eta * (delta - epsilon).powi(2) - temperature,
        };
        let factor = self.coefficient * exponent.exp();
        std::array::from_fn(|k| {
            self.polynomial[k][..self.lengths[k]]
                .iter()
                .rev()
                .fold(0., |v, a| v * delta + a)
                * factor
        })
    }
}

impl PreparedMultiparameter<'_> {
    /// Absolute pressure allowance from cancellation of the individual EOS
    /// terms. This includes a conservative 256-epsilon arithmetic/exp guard.
    /// It can exceed the pressure itself in extreme extrapolations; it is not
    /// a fitted-model accuracy estimate. Ordinary queries need no such relaxation.
    pub fn pressure_roundoff(&self, rho: f64) -> Result<f64> {
        if !rho.is_finite() || rho < 0. {
            return Err(Error::InvalidInput("nonnegative finite density required"));
        }
        let delta = rho / self.eos.rc();
        let mut scale = 1.;
        for (i, p) in self.eos.record.pol_terms.iter().enumerate() {
            scale += f64::from(p.d).abs() * self.tau_cache[i].abs() * delta.powi(p.d);
        }
        for (i, p) in self.eos.record.exp_terms.iter().enumerate() {
            let dl = delta.powi(p.l);
            scale += (f64::from(p.d).abs() + f64::from(p.l) * dl)
                * self.tau_exp_cache[i].abs()
                * delta.powi(p.d)
                * (-dl).exp();
        }
        for (i, p) in self.eos.record.expexp_terms.iter().enumerate() {
            let term = self.tau_ee_cache[i].abs()
                * delta.powi(p.d)
                * (-p.eta * (delta - p.epsilon).powi(2) - self.gaussian_temperature[i]).exp();
            scale += (f64::from(p.d).abs() + 2. * p.eta.abs() * delta * (delta + p.epsilon.abs()))
                * term;
        }
        let allowance =
            256. * f64::EPSILON * rho * self.eos.gas_constant() * self.temperature * scale;
        if !allowance.is_finite() {
            return Err(Error::InvalidInput(
                "nonfinite pressure cancellation allowance",
            ));
        }
        Ok(allowance)
    }

    fn visit_density_terms(&self, mut visit: impl FnMut(DensityTerm) -> Result<()>) -> Result<()> {
        for (i, t) in self.eos.record.pol_terms.iter().enumerate() {
            visit(DensityTerm::new(
                t.d,
                DensityExponent::Plain,
                self.tau_cache[i],
                1. - t.t,
                0,
            )?)?;
        }
        for (i, t) in self.eos.record.exp_terms.iter().enumerate() {
            visit(DensityTerm::new(
                t.d,
                DensityExponent::Power(t.l),
                self.tau_exp_cache[i],
                1. - t.t,
                1,
            )?)?;
        }
        for (i, t) in self.eos.record.expexp_terms.iter().enumerate() {
            visit(DensityTerm::new(
                t.d,
                DensityExponent::Gaussian {
                    eta: t.eta,
                    epsilon: t.epsilon,
                    temperature: self.gaussian_temperature[i],
                },
                self.tau_ee_cache[i],
                1. - t.t + self.beta_tau_diff[i],
                2,
            )?)?;
        }
        Ok(())
    }

    /// Analytic pressure derivatives through third density order and mixed
    /// T/rho/rho order. This query uses fixed-size polynomial storage, including
    /// the exact zero-density limit; it does not use numerical differences.
    pub fn pressure_density_jet(&self, rho: f64) -> Result<MultiparameterPressureJet> {
        if !rho.is_finite() || rho < 0. {
            return Err(Error::InvalidInput(
                "nonnegative finite molar density required",
            ));
        }
        let delta = rho / self.eos.rc();
        let mut family = [[0.; 7]; 3];
        self.visit_density_terms(|term| {
            let q = term.evaluate(delta);
            for (sum, &value) in family[term.family][..4].iter_mut().zip(&q) {
                *sum += value;
            }
            for (sum, &value) in family[term.family][4..].iter_mut().zip(&q) {
                *sum += term.temperature_factor * value;
            }
            Ok(())
        })?;
        let total: [f64; 7] = std::array::from_fn(|i| (family[0][i] + family[1][i]) + family[2][i]);
        let r = self.eos.gas_constant();
        let rt = r * self.temperature;
        let rc = self.eos.rc();
        // Preserve the established mp_pressure summation for its original
        // value/first derivatives; the additional derivatives use the analytic
        // exponential-polynomial recurrence above.
        let a = self.residual_at_delta(delta)?;
        let jet = MultiparameterPressureJet {
            value: rho * rt * (1. + a.d_delta),
            drho: rt * (1. + 2. * a.d_delta + a.d_deltadelta),
            drhorho: rt / rc * total[2],
            drhorhorho: rt / (rc * rc) * total[3],
            dt: rho * r * (1. + a.d_delta - a.d_deltatau),
            dtrho: r * (1. + total[5]),
            dtrhorho: r / rc * total[6],
        };
        if [
            jet.value,
            jet.drho,
            jet.drhorho,
            jet.drhorhorho,
            jet.dt,
            jet.dtrho,
            jet.dtrhorho,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(Error::InvalidInput(
                "nonfinite multiparameter pressure derivatives",
            ));
        }
        Ok(jet)
    }
}

impl Multiparameter {
    /// Refine the local fluid critical point from stored reducing T and density.
    /// This solves P_rho=P_rhorho=0, rather than returning table metadata.
    /// If that guess fails or reaches a nonordinary stationary point, guesses
    /// at 0.95 and 1.05 reducing density are tried and the higher-temperature
    /// ordinary result is selected. This is a local candidate policy.
    ///
    /// ```
    /// use thermopack_rs::Multiparameter;
    /// let model = Multiparameter::new("R134A")?;
    /// let critical = model.critical_point()?;
    /// let jet = model.prepare(critical.temperature)?
    ///     .pressure_density_jet(1.0 / critical.volume)?;
    /// assert!(jet.drho.abs() < 1e-6);
    /// # Ok::<(), thermopack_rs::Error>(())
    /// ```
    pub fn critical_point(&self) -> Result<PureCriticalPoint> {
        if let Ok(point) = self.critical_point_from(self.tc(), self.rc()) {
            return Ok(point);
        }
        // The full LJ fit has a nonordinary stationary point near its stored
        // reducing density and two neighboring ordinary candidates. Try both
        // sides and choose the higher-temperature local ordinary candidate.
        let mut best: Option<PureCriticalPoint> = None;
        for ratio in [0.95, 1.05] {
            if let Ok(point) = self.critical_point_from(self.tc(), ratio * self.rc())
                && best.is_none_or(|previous| point.temperature > previous.temperature)
            {
                best = Some(point);
            }
        }
        best.ok_or(Error::NotConverged(
            "multiparameter local ordinary critical point",
        ))
    }

    /// Local critical solve from explicit T and molar density. A positive
    /// pressure and ordinary positive third density derivative are required.
    /// This does not select globally among critical points of an extrapolated EOS.
    pub fn critical_point_from(&self, mut t: f64, mut rho: f64) -> Result<PureCriticalPoint> {
        positive(t, "critical initial temperature")?;
        positive(rho, "critical initial density")?;
        for _ in 0..80 {
            let q = self.prepare(t)?.pressure_density_jet(rho)?;
            let scale = self.gas_constant() * t;
            let norm = (q.drho / scale).hypot(rho * q.drhorho / scale);
            if norm < 2e-11 {
                if q.value <= 0. || q.drhorhorho <= 0. {
                    return Err(Error::NotConverged(
                        "nonordinary multiparameter stationary pressure point",
                    ));
                }
                return Ok(PureCriticalPoint {
                    temperature: t,
                    pressure: q.value,
                    volume: 1. / rho,
                    compressibility: q.value / (rho * scale),
                });
            }
            let [a, b, c, d] = [
                t * q.dtrho,
                rho * q.drhorho,
                t * q.dtrhorho,
                rho * q.drhorhorho,
            ];
            let det = a * d - b * c;
            if !det.is_finite() || det == 0. {
                break;
            }
            let st = (-q.drho * d + b * q.drhorho) / det;
            let sr = (-a * q.drhorho + q.drho * c) / det;
            let mut length = 1_f64.min(0.3 / st.abs().max(sr.abs()));
            let mut accepted = None;
            for _ in 0..24 {
                let nt = t * (length * st).exp();
                let nr = rho * (length * sr).exp();
                if let Ok(p) = self.prepare(nt)
                    && let Ok(nq) = p.pressure_density_jet(nr)
                    && (nq.drho / (self.gas_constant() * nt))
                        .hypot(nr * nq.drhorho / (self.gas_constant() * nt))
                        < norm
                {
                    accepted = Some((nt, nr));
                    break;
                }
                length *= 0.5;
            }
            (t, rho) =
                accepted.ok_or(Error::NotConverged("multiparameter critical backtracking"))?;
        }
        Err(Error::NotConverged("multiparameter critical equations"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn models() -> Vec<Multiparameter> {
        let mut m: Vec<_> = ["C3", "R134A", "N-H2", "O-H2", "P-H2"]
            .iter()
            .map(|c| Multiparameter::new(c).unwrap())
            .collect();
        for kind in [LennardJonesKind::Full, LennardJonesKind::TruncatedShifted] {
            m.push(Multiparameter::lennard_jones("AR", kind).unwrap());
        }
        m.push(Multiparameter::reduced_lennard_jones_truncated_shifted());
        m
    }
    #[test]
    fn analytic_density_jet_matches_pressure_and_independent_differences() {
        for m in models() {
            for tr in [0.7, 1., 1.3, 2.] {
                let t = tr * m.tc();
                let p = m.prepare(t).unwrap();
                for delta in [0., 1e-5, 0.2, 0.9, 1., 1.2, 2., 4.] {
                    let rho = delta * m.rc();
                    let q = p.pressure_density_jet(rho).unwrap();
                    let base = p.pressure_density(rho).unwrap();
                    assert_eq!(q.value, base.0);
                    assert_eq!(q.drho, base.1);
                    let ht = t * 1e-5;
                    let a = m
                        .prepare(t + ht)
                        .unwrap()
                        .pressure_density_jet(rho)
                        .unwrap();
                    let b = m
                        .prepare(t - ht)
                        .unwrap()
                        .pressure_density_jet(rho)
                        .unwrap();
                    for (actual, finite, scale) in [
                        (
                            q.dt,
                            (a.value - b.value) / (2. * ht),
                            m.gas_constant() * m.rc(),
                        ),
                        (q.dtrho, (a.drho - b.drho) / (2. * ht), m.gas_constant()),
                        (
                            q.dtrhorho,
                            (a.drhorho - b.drhorho) / (2. * ht),
                            m.gas_constant() / m.rc(),
                        ),
                    ] {
                        assert!(
                            (actual - finite).abs()
                                < 2e-7 * scale + 2e-6 * actual.abs().max(finite.abs()),
                            "{} {tr} {delta}: {actual} != {finite}",
                            m.record.ident
                        );
                    }
                    if rho > 0. {
                        let h = (1e-5 * rho).max(1e-6 * m.rc());
                        let a = p.pressure_density_jet(rho + h).unwrap();
                        let b = p.pressure_density_jet(rho - h).unwrap();
                        for (actual, finite, scale) in [
                            (
                                q.drhorho,
                                (a.drho - b.drho) / (2. * h),
                                m.gas_constant() * t / m.rc(),
                            ),
                            (
                                q.drhorhorho,
                                (a.drhorho - b.drhorho) / (2. * h),
                                m.gas_constant() * t / m.rc().powi(2),
                            ),
                        ] {
                            assert!(
                                (actual - finite).abs()
                                    < 3e-7 * scale + 3e-6 * actual.abs().max(finite.abs()),
                                "{} {tr} {delta}: {actual} != {finite}",
                                m.record.ident
                            );
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn local_critical_points_close_and_have_ordinary_curvature() {
        for m in models() {
            let c = m.critical_point().unwrap_or_else(|e| {
                panic!(
                    "{} {:?}: {e}; {:?}",
                    m.record.ident,
                    m.lj,
                    m.prepare(m.tc())
                        .unwrap()
                        .pressure_density_jet(m.rc())
                        .unwrap()
                )
            });
            let rho = 1. / c.volume;
            let q = m
                .prepare(c.temperature)
                .unwrap()
                .pressure_density_jet(rho)
                .unwrap();
            assert!(q.drho.abs() < 3e-11 * m.gas_constant() * c.temperature);
            assert!((rho * q.drhorho).abs() < 3e-11 * m.gas_constant() * c.temperature);
            assert!(q.drhorhorho > 0. && q.value > 0.);
            assert!((c.temperature / m.tc() - 1.).abs() < 0.02);
        }
    }

    #[test]
    fn full_lj_distinguishes_neighboring_ordinary_critical_candidates() {
        let m = Multiparameter::lennard_jones("AR", LennardJonesKind::Full).unwrap();
        assert!(m.critical_point_from(m.tc(), m.rc()).is_err());
        let low = m.critical_point_from(m.tc(), 0.95 * m.rc()).unwrap();
        let high = m.critical_point_from(m.tc(), 1.05 * m.rc()).unwrap();
        assert!(low.temperature < high.temperature && low.volume > high.volume);
        assert!((m.critical_point().unwrap().temperature - high.temperature).abs() < 1e-8);
    }
}
