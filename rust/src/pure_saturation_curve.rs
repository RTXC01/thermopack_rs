//! Pure-fluid saturation grids, including the critical endpoint.
//!
//! The source `saturation_curve::singleCompSaturation` distributes either
//! temperature or pressure to the database critical point. The physical API
//! instead terminates at the critical point of the active EOS.
//!
//! ```
//! use thermopack_rs::{Cubic, CubicKind};
//! use thermopack_rs::pure_saturation_curve::{PureSaturationOptions, PureSaturationEnd};
//! let model = Cubic::new(CubicKind::PengRobinson, &["C1"])?;
//! let curve = model.pure_saturation_curve(0, PureSaturationOptions::default())?;
//! assert_eq!(curve.end, PureSaturationEnd::CriticalPoint);
//! assert!(curve.failure.is_none());
//! # Ok::<(), thermopack_rs::Error>(())
//! ```
use crate::{Cubic, Error, Phase, R, Result, Saturation, error::positive};

#[derive(Debug, Clone, Copy)]
pub enum PureSaturationStart {
    Temperature(f64),
    Pressure(f64),
}

#[derive(Debug, Clone, Copy)]
pub struct PureSaturationOptions {
    pub start: PureSaturationStart,
    /// Requested linear pressure spacing (Pa). The point limit can enlarge it.
    /// Ignored by the logarithmic grid, whose requested log(P) spacing is 0.02.
    pub maximum_pressure_step: f64,
    /// The grid is uniform in the specified coordinate, or its logarithm.
    pub logarithmic: bool,
    /// Resamples the entire interval, rather than returning a bounded prefix.
    pub maximum_points: usize,
}
impl Default for PureSaturationOptions {
    fn default() -> Self {
        Self {
            start: PureSaturationStart::Pressure(1e5),
            maximum_pressure_step: 2e4,
            logarithmic: false,
            maximum_points: 100,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PureSaturationEnd {
    CriticalPoint,
    SourceCriticalParameters,
    NoConvergence,
}

#[derive(Debug, Clone)]
pub struct PureSaturationCurve {
    pub points: Vec<Saturation>,
    /// Zero-based index in the original model.
    pub component: usize,
    pub end: PureSaturationEnd,
    /// An accepted prefix is retained if an interior corrector fails.
    pub failure: Option<Error>,
}

struct Evaluation {
    point: Saturation,
    residual: f64,
    dt: f64,
    dp: f64,
}

struct Problem<'a> {
    model: &'a Cubic,
    z: Vec<f64>,
    component: usize,
    baseline: bool,
    phi_l: Vec<f64>,
    phi_v: Vec<f64>,
}
impl<'a> Problem<'a> {
    fn new(model: &'a Cubic, component: usize, baseline: bool) -> Result<Self> {
        if component >= model.components().len() {
            return Err(Error::InvalidInput("pure saturation component index"));
        }
        let mut z = vec![0.; model.components().len()];
        z[component] = 1.;
        let n = z.len();
        Ok(Self {
            model,
            z,
            component,
            baseline,
            phi_l: vec![0.; n],
            phi_v: vec![0.; n],
        })
    }
    fn evaluate(&mut self, t: f64, p: f64) -> Result<Evaluation> {
        let state = self.model.prepare(t, &self.z)?;
        let (vl, vv) = if self.baseline {
            (
                state.baseline_volume(p, Phase::Liquid)?,
                state.baseline_volume(p, Phase::Vapor)?,
            )
        } else {
            let (roots, count) = state.volumes(p)?;
            if count < 2 || roots[count - 1] - roots[0] <= 1e-8 * roots[count - 1] {
                return Err(Error::NotConverged("distinct pure saturation phases"));
            }
            (roots[0], roots[count - 1])
        };
        state.ln_fugacity_tp_into(p, vl, &mut self.phi_l)?;
        state.ln_fugacity_tp_into(p, vv, &mut self.phi_v)?;
        let hl = state.residual(vl)?.enthalpy;
        let hv = state.residual(vv)?.enthalpy;
        let dt = (hv - hl) / (R * t * t);
        let dp = (vl - vv) / (R * t);
        let dpdt = if dp == 0. {
            state.pressure(vl)?.dt
        } else {
            -dt / dp
        };
        let residual = self.phi_l[self.component] - self.phi_v[self.component];
        if [residual, dt, dp, dpdt].iter().any(|v| !v.is_finite()) {
            return Err(Error::NotConverged("nonfinite pure saturation equations"));
        }
        Ok(Evaluation {
            point: Saturation {
                temperature: t,
                pressure: p,
                liquid_volume: vl,
                vapor_volume: vv,
                dpdt,
            },
            residual,
            dt,
            dp,
        })
    }
    fn correct(
        &mut self,
        fixed: PureSaturationStart,
        guess: f64,
        low: f64,
        high: f64,
    ) -> Result<Saturation> {
        let mut q = guess.clamp(low, high).ln();
        let (low, high) = (low.ln(), high.ln());
        let tolerance = if self.baseline { 1e-10 } else { 1e-12 };
        for _ in 0..100 {
            let (t, p) = conditions(fixed, q.exp());
            let current = self.evaluate(t, p)?;
            // A coarse predictor can leave the two-root region. Equal phase
            // roots give a trivial zero fugacity difference, not coexistence.
            if current.point.vapor_volume - current.point.liquid_volume
                <= 1e-8 * current.point.vapor_volume
            {
                return Err(Error::NotConverged("coalesced pure curve trial"));
            }
            if current.residual.abs() < tolerance {
                return Ok(current.point);
            }
            let derivative = q.exp()
                * match fixed {
                    PureSaturationStart::Temperature(_) => current.dp,
                    PureSaturationStart::Pressure(_) => current.dt,
                };
            if derivative == 0. || !derivative.is_finite() {
                return Err(Error::NotConverged("pure saturation curve Jacobian"));
            }
            let step = (-current.residual / derivative).clamp(-0.2, 0.2);
            let mut fraction = 1.;
            let mut accepted = None;
            for _ in 0..20 {
                let next = (q + fraction * step).clamp(low, high);
                let (t, p) = conditions(fixed, next.exp());
                if let Ok(trial) = self.evaluate(t, p)
                    && trial.residual.abs() < current.residual.abs()
                {
                    accepted = Some(next);
                    break;
                }
                fraction *= 0.5;
            }
            q = accepted.ok_or(Error::NotConverged("pure saturation curve line search"))?;
        }
        Err(Error::NotConverged("pure saturation curve corrector"))
    }
    fn cold(&self, fixed: PureSaturationStart) -> Result<Saturation> {
        match fixed {
            PureSaturationStart::Temperature(t) => {
                self.model.pseudo_saturation_pressure(t, &self.z)
            }
            PureSaturationStart::Pressure(p) => {
                self.model.pseudo_saturation_temperature(p, &self.z)
            }
        }
    }
    fn endpoint(&mut self, endpoint: crate::critical::CriticalPoint) -> Result<Saturation> {
        if self.baseline {
            return self
                .evaluate(endpoint.temperature, endpoint.pressure)
                .map(|p| p.point);
        }
        let state = self.model.prepare(endpoint.temperature, &self.z)?;
        Ok(Saturation {
            temperature: endpoint.temperature,
            pressure: endpoint.pressure,
            liquid_volume: endpoint.volume,
            vapor_volume: endpoint.volume,
            dpdt: state.pressure(endpoint.volume)?.dt,
        })
    }
}

fn conditions(fixed: PureSaturationStart, variable: f64) -> (f64, f64) {
    match fixed {
        PureSaturationStart::Temperature(t) => (t, variable),
        PureSaturationStart::Pressure(p) => (variable, p),
    }
}

impl Cubic {
    /// Trace a selected pure component, retaining the original model and its
    /// parameters. The final node is the critical point of the active EOS.
    /// A start at that point, within 32 machine epsilons in the specified
    /// coordinate, returns one critical node.
    pub fn pure_saturation_curve(
        &self,
        component: usize,
        options: PureSaturationOptions,
    ) -> Result<PureSaturationCurve> {
        self.trace_pure_saturation_curve(component, options, false)
    }
    /// Retain the fork's grid, database-critical endpoint and raw cubic volume
    /// convention. A database endpoint need not be critical for modified EOS
    /// parameters. Raw iteration identity requires separate qualification.
    /// A one-point buffer retains only the endpoint. Critical-temperature
    /// starts retain the source's coarse metastable initialization estimate.
    pub fn baseline_pure_saturation_curve(
        &self,
        component: usize,
        options: PureSaturationOptions,
    ) -> Result<PureSaturationCurve> {
        self.trace_pure_saturation_curve(component, options, true)
    }
    fn trace_pure_saturation_curve(
        &self,
        component: usize,
        options: PureSaturationOptions,
        baseline: bool,
    ) -> Result<PureSaturationCurve> {
        if options.maximum_points == 0 || !baseline && options.maximum_points < 2 {
            return Err(Error::InvalidInput(
                "pure saturation grid needs at least two points",
            ));
        }
        if !options.logarithmic {
            positive(
                options.maximum_pressure_step,
                "pure saturation pressure spacing",
            )?;
        }
        let initial = match options.start {
            PureSaturationStart::Temperature(t) => t,
            PureSaturationStart::Pressure(p) => p,
        };
        positive(initial, "pure saturation start")?;
        let mut problem = Problem::new(self, component, baseline)?;
        let endpoint = if baseline {
            let record = &self.components()[component];
            crate::critical::CriticalPoint {
                temperature: record.critical_temperature,
                pressure: record.critical_pressure,
                volume: f64::NAN,
            }
        } else {
            crate::critical::critical_point(self, &problem.z)?
        };
        if !baseline {
            let target = match options.start {
                PureSaturationStart::Temperature(_) => endpoint.temperature,
                PureSaturationStart::Pressure(_) => endpoint.pressure,
            };
            if (initial - target).abs() <= 32. * f64::EPSILON * target {
                return Ok(PureSaturationCurve {
                    points: vec![problem.endpoint(endpoint)?],
                    component,
                    end: PureSaturationEnd::CriticalPoint,
                    failure: None,
                });
            }
        }
        // At a critical start the source can accept its coarse metastable
        // PureSat estimate once the ordinary phase roots have coalesced.
        let critical_temperature_start = baseline
            && matches!(options.start,PureSaturationStart::Temperature(t) if t==endpoint.temperature);
        let critical_pressure_start = baseline
            && matches!(options.start,PureSaturationStart::Pressure(p) if p==endpoint.pressure);
        let start = if critical_temperature_start {
            let pressure = self
                .baseline_pure_critical_pressure_seed(component)?
                .ln()
                .exp();
            problem.evaluate(endpoint.temperature, pressure)?.point
        } else if critical_pressure_start {
            problem
                .evaluate(endpoint.temperature.ln().exp(), endpoint.pressure)?
                .point
        } else {
            problem.cold(options.start)?
        };
        if !critical_temperature_start
            && !critical_pressure_start
            && (start.temperature >= endpoint.temperature || start.pressure >= endpoint.pressure)
        {
            return Err(Error::InvalidInput(
                "pure saturation start below critical endpoint required",
            ));
        }
        let requested_intervals = if options.logarithmic {
            (endpoint.pressure.ln() - start.pressure.ln()) / 0.02
        } else {
            (endpoint.pressure - start.pressure) / options.maximum_pressure_step
        };
        let count = if options.maximum_points == 1 {
            1
        } else if requested_intervals >= (options.maximum_points - 2) as f64 {
            options.maximum_points
        } else {
            requested_intervals.round() as usize + 2
        };
        let by_pressure = matches!(options.start, PureSaturationStart::Pressure(_));
        let begin = if by_pressure {
            start.pressure
        } else {
            start.temperature
        };
        let end = if by_pressure {
            endpoint.pressure
        } else {
            endpoint.temperature
        };
        let grid = |i: usize| {
            if options.logarithmic {
                (begin.ln() + (end.ln() - begin.ln()) / (count - 1) as f64 * i as f64).exp()
            } else {
                begin + (end - begin) / (count - 1) as f64 * i as f64
            }
        };
        let (low, high) = if by_pressure {
            (start.temperature, endpoint.temperature)
        } else {
            (start.pressure, endpoint.pressure)
        };
        let first = if baseline {
            problem.evaluate(start.temperature, start.pressure)?.point
        } else {
            start
        };
        let mut points = if count == 1 { Vec::new() } else { vec![first] };
        for i in 1..count - 1 {
            let previous = points.last().unwrap();
            let value = grid(i);
            let (fixed, guess) = if by_pressure {
                (
                    PureSaturationStart::Pressure(value),
                    previous.temperature + (value - previous.pressure) / previous.dpdt,
                )
            } else {
                (
                    PureSaturationStart::Temperature(value),
                    previous.pressure + previous.dpdt * (value - previous.temperature),
                )
            };
            match problem
                .correct(fixed, guess, low, high)
                .or_else(|_| problem.cold(fixed))
                .and_then(|point| {
                    if baseline {
                        problem
                            .evaluate(point.temperature, point.pressure)
                            .map(|p| p.point)
                    } else {
                        Ok(point)
                    }
                }) {
                Ok(point) => points.push(point),
                Err(error) => {
                    return Ok(PureSaturationCurve {
                        points,
                        component,
                        end: PureSaturationEnd::NoConvergence,
                        failure: Some(error),
                    });
                }
            }
        }
        let final_point = problem.endpoint(endpoint);
        let final_point = match final_point {
            Ok(point) => point,
            Err(error) => {
                return Ok(PureSaturationCurve {
                    points,
                    component,
                    end: PureSaturationEnd::NoConvergence,
                    failure: Some(error),
                });
            }
        };
        points.push(final_point);
        Ok(PureSaturationCurve {
            points,
            component,
            end: if baseline {
                PureSaturationEnd::SourceCriticalParameters
            } else {
                PureSaturationEnd::CriticalPoint
            },
            failure: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CubicKind, CubicOptions, MixingRule};
    #[test]
    fn analytic_saturation_derivatives_match_independent_perturbations() {
        for (kind, names, index, mixing_rule, translated) in [
            (
                CubicKind::PengRobinson,
                vec!["C1"],
                0,
                MixingRule::VdW,
                false,
            ),
            (
                CubicKind::SoaveRedlichKwong,
                vec!["CO2"],
                0,
                MixingRule::VdW,
                true,
            ),
            (
                CubicKind::RedlichKwong,
                vec!["C3"],
                0,
                MixingRule::VdW,
                false,
            ),
            (
                CubicKind::VanDerWaals,
                vec!["C1"],
                0,
                MixingRule::VdW,
                false,
            ),
            (
                CubicKind::PengRobinson,
                vec!["C1", "C2", "C3"],
                1,
                MixingRule::Unifac,
                false,
            ),
            (
                CubicKind::PengRobinson,
                vec!["C1", "C2", "C3"],
                2,
                MixingRule::HuronVidal2,
                false,
            ),
        ] {
            let model = Cubic::with_options(
                kind,
                &names,
                CubicOptions {
                    mixing_rule,
                    volume_translation: translated,
                    ..Default::default()
                },
            )
            .unwrap();
            let mut problem = Problem::new(&model, index, false).unwrap();
            let sat = problem
                .cold(PureSaturationStart::Temperature(
                    0.8 * model.components()[index].critical_temperature,
                ))
                .unwrap();
            let t = sat.temperature;
            let p = sat.pressure;
            let dt = 2e-6 * t;
            let dp = 2e-6 * p;
            let analytic = problem.evaluate(t, p).unwrap();
            let ft = (problem.evaluate(t + dt, p).unwrap().residual
                - problem.evaluate(t - dt, p).unwrap().residual)
                / (2. * dt);
            let fp = (problem.evaluate(t, p + dp).unwrap().residual
                - problem.evaluate(t, p - dp).unwrap().residual)
                / (2. * dp);
            assert!(
                (analytic.dt - ft).abs() < 1e-9 + 1e-6 * ft.abs(),
                "{names:?}: {} vs {ft}",
                analytic.dt
            );
            assert!(
                (analytic.dp - fp).abs() < 1e-15 + 1e-6 * fp.abs(),
                "{names:?}: {} vs {fp}",
                analytic.dp
            );
        }
    }
}
