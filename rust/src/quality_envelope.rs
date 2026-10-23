//! Continuation of the fixed-quality equations through log-coordinate turns.
use super::*;
use crate::SaturationKind;
use crate::envelope::{EnvelopeEnd, EnvelopeOptions};

/// A local phase-fraction locus with explicit partial-curve termination.
///
/// ```
/// use thermopack_rs::{Cubic, CubicKind};
/// use thermopack_rs::envelope::{EnvelopeOptions, EnvelopeEnd};
/// let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"])?;
/// let curve = model.quality_envelope(&[0.5, 0.5], 0.25, EnvelopeOptions::default())?;
/// assert_eq!(curve.end, EnvelopeEnd::InitialPressure);
/// assert_eq!(curve.points[0].vapor_fraction, 0.25);
/// # Ok::<(), thermopack_rs::Error>(())
/// ```
#[derive(Debug, Clone)]
pub struct QualityEnvelope {
    pub points: Vec<QualityState>,
    pub end: EnvelopeEnd,
    pub failure: Option<Error>,
}
impl QualityEnvelope {
    fn failed(points: Vec<QualityState>, error: Error) -> Self {
        Self {
            points,
            end: EnvelopeEnd::NoConvergence,
            failure: Some(error),
        }
    }
}
fn reported(mut p: QualityState, baseline: bool) -> QualityState {
    if !baseline && p.liquid_volume > p.vapor_volume {
        std::mem::swap(&mut p.liquid, &mut p.vapor);
        std::mem::swap(&mut p.liquid_volume, &mut p.vapor_volume);
        p.vapor_fraction = 1. - p.vapor_fraction;
    }
    p
}
impl Problem<'_> {
    fn tangent(&mut self, u: &[f64], s: usize) -> Result<Vec<f64>> {
        let d = self.active.len() + 2;
        let mut e = Evaluation::new(d);
        self.evaluate(u, s, u[s], &mut e, true)?;
        let mut rhs = vec![0.; d];
        rhs[d - 1] = 1.;
        solve_linear_in_place(&mut e.j, &mut rhs)?;
        Ok(rhs)
    }
    fn flip(&mut self, u: &mut [f64]) {
        self.beta = 1. - self.beta;
        for v in &mut u[..self.active.len()] {
            *v = -*v;
        }
    }
    fn single_roots(&mut self, u: &[f64]) -> Result<bool> {
        let m = self.active.len();
        let mut e = Evaluation::new(m + 2);
        self.evaluate(u, m, u[m], &mut e, false)?;
        let t = u[m].exp();
        let p = u[m + 1].exp();
        // The source first tests the overall feed, then its incipient side.
        let other = if self.beta > 0.5 { &self.x } else { &self.y };
        Ok(self.model.prepare(t, &self.z)?.volumes(p)?.1 == 1
            && self.model.prepare(t, other)?.volumes(p)?.1 == 1)
    }
}
impl Cubic {
    /// Trace a local fixed-quality locus using adaptive T/P/log-K coordinates.
    /// Density ordering defines reported vapor/liquid labels; swapping the
    /// formulation complements the reported fraction while preserving feed.
    /// Global multiphase stability is not established.
    pub fn quality_envelope(
        &self,
        z: &[f64],
        beta: f64,
        options: EnvelopeOptions,
    ) -> Result<QualityEnvelope> {
        self.trace_quality_envelope(z, beta, options, false)
    }
    /// Preserve source phase labels, quality complement and pressure-cap
    /// conventions. Exact adaptive-grid/failure identity needs qualification.
    pub fn baseline_quality_envelope(
        &self,
        z: &[f64],
        beta: f64,
        options: EnvelopeOptions,
    ) -> Result<QualityEnvelope> {
        self.trace_quality_envelope(z, beta, options, true)
    }
    fn trace_quality_envelope(
        &self,
        z: &[f64],
        beta: f64,
        mut options: EnvelopeOptions,
        baseline: bool,
    ) -> Result<QualityEnvelope> {
        let mut problem = Problem::new(self, z, beta)?;
        let source_pure = baseline
            && problem.z.iter().copied().fold(0., f64::max) == problem.z.iter().sum::<f64>();
        if problem.active.len() == 1 || source_pure || beta == 0. || beta == 1. {
            options.initial_kind = if beta < 0.5 {
                SaturationKind::Bubble
            } else {
                SaturationKind::Dew
            };
            let curve = if baseline {
                self.baseline_phase_envelope(z, options)?
            } else {
                self.phase_envelope(z, options)?
            };
            let pure = problem.active.len() == 1 || source_pure;
            let points = curve
                .points
                .into_iter()
                .map(|p| {
                    let b = if pure {
                        if baseline { 0. } else { beta }
                    } else {
                        f64::from(p.kind == SaturationKind::Dew)
                    };
                    let material_balance_residual = problem
                        .z
                        .iter()
                        .zip(&p.liquid)
                        .zip(&p.vapor)
                        .map(|((&z, &x), &y)| ((1. - b) * x + b * y - z).abs())
                        .fold(0., f64::max);
                    QualityState {
                        temperature: p.temperature,
                        pressure: p.pressure,
                        vapor_fraction: b,
                        feed: problem.z.clone(),
                        liquid: p.liquid,
                        vapor: p.vapor,
                        liquid_volume: p.liquid_volume,
                        vapor_volume: p.vapor_volume,
                        dpdt: p.dpdt,
                        iterations: p.iterations,
                        fugacity_residual: p.fugacity_residual,
                        material_balance_residual,
                    }
                })
                .collect();
            return Ok(QualityEnvelope {
                points,
                end: curve.end,
                failure: curve.failure,
            });
        }
        positive(
            options.maximum_pressure,
            "quality envelope maximum pressure",
        )?;
        positive(options.maximum_step, "quality envelope maximum step")?;
        let minimum = options
            .minimum_temperature
            .unwrap_or(self.temperature_bounds().0);
        positive(minimum, "quality envelope minimum temperature")?;
        let minimum = minimum.max(self.temperature_bounds().0);
        if options.maximum_points == 0 {
            return Err(Error::InvalidInput("quality envelope point limit"));
        }
        let start = self.quality_cold(z, beta, options.start)?;
        if !baseline && (start.pressure >= options.maximum_pressure || start.temperature < minimum)
        {
            return Err(Error::InvalidInput(
                "quality envelope initial state outside limits",
            ));
        }
        let mut u = problem.coordinates(&start.guess())?;
        let (mut s, spec) = problem.constraint(options.start)?;
        problem.baseline = baseline;
        let start = if baseline {
            let (q, it) = problem.solve(u, s, spec)?;
            u = q;
            problem.point(&u, s, it)?
        } else {
            start
        };
        let initial_pressure = start.pressure;
        let m = problem.active.len();
        let d = m + 2;
        let mut previous = u.clone();
        let mut direction = 1.;
        let minimum_step = (0.25 * options.maximum_step).min(0.05);
        let mut step = minimum_step;
        let mut points = vec![reported(start, baseline)];
        let mut switched = false;
        let mut switch_pending = false;
        while points.len() < options.maximum_points {
            let mut tangent = match problem.tangent(&u, s) {
                Ok(t) => t,
                Err(e) => return Ok(QualityEnvelope::failed(points, e)),
            };
            let sensitive = (0..d)
                .max_by(|&a, &b| tangent[a].abs().total_cmp(&tangent[b].abs()))
                .unwrap();
            if points.len() > 1 && sensitive != s {
                s = sensitive;
                direction = if u[s] >= previous[s] { 1. } else { -1. };
                let scale = tangent[s];
                for v in &mut tangent {
                    *v /= scale;
                }
            }
            let mut old = u.clone();
            let mut accepted = None;
            let attempts: &[f64] = &[1., 0.5, 2.];
            for &multiplier in attempts {
                let mut length = direction * step * multiplier;
                if !baseline && s < m && (u[s] + length).abs() < 1e-4 {
                    length = direction * (u[s].abs() + 1e-4);
                }
                let mut trial: Vec<_> = u
                    .iter()
                    .zip(&tangent)
                    .map(|(u, t)| u + t * length)
                    .collect();
                if s == m + 1 && trial[s] > options.maximum_pressure.ln() {
                    length = options.maximum_pressure.ln() - u[s];
                    trial = u
                        .iter()
                        .zip(&tangent)
                        .map(|(u, t)| u + t * length)
                        .collect();
                }
                if !baseline {
                    let bound = if trial[m] < minimum.ln() && trial[m] < u[m] {
                        Some((minimum, EnvelopeEnd::MinimumTemperature))
                    } else if trial[m] > self.temperature_bounds().1.ln() && trial[m] > u[m] {
                        Some((self.temperature_bounds().1, EnvelopeEnd::MaximumTemperature))
                    } else {
                        None
                    };
                    if let Some((temperature, mut end)) = bound
                        && let Ok((mut q, mut it)) =
                            problem.solve(trial.clone(), m, temperature.ln())
                    {
                        let pressure = q[m + 1].exp();
                        let cap = if pressure > options.maximum_pressure {
                            Some((options.maximum_pressure, EnvelopeEnd::MaximumPressure))
                        } else if pressure < initial_pressure {
                            Some((initial_pressure, EnvelopeEnd::InitialPressure))
                        } else {
                            None
                        };
                        let mut valid = true;
                        let mut fixed = m;
                        if let Some((p, reason)) = cap {
                            match problem.solve(q.clone(), m + 1, p.ln()) {
                                Ok((next, count)) => {
                                    q = next;
                                    it = count;
                                    end = reason;
                                    fixed = m + 1;
                                }
                                Err(_) => valid = false,
                            }
                        }
                        if valid
                            && q[m].exp() >= minimum * (1. - 1e-10)
                            && q[m].exp() <= self.temperature_bounds().1 * (1. + 1e-10)
                            && q[m + 1].exp() >= initial_pressure * (1. - 1e-10)
                            && q[m + 1].exp() <= options.maximum_pressure * (1. + 1e-10)
                            && (q[m].exp() - u[m].exp()).abs() <= 25.
                            && (q[m + 1].exp() - u[m + 1].exp()).abs() <= 1e6
                        {
                            return Ok(match problem.point(&q, fixed, it) {
                                Ok(p) => {
                                    points.push(reported(p, false));
                                    QualityEnvelope {
                                        points,
                                        end,
                                        failure: None,
                                    }
                                }
                                Err(e) => QualityEnvelope::failed(points, e),
                            });
                        }
                    }
                }
                let spec = trial[s];
                if let Ok((q, it)) = problem.solve(trial, s, spec)
                    && (q[m].exp() - u[m].exp()).abs() <= 25.
                    && (q[m + 1].exp() - u[m + 1].exp()).abs() <= 1e6
                {
                    accepted = Some((q, it));
                    break;
                }
            }
            if accepted.is_none() {
                // Source critical-crossing recovery: lock a nontrivial K,
                // hold the predicted T/P and try the other formulation.
                if s >= m {
                    let next_s = (0..m)
                        .max_by(|&a, &b| tangent[a].abs().total_cmp(&tangent[b].abs()))
                        .unwrap();
                    direction = (direction * tangent[s] * tangent[next_s]).signum();
                    s = next_s;
                    let scale = tangent[s];
                    for v in &mut tangent {
                        *v /= scale;
                    }
                    tangent[m] = 0.;
                    tangent[m + 1] = 0.;
                }
                let length = if (u[s] + direction * step).abs() < 1e-4 {
                    u[s].abs() + 1e-4
                } else {
                    step
                };
                let mut trial: Vec<_> = u
                    .iter()
                    .zip(&tangent)
                    .map(|(u, t)| u + t * direction * length)
                    .collect();
                if !switched {
                    problem.flip(&mut trial);
                }
                let spec = trial[s];
                if let Ok((q, it)) = problem.solve(trial.clone(), s, spec)
                    && (q[m].exp() - u[m].exp()).abs() <= 25.
                    && (q[m + 1].exp() - u[m + 1].exp()).abs() <= 1e6
                {
                    if !switched {
                        for v in &mut old[..m] {
                            *v = -*v;
                        }
                        for v in &mut u[..m] {
                            *v = -*v;
                        }
                        if s < m {
                            direction = -direction;
                        }
                        switched = true;
                        switch_pending = false;
                    }
                    accepted = Some((q, it));
                } else if !switched {
                    problem.flip(&mut trial);
                }
            }
            let Some((mut next, iterations)) = accepted else {
                return Ok(QualityEnvelope::failed(
                    points,
                    Error::NotConverged("quality envelope continuation step"),
                ));
            };
            let boundary = if next[m] < u[m] && next[m].exp() < minimum + 0.01 {
                Some((m, minimum.ln(), EnvelopeEnd::MinimumTemperature))
            } else if !baseline && next[m].exp() >= self.temperature_bounds().1 * (1. - 1e-12) {
                Some((
                    m,
                    self.temperature_bounds().1.ln(),
                    EnvelopeEnd::MaximumTemperature,
                ))
            } else if next[m + 1].exp() < initial_pressure {
                Some((m + 1, initial_pressure.ln(), EnvelopeEnd::InitialPressure))
            } else if next[m + 1].exp() >= options.maximum_pressure
                || !baseline && next[m + 1] >= options.maximum_pressure.ln() - 1e-12
            {
                Some((
                    m + 1,
                    options.maximum_pressure.ln(),
                    EnvelopeEnd::MaximumPressure,
                ))
            } else {
                None
            };
            if let Some((fixed, spec, end)) = boundary {
                let corrected = problem
                    .solve(next, fixed, spec)
                    .and_then(|(q, it)| problem.point(&q, fixed, it));
                return Ok(match corrected {
                    Ok(p) => {
                        points.push(reported(p, baseline));
                        QualityEnvelope {
                            points,
                            end,
                            failure: None,
                        }
                    }
                    Err(e) => QualityEnvelope::failed(points, e),
                });
            }
            match problem.point(&next, s, iterations) {
                Ok(p) => points.push(reported(p, baseline)),
                Err(e) => return Ok(QualityEnvelope::failed(points, e)),
            }
            previous = old;
            if !switched && switch_pending {
                let original = next.clone();
                problem.flip(&mut next);
                let spec = next[s];
                if let Ok((q, _)) = problem.solve(next.clone(), s, spec)
                    && (q[m] - original[m]).abs() < 1e-5
                    && (q[m + 1] - original[m + 1]).abs() < 1e-5
                {
                    next = q;
                    for v in &mut previous[..m] {
                        *v = -*v;
                    }
                    if s < m {
                        direction = -direction;
                    }
                    switched = true;
                } else {
                    problem.flip(&mut next);
                }
            }
            if !switched {
                switch_pending = match problem.single_roots(&next) {
                    Ok(v) => v,
                    Err(e) => return Ok(QualityEnvelope::failed(points, e)),
                };
            }
            u = next;
            if iterations < 3 {
                step *= 1.2;
            } else if iterations > 5 {
                step /= 1.2;
            }
            step = step.clamp(minimum_step, options.maximum_step);
        }
        Ok(QualityEnvelope {
            points,
            end: EnvelopeEnd::PointLimit,
            failure: None,
        })
    }
}
