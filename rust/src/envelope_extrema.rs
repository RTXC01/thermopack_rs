//! Local pressure/temperature extrema of a mixture bubble/dew envelope.
use super::*;

/// The stationary coordinate on a local saturation branch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnvelopeExtremum {
    Pressure,
    Temperature,
}

#[derive(Debug, Clone)]
pub struct EnvelopeExtremumPoint {
    pub kind: EnvelopeExtremum,
    pub state: MixtureSaturation,
    /// Dimensionless weighted fugacity derivative: T*d(ln phi)/dT for a
    /// pressure extremum, P*d(ln phi)/dP for a temperature extremum.
    pub stationarity_residual: f64,
}

#[derive(Debug, Clone)]
pub struct EnvelopeExtremumFailure {
    pub kind: EnvelopeExtremum,
    pub index: usize,
    pub error: Error,
}

/// Largest successfully refined interior maxima detected in the supplied
/// sampled trace. Missing candidates and failed refinements remain explicit.
/// A partial or metastable trace does not establish global envelope maxima.
///
/// ```
/// use thermopack_rs::{Cubic, CubicKind};
/// use thermopack_rs::envelope::EnvelopeOptions;
/// let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"])?;
/// let curve = model.phase_envelope(&[0.5, 0.5], EnvelopeOptions::default())?;
/// let maxima = model.envelope_extrema(&[0.5, 0.5], &curve)?;
/// assert!(maxima.failures.is_empty());
/// assert!(maxima.pressure.unwrap().stationarity_residual.abs() < 1e-9);
/// # Ok::<(), thermopack_rs::Error>(())
/// ```
#[derive(Debug, Clone, Default)]
pub struct EnvelopeExtrema {
    pub pressure: Option<EnvelopeExtremumPoint>,
    pub temperature: Option<EnvelopeExtremumPoint>,
    pub failures: Vec<EnvelopeExtremumFailure>,
}

struct ExtremumSeed {
    point: Vec<f64>,
    tangent: Vec<f64>,
}

const EXTREMUM_SEED_CACHE_CAPACITY: usize = 8;

impl Problem<'_> {
    fn bracket_extremum(
        &mut self,
        initial: &[f64],
        kind: EnvelopeExtremum,
    ) -> Result<EnvelopeExtremumPoint> {
        let m = self.active.len();
        let d = m + 2;
        let s = if kind == EnvelopeExtremum::Pressure {
            m
        } else {
            m + 1
        };
        let old_tolerance = self.tolerance;
        self.tolerance = 1e-12;
        let attempt = (|| {
            let (center, mut total) = self.solve(initial.to_vec(), s, initial[s])?;
            let tangent = self.tangent(&center, s)?;
            let origin = center[s];
            // Local to this model/feed/formulation/specification. Failed or
            // unresolved corrections never enter the bounded seed cache.
            let mut seeds = Vec::<ExtremumSeed>::with_capacity(EXTREMUM_SEED_CACHE_CAPACITY);
            let mut evaluate = |value: f64| -> Result<(Vec<f64>, f64)> {
                let nearby: Option<Vec<_>> = seeds
                    .iter()
                    .min_by(|a, b| {
                        (a.point[s] - value)
                            .abs()
                            .total_cmp(&(b.point[s] - value).abs())
                    })
                    .filter(|p| p.point != center)
                    .map(|p| {
                        p.point
                            .iter()
                            .zip(&p.tangent)
                            .map(|(u, t)| u + t * (value - p.point[s]))
                            .collect()
                    });
                let mut last = Error::NotConverged("envelope extremum inner correction");
                for inverted in [false, true] {
                    // Try nearby continuation first, then retain both original
                    // predictors and their ratio inversions as fallbacks.
                    for predictor in 0..3 {
                        let mut seed = match predictor {
                            0 if !inverted => match &nearby {
                                Some(p) => p.clone(),
                                None => continue,
                            },
                            0 => continue,
                            1 => center
                                .iter()
                                .zip(&tangent)
                                .map(|(u, t)| u + t * (value - origin))
                                .collect(),
                            _ => center.clone(),
                        };
                        if inverted {
                            for k in &mut seed[..m] {
                                *k = -*k;
                            }
                        }
                        match self.solve(seed, s, value) {
                            Ok((u, it)) => {
                                total += it;
                                let mut e = Evaluation::new(d);
                                match self.extremum_equations(&u, kind, &mut e) {
                                    Ok(()) => {
                                        let stationarity = e.f[d - 1];
                                        if stationarity.is_finite()
                                            && let Ok(tangent) = self.tangent(&u, s)
                                            && tangent.iter().all(|v| v.is_finite())
                                        {
                                            if seeds.len() == EXTREMUM_SEED_CACHE_CAPACITY {
                                                seeds.remove(0);
                                            }
                                            seeds.push(ExtremumSeed {
                                                point: u.clone(),
                                                tangent,
                                            });
                                        }
                                        return Ok((u, stationarity));
                                    }
                                    Err(error) => last = error,
                                }
                            }
                            Err(error) => last = error,
                        }
                    }
                }
                Err(last)
            };
            let (center, fc) = evaluate(origin)?;
            let mut solution = if fc.abs() < 1e-9 {
                Some((center, fc))
            } else {
                None
            };
            let mut bracket = None;
            if solution.is_none() {
                for radius in [0.001, 0.002, 0.005, 0.01, 0.02, 0.04, 0.08] {
                    for value in [origin - radius, origin + radius] {
                        if let Ok((u, f)) = evaluate(value) {
                            if f.abs() < 1e-9 {
                                solution = Some((u, f));
                                break;
                            }
                            if f * fc < 0. {
                                bracket = Some(if value < origin {
                                    (value, origin, f, fc)
                                } else {
                                    (origin, value, fc, f)
                                });
                                break;
                            }
                        }
                    }
                    if solution.is_some() || bracket.is_some() {
                        break;
                    }
                }
            }
            if solution.is_none() {
                let (mut lo, mut hi, mut fl, mut fh) =
                    bracket.ok_or(Error::NotConverged("envelope extremum bracket"))?;
                for _ in 0..80 {
                    let secant = (lo * fh - hi * fl) / (fh - fl);
                    let value = if secant > lo + 0.1 * (hi - lo) && secant < hi - 0.1 * (hi - lo) {
                        secant
                    } else {
                        0.5 * (lo + hi)
                    };
                    let (u, f) = evaluate(value)?;
                    if f.abs() < 1e-9 {
                        solution = Some((u, f));
                        break;
                    }
                    if f * fl > 0. {
                        lo = value;
                        fl = f;
                    } else {
                        hi = value;
                        fh = f;
                    }
                    if hi - lo < 8. * f64::EPSILON * value.abs().max(1.) {
                        break;
                    }
                }
            }
            let (u, f) =
                solution.ok_or(Error::NotConverged("envelope extremum scalar refinement"))?;
            let state = self.point(&u, s, total)?;
            let contrast = state
                .liquid
                .iter()
                .zip(&state.vapor)
                .map(|(x, y)| {
                    let d = x - y;
                    d * d
                })
                .sum::<f64>();
            Ok(EnvelopeExtremumPoint {
                kind,
                state,
                stationarity_residual: f * contrast,
            })
        })();
        self.tolerance = old_tolerance;
        attempt
    }

    fn extremum_equations(
        &mut self,
        u: &[f64],
        kind: EnvelopeExtremum,
        out: &mut Evaluation,
    ) -> Result<()> {
        let m = self.active.len();
        self.evaluate(u, m, u[m], out, true)?;
        let (parent, other, scale) = if kind == EnvelopeExtremum::Pressure {
            (&self.parent.dt, &self.other.dt, u[m].exp())
        } else {
            (&self.parent.dp, &self.other.dp, u[m + 1].exp())
        };
        out.f[m + 1] = self
            .active
            .iter()
            .map(|&i| self.weights[i] * (other[i] - parent[i]))
            .sum::<f64>()
            * scale;
        let contrast = self
            .composition
            .iter()
            .zip(&self.z)
            .map(|(x, z)| {
                let d = x - z;
                d * d
            })
            .sum::<f64>();
        if contrast < 1e-12 {
            return Err(Error::NotConverged("coincident envelope extremum phases"));
        }
        let rounding = 16.
            * f64::EPSILON
            * scale
            * self
                .active
                .iter()
                .map(|&i| self.weights[i].abs() * (other[i].abs() + parent[i].abs()))
                .sum::<f64>();
        if rounding / contrast > 1e-8 {
            return Err(Error::NotConverged(
                "unresolved envelope extremum phase contrast",
            ));
        }
        // The unscaled stationarity equation vanishes quadratically as the
        // phases coalesce, even away from a genuine P/T extremum. Scaling by
        // the squared composition contrast prevents false critical roots.
        out.f[m + 1] /= contrast;
        out.norm = out.f.iter().map(|v| v.abs()).fold(0., f64::max);
        if !out.norm.is_finite() {
            return Err(Error::NotConverged("nonfinite envelope extremum equations"));
        }
        Ok(())
    }

    fn extremum_jacobian(
        &mut self,
        u: &[f64],
        kind: EnvelopeExtremum,
        out: &mut Evaluation,
    ) -> Result<()> {
        self.extremum_equations(u, kind, out)?;
        let d = self.dim();
        let mut trial = u.to_vec();
        let mut e = Evaluation::new(d);
        // Only the stationarity row needs second fugacity derivatives.
        // The coexistence and material-balance rows remain analytic.
        for j in 0..d {
            let h = 2e-6;
            trial[j] = u[j] + h;
            self.extremum_equations(&trial, kind, &mut e)?;
            let plus = e.f[d - 1];
            trial[j] = u[j] - h;
            self.extremum_equations(&trial, kind, &mut e)?;
            out.j[(d - 1) * d + j] = (plus - e.f[d - 1]) / (2. * h);
            trial[j] = u[j];
        }
        Ok(())
    }
}

impl Cubic {
    /// Refine a stationary P or T on a mixture bubble/dew branch from a warm
    /// guess. Brackets stationarity along corrected coexistence states, with
    /// an augmented Newton fallback using a central-difference stationarity
    /// row. A stationary point is not necessarily a maximum.
    /// Critical/pure degeneracies and global stability are separate problems.
    pub fn refine_envelope_extremum(
        &self,
        z: &[f64],
        guess: &MixtureSaturation,
        kind: EnvelopeExtremum,
    ) -> Result<EnvelopeExtremumPoint> {
        let mut problem = Problem::new(self, z, guess.kind)?;
        let d = problem.dim();
        let m = d - 2;
        let mut u = problem.coordinates(guess)?;
        if let Ok(point) = problem.bracket_extremum(&u, kind) {
            return Ok(point);
        }
        let initial_direction = u[..m].to_vec();
        let mut e = Evaluation::new(d);
        let mut trial = Evaluation::new(d);
        for iteration in 0..60 {
            problem.extremum_jacobian(&u, kind, &mut e)?;
            if e.norm < 1e-10 {
                if u[..m].iter().map(|v| v.abs()).fold(0., f64::max) < 1e-6 {
                    return Err(Error::NotConverged("coincident envelope extremum"));
                }
                let s = if kind == EnvelopeExtremum::Pressure {
                    m
                } else {
                    m + 1
                };
                let state = problem.point(&u, s, iteration)?;
                return Ok(EnvelopeExtremumPoint {
                    kind,
                    stationarity_residual: e.f[d - 1]
                        * state
                            .liquid
                            .iter()
                            .zip(&state.vapor)
                            .map(|(x, y)| {
                                let d = x - y;
                                d * d
                            })
                            .sum::<f64>(),
                    state,
                });
            }
            let mut step: Vec<_> = e.f.iter().map(|v| -v).collect();
            solve_linear_in_place(&mut e.j, &mut step)?;
            let largest = step
                .iter()
                .enumerate()
                .map(|(i, v)| v.abs() / if i < m { 2. } else { 0.1 })
                .fold(1., f64::max);
            let mut fraction = 1. / largest;
            let mut accepted = None;
            for _ in 0..16 {
                let q: Vec<_> = u
                    .iter()
                    .zip(&step)
                    .map(|(u, dx)| u + fraction * dx)
                    .collect();
                if q[..m]
                    .iter()
                    .zip(&initial_direction)
                    .map(|(a, b)| a * b)
                    .sum::<f64>()
                    > 0.
                    && problem.extremum_equations(&q, kind, &mut trial).is_ok()
                    && trial.norm < e.norm
                {
                    accepted = Some(q);
                    break;
                }
                fraction *= 0.5;
            }
            u = accepted.ok_or(Error::NotConverged("envelope extremum backtracking"))?;
        }
        Err(Error::NotConverged("envelope extremum iteration"))
    }

    /// Detect sampled interior P/T peaks and refine them to stationary
    /// coexistence states. Reports every failed candidate; does not replace
    /// the original trace or treat endpoint limits as stationary maxima.
    pub fn envelope_extrema(&self, z: &[f64], curve: &PhaseEnvelope) -> Result<EnvelopeExtrema> {
        self.prepare(300., z)?;
        let mut result = EnvelopeExtrema::default();
        for (i, points) in curve.points.windows(3).enumerate() {
            let [a, b, c] = points else { unreachable!() };
            for kind in [EnvelopeExtremum::Pressure, EnvelopeExtremum::Temperature] {
                let value = |p: &MixtureSaturation| {
                    if kind == EnvelopeExtremum::Pressure {
                        p.pressure
                    } else {
                        p.temperature
                    }
                };
                if value(b) <= value(a) || value(b) < value(c) {
                    continue;
                }
                let attempt = self
                    .refine_envelope_extremum(z, b, kind)
                    .or_else(|_| self.refine_envelope_extremum(z, a, kind))
                    .or_else(|_| self.refine_envelope_extremum(z, c, kind))
                    .and_then(|p| {
                        let coordinate = |p: &MixtureSaturation| {
                            if kind == EnvelopeExtremum::Pressure {
                                p.temperature
                            } else {
                                p.pressure
                            }
                        };
                        let lo = coordinate(a).min(coordinate(b)).min(coordinate(c));
                        let hi = coordinate(a).max(coordinate(b)).max(coordinate(c));
                        // A formulation switch can jump across both nearby turns.
                        // Permit one neighboring cell beyond that switched sample.
                        let margin = if a.kind != b.kind || b.kind != c.kind {
                            hi - lo
                        } else {
                            0.
                        };
                        if coordinate(&p.state) < lo - margin || coordinate(&p.state) > hi + margin
                        {
                            Err(Error::NotConverged(
                                "envelope extremum left sampled bracket",
                            ))
                        } else {
                            for side in [-1., 1.] {
                                let constraint = if kind == EnvelopeExtremum::Pressure {
                                    EnvelopeConstraint::Temperature(
                                        p.state.temperature * (1. + side * 5e-5),
                                    )
                                } else {
                                    EnvelopeConstraint::Pressure(
                                        p.state.pressure * (1. + side * 5e-5),
                                    )
                                };
                                let neighbor = self.refine_envelope_state(&p.state, constraint)?;
                                if value(&neighbor) >= value(&p.state) {
                                    return Err(Error::NotConverged(
                                        "stationary envelope candidate is not a maximum",
                                    ));
                                }
                            }
                            Ok(p)
                        }
                    });
                match attempt {
                    Ok(p) => {
                        let entry = if kind == EnvelopeExtremum::Pressure {
                            &mut result.pressure
                        } else {
                            &mut result.temperature
                        };
                        if entry
                            .as_ref()
                            .is_none_or(|old| value(&old.state) < value(&p.state))
                        {
                            *entry = Some(p);
                        }
                    }
                    Err(error) => result.failures.push(EnvelopeExtremumFailure {
                        kind,
                        index: i + 1,
                        error,
                    }),
                }
            }
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CubicKind, CubicOptions, MixingRule};
    #[test]
    fn scaled_stationarity_augmented_jacobian_differentiates_all_rows_off_root() {
        for (rule, translated, quantum) in [
            (MixingRule::VdW, false, false),
            (MixingRule::Unifac, false, false),
            (MixingRule::HuronVidal1, false, false),
            (MixingRule::VdW, true, false),
            (MixingRule::VdW, false, true),
        ] {
            let model = if quantum {
                Cubic::quantum(&["H2", "NE"]).unwrap()
            } else {
                Cubic::with_options(
                    CubicKind::PengRobinson,
                    &["C1", "C3"],
                    CubicOptions {
                        mixing_rule: rule,
                        volume_translation: translated,
                        ..Default::default()
                    },
                )
                .unwrap()
            };
            for side in [SaturationKind::Bubble, SaturationKind::Dew] {
                let mut problem = Problem::new(&model, &[0.5, 0.5], side).unwrap();
                let u = vec![
                    0.7,
                    -0.4,
                    if quantum { 30_f64.ln() } else { 250_f64.ln() },
                    if quantum { 1e5_f64.ln() } else { 1e6_f64.ln() },
                ];
                for kind in [EnvelopeExtremum::Pressure, EnvelopeExtremum::Temperature] {
                    let mut e = Evaluation::new(4);
                    problem.extremum_jacobian(&u, kind, &mut e).unwrap();
                    for j in 0..4 {
                        let h = 7e-6;
                        let mut plus = u.clone();
                        plus[j] += h;
                        let mut minus = u.clone();
                        minus[j] -= h;
                        let mut a = Evaluation::new(4);
                        let mut b = Evaluation::new(4);
                        problem.extremum_equations(&plus, kind, &mut a).unwrap();
                        problem.extremum_equations(&minus, kind, &mut b).unwrap();
                        for i in 0..4 {
                            let finite = (a.f[i] - b.f[i]) / (2. * h);
                            let actual = e.j[4 * i + j];
                            assert!(
                                (finite - actual).abs()
                                    < 1e-5 + 3e-5 * finite.abs().max(actual.abs()),
                                "{rule:?} {translated} {quantum} {side:?} {kind:?} ({i},{j}): {actual} vs {finite}"
                            );
                        }
                    }
                }
            }
        }
    }
}
