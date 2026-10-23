//! Mixture VLE continuation in log(K), log(T), log(P) coordinates.
//! Equations and coordinate switching follow `saturation`/`saturation_curve`.
//!
//! ```
//! use thermopack_rs::{Cubic, CubicKind};
//! use thermopack_rs::envelope::{EnvelopeOptions, EnvelopeEnd};
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"])?;
//! let curve = model.phase_envelope(&[0.5, 0.5], EnvelopeOptions::default())?;
//! assert_eq!(curve.end, EnvelopeEnd::InitialPressure);
//! assert!(curve.points.iter().all(|p| p.liquid_volume <= p.vapor_volume));
//! # Ok(())
//! # }
//! ```
use crate::{
    Cubic, Error, MixtureSaturation, Phase, PreparedState, Result, SaturationKind, error::positive,
    numerics::solve_linear_in_place,
};

#[derive(Debug, Clone, Copy)]
pub enum EnvelopeConstraint {
    Temperature(f64),
    Pressure(f64),
    LogK { component: usize, value: f64 },
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CubicKind, CubicOptions, MixingRule};
    #[test]
    fn phase_jet_cache_preserves_values_derivatives_and_root_counts_after_trials() {
        let model = Cubic::with_options(
            CubicKind::PengRobinson,
            &["C1", "C3"],
            CubicOptions {
                mixing_rule: MixingRule::Unifac,
                ..Default::default()
            },
        )
        .unwrap();
        let mut jet = Jet::new(&model, 2);
        for (t, p, z, phase) in [
            (250., 1e5, [0.5, 0.5], Phase::Liquid),
            (250., 1e5, [0.5, 0.5], Phase::Vapor),
            (250., 2e6, [0.5, 0.5], Phase::Liquid),
            (250., 2e6, [0.2, 0.8], Phase::Vapor),
            (260., 2e6, [0.2, 0.8], Phase::Liquid),
            (250., 1e5, [0.5, 0.5], Phase::Vapor),
        ] {
            let v = jet.fill(t, p, &z, phase, true).unwrap();
            let fresh = model.prepare(t, &z).unwrap();
            let expected = fresh.volume(p, phase).unwrap();
            assert_eq!(v, expected);
            assert_eq!(jet.roots(t, p, &z).unwrap().1, fresh.volumes(p).unwrap().1);
            let mut phi = vec![0.; 2];
            let mut dt = vec![0.; 2];
            let mut dp = vec![0.; 2];
            let mut dn = vec![0.; 4];
            fresh.ln_fugacity_tp_into(p, v, &mut phi).unwrap();
            fresh
                .fugacity_derivatives_into(v, &mut dt, &mut dp, &mut dn)
                .unwrap();
            assert_eq!(jet.phi, phi);
            assert_eq!(jet.dt, dt);
            assert_eq!(jet.dp, dp);
            assert_eq!(jet.dn, dn);
            assert!(jet.fill(-1., p, &z, phase, true).is_err());
            assert!(jet.fill(t, -1., &z, phase, true).is_err());
        }
    }
    #[test]
    fn quantum_envelope_jacobian_matches_independent_differences() {
        let model = Cubic::quantum(&["H2", "NE"]).unwrap();
        let z = [0.5, 0.5];
        let seed = model.dew_temperature_with_guess(1e6, &z, 33.).unwrap();
        let mut problem = Problem::new(&model, &z, SaturationKind::Dew).unwrap();
        let u = problem.coordinates(&seed).unwrap();
        let d = problem.dim();
        let s = d - 1;
        let mut base = Evaluation::new(d);
        problem.evaluate(&u, s, u[s], &mut base, true).unwrap();
        for j in 0..d {
            let h = 2e-6;
            let mut plus = u.clone();
            let mut minus = u.clone();
            plus[j] += h;
            minus[j] -= h;
            let mut a = Evaluation::new(d);
            let mut b = Evaluation::new(d);
            problem.evaluate(&plus, s, u[s], &mut a, false).unwrap();
            problem.evaluate(&minus, s, u[s], &mut b, false).unwrap();
            for i in 0..d {
                let finite = (a.f[i] - b.f[i]) / (2. * h);
                let analytic = base.j[i * d + j];
                assert!(
                    (finite - analytic).abs() < 1e-6 + 1e-5 * finite.abs().max(analytic.abs()),
                    "{i},{j}: {analytic} vs {finite}"
                );
            }
        }
    }
    #[test]
    fn tangent_cache_tracks_coordinates_specification_and_formulation() {
        let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
        let z = [0.5, 0.5];
        let seed = model.bubble_pressure(250., &z).unwrap();
        let mut problem = Problem::new(&model, &z, SaturationKind::Bubble).unwrap();
        let mut q = problem.coordinates(&seed).unwrap();
        let m = problem.active.len();
        problem.tangent(&q, m + 1).unwrap();
        q[m] += 0.001;
        for coordinate in [0, m, m + 1] {
            let actual = problem.tangent(&q, coordinate).unwrap();
            let expected = Problem::new(&model, &z, SaturationKind::Bubble)
                .unwrap()
                .tangent(&q, coordinate)
                .unwrap();
            for (a, b) in actual.iter().zip(expected) {
                assert!((a - b).abs() < 1e-13);
            }
        }
        problem.flip(&mut q);
        assert!(problem.tangent_cache.is_none());
        let actual = problem.tangent(&q, 0).unwrap();
        let expected = Problem::new(&model, &z, SaturationKind::Dew)
            .unwrap()
            .tangent(&q, 0)
            .unwrap();
        for (a, b) in actual.iter().zip(expected) {
            assert!((a - b).abs() < 1e-13);
        }
    }
    #[test]
    fn augmented_envelope_jacobians_match_independent_differences() {
        for (kind, names, z, mixing_rule, translated) in [
            (
                CubicKind::PengRobinson,
                vec!["C1", "C2"],
                vec![0.5, 0.5],
                MixingRule::VdW,
                false,
            ),
            (
                CubicKind::SoaveRedlichKwong,
                vec!["C1", "C2"],
                vec![0.5, 0.5],
                MixingRule::VdW,
                true,
            ),
            (
                CubicKind::PengRobinson,
                vec!["C1", "C3"],
                vec![0.5, 0.5],
                MixingRule::Unifac,
                false,
            ),
            (
                CubicKind::PengRobinson,
                vec!["C1", "C2", "C3"],
                vec![0.5, 0.5, 0.],
                MixingRule::VdW,
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
            for saturation_kind in [SaturationKind::Bubble, SaturationKind::Dew] {
                let seed = if saturation_kind == SaturationKind::Bubble {
                    model.bubble_pressure(250., &z)
                } else {
                    model.dew_pressure(250., &z)
                }
                .unwrap();
                let mut problem = Problem::new(&model, &z, saturation_kind).unwrap();
                let u = problem.coordinates(&seed).unwrap();
                let d = problem.dim();
                let s = 0;
                let mut base = Evaluation::new(d);
                problem.evaluate(&u, s, u[s], &mut base, true).unwrap();
                for j in 0..d {
                    let h = 2e-6;
                    let mut plus = u.clone();
                    let mut minus = u.clone();
                    plus[j] += h;
                    minus[j] -= h;
                    let mut a = Evaluation::new(d);
                    let mut b = Evaluation::new(d);
                    problem.evaluate(&plus, s, u[s], &mut a, false).unwrap();
                    problem.evaluate(&minus, s, u[s], &mut b, false).unwrap();
                    for i in 0..d {
                        let finite = (a.f[i] - b.f[i]) / (2. * h);
                        let analytic = base.j[i * d + j];
                        assert!(
                            (finite - analytic).abs()
                                < 1e-6 + 1e-5 * finite.abs().max(analytic.abs()),
                            "{names:?} {saturation_kind:?} {i},{j}: {analytic} vs {finite}"
                        );
                    }
                }
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct EnvelopeOptions {
    /// Starting pressure (Pa) or temperature (K). LogK is a corrector
    /// constraint, not an automatic-initialization option.
    pub start: EnvelopeConstraint,
    pub initial_kind: SaturationKind,
    pub maximum_pressure: f64,
    /// Defaults to the model's minimum temperature, in K.
    pub minimum_temperature: Option<f64>,
    /// Maximum change of the selected logarithmic continuation coordinate.
    pub maximum_step: f64,
    pub maximum_points: usize,
}
impl Default for EnvelopeOptions {
    fn default() -> Self {
        Self {
            start: EnvelopeConstraint::Pressure(1e5),
            initial_kind: SaturationKind::Dew,
            maximum_pressure: 1.5e7,
            minimum_temperature: None,
            maximum_step: 0.25,
            maximum_points: 1000,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnvelopeEnd {
    InitialPressure,
    MaximumPressure,
    MinimumTemperature,
    PointLimit,
    NoConvergence,
}

#[derive(Debug, Clone)]
pub struct PhaseEnvelope {
    pub points: Vec<MixtureSaturation>,
    pub end: EnvelopeEnd,
    pub failure: Option<Error>,
}

fn reported_point(mut point: MixtureSaturation, source_conventions: bool) -> MixtureSaturation {
    if !source_conventions && point.liquid_volume > point.vapor_volume {
        std::mem::swap(&mut point.liquid, &mut point.vapor);
        std::mem::swap(&mut point.liquid_volume, &mut point.vapor_volume);
        point.kind = if point.kind == SaturationKind::Bubble {
            SaturationKind::Dew
        } else {
            SaturationKind::Bubble
        };
    }
    point
}

struct Jet<'a> {
    model: &'a Cubic,
    prepared: Option<PreparedState<'a>>,
    roots: Option<(f64, [f64; 3], usize)>,
    phi: Vec<f64>,
    dt: Vec<f64>,
    dp: Vec<f64>,
    dn: Vec<f64>,
}
impl<'a> Jet<'a> {
    fn new(model: &'a Cubic, n: usize) -> Self {
        Self {
            model,
            prepared: None,
            roots: None,
            phi: vec![0.; n],
            dt: vec![0.; n],
            dp: vec![0.; n],
            dn: vec![0.; n * n],
        }
    }
    fn roots(&mut self, t: f64, p: f64, z: &[f64]) -> Result<([f64; 3], usize)> {
        // Preparation depends on T and composition, not P. Each phase jet
        // belongs to one immutable model for the duration of this trace.
        if !self
            .prepared
            .as_ref()
            .is_some_and(|s| s.temperature() == t && s.composition() == z)
        {
            self.prepared = Some(self.model.prepare(t, z)?);
            self.roots = None;
        }
        if let Some((pressure, roots, count)) = self.roots
            && pressure == p
        {
            return Ok((roots, count));
        }
        let (roots, count) = self.prepared.as_ref().unwrap().volumes(p)?;
        self.roots = Some((p, roots, count));
        Ok((roots, count))
    }
    fn fill(&mut self, t: f64, p: f64, z: &[f64], phase: Phase, derivatives: bool) -> Result<f64> {
        let (roots, count) = self.roots(t, p, z)?;
        let state = self.prepared.as_ref().unwrap();
        let volume = match phase {
            Phase::Liquid => roots[0],
            Phase::Vapor => roots[count - 1],
            Phase::MinimumGibbs => return Err(Error::InvalidInput("envelope phase selection")),
        };
        state.ln_fugacity_tp_into(p, volume, &mut self.phi)?;
        if derivatives {
            state.fugacity_derivatives_into(volume, &mut self.dt, &mut self.dp, &mut self.dn)?;
        }
        Ok(volume)
    }
}
struct Evaluation {
    f: Vec<f64>,
    j: Vec<f64>,
    norm: f64,
    vl: f64,
    vv: f64,
}
impl Evaluation {
    fn new(n: usize) -> Self {
        Self {
            f: vec![0.; n],
            j: vec![0.; n * n],
            norm: 0.,
            vl: 0.,
            vv: 0.,
        }
    }
}
struct Problem<'a> {
    model: &'a Cubic,
    z: Vec<f64>,
    active: Vec<usize>,
    kind: SaturationKind,
    parent: Jet<'a>,
    other: Jet<'a>,
    weights: Vec<f64>,
    composition: Vec<f64>,
    tolerance: f64,
    tangent_cache: Option<(Vec<f64>, usize, SaturationKind, Vec<f64>)>,
}
impl<'a> Problem<'a> {
    fn new(model: &'a Cubic, z: &[f64], kind: SaturationKind) -> Result<Self> {
        model.prepare(300., z)?;
        let sum = z.iter().sum::<f64>();
        let z: Vec<_> = z.iter().map(|v| v / sum).collect();
        let active: Vec<_> = z
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| (v > 0.).then_some(i))
            .collect();
        if active.len() < 2 {
            return Err(Error::InvalidInput(
                "mixture envelope requires at least two active species",
            ));
        }
        let n = z.len();
        Ok(Self {
            model,
            z,
            active,
            kind,
            parent: Jet::new(model, n),
            other: Jet::new(model, n),
            weights: vec![0.; n],
            composition: vec![0.; n],
            tolerance: 1e-9,
            tangent_cache: None,
        })
    }
    fn dim(&self) -> usize {
        self.active.len() + 2
    }
    fn sign(&self) -> f64 {
        if self.kind == SaturationKind::Bubble {
            1.
        } else {
            -1.
        }
    }
    fn coordinates(&self, p: &MixtureSaturation) -> Result<Vec<f64>> {
        positive(p.temperature, "envelope seed temperature")?;
        positive(p.pressure, "envelope seed pressure")?;
        if p.liquid.len() != self.z.len()
            || p.vapor.len() != self.z.len()
            || p.liquid
                .iter()
                .chain(&p.vapor)
                .any(|v| !v.is_finite() || *v < 0.)
            || (p.liquid.iter().sum::<f64>() - 1.).abs() > 1e-9
            || (p.vapor.iter().sum::<f64>() - 1.).abs() > 1e-9
            || p.liquid
                .iter()
                .zip(&p.vapor)
                .zip(&self.z)
                .any(|((&x, &y), &z)| (z == 0.) != (x == 0.) || (z == 0.) != (y == 0.))
        {
            return Err(Error::InvalidInput(
                "normalized matching envelope phase compositions required",
            ));
        }
        let mut u: Vec<_> = self
            .active
            .iter()
            .map(|&i| (p.vapor[i] / p.liquid[i]).ln())
            .collect();
        u.extend([p.temperature.ln(), p.pressure.ln()]);
        Ok(u)
    }
    fn constraint(&self, c: EnvelopeConstraint) -> Result<(usize, f64)> {
        let m = self.active.len();
        match c {
            EnvelopeConstraint::Temperature(t) => {
                positive(t, "envelope temperature specification")?;
                Ok((m, t.ln()))
            }
            EnvelopeConstraint::Pressure(p) => {
                positive(p, "envelope pressure specification")?;
                Ok((m + 1, p.ln()))
            }
            EnvelopeConstraint::LogK { component, value } => {
                if !value.is_finite() {
                    return Err(Error::InvalidInput("finite envelope log K required"));
                }
                let i = self
                    .active
                    .iter()
                    .position(|&i| i == component)
                    .ok_or(Error::InvalidInput("active envelope component required"))?;
                Ok((i, value))
            }
        }
    }
    fn evaluate(
        &mut self,
        u: &[f64],
        s: usize,
        spec: f64,
        out: &mut Evaluation,
        derivatives: bool,
    ) -> Result<()> {
        let m = self.active.len();
        let d = m + 2;
        let sign = self.sign();
        let t = u[m].exp();
        let p = u[m + 1].exp();
        positive(t, "envelope iterate temperature")?;
        positive(p, "envelope iterate pressure")?;
        if t < self.model.temperature_bounds().0 * (1.0 - 1e-12)
            || t > self.model.temperature_bounds().1 * (1.0 + 1e-12)
            || u[..m].iter().any(|v| !v.is_finite() || v.abs() > 600.)
        {
            return Err(Error::NotConverged("envelope iterate bounds"));
        }
        for (k, &i) in self.active.iter().enumerate() {
            self.weights[i] = self.z[i] * (sign * u[k]).exp();
        }
        let sum = self.weights.iter().sum::<f64>();
        positive(sum, "envelope normalization")?;
        for (value, weight) in self.composition.iter_mut().zip(&self.weights) {
            *value = weight / sum;
        }
        let (fixed, trial) = if self.kind == SaturationKind::Bubble {
            (Phase::Liquid, Phase::Vapor)
        } else {
            (Phase::Vapor, Phase::Liquid)
        };
        let vf = self.parent.fill(t, p, &self.z, fixed, derivatives)?;
        let vo = self
            .other
            .fill(t, p, &self.composition, trial, derivatives)?;
        if self.kind == SaturationKind::Bubble {
            out.vl = vf;
            out.vv = vo;
        } else {
            out.vl = vo;
            out.vv = vf;
        }
        if derivatives {
            out.j.fill(0.);
        }
        for (r, &i) in self.active.iter().enumerate() {
            out.f[r] = self.parent.phi[i] - self.other.phi[i] - sign * u[r];
            if derivatives {
                for (c, &j) in self.active.iter().enumerate() {
                    out.j[r * d + c] = -sign
                        * (f64::from(r == c)
                            + self.other.dn[i * self.z.len() + j] * self.composition[j]);
                }
                out.j[r * d + m] = t * (self.parent.dt[i] - self.other.dt[i]);
                out.j[r * d + m + 1] = p * (self.parent.dp[i] - self.other.dp[i]);
            }
        }
        out.f[m] = sum - 1.;
        out.f[m + 1] = u[s] - spec;
        if derivatives {
            for (k, &i) in self.active.iter().enumerate() {
                out.j[m * d + k] = sign * self.weights[i];
            }
            out.j[(m + 1) * d + s] = 1.;
        }
        out.norm = out.f.iter().map(|f| f.abs()).fold(0., f64::max);
        if !out.norm.is_finite() {
            return Err(Error::NotConverged("nonfinite envelope equations"));
        }
        Ok(())
    }
    fn solve(&mut self, mut u: Vec<f64>, s: usize, spec: f64) -> Result<(Vec<f64>, usize)> {
        let d = self.dim();
        let m = d - 2;
        u[s] = spec;
        let mut current = Evaluation::new(d);
        let mut trial = Evaluation::new(d);
        for iteration in 0..60 {
            self.evaluate(&u, s, spec, &mut current, true)?;
            if current.norm < self.tolerance {
                if u[..m].iter().map(|v| v.abs()).fold(0., f64::max) < 1e-7 {
                    return Err(Error::NotConverged("trivial envelope root"));
                }
                let mut tangent = vec![0.; d];
                tangent[d - 1] = 1.;
                if solve_linear_in_place(&mut current.j, &mut tangent).is_ok() {
                    self.tangent_cache = Some((u.clone(), s, self.kind, tangent));
                }
                return Ok((u, iteration));
            }
            let mut step: Vec<_> = current.f.iter().map(|v| -v).collect();
            solve_linear_in_place(&mut current.j, &mut step)?;
            let largest = step
                .iter()
                .enumerate()
                .map(|(i, v)| v.abs() / if i < m { 2. } else { 0.2 })
                .fold(1., f64::max);
            for v in &mut step {
                *v /= largest;
            }
            let mut fraction = 1.;
            let mut accepted = None;
            for _ in 0..14 {
                let mut q: Vec<_> = u.iter().zip(&step).map(|(u, v)| u + fraction * v).collect();
                q[s] = spec;
                if self.evaluate(&q, s, spec, &mut trial, false).is_ok()
                    && trial.norm < current.norm
                {
                    accepted = Some(q);
                    break;
                }
                fraction *= 0.5;
            }
            u = accepted.ok_or(Error::NotConverged("envelope Newton backtracking"))?;
        }
        Err(Error::NotConverged("envelope Newton iteration"))
    }
    fn tangent(&mut self, u: &[f64], s: usize) -> Result<Vec<f64>> {
        if let Some((point, coordinate, kind, tangent)) = &self.tangent_cache
            && point == u
            && *coordinate == s
            && *kind == self.kind
        {
            return Ok(tangent.clone());
        }
        let d = self.dim();
        let mut e = Evaluation::new(d);
        self.evaluate(u, s, u[s], &mut e, true)?;
        let mut rhs = vec![0.; d];
        rhs[d - 1] = 1.;
        solve_linear_in_place(&mut e.j, &mut rhs)?;
        self.tangent_cache = Some((u.to_vec(), s, self.kind, rhs.clone()));
        Ok(rhs)
    }
    fn point(&mut self, u: &[f64], s: usize, iterations: usize) -> Result<MixtureSaturation> {
        let m = self.active.len();
        let mut e = Evaluation::new(m + 2);
        self.evaluate(u, s, u[s], &mut e, false)?;
        let tangent = self.tangent(u, s)?;
        let t = u[m].exp();
        let p = u[m + 1].exp();
        let (liquid, vapor) = if self.kind == SaturationKind::Bubble {
            (self.z.clone(), self.composition.clone())
        } else {
            (self.composition.clone(), self.z.clone())
        };
        Ok(MixtureSaturation {
            kind: self.kind,
            temperature: t,
            pressure: p,
            liquid,
            vapor,
            liquid_volume: e.vl,
            vapor_volume: e.vv,
            dpdt: p * tangent[m + 1] / (t * tangent[m]),
            iterations,
            fugacity_residual: e.norm,
        })
    }
    fn single_roots(&mut self, u: &[f64]) -> Result<bool> {
        let m = self.active.len();
        let t = u[m].exp();
        let p = u[m + 1].exp();
        let mut other = self.z.clone();
        for (k, &i) in self.active.iter().enumerate() {
            other[i] *= (self.sign() * u[k]).exp();
        }
        let sum = other.iter().sum::<f64>();
        for v in &mut other {
            *v /= sum;
        }
        Ok(self.parent.roots(t, p, &self.z)?.1 == 1 && self.other.roots(t, p, &other)?.1 == 1)
    }
    fn flip(&mut self, u: &mut [f64]) {
        self.tangent_cache = None;
        self.kind = if self.kind == SaturationKind::Bubble {
            SaturationKind::Dew
        } else {
            SaturationKind::Bubble
        };
        for v in &mut u[..self.active.len()] {
            *v = -*v;
        }
    }
}

impl Cubic {
    /// Correct a warm mixture VLE state at specified T, P or one log K value.
    /// The seed's specified-phase composition remains fixed. Allows pressure
    /// and temperature turning points without a cold saturation restart.
    pub fn refine_envelope_state(
        &self,
        seed: &MixtureSaturation,
        constraint: EnvelopeConstraint,
    ) -> Result<MixtureSaturation> {
        let z = if seed.kind == SaturationKind::Bubble {
            &seed.liquid
        } else {
            &seed.vapor
        };
        let mut problem = Problem::new(self, z, seed.kind)?;
        let u = problem.coordinates(seed)?;
        let (s, spec) = problem.constraint(constraint)?;
        let tangent = problem.tangent(&u, s)?;
        let predicted = u
            .iter()
            .zip(tangent)
            .map(|(v, t)| v + t * (spec - u[s]))
            .collect();
        let (u, iterations) = problem
            .solve(predicted, s, spec)
            .or_else(|_| problem.solve(u, s, spec))?;
        problem.point(&u, s, iterations)
    }
    /// Follow a mixture VLE envelope with adaptive most-sensitive log-variable
    /// prediction/correction. Partial curves carry an explicit termination.
    /// This is local VLE continuation; global VLLE stability is not established.
    pub fn phase_envelope(&self, z: &[f64], options: EnvelopeOptions) -> Result<PhaseEnvelope> {
        self.trace_phase_envelope(z, options, false)
    }
    /// Source convergence/pressure-limit conventions. Adaptive-grid and
    /// failure-trajectory identity require separate qualification.
    pub fn baseline_phase_envelope(
        &self,
        z: &[f64],
        options: EnvelopeOptions,
    ) -> Result<PhaseEnvelope> {
        self.trace_phase_envelope(z, options, true)
    }
    fn trace_phase_envelope(
        &self,
        z: &[f64],
        options: EnvelopeOptions,
        source_conventions: bool,
    ) -> Result<PhaseEnvelope> {
        positive(options.maximum_pressure, "envelope maximum pressure")?;
        positive(options.maximum_step, "envelope maximum step")?;
        let requested_minimum_temperature = options
            .minimum_temperature
            .unwrap_or(self.temperature_bounds().0);
        positive(
            requested_minimum_temperature,
            "envelope minimum temperature",
        )?;
        let minimum_temperature = requested_minimum_temperature.max(self.temperature_bounds().0);
        if options.maximum_points == 0 {
            return Err(Error::InvalidInput("envelope tracing limits"));
        }
        let mut problem = Problem::new(self, z, options.initial_kind)?;
        if source_conventions {
            problem.tolerance = 1e-8;
        }
        // Start near congruent coexistence before the mixture correction, as
        // in the fork's PureSat initialization. A Wilson-first Newton solve
        // can instead enter a separate dense-fluid branch (e.g. QC H2/Ne).
        let pseudo_seed = match options.start {
            EnvelopeConstraint::Pressure(p) => self
                .pseudo_saturation_temperature(p, &problem.z)
                .map(|s| s.temperature)
                .or_else(|_| self.pseudo_critical(&problem.z).map(|c| c.temperature)),
            EnvelopeConstraint::Temperature(t) => self
                .pseudo_saturation_pressure(t, &problem.z)
                .map(|s| s.pressure)
                .or_else(|_| self.pseudo_critical(&problem.z).map(|c| c.pressure)),
            _ => Err(Error::InvalidInput(
                "temperature or pressure envelope start required",
            )),
        }
        .ok();
        let start = match (options.start, options.initial_kind) {
            (EnvelopeConstraint::Pressure(p), SaturationKind::Dew) => pseudo_seed.map_or_else(
                || self.dew_temperature(p, &problem.z),
                |guess| {
                    self.dew_temperature_with_guess(p, &problem.z, guess)
                        .or_else(|_| self.dew_temperature(p, &problem.z))
                },
            )?,
            (EnvelopeConstraint::Pressure(p), SaturationKind::Bubble) => pseudo_seed.map_or_else(
                || self.bubble_temperature(p, &problem.z),
                |guess| {
                    self.bubble_temperature_with_guess(p, &problem.z, guess)
                        .or_else(|_| self.bubble_temperature(p, &problem.z))
                },
            )?,
            (EnvelopeConstraint::Temperature(t), SaturationKind::Dew) => pseudo_seed.map_or_else(
                || self.dew_pressure(t, &problem.z),
                |guess| {
                    self.dew_pressure_with_guess(t, &problem.z, guess)
                        .or_else(|_| self.dew_pressure(t, &problem.z))
                },
            )?,
            (EnvelopeConstraint::Temperature(t), SaturationKind::Bubble) => pseudo_seed
                .map_or_else(
                    || self.bubble_pressure(t, &problem.z),
                    |guess| {
                        self.bubble_pressure_with_guess(t, &problem.z, guess)
                            .or_else(|_| self.bubble_pressure(t, &problem.z))
                    },
                )?,
            _ => {
                return Err(Error::InvalidInput(
                    "start envelope at temperature or pressure",
                ));
            }
        };
        if start.pressure >= options.maximum_pressure || start.temperature < minimum_temperature {
            return Err(Error::InvalidInput("envelope initial state outside limits"));
        }
        let minimum_pressure = start.pressure;
        let m = problem.active.len();
        let d = m + 2;
        let mut u = problem.coordinates(&start)?;
        let (mut s, _) = problem.constraint(options.start)?;
        let mut previous = u.clone();
        let mut direction = 1.;
        let minimum_step = (0.25 * options.maximum_step).min(0.05);
        let mut step = minimum_step;
        let mut points = vec![reported_point(start, source_conventions)];
        let mut switched = false;
        let mut switch_pending = false;
        while points.len() < options.maximum_points {
            let mut tangent = match problem.tangent(&u, s) {
                Ok(t) => t,
                Err(e) => {
                    return Ok(PhaseEnvelope {
                        points,
                        end: EnvelopeEnd::NoConvergence,
                        failure: Some(e),
                    });
                }
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
            let old = u.clone();
            let mut accepted = None;
            for multiplier in [1., 0.5, 2., 0.25, 0.125] {
                let mut length = direction * step * multiplier;
                if s < m && (u[s] + length).abs() < 1e-4 {
                    length = direction * (u[s].abs() + 1e-4);
                }
                let mut trial: Vec<_> = u
                    .iter()
                    .zip(&tangent)
                    .map(|(u, t)| u + t * length)
                    .collect();
                if s == m + 1 && trial[s] > options.maximum_pressure.ln() {
                    let length = options.maximum_pressure.ln() - u[s];
                    trial = u
                        .iter()
                        .zip(&tangent)
                        .map(|(u, t)| u + t * length)
                        .collect();
                }
                if !source_conventions
                    && trial[m] < minimum_temperature.ln()
                    && trial[m] < u[m]
                    && let Ok((q, it)) = problem.solve(trial.clone(), m, minimum_temperature.ln())
                    && (q[m].exp() - u[m].exp()).abs() <= 25.
                    && (q[m + 1].exp() - u[m + 1].exp()).abs() <= 1e6
                {
                    let point = match problem.point(&q, m, it) {
                        Ok(p) => p,
                        Err(e) => {
                            return Ok(PhaseEnvelope {
                                points,
                                end: EnvelopeEnd::NoConvergence,
                                failure: Some(e),
                            });
                        }
                    };
                    points.push(reported_point(point, source_conventions));
                    return Ok(PhaseEnvelope {
                        points,
                        end: EnvelopeEnd::MinimumTemperature,
                        failure: None,
                    });
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
            let Some((mut next, iterations)) = accepted else {
                return Ok(PhaseEnvelope {
                    points,
                    end: EnvelopeEnd::NoConvergence,
                    failure: Some(Error::NotConverged("envelope predictor/corrector step")),
                });
            };
            let boundary = if next[m] < u[m] && next[m].exp() < minimum_temperature + 0.01 {
                Some((m, minimum_temperature.ln(), EnvelopeEnd::MinimumTemperature))
            } else if next[m + 1].exp() < minimum_pressure {
                Some((m + 1, minimum_pressure.ln(), EnvelopeEnd::InitialPressure))
            } else if next[m + 1].exp() >= options.maximum_pressure
                || !source_conventions && next[m + 1] >= options.maximum_pressure.ln() - 1e-12
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
                match problem
                    .solve(next, fixed, spec)
                    .and_then(|(q, it)| problem.point(&q, fixed, it))
                {
                    Ok(point) => {
                        points.push(reported_point(point, source_conventions));
                        return Ok(PhaseEnvelope {
                            points,
                            end,
                            failure: None,
                        });
                    }
                    Err(e) => {
                        return Ok(PhaseEnvelope {
                            points,
                            end: EnvelopeEnd::NoConvergence,
                            failure: Some(e),
                        });
                    }
                }
            }
            let point = match problem.point(&next, s, iterations) {
                Ok(p) => p,
                Err(e) => {
                    return Ok(PhaseEnvelope {
                        points,
                        end: EnvelopeEnd::NoConvergence,
                        failure: Some(e),
                    });
                }
            };
            points.push(reported_point(point, source_conventions));
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
                    Ok(value) => value,
                    Err(e) => {
                        return Ok(PhaseEnvelope {
                            points,
                            end: EnvelopeEnd::NoConvergence,
                            failure: Some(e),
                        });
                    }
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
        Ok(PhaseEnvelope {
            points,
            end: EnvelopeEnd::PointLimit,
            failure: None,
        })
    }
}
