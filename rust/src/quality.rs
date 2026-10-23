//! Local two-phase states at a specified molar vapor fraction.
//!
//! The fork's `saturation::sat_newton` equations use log K, log T and log P.
//! Both phase compositions respond to K at interior quality. The analytic
//! Jacobian includes their normalization away from the material-balance root.
//! Global stability and critical/coincident-phase continuation are separate.
//!
//! ```
//! use thermopack_rs::{Cubic, CubicKind};
//! let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"])?;
//! let state = model.quality_temperature(1e6, &[0.5, 0.5], 0.4)?;
//! assert!(state.fugacity_residual < 1e-9);
//! assert_eq!(state.vapor_fraction, 0.4);
//! # Ok::<(), thermopack_rs::Error>(())
//! ```
use crate::{
    Cubic, Error, Phase, Result,
    envelope::{EnvelopeConstraint, Jet},
    error::positive,
    numerics::solve_linear_in_place,
};

#[path = "quality_envelope.rs"]
mod curve;
pub use curve::QualityEnvelope;
#[path = "quality_caloric.rs"]
mod caloric;
pub use caloric::{QualityCaloricProperties, QualitySpecification};

/// Initial temperature, pressure and vapor/liquid mole-fraction ratios.
/// Ratios for absent feed species are ignored. A guess need not close the
/// material balance at the requested quality.
#[derive(Debug, Clone)]
pub struct QualityGuess {
    pub temperature: f64,
    pub pressure: f64,
    pub ratios: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct QualityState {
    pub temperature: f64,
    pub pressure: f64,
    pub vapor_fraction: f64,
    /// Original normalized feed, preserved exactly across warm corrections.
    pub feed: Vec<f64>,
    pub liquid: Vec<f64>,
    pub vapor: Vec<f64>,
    pub liquid_volume: f64,
    pub vapor_volume: f64,
    /// Local dP/dT with feed and molar vapor fraction held fixed.
    /// Infinite at a temperature turning point.
    pub dpdt: f64,
    pub iterations: usize,
    pub fugacity_residual: f64,
    pub material_balance_residual: f64,
}
impl QualityState {
    pub fn guess(&self) -> QualityGuess {
        QualityGuess {
            temperature: self.temperature,
            pressure: self.pressure,
            ratios: self
                .liquid
                .iter()
                .zip(&self.vapor)
                .map(|(&x, &y)| if x > 0. { y / x } else { 1. })
                .collect(),
        }
    }
}

fn check_quality(beta: f64) -> Result<()> {
    if !beta.is_finite() || !(0.0..=1.0).contains(&beta) {
        return Err(Error::InvalidInput(
            "molar vapor fraction in [0,1] required",
        ));
    }
    Ok(())
}
struct Evaluation {
    f: Vec<f64>,
    j: Vec<f64>,
    norm: f64,
    vl: f64,
    vv: f64,
}
impl Evaluation {
    fn new(d: usize) -> Self {
        Self {
            f: vec![0.; d],
            j: vec![0.; d * d],
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
    beta: f64,
    left: Jet<'a>,
    right: Jet<'a>,
    x: Vec<f64>,
    y: Vec<f64>,
    weights: Vec<f64>,
    baseline: bool,
}
impl<'a> Problem<'a> {
    fn new(model: &'a Cubic, z: &[f64], beta: f64) -> Result<Self> {
        check_quality(beta)?;
        let s = model.prepare(300., z)?;
        let z = s.composition().to_vec();
        let n = z.len();
        let active = z
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| (v > 0.).then_some(i))
            .collect();
        Ok(Self {
            model,
            z,
            active,
            beta,
            left: Jet::new(model, n),
            right: Jet::new(model, n),
            x: vec![0.; n],
            y: vec![0.; n],
            weights: vec![0.; n],
            baseline: false,
        })
    }
    fn coordinates(&self, g: &QualityGuess) -> Result<Vec<f64>> {
        positive(g.temperature, "quality seed temperature")?;
        positive(g.pressure, "quality seed pressure")?;
        if g.ratios.len() != self.z.len() {
            return Err(Error::InvalidInput("quality ratio length"));
        }
        let mut u = Vec::with_capacity(self.active.len() + 2);
        for &i in &self.active {
            positive(g.ratios[i], "quality seed ratio")?;
            u.push(g.ratios[i].ln());
        }
        u.extend([g.temperature.ln(), g.pressure.ln()]);
        Ok(u)
    }
    fn constraint(&self, c: impl Into<QualitySpecification>) -> Result<(usize, f64)> {
        let m = self.active.len();
        Ok(match c.into() {
            QualitySpecification::Temperature(t) => {
                positive(t, "quality temperature")?;
                (m, t.ln())
            }
            QualitySpecification::Pressure(p) => {
                positive(p, "quality pressure")?;
                (m + 1, p.ln())
            }
            QualitySpecification::LogK { component, value } => {
                if !value.is_finite() {
                    return Err(Error::InvalidInput("finite quality log K required"));
                }
                let i = self
                    .active
                    .iter()
                    .position(|&i| i == component)
                    .ok_or(Error::InvalidInput("active quality component required"))?;
                (i, value)
            }
            QualitySpecification::Entropy(value) => {
                if !value.is_finite() {
                    return Err(Error::InvalidInput("finite quality entropy required"));
                }
                (m + 2, value)
            }
            QualitySpecification::Enthalpy(value) => {
                if !value.is_finite() {
                    return Err(Error::InvalidInput("finite quality enthalpy required"));
                }
                (m + 3, value)
            }
        })
    }
    fn evaluate(
        &mut self,
        u: &[f64],
        s: usize,
        spec: f64,
        e: &mut Evaluation,
        derivatives: bool,
    ) -> Result<()> {
        let m = self.active.len();
        let d = m + 2;
        let n = self.z.len();
        let t = u[m].exp();
        let p = u[m + 1].exp();
        positive(t, "quality iterate temperature")?;
        positive(p, "quality iterate pressure")?;
        let (lo, hi) = self.model.temperature_bounds();
        if t < lo * (1. - 1e-12)
            || t > hi * (1. + 1e-12)
            || u[..m].iter().any(|v| !v.is_finite() || v.abs() > 600.)
        {
            return Err(Error::NotConverged("quality iterate bounds"));
        }
        for (k, &i) in self.active.iter().enumerate() {
            let ratio = u[k].exp();
            let denom = 1. - self.beta + self.beta * ratio;
            self.x[i] = self.z[i] / denom;
            self.y[i] = ratio * self.x[i];
            self.weights[i] = self.x[i] * self.y[i] / self.z[i];
        }
        let sx = self.x.iter().sum::<f64>();
        let sy = self.y.iter().sum::<f64>();
        positive(sx, "quality liquid normalization")?;
        positive(sy, "quality vapor normalization")?;
        e.f[m] = self.y.iter().zip(&self.x).map(|(y, x)| y - x).sum();
        for x in &mut self.x {
            *x /= sx;
        }
        for y in &mut self.y {
            *y /= sy;
        }
        if self.baseline {
            e.vl = self
                .left
                .fill_baseline(t, p, &self.x, Phase::Liquid, derivatives)?;
            e.vv = self
                .right
                .fill_baseline(t, p, &self.y, Phase::Vapor, derivatives)?;
        } else {
            e.vl = self.left.fill(t, p, &self.x, Phase::Liquid, derivatives)?;
            e.vv = self.right.fill(t, p, &self.y, Phase::Vapor, derivatives)?;
        }
        if derivatives {
            e.j.fill(0.);
        }
        // eos::thermo normalizes its input but returns derivatives at unit
        // total without rescaling them to the original phase amounts. Keep
        // that source Jacobian convention only on the baseline path.
        let (dx_scale, dy_scale) = if self.baseline { (1., 1.) } else { (sx, sy) };
        for (r, &i) in self.active.iter().enumerate() {
            e.f[r] = u[r] + self.right.phi[i] - self.left.phi[i];
            if derivatives {
                for (c, &j) in self.active.iter().enumerate() {
                    e.j[r * d + c] = f64::from(r == c)
                        + self.weights[j]
                            * ((1. - self.beta) * self.right.dn[i * n + j] / dy_scale
                                + self.beta * self.left.dn[i * n + j] / dx_scale);
                }
                e.j[r * d + m] = t * (self.right.dt[i] - self.left.dt[i]);
                e.j[r * d + m + 1] = p * (self.right.dp[i] - self.left.dp[i]);
            }
        }
        if s < d {
            e.f[m + 1] = u[s] - spec;
        } else {
            let entropy = s == d;
            let l = caloric::phase_jet(
                self.model,
                &self.left,
                e.vl,
                p,
                entropy,
                self.baseline,
                derivatives,
            )?;
            let v = caloric::phase_jet(
                self.model,
                &self.right,
                e.vv,
                p,
                entropy,
                self.baseline,
                derivatives,
            )?;
            let scale = spec.abs().max(1.);
            e.f[m + 1] = ((1. - self.beta) * l.value + self.beta * v.value - spec) / scale;
            if derivatives {
                for (k, &i) in self.active.iter().enumerate() {
                    e.j[(m + 1) * d + k] = self.beta * (1. - self.beta) * self.weights[i] / scale
                        * if self.baseline {
                            v.partial[i] + l.partial[i]
                        } else {
                            (v.partial[i] - v.value) / sy - (l.partial[i] - l.value) / sx
                        };
                }
                e.j[(m + 1) * d + m] = t * ((1. - self.beta) * l.dt + self.beta * v.dt) / scale;
                e.j[(m + 1) * d + m + 1] = p * ((1. - self.beta) * l.dp + self.beta * v.dp) / scale;
            }
        }
        if derivatives {
            for (k, &i) in self.active.iter().enumerate() {
                e.j[m * d + k] = self.weights[i];
            }
            if s < d {
                e.j[(m + 1) * d + s] = 1.;
            }
        }
        e.norm = e.f.iter().map(|v| v.abs()).fold(0., f64::max);
        if e.f.iter().any(|v| !v.is_finite()) {
            return Err(Error::NotConverged("nonfinite quality equations"));
        }
        Ok(())
    }
    fn solve(&mut self, mut u: Vec<f64>, s: usize, spec: f64) -> Result<(Vec<f64>, usize)> {
        if self.baseline {
            return self.solve_baseline(u, s, spec);
        }
        let d = self.active.len() + 2;
        let m = d - 2;
        if s < d {
            u[s] = spec;
        }
        let mut e = Evaluation::new(d);
        let mut trial = Evaluation::new(d);
        let mut q = u.clone();
        let mut step = vec![0.; d];
        for iteration in 0..60 {
            self.evaluate(&u, s, spec, &mut e, true)?;
            if e.norm < 1e-10 {
                if u[..m].iter().map(|v| v.abs()).fold(0., f64::max) < 1e-7 {
                    return Err(Error::NotConverged("trivial quality root"));
                }
                return Ok((u, iteration));
            }
            for (v, f) in step.iter_mut().zip(&e.f) {
                *v = -f;
            }
            solve_linear_in_place(&mut e.j, &mut step)?;
            let scale = step
                .iter()
                .enumerate()
                .map(|(i, v)| v.abs() / if i < m { 2. } else { 0.2 })
                .fold(1., f64::max);
            let mut fraction = 1. / scale;
            let mut accepted = false;
            for _ in 0..16 {
                for i in 0..d {
                    q[i] = u[i] + fraction * step[i];
                }
                if s < d {
                    q[s] = spec;
                }
                if self.evaluate(&q, s, spec, &mut trial, false).is_ok() && trial.norm < e.norm {
                    accepted = true;
                    break;
                }
                fraction *= 0.5;
            }
            if !accepted {
                return Err(Error::NotConverged("quality Newton backtracking"));
            }
            std::mem::swap(&mut u, &mut q);
        }
        Err(Error::NotConverged("quality Newton iteration"))
    }
    // sat_newton uses an unscaled Newton step and an L2 test after the
    // update, even if the initial residual is already below its tolerance.
    fn solve_baseline(
        &mut self,
        mut u: Vec<f64>,
        s: usize,
        spec: f64,
    ) -> Result<(Vec<f64>, usize)> {
        let d = self.active.len() + 2;
        let mut e = Evaluation::new(d);
        let mut trial = Evaluation::new(d);
        let mut q = u.clone();
        let mut step = vec![0.; d];
        let trivial = |u: &[f64]| u[..d - 2].iter().all(|v| (v.exp() - 1.).abs() < 1e-8);
        let (tmin, tmax) = self.model.temperature_bounds();
        let epsilon = 5000. * f64::EPSILON;
        let at_limit = |u: &[f64]| {
            u[d - 2] <= tmin.ln() + epsilon
                || u[d - 2] >= tmax.ln() - epsilon
                || u[d - 1] <= 10.0_f64.ln() + epsilon
                || u[d - 1] >= 1e8_f64.ln() - epsilon
                || u[..d - 2]
                    .iter()
                    .any(|&v| v <= f64::MIN_POSITIVE.ln() + epsilon || v >= f64::MAX.ln() - epsilon)
        };
        if at_limit(&u) {
            return Err(Error::NotConverged("source quality iterate at bounds"));
        }
        if trivial(&u) {
            return Err(Error::NotConverged("source trivial quality root"));
        }
        for iteration in 1..500 {
            self.evaluate(&u, s, spec, &mut e, true)?;
            for (v, f) in step.iter_mut().zip(&e.f) {
                *v = -f;
            }
            solve_linear_in_place(&mut e.j, &mut step)?;
            let old = e.f.iter().map(|f| f * f).sum::<f64>().sqrt();
            let mut fraction = 1.;
            let mut norm = f64::INFINITY;
            for _ in 0..=10 {
                for i in 0..d {
                    q[i] = u[i] + fraction * step[i];
                }
                // setXv clips to these bounds and sat_var_at_limit reports
                // ierr=-1 immediately, rather than backtracking into another
                // formal high-pressure solution beyond the source domain.
                if at_limit(&q) {
                    return Err(Error::NotConverged("source quality iterate at bounds"));
                }
                // sat_fun_newton zeroes its equations at K=1; the caller
                // immediately reports ierr=-2. It does not keep searching
                // through that trial with progressively smaller line steps.
                if trivial(&q) {
                    return Err(Error::NotConverged("source trivial quality root"));
                }
                norm = if self.evaluate(&q, s, spec, &mut trial, false).is_ok() {
                    trial.f.iter().map(|f| f * f).sum::<f64>().sqrt()
                } else {
                    f64::INFINITY
                };
                if norm < old {
                    break;
                }
                fraction *= 0.5;
            }
            if !norm.is_finite() {
                return Err(Error::NotConverged("source quality line search"));
            }
            std::mem::swap(&mut u, &mut q);
            if norm < 1e-8 {
                if u[..d - 2].iter().all(|v| (v.exp() - 1.).abs() < 1e-8) {
                    return Err(Error::NotConverged("source trivial quality root"));
                }
                return Ok((u, iteration));
            }
        }
        Err(Error::NotConverged("source quality Newton iteration"))
    }
    fn point(&mut self, u: &[f64], s: usize, iterations: usize) -> Result<QualityState> {
        self.point_at(u, s, u[s], iterations)
    }
    fn point_at(
        &mut self,
        u: &[f64],
        s: usize,
        spec: f64,
        iterations: usize,
    ) -> Result<QualityState> {
        let m = self.active.len();
        let d = m + 2;
        let mut e = Evaluation::new(d);
        self.evaluate(u, s, spec, &mut e, true)?;
        let mut tangent = vec![0.; d];
        tangent[d - 1] = 1.;
        solve_linear_in_place(&mut e.j, &mut tangent)?;
        let t = u[m].exp();
        let p = u[m + 1].exp();
        let material_balance_residual = self
            .z
            .iter()
            .zip(&self.x)
            .zip(&self.y)
            .map(|((&z, &x), &y)| ((1. - self.beta) * x + self.beta * y - z).abs())
            .fold(0., f64::max);
        Ok(QualityState {
            temperature: t,
            pressure: p,
            vapor_fraction: self.beta,
            feed: self.z.clone(),
            liquid: self.x.clone(),
            vapor: self.y.clone(),
            liquid_volume: e.vl,
            vapor_volume: e.vv,
            dpdt: p * tangent[m + 1] / (t * tangent[m]),
            iterations,
            fugacity_residual: e.f[..m].iter().map(|v| v.abs()).fold(0., f64::max),
            material_balance_residual,
        })
    }
}

impl Cubic {
    /// Follow the fork's raw-root, post-update L2 tolerance conventions from
    /// a supplied quality guess. Ill-conditioned near-trivial states can
    /// differ from a tighter physical correction; no global VLE claim is made.
    pub fn baseline_quality_state_with_guess(
        &self,
        z: &[f64],
        beta: f64,
        constraint: EnvelopeConstraint,
        guess: &QualityGuess,
    ) -> Result<QualityState> {
        self.baseline_quality_state_with_specification(z, beta, constraint.into(), guess)
    }
    /// Source sat_newton correction including entropy/enthalpy specifications.
    /// Preserves the fork's caloric Jacobian, including its documented sign
    /// and heat-capacity conventions. Global stability is not tested.
    pub fn baseline_quality_state_with_specification(
        &self,
        z: &[f64],
        beta: f64,
        specification: QualitySpecification,
        guess: &QualityGuess,
    ) -> Result<QualityState> {
        let mut problem = Problem::new(self, z, beta)?;
        problem.baseline = true;
        if problem.active.len() < 2 {
            return Err(Error::InvalidInput(
                "source quality corrector requires mixture",
            ));
        }
        let u = problem.coordinates(guess)?;
        let (s, spec) = problem.constraint(specification)?;
        let (u, it) = problem.solve(u, s, spec)?;
        problem.point_at(&u, s, spec, it)
    }
    /// Local coexistence at specified T, feed and molar vapor fraction.
    pub fn quality_pressure(&self, t: f64, z: &[f64], beta: f64) -> Result<QualityState> {
        self.quality_cold(z, beta, EnvelopeConstraint::Temperature(t))
    }
    /// Local coexistence at specified P, feed and molar vapor fraction.
    pub fn quality_temperature(&self, p: f64, z: &[f64], beta: f64) -> Result<QualityState> {
        self.quality_cold(z, beta, EnvelopeConstraint::Pressure(p))
    }
    /// Correct a supplied log-ratio seed at fixed T, P or one log K.
    /// This preserves access to local/retrograde branches; it does not select
    /// the globally stable split. Absent species remain exactly zero.
    pub fn quality_state_with_guess(
        &self,
        z: &[f64],
        beta: f64,
        constraint: EnvelopeConstraint,
        guess: &QualityGuess,
    ) -> Result<QualityState> {
        self.quality_state_with_specification(z, beta, constraint.into(), guess)
    }
    /// Local correction with a T/P/log-K or molar entropy/enthalpy target.
    /// Caloric targets use J/(mol K) and J/mol, including model reference levels.
    pub fn quality_state_with_specification(
        &self,
        z: &[f64],
        beta: f64,
        specification: QualitySpecification,
        guess: &QualityGuess,
    ) -> Result<QualityState> {
        let mut problem = Problem::new(self, z, beta)?;
        let u = problem.coordinates(guess)?;
        let (s, spec) = problem.constraint(specification)?;
        if problem.active.len() == 1 {
            return match specification {
                QualitySpecification::Temperature(t) => {
                    self.quality_pure(&problem.z, beta, EnvelopeConstraint::Temperature(t))
                }
                QualitySpecification::Pressure(p) => {
                    self.quality_pure(&problem.z, beta, EnvelopeConstraint::Pressure(p))
                }
                QualitySpecification::Entropy(value) => {
                    self.pure_quality_caloric(&problem.z, beta, value, true, guess)
                }
                QualitySpecification::Enthalpy(value) => {
                    self.pure_quality_caloric(&problem.z, beta, value, false, guess)
                }
                _ => Err(Error::InvalidInput(
                    "pure quality requires T/P/H/S specification",
                )),
            };
        }
        let (u, it) = problem.solve(u, s, spec)?;
        problem.point_at(&u, s, spec, it)
    }
    /// Warm continuation at the same feed and quality as a previous state.
    /// The stored feed is preserved, avoiding cumulative material-balance drift.
    pub fn refine_quality_state(
        &self,
        previous: &QualityState,
        constraint: EnvelopeConstraint,
    ) -> Result<QualityState> {
        check_quality(previous.vapor_fraction)?;
        if previous.liquid.len() != previous.vapor.len() {
            return Err(Error::InvalidInput("quality phase composition lengths"));
        }
        self.prepare(previous.temperature, &previous.liquid)?;
        self.prepare(previous.temperature, &previous.vapor)?;
        let beta = previous.vapor_fraction;
        self.quality_state_with_guess(&previous.feed, beta, constraint, &previous.guess())
    }
    fn quality_pure(&self, z: &[f64], beta: f64, c: EnvelopeConstraint) -> Result<QualityState> {
        let sat = match c {
            EnvelopeConstraint::Temperature(t) => self.pseudo_saturation_pressure(t, z)?,
            EnvelopeConstraint::Pressure(p) => self.pseudo_saturation_temperature(p, z)?,
            _ => {
                return Err(Error::InvalidInput(
                    "pure quality requires T or P specification",
                ));
            }
        };
        let state = self.prepare(sat.temperature, z)?;
        let mut fl = vec![0.; z.len()];
        let mut fv = fl.clone();
        state.ln_fugacity_tp_into(sat.pressure, sat.liquid_volume, &mut fl)?;
        state.ln_fugacity_tp_into(sat.pressure, sat.vapor_volume, &mut fv)?;
        let fugacity_residual = z
            .iter()
            .enumerate()
            .filter(|(_, x)| **x > 0.)
            .map(|(i, _)| (fl[i] - fv[i]).abs())
            .fold(0., f64::max);
        Ok(QualityState {
            temperature: sat.temperature,
            pressure: sat.pressure,
            vapor_fraction: beta,
            feed: z.to_vec(),
            liquid: z.to_vec(),
            vapor: z.to_vec(),
            liquid_volume: sat.liquid_volume,
            vapor_volume: sat.vapor_volume,
            dpdt: sat.dpdt,
            iterations: 0,
            fugacity_residual,
            material_balance_residual: 0.,
        })
    }
    fn quality_cold(&self, z: &[f64], beta: f64, c: EnvelopeConstraint) -> Result<QualityState> {
        let mut problem = Problem::new(self, z, beta)?;
        let (s, spec) = problem.constraint(c)?;
        if problem.active.len() == 1 {
            return self.quality_pure(&problem.z, beta, c);
        }
        let bubble = beta < 0.5;
        // The fork initializes its boundary through PureSat. This matters
        // when quantum mixtures admit more than one dense-fluid branch.
        let pseudo = match c {
            EnvelopeConstraint::Temperature(t) => self
                .pseudo_saturation_pressure(t, &problem.z)
                .map(|p| p.pressure),
            EnvelopeConstraint::Pressure(p) => self
                .pseudo_saturation_temperature(p, &problem.z)
                .map(|p| p.temperature),
            _ => Err(Error::InvalidInput(
                "quality initialization requires T or P",
            )),
        }
        .ok();
        let seed = match (c, bubble) {
            (EnvelopeConstraint::Temperature(t), true) => pseudo.map_or_else(
                || self.bubble_pressure(t, &problem.z),
                |g| {
                    self.bubble_pressure_with_guess(t, &problem.z, g)
                        .or_else(|_| self.bubble_pressure(t, &problem.z))
                },
            )?,
            (EnvelopeConstraint::Temperature(t), false) => pseudo.map_or_else(
                || self.dew_pressure(t, &problem.z),
                |g| {
                    self.dew_pressure_with_guess(t, &problem.z, g)
                        .or_else(|_| self.dew_pressure(t, &problem.z))
                },
            )?,
            (EnvelopeConstraint::Pressure(p), true) => pseudo.map_or_else(
                || self.bubble_temperature(p, &problem.z),
                |g| {
                    self.bubble_temperature_with_guess(p, &problem.z, g)
                        .or_else(|_| self.bubble_temperature(p, &problem.z))
                },
            )?,
            (EnvelopeConstraint::Pressure(p), false) => pseudo.map_or_else(
                || self.dew_temperature(p, &problem.z),
                |g| {
                    self.dew_temperature_with_guess(p, &problem.z, g)
                        .or_else(|_| self.dew_temperature(p, &problem.z))
                },
            )?,
            _ => {
                return Err(Error::InvalidInput(
                    "quality initialization requires T or P",
                ));
            }
        };
        let guess = QualityGuess {
            temperature: seed.temperature,
            pressure: seed.pressure,
            ratios: seed
                .liquid
                .iter()
                .zip(&seed.vapor)
                .map(|(&x, &y)| if x > 0. { y / x } else { 1. })
                .collect(),
        };
        let initial = problem.coordinates(&guess)?;
        if let Ok((u, it)) = problem.solve(initial.clone(), s, spec) {
            return problem.point(&u, s, it);
        }
        // Continue quality from the nearest boundary when the direct Newton
        // jump is too large. The specified T/P and original feed stay fixed.
        let origin = if bubble { 0. } else { 1. };
        let mut u = initial;
        let mut iterations = 0;
        for i in 1..=16 {
            problem.beta = origin + (beta - origin) * (i as f64 / 16.);
            let (next, it) = problem.solve(u, s, spec)?;
            u = next;
            iterations += it;
        }
        problem.beta = beta;
        problem.point(&u, s, iterations)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CubicKind, CubicOptions, MixingRule};
    #[test]
    fn caloric_quality_source_rows_and_independent_source_differences_are_retained() {
        let mut count = 0;
        for row in include_str!("../../tests/fixtures/quality_caloric_jacobian.csv")
            .lines()
            .skip(1)
        {
            let c: Vec<_> = row.split(';').collect();
            let f: Vec<_> = c[0].split('|').collect();
            let names: Vec<_> = f[2].split(',').collect();
            let model = if f[3] == "quantum" {
                Cubic::quantum(&names).unwrap()
            } else {
                Cubic::with_options(
                    match f[1] {
                        "SRK" => CubicKind::SoaveRedlichKwong,
                        "RK" => CubicKind::RedlichKwong,
                        "VDW" => CubicKind::VanDerWaals,
                        _ => CubicKind::PengRobinson,
                    },
                    &names,
                    CubicOptions {
                        volume_translation: f[3] == "translated",
                        mixing_rule: match f[3] {
                            "unifac" => MixingRule::Unifac,
                            "hv1" => MixingRule::HuronVidal1,
                            "hv2" => MixingRule::HuronVidal2,
                            "nrtl" => MixingRule::Nrtl,
                            _ => MixingRule::VdW,
                        },
                        ..Default::default()
                    },
                )
                .unwrap()
            };
            let z: Vec<f64> = f[4].split(',').map(|v| v.parse().unwrap()).collect();
            let beta = f[5].parse().unwrap();
            let o: Vec<_> = f[6].split(':').collect();
            let target = o[1].parse().unwrap();
            let constraint = if o[0] == "H" {
                QualitySpecification::Enthalpy(target)
            } else {
                QualitySpecification::Entropy(target)
            };
            let g: Vec<f64> = f[7].split(',').map(|v| v.parse().unwrap()).collect();
            let guess = QualityGuess {
                temperature: g[0],
                pressure: g[1],
                ratios: g[2..].to_vec(),
            };
            let mut problem = Problem::new(&model, &z, beta).unwrap();
            let d = problem.active.len() + 2;
            let u = problem.coordinates(&guess).unwrap();
            let (s, spec) = problem.constraint(constraint).unwrap();
            let source: Vec<f64> = c[1].split(',').map(|v| v.parse().unwrap()).collect();
            let fun: Vec<f64> = c[2].split(',').map(|v| v.parse().unwrap()).collect();
            let fd: Vec<f64> = c[3].split(',').map(|v| v.parse().unwrap()).collect();
            let indices: Vec<_> = problem
                .active
                .iter()
                .copied()
                .chain([z.len(), z.len() + 1])
                .collect();
            let mut e = Evaluation::new(d);
            problem.baseline = true;
            problem.evaluate(&u, s, spec, &mut e, true).unwrap();
            for (i, &source_i) in indices.iter().enumerate() {
                let a = e.f[i];
                let b = fun[source_i];
                assert!(
                    (a - b).abs() < 2e-8 + 2e-7 * a.abs().max(b.abs()),
                    "{} F{i}: {a} vs {b}",
                    c[0]
                );
                let a = e.j[(d - 1) * d + i];
                let b = source[source_i];
                assert!(
                    (a - b).abs() < 1e-8 + 2e-7 * a.abs().max(b.abs()),
                    "{} J{i}: {a} vs {b}",
                    c[0]
                );
            }
            problem.baseline = false;
            problem.evaluate(&u, s, spec, &mut e, true).unwrap();
            for (i, &source_i) in indices.iter().enumerate() {
                let a = e.j[(d - 1) * d + i];
                let b = fd[source_i];
                assert!(
                    (a - b).abs() < 3e-6 + 3e-5 * a.abs().max(b.abs()),
                    "{} finite {i}: {a} vs {b}",
                    c[0]
                );
            }
            count += 1;
        }
        assert_eq!(count, 280);
    }
    #[test]
    fn caloric_quality_rows_differentiate_normalized_phases_off_the_rr_root() {
        for (names, rule, translated, quantum) in [
            (vec!["C1", "C2"], MixingRule::VdW, false, false),
            (vec!["C1", "C3"], MixingRule::Unifac, false, false),
            (vec!["C1", "H2O"], MixingRule::HuronVidal1, false, false),
            (vec!["C1", "C2", "C3"], MixingRule::VdW, true, false),
            (vec!["H2", "NE"], MixingRule::VdW, false, true),
        ] {
            let model = if quantum {
                Cubic::quantum(&names).unwrap()
            } else {
                Cubic::with_options(
                    CubicKind::PengRobinson,
                    &names,
                    CubicOptions {
                        mixing_rule: rule,
                        volume_translation: translated,
                        ..Default::default()
                    },
                )
                .unwrap()
            };
            let mut z = vec![0.5, 0.5];
            if names.len() == 3 {
                z.push(0.);
            }
            for beta in [0., 0.1, 0.4, 0.9, 1.] {
                let mut problem = Problem::new(&model, &z, beta).unwrap();
                let m = problem.active.len();
                let d = m + 2;
                let mut u = vec![0.7, -0.4];
                u.extend([
                    if quantum { 30_f64.ln() } else { 250_f64.ln() },
                    if quantum { 1e5_f64.ln() } else { 1e6_f64.ln() },
                ]);
                for constraint in [
                    QualitySpecification::Enthalpy(-1e5),
                    QualitySpecification::Entropy(100.),
                ] {
                    let (s, spec) = problem.constraint(constraint).unwrap();
                    let mut base = Evaluation::new(d);
                    problem.evaluate(&u, s, spec, &mut base, true).unwrap();
                    for j in 0..d {
                        let h = 2e-6;
                        let mut plus = u.clone();
                        let mut minus = u.clone();
                        plus[j] += h;
                        minus[j] -= h;
                        let mut a = Evaluation::new(d);
                        let mut b = Evaluation::new(d);
                        problem.evaluate(&plus, s, spec, &mut a, false).unwrap();
                        problem.evaluate(&minus, s, spec, &mut b, false).unwrap();
                        let fd = (a.f[d - 1] - b.f[d - 1]) / (2. * h);
                        let exact = base.j[(d - 1) * d + j];
                        assert!(
                            (fd - exact).abs() < 2e-6 + 2e-5 * fd.abs().max(exact.abs()),
                            "{names:?} {beta} {constraint:?} {j}: {exact} vs {fd}"
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn normalized_quality_jacobian_matches_differences_away_from_material_balance() {
        for (names, feed, rule, translated) in [
            (vec!["C1", "C2"], vec![0.5, 0.5], MixingRule::VdW, false),
            (vec!["C1", "C3"], vec![0.5, 0.5], MixingRule::Unifac, false),
            (
                vec!["C1", "C2", "C3"],
                vec![0.2, 0.8, 0.],
                MixingRule::VdW,
                true,
            ),
            (
                vec!["C1", "H2O"],
                vec![0.8, 0.2],
                MixingRule::HuronVidal1,
                false,
            ),
        ] {
            let model = Cubic::with_options(
                CubicKind::PengRobinson,
                &names,
                CubicOptions {
                    mixing_rule: rule,
                    volume_translation: translated,
                    ..Default::default()
                },
            )
            .unwrap();
            for beta in [0., 0.05, 0.4, 0.5, 0.9, 1.] {
                let mut problem = Problem::new(&model, &feed, beta).unwrap();
                let m = problem.active.len();
                let d = m + 2;
                let mut u: Vec<_> = (0..m)
                    .map(|i| if i % 2 == 0 { 0.7 } else { -0.4 })
                    .collect();
                u.extend([250.0_f64.ln(), 1e6_f64.ln()]);
                let mut base = Evaluation::new(d);
                problem.evaluate(&u, m, u[m], &mut base, true).unwrap();
                for j in 0..d {
                    let h = 2e-6;
                    let mut up = u.clone();
                    let mut um = u.clone();
                    up[j] += h;
                    um[j] -= h;
                    let mut a = Evaluation::new(d);
                    let mut b = Evaluation::new(d);
                    problem.evaluate(&up, m, u[m], &mut a, false).unwrap();
                    problem.evaluate(&um, m, u[m], &mut b, false).unwrap();
                    for i in 0..d {
                        let fd = (a.f[i] - b.f[i]) / (2. * h);
                        let exact = base.j[i * d + j];
                        assert!(
                            (fd - exact).abs() < 2e-6 + 2e-5 * fd.abs().max(exact.abs()),
                            "{names:?} beta={beta} ({i},{j}) {exact} vs {fd}"
                        );
                    }
                }
            }
        }
    }
}
