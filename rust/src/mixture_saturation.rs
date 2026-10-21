//! Mixture bubble/dew states from equality of component fugacities.
//! Analytic T/P/mole-number derivatives drive a safeguarded log-variable
//! Newton solve. Zero fixed-phase components remain exactly absent.
//! SI units: K, Pa and m3/mol. These are separated VLE branch calculations;
//! a critical point or global multiphase stability is not implied.
//!
//! ```
//! use thermopack_rs::{Cubic,CubicKind};
//! let model=Cubic::new(CubicKind::PengRobinson,&["C1","C2"])?;
//! let bubble=model.bubble_pressure(180.,&[0.5,0.5])?;
//! assert!(bubble.vapor[0]>bubble.liquid[0]);
//! assert!(bubble.fugacity_residual<1e-9);
//! # Ok::<(),thermopack_rs::Error>(())
//! ```
use crate::{Cubic, Error, Phase, PreparedState, Result, error::positive, numerics};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SaturationKind {
    Bubble,
    Dew,
}

#[derive(Debug, Clone)]
pub struct MixtureSaturation {
    pub kind: SaturationKind,
    pub temperature: f64,
    pub pressure: f64,
    pub liquid: Vec<f64>,
    pub vapor: Vec<f64>,
    pub liquid_volume: f64,
    pub vapor_volume: f64,
    /// dP/dT along this branch with the specified phase composition fixed.
    /// For a mixture this includes the incipient phase's composition response.
    pub dpdt: f64,
    pub iterations: usize,
    pub fugacity_residual: f64,
}

#[derive(Clone, Copy)]
enum Fixed {
    Temperature(f64),
    Pressure(f64),
}
impl Fixed {
    fn conditions(self, variable: f64) -> (f64, f64) {
        match self {
            Self::Temperature(t) => (t, variable),
            Self::Pressure(p) => (variable, p),
        }
    }
}

impl Cubic {
    /// Continue a separated bubble/dew branch from a previous state, keeping
    /// its specified phase composition fixed. Uses its incipient composition
    /// as well as its pressure, so retrograde branches remain accessible.
    pub fn refine_saturation_pressure(
        &self,
        t: f64,
        previous: &MixtureSaturation,
    ) -> Result<MixtureSaturation> {
        self.refine_mixture_saturation(Fixed::Temperature(t), previous)
    }
    pub fn refine_saturation_temperature(
        &self,
        p: f64,
        previous: &MixtureSaturation,
    ) -> Result<MixtureSaturation> {
        self.refine_mixture_saturation(Fixed::Pressure(p), previous)
    }
    fn refine_mixture_saturation(
        &self,
        fixed: Fixed,
        previous: &MixtureSaturation,
    ) -> Result<MixtureSaturation> {
        let (parent, trial) = if previous.kind == SaturationKind::Bubble {
            (&previous.liquid, &previous.vapor)
        } else {
            (&previous.vapor, &previous.liquid)
        };
        let initial = match fixed {
            Fixed::Temperature(_) => previous.pressure,
            Fixed::Pressure(_) => previous.temperature,
        };
        let (t, p) = fixed.conditions(initial);
        positive(t, "saturation temperature")?;
        positive(p, "saturation pressure")?;
        let prepared = self.prepare(t, parent)?;
        if trial.len() != parent.len()
            || trial.iter().any(|v| !v.is_finite() || *v < 0.)
            || (trial.iter().sum::<f64>() - 1.).abs() > 1e-10
            || trial
                .iter()
                .zip(parent)
                .any(|(a, b)| (*a == 0.) != (*b == 0.))
        {
            return Err(Error::InvalidInput(
                "normalized incipient seed with matching zero components required",
            ));
        }
        let active: Vec<_> = parent
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| if v > 0. { Some(i) } else { None })
            .collect();
        let prepared = if matches!(fixed, Fixed::Temperature(_)) {
            Some(prepared)
        } else {
            None
        };
        let mut problem = Problem::new(self, parent, &active, previous.kind, fixed, prepared);
        problem.select_onset = false;
        problem.solve_seeded(initial, Some(trial))
    }
    pub fn bubble_pressure(&self, t: f64, liquid: &[f64]) -> Result<MixtureSaturation> {
        self.mixture_saturation(Fixed::Temperature(t), liquid, SaturationKind::Bubble, None)
    }
    pub fn dew_pressure(&self, t: f64, vapor: &[f64]) -> Result<MixtureSaturation> {
        self.mixture_saturation(Fixed::Temperature(t), vapor, SaturationKind::Dew, None)
    }
    pub fn bubble_temperature(&self, p: f64, liquid: &[f64]) -> Result<MixtureSaturation> {
        self.mixture_saturation(Fixed::Pressure(p), liquid, SaturationKind::Bubble, None)
    }
    pub fn dew_temperature(&self, p: f64, vapor: &[f64]) -> Result<MixtureSaturation> {
        self.mixture_saturation(Fixed::Pressure(p), vapor, SaturationKind::Dew, None)
    }
    /// Explicit pressure seed in Pa; other automatic seeds remain fallbacks.
    pub fn bubble_pressure_with_guess(
        &self,
        t: f64,
        liquid: &[f64],
        p: f64,
    ) -> Result<MixtureSaturation> {
        self.mixture_saturation(
            Fixed::Temperature(t),
            liquid,
            SaturationKind::Bubble,
            Some(p),
        )
    }
    pub fn dew_pressure_with_guess(
        &self,
        t: f64,
        vapor: &[f64],
        p: f64,
    ) -> Result<MixtureSaturation> {
        self.mixture_saturation(Fixed::Temperature(t), vapor, SaturationKind::Dew, Some(p))
    }
    /// Explicit temperature seed in K; other automatic seeds remain fallbacks.
    pub fn bubble_temperature_with_guess(
        &self,
        p: f64,
        liquid: &[f64],
        t: f64,
    ) -> Result<MixtureSaturation> {
        self.mixture_saturation(Fixed::Pressure(p), liquid, SaturationKind::Bubble, Some(t))
    }
    pub fn dew_temperature_with_guess(
        &self,
        p: f64,
        vapor: &[f64],
        t: f64,
    ) -> Result<MixtureSaturation> {
        self.mixture_saturation(Fixed::Pressure(p), vapor, SaturationKind::Dew, Some(t))
    }

    fn mixture_saturation(
        &self,
        fixed: Fixed,
        z: &[f64],
        kind: SaturationKind,
        guess: Option<f64>,
    ) -> Result<MixtureSaturation> {
        positive(
            match fixed {
                Fixed::Temperature(t) => t,
                Fixed::Pressure(p) => p,
            },
            "saturation specification",
        )?;
        if let Some(g) = guess {
            positive(g, "saturation initial guess")?;
        }
        if z.len() != self.components().len()
            || z.iter().any(|v| !v.is_finite() || *v < 0.)
            || (z.iter().sum::<f64>() - 1.).abs() > 1e-10
        {
            return Err(Error::InvalidInput(
                "normalized saturation phase composition required",
            ));
        }
        let total = z.iter().sum::<f64>();
        let z: Vec<_> = z.iter().map(|v| v / total).collect();
        let active: Vec<_> = z
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| if v > 0. { Some(i) } else { None })
            .collect();
        if active.len() == 1 {
            let s = match fixed {
                Fixed::Temperature(t) => self.pseudo_saturation_pressure(t, &z)?,
                Fixed::Pressure(p) => self.pseudo_saturation_temperature(p, &z)?,
            };
            let prepared = self.prepare(s.temperature, &z)?;
            let mut fl = vec![0.; z.len()];
            let mut fv = fl.clone();
            prepared.ln_fugacity_tp_into(s.pressure, s.liquid_volume, &mut fl)?;
            prepared.ln_fugacity_tp_into(s.pressure, s.vapor_volume, &mut fv)?;
            let fugacity_residual = (fl[active[0]] - fv[active[0]]).abs();
            return Ok(MixtureSaturation {
                kind,
                temperature: s.temperature,
                pressure: s.pressure,
                liquid: z.clone(),
                vapor: z,
                liquid_volume: s.liquid_volume,
                vapor_volume: s.vapor_volume,
                dpdt: s.dpdt,
                iterations: 0,
                fugacity_residual,
            });
        }
        let automatic = match fixed {
            Fixed::Temperature(t) => self.wilson_saturation_pressure(t, &z, kind),
            Fixed::Pressure(p) => self.wilson_saturation_temperature(p, &z, kind),
        };
        let seed = match (automatic, guess) {
            (Ok(s), _) => s,
            (Err(_), Some(g)) => g,
            (Err(e), None) => return Err(e),
        };
        let factors: &[f64] = match fixed {
            Fixed::Temperature(_) => &[1., 0.3, 3., 0.03, 10.],
            Fixed::Pressure(_) => &[1., 0.8, 1.2, 0.5, 1.5],
        };
        let mut seeds = Vec::with_capacity(7);
        if let Some(g) = guess {
            seeds.push(g);
        }
        seeds.extend(factors.iter().map(|f| seed * f));
        let fixed_prepared = match fixed {
            Fixed::Temperature(t) => Some(self.prepare(t, &z)?),
            Fixed::Pressure(_) => None,
        };
        let mut problem = Problem::new(self, &z, &active, kind, fixed, fixed_prepared);
        problem.select_onset = guess.is_none();
        for initial in seeds {
            if let Ok(s) = problem.solve(initial) {
                return Ok(s);
            }
        }
        // The source safe_* wrappers use pseudo-pure saturation as a seed.
        // Keep this more expensive initialization as a fallback after Wilson.
        let pseudo = match fixed {
            Fixed::Temperature(t) => self.pseudo_saturation_pressure(t, &z).map(|s| s.pressure),
            Fixed::Pressure(p) => self
                .pseudo_saturation_temperature(p, &z)
                .map(|s| s.temperature),
        };
        if let Ok(seed) = pseudo
            && let Ok(s) = problem.solve(seed)
        {
            return Ok(s);
        }
        if let Fixed::Pressure(p) = fixed {
            // The incipient composition can have a small Newton basin at
            // high pressure. Seed it from a pressure boundary at a nearby T,
            // then impose the exact requested pressure with the coupled solve.
            for factor in [
                1., 0.98, 1.02, 0.95, 1.05, 0.9, 1.1, 0.8, 1.2, 0.7, 1.3, 0.6, 1.4,
            ] {
                if let Ok(mut boundary) =
                    self.mixture_saturation(Fixed::Temperature(seed * factor), &z, kind, Some(p))
                {
                    for _ in 0..30 {
                        let error = (boundary.pressure / p).ln().abs();
                        if error < 0.01 {
                            let trial = if kind == SaturationKind::Bubble {
                                &boundary.vapor
                            } else {
                                &boundary.liquid
                            };
                            if let Ok(s) = problem.solve_seeded(boundary.temperature, Some(trial)) {
                                return Ok(s);
                            }
                        }
                        let dt = ((p - boundary.pressure) / boundary.dpdt)
                            .clamp(-0.1 * boundary.temperature, 0.1 * boundary.temperature);
                        if !dt.is_finite() || dt == 0. {
                            break;
                        }
                        let mut fraction = 1.;
                        let mut next = None;
                        for _ in 0..14 {
                            let t = boundary.temperature + fraction * dt;
                            if let Ok(candidate) = self.mixture_saturation(
                                Fixed::Temperature(t),
                                &z,
                                kind,
                                Some(boundary.pressure),
                            ) && (candidate.pressure / p).ln().abs() < error
                            {
                                next = Some(candidate);
                                break;
                            }
                            fraction *= 0.5;
                        }
                        if let Some(candidate) = next {
                            boundary = candidate;
                        } else {
                            break;
                        }
                    }
                }
            }
        }
        Err(Error::NotConverged("mixture bubble/dew saturation"))
    }

    fn wilson_log_pressure(&self, t: f64, i: usize) -> f64 {
        let c = &self.components()[i];
        c.critical_pressure.ln()
            + 5.373 * (1. + c.acentric_factor) * (1. - c.critical_temperature / t)
    }
    fn wilson_saturation_pressure(&self, t: f64, z: &[f64], kind: SaturationKind) -> Result<f64> {
        let sign = if kind == SaturationKind::Bubble {
            1.
        } else {
            -1.
        };
        let logs: Vec<_> = z
            .iter()
            .enumerate()
            .filter(|(_, v)| **v > 0.)
            .map(|(i, v)| v.ln() + sign * self.wilson_log_pressure(t, i))
            .collect();
        let p = (sign * log_sum(&logs)).exp();
        positive(p, "finite Wilson saturation pressure")?;
        Ok(p)
    }
    fn wilson_saturation_temperature(
        &self,
        p: f64,
        z: &[f64],
        kind: SaturationKind,
    ) -> Result<f64> {
        let mut low = self
            .components()
            .iter()
            .map(|c| c.critical_temperature)
            .fold(f64::INFINITY, f64::min)
            * 0.05;
        let mut high = self
            .components()
            .iter()
            .map(|c| c.critical_temperature)
            .fold(0., f64::max)
            * 4.;
        let value = |t| {
            self.wilson_saturation_pressure(t, z, kind)
                .map(|q| (q / p).ln())
        };
        // Both Wilson bubble and dew pressure are increasing with T.
        for _ in 0..12 {
            if value(low).is_ok_and(|v| v < 0.) && value(high).is_ok_and(|v| v > 0.) {
                break;
            }
            low *= 0.5;
            high *= 2.;
        }
        if !value(low).is_ok_and(|v| v < 0.) || !value(high).is_ok_and(|v| v > 0.) {
            return Err(Error::NoPhysicalRoot);
        }
        for _ in 0..80 {
            let mid = (low * high).sqrt();
            if value(mid)? > 0. {
                high = mid;
            } else {
                low = mid;
            }
        }
        Ok((low * high).sqrt())
    }
}

fn log_sum(logs: &[f64]) -> f64 {
    let max = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    max + logs.iter().map(|v| (v - max).exp()).sum::<f64>().ln()
}

struct Evaluation {
    residual: Vec<f64>,
    jacobian: Vec<f64>,
    other: Vec<f64>,
    norm: f64,
    vl: f64,
    vv: f64,
    sum: f64,
}
impl Evaluation {
    fn new(n: usize) -> Self {
        Self {
            residual: vec![0.; n],
            jacobian: vec![0.; n * n],
            other: vec![0.; n],
            norm: 0.,
            vl: 0.,
            vv: 0.,
            sum: 0.,
        }
    }
}
struct Problem<'a> {
    select_onset: bool,
    model: &'a Cubic,
    z: &'a [f64],
    active: &'a [usize],
    kind: SaturationKind,
    fixed: Fixed,
    fixed_prepared: Option<PreparedState<'a>>,
    weights: Vec<f64>,
    composition: Vec<f64>,
    phi_f: Vec<f64>,
    phi_t: Vec<f64>,
    dt_f: Vec<f64>,
    dp_f: Vec<f64>,
    dt_t: Vec<f64>,
    dp_t: Vec<f64>,
    dn_f: Vec<f64>,
    dn_t: Vec<f64>,
}
impl<'a> Problem<'a> {
    fn new(
        model: &'a Cubic,
        z: &'a [f64],
        active: &'a [usize],
        kind: SaturationKind,
        fixed: Fixed,
        fixed_prepared: Option<PreparedState<'a>>,
    ) -> Self {
        let n = z.len();
        Self {
            select_onset: true,
            model,
            z,
            active,
            kind,
            fixed,
            fixed_prepared,
            weights: vec![0.; n],
            composition: vec![0.; n],
            phi_f: vec![0.; n],
            phi_t: vec![0.; n],
            dt_f: vec![0.; n],
            dp_f: vec![0.; n],
            dt_t: vec![0.; n],
            dp_t: vec![0.; n],
            dn_f: vec![0.; n * n],
            dn_t: vec![0.; n * n],
        }
    }
    fn evaluate(&mut self, u: &[f64], out: &mut Evaluation, derivatives: bool) -> Result<()> {
        let m = self.active.len();
        let dim = m + 1;
        let variable = u[m].exp();
        positive(variable, "finite saturation iterate")?;
        let (t, p) = self.fixed.conditions(variable);
        for (k, &i) in self.active.iter().enumerate() {
            self.weights[i] = u[k].exp();
            positive(self.weights[i], "finite incipient phase amount")?;
        }
        let sum = self.weights.iter().sum::<f64>();
        positive(sum, "incipient phase normalization")?;
        for (i, v) in self.weights.iter().enumerate() {
            self.composition[i] = v / sum;
        }
        let dynamic;
        let f = if let Some(prepared) = self.fixed_prepared.as_ref() {
            prepared
        } else {
            dynamic = self.model.prepare(t, self.z)?;
            &dynamic
        };
        let trial = self.model.prepare(t, &self.composition)?;
        let (pf, pt) = if self.kind == SaturationKind::Bubble {
            (Phase::Liquid, Phase::Vapor)
        } else {
            (Phase::Vapor, Phase::Liquid)
        };
        let vf = f.volume(p, pf)?;
        let vt = trial.volume(p, pt)?;
        f.ln_fugacity_tp_into(p, vf, &mut self.phi_f)?;
        trial.ln_fugacity_tp_into(p, vt, &mut self.phi_t)?;
        if derivatives {
            f.fugacity_derivatives_into(vf, &mut self.dt_f, &mut self.dp_f, &mut self.dn_f)?;
            trial.fugacity_derivatives_into(vt, &mut self.dt_t, &mut self.dp_t, &mut self.dn_t)?;
            out.jacobian.fill(0.);
            out.other.fill(0.);
        }
        for (row, &i) in self.active.iter().enumerate() {
            out.residual[row] = self.z[i].ln() + self.phi_f[i] - u[row] - self.phi_t[i];
            if derivatives {
                for (col, &j) in self.active.iter().enumerate() {
                    out.jacobian[row * dim + col] = -f64::from(row == col)
                        - self.dn_t[i * self.z.len() + j] * self.composition[j];
                }
                match self.fixed {
                    Fixed::Temperature(_) => {
                        out.jacobian[row * dim + m] = p * (self.dp_f[i] - self.dp_t[i]);
                        out.other[row] = self.dt_f[i] - self.dt_t[i];
                    }
                    Fixed::Pressure(_) => {
                        out.jacobian[row * dim + m] = t * (self.dt_f[i] - self.dt_t[i]);
                        out.other[row] = self.dp_f[i] - self.dp_t[i];
                    }
                }
            }
        }
        out.residual[m] = sum - 1.;
        if derivatives {
            for (col, &i) in self.active.iter().enumerate() {
                out.jacobian[m * dim + col] = self.weights[i];
            }
        }
        out.norm = out.residual.iter().map(|r| r.abs()).fold(0., f64::max);
        if !out.norm.is_finite() {
            return Err(Error::NotConverged("non-finite saturation residual"));
        }
        (out.vl, out.vv) = if self.kind == SaturationKind::Bubble {
            (vf, vt)
        } else {
            (vt, vf)
        };
        out.sum = sum;
        Ok(())
    }
    fn solve(&mut self, initial: f64) -> Result<MixtureSaturation> {
        self.solve_seeded(initial, None)
    }
    fn solve_seeded(
        &mut self,
        initial: f64,
        incipient: Option<&[f64]>,
    ) -> Result<MixtureSaturation> {
        let m = self.active.len();
        let dim = m + 1;
        let (t, p) = self.fixed.conditions(initial);
        let sign = if self.kind == SaturationKind::Bubble {
            1.
        } else {
            -1.
        };
        let mut u: Vec<_> = self
            .active
            .iter()
            .map(|&i| {
                if let Some(composition) = incipient {
                    composition[i].ln()
                } else {
                    self.z[i].ln() + sign * (self.model.wilson_log_pressure(t, i) - p.ln())
                }
            })
            .collect();
        let sum = log_sum(&u);
        for v in &mut u {
            *v -= sum;
        }
        u.push(initial.ln());
        let mut current = Evaluation::new(dim);
        let mut trial = Evaluation::new(dim);
        let mut matrix = vec![0.; dim * dim];
        let mut step = vec![0.; dim];
        let mut next = vec![0.; dim];
        for iteration in 0..100 {
            self.evaluate(&u, &mut current, true)?;
            matrix.copy_from_slice(&current.jacobian);
            for (s, r) in step.iter_mut().zip(&current.residual) {
                *s = -r;
            }
            numerics::solve_linear_in_place(&mut matrix, &mut step)?;
            let correction = step.iter().map(|v| v.abs()).fold(0., f64::max);
            if current.norm < 1e-10 && correction < 1e-8 {
                let (t, p) = self.fixed.conditions(u[m].exp());
                if current.vv <= current.vl * (1. + 1e-7) {
                    return Err(Error::NoPhysicalRoot);
                }
                let phase_gap = (current.vv - current.vl) / current.vv;
                let composition_gap = self
                    .z
                    .iter()
                    .zip(&self.composition)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0., f64::max);
                if phase_gap < 1e-4 && composition_gap < 1e-4 {
                    return Err(Error::NotConverged(
                        "coalescing saturation phases require critical-point treatment",
                    ));
                }
                if self.select_onset {
                    // At the first onset, trial-phase TPD decreases on moving
                    // away from the stable parent phase. Its parameter
                    // derivative is the trial-weighted fugacity difference;
                    // Gibbs-Duhem removes the incipient composition response.
                    let slope = self
                        .active
                        .iter()
                        .map(|&i| {
                            self.composition[i]
                                * match self.fixed {
                                    Fixed::Temperature(_) => p * (self.dp_f[i] - self.dp_t[i]),
                                    Fixed::Pressure(_) => t * (self.dt_f[i] - self.dt_t[i]),
                                }
                        })
                        .sum::<f64>();
                    let expected = match (self.kind, self.fixed) {
                        (SaturationKind::Bubble, Fixed::Temperature(_))
                        | (SaturationKind::Dew, Fixed::Pressure(_)) => -1.,
                        _ => 1.,
                    };
                    if expected * slope <= 1e-12 {
                        return Err(Error::NoPhysicalRoot);
                    }
                }
                matrix.copy_from_slice(&current.jacobian);
                for (s, r) in step.iter_mut().zip(&current.other) {
                    *s = -r;
                }
                numerics::solve_linear_in_place(&mut matrix, &mut step)?;
                let (t, p) = self.fixed.conditions(u[m].exp());
                let dpdt = match self.fixed {
                    Fixed::Temperature(_) => p * step[m],
                    Fixed::Pressure(_) => 1. / (t * step[m]),
                };
                if !dpdt.is_finite() {
                    return Err(Error::NotConverged("singular saturation tangent"));
                }
                let incipient = self.composition.clone();
                let (liquid, vapor) = if self.kind == SaturationKind::Bubble {
                    (self.z.to_vec(), incipient)
                } else {
                    (incipient, self.z.to_vec())
                };
                let fugacity_residual = current.residual[..m]
                    .iter()
                    .map(|r| (r + current.sum.ln()).abs())
                    .fold(0., f64::max);
                return Ok(MixtureSaturation {
                    kind: self.kind,
                    temperature: t,
                    pressure: p,
                    liquid,
                    vapor,
                    liquid_volume: current.vl,
                    vapor_volume: current.vv,
                    dpdt,
                    iterations: iteration,
                    fugacity_residual,
                });
            }
            let scale = step.iter().map(|v| v.abs()).fold(0., f64::max);
            if scale > 2. {
                for s in &mut step {
                    *s *= 2. / scale;
                }
            }
            let mut fraction = 1.;
            let mut accepted = false;
            for _ in 0..24 {
                for ((v, a), d) in next.iter_mut().zip(&u).zip(&step) {
                    *v = a + fraction * d;
                }
                if self.evaluate(&next, &mut trial, false).is_ok()
                    && trial.norm < current.norm * (1. - 1e-4 * fraction)
                {
                    u.copy_from_slice(&next);
                    accepted = true;
                    break;
                }
                fraction *= 0.5;
            }
            if !accepted {
                return Err(Error::NotConverged("saturation Newton line search"));
            }
        }
        Err(Error::NotConverged("saturation Newton iteration limit"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn analytic_saturation_jacobians_match_independent_log_variable_differences() {
        let model = Cubic::new(crate::CubicKind::PengRobinson, &["C1", "C2", "C3"]).unwrap();
        for kind in [SaturationKind::Bubble, SaturationKind::Dew] {
            let z = if kind == SaturationKind::Bubble {
                [0.5, 0.5, 0.]
            } else {
                [0.9, 0.1, 0.]
            };
            let trial = if kind == SaturationKind::Bubble {
                [0.9_f64, 0.1]
            } else {
                [0.5_f64, 0.5]
            };
            for fixed in [Fixed::Temperature(180.), Fixed::Pressure(1e6)] {
                let fixed_prepared = match fixed {
                    Fixed::Temperature(t) => Some(model.prepare(t, &z).unwrap()),
                    _ => None,
                };
                let mut problem = Problem::new(&model, &z, &[0, 1], kind, fixed, fixed_prepared);
                let variable: f64 = match fixed {
                    Fixed::Temperature(_) => 1e6,
                    Fixed::Pressure(_) => 180.,
                };
                let u = [trial[0].ln(), trial[1].ln(), variable.ln()];
                let mut center = Evaluation::new(3);
                problem.evaluate(&u, &mut center, true).unwrap();
                for j in 0..3 {
                    let mut hi = u;
                    let mut lo = u;
                    let step = 1e-5;
                    hi[j] += step;
                    lo[j] -= step;
                    let mut a = Evaluation::new(3);
                    let mut b = Evaluation::new(3);
                    problem.evaluate(&hi, &mut a, false).unwrap();
                    problem.evaluate(&lo, &mut b, false).unwrap();
                    for i in 0..3 {
                        let finite = (a.residual[i] - b.residual[i]) / (2. * step);
                        let exact = center.jacobian[i * 3 + j];
                        assert!(
                            (finite - exact).abs() < 1e-7 + 1e-7 * exact.abs(),
                            "{kind:?} {i},{j}: {finite} != {exact}"
                        );
                    }
                }
            }
        }
    }
}
