//! Native VLE stability and TP flashes, following stability.f90/tp_solver.f90.
//! Negative tangent-plane trials seed a Gibbs minimization with analytic
//! mole-number Hessians. Zero-feed components remain exactly absent.
use crate::error::positive;
use crate::{Cubic, Error, Phase, Result, numerics};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum FlashPhase {
    TwoPhase = 0,
    Liquid = 1,
    Vapor = 2,
    Single = 4,
}

#[derive(Debug, Clone, Copy)]
pub struct FlashOptions {
    pub max_iterations: usize,
    pub stability_iterations: usize,
    pub fugacity_tolerance: f64,
}
impl Default for FlashOptions {
    fn default() -> Self {
        Self {
            max_iterations: 200,
            stability_iterations: 200,
            fugacity_tolerance: 1e-10,
        }
    }
}

#[derive(Debug, Clone)]
pub struct StabilityTrial {
    /// Michelsen modified TPD objective before normalizing trial amounts.
    /// Negative values prove instability; composition below is normalized.
    pub tpd: f64,
    pub composition: Vec<f64>,
    pub ln_fugacity: Vec<f64>,
    pub iterations: usize,
}

#[derive(Debug, Clone)]
pub struct Stability {
    /// Stability against the liquid-like and vapor-like searches used by TP.
    /// This does not qualify a complete multiphase global-stability search.
    pub stable: bool,
    pub feed_phase: FlashPhase,
    pub liquid_trial: StabilityTrial,
    pub vapor_trial: StabilityTrial,
}

#[derive(Debug, Clone)]
pub struct FlashResult {
    pub temperature: f64,
    pub pressure: f64,
    pub feed: Vec<f64>,
    pub phase: FlashPhase,
    pub liquid: Vec<f64>,
    pub vapor: Vec<f64>,
    /// Preserves the fork's -1 sentinel when phase is Single (undefined type).
    pub vapor_fraction: f64,
    pub liquid_fraction: f64,
    pub liquid_volume: Option<f64>,
    pub vapor_volume: Option<f64>,
    /// Primary split/temperature iterations; internal initial estimates and
    /// separate stability searches are excluded.
    pub iterations: usize,
    pub fugacity_residual: f64,
}

struct Evaluation {
    phi: Vec<f64>,
    dn: Vec<f64>,
    volume: f64,
}

struct SplitEvaluation {
    gibbs: f64,
    vapor_fraction: f64,
    liquid: Vec<f64>,
    vapor: Vec<f64>,
    liquid_properties: Evaluation,
    vapor_properties: Evaluation,
}

fn evaluate(
    model: &Cubic,
    t: f64,
    p: f64,
    x: &[f64],
    phase: Phase,
    jacobian: bool,
) -> Result<Evaluation> {
    let state = model.prepare(t, x)?;
    let volume = state.volume(p, phase)?;
    let mut phi = vec![0.0; x.len()];
    state.ln_fugacity_tp_into(p, volume, &mut phi)?;
    let mut dn = Vec::new();
    if jacobian {
        let mut dt = vec![0.0; x.len()];
        let mut dp = vec![0.0; x.len()];
        dn.resize(x.len() * x.len(), 0.0);
        state.fugacity_derivatives_into(volume, &mut dt, &mut dp, &mut dn)?;
    }
    Ok(Evaluation { phi, dn, volume })
}

fn normalize(values: &mut [f64]) -> Result<f64> {
    let sum: f64 = values.iter().sum();
    positive(sum, "positive finite total phase amount")?;
    for value in values {
        *value /= sum;
    }
    Ok(sum)
}

fn gibbs(x: &[f64], phi: &[f64]) -> f64 {
    x.iter()
        .zip(phi)
        .filter(|(x, _)| **x > 0.0)
        .map(|(x, f)| x * (x.ln() + f))
        .sum()
}

impl Cubic {
    /// Wilson equilibrium estimates and their standard critical-data formula.
    pub fn wilson_k(&self, t: f64, p: f64) -> Result<Vec<f64>> {
        positive(t, "temperature")?;
        positive(p, "pressure")?;
        self.components()
            .iter()
            .map(|c| {
                let k = (c.critical_pressure / p)
                    * (5.373 * (1.0 + c.acentric_factor) * (1.0 - c.critical_temperature / t))
                        .exp();
                if !k.is_finite() || k < 0. {
                    return Err(Error::InvalidInput(
                        "Wilson K outside finite nonnegative range",
                    ));
                }
                Ok(k)
            })
            .collect()
    }

    pub fn stability(&self, t: f64, p: f64, z: &[f64]) -> Result<Stability> {
        self.stability_with_options(t, p, z, FlashOptions::default())
            .map(|s| s.0)
    }

    fn stability_with_options(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        options: FlashOptions,
    ) -> Result<(Stability, Evaluation)> {
        positive(p, "pressure")?;
        let state = self.prepare(t, z)?;
        let (roots, count) = state.volumes(p)?;
        let mut best = roots[0];
        let mut phase = FlashPhase::Single;
        if count > 1 {
            let l = state.residual(roots[0])?;
            let v = state.residual(roots[count - 1])?;
            let gl = l.enthalpy - t * l.entropy;
            let gv = v.enthalpy - t * v.entropy;
            if gl < gv {
                phase = FlashPhase::Liquid;
            } else {
                phase = FlashPhase::Vapor;
                best = roots[count - 1];
            }
        }
        let mut feed_phi = vec![0.0; z.len()];
        state.ln_fugacity_tp_into(p, best, &mut feed_phi)?;
        let feed = Evaluation {
            phi: feed_phi,
            dn: Vec::new(),
            volume: best,
        };
        let mut k = self.wilson_k(t, p)?;
        for (i, c) in self.components().iter().enumerate() {
            if c.ident == "HE" {
                k[i] = (226970.0 / p) * (5.373 * (1.0 - 5.19 / t)).exp();
            }
        }
        // A Wilson estimate can underflow to zero at valid low temperatures.
        // Keep this numerical seed positive; the public estimate retains zero.
        for value in &mut k {
            *value = value.clamp(1e-300, 1e300);
        }
        let liquid = self.stability_trial(t, p, z, &feed.phi, &k, Phase::Liquid, options)?;
        let vapor = self.stability_trial(t, p, z, &feed.phi, &k, Phase::Vapor, options)?;
        let limit = -1000.0 * f64::EPSILON;
        let stable = liquid.tpd >= limit && vapor.tpd >= limit;
        Ok((
            Stability {
                stable,
                feed_phase: phase,
                liquid_trial: liquid,
                vapor_trial: vapor,
            },
            feed,
        ))
    }

    #[allow(clippy::too_many_arguments)]
    fn stability_trial(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        feed_phi: &[f64],
        k: &[f64],
        phase: Phase,
        options: FlashOptions,
    ) -> Result<StabilityTrial> {
        let n = z.len();
        let active: Vec<_> = (0..n).filter(|&i| z[i] > 0.0).collect();
        let d: Vec<_> = z
            .iter()
            .zip(feed_phi)
            .map(|(z, f)| if *z > 0.0 { z.ln() + f } else { 0.0 })
            .collect();
        let mut w: Vec<_> = z
            .iter()
            .zip(k)
            .map(|(z, k)| if phase == Phase::Vapor { z * k } else { z / k })
            .collect();
        let trial_evaluate = |w: &[f64], jacobian: bool| -> Result<(f64, Vec<f64>, Evaluation)> {
            let mut x = w.to_vec();
            normalize(&mut x)?;
            let e = evaluate(self, t, p, &x, phase, jacobian)?;
            let mut objective = 1.0;
            let mut residual = vec![0.0; n];
            for &i in &active {
                residual[i] = w[i].ln() + e.phi[i] - d[i];
                objective += w[i] * (residual[i] - 1.0);
            }
            Ok((objective, residual, e))
        };
        for iteration in 0..options.stability_iterations {
            let (objective, residual, e) = trial_evaluate(&w, iteration >= 3)?;
            let error = active
                .iter()
                .map(|&i| residual[i].abs())
                .fold(0.0, f64::max);
            if objective < -0.001 || error < 1e-10 {
                let mut composition = w.clone();
                normalize(&mut composition)?;
                return Ok(StabilityTrial {
                    tpd: objective,
                    composition,
                    ln_fugacity: e.phi,
                    iterations: iteration + 1,
                });
            }
            if iteration < 3 {
                for &i in &active {
                    w[i] = (d[i] - e.phi[i]).exp();
                }
                continue;
            }
            let m = active.len();
            let total: f64 = w.iter().sum();
            let gradient: Vec<_> = active.iter().map(|&i| w[i].sqrt() * residual[i]).collect();
            let mut hessian = vec![0.0; m * m];
            for (ii, &i) in active.iter().enumerate() {
                for (jj, &j) in active.iter().enumerate() {
                    hessian[ii * m + jj] = e.dn[i * n + j] * w[i].sqrt() * w[j].sqrt() / total;
                }
                hessian[ii * m + ii] += 1.0;
            }
            let step = numerics::modified_newton_direction(&hessian, &gradient);
            let slope: f64 = gradient.iter().zip(&step).map(|(g, s)| g * s).sum();
            let mut length = 1.0_f64;
            for (ii, &i) in active.iter().enumerate() {
                if step[ii] < 0.0 {
                    length = length.min(-0.99 * 2.0 * w[i].sqrt() / step[ii]);
                }
            }
            let mut accepted = false;
            for _ in 0..40 {
                let mut candidate = vec![0.0; n];
                for (ii, &i) in active.iter().enumerate() {
                    candidate[i] = (2.0 * w[i].sqrt() + length * step[ii]).powi(2) / 4.0;
                }
                if let Ok((value, _, _)) = trial_evaluate(&candidate, false)
                    && value <= objective + 1e-4 * length * slope + 4e-15 * objective.abs().max(1.0)
                {
                    w = candidate;
                    accepted = true;
                    break;
                }
                length *= 0.5;
            }
            if !accepted {
                return Err(Error::NotConverged("tangent-plane line search"));
            }
        }
        Err(Error::NotConverged("tangent-plane minimization"))
    }

    pub fn tp_flash(&self, t: f64, p: f64, z: &[f64]) -> Result<FlashResult> {
        self.tp_flash_with_options(t, p, z, FlashOptions::default())
    }

    pub fn tp_flash_with_options(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        options: FlashOptions,
    ) -> Result<FlashResult> {
        if options.max_iterations == 0
            || options.stability_iterations == 0
            || !options.fugacity_tolerance.is_finite()
            || options.fugacity_tolerance <= 0.0
        {
            return Err(Error::InvalidInput(
                "positive flash iteration limits and tolerance required",
            ));
        }
        let (stability, feed) = self.stability_with_options(t, p, z, options)?;
        if stability.stable {
            let phase = stability.feed_phase;
            return Ok(FlashResult {
                temperature: t,
                pressure: p,
                feed: z.to_vec(),
                phase,
                liquid: if phase == FlashPhase::Vapor {
                    vec![0.0; z.len()]
                } else {
                    z.to_vec()
                },
                vapor: if phase == FlashPhase::Vapor {
                    z.to_vec()
                } else {
                    vec![0.0; z.len()]
                },
                vapor_fraction: match phase {
                    FlashPhase::Vapor => 1.0,
                    FlashPhase::Liquid => 0.0,
                    _ => -1.0,
                },
                liquid_fraction: match phase {
                    FlashPhase::Vapor => 0.0,
                    FlashPhase::Liquid => 1.0,
                    _ => -1.0,
                },
                liquid_volume: if phase == FlashPhase::Vapor {
                    None
                } else {
                    Some(feed.volume)
                },
                vapor_volume: if phase == FlashPhase::Vapor {
                    Some(feed.volume)
                } else {
                    None
                },
                iterations: 0,
                fugacity_residual: 0.0,
            });
        }
        let limit = -1000.0 * f64::EPSILON;
        let liquid_unstable = stability.liquid_trial.tpd < limit;
        let vapor_unstable = stability.vapor_trial.tpd < limit;
        let mut log_k: Vec<_> = (0..z.len())
            .map(|i| {
                let l = if liquid_unstable {
                    stability.liquid_trial.ln_fugacity[i]
                } else {
                    feed.phi[i]
                };
                let v = if vapor_unstable {
                    stability.vapor_trial.ln_fugacity[i]
                } else {
                    feed.phi[i]
                };
                l - v
            })
            .collect();
        let mut allocation = vec![0.0; z.len()];
        // Successive substitution supplies a feasible split before Newton.
        for iteration in 0..20 {
            let k: Vec<_> = log_k.iter().map(|k| k.clamp(-700.0, 700.0).exp()).collect();
            let beta = numerics::rachford_rice(z, &k)?;
            let mut x: Vec<_> = z
                .iter()
                .zip(&k)
                .map(|(z, k)| z / ((1.0 - beta) + beta * k))
                .collect();
            let mut y: Vec<_> = x.iter().zip(&k).map(|(x, k)| x * k).collect();
            normalize(&mut x)?;
            normalize(&mut y)?;
            let l = evaluate(self, t, p, &x, Phase::Liquid, false)?;
            let v = evaluate(self, t, p, &y, Phase::Vapor, false)?;
            let mut error = 0.0_f64;
            for i in 0..z.len() {
                let new = l.phi[i] - v.phi[i];
                if z[i] > 0.0 {
                    error = error.max((new - log_k[i]).abs());
                }
                log_k[i] = new;
                allocation[i] = beta * y[i];
            }
            if beta > 1e-12 && beta < 1.0 - 1e-12 && error < options.fugacity_tolerance {
                let g = (1.0 - beta) * gibbs(&x, &l.phi) + beta * gibbs(&y, &v.phi);
                if g <= gibbs(z, &feed.phi) + 1e-12 {
                    return Ok(two_phase(
                        t,
                        p,
                        z,
                        x,
                        y,
                        beta,
                        l.volume,
                        v.volume,
                        iteration + 1,
                        error,
                    ));
                }
            }
            if iteration == 19 && (beta <= 0.0 || beta >= 1.0) {
                for i in 0..z.len() {
                    allocation[i] = if vapor_unstable {
                        (0.001 * stability.vapor_trial.composition[i])
                            .min(0.9 * z[i])
                            .max(1e-8 * z[i])
                    } else {
                        z[i] - (0.001 * stability.liquid_trial.composition[i])
                            .min(0.9 * z[i])
                            .max(1e-8 * z[i])
                    };
                }
            }
        }
        self.flash_newton(t, p, z, allocation, gibbs(z, &feed.phi), options)
    }

    fn flash_newton(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        mut amounts: Vec<f64>,
        g_feed: f64,
        options: FlashOptions,
    ) -> Result<FlashResult> {
        let n = z.len();
        let active: Vec<_> = (0..n).filter(|&i| z[i] > 0.0).collect();
        let m = active.len();
        let split_evaluate = |amounts: &[f64], jacobian: bool| -> Result<SplitEvaluation> {
            let mut x: Vec<_> = z.iter().zip(amounts).map(|(z, v)| z - v).collect();
            let mut y = amounts.to_vec();
            let bl = normalize(&mut x)?;
            let bv = normalize(&mut y)?;
            let l = evaluate(self, t, p, &x, Phase::Liquid, jacobian)?;
            let v = evaluate(self, t, p, &y, Phase::Vapor, jacobian)?;
            let g = bl * gibbs(&x, &l.phi) + bv * gibbs(&y, &v.phi);
            Ok(SplitEvaluation {
                gibbs: g,
                vapor_fraction: bv,
                liquid: x,
                vapor: y,
                liquid_properties: l,
                vapor_properties: v,
            })
        };
        for iteration in 0..options.max_iterations {
            let SplitEvaluation {
                gibbs: g,
                vapor_fraction: beta,
                liquid: x,
                vapor: y,
                liquid_properties: l,
                vapor_properties: v,
            } = split_evaluate(&amounts, true)?;
            let gradient: Vec<_> = active
                .iter()
                .map(|&i| y[i].ln() + v.phi[i] - x[i].ln() - l.phi[i])
                .collect();
            let error = gradient.iter().map(|g| g.abs()).fold(0.0, f64::max);
            if error < options.fugacity_tolerance && g > g_feed + 1e-12 {
                return Err(Error::NotConverged(
                    "TP split has higher Gibbs energy than the feed",
                ));
            }
            if error < options.fugacity_tolerance && g <= g_feed + 1e-12 {
                return Ok(two_phase(
                    t,
                    p,
                    z,
                    x,
                    y,
                    beta,
                    l.volume,
                    v.volume,
                    iteration + 21,
                    error,
                ));
            }
            let mut hessian = vec![0.0; m * m];
            for (ii, &i) in active.iter().enumerate() {
                for (jj, &j) in active.iter().enumerate() {
                    hessian[ii * m + jj] =
                        (v.dn[i * n + j] - 1.0) / beta + (l.dn[i * n + j] - 1.0) / (1.0 - beta);
                }
                hessian[ii * m + ii] += 1.0 / amounts[i] + 1.0 / (z[i] - amounts[i]);
            }
            let step = numerics::modified_newton_direction(&hessian, &gradient);
            let slope: f64 = gradient.iter().zip(&step).map(|(g, s)| g * s).sum();
            let mut length = 1.0_f64;
            for (ii, &i) in active.iter().enumerate() {
                if step[ii] < 0.0 {
                    length = length.min(-0.99 * amounts[i] / step[ii]);
                } else if step[ii] > 0.0 {
                    length = length.min(0.99 * (z[i] - amounts[i]) / step[ii]);
                }
            }
            let mut accepted = false;
            for _ in 0..40 {
                let mut candidate = amounts.clone();
                for (ii, &i) in active.iter().enumerate() {
                    candidate[i] += length * step[ii];
                }
                if let Ok(candidate_split) = split_evaluate(&candidate, false)
                    && candidate_split.gibbs <= g + 1e-4 * length * slope + 4e-15 * g.abs().max(1.0)
                {
                    amounts = candidate;
                    accepted = true;
                    break;
                }
                length *= 0.5;
            }
            if !accepted {
                return Err(Error::NotConverged("TP flash Gibbs line search"));
            }
        }
        Err(Error::NotConverged("TP flash Gibbs minimization"))
    }
}

#[allow(clippy::too_many_arguments)]
fn two_phase(
    t: f64,
    p: f64,
    z: &[f64],
    x: Vec<f64>,
    y: Vec<f64>,
    beta: f64,
    vl: f64,
    vv: f64,
    iterations: usize,
    error: f64,
) -> FlashResult {
    FlashResult {
        temperature: t,
        pressure: p,
        feed: z.to_vec(),
        phase: FlashPhase::TwoPhase,
        liquid: x,
        vapor: y,
        vapor_fraction: beta,
        liquid_fraction: 1.0 - beta,
        liquid_volume: Some(vl),
        vapor_volume: Some(vv),
        iterations,
        fugacity_residual: error,
    }
}
