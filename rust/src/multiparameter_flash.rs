//! Pure-fluid TP and caloric flashes for multiparameter Helmholtz models.
//!
//! TP phase selection compares the two source phase-seeded density solutions
//! when available. This deliberately does not claim a global search over all
//! stable density roots or solid phases. Caloric flashes use these homogeneous
//! branches and the model's fluid coexistence curve.
use super::*;
use crate::caloric::{CaloricFlashOptions, EquilibriumProperties};
use crate::{Error, FlashPhase, FlashResult, Phase};

pub(super) fn homogeneous(t: f64, p: f64, v: f64, phase: FlashPhase) -> FlashResult {
    let liquid = phase != FlashPhase::Vapor;
    FlashResult {
        temperature: t,
        pressure: p,
        feed: vec![1.],
        phase,
        liquid: vec![if liquid { 1. } else { 0. }],
        vapor: vec![if liquid { 0. } else { 1. }],
        vapor_fraction: if phase == FlashPhase::Single {
            -1.
        } else if liquid {
            0.
        } else {
            1.
        },
        liquid_fraction: if phase == FlashPhase::Single {
            -1.
        } else if liquid {
            1.
        } else {
            0.
        },
        liquid_volume: if liquid { Some(v) } else { None },
        vapor_volume: if liquid { None } else { Some(v) },
        iterations: 0,
        fugacity_residual: 0.,
    }
}

pub(super) fn split(sat: crate::Saturation, beta: f64) -> FlashResult {
    FlashResult {
        temperature: sat.temperature,
        pressure: sat.pressure,
        feed: vec![1.],
        phase: FlashPhase::TwoPhase,
        liquid: vec![1.],
        vapor: vec![1.],
        vapor_fraction: beta,
        liquid_fraction: 1. - beta,
        liquid_volume: Some(sat.liquid_volume),
        vapor_volume: Some(sat.vapor_volume),
        iterations: 0,
        fugacity_residual: 0.,
    }
}

pub(super) fn properties(
    prep: &PreparedMultiparameter<'_>,
    v: f64,
) -> Result<EquilibriumProperties> {
    // Evaluate each Helmholtz jet once, retaining finite caloric properties
    // even when the critical Cp denominator vanishes.
    let ideal = prep.ideal_jet(v)?;
    let res = prep.residual_jet(v)?;
    let r = prep.eos.gas_constant;
    let rt = r * prep.temperature;
    let p = rt / v * (1. + res.d_delta);
    let internal_energy = rt * (ideal.d_tau + res.d_tau);
    let enthalpy = internal_energy + p * v;
    let entropy = r * (ideal.d_tau + res.d_tau - ideal.alpha - res.alpha);
    let cv = -r * (ideal.d_tautau + res.d_tautau);
    let cp = cv
        + r * (1. + res.d_delta - res.d_deltatau).powi(2)
            / (1. + 2. * res.d_delta + res.d_deltadelta);
    let heat_capacity = (cp.is_finite() && cp > 0.).then_some(cp);
    if ![p, enthalpy, entropy, internal_energy]
        .iter()
        .all(|x| x.is_finite())
    {
        return Err(Error::InvalidInput(
            "non-finite multiparameter caloric state",
        ));
    }
    Ok(EquilibriumProperties {
        volume: v,
        enthalpy,
        entropy,
        internal_energy,
        heat_capacity,
    })
}

/// Local isobaric root continuation in log volume, initialized from the
/// matching coexistence endpoint. Every accepted step remains on a closed,
/// mechanically stable root; phase-seeded iteration is the fallback.
fn continue_volume(prep: &PreparedMultiparameter<'_>, p: f64, seed: f64) -> Result<f64> {
    positive(seed, "coexistence volume seed")?;
    let mut x = seed.ln();
    for _ in 0..60 {
        let v = x.exp();
        let state = prep.pressure(v)?;
        let pressure_scale = p.abs() + prep.eos.gas_constant * prep.temperature / v;
        let floor = 8. * f64::EPSILON * pressure_scale;
        let residual = state.value - p;
        if residual.abs() <= floor && state.dv < 0. {
            return Ok(v);
        }
        let derivative = state.dv * v;
        if !derivative.is_finite() || derivative >= 0. {
            return Err(Error::NotConverged("stable local isobaric volume root"));
        }
        let step = (-residual / derivative).clamp(-0.5, 0.5);
        let old_error = residual.abs();
        let mut fraction = 1.;
        let mut accepted = None;
        for _ in 0..24 {
            let trial_x = x + fraction * step;
            let trial_v = trial_x.exp();
            if let Ok(trial) = prep.pressure(trial_v)
                && trial.dv < 0.
                && (trial.value - p).abs() < old_error
            {
                accepted = Some(trial_x);
                break;
            }
            fraction *= 0.5;
        }
        if let Some(next) = accepted {
            x = next;
        } else if old_error <= 32. * f64::EPSILON * pressure_scale && state.dv < 0. {
            // Close to criticality, caloric properties can amplify small
            // density errors. Continue polishing beyond the usual density
            // tolerance, then accept only the best representable pressure
            // closure when no log-volume trial improves it.
            return Ok(v);
        } else {
            return Err(Error::NotConverged("backtracked isobaric volume root"));
        }
    }
    Err(Error::NotConverged("local isobaric volume root"))
}

#[derive(Clone, Copy)]
struct CaloricJacobian {
    pressure: f64,
    pressure_ln_t: f64,
    pressure_ln_rho: f64,
    caloric: f64,
    caloric_ln_t: f64,
    caloric_ln_rho: f64,
    mechanical: f64,
}

fn caloric_jacobian(
    prep: &PreparedMultiparameter<'_>,
    rho: f64,
    enthalpy: bool,
) -> Result<CaloricJacobian> {
    let v = 1. / rho;
    let ideal = prep.ideal_jet(v)?;
    let res = prep.residual_jet(v)?;
    let r = prep.eos.gas_constant;
    let t = prep.temperature;
    let mechanical = 1. + 2. * res.d_delta + res.d_deltadelta;
    let b = 1. + res.d_delta - res.d_deltatau;
    let pressure = rho * r * t * (1. + res.d_delta);
    let cv = -r * (ideal.d_tautau + res.d_tautau);
    let (caloric, caloric_ln_t, caloric_ln_rho) = if enthalpy {
        let u = r * t * (ideal.d_tau + res.d_tau);
        (
            u + pressure / rho,
            t * (cv + r * b),
            r * t * (res.d_deltatau + res.d_delta + res.d_deltadelta),
        )
    } else {
        (
            r * (ideal.d_tau + res.d_tau - ideal.alpha - res.alpha),
            cv,
            -r * b,
        )
    };
    Ok(CaloricJacobian {
        pressure,
        pressure_ln_t: rho * r * t * b,
        pressure_ln_rho: rho * r * t * mechanical,
        caloric,
        caloric_ln_t,
        caloric_ln_rho,
        mechanical,
    })
}

/// Coupled local pressure/caloric correction in log T and log density. This
/// resolves near-critical endpoints where the isobaric scalar slope becomes
/// too large for a temperature-only iteration to close the source tolerance.
struct CaloricSpecification {
    pressure: f64,
    value: f64,
    enthalpy: bool,
}

fn coupled_caloric(
    model: &Multiparameter,
    specification: CaloricSpecification,
    options: CaloricFlashOptions,
    branch: FlashPhase,
    mut t: f64,
    mut rho: f64,
    scalar_iterations: usize,
) -> Result<FlashResult> {
    let CaloricSpecification {
        pressure: p_spec,
        value: target,
        enthalpy,
    } = specification;
    let r = model.gas_constant();
    let caloric_scale = |t: f64| r * if enthalpy { t } else { 1. };
    let norm = |p: f64, q: f64, t: f64, rho: f64| {
        ((p - p_spec) / (p_spec.abs() + r * t * rho))
            .abs()
            .max((q - target).abs() / caloric_scale(t))
    };
    for iteration in 1..=40.min(options.max_iterations) {
        if !t.is_finite()
            || t < options.minimum_temperature
            || t > options.maximum_temperature
            || !rho.is_finite()
            || rho <= 0.
        {
            return Err(Error::NotConverged("coupled caloric trial outside bounds"));
        }
        let prep = model.prepare(t)?;
        let v = 1. / rho;
        let jac = caloric_jacobian(&prep, rho, enthalpy)?;
        let pressure = jac.pressure;
        let actual = jac.caloric;
        let p_error = (pressure - p_spec).abs() / (p_spec.abs() + r * t * rho);
        let caloric_error = (actual - target).abs() / caloric_scale(t);
        if p_error <= 32. * f64::EPSILON
            && caloric_error <= options.scaled_tolerance
            && jac.mechanical > 0.
        {
            let mut flash = homogeneous(t, p_spec, v, branch);
            flash.liquid = vec![1.];
            flash.vapor = vec![1.];
            flash.iterations = scalar_iterations + iteration;
            return Ok(flash);
        }
        if !(pressure.is_finite()
            && pressure > 0.
            && jac.mechanical.is_finite()
            && jac.mechanical > 0.)
        {
            return Err(Error::NotConverged(
                "coupled caloric state is not mechanically stable",
            ));
        }
        let a = jac.pressure_ln_t / p_spec;
        let bb = jac.pressure_ln_rho / p_spec;
        let c = jac.caloric_ln_t / caloric_scale(t);
        let dd = jac.caloric_ln_rho / caloric_scale(t);
        let det = a * dd - bb * c;
        let det_scale = (a * dd).abs() + (bb * c).abs();
        if !det.is_finite() || det.abs() <= 64. * f64::EPSILON * det_scale.max(f64::MIN_POSITIVE) {
            return Err(Error::NotConverged("singular coupled caloric Jacobian"));
        }
        let rp = pressure / p_spec - 1.;
        let rc = (actual - target) / caloric_scale(t);
        let dt = ((-rp * dd + bb * rc) / det).clamp(-0.25, 0.25);
        let drho = ((c * rp - a * rc) / det).clamp(-0.5, 0.5);
        let current_norm = norm(pressure, actual, t, rho);
        let mut fraction = 1.;
        let mut accepted = None;
        for _ in 0..28 {
            let trial_t = t * (fraction * dt).exp();
            let trial_rho = rho * (fraction * drho).exp();
            if trial_t >= options.minimum_temperature && trial_t <= options.maximum_temperature {
                let trial_prep = model.prepare(trial_t)?;
                if let Ok(trial_jac) = caloric_jacobian(&trial_prep, trial_rho, enthalpy)
                    && trial_jac.pressure > 0.
                    && trial_jac.mechanical > 0.
                    && norm(trial_jac.pressure, trial_jac.caloric, trial_t, trial_rho)
                        < current_norm
                {
                    accepted = Some((trial_t, trial_rho));
                    break;
                }
            }
            fraction *= 0.5;
        }
        (t, rho) = accepted.ok_or(Error::NotConverged("coupled caloric backtracking"))?;
    }
    Err(Error::NotConverged("coupled caloric flash"))
}

impl Multiparameter {
    /// Minimum dimensionless-Gibbs choice among liquid/vapor source density
    /// iterations that close on mechanically stable roots.
    /// A single available source phase-seeded root is returned as `Single`.
    pub fn tp_flash(&self, t: f64, p: f64) -> Result<FlashResult> {
        positive(t, "temperature")?;
        positive(p, "pressure")?;
        let prep = self.prepare(t)?;
        let mut candidates: [Option<(f64, FlashPhase, f64)>; 2] = [None, None];
        for (i, (phase, label)) in [
            (Phase::Liquid, FlashPhase::Liquid),
            (Phase::Vapor, FlashPhase::Vapor),
        ]
        .into_iter()
        .enumerate()
        {
            if let Ok(v) = prep.volume(p, phase)
                && let Ok(residual) = prep.residual_jet(v)
            {
                // The temperature-only ideal contribution cancels between
                // phases. Remove it before comparing nearly equal roots.
                let delta = 1. / (v * self.rc());
                let g = delta.ln() + residual.alpha + p * v / (self.gas_constant() * t);
                if g.is_finite() {
                    candidates[i] = Some((v, label, g));
                }
            }
        }
        let mut roots = [None, None];
        let mut n = 0;
        for c in candidates.into_iter().flatten() {
            if !roots[..n]
                .iter()
                .flatten()
                .any(|(v, _, _): &(f64, FlashPhase, f64)| {
                    ((v - c.0).abs() / v.abs().max(c.0.abs())) < 1e-8
                })
            {
                roots[n] = Some(c);
                n += 1;
            }
        }
        if n == 0 {
            return Err(Error::NoPhysicalRoot);
        }
        // At coexistence T/P do not specify a quality. Floating-point Gibbs
        // differences can choose a different endpoint from the source's
        // general stability trajectory; use saturated_flash for fixed quality.
        let selected = roots[..n]
            .iter()
            .flatten()
            .copied()
            .min_by(|a, b| a.2.total_cmp(&b.2))
            .ok_or(Error::NoPhysicalRoot)?;
        Ok(homogeneous(
            t,
            p,
            selected.0,
            if n == 1 {
                FlashPhase::Single
            } else {
                selected.1
            },
        ))
    }

    /// Fluid coexistence state at an explicitly supplied vapor fraction.
    pub fn saturated_flash(&self, t: f64, vapor_fraction: f64) -> Result<FlashResult> {
        if !vapor_fraction.is_finite() || !(0. ..=1.).contains(&vapor_fraction) {
            return Err(Error::InvalidInput("physical vapor fraction required"));
        }
        Ok(split(self.saturation_pressure(t)?, vapor_fraction))
    }

    /// Recompute phase states from this model; cached flash volumes are ignored.
    pub fn equilibrium_properties(&self, flash: &FlashResult) -> Result<EquilibriumProperties> {
        if flash.feed != [1.] {
            return Err(Error::InvalidInput(
                "multiparameter flash requires pure normalized feed",
            ));
        }
        positive(flash.temperature, "temperature")?;
        positive(flash.pressure, "pressure")?;
        let prep = self.prepare(flash.temperature)?;
        if flash.phase != FlashPhase::TwoPhase {
            let v = match flash.phase {
                FlashPhase::Liquid => prep.volume(flash.pressure, Phase::Liquid)?,
                FlashPhase::Vapor => prep.volume(flash.pressure, Phase::Vapor)?,
                // SINGLEPH properties use the source liquid iterator. The
                // vapor fallback supports a sole vapor candidate if that
                // liquid-seeded iteration failed during the TP calculation.
                _ => prep
                    .volume(flash.pressure, Phase::Liquid)
                    .or_else(|_| prep.volume(flash.pressure, Phase::Vapor))?,
            };
            return properties(&prep, v);
        }
        let beta = flash.vapor_fraction;
        if !beta.is_finite()
            || !(0. ..=1.).contains(&beta)
            || !flash.liquid_fraction.is_finite()
            || (flash.liquid_fraction + beta - 1.).abs() > 1e-12
            || flash.liquid != [1.]
            || flash.vapor != [1.]
        {
            return Err(Error::InvalidInput(
                "physical pure two-phase amounts required",
            ));
        }
        let sat = self.saturation_pressure(flash.temperature)?;
        if (flash.pressure - sat.pressure).abs() > 1e-8 * flash.pressure {
            return Err(Error::InvalidInput(
                "pure split is off the coexistence curve",
            ));
        }
        let l = properties(&prep, sat.liquid_volume)?;
        let v = properties(&prep, sat.vapor_volume)?;
        let mix = |a, b| (1. - beta) * a + beta * b;
        Ok(EquilibriumProperties {
            volume: mix(l.volume, v.volume),
            enthalpy: mix(l.enthalpy, v.enthalpy),
            entropy: mix(l.entropy, v.entropy),
            internal_energy: mix(l.internal_energy, v.internal_energy),
            heat_capacity: None,
        })
    }

    /// Solve a pure-fluid pressure/enthalpy specification. Homogeneous states
    /// use the source 80..999 K bounds; reduced LJTS uses 0.3..3 times Tc.
    pub fn ph_flash(&self, p: f64, h: f64) -> Result<FlashResult> {
        self.ph_flash_with_options(p, h, self.default_caloric_options())
    }
    /// Solve a pure-fluid pressure/entropy specification, with the same
    /// model-dependent default bounds as [`Self::ph_flash`].
    pub fn ps_flash(&self, p: f64, s: f64) -> Result<FlashResult> {
        self.ps_flash_with_options(p, s, self.default_caloric_options())
    }
    fn default_caloric_options(&self) -> CaloricFlashOptions {
        if self.gas_constant() == 1. {
            CaloricFlashOptions {
                minimum_temperature: 0.3 * self.tc(),
                maximum_temperature: 3. * self.tc(),
                ..Default::default()
            }
        } else {
            CaloricFlashOptions::default()
        }
    }
    /// Solve pressure/enthalpy using explicit temperature bounds and guess.
    pub fn ph_flash_with_options(
        &self,
        p: f64,
        h: f64,
        o: CaloricFlashOptions,
    ) -> Result<FlashResult> {
        self.caloric_flash(p, h, true, o)
    }
    /// Solve pressure/entropy using explicit temperature bounds and guess.
    pub fn ps_flash_with_options(
        &self,
        p: f64,
        s: f64,
        o: CaloricFlashOptions,
    ) -> Result<FlashResult> {
        self.caloric_flash(p, s, false, o)
    }

    fn caloric_flash(
        &self,
        p: f64,
        target: f64,
        enthalpy: bool,
        o: CaloricFlashOptions,
    ) -> Result<FlashResult> {
        positive(p, "pressure")?;
        positive(o.minimum_temperature, "minimum temperature")?;
        positive(o.maximum_temperature, "maximum temperature")?;
        positive(o.scaled_tolerance, "caloric residual tolerance")?;
        if !target.is_finite()
            || o.minimum_temperature >= o.maximum_temperature
            || o.max_iterations == 0
            || o.temperature_guess
                .is_some_and(|t| !t.is_finite() || t <= 0.)
        {
            return Err(Error::InvalidInput(
                "finite caloric specification and valid temperature bounds required",
            ));
        }
        let value = |q: &EquilibriumProperties| if enthalpy { q.enthalpy } else { q.entropy };
        let critical_pressure = if self.pc() > 0. {
            self.pc()
        } else {
            self.prepare(self.tc())?.pressure(1. / self.rc())?.value
        };
        let mut branch = FlashPhase::Single;
        let mut saturation_anchor = None;
        if p < critical_pressure
            && let Ok(sat) = self.saturation_temperature(p)
        {
            saturation_anchor = Some(sat);
            let prep = self.prepare(sat.temperature)?;
            let l = properties(&prep, sat.liquid_volume)?;
            let v = properties(&prep, sat.vapor_volume)?;
            let (a, b) = (value(&l), value(&v));
            let tolerance = o.scaled_tolerance
                * self.gas_constant()
                * if enthalpy { sat.temperature } else { 1. };
            let beta = if b > a && target >= a && target <= b {
                Some((target - a) / (b - a))
            } else if (target - a).abs() <= tolerance {
                Some(0.)
            } else if (target - b).abs() <= tolerance {
                Some(1.)
            } else {
                None
            };
            if let Some(beta) = beta {
                let mixed = (1. - beta) * a + beta * b;
                if (mixed - target).abs() > tolerance {
                    return Err(Error::NotConverged("latent caloric endpoint closure"));
                }
                let mut result = split(sat, beta);
                result.pressure = p;
                return Ok(result);
            }
            if target < a {
                branch = FlashPhase::Liquid;
            } else {
                branch = FlashPhase::Vapor;
            }
        }
        if o.minimum_temperature >= o.maximum_temperature {
            return Err(Error::InvalidInput(
                "caloric specification outside temperature bounds",
            ));
        }
        let calc = |t: f64| -> Result<(f64, EquilibriumProperties)> {
            let prep = self.prepare(t)?;
            // The source pure caloric solver uses the liquid density iterator
            // for its SINGLEPH path, including supercritical specifications.
            let ph = match branch {
                FlashPhase::Vapor => Phase::Vapor,
                _ => Phase::Liquid,
            };
            let vol = if let Some(sat) = saturation_anchor {
                let seed = if branch == FlashPhase::Vapor {
                    sat.vapor_volume
                } else {
                    sat.liquid_volume
                };
                if (t - sat.temperature).abs() <= 8. * f64::EPSILON * sat.temperature {
                    seed
                } else {
                    continue_volume(&prep, p, seed).or_else(|_| prep.volume(p, ph))?
                }
            } else {
                prep.volume(p, ph)?
            };
            let q = properties(&prep, vol)?;
            Ok((vol, q))
        };
        let scale = |t: f64| self.gas_constant() * if enthalpy { t } else { 1. };
        let mut t = o
            .temperature_guess
            .unwrap_or(self.tc())
            .clamp(o.minimum_temperature, o.maximum_temperature);
        let fallback = |temperature: f64, volume: f64, iterations: usize| {
            coupled_caloric(
                self,
                CaloricSpecification {
                    pressure: p,
                    value: target,
                    enthalpy,
                },
                o,
                branch,
                temperature,
                1. / volume,
                iterations,
            )
        };
        let (mut volume, mut q) = match calc(t) {
            Ok(state) => state,
            Err(e) => {
                let seed = if let Some(sat) = saturation_anchor {
                    if branch == FlashPhase::Vapor {
                        sat.vapor_volume
                    } else {
                        sat.liquid_volume
                    }
                } else {
                    let prep = self.prepare(t)?;
                    prep.volume(
                        p,
                        if branch == FlashPhase::Vapor {
                            Phase::Vapor
                        } else {
                            Phase::Liquid
                        },
                    )?
                };
                if !matches!(&e, Error::NoPhysicalRoot | Error::NotConverged(_)) {
                    return Err(e);
                }
                return fallback(t, seed, 0);
            }
        };
        let mut residual = value(&q) - target;
        let mut low = t;
        let mut high = t;
        let mut bracketed = false;
        let mut step = 0.1 * t;
        let mut expansion_iterations = 0;
        for iteration_count in 0..32 {
            expansion_iterations = iteration_count + 1;
            if residual.abs() <= o.scaled_tolerance * scale(t) {
                let mut flash = homogeneous(t, p, volume, branch);
                flash.liquid = vec![1.];
                flash.vapor = vec![1.];
                flash.iterations = expansion_iterations;
                return Ok(flash);
            }
            let direction = if residual > 0. { -1. } else { 1. };
            let mut bounded = if direction < 0. {
                step.min(0.3 * t)
            } else {
                step
            };
            if let Some(cp) = q.heat_capacity {
                let slope = if enthalpy { cp } else { cp / t };
                if slope > 0. {
                    bounded = bounded.min(1.25 * (residual / slope).abs());
                }
            }
            let next =
                (t + direction * bounded).clamp(o.minimum_temperature, o.maximum_temperature);
            if next == t {
                return fallback(t, volume, expansion_iterations);
            }
            let mut trial = next;
            let mut evaluated = None;
            let mut last_error = None;
            for _ in 0..16 {
                match calc(trial) {
                    Ok(result) => {
                        evaluated = Some(result);
                        break;
                    }
                    Err(e @ (Error::NotConverged(_) | Error::NoPhysicalRoot)) => {
                        last_error = Some(e)
                    }
                    Err(e) => return Err(e),
                }
                trial = 0.5 * (trial + t);
            }
            let (nv, nq) = match evaluated {
                Some(q) => q,
                None => {
                    let _ = last_error;
                    return fallback(t, volume, expansion_iterations);
                }
            };
            let next = trial;
            let nr = value(&nq) - target;
            if nr.signum() != residual.signum() {
                low = t.min(next);
                high = t.max(next);
                volume = nv;
                q = nq;
                residual = nr;
                t = next;
                bracketed = true;
                break;
            }
            step = (next - t).abs() * 2.;
            t = next;
            volume = nv;
            q = nq;
            residual = nr;
        }
        if !bracketed {
            return fallback(t, volume, expansion_iterations);
        }
        for iteration_count in 1..=o.max_iterations {
            if residual.abs() <= o.scaled_tolerance * scale(t) {
                let mut flash = homogeneous(t, p, volume, branch);
                flash.liquid = vec![1.];
                flash.vapor = vec![1.];
                flash.iterations = expansion_iterations + iteration_count;
                return Ok(flash);
            }
            if residual < 0. {
                low = low.max(t);
            } else {
                high = high.min(t);
            }
            if high - low <= 8. * f64::EPSILON * t {
                return fallback(t, volume, expansion_iterations + iteration_count);
            }
            let slope = q
                .heat_capacity
                .map(|c| if enthalpy { c } else { c / t })
                .filter(|s| s.is_finite() && *s > 0.);
            let newton = slope.map(|s| t - residual / s).unwrap_or(f64::NAN);
            let next = if newton > low && newton < high {
                newton
            } else {
                (low + high) * 0.5
            };
            // If the local EOS branch cannot be evaluated, approach from the opposite bracket endpoint.
            let (nv, nq, nt) = match calc(next) {
                Ok((v, q)) => (v, q, next),
                Err(Error::NotConverged(_) | Error::NoPhysicalRoot) => {
                    let fallback_t = if residual < 0. { high } else { low };
                    match calc(fallback_t) {
                        Ok((v, q)) => (v, q, fallback_t),
                        Err(_) => {
                            return fallback(t, volume, expansion_iterations + iteration_count);
                        }
                    }
                }
                Err(e) => return Err(e),
            };
            volume = nv;
            q = nq;
            t = nt;
            residual = value(&q) - target;
        }
        fallback(t, volume, expansion_iterations + o.max_iterations)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64) {
        let scale = actual.abs().max(expected.abs()).max(1e-8);
        assert!(
            (actual - expected).abs() <= 3e-5 * scale + 2e-8,
            "finite difference {actual:e}, analytic {expected:e}"
        );
    }

    #[test]
    fn coupled_caloric_jacobian_matches_log_coordinate_differences() {
        let mut models = vec![
            Multiparameter::new("C3").unwrap(),
            Multiparameter::new("R134A").unwrap(),
            Multiparameter::reduced_lennard_jones_truncated_shifted(),
        ];
        let step: f64 = 2e-5;
        for model in models.drain(..) {
            let t = 1.2 * model.tc();
            let rho = 0.2 * model.rc();
            let prep = model.prepare(t).unwrap();
            for enthalpy in [true, false] {
                let j = caloric_jacobian(&prep, rho, enthalpy).unwrap();
                assert!(j.mechanical > 0., "test point must be mechanically stable");
                let plus_t = model.prepare(t * step.exp()).unwrap();
                let minus_t = model.prepare(t * (-step).exp()).unwrap();
                let jt_plus = caloric_jacobian(&plus_t, rho, enthalpy).unwrap();
                let jt_minus = caloric_jacobian(&minus_t, rho, enthalpy).unwrap();
                close(
                    (jt_plus.pressure - jt_minus.pressure) / (2. * step),
                    j.pressure_ln_t,
                );
                close(
                    (jt_plus.caloric - jt_minus.caloric) / (2. * step),
                    j.caloric_ln_t,
                );

                let jrho_plus = caloric_jacobian(&prep, rho * step.exp(), enthalpy).unwrap();
                let jrho_minus = caloric_jacobian(&prep, rho * (-step).exp(), enthalpy).unwrap();
                close(
                    (jrho_plus.pressure - jrho_minus.pressure) / (2. * step),
                    j.pressure_ln_rho,
                );
                close(
                    (jrho_plus.caloric - jrho_minus.caloric) / (2. * step),
                    j.caloric_ln_rho,
                );
            }
        }
    }
}
