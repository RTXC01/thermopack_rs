//! Coupled UV/SV equations from uv_solver/sv_solver/state_functions. Phase
//! amounts are updated through the smaller amount to preserve trace phases.
use super::*;
use crate::TvProperties;

#[derive(Debug, Clone, Copy)]
/// Bounds and stopping criteria for a molar UV/SV equilibrium flash.
/// Use `Cubic::volume_flash_options` to inherit model temperature bounds.
pub struct VolumeFlashOptions {
    /// Positive starting temperature (K); clamped to the search bounds.
    pub temperature_guess: Option<f64>,
    /// Positive starting pressure (Pa); clamped to the search bounds.
    pub pressure_guess: Option<f64>,
    /// Lower search temperature (K).
    pub minimum_temperature: f64,
    /// Upper search temperature (K).
    pub maximum_temperature: f64,
    /// Lower search pressure (Pa).
    pub minimum_pressure: f64,
    /// Upper search pressure (Pa).
    pub maximum_pressure: f64,
    /// Maximum coupled iterations for each initialization attempt.
    pub max_iterations: usize,
    /// Dimensionless chemical-potential, caloric and volume tolerance.
    pub scaled_tolerance: f64,
    /// Initialization and final equilibrium branch verification options.
    pub tp: FlashOptions,
}
impl Default for VolumeFlashOptions {
    fn default() -> Self {
        Self {
            temperature_guess: None,
            pressure_guess: None,
            minimum_temperature: 80.,
            maximum_temperature: 999.,
            minimum_pressure: 10.,
            maximum_pressure: 1e8,
            max_iterations: 100,
            scaled_tolerance: 1e-12,
            tp: FlashOptions::default(),
        }
    }
}

#[derive(Clone, Copy)]
enum Spec {
    Energy,
    Entropy,
}
impl Spec {
    fn value(self, p: &TvProperties) -> f64 {
        match self {
            Self::Energy => p.internal_energy,
            Self::Entropy => p.entropy,
        }
    }
    fn scale(self, t: f64) -> f64 {
        match self {
            Self::Energy => R * t,
            Self::Entropy => R,
        }
    }
    fn closure(self, p: &EquilibriumProperties) -> f64 {
        match self {
            Self::Energy => p.internal_energy,
            Self::Entropy => p.entropy,
        }
    }
}

struct Evaluation {
    residual: Vec<f64>,
    matrix: Vec<f64>,
    flash: FlashResult,
    closure: f64,
}
struct Amounts {
    vapor: Vec<f64>,
    liquid: Vec<f64>,
}

impl Cubic {
    pub fn volume_flash_options(&self) -> VolumeFlashOptions {
        let (minimum_temperature, maximum_temperature) = self.temperature_bounds();
        VolumeFlashOptions {
            minimum_temperature,
            maximum_temperature,
            ..Default::default()
        }
    }
    /// Recover equilibrium from molar U (J/mol), V (m3/mol) and normalized z.
    pub fn uv_flash(&self, z: &[f64], energy: f64, volume: f64) -> Result<FlashResult> {
        self.uv_flash_with_options(z, energy, volume, self.volume_flash_options())
    }
    pub fn uv_flash_with_options(
        &self,
        z: &[f64],
        energy: f64,
        volume: f64,
        options: VolumeFlashOptions,
    ) -> Result<FlashResult> {
        self.volume_flash(z, energy, volume, Spec::Energy, options)
    }
    /// Recover equilibrium from molar S (J/(mol K)), V (m3/mol), normalized z.
    pub fn sv_flash(&self, z: &[f64], entropy: f64, volume: f64) -> Result<FlashResult> {
        self.sv_flash_with_options(z, entropy, volume, self.volume_flash_options())
    }
    pub fn sv_flash_with_options(
        &self,
        z: &[f64],
        entropy: f64,
        volume: f64,
        options: VolumeFlashOptions,
    ) -> Result<FlashResult> {
        self.volume_flash(z, entropy, volume, Spec::Entropy, options)
    }

    fn volume_flash(
        &self,
        z: &[f64],
        target: f64,
        volume: f64,
        spec: Spec,
        options: VolumeFlashOptions,
    ) -> Result<FlashResult> {
        positive(volume, "molar volume")?;
        positive(options.minimum_temperature, "minimum temperature")?;
        positive(options.maximum_temperature, "maximum temperature")?;
        positive(options.minimum_pressure, "minimum pressure")?;
        positive(options.maximum_pressure, "maximum pressure")?;
        positive(options.scaled_tolerance, "volume flash tolerance")?;
        if !target.is_finite()
            || options.minimum_temperature >= options.maximum_temperature
            || options.minimum_pressure >= options.maximum_pressure
            || options.max_iterations == 0
            || options
                .temperature_guess
                .is_some_and(|t| !t.is_finite() || t <= 0.)
            || options
                .pressure_guess
                .is_some_and(|p| !p.is_finite() || p <= 0.)
            || options.tp.max_iterations == 0
            || options.tp.stability_iterations == 0
            || !options.tp.fugacity_tolerance.is_finite()
            || options.tp.fugacity_tolerance <= 0.
        {
            return Err(Error::InvalidInput(
                "finite volume specification and valid flash bounds/guesses required",
            ));
        }
        let critical = self.pseudo_critical(z)?;
        let initial = options
            .temperature_guess
            .unwrap_or(critical.temperature)
            .clamp(options.minimum_temperature, options.maximum_temperature);
        let mut homogeneous = initial;
        for _ in 0..24 {
            let Ok(state) = self
                .prepare(homogeneous, z)
                .and_then(|s| s.caloric_tv(volume))
            else {
                break;
            };
            let residual = target - spec.value(&state);
            if residual.abs() < options.scaled_tolerance * spec.scale(homogeneous) {
                break;
            }
            let derivative = match spec {
                Spec::Energy => state.heat_capacity_v,
                Spec::Entropy => state.heat_capacity_v / homogeneous,
            };
            if !derivative.is_finite() || derivative == 0. {
                break;
            }
            let change = (residual / derivative).clamp(-0.3 * homogeneous, 0.3 * homogeneous);
            let next = (homogeneous + change)
                .clamp(options.minimum_temperature, options.maximum_temperature);
            if next == homogeneous {
                break;
            }
            homogeneous = next;
        }
        let homop = self
            .prepare(homogeneous, z)
            .and_then(|s| s.pressure(volume))
            .map(|p| p.value)
            .unwrap_or(0.);
        let mut seeds = vec![(
            initial,
            Some(
                options
                    .pressure_guess
                    .unwrap_or(R * initial / volume)
                    .clamp(options.minimum_pressure, options.maximum_pressure),
            ),
        )];
        if homop >= options.minimum_pressure && homop <= options.maximum_pressure {
            seeds.insert(0, (homogeneous, Some(homop)));
        }
        // Coexistence seeds include targets whose hypothetical homogeneous
        // state lies inside a negative-pressure loop, including pure latent.
        for t in [
            homogeneous,
            initial,
            0.8 * critical.temperature,
            0.6 * critical.temperature,
            0.95 * critical.temperature,
            0.4 * critical.temperature,
        ] {
            let t = t.clamp(options.minimum_temperature, options.maximum_temperature);
            seeds.push((t, None));
        }
        seeds.push((
            initial,
            Some(
                (0.1 * critical.pressure).clamp(options.minimum_pressure, options.maximum_pressure),
            ),
        ));
        for (t, p) in seeds {
            let p = if let Some(p) = p {
                p
            } else if let Ok(sat) = self.pseudo_saturation_pressure(t, z)
                && sat.pressure >= options.minimum_pressure
                && sat.pressure <= options.maximum_pressure
            {
                sat.pressure
            } else {
                continue;
            };
            if let Ok(seed) = self.tp_flash_with_options(t, p, z, options.tp)
                && let Ok(candidate) = coupled(self, z, target, volume, spec, seed, options)
                && let Ok(result) = accept(self, candidate, target, volume, spec, options)
            {
                return Ok(result);
            }
            // Pure TP cannot determine quality. A finite amount initializer
            // also lets the coupled mixture equations find a separated state.
            let pure = z.iter().filter(|z| **z > 0.).count() == 1;
            let (x, y, bv) = if pure {
                (z.to_vec(), z.to_vec(), 0.5)
            } else {
                let k = self.wilson_k(t, p)?;
                let mut vapor: Vec<_> = z
                    .iter()
                    .zip(&k)
                    .map(|(z, k)| 0.5 * z * k / (0.5 + 0.5 * k))
                    .collect();
                for (v, &z) in vapor.iter_mut().zip(z) {
                    *v = v.clamp(1e-6 * z, (1. - 1e-6) * z);
                }
                let bv: f64 = vapor.iter().sum();
                let liquid: Vec<_> = z
                    .iter()
                    .zip(&vapor)
                    .map(|(z, v)| (z - v) / (1. - bv))
                    .collect();
                (liquid, vapor.iter().map(|v| v / bv).collect(), bv)
            };
            let seed = FlashResult {
                temperature: t,
                pressure: p,
                feed: z.to_vec(),
                phase: FlashPhase::TwoPhase,
                liquid: x,
                vapor: y,
                vapor_fraction: bv,
                liquid_fraction: 1. - bv,
                liquid_volume: None,
                vapor_volume: None,
                iterations: 0,
                fugacity_residual: 0.,
            };
            if let Ok(candidate) = coupled(self, z, target, volume, spec, seed, options)
                && let Ok(result) = accept(self, candidate, target, volume, spec, options)
            {
                return Ok(result);
            }
        }
        Err(Error::NotConverged(
            "UV/SV coupled phase and caloric equations",
        ))
    }
}

#[allow(clippy::too_many_arguments)]
fn evaluate(
    model: &Cubic,
    t: f64,
    p: f64,
    z: &[f64],
    amounts: Option<&Amounts>,
    root: Phase,
    target: f64,
    volume: f64,
    spec: Spec,
) -> Result<Evaluation> {
    let n = z.len();
    let dim = if amounts.is_some() { n + 2 } else { 2 };
    let mut matrix = vec![0.; dim * dim];
    let mut residual = vec![0.; dim];
    let mut flash = FlashResult {
        temperature: t,
        pressure: p,
        feed: z.to_vec(),
        phase: FlashPhase::Single,
        liquid: z.to_vec(),
        vapor: vec![0.; n],
        vapor_fraction: -1.,
        liquid_fraction: -1.,
        liquid_volume: None,
        vapor_volume: None,
        iterations: 0,
        fugacity_residual: 0.,
    };
    let (h, s, v, cp, vt, vp) = if let Some(amounts) = amounts {
        let bv: f64 = amounts.vapor.iter().sum();
        let bl: f64 = amounts.liquid.iter().sum();
        positive(bv, "vapor amount")?;
        positive(bl, "liquid amount")?;
        let x: Vec<_> = amounts.liquid.iter().map(|n| n / bl).collect();
        let y: Vec<_> = amounts.vapor.iter().map(|n| n / bv).collect();
        let liquid = model.prepare(t, &x)?;
        let vapor = model.prepare(t, &y)?;
        let vl = liquid.volume(p, Phase::Liquid)?;
        let vv = vapor.volume(p, Phase::Vapor)?;
        let (hl, sl, cpl) = phase_caloric(&liquid, vl)?;
        let (hv, sv, cpv) = phase_caloric(&vapor, vv)?;
        let pl = liquid.pressure(vl)?;
        let pv = vapor.pressure(vv)?;
        let mut fl = vec![0.; n];
        let mut fv = vec![0.; n];
        liquid.ln_fugacity_tp_into(p, vl, &mut fl)?;
        vapor.ln_fugacity_tp_into(p, vv, &mut fv)?;
        let mut dtl = vec![0.; n];
        let mut dtv = vec![0.; n];
        let mut dpl = vec![0.; n];
        let mut dpv = vec![0.; n];
        let mut dnl = vec![0.; n * n];
        let mut dnv = vec![0.; n * n];
        liquid.fugacity_derivatives_into(vl, &mut dtl, &mut dpl, &mut dnl)?;
        vapor.fugacity_derivatives_into(vv, &mut dtv, &mut dpv, &mut dnv)?;
        for i in 0..n {
            if z[i] > 0. {
                residual[i] = y[i].max(f64::EPSILON.powi(5)).ln()
                    - x[i].max(f64::EPSILON.powi(5)).ln()
                    + fv[i]
                    - fl[i];
                matrix[i * dim + n] = t * (dtv[i] - dtl[i]);
                matrix[n * dim + i] = matrix[i * dim + n];
                matrix[i * dim + n + 1] = p * (dpv[i] - dpl[i]);
                matrix[(n + 1) * dim + i] = matrix[i * dim + n + 1];
                for j in 0..n {
                    if z[j] > 0. {
                        matrix[i * dim + j] =
                            (dnv[i * n + j] - 1.) / bv + (dnl[i * n + j] - 1.) / bl;
                    }
                }
                matrix[i * dim + i] += 1. / amounts.vapor[i] + 1. / amounts.liquid[i];
            } else {
                matrix[i * dim + i] = 1.;
            }
        }
        flash.phase = FlashPhase::TwoPhase;
        flash.liquid = x;
        flash.vapor = y;
        flash.vapor_fraction = bv;
        flash.liquid_fraction = bl;
        flash.liquid_volume = Some(vl);
        flash.vapor_volume = Some(vv);
        flash.fugacity_residual = residual[..n].iter().map(|f| f.abs()).fold(0., f64::max);
        (
            bl * hl + bv * hv,
            bl * sl + bv * sv,
            bl * vl + bv * vv,
            bl * cpl + bv * cpv,
            -bl * pl.dt / pl.dv - bv * pv.dt / pv.dv,
            bl / pl.dv + bv / pv.dv,
        )
    } else {
        let state = model.prepare(t, z)?;
        let v = state.volume(p, root)?;
        let (h, s, cp) = phase_caloric(&state, v)?;
        let pressure = state.pressure(v)?;
        flash.phase = if root == Phase::Vapor {
            FlashPhase::Vapor
        } else {
            FlashPhase::Liquid
        };
        if root == Phase::Vapor {
            flash.vapor = z.to_vec();
            flash.liquid.fill(0.);
            flash.vapor_fraction = 1.;
            flash.liquid_fraction = 0.;
            flash.vapor_volume = Some(v);
        } else {
            flash.vapor_fraction = 0.;
            flash.liquid_fraction = 1.;
            flash.liquid_volume = Some(v);
        }
        (h, s, v, cp, -pressure.dt / pressure.dv, 1. / pressure.dv)
    };
    let it = dim - 2;
    let ip = dim - 1;
    residual[it] = match spec {
        Spec::Energy => (target + p * volume - h) / (R * t),
        Spec::Entropy => (target - s) / R,
    };
    residual[ip] = p * (v - volume) / (R * t);
    matrix[it * dim + it] = -cp / R;
    matrix[it * dim + ip] = p * vt / R;
    matrix[ip * dim + it] = p * vt / R;
    matrix[ip * dim + ip] = p * p * vp / (R * t);
    let value = match spec {
        Spec::Energy => h - p * v,
        Spec::Entropy => s,
    };
    let closure = ((value - target) / spec.scale(t))
        .abs()
        .max(((v - volume) / volume).abs());
    Ok(Evaluation {
        residual,
        matrix,
        flash,
        closure,
    })
}

#[allow(clippy::too_many_arguments)]
fn coupled(
    model: &Cubic,
    z: &[f64],
    target: f64,
    volume: f64,
    spec: Spec,
    seed: FlashResult,
    options: VolumeFlashOptions,
) -> Result<FlashResult> {
    let mut t = seed.temperature;
    let mut p = seed.pressure;
    let root = if seed.phase == FlashPhase::Vapor {
        Phase::Vapor
    } else {
        Phase::Liquid
    };
    let mut amounts = if seed.phase == FlashPhase::TwoPhase {
        Some(Amounts {
            vapor: seed.vapor.iter().map(|y| y * seed.vapor_fraction).collect(),
            liquid: seed
                .liquid
                .iter()
                .map(|x| x * seed.liquid_fraction)
                .collect(),
        })
    } else {
        None
    };
    for iteration in 0..options.max_iterations {
        let mut current = evaluate(model, t, p, z, amounts.as_ref(), root, target, volume, spec)?;
        let error = current.residual.iter().map(|f| f.abs()).fold(0., f64::max);
        if error <= options.scaled_tolerance && current.closure <= options.scaled_tolerance * 10. {
            current.flash.iterations = iteration + 1;
            return Ok(current.flash);
        }
        let dim = current.residual.len();
        let mut step: Vec<_> = current.residual.iter().map(|f| -f).collect();
        let mut pivots = vec![0; dim];
        let mut work = vec![0.; dim];
        if !numerics::solve_symmetric_in_place(
            &mut current.matrix,
            &mut step,
            &mut pivots,
            &mut work,
        )? {
            return Err(Error::NotConverged("singular UV/SV state matrix"));
        }
        let mut scale = 1.0_f64;
        for (position, low, high, value) in [
            (
                dim - 2,
                options.minimum_temperature.ln(),
                options.maximum_temperature.ln(),
                t.ln(),
            ),
            (
                dim - 1,
                options.minimum_pressure.ln(),
                options.maximum_pressure.ln(),
                p.ln(),
            ),
        ] {
            let delta = step[position];
            if value + delta < low {
                scale = scale.min((low - value) / delta);
            }
            if value + delta > high {
                scale = scale.min((high - value) / delta);
            }
        }
        if let Some(old) = &amounts {
            for (i, &delta) in step[..z.len()].iter().enumerate() {
                if delta < 0. {
                    scale = scale.min(-old.vapor[i] / delta);
                }
                if delta > 0. {
                    scale = scale.min(old.liquid[i] / delta);
                }
            }
        }
        if scale < 1. {
            scale *= 1. - 1e-12;
        }
        let norm: f64 = current.residual.iter().map(|f| f * f).sum();
        let mut accepted = false;
        for _ in 0..24 {
            let next_t = (t.ln() + scale * step[dim - 2]).exp();
            let next_p = (p.ln() + scale * step[dim - 1]).exp();
            let next_amounts = amounts.as_ref().map(|old| {
                let mut new = Amounts {
                    vapor: old.vapor.clone(),
                    liquid: old.liquid.clone(),
                };
                for i in 0..z.len() {
                    if old.vapor[i] >= old.liquid[i] {
                        new.liquid[i] -= scale * step[i];
                        new.vapor[i] = z[i] - new.liquid[i];
                    } else {
                        new.vapor[i] += scale * step[i];
                        new.liquid[i] = z[i] - new.vapor[i];
                    }
                }
                new
            });
            if let Ok(candidate) = evaluate(
                model,
                next_t,
                next_p,
                z,
                next_amounts.as_ref(),
                root,
                target,
                volume,
                spec,
            ) && candidate.residual.iter().map(|f| f * f).sum::<f64>() < norm
            {
                t = next_t;
                p = next_p;
                amounts = next_amounts;
                accepted = true;
                break;
            }
            scale *= 0.5;
        }
        if !accepted {
            return Err(Error::NotConverged("UV/SV state line search"));
        }
    }
    Err(Error::NotConverged("UV/SV state iteration limit"))
}

fn accept(
    model: &Cubic,
    candidate: FlashResult,
    target: f64,
    volume: f64,
    spec: Spec,
    options: VolumeFlashOptions,
) -> Result<FlashResult> {
    let props = model.equilibrium_properties(&candidate)?;
    if (spec.closure(&props) - target).abs()
        > options.scaled_tolerance * spec.scale(candidate.temperature) * 10.
        || (props.volume - volume).abs() > options.scaled_tolerance * volume * 10.
    {
        return Err(Error::NotConverged(
            "UV/SV independent specification closure",
        ));
    }
    // The coupled chemical-potential equation fixes pure coexistence. A TP
    // flash at that T,P cannot determine quality and must not replace it.
    if candidate.phase == FlashPhase::TwoPhase
        && candidate.feed.iter().filter(|z| **z > 0.).count() == 1
    {
        let vl = candidate.liquid_volume.ok_or(Error::NoPhysicalRoot)?;
        let vv = candidate.vapor_volume.ok_or(Error::NoPhysicalRoot)?;
        if (vv - vl).abs() > 1e-8 * vv.abs().max(vl.abs()) {
            return Ok(candidate);
        }
    }
    let stable = model.tp_flash_with_options(
        candidate.temperature,
        candidate.pressure,
        &candidate.feed,
        options.tp,
    )?;
    let stable_props = model.equilibrium_properties(&stable)?;
    // TP and the coupled solve have independent stopping criteria. Retain
    // the tightly closed UV/SV result; TP only verifies its equilibrium
    // branch at the precision supported by its chemical-potential solve.
    let stability_tolerance =
        (options.scaled_tolerance * 10.).max(options.tp.fugacity_tolerance * 100.);
    if (spec.closure(&stable_props) - target).abs()
        > stability_tolerance * spec.scale(candidate.temperature)
        || (stable_props.volume - volume).abs() > stability_tolerance * volume
    {
        return Err(Error::NotConverged("UV/SV phase stability closure"));
    }
    if dimensionless_gibbs(model, &candidate)?
        > dimensionless_gibbs(model, &stable)? + stability_tolerance
    {
        return Err(Error::NotConverged("UV/SV higher Gibbs energy branch"));
    }
    if stable.phase != FlashPhase::TwoPhase {
        if (spec.closure(&stable_props) - target).abs()
            > options.scaled_tolerance * spec.scale(candidate.temperature) * 10.
            || (stable_props.volume - volume).abs() > options.scaled_tolerance * volume * 10.
        {
            return Err(Error::NotConverged("UV/SV homogeneous branch closure"));
        }
        // Keep native TP phase labels and homogeneous fraction sentinels.
        // There are no uncertain phase compositions to replace in this case.
        return Ok(FlashResult {
            iterations: candidate.iterations,
            ..stable
        });
    }
    if candidate.phase != FlashPhase::TwoPhase {
        return Err(Error::NotConverged("UV/SV missing stable phase split"));
    }
    Ok(candidate)
}

fn dimensionless_gibbs(model: &Cubic, flash: &FlashResult) -> Result<f64> {
    let phase_value = |composition: &[f64], root| -> Result<f64> {
        let state = model.prepare(flash.temperature, composition)?;
        let volume = state.volume(flash.pressure, root)?;
        let mut phi = vec![0.; composition.len()];
        state.ln_fugacity_tp_into(flash.pressure, volume, &mut phi)?;
        Ok(composition
            .iter()
            .zip(phi)
            .filter(|(z, _)| **z > 0.)
            .map(|(z, phi)| z * (z.ln() + phi))
            .sum())
    };
    if flash.phase == FlashPhase::TwoPhase {
        Ok(
            flash.liquid_fraction * phase_value(&flash.liquid, Phase::Liquid)?
                + flash.vapor_fraction * phase_value(&flash.vapor, Phase::Vapor)?,
        )
    } else {
        phase_value(
            &flash.feed,
            if flash.phase == FlashPhase::Vapor {
                Phase::Vapor
            } else {
                Phase::Liquid
            },
        )
    }
}
