//! Equilibrium caloric properties and pressure/enthalpy or pressure/entropy
//! flashes. Uses analytic equilibrium heat capacity with safeguarded Newton.
//!
//! ```
//! use thermopack_rs::{Cubic, CubicKind};
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2", "C3"])?;
//! let z = [0.5, 0.3, 0.2];
//! let initial = model.tp_flash(200.0, 1e6, &z)?;
//! let properties = model.equilibrium_properties(&initial)?;
//! let recovered = model.ph_flash(1e6, &z, properties.enthalpy)?;
//! assert!((recovered.temperature - 200.0).abs() < 1e-6);
//! # Ok(())
//! # }
//! ```
use crate::error::positive;
use crate::{
    Cubic, Error, FlashOptions, FlashPhase, FlashResult, Phase, PreparedState, R, Result, numerics,
};

#[derive(Debug, Clone, Copy)]
pub struct EquilibriumProperties {
    pub volume: f64,
    pub enthalpy: f64,
    pub entropy: f64,
    pub internal_energy: f64,
    /// dH/dT along equilibrium at fixed P. None for pure two-phase states
    /// (latent heat at fixed T) or a numerically singular response matrix.
    pub heat_capacity: Option<f64>,
}

#[derive(Debug, Clone, Copy)]
pub struct CaloricFlashOptions {
    pub temperature_guess: Option<f64>,
    pub minimum_temperature: f64,
    pub maximum_temperature: f64,
    pub max_iterations: usize,
    /// Absolute tolerance for residual H/(RT) or S/R, as in the source solvers.
    pub scaled_tolerance: f64,
    pub tp: FlashOptions,
}
impl Default for CaloricFlashOptions {
    fn default() -> Self {
        Self {
            temperature_guess: None,
            minimum_temperature: 80.0,
            maximum_temperature: 999.0,
            max_iterations: 100,
            scaled_tolerance: 1e-10,
            tp: FlashOptions::default(),
        }
    }
}

#[derive(Clone, Copy)]
enum Mode {
    Enthalpy,
    Entropy,
}
impl Mode {
    fn value(self, p: &EquilibriumProperties) -> f64 {
        match self {
            Self::Enthalpy => p.enthalpy,
            Self::Entropy => p.entropy,
        }
    }
    fn scale(self, t: f64) -> f64 {
        match self {
            Self::Enthalpy => R * t,
            Self::Entropy => R,
        }
    }
    fn slope(self, p: &EquilibriumProperties, t: f64) -> Option<f64> {
        p.heat_capacity.map(|cp| match self {
            Self::Enthalpy => cp,
            Self::Entropy => cp / t,
        })
    }
}

impl Cubic {
    /// Properties of the equilibrium result, including the change of phase
    /// amounts/compositions in the equilibrium heat capacity.
    pub fn equilibrium_properties(&self, flash: &FlashResult) -> Result<EquilibriumProperties> {
        let t = flash.temperature;
        let p = flash.pressure;
        positive(p, "pressure")?;
        if flash.phase != FlashPhase::TwoPhase {
            let state = self.prepare(t, &flash.feed)?;
            let v = state.volume(p, Phase::MinimumGibbs)?;
            let (h, s, cp) = phase_caloric(&state, v)?;
            return Ok(EquilibriumProperties {
                volume: v,
                enthalpy: h,
                entropy: s,
                internal_energy: h - p * v,
                heat_capacity: if cp.is_finite() && cp > 0.0 {
                    Some(cp)
                } else {
                    None
                },
            });
        }
        let beta = flash.vapor_fraction;
        if !beta.is_finite() || !(0.0..=1.0).contains(&beta) {
            return Err(Error::InvalidInput("physical two-phase fraction required"));
        }
        let l = self.prepare(t, &flash.liquid)?;
        let v = self.prepare(t, &flash.vapor)?;
        let vl = l.volume(p, Phase::Liquid)?;
        let vv = v.volume(p, Phase::Vapor)?;
        let (hl, sl, cpl) = phase_caloric(&l, vl)?;
        let (hv, sv, cpv) = phase_caloric(&v, vv)?;
        let h = (1.0 - beta) * hl + beta * hv;
        let s = (1.0 - beta) * sl + beta * sv;
        let volume = (1.0 - beta) * vl + beta * vv;
        let active: Vec<_> = flash
            .feed
            .iter()
            .enumerate()
            .filter(|(_, z)| **z > 0.0)
            .map(|(i, _)| i)
            .collect();
        let capacity = if active.len() > 1 && beta > 0.0 && beta < 1.0 {
            let n = flash.feed.len();
            let mut dtl = vec![0.0; n];
            let mut dtv = vec![0.0; n];
            let mut dp = vec![0.0; n];
            let mut dnl = vec![0.0; n * n];
            let mut dnv = vec![0.0; n * n];
            l.fugacity_derivatives_into(vl, &mut dtl, &mut dp, &mut dnl)?;
            v.fugacity_derivatives_into(vv, &mut dtv, &mut dp, &mut dnv)?;
            let m = active.len();
            let mut matrix = vec![0.0; m * m];
            let delta: Vec<_> = active.iter().map(|&i| dtv[i] - dtl[i]).collect();
            for (ii, &i) in active.iter().enumerate() {
                for (jj, &j) in active.iter().enumerate() {
                    matrix[ii * m + jj] =
                        (dnv[i * n + j] - 1.0) / beta + (dnl[i * n + j] - 1.0) / (1.0 - beta);
                }
                matrix[ii * m + ii] +=
                    1.0 / (beta * flash.vapor[i]) + 1.0 / ((1.0 - beta) * flash.liquid[i]);
            }
            let mut response = delta.clone();
            if numerics::solve_linear_in_place(&mut matrix, &mut response).is_ok() {
                let frozen = (1.0 - beta) * cpl + beta * cpv;
                let cp = frozen
                    + R * t * t * delta.iter().zip(&response).map(|(a, b)| a * b).sum::<f64>();
                if cp.is_finite() && cp > 0.0 {
                    Some(cp)
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        };
        Ok(EquilibriumProperties {
            volume,
            enthalpy: h,
            entropy: s,
            internal_energy: h - p * volume,
            heat_capacity: capacity,
        })
    }

    /// Recover T and phase split from P (Pa), normalized z, H (J/mol).
    /// Default temperature bounds follow the source: 80 to 999 K.
    pub fn ph_flash(&self, p: f64, z: &[f64], enthalpy: f64) -> Result<FlashResult> {
        self.ph_flash_with_options(p, z, enthalpy, CaloricFlashOptions::default())
    }

    pub fn ph_flash_with_options(
        &self,
        p: f64,
        z: &[f64],
        enthalpy: f64,
        options: CaloricFlashOptions,
    ) -> Result<FlashResult> {
        self.caloric_flash(p, z, enthalpy, Mode::Enthalpy, options)
    }

    /// Recover T and phase split from P (Pa), normalized z, S (J/(mol K)).
    pub fn ps_flash(&self, p: f64, z: &[f64], entropy: f64) -> Result<FlashResult> {
        self.ps_flash_with_options(p, z, entropy, CaloricFlashOptions::default())
    }

    pub fn ps_flash_with_options(
        &self,
        p: f64,
        z: &[f64],
        entropy: f64,
        options: CaloricFlashOptions,
    ) -> Result<FlashResult> {
        self.caloric_flash(p, z, entropy, Mode::Entropy, options)
    }

    fn caloric_flash(
        &self,
        p: f64,
        z: &[f64],
        target: f64,
        mode: Mode,
        options: CaloricFlashOptions,
    ) -> Result<FlashResult> {
        positive(p, "pressure")?;
        positive(options.minimum_temperature, "minimum temperature")?;
        positive(options.maximum_temperature, "maximum temperature")?;
        positive(options.scaled_tolerance, "caloric residual tolerance")?;
        if !target.is_finite()
            || options.minimum_temperature >= options.maximum_temperature
            || options.max_iterations == 0
            || options
                .temperature_guess
                .is_some_and(|t| !t.is_finite() || t <= 0.0)
        {
            return Err(Error::InvalidInput(
                "finite caloric specification and valid temperature bounds required",
            ));
        }
        let critical = self.pseudo_critical(z)?;
        let initial = options
            .temperature_guess
            .unwrap_or(critical.temperature)
            .clamp(options.minimum_temperature, options.maximum_temperature);
        let active: Vec<_> = z
            .iter()
            .enumerate()
            .filter(|(_, z)| **z > 0.0)
            .map(|(i, _)| i)
            .collect();
        if active.len() == 1 {
            return self.pure_caloric_flash(p, z, target, mode, options, active[0]);
        }
        // Solve the homogeneous caloric problem first. It is cheap and is
        // accepted only after a TP stability/split check at that temperature.
        let homogeneous = |t| -> Result<(FlashResult, EquilibriumProperties)> {
            let flash = self.homogeneous_flash(t, p, z, FlashPhase::Single)?;
            let properties = self.equilibrium_properties(&flash)?;
            Ok((flash, properties))
        };
        // A minimum-Gibbs homogeneous branch jumps at fixed-composition
        // coexistence. A target in that jump has no homogeneous root; do not
        // spend the full equilibrium iteration budget bisecting it.
        let mut homogeneous_options = options;
        homogeneous_options.max_iterations = homogeneous_options.max_iterations.min(16);
        if let Ok(mut candidate) =
            solve_scalar(target, mode, initial, homogeneous_options, &homogeneous)
        {
            let stable = self.tp_flash_with_options(candidate.temperature, p, z, options.tp)?;
            let props = self.equilibrium_properties(&stable)?;
            if (mode.value(&props) - target).abs()
                <= options.scaled_tolerance * mode.scale(candidate.temperature)
            {
                candidate = stable;
                return Ok(self.caloric_phase_convention(
                    candidate,
                    critical.pressure,
                    critical.volume,
                    mode,
                ));
            }
        }
        let equilibrium = |t| -> Result<(FlashResult, EquilibriumProperties)> {
            let flash = self.tp_flash_with_options(t, p, z, options.tp)?;
            let properties = self.equilibrium_properties(&flash)?;
            Ok((flash, properties))
        };
        let flash = solve_scalar(target, mode, initial, options, &equilibrium)?;
        Ok(self.caloric_phase_convention(flash, critical.pressure, critical.volume, mode))
    }

    fn homogeneous_flash(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        phase: FlashPhase,
    ) -> Result<FlashResult> {
        let state = self.prepare(t, z)?;
        let root = match phase {
            FlashPhase::Liquid => Phase::Liquid,
            FlashPhase::Vapor => Phase::Vapor,
            _ => Phase::MinimumGibbs,
        };
        let volume = state.volume(p, root)?;
        Ok(FlashResult {
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
                Some(volume)
            },
            vapor_volume: if phase == FlashPhase::Vapor {
                Some(volume)
            } else {
                None
            },
            iterations: 0,
            fugacity_residual: 0.0,
        })
    }

    fn caloric_phase_convention(
        &self,
        mut flash: FlashResult,
        pc: f64,
        vc: f64,
        mode: Mode,
    ) -> FlashResult {
        if matches!(mode, Mode::Entropy) {
            return flash;
        }
        if flash.phase == FlashPhase::Single && flash.pressure < pc {
            let volume = flash
                .liquid_volume
                .or(flash.vapor_volume)
                .unwrap_or(f64::NAN);
            if volume > vc {
                flash.phase = FlashPhase::Vapor;
                flash.vapor_fraction = 1.0;
                flash.liquid_fraction = 0.0;
                flash.vapor = flash.feed.clone();
                flash.liquid.fill(0.0);
                flash.vapor_volume = Some(volume);
                flash.liquid_volume = None;
            } else {
                flash.phase = FlashPhase::Liquid;
                flash.vapor_fraction = 0.0;
                flash.liquid_fraction = 1.0;
            }
        }
        if flash.phase != FlashPhase::TwoPhase {
            // Match the PH wrapper's single-phase output initialization.
            flash.liquid.fill(0.0);
            flash.vapor.fill(0.0);
            if flash.phase == FlashPhase::Single {
                flash.liquid_fraction = 0.0;
            }
        }
        flash
    }

    #[allow(clippy::too_many_arguments)]
    fn pure_caloric_flash(
        &self,
        p: f64,
        z: &[f64],
        target: f64,
        mode: Mode,
        mut options: CaloricFlashOptions,
        index: usize,
    ) -> Result<FlashResult> {
        let mut pure = Cubic::new(self.kind(), &[self.components()[index].ident])?;
        let alpha = self.alpha(index)?;
        pure.set_alpha(
            0,
            alpha.correlation,
            &alpha.parameters[..alpha.correlation.parameter_count()],
        )?;
        let critical = pure.pseudo_critical(&[1.0])?;
        let mut phase = FlashPhase::Single;
        if p < critical.pressure {
            let sat = pure.saturation_temperature(p)?;
            let state = self.prepare(sat.temperature, z)?;
            let liquid = match mode {
                Mode::Enthalpy => state.enthalpy(sat.liquid_volume)?,
                Mode::Entropy => state.entropy(sat.liquid_volume)?,
            };
            let vapor = match mode {
                Mode::Enthalpy => state.enthalpy(sat.vapor_volume)?,
                Mode::Entropy => state.entropy(sat.vapor_volume)?,
            };
            if target >= liquid && target <= vapor {
                if sat.temperature < options.minimum_temperature
                    || sat.temperature > options.maximum_temperature
                {
                    return Err(Error::InvalidInput(
                        "caloric specification outside temperature bounds",
                    ));
                }
                let beta = (target - liquid) / (vapor - liquid);
                return Ok(FlashResult {
                    temperature: sat.temperature,
                    pressure: p,
                    feed: z.to_vec(),
                    phase: FlashPhase::TwoPhase,
                    liquid: z.to_vec(),
                    vapor: z.to_vec(),
                    vapor_fraction: beta,
                    liquid_fraction: 1.0 - beta,
                    liquid_volume: Some(sat.liquid_volume),
                    vapor_volume: Some(sat.vapor_volume),
                    iterations: 0,
                    fugacity_residual: 0.0,
                });
            }
            if target < liquid {
                phase = FlashPhase::Liquid;
                options.maximum_temperature = options.maximum_temperature.min(sat.temperature);
            } else {
                phase = FlashPhase::Vapor;
                options.minimum_temperature = options.minimum_temperature.max(sat.temperature);
            }
        }
        if options.minimum_temperature >= options.maximum_temperature {
            return Err(Error::InvalidInput(
                "caloric specification outside temperature bounds",
            ));
        }
        let initial = options
            .temperature_guess
            .unwrap_or(critical.temperature)
            .clamp(options.minimum_temperature, options.maximum_temperature);
        let evaluate = |t| -> Result<(FlashResult, EquilibriumProperties)> {
            let flash = self.homogeneous_flash(t, p, z, phase)?;
            let state = self.prepare(t, z)?;
            let volume = if phase == FlashPhase::Vapor {
                flash.vapor_volume
            } else {
                flash.liquid_volume
            }
            .ok_or(Error::NoPhysicalRoot)?;
            let (h, s, cp) = phase_caloric(&state, volume)?;
            Ok((
                flash,
                EquilibriumProperties {
                    volume,
                    enthalpy: h,
                    entropy: s,
                    internal_energy: h - p * volume,
                    heat_capacity: if cp.is_finite() && cp > 0.0 {
                        Some(cp)
                    } else {
                        None
                    },
                },
            ))
        };
        let mut flash = solve_scalar(target, mode, initial, options, &evaluate)?;
        // Pure PH/PS interfaces return z in both composition arrays.
        flash.liquid = z.to_vec();
        flash.vapor = z.to_vec();
        Ok(flash)
    }
}

// Evaluate shared ideal integrals and residual derivatives once per phase.
fn phase_caloric(state: &PreparedState<'_>, volume: f64) -> Result<(f64, f64, f64)> {
    let residual = state.residual(volume)?;
    let ideal = state.ideal(state.pressure(volume)?.value)?;
    Ok((
        ideal.enthalpy + residual.enthalpy,
        ideal.entropy + residual.entropy,
        ideal.enthalpy_dt + residual.cp,
    ))
}

fn solve_scalar<F>(
    target: f64,
    mode: Mode,
    initial: f64,
    options: CaloricFlashOptions,
    evaluate: &F,
) -> Result<FlashResult>
where
    F: Fn(f64) -> Result<(FlashResult, EquilibriumProperties)>,
{
    let (mut flash, mut properties) = evaluate(initial)?;
    let mut residual = mode.value(&properties) - target;
    let mut low = initial;
    let mut high = initial;
    let mut step = (0.1 * initial).max(5.0);
    let direction = if residual > 0.0 { -1.0 } else { 1.0 };
    let mut bracketed = false;
    for _ in 0..32 {
        if residual.abs() <= options.scaled_tolerance * mode.scale(flash.temperature) {
            return Ok(flash);
        }
        let next = (flash.temperature + direction * step)
            .clamp(options.minimum_temperature, options.maximum_temperature);
        if next == flash.temperature {
            return Err(Error::InvalidInput(
                "caloric specification outside temperature bounds",
            ));
        }
        let (new_flash, new_properties) = evaluate_nearby(next, flash.temperature, evaluate)?;
        let next = new_flash.temperature;
        let new_residual = mode.value(&new_properties) - target;
        if new_residual.signum() != residual.signum() {
            low = flash.temperature.min(next);
            high = flash.temperature.max(next);
            flash = new_flash;
            properties = new_properties;
            residual = new_residual;
            bracketed = true;
            break;
        }
        step = (next - flash.temperature).abs() * 2.0;
        flash = new_flash;
        properties = new_properties;
        residual = new_residual;
    }
    if !bracketed {
        return Err(Error::NotConverged("caloric temperature bracket"));
    }
    for iteration in 0..options.max_iterations {
        if residual.abs() <= options.scaled_tolerance * mode.scale(flash.temperature) {
            flash.iterations += iteration + 1;
            return Ok(flash);
        }
        if residual < 0.0 {
            low = flash.temperature;
        } else {
            high = flash.temperature;
        }
        if high - low <= 8.0 * f64::EPSILON * flash.temperature {
            return Err(Error::NotConverged("caloric bracket has no closing root"));
        }
        let newton = mode
            .slope(&properties, flash.temperature)
            .map(|cp| flash.temperature - residual / cp)
            .unwrap_or(f64::NAN);
        let next = if newton > low && newton < high {
            newton
        } else {
            (low + high) * 0.5
        };
        // If a trial enters a region the TP solver cannot resolve, approach
        // from the opposite valid bracket endpoint instead of getting trapped
        // repeatedly next to the same unresolved interval.
        let anchor = if residual < 0.0 { high } else { low };
        let point = evaluate_nearby(next, anchor, evaluate)?;
        flash = point.0;
        properties = point.1;
        residual = mode.value(&properties) - target;
    }
    Err(Error::NotConverged("caloric temperature solve"))
}

fn evaluate_nearby<F>(
    mut candidate: f64,
    previous: f64,
    evaluate: &F,
) -> Result<(FlashResult, EquilibriumProperties)>
where
    F: Fn(f64) -> Result<(FlashResult, EquilibriumProperties)>,
{
    let mut last_error = Error::NotConverged("caloric trial state");
    for _ in 0..16 {
        match evaluate(candidate) {
            Ok(value) => return Ok(value),
            Err(error @ (Error::NotConverged(_) | Error::NoPhysicalRoot)) => last_error = error,
            Err(error) => return Err(error),
        }
        candidate = 0.5 * (candidate + previous);
        if candidate == previous {
            break;
        }
    }
    Err(last_error)
}
