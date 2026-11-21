//! Pure multiparameter equilibrium on a fixed-volume path.
use super::*;
use crate::caloric::{EquilibriumProperties, VolumeFlashOptions};
use crate::{FlashPhase, FlashResult};

/// Equilibrium at fixed temperature and molar volume. In the two-phase region
/// the reported heat capacity includes the quality change required to hold V.
#[derive(Debug, Clone)]
pub struct MultiparameterEquilibriumTv {
    pub flash: FlashResult,
    pub properties: EquilibriumProperties,
    pub heat_capacity_v: f64,
}

#[derive(Clone, Copy)]
enum Spec {
    Energy,
    Entropy,
}
impl Spec {
    fn get(self, q: &MultiparameterEquilibriumTv) -> f64 {
        match self {
            Self::Energy => q.properties.internal_energy,
            Self::Entropy => q.properties.entropy,
        }
    }
    fn scale(self, model: &Multiparameter, t: f64) -> f64 {
        model.gas_constant() * if matches!(self, Self::Energy) { t } else { 1. }
    }
    fn slope(self, state: &MultiparameterEquilibriumTv) -> f64 {
        match self {
            Self::Energy => state.heat_capacity_v,
            Self::Entropy => state.heat_capacity_v / state.flash.temperature,
        }
    }
}

impl Multiparameter {
    /// Generic source UV/SV solution bounds. The component triple temperature
    /// only adjusts UV initial guesses; callers may supply narrower bounds.
    pub fn volume_flash_options(&self) -> VolumeFlashOptions {
        if self.gas_constant() == 1. {
            let (_, maximum_pressure) = self.lennard_jones_limits().unwrap_or((11., 70.));
            return VolumeFlashOptions {
                minimum_temperature: 0.3 * self.tc(),
                maximum_temperature: 3. * self.tc(),
                minimum_pressure: 1e-8,
                maximum_pressure,
                ..Default::default()
            };
        }
        VolumeFlashOptions {
            minimum_temperature: 80.,
            maximum_temperature: 999.,
            minimum_pressure: 10.,
            maximum_pressure: 1e8,
            ..Default::default()
        }
    }

    /// Evaluate pure liquid/vapor branch equilibrium at (T,V).
    pub fn equilibrium_tv(&self, t: f64, volume: f64) -> Result<MultiparameterEquilibriumTv> {
        positive(volume, "molar volume")?;
        let prep = self.prepare(t)?;
        if t < self.tc() {
            let sat = self.saturation_pressure(t)?;
            if volume >= sat.liquid_volume && volume <= sat.vapor_volume {
                let gap = sat.vapor_volume - sat.liquid_volume;
                if !gap.is_finite() || gap <= 0. {
                    return Err(Error::NotConverged("multiparameter coexistence volume gap"));
                }
                let beta = ((volume - sat.liquid_volume) / (sat.vapor_volume - sat.liquid_volume))
                    .clamp(0., 1.);
                let l = flash::properties(&prep, sat.liquid_volume)?;
                let v = flash::properties(&prep, sat.vapor_volume)?;
                let response = |vol: f64| -> Result<f64> {
                    let pressure = prep.pressure(vol)?;
                    let pv = pressure.dv;
                    if !pv.is_finite() || pv >= 0. {
                        return Err(Error::NoPhysicalRoot);
                    }
                    // dP_sat/dT is Clapeyron-exact for the coexistence pair.
                    let pt = pressure.dt;
                    let cv = lcv(self, &prep, vol)?;
                    let dpt = pt - sat.dpdt;
                    Ok(cv - t * (dpt * dpt) / pv)
                };
                let cv =
                    (1. - beta) * response(sat.liquid_volume)? + beta * response(sat.vapor_volume)?;
                positive(cv, "positive multiparameter equilibrium Cv")?;
                let mix = |a: f64, b: f64| (1. - beta) * a + beta * b;
                let props = EquilibriumProperties {
                    volume,
                    enthalpy: mix(l.enthalpy, v.enthalpy),
                    entropy: mix(l.entropy, v.entropy),
                    internal_energy: mix(l.internal_energy, v.internal_energy),
                    heat_capacity: None,
                };
                return Ok(MultiparameterEquilibriumTv {
                    flash: flash::split(sat, beta),
                    properties: props,
                    heat_capacity_v: cv,
                });
            }
        }
        let props = flash::properties(&prep, volume)?;
        let pressure = prep.pressure(volume)?.value;
        positive(pressure, "positive multiparameter pressure")?;
        let pressure_state = prep.pressure(volume)?;
        let stiffness_tol = 32. * f64::EPSILON * self.gas_constant() * t;
        if pressure_state.dv >= -stiffness_tol / (volume * volume) {
            return Err(Error::NoPhysicalRoot);
        }
        // Outside the coexistence interval accept a mechanically stable pure
        // branch; this is not an exhaustive global-stability proof.
        let phase = if t >= self.tc() {
            FlashPhase::Single
        } else if volume < 1. / self.rc() {
            FlashPhase::Liquid
        } else {
            FlashPhase::Vapor
        };
        let flash = flash::homogeneous(t, pressure, volume, phase);
        let cv = lcv(self, &prep, volume)?;
        positive(cv, "positive multiparameter Cv")?;
        Ok(MultiparameterEquilibriumTv {
            flash,
            properties: props,
            heat_capacity_v: cv,
        })
    }

    pub fn uv_flash(&self, energy: f64, volume: f64) -> Result<FlashResult> {
        self.uv_flash_with_options(energy, volume, self.volume_flash_options())
    }
    pub fn sv_flash(&self, entropy: f64, volume: f64) -> Result<FlashResult> {
        self.sv_flash_with_options(entropy, volume, self.volume_flash_options())
    }
    pub fn uv_flash_with_options(
        &self,
        energy: f64,
        volume: f64,
        options: VolumeFlashOptions,
    ) -> Result<FlashResult> {
        let mut options = options;
        if self.gas_constant() != 1. {
            // twoPhaseUVsingleComp uses component Ttr to repair an invalid
            // initial guess. Its inner singleCompUV_Tv/singleCompUV solvers
            // bound the solution by tpTmin/tpTmax, allowing roots below Ttr.
            let seed_min = options
                .minimum_temperature
                .max(self.source_triple_temperature());
            if seed_min < options.maximum_temperature
                && options.temperature_guess.is_some_and(|guess| {
                    guess.is_finite()
                        && guess > 0.
                        && (guess < seed_min || guess > options.maximum_temperature)
                })
            {
                options.temperature_guess = Some(0.5 * (seed_min + options.maximum_temperature));
            }
        }
        self.volume_flash(energy, volume, Spec::Energy, options)
    }
    pub fn sv_flash_with_options(
        &self,
        entropy: f64,
        volume: f64,
        options: VolumeFlashOptions,
    ) -> Result<FlashResult> {
        self.volume_flash(entropy, volume, Spec::Entropy, options)
    }

    fn volume_flash(
        &self,
        target: f64,
        volume: f64,
        spec: Spec,
        o: VolumeFlashOptions,
    ) -> Result<FlashResult> {
        positive(volume, "molar volume")?;
        if !target.is_finite()
            || !o.minimum_temperature.is_finite()
            || !o.maximum_temperature.is_finite()
            || !o.minimum_pressure.is_finite()
            || !o.maximum_pressure.is_finite()
            || o.minimum_temperature <= 0.
            || o.minimum_temperature >= o.maximum_temperature
            || o.max_iterations == 0
            || !o.scaled_tolerance.is_finite()
            || o.scaled_tolerance <= 0.
            || o.minimum_pressure <= 0.
            || o.minimum_pressure >= o.maximum_pressure
            || o.temperature_guess
                .is_some_and(|x| !x.is_finite() || x <= 0.)
            || o.pressure_guess.is_some_and(|x| !x.is_finite() || x <= 0.)
        {
            return Err(Error::InvalidInput(
                "finite UV/SV specification and valid bounds required",
            ));
        }
        let lo = o.minimum_temperature;
        let hi = o.maximum_temperature;
        // A pressure seed is validated but unnecessary on a fixed-volume
        // path: the native EOS supplies P(T,V) directly for every trial.
        let initial = o.temperature_guess.unwrap_or(self.tc()).clamp(lo, hi);
        let mut state = self
            .equilibrium_tv(initial, volume)
            .or_else(|_| self.equilibrium_tv(self.tc().clamp(lo, hi), volume))
            .or_else(|_| self.equilibrium_tv(hi, volume))?;
        validate_trial(&state, target, spec)?;
        let mut residual = spec.get(&state) - target;
        let mut low = state.flash.temperature;
        let mut high = low;
        let direction = if residual > 0. { -1. } else { 1. };
        let mut step = 0.1 * state.flash.temperature;
        let mut bracketed = false;

        // Expand from the supplied (or critical) seed. The analytic equilibrium
        // Cv gives dU/dT|V; dS/dT|V is Cv/T. Each trial stays on a local valid
        // branch through bounded temperature backtracking. Pressure limits apply
        // to the accepted root, including when the initial trial exceeds them.
        for _ in 0..48 {
            if closes(residual, spec, self, state.flash.temperature, &o) {
                return accept_flash(state.flash, residual, spec, self, &o, 0);
            }
            let slope = spec.slope(&state);
            if !slope.is_finite() || slope <= 0. {
                return Err(Error::NotConverged(
                    "nonpositive MEOS isochoric caloric slope",
                ));
            }
            let mut bounded = step.min(1.25 * (residual / slope).abs());
            if direction < 0. {
                bounded = bounded.min(0.3 * state.flash.temperature);
            }
            let next = (state.flash.temperature + direction * bounded).clamp(lo, hi);
            if next == state.flash.temperature {
                return Err(Error::NotConverged(
                    "MEOS caloric specification outside temperature bounds",
                ));
            }
            let trial = volume_trial(self, next, state.flash.temperature, volume, target, spec)?;
            let next_residual = spec.get(&trial) - target;
            if next_residual.signum() != residual.signum() {
                low = state.flash.temperature.min(trial.flash.temperature);
                high = state.flash.temperature.max(trial.flash.temperature);
                state = trial;
                residual = next_residual;
                bracketed = true;
                break;
            }
            step = (trial.flash.temperature - state.flash.temperature).abs() * 2.;
            state = trial;
            residual = next_residual;
        }
        if !bracketed {
            return Err(Error::NotConverged(
                "MEOS isochoric caloric temperature bracket",
            ));
        }

        for iteration in 1..=o.max_iterations {
            if closes(residual, spec, self, state.flash.temperature, &o) {
                return accept_flash(state.flash, residual, spec, self, &o, iteration);
            }
            if residual < 0. {
                low = state.flash.temperature;
            } else {
                high = state.flash.temperature;
            }
            let t = state.flash.temperature;
            if high - low <= 8. * f64::EPSILON * t {
                return Err(Error::NotConverged(
                    "MEOS caloric temperature precision without closure",
                ));
            }
            let slope = spec.slope(&state);
            let newton = if slope.is_finite() && slope > 0. {
                t - residual / slope
            } else {
                f64::NAN
            };
            let next = if newton > low && newton < high {
                newton
            } else {
                low + 0.5 * (high - low)
            };
            let anchor = if residual < 0. { high } else { low };
            let trial = volume_trial(self, next, anchor, volume, target, spec)?;
            let next_residual = spec.get(&trial) - target;
            if trial.flash.temperature == t && !closes(next_residual, spec, self, t, &o) {
                return Err(Error::NotConverged(
                    "MEOS caloric temperature precision without closure",
                ));
            }
            state = trial;
            residual = next_residual;
        }
        Err(Error::NotConverged(
            "MEOS isochoric caloric temperature solve",
        ))
    }

    fn source_triple_temperature(&self) -> f64 {
        let component = self
            .lennard_jones_parameters()
            .map_or(self.record.ident, |p| p.component);
        crate::component::component(component)
            .map_or(0., |record| record.source_triple_temperature())
    }
}

fn closes(r: f64, spec: Spec, model: &Multiparameter, t: f64, o: &VolumeFlashOptions) -> bool {
    r.is_finite() && r.abs() <= o.scaled_tolerance * spec.scale(model, t)
}

fn validate_trial(q: &MultiparameterEquilibriumTv, target: f64, spec: Spec) -> Result<()> {
    // Pressure bounds constrain the accepted solution. A cold dense target
    // can be valid even when the initial critical-temperature trial exceeds
    // Pmax. The isochoric continuation must be allowed to descend to it.
    if !spec.get(q).is_finite() || !target.is_finite() {
        return Err(Error::InvalidInput("nonfinite MEOS UV/SV caloric state"));
    }
    Ok(())
}

fn accept_flash(
    mut flash: FlashResult,
    residual: f64,
    spec: Spec,
    model: &Multiparameter,
    o: &VolumeFlashOptions,
    iterations: usize,
) -> Result<FlashResult> {
    let scale = spec.scale(model, flash.temperature);
    if !residual.is_finite() || residual.abs() > o.scaled_tolerance * scale {
        return Err(Error::NotConverged("MEOS UV/SV caloric closure"));
    }
    if flash.pressure < o.minimum_pressure || flash.pressure > o.maximum_pressure {
        return Err(Error::InvalidInput(
            "MEOS UV/SV state pressure outside bounds",
        ));
    }
    flash.iterations = iterations;
    Ok(flash)
}

fn volume_trial(
    model: &Multiparameter,
    mut t: f64,
    anchor: f64,
    volume: f64,
    target: f64,
    spec: Spec,
) -> Result<MultiparameterEquilibriumTv> {
    let mut error = Error::NotConverged("MEOS UV/SV isochoric trial");
    for _ in 0..20 {
        match model.equilibrium_tv(t, volume) {
            Ok(q) => {
                validate_trial(&q, target, spec)?;
                return Ok(q);
            }
            Err(e) => error = e,
        }
        let next = 0.5 * (t + anchor);
        if next == t || next == anchor {
            break;
        }
        t = next;
    }
    Err(error)
}

fn lcv(model: &Multiparameter, prep: &PreparedMultiparameter<'_>, volume: f64) -> Result<f64> {
    let i = prep.ideal_jet(volume)?;
    let r = prep.residual_jet(volume)?;
    let cv = -model.gas_constant() * (i.d_tautau + r.d_tautau);
    if cv.is_finite() {
        Ok(cv)
    } else {
        Err(Error::InvalidInput("nonfinite multiparameter Cv"))
    }
}
