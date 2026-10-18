//! Pure PeTS flashes. No composition iteration is needed for a single species.
//!
//! ```
//! use thermopack_rs::Pets;
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let model = Pets::new()?;
//! let initial = model.tp_flash(300.0, 1e5)?;
//! let state = model.equilibrium_properties(&initial)?;
//! let recovered = model.ph_flash(1e5, state.enthalpy)?;
//! assert!((recovered.temperature - 300.0).abs() < 1e-6);
//! let latent = model.saturated_flash(100.0, 0.4)?;
//! let properties = model.equilibrium_properties(&latent)?;
//! let quality = model.ps_flash(latent.pressure, properties.entropy)?;
//! assert!((quality.vapor_fraction - 0.4).abs() < 1e-8);
//! # Ok(())
//! # }
//! ```
use super::*;
use crate::caloric::{Mode, solve_scalar};
use crate::{CaloricFlashOptions, EquilibriumProperties, FlashPhase, FlashResult, Phase};

fn homogeneous(t: f64, p: f64, v: f64, phase: FlashPhase) -> FlashResult {
    let vapor = phase == FlashPhase::Vapor;
    FlashResult {
        temperature: t,
        pressure: p,
        feed: vec![1.],
        phase,
        liquid: vec![if vapor { 0. } else { 1. }],
        vapor: vec![if vapor { 1. } else { 0. }],
        vapor_fraction: if phase == FlashPhase::Single {
            -1.
        } else if vapor {
            1.
        } else {
            0.
        },
        liquid_fraction: if phase == FlashPhase::Single {
            -1.
        } else if vapor {
            0.
        } else {
            1.
        },
        liquid_volume: if vapor { None } else { Some(v) },
        vapor_volume: if vapor { Some(v) } else { None },
        iterations: 0,
        fugacity_residual: 0.,
    }
}

fn split(sat: crate::Saturation, beta: f64) -> FlashResult {
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

fn properties(v: f64, state: PetsTvProperties) -> EquilibriumProperties {
    EquilibriumProperties {
        volume: v,
        enthalpy: state.enthalpy,
        entropy: state.entropy,
        internal_energy: state.internal_energy,
        heat_capacity: state.heat_capacity_p.filter(|cp| *cp > 0.),
    }
}

impl Pets {
    /// Minimum-Gibbs homogeneous pure-fluid state. TP alone does not determine
    /// vapor quality at coexistence; use `saturated_flash` to specify quality.
    /// A single available pressure root uses the `Single` phase sentinel.
    pub fn tp_flash(&self, t: f64, p: f64) -> Result<FlashResult> {
        let eq = self.prepare_equilibrium(t)?;
        let roots = eq.volumes(p)?;
        let v = eq.minimum_gibbs_volume(p, &roots)?;
        let phase = if roots.len() == 1 {
            FlashPhase::Single
        } else if v == roots[roots.len() - 1] {
            FlashPhase::Vapor
        } else {
            FlashPhase::Liquid
        };
        Ok(homogeneous(t, p, v, phase))
    }

    /// Pure coexistence state at the specified molar vapor fraction, including
    /// the saturated liquid/vapor endpoints. No solid-phase comparison is made.
    pub fn saturated_flash(&self, t: f64, vapor_fraction: f64) -> Result<FlashResult> {
        if !vapor_fraction.is_finite() || !(0. ..=1.).contains(&vapor_fraction) {
            return Err(Error::InvalidInput("physical vapor fraction required"));
        }
        Ok(split(self.saturation_pressure(t)?, vapor_fraction))
    }

    /// Molar V/H/S/U and equilibrium isobaric heat capacity. Phase volumes are
    /// recomputed for this model; cached volumes in `FlashResult` are ignored.
    /// A two-phase input must lie on this model's coexistence curve.
    pub fn equilibrium_properties(&self, flash: &FlashResult) -> Result<EquilibriumProperties> {
        if flash.feed != [1.] {
            return Err(Error::InvalidInput("PeTS requires pure normalized feed"));
        }
        let t = flash.temperature;
        let p = flash.pressure;
        positive(p, "pressure")?;
        let eq = self.prepare_equilibrium(t)?;
        if flash.phase != FlashPhase::TwoPhase {
            let v = eq.volume(p, Phase::MinimumGibbs)?;
            return Ok(properties(v, self.caloric_tv(t, v, 1.)?));
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
        let sat = eq.saturation_pressure()?;
        if (p - sat.pressure).abs() > 1e-8 * p {
            return Err(Error::InvalidInput(
                "pure split is off the coexistence curve",
            ));
        }
        let l = self.caloric_tv(t, sat.liquid_volume, 1.)?;
        let v = self.caloric_tv(t, sat.vapor_volume, 1.)?;
        let weighted = |l, v| (1. - beta) * l + beta * v;
        Ok(EquilibriumProperties {
            volume: weighted(sat.liquid_volume, sat.vapor_volume),
            enthalpy: weighted(l.enthalpy, v.enthalpy),
            entropy: weighted(l.entropy, v.entropy),
            internal_energy: weighted(l.internal_energy, v.internal_energy),
            heat_capacity: None,
        })
    }

    pub fn ph_flash(&self, p: f64, enthalpy: f64) -> Result<FlashResult> {
        self.ph_flash_with_options(p, enthalpy, CaloricFlashOptions::default())
    }
    pub fn ps_flash(&self, p: f64, entropy: f64) -> Result<FlashResult> {
        self.ps_flash_with_options(p, entropy, CaloricFlashOptions::default())
    }
    pub fn ph_flash_with_options(
        &self,
        p: f64,
        enthalpy: f64,
        options: CaloricFlashOptions,
    ) -> Result<FlashResult> {
        self.caloric_flash(p, enthalpy, Mode::Enthalpy, options)
    }
    pub fn ps_flash_with_options(
        &self,
        p: f64,
        entropy: f64,
        options: CaloricFlashOptions,
    ) -> Result<FlashResult> {
        self.caloric_flash(p, entropy, Mode::Entropy, options)
    }

    fn caloric_flash(
        &self,
        p: f64,
        target: f64,
        mode: Mode,
        mut options: CaloricFlashOptions,
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
                .is_some_and(|t| !t.is_finite() || t <= 0.)
        {
            return Err(Error::InvalidInput(
                "finite caloric specification and valid temperature bounds required",
            ));
        }
        let value = |state: &PetsTvProperties| match mode {
            Mode::Enthalpy => state.enthalpy,
            Mode::Entropy => state.entropy,
        };
        let c = self.critical_point()?;
        let mut phase = FlashPhase::Single;
        if p < c.pressure {
            let sat = self.saturation_temperature(p)?;
            let l = value(&self.caloric_tv(sat.temperature, sat.liquid_volume, 1.)?);
            let v = value(&self.caloric_tv(sat.temperature, sat.vapor_volume, 1.)?);
            if v <= l {
                return Err(Error::NotConverged(
                    "positive PeTS latent caloric difference",
                ));
            }
            if target >= l && target <= v {
                if sat.temperature < options.minimum_temperature
                    || sat.temperature > options.maximum_temperature
                {
                    return Err(Error::InvalidInput(
                        "caloric specification outside temperature bounds",
                    ));
                }
                return Ok(split(sat, (target - l) / (v - l)));
            }
            if target < l {
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
            .unwrap_or(c.temperature)
            .clamp(options.minimum_temperature, options.maximum_temperature);
        let evaluate = |t| -> Result<(FlashResult, EquilibriumProperties)> {
            let v = self.volume(t, p, Phase::MinimumGibbs)?;
            let state = self.caloric_tv(t, v, 1.)?;
            Ok((homogeneous(t, p, v, phase), properties(v, state)))
        };
        let mut result = solve_scalar(target, mode, initial, options, &evaluate)?;
        // Pure PH/PS source interfaces populate both composition arrays.
        result.liquid = vec![1.];
        result.vapor = vec![1.];
        Ok(result)
    }
}
