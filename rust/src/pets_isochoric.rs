//! Pure PeTS UV/SV recovery along the equilibrium isochore.
//!
//! ```
//! use thermopack_rs::Pets;
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let model = Pets::new()?;
//! let initial = model.equilibrium_tv(100.0, 1e-3)?;
//! let recovered = model.uv_flash(initial.properties.internal_energy, 1e-3)?;
//! assert!((recovered.temperature - 100.0).abs() < 1e-6);
//! # Ok(())
//! # }
//! ```
use super::*;
use crate::VolumeFlashOptions;

/// Pure equilibrium at specified temperature and molar volume. Volume fixes
/// vapor quality in coexistence. `heat_capacity_v` includes phase redistribution.
#[derive(Debug, Clone)]
pub struct PetsEquilibriumTv {
    pub flash: FlashResult,
    pub properties: EquilibriumProperties,
    /// dU/dT at constant total molar volume, in J/(mol K). At a saturated
    /// endpoint this is the derivative from the coexistence side.
    pub heat_capacity_v: f64,
}

#[derive(Clone, Copy)]
enum Specification {
    Energy,
    Entropy,
}
impl Specification {
    fn value(self, state: &PetsEquilibriumTv) -> f64 {
        match self {
            Self::Energy => state.properties.internal_energy,
            Self::Entropy => state.properties.entropy,
        }
    }
    fn scale(self, t: f64) -> f64 {
        match self {
            Self::Energy => R * t,
            Self::Entropy => R,
        }
    }
    fn slope(self, state: &PetsEquilibriumTv) -> f64 {
        match self {
            Self::Energy => state.heat_capacity_v,
            Self::Entropy => state.heat_capacity_v / state.flash.temperature,
        }
    }
}

impl Pets {
    /// Fluid equilibrium at positive T and molar V. Uses the same packing
    /// interval and branch-discovery limits as the pressure-root solver.
    pub fn equilibrium_tv(&self, t: f64, volume: f64) -> Result<PetsEquilibriumTv> {
        positive(volume, "molar volume")?;
        let prepared = self.prepare(t)?;
        if prepared.packing / volume > 0.75 {
            return Err(Error::InvalidInput(
                "PeTS equilibrium packing fraction exceeds 0.75",
            ));
        }
        let critical = self.critical_point()?;
        let eq = if t < critical.temperature {
            Some(prepared.equilibrium()?)
        } else {
            None
        };
        if let Some(eq) = &eq {
            let sat = eq.saturation_pressure()?;
            if volume >= sat.liquid_volume && volume <= sat.vapor_volume {
                let beta = (volume - sat.liquid_volume) / (sat.vapor_volume - sat.liquid_volume);
                let liquid = prepared.caloric_tv(sat.liquid_volume, 1.)?;
                let vapor = prepared.caloric_tv(sat.vapor_volume, 1.)?;
                let weight = |l, v| (1. - beta) * l + beta * v;
                let u = weight(liquid.internal_energy, vapor.internal_energy);
                let entropy = weight(liquid.entropy, vapor.entropy);
                // dV_i/dT on saturation=(P_sat' - P_T_i)/P_V_i.
                // Eliminating beta' using constant total V and Clapeyron gives
                // Cv_eq = weighted(Cv_i - T*(P_T_i-P_sat')^2/P_V_i).
                // This avoids subtracting caloric reference offsets.
                let response = |state: &PetsTvProperties| {
                    let dpt = state.pressure.dt - sat.dpdt;
                    state.heat_capacity_v - t * (dpt * dpt) / state.pressure.dv
                };
                let cv = weight(response(&liquid), response(&vapor));
                positive(cv, "positive finite PeTS equilibrium Cv")?;
                let properties = EquilibriumProperties {
                    volume,
                    enthalpy: u + sat.pressure * volume,
                    entropy,
                    internal_energy: u,
                    heat_capacity: None,
                };
                if !properties.enthalpy.is_finite() || !entropy.is_finite() || !u.is_finite() {
                    return Err(Error::InvalidInput(
                        "non-finite PeTS coexistence properties",
                    ));
                }
                return Ok(PetsEquilibriumTv {
                    flash: split(sat, beta),
                    properties,
                    heat_capacity_v: cv,
                });
            }
        }
        let state = prepared.caloric_tv(volume, 1.)?;
        positive(state.pressure.value, "positive PeTS equilibrium pressure")?;
        positive(state.heat_capacity_v, "positive finite PeTS equilibrium Cv")?;
        let stiffness_roundoff = 32. * f64::EPSILON * R * t / (volume * volume);
        if state.pressure.dv > stiffness_roundoff {
            return Err(Error::NoPhysicalRoot);
        }
        if let Some(eq) = eq {
            let v = eq.volume(state.pressure.value, Phase::MinimumGibbs)?;
            if (v - volume).abs() > 1e-7 * volume {
                return Err(Error::NotConverged(
                    "PeTS TV state is not the minimum-Gibbs branch",
                ));
            }
        }
        Ok(PetsEquilibriumTv {
            flash: homogeneous(t, state.pressure.value, volume, FlashPhase::Single),
            properties: properties(volume, state),
            heat_capacity_v: state.heat_capacity_v,
        })
    }

    /// Recover fluid equilibrium from molar U (J/mol) and V (m3/mol).
    pub fn uv_flash(&self, energy: f64, volume: f64) -> Result<FlashResult> {
        self.uv_flash_with_options(energy, volume, VolumeFlashOptions::default())
    }
    /// Recover fluid equilibrium from molar S (J/(mol K)) and V (m3/mol).
    pub fn sv_flash(&self, entropy: f64, volume: f64) -> Result<FlashResult> {
        self.sv_flash_with_options(entropy, volume, VolumeFlashOptions::default())
    }
    pub fn uv_flash_with_options(
        &self,
        energy: f64,
        volume: f64,
        options: VolumeFlashOptions,
    ) -> Result<FlashResult> {
        self.volume_flash(energy, volume, Specification::Energy, options)
    }
    pub fn sv_flash_with_options(
        &self,
        entropy: f64,
        volume: f64,
        options: VolumeFlashOptions,
    ) -> Result<FlashResult> {
        self.volume_flash(entropy, volume, Specification::Entropy, options)
    }

    fn volume_flash(
        &self,
        target: f64,
        volume: f64,
        spec: Specification,
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
                .is_some_and(|v| !v.is_finite() || v <= 0.)
            || options
                .pressure_guess
                .is_some_and(|v| !v.is_finite() || v <= 0.)
        {
            return Err(Error::InvalidInput(
                "finite volume specification and valid flash bounds/guesses required",
            ));
        }
        let critical = self.critical_point()?;
        let initial = options
            .temperature_guess
            .unwrap_or(critical.temperature)
            .clamp(options.minimum_temperature, options.maximum_temperature);
        let mut state = self
            .equilibrium_tv(initial, volume)
            .or_else(|_| {
                self.equilibrium_tv(
                    critical
                        .temperature
                        .clamp(options.minimum_temperature, options.maximum_temperature),
                    volume,
                )
            })
            .or_else(|_| self.equilibrium_tv(options.maximum_temperature, volume))?;
        let mut error = spec.value(&state) - target;
        let mut low = state.flash.temperature;
        let mut high = low;
        let direction = if error > 0. { -1. } else { 1. };
        let mut step = 0.1 * low;
        let mut bracketed = false;
        let accept = |mut state: PetsEquilibriumTv, iterations| -> Result<FlashResult> {
            if state.flash.pressure < options.minimum_pressure
                || state.flash.pressure > options.maximum_pressure
            {
                return Err(Error::InvalidInput(
                    "PeTS volume-flash pressure outside bounds",
                ));
            }
            state.flash.iterations = iterations;
            Ok(state.flash)
        };
        for _ in 0..40 {
            if error.abs() <= options.scaled_tolerance * spec.scale(state.flash.temperature) {
                return accept(state, 0);
            }
            let bounded_step = step.min(1.25 * (error / spec.slope(&state)).abs());
            let bounded_step = if direction < 0. {
                bounded_step.min(0.3 * state.flash.temperature)
            } else {
                bounded_step
            };
            let next = (state.flash.temperature + direction * bounded_step)
                .clamp(options.minimum_temperature, options.maximum_temperature);
            if next == state.flash.temperature {
                return Err(Error::InvalidInput(
                    "PeTS caloric specification outside temperature bounds",
                ));
            }
            let new = self.volume_trial(next, state.flash.temperature, volume)?;
            let new_error = spec.value(&new) - target;
            if new_error.signum() != error.signum() {
                low = state.flash.temperature.min(new.flash.temperature);
                high = state.flash.temperature.max(new.flash.temperature);
                state = new;
                error = new_error;
                bracketed = true;
                break;
            }
            step = (new.flash.temperature - state.flash.temperature).abs() * 2.;
            state = new;
            error = new_error;
        }
        if !bracketed {
            return Err(Error::NotConverged("PeTS volume-flash temperature bracket"));
        }
        for i in 0..options.max_iterations {
            if error.abs() <= options.scaled_tolerance * spec.scale(state.flash.temperature) {
                return accept(state, i + 1);
            }
            if error < 0. {
                low = state.flash.temperature;
            } else {
                high = state.flash.temperature;
            }
            if high - low <= 8. * f64::EPSILON * state.flash.temperature {
                return Err(Error::NotConverged(
                    "PeTS volume-flash bracket has no closing root",
                ));
            }
            let newton = state.flash.temperature - error / spec.slope(&state);
            let next = if newton > low && newton < high {
                newton
            } else {
                low + (high - low) * 0.5
            };
            let anchor = if error < 0. { high } else { low };
            state = self.volume_trial(next, anchor, volume)?;
            error = spec.value(&state) - target;
        }
        Err(Error::NotConverged("PeTS volume-flash temperature solve"))
    }

    fn volume_trial(&self, mut t: f64, anchor: f64, volume: f64) -> Result<PetsEquilibriumTv> {
        let mut error = Error::NotConverged("PeTS volume-flash trial");
        for _ in 0..16 {
            match self.equilibrium_tv(t, volume) {
                Ok(state) => return Ok(state),
                Err(e) => error = e,
            }
            t = (t + anchor) * 0.5;
            if t == anchor {
                break;
            }
        }
        Err(error)
    }
}
