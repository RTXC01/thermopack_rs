//! Residual TP ideal-mixture rule over pure NIST multiparameter models.
//!
//! This ports internal component kernels only. The original public
//! `NIST_MEOS_MIX` constructor is broken in this fork; this module does not
//! claim parity for total caloric properties, flashes, or global stability.
//! Ported from `src/multiparameter_idealmix.f90` and `multiparameter_base.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.

use crate::error::positive;
use crate::{Error, Multiparameter, Phase, R, Result};

/// Mixture assembled from the NIST pure-fluid multiparameter fits supported by
/// ThermoPack's `NIST_MEOS` allocator.
#[derive(Debug, Clone)]
pub struct MultiparameterIdealMixture {
    components: Vec<Multiparameter>,
}

/// Residual TP properties for one homogeneous ideal-mixture-rule state.
///
/// Composition coefficients are linear coefficients in `sum(x_i q_i)`;
/// they are not constrained mole-number derivatives. The source's fugacity
/// composition derivative matrix is identically zero for this mixing rule.
#[derive(Debug, Clone)]
pub struct MultiparameterIdealMixtureState {
    pub pressure: f64,
    pub temperature: f64,
    pub mole_fractions: Vec<f64>,
    pub component_volumes: Vec<f64>,
    pub component_compressibility_factors: Vec<f64>,
    pub compressibility_factor: f64,
    pub compressibility_factor_dt: f64,
    pub compressibility_factor_dp: f64,
    pub compressibility_factor_composition_coefficients: Vec<f64>,
    pub residual_enthalpy: f64,
    pub residual_enthalpy_dt: f64,
    pub residual_enthalpy_dp: f64,
    pub residual_enthalpy_composition_coefficients: Vec<f64>,
    pub residual_entropy: f64,
    pub residual_entropy_dt: f64,
    pub residual_entropy_dp: f64,
    pub residual_entropy_composition_coefficients: Vec<f64>,
    pub residual_gibbs_energy: f64,
    pub residual_gibbs_energy_dt: f64,
    pub residual_gibbs_energy_dp: f64,
    /// Source-compatible, buggy `calc_multiparameter_idealmix_Gres` dG/dP:
    /// this routine returns the weighted residual Gibbs sum.
    pub baseline_gibbs_dp: f64,
    pub residual_gibbs_composition_coefficients: Vec<f64>,
    pub ln_fugacity_coefficients: Vec<f64>,
    /// `exp(ln_fugacity_coefficients)` when every value is representable.
    pub fugacity_coefficients: Option<Vec<f64>>,
    pub ln_fugacity_coefficients_dt: Vec<f64>,
    pub ln_fugacity_coefficients_dp: Vec<f64>,
    pub ln_fugacity_coefficients_composition_derivatives: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Copy)]
struct ComponentState {
    volume: f64,
    z: f64,
    z_t: f64,
    z_p: f64,
    h: f64,
    h_t: f64,
    h_p: f64,
    s: f64,
    s_t: f64,
    s_p: f64,
    g: f64,
    g_t: f64,
    g_p: f64,
    lnphi: f64,
    lnphi_t: f64,
    lnphi_p: f64,
}

impl MultiparameterIdealMixture {
    /// Construct from the exact NIST components supported by the source
    /// allocator: `C3`, `R134A`, `N-H2`, `O-H2`, and `P-H2`.
    pub fn new(components: &[&str]) -> Result<Self> {
        if components.is_empty() {
            return Err(Error::InvalidInput(
                "mixture requires at least one component",
            ));
        }
        let mut models = Vec::with_capacity(components.len());
        for name in components {
            let model = Multiparameter::new(name)?;
            if !matches!(
                model.record().ident,
                "C3" | "R134A" | "N-H2" | "O-H2" | "P-H2"
            ) {
                return Err(Error::Unsupported(
                    "NIST ideal mixture supports C3, R134A, N-H2, O-H2, and P-H2",
                ));
            }
            models.push(model);
        }
        Ok(Self { components: models })
    }

    pub fn component_count(&self) -> usize {
        self.components.len()
    }

    /// Evaluate the residual TP ideal-mixture rule on a stable homogeneous
    /// liquid or vapor branch. Every constituent is evaluated at common T/P,
    /// including constituents with zero mole fraction, as in the source.
    pub fn residual_tp(
        &self,
        temperature: f64,
        pressure: f64,
        mole_fractions: &[f64],
        phase: Phase,
    ) -> Result<MultiparameterIdealMixtureState> {
        positive(temperature, "temperature")?;
        positive(pressure, "pressure")?;
        if phase == Phase::MinimumGibbs {
            return Err(Error::InvalidInput(
                "ideal mixture requires liquid or vapor phase",
            ));
        }
        if mole_fractions.len() != self.components.len() {
            return Err(Error::InvalidInput(
                "mole-fraction length does not match components",
            ));
        }
        if mole_fractions.iter().any(|x| !x.is_finite() || *x < 0.) {
            return Err(Error::InvalidInput(
                "mole fractions must be finite and nonnegative",
            ));
        }
        let total: f64 = mole_fractions.iter().sum();
        if !total.is_finite() || (total - 1.).abs() > 1e-12 {
            return Err(Error::InvalidInput("mole fractions must sum to one"));
        }

        let mut states = Vec::with_capacity(self.components.len());
        for model in &self.components {
            let prep = model.prepare(temperature)?;
            let volume = prep.volume(pressure, phase)?;
            let residual = prep.residual_jet(volume)?;
            let pressure_value = R * temperature / volume * (1. + residual.d_delta);
            // The original MEOS routine forms PT from the supplied pressure,
            // rather than the EOS pressure reconstructed from the root.
            let pressure_t = pressure / temperature - R * residual.d_deltatau / volume;
            let pressure_v = -R * temperature / (volume * volume)
                * (1. + 2. * residual.d_delta + residual.d_deltadelta);
            let z = pressure * volume / (R * temperature);
            let closure_tolerance =
                2e-8 * pressure + 512. * f64::EPSILON * R * temperature / volume;
            if !volume.is_finite()
                || volume <= 0.
                || !pressure_v.is_finite()
                || pressure_v >= 0.
                || !pressure_value.is_finite()
                || (pressure_value - pressure).abs() > closure_tolerance
            {
                return Err(Error::NotConverged(
                    "multiparameter ideal-mixture component root is not stable and pressure-closing",
                ));
            }
            let v_t = -pressure_t / pressure_v;
            let v_p = 1. / pressure_v;
            // Match source's input-pressure derivative form for Z.
            let z_t = z * (v_t / volume - 1. / temperature);
            let z_p = z * (1. / pressure + v_p / volume);
            let h = R * temperature * residual.d_tau + pressure * volume - R * temperature;
            let h_t = temperature * v_t * pressure_t - R * (residual.d_tautau + 1.);
            let h_p = volume - temperature * v_t;
            let s = R
                * (residual.d_tau - residual.alpha + (pressure * volume / (R * temperature)).ln());
            let s_t = h_t / temperature;
            let s_p = -v_t + R / pressure;
            let g = R * temperature * (residual.alpha + z - 1. - z.ln());
            let g_t = -s;
            let g_p = volume - R * temperature / pressure;
            let lnphi = residual.alpha + residual.d_delta - z.ln();
            let lnphi_t = -h / (R * temperature * temperature);
            let lnphi_p = volume / (R * temperature) - 1. / pressure;
            let values = [
                z, z_t, z_p, h, h_t, h_p, s, s_t, s_p, g, g_t, g_p, lnphi, lnphi_t, lnphi_p,
            ];
            if values.iter().any(|v| !v.is_finite()) {
                return Err(Error::InvalidInput("non-finite ideal-mixture property"));
            }
            states.push(ComponentState {
                volume,
                z,
                z_t,
                z_p,
                h,
                h_t,
                h_p,
                s,
                s_t,
                s_p,
                g,
                g_t,
                g_p,
                lnphi,
                lnphi_t,
                lnphi_p,
            });
        }

        let weighted = |f: fn(&ComponentState) -> f64| -> f64 {
            states
                .iter()
                .zip(mole_fractions)
                .map(|(s, x)| x * f(s))
                .sum()
        };
        let component_volumes = states.iter().map(|s| s.volume).collect();
        let component_z: Vec<_> = states.iter().map(|s| s.z).collect();
        let lnphi: Vec<_> = states.iter().map(|s| s.lnphi).collect();
        let fugacity_coefficients = lnphi.iter().map(|v| v.exp()).collect::<Vec<_>>();
        let fugacity_coefficients = fugacity_coefficients
            .iter()
            .all(|v| v.is_finite())
            .then_some(fugacity_coefficients);
        let n = states.len();
        let output = MultiparameterIdealMixtureState {
            pressure,
            temperature,
            mole_fractions: mole_fractions.to_vec(),
            component_volumes,
            component_compressibility_factors: component_z.clone(),
            compressibility_factor: weighted(|s| s.z),
            compressibility_factor_dt: weighted(|s| s.z_t),
            compressibility_factor_dp: weighted(|s| s.z_p),
            compressibility_factor_composition_coefficients: component_z,
            residual_enthalpy: weighted(|s| s.h),
            residual_enthalpy_dt: weighted(|s| s.h_t),
            residual_enthalpy_dp: weighted(|s| s.h_p),
            residual_enthalpy_composition_coefficients: states.iter().map(|s| s.h).collect(),
            residual_entropy: weighted(|s| s.s),
            residual_entropy_dt: weighted(|s| s.s_t),
            residual_entropy_dp: weighted(|s| s.s_p),
            residual_entropy_composition_coefficients: states.iter().map(|s| s.s).collect(),
            residual_gibbs_energy: weighted(|s| s.g),
            residual_gibbs_energy_dt: weighted(|s| s.g_t),
            residual_gibbs_energy_dp: weighted(|s| s.g_p),
            baseline_gibbs_dp: weighted(|s| s.g),
            residual_gibbs_composition_coefficients: states.iter().map(|s| s.g).collect(),
            ln_fugacity_coefficients: lnphi,
            fugacity_coefficients,
            ln_fugacity_coefficients_dt: states.iter().map(|s| s.lnphi_t).collect(),
            ln_fugacity_coefficients_dp: states.iter().map(|s| s.lnphi_p).collect(),
            ln_fugacity_coefficients_composition_derivatives: vec![vec![0.; n]; n],
        };
        let scalars = [
            output.compressibility_factor,
            output.compressibility_factor_dt,
            output.compressibility_factor_dp,
            output.residual_enthalpy,
            output.residual_enthalpy_dt,
            output.residual_enthalpy_dp,
            output.residual_entropy,
            output.residual_entropy_dt,
            output.residual_entropy_dp,
            output.residual_gibbs_energy,
            output.residual_gibbs_energy_dt,
            output.residual_gibbs_energy_dp,
            output.baseline_gibbs_dp,
        ];
        let vectors_are_finite = output
            .component_volumes
            .iter()
            .chain(&output.component_compressibility_factors)
            .chain(&output.compressibility_factor_composition_coefficients)
            .chain(&output.residual_enthalpy_composition_coefficients)
            .chain(&output.residual_entropy_composition_coefficients)
            .chain(&output.residual_gibbs_composition_coefficients)
            .chain(&output.ln_fugacity_coefficients)
            .chain(&output.ln_fugacity_coefficients_dt)
            .chain(&output.ln_fugacity_coefficients_dp)
            .all(|v| v.is_finite());
        if scalars.iter().any(|v| !v.is_finite()) || !vectors_are_finite {
            return Err(Error::InvalidInput(
                "non-finite weighted ideal-mixture property",
            ));
        }
        Ok(output)
    }
}
