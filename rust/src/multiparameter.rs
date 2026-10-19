//! High-precision NIST Multiparameter Equations of State.
//! Implements fundamental reduced Helmholtz energy formulations for Propane (C3),
//! R134a, Normal Hydrogen (N-H2), Ortho Hydrogen (O-H2), and Para Hydrogen (P-H2).
//! Ported from `src/multiparameter_base.f90` and `src/multiparameter_*.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.

use crate::error::positive;
use crate::multiparameter_data::{
    AncillarySatDensity, IdealGasHelmholtz, MultiparameterRecord, lookup_multiparameter_record,
};
use crate::{Error, Phase, Pressure, R, Result};

/// Reduced Helmholtz jet containing the value and scaled derivatives.
/// Follows the notation of `multiparameter_base.f90`.
#[derive(Debug, Clone, Copy, Default)]
pub struct HelmholtzJet {
    /// Dimensionless Helmholtz energy: `alpha`
    pub alpha: f64,
    /// Dimensionless first derivative w.r.t tau: `tau * d_alpha / d_tau`
    pub d_tau: f64,
    /// Dimensionless second derivative w.r.t tau: `tau^2 * d2_alpha / d_tau^2`
    pub d_tautau: f64,
    /// Dimensionless first derivative w.r.t delta: `delta * d_alpha / d_delta`
    pub d_delta: f64,
    /// Dimensionless second derivative w.r.t delta: `delta^2 * d2_alpha / d_delta^2`
    pub d_deltadelta: f64,
    /// Dimensionless cross-derivative: `delta * tau * d2_alpha / (d_delta d_tau)`
    pub d_deltatau: f64,
}

/// NIST Multiparameter equation of state for a pure fluid.
#[derive(Debug, Clone, Copy)]
pub struct Multiparameter {
    record: &'static MultiparameterRecord,
}

impl Multiparameter {
    /// Create a new multiparameter EOS for the specified fluid name.
    /// Supports `"C3"`, `"Propane"`, `"R134A"`, `"N-H2"`, `"H2"`, `"O-H2"`, `"P-H2"`.
    pub fn new(fluid_name: &str) -> Result<Self> {
        let record = lookup_multiparameter_record(fluid_name)
            .ok_or_else(|| Error::UnknownComponent(fluid_name.to_owned()))?;
        Ok(Self { record })
    }

    /// Access the underlying fluid record metadata.
    pub fn record(&self) -> &'static MultiparameterRecord {
        self.record
    }

    /// Critical temperature (K).
    pub fn tc(&self) -> f64 {
        self.record.tc
    }

    /// Critical pressure (Pa).
    pub fn pc(&self) -> f64 {
        self.record.pc
    }

    /// Critical density (mol/m^3).
    pub fn rc(&self) -> f64 {
        self.record.rc
    }

    /// Molar mass (kg/mol).
    pub fn molar_mass(&self) -> f64 {
        self.record.molar_mass
    }

    /// Prepare state calculations at a specified temperature (K).
    /// Precomputes temperature-dependent powers and factors for allocation-free evaluation.
    pub fn prepare(&self, temperature: f64) -> Result<PreparedMultiparameter<'_>> {
        positive(temperature, "temperature")?;
        let tau = self.record.tc / temperature;
        let mut prep = PreparedMultiparameter {
            eos: self,
            temperature,
            tau,
            tau_cache: Vec::new(),
            tau_exp_cache: Vec::new(),
            tau_ee_cache: Vec::new(),
            beta_tau_diff: Vec::new(),
        };

        // Precompute temperature factors
        prep.tau_cache = self
            .record
            .pol_terms
            .iter()
            .map(|p| p.n * tau.powf(p.t))
            .collect();

        prep.tau_exp_cache = self
            .record
            .exp_terms
            .iter()
            .map(|e| e.n * tau.powf(e.t))
            .collect();

        prep.tau_ee_cache = self
            .record
            .expexp_terms
            .iter()
            .map(|ee| {
                let diff_t = tau - ee.gamma;
                let t_fac = (-ee.beta * diff_t.powi(2)).exp();
                ee.n * tau.powf(ee.t) * t_fac
            })
            .collect();

        prep.beta_tau_diff = self
            .record
            .expexp_terms
            .iter()
            .map(|ee| 2.0 * ee.beta * tau * (tau - ee.gamma))
            .collect();

        Ok(prep)
    }
}

/// Prepared multiparameter state at fixed temperature.
/// Executes all volume, pressure, caloric, and acoustic calculations with zero heap allocations.
#[derive(Debug, Clone)]
pub struct PreparedMultiparameter<'a> {
    eos: &'a Multiparameter,
    temperature: f64,
    tau: f64,
    tau_cache: Vec<f64>,
    tau_exp_cache: Vec<f64>,
    tau_ee_cache: Vec<f64>,
    beta_tau_diff: Vec<f64>,
}

impl<'a> PreparedMultiparameter<'a> {
    /// Return the prepared temperature (K).
    pub fn temperature(&self) -> f64 {
        self.temperature
    }

    /// Return the reduced temperature ratio `tau = Tc / T`.
    pub fn tau(&self) -> f64 {
        self.tau
    }

    /// Calculate the ideal gas reduced Helmholtz energy jet at molar volume `v` (m^3/mol).
    pub fn ideal_jet(&self, v: f64) -> Result<HelmholtzJet> {
        positive(v, "volume")?;
        let rho = 1.0 / v;
        let delta = rho / self.eos.record.rc;
        let tau = self.tau;

        match self.eos.record.ideal {
            IdealGasHelmholtz::Einstein { c_tau, a, b, v: vi } => {
                let mut jet = HelmholtzJet {
                    alpha: delta.ln() + c_tau * tau.ln() + a[0] + a[1] * tau,
                    d_delta: 1.0,
                    d_deltadelta: -1.0,
                    d_tau: c_tau + a[1] * tau,
                    d_tautau: -c_tau,
                    d_deltatau: 0.0,
                };
                for (&bi, &v_coef) in b.iter().zip(vi.iter()) {
                    let exp_term = (-bi * tau).exp();
                    let one_minus_exp = 1.0 - exp_term;
                    jet.alpha += v_coef * one_minus_exp.ln();
                    let denom = (bi * tau).exp() - 1.0;
                    jet.d_tau += tau * v_coef * bi / denom;
                    let denom2 = denom * denom;
                    jet.d_tautau -= tau * tau * v_coef * bi * bi * (bi * tau).exp() / denom2;
                }
                Ok(jet)
            }
            IdealGasHelmholtz::PowerLaw { a } => {
                let tau_m05 = tau.powf(-0.5);
                let tau_m075 = tau.powf(-0.75);
                let jet = HelmholtzJet {
                    alpha: a[0]
                        + a[1] * tau
                        + a[2] * tau.ln()
                        + delta.ln()
                        + a[3] * tau_m05
                        + a[4] * tau_m075,
                    d_delta: 1.0,
                    d_deltadelta: -1.0,
                    d_tau: a[1] * tau + a[2] - 0.5 * a[3] * tau_m05 - 0.75 * a[4] * tau_m075,
                    d_tautau: -a[2] + 0.75 * a[3] * tau_m05 + 1.3125 * a[4] * tau_m075,
                    d_deltatau: 0.0,
                };
                Ok(jet)
            }
        }
    }

    /// Calculate the residual reduced Helmholtz energy jet at molar volume `v` (m^3/mol).
    pub fn residual_jet(&self, v: f64) -> Result<HelmholtzJet> {
        positive(v, "volume")?;
        let rho = 1.0 / v;
        let delta = rho / self.eos.record.rc;

        let mut jet = HelmholtzJet::default();

        // 1. Polynomial terms
        for (i, p) in self.eos.record.pol_terms.iter().enumerate() {
            let term = self.tau_cache[i] * delta.powi(p.d);
            jet.alpha += term;
            jet.d_delta += p.d as f64 * term;
            jet.d_deltadelta += (p.d * (p.d - 1)) as f64 * term;
            jet.d_tau += p.t * term;
            jet.d_tautau += p.t * (p.t - 1.0) * term;
            jet.d_deltatau += (p.d as f64) * p.t * term;
        }

        // 2. Exponential terms
        for (i, e) in self.eos.record.exp_terms.iter().enumerate() {
            let delta_l = delta.powi(e.l);
            let term = self.tau_exp_cache[i] * delta.powi(e.d) * (-delta_l).exp();
            let d_eff = e.d as f64 - (e.l as f64) * delta_l;
            jet.alpha += term;
            jet.d_delta += d_eff * term;
            jet.d_deltadelta += (d_eff * (d_eff - 1.0) - (e.l * e.l) as f64 * delta_l) * term;
            jet.d_tau += e.t * term;
            jet.d_tautau += e.t * (e.t - 1.0) * term;
            jet.d_deltatau += d_eff * e.t * term;
        }

        // 3. Double-exponential (Gaussian bell) terms
        for (i, ee) in self.eos.record.expexp_terms.iter().enumerate() {
            let diff_d = delta - ee.epsilon;
            let d_fac = (-ee.eta * diff_d.powi(2)).exp();
            let term = self.tau_ee_cache[i] * delta.powi(ee.d) * d_fac;
            let th_delta = ee.d as f64 - 2.0 * ee.eta * delta * diff_d;
            let th_tau = ee.t - self.beta_tau_diff[i];
            jet.alpha += term;
            jet.d_delta += th_delta * term;
            jet.d_deltadelta +=
                (th_delta.powi(2) - ee.d as f64 - 2.0 * ee.eta * delta.powi(2)) * term;
            jet.d_tau += th_tau * term;
            jet.d_tautau += (th_tau.powi(2) - ee.t - 2.0 * ee.beta * self.tau.powi(2)) * term;
            jet.d_deltatau += th_delta * th_tau * term;
        }

        Ok(jet)
    }

    /// Calculate pressure (Pa) and derivatives dP/dT and dP/dv at molar volume `v` (m^3/mol).
    pub fn pressure(&self, v: f64) -> Result<Pressure> {
        positive(v, "volume")?;
        let rho = 1.0 / v;
        let res = self.residual_jet(v)?;
        let t = self.temperature;

        let p_val = rho * R * t * (1.0 + res.d_delta);
        let dp_drho = R * t * (1.0 + 2.0 * res.d_delta + res.d_deltadelta);
        let dp_dv = -rho * rho * dp_drho;
        let dp_dt = rho * R * (1.0 + res.d_delta - res.d_deltatau);

        Ok(Pressure {
            value: p_val,
            dt: dp_dt,
            dv: dp_dv,
        })
    }

    /// Molar internal energy U [J/mol].
    pub fn internal_energy(&self, v: f64) -> Result<f64> {
        let ideal = self.ideal_jet(v)?;
        let res = self.residual_jet(v)?;
        Ok(R * self.temperature * (ideal.d_tau + res.d_tau))
    }

    /// Molar entropy S [J/(mol K)].
    pub fn entropy(&self, v: f64) -> Result<f64> {
        let ideal = self.ideal_jet(v)?;
        let res = self.residual_jet(v)?;
        Ok(R * ((ideal.d_tau + res.d_tau) - (ideal.alpha + res.alpha)))
    }

    /// Molar enthalpy H [J/mol].
    pub fn enthalpy(&self, v: f64) -> Result<f64> {
        let u = self.internal_energy(v)?;
        let p = self.pressure(v)?.value;
        Ok(u + p * v)
    }

    /// Isochoric heat capacity Cv [J/(mol K)].
    pub fn heat_capacity_v(&self, v: f64) -> Result<f64> {
        let ideal = self.ideal_jet(v)?;
        let res = self.residual_jet(v)?;
        Ok(-R * (ideal.d_tautau + res.d_tautau))
    }

    /// Isobaric heat capacity Cp [J/(mol K)].
    pub fn heat_capacity_p(&self, v: f64) -> Result<f64> {
        let cv = self.heat_capacity_v(v)?;
        let res = self.residual_jet(v)?;
        let num = (1.0 + res.d_delta - res.d_deltatau).powi(2);
        let den = 1.0 + 2.0 * res.d_delta + res.d_deltadelta;
        Ok(cv + R * num / den)
    }

    /// Speed of sound w [m/s].
    pub fn speed_of_sound(&self, v: f64) -> Result<f64> {
        let ideal = self.ideal_jet(v)?;
        let res = self.residual_jet(v)?;
        let num = (1.0 + res.d_delta - res.d_deltatau).powi(2);
        let den_caloric = ideal.d_tautau + res.d_tautau;
        let adim = 1.0 + 2.0 * res.d_delta + res.d_deltadelta - num / den_caloric;
        let w2 = adim * R * self.temperature / self.eos.record.molar_mass;
        if w2 < 0.0 {
            return Err(Error::NoPhysicalRoot);
        }
        Ok(w2.sqrt())
    }

    /// Natural logarithm of fugacity coefficient `ln(phi)`.
    pub fn ln_fugacity(&self, v: f64) -> Result<f64> {
        let _rho = 1.0 / v;
        let res = self.residual_jet(v)?;
        let z = 1.0 + res.d_delta;
        if z <= 0.0 {
            return Err(Error::NoPhysicalRoot);
        }
        Ok(res.alpha + res.d_delta - z.ln())
    }

    /// Estimate saturation reduced density `delta_sat = rho / rho_c` using ancillary equations.
    pub fn sat_delta_estimate(&self, phase: Phase) -> Result<f64> {
        let theta = 1.0 - 1.0 / self.tau;
        let th_pos = theta.max(0.0);

        match self.eos.record.ancillary {
            AncillarySatDensity::LemmonSpan {
                n_liq,
                expo_liq,
                n_vap,
                expo_vap,
            } => match phase {
                Phase::Liquid => {
                    let mut sum_val = 0.0;
                    for (&n, &exp) in n_liq.iter().zip(expo_liq.iter()) {
                        sum_val += n * th_pos.powf(exp);
                    }
                    Ok(1.0 + sum_val)
                }
                Phase::Vapor => {
                    let mut sum_val = 0.0;
                    for (&n, &exp) in n_vap.iter().zip(expo_vap.iter()) {
                        sum_val += n * th_pos.powf(exp);
                    }
                    Ok(sum_val.exp())
                }
                Phase::MinimumGibbs => Ok(1.0),
            },
            AncillarySatDensity::TillnerRoth {
                n_liq,
                expo_liq,
                n_vap,
                expo_vap,
            } => {
                let scale = 1.0 / (self.eos.record.molar_mass * self.eos.record.rc);
                match phase {
                    Phase::Liquid => {
                        let mut sum_val = 0.0;
                        for (&n, &exp) in n_liq.iter().zip(expo_liq.iter()) {
                            sum_val += n * th_pos.powf(exp);
                        }
                        Ok(sum_val * scale)
                    }
                    Phase::Vapor => {
                        let mut sum_val = 0.0;
                        for (&n, &exp) in n_vap.iter().zip(expo_vap.iter()) {
                            sum_val += n * th_pos.powf(exp);
                        }
                        Ok(sum_val.exp() * 516.86 * scale)
                    }
                    Phase::MinimumGibbs => Ok(1.0),
                }
            }
            AncillarySatDensity::TriplePointRatio {
                liq_ratio,
                vap_ratio,
            } => match phase {
                Phase::Liquid => Ok(liq_ratio),
                Phase::Vapor => Ok(vap_ratio),
                Phase::MinimumGibbs => Ok(1.0),
            },
        }
    }

    /// Robust 1D Newton solver for molar volume `v` (m^3/mol) at specified pressure `p_spec` (Pa) and phase.
    pub fn volume(&self, p_spec: f64, phase: Phase) -> Result<f64> {
        positive(p_spec, "pressure")?;
        let rc = self.eos.record.rc;

        // Determine initial density estimate
        let mut rho = if self.temperature > self.eos.record.tc {
            // Supercritical: ideal gas
            p_spec / (R * self.temperature)
        } else {
            match phase {
                Phase::Liquid => {
                    let est_delta = self.sat_delta_estimate(Phase::Liquid)?;
                    1.01 * rc * est_delta
                }
                Phase::Vapor | Phase::MinimumGibbs => p_spec / (R * self.temperature),
            }
        };

        let max_iter = 100;
        let releps_p = 1e-8;
        let releps_rho = 1e-8;

        for _ in 0..max_iter {
            let v = 1.0 / rho;
            let p_res = self.pressure(v)?;
            let p_curr = p_res.value;
            let dp_drho = -v * v * p_res.dv;

            let diff_p = p_curr - p_spec;
            if diff_p.abs() < releps_p * p_spec + 1e-4 {
                return Ok(v);
            }

            if dp_drho <= 0.0 {
                // If derivative is non-positive, shift density
                rho = match phase {
                    Phase::Liquid => rho * 1.05,
                    _ => rho * 0.95,
                };
                continue;
            }

            let delta_rho = diff_p / dp_drho;
            let mut step = delta_rho;
            if step > 0.5 * rho {
                step = 0.5 * rho;
            } else if step < -0.5 * rho {
                step = -0.5 * rho;
            }

            let rho_next = rho - step;
            if rho_next <= 0.0 {
                rho *= 0.5;
                continue;
            }

            if (rho_next - rho).abs() < releps_rho * rho {
                return Ok(1.0 / rho_next);
            }

            rho = rho_next;
        }

        Ok(1.0 / rho)
    }
}
