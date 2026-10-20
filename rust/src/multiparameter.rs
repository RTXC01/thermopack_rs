//! Pure NIST and Lennard-Jones multiparameter equations of state.
//! Five NIST fluids, full/truncated-shifted LJ for four noble gases, and reduced
//! LJTS share native Helmholtz jets and allocation-free temperature preparation.
//! Ported from `src/multiparameter_base.f90` and `src/multiparameter_*.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.

use crate::error::positive;
use crate::multiparameter_data::{
    AncillarySatDensity, IdealGasHelmholtz, MultiparameterRecord, lookup_multiparameter_record,
};
use crate::{Error, Phase, Pressure, R, Result};

#[path = "lj_meos_data.rs"]
mod lj_data;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LennardJonesKind {
    Full,
    /// Source LJTS truncated-shifted fit.
    TruncatedShifted,
}

#[derive(Debug, Clone, Copy)]
pub struct LennardJonesParameters {
    pub component: &'static str,
    pub kind: LennardJonesKind,
    pub sigma: f64,
    pub epsilon_over_k: f64,
    pub maximum_temperature: f64,
    pub maximum_pressure: f64,
    /// The original constructor stores g/mol in a field documented as kg/mol.
    pub baseline_molar_mass: f64,
}

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

impl HelmholtzJet {
    fn add(self, b: Self) -> Self {
        Self {
            alpha: self.alpha + b.alpha,
            d_tau: self.d_tau + b.d_tau,
            d_tautau: self.d_tautau + b.d_tautau,
            d_delta: self.d_delta + b.d_delta,
            d_deltadelta: self.d_deltadelta + b.d_deltadelta,
            d_deltatau: self.d_deltatau + b.d_deltatau,
        }
    }
    fn finite(self) -> bool {
        [
            self.alpha,
            self.d_tau,
            self.d_tautau,
            self.d_delta,
            self.d_deltadelta,
            self.d_deltatau,
        ]
        .iter()
        .all(|v| v.is_finite())
    }
}

#[derive(Debug, Clone, Copy)]
pub struct MultiparameterState {
    pub ideal: HelmholtzJet,
    pub residual: HelmholtzJet,
    pub pressure: Pressure,
    pub internal_energy: f64,
    pub enthalpy: f64,
    pub entropy: f64,
    pub helmholtz_energy: f64,
    pub gibbs_energy: f64,
    pub heat_capacity_v: f64,
    pub heat_capacity_p: f64,
    /// Absent when Z<=0 (the source logarithm is non-finite).
    pub ln_fugacity: Option<f64>,
    /// Absent when the squared sound speed is negative or non-finite.
    pub speed_of_sound: Option<f64>,
    /// Source meos method, including LJ constructor's mass-unit convention.
    pub baseline_speed_of_sound: Option<f64>,
}

/// NIST Multiparameter equation of state for a pure fluid.
#[derive(Debug, Clone, Copy)]
pub struct Multiparameter {
    record: MultiparameterRecord,
    gas_constant: f64,
    lj: Option<LennardJonesParameters>,
}

impl Multiparameter {
    /// Create a new multiparameter EOS for the specified fluid name.
    /// Supports `"C3"`, `"Propane"`, `"R134A"`, `"N-H2"`, `"H2"`, `"O-H2"`, `"P-H2"`.
    pub fn new(fluid_name: &str) -> Result<Self> {
        let record = lookup_multiparameter_record(fluid_name)
            .ok_or_else(|| Error::UnknownComponent(fluid_name.to_owned()))?;
        Ok(Self {
            record: *record,
            gas_constant: R,
            lj: None,
        })
    }

    /// Access the underlying fluid record metadata.
    pub fn record(&self) -> &MultiparameterRecord {
        &self.record
    }

    pub fn gas_constant(&self) -> f64 {
        self.gas_constant
    }

    pub fn lennard_jones_parameters(&self) -> Option<LennardJonesParameters> {
        self.lj
    }

    /// Source meos field; LJ noble-gas constructors leave it in g/mol.
    pub fn baseline_molar_mass(&self) -> f64 {
        self.lj
            .map_or(self.record.molar_mass, |p| p.baseline_molar_mass)
    }

    pub fn lennard_jones_limits(&self) -> Option<(f64, f64)> {
        self.lj
            .map(|p| (p.maximum_temperature, p.maximum_pressure))
            .or({
                if self.gas_constant == 1. {
                    Some(lj_data::LJTS_LIMITS)
                } else {
                    None
                }
            })
    }

    /// Source LJ_MEOS / LJTS_MEOS constructor for NE, AR, KR and XE.
    /// Scaling and reference offsets are owned by this model.
    pub fn lennard_jones(component: &str, kind: LennardJonesKind) -> Result<Self> {
        let component = crate::component::component(component)?;
        let truncated = kind == LennardJonesKind::TruncatedShifted;
        let parameters = lj_data::GASES
            .iter()
            .find(|p| p.component == component.ident && p.truncated == truncated)
            .ok_or_else(|| Error::UnknownComponent(component.ident.to_owned()))?;
        let mut eos = Self::lj_template(kind);
        eos.record.tc *= parameters.epsilon;
        eos.record.rc /= 6.02214076e23 * parameters.sigma.powi(3);
        eos.record.t_triple *= parameters.epsilon;
        eos.record.rgas_fit = R;
        eos.record.molar_mass = component.molar_mass * 1e-3;
        eos.gas_constant = R;
        let (max_t, max_p) = if truncated {
            lj_data::LJTS_LIMITS
        } else {
            lj_data::LJ_LIMITS
        };
        eos.lj = Some(LennardJonesParameters {
            component: component.ident,
            kind,
            sigma: parameters.sigma,
            epsilon_over_k: parameters.epsilon,
            maximum_temperature: max_t * parameters.epsilon,
            maximum_pressure: max_p * parameters.epsilon * 1.380649e-23 / parameters.sigma.powi(3),
            baseline_molar_mass: component.molar_mass,
        });
        eos.record.pc = eos
            .prepare(eos.record.tc)?
            .pressure(1. / eos.record.rc)?
            .value;
        let (liq, vap) = {
            let triple = eos.prepare(eos.record.t_triple)?;
            (
                triple.sat_delta_estimate(Phase::Liquid)?,
                triple.sat_delta_estimate(Phase::Vapor)?,
            )
        };
        eos.record.rholiq_triple = eos.record.rc * liq;
        eos.record.rhovap_triple = eos.record.rc * vap;
        Ok(eos)
    }

    /// Original public init_LJTS, before noble-gas SI scaling: R=1, M=1,
    /// Tc=1.086 and rhoc=0.319. Inputs/outputs are in the source reduced units.
    pub fn reduced_lennard_jones_truncated_shifted() -> Self {
        Self::lj_template(LennardJonesKind::TruncatedShifted)
    }

    fn lj_template(kind: LennardJonesKind) -> Self {
        let mut record = match kind {
            LennardJonesKind::Full => lj_data::LJ_TEMPLATE,
            LennardJonesKind::TruncatedShifted => lj_data::LJTS_TEMPLATE,
        };
        if kind == LennardJonesKind::TruncatedShifted {
            // Original init_LJTS entropy alignment at T*=0.8, p*=0.001.
            let delta = 0.001 / (0.8 * record.rc);
            let tau = record.tc / 0.8;
            if let IdealGasHelmholtz::Einstein { ref mut a, .. } = record.ideal {
                let alpha = delta.ln() + 1.5 * tau.ln() + a[1] * tau + a[0];
                a[0] = -(alpha - (1.5 + a[1] * tau));
            }
        }
        Self {
            record,
            gas_constant: 1.,
            lj: None,
        }
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
        if self.record.pol_terms.len() > 8
            || self.record.exp_terms.len() > 13
            || self.record.expexp_terms.len() > 11
        {
            return Err(Error::InvalidInput(
                "multiparameter fit exceeds inline cache capacity",
            ));
        }
        let tau = self.record.tc / temperature;
        let mut prep = PreparedMultiparameter {
            eos: self,
            temperature,
            tau,
            tau_cache: [0.; 8],
            tau_exp_cache: [0.; 13],
            tau_ee_cache: [0.; 11],
            gaussian_temperature: [0.; 11],
            beta_tau_diff: [0.; 11],
            ideal_cache: HelmholtzJet::default(),
        };
        for (out, p) in prep.tau_cache.iter_mut().zip(self.record.pol_terms) {
            *out = p.n * tau.powf(p.t);
        }
        for (out, p) in prep.tau_exp_cache.iter_mut().zip(self.record.exp_terms) {
            *out = p.n * tau.powf(p.t);
        }
        for (i, p) in self.record.expexp_terms.iter().enumerate() {
            prep.tau_ee_cache[i] = p.n * tau.powf(p.t);
            prep.gaussian_temperature[i] = p.beta * (tau - p.gamma).powi(2);
            prep.beta_tau_diff[i] = 2. * p.beta * tau * (tau - p.gamma);
        }
        prep.ideal_cache = prep.ideal_at_delta(1.)?;
        if !tau.is_finite()
            || tau <= 0.
            || prep
                .tau_cache
                .iter()
                .chain(&prep.tau_exp_cache)
                .chain(&prep.tau_ee_cache)
                .chain(&prep.gaussian_temperature)
                .chain(&prep.beta_tau_diff)
                .any(|v| !v.is_finite())
            || !prep.ideal_cache.finite()
        {
            return Err(Error::InvalidInput(
                "non-finite multiparameter temperature functions",
            ));
        }
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
    tau_cache: [f64; 8],
    tau_exp_cache: [f64; 13],
    tau_ee_cache: [f64; 11],
    gaussian_temperature: [f64; 11],
    beta_tau_diff: [f64; 11],
    ideal_cache: HelmholtzJet,
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
        let mut jet = self.ideal_cache;
        jet.alpha += (1. / v / self.eos.record.rc).ln();
        if !jet.finite() {
            return Err(Error::InvalidInput("non-finite ideal Helmholtz jet"));
        }
        Ok(jet)
    }

    fn ideal_at_delta(&self, delta: f64) -> Result<HelmholtzJet> {
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
                    jet.d_tau += tau * v_coef * bi * exp_term / one_minus_exp;
                    jet.d_tautau -=
                        tau * tau * v_coef * bi * bi * exp_term / (one_minus_exp * one_minus_exp);
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
        let delta = 1. / (v * self.eos.record.rc);

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

        let polynomial = jet;
        jet = HelmholtzJet::default();
        // 2. Exponential terms
        for (i, e) in self.eos.record.exp_terms.iter().enumerate() {
            let delta_l = delta.powi(e.l);
            let term = self.tau_exp_cache[i] * delta.powi(e.d) * (-delta_l).exp();
            let d_eff = e.d as f64 - (e.l as f64) * delta_l;
            jet.alpha += term;
            jet.d_delta += d_eff * term;
            jet.d_deltadelta += (d_eff * (e.d as f64 - 1. - e.l as f64 * delta_l)
                - (e.l * e.l) as f64 * delta_l)
                * term;
            jet.d_tau += e.t * term;
            jet.d_tautau += e.t * (e.t - 1.0) * term;
            jet.d_deltatau += d_eff * e.t * term;
        }

        let exponential = jet;
        jet = HelmholtzJet::default();
        // 3. Double-exponential (Gaussian bell) terms
        for (i, ee) in self.eos.record.expexp_terms.iter().enumerate() {
            let diff_d = delta - ee.epsilon;
            let d_fac = (-ee.eta * diff_d.powi(2) - self.gaussian_temperature[i]).exp();
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

        let jet = polynomial.add(exponential).add(jet);
        if !jet.finite() {
            return Err(Error::InvalidInput("non-finite residual Helmholtz jet"));
        }
        Ok(jet)
    }

    /// Share both Helmholtz jets across all homogeneous state properties.
    /// Unstable/extrapolated states are evaluated; no phase selection is implied.
    pub fn state(&self, v: f64) -> Result<MultiparameterState> {
        let ideal = self.ideal_jet(v)?;
        let residual = self.residual_jet(v)?;
        let r = self.eos.gas_constant;
        let t = self.temperature;
        let rho = 1. / v;
        let z = 1. + residual.d_delta;
        let thermal = 1. + residual.d_delta - residual.d_deltatau;
        let mechanical = 1. + 2. * residual.d_delta + residual.d_deltadelta;
        let caloric = ideal.d_tautau + residual.d_tautau;
        let pressure = Pressure {
            value: rho * r * t * z,
            dt: rho * r * thermal,
            dv: -rho * rho * r * t * mechanical,
        };
        let u = r * t * (ideal.d_tau + residual.d_tau);
        let s = r * (ideal.d_tau + residual.d_tau - ideal.alpha - residual.alpha);
        let a = r * t * (ideal.alpha + residual.alpha);
        let cv = -r * caloric;
        let cp = cv + r * thermal * thermal / mechanical;
        let sound_square =
            (mechanical - thermal * thermal / caloric) * r * t / self.eos.record.molar_mass;
        let sound = if sound_square >= 0. && sound_square.is_finite() {
            Some(sound_square.sqrt())
        } else {
            None
        };
        let baseline_sound = if let Some(p) = self.eos.lj {
            let w2 = (mechanical - thermal * thermal / caloric) * r * t / p.baseline_molar_mass;
            if w2 >= 0. && w2.is_finite() {
                Some(w2.sqrt())
            } else {
                None
            }
        } else {
            sound
        };
        let state = MultiparameterState {
            ideal,
            residual,
            pressure,
            internal_energy: u,
            enthalpy: u + pressure.value * v,
            entropy: s,
            helmholtz_energy: a,
            gibbs_energy: a + pressure.value * v,
            heat_capacity_v: cv,
            heat_capacity_p: cp,
            ln_fugacity: if z > 0. {
                Some(residual.alpha + residual.d_delta - z.ln())
            } else {
                None
            },
            speed_of_sound: sound,
            baseline_speed_of_sound: baseline_sound,
        };
        if [
            pressure.value,
            pressure.dt,
            pressure.dv,
            u,
            s,
            a,
            cv,
            cp,
            state.enthalpy,
            state.gibbs_energy,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(Error::InvalidInput("non-finite multiparameter state"));
        }
        Ok(state)
    }

    /// Reduced extensive residual F=Ares/(RT) and all first/second T/V/n derivatives.
    pub fn residual_helmholtz(
        &self,
        volume: f64,
        moles: f64,
    ) -> Result<crate::helmholtz::PureHelmholtz> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        let v = volume / moles;
        self.extensive_jet(self.residual_jet(v)?, volume, moles)
    }

    /// Reduced extensive ideal F=Aideal/(RT) and all first/second T/V/n derivatives.
    pub fn ideal_helmholtz(
        &self,
        volume: f64,
        moles: f64,
    ) -> Result<crate::helmholtz::PureHelmholtz> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        let v = volume / moles;
        self.extensive_jet(self.ideal_jet(v)?, volume, moles)
    }

    fn extensive_jet(
        &self,
        a: HelmholtzJet,
        volume: f64,
        moles: f64,
    ) -> Result<crate::helmholtz::PureHelmholtz> {
        let rho = moles / volume;
        crate::helmholtz::DensityHelmholtz {
            value: a.alpha,
            dt: -a.d_tau / self.temperature,
            drho: a.d_delta / rho,
            dtt: (a.d_tautau + 2. * a.d_tau) / (self.temperature * self.temperature),
            dtrho: -a.d_deltatau / (self.temperature * rho),
            drhorho: a.d_deltadelta / (rho * rho),
        }
        .extensive(volume, moles)
    }

    /// Source meos sound-speed convention, including its LJ mass-unit mismatch.
    pub fn baseline_speed_of_sound(&self, v: f64) -> Result<f64> {
        Ok(self.speed_of_sound(v)?
            * (self.eos.molar_mass() / self.eos.baseline_molar_mass()).sqrt())
    }

    /// Calculate pressure (Pa) and derivatives dP/dT and dP/dv at molar volume `v` (m^3/mol).
    pub fn pressure(&self, v: f64) -> Result<Pressure> {
        positive(v, "volume")?;
        let rho = 1.0 / v;
        let res = self.residual_jet(v)?;
        let t = self.temperature;

        let p_val = rho * self.eos.gas_constant * t * (1.0 + res.d_delta);
        let dp_drho = self.eos.gas_constant * t * (1.0 + 2.0 * res.d_delta + res.d_deltadelta);
        let dp_dv = -rho * rho * dp_drho;
        let dp_dt = rho * self.eos.gas_constant * (1.0 + res.d_delta - res.d_deltatau);

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
        Ok(self.eos.gas_constant * self.temperature * (ideal.d_tau + res.d_tau))
    }

    /// Molar entropy S [J/(mol K)].
    pub fn entropy(&self, v: f64) -> Result<f64> {
        let ideal = self.ideal_jet(v)?;
        let res = self.residual_jet(v)?;
        Ok(self.eos.gas_constant * ((ideal.d_tau + res.d_tau) - (ideal.alpha + res.alpha)))
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
        Ok(-self.eos.gas_constant * (ideal.d_tautau + res.d_tautau))
    }

    /// Isobaric heat capacity Cp [J/(mol K)].
    pub fn heat_capacity_p(&self, v: f64) -> Result<f64> {
        let cv = self.heat_capacity_v(v)?;
        let res = self.residual_jet(v)?;
        let num = (1.0 + res.d_delta - res.d_deltatau).powi(2);
        let den = 1.0 + 2.0 * res.d_delta + res.d_deltadelta;
        Ok(cv + self.eos.gas_constant * num / den)
    }

    /// Speed of sound w [m/s].
    pub fn speed_of_sound(&self, v: f64) -> Result<f64> {
        let ideal = self.ideal_jet(v)?;
        let res = self.residual_jet(v)?;
        let num = (1.0 + res.d_delta - res.d_deltatau).powi(2);
        let den_caloric = ideal.d_tautau + res.d_tautau;
        let adim = 1.0 + 2.0 * res.d_delta + res.d_deltadelta - num / den_caloric;
        let w2 = adim * self.eos.gas_constant * self.temperature / self.eos.record.molar_mass;
        if w2 < 0.0 || !w2.is_finite() {
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
        if phase == Phase::MinimumGibbs {
            return Err(Error::InvalidInput(
                "ancillary density requires liquid or vapor",
            ));
        }
        if theta < 0.
            && matches!(
                self.eos.record.ancillary,
                AncillarySatDensity::LemmonSpan { .. } | AncillarySatDensity::TillnerRoth { .. }
            )
        {
            if self.eos.lj.is_some() || self.eos.gas_constant == 1. {
                return Ok(1.);
            }
            return Err(Error::InvalidInput(
                "non-finite supercritical ancillary density",
            ));
        }
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
            p_spec / (self.eos.gas_constant * self.temperature)
        } else {
            match phase {
                Phase::Liquid => {
                    let est_delta = self.sat_delta_estimate(Phase::Liquid)?;
                    1.01 * rc * est_delta
                }
                Phase::Vapor | Phase::MinimumGibbs => {
                    p_spec / (self.eos.gas_constant * self.temperature)
                }
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
