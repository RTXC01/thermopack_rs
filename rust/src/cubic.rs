//! Conventional cubic equations, classic alpha and van der Waals mixing.
//! Ported from cbmix.f90, cbalpha.f90, cubic.f90 and cbhelm.f90.
use crate::component::interaction_with_reference;
use crate::error::positive;
use crate::{Component, Error, R, Result, component, numerics};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CubicKind {
    PengRobinson,
    SoaveRedlichKwong,
    RedlichKwong,
    VanDerWaals,
}

impl CubicKind {
    pub fn name(self) -> &'static str {
        match self {
            Self::PengRobinson => "PR",
            Self::SoaveRedlichKwong => "SRK",
            Self::RedlichKwong => "RK",
            Self::VanDerWaals => "VDW",
        }
    }
    fn deltas(self) -> (f64, f64) {
        match self {
            Self::PengRobinson => (1.0 + 2.0_f64.sqrt(), 1.0 - 2.0_f64.sqrt()),
            Self::SoaveRedlichKwong | Self::RedlichKwong => (1.0, 0.0),
            Self::VanDerWaals => (0.0, 0.0),
        }
    }
}

/// Root choice. MinimumGibbs compares roots at the supplied composition;
/// it does not test stability against phases with a different composition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    Liquid,
    Vapor,
    MinimumGibbs,
}

#[derive(Debug, Clone)]
pub struct Cubic {
    kind: CubicKind,
    components: Vec<Component>,
    ideal_gases: Vec<Result<crate::ideal::IdealGas>>,
    a0: Vec<f64>,
    bi: Vec<f64>,
    alpha: Vec<crate::Alpha>,
    beta: Vec<crate::Beta>,
    variable_beta: bool,
    kij: Vec<f64>,
    lij: Vec<f64>,
    nonzero_kij: Vec<(usize, usize, f64)>,
    translation_parameters: Option<Vec<crate::VolumeTranslation>>,
    translation_enabled: bool,
    pub(crate) caloric_minimum: f64,
    pub(crate) caloric_maximum: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct CubicOptions<'a> {
    pub alpha: Option<crate::AlphaCorrelation>,
    pub parameter_reference: &'a str,
    pub volume_translation: bool,
}
impl Default for CubicOptions<'_> {
    fn default() -> Self {
        Self {
            alpha: None,
            parameter_reference: "DEFAULT",
            volume_translation: false,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Pressure {
    /// Pressure (Pa).
    pub value: f64,
    /// dP/dT at fixed molar volume [Pa/K].
    pub dt: f64,
    /// dP/dv at fixed temperature [Pa mol/m3].
    pub dv: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Residual {
    pub compressibility: f64,
    /// Residual Helmholtz at fixed T, v [J/mol].
    pub helmholtz_tv: f64,
    /// Residual enthalpy relative to ideal gas at T, P [J/mol].
    pub enthalpy: f64,
    /// Residual entropy relative to ideal gas at T, P [J/(mol K)].
    pub entropy: f64,
    pub cv: f64,
    pub cp: f64,
}

/// Molar caloric properties at fixed T,V, independently of pressure sign.
#[derive(Debug, Clone, Copy)]
pub struct TvProperties {
    /// Pressure (Pa), including negative homogeneous values.
    pub pressure: f64,
    /// Molar internal energy [J/mol].
    pub internal_energy: f64,
    /// Molar entropy [J/(mol K)].
    pub entropy: f64,
    /// Physical derivative dU/dT at fixed volume [J/(mol K)].
    pub heat_capacity_v: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Saturation {
    pub temperature: f64,
    pub pressure: f64,
    pub liquid_volume: f64,
    pub vapor_volume: f64,
    /// dPsat/dT [Pa/K] from the Clapeyron equation.
    pub dpdt: f64,
}

/// Critical point of the EOS at fixed composition (not a mixture critical point).
#[derive(Debug, Clone, Copy)]
pub struct PseudoCritical {
    pub temperature: f64,
    pub pressure: f64,
    pub volume: f64,
    pub compressibility: f64,
}

/// Cached temperature and composition parameters. Volume/pressure scans
/// perform no allocations, database lookups or component-pair mixing.
#[derive(Debug, Clone)]
pub struct PreparedState<'a> {
    model: &'a Cubic,
    temperature: f64,
    composition: Vec<f64>,
    a: f64,
    at: f64,
    att: f64,
    b: f64,
    a_partial: Vec<f64>,
    b_partial: Vec<f64>,
    factors: Vec<(f64, f64, f64)>,
    beta: Option<crate::beta::PreparedBeta>,
    translation: Option<crate::translation::PreparedTranslation>,
}

impl Cubic {
    /// Pure-component vapor-liquid equilibrium at fixed temperature.
    /// Brackets coexistence between the EOS spinodals and equates fugacities.
    pub fn saturation_pressure(&self, temperature: f64) -> Result<Saturation> {
        if self.components.len() != 1 {
            return Err(Error::Unsupported("mixture saturation pressure"));
        }
        self.pseudo_saturation_pressure(temperature, &[1.0])
    }

    /// Coexistence with both phase compositions fixed to z, as in PureSat.
    /// For mixtures this is a numerical initial estimate, not VLE saturation.
    pub fn pseudo_saturation_pressure(&self, temperature: f64, z: &[f64]) -> Result<Saturation> {
        let state = self.prepare(temperature, z)?;
        let (_, _, qc) = omega_with_critical_ratio(self.kind);
        let center = qc * state.b;
        if state.eos_pressure(center)?.dv <= 0.0 {
            return Err(Error::InvalidInput(
                "temperature at or above EOS critical temperature",
            ));
        }
        let bisect_spinodal = |mut low: f64, mut high: f64| -> Result<f64> {
            let negative_at_low = state.eos_pressure(low)?.dv < 0.0;
            for _ in 0..80 {
                let mid = (low + high) * 0.5;
                if (state.eos_pressure(mid)?.dv < 0.0) == negative_at_low {
                    low = mid;
                } else {
                    high = mid;
                }
            }
            Ok((low + high) * 0.5)
        };
        let liquid_spinodal = bisect_spinodal(state.b * (1.0 + 1e-8), center)?;
        let mut high_volume = center * 2.0;
        for _ in 0..60 {
            if state.eos_pressure(high_volume)?.dv < 0.0 {
                break;
            }
            high_volume *= 2.0;
        }
        let vapor_spinodal = bisect_spinodal(center, high_volume)?;
        let min_pressure = state.eos_pressure(liquid_spinodal)?.value;
        let max_pressure = state.eos_pressure(vapor_spinodal)?.value;
        let margin = 1e-4 * (max_pressure - min_pressure.max(0.0));
        let mut low = (min_pressure + margin).max(max_pressure * 1e-9).ln();
        let mut high = (max_pressure - margin).ln();
        let evaluate = |log_pressure: f64| -> Result<(f64, f64, f64)> {
            let (volumes, count) = state.eos_volumes(log_pressure.exp())?;
            if count != 2 {
                return Err(Error::NotConverged("distinct saturation roots"));
            }
            let difference = log_pressure.exp() * (volumes[0] - volumes[1]) / (R * temperature)
                - ((volumes[0] - state.b) / (volumes[1] - state.b)).ln()
                - state.a
                    * (state.attraction_integral(volumes[0])
                        - state.attraction_integral(volumes[1]))
                    / (R * temperature);
            Ok((difference, volumes[0], volumes[1]))
        };
        if evaluate(low)?.0 * evaluate(high)?.0 > 0.0 {
            return Err(Error::NotConverged("saturation fugacity bracket"));
        }
        let mut low_f = evaluate(low)?.0;
        for _ in 0..100 {
            let mid = (low + high) * 0.5;
            let (f, vl, vv) = evaluate(mid)?;
            if f.abs() < 1e-13 || high - low < 1e-14 {
                let shift = state.volume_translation().value;
                let (vl, vv) = (vl - shift, vv - shift);
                let latent_heat = state.residual(vv)?.enthalpy - state.residual(vl)?.enthalpy;
                return Ok(Saturation {
                    temperature,
                    pressure: mid.exp(),
                    liquid_volume: vl,
                    vapor_volume: vv,
                    dpdt: latent_heat / (temperature * (vv - vl)),
                });
            }
            if (f > 0.0) == (low_f > 0.0) {
                low = mid;
                low_f = f;
            } else {
                high = mid;
            }
        }
        Err(Error::NotConverged("saturation pressure"))
    }

    /// Inverse pure-fluid saturation pressure, including both phase volumes.
    pub fn saturation_temperature(&self, pressure: f64) -> Result<Saturation> {
        if self.components.len() != 1 {
            return Err(Error::Unsupported("mixture saturation temperature"));
        }
        self.pseudo_saturation_temperature(pressure, &[1.0])
    }

    /// Fixed-composition saturation estimate; mixtures are not in true VLE.
    pub fn pseudo_saturation_temperature(&self, pressure: f64, z: &[f64]) -> Result<Saturation> {
        positive(pressure, "saturation pressure")?;
        let critical = self.pseudo_critical(z)?;
        if pressure >= critical.pressure {
            return Err(Error::InvalidInput(
                "pressure at or above fixed-composition critical pressure",
            ));
        }
        let w: f64 = z
            .iter()
            .zip(&self.components)
            .map(|(x, c)| x * c.acentric_factor)
            .sum();
        let mut low = 0.15 * critical.temperature;
        let mut high = critical.temperature * (1.0 - 1e-9);
        let mut t = (critical.temperature
            / (1.0 - (pressure / critical.pressure).ln() / (5.373 * (1.0 + w))))
            .clamp(0.3 * critical.temperature, 0.999 * critical.temperature);
        for _ in 0..80 {
            let sat = self.pseudo_saturation_pressure(t, z)?;
            let residual = sat.pressure - pressure;
            if residual.abs() <= 1e-10 * pressure || high - low < 1e-10 * critical.temperature {
                return Ok(sat);
            }
            if residual < 0.0 {
                low = t;
            } else {
                high = t;
            }
            let next = t - residual / sat.dpdt;
            t = if next > low && next < high {
                next
            } else {
                (low + high) * 0.5
            };
        }
        Err(Error::NotConverged("saturation temperature"))
    }

    /// Port of cbCalcPseudo for the native conventional cubic family.
    pub fn pseudo_critical(&self, z: &[f64]) -> Result<PseudoCritical> {
        let mut t: f64 = z
            .iter()
            .zip(&self.components)
            .map(|(x, c)| x * c.critical_temperature)
            .sum();
        let (oa, ob, qc) = omega_with_critical_ratio(self.kind);
        for _ in 0..80 {
            let state = self.prepare(t, z)?;
            let f = state.a - (oa / ob) * R * t * state.b;
            let df =
                state.at - (oa / ob) * R * (state.b + t * state.beta.as_ref().map_or(0., |b| b.dt));
            let arg = 1.0 - 2.0 * f * state.att / (df * df);
            let mut step = -2.0 * f / (df * (1.0 + if arg > 0.0 { arg.sqrt() } else { 1.0 }));
            step = step.clamp(-0.5 * t, 0.5 * t);
            t += step;
            if step.abs() < 1e-10 {
                let state = self.prepare(t, z)?;
                let p = ob * R * t / state.b;
                let volume = qc * state.b;
                return Ok(PseudoCritical {
                    temperature: t,
                    pressure: p,
                    volume,
                    compressibility: p * volume / (R * t),
                });
            }
        }
        Err(Error::NotConverged("fixed-composition critical point"))
    }
    pub fn new(kind: CubicKind, names: &[&str]) -> Result<Self> {
        Self::with_options(kind, names, CubicOptions::default())
    }

    /// Match the qcubic helper: PR, combined parameter priorities and 2 K bound.
    pub fn quantum(names: &[&str]) -> Result<Self> {
        let mut model = Self::with_options(
            CubicKind::PengRobinson,
            names,
            CubicOptions {
                parameter_reference: "QuantumCubic/tcPR/DEFAULT",
                ..Default::default()
            },
        )?;
        model.caloric_minimum = 2.;
        Ok(model)
    }

    pub fn temperature_bounds(&self) -> (f64, f64) {
        (self.caloric_minimum, self.caloric_maximum)
    }

    pub fn set_temperature_bounds(&mut self, minimum: f64, maximum: f64) -> Result<()> {
        positive(minimum, "minimum temperature")?;
        positive(maximum, "maximum temperature")?;
        if minimum >= maximum {
            return Err(Error::InvalidInput(
                "minimum temperature must be below maximum",
            ));
        }
        self.caloric_minimum = minimum;
        self.caloric_maximum = maximum;
        Ok(())
    }

    /// Public cubic initialization, including tcPR's forced TWU/translation.
    pub fn with_options(
        kind: CubicKind,
        names: &[&str],
        options: CubicOptions<'_>,
    ) -> Result<Self> {
        // eoslibinit stores the public constructor reference in len=100.
        let mut end = options.parameter_reference.len().min(100);
        while !options.parameter_reference.is_char_boundary(end) {
            end -= 1;
        }
        let parameter_reference = &options.parameter_reference[..end];
        let reference = parameter_reference.to_ascii_uppercase();
        let quantum = reference.contains("QUANTUMCUBIC");
        let tcpr = reference.contains("TCPR");
        let requested = if tcpr || quantum {
            crate::AlphaCorrelation::Twu
        } else {
            options.alpha.unwrap_or(if kind == CubicKind::RedlichKwong {
                crate::AlphaCorrelation::RedlichKwong
            } else {
                crate::AlphaCorrelation::Classic
            })
        };
        let mut model = Self::new_with_references(
            kind,
            names,
            requested,
            parameter_reference,
            parameter_reference,
            quantum,
        )?;
        if tcpr || quantum || options.volume_translation {
            model.enable_volume_translation(parameter_reference)?;
        }
        Ok(model)
    }

    /// Initialize alpha as in tpInitAlphaCorr. The reference selects alpha
    /// records only; other model parameters retain their default selection.
    pub fn with_alpha(
        kind: CubicKind,
        names: &[&str],
        requested: crate::AlphaCorrelation,
        alpha_reference: &str,
    ) -> Result<Self> {
        Self::new_with_references(kind, names, requested, alpha_reference, "DEFAULT", false)
    }

    fn new_with_references(
        kind: CubicKind,
        names: &[&str],
        requested: crate::AlphaCorrelation,
        alpha_reference: &str,
        parameter_reference: &str,
        quantum: bool,
    ) -> Result<Self> {
        if names.is_empty() {
            return Err(Error::InvalidInput("at least one component required"));
        }
        let mut components: Vec<_> = names
            .iter()
            .map(|name| component(name).copied())
            .collect::<Result<_>>()?;
        if quantum {
            for component in &mut components {
                crate::beta::quantum_component(component);
            }
        }
        let (omega_a, omega_b) = omega(kind);
        let ideal_gases = components
            .iter()
            .map(|c| crate::ideal::IdealGas::new(&crate::component::ideal_component(c)))
            .collect();
        let alpha = components
            .iter()
            .map(|c| crate::Alpha::initialize(kind, c, requested, alpha_reference))
            .collect::<Result<_>>()?;
        let mut model = Self {
            kind,
            ideal_gases,
            a0: Vec::with_capacity(names.len()),
            bi: Vec::with_capacity(names.len()),
            alpha,
            variable_beta: quantum,
            beta: components
                .iter()
                .map(|c| {
                    if quantum {
                        crate::Beta::quantum(c)
                    } else {
                        crate::Beta::default()
                    }
                })
                .collect(),
            kij: vec![0.0; names.len().pow(2)],
            lij: vec![0.0; names.len().pow(2)],
            nonzero_kij: Vec::new(),
            components,
            translation_parameters: None,
            translation_enabled: false,
            caloric_minimum: 80.,
            caloric_maximum: 999.,
        };
        for c in &model.components {
            positive(c.critical_temperature, "critical temperature")?;
            positive(c.critical_pressure, "critical pressure")?;
            model
                .a0
                .push(omega_a * (R * c.critical_temperature).powi(2) / c.critical_pressure);
            model
                .bi
                .push(omega_b * R * c.critical_temperature / c.critical_pressure);
        }
        for i in 0..names.len() {
            for j in i + 1..names.len() {
                let a = model.components[i].ident;
                let b = model.components[j].ident;
                let database_eos = if kind == CubicKind::RedlichKwong {
                    "SRK"
                } else {
                    kind.name()
                };
                let k = interaction_with_reference(database_eos, a, b, false, parameter_reference);
                let l = interaction_with_reference(database_eos, a, b, true, parameter_reference);
                model.kij[i * names.len() + j] = k;
                model.kij[j * names.len() + i] = k;
                model.lij[i * names.len() + j] = l;
                model.lij[j * names.len() + i] = l;
            }
        }
        model.rebuild_pairs();
        Ok(model)
    }

    pub fn kind(&self) -> CubicKind {
        self.kind
    }

    /// Select Peneloux parameters from the source database/Rackett fallback.
    /// This is the lower-level initializer and honors the reference argument.
    pub fn enable_volume_translation(&mut self, reference: &str) -> Result<()> {
        let parameters = self
            .components
            .iter()
            .map(|c| crate::VolumeTranslation::initialize(self.kind, c, reference))
            .collect::<Result<Vec<_>>>()?;
        self.translation_parameters = Some(parameters);
        self.translation_enabled = true;
        Ok(())
    }
    /// Match init_peneloux_volume_translation, whose source wrapper ignores ref.
    pub fn init_peneloux_volume_translation(&mut self, _parameter_reference: &str) -> Result<()> {
        self.enable_volume_translation("DEFAULT")
    }
    pub fn disable_volume_translation(&mut self) {
        self.translation_enabled = false;
    }
    pub fn volume_translation_enabled(&self) -> bool {
        self.translation_enabled
    }
    pub fn set_volume_translation_enabled(&mut self, enabled: bool) {
        self.translation_enabled = enabled;
    }
    pub fn volume_translation(&self, index: usize) -> Result<crate::VolumeTranslation> {
        if index >= self.components.len() {
            return Err(Error::InvalidInput("component index out of range"));
        }
        Ok(self
            .translation_parameters
            .as_ref()
            .map_or(crate::VolumeTranslation::default(), |values| values[index]))
    }
    /// Set ci parameters without changing the enabled flag, as in set_ci.
    pub fn set_volume_translation(
        &mut self,
        index: usize,
        parameters: crate::VolumeTranslation,
    ) -> Result<()> {
        if index >= self.components.len() {
            return Err(Error::InvalidInput("component index out of range"));
        }
        crate::VolumeTranslation::new(parameters.kind, parameters.coefficients)?;
        self.translation_parameters.get_or_insert_with(|| {
            vec![crate::VolumeTranslation::default(); self.components.len()]
        })[index] = parameters;
        Ok(())
    }
    pub fn alpha(&self, index: usize) -> Result<&crate::Alpha> {
        self.alpha
            .get(index)
            .ok_or(Error::InvalidInput("component index out of range"))
    }

    /// Set one component's effective alpha with explicit coefficients.
    pub fn set_alpha(
        &mut self,
        index: usize,
        correlation: crate::AlphaCorrelation,
        parameters: &[f64],
    ) -> Result<()> {
        let destination = self
            .alpha
            .get_mut(index)
            .ok_or(Error::InvalidInput("component index out of range"))?;
        *destination = crate::Alpha::new(correlation, parameters)?;
        Ok(())
    }
    /// Database metadata. Use ideal_gas for this model's active caloric parameters.
    pub fn components(&self) -> &[Component] {
        &self.components
    }

    pub fn beta(&self, index: usize) -> Result<crate::Beta> {
        self.beta
            .get(index)
            .copied()
            .ok_or(Error::InvalidInput("component index out of range"))
    }
    pub fn set_beta(&mut self, index: usize, beta: crate::Beta) -> Result<()> {
        crate::Beta::new(beta.correlation, beta.parameters)?;
        *self
            .beta
            .get_mut(index)
            .ok_or(Error::InvalidInput("component index out of range"))? = beta;
        self.variable_beta = self
            .beta
            .iter()
            .any(|b| b.correlation != crate::BetaCorrelation::Classic);
        Ok(())
    }

    pub(crate) fn pure_component(&self, index: usize) -> Result<Self> {
        let mut pure = Self::new(self.kind, &[self.components[index].ident])?;
        pure.components[0] = self.components[index];
        pure.a0[0] = self.a0[index];
        pure.bi[0] = self.bi[index];
        pure.alpha[0] = self.alpha[index];
        pure.beta[0] = self.beta[index];
        pure.variable_beta = pure.beta[0].correlation != crate::BetaCorrelation::Classic;
        pure.ideal_gases[0] = self.ideal_gases[index].clone();
        pure.set_volume_translation(0, self.volume_translation(index)?)?;
        pure.translation_enabled = self.translation_enabled;
        pure.caloric_minimum = self.caloric_minimum;
        pure.caloric_maximum = self.caloric_maximum;
        Ok(pure)
    }

    /// Active ideal parameters; component indices are zero based.
    pub fn ideal_gas(&self, index: usize) -> Result<&crate::ideal::IdealGas> {
        self.ideal_gases
            .get(index)
            .ok_or(Error::InvalidInput("component index out of range"))?
            .as_ref()
            .map_err(Clone::clone)
    }

    /// Configure this model's ideal parameters without changing the database
    /// or any other model instance. Existing prepared states borrow the model.
    pub fn ideal_gas_mut(&mut self, index: usize) -> Result<&mut crate::ideal::IdealGas> {
        self.ideal_gases
            .get_mut(index)
            .ok_or(Error::InvalidInput("component index out of range"))?
            .as_mut()
            .map_err(|e| e.clone())
    }

    pub fn kij(&self, i: usize, j: usize) -> Result<f64> {
        self.validate_pair(i, j)?;
        Ok(self.kij[i * self.components.len() + j])
    }

    pub fn set_kij(&mut self, i: usize, j: usize, value: f64) -> Result<()> {
        self.validate_pair(i, j)?;
        if !value.is_finite() || (i == j && value != 0.0) {
            return Err(Error::InvalidInput(
                "finite kij with zero diagonal required",
            ));
        }
        let n = self.components.len();
        self.kij[i * n + j] = value;
        self.kij[j * n + i] = value;
        self.rebuild_pairs();
        Ok(())
    }

    pub fn set_lij(&mut self, i: usize, j: usize, value: f64) -> Result<()> {
        self.validate_pair(i, j)?;
        if !value.is_finite() || (i == j && value != 0.0) {
            return Err(Error::InvalidInput(
                "finite lij with zero diagonal required",
            ));
        }
        let n = self.components.len();
        self.lij[i * n + j] = value;
        self.lij[j * n + i] = value;
        Ok(())
    }

    fn validate_pair(&self, i: usize, j: usize) -> Result<()> {
        if i < self.components.len() && j < self.components.len() {
            Ok(())
        } else {
            Err(Error::InvalidInput("component index out of range"))
        }
    }

    fn rebuild_pairs(&mut self) {
        self.nonzero_kij.clear();
        let n = self.components.len();
        for i in 0..n {
            for j in i + 1..n {
                let k = self.kij[i * n + j];
                if k != 0.0 {
                    self.nonzero_kij.push((i, j, k));
                }
            }
        }
    }

    pub fn prepare(&self, temperature: f64, composition: &[f64]) -> Result<PreparedState<'_>> {
        positive(temperature, "temperature")?;
        if composition.len() != self.components.len()
            || composition.iter().any(|x| !x.is_finite() || *x < 0.0)
            || (composition.iter().sum::<f64>() - 1.0).abs() > 1e-10
        {
            return Err(Error::InvalidInput(
                "normalized composition matching the model required",
            ));
        }
        let n = self.components.len();
        let mut factors = Vec::with_capacity(n);
        let (mut sum, mut sum_t, mut sum_tt) = (0.0, 0.0, 0.0);
        for (i, &x) in composition.iter().enumerate() {
            let tc = self.components[i].critical_temperature;
            let crate::AlphaValue {
                value: alpha,
                dt: alpha_t,
                dtt: alpha_tt,
            } = self.alpha[i].evaluate(temperature, tc)?;
            let f = (self.a0[i] * alpha.max(0.0)).sqrt();
            let ft = if f > 0.0 {
                0.5 * self.a0[i] * alpha_t / f
            } else {
                0.0
            };
            let ftt = if f > 0.0 {
                (0.5 * self.a0[i] * alpha_tt - ft * ft) / f
            } else {
                0.0
            };
            factors.push((f, ft, ftt));
            sum += x * f;
            sum_t += x * ft;
            sum_tt += x * ftt;
        }
        let mut a = sum * sum;
        let mut at = 2.0 * sum * sum_t;
        let mut att = 2.0 * (sum_t * sum_t + sum * sum_tt);
        let mut a_partial: Vec<_> = factors.iter().map(|f| f.0 * sum).collect();
        for &(i, j, k) in &self.nonzero_kij {
            let (fi, fit, fitt) = factors[i];
            let (fj, fjt, fjtt) = factors[j];
            let pair = k * fi * fj;
            let xx = 2.0 * composition[i] * composition[j];
            a -= xx * pair;
            at -= xx * k * (fit * fj + fi * fjt);
            att -= xx * k * (fitt * fj + 2.0 * fit * fjt + fi * fjtt);
            a_partial[i] -= composition[j] * pair;
            a_partial[j] -= composition[i] * pair;
        }
        let mut b = 0.0;
        let mut beta = if self.variable_beta {
            let components = self
                .beta
                .iter()
                .zip(&self.components)
                .zip(&self.bi)
                .map(|((beta, component), &base)| {
                    let value = beta.evaluate(temperature, component.critical_temperature)?;
                    Ok(crate::BetaValue {
                        value: base * value.value,
                        dt: base * value.dt,
                        dtt: base * value.dtt,
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            Some(crate::beta::PreparedBeta {
                partial_t: components.iter().map(|c| c.dt).collect(),
                components,
                dt: 0.,
                dtt: 0.,
            })
        } else {
            None
        };
        let mut b_partial: Vec<_> = match &beta {
            Some(beta) => beta.components.iter().map(|b| b.value).collect(),
            None => self.bi.clone(),
        };
        if self.lij.iter().all(|l| *l == 0.0) {
            b = composition
                .iter()
                .zip(&b_partial)
                .map(|(x, bi)| x * bi)
                .sum();
            if let Some(beta) = &mut beta {
                beta.dt = composition
                    .iter()
                    .zip(&beta.components)
                    .map(|(x, b)| x * b.dt)
                    .sum();
                beta.dtt = composition
                    .iter()
                    .zip(&beta.components)
                    .map(|(x, b)| x * b.dtt)
                    .sum();
            }
        } else {
            for i in 0..n {
                b_partial[i] = 0.0;
                let bii = beta.as_ref().map_or(self.bi[i], |b| b.components[i].value);
                let mut bt_partial = 0.;
                let mut btt_partial = 0.;
                for (j, &xj) in composition.iter().enumerate() {
                    let bjj = beta.as_ref().map_or(self.bi[j], |b| b.components[j].value);
                    b_partial[i] += xj * (1.0 - self.lij[i * n + j]) * (bii + bjj) / 2.0;
                    if let Some(beta) = &beta {
                        let weight = xj * (1. - self.lij[i * n + j]) / 2.;
                        bt_partial += weight * (beta.components[i].dt + beta.components[j].dt);
                        btt_partial += weight * (beta.components[i].dtt + beta.components[j].dtt);
                    }
                }
                b += composition[i] * b_partial[i];
                if let Some(beta) = &mut beta {
                    beta.partial_t[i] = 2. * bt_partial;
                    beta.dt += composition[i] * bt_partial;
                    beta.dtt += composition[i] * btt_partial;
                }
            }
            for bi in &mut b_partial {
                *bi = 2.0 * *bi - b;
            }
            if let Some(beta) = &mut beta {
                for bit in &mut beta.partial_t {
                    *bit -= beta.dt;
                }
            }
        }
        positive(b, "mixture covolume")?;
        if !a.is_finite() || !at.is_finite() || !att.is_finite() {
            return Err(Error::InvalidInput("temperature outside numerical range"));
        }
        let translation = if self.translation_enabled {
            // Borrow active coefficients in the hot loop; the public getter
            // returns a full parameter record, including provenance metadata.
            let components = match &self.translation_parameters {
                Some(parameters) => parameters
                    .iter()
                    .map(|p| p.evaluate(temperature))
                    .collect::<Result<Vec<_>>>()?,
                None => vec![
                    crate::TranslationValue {
                        value: 0.,
                        dt: 0.,
                        dtt: 0.
                    };
                    n
                ],
            };
            let mut prepared = crate::translation::PreparedTranslation {
                c: 0.,
                ct: 0.,
                ctt: 0.,
                components,
            };
            for (&x, value) in composition.iter().zip(&prepared.components) {
                prepared.c += x * value.value;
                prepared.ct += x * value.dt;
                prepared.ctt += x * value.dtt;
            }
            Some(prepared)
        } else {
            None
        };
        Ok(PreparedState {
            model: self,
            temperature,
            composition: composition.to_vec(),
            a,
            at,
            att,
            b,
            a_partial,
            b_partial,
            factors,
            beta,
            translation,
        })
    }
}

impl PreparedState<'_> {
    /// Analytic derivatives of log fugacity coefficients at constant T/P.
    /// `dt` is 1/K, `dp` is 1/Pa, and `dn` is a row-major n*n matrix
    /// of derivatives with respect to mole numbers at total amount 1 mol.
    /// All output buffers are caller-owned and must have the exact dimensions.
    pub fn fugacity_derivatives_into(
        &self,
        volume: f64,
        dt: &mut [f64],
        dp: &mut [f64],
        dn: &mut [f64],
    ) -> Result<()> {
        let n = self.composition.len();
        if dt.len() != n || dp.len() != n || dn.len() != n * n {
            return Err(Error::InvalidInput("fugacity derivative buffer dimensions"));
        }
        self.volume_check(volume)?;
        let volume = self.eos_volume(volume);
        let p = self.eos_pressure(volume)?;
        positive(p.value, "pressure for fugacity derivatives")?;
        if p.dv == 0.0 {
            return Err(Error::InvalidInput("singular volume derivative"));
        }
        let t = self.temperature;
        let rt = R * t;
        let free = volume - self.b;
        let (d1, d2) = self.model.kind.deltas();
        let denom = (volume + d1 * self.b) * (volume + d2 * self.b);
        let denom_b = (d1 + d2) * volume + 2.0 * d1 * d2 * self.b;
        let (integral, ib, ibb) = self.attraction_integral_derivatives(volume);
        let pressure_n = |i: usize| {
            rt / free + rt * self.b_partial[i] / free.powi(2) - 2.0 * self.a_partial[i] / denom
                + self.a * denom_b * self.b_partial[i] / denom.powi(2)
        };
        for i in 0..n {
            let mut ai_t = 0.0;
            for j in 0..n {
                ai_t += self.composition[j]
                    * (1.0 - self.model.kij[i * n + j])
                    * (self.factors[i].1 * self.factors[j].0
                        + self.factors[i].0 * self.factors[j].1);
            }
            let bi = self.b_partial[i];
            let di = 2.0 * self.a_partial[i];
            let pi = pressure_n(i);
            let mut fnt =
                -(2.0 * ai_t - di / t) * integral / rt - (self.at - self.a / t) * ib * bi / rt;
            if let Some(beta) = &self.beta {
                let bt = beta.dt;
                let bit = beta.partial_t[i];
                fnt += (bt + bit) / free + bi * bt / (free * free)
                    - (di * ib * bt + self.a * ibb * bt * bi + self.a * ib * bit) / rt;
            }
            dt[i] = fnt + pi * p.dt / (rt * p.dv) + 1.0 / t;
            dp[i] = -pi / (rt * p.dv) - 1.0 / p.value;
            for j in 0..n {
                let bj = self.b_partial[j];
                let dj = 2.0 * self.a_partial[j];
                let bii = self
                    .beta
                    .as_ref()
                    .map_or(self.model.bi[i], |b| b.components[i].value);
                let bjj = self
                    .beta
                    .as_ref()
                    .map_or(self.model.bi[j], |b| b.components[j].value);
                let bij = (1.0 - self.model.lij[i * n + j]) * (bii + bjj) / 2.0;
                let bnn = 2.0 * bij - bi - bj;
                let dnn =
                    2.0 * (1.0 - self.model.kij[i * n + j]) * self.factors[i].0 * self.factors[j].0;
                let fnn = (bi + bj + bnn) / free + bi * bj / free.powi(2)
                    - (dnn * integral
                        + (di * bj + dj * bi) * ib
                        + self.a * ibb * bi * bj
                        + self.a * ib * bnn)
                        / rt;
                dn[i * n + j] = fnn + 1.0 + pi * pressure_n(j) / (rt * p.dv);
            }
        }
        if let Some(shift) = &self.translation {
            for (i, component) in shift.components.iter().enumerate() {
                dt[i] += p.value / R * (component.value / (t * t) - component.dt / t);
                dp[i] -= component.value / rt;
            }
        }
        Ok(())
    }

    fn attraction_integral_derivatives(&self, volume: f64) -> (f64, f64, f64) {
        let (d1, d2) = self.model.kind.deltas();
        if d1 == d2 {
            return (1.0 / volume, 0.0, 0.0);
        }
        let ratio = self.b / volume;
        if ratio < 0.01 {
            let (mut i, mut ib, mut ibb) = (0.0, 0.0, 0.0);
            for k in 0_i32..12 {
                let coefficient = (-1.0_f64).powi(k) * (d1.powi(k + 1) - d2.powi(k + 1))
                    / ((d1 - d2) * (k + 1) as f64);
                i += coefficient * ratio.powi(k) / volume;
                if k >= 1 {
                    ib += coefficient * k as f64 * ratio.powi(k - 1) / volume.powi(2);
                }
                if k >= 2 {
                    ibb += coefficient * (k * (k - 1)) as f64 * ratio.powi(k - 2) / volume.powi(3);
                }
            }
            (i, ib, ibb)
        } else {
            let integral = self.attraction_integral(volume);
            let log_b = d1 / (volume + d1 * self.b) - d2 / (volume + d2 * self.b);
            let log_bb =
                -(d1 / (volume + d1 * self.b)).powi(2) + (d2 / (volume + d2 * self.b)).powi(2);
            let ib = (log_b / (d1 - d2) - integral) / self.b;
            let ibb = (log_bb / (d1 - d2) - 2.0 * ib) / self.b;
            (integral, ib, ibb)
        }
    }
    /// Ideal mixture enthalpy, entropy (including entropy of mixing) and Cp.
    pub fn ideal(&self, pressure: f64) -> Result<crate::ideal::IdealProperties> {
        positive(pressure, "pressure")?;
        let mut result = crate::ideal::IdealProperties {
            cp: 0.0,
            enthalpy_dt: 0.0,
            enthalpy: 0.0,
            entropy: 0.0,
        };
        for (&x, gas) in self.composition.iter().zip(&self.model.ideal_gases) {
            if x == 0.0 {
                continue;
            }
            let pure = gas
                .as_ref()
                .map_err(Clone::clone)?
                .properties(self.temperature, pressure)?;
            result.cp += x * pure.cp;
            result.enthalpy_dt += x * pure.enthalpy_dt;
            result.enthalpy += x * pure.enthalpy;
            result.entropy += x * (pure.entropy - R * x.ln());
        }
        Ok(result)
    }

    pub fn enthalpy(&self, volume: f64) -> Result<f64> {
        Ok(self.ideal(self.pressure(volume)?.value)?.enthalpy + self.residual(volume)?.enthalpy)
    }

    pub fn entropy(&self, volume: f64) -> Result<f64> {
        Ok(self.ideal(self.pressure(volume)?.value)?.entropy + self.residual(volume)?.entropy)
    }

    /// Molar U (J/mol), S (J/(mol K)), Cv and pressure at the supplied molar
    /// volume. Heat capacities differentiate actual native caloric values.
    pub fn caloric_tv(&self, volume: f64) -> Result<TvProperties> {
        self.volume_check(volume)?;
        let v = self.eos_volume(volume);
        let t = self.temperature;
        let ideal = self.ideal(1e5)?;
        let pressure = self.eos_pressure(v)?;
        let integral = self.attraction_integral(v);
        let (_, ab) = self.eos_thermal_terms(v, integral);
        let beta_t = self.beta.as_ref().map_or(0., |b| b.dt);
        let translation_t = self.translation.as_ref().map_or(0., |c| c.ct);
        let u = ideal.enthalpy - R * t + (t * self.at - self.a) * integral - t * ab * beta_t
            + pressure.value * t * translation_t;
        let s = ideal.entropy + R * ((v - self.b).ln() - (R * t / 1e5).ln()) + self.at * integral
            - ab * beta_t
            + pressure.value * translation_t;
        let cv = ideal.enthalpy_dt - R + self.residual_cv_tv(volume)?;
        if [pressure.value, u, s, cv]
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(Error::InvalidInput("non-finite TV caloric properties"));
        }
        Ok(TvProperties {
            pressure: pressure.value,
            internal_energy: u,
            entropy: s,
            heat_capacity_v: cv,
        })
    }

    /// Single-phase heat capacities (Cp, Cv) [J/(mol K)].
    pub fn heat_capacities(&self, volume: f64) -> Result<(f64, f64)> {
        let ideal_cp = self.ideal(self.pressure(volume)?.value)?.enthalpy_dt;
        let residual = self.residual(volume)?;
        Ok((ideal_cp + residual.cp, ideal_cp - R + residual.cv))
    }

    /// Compatibility with this fork's `enthalpy(..., dhdt=True)`.
    /// Preserves the vdW F_TT expression and type-10 ideal-H2 Cp discrepancy.
    /// Use `heat_capacities` for a derivative consistent with enthalpy values.
    /// See docs/rust/baseline-defects.md for the exact discrepancy.
    pub fn baseline_enthalpy_derivative_tp(&self, volume: f64) -> Result<f64> {
        let ideal = self.ideal(self.pressure(volume)?.value)?;
        let cp = ideal.cp + self.residual(volume)?.cp;
        if self.model.kind == CubicKind::VanDerWaals {
            let volume = self.eos_volume(volume);
            Ok(
                cp + 2.0 * self.a * self.temperature / (1000.0 * volume * volume)
                    - 2.0 * self.a / (self.temperature * volume),
            )
        } else {
            Ok(cp)
        }
    }

    /// Single-phase equilibrium sound speed [m/s].
    pub fn speed_of_sound(&self, volume: f64) -> Result<f64> {
        let cv = self.ideal(1e5)?.enthalpy_dt - R + self.residual_cv_tv(volume)?;
        let mass: f64 = self
            .composition
            .iter()
            .zip(&self.model.components)
            .map(|(x, c)| x * c.molar_mass * 1e-3)
            .sum();
        let pressure = self.pressure(volume)?;
        let squared =
            -volume * volume * (pressure.dv - self.temperature * pressure.dt.powi(2) / cv) / mass;
        positive(squared, "sound speed squared")?;
        Ok(squared.sqrt())
    }

    /// Source TV sound speed, retaining reported H2/vdW Cv and zero on an
    /// invalid isentropic stiffness. Uses molar volume, independent of total n.
    pub fn baseline_speed_of_sound(&self, volume: f64) -> Result<f64> {
        let mut cv = self.ideal(1e5)?.cp - R + self.residual_cv_tv(volume)?;
        if self.model.kind == CubicKind::VanDerWaals {
            let eos_volume = self.eos_volume(volume);
            cv += 2. * self.a * self.temperature / (1000. * eos_volume * eos_volume)
                - 2. * self.a / (self.temperature * eos_volume);
        }
        let pressure = self.pressure(volume)?;
        let derivative = pressure.dv - self.temperature * pressure.dt * pressure.dt / cv;
        if derivative < 0. {
            let mass: f64 = self
                .composition
                .iter()
                .zip(&self.model.components)
                .map(|(x, c)| x * c.molar_mass * 1e-3)
                .sum();
            Ok((-volume * volume * derivative / mass).sqrt())
        } else {
            Ok(0.)
        }
    }
    pub fn temperature(&self) -> f64 {
        self.temperature
    }
    pub fn composition(&self) -> &[f64] {
        &self.composition
    }
    pub fn covolume(&self) -> f64 {
        self.b - self.translation.as_ref().map_or(0., |c| c.c)
    }

    /// Untranslated repulsive covolume (m3/mol).
    pub fn eos_covolume(&self) -> f64 {
        self.b
    }
    pub fn volume_translation(&self) -> crate::TranslationValue {
        self.translation.as_ref().map_or(
            crate::TranslationValue {
                value: 0.,
                dt: 0.,
                dtt: 0.,
            },
            |c| crate::TranslationValue {
                value: c.c,
                dt: c.ct,
                dtt: c.ctt,
            },
        )
    }
    #[inline]
    fn eos_volume(&self, volume: f64) -> f64 {
        volume + self.translation.as_ref().map_or(0., |c| c.c)
    }

    #[inline]
    fn volume_check(&self, volume: f64) -> Result<()> {
        if volume.is_finite() && volume > 0. && self.eos_volume(volume) > self.b {
            Ok(())
        } else {
            Err(Error::InvalidInput("molar volume must exceed covolume"))
        }
    }

    #[inline]
    pub fn pressure(&self, volume: f64) -> Result<Pressure> {
        self.volume_check(volume)?;
        let mut pressure = self.eos_pressure(self.eos_volume(volume))?;
        if let Some(shift) = &self.translation {
            pressure.dt += pressure.dv * shift.ct;
        }
        Ok(pressure)
    }

    #[inline]
    fn eos_pressure(&self, volume: f64) -> Result<Pressure> {
        if !volume.is_finite() || volume <= self.b {
            return Err(Error::InvalidInput("EOS volume must exceed covolume"));
        }
        let (d1, d2) = self.model.kind.deltas();
        let denom = (volume + d1 * self.b) * (volume + d2 * self.b);
        let free = volume - self.b;
        let beta_dt = self.beta.as_ref().map_or(0., |beta| {
            let denom_b = (d1 + d2) * volume + 2. * d1 * d2 * self.b;
            (R * self.temperature / (free * free) + self.a * denom_b / (denom * denom)) * beta.dt
        });
        Ok(Pressure {
            value: R * self.temperature / free - self.a / denom,
            dt: R / free - self.at / denom + beta_dt,
            dv: -R * self.temperature / (free * free)
                + self.a * (2.0 * volume + (d1 + d2) * self.b) / (denom * denom),
        })
    }

    fn attraction_integral(&self, volume: f64) -> f64 {
        let (d1, d2) = self.model.kind.deltas();
        if d1 == d2 {
            1.0 / volume
        } else {
            ((d1 - d2) * self.b / (volume + d2 * self.b)).ln_1p() / ((d1 - d2) * self.b)
        }
    }

    // Untranslated Cv and A_b for temperature-dependent covolume. Shared by
    // TP caloric properties and the TV sound path, including negative P.
    fn eos_thermal_terms(&self, volume: f64, integral: f64) -> (f64, f64) {
        let t = self.temperature;
        let mut cv = t * self.att * integral;
        let mut ab = 0.;
        if let Some(beta) = &self.beta {
            let rt = R * t;
            let free = volume - self.b;
            let (_, ib, ibb) = self.attraction_integral_derivatives(volume);
            ab = rt / free - self.a * ib;
            cv -= t
                * (2. * (R / free - self.at * ib) * beta.dt
                    + (rt / (free * free) - self.a * ibb) * beta.dt * beta.dt
                    + ab * beta.dtt);
        }
        (cv, ab)
    }

    fn residual_cv_tv(&self, actual_volume: f64) -> Result<f64> {
        self.volume_check(actual_volume)?;
        let volume = self.eos_volume(actual_volume);
        let (mut cv, _) = self.eos_thermal_terms(volume, self.attraction_integral(volume));
        if let Some(shift) = &self.translation {
            let pressure = self.eos_pressure(volume)?;
            let t = self.temperature;
            cv += 2. * t * pressure.dt * shift.ct
                + t * pressure.dv * shift.ct * shift.ct
                + pressure.value * t * shift.ctt;
        }
        Ok(cv)
    }

    pub fn residual(&self, volume: f64) -> Result<Residual> {
        self.volume_check(volume)?;
        let actual_volume = volume;
        let volume = self.eos_volume(volume);
        let pressure = self.eos_pressure(volume)?;
        positive(pressure.value, "pressure for T-P residual properties")?;
        let z = pressure.value * volume / (R * self.temperature);
        let integral = self.attraction_integral(volume);
        let (cv, beta_ab) = self.eos_thermal_terms(volume, integral);
        let mut residual = Residual {
            compressibility: z,
            helmholtz_tv: -R * self.temperature * (-self.b / volume).ln_1p() - self.a * integral,
            enthalpy: R * self.temperature * (z - 1.0)
                + (self.temperature * self.at - self.a) * integral,
            entropy: R * (z * (1.0 - self.b / volume)).ln() + self.at * integral,
            cv,
            cp: cv - self.temperature * pressure.dt.powi(2) / pressure.dv - R,
        };
        if let Some(beta) = &self.beta {
            let t = self.temperature;
            residual.enthalpy -= t * beta_ab * beta.dt;
            residual.entropy -= beta_ab * beta.dt;
        }
        if let Some(shift) = &self.translation {
            let t = self.temperature;
            residual.compressibility = pressure.value * actual_volume / (R * t);
            residual.helmholtz_tv += R * t * (-shift.c / volume).ln_1p();
            residual.enthalpy += pressure.value * (t * shift.ct - shift.c);
            residual.entropy += pressure.value * shift.ct;
            residual.cv += 2. * t * pressure.dt * shift.ct
                + t * pressure.dv * shift.ct * shift.ct
                + pressure.value * t * shift.ctt;
            residual.cp += pressure.value * t * shift.ctt;
        }
        Ok(residual)
    }

    /// Write log fugacity coefficients into a caller-owned buffer.
    pub fn ln_fugacity_into(&self, volume: f64, output: &mut [f64]) -> Result<()> {
        self.ln_fugacity_tp_into(self.pressure(volume)?.value, volume, output)
    }

    /// Fugacity at an imposed TP pressure and its solved molar volume.
    /// Using the specified pressure avoids reconstructing it by subtracting
    /// large repulsive/attractive terms in dense, low-pressure liquids.
    pub fn ln_fugacity_tp_into(&self, p: f64, volume: f64, output: &mut [f64]) -> Result<()> {
        if output.len() != self.composition.len() {
            return Err(Error::InvalidInput("fugacity buffer length"));
        }
        self.volume_check(volume)?;
        positive(p, "pressure for fugacity")?;
        let volume = self.eos_volume(volume);
        let rt = R * self.temperature;
        let z = p * volume / rt;
        let repulsion = (z * (1.0 - self.b / volume)).ln();
        let integral = self.attraction_integral(volume) / rt;
        for (i, phi) in output.iter_mut().enumerate() {
            let br = self.b_partial[i] / self.b;
            *phi = br * (z - 1.0) - repulsion - (2.0 * self.a_partial[i] - self.a * br) * integral;
            if let Some(shift) = &self.translation {
                *phi -= p * shift.components[i].value / rt;
            }
        }
        Ok(())
    }

    /// Physical roots at positive pressure, ascending by volume.
    pub fn volumes(&self, pressure: f64) -> Result<([f64; 3], usize)> {
        let (mut roots, count) = self.eos_volumes(pressure)?;
        if let Some(shift) = &self.translation {
            let mut kept = 0;
            for i in 0..count {
                let physical = roots[i] - shift.c;
                if physical > 0. {
                    roots[kept] = physical;
                    kept += 1;
                }
            }
            if kept == 0 {
                return Err(Error::NoPhysicalRoot);
            }
            return Ok((roots, kept));
        }
        Ok((roots, count))
    }

    fn eos_volumes(&self, pressure: f64) -> Result<([f64; 3], usize)> {
        positive(pressure, "pressure")?;
        let rt = R * self.temperature;
        let a = self.a * pressure / (rt * rt);
        let b = self.b * pressure / rt;
        let (d1, d2) = self.model.kind.deltas();
        let u = d1 + d2;
        let w = d1 * d2;
        let (zroots, nz) = numerics::cubic_roots(
            (u - 1.0) * b - 1.0,
            a + w * b * b - u * b * (b + 1.0),
            -a * b - w * b * b * (b + 1.0),
        );
        let mut roots = [f64::NAN; 3];
        let mut count = 0;
        for &z in &zroots[..nz] {
            let mut v = z * rt / pressure;
            if !v.is_finite() || v <= self.b {
                continue;
            }
            for _ in 0..4 {
                let p = self.eos_pressure(v)?;
                if p.dv == 0.0 {
                    break;
                }
                let step = (p.value - pressure) / p.dv;
                let candidate = v - step;
                if candidate > self.b && step.abs() < 0.1 * (v - self.b) {
                    v = candidate;
                } else {
                    break;
                }
            }
            let p = self.eos_pressure(v)?;
            if p.dv > 0.0 || (p.value - pressure).abs() > 1e-7 * pressure + 1e-3 {
                continue;
            }
            if count > 0 && (v - roots[count - 1]).abs() < 1e-12 * v {
                continue;
            }
            roots[count] = v;
            count += 1;
        }
        if count == 0 {
            Err(Error::NoPhysicalRoot)
        } else {
            Ok((roots, count))
        }
    }

    pub fn volume(&self, pressure: f64, phase: Phase) -> Result<f64> {
        let (roots, count) = self.volumes(pressure)?;
        match phase {
            Phase::Liquid => Ok(roots[0]),
            Phase::Vapor => Ok(roots[count - 1]),
            Phase::MinimumGibbs => {
                let mut best = roots[0];
                let mut gibbs = f64::INFINITY;
                for &v in &roots[..count] {
                    let residual = self.residual(v)?;
                    let g = residual.enthalpy - self.temperature * residual.entropy;
                    if g < gibbs {
                        gibbs = g;
                        best = v;
                    }
                }
                Ok(best)
            }
        }
    }
}

// Preserve the source iteration, including evaluating omega from the final
// iteration's n/n1/n2 (cbCalcOmegaZc), rather than rounded textbook constants.
fn omega(kind: CubicKind) -> (f64, f64) {
    let (a, b, _) = omega_with_critical_ratio(kind);
    (a, b)
}

fn omega_with_critical_ratio(kind: CubicKind) -> (f64, f64, f64) {
    let (d1, d2) = kind.deltas();
    let (m1, m2) = (-d2, -d1);
    let mut zcomb = 5.0;
    let (mut n, mut n1, mut n2) = (0.0, 0.0, 0.0);
    for _ in 0..9 {
        n = zcomb - 1.0;
        n1 = zcomb - m1;
        n2 = zcomb - m2;
        let f = n * n1 * n1 + n * n1 * n2 + n * n2 * n2 - n1 * n1 * n2 - n1 * n2 * n2;
        let fd = 3.0 * (n * n1 + n * n2 - n1 * n2);
        let fdd = 6.0 * n;
        let arg = 1.0 - 2.0 * f * fdd / (fd * fd);
        let arg = if arg > 0.0 { arg.sqrt() } else { 1.0 };
        let dz = -2.0 * f / (fd * (1.0 + arg));
        zcomb += dz;
        if dz.abs() <= 1e-10 {
            break;
        }
    }
    let arg = n1 * n2 / (n * n1 + n * n2);
    let b = (1.0 - arg) / (zcomb - 1.0);
    (b * (1.0 / n - b) * n1 * n2, b, zcomb)
}
