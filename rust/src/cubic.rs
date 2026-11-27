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
    b_exponent: f64,
    nonzero_kij: Vec<(usize, usize, f64)>,
    translation_parameters: Option<Vec<crate::VolumeTranslation>>,
    translation_enabled: bool,
    pub(crate) caloric_minimum: f64,
    pub(crate) caloric_maximum: f64,
    pub mixing_rule: crate::excess_gibbs::MixingRule,
    pub excess_gibbs: Option<crate::excess_gibbs::ExcessGibbsMixture>,
    /// The binary parameters of the Wong-Sandler mixing rule.
    pub wong_sandler: Option<crate::wong_sandler::WongSandlerData>,
}

#[derive(Debug, Clone, Copy)]
pub struct CubicOptions<'a> {
    pub alpha: Option<crate::AlphaCorrelation>,
    pub parameter_reference: &'a str,
    pub volume_translation: bool,
    pub mixing_rule: crate::excess_gibbs::MixingRule,
    /// The binary parameters of the Wong-Sandler mixing rule.
    pub wong_sandler: Option<&'a crate::wong_sandler::WongSandlerData>,
}
impl Default for CubicOptions<'_> {
    fn default() -> Self {
        Self {
            alpha: None,
            parameter_reference: "DEFAULT",
            volume_translation: false,
            mixing_rule: crate::excess_gibbs::MixingRule::VdW,
            wong_sandler: None,
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

/// Complete 10-block reduced residual Helmholtz energy `F = A^res / (R T)` and
/// all first and second `(T, V, n)` derivatives for a cubic equation of state.
#[derive(Debug, Clone)]
pub struct CubicHelmholtz {
    /// Dimensionless reduced residual Helmholtz energy `F = A^res / (R T)`.
    pub f: f64,
    /// `dF / dT` [1/K].
    pub f_t: f64,
    /// `dF / dV` [1/m^3].
    pub f_v: f64,
    /// `dF / dn_i` [1/mol].
    pub f_n: Vec<f64>,
    /// `d^2F / dT^2` [1/K^2].
    pub f_tt: f64,
    /// `d^2F / (dT dV)` [1/(K m^3)].
    pub f_tv: f64,
    /// `d^2F / (dT dn_i)` [1/(K mol)].
    pub f_tn: Vec<f64>,
    /// `d^2F / dV^2` [1/m^6].
    pub f_vv: f64,
    /// `d^2F / (dV dn_i)` [1/(m^3 mol)].
    pub f_vn: Vec<f64>,
    /// `d^2F / (dn_i dn_j)` [1/mol^2], row-major `nc * nc`.
    pub f_nn: Vec<f64>,
}

/// Cached temperature and composition parameters. Volume/pressure scans
/// perform no allocations, database lookups or component-pair mixing.
#[derive(Debug, Clone)]
pub struct PreparedState<'a> {
    model: &'a Cubic,
    temperature: f64,
    /// Component buffers: the composition, the partial attraction parameters,
    /// the partial covolumes and the three square-root factors share one
    /// allocation.
    buffer: Vec<f64>,
    a: f64,
    at: f64,
    att: f64,
    b: f64,
    ai_t: Option<Vec<f64>>,
    /// GE attraction Hessian, replacing the classical pair Hessian.
    ge_a_nn: Option<Vec<f64>>,
    beta: Option<crate::beta::PreparedBeta>,
    translation: Option<crate::translation::PreparedTranslation>,
    /// The second mole derivatives of the mixture coefficients of the
    /// Wong-Sandler mixing rule.
    wong_sandler: Option<crate::wong_sandler::PreparedWongSandler>,
}

impl PreparedState<'_> {
    /// The fork's gflag=1 cubic compressibility iteration, including its raw
    /// last iterate at a nearly multiple root. This compatibility result is
    /// not filtered by mechanical stability or pressure closure; use `volume`
    /// for physical root selection. Minimum-Gibbs selection is not this path.
    pub fn baseline_volume(&self, pressure: f64, phase: Phase) -> Result<f64> {
        self.baseline_volume_with_extrema(pressure, phase, false)
    }
    pub(crate) fn baseline_volume_with_extrema(
        &self,
        pressure: f64,
        phase: Phase,
        meta: bool,
    ) -> Result<f64> {
        positive(pressure, "cubic pressure")?;
        if phase == Phase::MinimumGibbs {
            return Err(Error::InvalidInput("baseline cubic phase"));
        }
        let t = self.temperature;
        let a = self.a * pressure / (t * R * t * R);
        let b = self.b * pressure / (t * R);
        let (d1, d2) = self.model.kind.deltas();
        let (m1, m2) = (-d2, -d1);
        let pp = b * (-m1 - m2 - 1.0) - 1.0;
        let qq = b * ((m1 + m2) * (b + 1.0) + b * m1 * m2) + a;
        let rr = -(b * b * b + b * b) * (m1 * m2) - a * b;
        let (_, ob, qc) = omega_with_critical_ratio(self.model.kind);
        let zc = ob * qc;
        let f3 = rr + zc * (qq + zc * (pp + zc));
        let p2 = 2.0 * pp;
        let discriminant = p2 * p2 - 12.0 * qq;
        let mut liquid = f3 > 0.0;
        let mut extremum = None;
        if discriminant > 0.0 {
            let z1 = (-2.0 * pp - discriminant.sqrt()) / 6.0;
            let z2 = (-2.0 * pp + discriminant.sqrt()) / 6.0;
            let f1 = rr + z1 * (qq + z1 * (pp + z1));
            let f2 = rr + z2 * (qq + z2 * (pp + z2));
            if f1 * f2 <= 0.0 {
                liquid = phase == Phase::Liquid;
            } else if meta {
                if liquid && phase == Phase::Vapor {
                    extremum = Some(z2.max(b * 1.001));
                }
                if !liquid && phase == Phase::Liquid {
                    extremum = Some(z1.max(b * 1.001));
                }
            }
        }
        let solve = |mut z: f64| {
            let eps = 10.0 * f64::EPSILON;
            let mut dz = 1e10 * eps;
            let mut f = rr + z * (qq + z * (pp + z));
            let mut iteration = 0;
            while iteration < 50 && (dz / z).abs() >= eps {
                let df = qq + z * (2.0 * pp + 3.0 * z);
                let ddf = 2.0 * pp + 6.0 * z;
                let arg = 1.0 - 2.0 * f * ddf / (df * df);
                let arg = if arg < 0.0 { 1.0 } else { arg.sqrt() };
                dz = -(f / df) * 2.0 / (1.0 + arg);
                let mut next = z + dz;
                if (dz / z).abs() < eps {
                    z = next;
                } else {
                    let mut fn_ = rr + next * (qq + next * (pp + next));
                    let mut relax = 0;
                    while fn_.abs() >= f.abs() && relax < 15 {
                        relax += 1;
                        dz *= 0.5;
                        next = z + dz;
                        fn_ = rr + next * (qq + next * (pp + next));
                    }
                    z = next;
                    f = fn_;
                    iteration += 1;
                }
            }
            z
        };
        let mut z =
            extremum.unwrap_or_else(|| solve(if liquid { b } else { 2.0_f64.max(b * 1.001) }));
        if extremum.is_none() && liquid && z < b {
            z = solve(2.0_f64.max(b * 1.001));
        }
        let volume = z * R * t / pressure - self.volume_translation().value;
        if !volume.is_finite() || volume <= 0.0 {
            return Err(Error::NotConverged("baseline cubic volume"));
        }
        Ok(volume)
    }

    // Analytic third composition derivative for classic quadratic mixing.
    // Only a scalar direction is differentiated; no n*n*n tensor is built.
    pub(crate) fn critical_third_directional(&self, volume: f64, direction: &[f64]) -> Option<f64> {
        use crate::hyperdual::HyperDual;
        if self.model.excess_gibbs.is_some() {
            return None;
        }
        let n = self.composition_slice().len();
        let mut s = HyperDual::with_order(0.0, 3);
        s.f1 = 1.0;
        s.f2 = 1.0;
        s.f3 = 1.0;
        let real = |v| HyperDual::with_order(v, 3);
        let amount = real(1.0) + s * direction.iter().sum::<f64>();
        let mut d1 = 0.0;
        let mut d2 = 0.0;
        let mut b1 = 0.0;
        let mut b2 = 0.0;
        let mut shift = 0.0;
        let mut ideal = real(0.0);
        for i in 0..n {
            let zi = self.composition_slice()[i];
            if zi > 0.0 {
                let ni = real(zi) + s * direction[i];
                ideal = ideal + ni * ni.ln();
            }
            d1 += 2.0 * direction[i] * self.a_partial_slice()[i];
            b1 += direction[i] * (self.b_partial_slice()[i] + self.b);
            shift += direction[i]
                * self
                    .translation
                    .as_ref()
                    .map_or(0.0, |c| c.components[i].value);
            let bi = self
                .beta
                .as_ref()
                .map_or(self.model.bi[i], |b| b.components[i].value);
            for (j, &dj) in direction.iter().enumerate() {
                let bj = self
                    .beta
                    .as_ref()
                    .map_or(self.model.bi[j], |b| b.components[j].value);
                b2 += direction[i] * dj * self.model.pair_covolume(i, j, bi, bj);
                d2 += direction[i]
                    * dj
                    * (1.0 - self.model.vdw_kij(i, j))
                    * self.factors().0[i]
                    * self.factors().0[j];
            }
        }
        let attraction = real(self.a) + s * d1 + s * s * d2;
        let covolume = (real(self.b) + s * b1 + s * s * b2) / amount;
        let v = real(self.eos_volume(volume)) + s * shift;
        let (delta1, delta2) = self.model.kind.deltas();
        let integral = if delta1 == delta2 {
            real(1.0) / v
        } else {
            ((v + covolume * delta1) / (v + covolume * delta2)).ln()
                / (covolume * (delta1 - delta2))
        };
        Some(
            (ideal - amount * (v - covolume).ln() - attraction * integral / (R * self.temperature))
                .f123,
        )
    }

    /// Log fugacities (Pa) and their mole derivatives at fixed T and total V.
    ///
    /// The state contains one mole, so `volume` is both its molar and total
    /// volume. `dn` is a row-major n*n Hessian of A/(RT). Unlike TP derivatives,
    /// this expression remains regular at a mechanical spinodal and permits
    /// negative pressure. Absent species have log fugacity -infinity and an
    /// infinite ideal diagonal; their off-diagonal entries remain finite.
    /// Classic mixing and Wong-Sandler are supported; other activity rules
    /// still require their attraction composition Hessian.
    pub fn log_fugacity_tv_into(
        &self,
        volume: f64,
        log_fugacity: &mut [f64],
        dn: &mut [f64],
    ) -> Result<()> {
        let n = self.composition_slice().len();
        if log_fugacity.len() != n || dn.len() != n * n {
            return Err(Error::InvalidInput("TV fugacity buffer dimensions"));
        }
        if self.model.excess_gibbs.is_some() && self.wong_sandler.is_none() {
            return Err(Error::Unsupported(
                "TV attraction Hessian for this activity rule",
            ));
        }
        self.volume_check(volume)?;
        let v = self.eos_volume(volume);
        let p = self.eos_pressure(v)?;
        let rt = R * self.temperature;
        let free = v - self.b;
        let (d1, d2) = self.model.kind.deltas();
        let denom = (v + d1 * self.b) * (v + d2 * self.b);
        let denom_b = (d1 + d2) * v + 2.0 * d1 * d2 * self.b;
        let (integral, ib, ibb) = self.attraction_integral_derivatives(v);
        let pressure_n = |i: usize| {
            rt / free + rt * self.b_partial_slice()[i] / (free * free)
                - 2.0 * self.a_partial_slice()[i] / denom
                + self.a * denom_b * self.b_partial_slice()[i] / (denom * denom)
        };
        let shift = |i: usize| {
            self.translation
                .as_ref()
                .map_or(0.0, |c| c.components[i].value)
        };
        for i in 0..n {
            let bi = self.b_partial_slice()[i];
            let di = 2.0 * self.a_partial_slice()[i];
            log_fugacity[i] = (self.composition_slice()[i] * rt / free).ln() + bi / free
                - (di * integral + self.a * ib * bi) / rt
                - p.value * shift(i) / rt;
            for j in 0..n {
                let bj = self.b_partial_slice()[j];
                let dj = 2.0 * self.a_partial_slice()[j];
                let (bnn, dnn) = if let Some(ws) = &self.wong_sandler {
                    (ws.b_ij[i * n + j], ws.a_ij[i * n + j])
                } else {
                    let bii = self
                        .beta
                        .as_ref()
                        .map_or(self.model.bi[i], |b| b.components[i].value);
                    let bjj = self
                        .beta
                        .as_ref()
                        .map_or(self.model.bi[j], |b| b.components[j].value);
                    let bij = self.model.pair_covolume(i, j, bii, bjj);
                    (
                        2.0 * bij - bi - bj,
                        2.0 * (1.0 - self.model.vdw_kij(i, j))
                            * self.factors().0[i]
                            * self.factors().0[j],
                    )
                };
                dn[i * n + j] = (bi + bj + bnn) / free + bi * bj / (free * free)
                    - (dnn * integral
                        + (di * bj + dj * bi) * ib
                        + self.a * ibb * bi * bj
                        + self.a * ib * bnn)
                        / rt
                    - (shift(i) * pressure_n(j)
                        + shift(j) * pressure_n(i)
                        + shift(i) * shift(j) * p.dv)
                        / rt;
                if i == j {
                    dn[i * n + j] += 1.0 / self.composition_slice()[i];
                }
            }
        }
        Ok(())
    }

    /// Number of components of the model.
    #[inline]
    fn count(&self) -> usize {
        self.model.components.len()
    }

    /// Mole fractions of the prepared state.
    #[inline]
    fn composition_slice(&self) -> &[f64] {
        &self.buffer[..self.count()]
    }

    /// Partial attraction parameters of the prepared state.
    #[inline]
    fn a_partial_slice(&self) -> &[f64] {
        let n = self.count();
        &self.buffer[n..2 * n]
    }

    /// Partial covolumes of the prepared state.
    #[inline]
    fn b_partial_slice(&self) -> &[f64] {
        let n = self.count();
        &self.buffer[2 * n..3 * n]
    }

    /// Square-root factors of the prepared state.
    #[inline]
    fn factors(&self) -> (&[f64], &[f64], &[f64]) {
        let n = self.count();
        (
            &self.buffer[3 * n..4 * n],
            &self.buffer[4 * n..5 * n],
            &self.buffer[5 * n..6 * n],
        )
    }
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
    // Source-order PureSat initialization at a pure database critical T.
    // The L/mol arithmetic of cbCalcPseudo controls the T > Tpc branch;
    // reducing it algebraically to SI can change that comparison by one ulp.
    pub(crate) fn baseline_pure_critical_pressure_seed(&self, component: usize) -> Result<f64> {
        let c = &self.components[component];
        let (oa, ob, qc) = source_order_omega(self.kind);
        let kr = 1000. * R;
        let krt = kr * c.critical_temperature;
        let a0 = oa * (krt * krt) / c.critical_pressure;
        let b0 = ob * kr * c.critical_temperature / c.critical_pressure;
        let coefficients = |t: f64| -> Result<(f64, f64, f64, f64, f64)> {
            let alpha = self.alpha[component].evaluate(t, c.critical_temperature)?;
            let beta = self.beta[component].evaluate(t, c.critical_temperature)?;
            let f = (a0 * alpha.value).sqrt();
            let ft = (0.5 * a0 / f) * alpha.dt;
            let ftt =
                -0.25 * a0 / f * (alpha.dt * alpha.dt - 2. * alpha.value * alpha.dtt) / alpha.value;
            Ok((
                f * f,
                2. * f * ft,
                2. * (ft * ft + f * ftt),
                b0 * beta.value,
                b0 * beta.dt,
            ))
        };
        let mut tpc = c.critical_temperature;
        let mut last = None;
        let mut converged = false;
        for iteration in 0..20 {
            let (a, at, att, b, bt) = coefficients(tpc)?;
            let f = a - oa / ob * kr * tpc * b;
            let df = at - oa / ob * kr * (b + tpc * bt);
            let arg = 1. - 2. * f * att / (df * df);
            let arg = if arg < 0. { 1. } else { arg.sqrt() };
            let step = (-f) / df * 2. / (1. + arg);
            tpc += step;
            last = Some((a, b));
            if step.abs() <= 1e-10 {
                converged = iteration < 19;
                break;
            }
        }
        if !converged {
            return Err(Error::NotConverged("source pure pseudo critical"));
        }
        let (a, b) = last.ok_or(Error::NotConverged("source pure pseudo critical"))?;
        let v = (ob * qc) * b / ob;
        let (d1, d2) = self.kind.deltas();
        let ppc = kr * tpc / (v - b) - a / ((v + d2 * b) * (v + d1 * b));
        positive(ppc, "source pure pseudo critical pressure")?;
        let t = c.critical_temperature;
        if t > tpc {
            return Ok(ppc);
        }
        let mut z = vec![0.; self.components.len()];
        z[component] = 1.;
        let state = self.prepare(t, &z)?;
        let evaluate = |p: f64| -> Result<(f64, f64)> {
            let mut values = [0.; 2];
            let mut volumes = [0.; 2];
            for (i, phase) in [Phase::Liquid, Phase::Vapor].into_iter().enumerate() {
                let physical = state.baseline_volume_with_extrema(p, phase, true)?;
                let v = state.eos_volume(physical);
                let rt = R * t;
                let zz = p * v / rt;
                let integral = state.attraction_integral(v);
                let f = -(1. - state.b / v).ln() - state.a * integral / rt;
                values[i] = rt * f + p * v - rt * (1. + zz.ln());
                volumes[i] = v;
            }
            let result = (values[1] - values[0], volumes[1] - volumes[0]);
            if !result.0.is_finite() || !result.1.is_finite() {
                return Err(Error::NotConverged("source PureSat seed equations"));
            }
            Ok(result)
        };
        let low = 10.;
        let high = (1.5 * ppc * t / tpc).min(ppc);
        if high < low {
            return Err(Error::NotConverged("source PureSat pressure bounds"));
        }
        let mut p = (0.1 * ppc * t / tpc).max(low);
        if evaluate(p)?.1.abs() < 1e-12 {
            p = low + 1e-12;
            if evaluate(p)?.1.abs() < 1e-12 {
                p = ppc * t / tpc;
            }
        }
        for _ in 0..30 {
            let (f, df) = evaluate(p)?;
            if df == 0. {
                return Err(Error::NotConverged("source PureSat seed Jacobian"));
            }
            let mut step = -f / df;
            if step.abs() > 1e-9 {
                if p + step < low {
                    step = low - p;
                }
                if p + step > high {
                    step = high - p;
                }
            }
            let mut fraction = 1.;
            let mut next = (p + step).clamp(low, high);
            let mut residual = evaluate(next)?.0;
            for _ in 0..3 {
                if residual.abs() < f.abs() {
                    break;
                }
                fraction *= 0.5;
                next = (p + fraction * step).clamp(low, high);
                residual = evaluate(next)?.0;
            }
            p = next;
            if residual.abs() < 1e-4 {
                return Ok(p);
            }
        }
        Err(Error::NotConverged(
            "source PureSat critical initialization",
        ))
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
        if options.wong_sandler.is_some() && !options.mixing_rule.is_wong_sandler() {
            return Err(Error::InvalidInput(
                "the Wong-Sandler parameters need the Wong-Sandler mixing rule",
            ));
        }
        if options.mixing_rule.is_excess_gibbs() {
            if matches!(
                options.mixing_rule,
                crate::excess_gibbs::MixingRule::HuronVidal1
                    | crate::excess_gibbs::MixingRule::HuronVidal2
                    | crate::excess_gibbs::MixingRule::Nrtl
                    | crate::excess_gibbs::MixingRule::HvWongSandler
            ) {
                // cbSelect's GE path requests "Classic", whereas this fork's
                // binary records are labelled "vdW" (and cbMixHVWongSandler has
                // no branch in tpSelectInteractionParameters). Its lookup
                // returns zero and leaves covolume interactions at zero.
                model.kij.fill(0.0);
                model.lij.fill(0.0);
                model.rebuild_pairs();
            }
            let comp_idents: Vec<&str> = model.components.iter().map(|c| c.ident).collect();
            let database_eos = if kind == CubicKind::RedlichKwong {
                "SRK"
            } else {
                kind.name()
            };
            model.excess_gibbs = Some(crate::excess_gibbs::ExcessGibbsMixture::new(
                &comp_idents,
                options.mixing_rule,
                database_eos,
                Some(parameter_reference),
                &model.kij,
            )?);
            if let Some(data) = options.wong_sandler {
                model.excess_gibbs.as_mut().unwrap().wong_sandler = Some(data.clone());
                model.wong_sandler = Some(data.clone());
            }
            model.mixing_rule = options.mixing_rule;
        }
        Ok(model)
    }

    /// Construct a cubic base model for CPA with CPA-specific excess Gibbs
    /// initialization (`"CPA-SRK"` / `"CPA-PR"`).
    pub(crate) fn new_cpa_base(
        kind: CubicKind,
        names: &[&str],
        mixing_rule: crate::excess_gibbs::MixingRule,
        alpha: crate::AlphaCorrelation,
        parameter_reference: &str,
        wong_sandler: Option<&crate::wong_sandler::WongSandlerData>,
    ) -> Result<Self> {
        let mut model = Self::new_with_references(
            kind,
            names,
            alpha,
            parameter_reference,
            parameter_reference,
            false,
        )?;
        model.kij.fill(0.0);
        model.lij.fill(0.0);
        model.rebuild_pairs();
        if wong_sandler.is_some() && !mixing_rule.is_wong_sandler() {
            return Err(Error::InvalidInput(
                "the Wong-Sandler parameters need a Wong-Sandler mixing rule",
            ));
        }
        if mixing_rule.is_excess_gibbs() {
            let comp_idents: Vec<&str> = model.components.iter().map(|c| c.ident).collect();
            let cpa_eos_name = match kind {
                CubicKind::SoaveRedlichKwong => "CPA-SRK",
                CubicKind::PengRobinson => "CPA-PR",
                _ => kind.name(),
            };
            model.excess_gibbs = Some(crate::excess_gibbs::ExcessGibbsMixture::new(
                &comp_idents,
                mixing_rule,
                cpa_eos_name,
                Some(parameter_reference),
                &model.kij,
            )?);
            if mixing_rule == crate::excess_gibbs::MixingRule::WongSandler {
                let ws_data = wong_sandler
                    .cloned()
                    .unwrap_or_else(|| crate::wong_sandler::WongSandlerData::new(names.len()));
                model.excess_gibbs.as_mut().unwrap().wong_sandler = Some(ws_data.clone());
                model.wong_sandler = Some(ws_data);
            } else if let Some(data) = wong_sandler {
                model.excess_gibbs.as_mut().unwrap().wong_sandler = Some(data.clone());
                model.wong_sandler = Some(data.clone());
            }
            model.mixing_rule = mixing_rule;
        }
        Ok(model)
    }

    /// The pure component attraction parameters of the model.
    ///
    /// The values are the products of the critical attraction and the alpha
    /// function at the given temperature.
    pub fn pure_attraction(&self, temperature: f64) -> Result<Vec<f64>> {
        let count = self.components.len();
        let mut values = Vec::with_capacity(count);
        for i in 0..count {
            let critical_temperature = self.components[i].critical_temperature;
            let crate::AlphaValue { value: alpha, .. } =
                self.alpha[i].evaluate(temperature, critical_temperature)?;
            values.push(self.a0[i] * alpha);
        }
        Ok(values)
    }

    /// The excess Gibbs energy of the mixture at the infinite pressure limit.
    ///
    /// The values are the ones of the original excess_gibbs::gEinf.
    pub fn excess_gibbs_jet(
        &self,
        temperature: f64,
        composition: &[f64],
    ) -> Result<crate::excess_gibbs::ExcessGibbsJet> {
        let count = self.components.len();
        if composition.len() != count {
            return Err(Error::InvalidInput("composition mismatch"));
        }
        let mixture = self.excess_gibbs.as_ref().ok_or(Error::Unsupported(
            "the excess Gibbs model needs the mixing rule of the model",
        ))?;
        let mut a_pure = Vec::with_capacity(count);
        let mut da_dt_pure = Vec::with_capacity(count);
        let mut d2a_dt2_pure = Vec::with_capacity(count);
        for i in 0..count {
            let critical_temperature = self.components[i].critical_temperature;
            let crate::AlphaValue {
                value: alpha,
                dt: alpha_t,
                dtt: alpha_tt,
            } = self.alpha[i].evaluate(temperature, critical_temperature)?;
            a_pure.push(self.a0[i] * alpha);
            da_dt_pure.push(self.a0[i] * alpha_t);
            d2a_dt2_pure.push(self.a0[i] * alpha_tt);
        }
        crate::excess_gibbs::excess_gibbs_jet(
            mixture,
            self.kind,
            temperature,
            composition,
            &self.bi,
            &a_pure,
            &da_dt_pure,
            &d2a_dt2_pure,
        )
    }

    /// The mixture coefficients of the Wong-Sandler mixing rule.
    ///
    /// The rule replaces the binary attraction and the binary covolume of the
    /// cubic equation of state with its own coefficients. The excess Gibbs
    /// energy at the infinite pressure limit comes from the excess Gibbs
    /// model of the mixture, which the Wong-Sandler rule shares with the
    /// Huron-Vidal and the NRTL models.
    pub fn wong_sandler_coefficients(
        &self,
        temperature: f64,
        composition: &[f64],
    ) -> Result<crate::wong_sandler::WongSandlerCoefficients> {
        let count = self.components.len();
        if composition.len() != count {
            return Err(Error::InvalidInput("composition mismatch"));
        }
        let mixture = self.excess_gibbs.as_ref().ok_or(Error::Unsupported(
            "the Wong-Sandler mixing rule needs the excess Gibbs model",
        ))?;
        let mut a_pure = Vec::with_capacity(count);
        let mut da_dt_pure = Vec::with_capacity(count);
        let mut d2a_dt2_pure = Vec::with_capacity(count);
        for i in 0..count {
            let critical_temperature = self.components[i].critical_temperature;
            let crate::AlphaValue {
                value: alpha,
                dt: alpha_t,
                dtt: alpha_tt,
            } = self.alpha[i].evaluate(temperature, critical_temperature)?;
            a_pure.push(self.a0[i] * alpha);
            da_dt_pure.push(self.a0[i] * alpha_t);
            d2a_dt2_pure.push(self.a0[i] * alpha_tt);
        }
        let jet = crate::excess_gibbs::excess_gibbs_jet(
            mixture,
            self.kind,
            temperature,
            composition,
            &self.bi,
            &a_pure,
            &da_dt_pure,
            &d2a_dt2_pure,
        )?;
        let limit_c = crate::excess_gibbs::excess_gibbs_limit(
            self.kind,
            crate::excess_gibbs::GeLimit::Infinite,
        );
        if self.mixing_rule == crate::excess_gibbs::MixingRule::HvWongSandler {
            crate::wong_sandler::hv_wong_sandler_mix(
                &self.kij,
                &self.lij,
                temperature,
                composition,
                &a_pure,
                &da_dt_pure,
                &d2a_dt2_pure,
                &self.bi,
                limit_c,
                &jet,
            )
        } else {
            let data = self.wong_sandler.as_ref().ok_or(Error::Unsupported(
                "the Wong-Sandler mixing rule needs the binary parameters",
            ))?;
            crate::wong_sandler::wong_sandler_mix(
                data,
                temperature,
                composition,
                &a_pure,
                &da_dt_pure,
                &d2a_dt2_pure,
                &self.bi,
                limit_c,
                &jet,
            )
        }
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
        let components: Vec<_> = names
            .iter()
            .map(|name| component(name).copied())
            .collect::<Result<_>>()?;
        Self::from_components(
            kind,
            components,
            requested,
            alpha_reference,
            parameter_reference,
            quantum,
        )
    }

    /// Builds a cubic model with pseudo components.
    ///
    /// This is the native form of the source `init_cubic_pseudo`. The
    /// components named `PSEUDO` take the critical temperature, the critical
    /// pressure, the acentric factor and the molar mass of the corresponding
    /// list entries, and the remaining components keep the values of the
    /// source database. The lists are indexed by the position of the
    /// component, as in the source.
    pub fn new_pseudo(
        kind: CubicKind,
        names: &[&str],
        critical_temperatures: &[f64],
        critical_pressures: &[f64],
        acentric_factors: &[f64],
        molar_masses: Option<&[f64]>,
    ) -> Result<Self> {
        if names.is_empty() {
            return Err(Error::InvalidInput("at least one component required"));
        }
        if critical_temperatures.len() != names.len()
            || critical_pressures.len() != names.len()
            || acentric_factors.len() != names.len()
        {
            return Err(Error::InvalidInput(
                "the pseudo component lists must match the component list",
            ));
        }
        if let Some(masses) = molar_masses
            && masses.len() != names.len()
        {
            return Err(Error::InvalidInput(
                "the molar mass list must match the component list",
            ));
        }
        let mut components: Vec<_> = names
            .iter()
            .map(|name| component(name).copied())
            .collect::<Result<_>>()?;
        for (index, entry) in components.iter_mut().enumerate() {
            if !entry.ident.eq_ignore_ascii_case("PSEUDO") {
                continue;
            }
            entry.critical_temperature = critical_temperatures[index];
            entry.critical_pressure = critical_pressures[index];
            entry.acentric_factor = acentric_factors[index];
            if let Some(masses) = molar_masses {
                entry.molar_mass = masses[index];
            }
        }
        Self::from_components(
            kind,
            components,
            crate::AlphaCorrelation::Classic,
            "DEFAULT",
            "DEFAULT",
            false,
        )
    }

    fn from_components(
        kind: CubicKind,
        mut components: Vec<crate::component::Component>,
        requested: crate::AlphaCorrelation,
        alpha_reference: &str,
        parameter_reference: &str,
        quantum: bool,
    ) -> Result<Self> {
        let count = components.len();
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
            a0: Vec::with_capacity(count),
            bi: Vec::with_capacity(count),
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
            kij: vec![0.0; count.pow(2)],
            lij: vec![0.0; count.pow(2)],
            b_exponent: 1.0,
            nonzero_kij: Vec::new(),
            components,
            translation_parameters: None,
            translation_enabled: false,
            caloric_minimum: 80.,
            caloric_maximum: 999.,
            mixing_rule: crate::excess_gibbs::MixingRule::VdW,
            excess_gibbs: None,
            wong_sandler: None,
        };
        for c in &model.components {
            positive(c.critical_temperature, "critical temperature")?;
            positive(c.critical_pressure, "critical pressure")?;
            let rt_c = R * c.critical_temperature;
            model.a0.push(omega_a * (rt_c * rt_c) / c.critical_pressure);
            model
                .bi
                .push(omega_b * R * c.critical_temperature / c.critical_pressure);
        }
        for i in 0..count {
            for j in i + 1..count {
                let a = model.components[i].ident;
                let b = model.components[j].ident;
                let database_eos = if kind == CubicKind::RedlichKwong {
                    "SRK"
                } else {
                    kind.name()
                };
                let k = interaction_with_reference(database_eos, a, b, false, parameter_reference);
                let l = interaction_with_reference(database_eos, a, b, true, parameter_reference);
                model.kij[i * count + j] = k;
                model.kij[j * count + i] = k;
                model.lij[i * count + j] = l;
                model.lij[j * count + i] = l;
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
        pure.b_exponent = self.b_exponent;
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

    /// The matrix of the third virial coefficients of a binary model.
    ///
    /// This is the native form of the original `binaryThirdVirialCoeffMatrix`,
    /// which is only defined for a binary system. The pure coefficients come
    /// from the pure fluids. The cross coefficient of the pair follows from
    /// the two unsymmetric mixtures of the pair, with the linear system of the
    /// original. The matrix is stored row by row.
    pub fn binary_third_virial_coefficient_matrix(&self, temperature: f64) -> Result<Vec<f64>> {
        let count = self.components().len();
        if count != 2 {
            return Err(Error::InvalidInput(
                "the third virial coefficient matrix is only defined for a binary model",
            ));
        }
        let mut matrix = vec![0.0; 4];
        for i in 0..count {
            let mut composition = vec![0.0; count];
            composition[i] = 1.0;
            let state = self.prepare(temperature, &composition)?;
            matrix[i * count + i] = state.virial_coefficients()?.1;
        }
        let share = 0.75;
        // The two mixtures of the pair with the right hand side of the system
        // of the original.
        let mut mixtures = [(0.0, 0.0, 0.0); 2];
        for (index, (first, second)) in [(share, 1.0 - share), (1.0 - share, share)]
            .into_iter()
            .enumerate()
        {
            let mut composition = vec![0.0; count];
            composition[0] = first;
            composition[1] = second;
            let state = self.prepare(temperature, &composition)?;
            let mixed = state.virial_coefficients()?.1;
            let right =
                -mixed + matrix[0] * first * first * first + matrix[3] * second * second * second;
            mixtures[index] = (first, second, right);
        }
        let (one, two, three) = mixtures[0];
        let (four, five, six) = mixtures[1];
        let a11 = 3.0 * one * one * two;
        let a12 = 3.0 * two * two * one;
        let a21 = 3.0 * four * four * five;
        let a22 = 3.0 * (1.0 - share) * five * five * four;
        let determinant = a11 * a22 - a12 * a21;
        if determinant == 0.0 || !determinant.is_finite() {
            return Err(Error::NotConverged(
                "the third virial coefficient system of the original",
            ));
        }
        matrix[1] = (three * a22 - a12 * six) / determinant;
        matrix[2] = (a11 * six - three * a21) / determinant;
        Ok(matrix)
    }

    /// The matrix of the second virial coefficients of the model.
    ///
    /// The pure coefficients and the cross coefficients follow the original
    /// `secondVirialCoeffMatrix`. The cross coefficient of a pair follows from
    /// the virial coefficient of the equal mixture of the pair and the pure
    /// coefficients, with the equal mixture as the reference. The matrix is
    /// stored row by row.
    pub fn second_virial_coefficient_matrix(&self, temperature: f64) -> Result<Vec<f64>> {
        let count = self.components().len();
        let mut matrix = vec![0.0; count * count];
        for i in 0..count {
            let mut composition = vec![0.0; count];
            composition[i] = 1.0;
            let state = self.prepare(temperature, &composition)?;
            matrix[i * count + i] = state.virial_coefficients()?.0;
        }
        let share = 0.5;
        for i in 0..count {
            for j in i + 1..count {
                let mut composition = vec![0.0; count];
                composition[i] = share;
                composition[j] = 1.0 - share;
                let state = self.prepare(temperature, &composition)?;
                let mixed = state.virial_coefficients()?.0;
                let cross = (mixed
                    - matrix[i * count + i] * share * share
                    - matrix[j * count + j] * (1.0 - share) * (1.0 - share))
                    / (2.0 * share * (1.0 - share));
                matrix[i * count + j] = cross;
                matrix[j * count + i] = cross;
            }
        }
        Ok(matrix)
    }

    #[inline]
    fn vdw_kij(&self, i: usize, j: usize) -> f64 {
        let (lo, hi) = if i <= j { (i, j) } else { (j, i) };
        self.kij[lo * self.components.len() + hi]
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
        if let Some(ge) = &mut self.excess_gibbs {
            ge.pairs[i * n + j].kij = value;
            ge.pairs[j * n + i].kij = value;
        }
        self.rebuild_pairs();
        Ok(())
    }

    /// Set a single `(i, j)` entry of `kij` without modifying `(j, i)`,
    /// mirroring `thermopack_setkij`.
    pub fn set_kij_entry(&mut self, i: usize, j: usize, value: f64) -> Result<()> {
        self.validate_pair(i, j)?;
        if !value.is_finite() || (i == j && value != 0.0) {
            return Err(Error::InvalidInput(
                "finite kij with zero diagonal required",
            ));
        }
        let n = self.components.len();
        self.kij[i * n + j] = value;
        if let Some(ge) = &mut self.excess_gibbs {
            ge.pairs[i * n + j].kij = value;
        }
        self.rebuild_pairs();
        Ok(())
    }

    /// Get binary covolume interaction parameter `lij` for pair `(i, j)`
    /// (`thermopack_getlij`).
    pub fn lij(&self, i: usize, j: usize) -> Result<f64> {
        self.validate_pair(i, j)?;
        Ok(self.lij[i * self.components.len() + j])
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

    /// Power-mean covolume mixing exponent `s` (`1.0` = linear, `4/3` = VTPR, `2.0` = UMR).
    #[inline]
    pub fn b_exponent(&self) -> f64 {
        self.b_exponent
    }

    /// Set power-mean covolume mixing exponent `s > 0` (`cbCalcLowcasebij`).
    pub fn set_b_exponent(&mut self, exponent: f64) -> Result<()> {
        positive(exponent, "covolume mixing exponent")?;
        self.b_exponent = exponent;
        Ok(())
    }

    #[inline]
    fn has_power_covolume(&self) -> bool {
        (self.b_exponent - 1.0).abs() >= 1.0e-15
    }

    #[inline]
    fn uses_quadratic_covolume(&self) -> bool {
        self.has_power_covolume() || self.lij.iter().any(|&l| l != 0.0)
    }

    #[inline]
    fn pair_covolume(&self, i: usize, j: usize, bii: f64, bjj: f64) -> f64 {
        let n = self.components.len();
        let lij_fac = 1.0 - self.lij[i * n + j];
        if !self.has_power_covolume() {
            lij_fac * 0.5 * (bii + bjj)
        } else {
            let s = self.b_exponent;
            let inv_s = 1.0 / s;
            let g = 0.5 * (bii.powf(inv_s) + bjj.powf(inv_s));
            lij_fac * g.powf(s)
        }
    }

    #[inline]
    fn pair_covolume_with_t(
        &self,
        i: usize,
        j: usize,
        bi: crate::BetaValue,
        bj: crate::BetaValue,
    ) -> (f64, f64, f64) {
        let n = self.components.len();
        let lij_fac = 1.0 - self.lij[i * n + j];
        if !self.has_power_covolume() {
            (
                lij_fac * 0.5 * (bi.value + bj.value),
                lij_fac * 0.5 * (bi.dt + bj.dt),
                lij_fac * 0.5 * (bi.dtt + bj.dtt),
            )
        } else {
            let s = self.b_exponent;
            let p = 1.0 / s;
            let ui = bi.value.powf(p);
            let uj = bj.value.powf(p);
            let ui_t = p * ui * bi.dt / bi.value;
            let uj_t = p * uj * bj.dt / bj.value;
            let ui_tt = (p * ui / bi.value) * (bi.dtt + (p - 1.0) * bi.dt * bi.dt / bi.value);
            let uj_tt = (p * uj / bj.value) * (bj.dtt + (p - 1.0) * bj.dt * bj.dt / bj.value);
            let g = 0.5 * (ui + uj);
            let g_t = 0.5 * (ui_t + uj_t);
            let g_tt = 0.5 * (ui_tt + uj_tt);
            let bij = lij_fac * g.powf(s);
            let bij_t = s * bij * g_t / g;
            let bij_tt = (s * bij / g) * (g_tt + (s - 1.0) * g_t * g_t / g);
            (bij, bij_t, bij_tt)
        }
    }

    /// Power-mean pair covolume `lowcase_bij(i, j) = (0.5 * (b_i^(1/s) + b_j^(1/s)))^s`
    /// in SI units (`m^3/mol`), matching `cbmix::cbCalcLowcasebij`.
    pub fn lowcase_bij(&self, i: usize, j: usize) -> Result<f64> {
        self.validate_pair(i, j)?;
        if !self.has_power_covolume() {
            Ok(0.5 * (self.bi[i] + self.bi[j]))
        } else {
            let s = self.b_exponent;
            let inv_s = 1.0 / s;
            Ok((0.5 * (self.bi[i].powf(inv_s) + self.bi[j].powf(inv_s))).powf(s))
        }
    }

    /// Mixture molar covolume `Bmix(z)` [m^3/mol] evaluated via `lowcase_bij`,
    /// mirroring `cbmix::calcBmixCubic`.
    pub fn calc_bmix_cubic(&self, z: &[f64]) -> Result<f64> {
        let n = self.components.len();
        if z.len() != n || z.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(Error::InvalidInput("non-negative composition required"));
        }
        let sumn: f64 = z.iter().copied().sum();
        positive(sumn, "total moles")?;
        if !self.has_power_covolume() {
            let sum_b: f64 = z.iter().zip(&self.bi).map(|(&zi, &bi)| zi * bi).sum();
            return Ok(sum_b / sumn);
        }
        let mut bmix = 0.0;
        for (i, &zi) in z.iter().enumerate() {
            let mut bi = 0.0;
            for (j, &zj) in z.iter().enumerate() {
                bi += zj * self.lowcase_bij(i, j)?;
            }
            bmix += zi * bi;
        }
        Ok(bmix / (sumn * sumn))
    }

    pub(crate) fn validate_pair_indices(&self, i: usize, j: usize) -> Result<()> {
        self.validate_pair(i, j)
    }

    pub(crate) fn set_alpha_direct(&mut self, index: usize, alpha: crate::Alpha) -> Result<()> {
        let destination = self
            .alpha
            .get_mut(index)
            .ok_or(Error::InvalidInput("component index out of range"))?;
        *destination = alpha;
        Ok(())
    }

    pub(crate) fn set_component_criticals_internal(
        &mut self,
        i: usize,
        tc: f64,
        pc: f64,
        acf: f64,
        a0: f64,
        bi: f64,
    ) -> Result<()> {
        let comp = self
            .components
            .get_mut(i)
            .ok_or(Error::InvalidInput("component index out of range"))?;
        comp.critical_temperature = tc;
        comp.critical_pressure = pc;
        comp.acentric_factor = acf;
        self.a0[i] = a0;
        self.bi[i] = bi;
        Ok(())
    }

    pub(crate) fn set_pure_ab_override(&mut self, i: usize, a0: f64, bi: f64) -> Result<()> {
        positive(a0, "pure attraction parameter a0")?;
        positive(bi, "pure covolume parameter b")?;
        if i >= self.components.len() {
            return Err(Error::InvalidInput("component index out of range"));
        }
        self.a0[i] = a0;
        self.bi[i] = bi;
        Ok(())
    }

    pub(crate) fn pure_a0(&self) -> &[f64] {
        &self.a0
    }

    pub(crate) fn pure_bi(&self) -> &[f64] {
        &self.bi
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

    /// Value-only molar mixture attraction `a(T, x)` (`[Pa m^6 / mol^2]`) and
    /// covolume `b(T, x)` (`[m^3 / mol]`) without computing or allocating
    /// derivative buffers or composition Hessians.
    pub(crate) fn mixture_ab_values(
        &self,
        temperature: f64,
        composition: &[f64],
    ) -> Result<(f64, f64)> {
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
        if self.mixing_rule.is_wong_sandler() {
            let mixture = self.excess_gibbs.as_ref().ok_or(Error::Unsupported(
                "the Wong-Sandler mixing rule needs the excess Gibbs model",
            ))?;
            let mut a_pure = Vec::with_capacity(n);
            for i in 0..n {
                let tc = self.components[i].critical_temperature;
                let alpha = self.alpha[i].evaluate(temperature, tc)?.value;
                let ai = self.a0[i] * alpha;
                if !ai.is_finite() {
                    return Err(Error::InvalidInput("temperature outside numerical range"));
                }
                a_pure.push(ai);
            }
            let g_ex = crate::excess_gibbs::excess_gibbs_value(
                mixture,
                self.kind,
                temperature,
                composition,
                &self.bi,
                &a_pure,
            )?;
            let limit_c = crate::excess_gibbs::excess_gibbs_limit(
                self.kind,
                crate::excess_gibbs::GeLimit::Infinite,
            );
            let (a, b) = if self.mixing_rule == crate::excess_gibbs::MixingRule::HvWongSandler {
                crate::wong_sandler::hv_wong_sandler_ab_values(
                    &self.kij,
                    &self.lij,
                    temperature,
                    composition,
                    &a_pure,
                    &self.bi,
                    limit_c,
                    g_ex,
                )?
            } else {
                let data = self.wong_sandler.as_ref().ok_or(Error::Unsupported(
                    "the Wong-Sandler mixing rule needs the binary parameters",
                ))?;
                crate::wong_sandler::wong_sandler_ab_values(
                    data,
                    temperature,
                    composition,
                    &a_pure,
                    &self.bi,
                    limit_c,
                    g_ex,
                )?
            };
            positive(b, "mixture covolume")?;
            if !a.is_finite() {
                return Err(Error::InvalidInput("temperature outside numerical range"));
            }
            return Ok((a, b));
        }

        // Evaluate scalar covolume b(T, x) without derivative vectors.
        let mut b_pure = Vec::with_capacity(n);
        if self.variable_beta {
            for i in 0..n {
                let tc = self.components[i].critical_temperature;
                let bv = self.beta[i].evaluate(temperature, tc)?.value;
                b_pure.push(self.bi[i] * bv);
            }
        } else {
            b_pure.extend_from_slice(&self.bi);
        }
        let b = if !self.uses_quadratic_covolume() {
            composition
                .iter()
                .zip(b_pure.iter())
                .map(|(x, bi)| x * bi)
                .sum()
        } else {
            let mut b_sum = 0.0;
            for i in 0..n {
                let mut b_i = 0.0;
                for (j, &xj) in composition.iter().enumerate() {
                    b_i += xj * self.pair_covolume(i, j, b_pure[i], b_pure[j]);
                }
                b_sum += composition[i] * b_i;
            }
            b_sum
        };
        positive(b, "mixture covolume")?;

        if let Some(ge) = &self.excess_gibbs {
            let mut a_pure = Vec::with_capacity(n);
            for i in 0..n {
                let tc = self.components[i].critical_temperature;
                let alpha = self.alpha[i].evaluate(temperature, tc)?.value;
                let ai = self.a0[i] * alpha;
                if !ai.is_finite() {
                    return Err(Error::InvalidInput("temperature outside numerical range"));
                }
                a_pure.push(ai);
            }
            let a = crate::excess_gibbs::excess_gibbs_cubic_a_value(
                ge,
                self.kind,
                temperature,
                composition,
                &self.bi,
                &a_pure,
                b,
            )?;
            if !a.is_finite() {
                return Err(Error::InvalidInput("temperature outside numerical range"));
            }
            return Ok((a, b));
        }

        let mut f = Vec::with_capacity(n);
        let mut sum = 0.0;
        for (i, &x) in composition.iter().enumerate() {
            let tc = self.components[i].critical_temperature;
            let alpha = self.alpha[i].evaluate(temperature, tc)?.value;
            let fi = (self.a0[i] * alpha.max(0.0)).sqrt();
            f.push(fi);
            sum += x * fi;
        }
        let mut a = sum * sum;
        for &(i, j, k) in &self.nonzero_kij {
            a -= 2.0 * composition[i] * composition[j] * k * f[i] * f[j];
        }
        if !a.is_finite() {
            return Err(Error::InvalidInput("temperature outside numerical range"));
        }
        Ok((a, b))
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
        // The composition, the partial attraction parameters, the partial
        // covolumes and the three square-root factors share one allocation.
        let mut buffer = vec![0.0; 6 * n];
        buffer[..n].copy_from_slice(composition);
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
            buffer[3 * n + i] = f;
            buffer[4 * n + i] = ft;
            buffer[5 * n + i] = ftt;
            sum += x * f;
            sum_t += x * ft;
            sum_tt += x * ftt;
        }
        let mut a = sum * sum;
        let mut at = 2.0 * sum * sum_t;
        let mut att = 2.0 * (sum_t * sum_t + sum * sum_tt);
        let (head, factors) = buffer.split_at_mut(3 * n);
        let (_, tail) = head.split_at_mut(n);
        let (a_partial, b_partial) = tail.split_at_mut(n);
        let (f, ft, ftt) = (&factors[..n], &factors[n..2 * n], &factors[2 * n..3 * n]);
        for i in 0..n {
            a_partial[i] = f[i] * sum;
        }
        for &(i, j, k) in &self.nonzero_kij {
            let (fi, fit, fitt) = (f[i], ft[i], ftt[i]);
            let (fj, fjt, fjtt) = (f[j], ft[j], ftt[j]);
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
        match &beta {
            Some(beta) => {
                for (slot, component) in b_partial.iter_mut().zip(beta.components.iter()) {
                    *slot = component.value;
                }
            }
            None => b_partial.copy_from_slice(&self.bi),
        }
        if !self.uses_quadratic_covolume() {
            b = composition
                .iter()
                .zip(b_partial.iter())
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
                    if let Some(beta) = &beta {
                        let (bij, bij_t, bij_tt) =
                            self.pair_covolume_with_t(i, j, beta.components[i], beta.components[j]);
                        b_partial[i] += xj * bij;
                        bt_partial += xj * bij_t;
                        btt_partial += xj * bij_tt;
                    } else {
                        b_partial[i] += xj * self.pair_covolume(i, j, bii, self.bi[j]);
                    }
                }
                b += composition[i] * b_partial[i];
                if let Some(beta) = &mut beta {
                    beta.partial_t[i] = 2. * bt_partial;
                    beta.dt += composition[i] * bt_partial;
                    beta.dtt += composition[i] * btt_partial;
                }
            }
            for bi in b_partial.iter_mut() {
                *bi = 2.0 * *bi - b;
            }
            if let Some(beta) = &mut beta {
                for bit in &mut beta.partial_t {
                    *bit -= beta.dt;
                }
            }
        }
        positive(b, "mixture covolume")?;

        let wong_sandler = self.mixing_rule.is_wong_sandler();
        let mut ai_t = None;
        let mut ge_a_nn = None;
        let mut prepared_wong_sandler = None;
        if wong_sandler {
            // The Wong-Sandler rule replaces the attraction and the covolume
            // of the mixture with its own coefficients. The partial
            // attraction parameters of the native library are the half of the
            // partial derivatives of the original.
            let coefficients = self.wong_sandler_coefficients(temperature, composition)?;
            a = coefficients.a;
            at = coefficients.a_t;
            att = coefficients.a_tt;
            for (slot, &ai) in buffer[n..2 * n].iter_mut().zip(coefficients.a_i().iter()) {
                *slot = 0.5 * ai;
            }
            ai_t = Some(coefficients.a_it().iter().map(|&ait| 0.5 * ait).collect());
            b = coefficients.b;
            for (slot, &bi) in buffer[2 * n..3 * n]
                .iter_mut()
                .zip(coefficients.b_i().iter())
            {
                *slot = bi;
            }
            // The covolume of the rule depends on the temperature. The
            // carrier of the temperature dependent covolume of the classic
            // rules holds the derivatives of the covolume and of its partial
            // derivatives.
            beta = Some(crate::beta::PreparedBeta {
                partial_t: coefficients.b_it().to_vec(),
                components: coefficients
                    .b_i()
                    .iter()
                    .enumerate()
                    .map(|(index, &value)| crate::BetaValue {
                        value,
                        dt: coefficients.b_it()[index],
                        dtt: 0.0,
                    })
                    .collect(),
                dt: coefficients.b_t,
                dtt: coefficients.b_tt,
            });
            prepared_wong_sandler = Some(crate::wong_sandler::PreparedWongSandler {
                a_ij: coefficients.a_ij().to_vec(),
                b_ij: coefficients.b_ij().to_vec(),
            });
        } else if let Some(ge) = &self.excess_gibbs {
            let mut a_pure = Vec::with_capacity(n);
            let mut da_pure = Vec::with_capacity(n);
            let mut d2a_pure = Vec::with_capacity(n);
            for i in 0..n {
                let tc = self.components[i].critical_temperature;
                let crate::AlphaValue {
                    value: alpha,
                    dt: alpha_t,
                    dtt: alpha_tt,
                } = self.alpha[i].evaluate(temperature, tc)?;
                a_pure.push(self.a0[i] * alpha);
                da_pure.push(self.a0[i] * alpha_t);
                d2a_pure.push(self.a0[i] * alpha_tt);
            }
            let mut b_second = vec![0.; n * n];
            for i in 0..n {
                for j in 0..n {
                    let bij = self.pair_covolume(i, j, self.bi[i], self.bi[j]);
                    b_second[i * n + j] = 2.0 * bij - buffer[2 * n + i] - buffer[2 * n + j];
                }
            }
            let jet = crate::excess_gibbs::excess_gibbs_cubic_mix_with_covolume(
                ge,
                self.kind,
                temperature,
                composition,
                &self.bi,
                &a_pure,
                &da_pure,
                &d2a_pure,
                b,
                &buffer[2 * n..3 * n],
                &b_second,
            )?;
            for (slot, &ai) in buffer[n..2 * n].iter_mut().zip(jet.da_dni.iter()) {
                *slot = 0.5 * ai;
            }
            ai_t = Some(jet.d2a_dnid_t.iter().map(|&ait| 0.5 * ait).collect());
            ge_a_nn = Some(jet.d2a_dn_dn);
            a = jet.a_mix;
            at = jet.da_dt;
            att = jet.d2a_dt2;
        }

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
            buffer,
            a,
            at,
            att,
            b,
            ai_t,
            ge_a_nn,
            beta,
            translation,
            wong_sandler: prepared_wong_sandler,
        })
    }
    /// Complete 10-block reduced residual Helmholtz energy `F = A^res / (R T)`
    /// and all first and second `(T, V, n)` derivatives at `(temperature, volume, moles)`,
    /// evaluated in unit-mole coordinates with exact Euler scaling.
    pub fn residual_helmholtz_tv(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<CubicHelmholtz> {
        positive(temperature, "temperature")?;
        positive(volume, "volume")?;
        let nc = self.components.len();
        if moles.len() != nc || moles.iter().any(|&v| !v.is_finite() || v < 0.0) {
            return Err(Error::InvalidInput(
                "finite non-negative mole numbers matching the model required",
            ));
        }
        let total: f64 = moles.iter().copied().sum();
        positive(total, "total moles")?;
        let v_unit = volume / total;
        positive(v_unit, "molar volume")?;
        let x: Vec<f64> = moles.iter().map(|&v| v / total).collect();
        let state = self.prepare(temperature, &x)?;
        let mut h = state.residual_helmholtz_tv_unit(v_unit)?;
        if total != 1.0 {
            let inv_n = 1.0 / total;
            h.f *= total;
            h.f_t *= total;
            h.f_tt *= total;
            h.f_vv *= inv_n;
            for v in &mut h.f_vn {
                *v *= inv_n;
            }
            for v in &mut h.f_nn {
                *v *= inv_n;
            }
        }
        let all_finite = [h.f, h.f_t, h.f_v, h.f_tt, h.f_tv, h.f_vv]
            .iter()
            .all(|v| v.is_finite())
            && h.f_n.iter().all(|v| v.is_finite())
            && h.f_tn.iter().all(|v| v.is_finite())
            && h.f_vn.iter().all(|v| v.is_finite())
            && h.f_nn.iter().all(|v| v.is_finite());
        if !all_finite {
            return Err(Error::InvalidInput(
                "non-finite cubic residual Helmholtz derivative",
            ));
        }
        Ok(h)
    }
}

#[derive(Debug, Clone)]
pub(crate) struct CubicMixtureJet {
    pub(crate) b_t: f64,
    pub(crate) b_tt: f64,
    pub(crate) b_i: Vec<f64>,
    pub(crate) b_it: Vec<f64>,
    pub(crate) b_ij: Vec<f64>,
    pub(crate) d: f64,
    pub(crate) d_t: f64,
    pub(crate) d_tt: f64,
    pub(crate) d_i: Vec<f64>,
    pub(crate) d_it: Vec<f64>,
    pub(crate) d_ij: Vec<f64>,
}

impl PreparedState<'_> {
    /// Molar mixture attraction `a(T, x)` (`[Pa m^6 / mol^2]`) and covolume
    /// `b(T, x)` (`[m^3 / mol]`) without computing derivative Hessians.
    #[inline]
    pub(crate) fn mixture_ab(&self) -> (f64, f64) {
        (self.a, self.b)
    }

    #[inline]
    pub(crate) fn mixture_a_derivatives(&self) -> (f64, f64, f64) {
        (self.a, self.at, self.att)
    }

    /// Allocation-free scalar pressure jet `(f, f_v, f_vv, f_tv)` on a unit-mole
    /// basis (`N = 1`, molar volume `v_unit [m^3/mol]`).
    pub(crate) fn residual_helmholtz_pressure_jet_unit(
        &self,
        v_unit: f64,
    ) -> Result<(f64, f64, f64, f64)> {
        self.volume_check(v_unit)?;
        let t = self.temperature;
        let rt = R * t;
        let inv_rt = 1.0 / rt;
        let inv_t = 1.0 / t;
        let b_t = self.beta.as_ref().map_or(0.0, |beta| beta.dt);
        let (c, c_t) = self
            .translation
            .as_ref()
            .map_or((0.0, 0.0), |tr| (tr.c, tr.ct));
        let v_eos = v_unit + c;
        let free = v_eos - self.b;
        let inv_free2 = 1.0 / (free * free);
        let w_t = c_t - b_t;

        let f_rep = -((c - self.b) / v_unit).ln_1p();
        let f_rep_v = (c - self.b) / (v_unit * free);
        let f_rep_tv = inv_free2 * w_t;
        let f_rep_vv = (self.b - c) * (2.0 * v_unit + c - self.b) / (v_unit * v_unit * free * free);

        let (d1, d2) = self.model.kind.deltas();
        let (i_val, _, _) = self.attraction_integral_derivatives(v_eos);
        let denom = (v_eos + d1 * self.b) * (v_eos + d2 * self.b);
        let denom2 = denom * denom;
        let i_ve = -1.0 / denom;
        let i_veve = (2.0 * v_eos + (d1 + d2) * self.b) / denom2;
        let i_veb = ((d1 + d2) * v_eos + 2.0 * d1 * d2 * self.b) / denom2;

        let alpha = self.a * inv_rt;
        let alpha_t = (self.at - self.a * inv_t) * inv_rt;
        let i_tv = i_veve * c_t + i_veb * b_t;

        Ok((
            f_rep - alpha * i_val,
            f_rep_v - alpha * i_ve,
            f_rep_vv - alpha * i_veve,
            f_rep_tv - (alpha_t * i_ve + alpha * i_tv),
        ))
    }

    /// Complete 10-block reduced residual Helmholtz energy `F = A^res / (R T)`
    /// and all first and second `(T, V, n)` derivatives on a unit-mole basis
    /// (`N = 1`, molar volume `v_unit [m^3/mol]`).
    pub fn residual_helmholtz_tv_unit(&self, v_unit: f64) -> Result<CubicHelmholtz> {
        self.volume_check(v_unit)?;
        let nc = self.count();
        let t = self.temperature;
        let rt = R * t;
        let inv_rt = 1.0 / rt;
        let inv_t = 1.0 / t;

        let mix = self.mixture_parameter_jet(1.0)?;
        let (c, c_t, c_tt) = self
            .translation
            .as_ref()
            .map_or((0.0, 0.0, 0.0), |tr| (tr.c, tr.ct, tr.ctt));
        let v_eos = v_unit + c;
        let free = v_eos - self.b;
        let inv_free = 1.0 / free;
        let inv_free2 = inv_free * inv_free;

        // Repulsive + volume-shift term: F_rep = -N ln((V + C - B) / V)
        let w_t = c_t - mix.b_t;
        let w_tt = c_tt - mix.b_tt;
        let f_rep = -((c - self.b) / v_unit).ln_1p();
        let f_rep_t = -inv_free * w_t;
        let f_rep_v = (c - self.b) / (v_unit * free);
        let f_rep_tt = inv_free2 * w_t * w_t - inv_free * w_tt;
        let f_rep_tv = inv_free2 * w_t;
        let f_rep_vv = (self.b - c) * (2.0 * v_unit + c - self.b) / (v_unit * v_unit * free * free);

        // Attraction term: F_att = -(D / (R T)) * I(V_eos, B)
        let (d1, d2) = self.model.kind.deltas();
        let (i_val, i_b, i_bb) = self.attraction_integral_derivatives(v_eos);
        let denom = (v_eos + d1 * self.b) * (v_eos + d2 * self.b);
        let denom2 = denom * denom;
        let i_ve = -1.0 / denom;
        let i_veve = (2.0 * v_eos + (d1 + d2) * self.b) / denom2;
        let i_veb = ((d1 + d2) * v_eos + 2.0 * d1 * d2 * self.b) / denom2;

        let alpha = mix.d * inv_rt;
        let alpha_t = (mix.d_t - mix.d * inv_t) * inv_rt;
        let alpha_tt = (mix.d_tt - 2.0 * mix.d_t * inv_t + 2.0 * mix.d * inv_t * inv_t) * inv_rt;

        let i_t = i_ve * c_t + i_b * mix.b_t;
        let i_tv = i_veve * c_t + i_veb * mix.b_t;
        let i_tt = i_veve * c_t * c_t
            + 2.0 * i_veb * c_t * mix.b_t
            + i_bb * mix.b_t * mix.b_t
            + i_ve * c_tt
            + i_b * mix.b_tt;

        let f = f_rep - alpha * i_val;
        let f_t = f_rep_t - (alpha_t * i_val + alpha * i_t);
        let f_v = f_rep_v - alpha * i_ve;
        let f_tt = f_rep_tt - (alpha_tt * i_val + 2.0 * alpha_t * i_t + alpha * i_tt);
        let f_tv = f_rep_tv - (alpha_t * i_ve + alpha * i_tv);
        let f_vv = f_rep_vv - alpha * i_veve;

        let mut f_n = vec![0.0; nc];
        let mut f_tn = vec![0.0; nc];
        let mut f_vn = vec![0.0; nc];
        let mut f_nn = vec![0.0; nc * nc];

        let mut w_i_vec = vec![0.0; nc];
        let mut alpha_i_vec = vec![0.0; nc];
        let mut i_i_vec = vec![0.0; nc];

        for i in 0..nc {
            let (ci, cit) = self.translation.as_ref().map_or((0.0, 0.0), |tr| {
                (tr.components[i].value, tr.components[i].dt)
            });
            let bi = mix.b_i[i];
            let bit = mix.b_it[i];
            let wi = ci - bi;
            let wit = cit - bit;
            w_i_vec[i] = wi;

            let ai = mix.d_i[i] * inv_rt;
            let ait = (mix.d_it[i] - mix.d_i[i] * inv_t) * inv_rt;
            alpha_i_vec[i] = ai;

            let ii = i_ve * ci + i_b * bi;
            let ivi = i_veve * ci + i_veb * bi;
            let iti = i_veve * c_t * ci
                + i_veb * (c_t * bi + ci * mix.b_t)
                + i_bb * mix.b_t * bi
                + i_ve * cit
                + i_b * bit;
            i_i_vec[i] = ii;

            f_n[i] = (f_rep - inv_free * wi) - (ai * i_val + alpha * ii);
            f_tn[i] = (f_rep_t + inv_free2 * w_t * wi - inv_free * wit)
                - (ait * i_val + ai * i_t + alpha_t * ii + alpha * iti);
            f_vn[i] = (f_rep_v + inv_free2 * wi) - (ai * i_ve + alpha * ivi);
        }

        for i in 0..nc {
            let ci = self
                .translation
                .as_ref()
                .map_or(0.0, |tr| tr.components[i].value);
            let bi = mix.b_i[i];
            let wi = w_i_vec[i];
            let ai = alpha_i_vec[i];
            let ii = i_i_vec[i];
            for j in 0..nc {
                let cj = self
                    .translation
                    .as_ref()
                    .map_or(0.0, |tr| tr.components[j].value);
                let bj = mix.b_i[j];
                let wj = w_i_vec[j];
                let bij = mix.b_ij[i * nc + j];
                let aij = mix.d_ij[i * nc + j] * inv_rt;
                let iij =
                    i_veve * ci * cj + i_veb * (ci * bj + cj * bi) + i_bb * bi * bj + i_b * bij;
                let rep_ij = -inv_free * (wi + wj - bij) + inv_free2 * wi * wj;
                let att_ij = -(aij * i_val + ai * i_i_vec[j] + alpha_i_vec[j] * ii + alpha * iij);
                f_nn[i * nc + j] = rep_ij + att_ij;
            }
        }

        Ok(CubicHelmholtz {
            f,
            f_t,
            f_v,
            f_n,
            f_tt,
            f_tv,
            f_tn,
            f_vv,
            f_vn,
            f_nn,
        })
    }

    /// Extensive mixture parameters `B(T, n)` (`[m^3]`) and `D(T, n)` (`[Pa m^6]`)
    /// with full first and second derivatives with respect to `(T, n)` for total
    /// mole number `moles`, preserving the exact scaling of `cbmix::cbCalcMixtureParams`.
    pub(crate) fn mixture_parameter_jet(&self, moles: f64) -> Result<CubicMixtureJet> {
        positive(moles, "total moles")?;
        let n = self.composition_slice().len();
        let n2 = moles * moles;
        let inv_n = 1.0 / moles;

        let b_t = moles * self.beta.as_ref().map_or(0.0, |beta| beta.dt);
        let b_tt = moles * self.beta.as_ref().map_or(0.0, |beta| beta.dtt);
        let b_i = self.b_partial_slice().to_vec();
        let b_it = if let Some(beta) = &self.beta {
            beta.partial_t.clone()
        } else {
            vec![0.0; n]
        };

        let mut b_ij = vec![0.0; n * n];
        if let Some(ws) = &self.wong_sandler {
            for (dst, &src) in b_ij.iter_mut().zip(ws.b_ij.iter()) {
                *dst = src * inv_n;
            }
        } else if self.model.uses_quadratic_covolume() {
            for i in 0..n {
                let bii = self
                    .beta
                    .as_ref()
                    .map_or(self.model.bi[i], |beta| beta.components[i].value);
                for j in 0..n {
                    let bjj = self
                        .beta
                        .as_ref()
                        .map_or(self.model.bi[j], |beta| beta.components[j].value);
                    let bij = self.model.pair_covolume(i, j, bii, bjj);
                    b_ij[i * n + j] = (2.0 * bij - b_i[i] - b_i[j]) * inv_n;
                }
            }
        }

        let d = n2 * self.a;
        let d_t = n2 * self.at;
        let d_tt = n2 * self.att;
        let d_i: Vec<f64> = self
            .a_partial_slice()
            .iter()
            .map(|&ai| 2.0 * moles * ai)
            .collect();
        let d_it: Vec<f64> = if let Some(ait) = &self.ai_t {
            ait.iter().map(|&v| 2.0 * moles * v).collect()
        } else {
            let (f, ft, _) = self.factors();
            let comp = self.composition_slice();
            let mut out = vec![0.0; n];
            for i in 0..n {
                let mut sum = 0.0;
                for j in 0..n {
                    sum +=
                        comp[j] * (1.0 - self.model.vdw_kij(i, j)) * (ft[i] * f[j] + f[i] * ft[j]);
                }
                out[i] = 2.0 * moles * sum;
            }
            out
        };
        let d_ij: Vec<f64> = if let Some(ws) = &self.wong_sandler {
            ws.a_ij.clone()
        } else if let Some(ge_nn) = &self.ge_a_nn {
            ge_nn.clone()
        } else {
            let (f, _, _) = self.factors();
            let mut out = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..n {
                    out[i * n + j] = 2.0 * (1.0 - self.model.vdw_kij(i, j)) * f[i] * f[j];
                }
            }
            out
        };

        Ok(CubicMixtureJet {
            b_t,
            b_tt,
            b_i,
            b_it,
            b_ij,
            d,
            d_t,
            d_tt,
            d_i,
            d_it,
            d_ij,
        })
    }

    /// Pure-fluid reduced residual Helmholtz energy `alpr = A_res / (R T)` and its
    /// first and second derivatives `[alpr, alpr_T, alpr_v, alpr_TT, alpr_Tv, alpr_vv]`
    /// at molar volume `volume` (`[m^3/mol]`) for an untranslated, constant-covolume
    /// reference cubic EOS.
    pub(crate) fn pure_residual_helmholtz_tv(&self, volume: f64) -> Result<[f64; 6]> {
        self.volume_check(volume)?;
        let t = self.temperature;
        let rt = R * t;
        let free = volume - self.b;
        let (d1, d2) = self.model.kind.deltas();
        let integral = self.attraction_integral(volume);
        let denom = (volume + d1 * self.b) * (volume + d2 * self.b);
        let integral_v = -1.0 / denom;
        let integral_vv = (2.0 * volume + (d1 + d2) * self.b) / (denom * denom);

        let u = self.a / rt;
        let u_t = (self.at - self.a / t) / rt;
        let u_tt = (self.att - 2.0 * self.at / t + 2.0 * self.a / (t * t)) / rt;

        let alpr = -(-self.b / volume).ln_1p() - u * integral;
        let alpr_t = -u_t * integral;
        let alpr_tt = -u_tt * integral;
        let alpr_v = -self.b / (volume * free) - u * integral_v;
        let alpr_tv = -u_t * integral_v;
        let alpr_vv =
            self.b * (2.0 * volume - self.b) / (volume * volume * free * free) - u * integral_vv;
        let out = [alpr, alpr_t, alpr_v, alpr_tt, alpr_tv, alpr_vv];
        if out.iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput(
                "non-finite cubic reference Helmholtz jet",
            ));
        }
        Ok(out)
    }

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
        let n = self.composition_slice().len();
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
            rt / free + rt * self.b_partial_slice()[i] / (free * free)
                - 2.0 * self.a_partial_slice()[i] / denom
                + self.a * denom_b * self.b_partial_slice()[i] / (denom * denom)
        };
        for i in 0..n {
            let ai_t = if let Some(ait) = &self.ai_t {
                ait[i]
            } else {
                let mut s = 0.0;
                let (f, ft, _) = self.factors();
                for j in 0..n {
                    s += self.composition_slice()[j]
                        * (1.0 - self.model.vdw_kij(i, j))
                        * (ft[i] * f[j] + f[i] * ft[j]);
                }
                s
            };
            let bi = self.b_partial_slice()[i];
            let di = 2.0 * self.a_partial_slice()[i];
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
                let bj = self.b_partial_slice()[j];
                let dj = 2.0 * self.a_partial_slice()[j];
                let (bnn, dnn) = if let Some(ws) = &self.wong_sandler {
                    // The rule delivers the second mole derivatives of the
                    // attraction and of the covolume.
                    (ws.b_ij[i * n + j], ws.a_ij[i * n + j])
                } else {
                    let bii = self
                        .beta
                        .as_ref()
                        .map_or(self.model.bi[i], |b| b.components[i].value);
                    let bjj = self
                        .beta
                        .as_ref()
                        .map_or(self.model.bi[j], |b| b.components[j].value);
                    let bij = self.model.pair_covolume(i, j, bii, bjj);
                    (
                        2.0 * bij - bi - bj,
                        self.ge_a_nn.as_ref().map_or_else(
                            || {
                                2.0 * (1.0 - self.model.vdw_kij(i, j))
                                    * self.factors().0[i]
                                    * self.factors().0[j]
                            },
                            |h| h[i * n + j],
                        ),
                    )
                };
                let fnn = (bi + bj + bnn) / free + bi * bj / (free * free)
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
            let v2 = volume * volume;
            let v3 = v2 * volume;
            for k in 0_i32..12 {
                let sign = if k & 1 == 0 { 1.0 } else { -1.0 };
                let coefficient = sign
                    * (crate::numerics::ipow(d1, k + 1) - crate::numerics::ipow(d2, k + 1))
                    / ((d1 - d2) * (k + 1) as f64);
                i += coefficient * crate::numerics::ipow(ratio, k) / volume;
                if k >= 1 {
                    ib += coefficient * k as f64 * crate::numerics::ipow(ratio, k - 1) / v2;
                }
                if k >= 2 {
                    ibb += coefficient * (k * (k - 1)) as f64 * crate::numerics::ipow(ratio, k - 2)
                        / v3;
                }
            }
            (i, ib, ibb)
        } else {
            let integral = self.attraction_integral(volume);
            let t1 = d1 / (volume + d1 * self.b);
            let t2 = d2 / (volume + d2 * self.b);
            let log_b = t1 - t2;
            let log_bb = -(t1 * t1) + (t2 * t2);
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
        for (&x, gas) in self.composition_slice().iter().zip(&self.model.ideal_gases) {
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
            .composition_slice()
            .iter()
            .zip(&self.model.components)
            .map(|(x, c)| x * c.molar_mass * 1e-3)
            .sum();
        let pressure = self.pressure(volume)?;
        let squared =
            -volume * volume * (pressure.dv - self.temperature * (pressure.dt * pressure.dt) / cv)
                / mass;
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
                .composition_slice()
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
        self.composition_slice()
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

    /// The second volume derivative of the pressure at constant temperature
    /// in pascal mole squared per cubic metre. The volume translation is a
    /// constant shift and does not change it.
    pub fn pressure_second_derivative(&self, volume: f64) -> Result<f64> {
        self.volume_check(volume)?;
        let (d1, d2) = self.model.kind.deltas();
        let eos = self.eos_volume(volume);
        let free = eos - self.b;
        let denom = (eos + d1 * self.b) * (eos + d2 * self.b);
        let slope = 2.0 * eos + (d1 + d2) * self.b;
        Ok(2.0 * R * self.temperature / (free * free * free)
            + self.a * (2.0 / (denom * denom) - 2.0 * slope * slope / (denom * denom * denom)))
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

    /// The second and third virial coefficients of the state.
    ///
    /// The coefficients follow the definition of the original, the pressure
    /// expansion `P/(R T) = rho + B rho^2 + C rho^3` at the zero density. The
    /// original estimates them from the residual Helmholtz energy at a large
    /// molar volume, and the native model uses the same estimate.
    pub fn virial_coefficients(&self) -> Result<(f64, f64)> {
        // The large molar volume of the original estimate.
        const LARGE_VOLUME: f64 = 1.0e7;
        let volume = self.eos_volume(LARGE_VOLUME);
        let (d1, d2) = self.model.kind.deltas();
        let free = volume - self.b;
        let denom = (volume + d1 * self.b) * (volume + d2 * self.b);
        let slope = 2.0 * volume + (d1 + d2) * self.b;
        let thermal = R * self.temperature;
        // The reduced residual Helmholtz energy and its volume derivatives.
        let first = (-thermal * self.b / (volume * free) + self.a / denom) / thermal;
        let second = (thermal * self.b * (2.0 * volume - self.b) / (volume * volume * free * free)
            - self.a * slope / (denom * denom))
            / thermal;
        let lv2 = LARGE_VOLUME * LARGE_VOLUME;
        Ok((
            -lv2 * first,
            (lv2 * lv2) * (second + 2.0 * first / LARGE_VOLUME),
        ))
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
            cp: cv - self.temperature * (pressure.dt * pressure.dt) / pressure.dv - R,
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
        if output.len() != self.composition_slice().len() {
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
            let br = self.b_partial_slice()[i] / self.b;
            *phi = br * (z - 1.0)
                - repulsion
                - (2.0 * self.a_partial_slice()[i] - self.a * br) * integral;
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
pub(crate) fn omega_ab(kind: CubicKind) -> (f64, f64) {
    omega(kind)
}

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

// Preserve Fortran's exponentiation grouping as well as its iteration order.
// This compatibility helper is kept separate from the established SI kernels.
fn source_order_omega(kind: CubicKind) -> (f64, f64, f64) {
    let (d1, d2) = kind.deltas();
    let (m1, m2) = (-d2, -d1);
    let mut zcomb = 5.;
    let (mut n, mut n1, mut n2) = (0., 0., 0.);
    for _ in 0..9 {
        n = zcomb - 1.;
        n1 = zcomb - m1;
        n2 = zcomb - m2;
        let f = n * (n1 * n1) + n * n1 * n2 + n * (n2 * n2) - (n1 * n1) * n2 - n1 * (n2 * n2);
        let fd = 3. * (n * n1 + n * n2 - n1 * n2);
        let fdd = 6. * n;
        let arg = 1. - 2. * f * fdd / (fd * fd);
        let arg = if arg > 0. { arg.sqrt() } else { 1. };
        let step = (-2.) * f / (fd * (1. + arg));
        zcomb += step;
        if step.abs() <= 1e-10 {
            break;
        }
    }
    let arg = n1 * n2 / (n * n1 + n * n2);
    let b = (1. - arg) / (zcomb - 1.);
    (b * (1. / n - b) * n1 * n2, b, zcomb)
}

#[cfg(test)]
mod mixture_ab_value_tests {
    use super::*;
    use crate::excess_gibbs::MixingRule;
    use crate::wong_sandler::{Fraction, WongSandlerData};

    fn assert_ab_matches(label: &str, model: &Cubic, t: f64, x: &[f64]) -> Result<()> {
        let (a_val, b_val) = model.mixture_ab_values(t, x)?;
        let (a_prep, b_prep) = model.prepare(t, x)?.mixture_ab();
        let a_tol = 1e-13 * a_prep.abs().max(1.0);
        let b_tol = 1e-13 * b_prep.abs().max(1e-12);
        assert!(
            (a_val - a_prep).abs() <= a_tol,
            "{label} a mismatch: value={a_val:.16e}, prepare={a_prep:.16e}"
        );
        assert!(
            (b_val - b_prep).abs() <= b_tol,
            "{label} b mismatch: value={b_val:.16e}, prepare={b_prep:.16e}"
        );
        Ok(())
    }

    #[test]
    fn mixture_ab_values_matches_prepare_across_all_mixing_rules() -> Result<()> {
        // 1. VdW (with nonzero kij and lij)
        let mut vdw = Cubic::new(CubicKind::SoaveRedlichKwong, &["C1", "C3"])?;
        vdw.set_kij(0, 1, 0.04)?;
        vdw.set_lij(0, 1, 0.02)?;
        assert_ab_matches("VdW+lij", &vdw, 310.0, &[0.4, 0.6])?;

        // 2. HV1
        let hv1 = Cubic::with_options(
            CubicKind::SoaveRedlichKwong,
            &["C1", "H2O"],
            CubicOptions {
                mixing_rule: MixingRule::HuronVidal1,
                ..Default::default()
            },
        )?;
        assert_ab_matches("HV1", &hv1, 350.0, &[0.4, 0.6])?;

        // 3. HV2
        let hv2 = Cubic::with_options(
            CubicKind::PengRobinson,
            &["C1", "H2O"],
            CubicOptions {
                mixing_rule: MixingRule::HuronVidal2,
                ..Default::default()
            },
        )?;
        assert_ab_matches("HV2", &hv2, 400.0, &[0.3, 0.7])?;

        // 4. NRTL
        let nrtl = Cubic::with_options(
            CubicKind::PengRobinson,
            &["H2S", "C3"],
            CubicOptions {
                mixing_rule: MixingRule::Nrtl,
                parameter_reference: "Dicko2012",
                ..Default::default()
            },
        )?;
        assert_ab_matches("NRTL", &nrtl, 360.0, &[0.45, 0.55])?;

        // 5. UNIFAC
        let unifac = Cubic::with_options(
            CubicKind::SoaveRedlichKwong,
            &["C1", "C3"],
            CubicOptions {
                mixing_rule: MixingRule::Unifac,
                ..Default::default()
            },
        )?;
        assert_ab_matches("UNIFAC", &unifac, 320.0, &[0.3, 0.7])?;

        // 6. Wong-Sandler
        let mut ws = WongSandlerData::new(2);
        ws.set_kij(0, 1, Fraction::constant(0.1))?;
        ws.set_kij(1, 0, Fraction::constant(0.1))?;
        ws.set_tauij(0, 1, Fraction::constant(0.25))?;
        ws.set_tauij(1, 0, Fraction::constant(-0.15))?;
        ws.set_alphaij(0, 1, 0.3)?;
        ws.set_alphaij(1, 0, 0.3)?;
        let ws_model = Cubic::with_options(
            CubicKind::PengRobinson,
            &["CO2", "C1"],
            CubicOptions {
                mixing_rule: MixingRule::WongSandler,
                wong_sandler: Some(&ws),
                ..Default::default()
            },
        )?;
        assert_ab_matches("WongSandler", &ws_model, 320.0, &[0.5, 0.5])?;

        // 7. Overflowing Wong-Sandler activity sum with zero-mole component must
        // propagate an error instead of silently skipping a NaN row.
        let mut ws_overflow = WongSandlerData::new(2);
        ws_overflow.set_tauij(0, 1, Fraction::constant(-1000.0))?;
        ws_overflow.set_alphaij(0, 1, 1.0)?;
        let ws_overflow_model = Cubic::with_options(
            CubicKind::PengRobinson,
            &["CO2", "C1"],
            CubicOptions {
                mixing_rule: MixingRule::WongSandler,
                wong_sandler: Some(&ws_overflow),
                ..Default::default()
            },
        )?;
        assert!(ws_overflow_model.prepare(320.0, &[1.0, 0.0]).is_err());
        assert!(
            ws_overflow_model
                .mixture_ab_values(320.0, &[1.0, 0.0])
                .is_err()
        );

        Ok(())
    }
}
