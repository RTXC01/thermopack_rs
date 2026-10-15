//! Conventional cubic equations, classic alpha and van der Waals mixing.
//! Ported from cbmix.f90, cbalpha.f90, cubic.f90 and cbhelm.f90.
use crate::component::interaction;
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
    components: Vec<&'static Component>,
    a0: Vec<f64>,
    bi: Vec<f64>,
    kappa: Vec<f64>,
    kij: Vec<f64>,
    lij: Vec<f64>,
    nonzero_kij: Vec<(usize, usize, f64)>,
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
        if state.pressure(center)?.dv <= 0.0 {
            return Err(Error::InvalidInput(
                "temperature at or above EOS critical temperature",
            ));
        }
        let bisect_spinodal = |mut low: f64, mut high: f64| -> Result<f64> {
            let negative_at_low = state.pressure(low)?.dv < 0.0;
            for _ in 0..80 {
                let mid = (low + high) * 0.5;
                if (state.pressure(mid)?.dv < 0.0) == negative_at_low {
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
            if state.pressure(high_volume)?.dv < 0.0 {
                break;
            }
            high_volume *= 2.0;
        }
        let vapor_spinodal = bisect_spinodal(center, high_volume)?;
        let min_pressure = state.pressure(liquid_spinodal)?.value;
        let max_pressure = state.pressure(vapor_spinodal)?.value;
        let margin = 1e-4 * (max_pressure - min_pressure.max(0.0));
        let mut low = (min_pressure + margin).max(max_pressure * 1e-9).ln();
        let mut high = (max_pressure - margin).ln();
        let evaluate = |log_pressure: f64| -> Result<(f64, f64, f64)> {
            let (volumes, count) = state.volumes(log_pressure.exp())?;
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
            let df = state.at - (oa / ob) * R * state.b;
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
        if names.is_empty() {
            return Err(Error::InvalidInput("at least one component required"));
        }
        let components: Vec<_> = names
            .iter()
            .map(|name| component(name))
            .collect::<Result<_>>()?;
        let (omega_a, omega_b) = omega(kind);
        let mut model = Self {
            kind,
            a0: Vec::with_capacity(names.len()),
            bi: Vec::with_capacity(names.len()),
            kappa: Vec::with_capacity(names.len()),
            kij: vec![0.0; names.len().pow(2)],
            lij: vec![0.0; names.len().pow(2)],
            nonzero_kij: Vec::new(),
            components,
        };
        for c in &model.components {
            positive(c.critical_temperature, "critical temperature")?;
            positive(c.critical_pressure, "critical pressure")?;
            let w = c.acentric_factor;
            model
                .a0
                .push(omega_a * (R * c.critical_temperature).powi(2) / c.critical_pressure);
            model
                .bi
                .push(omega_b * R * c.critical_temperature / c.critical_pressure);
            model.kappa.push(match kind {
                CubicKind::PengRobinson => 0.37464 + 1.54226 * w - 0.26992 * w * w,
                CubicKind::SoaveRedlichKwong => 0.48 + 1.574 * w - 0.176 * w * w,
                _ => 0.0,
            });
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
                let k = interaction(database_eos, a, b, false);
                let l = interaction(database_eos, a, b, true);
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
    pub fn components(&self) -> &[&'static Component] {
        &self.components
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
            let tr = temperature / tc;
            let k = self.kappa[i];
            let (alpha, alpha_t, alpha_tt) = match self.kind {
                CubicKind::VanDerWaals => (1.0, 0.0, 0.0),
                CubicKind::RedlichKwong => {
                    let alpha = 1.0 / tr.sqrt();
                    (
                        alpha,
                        -0.5 * alpha / temperature,
                        0.75 * alpha / temperature.powi(2),
                    )
                }
                _ => {
                    let yy = -2.0 * k * (k + 1.0);
                    (
                        1.0 + k * k * (tr - 1.0) + yy * (tr.sqrt() - 1.0),
                        k * k / tc + 0.5 * yy / (temperature * tc).sqrt(),
                        -yy / (4.0 * temperature * (temperature * tc).sqrt()),
                    )
                }
            };
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
        let mut b_partial = self.bi.clone();
        if self.lij.iter().all(|l| *l == 0.0) {
            b = composition.iter().zip(&self.bi).map(|(x, bi)| x * bi).sum();
        } else {
            for i in 0..n {
                b_partial[i] = 0.0;
                for (j, &xj) in composition.iter().enumerate() {
                    b_partial[i] +=
                        xj * (1.0 - self.lij[i * n + j]) * (self.bi[i] + self.bi[j]) / 2.0;
                }
                b += composition[i] * b_partial[i];
            }
            for bi in &mut b_partial {
                *bi = 2.0 * *bi - b;
            }
        }
        positive(b, "mixture covolume")?;
        if !a.is_finite() || !at.is_finite() || !att.is_finite() {
            return Err(Error::InvalidInput("temperature outside numerical range"));
        }
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
        let p = self.pressure(volume)?;
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
            let fnt =
                -(2.0 * ai_t - di / t) * integral / rt - (self.at - self.a / t) * ib * bi / rt;
            dt[i] = fnt + pi * p.dt / (rt * p.dv) + 1.0 / t;
            dp[i] = -pi / (rt * p.dv) - 1.0 / p.value;
            for j in 0..n {
                let bj = self.b_partial[j];
                let dj = 2.0 * self.a_partial[j];
                let bij =
                    (1.0 - self.model.lij[i * n + j]) * (self.model.bi[i] + self.model.bi[j]) / 2.0;
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
            enthalpy: 0.0,
            entropy: 0.0,
        };
        for (&x, c) in self.composition.iter().zip(&self.model.components) {
            if x == 0.0 {
                continue;
            }
            let pure = c.ideal(self.temperature, pressure)?;
            result.cp += x * pure.cp;
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

    /// Single-phase heat capacities (Cp, Cv) [J/(mol K)].
    pub fn heat_capacities(&self, volume: f64) -> Result<(f64, f64)> {
        let ideal_cp = self.ideal(self.pressure(volume)?.value)?.cp;
        let residual = self.residual(volume)?;
        Ok((ideal_cp + residual.cp, ideal_cp - R + residual.cv))
    }

    /// Compatibility with this fork's `enthalpy(..., dhdt=True)`.
    /// For vdW only, preserves the erroneous F_TT expression in cubic.f90.
    /// Use `heat_capacities` for a derivative consistent with enthalpy values.
    /// See docs/rust/baseline-defects.md for the exact discrepancy.
    pub fn baseline_enthalpy_derivative_tp(&self, volume: f64) -> Result<f64> {
        let cp = self.heat_capacities(volume)?.0;
        if self.model.kind == CubicKind::VanDerWaals {
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
        let (cp, cv) = self.heat_capacities(volume)?;
        let mass: f64 = self
            .composition
            .iter()
            .zip(&self.model.components)
            .map(|(x, c)| x * c.molar_mass * 1e-3)
            .sum();
        let squared = -volume * volume * self.pressure(volume)?.dv * cp / (cv * mass);
        positive(squared, "sound speed squared")?;
        Ok(squared.sqrt())
    }
    pub fn temperature(&self) -> f64 {
        self.temperature
    }
    pub fn composition(&self) -> &[f64] {
        &self.composition
    }
    pub fn covolume(&self) -> f64 {
        self.b
    }

    fn volume_check(&self, volume: f64) -> Result<()> {
        if volume.is_finite() && volume > self.b {
            Ok(())
        } else {
            Err(Error::InvalidInput("molar volume must exceed covolume"))
        }
    }

    pub fn pressure(&self, volume: f64) -> Result<Pressure> {
        self.volume_check(volume)?;
        let (d1, d2) = self.model.kind.deltas();
        let denom = (volume + d1 * self.b) * (volume + d2 * self.b);
        let free = volume - self.b;
        Ok(Pressure {
            value: R * self.temperature / free - self.a / denom,
            dt: R / free - self.at / denom,
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

    pub fn residual(&self, volume: f64) -> Result<Residual> {
        let pressure = self.pressure(volume)?;
        positive(pressure.value, "pressure for T-P residual properties")?;
        let z = pressure.value * volume / (R * self.temperature);
        let integral = self.attraction_integral(volume);
        let cv = self.temperature * self.att * integral;
        Ok(Residual {
            compressibility: z,
            helmholtz_tv: -R * self.temperature * (-self.b / volume).ln_1p() - self.a * integral,
            enthalpy: R * self.temperature * (z - 1.0)
                + (self.temperature * self.at - self.a) * integral,
            entropy: R * (z * (1.0 - self.b / volume)).ln() + self.at * integral,
            cv,
            cp: cv - self.temperature * pressure.dt.powi(2) / pressure.dv - R,
        })
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
        let rt = R * self.temperature;
        let z = p * volume / rt;
        let repulsion = (z * (1.0 - self.b / volume)).ln();
        let integral = self.attraction_integral(volume) / rt;
        for (i, phi) in output.iter_mut().enumerate() {
            let br = self.b_partial[i] / self.b;
            *phi = br * (z - 1.0) - repulsion - (2.0 * self.a_partial[i] - self.a * br) * integral;
        }
        Ok(())
    }

    /// Physical roots at positive pressure, ascending by volume.
    pub fn volumes(&self, pressure: f64) -> Result<([f64; 3], usize)> {
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
                let p = self.pressure(v)?;
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
            let p = self.pressure(v)?;
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
