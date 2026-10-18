//! Pure PeTS root, critical and coexistence calculations.
use super::*;
use crate::{Phase, Saturation};
use std::sync::OnceLock;

/// True pure-fluid critical point, with molar volume and compressibility.
#[derive(Debug, Clone, Copy)]
pub struct PureCriticalPoint {
    pub temperature: f64,
    pub pressure: f64,
    pub volume: f64,
    pub compressibility: f64,
}

/// The source nonassociation volume iterator always returns its finite last
/// iterate. These flags expose the two source stopping conditions separately.
#[derive(Debug, Clone, Copy)]
pub struct PetsVolumeIteration {
    pub volume: f64,
    pub iterations: usize,
    pub volume_converged: bool,
    pub pressure_converged: bool,
    pub relative_pressure_error: f64,
}

#[derive(Debug, Clone, Copy)]
struct ReducedCritical {
    temperature: f64,
    eta: f64,
    pressure: f64,
}

/// Temperature-prepared pure-fluid branches. Spinodals are found once and
/// reused across pressures; each query returns mechanically stable roots.
#[derive(Debug, Clone)]
pub struct PetsEquilibrium<'a> {
    prepared: PreparedPets<'a>,
    boundaries: Vec<f64>,
}
static CRITICAL: OnceLock<Result<ReducedCritical>> = OnceLock::new();

// Taylor coefficients in packing fraction, through fourth order. Fixed-size
// scalar arithmetic supplies analytic density derivatives for critical work.
#[derive(Clone, Copy)]
struct Taylor([f64; 5]);
impl Taylor {
    fn constant(value: f64) -> Self {
        Self([value, 0., 0., 0., 0.])
    }
    fn variable(value: f64) -> Self {
        Self([value, 1., 0., 0., 0.])
    }
    fn scale(self, value: f64) -> Self {
        Self(self.0.map(|x| x * value))
    }
    fn inverse(self) -> Self {
        let mut v = [0.; 5];
        v[0] = 1. / self.0[0];
        for k in 1..5 {
            let mut sum = 0.;
            for j in 1..=k {
                sum += self.0[j] * v[k - j];
            }
            v[k] = -sum / self.0[0];
        }
        Self(v)
    }
}
impl std::ops::Add for Taylor {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        let mut r = self;
        for i in 0..5 {
            r.0[i] += b.0[i];
        }
        r
    }
}
impl std::ops::Sub for Taylor {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + b.scale(-1.)
    }
}
impl std::ops::Mul for Taylor {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        let mut r = [0.; 5];
        for (k, v) in r.iter_mut().enumerate() {
            for j in 0..=k {
                *v += self.0[j] * b.0[k - j];
            }
        }
        Self(r)
    }
}

fn pressure_series(tr: f64, eta: f64) -> [f64; 4] {
    let e = Taylor::variable(eta);
    let free = Taylor::constant(1.) - e;
    let f2 = free * free;
    let f4 = f2 * f2;
    let hs = (e.scale(4.) - (e * e).scale(3.)) * f2.inverse();
    let c = f4 * (f4 + e.scale(8.) - (e * e).scale(2.)).inverse();
    let polynomial = |coefficients: &[f64; 7]| {
        let mut v = Taylor::constant(coefficients[6]);
        for &a in coefficients[..6].iter().rev() {
            v = v * e + Taylor::constant(a);
        }
        v
    };
    let diameter_ratio = 1. - 0.127112544 * (-3.052785558 / tr).exp();
    let q = 1. / tr;
    let alpha = hs
        - (e * (polynomial(&A).scale(2. * q) + (c * polynomial(&B)).scale(q * q)))
            .scale(6. / diameter_ratio.powi(3));
    let ar = alpha.0[1];
    let arr = 2. * alpha.0[2];
    let arrr = 6. * alpha.0[3];
    let arrrr = 24. * alpha.0[4];
    [
        eta + eta * eta * ar,
        1. + 2. * eta * ar + eta * eta * arr,
        2. * ar + 4. * eta * arr + eta * eta * arrr,
        6. * arr + 6. * eta * arrr + eta * eta * arrrr,
    ]
}

fn reduced_critical() -> Result<ReducedCritical> {
    CRITICAL
        .get_or_init(|| {
            let (mut t, mut e) = (1.1, 0.15);
            for _ in 0..50 {
                let f = pressure_series(t, e);
                if f[1].abs() < 2e-13 && (e * f[2]).abs() < 2e-13 {
                    return Ok(ReducedCritical {
                        temperature: t,
                        eta: e,
                        pressure: f[0],
                    });
                }
                // Only the temperature column is differenced; the density
                // derivatives, including the curvature column, are analytic.
                let dt = t * 1e-5;
                let lo = pressure_series(t - dt, e);
                let hi = pressure_series(t + dt, e);
                let a = (hi[1] - lo[1]) / (2. * dt);
                let b = f[2];
                let c = (hi[2] - lo[2]) / (2. * dt);
                let d = f[3];
                let determinant = a * d - b * c;
                if !determinant.is_finite() || determinant == 0. {
                    break;
                }
                let step_t = (-f[1] * d + b * f[2]) / determinant;
                let step_e = (-a * f[2] + f[1] * c) / determinant;
                let norm = f[1] * f[1] + (e * f[2]).powi(2);
                let mut length = 1.;
                let mut accepted = false;
                for _ in 0..20 {
                    let nt = t + length * step_t;
                    let ne = e + length * step_e;
                    if (0.2..4.).contains(&nt) && (0.001..0.7).contains(&ne) {
                        let nf = pressure_series(nt, ne);
                        if nf[1] * nf[1] + (ne * nf[2]).powi(2) < norm {
                            t = nt;
                            e = ne;
                            accepted = true;
                            break;
                        }
                    }
                    length *= 0.5;
                }
                if !accepted {
                    break;
                }
            }
            Err(Error::NotConverged("PeTS reduced critical point"))
        })
        .clone()
}

impl Pets {
    pub fn prepare_equilibrium(&self, t: f64) -> Result<PetsEquilibrium<'_>> {
        self.prepare(t)?.equilibrium()
    }
    pub fn volume(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        self.prepare_equilibrium(t)?.volume(p, phase)
    }
    pub fn saturation_pressure(&self, t: f64) -> Result<Saturation> {
        self.prepare_equilibrium(t)?.saturation_pressure()
    }
    pub fn saturation_temperature(&self, p: f64) -> Result<Saturation> {
        positive(p, "saturation pressure")?;
        let critical = self.critical_point()?;
        if p >= critical.pressure {
            return Err(Error::InvalidInput(
                "pressure at or above PeTS critical pressure",
            ));
        }
        let mut low = 0.5 * critical.temperature;
        let mut low_sat = self.saturation_pressure(low)?;
        for _ in 0..20 {
            if low_sat.pressure <= p {
                break;
            }
            low *= 0.8;
            low_sat = self.saturation_pressure(low)?;
        }
        if low_sat.pressure > p {
            return Err(Error::NotConverged("PeTS saturation temperature bracket"));
        }
        let mut high = critical.temperature * (1. - 1e-8);
        let mut t =
            (critical.temperature / (1. - (p / critical.pressure).ln() / 7.)).clamp(low, high);
        for _ in 0..80 {
            let sat = self.saturation_pressure(t)?;
            let error = sat.pressure - p;
            if error.abs() <= 1e-10 * p {
                return Ok(sat);
            }
            if error < 0. {
                low = t;
            } else {
                high = t;
            }
            let next = t - error / sat.dpdt;
            t = if next > low && next < high {
                next
            } else {
                low + (high - low) * 0.5
            };
        }
        Err(Error::NotConverged("PeTS saturation temperature"))
    }
    pub fn critical_point(&self) -> Result<PureCriticalPoint> {
        let critical = reduced_critical()?;
        let t = critical.temperature * self.epsilon_k;
        let prepared = self.prepare(t)?;
        let result = PureCriticalPoint {
            temperature: t,
            volume: prepared.packing / critical.eta,
            pressure: R * t / prepared.packing * critical.pressure,
            compressibility: critical.pressure / critical.eta,
        };
        if [
            result.temperature,
            result.volume,
            result.pressure,
            result.compressibility,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.)
        {
            return Err(Error::InvalidInput("non-finite PeTS critical point"));
        }
        Ok(result)
    }
    pub fn baseline_volume_iteration(
        &self,
        t: f64,
        p: f64,
        phase: Phase,
    ) -> Result<PetsVolumeIteration> {
        self.prepare(t)?.baseline_volume_iteration(p, phase)
    }
}

impl<'a> PreparedPets<'a> {
    fn pressure_eta(&self, eta: f64) -> Result<(f64, f64, f64)> {
        let rho = eta / self.packing;
        let contributions = self.contributions(rho)?;
        let a = contributions.total;
        let rt = R * self.t;
        let pressure = rt * rho * (1. + rho * a.drho);
        let derivative = rt / self.packing * (1. + 2. * rho * a.drho + rho * rho * a.drhorho);
        let scale = rt
            * rho
            * (1.
                + rho
                    * (contributions.hard_sphere.drho.abs() + contributions.dispersion.drho.abs()));
        if [pressure, derivative, scale].iter().any(|x| !x.is_finite()) {
            return Err(Error::InvalidInput("non-finite PeTS density pressure"));
        }
        Ok((pressure, derivative, scale))
    }
    pub fn equilibrium(&self) -> Result<PetsEquilibrium<'a>> {
        let critical = reduced_critical()?;
        let mut grid = Vec::with_capacity(390);
        grid.push(0.);
        grid.push(critical.eta);
        for i in 0..=128 {
            grid.push(10_f64.powf(-14. + 13. * i as f64 / 128.));
        }
        for i in 1..=256 {
            grid.push(0.75 * i as f64 / 256.);
        }
        grid.sort_by(f64::total_cmp);
        grid.dedup();
        let mut boundaries = vec![0.];
        let mut left = grid[0];
        let mut derivative = self.pressure_eta(left)?.1;
        for right in grid.into_iter().skip(1) {
            let next = self.pressure_eta(right)?.1;
            if (derivative < 0.) != (next < 0.) {
                let (mut lo, mut hi) = (left, right);
                let negative = derivative < 0.;
                for _ in 0..80 {
                    let mid = lo + (hi - lo) * 0.5;
                    if mid == lo || mid == hi {
                        break;
                    }
                    if (self.pressure_eta(mid)?.1 < 0.) == negative {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                boundaries.push(lo + (hi - lo) * 0.5);
            }
            left = right;
            derivative = next;
        }
        boundaries.push(0.75);
        Ok(PetsEquilibrium {
            prepared: *self,
            boundaries,
        })
    }
    pub fn volume(&self, p: f64, phase: Phase) -> Result<f64> {
        self.equilibrium()?.volume(p, phase)
    }
    pub fn baseline_volume_iteration(
        &self,
        p_spec: f64,
        phase: Phase,
    ) -> Result<PetsVolumeIteration> {
        positive(p_spec, "pressure")?;
        if phase == Phase::MinimumGibbs {
            return Err(Error::InvalidInput(
                "baseline SAFT volume requires liquid or vapor seed",
            ));
        }
        let mut zeta = if phase == Phase::Vapor { 1e-10 } else { 0.5 };
        let (mut low, mut high) = (1e-10, 0.75);
        let mut v = self.packing / zeta;
        let mut pressure = self.pressure_tv(v, 1.)?;
        let mut f = (1. - zeta) * (pressure.value - p_spec);
        let mut iter = 0;
        loop {
            let derivative = p_spec - pressure.value - ((1. - zeta) / zeta) * pressure.dv * v;
            zeta -= f / derivative;
            if zeta > high || zeta < low {
                zeta = (low + high) / 2.;
            }
            let old = v;
            v = self.packing / zeta;
            pressure = self.pressure_tv(v, 1.)?;
            f = (1. - zeta) * (pressure.value - p_spec);
            let v_converged = (v - old).abs() < 1e5 * f64::EPSILON * old;
            let p_converged = (pressure.value - p_spec).abs() < 1e8 * f64::EPSILON * p_spec;
            if v_converged && p_converged {
                return Ok(PetsVolumeIteration {
                    volume: v,
                    iterations: iter + 1,
                    volume_converged: true,
                    pressure_converged: true,
                    relative_pressure_error: (pressure.value - p_spec) / p_spec,
                });
            }
            iter += 1;
            if v_converged || iter > 40 {
                return Ok(PetsVolumeIteration {
                    volume: v,
                    iterations: iter,
                    volume_converged: v_converged,
                    pressure_converged: p_converged,
                    relative_pressure_error: (pressure.value - p_spec) / p_spec,
                });
            }
            if f > 0. {
                high = zeta;
            } else if f < 0. {
                low = zeta;
            }
        }
    }
}

impl PetsEquilibrium<'_> {
    /// Molar spinodal volumes, sorted ascending. No critical/end-point
    /// spinodal is invented for a monotone isotherm.
    pub fn spinodal_volumes(&self) -> Vec<f64> {
        self.boundaries[1..self.boundaries.len() - 1]
            .iter()
            .rev()
            .map(|e| self.prepared.packing / e)
            .collect()
    }
    fn root_eta(&self, p: f64, mut lo: f64, mut hi: f64) -> Result<f64> {
        let s = &self.prepared;
        let ideal = p * s.packing / (R * s.t);
        let mut e = if ideal > lo && ideal < hi {
            ideal
        } else {
            lo + (hi - lo) * 0.5
        };
        for _ in 0..120 {
            let (pressure, derivative, scale) = s.pressure_eta(e)?;
            let error = pressure - p;
            let tolerance = 2e-12 * p + 32. * f64::EPSILON * scale;
            let step = error / derivative;
            if error.abs() <= tolerance
                && (step.abs() <= 16. * f64::EPSILON * e || hi - lo <= 32. * f64::EPSILON * e)
            {
                return Ok(e);
            }
            if error > 0. {
                hi = e;
            } else {
                lo = e;
            }
            let next = e - step;
            let next = if next > lo && next < hi && derivative > 0. {
                next
            } else {
                lo + (hi - lo) * 0.5
            };
            if next == e || next == lo || next == hi {
                if error.abs() <= tolerance {
                    return Ok(e);
                }
                break;
            }
            e = next;
        }
        Err(Error::NotConverged("PeTS pressure branch"))
    }
    pub fn volumes(&self, p: f64) -> Result<Vec<f64>> {
        positive(p, "pressure")?;
        let reduced = reduced_critical()?;
        if (self.prepared.t / self.prepared.model.epsilon_k - reduced.temperature).abs()
            <= 8. * f64::EPSILON * reduced.temperature
        {
            let critical = self.prepared.model.critical_point()?;
            if (p - critical.pressure).abs() <= 64. * f64::EPSILON * critical.pressure {
                return Ok(vec![critical.volume]);
            }
        }
        let mut roots = Vec::new();
        let s = &self.prepared;
        for interval in self.boundaries.windows(2) {
            let (lo, hi) = (interval[0], interval[1]);
            if s.pressure_eta(lo + (hi - lo) * 0.5)?.1 < 0. {
                continue;
            }
            let (pl, _, sl) = s.pressure_eta(lo)?;
            let (ph, _, sh) = s.pressure_eta(hi)?;
            let tolerance = 32. * f64::EPSILON * sl.max(sh);
            if p < pl - tolerance || p > ph + tolerance {
                continue;
            }
            let e = self.root_eta(p, lo, hi)?;
            let v = s.packing / e;
            if v.is_finite() && v > 0. {
                roots.push(v);
            }
        }
        roots.sort_by(f64::total_cmp);
        roots.dedup_by(|a, b| (*a - *b).abs() <= 32. * f64::EPSILON * a.abs().max(b.abs()));
        if roots.is_empty() {
            return Err(Error::NoPhysicalRoot);
        }
        Ok(roots)
    }
    pub fn volume(&self, p: f64, phase: Phase) -> Result<f64> {
        let roots = self.volumes(p)?;
        match phase {
            Phase::Liquid => Ok(roots[0]),
            Phase::Vapor => Ok(roots[roots.len() - 1]),
            Phase::MinimumGibbs => {
                let mut result = roots[0];
                let mut best = f64::INFINITY;
                for v in roots {
                    let a = self.prepared.residual_density(1. / v)?.value;
                    let g = a + p * v / (R * self.prepared.t) - v.ln();
                    if g < best {
                        best = g;
                        result = v;
                    }
                }
                Ok(result)
            }
        }
    }
    pub fn saturation_pressure(&self) -> Result<Saturation> {
        let s = &self.prepared;
        let critical = reduced_critical()?;
        if s.t >= critical.temperature * s.model.epsilon_k || self.boundaries.len() < 4 {
            return Err(Error::InvalidInput(
                "temperature at or above PeTS critical temperature",
            ));
        }
        let pmax = s.pressure_eta(self.boundaries[1])?.0;
        let pmin = s.pressure_eta(self.boundaries[2])?.0;
        if pmax <= 0. || pmax <= pmin {
            return Err(Error::NotConverged("PeTS coexistence pressure bounds"));
        }
        let evaluate = |logp: f64| -> Result<(f64, f64, f64)> {
            let p = logp.exp();
            let roots = self.volumes(p)?;
            if roots.len() < 2 {
                return Err(Error::NotConverged("distinct PeTS coexistence roots"));
            }
            let vl = roots[0];
            let vv = roots[roots.len() - 1];
            let al = s.residual_density(1. / vl)?.value;
            let av = s.residual_density(1. / vv)?.value;
            Ok((
                al - av + p * (vl - vv) / (R * s.t) + vv.ln() - vl.ln(),
                vl,
                vv,
            ))
        };
        let margin = 1e-7 * (pmax - pmin.max(0.));
        let mut high = (pmax - margin).ln();
        let mut low = if pmin > 0. {
            (pmin + margin).ln()
        } else {
            high - 30.
        };
        let mut f_low = evaluate(low)?.0;
        let f_high = evaluate(high)?.0;
        while f_low < 0. && pmin <= 0. && low > -680. {
            low -= 30.;
            f_low = evaluate(low)?.0;
        }
        if f_low < 0. || f_high > 0. {
            return Err(Error::NotConverged(
                "PeTS coexistence chemical potential bracket",
            ));
        }
        let mut logp = low + (high - low) * 0.5;
        for _ in 0..120 {
            let (f, vl, vv) = evaluate(logp)?;
            let p = logp.exp();
            let derivative = p * (vl - vv) / (R * s.t);
            if (f.abs() < 5e-13 && (f / derivative).abs() < 1e-11) || high - low < 2e-13 {
                let fl = s.residual_helmholtz(vl, 1.)?;
                let fv = s.residual_helmholtz(vv, 1.)?;
                let latent_heat = -R * s.t * s.t * (fv.dt - fl.dt) + p * (vv - vl);
                return Ok(Saturation {
                    temperature: s.t,
                    pressure: p,
                    liquid_volume: vl,
                    vapor_volume: vv,
                    dpdt: latent_heat / (s.t * (vv - vl)),
                });
            }
            if f > 0. {
                low = logp;
            } else {
                high = logp;
            }
            let next = logp - f / derivative;
            logp = if next > low && next < high {
                next
            } else {
                low + (high - low) * 0.5
            };
        }
        Err(Error::NotConverged("PeTS saturation pressure"))
    }
}
