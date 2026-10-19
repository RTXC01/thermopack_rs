//! Source-order MBWR density iteration, including phase fallbacks and raw
//! last iterates. This is distinct from bounded mechanically stable root sets.
//!
//! ```
//! use thermopack_rs::{Mbwr, MbwrKind, Phase};
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let model = Mbwr::new(MbwrKind::Mbwr32, "C1")?;
//! let initialization = model.source_density_parameters()?;
//! let prepared = model.prepare_baseline_density(150., initialization)?;
//! let result = prepared.density(1e6, Phase::Vapor)?;
//! assert!(result.pressure_converged && result.mechanically_stable);
//! assert!(result.density > 0.);
//! # Ok(())
//! # }
//! ```
use super::*;
use crate::Phase;

const KR: f64 = 1000. * R;
const RHO_EPS: f64 = 5000. * f64::EPSILON;
const P_EPS: f64 = 1e-6;

/// Explicit initialization data used by the original density routine. Units
/// are K, Pa and mol/m3. The source initializer's coarse critical search is
/// separate; its table values, estimates or user parameters can be passed here.
#[derive(Debug, Clone, Copy)]
pub struct MbwrDensityParameters {
    pub critical_temperature: f64,
    pub critical_pressure: f64,
    pub critical_density: f64,
    pub acentric_factor: f64,
}
impl MbwrDensityParameters {
    pub fn new(tc: f64, pc: f64, rhoc: f64, acentric_factor: f64) -> Result<Self> {
        let p = Self {
            critical_temperature: tc,
            critical_pressure: pc,
            critical_density: rhoc,
            acentric_factor,
        };
        p.validate()?;
        Ok(p)
    }
    fn validate(self) -> Result<()> {
        positive(
            self.critical_temperature,
            "MBWR density critical temperature",
        )?;
        positive(self.critical_pressure, "MBWR density critical pressure")?;
        positive(self.critical_density, "MBWR density critical density")?;
        if !self.acentric_factor.is_finite() {
            return Err(Error::InvalidInput(
                "finite MBWR density acentric factor required",
            ));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MbwrDensityPath {
    RequestedPhaseNewton,
    CurvatureSwitchNewton,
    OtherPhaseNewton,
    BareNewtonFallback,
    LegacyTplibFallback,
}
/// Finite original last iterate and phase flag. A finite return is not
/// necessarily a converged, positive, or mechanically stable density.
#[derive(Debug, Clone, Copy)]
pub struct MbwrDensityIteration {
    pub density: f64,
    pub phase: Phase,
    pub pressure: f64,
    pub pressure_density_derivative: f64,
    pub relative_pressure_error: f64,
    pub pressure_converged: bool,
    pub mechanically_stable: bool,
    pub path: MbwrDensityPath,
}
/// Owned fixed-temperature coefficients in original MBWR_coef/makeParam order.
/// Density queries do not allocate or change a global active model.
#[derive(Debug, Clone)]
pub struct PreparedMbwrDensity {
    t: f64,
    parameters: [f64; 16],
    polynomial_count: usize,
    exponential_count: usize,
    critical: MbwrDensityParameters,
    b_srk: f64,
    m_srk: f64,
    a0_srk: f64,
}

impl Mbwr {
    /// Original default critical/acentric initialization, cached per model.
    /// Missing critical tables use the source's coarse grid search.
    pub fn source_density_parameters(&self) -> Result<MbwrDensityParameters> {
        Ok(self.source_initialization()?.parameters)
    }
    pub fn prepare_baseline_density(
        &self,
        t: f64,
        critical: MbwrDensityParameters,
    ) -> Result<PreparedMbwrDensity> {
        positive(t, "temperature")?;
        critical.validate()?;
        let (a, polynomial_count, exponential_count) = self.source_pressure_coefficients(t)?;
        let tc = critical.critical_temperature;
        let pc = critical.critical_pressure;
        let w = critical.acentric_factor;
        let result = PreparedMbwrDensity {
            t,
            parameters: a,
            polynomial_count,
            exponential_count,
            critical,
            b_srk: 0.08664 * KR * tc / pc,
            m_srk: 0.480 + 1.547 * w - 0.176 * w.powi(2),
            a0_srk: 0.42747 * (KR * tc).powi(2) / pc,
        };
        if [result.b_srk, result.m_srk, result.a0_srk]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(Error::InvalidInput("non-finite MBWR density preparation"));
        }
        Ok(result)
    }
    pub(super) fn source_pressure_coefficients(&self, t: f64) -> Result<([f64; 16], usize, usize)> {
        positive(t, "temperature")?;
        let c = &self.coefficients;
        let t2 = t * t;
        let t3 = t2 * t;
        let mut a = [0.; 16];
        a[0] = self.gamma;
        a[1] = KR * t;
        let (polynomial_count, exponential_count) = if self.kind == MbwrKind::Bender19 {
            // Original MBWR_coef(1), including its subtraction/division tree.
            a[2] = c[0] * t - c[1] - c[2] / t - c[3] / t2 - c[4] / t3;
            a[3] = c[5] * t + c[6] + c[7] / t;
            a[4] = c[8] * t + c[9];
            a[5] = c[10] * t + c[11];
            a[6] = c[12];
            a[7] = (c[13] + c[14] / t + c[15] / t2) / t2;
            a[8] = (c[16] + c[17] / t + c[18] / t2) / t2;
            for v in &mut a[2..9] {
                *v *= 1000.;
            }
            (6, 2)
        } else {
            a[2] = c[0] * t + c[1] * t.sqrt() + c[2] + c[3] / t + c[4] / t2;
            a[3] = c[5] * t + c[6] + c[7] / t + c[8] / t2;
            a[4] = c[9] * t + c[10] + c[11] / t;
            a[5] = c[12];
            a[6] = c[13] / t + c[14] / t2;
            a[7] = c[15] / t;
            a[8] = c[16] / t + c[17] / t2;
            a[9] = c[18] / t2;
            a[10] = (c[19] + c[20] / t) / t2;
            a[11] = (c[21] + c[22] / t2) / t2;
            a[12] = (c[23] + c[24] / t) / t2;
            a[13] = (c[25] + c[26] / t2) / t2;
            a[14] = (c[27] + c[28] / t) / t2;
            a[15] = (c[29] + c[30] / t + c[31] / t2) / t2;
            (9, 6)
        };
        if a.iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput("non-finite MBWR pressure coefficients"));
        }
        Ok((a, polynomial_count, exponential_count))
    }
    pub fn baseline_density(
        &self,
        t: f64,
        p: f64,
        phase: Phase,
        critical: MbwrDensityParameters,
    ) -> Result<MbwrDensityIteration> {
        self.prepare_baseline_density(t, critical)?
            .density(p, phase)
    }
}

impl PreparedMbwrDensity {
    /// Gamma followed by original BP/BE coefficients, with backend density in
    /// mol/L. These parameters are read-only and include the ideal RT term.
    pub fn source_parameters(&self) -> &[f64] {
        &self.parameters[..1 + self.polynomial_count + self.exponential_count]
    }
    // Direct loop order from MBWR_pressure; negative intermediate densities
    // are deliberately retained for the legacy iterator. Backend units mol/L.
    #[inline]
    fn pressure<const SECOND: bool>(&self, rho: f64) -> [f64; 3] {
        if self.polynomial_count == 6 {
            Self::pressure_shape::<6, 2, true, SECOND>(&self.parameters, rho)
        } else {
            Self::pressure_shape::<9, 6, true, SECOND>(&self.parameters, rho)
        }
    }
    pub(super) fn pressure_value_backend(
        parameters: &[f64; 16],
        polynomial_count: usize,
        rho: f64,
    ) -> f64 {
        if polynomial_count == 6 {
            Self::pressure_shape::<6, 2, false, false>(parameters, rho)[0]
        } else {
            Self::pressure_shape::<9, 6, false, false>(parameters, rho)[0]
        }
    }
    // The source-indexed loop also measured faster than the iterator rewrite
    // on the recorded Windows GNU workload; keep its fixed family bounds.
    #[allow(clippy::needless_range_loop)]
    fn pressure_shape<
        const POLYNOMIAL: usize,
        const EXPONENTIAL: usize,
        const FIRST: bool,
        const SECOND: bool,
    >(
        parameters: &[f64; 16],
        rho: f64,
    ) -> [f64; 3] {
        let gamma = parameters[0];
        let r2 = rho * rho;
        // Positive integer-power trees used by the original __powidf2 loop.
        // Sharing them avoids repeated runtime powi calls without reordering
        // the polynomial/Gaussian sums or their rounding-sensitive products.
        let r3 = rho * r2;
        let r4 = r2 * r2;
        let r5 = rho * r4;
        let r8 = r4 * r4;
        let powers = [
            1.,
            rho,
            r2,
            r3,
            r4,
            r5,
            r2 * r4,
            r3 * r4,
            r8,
            rho * r8,
            r2 * r8,
            r3 * r8,
            r4 * r8,
            r5 * r8,
        ];
        let ex = (-gamma * r2).exp();
        let mut poly = [0.; 3];
        let mut expo = [0.; 3];
        for n in 1..=POLYNOMIAL {
            let term = parameters[n] * powers[n];
            poly[0] += term;
            if FIRST {
                poly[1] += term * n as f64;
            }
            if SECOND {
                poly[2] += term * (n * (n - 1)) as f64;
            }
        }
        for m in 1..=EXPONENTIAL {
            let term = parameters[m + POLYNOMIAL] * powers[2 * m + 1];
            expo[0] += term;
            if FIRST {
                expo[1] += term * ((2 * m + 1) as f64 / rho - 2. * gamma * rho);
            }
            if SECOND {
                expo[2] += term
                    * (2.
                        * ((m * (2 * m + 1)) as f64 / r2
                            - gamma * ((4 * m + 3) as f64 - 2. * gamma * r2)));
            }
        }
        [
            poly[0] + ex * expo[0],
            if FIRST {
                poly[1] / rho + ex * expo[1]
            } else {
                0.
            },
            if SECOND {
                poly[2] / r2 + ex * expo[2]
            } else {
                0.
            },
        ]
    }
    fn liquid_guess(&self, p: f64) -> f64 {
        let tr = self.t / self.critical.critical_temperature;
        if tr > 1.05 {
            return p / (KR * self.t);
        }
        let pr = p / self.critical.critical_pressure;
        if tr > 0.95 && pr > 0.95 && tr < 1.01 && pr < 1.05 {
            return self.critical.critical_density * 1e-3 * 1.5;
        }
        let a = self.a0_srk * (1. + self.m_srk * (1. - tr.sqrt())).powi(2);
        let big_a = a * p / (KR * self.t).powi(2);
        let big_b = self.b_srk * p / (KR * self.t);
        let r = (big_a - big_b - big_b.powi(2)) - 1. / 3.;
        let q = -2. / 27. + 1. / 3. * (big_a - big_b - big_b.powi(2)) - big_a * big_b;
        let discr = (r / 3.).powi(3) + (q / 2.).powi(2);
        let cube = |v: f64| v.abs().powf(1. / 3.).copysign(v);
        let y = if discr >= 1e-11 {
            let u = cube(-q / 2. + discr.sqrt());
            if u.abs() > 1e-16 {
                u - r / (3. * u)
            } else {
                return p / (KR * self.t);
            }
        } else if discr >= 0. {
            let u = (-q / 2.).min(cube(-q / 2. + discr.sqrt()));
            if -u < -1. / 3. { u - r / (3. * u) } else { -u }
        } else {
            let theta = (-r.powi(3) / 27.).sqrt();
            let phi = (-q / (2. * theta)).acos();
            let pi = 4. * 1_f64.atan();
            let y1 = (phi / 3.).cos();
            let y2 = ((phi + 2. * pi) / 3.).cos();
            let y3 = ((phi + 4. * pi) / 3.).cos();
            2. * theta.powf(1. / 3.) * y1.min(y2).min(y3)
        };
        let rho = p / (KR * self.t * (y + 1. / 3.));
        if rho <= 1e-14 || rho > 1e6 {
            p / (KR * self.t)
        } else {
            rho
        }
    }
    fn newton(&self, p: f64, initial: f64, phase: Phase, known: Option<[f64; 3]>) -> f64 {
        let mut q = known.unwrap_or_else(|| self.pressure::<false>(initial));
        let mut rho = initial;
        let mut old = 1.;
        let mut dev: f64 = 5e10;
        let mut previous = [0.; 3];
        let tr = self.t / self.critical.critical_temperature;
        let mut iter = 0;
        while dev.abs() / p > P_EPS || (old - rho).abs() / old > RHO_EPS {
            if iter > 0 {
                let drho = rho - old;
                if !((drho / old).abs() < 1e-9 || tr > 1.01) {
                    let missing = match phase {
                        Phase::Vapor => {
                            q[1] < 0.
                                || q[0] < 0.999 * previous[0]
                                || (q[0] - previous[0]) / drho < 0.9999 * q[1]
                        }
                        Phase::Liquid => {
                            q[1] < 0.
                                || q[0] > 1.001 * previous[0]
                                || (q[0] - previous[0]) / drho < 0.9999 * q[1]
                        }
                        Phase::MinimumGibbs => unreachable!(),
                    };
                    if missing {
                        return -1.;
                    }
                }
            }
            old = rho;
            previous = q;
            dev = q[0] - p;
            rho = old - dev / q[1];
            iter += 1;
            if iter == 15 || rho < 0. {
                return -1.;
            }
            q = self.pressure::<false>(rho);
        }
        rho
    }
    fn bare_newton(&self, p: f64, initial: f64) -> f64 {
        let mut rho = initial;
        let mut old = 1.;
        let mut dev: f64 = 5e10;
        let mut iter = 0;
        while dev.abs() / p > P_EPS || (old - rho).abs() / old > RHO_EPS {
            old = rho;
            let q = self.pressure::<false>(old);
            dev = q[0] - p;
            rho = old - dev / q[1];
            iter += 1;
            if q[1] < 0. || iter > 20 || rho < 0. {
                return -1.;
            }
        }
        rho
    }
    fn legacy(&self, p: f64, phase: Phase) -> f64 {
        let starts = [self.liquid_guess(p) * 1.4, 2.5e-4];
        let index = if phase == Phase::Liquid { 0 } else { 1 };
        let mut rho = starts[index];
        // The original evaluates liquid-start pressure even for a vapor D.
        let mut q = self.pressure::<true>(starts[0]);
        let mut dp = p - q[0];
        let mut new_iteration = true;
        let mut iteration = 0;
        let mut root = 0;
        while new_iteration {
            iteration += 1;
            let dps = dp;
            let ds = rho;
            let mut sq = (q[1] / q[2]).powi(2) + 2. * dp / q[2];
            let pdiv = -q[1] / q[2];
            let sign = if pdiv > 0. { -1. } else { 1. };
            if sq < 0. {
                sq = 0.;
            }
            let mut dd = pdiv + sign * sq.sqrt();
            if dd > 2. {
                dd = 2.;
            }
            rho += dd;
            let mut reduction = true;
            let mut count = 0;
            while reduction {
                count += 1;
                q = self.pressure::<true>(rho);
                dp = p - q[0];
                if dd.abs() <= 1e-8 {
                    reduction = false;
                }
                if reduction {
                    if dp.abs() <= dps.abs() {
                        reduction = false;
                    } else {
                        if 1. - dp / dps <= 1e-2 {
                            dd /= 2.;
                        } else if dp / dps - 1. >= 1e-2 {
                            dd = -dd;
                        } else {
                            reduction = false;
                        }
                        rho = ds + dd;
                    }
                }
                if count > 20 {
                    reduction = false;
                    new_iteration = false;
                }
                q = self.pressure::<true>(rho);
            }
            if dd.abs() <= 1e-8 {
                root += 1;
                let try_other = rho < 0. || q[1].abs() < 1e-2 || (dp / p).abs() > 1e-2;
                if try_other && root < 2 {
                    new_iteration = true;
                    rho = starts[1 - index];
                    q = self.pressure::<true>(rho);
                    dp = p - q[0];
                } else {
                    return rho;
                }
                // GUNDER=1 always returns directly for a wrong-curvature
                // phase root. The requested phase flag is retained by caller.
            }
            if iteration > 50 {
                new_iteration = false;
            }
        }
        rho
    }
    pub fn density(&self, p: f64, requested: Phase) -> Result<MbwrDensityIteration> {
        positive(p, "pressure")?;
        if requested == Phase::MinimumGibbs {
            return Err(Error::InvalidInput(
                "baseline MBWR density requires liquid or vapor seed",
            ));
        }
        let mut phase = requested;
        let mut rho = f64::NAN;
        let vapor = if p < 100. { 1e-12 } else { 1e-6 };
        let mut path = MbwrDensityPath::RequestedPhaseNewton;
        if phase == Phase::Liquid {
            let guess = self.liquid_guess(p) * 1.35;
            let q = self.pressure::<true>(guess);
            if q[2] < 0. || q[1] < 0. {
                phase = Phase::Vapor;
                path = MbwrDensityPath::CurvatureSwitchNewton;
            } else {
                rho = self.newton(p, guess, phase, Some(q));
            }
        }
        if phase == Phase::Vapor {
            rho = self.newton(p, vapor, phase, None);
        }
        if rho < 0. {
            phase = if phase == Phase::Liquid {
                Phase::Vapor
            } else {
                Phase::Liquid
            };
            let guess = if phase == Phase::Liquid {
                self.liquid_guess(p) * 1.35
            } else {
                vapor
            };
            rho = self.newton(p, guess, phase, None);
            path = MbwrDensityPath::OtherPhaseNewton;
        }
        if rho <= 0. || rho.is_nan() {
            phase = Phase::Liquid;
            path = MbwrDensityPath::BareNewtonFallback;
            rho = self.bare_newton(p, self.liquid_guess(p + 1e5) * 1.35 * 0.8);
            if rho < 0. {
                rho = self.legacy(p, requested);
                phase = requested;
                path = MbwrDensityPath::LegacyTplibFallback;
            }
        }
        let q = self.pressure::<false>(rho);
        let relative = (q[0] - p) / p;
        let density = rho * 1000.;
        let pressure_density_derivative = q[1] * 1e-3;
        if [density, q[0], pressure_density_derivative, relative]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(Error::NotConverged(
                "non-finite source MBWR density iterate",
            ));
        }
        Ok(MbwrDensityIteration {
            density,
            phase,
            pressure: q[0],
            pressure_density_derivative,
            relative_pressure_error: relative,
            pressure_converged: relative.abs() <= P_EPS,
            mechanically_stable: rho > 0. && q[1] >= 0.,
            path,
        })
    }
}
