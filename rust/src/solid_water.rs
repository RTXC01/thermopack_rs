//! Ice Ih Gibbs energy and its analytic derivatives from `h2o_gibbs.f90`.
//! Source reference shifts are owned per model. Values use SI molar units.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! ```
//! use thermopack_rs::SolidWater;
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let mut ice = SolidWater::new();
//! ice.initialize(100., -10000.)?;
//! let state = ice.prepare(250.)?.state(1e5)?;
//! assert!(state.volume > 0. && state.heat_capacity_p > 0.);
//! # Ok(())
//! # }
//! ```
use crate::{Error, Result};

const TT: f64 = 273.160;
const PT: f64 = 611.657;
const P0: f64 = 101325.;
const T0: f64 = 273.152519;
const MW: f64 = 18.01528e-3;
const G0: [f64; 5] = [
    -632020.233449497,
    0.655022213658955,
    -1.893699293261e-8,
    3.39746123271053e-15,
    -5.56464869058991e-22,
];

// Only the complex operations required by the two fixed analytic terms.
// Principal log uses atan2; reciprocal uses a scaled division to avoid norm^2
// overflow. No general-purpose complex dependency or heap storage is needed.
#[derive(Debug, Clone, Copy, Default)]
struct C {
    re: f64,
    im: f64,
}
impl C {
    const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    fn scale(self, a: f64) -> Self {
        Self::new(self.re * a, self.im * a)
    }
    fn shift(self, a: f64) -> Self {
        Self::new(self.re + a, self.im)
    }
    fn log(self) -> Self {
        Self::new(self.re.hypot(self.im).ln(), self.im.atan2(self.re))
    }
    fn divide_real(a: f64, b: Self) -> Self {
        if b.re.abs() >= b.im.abs() {
            let ratio = b.im / b.re;
            let den = b.re + b.im * ratio;
            Self::new(a / den, -a * ratio / den)
        } else {
            let ratio = b.re / b.im;
            let den = b.re * ratio + b.im;
            Self::new(a * ratio / den, -a / den)
        }
    }
}
impl std::ops::Add for C {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self::new(self.re + b.re, self.im + b.im)
    }
}
impl std::ops::Sub for C {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self::new(self.re - b.re, self.im - b.im)
    }
}
impl std::ops::Mul for C {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self::new(
            self.re * b.re - self.im * b.im,
            self.re * b.im + self.im * b.re,
        )
    }
}
const TERMS: [C; 2] = [
    C::new(3.68017112855051e-2, 5.10878114959572e-2),
    C::new(0.337315741065416, 0.335449415919309),
];
const R1: C = C::new(44.7050716285388, 65.6876847463481);
const R2: [C; 3] = [
    C::new(-72.597457432922, -78.100842711287),
    C::new(-5.57107698030123e-5, 4.64578634580806e-5),
    C::new(2.34801409215913e-11, -2.85651142904972e-11),
];

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolidWater {
    g00: f64,
    s0: f64,
}
impl Default for SolidWater {
    fn default() -> Self {
        Self {
            g00: G0[0],
            s0: -3327.33756492168,
        }
    }
}

/// Molar Gibbs energy (J/mol), T derivatives in K, P derivatives in Pa.
#[derive(Debug, Clone, Copy)]
pub struct SolidWaterGibbs {
    pub value: f64,
    pub dt: f64,
    pub dp: f64,
    pub dtt: f64,
    pub dtp: f64,
    pub dpp: f64,
}
#[derive(Debug, Clone, Copy)]
pub struct SolidWaterState {
    pub gibbs: SolidWaterGibbs,
    pub volume: f64,
    pub entropy: f64,
    pub enthalpy: f64,
    pub internal_energy: f64,
    pub helmholtz_energy: f64,
    /// Physical Cp = -T G_TT in J/(mol K).
    pub heat_capacity_p: f64,
    /// Original sh2o_heat_capacity multiplies the molar result by MW again.
    pub baseline_heat_capacity_p: f64,
}
#[derive(Debug, Clone, Copy)]
pub struct PreparedSolidWater<'a> {
    model: &'a SolidWater,
    t: f64,
    b: [C; 2],
    bt: [C; 2],
    btt: [C; 2],
}

impl SolidWater {
    pub const REFERENCE_TEMPERATURE: f64 = T0;
    pub const REFERENCE_PRESSURE: f64 = P0;
    pub const MELTING_ENTROPY: f64 = MW * 333.426516 * 1e3 / T0;
    pub fn new() -> Self {
        Self::default()
    }
    /// Original sho2_init order: entropy shift, then Gibbs shift at T0/P0.
    /// A rejected non-finite reference leaves the model unchanged.
    pub fn initialize(&mut self, liquid_entropy: f64, liquid_gibbs: f64) -> Result<()> {
        if !liquid_entropy.is_finite() || !liquid_gibbs.is_finite() {
            return Err(Error::InvalidInput("finite ice reference values required"));
        }
        let mut next = *self;
        let entropy = self.entropy(T0, P0)?;
        let target = liquid_entropy - Self::MELTING_ENTROPY;
        next.s0 += (target - entropy) / MW;
        let current = next.gibbs(T0, P0)?;
        next.g00 += (liquid_gibbs - current) / MW;
        if !next.s0.is_finite() || !next.g00.is_finite() {
            return Err(Error::InvalidInput("non-finite ice reference shift"));
        }
        *self = next;
        Ok(())
    }
    /// Prepare the two complex temperature functions once for pressure scans.
    /// Zero temperature is allowed; negative/non-finite T is rejected.
    pub fn prepare(&self, t: f64) -> Result<PreparedSolidWater<'_>> {
        if !t.is_finite() || t < 0. {
            return Err(Error::InvalidInput(
                "nonnegative finite ice temperature required",
            ));
        }
        let tau = t / TT;
        let mut b = [C::default(); 2];
        let mut bt = b;
        let mut btt = b;
        for (k, &a) in TERMS.iter().enumerate() {
            let minus = a.shift(-tau);
            let plus = a.shift(tau);
            let lm = minus.log();
            let lp = plus.log();
            b[k] = minus * lm + plus * lp - a.scale(2.) * a.log() - C::divide_real(tau * tau, a);
            bt[k] = lp - lm - C::divide_real(2. * tau, a);
            btt[k] = C::divide_real(1., minus) + C::divide_real(1., plus) - C::divide_real(2., a);
        }
        if b.iter()
            .chain(&bt)
            .chain(&btt)
            .any(|v| !v.re.is_finite() || !v.im.is_finite())
        {
            return Err(Error::InvalidInput("non-finite ice temperature functions"));
        }
        Ok(PreparedSolidWater {
            model: self,
            t,
            b,
            bt,
            btt,
        })
    }
    pub fn state(&self, t: f64, p: f64) -> Result<SolidWaterState> {
        self.prepare(t)?.state(p)
    }
    pub fn gibbs(&self, t: f64, p: f64) -> Result<f64> {
        Ok(self.prepare(t)?.gibbs(p)?.value)
    }
    pub fn dg_dt(&self, t: f64, p: f64) -> Result<f64> {
        Ok(self.prepare(t)?.gibbs(p)?.dt)
    }
    pub fn dg_dp(&self, t: f64, p: f64) -> Result<f64> {
        Ok(self.prepare(t)?.gibbs(p)?.dp)
    }
    pub fn d2g_dt2(&self, t: f64, p: f64) -> Result<f64> {
        Ok(self.prepare(t)?.gibbs(p)?.dtt)
    }
    pub fn d2g_dtdp(&self, t: f64, p: f64) -> Result<f64> {
        Ok(self.prepare(t)?.gibbs(p)?.dtp)
    }
    pub fn d2g_dp2(&self, t: f64, p: f64) -> Result<f64> {
        Ok(self.prepare(t)?.gibbs(p)?.dpp)
    }
    pub fn specific_volume(&self, t: f64, p: f64) -> Result<f64> {
        self.dg_dp(t, p)
    }
    pub fn entropy(&self, t: f64, p: f64) -> Result<f64> {
        Ok(-self.dg_dt(t, p)?)
    }
    pub fn enthalpy(&self, t: f64, p: f64) -> Result<f64> {
        Ok(self.state(t, p)?.enthalpy)
    }
    pub fn internal_energy(&self, t: f64, p: f64) -> Result<f64> {
        Ok(self.state(t, p)?.internal_energy)
    }
    pub fn helmholtz(&self, t: f64, p: f64) -> Result<f64> {
        Ok(self.state(t, p)?.helmholtz_energy)
    }
    pub fn heat_capacity(&self, t: f64, p: f64) -> Result<f64> {
        Ok(-t * self.d2g_dt2(t, p)?)
    }
    pub fn baseline_heat_capacity(&self, t: f64, p: f64) -> Result<f64> {
        Ok(self.heat_capacity(t, p)? * MW)
    }
}
impl PreparedSolidWater<'_> {
    pub fn gibbs(&self, p: f64) -> Result<SolidWaterGibbs> {
        if !p.is_finite() {
            return Err(Error::InvalidInput("finite ice pressure required"));
        }
        // Keep the original reduced-pressure subtraction and summation order.
        let x = p / PT - P0 / PT;
        let powers = [1., x, x * x, x * (x * x), (x * x) * (x * x)];
        let mut g0 = 0.;
        let mut gp = 0.;
        let mut gpp = 0.;
        let mut coefficients = G0;
        coefficients[0] = self.model.g00;
        for (k, &c) in coefficients.iter().enumerate() {
            g0 += c * powers[k];
            if k >= 1 {
                gp += c * k as f64 * powers[k - 1] / PT;
            }
            if k >= 2 {
                gpp += c * k as f64 * (k - 1) as f64 * powers[k - 2] / (PT * PT);
            }
        }
        let r2 = R2[0] + R2[1].scale(x) + R2[2].scale(powers[2]);
        let rp = R2[1].scale(1. / PT) + R2[2].scale(2. * x / PT);
        let rpp = R2[2].scale(2. / (PT * PT));
        let g = R1 * self.b[0] + r2 * self.b[1];
        let gt = R1 * self.bt[0] + r2 * self.bt[1];
        let gtt = R1 * self.btt[0] + r2 * self.btt[1];
        let result = SolidWaterGibbs {
            value: (g0 - self.model.s0 * self.t + TT * g.re) * MW,
            dt: (-self.model.s0 + gt.re) * MW,
            dtt: gtt.re / TT * MW,
            dp: (gp + TT * (rp * self.b[1]).re) * MW,
            dtp: (rp * self.bt[1]).re * MW,
            dpp: (gpp + TT * (rpp * self.b[1]).re) * MW,
        };
        if [
            result.value,
            result.dt,
            result.dp,
            result.dtt,
            result.dtp,
            result.dpp,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(Error::InvalidInput("non-finite ice Gibbs derivatives"));
        }
        Ok(result)
    }
    pub fn state(&self, p: f64) -> Result<SolidWaterState> {
        let g = self.gibbs(p)?;
        let entropy = -g.dt;
        let h = g.value + self.t * entropy;
        let cp = -self.t * g.dtt;
        let state = SolidWaterState {
            gibbs: g,
            volume: g.dp,
            entropy,
            enthalpy: h,
            internal_energy: h - p * g.dp,
            helmholtz_energy: g.value - p * g.dp,
            heat_capacity_p: cp,
            baseline_heat_capacity_p: cp * MW,
        };
        if [
            h,
            state.internal_energy,
            state.helmholtz_energy,
            cp,
            state.baseline_heat_capacity_p,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(Error::InvalidInput("non-finite ice state"));
        }
        Ok(state)
    }
}
