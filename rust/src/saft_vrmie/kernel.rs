//! Directional evaluation of the non-associating residual terms.
//! All extensive coordinates are seeded at unit total moles before scaling.
use super::*;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

const STACK_NC: usize = 8;

enum SmallBuf<T: Copy> {
    Stack([T; STACK_NC], usize),
    Heap(Vec<T>),
}

impl<T: Copy> SmallBuf<T> {
    #[inline]
    fn from_elem(val: T, len: usize) -> Self {
        if len <= STACK_NC {
            Self::Stack([val; STACK_NC], len)
        } else {
            Self::Heap(vec![val; len])
        }
    }

    #[inline]
    fn as_slice(&self) -> &[T] {
        match self {
            Self::Stack(arr, len) => &arr[..*len],
            Self::Heap(vec) => vec.as_slice(),
        }
    }

    #[inline]
    fn as_mut_slice(&mut self) -> &mut [T] {
        match self {
            Self::Stack(arr, len) => &mut arr[..*len],
            Self::Heap(vec) => vec.as_mut_slice(),
        }
    }
}

enum SmallMat<T: Copy> {
    Stack([T; STACK_NC * STACK_NC], usize),
    Heap(Vec<T>),
}

impl<T: Copy> SmallMat<T> {
    #[inline]
    fn from_elem(val: T, len: usize) -> Self {
        if len <= STACK_NC * STACK_NC {
            Self::Stack([val; STACK_NC * STACK_NC], len)
        } else {
            Self::Heap(vec![val; len])
        }
    }

    #[inline]
    fn as_slice(&self) -> &[T] {
        match self {
            Self::Stack(arr, len) => &arr[..*len],
            Self::Heap(vec) => vec.as_slice(),
        }
    }

    #[inline]
    fn as_mut_slice(&mut self) -> &mut [T] {
        match self {
            Self::Stack(arr, len) => &mut arr[..*len],
            Self::Heap(vec) => vec.as_mut_slice(),
        }
    }
}

/// First-order dual number (value + one directional derivative).
#[derive(Debug, Clone, Copy)]
pub(super) struct Dual1 {
    pub f0: f64,
    pub f1: f64,
}

impl Add for Dual1 {
    type Output = Self;
    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        Self {
            f0: self.f0 + rhs.f0,
            f1: self.f1 + rhs.f1,
        }
    }
}

impl AddAssign for Dual1 {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) {
        self.f0 += rhs.f0;
        self.f1 += rhs.f1;
    }
}

impl Sub for Dual1 {
    type Output = Self;
    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        Self {
            f0: self.f0 - rhs.f0,
            f1: self.f1 - rhs.f1,
        }
    }
}

impl SubAssign for Dual1 {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Self) {
        self.f0 -= rhs.f0;
        self.f1 -= rhs.f1;
    }
}

impl Neg for Dual1 {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        Self {
            f0: -self.f0,
            f1: -self.f1,
        }
    }
}

impl Mul for Dual1 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: Self) -> Self {
        Self {
            f0: self.f0 * rhs.f0,
            f1: self.f0 * rhs.f1 + self.f1 * rhs.f0,
        }
    }
}

impl Div for Dual1 {
    type Output = Self;
    #[inline(always)]
    fn div(self, rhs: Self) -> Self {
        let inv = 1.0 / rhs.f0;
        let q = self.f0 * inv;
        Self {
            f0: q,
            f1: (self.f1 - q * rhs.f1) * inv,
        }
    }
}

/// Second-order hyperdual number (value + two first derivatives + mixed second derivative).
#[derive(Debug, Clone, Copy)]
pub(super) struct Dual2 {
    pub f0: f64,
    pub f1: f64,
    pub f2: f64,
    pub f12: f64,
}

impl Add for Dual2 {
    type Output = Self;
    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        Self {
            f0: self.f0 + rhs.f0,
            f1: self.f1 + rhs.f1,
            f2: self.f2 + rhs.f2,
            f12: self.f12 + rhs.f12,
        }
    }
}

impl AddAssign for Dual2 {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) {
        self.f0 += rhs.f0;
        self.f1 += rhs.f1;
        self.f2 += rhs.f2;
        self.f12 += rhs.f12;
    }
}

impl Sub for Dual2 {
    type Output = Self;
    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        Self {
            f0: self.f0 - rhs.f0,
            f1: self.f1 - rhs.f1,
            f2: self.f2 - rhs.f2,
            f12: self.f12 - rhs.f12,
        }
    }
}

impl SubAssign for Dual2 {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Self) {
        self.f0 -= rhs.f0;
        self.f1 -= rhs.f1;
        self.f2 -= rhs.f2;
        self.f12 -= rhs.f12;
    }
}

impl Neg for Dual2 {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        Self {
            f0: -self.f0,
            f1: -self.f1,
            f2: -self.f2,
            f12: -self.f12,
        }
    }
}

impl Mul for Dual2 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: Self) -> Self {
        Self {
            f0: self.f0 * rhs.f0,
            f1: self.f0 * rhs.f1 + self.f1 * rhs.f0,
            f2: self.f0 * rhs.f2 + self.f2 * rhs.f0,
            f12: self.f0 * rhs.f12 + self.f1 * rhs.f2 + self.f2 * rhs.f1 + self.f12 * rhs.f0,
        }
    }
}

impl Div for Dual2 {
    type Output = Self;
    #[inline(always)]
    fn div(self, rhs: Self) -> Self {
        let inv = 1.0 / rhs.f0;
        let q0 = self.f0 * inv;
        let q1 = (self.f1 - q0 * rhs.f1) * inv;
        let q2 = (self.f2 - q0 * rhs.f2) * inv;
        let q12 = (self.f12 - q1 * rhs.f2 - q2 * rhs.f1 - q0 * rhs.f12) * inv;
        Self {
            f0: q0,
            f1: q1,
            f2: q2,
            f12: q12,
        }
    }
}

pub(super) trait Scalar:
    Copy
    + Add<Output = Self>
    + AddAssign
    + Sub<Output = Self>
    + SubAssign
    + Neg<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
{
    fn real(v: f64) -> Self;
    fn seed(v: f64, d1: bool, d2: bool) -> Self;
    fn val(self) -> f64;
    fn add_f(self, rhs: f64) -> Self;
    fn mul_f(self, rhs: f64) -> Self;
    fn sub_from_f(lhs: f64, rhs: Self) -> Self;
    fn div_from_f(lhs: f64, rhs: Self) -> Self;
    fn pow2(self) -> Self;
    fn pow3(self) -> Self;
    fn pow4(self) -> Self;
    fn pow7(self) -> Self;
    fn pow_real(self, exp: f64) -> Self;
    fn exp(self) -> Self;
    fn ln(self) -> Self;
}

impl Scalar for f64 {
    #[inline(always)]
    fn real(v: f64) -> Self {
        v
    }
    #[inline(always)]
    fn seed(v: f64, _d1: bool, _d2: bool) -> Self {
        v
    }
    #[inline(always)]
    fn val(self) -> f64 {
        self
    }
    #[inline(always)]
    fn add_f(self, rhs: f64) -> Self {
        self + rhs
    }
    #[inline(always)]
    fn mul_f(self, rhs: f64) -> Self {
        self * rhs
    }
    #[inline(always)]
    fn sub_from_f(lhs: f64, rhs: Self) -> Self {
        lhs - rhs
    }
    #[inline(always)]
    fn div_from_f(lhs: f64, rhs: Self) -> Self {
        lhs / rhs
    }
    #[inline(always)]
    fn pow2(self) -> Self {
        self * self
    }
    #[inline(always)]
    fn pow3(self) -> Self {
        self * self * self
    }
    #[inline(always)]
    fn pow4(self) -> Self {
        let s2 = self * self;
        s2 * s2
    }
    #[inline(always)]
    fn pow7(self) -> Self {
        let s2 = self * self;
        let s4 = s2 * s2;
        s4 * s2 * self
    }
    #[inline(always)]
    fn pow_real(self, exp: f64) -> Self {
        self.powf(exp)
    }
    #[inline(always)]
    fn exp(self) -> Self {
        f64::exp(self)
    }
    #[inline(always)]
    fn ln(self) -> Self {
        f64::ln(self)
    }
}

impl Scalar for Dual1 {
    #[inline(always)]
    fn real(v: f64) -> Self {
        Self { f0: v, f1: 0.0 }
    }
    #[inline(always)]
    fn seed(v: f64, d1: bool, _d2: bool) -> Self {
        Self {
            f0: v,
            f1: f64::from(d1),
        }
    }
    #[inline(always)]
    fn val(self) -> f64 {
        self.f0
    }
    #[inline(always)]
    fn add_f(self, rhs: f64) -> Self {
        Self {
            f0: self.f0 + rhs,
            f1: self.f1,
        }
    }
    #[inline(always)]
    fn mul_f(self, rhs: f64) -> Self {
        Self {
            f0: self.f0 * rhs,
            f1: self.f1 * rhs,
        }
    }
    #[inline(always)]
    fn sub_from_f(lhs: f64, rhs: Self) -> Self {
        Self {
            f0: lhs - rhs.f0,
            f1: -rhs.f1,
        }
    }
    #[inline(always)]
    fn div_from_f(lhs: f64, rhs: Self) -> Self {
        let inv = 1.0 / rhs.f0;
        let q = lhs * inv;
        Self {
            f0: q,
            f1: -q * rhs.f1 * inv,
        }
    }
    #[inline(always)]
    fn pow2(self) -> Self {
        self * self
    }
    #[inline(always)]
    fn pow3(self) -> Self {
        self * self * self
    }
    #[inline(always)]
    fn pow4(self) -> Self {
        let s2 = self * self;
        s2 * s2
    }
    #[inline(always)]
    fn pow7(self) -> Self {
        let s2 = self * self;
        let s4 = s2 * s2;
        s4 * s2 * self
    }
    #[inline(always)]
    fn pow_real(self, exp: f64) -> Self {
        let p0 = self.f0.powf(exp);
        let d1 = if self.f0 != 0.0 {
            exp * p0 / self.f0
        } else {
            0.0
        };
        Self {
            f0: p0,
            f1: self.f1 * d1,
        }
    }
    #[inline(always)]
    fn exp(self) -> Self {
        let v = self.f0.exp();
        Self {
            f0: v,
            f1: v * self.f1,
        }
    }
    #[inline(always)]
    fn ln(self) -> Self {
        Self {
            f0: self.f0.ln(),
            f1: self.f1 / self.f0,
        }
    }
}

impl Scalar for Dual2 {
    #[inline(always)]
    fn real(v: f64) -> Self {
        Self {
            f0: v,
            f1: 0.0,
            f2: 0.0,
            f12: 0.0,
        }
    }
    #[inline(always)]
    fn seed(v: f64, d1: bool, d2: bool) -> Self {
        Self {
            f0: v,
            f1: f64::from(d1),
            f2: f64::from(d2),
            f12: 0.0,
        }
    }
    #[inline(always)]
    fn val(self) -> f64 {
        self.f0
    }
    #[inline(always)]
    fn add_f(self, rhs: f64) -> Self {
        Self {
            f0: self.f0 + rhs,
            ..self
        }
    }
    #[inline(always)]
    fn mul_f(self, rhs: f64) -> Self {
        Self {
            f0: self.f0 * rhs,
            f1: self.f1 * rhs,
            f2: self.f2 * rhs,
            f12: self.f12 * rhs,
        }
    }
    #[inline(always)]
    fn sub_from_f(lhs: f64, rhs: Self) -> Self {
        Self {
            f0: lhs - rhs.f0,
            f1: -rhs.f1,
            f2: -rhs.f2,
            f12: -rhs.f12,
        }
    }
    #[inline(always)]
    fn div_from_f(lhs: f64, rhs: Self) -> Self {
        let inv = 1.0 / rhs.f0;
        let q0 = lhs * inv;
        let q1 = -q0 * rhs.f1 * inv;
        let q2 = -q0 * rhs.f2 * inv;
        let q12 = -(q1 * rhs.f2 + q2 * rhs.f1 + q0 * rhs.f12) * inv;
        Self {
            f0: q0,
            f1: q1,
            f2: q2,
            f12: q12,
        }
    }
    #[inline(always)]
    fn pow2(self) -> Self {
        self * self
    }
    #[inline(always)]
    fn pow3(self) -> Self {
        self * self * self
    }
    #[inline(always)]
    fn pow4(self) -> Self {
        let s2 = self * self;
        s2 * s2
    }
    #[inline(always)]
    fn pow7(self) -> Self {
        let s2 = self * self;
        let s4 = s2 * s2;
        s4 * s2 * self
    }
    #[inline(always)]
    fn pow_real(self, exp: f64) -> Self {
        let p0 = self.f0.powf(exp);
        if self.f0 != 0.0 {
            let inv = 1.0 / self.f0;
            let d1 = exp * p0 * inv;
            let d2 = (exp - 1.0) * d1 * inv;
            Self {
                f0: p0,
                f1: self.f1 * d1,
                f2: self.f2 * d1,
                f12: self.f12 * d1 + self.f1 * self.f2 * d2,
            }
        } else {
            Self::real(p0)
        }
    }
    #[inline(always)]
    fn exp(self) -> Self {
        let v = self.f0.exp();
        Self {
            f0: v,
            f1: v * self.f1,
            f2: v * self.f2,
            f12: v * (self.f12 + self.f1 * self.f2),
        }
    }
    #[inline(always)]
    fn ln(self) -> Self {
        let inv = 1.0 / self.f0;
        let f1 = self.f1 * inv;
        let f2 = self.f2 * inv;
        Self {
            f0: self.f0.ln(),
            f1,
            f2,
            f12: self.f12 * inv - f1 * f2,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum EvalMode {
    Total,
    Dispersion,
    HardSphere,
}

pub(super) struct Terms<T> {
    pub hs: T,
    pub a1: T,
    pub a2: T,
    pub a3: T,
    pub chain: T,
    pub total: T,
    pub dispersion: T,
    pub hard_sphere: T,
}

impl PreparedSaftVrMie<'_> {
    #[inline]
    pub(super) fn eval_scalar(&self, rho_value: f64) -> Result<Terms<f64>> {
        self.kernel_typed::<f64>(rho_value, None, None, EvalMode::Total)
    }

    #[inline]
    pub(super) fn eval_first(&self, rho_value: f64, dir1: usize) -> Result<Terms<Dual1>> {
        self.kernel_typed::<Dual1>(rho_value, Some(dir1), None, EvalMode::Total)
    }

    #[inline]
    pub(super) fn eval_second(
        &self,
        rho_value: f64,
        dir1: usize,
        dir2: usize,
        mode: EvalMode,
    ) -> Result<Terms<Dual2>> {
        self.kernel_typed::<Dual2>(rho_value, Some(dir1), Some(dir2), mode)
    }

    fn kernel_typed<T: Scalar>(
        &self,
        rho_value: f64,
        dir1: Option<usize>,
        dir2: Option<usize>,
        mode: EvalMode,
    ) -> Result<Terms<T>> {
        let nc = self.model.num_components();
        let seed = |value, index| T::seed(value, dir1 == Some(index), dir2 == Some(index));
        let zero = T::real(0.0);
        let t = seed(self.temperature, 0);
        let v = seed(1.0 / rho_value, 1);

        let has_mole_dir = dir1.is_some_and(|d| d >= 2) || dir2.is_some_and(|d| d >= 2);
        let mut x_buf = SmallBuf::from_elem(zero, nc);
        let mut xs_buf = SmallBuf::from_elem(zero, nc);
        let (sum, segments) = if has_mole_dir {
            let mut sum = zero;
            for (i, &n) in self.mole_fractions.iter().enumerate() {
                let ni = seed(n, 2 + i);
                x_buf.as_mut_slice()[i] = ni;
                sum += ni;
            }
            let mut segments = zero;
            for (x_slot, comp) in x_buf.as_mut_slice().iter_mut().zip(&self.model.components) {
                *x_slot = *x_slot / sum;
                segments += x_slot.mul_f(comp.m);
            }
            for ((xs_slot, &xi), comp) in xs_buf
                .as_mut_slice()
                .iter_mut()
                .zip(x_buf.as_slice())
                .zip(&self.model.components)
            {
                *xs_slot = xi.mul_f(comp.m) / segments;
            }
            (sum, segments)
        } else {
            for (x_slot, &xi) in x_buf.as_mut_slice().iter_mut().zip(&self.mole_fractions) {
                *x_slot = T::real(xi);
            }
            for (xs_slot, &xsi) in xs_buf.as_mut_slice().iter_mut().zip(&self.xs) {
                *xs_slot = T::real(xsi);
            }
            (T::real(1.0), T::real(self.m_bar))
        };
        let x = x_buf.as_slice();
        let xs = xs_buf.as_slice();
        let rho = sum / v;

        let has_temp_dir = dir1 == Some(0) || dir2 == Some(0);
        let dt = t.add_f(-self.temperature);
        let mut d_buf = SmallMat::from_elem(zero, nc * nc);
        if has_temp_dir {
            for (k, d_slot) in d_buf.as_mut_slice().iter_mut().enumerate() {
                *d_slot = (dt.mul_f(0.5 * self.d_tt[k]).add_f(self.d_t[k]) * dt).add_f(self.d[k]);
            }
        } else {
            for (d_slot, &dk) in d_buf.as_mut_slice().iter_mut().zip(&self.d) {
                *d_slot = T::real(dk);
            }
        }
        let d = d_buf.as_slice();
        let lift = |p: [f64; 3]| {
            if has_temp_dir {
                (dt.mul_f(0.5 * p[2]).add_f(p[1]) * dt).add_f(p[0])
            } else {
                T::real(p[0])
            }
        };

        let single = self.model.reference == HardSphereReference::SingleComponent;
        let options = self.model.options;
        let terms_cfg = self.model.terms;
        let packing_prefactor = (rho * segments).mul_f(PI * N_A / 6.0);

        let (eta, zeta_bar, mean_d) = if !has_temp_dir && !has_mole_dir {
            let eta = rho.mul_f(self.conv_eta);
            let zeta_bar = rho.mul_f(self.conv_zeta);
            let mean_d = if single {
                let mut md = 0.0;
                for i in 0..nc {
                    for j in 0..nc {
                        md += self.xs[i] * self.xs[j] * self.d[i * nc + j];
                    }
                }
                T::real(md)
            } else {
                zero
            };
            (eta, zeta_bar, mean_d)
        } else {
            let mut mean_d = zero;
            let mut mean_sigma = zero;
            let mut eta_acc = zero;
            let mut zeta_acc = zero;
            for i in 0..nc {
                for j in i..nc {
                    let k = i * nc + j;
                    let weight = if i == j {
                        xs[i] * xs[i]
                    } else {
                        (xs[i] * xs[j]).mul_f(2.0)
                    };
                    let sigma = lift(self.effective[k].0);
                    match options.packing {
                        PackingRule::Lafitte => {
                            eta_acc += weight * d[k].pow3();
                            zeta_acc += weight * sigma.pow3();
                        }
                        PackingRule::Linear if i == j => {
                            eta_acc += xs[i] * d[k].pow3();
                            zeta_acc += xs[i] * sigma.pow3();
                        }
                        PackingRule::Leonard => {
                            mean_d += weight * d[k];
                            mean_sigma += weight * sigma;
                        }
                        _ => {}
                    }
                }
            }
            if options.packing == PackingRule::Leonard {
                eta_acc = mean_d.pow3();
                zeta_acc = mean_sigma.pow3();
            }
            (
                eta_acc * packing_prefactor,
                zeta_acc * packing_prefactor,
                mean_d,
            )
        };

        if !(0.0..1.0).contains(&eta.val()) {
            return Err(Error::InvalidInput(
                "density exceeds maximum packing fraction",
            ));
        }

        let need_hs = match mode {
            EvalMode::Total => terms_cfg.hard_sphere || options.extra_hard_sphere,
            EvalMode::HardSphere => terms_cfg.hard_sphere || options.extra_hard_sphere,
            EvalMode::Dispersion => false,
        };
        let need_moments =
            need_hs || (mode != EvalMode::HardSphere && options.bmcsl_compressibility);
        let (z3, om3, k_mix, l_mix, m2) = if need_moments {
            let mut m1 = zero;
            let mut m2 = zero;
            let mut m3 = zero;
            for (i, &weight) in xs.iter().enumerate() {
                let di = d[i * nc + i].mul_f(1e10);
                let di2 = di * di;
                m1 += weight * di;
                m2 += weight * di2;
                m3 += weight * (di2 * di);
            }
            let z3 = (packing_prefactor * m3).mul_f(1e-30);
            let om3 = T::sub_from_f(1.0, z3);
            let k_mix = (m1 * m2) / m3;
            let l_mix = m2.pow3() / m3.pow2();
            (z3, om3, k_mix, l_mix, m2)
        } else {
            (zero, T::real(1.0), zero, zero, zero)
        };

        // 1. Hard-sphere term
        let a_hs = if need_hs && terms_cfg.hard_sphere {
            if self.model.reference == HardSphereReference::NonAdditive {
                santos(x, d, rho, segments, options.kolafa_pure_reference)
            } else if single {
                pure_hs(eta, options.kolafa_pure_reference)
            } else if nc == 1 {
                pure_hs(z3, false)
            } else {
                (k_mix * z3).mul_f(3.0) / om3
                    + (l_mix * z3) / om3.pow2()
                    + l_mix.add_f(-1.0) * om3.ln()
            }
        } else {
            zero
        };

        let mut extra = zero;
        if need_hs && options.extra_hard_sphere && nc > 1 {
            let z2 = (packing_prefactor * m2).mul_f(1e-20);
            let inv_om3 = T::div_from_f(1.0, om3);
            for i in 0..nc {
                let di = d[i * nc + i];
                for j in i + 1..nc {
                    let dj = d[j * nc + j];
                    let additive = (di + dj).mul_f(0.5);
                    let mu = (di * dj) / (di + dj);
                    let h = z2 * mu * inv_om3;
                    let g = (h.mul_f(3.0) + h.pow2().mul_f(2.0)).add_f(1.0) * inv_om3;
                    extra += (rho
                        * x[i]
                        * x[j]
                        * g
                        * additive.pow2()
                        * (lift(self.exact_d[i * nc + j]) - additive))
                        .mul_f(4.0 * PI * N_A);
                }
            }
        }

        let hard_sphere = segments * a_hs + extra;
        if mode == EvalMode::HardSphere {
            return Ok(Terms {
                hs: a_hs,
                a1: zero,
                a2: zero,
                a3: zero,
                chain: zero,
                dispersion: zero,
                hard_sphere,
                total: sum * hard_sphere,
            });
        }

        let beta = T::div_from_f(1.0, t);

        // 2. Dispersion A1, A2, A3
        let mut a1_mix = zero;
        let mut a2_mix = zero;
        let mut a3_mix = zero;

        let eta2 = eta * eta;
        let eta3 = eta2 * eta;
        let eta4 = eta2 * eta2;
        let om_eta = T::sub_from_f(1.0, eta);
        let om_eta2 = om_eta * om_eta;
        let denum_khs = (eta.mul_f(4.0) + eta2.mul_f(4.0) - eta3.mul_f(4.0) + eta4).add_f(1.0);
        let khs = if options.bmcsl_compressibility {
            om3.pow4()
                / (z3
                    * (k_mix.mul_f(6.0).add_f(-2.0)
                        + z3 * (l_mix.mul_f(9.0) - k_mix.mul_f(6.0)).add_f(1.0)
                        + z3.pow2() * l_mix * z3.add_f(-4.0)))
                .add_f(1.0)
        } else {
            (om_eta2 * om_eta2) / denum_khs
        };

        let need_disp = terms_cfg.first || terms_cfg.second || terms_cfg.third;
        if need_disp {
            let zb2 = zeta_bar * zeta_bar;
            let zb4 = zb2 * zb2;
            let zb7 = zeta_bar.pow7();
            for i in 0..nc {
                for j in i..nc {
                    let idx = i * nc + j;
                    let la = self.model.lambda_a_ij[idx];
                    let lr = self.model.lambda_r_ij[idx];
                    let eps = self.model.eps_divk_ij[idx];
                    let c_fac = self.model.c_ij[idx];
                    let dij = if single { mean_d } else { d[idx] };
                    let sigma = self.model.sigma_ij[idx];
                    let x0 = T::div_from_f(sigma, dij);
                    let sigma_eff = lift(self.effective[idx].0);
                    let eps_eff = lift(self.effective[idx].1);
                    let y = sigma_eff / dij;
                    let y_gt_1 = y.val() > 1.0;
                    let integrate = |lambda| {
                        let b = if y_gt_1 {
                            calc_b_tilde(y, eta, lambda, eps)
                        } else {
                            zero
                        };
                        x0.pow_real(lambda) * (calc_a1_sutherland(eta, lambda, eps) + b)
                    };
                    let order = self.model.components[i].fh_order as usize;
                    let dfh = beta.mul_f(self.model.d_feynhibbs_param_ij[idx] / (sigma * sigma));
                    let qa = quantum_factors(la, dfh, options.quantum_variant);
                    let qr = quantum_factors(lr, dfh, options.quantum_variant);
                    let mut a1_tilde = zero;
                    let mut a22_tilde = zero;
                    let mut alpha = zero;
                    let ratio = T::div_from_f(sigma, sigma_eff);
                    let need_alpha = order > 0 && (terms_cfg.second || terms_cfg.third);
                    for k in 0..=order {
                        let a = la + 2.0 * k as f64;
                        let r = lr + 2.0 * k as f64;
                        if terms_cfg.first {
                            a1_tilde += (qa[k] * integrate(a) - qr[k] * integrate(r)).mul_f(c_fac);
                        }
                        if need_alpha {
                            alpha += qa[k] * ratio.pow_real(a).mul_f(1.0 / (a - 3.0))
                                - qr[k] * ratio.pow_real(r).mul_f(1.0 / (r - 3.0));
                        }
                        if terms_cfg.second && (options.quantum_correct_a2 || k == 0) {
                            let a2_order = if options.quantum_correct_a2 { order } else { 0 };
                            for l in 0..=a2_order {
                                let offset = 2.0 * l as f64;
                                a22_tilde += qa[k] * qa[l] * integrate(a + la + offset)
                                    + qr[k] * qr[l] * integrate(r + lr + offset)
                                    - (qa[k] * qr[l] * integrate(a + lr + offset)).mul_f(2.0);
                            }
                        }
                    }

                    let weight = if i == j {
                        xs[i] * xs[i]
                    } else {
                        (xs[i] * xs[j]).mul_f(2.0)
                    };
                    let eta_ij = if terms_cfg.first || terms_cfg.second {
                        packing_prefactor * dij.pow3()
                    } else {
                        zero
                    };
                    if terms_cfg.first {
                        a1_mix += weight * (eta_ij * a1_tilde);
                    }
                    if terms_cfg.second || terms_cfg.third {
                        let f_alpha = if order == 0 {
                            self.model.f_alpha_ij[idx].map(T::real)
                        } else {
                            alpha_functions(alpha.mul_f(c_fac * eps) / eps_eff)
                        };
                        if terms_cfg.second {
                            let chi = zeta_bar * (f_alpha[0] + f_alpha[1] * zb4 + f_alpha[2] * zb7);
                            let a2_ij = (khs * chi.add_f(1.0) * eta_ij * a22_tilde)
                                .mul_f(0.5 * eps * c_fac * c_fac);
                            a2_mix += weight * a2_ij;
                        }
                        if terms_cfg.third {
                            let a3_packing = if options.pair_prefactor_a3 {
                                packing_prefactor * sigma_eff.pow3()
                            } else {
                                zeta_bar
                            };
                            let a3_ij = -(eps_eff.pow3()
                                * f_alpha[3]
                                * a3_packing
                                * (f_alpha[4] * zeta_bar + f_alpha[5] * zb2).exp());
                            a3_mix += weight * a3_ij;
                        }
                    }
                }
            }
        }

        let mut dispersion = segments * beta * (a1_mix + beta * (a2_mix + beta * a3_mix));
        if let Some((rc, shift)) = self.model.cutoff {
            let la = self.model.lambda_a_ij[0];
            let lr = self.model.lambda_r_ij[0];
            let sigma = self.model.sigma_ij[0];
            let dfh = beta.mul_f(self.model.d_feynhibbs_param_ij[0] / (sigma * sigma));
            let order = self.model.components[0].fh_order as usize;
            let tail = |lambda: f64| {
                let qf = quantum_factors(lambda, dfh, options.quantum_variant);
                let mut acc = zero;
                for (k, &q) in qf.iter().take(1 + order).enumerate() {
                    let exponent = lambda + 2.0 * k as f64 - 3.0;
                    let factor =
                        rc.powf(-exponent) * (1.0 / exponent + if shift { 1.0 / 3.0 } else { 0.0 });
                    acc += q.mul_f(factor);
                }
                acc
            };
            let pref =
                2.0 * PI * N_A * self.model.c_ij[0] * sigma.powi(3) * self.model.eps_divk_ij[0];
            dispersion -= (rho * beta * (tail(lr) - tail(la))).mul_f(pref);
        }

        if mode == EvalMode::Dispersion {
            return Ok(Terms {
                hs: zero,
                a1: a1_mix,
                a2: a2_mix,
                a3: a3_mix,
                chain: zero,
                dispersion,
                hard_sphere: zero,
                total: sum * dispersion,
            });
        }

        // 3. Chain term
        let mut a_chain = zero;
        if terms_cfg.chain {
            let d2 = om_eta2;
            let d3 = d2 * om_eta;
            let k0 = -om_eta.ln()
                + (eta.mul_f(42.0) - eta2.mul_f(39.0) + eta3.mul_f(9.0) - eta4.mul_f(2.0))
                    / d3.mul_f(6.0);
            let k1 = (eta4 + eta2.mul_f(6.0) - eta.mul_f(12.0)) / d3.mul_f(2.0);
            let k2 = eta2.mul_f(-3.0) / d2.mul_f(8.0);
            let k3 = (-eta4 + eta2.mul_f(3.0) + eta.mul_f(3.0)) / d3.mul_f(6.0);
            let num_khs = d2 * d2;
            let num_khs_d = d3.mul_f(-4.0);
            let den_khs_d = (eta.mul_f(8.0) - eta2.mul_f(12.0) + eta3.mul_f(4.0)).add_f(4.0);
            let khs_e = (num_khs_d * denum_khs - num_khs * den_khs_d) / denum_khs.pow2();
            let gam_exp = (zeta_bar.mul_f(-6.7) - (zeta_bar * zeta_bar).mul_f(8.0)).exp();

            for i in 0..nc {
                let mi = self.model.components[i].m;
                if (mi - 1.0).abs() > 1e-12 {
                    let idx = i * nc + i;
                    let la = self.model.lambda_a_ij[idx];
                    let lr = self.model.lambda_r_ij[idx];
                    let eps = self.model.eps_divk_ij[idx];
                    let c_fac = self.model.c_ij[idx];
                    let x0 = T::div_from_f(
                        self.model.sigma_ij[idx],
                        if single { mean_d } else { d[idx] },
                    );
                    let alpha = self.model.alpha_ij[idx];
                    let x0_gt_1 = x0.val() > 1.0;

                    // gd
                    let x0_2 = x0 * x0;
                    let exparg = k0 + k1 * x0 + k2 * x0_2 + k3 * (x0_2 * x0);
                    let gd = exparg.exp();

                    // g1
                    let g1 = if terms_cfg.first {
                        let (asr, asr_e) = calc_a1_sutherland_with_deriv(eta, lr, eps);
                        let (asa, asa_e) = calc_a1_sutherland_with_deriv(eta, la, eps);
                        let (br, br_e, ba, ba_e) = if x0_gt_1 {
                            let (br, br_e) = calc_b_tilde_with_deriv(x0, eta, lr, eps);
                            let (ba, ba_e) = calc_b_tilde_with_deriv(x0, eta, la, eps);
                            (br, br_e, ba, ba_e)
                        } else {
                            (zero, zero, zero, zero)
                        };
                        let x0_la = x0.pow_real(la);
                        let x0_lr = x0.pow_real(lr);
                        let term_a = x0_la * (asa + ba);
                        let term_r = x0_lr * (asr + br);
                        let a1_tilde = (term_a - term_r).mul_f(c_fac);
                        let a1_tilde_e =
                            (x0_la * (asa_e + ba_e) - x0_lr * (asr_e + br_e)).mul_f(c_fac);
                        let g11 = (a1_tilde + eta * a1_tilde_e).mul_f(1.0 / (4.0 * eps));
                        let g12 = (term_r.mul_f(lr) - term_a.mul_f(la)).mul_f(c_fac / (12.0 * eps));
                        g11 + g12
                    } else {
                        zero
                    };

                    // g2
                    let g2 = if terms_cfg.second {
                        let lam_2r = 2.0 * lr;
                        let lam_2a = 2.0 * la;
                        let lam_ar = la + lr;
                        let (as_2r, as_2r_e) = calc_a1_sutherland_with_deriv(eta, lam_2r, eps);
                        let (as_2a, as_2a_e) = calc_a1_sutherland_with_deriv(eta, lam_2a, eps);
                        let (as_ar, as_ar_e) = calc_a1_sutherland_with_deriv(eta, lam_ar, eps);
                        let (b_2r, b_2r_e, b_2a, b_2a_e, b_ar, b_ar_e) = if x0_gt_1 {
                            let (b_2r, b_2r_e) = calc_b_tilde_with_deriv(x0, eta, lam_2r, eps);
                            let (b_2a, b_2a_e) = calc_b_tilde_with_deriv(x0, eta, lam_2a, eps);
                            let (b_ar, b_ar_e) = calc_b_tilde_with_deriv(x0, eta, lam_ar, eps);
                            (b_2r, b_2r_e, b_2a, b_2a_e, b_ar, b_ar_e)
                        } else {
                            (zero, zero, zero, zero, zero, zero)
                        };
                        let x0_2a = x0.pow_real(lam_2a);
                        let x0_ar = x0.pow_real(lam_ar);
                        let x0_2r = x0.pow_real(lam_2r);
                        let term_2a = x0_2a * (as_2a + b_2a);
                        let term_ar = x0_ar * (as_ar + b_ar);
                        let term_2r = x0_2r * (as_2r + b_2r);
                        let a22_tilde = term_2a - term_ar.mul_f(2.0) + term_2r;
                        let a22_tilde_e = x0_2a * (as_2a_e + b_2a_e)
                            - (x0_ar * (as_ar_e + b_ar_e)).mul_f(2.0)
                            + x0_2r * (as_2r_e + b_2r_e);
                        let c_fac2 = c_fac * c_fac;
                        let g21 = (khs * (a22_tilde + eta * a22_tilde_e) + eta * khs_e * a22_tilde)
                            .mul_f(c_fac2 / (8.0 * eps));
                        let g22 = (khs
                            * (term_ar.mul_f(lam_ar) - term_2a.mul_f(la) - term_2r.mul_f(lr)))
                        .mul_f(c_fac2 / (12.0 * eps));

                        let theta = beta.mul_f(eps).exp().add_f(-1.0);
                        let tanh_factor = 10.0 * (1.0 - (10.0 * (0.57 - alpha)).tanh());
                        let gam_c = (theta * zeta_bar * gam_exp).mul_f(tanh_factor);
                        gam_c.add_f(1.0) * (g21 + g22)
                    } else {
                        zero
                    };

                    let beta_eps = beta.mul_f(eps);
                    let w = (beta_eps * g1 + beta_eps.pow2() * g2) / gd;
                    let log_g_contact = exparg + w;

                    a_chain -= (x[i] * log_g_contact).mul_f(mi - 1.0);
                }
            }
        }

        Ok(Terms {
            hs: a_hs,
            a1: a1_mix,
            a2: a2_mix,
            a3: a3_mix,
            chain: a_chain,
            dispersion,
            hard_sphere,
            total: sum * (hard_sphere + dispersion + a_chain),
        })
    }
}

#[inline]
fn eff_eta_coeffs(lam: f64) -> [f64; 4] {
    let inv = 1.0 / lam;
    let inv2 = inv * inv;
    let inv3 = inv2 * inv;
    let inv_lam = [1.0, inv, inv2, inv3];
    let mut c_vec = [0.0; 4];
    for (j, cj) in c_vec.iter_mut().enumerate() {
        *cj = LAM_COEFF[0][j] * inv_lam[0]
            + LAM_COEFF[1][j] * inv_lam[1]
            + LAM_COEFF[2][j] * inv_lam[2]
            + LAM_COEFF[3][j] * inv_lam[3];
    }
    c_vec
}

#[inline]
fn calc_eff_eta<T: Scalar>(eta: T, lam: f64) -> T {
    let c = eff_eta_coeffs(lam);
    let p = eta.mul_f(c[3]).add_f(c[2]);
    let p = (eta * p).add_f(c[1]);
    let p = (eta * p).add_f(c[0]);
    eta * p
}

#[inline]
fn calc_eff_eta_with_deriv<T: Scalar>(eta: T, lam: f64) -> (T, T) {
    let c = eff_eta_coeffs(lam);
    let p = eta.mul_f(c[3]).add_f(c[2]);
    let p = (eta * p).add_f(c[1]);
    let ef = eta * (eta * p).add_f(c[0]);
    let dp = eta.mul_f(4.0 * c[3]).add_f(3.0 * c[2]);
    let dp = (eta * dp).add_f(2.0 * c[1]);
    let ef_e = (eta * dp).add_f(c[0]);
    (ef, ef_e)
}

#[inline]
fn calc_a1_sutherland<T: Scalar>(eta: T, lam: f64, eps_val: f64) -> T {
    let ef = calc_eff_eta(eta, lam);
    let pref = -12.0 * eps_val / (lam - 3.0);
    let om = T::sub_from_f(1.0, ef);
    let denum = om.pow3();
    T::sub_from_f(1.0, ef.mul_f(0.5)).mul_f(pref) / denum
}

#[inline]
fn calc_a1_sutherland_with_deriv<T: Scalar>(eta: T, lam: f64, eps_val: f64) -> (T, T) {
    let (ef, ef_e) = calc_eff_eta_with_deriv(eta, lam);
    let pref = 6.0 * eps_val / (lam - 3.0);
    let om = T::sub_from_f(1.0, ef);
    let denum = om.pow3();
    let as_val = T::sub_from_f(1.0, ef.mul_f(0.5)).mul_f(-2.0 * pref) / denum;
    let temp1 = T::sub_from_f(5.0, ef.mul_f(2.0)).mul_f(-pref) / (denum * om);
    let as_e = temp1 * ef_e;
    (as_val, as_e)
}

#[inline]
fn calc_b_tilde<T: Scalar>(x0: T, eta: T, lam: f64, eps_val: f64) -> T {
    let p3 = x0.pow_real(3.0 - lam);
    let p4 = p3 * x0;
    let i_term = p3.add_f(-1.0).mul_f(-1.0 / (lam - 3.0));
    let j_term = (p4.mul_f(lam - 3.0) - p3.mul_f(lam - 4.0))
        .add_f(-1.0)
        .mul_f(-1.0 / ((lam - 3.0) * (lam - 4.0)));
    let om = T::sub_from_f(1.0, eta);
    let inv_d3 = T::div_from_f(1.0, om.pow3());
    let keta10 = T::sub_from_f(2.0, eta) * inv_d3;
    let keta20 = (eta * eta.add_f(1.0)).mul_f(-9.0) * inv_d3;
    (keta10 * i_term + keta20 * j_term).mul_f(6.0 * eps_val)
}

#[inline]
fn calc_b_tilde_with_deriv<T: Scalar>(x0: T, eta: T, lam: f64, eps_val: f64) -> (T, T) {
    let p3 = x0.pow_real(3.0 - lam);
    let p4 = p3 * x0;
    let i_term = p3.add_f(-1.0).mul_f(-1.0 / (lam - 3.0));
    let j_term = (p4.mul_f(lam - 3.0) - p3.mul_f(lam - 4.0))
        .add_f(-1.0)
        .mul_f(-1.0 / ((lam - 3.0) * (lam - 4.0)));
    let om = T::sub_from_f(1.0, eta);
    let inv_om = T::div_from_f(1.0, om);
    let inv_d3 = inv_om.pow3();
    let inv_d4 = inv_d3 * inv_om;
    let keta10 = T::sub_from_f(2.0, eta) * inv_d3;
    let keta11 = T::sub_from_f(5.0, eta.mul_f(2.0)) * inv_d4;
    let keta20 = (eta * eta.add_f(1.0)).mul_f(-9.0) * inv_d3;
    let keta21 = (eta.pow2() + eta.mul_f(4.0)).add_f(1.0).mul_f(-9.0) * inv_d4;
    let b_val = (keta10 * i_term + keta20 * j_term).mul_f(6.0 * eps_val);
    let b_e = (keta11 * i_term + keta21 * j_term).mul_f(6.0 * eps_val);
    (b_val, b_e)
}

#[inline]
fn quantum_factors<T: Scalar>(lambda: f64, dfh: T, variant: QuantumVariant) -> [T; 3] {
    let first = dfh.mul_f(lambda * (lambda - 1.0));
    [
        T::real(1.0),
        first,
        (dfh * first).mul_f(variant.second_factor() * (lambda + 1.0) * (lambda + 2.0)),
    ]
}

#[inline]
fn alpha_functions<T: Scalar>(alpha: T) -> [T; 6] {
    std::array::from_fn(|k| {
        let num = alpha.mul_f(PHI_MAT[3][k]).add_f(PHI_MAT[2][k]);
        let num = (alpha * num).add_f(PHI_MAT[1][k]);
        let num = (alpha * num).add_f(PHI_MAT[0][k]);
        let den = alpha.mul_f(PHI_MAT[6][k]).add_f(PHI_MAT[5][k]);
        let den = (alpha * den).add_f(PHI_MAT[4][k]);
        let den = (alpha * den).add_f(1.0);
        num / den
    })
}

#[inline]
fn pure_hs<T: Scalar>(eta: T, kolafa: bool) -> T {
    let om = T::sub_from_f(1.0, eta);
    if kolafa {
        (T::div_from_f(2.5, om.pow2())
            + T::div_from_f(10.0, om)
            + eta.mul_f(2.0)
            + om.ln().mul_f(5.0))
        .mul_f(1.0 / 3.0)
    } else {
        (eta * T::sub_from_f(4.0, eta.mul_f(3.0))) / om.pow2()
    }
}

// Santos second and third virial coefficients in units of angstrom^3
// and angstrom^6. The cyclic B_ijk terms give identical sums over x_i x_j x_k.
fn santos<T: Scalar>(x: &[T], d: &[T], rho: T, segments: T, kolafa: bool) -> T {
    let nc = x.len();
    let zero = T::real(0.0);
    let diameter = |i: usize, j: usize| d[i * nc + j].mul_f(1e10);
    let mut moment = zero;
    let mut b2 = zero;
    let mut b3 = zero;
    for (i, &xi) in x.iter().enumerate() {
        moment += xi * diameter(i, i).pow3();
        for (j, &xj) in x.iter().enumerate() {
            let dij = diameter(i, j);
            let w = (xi * xj * dij.pow3()).mul_f(4.0);
            b2 += w;
            for (k, &xk) in x.iter().enumerate() {
                let dik = diameter(i, k);
                let djk = diameter(j, k);
                let clip = |v: T| if v.val() > 0.0 { v } else { zero };
                let dkij = clip(dik + djk - dij);
                let ckij = dkij.pow2()
                    * (dkij + (clip(dij + dik - djk) * clip(dij + djk - dik)).mul_f(1.5) / dij);
                b3 += w * xk * ckij;
            }
        }
    }
    let eta = (rho * moment).mul_f(PI * N_A * 1e-30 / 6.0);
    let om = T::sub_from_f(1.0, eta);
    let b2 = b2 / moment;
    let b3 = b3 / moment.pow2();
    (-((b2.mul_f(10.0) - b3.mul_f(4.0)) * om.ln()) + (b3 - b2) * pure_hs(eta, kolafa))
        / segments.mul_f(6.0)
}
