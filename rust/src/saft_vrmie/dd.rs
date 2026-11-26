//! Pure IEEE-754 double-double (`hi + lo`) arithmetic using Veltkamp-Dekker
//! error-free transformations (zero `unsafe`, zero software `fma()` DLL calls
//! on baseline `x86_64` targets, and CRT-independent transcendental evaluation).

/// Double-double number representing `self.0 + self.1` with `|self.1| <= 0.5 * ulp(self.0)`.
#[derive(Debug, Clone, Copy)]
pub(super) struct Dd(pub(super) f64, pub(super) f64);

const LN2_DD: Dd = Dd(std::f64::consts::LN_2, 2.319_046_813_846_299_6e-17);
const INV_LN2: f64 = std::f64::consts::LOG2_E;
#[allow(clippy::excessive_precision)]
pub(super) const PI_NA_DIV6_DD: Dd = Dd(3.153_185_528_416_608_99e23, 4.283_795_407_596_530_39e6);

/// Reciprocal factorials `1 / k!` for `k = 2..=12` in 106-bit double-double precision.
#[allow(clippy::excessive_precision)]
const INV_FACT: [Dd; 11] = [
    Dd(5.000_000_000_000_000_00e-1, 0.0),
    Dd(1.666_666_666_666_666_57e-1, 9.251_858_538_542_970_66e-18),
    Dd(4.166_666_666_666_666_44e-2, 2.312_964_634_635_742_66e-18),
    Dd(8.333_333_333_333_333_22e-3, 1.156_482_317_317_871_38e-19),
    Dd(1.388_888_888_888_888_94e-3, -5.300_543_954_373_577_06e-20),
    Dd(1.984_126_984_126_984_13e-4, 1.720_955_829_342_070_53e-22),
    Dd(2.480_158_730_158_730_16e-5, 2.151_194_786_677_588_16e-23),
    Dd(2.755_731_922_398_589_25e-6, -1.858_393_274_046_472_08e-22),
    Dd(2.755_731_922_398_588_83e-7, 2.376_771_462_225_029_73e-23),
    Dd(2.505_210_838_544_172_02e-8, -1.448_814_070_935_911_97e-24),
    Dd(2.087_675_698_786_810_02e-9, -1.207_345_059_113_259_97e-25),
];

impl Dd {
    pub(super) const ZERO: Self = Self(0.0, 0.0);
    pub(super) const ONE: Self = Self(1.0, 0.0);

    #[inline(always)]
    pub(super) const fn from_f64(x: f64) -> Self {
        Self(x, 0.0)
    }

    #[inline(always)]
    pub(super) fn is_finite(self) -> bool {
        self.0.is_finite() && self.1.is_finite()
    }

    #[inline(always)]
    pub(super) fn quick_two_sum(a: f64, b: f64) -> Self {
        let s = a + b;
        let e = b - (s - a);
        Self(s, e)
    }

    #[inline(always)]
    pub(super) fn two_sum(a: f64, b: f64) -> Self {
        let s = a + b;
        let v = s - a;
        let e = (a - (s - v)) + (b - v);
        Self(s, e)
    }

    #[inline(always)]
    pub(super) fn two_prod(a: f64, b: f64) -> Self {
        const SPLIT: f64 = 134_217_729.0; // 2^27 + 1
        let p = a * b;
        let ca = SPLIT * a;
        let cb = SPLIT * b;
        if p.is_finite() && ca.is_finite() && cb.is_finite() {
            let a_hi = ca - (ca - a);
            let a_lo = a - a_hi;
            let b_hi = cb - (cb - b);
            let b_lo = b - b_hi;
            let err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
            Self(p, err)
        } else if p.is_finite() {
            Self::two_prod_prescaled(a, b, p, ca, cb)
        } else {
            Self(p, 0.0)
        }
    }

    #[cold]
    #[inline(never)]
    fn two_prod_prescaled(a: f64, b: f64, p: f64, ca: f64, cb: f64) -> Self {
        const SPLIT: f64 = 134_217_729.0; // 2^27 + 1
        const SCALE_UP: f64 = 268_435_456.0; // 2^28
        const SCALE_DOWN: f64 = 1.0 / SCALE_UP; // 2^-28 (exact)
        let (sa, fa) = if ca.is_finite() {
            (a, 1.0)
        } else {
            (a * SCALE_DOWN, SCALE_UP)
        };
        let (sb, fb) = if cb.is_finite() {
            (b, 1.0)
        } else {
            (b * SCALE_DOWN, SCALE_UP)
        };
        let ca2 = SPLIT * sa;
        let cb2 = SPLIT * sb;
        if ca2.is_finite() && cb2.is_finite() {
            let ps = sa * sb;
            let a_hi = ca2 - (ca2 - sa);
            let a_lo = sa - a_hi;
            let b_hi = cb2 - (cb2 - sb);
            let b_lo = sb - b_hi;
            let err_s = ((a_hi * b_hi - ps) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
            let err = (err_s * fa) * fb;
            if err.is_finite() {
                return Self(p, err);
            }
        }
        Self(p, 0.0)
    }

    #[inline(always)]
    pub(super) fn neg(self) -> Self {
        Self(-self.0, -self.1)
    }

    #[inline(always)]
    pub(super) fn add(self, rhs: Self) -> Self {
        let s = Self::two_sum(self.0, rhs.0);
        let t = Self::two_sum(self.1, rhs.1);
        let c = Self::quick_two_sum(s.0, s.1 + t.0);
        Self::quick_two_sum(c.0, c.1 + t.1)
    }

    #[inline(always)]
    pub(super) fn sub(self, rhs: Self) -> Self {
        self.add(rhs.neg())
    }

    #[inline(always)]
    pub(super) fn add_f64(self, rhs: f64) -> Self {
        let s = Self::two_sum(self.0, rhs);
        Self::quick_two_sum(s.0, s.1 + self.1)
    }

    #[inline(always)]
    pub(super) fn mul(self, rhs: Self) -> Self {
        let p = Self::two_prod(self.0, rhs.0);
        Self::quick_two_sum(p.0, p.1 + (self.0 * rhs.1 + self.1 * rhs.0))
    }

    #[inline(always)]
    pub(super) fn mul_f64(self, rhs: f64) -> Self {
        let p = Self::two_prod(self.0, rhs);
        Self::quick_two_sum(p.0, p.1 + self.1 * rhs)
    }

    #[inline(always)]
    pub(super) fn scale2(self) -> Self {
        Self(2.0 * self.0, 2.0 * self.1)
    }

    #[inline]
    pub(super) fn div(self, rhs: Self) -> Self {
        if rhs.0 == 0.0 {
            return Self::from_f64(self.0 / rhs.0);
        }
        let q1 = self.0 / rhs.0;
        let r1 = self.sub(rhs.mul_f64(q1));
        let q2 = r1.0 / rhs.0;
        let r2 = r1.sub(rhs.mul_f64(q2));
        let q3 = r2.0 / rhs.0;
        let s = Self::quick_two_sum(q1, q2);
        s.add_f64(q3)
    }

    #[inline]
    pub(super) fn sqrt(self) -> Self {
        if self.0 <= 0.0 || !self.0.is_finite() {
            return Self::from_f64(self.0.sqrt());
        }
        let x = Self::from_f64(self.0.sqrt());
        x.add(self.div(x)).mul_f64(0.5)
    }

    #[inline]
    pub(super) fn scale_pow2(self, k: i32) -> Self {
        if k == 0 {
            return self;
        }
        if (-1022..=1023).contains(&k) {
            let factor = f64::from_bits(((k + 1023) as u64) << 52);
            let hi = self.0 * factor;
            if hi.is_finite() {
                Self(hi, self.1 * factor)
            } else {
                Self(hi, 0.0)
            }
        } else {
            let half = k / 2;
            self.scale_pow2(half).scale_pow2(k - half)
        }
    }

    /// Return `(mantissa, k)` such that `exp(self) == mantissa * 2^k` with
    /// `mantissa` in `[sqrt(0.5), sqrt(2)]`, for finite `self` with `|self| < 1e6`.
    #[inline]
    pub(super) fn exp_mantissa_and_pow2(self) -> (Self, i32) {
        let k_f = (self.0 * INV_LN2).round();
        let k = k_f as i32;
        let r = self.sub(LN2_DD.mul_f64(k_f));
        // Reduce further by 2^-4 so |s| <= ln(2)/32 < 0.0217 -> s^13/13! < 4e-32.
        let s = Self(r.0 * 0.0625, r.1 * 0.0625);
        let mut poly = INV_FACT[10];
        for idx in (0..10).rev() {
            poly = poly.mul(s).add(INV_FACT[idx]);
        }
        let mut e = Self::ONE.add(s.add(s.mul(s).mul(poly)));
        for _ in 0..4 {
            e = e.mul(e);
        }
        (e, k)
    }

    /// CRT-independent 106-bit `exp(self)` via `2^k * exp(r / 16)^16` with degree-12 Taylor.
    #[inline]
    pub(super) fn exp(self) -> Self {
        if !self.0.is_finite() {
            return Self::from_f64(self.0.exp());
        }
        // Leave the boundary near ln(f64::MAX) to the final binary scaling.
        if self.0 >= 710.0 {
            return Self::from_f64(f64::INFINITY);
        }
        if self.0 < -745.0 {
            return Self::ZERO;
        }
        let (e, k) = self.exp_mantissa_and_pow2();
        e.scale_pow2(k)
    }

    /// CRT-independent 106-bit `exp(self) - 1`.
    #[inline]
    pub(super) fn exp_m1(self) -> Self {
        if self.0.abs() < 0.025 {
            let s = self;
            let mut poly = INV_FACT[10];
            for idx in (0..10).rev() {
                poly = poly.mul(s).add(INV_FACT[idx]);
            }
            s.add(s.mul(s).mul(poly))
        } else {
            let exp = self.exp();
            if exp.0.is_finite() {
                exp.sub(Self::ONE)
            } else {
                exp
            }
        }
    }

    /// CRT-independent 106-bit `ln(self)` via 2 Newton steps on `exp(y) = self`.
    #[inline]
    pub(super) fn ln(self) -> Self {
        if self.0 <= 0.0 || !self.0.is_finite() {
            return Self::from_f64(self.0.ln());
        }
        let mut y = Self::from_f64(self.0.ln());
        let e1 = y.exp();
        y = y.add(self.sub(e1).div(e1));
        let e2 = y.exp();
        y.add(self.sub(e2).div(e2))
    }

    /// Exact binary exponentiation for integer exponents in 106-bit precision.
    #[inline]
    pub(super) fn powi(self, n: i32) -> Self {
        if n == 0 {
            return Self::ONE;
        }
        let mut base = if n < 0 { Self::ONE.div(self) } else { self };
        let mut exp = n.unsigned_abs();
        let mut acc = Self::ONE;
        while exp > 0 {
            if exp & 1 == 1 {
                acc = acc.mul(base);
            }
            exp >>= 1;
            if exp > 0 {
                base = base.mul(base);
            }
        }
        acc
    }

    /// CRT-independent 106-bit `self^exp` (binary exponentiation for small exact integers,
    /// `exp(exp * ln(self))` otherwise).
    #[inline]
    pub(super) fn pow_real(self, exp: f64) -> Self {
        if (-64.0..=64.0).contains(&exp) && exp == (exp as i32) as f64 {
            return self.powi(exp as i32);
        }
        self.ln().mul_f64(exp).exp()
    }

    /// CRT-independent 106-bit `tanh(self)`.
    #[inline]
    pub(super) fn tanh(self) -> Self {
        if self.0 > 40.0 {
            return Self::ONE;
        }
        if self.0 < -40.0 {
            return Self::ONE.neg();
        }
        let e2 = self.scale2().exp();
        e2.sub(Self::ONE).div(e2.add(Self::ONE))
    }
}
