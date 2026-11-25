//! Dufal association integral and its TVN derivatives at unit total moles.
use super::{N_A, PreparedSaftVrMie, SaftVrMie};
use crate::Result;
use crate::association::{AssocRdfFirst, AssocRdfSecond, AssocRdfValue, AssociationRdfProvider};
use crate::rdf::RdfDerivatives;

/// Double-double (`hi + lo`) error-free arithmetic using pure IEEE-754
/// Veltkamp-Dekker splitting (avoiding external software `fma` calls on
/// baseline `x86_64` targets) for cancellation-sensitive evaluations of the
/// degree-10 Dufal association polynomial.
#[derive(Debug, Clone, Copy)]
pub(super) struct Dd(pub(super) f64, pub(super) f64);

impl Dd {
    const ZERO: Self = Self(0.0, 0.0);

    #[inline(always)]
    fn quick_two_sum(a: f64, b: f64) -> Self {
        let s = a + b;
        let e = b - (s - a);
        Self(s, e)
    }

    #[inline(always)]
    fn two_sum(a: f64, b: f64) -> Self {
        let s = a + b;
        let v = s - a;
        let e = (a - (s - v)) + (b - v);
        Self(s, e)
    }

    #[inline(always)]
    fn two_prod(a: f64, b: f64) -> Self {
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
        } else {
            Self(p, 0.0)
        }
    }

    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        let s = Self::two_sum(self.0, rhs.0);
        Self::quick_two_sum(s.0, s.1 + self.1 + rhs.1)
    }

    #[inline(always)]
    fn add_f64(self, rhs: f64) -> Self {
        let s = Self::two_sum(self.0, rhs);
        Self::quick_two_sum(s.0, s.1 + self.1)
    }

    #[inline(always)]
    fn mul_f64(self, rhs: f64) -> Self {
        let p = Self::two_prod(self.0, rhs);
        Self::quick_two_sum(p.0, p.1 + self.1 * rhs)
    }

    #[inline(always)]
    fn scale2(self) -> Self {
        Self(2.0 * self.0, 2.0 * self.1)
    }
}

/// Temperature-dependent outer coefficients `b_d(T)`, `db_d/dt`, `d^2b_d/dt^2`
/// (`d = 0..=10`, stored in Horner order `d = 10` down to `0`) precomputed once
/// in double-double precision per associating pair on [`PreparedSaftVrMie`].
#[derive(Debug, Clone, Copy)]
pub(super) struct AssocTempCoeffs {
    b: [Dd; 11],
    bt: [Dd; 11],
    btt: [Dd; 11],
}

impl AssocTempCoeffs {
    pub(super) fn compute(model: &SaftVrMie, pair_idx: usize, temperature: f64) -> Self {
        let eps = model.eps_divk_ij[pair_idx];
        let t = temperature / eps;
        let mut coef = model.integral_coefficients[pair_idx].iter().rev();
        let mut b_out = [Dd::ZERO; 11];
        let mut bt_out = [Dd::ZERO; 11];
        let mut btt_out = [Dd::ZERO; 11];
        for degree in 0..=10 {
            let (mut b, mut bt, mut btt) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
            for _ in 0..=degree {
                btt = btt.mul_f64(t).add(bt.scale2());
                bt = bt.mul_f64(t).add(b);
                b = b.mul_f64(t).add_f64(*coef.next().unwrap());
            }
            b_out[degree] = b;
            bt_out[degree] = bt;
            btt_out[degree] = btt;
        }
        Self {
            b: b_out,
            bt: bt_out,
            btt: btt_out,
        }
    }
}

pub(super) struct MieRdf<'a> {
    state: &'a PreparedSaftVrMie<'a>,
    density: f64,
    reduced: f64,
}

impl<'a> MieRdf<'a> {
    #[inline]
    pub(super) fn new(state: &'a PreparedSaftVrMie<'a>, density: f64) -> Self {
        Self {
            state,
            density,
            reduced: state.assoc_reduced_factor * density,
        }
    }

    /// Evaluate the degree-10 Dufal polynomial `g(r, t)` and requested derivatives.
    ///
    /// Uses precomputed temperature polynomials `b_d(t), bt_d(t), btt_d(t)` with an
    /// 11-step standard `f64` Horner fast path and a running cancellation monitor;
    /// when catastrophic cancellation is detected (`abs_sum > 1024 * |val|`, where
    /// more than 10 bits of mantissa could be lost on dense cold liquids), it re-evaluates
    /// only the 11-step density polynomial in Dekker double-double precision so
    /// dense cold liquids remain smooth without slowing routine states.
    pub(super) fn integral<const ORDER: u8>(&self, i: usize, j: usize) -> [f64; 6] {
        let model = self.state.model;
        let idx = i * model.num_components() + j;
        let eps = model.eps_divk_ij[idx];
        let r = self.reduced;
        let owned;
        let coeffs = match self.state.assoc_t_coeffs.get(idx).and_then(|c| c.as_ref()) {
            Some(c) => c,
            None => {
                owned = AssocTempCoeffs::compute(model, idx, self.state.temperature);
                &owned
            }
        };

        // Fast 11-step f64 Horner path with cheap running absolute-sum cancellation monitor.
        let ar = r.abs();
        let (mut g, mut gr, mut gt, mut grr, mut grt, mut gtt) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let (mut abs_g, mut abs_gr, mut abs_gt) = (0.0, 0.0, 0.0);
        for d in 0..=10 {
            let bd = coeffs.b[d].0;
            if ORDER >= 2 {
                grr = grr * r + 2.0 * gr;
                grt = grt * r + gt;
                gtt = gtt * r + coeffs.btt[d].0;
            }
            if ORDER >= 1 {
                let btd = coeffs.bt[d].0;
                gr = gr * r + g;
                gt = gt * r + btd;
                abs_gr = abs_gr * ar + abs_g;
                abs_gt = abs_gt * ar + btd.abs();
            }
            g = g * r + bd;
            abs_g = abs_g * ar + bd.abs();
        }

        const CANCEL_GATE: f64 = 1024.0;
        let cancelled = abs_g > CANCEL_GATE * g.abs()
            || (ORDER >= 1 && (abs_gr > CANCEL_GATE * gr.abs() || abs_gt > CANCEL_GATE * gt.abs()));
        if !cancelled {
            return [g, gr, gt / eps, grr, grt / eps, gtt / (eps * eps)];
        }

        self.integral_dd::<ORDER>(coeffs, eps, r)
    }

    #[inline(never)]
    fn integral_dd<const ORDER: u8>(&self, coeffs: &AssocTempCoeffs, eps: f64, r: f64) -> [f64; 6] {
        let (mut g, mut gr, mut gt, mut grr, mut grt, mut gtt) =
            (Dd::ZERO, Dd::ZERO, Dd::ZERO, Dd::ZERO, Dd::ZERO, Dd::ZERO);
        for d in 0..=10 {
            if ORDER >= 2 {
                grr = grr.mul_f64(r).add(gr.scale2());
                grt = grt.mul_f64(r).add(gt);
                gtt = gtt.mul_f64(r).add(coeffs.btt[d]);
            }
            if ORDER >= 1 {
                gr = gr.mul_f64(r).add(g);
                gt = gt.mul_f64(r).add(coeffs.bt[d]);
            }
            g = g.mul_f64(r).add(coeffs.b[d]);
        }
        [
            g.0,
            gr.0,
            gt.0 / eps,
            grr.0,
            grt.0 / eps,
            gtt.0 / (eps * eps),
        ]
    }

    #[cfg(test)]
    pub(super) fn integral_high_precision<const ORDER: u8>(&self, i: usize, j: usize) -> [f64; 6] {
        let model = self.state.model;
        let idx = i * model.num_components() + j;
        let eps = model.eps_divk_ij[idx];
        let coeffs = AssocTempCoeffs::compute(model, idx, self.state.temperature);
        self.integral_dd::<ORDER>(&coeffs, eps, self.reduced)
    }

    fn reduced_moles(&self) -> Vec<f64> {
        let p = self.state;
        let nc = p.model.num_components();
        (0..nc)
            .map(|k| {
                let cross: f64 =
                    p.xs.iter()
                        .enumerate()
                        .map(|(l, x)| x * p.model.sigma_ij[k * nc + l].powi(3))
                        .sum();
                p.model.components[k].m * (2. * N_A * self.density * cross - self.reduced / p.m_bar)
            })
            .collect()
    }
}

impl AssociationRdfProvider for MieRdf<'_> {
    fn pair_dependent(&self) -> bool {
        true
    }

    fn rdf_value(&self, i: usize, j: usize) -> Result<AssocRdfValue> {
        Ok(AssocRdfValue {
            g: self.integral::<0>(i, j)[0],
            covol: N_A,
        })
    }

    fn rdf_first(&self, i: usize, j: usize) -> Result<AssocRdfFirst> {
        let [g, gr, gt, ..] = self.integral::<1>(i, j);
        Ok(AssocRdfFirst {
            g,
            g_t: gt,
            g_v: -gr * self.reduced * self.density,
            g_n: self.reduced_moles().into_iter().map(|v| gr * v).collect(),
            covol: N_A,
            covol_t: 0.,
        })
    }

    fn rdf_second(&self, i: usize, j: usize) -> Result<AssocRdfSecond> {
        let p = self.state;
        let nc = p.model.num_components();
        let [g, gr, gt, grr, grt, gtt] = self.integral::<2>(i, j);
        let rn = self.reduced_moles();
        let rv = -self.reduced * self.density;
        let mut nn = Vec::with_capacity(nc * nc);
        for k in 0..nc {
            let mk = p.model.components[k].m;
            for l in 0..nc {
                let ml = p.model.components[l].m;
                let rkl =
                    (2. * N_A * self.density * mk * ml * p.model.sigma_ij[k * nc + l].powi(3)
                        - mk * rn[l]
                        - ml * rn[k])
                        / p.m_bar;
                nn.push(grr * rn[k] * rn[l] + gr * rkl);
            }
        }
        Ok(AssocRdfSecond {
            covol: N_A,
            covol_t: 0.,
            covol_tt: 0.,
            rdf: RdfDerivatives {
                g,
                temperature: gt,
                volume: gr * rv,
                moles: rn.iter().map(|d| gr * d).collect(),
                temperature_second: gtt,
                temperature_volume: grt * rv,
                temperature_moles: rn.iter().map(|d| grt * d).collect(),
                volume_second: grr * rv * rv - 2. * gr * rv * self.density,
                volume_moles: rn
                    .iter()
                    .map(|d| (grr * rv - gr * self.density) * d)
                    .collect(),
                moles_second: nn,
            },
        })
    }
}
