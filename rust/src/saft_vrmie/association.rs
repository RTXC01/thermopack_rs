//! Dufal association integral and its TVN derivatives at unit total moles.
use super::dd::Dd;
use super::{N_A, PreparedSaftVrMie, RGAS, SaftVrMie};
use crate::association::{
    AssocRdfFirst, AssocRdfSecond, AssocRdfValue, AssociationRdfProvider, DeltaCombRule,
};
use crate::rdf::RdfDerivatives;
use crate::{Error, Result};

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
        let t = Dd::from_f64(temperature).div(Dd::from_f64(eps));
        let mut coef = model.integral_coefficients[pair_idx].iter().rev();
        let mut b_out = [Dd::ZERO; 11];
        let mut bt_out = [Dd::ZERO; 11];
        let mut btt_out = [Dd::ZERO; 11];
        for degree in 0..=10 {
            let (mut b, mut bt, mut btt) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
            for _ in 0..=degree {
                btt = btt.mul(t).add(bt.scale2());
                bt = bt.mul(t).add(b);
                b = b.mul(t).add_f64(*coef.next().unwrap());
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

        let r_dd = Dd::two_prod(self.state.assoc_reduced_factor, self.density);
        self.integral_dd::<ORDER>(coeffs, eps, r_dd)
    }

    #[inline(never)]
    fn integral_dd<const ORDER: u8>(&self, coeffs: &AssocTempCoeffs, eps: f64, r: Dd) -> [f64; 6] {
        let (mut g, mut gr, mut gt, mut grr, mut grt, mut gtt) =
            (Dd::ZERO, Dd::ZERO, Dd::ZERO, Dd::ZERO, Dd::ZERO, Dd::ZERO);
        for d in 0..=10 {
            if ORDER >= 2 {
                grr = grr.mul(r).add(gr.scale2());
                grt = grt.mul(r).add(gt);
                gtt = gtt.mul(r).add(coeffs.btt[d]);
            }
            if ORDER >= 1 {
                gr = gr.mul(r).add(g);
                gt = gt.mul(r).add(coeffs.bt[d]);
            }
            g = g.mul(r).add(coeffs.b[d]);
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

    #[inline]
    fn integral_g_gr_dd(coeffs: &AssocTempCoeffs, r: Dd) -> (Dd, Dd) {
        let (mut g, mut gr) = (Dd::ZERO, Dd::ZERO);
        for &bd in &coeffs.b {
            gr = gr.mul(r).add(g);
            g = g.mul(r).add(bd);
        }
        (g, gr)
    }

    #[cfg(test)]
    pub(super) fn integral_high_precision<const ORDER: u8>(&self, i: usize, j: usize) -> [f64; 6] {
        let model = self.state.model;
        let idx = i * model.num_components() + j;
        let eps = model.eps_divk_ij[idx];
        let coeffs = AssocTempCoeffs::compute(model, idx, self.state.temperature);
        let r_dd = Dd::two_prod(self.state.assoc_reduced_factor, self.density);
        self.integral_dd::<ORDER>(&coeffs, eps, r_dd)
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

/// Evaluate the unit-mole association volume derivative `dF_assoc / dV` in
/// 106-bit double-double precision using the Wertheim/Michelsen envelope
/// theorem (`dF_assoc/dV = (dQ/dV)_X`), eliminating `r = factor * rho`
/// quantization and `f64` site-fraction/summation rounding on cancellation-
/// sensitive dense liquid states.
#[inline(never)]
pub(super) fn assoc_f_v_dd(
    state: &PreparedSaftVrMie<'_>,
    rho: f64,
    x_k_hint: &[f64],
) -> Result<Dd> {
    let assoc = &state.model.association;
    let ns = assoc.num_sites;
    if ns == 0 {
        return Ok(Dd::ZERO);
    }
    let nc = state.model.num_components();
    let rho_dd = Dd::from_f64(rho);
    let mut sum_sig3 = Dd::ZERO;
    for i in 0..nc {
        let xsi = Dd::from_f64(state.xs[i]);
        for j in 0..nc {
            let s = Dd::from_f64(state.model.sigma_ij[i * nc + j]);
            let s3 = s.mul(s).mul(s);
            sum_sig3 = sum_sig3.add(xsi.mul_f64(state.xs[j]).mul(s3));
        }
    }
    let r_dd = sum_sig3.mul_f64(N_A).mul_f64(state.m_bar).mul(rho_dd);
    let rt_dd = Dd::two_prod(RGAS, state.temperature);

    // Precompute pair (g, gr) in double-double precision for active associating pairs.
    let mut g_pair = vec![(Dd::ZERO, Dd::ZERO); nc * nc];
    for i in 0..nc {
        if assoc.comp_vs_sites[i].is_none() {
            continue;
        }
        for j in i..nc {
            if assoc.comp_vs_sites[j].is_none() {
                continue;
            }
            let idx = i * nc + j;
            let owned;
            let coeffs = match state.assoc_t_coeffs.get(idx).and_then(|c| c.as_ref()) {
                Some(c) => c,
                None => {
                    owned = AssocTempCoeffs::compute(state.model, idx, state.temperature);
                    &owned
                }
            };
            let pair_val = MieRdf::integral_g_gr_dd(coeffs, r_dd);
            g_pair[idx] = pair_val;
            g_pair[j * nc + i] = pair_val;
        }
    }

    let mut delta = vec![Dd::ZERO; ns * ns];
    // delta_rho_minus_v[kl] = Delta_kl * rho - d(Delta_kl)/dV
    // For Standard: (g + r * gr) * rho * beta * N_A * (exp(eps/(RT)) - 1)
    let mut delta_v = vec![Dd::ZERO; ns * ns];

    for k in 0..ns {
        let ic = assoc.site_to_comp[k];
        for l in k..ns {
            let jc = assoc.site_to_comp[l];
            if assoc.delta_comb_rule == DeltaCombRule::Elliot && ic != jc {
                continue;
            }
            let kl = k * ns + l;
            let lk = l * ns + k;
            let beta = assoc.beta_kl[kl];
            if beta == 0.0 {
                continue;
            }
            let eps = assoc.eps_kl[kl];
            let h = if eps > 0.0 {
                let arg = Dd::from_f64(eps).div(rt_dd);
                let bn = Dd::two_prod(beta, N_A);
                if arg.0 > 100.0 && arg.0 < 1.0e6 {
                    // For arg > 100, exp(-arg) < 3.8e-44 << 2^-106 so exp(arg) - 1 == exp(arg)
                    // to full 106-bit precision; multiplying `beta * N_A` by the O(1) mantissa
                    // before binary `2^k` scaling avoids intermediate overflow near `ln(f64::MAX)`.
                    let (mant, k) = arg.exp_mantissa_and_pow2();
                    bn.mul(mant).scale_pow2(k)
                } else {
                    bn.mul(arg.exp_m1())
                }
            } else {
                Dd::ZERO
            };
            let (g, gr) = g_pair[ic * nc + jc];
            let d = g.mul(h);
            let dv = gr.mul(r_dd).mul(rho_dd).mul(h).neg();
            if !d.0.is_finite() || !dv.0.is_finite() {
                return Err(Error::InvalidInput(
                    "non-finite association strength in high-precision evaluation",
                ));
            }
            delta[kl] = d;
            delta[lk] = d;
            delta_v[kl] = dv;
            delta_v[lk] = dv;
        }
    }

    if assoc.delta_comb_rule == DeltaCombRule::Elliot {
        for ic in 0..nc {
            let Some((k1, k2)) = assoc.comp_vs_sites[ic] else {
                continue;
            };
            for jc in (ic + 1)..nc {
                let Some((l1, l2)) = assoc.comp_vs_sites[jc] else {
                    continue;
                };
                let ka = k1 * ns + k2;
                let lb = l1 * ns + l2;
                let da = delta[ka];
                let db = delta[lb];
                if da.0 < 0.0 || db.0 < 0.0 {
                    return Err(Error::InvalidInput("negative Delta in Elliot rule"));
                }
                let d_ab = da.sqrt().mul(db.sqrt());
                let dv_ab = if d_ab.0 != 0.0 && da.0 > 0.0 && db.0 > 0.0 {
                    d_ab.mul(delta_v[ka].div(da).add(delta_v[lb].div(db)))
                        .mul_f64(0.5)
                } else {
                    Dd::ZERO
                };
                if !d_ab.0.is_finite() || !dv_ab.0.is_finite() {
                    return Err(Error::InvalidInput(
                        "non-finite Elliot association strength in high-precision evaluation",
                    ));
                }
                for k in k1..=k2 {
                    for l in l1..=l2 {
                        let kl = k * ns + l;
                        let lk = l * ns + k;
                        if assoc.beta_kl[kl] <= 0.0 {
                            continue;
                        }
                        delta[kl] = d_ab;
                        delta[lk] = d_ab;
                        delta_v[kl] = dv_ab;
                        delta_v[lk] = dv_ab;
                    }
                }
            }
        }
    }

    // Refine X_k from the 53-bit f64 hint to 106-bit double-double via damped
    // log-space Newton (`X_k = exp(u_k)`) on the O(1)-scaled site residual
    //   r_k(u) = 1 - X_k - X_k * sum_l (rho * m_l * Delta_kl) * X_l = 0.
    // Log-space steps stay strictly positive, avoid `1/X_k^2` scaling blow-up,
    // and move along constant-product manifolds (`X_A * X_B = const`) when an
    // extreme `f64` hint exhibits bipartite symmetry breaking.
    let m_u: Vec<Dd> = assoc
        .site_to_comp
        .iter()
        .map(|&ic| Dd::from_f64(state.mole_fractions[ic]))
        .collect();
    if x_k_hint.len() != ns {
        return Err(Error::InvalidInput(
            "association site fraction hint mismatch",
        ));
    }
    let mut x_k = Vec::with_capacity(ns);
    for &v in x_k_hint {
        if !v.is_finite() || v <= 0.0 {
            return Err(Error::InvalidInput(
                "invalid association site fraction hint",
            ));
        }
        x_k.push(Dd::from_f64(v));
    }

    let mut jac = vec![Dd::ZERO; ns * ns];
    let mut rhs = vec![Dd::ZERO; ns];
    let mut last_max_res = f64::INFINITY;
    for _ in 0..24 {
        let mut max_res = 0.0_f64;
        for k in 0..ns {
            let xk = x_k[k];
            let mut sum_cx = Dd::ZERO;
            for l in 0..ns {
                let cx_l = delta[k * ns + l].mul(x_k[l]).mul(rho_dd).mul(m_u[l]);
                sum_cx = sum_cx.add(cx_l);
                jac[k * ns + l] = xk.mul(cx_l).neg();
            }
            let x_sum = xk.mul(sum_cx);
            let res = Dd::ONE.sub(xk).sub(x_sum);
            max_res = max_res.max(res.0.abs());
            rhs[k] = res.neg();
            jac[k * ns + k] = jac[k * ns + k].sub(xk.add(x_sum));
        }
        last_max_res = max_res;
        if !max_res.is_finite() {
            return Err(Error::NotConverged(
                "non-finite residual in high-precision association site solver",
            ));
        }
        if max_res <= 1.0e-28 {
            break;
        }
        let Some(du) = solve_linear_dd(ns, &mut jac, &mut rhs) else {
            return Err(Error::NotConverged(
                "singular Jacobian in high-precision association site solver",
            ));
        };
        let mut max_du = 0.0_f64;
        for d in &du {
            if !d.0.is_finite() {
                return Err(Error::NotConverged(
                    "non-finite step in high-precision association site solver",
                ));
            }
            max_du = max_du.max(d.0.abs());
        }
        let scale = if max_du > 2.0 { 2.0 / max_du } else { 1.0 };
        for k in 0..ns {
            let step = if scale < 1.0 {
                du[k].mul_f64(scale)
            } else {
                du[k]
            };
            let next = x_k[k].mul(step.exp());
            if !next.0.is_finite() || next.0 <= 0.0 {
                return Err(Error::NotConverged(
                    "invalid site fraction update in high-precision association solver",
                ));
            }
            x_k[k] = if next.0 > 1.0 { Dd::ONE } else { next };
        }
        if max_du <= 1.0e-28 {
            break;
        }
    }
    if !last_max_res.is_finite() || last_max_res > 1.0e-12 {
        return Err(Error::NotConverged(
            "high-precision association site solver did not converge",
        ));
    }

    let mut f_v = Dd::ZERO;
    for k in 0..ns {
        let mk_xk = m_u[k].mul(x_k[k]);
        if mk_xk.0 == 0.0 {
            continue;
        }
        let half_mk_xk_rho = mk_xk.mul(rho_dd).mul_f64(0.5);
        for l in 0..ns {
            let kl = k * ns + l;
            if delta[kl].0 == 0.0 && delta_v[kl].0 == 0.0 {
                continue;
            }
            let ml_xl = m_u[l].mul(x_k[l]);
            let term_x = delta[kl].mul(ml_xl).mul(rho_dd).sub(delta_v[kl].mul(ml_xl));
            f_v = f_v.add(half_mk_xk_rho.mul(term_x));
        }
    }
    if !f_v.0.is_finite() {
        return Err(Error::NotConverged(
            "non-finite high-precision association volume derivative",
        ));
    }
    Ok(f_v)
}

fn solve_linear_dd(n: usize, a: &mut [Dd], b: &mut [Dd]) -> Option<Vec<Dd>> {
    for col in 0..n {
        let mut piv = col;
        let mut max_val = a[col * n + col].0.abs();
        for row in (col + 1)..n {
            let v = a[row * n + col].0.abs();
            if v > max_val {
                max_val = v;
                piv = row;
            }
        }
        if max_val == 0.0 || !max_val.is_finite() {
            return None;
        }
        if piv != col {
            for j in col..n {
                a.swap(col * n + j, piv * n + j);
            }
            b.swap(col, piv);
        }
        let diag = a[col * n + col];
        for row in (col + 1)..n {
            let fac = a[row * n + col].div(diag);
            if fac.0 == 0.0 && fac.1 == 0.0 {
                continue;
            }
            for j in (col + 1)..n {
                a[row * n + j] = a[row * n + j].sub(fac.mul(a[col * n + j]));
            }
            b[row] = b[row].sub(fac.mul(b[col]));
        }
    }
    let mut x = vec![Dd::ZERO; n];
    for row in (0..n).rev() {
        let mut s = b[row];
        for j in (row + 1)..n {
            s = s.sub(a[row * n + j].mul(x[j]));
        }
        x[row] = s.div(a[row * n + row]);
    }
    Some(x)
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
