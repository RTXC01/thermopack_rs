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
            if !d.is_finite() || !dv.is_finite() {
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
                if !d_ab.is_finite() || !dv_ab.is_finite() {
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

    // Refine X_k to 106-bit double-double precision.
    //
    // First compute a conservative exact equitable partition (1-WL color refinement
    // using canonical sorted nonzero coefficient multisets) of the scaled site
    // interaction matrix `C_{kl} = rho * m_l * Delta_{kl}`: whenever every site in
    // class `a` has the exact same multiset of nonzero weights into every class `b`,
    // their exact real row sums `W_{ab}` are guaranteed identical without rounding
    // ambiguity, so any positive solution `Y` of `1/Y_a = 1 + sum_b W_{ab} Y_b`
    // expands to the unique Wertheim equilibrium `X_k = Y_{class(k)}`.
    //
    // This merges equivalent donor/acceptor classes ONLY when the interaction graph
    // and composition mathematically justify it (e.g. pure 4C/2B/4A -> scalar
    // quadratic, symmetric 4C+2B -> strictly diagonally dominant 2x2, asymmetric
    // 4C+1ea or 4C+3B -> donor and acceptor classes stay distinct), eliminating the
    // finite-precision singular constant-product nullspace when `X_k < 2^-106`.
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
    for &v in x_k_hint {
        if !v.is_finite() || v <= 0.0 {
            return Err(Error::InvalidInput(
                "invalid association site fraction hint",
            ));
        }
    }

    let mut c_mat = vec![Dd::ZERO; ns * ns];
    for k in 0..ns {
        for l in 0..ns {
            let d_kl = delta[k * ns + l];
            if d_kl.0 == 0.0 && d_kl.1 == 0.0 {
                continue;
            }
            let c_kl = rho_dd.mul(m_u[l]).mul(d_kl);
            if !c_kl.is_finite() || c_kl.0 < 0.0 || (c_kl.0 == 0.0 && c_kl.1 < 0.0) {
                return Err(Error::InvalidInput(
                    "invalid scaled association matrix in high-precision evaluation",
                ));
            }
            c_mat[k * ns + l] = c_kl;
        }
    }

    let part = equitable_site_partition(ns, &c_mat)?;
    let ng = part.num_classes;
    let w_mat = &part.w_mat;
    let mut y_g = vec![Dd::ONE; ng];
    if ng == 1 {
        let w = w_mat[0];
        let y_bound = solve_positive_quadratic_dd(Dd::ONE, w)?;
        let res = Dd::ONE.sub(y_bound).sub(w.mul(y_bound).mul(y_bound));
        if !res.is_finite() || res.0.abs() > 1.0e-12 {
            return Err(Error::NotConverged(
                "symmetric site fraction residual failed tolerance in high-precision association solver",
            ));
        }
        y_g[0] = y_bound;
    } else if ng == 2 && is_dd_zero(w_mat[0]) && is_dd_zero(w_mat[3]) {
        // Exact 2-class bipartite reduction (`1/Y_0 = 1 + W_01 Y_1`, `1/Y_1 = 1 + W_10 Y_0`).
        // Cancel identical terms across the two class-aggregated multisets before forming
        // `D = W_01 - W_10` so trace couplings (e.g. `1e-40`) are retained without O(1)
        // cancellation even when `Y_0, Y_1 << 2^-106`.
        let (y0, y1) = solve_bipartite_two_class_dd(&part)?;
        y_g[0] = y0;
        y_g[1] = y1;
    } else if let Some(hub) = find_star_bipartite_hub(ng, w_mat) {
        // Exact star-bipartite reduction (one hub class coupled to >=2 non-interacting leaf
        // classes, e.g. 4C + 1ea with distinct cross-association strength): eliminate all leaf
        // classes algebraically and solve the strictly monotone 1D hub equation in DD log-space.
        y_g = solve_star_bipartite_dd(&part, hub)?;
    } else {
        // Initialize each class from the geometric mean of its member hints
        // (preserving exact single-site hints and immediately restoring the
        // symmetric constant-product midpoint when a bipartite f64 hint drifted).
        for (a, slot) in y_g.iter_mut().enumerate() {
            let mut count = 0usize;
            let mut first = 0usize;
            let mut log_sum = 0.0_f64;
            for (k, &cls) in part.class_of.iter().enumerate() {
                if cls == a {
                    if count == 0 {
                        first = k;
                    }
                    log_sum += x_k_hint[k].ln();
                    count += 1;
                }
            }
            let init = if count == 1 {
                x_k_hint[first]
            } else {
                let gm = (log_sum / count as f64).exp();
                if gm.is_finite() && gm > 0.0 {
                    gm.min(1.0)
                } else {
                    x_k_hint[first]
                }
            };
            *slot = Dd::from_f64(init);
        }

        let mut jac = vec![Dd::ZERO; ng * ng];
        let mut rhs = vec![Dd::ZERO; ng];
        let mut last_max_res = f64::INFINITY;
        for _ in 0..24 {
            let mut max_res = 0.0_f64;
            for a in 0..ng {
                let ya = y_g[a];
                let mut sum_wy = Dd::ZERO;
                for b in 0..ng {
                    let wy_b = w_mat[a * ng + b].mul(y_g[b]);
                    sum_wy = sum_wy.add(wy_b);
                    jac[a * ng + b] = ya.mul(wy_b).neg();
                }
                let y_sum = ya.mul(sum_wy);
                let res = Dd::ONE.sub(ya).sub(y_sum);
                if !res.is_finite() {
                    return Err(Error::NotConverged(
                        "non-finite residual in high-precision association site solver",
                    ));
                }
                max_res = max_res.max(res.0.abs());
                rhs[a] = res.neg();
                jac[a * ng + a] = jac[a * ng + a].sub(ya.add(y_sum));
            }
            last_max_res = max_res;
            if max_res <= 1.0e-28 {
                break;
            }
            let Some(du) = solve_linear_dd(ng, &mut jac, &mut rhs) else {
                return Err(Error::NotConverged(
                    "singular Jacobian in high-precision association site solver",
                ));
            };
            let mut max_du = 0.0_f64;
            for d in &du {
                if !d.is_finite() {
                    return Err(Error::NotConverged(
                        "non-finite step in high-precision association site solver",
                    ));
                }
                max_du = max_du.max(d.0.abs());
            }
            let scale = if max_du > 2.0 { 2.0 / max_du } else { 1.0 };
            for a in 0..ng {
                let step = if scale < 1.0 {
                    du[a].mul_f64(scale)
                } else {
                    du[a]
                };
                let next = y_g[a].mul(step.exp());
                if !next.is_finite() || next.0 <= 0.0 {
                    return Err(Error::NotConverged(
                        "invalid site fraction update in high-precision association solver",
                    ));
                }
                y_g[a] = if next.0 > 1.0 { Dd::ONE } else { next };
            }
            if max_du <= 1.0e-28 && max_res <= 1.0e-20 {
                break;
            }
        }
        if !last_max_res.is_finite() || last_max_res > 1.0e-12 {
            return Err(Error::NotConverged(
                "high-precision association site solver did not converge",
            ));
        }
    }

    let x_k: Vec<Dd> = part.class_of.iter().map(|&a| y_g[a]).collect();

    // Certify every site fraction bound and every original ns-site Wertheim
    // residual before forming dF_assoc/dV.
    for k in 0..ns {
        let xk = x_k[k];
        if !xk.is_finite() || xk.0 <= 0.0 || xk.0 > 1.0 {
            return Err(Error::NotConverged(
                "site fraction out of domain in high-precision association solver",
            ));
        }
        let mut sum_cx = Dd::ZERO;
        for l in 0..ns {
            sum_cx = sum_cx.add(c_mat[k * ns + l].mul(x_k[l]));
        }
        let res = Dd::ONE.sub(xk).sub(xk.mul(sum_cx));
        if !res.is_finite() || res.0.abs() > 1.0e-12 {
            return Err(Error::NotConverged(
                "high-precision association site residual exceeds tolerance",
            ));
        }
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
    if !f_v.is_finite() {
        return Err(Error::NotConverged(
            "non-finite high-precision association volume derivative",
        ));
    }
    Ok(f_v)
}

struct SitePartition {
    ns: usize,
    class_of: Vec<usize>,
    class_rep: Vec<usize>,
    num_classes: usize,
    w_mat: Vec<Dd>,
    norm_c: Vec<Dd>,
}

#[inline(always)]
fn is_dd_zero(v: Dd) -> bool {
    v.0 == 0.0 && v.1 == 0.0
}

#[inline(always)]
fn canonical_f64_bits(x: f64) -> u64 {
    if x == 0.0 { 0u64 } else { x.to_bits() }
}

/// Solve `w * y^2 + b * y - 1 = 0` for the unique positive root `y in (0, 1]`
/// given `b >= 1` and `w >= 0` in 106-bit double-double precision, using exact
/// binary scaling (`2^p`) so `b^2 + 4*w` never overflows `f64::MAX`.
fn solve_positive_quadratic_dd(b: Dd, w: Dd) -> Result<Dd> {
    if !b.is_finite()
        || !w.is_finite()
        || b.0 < 1.0 - 1.0e-15
        || w.0 < 0.0
        || (w.0 == 0.0 && w.1 < 0.0)
    {
        return Err(Error::InvalidInput(
            "invalid quadratic coefficients in high-precision association solver",
        ));
    }
    if is_dd_zero(w) {
        let y = Dd::ONE.div(b);
        if !y.is_finite() || y.0 <= 0.0 || y.0 > 1.0 + 1.0e-28 {
            return Err(Error::NotConverged(
                "out-of-domain linear site fraction in high-precision association solver",
            ));
        }
        return Ok(if y.0 > 1.0 { Dd::ONE } else { y });
    }
    let scale_ref = b.0.max(w.0.sqrt());
    let disc_sqrt = if scale_ref > 1.0e140 {
        let p = ((scale_ref.to_bits() >> 52) as i32) - 1023;
        let bs = b.scale_pow2(-p);
        let ws = w.scale_pow2(-2 * p);
        bs.mul(bs).add(ws.mul_f64(4.0)).sqrt().scale_pow2(p)
    } else {
        b.mul(b).add(w.mul_f64(4.0)).sqrt()
    };
    let y = Dd::from_f64(2.0).div(b.add(disc_sqrt));
    if !y.is_finite() || y.0 <= 0.0 || y.0 > 1.0 + 1.0e-28 {
        return Err(Error::NotConverged(
            "out-of-domain quadratic site fraction in high-precision association solver",
        ));
    }
    Ok(if y.0 > 1.0 { Dd::ONE } else { y })
}

/// Compute `W_{a, b} - W_{b, a}` by cancelling bitwise-identical nonzero entries
/// across the two sorted class-row multisets before summing the remaining terms,
/// preserving tiny positive trace couplings without subtracting rounded `O(1)` sums.
fn class_weight_diff_dd(part: &SitePartition, a: usize, b: usize) -> Result<Dd> {
    let collect_sorted = |row_cls: usize, col_cls: usize| -> Vec<Dd> {
        let rep = part.class_rep[row_cls];
        let mut terms = Vec::with_capacity(part.ns);
        for (l, &cls) in part.class_of.iter().enumerate() {
            if cls == col_cls {
                let c_kl = part.norm_c[rep * part.ns + l];
                if !is_dd_zero(c_kl) {
                    terms.push(c_kl);
                }
            }
        }
        terms.sort_unstable_by_key(|v| (canonical_f64_bits(v.0), canonical_f64_bits(v.1)));
        terms
    };

    let a_terms = collect_sorted(a, b);
    let b_terms = collect_sorted(b, a);
    let mut only_a = Vec::with_capacity(a_terms.len());
    let mut only_b = Vec::with_capacity(b_terms.len());
    let (mut i, mut j) = (0usize, 0usize);
    while i < a_terms.len() && j < b_terms.len() {
        let ka = (
            canonical_f64_bits(a_terms[i].0),
            canonical_f64_bits(a_terms[i].1),
        );
        let kb = (
            canonical_f64_bits(b_terms[j].0),
            canonical_f64_bits(b_terms[j].1),
        );
        match ka.cmp(&kb) {
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => {
                only_a.push(a_terms[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                only_b.push(b_terms[j]);
                j += 1;
            }
        }
    }
    only_a.extend_from_slice(&a_terms[i..]);
    only_b.extend_from_slice(&b_terms[j..]);

    let mut pos = Dd::ZERO;
    for &v in &only_a {
        pos = pos.add(v);
    }
    let mut neg = Dd::ZERO;
    for &v in &only_b {
        neg = neg.add(v);
    }
    let diff = pos.sub(neg);
    if !diff.is_finite() {
        return Err(Error::InvalidInput(
            "non-finite bipartite balance difference in high-precision association solver",
        ));
    }
    Ok(diff)
}

/// Solve an exact 2-class bipartite system (`W_{00} = W_{11} = 0`) via exact
/// multiset-difference elimination into a single positive quadratic.
fn solve_bipartite_two_class_dd(part: &SitePartition) -> Result<(Dd, Dd)> {
    let w01 = part.w_mat[1];
    let w10 = part.w_mat[2];
    let diff = class_weight_diff_dd(part, 0, 1)?;
    let (y0, y1) = if diff.0 > 0.0 || (diff.0 == 0.0 && diff.1 >= 0.0) {
        // W_01 >= W_10: W_10 * Y_0^2 + (1 + (W_01 - W_10)) * Y_0 - 1 = 0
        let y0 = solve_positive_quadratic_dd(Dd::ONE.add(diff), w10)?;
        let y1 = Dd::ONE.div(Dd::ONE.add(w10.mul(y0)));
        (y0, y1)
    } else {
        // W_10 > W_01: W_01 * Y_1^2 + (1 + (W_10 - W_01)) * Y_1 - 1 = 0
        let y1 = solve_positive_quadratic_dd(Dd::ONE.add(diff.neg()), w01)?;
        let y0 = Dd::ONE.div(Dd::ONE.add(w01.mul(y1)));
        (y0, y1)
    };
    if !y0.is_finite()
        || !y1.is_finite()
        || y0.0 <= 0.0
        || y1.0 <= 0.0
        || y0.0 > 1.0 + 1.0e-28
        || y1.0 > 1.0 + 1.0e-28
    {
        return Err(Error::NotConverged(
            "out-of-domain bipartite site fraction in high-precision association solver",
        ));
    }
    Ok((
        if y0.0 > 1.0 { Dd::ONE } else { y0 },
        if y1.0 > 1.0 { Dd::ONE } else { y1 },
    ))
}

/// Return `Some(hub)` if the `ng x ng` reduced matrix (`ng >= 3`) is star-bipartite:
/// `W_{hub, hub} == 0` and all leaf-leaf entries `W_{a, b} == 0` for `a != hub, b != hub`.
fn find_star_bipartite_hub(ng: usize, w_mat: &[Dd]) -> Option<usize> {
    if ng < 3 {
        return None;
    }
    'hub: for h in 0..ng {
        if !is_dd_zero(w_mat[h * ng + h]) {
            continue;
        }
        for a in 0..ng {
            if a == h {
                continue;
            }
            for b in 0..ng {
                if b == h {
                    continue;
                }
                if !is_dd_zero(w_mat[a * ng + b]) {
                    continue 'hub;
                }
            }
        }
        return Some(h);
    }
    None
}

/// Solve an exact star-bipartite system (`ng >= 3`) by eliminating all leaf classes
/// `Y_b = 1 / (1 + W_{bh} Y_h)` and solving the strictly decreasing 1D hub equation
/// `(1 + (W_{ph} - W_{hp}) Y_h) / (1 + W_{ph} Y_h) - Y_h - sum_{b != h, p} W_{hb} Y_h / (1 + W_{bh} Y_h) = 0`
/// via safeguarded DD log-space Newton-bisection.
fn solve_star_bipartite_dd(part: &SitePartition, hub: usize) -> Result<Vec<Dd>> {
    let ng = part.num_classes;
    let w_mat = &part.w_mat;
    let mut dom = usize::MAX;
    let mut max_w = -1.0_f64;
    for b in 0..ng {
        if b == hub {
            continue;
        }
        let w = w_mat[hub * ng + b].0;
        if w > max_w {
            max_w = w;
            dom = b;
        }
    }
    if dom == usize::MAX {
        return Err(Error::NotConverged(
            "missing leaf class in star-bipartite association solver",
        ));
    }

    let w_hp = w_mat[hub * ng + dom];
    let w_ph = w_mat[dom * ng + hub];
    let d_p = class_weight_diff_dd(part, dom, hub)?; // W_{ph} - W_{hp}

    let eval_g_and_dg_du = |y: Dd| -> Result<(Dd, Dd, f64)> {
        let den_p = Dd::ONE.add(w_ph.mul(y));
        let term_p = Dd::ONE.add(d_p.mul(y)).div(den_p);
        let mut g = term_p.sub(y);
        let mut dg_du = w_hp.mul(y).div(den_p.mul(den_p)).add(y).neg();
        let mut scale = term_p.0.abs().max(y.0.abs());
        for b in 0..ng {
            if b == hub || b == dom {
                continue;
            }
            let w_hb = w_mat[hub * ng + b];
            let w_bh = w_mat[b * ng + hub];
            if is_dd_zero(w_hb) {
                continue;
            }
            let den_b = Dd::ONE.add(w_bh.mul(y));
            let term_b = w_hb.mul(y).div(den_b);
            g = g.sub(term_b);
            dg_du = dg_du.sub(term_b.div(den_b));
            scale = scale.max(term_b.0.abs());
        }
        if !g.is_finite() || !dg_du.is_finite() || !scale.is_finite() {
            return Err(Error::NotConverged(
                "non-finite star-bipartite residual in high-precision association solver",
            ));
        }
        Ok((g, dg_du, scale))
    };

    // Upper bound on y_h: strictly <= 1, and if d_p < 0 then 1 + d_p * y_h > 0 requires y_h < 1 / (-d_p).
    let y_hi = if d_p.0 < 0.0 {
        let cap = Dd::ONE.div(d_p.neg());
        if cap.0 < 1.0 { cap } else { Dd::ONE }
    } else {
        Dd::ONE
    };
    let mut u_hi = y_hi.ln();
    let (g_hi, _, _) = eval_g_and_dg_du(y_hi)?;
    if g_hi.0 > 0.0 {
        return Err(Error::NotConverged(
            "failed to bracket star-bipartite hub root from above",
        ));
    }

    // Start from the 2-class estimate against the dominant leaf, then lower u_lo until g > 0.
    let y_est = solve_positive_quadratic_dd(
        if d_p.0 < 0.0 {
            Dd::ONE.add(d_p.neg())
        } else {
            Dd::ONE
        },
        Dd::from_f64(w_ph.0.max(w_hp.0)),
    )?;
    let mut u_lo = if y_est.0 < y_hi.0 {
        y_est.ln()
    } else {
        u_hi.sub(Dd::ONE)
    };
    for _ in 0..80 {
        let (g_lo, _, _) = eval_g_and_dg_du(u_lo.exp())?;
        if g_lo.0 >= 0.0 {
            break;
        }
        u_lo = u_lo.sub(Dd::from_f64(10.0));
        if u_lo.0 < -708.0 {
            return Err(Error::NotConverged(
                "star-bipartite hub site fraction underflowed",
            ));
        }
    }

    let mut u = u_lo.add(u_hi).mul_f64(0.5);
    for _ in 0..80 {
        let y = u.exp();
        let (g, dg_du, scale) = eval_g_and_dg_du(y)?;
        if g.0.abs() <= 1.0e-28 * scale.max(1.0e-300) {
            break;
        }
        if g.0 > 0.0 {
            u_lo = u;
        } else {
            u_hi = u;
        }
        let width = u_hi.sub(u_lo);
        if width.0 <= 1.0e-28 * u.0.abs().max(1.0) {
            u = u_lo.add(u_hi).mul_f64(0.5);
            break;
        }
        let next = if dg_du.0 < 0.0 {
            let cand = u.sub(g.div(dg_du));
            if cand.is_finite() && cand.0 > u_lo.0 && cand.0 < u_hi.0 {
                cand
            } else {
                u_lo.add(u_hi).mul_f64(0.5)
            }
        } else {
            u_lo.add(u_hi).mul_f64(0.5)
        };
        u = next;
    }

    let y_h = u.exp();
    if !y_h.is_finite() || y_h.0 <= 0.0 || y_h.0 > 1.0 + 1.0e-28 {
        return Err(Error::NotConverged(
            "out-of-domain star-bipartite hub site fraction",
        ));
    }
    let y_h = if y_h.0 > 1.0 { Dd::ONE } else { y_h };
    let mut y_g = vec![Dd::ONE; ng];
    y_g[hub] = y_h;
    for (b, slot) in y_g.iter_mut().enumerate() {
        if b == hub {
            continue;
        }
        let yb = Dd::ONE.div(Dd::ONE.add(w_mat[b * ng + hub].mul(y_h)));
        if !yb.is_finite() || yb.0 <= 0.0 || yb.0 > 1.0 + 1.0e-28 {
            return Err(Error::NotConverged(
                "out-of-domain star-bipartite leaf site fraction",
            ));
        }
        *slot = if yb.0 > 1.0 { Dd::ONE } else { yb };
    }
    Ok(y_g)
}

/// Compute a conservative exact equitable partition of the `ns x ns` scaled
/// interaction matrix `c_mat` (`C_{kl} = rho * m_l * Delta_{kl}`), returning
/// `SitePartition` where `w_mat[a * ng + b]` is the class-aggregated weight
/// `sum_{l in class b} C_{k, l}` shared identically by every site `k` in class `a`.
///
/// Refinement keys preserve the prior class `(class_of[k], sorted_nonzero_triples)`
/// where each triple is `(class_of[l], hi_bits, lo_bits)` of the canonicalized
/// nonzero entry `C_{kl}`. Two sites remain in the same class only when their
/// prior class and sorted multisets of nonzero weights into every target class
/// match identically. This guarantees exact real-sum equality with zero rounding
/// false merges (and ensures each pass can only split existing classes), while
/// deliberately leaving distinct-multiplicity accidental sum coincidences (such
/// as `[2*w]` vs `[w, w]`) in separate classes for the multi-class Newton solver.
fn equitable_site_partition(ns: usize, c_mat: &[Dd]) -> Result<SitePartition> {
    type SitePartitionKey = (usize, Vec<(usize, u64, u64)>);
    let mut class_of = vec![0usize; ns];
    let mut num_classes = 1usize;
    let mut keys: Vec<SitePartitionKey> =
        (0..ns).map(|_| (0usize, Vec::with_capacity(ns))).collect();

    let mut norm_c = vec![Dd::ZERO; ns * ns];
    for (slot, &entry) in norm_c.iter_mut().zip(c_mat.iter()) {
        if !entry.is_finite() {
            return Err(Error::InvalidInput(
                "non-finite entry in association site partition",
            ));
        }
        let norm = Dd::two_sum(entry.0, entry.1);
        if !norm.is_finite() {
            return Err(Error::InvalidInput(
                "non-finite normalized entry in association site partition",
            ));
        }
        *slot = norm;
    }

    for _ in 0..ns {
        for (k, (prior_cls, triples)) in keys.iter_mut().enumerate() {
            *prior_cls = class_of[k];
            triples.clear();
            for (l, &target_cls) in class_of.iter().enumerate() {
                let c_kl = norm_c[k * ns + l];
                if !is_dd_zero(c_kl) {
                    triples.push((
                        target_cls,
                        canonical_f64_bits(c_kl.0),
                        canonical_f64_bits(c_kl.1),
                    ));
                }
            }
            triples.sort_unstable();
        }

        let mut reps: Vec<usize> = Vec::with_capacity(ns);
        let mut next_class_of = vec![0usize; ns];
        for k in 0..ns {
            if let Some(pos) = reps.iter().position(|&rep| keys[rep] == keys[k]) {
                next_class_of[k] = pos;
            } else {
                next_class_of[k] = reps.len();
                reps.push(k);
            }
        }

        let next_num = reps.len();
        class_of = next_class_of;
        if next_num == num_classes {
            break;
        }
        num_classes = next_num;
    }

    let mut class_rep = vec![usize::MAX; num_classes];
    for (k, &cls) in class_of.iter().enumerate() {
        if class_rep[cls] == usize::MAX {
            class_rep[cls] = k;
        }
    }

    let mut scratch = Vec::with_capacity(ns);
    let mut w_mat = vec![Dd::ZERO; num_classes * num_classes];
    for a in 0..num_classes {
        let rep = class_rep[a];
        for b in 0..num_classes {
            scratch.clear();
            for (l, &cls) in class_of.iter().enumerate() {
                if cls == b {
                    let c_kl = norm_c[rep * ns + l];
                    if !is_dd_zero(c_kl) {
                        scratch.push(c_kl);
                    }
                }
            }
            scratch.sort_unstable_by_key(|v| (canonical_f64_bits(v.0), canonical_f64_bits(v.1)));
            let mut acc = Dd::ZERO;
            for &v in &scratch {
                acc = acc.add(v);
            }
            if !acc.is_finite() {
                return Err(Error::InvalidInput(
                    "non-finite aggregated association weight in high-precision evaluation",
                ));
            }
            w_mat[a * num_classes + b] = acc;
        }
    }

    Ok(SitePartition {
        ns,
        class_of,
        class_rep,
        num_classes,
        w_mat,
        norm_c,
    })
}

fn solve_linear_dd(n: usize, a: &mut [Dd], b: &mut [Dd]) -> Option<Vec<Dd>> {
    for col in 0..n {
        let mut piv = col;
        let mut max_val = 0.0_f64;
        for row in col..n {
            let entry = a[row * n + col];
            if !entry.is_finite() {
                return None;
            }
            let v = entry.0.abs();
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
            if !fac.is_finite() {
                return None;
            }
            if fac.0 == 0.0 && fac.1 == 0.0 {
                continue;
            }
            for j in (col + 1)..n {
                let updated = a[row * n + j].sub(fac.mul(a[col * n + j]));
                if !updated.is_finite() {
                    return None;
                }
                a[row * n + j] = updated;
            }
            let updated_b = b[row].sub(fac.mul(b[col]));
            if !updated_b.is_finite() {
                return None;
            }
            b[row] = updated_b;
        }
    }
    let mut x = vec![Dd::ZERO; n];
    for row in (0..n).rev() {
        let mut s = b[row];
        if !s.is_finite() {
            return None;
        }
        for j in (row + 1)..n {
            s = s.sub(a[row * n + j].mul(x[j]));
        }
        let val = s.div(a[row * n + row]);
        if !val.is_finite() {
            return None;
        }
        x[row] = val;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equitable_partition_merges_only_justified_site_classes() {
        // Pure 4C: 2 donors (0,1) interact with 2 acceptors (2,3) -> 1 class.
        let y = Dd::from_f64(1.8e70);
        let mut c4 = vec![Dd::ZERO; 16];
        for k in 0..2 {
            for l in 2..4 {
                c4[k * 4 + l] = y;
                c4[l * 4 + k] = y;
            }
        }
        let p4 = equitable_site_partition(4, &c4).unwrap();
        assert_eq!(p4.num_classes, 1);
        assert_eq!(p4.class_of, vec![0, 0, 0, 0]);
        assert!((p4.w_mat[0].0 - 3.6e70).abs() <= 1e55);

        // Pure 3B: 2 positive (0,1) interact with 1 negative (2) -> 2 classes.
        let mut c3 = vec![Dd::ZERO; 9];
        for k in 0..2 {
            c3[k * 3 + 2] = y;
            c3[2 * 3 + k] = y;
        }
        let p3 = equitable_site_partition(3, &c3).unwrap();
        assert_eq!(p3.num_classes, 2);
        assert_eq!(p3.class_of[0], p3.class_of[1]);
        assert_ne!(p3.class_of[0], p3.class_of[2]);

        // Asymmetric 4C (0,1 negative; 2,3 positive) + 1ea (4 positive):
        // negative 4C sites couple to site 4 while positive 4C sites do not -> 3 classes.
        let y12 = Dd::from_f64(5.0e10);
        let y21 = Dd::from_f64(2.0e10);
        let mut c5 = vec![Dd::ZERO; 25];
        for k in 0..2 {
            for l in 2..4 {
                c5[k * 5 + l] = y;
                c5[l * 5 + k] = y;
            }
            c5[k * 5 + 4] = y12;
            c5[4 * 5 + k] = y21;
        }
        let p5 = equitable_site_partition(5, &c5).unwrap();
        assert_eq!(p5.num_classes, 3);
        assert_eq!(p5.class_of[0], p5.class_of[1]);
        assert_eq!(p5.class_of[2], p5.class_of[3]);
        assert_ne!(p5.class_of[0], p5.class_of[2]);
        assert_ne!(p5.class_of[0], p5.class_of[4]);
        assert_ne!(p5.class_of[2], p5.class_of[4]);

        // Round 9 counterexample (PRRT_kwDOVJkJcc68MbDL):
        // C01 = C10 = C23 = C32 = Dd(1.0, 2^-54), C12 = C21 = Dd(2^-160, 0.0), others zero.
        // Although the rounded DD row sums of sites 0/3 and 1/2 are bitwise identical,
        // exact multiset partitioning must keep {0, 3} distinct from {1, 2}.
        let w_edge = Dd(1.0, 2.0_f64.powi(-54));
        let w_tiny = Dd(2.0_f64.powi(-160), 0.0);
        let mut c_tiny = vec![Dd::ZERO; 16];
        c_tiny[1] = w_edge;
        c_tiny[4] = w_edge;
        c_tiny[2 * 4 + 3] = w_edge;
        c_tiny[3 * 4 + 2] = w_edge;
        c_tiny[4 + 2] = w_tiny;
        c_tiny[2 * 4 + 1] = w_tiny;
        let pt = equitable_site_partition(4, &c_tiny).unwrap();
        assert_eq!(pt.num_classes, 2);
        assert_eq!(pt.class_of[0], pt.class_of[3]);
        assert_eq!(pt.class_of[1], pt.class_of[2]);
        assert_ne!(pt.class_of[0], pt.class_of[1]);

        // Symmetric 4C (0,1-; 2,3+) + 2B (4-; 5+) merges to 2 classes (one per component).
        let c11 = Dd::from_f64(1.0e30);
        let c12 = Dd::from_f64(3.0e15);
        let c21 = Dd::from_f64(7.0e15);
        let c22 = Dd::from_f64(2.0e10);
        let mut c6 = vec![Dd::ZERO; 36];
        for k in 0..2 {
            for l in 2..4 {
                c6[k * 6 + l] = c11;
                c6[l * 6 + k] = c11;
            }
            c6[k * 6 + 5] = c12;
            c6[5 * 6 + k] = c21;
        }
        for k in 2..4 {
            c6[k * 6 + 4] = c12;
            c6[4 * 6 + k] = c21;
        }
        c6[4 * 6 + 5] = c22;
        c6[5 * 6 + 4] = c22;
        let p6 = equitable_site_partition(6, &c6).unwrap();
        assert_eq!(p6.num_classes, 2);
        assert_eq!(p6.class_of[0], p6.class_of[1]);
        assert_eq!(p6.class_of[0], p6.class_of[2]);
        assert_eq!(p6.class_of[0], p6.class_of[3]);
        assert_eq!(p6.class_of[4], p6.class_of[5]);
        assert_ne!(p6.class_of[0], p6.class_of[4]);

        // Round 10 trace-asymmetric bipartite 4C + 1ea (equal self/cross Delta, x_trace = 1e-40):
        // preserves exact 2-class split ({0,1} vs {2,3,4}) and solves for X_- ~= 5.13e-43,
        // X_+ ~= 5.05e-41 without losing the 1e-40 trace balance to O(1) sum rounding.
        let y_r10 = Dd::from_f64(1.929_832_031_222_762e82);
        let y_trace = y_r10.mul_f64(1.0e-40);
        let mut c_r10 = vec![Dd::ZERO; 25];
        for k in 0..2 {
            for l in 2..4 {
                c_r10[k * 5 + l] = y_r10;
                c_r10[l * 5 + k] = y_r10;
            }
            c_r10[k * 5 + 4] = y_trace;
            c_r10[4 * 5 + k] = y_r10;
        }
        let p_r10 = equitable_site_partition(5, &c_r10).unwrap();
        assert_eq!(p_r10.num_classes, 2);
        assert_eq!(p_r10.class_of, vec![0, 0, 1, 1, 1]);
        let (x_minus, x_plus) = solve_bipartite_two_class_dd(&p_r10).unwrap();
        assert!((x_minus.0 - 5.129_181_126_449_266e-43).abs() <= 1.0e-55);
        assert!((x_plus.0 - 5.051_291_811_264_493e-41).abs() <= 1.0e-53);
    }

    #[test]
    fn dd_finiteness_and_linear_solver_reject_nan_parts() {
        assert!(!Dd(f64::NAN, 0.0).is_finite());
        assert!(!Dd(0.0, f64::NAN).is_finite());
        assert!(!Dd(1.0, f64::INFINITY).is_finite());
        assert!(Dd(1.0, -1.0e-17).is_finite());

        // Non-pivot row carrying a NaN must cause solve_linear_dd to reject rather than mask it.
        let mut a = vec![
            Dd::from_f64(2.0),
            Dd::ZERO,
            Dd(f64::NAN, 0.0),
            Dd::from_f64(3.0),
        ];
        let mut b = vec![Dd::ONE, Dd::ONE];
        assert!(solve_linear_dd(2, &mut a, &mut b).is_none());
    }
}
