//! Directional evaluation of the non-associating residual terms.
//! All extensive coordinates are seeded at unit total moles before scaling.
use super::*;
use crate::hyperdual::HyperDual as H;

pub(super) struct Terms {
    pub hs: H,
    pub a1: H,
    pub a2: H,
    pub a3: H,
    pub chain: H,
    pub total: H,
}

impl PreparedSaftVrMie<'_> {
    pub(super) fn kernel(
        &self,
        rho_value: f64,
        dir1: Option<usize>,
        dir2: Option<usize>,
    ) -> Result<Terms> {
        let nc = self.model.num_components();
        let order = dir1.is_some() as u8 + dir2.is_some() as u8;
        let seed = |value, index| {
            let mut z = H::with_order(value, order);
            z.f1 = f64::from(dir1 == Some(index));
            z.f2 = f64::from(dir2 == Some(index));
            z
        };
        let zero = H::with_order(0., order);
        let t = seed(self.temperature, 0);
        let v = seed(1. / rho_value, 1);
        let amounts: Vec<_> = self
            .mole_fractions
            .iter()
            .enumerate()
            .map(|(i, &n)| seed(n, 2 + i))
            .collect();
        let sum = amounts.iter().fold(zero, |acc, &n| acc + n);
        let x: Vec<_> = amounts.iter().map(|&n| n / sum).collect();
        let segments = x
            .iter()
            .zip(&self.model.components)
            .fold(zero, |a, (&x, c)| a + x * c.m);
        let xs: Vec<_> = x
            .iter()
            .zip(&self.model.components)
            .map(|(&x, c)| x * c.m / segments)
            .collect();
        let rho = sum / v;
        let dt = t - self.temperature;
        let d: Vec<_> = (0..nc * nc)
            .map(|k| self.d[k] + dt * (self.d_t[k] + dt * (0.5 * self.d_tt[k])))
            .collect();
        let mut eta = zero;
        let mut zeta_bar = zero;
        for i in 0..nc {
            for j in 0..nc {
                let k = i * nc + j;
                let w = xs[i] * xs[j] * (PI * N_A / 6.) * rho * segments;
                eta += w * d[k] * d[k] * d[k];
                zeta_bar += w * self.model.sigma_ij[k].powi(3);
            }
        }
        if !(0.0..1.0).contains(&eta.f0) {
            return Err(Error::InvalidInput(
                "density exceeds maximum packing fraction",
            ));
        }
        let beta = 1. / t;
        // 1. Hard-sphere term (BMCSL for mixture, Carnahan-Starling for pure fluid)
        let a_hs = if nc == 1 {
            let om = 1.0 - eta;
            (4.0 * eta - 3.0 * eta * eta) / (om * om)
        } else {
            let mut z0 = zero;
            let mut z1 = zero;
            let mut z2 = zero;
            let mut z3 = zero;
            for i in 0..nc {
                let di = d[i * nc + i];
                let factor = (PI * N_A / 6.0) * rho * self.model.components[i].m * x[i];
                let di2 = di * di;
                z0 += factor;
                z1 += factor * di;
                z2 += factor * di2;
                z3 += factor * (di2 * di);
            }
            let em3 = 1.0 - z3;
            let pref = 1.0 / z0;
            let z2_3 = z2 * z2 * z2;
            pref * (3.0 * z1 * z2 / em3
                + z2_3 / (z3 * em3 * em3)
                + (z2_3 / (z3 * z3) - z0) * em3.ln())
        };

        // 2. Dispersion A1, A2, A3
        let mut a1_mix = zero;
        let mut a2_mix = zero;
        let mut a3_mix = zero;

        // Carnahan-Starling isothermal compressibility Khs
        let eta2 = eta * eta;
        let eta3 = eta2 * eta;
        let eta4 = eta2 * eta2;
        let om_eta = 1.0 - eta;
        let om_eta2 = om_eta * om_eta;
        let denum_khs = 1.0 + 4.0 * eta + 4.0 * eta2 - 4.0 * eta3 + eta4;
        let khs = (om_eta2 * om_eta2) / denum_khs;

        for i in 0..nc {
            for j in 0..nc {
                let idx = i * nc + j;
                let la = self.model.lambda_a_ij[idx];
                let lr = self.model.lambda_r_ij[idx];
                let eps = self.model.eps_divk_ij[idx];
                let c_fac = self.model.c_ij[idx];
                let x0 = self.model.sigma_ij[idx] / d[idx];
                let f_alpha = &self.model.f_alpha_ij[idx];

                // a1_ij
                let asr = calc_a1_sutherland(eta, lr, eps);
                let asa = calc_a1_sutherland(eta, la, eps);
                let (br, ba) = if x0.f0 > 1.0 {
                    (
                        calc_b_tilde(x0, eta, lr, eps),
                        calc_b_tilde(x0, eta, la, eps),
                    )
                } else {
                    (zero, zero)
                };
                let a1_tilde =
                    c_fac * (x0.pow_real(la) * (asa + ba) - x0.pow_real(lr) * (asr + br));
                let dij = d[idx];
                let eta_ij = (PI * N_A / 6.0) * segments * rho * (dij * dij * dij);
                let a1_ij = eta_ij * a1_tilde;

                // a2_ij
                let lam_2r = 2.0 * lr;
                let lam_2a = 2.0 * la;
                let lam_ar = la + lr;
                let as_2r = calc_a1_sutherland(eta, lam_2r, eps);
                let as_2a = calc_a1_sutherland(eta, lam_2a, eps);
                let as_ar = calc_a1_sutherland(eta, lam_ar, eps);
                let (b_2r, b_2a, b_ar) = if x0.f0 > 1.0 {
                    (
                        calc_b_tilde(x0, eta, lam_2r, eps),
                        calc_b_tilde(x0, eta, lam_2a, eps),
                        calc_b_tilde(x0, eta, lam_ar, eps),
                    )
                } else {
                    (zero, zero, zero)
                };
                let a22_tilde = x0.pow_real(lam_2a) * (as_2a + b_2a)
                    - 2.0 * x0.pow_real(lam_ar) * (as_ar + b_ar)
                    + x0.pow_real(lam_2r) * (as_2r + b_2r);
                let zb2 = zeta_bar * zeta_bar;
                let chi = zeta_bar
                    * (f_alpha[0] + f_alpha[1] * (zb2 * zb2) + f_alpha[2] * zeta_bar.pow_int(7));
                let a2_ij = 0.5 * eps * (c_fac * c_fac) * khs * (1.0 + chi) * eta_ij * a22_tilde;

                // a3_ij
                let a3_ij = -(eps * eps * eps)
                    * f_alpha[3]
                    * zeta_bar
                    * (f_alpha[4] * zeta_bar + f_alpha[5] * zb2).exp();

                let weight = xs[i] * xs[j];
                a1_mix += weight * a1_ij;
                a2_mix += weight * a2_ij;
                a3_mix += weight * a3_ij;
            }
        }

        // 3. Chain term
        let mut a_chain = zero;
        for i in 0..nc {
            let mi = self.model.components[i].m;
            if (mi - 1.0).abs() > 1e-12 {
                let idx = i * nc + i;
                let la = self.model.lambda_a_ij[idx];
                let lr = self.model.lambda_r_ij[idx];
                let eps = self.model.eps_divk_ij[idx];
                let c_fac = self.model.c_ij[idx];
                let x0 = self.model.sigma_ij[idx] / d[idx];
                let alpha = self.model.alpha_ij[idx];

                // gd
                let one_m_eta = 1.0 - eta;
                let d2 = one_m_eta * one_m_eta;
                let d3 = d2 * one_m_eta;
                let k0 = -one_m_eta.ln()
                    + (42.0 * eta - 39.0 * eta2 + 9.0 * eta3 - 2.0 * eta4) / (6.0 * d3);
                let k1 = (eta4 + 6.0 * eta2 - 12.0 * eta) / (2.0 * d3);
                let k2 = (-3.0 * eta2) / (8.0 * d2);
                let k3 = (-eta4 + 3.0 * eta2 + 3.0 * eta) / (6.0 * d3);
                let x0_2 = x0 * x0;
                let exparg = k0 + k1 * x0 + k2 * x0_2 + k3 * (x0_2 * x0);
                let gd = exparg.exp();

                // g1
                let (asr, asr_e) = calc_a1_sutherland_with_deriv(eta, lr, eps);
                let (asa, asa_e) = calc_a1_sutherland_with_deriv(eta, la, eps);
                let (br, br_e, ba, ba_e) = if x0.f0 > 1.0 {
                    let (br, br_e) = calc_b_tilde_with_deriv(x0, eta, lr, eps);
                    let (ba, ba_e) = calc_b_tilde_with_deriv(x0, eta, la, eps);
                    (br, br_e, ba, ba_e)
                } else {
                    (zero, zero, zero, zero)
                };
                let a1_tilde =
                    c_fac * (x0.pow_real(la) * (asa + ba) - x0.pow_real(lr) * (asr + br));
                let a1_tilde_e =
                    c_fac * (x0.pow_real(la) * (asa_e + ba_e) - x0.pow_real(lr) * (asr_e + br_e));
                let g11 = (1.0 / (4.0 * eps)) * (a1_tilde + eta * a1_tilde_e);
                let g12 = (c_fac / (12.0 * eps))
                    * (-la * x0.pow_real(la) * (asa + ba) + lr * x0.pow_real(lr) * (asr + br));
                let g1 = g11 + g12;

                // g2
                let lam_2r = 2.0 * lr;
                let lam_2a = 2.0 * la;
                let lam_ar = la + lr;
                let (as_2r, as_2r_e) = calc_a1_sutherland_with_deriv(eta, lam_2r, eps);
                let (as_2a, as_2a_e) = calc_a1_sutherland_with_deriv(eta, lam_2a, eps);
                let (as_ar, as_ar_e) = calc_a1_sutherland_with_deriv(eta, lam_ar, eps);
                let (b_2r, b_2r_e, b_2a, b_2a_e, b_ar, b_ar_e) = if x0.f0 > 1.0 {
                    let (b_2r, b_2r_e) = calc_b_tilde_with_deriv(x0, eta, lam_2r, eps);
                    let (b_2a, b_2a_e) = calc_b_tilde_with_deriv(x0, eta, lam_2a, eps);
                    let (b_ar, b_ar_e) = calc_b_tilde_with_deriv(x0, eta, lam_ar, eps);
                    (b_2r, b_2r_e, b_2a, b_2a_e, b_ar, b_ar_e)
                } else {
                    (zero, zero, zero, zero, zero, zero)
                };
                let a22_tilde = x0.pow_real(lam_2a) * (as_2a + b_2a)
                    - 2.0 * x0.pow_real(lam_ar) * (as_ar + b_ar)
                    + x0.pow_real(lam_2r) * (as_2r + b_2r);
                let a22_tilde_e = x0.pow_real(lam_2a) * (as_2a_e + b_2a_e)
                    - 2.0 * x0.pow_real(lam_ar) * (as_ar_e + b_ar_e)
                    + x0.pow_real(lam_2r) * (as_2r_e + b_2r_e);
                let num_khs = d2 * d2;
                let num_khs_d = -4.0 * d3;
                let den_khs_d = 4.0 + 8.0 * eta - 12.0 * eta2 + 4.0 * eta3;
                let khs_e = (num_khs_d * denum_khs - num_khs * den_khs_d) / (denum_khs * denum_khs);
                let c_fac2 = c_fac * c_fac;
                let g21 = (c_fac2 / (8.0 * eps))
                    * (khs * (a22_tilde + eta * a22_tilde_e) + eta * khs_e * a22_tilde);
                let g22 = khs
                    * (c_fac2 / (12.0 * eps))
                    * (-la * x0.pow_real(lam_2a) * (as_2a + b_2a)
                        + lam_ar * x0.pow_real(lam_ar) * (as_ar + b_ar)
                        - lr * x0.pow_real(lam_2r) * (as_2r + b_2r));

                let theta = (eps / t).exp() - 1.0;
                let gam_c = theta
                    * 10.0
                    * (1.0 - (10.0 * (0.57 - alpha)).tanh())
                    * zeta_bar
                    * (-6.7 * zeta_bar - 8.0 * (zeta_bar * zeta_bar)).exp();
                let g2 = (1.0 + gam_c) * (g21 + g22);

                let beta_eps = eps / t;
                let w = (beta_eps * g1 + (beta_eps * beta_eps) * g2) / gd;
                let log_g_contact = exparg + w;

                a_chain -= x[i] * (mi - 1.0) * log_g_contact;
            }
        }

        let am = a_hs + beta * a1_mix + beta * beta * a2_mix + (beta * beta * beta) * a3_mix;
        Ok(Terms {
            hs: a_hs,
            a1: a1_mix,
            a2: a2_mix,
            a3: a3_mix,
            chain: a_chain,
            total: sum * (segments * am + a_chain),
        })
    }
}

fn calc_eff_eta(eta: H, lam: f64) -> H {
    let inv_lam = [1.0, 1.0 / lam, 1.0 / (lam * lam), 1.0 / (lam * lam * lam)];
    let mut c_vec = [0.0; 4];
    for j in 0..4 {
        for i in 0..4 {
            c_vec[j] += LAM_COEFF[i][j] * inv_lam[i];
        }
    }
    eta * (c_vec[0] + eta * (c_vec[1] + eta * (c_vec[2] + eta * c_vec[3])))
}

fn calc_eff_eta_with_deriv(eta: H, lam: f64) -> (H, H) {
    let inv_lam = [1.0, 1.0 / lam, 1.0 / (lam * lam), 1.0 / (lam * lam * lam)];
    let mut c_vec = [0.0; 4];
    for j in 0..4 {
        for i in 0..4 {
            c_vec[j] += LAM_COEFF[i][j] * inv_lam[i];
        }
    }
    let ef = eta * (c_vec[0] + eta * (c_vec[1] + eta * (c_vec[2] + eta * c_vec[3])));
    let ef_e = c_vec[0] + eta * (2.0 * c_vec[1] + eta * (3.0 * c_vec[2] + 4.0 * eta * c_vec[3]));
    (ef, ef_e)
}

fn calc_a1_sutherland(eta: H, lam: f64, eps_val: f64) -> H {
    let ef = calc_eff_eta(eta, lam);
    let pref = 6.0 * eps_val / (lam - 3.0);
    let om = 1.0 - ef;
    let denum = om * om * om;
    -2.0 * pref * (1.0 - 0.5 * ef) / denum
}

fn calc_a1_sutherland_with_deriv(eta: H, lam: f64, eps_val: f64) -> (H, H) {
    let (ef, ef_e) = calc_eff_eta_with_deriv(eta, lam);
    let pref = 6.0 * eps_val / (lam - 3.0);
    let om = 1.0 - ef;
    let denum = om * om * om;
    let as_val = -2.0 * pref * (1.0 - 0.5 * ef) / denum;
    let temp1 = -pref * (5.0 - 2.0 * ef) / denum / (1.0 - ef);
    let as_e = temp1 * ef_e;
    (as_val, as_e)
}

fn calc_b_tilde(x0: H, eta: H, lam: f64, eps_val: f64) -> H {
    let i_term = -(x0.pow_real(3.0 - lam) - 1.0) / (lam - 3.0);
    let j_term =
        -((lam - 3.0) * x0.pow_real(4.0 - lam) - (lam - 4.0) * x0.pow_real(3.0 - lam) - 1.0)
            / ((lam - 3.0) * (lam - 4.0));
    let om = 1.0 - eta;
    let d3 = om * om * om;
    let keta10 = (2.0 - eta) / d3;
    let keta20 = -9.0 * eta * (1.0 + eta) / d3;
    6.0 * eps_val * (keta10 * i_term + keta20 * j_term)
}

fn calc_b_tilde_with_deriv(x0: H, eta: H, lam: f64, eps_val: f64) -> (H, H) {
    let i_term = -(x0.pow_real(3.0 - lam) - 1.0) / (lam - 3.0);
    let j_term =
        -((lam - 3.0) * x0.pow_real(4.0 - lam) - (lam - 4.0) * x0.pow_real(3.0 - lam) - 1.0)
            / ((lam - 3.0) * (lam - 4.0));
    let om = 1.0 - eta;
    let d3 = om * om * om;
    let d4 = (1.0 - eta) * d3;
    let keta10 = (2.0 - eta) / d3;
    let keta11 = (5.0 - 2.0 * eta) / d4;
    let keta20 = -9.0 * eta * (1.0 + eta) / d3;
    let keta21 = -9.0 * (eta * eta + 4.0 * eta + 1.0) / d4;
    let b_val = 6.0 * eps_val * (keta10 * i_term + keta20 * j_term);
    let b_e = 6.0 * eps_val * (keta11 * i_term + keta21 * j_term);
    (b_val, b_e)
}
