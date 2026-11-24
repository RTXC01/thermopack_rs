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
    pub dispersion: H,
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
        let lift = |p: [f64; 3]| p[0] + dt * (p[1] + dt * (0.5 * p[2]));
        let single = self.model.reference == HardSphereReference::SingleComponent;
        let mut mean_d = zero;
        let mut mean_sigma = zero;
        let mut eta = zero;
        let mut zeta_bar = zero;
        for i in 0..nc {
            for j in 0..nc {
                let k = i * nc + j;
                let w = xs[i] * xs[j] * (PI * N_A / 6.) * rho * segments;
                eta += w * d[k] * d[k] * d[k];
                let sigma = lift(self.effective[k].0);
                zeta_bar += w * sigma * sigma * sigma;
                if single {
                    mean_d += xs[i] * xs[j] * d[k];
                    mean_sigma += xs[i] * xs[j] * sigma;
                }
            }
        }
        if single {
            eta = segments * rho * (PI * N_A / 6.) * mean_d.pow_int(3);
            zeta_bar = segments * rho * (PI * N_A / 6.) * mean_sigma.pow_int(3);
        }
        if !(0.0..1.0).contains(&eta.f0) {
            return Err(Error::InvalidInput(
                "density exceeds maximum packing fraction",
            ));
        }
        let beta = 1. / t;
        // 1. Hard-sphere term (BMCSL for mixture, Carnahan-Starling for pure fluid)
        let mut extra = zero;
        let mut a_hs = if self.model.reference == HardSphereReference::NonAdditive {
            santos(&x, &d, rho, segments)
        } else if nc == 1 || single {
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
            if matches!(
                self.model.reference,
                HardSphereReference::Additive | HardSphereReference::AdditiveExact
            ) {
                for i in 0..nc {
                    for j in i + 1..nc {
                        let di = d[i * nc + i];
                        let dj = d[j * nc + j];
                        let additive = 0.5 * (di + dj);
                        let mu = di * dj / (di + dj);
                        let h = z2 * mu / em3;
                        let g = (1. + 3. * h + 2. * h * h) / em3;
                        extra += 4.
                            * PI
                            * N_A
                            * rho
                            * x[i]
                            * x[j]
                            * g
                            * additive
                            * additive
                            * (lift(self.exact_d[i * nc + j]) - additive);
                    }
                }
            }
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
                let dij = if single { mean_d } else { d[idx] };
                let sigma = self.model.sigma_ij[idx];
                let x0 = sigma / dij;
                let sigma_eff = lift(self.effective[idx].0);
                let eps_eff = lift(self.effective[idx].1);
                let y = sigma_eff / dij;
                let integrate = |lambda| {
                    let b = if y.f0 > 1. {
                        calc_b_tilde(y, eta, lambda, eps)
                    } else {
                        zero
                    };
                    x0.pow_real(lambda) * (calc_a1_sutherland(eta, lambda, eps) + b)
                };
                let order = self.model.components[i].fh_order as usize;
                let dfh = self.model.d_feynhibbs_param_ij[idx] / (sigma * sigma) / t;
                let qa = quantum_factors(la, dfh);
                let qr = quantum_factors(lr, dfh);
                let mut a1_tilde = zero;
                let mut a22_tilde = zero;
                let mut alpha = zero;
                let ratio = sigma / sigma_eff;
                // Integrate each inverse power of the FH potential and of its square.
                for k in 0..=order {
                    let a = la + 2. * k as f64;
                    let r = lr + 2. * k as f64;
                    a1_tilde += c_fac * (qa[k] * integrate(a) - qr[k] * integrate(r));
                    if order > 0 {
                        alpha += qa[k] * ratio.pow_real(a) / (a - 3.)
                            - qr[k] * ratio.pow_real(r) / (r - 3.);
                    }
                    for l in 0..=order {
                        let offset = 2. * l as f64;
                        a22_tilde += qa[k] * qa[l] * integrate(a + la + offset)
                            + qr[k] * qr[l] * integrate(r + lr + offset)
                            - 2. * qa[k] * qr[l] * integrate(a + lr + offset);
                    }
                }
                let f_alpha = if order == 0 {
                    self.model.f_alpha_ij[idx].map(|v| zero + v)
                } else {
                    alpha_functions(alpha * (c_fac * eps) / eps_eff)
                };
                let eta_ij = (PI * N_A / 6.) * segments * rho * dij * dij * dij;
                let a1_ij = eta_ij * a1_tilde;
                let zb2 = zeta_bar * zeta_bar;
                let chi = zeta_bar
                    * (f_alpha[0] + f_alpha[1] * (zb2 * zb2) + f_alpha[2] * zeta_bar.pow_int(7));
                let a2_ij = 0.5 * eps * (c_fac * c_fac) * khs * (1.0 + chi) * eta_ij * a22_tilde;

                // a3_ij
                let a3_ij = -(eps_eff * eps_eff * eps_eff)
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
            if self.model.terms.chain && (mi - 1.0).abs() > 1e-12 {
                let idx = i * nc + i;
                let la = self.model.lambda_a_ij[idx];
                let lr = self.model.lambda_r_ij[idx];
                let eps = self.model.eps_divk_ij[idx];
                let c_fac = self.model.c_ij[idx];
                let x0 = self.model.sigma_ij[idx] / if single { mean_d } else { d[idx] };
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
                let g1 = if self.model.terms.first {
                    g11 + g12
                } else {
                    zero
                };

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
                let g2 = if self.model.terms.second {
                    (1.0 + gam_c) * (g21 + g22)
                } else {
                    zero
                };

                let beta_eps = eps / t;
                let w = (beta_eps * g1 + (beta_eps * beta_eps) * g2) / gd;
                let log_g_contact = exparg + w;

                a_chain -= x[i] * (mi - 1.0) * log_g_contact;
            }
        }

        let terms = self.model.terms;
        if !terms.hard_sphere {
            a_hs = zero;
        }
        if !terms.first {
            a1_mix = zero;
        }
        if !terms.second {
            a2_mix = zero;
        }
        if !terms.third {
            a3_mix = zero;
        }
        let mut dispersion = segments * beta * (a1_mix + beta * (a2_mix + beta * a3_mix));
        if let Some((rc, shift)) = self.model.cutoff {
            let la = self.model.lambda_a_ij[0];
            let lr = self.model.lambda_r_ij[0];
            let sigma = self.model.sigma_ij[0];
            let dfh = self.model.d_feynhibbs_param_ij[0] / sigma.powi(2) / t;
            let tail = |lambda: f64| {
                quantum_factors(lambda, dfh)
                    .iter()
                    .take(1 + self.model.components[0].fh_order as usize)
                    .enumerate()
                    .fold(zero, |sum, (k, &q)| {
                        let exponent = lambda + 2. * k as f64 - 3.;
                        sum + q
                            * rc.powf(-exponent)
                            * (1. / exponent + if shift { 1. / 3. } else { 0. })
                    })
            };
            dispersion -= 2.
                * PI
                * N_A
                * self.model.c_ij[0]
                * sigma.powi(3)
                * self.model.eps_divk_ij[0]
                * rho
                * beta
                * (tail(lr) - tail(la));
        }
        Ok(Terms {
            hs: a_hs,
            a1: a1_mix,
            a2: a2_mix,
            a3: a3_mix,
            chain: a_chain,
            dispersion,
            total: sum * (segments * a_hs + dispersion + a_chain + extra),
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

fn quantum_factors(lambda: f64, dfh: H) -> [H; 3] {
    let first = lambda * (lambda - 1.) * dfh;
    [
        H::with_order(1., dfh.order),
        first,
        0.5 * (lambda + 1.) * (lambda + 2.) * dfh * first,
    ]
}

fn alpha_functions(alpha: H) -> [H; 6] {
    std::array::from_fn(|k| {
        let num = PHI_MAT[0][k]
            + alpha * (PHI_MAT[1][k] + alpha * (PHI_MAT[2][k] + alpha * PHI_MAT[3][k]));
        let den = 1. + alpha * (PHI_MAT[4][k] + alpha * (PHI_MAT[5][k] + alpha * PHI_MAT[6][k]));
        num / den
    })
}

// Santos second and third virial coefficients in units of angstrom^3
// and angstrom^6. The cyclic B_ijk terms give identical sums over x_i x_j x_k.
fn santos(x: &[H], d: &[H], rho: H, segments: H) -> H {
    let nc = x.len();
    let zero = H::with_order(0., rho.order);
    let diameter = |i: usize, j: usize| d[i * nc + j] * 1e10;
    let mut moment = zero;
    let mut b2 = zero;
    let mut b3 = zero;
    for i in 0..nc {
        moment += x[i] * diameter(i, i).pow_int(3);
        for j in 0..nc {
            let dij = diameter(i, j);
            let w = x[i] * x[j] * dij.pow_int(3) * 4.;
            b2 += w;
            for k in 0..nc {
                let dik = diameter(i, k);
                let djk = diameter(j, k);
                let clip = |v: H| if v.f0 > 0. { v } else { zero };
                let dkij = clip(dik + djk - dij);
                let ckij = dkij
                    * dkij
                    * (dkij + 1.5 * clip(dij + dik - djk) * clip(dij + djk - dik) / dij);
                b3 += w * x[k] * ckij;
            }
        }
    }
    let eta = (PI * N_A * 1e-30 / 6.) * rho * moment;
    let om = 1. - eta;
    let b2 = b2 / moment;
    let b3 = b3 / (moment * moment);
    (-(10. * b2 - 4. * b3) * om.ln() + (b3 - b2) * (4. - 3. * eta) * eta / (om * om))
        / (6. * segments)
}
