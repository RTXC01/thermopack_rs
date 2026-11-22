//! Multipolar (dipole-dipole, quadrupole-quadrupole, dipole-quadrupole)
//! perturbation Helmholtz energy ported from ThermoPack's `src/multipol_var.f90`
//! and `src/multipol.f90`.
//!
//! References:
//! - Gross (2005), AIChE J. 51(9), 2556–2568 (`10.1002/aic.10502`)
//! - Gross & Vrabec (2006), AIChE J. 52(3), 1194–1204 (`10.1002/aic.10683`)
//! - Vrabec & Gross (2008), J. Phys. Chem. B 112(1), 51–60 (`10.1021/jp072619u`)

use std::f64::consts::PI;

use crate::error::positive;
use crate::hyperdual::HyperDual;
use crate::{Error, Result};

/// Boltzmann constant [J/K] matching `thermopack_constants::kB_const`.
const KB_CONST: f64 = 1.380_649e-23;
/// Avogadro constant [1/mol] matching `thermopack_constants::N_Avogadro`.
const N_AVOGADRO: f64 = crate::N_A;
/// Small threshold matching `numconstants::small`.
const SMALL: f64 = 1.0e-12;
/// Dipole-quadrupole three-body constant matching `multipol::ALPHA`.
const ALPHA_DQ: f64 = 1.193_735_0;

// Dipole-dipole constants (Gross & Vrabec 2006), indexed [row 0..3][col 0..5]
const AD: [[f64; 5]; 3] = [
    [
        0.30435038064,
        -0.13585877707,
        1.44933285154,
        0.35569769252,
        -2.06533084541,
    ],
    [
        0.95346405973,
        -1.83963831920,
        2.01311801180,
        -7.37249576667,
        8.23741345333,
    ],
    [
        -1.16100802773,
        4.52586067320,
        0.97512223853,
        -12.2810377713,
        5.93975747420,
    ],
];

const BD: [[f64; 5]; 3] = [
    [0.21879385627, -1.18964307357, 1.16268885692, 0.0, 0.0],
    [-0.58731641193, 1.24891317047, -0.50852797392, 0.0, 0.0],
    [3.48695755800, -14.9159739347, 15.3720218600, 0.0, 0.0],
];

const CD: [[f64; 4]; 3] = [
    [-0.06467735252, 0.19758818347, -0.80875619458, 0.69028490492],
    [
        -0.95208758351,
        2.99242575222,
        -2.38026356489,
        -0.27012609786,
    ],
    [-0.62609792333, 1.29246858189, 1.65427830900, -3.43967436378],
];

// Quadrupole-quadrupole constants (Gross 2005), indexed [row 0..3][col 0..5]
const AQ: [[f64; 5]; 3] = [
    [
        1.237830788,
        2.435503144,
        1.633090469,
        -1.611815241,
        6.977118504,
    ],
    [
        1.285410878,
        -11.46561451,
        22.08689285,
        7.46913832,
        -17.19777208,
    ],
    [
        1.794295401,
        0.769510293,
        7.264792255,
        94.48669892,
        -77.1484579,
    ],
];

const BQ: [[f64; 5]; 3] = [
    [0.454271755, -4.501626435, 3.585886783, 0.0, 0.0],
    [-0.813734006, 10.06402986, -10.87663092, 0.0, 0.0],
    [6.868267516, -5.173223765, -17.2402066, 0.0, 0.0],
];

const CQ: [[f64; 4]; 3] = [
    [-0.500043713, 6.531869153, -16.01477983, 14.42597018],
    [2.000209381, -6.78386584, 20.38324603, -10.89598394],
    [3.135827145, 7.247588801, 3.075947834, 0.0],
];

// Dipole-quadrupole constants (Vrabec & Gross 2008)
const ADQ: [[f64; 4]; 3] = [
    [0.697094963, -0.633554144, 2.945509028, -1.467027314],
    [-0.673459279, -1.425899106, 4.19441392, 1.0266216],
    [0.670340770, -4.338471826, 7.234168360, 0.0],
];

const BDQ: [[f64; 4]; 3] = [
    [-0.484038322, 1.970405465, -2.118572671, 0.0],
    [0.67651011, -3.013867512, 0.46742656, 0.0],
    [-1.167560146, 2.13488432, 0.0, 0.0],
];

const CDQ: [[f64; 3]; 2] = [
    [0.795009692, 3.386863396, 0.475106328],
    [-2.099579397, -5.941376392, -0.178820384],
];

/// Precomputed multipolar parameters mirroring `multipol_var::multipol_param`.
#[derive(Debug, Clone)]
pub struct MultipolParams {
    pub nc: usize,
    pub mu_star_2: Vec<f64>,
    pub q_star_2: Vec<f64>,
    pub m: Vec<f64>,
    pub eps_divk_ij: Vec<f64>,
    pub sigma_ij: Vec<f64>,
    pub sigma_ij_3: Vec<f64>,
    pub sigma_ij_5: Vec<f64>,
    // Precomputed pair/triple correlation coefficients:
    // pair index: (i * nc + j) * len + l
    a_ij_qq: Vec<[f64; 5]>,
    b_ij_qq: Vec<[f64; 5]>,
    c_ijk_qq: Vec<[f64; 4]>,
    a_ij_dd: Vec<[f64; 5]>,
    b_ij_dd: Vec<[f64; 5]>,
    c_ijk_dd: Vec<[f64; 4]>,
    a_ij_dq: Vec<[f64; 4]>,
    b_ij_dq: Vec<[f64; 4]>,
    c_ijk_dq: Vec<[f64; 3]>,
    pub mu_indices: Vec<usize>,
    pub q_indices: Vec<usize>,
    pub enable_qq: bool,
    pub enable_dd: bool,
    pub enable_dq: bool,
}

/// Temperature- and volume-derivative multipolar Helmholtz jet.
#[derive(Debug, Clone, Copy, Default)]
pub struct MultipolPressureJet {
    pub f: f64,
    pub f_t: f64,
    pub f_v: f64,
    pub f_tt: f64,
    pub f_tv: f64,
    pub f_vv: f64,
}

/// Complete 10-block multipolar reduced residual Helmholtz energy (`F = A_mp / (R T)`).
#[derive(Debug, Clone)]
pub struct MultipolHelmholtz {
    pub f: f64,
    pub f_t: f64,
    pub f_v: f64,
    pub f_n: Vec<f64>,
    pub f_tt: f64,
    pub f_tv: f64,
    pub f_tn: Vec<f64>,
    pub f_vv: f64,
    pub f_vn: Vec<f64>,
    /// Row-major `nc * nc` second mole derivative matrix.
    pub f_nn: Vec<f64>,
}

impl MultipolHelmholtz {
    /// All-zero multipolar Helmholtz bundle for `nc` components.
    pub fn zeros(nc: usize) -> Self {
        Self {
            f: 0.0,
            f_t: 0.0,
            f_v: 0.0,
            f_n: vec![0.0; nc],
            f_tt: 0.0,
            f_tv: 0.0,
            f_tn: vec![0.0; nc],
            f_vv: 0.0,
            f_vn: vec![0.0; nc],
            f_nn: vec![0.0; nc * nc],
        }
    }
}

impl MultipolParams {
    /// Initialize multipolar parameters from pure moments (`mu` `[D]`, `q` `[ÅD]`),
    /// segment numbers `m`, and full `nc x nc` pair matrices `sigma_ij` `[m]` and
    /// `eps_divk_ij` `[K]`.
    pub fn new(
        mu: &[f64],
        q: &[f64],
        m: &[f64],
        sigma_ij: &[f64],
        eps_divk_ij: &[f64],
    ) -> Result<Self> {
        let nc = m.len();
        if mu.len() != nc
            || q.len() != nc
            || sigma_ij.len() != nc * nc
            || eps_divk_ij.len() != nc * nc
        {
            return Err(Error::InvalidInput(
                "multipolar parameter dimension mismatch",
            ));
        }

        let mut sigma_ij_3 = Vec::with_capacity(nc * nc);
        let mut sigma_ij_5 = Vec::with_capacity(nc * nc);
        for &s in sigma_ij {
            let s2 = s * s;
            let s3 = s2 * s;
            sigma_ij_3.push(s3);
            sigma_ij_5.push(s3 * s2);
        }

        let mut mu_star_2 = Vec::with_capacity(nc);
        let mut q_star_2 = Vec::with_capacity(nc);
        let mut mu_indices = Vec::new();
        let mut q_indices = Vec::new();

        for i in 0..nc {
            let ii = i * nc + i;
            let ms2 =
                1.0e-49 * (mu[i] * mu[i]) / (m[i] * KB_CONST * eps_divk_ij[ii] * sigma_ij_3[ii]);
            let qs2 =
                1.0e-69 * (q[i] * q[i]) / (m[i] * KB_CONST * eps_divk_ij[ii] * sigma_ij_5[ii]);
            if !ms2.is_finite() || !qs2.is_finite() {
                return Err(Error::InvalidInput("non-finite reduced multipolar moments"));
            }
            mu_star_2.push(ms2);
            q_star_2.push(qs2);
            if mu[i].abs() > SMALL {
                mu_indices.push(i);
            }
            if q[i].abs() > SMALL {
                q_indices.push(i);
            }
        }

        let mut a_ij_qq = vec![[0.0; 5]; nc * nc];
        let mut b_ij_qq = vec![[0.0; 5]; nc * nc];
        let mut c_ijk_qq = vec![[0.0; 4]; nc * nc * nc];
        let mut a_ij_dd = vec![[0.0; 5]; nc * nc];
        let mut b_ij_dd = vec![[0.0; 5]; nc * nc];
        let mut c_ijk_dd = vec![[0.0; 4]; nc * nc * nc];
        let mut a_ij_dq = vec![[0.0; 4]; nc * nc];
        let mut b_ij_dq = vec![[0.0; 4]; nc * nc];
        let mut c_ijk_dq = vec![[0.0; 3]; nc * nc * nc];

        for i in 0..nc {
            let mi = m[i].min(2.0);
            for j in 0..nc {
                let mj = m[j].min(2.0);
                let mij = (mi * mj).sqrt();
                let h = (mij - 1.0) / mij;
                let g = h * (mij - 2.0) / mij;
                let ij = i * nc + j;
                let eps_ij = eps_divk_ij[ij];
                for l in 0..5 {
                    a_ij_qq[ij][l] = AQ[0][l] + AQ[1][l] * h + AQ[2][l] * g;
                    b_ij_qq[ij][l] = (BQ[0][l] + BQ[1][l] * h + BQ[2][l] * g) * eps_ij;
                    a_ij_dd[ij][l] = AD[0][l] + AD[1][l] * h + AD[2][l] * g;
                    b_ij_dd[ij][l] = (BD[0][l] + BD[1][l] * h + BD[2][l] * g) * eps_ij;
                    if l < 4 {
                        a_ij_dq[ij][l] = ADQ[0][l] + ADQ[1][l] * h + ADQ[2][l] * g;
                        b_ij_dq[ij][l] = (BDQ[0][l] + BDQ[1][l] * h + BDQ[2][l] * g) * eps_ij;
                    }
                }
                for (k, &m_k) in m.iter().enumerate() {
                    let mk = m_k.min(2.0);
                    let mijk = (mi * mj * mk).cbrt();
                    let hk = (mijk - 1.0) / mijk;
                    let gk = hk * (mijk - 2.0) / mijk;
                    let ijk = (i * nc + j) * nc + k;
                    for l in 0..4 {
                        c_ijk_qq[ijk][l] = CQ[0][l] + CQ[1][l] * hk + CQ[2][l] * gk;
                        c_ijk_dd[ijk][l] = CD[0][l] + CD[1][l] * hk + CD[2][l] * gk;
                        if l < 3 {
                            c_ijk_dq[ijk][l] = CDQ[0][l] + CDQ[1][l] * hk;
                        }
                    }
                }
            }
        }

        Ok(Self {
            nc,
            mu_star_2,
            q_star_2,
            m: m.to_vec(),
            eps_divk_ij: eps_divk_ij.to_vec(),
            sigma_ij: sigma_ij.to_vec(),
            sigma_ij_3,
            sigma_ij_5,
            a_ij_qq,
            b_ij_qq,
            c_ijk_qq,
            a_ij_dd,
            b_ij_dd,
            c_ijk_dd,
            a_ij_dq,
            b_ij_dq,
            c_ijk_dq,
            mu_indices,
            q_indices,
            enable_qq: true,
            enable_dd: true,
            enable_dq: true,
        })
    }

    /// Whether any polar interaction term is active.
    #[inline]
    pub fn is_active(&self) -> bool {
        (self.enable_dd && !self.mu_indices.is_empty())
            || (self.enable_qq && !self.q_indices.is_empty())
            || (self.enable_dq && !self.mu_indices.is_empty() && !self.q_indices.is_empty())
    }

    #[inline]
    fn hd_j2_5(t: HyperDual, eta: HyperDual, a: &[f64; 5], b: &[f64; 5]) -> HyperDual {
        let mut sum = HyperDual::real(0.0);
        let mut eta_l = HyperDual::real(1.0);
        for l in 0..5 {
            let coeff = HyperDual::real(a[l]) + HyperDual::real(b[l]) / t;
            sum = sum + coeff * eta_l;
            if l < 4 {
                eta_l = eta_l * eta;
            }
        }
        sum
    }

    #[inline]
    fn hd_j2_4(t: HyperDual, eta: HyperDual, a: &[f64; 4], b: &[f64; 4]) -> HyperDual {
        let mut sum = HyperDual::real(0.0);
        let mut eta_l = HyperDual::real(1.0);
        for l in 0..4 {
            let coeff = HyperDual::real(a[l]) + HyperDual::real(b[l]) / t;
            sum = sum + coeff * eta_l;
            if l < 3 {
                eta_l = eta_l * eta;
            }
        }
        sum
    }

    #[inline]
    fn hd_j3_4(eta: HyperDual, c: &[f64; 4]) -> HyperDual {
        let mut sum = HyperDual::real(0.0);
        let mut eta_l = HyperDual::real(1.0);
        for (l, &cl) in c.iter().enumerate() {
            sum = sum + HyperDual::real(cl) * eta_l;
            if l < 3 {
                eta_l = eta_l * eta;
            }
        }
        sum
    }

    #[inline]
    fn hd_j3_3(eta: HyperDual, c: &[f64; 3]) -> HyperDual {
        let mut sum = HyperDual::real(0.0);
        let mut eta_l = HyperDual::real(1.0);
        for (l, &cl) in c.iter().enumerate() {
            sum = sum + HyperDual::real(cl) * eta_l;
            if l < 2 {
                eta_l = eta_l * eta;
            }
        }
        sum
    }

    fn hyperdual_f_qq(
        &self,
        t: HyperDual,
        v: HyperDual,
        n: &[HyperDual],
        eta: HyperDual,
    ) -> HyperDual {
        let nc = self.nc;
        let t2 = t * t;
        let t3 = t2 * t;
        let v2 = v * v;
        let pre2 = HyperDual::real(-PI * 0.75 * 0.75 * N_AVOGADRO) / v / t2;
        let pre3 = HyperDual::real((PI * 0.75) * (PI * 0.75) * N_AVOGADRO * N_AVOGADRO) / v2 / t3;

        let mut a2 = HyperDual::real(0.0);
        let mut a3 = HyperDual::real(0.0);

        let mut add_a2 = |is: usize, js: usize, fac: f64| {
            let ij = is * nc + js;
            let ii = is * nc + is;
            let jj = js * nc + js;
            let j2 = Self::hd_j2_5(t, eta, &self.a_ij_qq[ij], &self.b_ij_qq[ij]);
            let s_ij = self.sigma_ij[ij];
            let coeff = fac
                * self.eps_divk_ij[ii]
                * self.eps_divk_ij[jj]
                * self.sigma_ij_5[ii]
                * self.sigma_ij_5[jj]
                / (self.sigma_ij_5[ij] * s_ij * s_ij)
                * self.q_star_2[is]
                * self.q_star_2[js];
            a2 = a2 + n[is] * n[js] * HyperDual::real(coeff) * j2;
        };

        let mut add_a3 = |is: usize, js: usize, ks: usize, fac: f64| {
            let ijk = (is * nc + js) * nc + ks;
            let ii = is * nc + is;
            let jj = js * nc + js;
            let kk = ks * nc + ks;
            let j3 = Self::hd_j3_4(eta, &self.c_ijk_qq[ijk]);
            let coeff = fac
                * self.eps_divk_ij[ii]
                * self.eps_divk_ij[jj]
                * self.eps_divk_ij[kk]
                * self.sigma_ij_5[ii]
                * self.sigma_ij_5[jj]
                * self.sigma_ij_5[kk]
                / (self.sigma_ij_3[is * nc + js]
                    * self.sigma_ij_3[is * nc + ks]
                    * self.sigma_ij_3[js * nc + ks])
                * self.q_star_2[is]
                * self.q_star_2[js]
                * self.q_star_2[ks];
            a3 = a3 + n[is] * n[js] * n[ks] * HyperDual::real(coeff) * j3;
        };

        let nq = self.q_indices.len();
        for ii in 0..nq {
            let i = self.q_indices[ii];
            add_a2(i, i, 1.0);
            add_a3(i, i, i, 1.0);
            for jj in (ii + 1)..nq {
                let j = self.q_indices[jj];
                add_a2(i, j, 2.0);
                add_a3(i, i, j, 3.0);
                add_a3(i, j, j, 3.0);
                for &k in &self.q_indices[(jj + 1)..nq] {
                    add_a3(i, j, k, 6.0);
                }
            }
        }

        a2 = pre2 * a2;
        a3 = pre3 * a3;
        if a2.f0 - a3.f0 != 0.0 {
            (a2 * a2) / (a2 - a3)
        } else {
            a2
        }
    }

    fn hyperdual_f_dd(
        &self,
        t: HyperDual,
        v: HyperDual,
        n: &[HyperDual],
        eta: HyperDual,
    ) -> HyperDual {
        let nc = self.nc;
        let t2 = t * t;
        let t3 = t2 * t;
        let v2 = v * v;
        let pre2 = HyperDual::real(-PI * N_AVOGADRO) / v / t2;
        let pre3 = HyperDual::real(-(4.0 / 3.0) * (PI * N_AVOGADRO) * (PI * N_AVOGADRO)) / v2 / t3;

        let mut a2 = HyperDual::real(0.0);
        let mut a3 = HyperDual::real(0.0);

        let mut add_a2 = |is: usize, js: usize, fac: f64| {
            let ij = is * nc + js;
            let ii = is * nc + is;
            let jj = js * nc + js;
            let j2 = Self::hd_j2_5(t, eta, &self.a_ij_dd[ij], &self.b_ij_dd[ij]);
            let coeff = fac
                * self.eps_divk_ij[ii]
                * self.eps_divk_ij[jj]
                * self.sigma_ij_3[ii]
                * self.sigma_ij_3[jj]
                / self.sigma_ij_3[ij]
                * self.mu_star_2[is]
                * self.mu_star_2[js];
            a2 = a2 + n[is] * n[js] * HyperDual::real(coeff) * j2;
        };

        let mut add_a3 = |is: usize, js: usize, ks: usize, fac: f64| {
            let ijk = (is * nc + js) * nc + ks;
            let ii = is * nc + is;
            let jj = js * nc + js;
            let kk = ks * nc + ks;
            let j3 = Self::hd_j3_4(eta, &self.c_ijk_dd[ijk]);
            let coeff = fac
                * self.eps_divk_ij[ii]
                * self.eps_divk_ij[jj]
                * self.eps_divk_ij[kk]
                * self.sigma_ij_3[ii]
                * self.sigma_ij_3[jj]
                * self.sigma_ij_3[kk]
                / (self.sigma_ij[is * nc + js]
                    * self.sigma_ij[is * nc + ks]
                    * self.sigma_ij[js * nc + ks])
                * self.mu_star_2[is]
                * self.mu_star_2[js]
                * self.mu_star_2[ks];
            a3 = a3 + n[is] * n[js] * n[ks] * HyperDual::real(coeff) * j3;
        };

        let nmu = self.mu_indices.len();
        for ii in 0..nmu {
            let i = self.mu_indices[ii];
            add_a2(i, i, 1.0);
            add_a3(i, i, i, 1.0);
            for jj in (ii + 1)..nmu {
                let j = self.mu_indices[jj];
                add_a2(i, j, 2.0);
                add_a3(i, i, j, 3.0);
                add_a3(i, j, j, 3.0);
                for &k in &self.mu_indices[(jj + 1)..nmu] {
                    add_a3(i, j, k, 6.0);
                }
            }
        }

        a2 = pre2 * a2;
        a3 = pre3 * a3;
        if a2.f0 - a3.f0 != 0.0 {
            (a2 * a2) / (a2 - a3)
        } else {
            a2
        }
    }

    fn hyperdual_f_dq(
        &self,
        t: HyperDual,
        v: HyperDual,
        n: &[HyperDual],
        eta: HyperDual,
    ) -> HyperDual {
        let nc = self.nc;
        let t2 = t * t;
        let t3 = t2 * t;
        let v2 = v * v;
        let pre2 = HyperDual::real(-(9.0 / 4.0) * PI * N_AVOGADRO) / v / t2;
        let pre3 = HyperDual::real(-(PI * N_AVOGADRO) * (PI * N_AVOGADRO)) / v2 / t3;

        let mut a2 = HyperDual::real(0.0);
        let mut a3 = HyperDual::real(0.0);

        for &i in &self.mu_indices {
            let ii = i * nc + i;
            for &j in &self.q_indices {
                let ij = i * nc + j;
                let jj = j * nc + j;
                let j2 = Self::hd_j2_4(t, eta, &self.a_ij_dq[ij], &self.b_ij_dq[ij]);
                let coeff2 = self.eps_divk_ij[ii]
                    * self.eps_divk_ij[jj]
                    * self.sigma_ij_3[ii]
                    * self.sigma_ij_5[jj]
                    / self.sigma_ij_5[ij]
                    * self.mu_star_2[i]
                    * self.q_star_2[j];
                a2 = a2 + n[i] * n[j] * HyperDual::real(coeff2) * j2;

                // add_to_a3_mu(i, k, j)
                for &k in &self.mu_indices {
                    let ijk = (i * nc + k) * nc + j;
                    let kk = k * nc + k;
                    let j3 = Self::hd_j3_3(eta, &self.c_ijk_dq[ijk]);
                    let s_num = self.sigma_ij[ii] * self.sigma_ij[kk] * self.sigma_ij[jj];
                    let s_num2 = s_num * s_num;
                    let s_den = self.sigma_ij[i * nc + k]
                        * self.sigma_ij[i * nc + j]
                        * self.sigma_ij[k * nc + j];
                    let coeff3 = self.eps_divk_ij[ii]
                        * self.eps_divk_ij[kk]
                        * self.eps_divk_ij[jj]
                        * (s_num2 * s_num2)
                        / (s_den * s_den)
                        * self.mu_star_2[k]
                        * self.mu_star_2[i]
                        * self.q_star_2[j];
                    a3 = a3 + n[i] * n[k] * n[j] * HyperDual::real(coeff3) * j3;
                }

                // add_to_a3_Q(i, j, k)
                for &k in &self.q_indices {
                    let ijk = (i * nc + j) * nc + k;
                    let kk = k * nc + k;
                    let j3 = Self::hd_j3_3(eta, &self.c_ijk_dq[ijk]);
                    let s_num = self.sigma_ij[ii] * self.sigma_ij[jj] * self.sigma_ij[kk];
                    let s_num2 = s_num * s_num;
                    let s_den = self.sigma_ij[i * nc + j]
                        * self.sigma_ij[i * nc + k]
                        * self.sigma_ij[j * nc + k];
                    let coeff3 = ALPHA_DQ
                        * self.eps_divk_ij[ii]
                        * self.eps_divk_ij[jj]
                        * self.eps_divk_ij[kk]
                        * (s_num2 * s_num2)
                        / (s_den * s_den)
                        * self.q_star_2[j]
                        * self.mu_star_2[i]
                        * self.q_star_2[k];
                    a3 = a3 + n[i] * n[j] * n[k] * HyperDual::real(coeff3) * j3;
                }
            }
        }

        a2 = pre2 * a2;
        a3 = pre3 * a3;
        if a2.f0 - a3.f0 != 0.0 {
            (a2 * a2) / (a2 - a3)
        } else {
            a2
        }
    }

    fn eval_hyperdual(&self, t: HyperDual, v: HyperDual, n: &[HyperDual]) -> HyperDual {
        let nc = self.nc;
        let mut eta_sum = HyperDual::real(0.0);
        for i in 0..nc {
            let ii = i * nc + i;
            let expo = (HyperDual::real(-3.0 * self.eps_divk_ij[ii]) / t).exp();
            let d = HyperDual::real(self.sigma_ij[ii])
                * (HyperDual::real(1.0) - HyperDual::real(0.12) * expo);
            let d3 = d * d * d;
            eta_sum = eta_sum + HyperDual::real(self.m[i]) * n[i] * d3;
        }
        let eta = (HyperDual::real(N_AVOGADRO * PI / 6.0) / v) * eta_sum;

        let mut f = HyperDual::real(0.0);
        if self.enable_dd && !self.mu_indices.is_empty() {
            f = f + self.hyperdual_f_dd(t, v, n, eta);
        }
        if self.enable_qq && !self.q_indices.is_empty() {
            f = f + self.hyperdual_f_qq(t, v, n, eta);
        }
        if self.enable_dq && !self.mu_indices.is_empty() && !self.q_indices.is_empty() {
            f = f + self.hyperdual_f_dq(t, v, n, eta);
        }
        f
    }

    /// Value-only multipolar reduced residual Helmholtz energy `F = A_mp / (R T)`.
    ///
    /// Evaluates in unit-mole `(v_u = V / N, x = n / N)` space and scales by `N`.
    pub fn fres_value(&self, t: f64, v: f64, n: &[f64]) -> Result<f64> {
        let sumn = validate_state(t, v, n, self.nc)?;
        if !self.is_active() {
            return Ok(0.0);
        }
        let v_u = v / sumn;
        let t_hd = HyperDual::with_order(t, 0);
        let v_hd = HyperDual::with_order(v_u, 0);
        let n_hd: Vec<HyperDual> = n
            .iter()
            .map(|&ni| HyperDual::with_order(ni / sumn, 0))
            .collect();
        let f_u = self.eval_hyperdual(t_hd, v_hd, &n_hd).f0;
        let f = sumn * f_u;
        if !f.is_finite() {
            return Err(Error::InvalidInput(
                "non-finite multipolar Helmholtz energy",
            ));
        }
        Ok(f)
    }

    /// Multipolar `(F, F_T, F_V, F_TT, F_TV, F_VV)` without computing `n^2` mole Hessians.
    pub fn fres_pressure_jet(&self, t: f64, v: f64, n: &[f64]) -> Result<MultipolPressureJet> {
        let sumn = validate_state(t, v, n, self.nc)?;
        if !self.is_active() {
            return Ok(MultipolPressureJet::default());
        }
        let v_u = v / sumn;
        let n_base: Vec<HyperDual> = n.iter().map(|&ni| HyperDual::real(ni / sumn)).collect();

        // Pass 1: (T, T) -> F, F_T, F_TT
        let r_tt = self.eval_hyperdual(
            HyperDual::second(t, 1.0, 1.0),
            HyperDual::real(v_u),
            &n_base,
        );
        // Pass 2: (V, V) -> F_V, F_VV
        let r_vv = self.eval_hyperdual(
            HyperDual::real(t),
            HyperDual::second(v_u, 1.0, 1.0),
            &n_base,
        );
        // Pass 3: (T, V) -> F_TV
        let r_tv = self.eval_hyperdual(
            HyperDual::second(t, 1.0, 0.0),
            HyperDual::second(v_u, 0.0, 1.0),
            &n_base,
        );

        let out = MultipolPressureJet {
            f: sumn * r_tt.f0,
            f_t: sumn * r_tt.f1,
            f_v: r_vv.f1,
            f_tt: sumn * r_tt.f12,
            f_tv: r_tv.f12,
            f_vv: r_vv.f12 / sumn,
        };
        if ![out.f, out.f_t, out.f_v, out.f_tt, out.f_tv, out.f_vv]
            .iter()
            .all(|val| val.is_finite())
        {
            return Err(Error::InvalidInput(
                "non-finite multipolar pressure-jet derivatives",
            ));
        }
        Ok(out)
    }

    /// Multipolar `F` and first mole derivatives `F_n` (for fast fugacity evaluation
    /// without computing the full `n^2` mole Hessian).
    pub fn fres_fn(&self, t: f64, v: f64, n: &[f64]) -> Result<Vec<f64>> {
        let sumn = validate_state(t, v, n, self.nc)?;
        if !self.is_active() {
            return Ok(vec![0.0; self.nc]);
        }
        let v_u = v / sumn;
        let t_hd = HyperDual::with_order(t, 1);
        let v_hd = HyperDual::with_order(v_u, 1);
        let mut n_hd: Vec<HyperDual> = n
            .iter()
            .map(|&ni| HyperDual::with_order(ni / sumn, 1))
            .collect();
        let mut f_n = vec![0.0; self.nc];
        for i in 0..self.nc {
            n_hd[i].f1 = 1.0;
            let res = self.eval_hyperdual(t_hd, v_hd, &n_hd);
            n_hd[i].f1 = 0.0;
            if !res.f1.is_finite() {
                return Err(Error::InvalidInput("non-finite multipolar mole derivative"));
            }
            f_n[i] = res.f1;
        }
        Ok(f_n)
    }

    /// Complete 10-block multipolar reduced residual Helmholtz energy (`F..F_nn`).
    ///
    /// Mirrors `multipol::add_hyperdual_fres_multipol` via `hyperdual_fres_wrapper`,
    /// evaluated in unit-mole `(v_u = V / N, x = n / N)` coordinates and Euler-scaled
    /// to `(V, n)` so `n_i n_j n_k / V^2` never underflows at tiny `N`.
    pub fn fres_all(&self, t: f64, v: f64, n: &[f64]) -> Result<MultipolHelmholtz> {
        let sumn = validate_state(t, v, n, self.nc)?;
        let nc = self.nc;
        if !self.is_active() {
            return Ok(MultipolHelmholtz::zeros(nc));
        }

        let v_u = v / sumn;
        let x: Vec<f64> = n.iter().map(|&ni| ni / sumn).collect();
        let mut n_hd: Vec<HyperDual> = x.iter().map(|&xi| HyperDual::real(xi)).collect();

        let r_tt = self.eval_hyperdual(HyperDual::second(t, 1.0, 1.0), HyperDual::real(v_u), &n_hd);
        let r_vv = self.eval_hyperdual(HyperDual::real(t), HyperDual::second(v_u, 1.0, 1.0), &n_hd);
        let r_tv = self.eval_hyperdual(
            HyperDual::second(t, 1.0, 0.0),
            HyperDual::second(v_u, 0.0, 1.0),
            &n_hd,
        );

        let mut f_n = vec![0.0; nc];
        let mut f_tn = vec![0.0; nc];
        let mut f_vn_u = vec![0.0; nc];
        let mut f_nn_u = vec![0.0; nc * nc];

        let t_e1 = HyperDual::second(t, 1.0, 0.0);
        let v_e1 = HyperDual::second(v_u, 1.0, 0.0);
        let t_0 = HyperDual::real(t);
        let v_0 = HyperDual::real(v_u);

        for i in 0..nc {
            n_hd[i].f2 = 1.0;
            let r_tn = self.eval_hyperdual(t_e1, v_0, &n_hd);
            let r_vn = self.eval_hyperdual(t_0, v_e1, &n_hd);
            f_n[i] = r_tn.f2;
            f_tn[i] = r_tn.f12;
            f_vn_u[i] = r_vn.f12;
            n_hd[i].f2 = 0.0;

            for j in i..nc {
                n_hd[i].f1 = 1.0;
                n_hd[j].f2 = 1.0;
                let r_nn = self.eval_hyperdual(t_0, v_0, &n_hd);
                n_hd[i].f1 = 0.0;
                n_hd[j].f2 = 0.0;
                f_nn_u[i * nc + j] = r_nn.f12;
                f_nn_u[j * nc + i] = r_nn.f12;
            }
        }

        let out = MultipolHelmholtz {
            f: sumn * r_tt.f0,
            f_t: sumn * r_tt.f1,
            f_v: r_vv.f1,
            f_n,
            f_tt: sumn * r_tt.f12,
            f_tv: r_tv.f12,
            f_tn,
            f_vv: r_vv.f12 / sumn,
            f_vn: f_vn_u.into_iter().map(|val| val / sumn).collect(),
            f_nn: f_nn_u.into_iter().map(|val| val / sumn).collect(),
        };

        if ![out.f, out.f_t, out.f_v, out.f_tt, out.f_tv, out.f_vv]
            .iter()
            .all(|val| val.is_finite())
            || !out.f_n.iter().all(|val| val.is_finite())
            || !out.f_tn.iter().all(|val| val.is_finite())
            || !out.f_vn.iter().all(|val| val.is_finite())
            || !out.f_nn.iter().all(|val| val.is_finite())
        {
            return Err(Error::InvalidInput(
                "non-finite multipolar Helmholtz derivatives",
            ));
        }

        Ok(out)
    }
}

fn validate_state(t: f64, v: f64, n: &[f64], expected_nc: usize) -> Result<f64> {
    positive(t, "temperature")?;
    positive(v, "volume")?;
    if n.len() != expected_nc {
        return Err(Error::InvalidInput("composition length mismatch"));
    }
    let mut sumn = 0.0;
    for &ni in n {
        if !ni.is_finite() || ni < 0.0 {
            return Err(Error::InvalidInput("non-negative mole numbers required"));
        }
        sumn += ni;
    }
    if !sumn.is_finite() || sumn <= 0.0 {
        return Err(Error::InvalidInput("positive total moles required"));
    }
    Ok(sumn)
}
