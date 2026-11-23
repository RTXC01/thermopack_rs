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

/// Coefficients of one family mole number in A2 and A3. Splitting A3 by
/// degree in this family keeps its quadratic contribution visible at trace.
#[derive(Clone, Copy)]
struct DqFamilyTerm {
    alpha: HyperDual,
    beta_linear: HyperDual,
    beta_quadratic: HyperDual,
}

impl DqFamilyTerm {
    fn beta(self) -> HyperDual {
        self.beta_linear + self.beta_quadratic
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

        if m.iter().any(|&mi| !mi.is_finite() || mi <= 0.0) {
            return Err(Error::InvalidInput(
                "segment numbers must be finite and strictly positive",
            ));
        }
        if mu.iter().any(|&v| !v.is_finite()) || q.iter().any(|&v| !v.is_finite()) {
            return Err(Error::InvalidInput("polar moments must be finite"));
        }
        if eps_divk_ij
            .iter()
            .any(|&e| !e.is_finite() || e <= 0.0 || !(e * e).is_finite())
        {
            return Err(Error::InvalidInput(
                "dispersion energy matrix entries must be finite, strictly positive, and have finite square",
            ));
        }

        let mut sigma_ij_3 = Vec::with_capacity(nc * nc);
        let mut sigma_ij_5 = Vec::with_capacity(nc * nc);
        for &s in sigma_ij {
            if !s.is_finite() || s <= 0.0 {
                return Err(Error::InvalidInput(
                    "segment diameter matrix entries must be finite and strictly positive",
                ));
            }
            let s2 = s * s;
            let s3 = s2 * s;
            let s5 = s3 * s2;
            if !s3.is_finite() || s3 <= 0.0 || !s5.is_finite() || s5 <= 0.0 {
                return Err(Error::InvalidInput(
                    "segment diameter powers overflowed or underflowed",
                ));
            }
            sigma_ij_3.push(s3);
            sigma_ij_5.push(s5);
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
                    if !b_ij_qq[ij][l].is_finite() || !b_ij_dd[ij][l].is_finite() {
                        return Err(Error::InvalidInput(
                            "non-finite multipolar correlation coefficients",
                        ));
                    }
                    if l < 4 {
                        a_ij_dq[ij][l] = ADQ[0][l] + ADQ[1][l] * h + ADQ[2][l] * g;
                        b_ij_dq[ij][l] = (BDQ[0][l] + BDQ[1][l] * h + BDQ[2][l] * g) * eps_ij;
                        if !b_ij_dq[ij][l].is_finite() {
                            return Err(Error::InvalidInput(
                                "non-finite multipolar correlation coefficients",
                            ));
                        }
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

        let out = Self {
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
        };
        out.validate_active_coefficients()?;
        Ok(out)
    }

    #[inline]
    fn dd_coeff2(&self, i: usize, j: usize, fac: f64) -> f64 {
        let nc = self.nc;
        let ii = i * nc + i;
        let jj = j * nc + j;
        let ij = i * nc + j;
        fac * (-PI)
            * ((N_AVOGADRO * self.sigma_ij_3[ii]) * (self.eps_divk_ij[ii] * self.mu_star_2[i]))
            * ((self.sigma_ij_3[jj] / self.sigma_ij_3[ij])
                * (self.eps_divk_ij[jj] * self.mu_star_2[j]))
    }

    #[inline]
    fn dd_coeff3(&self, i: usize, j: usize, k: usize, fac: f64) -> f64 {
        let nc = self.nc;
        let ii = i * nc + i;
        let jj = j * nc + j;
        let kk = k * nc + k;
        let s_k = self.sigma_ij[kk];
        let s_ratio = (s_k / self.sigma_ij[i * nc + j])
            * (s_k / self.sigma_ij[i * nc + k])
            * (s_k / self.sigma_ij[j * nc + k]);
        fac * (-(4.0 / 3.0) * PI * PI)
            * ((N_AVOGADRO * self.sigma_ij_3[ii]) * (self.eps_divk_ij[ii] * self.mu_star_2[i]))
            * ((N_AVOGADRO * self.sigma_ij_3[jj]) * (self.eps_divk_ij[jj] * self.mu_star_2[j]))
            * (s_ratio * (self.eps_divk_ij[kk] * self.mu_star_2[k]))
    }

    #[inline]
    fn qq_coeff2(&self, i: usize, j: usize, fac: f64) -> f64 {
        let nc = self.nc;
        let ii = i * nc + i;
        let jj = j * nc + j;
        let ij = i * nc + j;
        let r_i = self.sigma_ij[ii] / self.sigma_ij[ij];
        fac * (-PI * 0.75 * 0.75)
            * ((N_AVOGADRO * self.sigma_ij_3[ii])
                * (r_i * r_i)
                * (self.eps_divk_ij[ii] * self.q_star_2[i]))
            * ((self.sigma_ij_5[jj] / self.sigma_ij_5[ij])
                * (self.eps_divk_ij[jj] * self.q_star_2[j]))
    }

    #[inline]
    fn qq_coeff3(&self, i: usize, j: usize, k: usize, fac: f64) -> f64 {
        let nc = self.nc;
        let ii = i * nc + i;
        let jj = j * nc + j;
        let kk = k * nc + k;
        let s_i = self.sigma_ij[ii];
        let s_j = self.sigma_ij[jj];
        let s_k = self.sigma_ij[kk];
        let s_ij = self.sigma_ij[i * nc + j];
        let s_ik = self.sigma_ij[i * nc + k];
        let s_jk = self.sigma_ij[j * nc + k];
        let r_k_ik = s_k / s_ik;
        let r_k_jk = s_k / s_jk;
        fac * ((PI * 0.75) * (PI * 0.75))
            * ((N_AVOGADRO * self.sigma_ij_3[ii])
                * ((s_i / s_ij) * (s_i / s_ik))
                * (self.eps_divk_ij[ii] * self.q_star_2[i]))
            * ((N_AVOGADRO * self.sigma_ij_3[jj])
                * ((s_j / s_ij) * (s_j / s_jk))
                * (self.eps_divk_ij[jj] * self.q_star_2[j]))
            * (((s_k / s_ij) * (r_k_ik * r_k_ik) * (r_k_jk * r_k_jk))
                * (self.eps_divk_ij[kk] * self.q_star_2[k]))
    }

    #[inline]
    fn dq_coeff2(&self, i: usize, j: usize) -> f64 {
        let nc = self.nc;
        let ii = i * nc + i;
        let jj = j * nc + j;
        let ij = i * nc + j;
        (-(9.0 / 4.0) * PI)
            * ((N_AVOGADRO * self.sigma_ij_3[ii]) * (self.eps_divk_ij[ii] * self.mu_star_2[i]))
            * ((self.sigma_ij_5[jj] / self.sigma_ij_5[ij])
                * (self.eps_divk_ij[jj] * self.q_star_2[j]))
    }

    #[inline]
    fn dq_coeff3_mu(&self, i: usize, k: usize, j: usize) -> f64 {
        let nc = self.nc;
        let ii = i * nc + i;
        let kk = k * nc + k;
        let jj = j * nc + j;
        let s_i = self.sigma_ij[ii];
        let s_k = self.sigma_ij[kk];
        let s_j = self.sigma_ij[jj];
        let s_ik = self.sigma_ij[i * nc + k];
        let s_ij = self.sigma_ij[i * nc + j];
        let s_kj = self.sigma_ij[k * nc + j];
        let r_j_ij = s_j / s_ij;
        let r_j_kj = s_j / s_kj;
        (-PI * PI)
            * ((N_AVOGADRO * self.sigma_ij_3[ii])
                * (s_i / s_ik)
                * (self.eps_divk_ij[ii] * self.mu_star_2[i]))
            * ((N_AVOGADRO * self.sigma_ij_3[kk])
                * (s_k / s_ik)
                * (self.eps_divk_ij[kk] * self.mu_star_2[k]))
            * (((r_j_ij * r_j_ij) * (r_j_kj * r_j_kj)) * (self.eps_divk_ij[jj] * self.q_star_2[j]))
    }

    #[inline]
    fn dq_coeff3_q(&self, i: usize, j: usize, k: usize) -> f64 {
        let nc = self.nc;
        let ii = i * nc + i;
        let jj = j * nc + j;
        let kk = k * nc + k;
        let s_i = self.sigma_ij[ii];
        let s_j = self.sigma_ij[jj];
        let s_k = self.sigma_ij[kk];
        let s_ij = self.sigma_ij[i * nc + j];
        let s_ik = self.sigma_ij[i * nc + k];
        let s_jk = self.sigma_ij[j * nc + k];
        let r_k_ik = s_k / s_ik;
        let r_k_jk = s_k / s_jk;
        (-PI * PI * ALPHA_DQ)
            * ((N_AVOGADRO * self.sigma_ij_3[ii])
                * (s_i / s_ij)
                * (self.eps_divk_ij[ii] * self.mu_star_2[i]))
            * ((N_AVOGADRO * self.sigma_ij_3[jj])
                * (s_j / s_ij)
                * (self.eps_divk_ij[jj] * self.q_star_2[j]))
            * (((r_k_ik * r_k_ik) * (r_k_jk * r_k_jk)) * (self.eps_divk_ij[kk] * self.q_star_2[k]))
    }

    fn validate_active_coefficients(&self) -> Result<()> {
        let nmu = self.mu_indices.len();
        for ii in 0..nmu {
            let i = self.mu_indices[ii];
            if !self.dd_coeff2(i, i, 1.0).is_finite() || !self.dd_coeff3(i, i, i, 1.0).is_finite() {
                return Err(Error::InvalidInput("non-finite DD multipolar coefficient"));
            }
            for jj in (ii + 1)..nmu {
                let j = self.mu_indices[jj];
                if !self.dd_coeff2(i, j, 2.0).is_finite()
                    || !self.dd_coeff3(i, i, j, 3.0).is_finite()
                    || !self.dd_coeff3(i, j, j, 3.0).is_finite()
                {
                    return Err(Error::InvalidInput("non-finite DD multipolar coefficient"));
                }
                for &k in &self.mu_indices[(jj + 1)..nmu] {
                    if !self.dd_coeff3(i, j, k, 6.0).is_finite() {
                        return Err(Error::InvalidInput("non-finite DD multipolar coefficient"));
                    }
                }
            }
        }
        let nq = self.q_indices.len();
        for ii in 0..nq {
            let i = self.q_indices[ii];
            if !self.qq_coeff2(i, i, 1.0).is_finite() || !self.qq_coeff3(i, i, i, 1.0).is_finite() {
                return Err(Error::InvalidInput("non-finite QQ multipolar coefficient"));
            }
            for jj in (ii + 1)..nq {
                let j = self.q_indices[jj];
                if !self.qq_coeff2(i, j, 2.0).is_finite()
                    || !self.qq_coeff3(i, i, j, 3.0).is_finite()
                    || !self.qq_coeff3(i, j, j, 3.0).is_finite()
                {
                    return Err(Error::InvalidInput("non-finite QQ multipolar coefficient"));
                }
                for &k in &self.q_indices[(jj + 1)..nq] {
                    if !self.qq_coeff3(i, j, k, 6.0).is_finite() {
                        return Err(Error::InvalidInput("non-finite QQ multipolar coefficient"));
                    }
                }
            }
        }
        for &i in &self.mu_indices {
            for &j in &self.q_indices {
                if !self.dq_coeff2(i, j).is_finite() {
                    return Err(Error::InvalidInput("non-finite DQ multipolar coefficient"));
                }
                for &k in &self.mu_indices {
                    if !self.dq_coeff3_mu(i, k, j).is_finite() {
                        return Err(Error::InvalidInput("non-finite DQ multipolar coefficient"));
                    }
                }
                for &k in &self.q_indices {
                    if !self.dq_coeff3_q(i, j, k).is_finite() {
                        return Err(Error::InvalidInput("non-finite DQ multipolar coefficient"));
                    }
                }
            }
        }
        Ok(())
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

    /// Evaluate u^2/(u-v). The Hessian needs (u dv - v du)/(u-v)
    /// for each seed. Form these normalized products only at second order.
    fn pade_exact_with_wronskian(
        u: HyperDual,
        v: HyperDual,
        order: u8,
        normalized_wronskians: impl FnOnce(f64) -> (f64, f64),
    ) -> Result<HyperDual> {
        let w0 = u.f0 - v.f0;
        if !w0.is_finite() || w0 == 0.0 {
            return Err(Error::InvalidInput("singular multipolar Pade denominator"));
        }
        let r = u.f0 / w0;
        let f0 = u.f0 * r;
        if !r.is_finite() || !f0.is_finite() {
            return Err(Error::InvalidInput("non-finite multipolar Pade result"));
        }
        let (f1, f2) = if order >= 1 {
            let f1 = r * (2.0 * u.f1 - r * (u.f1 - v.f1));
            let f2 = r * (2.0 * u.f2 - r * (u.f2 - v.f2));
            if !f1.is_finite() || !f2.is_finite() {
                return Err(Error::InvalidInput("non-finite multipolar Pade result"));
            }
            (f1, f2)
        } else {
            (0.0, 0.0)
        };
        let f12 = if order >= 2 {
            let (t1, t2) = normalized_wronskians(w0);
            let f12 = r * (2.0 * u.f12 - r * (u.f12 - v.f12)) + 2.0 * t1 * (t2 / w0);
            if !f12.is_finite() {
                return Err(Error::InvalidInput("non-finite multipolar Pade result"));
            }
            f12
        } else {
            0.0
        };
        Ok(HyperDual {
            f0,
            f1,
            f2,
            f3: 0.0,
            f12,
            f13: 0.0,
            f23: 0.0,
            f123: 0.0,
            order,
        })
    }

    fn pade_exact(u: HyperDual, v: HyperDual, order: u8) -> Result<HyperDual> {
        Self::pade_exact_with_wronskian(u, v, order, |w| {
            let u_over_w = u.f0 / w;
            let v_over_w = v.f0 / w;
            (
                u_over_w * v.f1 - v_over_w * u.f1,
                u_over_w * v.f2 - v_over_w * u.f2,
            )
        })
    }

    fn pade_dd_qq(a2: HyperDual, a3: HyperDual, order: u8) -> Result<HyperDual> {
        // For DD and QQ when the polar family is absent or ultra-trace,
        // `a_2 = O(eps^2)` and `a_3 = O(eps^3)`, so `a_3 / a_2` is below f64
        // precision and `a_2` is exact through second derivatives without any
        // `w_0` underflow.
        if a2.f0.abs() < 1.0e-100 && a3.f0.abs() <= 1.0e-20 * a2.f0.abs() {
            return Ok(a2);
        }
        Self::pade_exact(a2, a3, order)
    }

    /// Test proportionality of A2 and the part of A3 linear in this family.
    /// The varying species' moment cancels in the coefficient ratios.
    /// For fixed i, coeff3_mu(i,i,j)/coeff2(i,j) varies as s_ij/s_jj.
    /// For fixed j, coeff3_q(i,j,j)/coeff2(i,j) varies as s_ii*s_ij.
    /// Matching each active cross diameter this way also gives the common
    /// scaling for triples containing two distinct opposite species.
    fn dq_linear_proportional(
        &self,
        a: usize,
        b: usize,
        quadrupoles: bool,
        opposite_indices: &[usize],
        n: &[HyperDual],
    ) -> bool {
        if a == b {
            return true;
        }
        if self.m[a].min(2.0) != self.m[b].min(2.0) {
            return false;
        }
        let nc = self.nc;
        let sa = self.sigma_ij[a * nc + a];
        let sb = self.sigma_ij[b * nc + b];
        opposite_indices
            .iter()
            .filter(|&&i| n[i].f0 > 0.0)
            .all(|&i| {
                let (ia, ib) = if quadrupoles {
                    (i * nc + a, i * nc + b)
                } else {
                    (a * nc + i, b * nc + i)
                };
                let (x, y) = if quadrupoles { (sb, sa) } else { (sa, sb) };
                let ca = self.sigma_ij[ia];
                let cb = self.sigma_ij[ib];
                let product = ca * x;
                // The constructor bounds diameter fifth powers, so both products
                // and their rounding residuals are in range. Compare both parts
                // to distinguish even ratios that round to the same float.
                self.eps_divk_ij[ia] == self.eps_divk_ij[ib]
                    && product == cb * y
                    && ca.mul_add(x, -product) == cb.mul_add(y, -product)
            })
    }

    /// Exact L'Hôpital DQ Padé boundary evaluator when `a_2(0) == 0 && a_3(0) == 0`.
    fn pade_dq_boundary(
        a2: HyperDual,
        a3: HyperDual,
        same_direction: bool,
        proportional: bool,
        order: u8,
    ) -> Result<HyperDual> {
        let f1_nz = a2.f1 != 0.0;
        let f2_nz = a2.f2 != 0.0;
        if !f1_nz && !f2_nz {
            return Ok(a2);
        }
        if f1_nz && !f2_nz {
            let d = a2.f1 - a3.f1;
            if !d.is_finite() || d == 0.0 {
                return Err(Error::InvalidInput(
                    "singular multipolar DQ Pade boundary denominator",
                ));
            }
            let r = a2.f1 / d;
            let f1 = a2.f1 * r;
            let f12 = if order >= 2 {
                2.0 * r * a2.f12 - (r * r) * (a2.f12 - a3.f12)
            } else {
                0.0
            };
            if !f1.is_finite() || !f12.is_finite() {
                return Err(Error::InvalidInput(
                    "non-finite multipolar DQ Pade boundary derivative",
                ));
            }
            return Ok(HyperDual {
                f0: 0.0,
                f1,
                f2: 0.0,
                f3: 0.0,
                f12,
                f13: 0.0,
                f23: 0.0,
                f123: 0.0,
                order,
            });
        }
        if !f1_nz && f2_nz {
            let d = a2.f2 - a3.f2;
            if !d.is_finite() || d == 0.0 {
                return Err(Error::InvalidInput(
                    "singular multipolar DQ Pade boundary denominator",
                ));
            }
            let r = a2.f2 / d;
            let f2 = a2.f2 * r;
            let f12 = if order >= 2 {
                2.0 * r * a2.f12 - (r * r) * (a2.f12 - a3.f12)
            } else {
                0.0
            };
            if !f2.is_finite() || !f12.is_finite() {
                return Err(Error::InvalidInput(
                    "non-finite multipolar DQ Pade boundary derivative",
                ));
            }
            return Ok(HyperDual {
                f0: 0.0,
                f1: 0.0,
                f2,
                f3: 0.0,
                f12,
                f13: 0.0,
                f23: 0.0,
                f123: 0.0,
                order,
            });
        }

        // Both `a2.f1 != 0.0` and `a2.f2 != 0.0`.
        let d1 = a2.f1 - a3.f1;
        let d2 = a2.f2 - a3.f2;
        if !d1.is_finite() || d1 == 0.0 || !d2.is_finite() || d2 == 0.0 {
            return Err(Error::InvalidInput(
                "singular multipolar DQ Pade boundary denominator",
            ));
        }
        // Only allow the mixed boundary derivative when the directions are
        // identical, structurally proportional, or `a_3` vanishes to higher
        // order (`a3.f1 == 0 && a3.f2 == 0`, i.e., opposite family also absent).
        if same_direction || proportional || (a3.f1 == 0.0 && a3.f2 == 0.0) {
            // Structural proportionality gives a common directional ratio.
            // Either nonzero direction therefore supplies the Hessian limit.
            let r = a2.f1 / d1;
            let f1 = a2.f1 * r;
            let f2 = a2.f2 * (a2.f2 / d2);
            let f12 = if order >= 2 {
                r * (2.0 * a2.f12 - r * (a2.f12 - a3.f12))
            } else {
                0.0
            };
            if !f1.is_finite() || !f2.is_finite() || !f12.is_finite() {
                return Err(Error::InvalidInput(
                    "non-finite multipolar DQ Pade boundary derivative",
                ));
            }
            Ok(HyperDual {
                f0: 0.0,
                f1,
                f2,
                f3: 0.0,
                f12,
                f13: 0.0,
                f23: 0.0,
                f123: 0.0,
                order,
            })
        } else {
            Err(Error::InvalidInput(
                "undefined multipolar DQ Pade mixed boundary derivative with multiple absent polar species",
            ))
        }
    }

    /// Evaluates `u^2 / (u - v)` for linear family sums `u = sum_a n[f_a] * alpha[a]`,
    /// `v = sum_a n[f_a] * beta[a]` using a cancellation-free Wronskian that cancels
    /// diagonal and structurally proportional species terms before subtraction.
    #[allow(clippy::too_many_arguments)]
    fn pade_dq_family_sum(
        &self,
        quadrupoles: bool,
        family_indices: &[usize],
        opposite_indices: &[usize],
        n: &[HyperDual],
        terms: &[DqFamilyTerm],
        same_direction: bool,
        order: u8,
    ) -> Result<HyperDual> {
        let mut u = HyperDual::with_order(0.0, order);
        let mut v = HyperDual::with_order(0.0, order);
        for (&i, term) in family_indices.iter().zip(terms) {
            u = u + n[i] * term.alpha;
            v = v + n[i] * term.beta();
        }
        if u.f0 == 0.0 && v.f0 == 0.0 {
            let mut proportional = true;
            'check: for &fa in family_indices {
                if n[fa].f1 != 0.0 {
                    for &fb in family_indices {
                        if n[fb].f2 != 0.0
                            && !self.dq_linear_proportional(
                                fa,
                                fb,
                                quadrupoles,
                                opposite_indices,
                                n,
                            )
                        {
                            proportional = false;
                            break 'check;
                        }
                    }
                }
            }
            return Self::pade_dq_boundary(u, v, same_direction, proportional, order);
        }
        Self::pade_exact_with_wronskian(u, v, order, |w| {
            let mut t1 = 0.0;
            let mut t2 = 0.0;
            let u_over_w = u.f0 / w;
            let v_over_w = v.f0 / w;
            for (a, (&fa, ta)) in family_indices.iter().zip(terms).enumerate() {
                let na = n[fa];
                let ba = ta.beta();
                for (&fb, tb) in family_indices[a + 1..].iter().zip(&terms[a + 1..]) {
                    let nb = n[fb];
                    // Cancel only the proven linear determinant. The quadratic
                    // terms can have different same-family cross diameters.
                    let (ba0, bb0) =
                        if self.dq_linear_proportional(fa, fb, quadrupoles, opposite_indices, n) {
                            (ta.beta_quadratic.f0, tb.beta_quadratic.f0)
                        } else {
                            (ba.f0, tb.beta().f0)
                        };
                    let cross = |dna, dnb| {
                        let weight = (na.f0 * dnb - nb.f0 * dna) / w;
                        (weight * ta.alpha.f0) * bb0 - (weight * tb.alpha.f0) * ba0
                    };
                    t1 += cross(na.f1, nb.f1);
                    t2 += cross(na.f2, nb.f2);
                }
                // Coefficient derivatives carry another family mole factor.
                t1 += u_over_w * (na.f0 * ba.f1) - v_over_w * (na.f0 * ta.alpha.f1);
                t2 += u_over_w * (na.f0 * ba.f2) - v_over_w * (na.f0 * ta.alpha.f2);
            }
            (t1, t2)
        })
    }

    fn hyperdual_f_qq(
        &self,
        t: HyperDual,
        v: HyperDual,
        n: &[HyperDual],
        eta: HyperDual,
    ) -> Result<HyperDual> {
        let nc = self.nc;
        let order = t.order.max(v.order);
        let t2 = t * t;
        let t3 = t2 * t;
        let v2 = v * v;
        let inv_vt2 = HyperDual::real(1.0) / v / t2;
        let inv_v2t3 = HyperDual::real(1.0) / v2 / t3;

        let mut a2 = HyperDual::with_order(0.0, order);
        let mut a3 = HyperDual::with_order(0.0, order);

        let mut add_a2 = |is: usize, js: usize, fac: f64| {
            let ij = is * nc + js;
            let j2 = Self::hd_j2_5(t, eta, &self.a_ij_qq[ij], &self.b_ij_qq[ij]);
            let coeff = self.qq_coeff2(is, js, fac);
            a2 = a2 + (HyperDual::real(coeff) * inv_vt2 * j2) * (n[is] * n[js]);
        };

        let mut add_a3 = |is: usize, js: usize, ks: usize, fac: f64| {
            let ijk = (is * nc + js) * nc + ks;
            let j3 = Self::hd_j3_4(eta, &self.c_ijk_qq[ijk]);
            let coeff = self.qq_coeff3(is, js, ks, fac);
            a3 = a3 + (HyperDual::real(coeff) * inv_v2t3 * j3) * (n[is] * n[js] * n[ks]);
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

        Self::pade_dd_qq(a2, a3, order)
    }

    fn hyperdual_f_dd(
        &self,
        t: HyperDual,
        v: HyperDual,
        n: &[HyperDual],
        eta: HyperDual,
    ) -> Result<HyperDual> {
        let nc = self.nc;
        let order = t.order.max(v.order);
        let t2 = t * t;
        let t3 = t2 * t;
        let v2 = v * v;
        let inv_vt2 = HyperDual::real(1.0) / v / t2;
        let inv_v2t3 = HyperDual::real(1.0) / v2 / t3;

        let mut a2 = HyperDual::with_order(0.0, order);
        let mut a3 = HyperDual::with_order(0.0, order);

        let mut add_a2 = |is: usize, js: usize, fac: f64| {
            let ij = is * nc + js;
            let j2 = Self::hd_j2_5(t, eta, &self.a_ij_dd[ij], &self.b_ij_dd[ij]);
            let coeff = self.dd_coeff2(is, js, fac);
            a2 = a2 + (HyperDual::real(coeff) * inv_vt2 * j2) * (n[is] * n[js]);
        };

        let mut add_a3 = |is: usize, js: usize, ks: usize, fac: f64| {
            let ijk = (is * nc + js) * nc + ks;
            let j3 = Self::hd_j3_4(eta, &self.c_ijk_dd[ijk]);
            let coeff = self.dd_coeff3(is, js, ks, fac);
            a3 = a3 + (HyperDual::real(coeff) * inv_v2t3 * j3) * (n[is] * n[js] * n[ks]);
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

        Self::pade_dd_qq(a2, a3, order)
    }

    fn hyperdual_f_dq(
        &self,
        t: HyperDual,
        v: HyperDual,
        n: &[HyperDual],
        eta: HyperDual,
        same_direction: bool,
    ) -> Result<HyperDual> {
        let nc = self.nc;
        let order = t.order.max(v.order);
        let t2 = t * t;
        let t3 = t2 * t;
        let v2 = v * v;
        let inv_vt2 = HyperDual::real(1.0) / v / t2;
        let inv_v2t3 = HyperDual::real(1.0) / v2 / t3;

        let dq_term2 = |i: usize, j: usize| -> HyperDual {
            let ij = i * nc + j;
            let j2 = Self::hd_j2_4(t, eta, &self.a_ij_dq[ij], &self.b_ij_dq[ij]);
            HyperDual::real(self.dq_coeff2(i, j)) * inv_vt2 * j2
        };

        let dq_term3_mu = |i: usize, k: usize, j: usize| -> HyperDual {
            let ijk = (i * nc + k) * nc + j;
            let j3 = Self::hd_j3_3(eta, &self.c_ijk_dq[ijk]);
            HyperDual::real(self.dq_coeff3_mu(i, k, j)) * inv_v2t3 * j3
        };

        let dq_term3_q = |i: usize, j: usize, k: usize| -> HyperDual {
            let ijk = (i * nc + j) * nc + k;
            let j3 = Self::hd_j3_3(eta, &self.c_ijk_dq[ijk]);
            HyperDual::real(self.dq_coeff3_q(i, j, k)) * inv_v2t3 * j3
        };

        // Single dipole and single quadrupole — structurally factor out
        // BOTH `n[m]` and `n[q]` so the Padé quotient is O(1) smooth everywhere.
        if self.mu_indices.len() == 1 && self.q_indices.len() == 1 {
            let m = self.mu_indices[0];
            let q = self.q_indices[0];
            let a2_bar = dq_term2(m, q);
            let a3_bar = n[m] * dq_term3_mu(m, m, q) + n[q] * dq_term3_q(m, q, q);
            let p_bar = Self::pade_exact(a2_bar, a3_bar, order)?;
            return Ok(n[m] * n[q] * p_bar);
        }

        // Expand in the trace family. A singleton opposite family has a
        // common mole factor in A2 and A3 which can be removed exactly.
        let s_mu: f64 = self.mu_indices.iter().map(|&i| n[i].f0).sum();
        let s_q: f64 = self.q_indices.iter().map(|&j| n[j].f0).sum();
        let quadrupoles = self.mu_indices.len() == 1 || (self.q_indices.len() > 1 && s_q <= s_mu);
        let (family, opposite) = if quadrupoles {
            (&self.q_indices, &self.mu_indices)
        } else {
            (&self.mu_indices, &self.q_indices)
        };
        let mut terms = Vec::with_capacity(family.len());
        for &a in family {
            let mut alpha = HyperDual::with_order(0.0, order);
            let mut beta_linear = HyperDual::with_order(0.0, order);
            let mut beta_quadratic = HyperDual::with_order(0.0, order);
            for &b in opposite {
                let weight = if opposite.len() == 1 {
                    HyperDual::real(1.0)
                } else {
                    n[b]
                };
                let (i, j) = if quadrupoles { (b, a) } else { (a, b) };
                alpha = alpha + weight * dq_term2(i, j);
                for &k in opposite {
                    let term = if quadrupoles {
                        dq_term3_mu(i, k, j)
                    } else {
                        dq_term3_q(i, j, k)
                    };
                    beta_linear = beta_linear + weight * n[k] * term;
                }
                for &k in family {
                    let term = if quadrupoles {
                        dq_term3_q(i, j, k)
                    } else {
                        dq_term3_mu(i, k, j)
                    };
                    beta_quadratic = beta_quadratic + weight * n[k] * term;
                }
            }
            terms.push(DqFamilyTerm {
                alpha,
                beta_linear,
                beta_quadratic,
            });
        }
        let pade = self.pade_dq_family_sum(
            quadrupoles,
            family,
            opposite,
            n,
            &terms,
            same_direction,
            order,
        )?;
        if opposite.len() == 1 {
            Ok(n[opposite[0]] * pade)
        } else {
            Ok(pade)
        }
    }

    fn eval_hyperdual(
        &self,
        t: HyperDual,
        v: HyperDual,
        n: &[HyperDual],
        same_direction: bool,
    ) -> Result<HyperDual> {
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
            f = f + self.hyperdual_f_dd(t, v, n, eta)?;
        }
        if self.enable_qq && !self.q_indices.is_empty() {
            f = f + self.hyperdual_f_qq(t, v, n, eta)?;
        }
        if self.enable_dq && !self.mu_indices.is_empty() && !self.q_indices.is_empty() {
            f = f + self.hyperdual_f_dq(t, v, n, eta, same_direction)?;
        }
        Ok(f)
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
        let f_u = self.eval_hyperdual(t_hd, v_hd, &n_hd, false)?.f0;
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
            true,
        )?;
        // Pass 2: (V, V) -> F_V, F_VV
        let r_vv = self.eval_hyperdual(
            HyperDual::real(t),
            HyperDual::second(v_u, 1.0, 1.0),
            &n_base,
            true,
        )?;
        // Pass 3: (T, V) -> F_TV
        let r_tv = self.eval_hyperdual(
            HyperDual::second(t, 1.0, 0.0),
            HyperDual::second(v_u, 0.0, 1.0),
            &n_base,
            false,
        )?;

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
            let res = self.eval_hyperdual(t_hd, v_hd, &n_hd, false)?;
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

        let r_tt = self.eval_hyperdual(
            HyperDual::second(t, 1.0, 1.0),
            HyperDual::real(v_u),
            &n_hd,
            true,
        )?;
        let r_vv = self.eval_hyperdual(
            HyperDual::real(t),
            HyperDual::second(v_u, 1.0, 1.0),
            &n_hd,
            true,
        )?;
        let r_tv = self.eval_hyperdual(
            HyperDual::second(t, 1.0, 0.0),
            HyperDual::second(v_u, 0.0, 1.0),
            &n_hd,
            false,
        )?;

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
            let r_tn = self.eval_hyperdual(t_e1, v_0, &n_hd, false)?;
            let r_vn = self.eval_hyperdual(t_0, v_e1, &n_hd, false)?;
            f_n[i] = r_tn.f2;
            f_tn[i] = r_tn.f12;
            f_vn_u[i] = r_vn.f12;
            n_hd[i].f2 = 0.0;

            for j in i..nc {
                n_hd[i].f1 = 1.0;
                n_hd[j].f2 = 1.0;
                let r_nn = self.eval_hyperdual(t_0, v_0, &n_hd, i == j)?;
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
