//! The pair correlation function of the additive hard-sphere model of
//! Boublik.
//!
//! The original `src/hardsphere_bmcsl.f90:calc_bmcsl_dgij_dzeta` evaluates the
//! reduced Helmholtz energy of a pair of components and its first and second
//! derivatives with respect to the packing fractions of the mixture. The
//! expression is the one of the paper of Boublik, "Hard-Sphere Equation of
//! State", J. Chem. Phys. 53, 471 (1970), as the original states.
//!
//! The module holds the three pieces of the original that the radial
//! distribution function of the PC-SAFT model needs:
//!
//! - the hard-sphere diameter of every component with its temperature
//!   derivatives, [`HsDiameters`],
//! - the packing fractions of the mixture with their derivatives,
//!   [`packing_fractions`],
//! - the reduced Helmholtz energy of the pair and the chain rule to the
//!   temperature, the volume and the mole numbers, [`pair_correlation`].

/// The pair correlation function of the additive hard-sphere model and its
/// derivatives with respect to the packing fractions.
///
/// The derivatives are with respect to the second and the third packing
/// fraction and with respect to the reduced diameter of the pair, in the order
/// of the original `dg_dzeta` and `d2g_dzeta2`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BmcslPairCorrelation {
    /// The pair correlation function.
    pub value: f64,
    /// The first derivatives with respect to the packing fractions.
    pub dzeta: [f64; 3],
    /// The second derivatives with respect to the packing fractions.
    pub d2zeta: [[f64; 3]; 3],
}

/// The pair correlation function of the additive hard-sphere model.
///
/// `zeta_2` is the second packing fraction of the mixture, `zeta_3` the third
/// one and `mu` the reduced diameter of the pair, the ratio of the product of
/// the hard-sphere diameters and their sum.
pub fn boublik_pair_correlation(zeta_2: f64, zeta_3: f64, mu: f64) -> BmcslPairCorrelation {
    let div_diff_1 = 1.0 / (1.0 - zeta_3);
    let div_diff_2 = div_diff_1 * div_diff_1;
    let div_diff_3 = div_diff_1 * div_diff_2;
    let div_diff_4 = div_diff_2 * div_diff_2;
    let mu_2 = mu * mu;
    let zeta_2_2 = zeta_2 * zeta_2;

    let value = div_diff_1 + 3.0 * zeta_2 * div_diff_2 * mu + 2.0 * zeta_2_2 * mu_2 * div_diff_3;

    // The first derivatives with respect to (zeta_2, zeta_3, mu).
    let dzeta = [
        3.0 * mu * div_diff_2 + 4.0 * div_diff_3 * mu_2 * zeta_2,
        div_diff_2 + 6.0 * zeta_2 * div_diff_3 * mu + 6.0 * zeta_2_2 * div_diff_4 * mu_2,
        3.0 * zeta_2 * div_diff_2 + 4.0 * zeta_2_2 * div_diff_3 * mu,
    ];

    // The second derivatives with respect to (zeta_2, zeta_3, mu).
    let mut d2zeta = [[0.0; 3]; 3];
    d2zeta[0][0] = 4.0 * div_diff_3 * mu_2;
    d2zeta[1][1] = 2.0 * div_diff_3
        + 18.0 * zeta_2 * mu * div_diff_4
        + 24.0 * zeta_2_2 * mu_2 * div_diff_3 * div_diff_2;
    d2zeta[2][2] = 4.0 * zeta_2_2 * div_diff_3;
    d2zeta[2][0] = 3.0 * div_diff_2 + 8.0 * zeta_2 * mu * div_diff_3;
    d2zeta[0][2] = d2zeta[2][0];
    d2zeta[0][1] = 6.0 * div_diff_3 * mu + 12.0 * zeta_2 * div_diff_4 * mu_2;
    d2zeta[1][0] = d2zeta[0][1];
    d2zeta[2][1] = 6.0 * zeta_2 * div_diff_3 + 12.0 * zeta_2_2 * div_diff_4 * mu;
    d2zeta[1][2] = d2zeta[2][1];

    BmcslPairCorrelation {
        value,
        dzeta,
        d2zeta,
    }
}

/// The hard-sphere diameter of every component of the mixture and its
/// temperature derivatives.
///
/// The original holds the diameter in the type `hs_diameter` of
/// `src/hardsphere_bmcsl.f90`.
#[derive(Debug, Clone)]
pub struct HsDiameters {
    /// The hard-sphere diameter of every component \[m\].
    pub d: Vec<f64>,
    /// The temperature derivative of the diameter \[m/K\].
    pub d_t: Vec<f64>,
    /// The second temperature derivative of the diameter \[m/K^2\].
    pub d_tt: Vec<f64>,
}

impl HsDiameters {
    /// The hard-sphere diameter of the pair of the components `i` and `j` and
    /// its temperature derivatives.
    ///
    /// The reduced diameter of the pair is the ratio of the product of the
    /// diameters and their sum, the quantity that the original stores in the
    /// fifth entry of the packing fractions.
    pub fn reduced_diameter(&self, i: usize, j: usize) -> ReducedDiameter {
        let d_i = self.d[i];
        let d_j = self.d[j];
        let sum = d_i + d_j;
        let mu = d_i * d_j / sum;
        let mu_t = (d_i * d_i * self.d_t[j] + d_j * d_j * self.d_t[i]) / (sum * sum);
        let first = 2.0 * sum * self.d_t[i] * self.d_t[j]
            + d_i * d_i * self.d_tt[j]
            + d_j * d_j * self.d_tt[i];
        let second =
            (d_i * d_i * self.d_t[j] + d_j * d_j * self.d_t[i]) * (self.d_t[i] + self.d_t[j]);
        ReducedDiameter {
            value: mu,
            d_t: mu_t,
            d_tt: first / (sum * sum) - 2.0 * second / (sum * sum * sum),
        }
    }
}

/// The reduced diameter of a pair of components and its temperature
/// derivatives.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReducedDiameter {
    /// The reduced diameter of the pair.
    pub value: f64,
    /// The temperature derivative of the reduced diameter.
    pub d_t: f64,
    /// The second temperature derivative of the reduced diameter.
    pub d_tt: f64,
}

/// The moments of the number density of the hard-sphere mixture and their
/// derivatives with respect to the temperature and the volume.
///
/// The moments are the ones of the original
/// `src/hardsphere_bmcsl.f90:calc_bmcsl_zeta_and_derivatives`, in the order
/// zeta_0 to zeta_3. The container holds no derivative with respect to the
/// mole numbers, so the evaluation needs no allocation.
#[derive(Debug, Clone)]
pub struct PackingMoments {
    /// The number of components of the mixture.
    pub nc: usize,
    /// The moments zeta_0 to zeta_3 of the number density.
    pub zet: [f64; 4],
    /// The temperature derivatives of the moments.
    pub zet_t: [f64; 4],
    /// The second temperature derivatives of the moments.
    pub zet_tt: [f64; 4],
    /// The volume derivatives of the moments.
    pub zet_v: [f64; 4],
    /// The second volume derivatives of the moments.
    pub zet_vv: [f64; 4],
    /// The temperature and volume derivatives of the moments.
    pub zet_tv: [f64; 4],
}

/// The packing fractions of the hard-sphere mixture and their derivatives.
#[derive(Debug, Clone)]
pub struct PackingFractions {
    /// The number of components of the mixture.
    pub nc: usize,
    /// The moments of the number density with their temperature and volume
    /// derivatives.
    pub moments: PackingMoments,
    /// The mole number derivatives of the moments.
    pub zet_n: Vec<[f64; 4]>,
    /// The volume and mole number derivatives of the moments.
    pub zet_vn: Vec<[f64; 4]>,
    /// The temperature and mole number derivatives of the moments.
    pub zet_tn: Vec<[f64; 4]>,
}

/// The packing fractions of the hard-sphere mixture.
///
/// `volume` is the volume of the mixture \[m^3\], `moles` the mole numbers
/// \[mol\], `diameters` the hard-sphere diameters with their temperature
/// derivatives and `segments` the number of segments of every component.
pub fn packing_moments(
    volume: f64,
    moles: &[f64],
    diameters: &HsDiameters,
    segments: &[f64],
) -> PackingMoments {
    let nc = moles.len();
    let prefactor = std::f64::consts::PI * crate::pc_saft::N_AVOGADRO / (6.0 * volume);

    let mut zet = [0.0; 4];
    let mut zet_t = [0.0; 4];
    let mut zet_tt = [0.0; 4];
    let mut zet_tv = [0.0; 4];

    for component in 0..nc {
        let m = segments[component];
        let d = diameters.d[component];
        let d_t = diameters.d_t[component];
        let d_tt = diameters.d_tt[component];
        let n = moles[component];
        for l in 0..4 {
            zet[l] += m * n * d.powi(l as i32);
            if l > 0 {
                let exponent = l as i32 - 1;
                let d_previous = d.powi(exponent);
                zet_t[l] += l as f64 * m * n * d_previous * d_t;
                zet_tt[l] += l as f64 * m * n * d_previous * d_tt;
                zet_tv[l] -= l as f64 * m * n * d_previous * d_t / volume;
                if l > 1 {
                    zet_tt[l] +=
                        l as f64 * (l as f64 - 1.0) * m * n * d.powi(exponent - 1) * d_t * d_t;
                }
            }
        }
    }

    let mut zet_v = [0.0; 4];
    let mut zet_vv = [0.0; 4];
    for l in 0..4 {
        zet[l] *= prefactor;
        zet_t[l] *= prefactor;
        zet_tt[l] *= prefactor;
        zet_tv[l] *= prefactor;
        zet_v[l] = -zet[l] / volume;
        zet_vv[l] = 2.0 * zet[l] / (volume * volume);
    }

    PackingMoments {
        nc,
        zet,
        zet_t,
        zet_tt,
        zet_v,
        zet_vv,
        zet_tv,
    }
}

/// The packing fractions of the hard-sphere mixture.
///
/// `volume` is the volume of the mixture \[m^3\], `moles` the mole numbers
/// \[mol\], `diameters` the hard-sphere diameters with their temperature
/// derivatives and `segments` the number of segments of every component.
pub fn packing_fractions(
    volume: f64,
    moles: &[f64],
    diameters: &HsDiameters,
    segments: &[f64],
) -> PackingFractions {
    let nc = moles.len();
    let moments = packing_moments(volume, moles, diameters, segments);
    let prefactor = std::f64::consts::PI * crate::pc_saft::N_AVOGADRO / (6.0 * volume);

    let mut zet_n = vec![[0.0; 4]; nc];
    let mut zet_vn = vec![[0.0; 4]; nc];
    let mut zet_tn = vec![[0.0; 4]; nc];
    for component in 0..nc {
        let m = segments[component];
        let d = diameters.d[component];
        let d_t = diameters.d_t[component];
        for l in 0..4 {
            let power = d.powi(l as i32);
            zet_n[component][l] = prefactor * m * power;
            zet_vn[component][l] = -prefactor * m * power / volume;
            if l > 0 {
                zet_tn[component][l] = prefactor * l as f64 * m * d.powi(l as i32 - 1) * d_t;
            }
        }
    }

    PackingFractions {
        nc,
        moments,
        zet_n,
        zet_vn,
        zet_tn,
    }
}

/// The pair correlation function of a pair of components and its derivatives
/// with respect to the temperature, the volume and the mole numbers.
///
/// The chain rule is the one of the original
/// `src/hardsphere_bmcsl.f90:calc_bmcsl_gij`.
#[derive(Debug, Clone)]
pub struct PairCorrelation {
    /// The pair correlation function.
    pub value: f64,
    /// The temperature derivative of the pair correlation function.
    pub d_t: f64,
    /// The volume derivative of the pair correlation function.
    pub d_v: f64,
    /// The mole number derivatives of the pair correlation function.
    pub d_n: Vec<f64>,
    /// The second temperature derivative of the pair correlation function.
    pub d_tt: f64,
    /// The temperature and volume derivative of the pair correlation function.
    pub d_tv: f64,
    /// The temperature and mole number derivatives of the pair correlation
    /// function.
    pub d_tn: Vec<f64>,
    /// The second volume derivative of the pair correlation function.
    pub d_vv: f64,
    /// The volume and mole number derivatives of the pair correlation function.
    pub d_vn: Vec<f64>,
    /// The second mole number derivatives of the pair correlation function.
    pub d_nn: Vec<Vec<f64>>,
}

/// The pair correlation function of the pair of the components `i` and `j`.
///
/// `packing` holds the packing fractions of the mixture and `diameters` the
/// hard-sphere diameters with their temperature derivatives.
pub fn pair_correlation(
    packing: &PackingFractions,
    diameters: &HsDiameters,
    i: usize,
    j: usize,
) -> PairCorrelation {
    let nc = packing.nc;
    let zeta_2 = packing.moments.zet[2];
    let zeta_3 = packing.moments.zet[3];
    let mu = diameters.reduced_diameter(i, j);

    let kernel = boublik_pair_correlation(zeta_2, zeta_3, mu.value);
    let dg = kernel.dzeta;
    let d2g = kernel.d2zeta;

    // The derivatives of the three variables of the kernel, in the order of
    // the packing fraction of the second order, the packing fraction of the
    // third order and the reduced diameter of the pair.
    let zet_t = [packing.moments.zet_t[2], packing.moments.zet_t[3], mu.d_t];
    let zet_tt = [
        packing.moments.zet_tt[2],
        packing.moments.zet_tt[3],
        mu.d_tt,
    ];
    let zet_v = [packing.moments.zet_v[2], packing.moments.zet_v[3], 0.0];
    let zet_vv = [packing.moments.zet_vv[2], packing.moments.zet_vv[3], 0.0];
    let zet_tv = [packing.moments.zet_tv[2], packing.moments.zet_tv[3], 0.0];

    let mut d_n = vec![0.0; nc];
    let mut d_tn = vec![0.0; nc];
    let mut d_vn = vec![0.0; nc];
    let mut d_nn = vec![vec![0.0; nc]; nc];
    for component in 0..nc {
        let zet_n = [
            packing.zet_n[component][2],
            packing.zet_n[component][3],
            0.0,
        ];
        let zet_vn = [
            packing.zet_vn[component][2],
            packing.zet_vn[component][3],
            0.0,
        ];
        let zet_tn = [
            packing.zet_tn[component][2],
            packing.zet_tn[component][3],
            0.0,
        ];

        d_n[component] = dot(&dg, &zet_n);
        d_vn[component] = dot(&dg, &zet_vn) + double_sum(&d2g, &zet_v, &zet_n);
        d_tn[component] = dot(&dg, &zet_tn) + double_sum(&d2g, &zet_t, &zet_n);
        for (other, entry) in d_nn[component].iter_mut().enumerate() {
            let zet_n_other = [packing.zet_n[other][2], packing.zet_n[other][3], 0.0];
            *entry = double_sum(&d2g, &zet_n, &zet_n_other);
        }
    }

    PairCorrelation {
        value: kernel.value,
        d_t: dot(&dg, &zet_t),
        d_v: dot(&dg, &zet_v),
        d_n,
        d_tt: dot(&dg, &zet_tt) + double_sum(&d2g, &zet_t, &zet_t),
        d_tv: dot(&dg, &zet_tv) + double_sum(&d2g, &zet_t, &zet_v),
        d_tn,
        d_vv: dot(&dg, &zet_vv) + double_sum(&d2g, &zet_v, &zet_v),
        d_vn,
        d_nn,
    }
}

/// The logarithm of the pair correlation function and its derivatives.
///
/// The expression follows the original
/// `src/hardsphere_bmcsl.f90:calc_bmcsl_lngij`.
#[derive(Debug, Clone)]
pub struct LogPairCorrelation {
    /// The logarithm of the pair correlation function.
    pub value: f64,
    /// The temperature derivative of the logarithm.
    pub d_t: f64,
    /// The volume derivative of the logarithm.
    pub d_v: f64,
    /// The mole number derivatives of the logarithm.
    pub d_n: Vec<f64>,
    /// The second temperature derivative of the logarithm.
    pub d_tt: f64,
    /// The temperature and volume derivative of the logarithm.
    pub d_tv: f64,
    /// The temperature and mole number derivatives of the logarithm.
    pub d_tn: Vec<f64>,
    /// The second volume derivative of the logarithm.
    pub d_vv: f64,
    /// The volume and mole number derivatives of the logarithm.
    pub d_vn: Vec<f64>,
    /// The second mole number derivatives of the logarithm.
    pub d_nn: Vec<Vec<f64>>,
}

impl PairCorrelation {
    /// The logarithm of the pair correlation function and its derivatives.
    pub fn logarithm(&self) -> LogPairCorrelation {
        let g = self.value;
        let g_t = self.d_t;
        let g_v = self.d_v;
        let count = self.d_n.len();
        let mut d_n = Vec::with_capacity(count);
        let mut d_tn = Vec::with_capacity(count);
        let mut d_vn = Vec::with_capacity(count);
        let mut d_nn = vec![vec![0.0; count]; count];
        for (component, row) in d_nn.iter_mut().enumerate() {
            let g_n = self.d_n[component];
            d_n.push(g_n / g);
            d_tn.push(self.d_tn[component] / g - g_t * g_n / (g * g));
            d_vn.push(self.d_vn[component] / g - g_v * g_n / (g * g));
            for (other, entry) in row.iter_mut().enumerate() {
                *entry = self.d_nn[component][other] / g - g_n * self.d_n[other] / (g * g);
            }
        }
        LogPairCorrelation {
            value: g.ln(),
            d_t: g_t / g,
            d_v: g_v / g,
            d_n,
            d_tt: self.d_tt / g - g_t * g_t / (g * g),
            d_tv: self.d_tv / g - g_t * g_v / (g * g),
            d_tn,
            d_vv: self.d_vv / g - g_v * g_v / (g * g),
            d_vn,
            d_nn,
        }
    }
}

/// The reduced Helmholtz energy of the additive hard-sphere model and its
/// derivatives with respect to the packing fractions.
///
/// The expression and its derivatives follow the original
/// `src/hardsphere_bmcsl.f90:calc_bmcsl_dalpha_dzeta`, which the original
/// keeps private.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HardSphereEnergy {
    /// The reduced Helmholtz energy of the mixture.
    pub value: f64,
    /// The first derivatives with respect to the packing fractions.
    pub dzeta: [f64; 4],
    /// The second derivatives with respect to the packing fractions.
    pub d2zeta: [[f64; 4]; 4],
}

/// The reduced Helmholtz energy of the additive hard-sphere model.
///
/// `zeta` holds the packing fractions of the mixture, from the zeroth to the
/// third order.
pub fn boublik_helmholtz_energy(zeta: [f64; 4]) -> HardSphereEnergy {
    let zeta_0 = zeta[0];
    let zeta_1 = zeta[1];
    let zeta_2 = zeta[2];
    let zeta_3 = zeta[3];

    let prefactor = 1.0 / zeta_0;
    let en_m_zeta_3 = 1.0 - zeta_3;
    let ln_term = en_m_zeta_3.abs().ln();

    let value = prefactor
        * (3.0 * zeta_1 * zeta_2 / en_m_zeta_3
            + zeta_2.powi(3) / (zeta_3 * en_m_zeta_3.powi(2))
            + (zeta_2.powi(3) / zeta_3.powi(2) - zeta_0) * ln_term);

    // The first order derivatives with respect to the packing fractions.
    let dzeta = [
        -prefactor * (value + ln_term),
        prefactor * 3.0 * zeta_2 / en_m_zeta_3,
        prefactor
            * (3.0 * zeta_1 / en_m_zeta_3
                + 3.0 * zeta_2.powi(2) / (zeta_3 * en_m_zeta_3.powi(2))
                + 3.0 * zeta_2.powi(2) / zeta_3.powi(2) * ln_term),
        prefactor
            * (3.0 * zeta_1 * zeta_2 / en_m_zeta_3.powi(2)
                + zeta_2.powi(3) * (3.0 * zeta_3 - 1.0) / (zeta_3.powi(2) * en_m_zeta_3.powi(3))
                - 2.0 * zeta_2.powi(3) / zeta_3.powi(3) * ln_term
                - (zeta_2.powi(3) - zeta_0 * zeta_3.powi(2)) / (zeta_3.powi(2) * en_m_zeta_3)),
    ];

    // The second order derivatives with respect to the packing fractions.
    let mut d2zeta = [[0.0; 4]; 4];
    d2zeta[0][0] = prefactor * prefactor * (2.0 * value + 2.0 * ln_term);
    d2zeta[1][1] = 0.0;
    d2zeta[2][2] =
        (prefactor * 6.0 * zeta_2 / zeta_3) * (1.0 / en_m_zeta_3.powi(2) + ln_term / zeta_3);
    d2zeta[3][3] = prefactor
        * (6.0 * zeta_1 * zeta_2 / en_m_zeta_3.powi(3)
            + 3.0 * zeta_2.powi(3) / (zeta_3.powi(2) * en_m_zeta_3.powi(3))
            + 3.0 * zeta_2.powi(3) * (3.0 * zeta_3 - 1.0) / (zeta_3.powi(2) * en_m_zeta_3.powi(4))
            - 2.0 * zeta_2.powi(3) * (3.0 * zeta_3 - 1.0) / (zeta_3.powi(3) * en_m_zeta_3.powi(3))
            + 2.0 * zeta_2.powi(3) / (zeta_3.powi(3) * en_m_zeta_3)
            + 6.0 * zeta_2.powi(3) * ln_term / zeta_3.powi(4)
            + 2.0 * zeta_0 * zeta_3 / (zeta_3.powi(2) * en_m_zeta_3)
            + 2.0 * (zeta_2.powi(3) - zeta_0 * zeta_3.powi(2)) / (zeta_3.powi(3) * en_m_zeta_3)
            - (zeta_2.powi(3) - zeta_0 * zeta_3.powi(2)) / (zeta_3.powi(2) * en_m_zeta_3.powi(2)));
    d2zeta[0][1] = -prefactor * dzeta[1];
    d2zeta[1][0] = d2zeta[0][1];
    d2zeta[0][2] = -prefactor * dzeta[2];
    d2zeta[2][0] = d2zeta[0][2];
    d2zeta[0][3] = -prefactor * (dzeta[3] - 1.0 / en_m_zeta_3);
    d2zeta[3][0] = d2zeta[0][3];
    d2zeta[1][2] = 3.0 * prefactor / en_m_zeta_3;
    d2zeta[2][1] = d2zeta[1][2];
    d2zeta[2][3] = prefactor
        * (3.0 * zeta_1 / en_m_zeta_3.powi(2)
            + 3.0 * zeta_2.powi(2) * (3.0 * zeta_3 - 1.0) / (zeta_3.powi(2) * en_m_zeta_3.powi(3))
            - 6.0 * zeta_2.powi(2) / zeta_3.powi(3) * ln_term
            - 3.0 * zeta_2.powi(2) / (zeta_3.powi(2) * en_m_zeta_3));
    d2zeta[3][2] = d2zeta[2][3];
    d2zeta[1][3] = dzeta[1] / en_m_zeta_3;
    d2zeta[3][1] = d2zeta[1][3];

    HardSphereEnergy {
        value,
        dzeta,
        d2zeta,
    }
}

/// The temperature derivatives of the reduced Helmholtz energy of the
/// hard-sphere mixture.
///
/// The function delivers the first and the second derivative with respect to
/// the temperature at a fixed volume and fixed mole numbers, and needs no
/// allocation.
pub fn hard_sphere_energy_t_derivatives(moments: &PackingMoments) -> (f64, f64) {
    let energy = boublik_helmholtz_energy(moments.zet);
    let d_t = dot4(&energy.dzeta, &moments.zet_t);
    let d_tt = dot4(&energy.dzeta, &moments.zet_tt)
        + double_sum4(&energy.d2zeta, &moments.zet_t, &moments.zet_t);
    (d_t, d_tt)
}

/// The pair correlation function of a pair of components and its temperature
/// derivatives.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairTemperatureDerivatives {
    /// The pair correlation function.
    pub value: f64,
    /// The temperature derivative of the pair correlation function.
    pub d_t: f64,
    /// The second temperature derivative of the pair correlation function.
    pub d_tt: f64,
    /// The logarithm of the pair correlation function.
    pub log_value: f64,
    /// The temperature derivative of the logarithm.
    pub log_d_t: f64,
    /// The second temperature derivative of the logarithm.
    pub log_d_tt: f64,
}

/// The pair correlation function of the pair of the components `i` and `j`
/// and its temperature derivatives.
///
/// The function delivers the first and the second derivative with respect to
/// the temperature at a fixed volume and fixed mole numbers, and needs no
/// allocation.
pub fn pair_correlation_t_derivatives(
    moments: &PackingMoments,
    diameters: &HsDiameters,
    i: usize,
    j: usize,
) -> PairTemperatureDerivatives {
    let mu = diameters.reduced_diameter(i, j);
    let kernel = boublik_pair_correlation(moments.zet[2], moments.zet[3], mu.value);
    let zet_t = [moments.zet_t[2], moments.zet_t[3], mu.d_t];
    let zet_tt = [moments.zet_tt[2], moments.zet_tt[3], mu.d_tt];

    let d_t = dot(&kernel.dzeta, &zet_t);
    let d_tt = dot(&kernel.dzeta, &zet_tt) + double_sum(&kernel.d2zeta, &zet_t, &zet_t);

    let g = kernel.value;
    PairTemperatureDerivatives {
        value: g,
        d_t,
        d_tt,
        log_value: g.ln(),
        log_d_t: d_t / g,
        log_d_tt: d_tt / g - d_t * d_t / (g * g),
    }
}

/// The reduced Helmholtz energy of the hard-sphere mixture and its
/// derivatives with respect to the temperature, the volume and the mole
/// numbers.
///
/// The chain rule is the one of the original
/// `src/hardsphere_bmcsl.f90:calc_bmcsl_helmholtzenergy`.
#[derive(Debug, Clone)]
pub struct HardSphereEnergyDerivatives {
    /// The reduced Helmholtz energy of the mixture.
    pub value: f64,
    /// The temperature derivative of the reduced Helmholtz energy.
    pub d_t: f64,
    /// The volume derivative of the reduced Helmholtz energy.
    pub d_v: f64,
    /// The mole number derivatives of the reduced Helmholtz energy.
    pub d_n: Vec<f64>,
    /// The second temperature derivative of the reduced Helmholtz energy.
    pub d_tt: f64,
    /// The temperature and volume derivative of the reduced Helmholtz energy.
    pub d_tv: f64,
    /// The temperature and mole number derivatives of the reduced Helmholtz
    /// energy.
    pub d_tn: Vec<f64>,
    /// The second volume derivative of the reduced Helmholtz energy.
    pub d_vv: f64,
    /// The volume and mole number derivatives of the reduced Helmholtz energy.
    pub d_vn: Vec<f64>,
    /// The second mole number derivatives of the reduced Helmholtz energy.
    pub d_nn: Vec<Vec<f64>>,
}

/// The reduced Helmholtz energy of the hard-sphere mixture.
///
/// `packing` holds the packing fractions of the mixture with their
/// derivatives.
pub fn hard_sphere_energy(packing: &PackingFractions) -> HardSphereEnergyDerivatives {
    let nc = packing.nc;
    let energy = boublik_helmholtz_energy(packing.moments.zet);
    let dalpha = energy.dzeta;
    let d2alpha = energy.d2zeta;

    let mut d_n = vec![0.0; nc];
    let mut d_tn = vec![0.0; nc];
    let mut d_vn = vec![0.0; nc];
    let mut d_nn = vec![vec![0.0; nc]; nc];
    for component in 0..nc {
        let zet_n = packing.zet_n[component];
        let zet_vn = packing.zet_vn[component];
        let zet_tn = packing.zet_tn[component];

        d_n[component] = dot4(&dalpha, &zet_n);
        d_vn[component] =
            dot4(&dalpha, &zet_vn) + double_sum4(&d2alpha, &packing.moments.zet_v, &zet_n);
        d_tn[component] =
            dot4(&dalpha, &zet_tn) + double_sum4(&d2alpha, &packing.moments.zet_t, &zet_n);
        for (other, entry) in d_nn[component].iter_mut().enumerate() {
            *entry = double_sum4(&d2alpha, &zet_n, &packing.zet_n[other]);
        }
    }

    HardSphereEnergyDerivatives {
        value: energy.value,
        d_t: dot4(&dalpha, &packing.moments.zet_t),
        d_v: dot4(&dalpha, &packing.moments.zet_v),
        d_n,
        d_tt: dot4(&dalpha, &packing.moments.zet_tt)
            + double_sum4(&d2alpha, &packing.moments.zet_t, &packing.moments.zet_t),
        d_tv: dot4(&dalpha, &packing.moments.zet_tv)
            + double_sum4(&d2alpha, &packing.moments.zet_t, &packing.moments.zet_v),
        d_tn,
        d_vv: dot4(&dalpha, &packing.moments.zet_vv)
            + double_sum4(&d2alpha, &packing.moments.zet_v, &packing.moments.zet_v),
        d_vn,
        d_nn,
    }
}

/// The scalar product of two vectors of the chain rule.
fn dot(left: &[f64; 3], right: &[f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

/// The scalar product of two vectors of the chain rule of the packing
/// fractions.
fn dot4(left: &[f64; 4], right: &[f64; 4]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2] + left[3] * right[3]
}

/// The double sum of the second derivatives of the chain rule.
fn double_sum(second: &[[f64; 3]; 3], left: &[f64; 3], right: &[f64; 3]) -> f64 {
    let mut sum = 0.0;
    for l in 0..3 {
        for k in 0..3 {
            sum += second[l][k] * left[l] * right[k];
        }
    }
    sum
}

/// The double sum of the second derivatives of the chain rule of the packing
/// fractions.
fn double_sum4(second: &[[f64; 4]; 4], left: &[f64; 4], right: &[f64; 4]) -> f64 {
    let mut sum = 0.0;
    for l in 0..4 {
        for k in 0..4 {
            sum += second[l][k] * left[l] * right[k];
        }
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The derivatives with respect to the packing fractions are checked
    /// against a central difference of the expression.
    #[test]
    fn the_derivatives_match_the_difference_quotient() {
        let states = [(0.05, 0.1, 0.3), (0.12, 0.25, 0.4), (0.3, 0.4, 0.25)];
        for (zeta_2, zeta_3, mu) in states {
            let value = boublik_pair_correlation(zeta_2, zeta_3, mu);
            let step = 1.0e-6;
            for index in 0..3 {
                let mut up = [zeta_2, zeta_3, mu];
                let mut down = [zeta_2, zeta_3, mu];
                up[index] += step;
                down[index] -= step;
                let quotient = (boublik_pair_correlation(up[0], up[1], up[2]).value
                    - boublik_pair_correlation(down[0], down[1], down[2]).value)
                    / (2.0 * step);
                let difference = (quotient - value.dzeta[index]).abs();
                assert!(
                    difference < 1.0e-6,
                    "the first derivative {index} differs by {difference}"
                );
                for second in 0..3 {
                    let mut up = [zeta_2, zeta_3, mu];
                    let mut down = [zeta_2, zeta_3, mu];
                    up[index] += step;
                    down[index] -= step;
                    let quotient = (boublik_pair_correlation(up[0], up[1], up[2]).dzeta[second]
                        - boublik_pair_correlation(down[0], down[1], down[2]).dzeta[second])
                        / (2.0 * step);
                    let difference = (quotient - value.d2zeta[index][second]).abs();
                    assert!(
                        difference < 1.0e-4,
                        "the second derivative {index},{second} differs by {difference}"
                    );
                }
            }
        }
    }

    /// The segment numbers of the mixture of the test.
    const SEGMENTS: [f64; 2] = [2.0, 1.0];

    /// The hard-sphere diameters of the mixture of the test.
    ///
    /// The diameters follow the expression of the PC-SAFT model of the
    /// original, so the temperature derivatives are the ones of
    /// `src/pc_saft_nonassoc.f90:calc_d`.
    fn diameters(temperature: f64) -> HsDiameters {
        let sigma = [3.5e-10, 2.6e-10];
        let eps = [90.0, 120.0];
        let mut d = Vec::new();
        let mut d_t = Vec::new();
        let mut d_tt = Vec::new();
        for component in 0..2 {
            let tinv = 1.0 / temperature;
            let expo = sigma[component] * (-3.0 * eps[component] * tinv).exp();
            d.push(sigma[component] - 0.12 * expo);
            d_t.push(-0.36 * eps[component] * tinv * tinv * expo);
            d_tt.push(
                expo * (0.72 * eps[component] - tinv * 1.08 * eps[component].powi(2))
                    * tinv.powi(3),
            );
        }
        HsDiameters { d, d_t, d_tt }
    }

    /// The pair correlation function of the pair of the test.
    fn correlation(temperature: f64, volume: f64, moles: &[f64], i: usize, j: usize) -> f64 {
        let diameters = diameters(temperature);
        let packing = packing_fractions(volume, moles, &diameters, &SEGMENTS);
        let mu = diameters.reduced_diameter(i, j);
        boublik_pair_correlation(packing.moments.zet[2], packing.moments.zet[3], mu.value).value
    }

    /// The reduced Helmholtz energy of the hard-sphere mixture of the test.
    fn energy(temperature: f64, volume: f64, moles: &[f64]) -> f64 {
        let diameters = diameters(temperature);
        let packing = packing_fractions(volume, moles, &diameters, &SEGMENTS);
        boublik_helmholtz_energy(packing.moments.zet).value
    }

    /// The chain rule of the reduced Helmholtz energy of the hard-sphere
    /// mixture is checked against a central difference of the expression.
    /// The chain rule of the reduced Helmholtz energy of the hard-sphere
    /// mixture is checked against a central difference of the expression.
    #[test]
    fn the_energy_chain_rule_matches_the_difference_quotient() {
        let temperature = 300.0;
        let volume = 1.0e-3;
        let moles = [0.05, 0.05];

        let diameters = diameters(temperature);
        let packing = packing_fractions(volume, &moles, &diameters, &SEGMENTS);
        let value = hard_sphere_energy(&packing);

        let step_t = temperature * 1.0e-5;
        let step_v = volume * 1.0e-5;
        let step_n = moles[0] * 1.0e-5;

        let quotient_t = (energy(temperature + step_t, volume, &moles)
            - energy(temperature - step_t, volume, &moles))
            / (2.0 * step_t);
        assert_close(quotient_t, value.d_t, 1.0e-6, "the temperature derivative");

        let quotient_v = (energy(temperature, volume + step_v, &moles)
            - energy(temperature, volume - step_v, &moles))
            / (2.0 * step_v);
        assert_close(quotient_v, value.d_v, 1.0e-6, "the volume derivative");

        for component in 0..2 {
            let mut up = moles;
            let mut down = moles;
            up[component] += step_n;
            down[component] -= step_n;
            let quotient = (energy(temperature, volume, &up) - energy(temperature, volume, &down))
                / (2.0 * step_n);
            assert_close(
                quotient,
                value.d_n[component],
                1.0e-6,
                "the mole number derivative",
            );
        }

        let step_t = temperature * 1.0e-3;
        let step_v = volume * 1.0e-3;
        let step_n = moles[0] * 1.0e-3;

        let quotient_tt = (energy(temperature + step_t, volume, &moles)
            - 2.0 * energy(temperature, volume, &moles)
            + energy(temperature - step_t, volume, &moles))
            / (step_t * step_t);
        assert_close(
            quotient_tt,
            value.d_tt,
            1.0e-4,
            "the second temperature derivative",
        );

        let quotient_vv = (energy(temperature, volume + step_v, &moles)
            - 2.0 * energy(temperature, volume, &moles)
            + energy(temperature, volume - step_v, &moles))
            / (step_v * step_v);
        assert_close(
            quotient_vv,
            value.d_vv,
            1.0e-4,
            "the second volume derivative",
        );

        for component in 0..2 {
            let mut up = moles;
            let mut down = moles;
            up[component] += step_n;
            down[component] -= step_n;
            let quotient_tn = (energy(temperature + step_t, volume, &up)
                - energy(temperature + step_t, volume, &down)
                - energy(temperature - step_t, volume, &up)
                + energy(temperature - step_t, volume, &down))
                / (4.0 * step_t * step_n);
            assert_close(
                quotient_tn,
                value.d_tn[component],
                1.0e-4,
                "the temperature and mole number derivative",
            );

            let quotient_vn = (energy(temperature, volume + step_v, &up)
                - energy(temperature, volume + step_v, &down)
                - energy(temperature, volume - step_v, &up)
                + energy(temperature, volume - step_v, &down))
                / (4.0 * step_v * step_n);
            assert_close(
                quotient_vn,
                value.d_vn[component],
                1.0e-4,
                "the volume and mole number derivative",
            );
        }
    }

    /// The chain rule of the pair correlation function is checked against a
    /// central difference of the expression.
    #[test]
    fn the_chain_rule_matches_the_difference_quotient() {
        let temperature = 300.0;
        let volume = 1.0e-3;
        let moles = [0.05, 0.05];
        let (i, j) = (0, 1);

        let diameters = diameters(temperature);
        let packing = packing_fractions(volume, &moles, &diameters, &SEGMENTS);
        let value = pair_correlation(&packing, &diameters, i, j);

        // The relative steps of the difference quotients. The second
        // derivatives of the volume and the mole numbers need a larger step,
        // because the packing fractions scale with the inverse of the volume.
        let step_t = temperature * 1.0e-5;
        let step_v = volume * 1.0e-5;
        let step_n = moles[0] * 1.0e-5;

        let quotient_t = (correlation(temperature + step_t, volume, &moles, i, j)
            - correlation(temperature - step_t, volume, &moles, i, j))
            / (2.0 * step_t);
        assert_close(quotient_t, value.d_t, 1.0e-6, "the temperature derivative");

        let quotient_v = (correlation(temperature, volume + step_v, &moles, i, j)
            - correlation(temperature, volume - step_v, &moles, i, j))
            / (2.0 * step_v);
        assert_close(quotient_v, value.d_v, 1.0e-6, "the volume derivative");

        for component in 0..2 {
            let mut up = moles;
            let mut down = moles;
            up[component] += step_n;
            down[component] -= step_n;
            let quotient = (correlation(temperature, volume, &up, i, j)
                - correlation(temperature, volume, &down, i, j))
                / (2.0 * step_n);
            assert_close(
                quotient,
                value.d_n[component],
                1.0e-6,
                "the mole number derivative",
            );
        }

        let step_t = temperature * 1.0e-3;
        let step_v = volume * 1.0e-3;
        let step_n = moles[0] * 1.0e-3;

        let quotient_tt = (correlation(temperature + step_t, volume, &moles, i, j)
            - 2.0 * correlation(temperature, volume, &moles, i, j)
            + correlation(temperature - step_t, volume, &moles, i, j))
            / (step_t * step_t);
        assert_close(
            quotient_tt,
            value.d_tt,
            1.0e-4,
            "the second temperature derivative",
        );

        let quotient_vv = (correlation(temperature, volume + step_v, &moles, i, j)
            - 2.0 * correlation(temperature, volume, &moles, i, j)
            + correlation(temperature, volume - step_v, &moles, i, j))
            / (step_v * step_v);
        assert_close(
            quotient_vv,
            value.d_vv,
            1.0e-4,
            "the second volume derivative",
        );

        let quotient_tv = (correlation(temperature + step_t, volume + step_v, &moles, i, j)
            - correlation(temperature + step_t, volume - step_v, &moles, i, j)
            - correlation(temperature - step_t, volume + step_v, &moles, i, j)
            + correlation(temperature - step_t, volume - step_v, &moles, i, j))
            / (4.0 * step_t * step_v);
        assert_close(
            quotient_tv,
            value.d_tv,
            1.0e-4,
            "the temperature and volume derivative",
        );

        for component in 0..2 {
            let mut up = moles;
            let mut down = moles;
            up[component] += step_n;
            down[component] -= step_n;
            let quotient_tn = (correlation(temperature + step_t, volume, &up, i, j)
                - correlation(temperature + step_t, volume, &down, i, j)
                - correlation(temperature - step_t, volume, &up, i, j)
                + correlation(temperature - step_t, volume, &down, i, j))
                / (4.0 * step_t * step_n);
            assert_close(
                quotient_tn,
                value.d_tn[component],
                1.0e-4,
                "the temperature and mole number derivative",
            );

            let quotient_vn = (correlation(temperature, volume + step_v, &up, i, j)
                - correlation(temperature, volume + step_v, &down, i, j)
                - correlation(temperature, volume - step_v, &up, i, j)
                + correlation(temperature, volume - step_v, &down, i, j))
                / (4.0 * step_v * step_n);
            assert_close(
                quotient_vn,
                value.d_vn[component],
                1.0e-4,
                "the volume and mole number derivative",
            );
        }

        let step_n = moles[0] * 1.0e-2;
        for component in 0..2 {
            let mut up = moles;
            let mut down = moles;
            up[component] += step_n;
            down[component] -= step_n;
            let quotient = (correlation(temperature, volume, &up, i, j)
                - 2.0 * correlation(temperature, volume, &moles, i, j)
                + correlation(temperature, volume, &down, i, j))
                / (step_n * step_n);
            assert_close(
                quotient,
                value.d_nn[component][component],
                1.0e-4,
                "the second mole number derivative",
            );
        }

        let step_n = moles[0] * 1.0e-3;
        let mut up = moles;
        let mut down = moles;
        up[0] += step_n;
        down[0] -= step_n;
        let mut up_other = moles;
        let mut down_other = moles;
        up_other[1] += step_n;
        down_other[1] -= step_n;
        let quotient = (correlation(temperature, volume, &[up[0], up_other[1]], i, j)
            - correlation(temperature, volume, &[up[0], down_other[1]], i, j)
            - correlation(temperature, volume, &[down[0], up_other[1]], i, j)
            + correlation(temperature, volume, &[down[0], down_other[1]], i, j))
            / (4.0 * step_n * step_n);
        assert_close(
            quotient,
            value.d_nn[0][1],
            1.0e-4,
            "the mixed second mole number derivative",
        );
    }

    /// The quotient and the derivative of the test agree relatively.
    fn assert_close(quotient: f64, derivative: f64, tolerance: f64, label: &str) {
        let scale = derivative.abs().max(1.0e-12);
        let difference = (quotient - derivative).abs() / scale;
        assert!(
            difference < tolerance,
            "{label} differs by the relative difference {difference}"
        );
    }
}
