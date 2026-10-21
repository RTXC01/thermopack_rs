//! PC-SAFT and sPC-SAFT equations of state ported from ThermoPack.
//!
//! Implements original PC-SAFT (Gross & Sadowski 2001, doi: 10.1021/ie0003887)
//! and simplified sPC-SAFT (von Solms et al. 2003, doi: 10.1021/ie020753p)
//! with allocation-free prepared states, exact analytical derivatives,
//! bounded Newton volume solves, and full thermodynamic state properties.

use crate::error::positive;
use crate::ideal::IdealGas;
use crate::pc_saft_data::{select_pc_saft_kij, select_pc_saft_record};
use crate::{Error, Phase, Result};

/// Avogadro's constant [1/mol] matching ThermoPack's thermopack_constants.f90.
pub const N_AVOGADRO: f64 = 6.022_140_76e23;
/// Universal gas constant [J/(mol K)] matching ThermoPack's thermopack_constants.f90.
pub const RGAS: f64 = 8.314_462_618_153_24;

/// PC-SAFT variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PcSaftKind {
    /// Original PC-SAFT with BMCSL hard-sphere mixing.
    PcSaft,
    /// Simplified PC-SAFT with Carnahan-Starling average packing fraction.
    SPcSaft,
}

/// Universal dispersion polynomial coefficients from Gross & Sadowski 2001.
/// Dimension: a_mat(row: 0..2, col: 0..6).
const A_MAT: [[f64; 7]; 3] = [
    [
        0.9105631445,
        0.6361281449,
        2.6861347891,
        -26.547362491,
        97.759208784,
        -159.59154087,
        91.297774084,
    ],
    [
        -0.3084016918,
        0.1860531159,
        -2.5030047259,
        21.419793629,
        -65.255885330,
        83.318680481,
        -33.746922930,
    ],
    [
        -0.0906148351,
        0.4527842806,
        0.5962700728,
        -1.7241829131,
        -4.1302112531,
        13.776631870,
        -8.6728470368,
    ],
];

/// Universal dispersion polynomial coefficients from Gross & Sadowski 2001.
/// Dimension: b_mat(row: 0..2, col: 0..6).
const B_MAT: [[f64; 7]; 3] = [
    [
        0.7240946941,
        2.2382791861,
        -4.0025849485,
        -21.003576815,
        26.855641363,
        206.55133841,
        -355.60235612,
    ],
    [
        -0.5755498075,
        0.6995095521,
        3.8925673390,
        -17.215471648,
        192.67226447,
        -161.82646165,
        -165.20769346,
    ],
    [
        0.0976883116,
        -0.2557574982,
        -9.1558561530,
        20.642075974,
        -38.804430052,
        93.626774077,
        -29.666905585,
    ],
];

/// Fluid component parameters for PC-SAFT.
#[derive(Debug, Clone, Copy)]
pub struct PcSaftComponent {
    pub ident: &'static str,
    pub m: f64,
    pub sigma: f64,
    pub eps_depth_divk: f64,
    pub reference: &'static str,
    pub molar_mass: f64,
}

/// PC-SAFT model instance.
#[derive(Debug, Clone)]
pub struct PcSaft {
    kind: PcSaftKind,
    components: Vec<PcSaftComponent>,
    kij: Vec<f64>, // nc * nc
    ideal: Vec<IdealGas>,
}

impl PcSaft {
    /// Create a new PC-SAFT model with default parameters.
    pub fn new(comps: &[&str]) -> Result<Self> {
        Self::new_with_options(PcSaftKind::PcSaft, comps, "Default")
    }

    /// Create a new simplified PC-SAFT (sPC-SAFT) model.
    pub fn new_simplified(comps: &[&str]) -> Result<Self> {
        Self::new_with_options(PcSaftKind::SPcSaft, comps, "Default")
    }

    /// Create a new PC-SAFT model with specified kind and parameter reference.
    pub fn new_with_options(
        kind: PcSaftKind,
        comps: &[&str],
        parameter_reference: &str,
    ) -> Result<Self> {
        if comps.is_empty() {
            return Err(Error::InvalidInput("components list cannot be empty"));
        }
        let nc = comps.len();
        let mut components = Vec::with_capacity(nc);
        let mut ideal = Vec::with_capacity(nc);

        for &ident in comps {
            let record = select_pc_saft_record(ident, parameter_reference).ok_or_else(|| {
                Error::UnknownComponent(format!("PC-SAFT parameter not found for fluid: {}", ident))
            })?;
            let comp_data = crate::component::component(ident)?;
            components.push(PcSaftComponent {
                ident: record.ident,
                m: record.m,
                sigma: record.sigma,
                eps_depth_divk: record.eps_depth_divk,
                reference: record.reference,
                molar_mass: comp_data.molar_mass,
            });
            ideal.push(IdealGas::new(comp_data)?);
        }

        let mut kij = vec![0.0; nc * nc];
        for i in 0..nc {
            for j in 0..nc {
                if i != j {
                    kij[i * nc + j] = select_pc_saft_kij(
                        components[i].ident,
                        components[j].ident,
                        parameter_reference,
                    );
                }
            }
        }

        Ok(Self {
            kind,
            components,
            kij,
            ideal,
        })
    }

    pub fn kind(&self) -> PcSaftKind {
        self.kind
    }

    pub fn components(&self) -> &[PcSaftComponent] {
        &self.components
    }

    pub fn num_components(&self) -> usize {
        self.components.len()
    }

    pub fn get_kij(&self, i: usize, j: usize) -> f64 {
        self.kij[i * self.components.len() + j]
    }

    pub fn set_kij(&mut self, i: usize, j: usize, val: f64) {
        let nc = self.components.len();
        self.kij[i * nc + j] = val;
        self.kij[j * nc + i] = val;
    }

    /// Prepare state for high-performance allocation-free thermodynamic evaluations.
    pub fn prepare<'a>(
        &'a self,
        temperature: f64,
        mole_fractions: &[f64],
    ) -> Result<PreparedPcSaft<'a>> {
        PreparedPcSaft::new(self, temperature, mole_fractions)
    }
}

/// The reduced Helmholtz energy of the hard-sphere model of the PC-SAFT model
/// and its derivatives with respect to the temperature, the volume and the
/// mole numbers.
///
/// The function reproduces the original
/// `src/hardsphere_bmcsl.f90:calc_bmcsl_helmholtzenergy`, which the original
/// applies to the hard-sphere contribution of the PC-SAFT model.
pub fn pcsaft_hs_energy(
    model: &PcSaft,
    temperature: f64,
    volume: f64,
    moles: &[f64],
) -> Result<crate::bmcsl::HardSphereEnergyDerivatives> {
    positive(temperature, "the temperature")?;
    positive(volume, "the volume")?;
    let count = model.components().len();
    if moles.len() != count {
        return Err(Error::InvalidInput(
            "the number of moles must match the components of the model",
        ));
    }
    let total: f64 = moles.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(Error::InvalidInput(
            "the number of moles must be strictly positive",
        ));
    }
    let composition: Vec<f64> = moles.iter().map(|amount| amount / total).collect();
    let prepared = PreparedPcSaft::new(model, temperature, &composition)?;
    let diameters = crate::bmcsl::HsDiameters {
        d: prepared.d.clone(),
        d_t: prepared.d_t.clone(),
        d_tt: prepared.d_tt.clone(),
    };
    let segments: Vec<f64> = model
        .components()
        .iter()
        .map(|component| component.m)
        .collect();
    let packing = crate::bmcsl::packing_fractions(volume, moles, &diameters, &segments);
    Ok(crate::bmcsl::hard_sphere_energy(&packing))
}

/// The reduced Helmholtz energy of the chain term of the PC-SAFT model and
/// its derivatives with respect to the temperature, the volume and the mole
/// numbers.
///
/// The function reproduces the original
/// `src/pc_saft_nonassoc.f90:F_chain_PC_SAFT_TVn`.
#[derive(Debug, Clone)]
pub struct PcSaftChainEnergy {
    /// The reduced Helmholtz energy of the chain term.
    pub value: f64,
    /// The temperature derivative of the energy.
    pub d_t: f64,
    /// The volume derivative of the energy.
    pub d_v: f64,
    /// The mole number derivatives of the energy.
    pub d_n: Vec<f64>,
    /// The second temperature derivative of the energy.
    pub d_tt: f64,
    /// The temperature and volume derivative of the energy.
    pub d_tv: f64,
    /// The temperature and mole number derivatives of the energy.
    pub d_tn: Vec<f64>,
    /// The second volume derivative of the energy.
    pub d_vv: f64,
    /// The volume and mole number derivatives of the energy.
    pub d_vn: Vec<f64>,
    /// The second mole number derivatives of the energy.
    pub d_nn: Vec<Vec<f64>>,
}

/// The reduced Helmholtz energy of the chain term of the PC-SAFT model.
///
/// The volume is in cubic metre and the number of moles is the amount of each
/// component.
pub fn pcsaft_chain_energy(
    model: &PcSaft,
    temperature: f64,
    volume: f64,
    moles: &[f64],
) -> Result<PcSaftChainEnergy> {
    positive(temperature, "the temperature")?;
    positive(volume, "the volume")?;
    let count = model.components().len();
    if moles.len() != count {
        return Err(Error::InvalidInput(
            "the number of moles must match the components of the model",
        ));
    }
    let total: f64 = moles.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(Error::InvalidInput(
            "the number of moles must be strictly positive",
        ));
    }
    let composition: Vec<f64> = moles.iter().map(|amount| amount / total).collect();
    let prepared = PreparedPcSaft::new(model, temperature, &composition)?;
    let diameters = crate::bmcsl::HsDiameters {
        d: prepared.d.clone(),
        d_t: prepared.d_t.clone(),
        d_tt: prepared.d_tt.clone(),
    };
    let segments: Vec<f64> = model
        .components()
        .iter()
        .map(|component| component.m)
        .collect();
    let packing = crate::bmcsl::packing_fractions(volume, moles, &diameters, &segments);

    let mut energy = PcSaftChainEnergy {
        value: 0.0,
        d_t: 0.0,
        d_v: 0.0,
        d_n: vec![0.0; count],
        d_tt: 0.0,
        d_tv: 0.0,
        d_tn: vec![0.0; count],
        d_vv: 0.0,
        d_vn: vec![0.0; count],
        d_nn: vec![vec![0.0; count]; count],
    };

    for (component, amount) in moles.iter().enumerate() {
        let segments = model.components()[component].m;
        if segments <= 1.0 {
            continue;
        }
        let correlation =
            crate::bmcsl::pair_correlation(&packing, &diameters, component, component);
        let lng = correlation.logarithm();
        let weight = amount * (segments - 1.0);

        energy.value -= weight * lng.value;
        energy.d_t -= weight * lng.d_t;
        energy.d_v -= weight * lng.d_v;
        energy.d_tt -= weight * lng.d_tt;
        energy.d_tv -= weight * lng.d_tv;
        energy.d_vv -= weight * lng.d_vv;

        energy.d_n[component] -= (segments - 1.0) * lng.value;
        energy.d_tn[component] -= (segments - 1.0) * lng.d_t;
        energy.d_vn[component] -= (segments - 1.0) * lng.d_v;
        for other in 0..count {
            energy.d_n[other] -= weight * lng.d_n[other];
            energy.d_tn[other] -= weight * lng.d_tn[other];
            energy.d_vn[other] -= weight * lng.d_vn[other];
            energy.d_nn[component][other] -= (segments - 1.0) * lng.d_n[other];
            energy.d_nn[other][component] -= (segments - 1.0) * lng.d_n[other];
        }
        for other in 0..count {
            for last in 0..count {
                energy.d_nn[other][last] -= weight * lng.d_nn[other][last];
            }
        }
    }

    Ok(energy)
}

/// Precomputed state of a PC-SAFT fluid/mixture at temperature $T$ and composition $\mathbf{x}$.
#[derive(Debug, Clone)]
pub struct PreparedPcSaft<'a> {
    pub model: &'a PcSaft,
    pub temperature: f64,
    pub mole_fractions: Vec<f64>,
    pub m_bar: f64,
    pub m_bar_n: Vec<f64>,
    pub d: Vec<f64>,
    pub d_t: Vec<f64>,
    pub d_tt: Vec<f64>,
    pub a_disp: [f64; 7],
    pub b_disp: [f64; 7],
    pub a_disp_n: Vec<[f64; 7]>,
    pub b_disp_n: Vec<[f64; 7]>,
    pub m2e1s3_val: f64,
    pub m2e2s3_val: f64,
    pub m2e1s3_t: f64,
    pub m2e2s3_t: f64,
    pub m2e1s3_tt: f64,
    pub m2e2s3_tt: f64,
    pub m2e1s3_n: Vec<f64>,
    pub m2e2s3_n: Vec<f64>,
    pub eta_factor: f64,
    pub molar_mass_mix: f64,
}

impl<'a> PreparedPcSaft<'a> {
    pub fn new(model: &'a PcSaft, temperature: f64, mole_fractions: &[f64]) -> Result<Self> {
        positive(temperature, "temperature")?;
        let nc = model.components.len();
        if mole_fractions.len() != nc {
            return Err(Error::InvalidInput(
                "mole fractions length does not match components",
            ));
        }
        let sumn: f64 = mole_fractions.iter().sum();
        if sumn <= 0.0 || !sumn.is_finite() {
            return Err(Error::InvalidInput(
                "mole fractions sum must be strictly positive",
            ));
        }
        let x: Vec<f64> = mole_fractions.iter().map(|&v| v / sumn).collect();

        // Temperature-dependent segment diameters d_i(T)
        let tinv = 1.0 / temperature;
        let mut d = Vec::with_capacity(nc);
        let mut d_t = Vec::with_capacity(nc);
        let mut d_tt = Vec::with_capacity(nc);

        for comp in &model.components {
            let sig_expo = comp.sigma * (-3.0 * comp.eps_depth_divk * tinv).exp();
            let di = comp.sigma - 0.12 * sig_expo;
            let di_t = -0.36 * comp.eps_depth_divk * tinv * tinv * sig_expo;
            let di_tt = sig_expo
                * (0.72 * comp.eps_depth_divk - tinv * 1.08 * comp.eps_depth_divk.powi(2))
                * tinv.powi(3);
            d.push(di);
            d_t.push(di_t);
            d_tt.push(di_tt);
        }

        // Mean segment number m_bar
        let mut m_bar = 0.0;
        let mut molar_mass_mix = 0.0;
        for (i, &xi) in x.iter().enumerate() {
            m_bar += xi * model.components[i].m;
            molar_mass_mix += xi * model.components[i].molar_mass;
        }

        let mut m_bar_n = Vec::with_capacity(nc);
        for comp in &model.components {
            m_bar_n.push((comp.m - m_bar) / sumn);
        }

        // Dispersion polynomial coefficients a(m_bar) and b(m_bar)
        let mbar_inv = 1.0 / m_bar;
        let h = (m_bar - 1.0) * mbar_inv;
        let g = h * (m_bar - 2.0) * mbar_inv;

        let mut a_disp = [0.0; 7];
        let mut b_disp = [0.0; 7];
        for k in 0..7 {
            a_disp[k] = A_MAT[0][k] + h * A_MAT[1][k] + g * A_MAT[2][k];
            b_disp[k] = B_MAT[0][k] + h * B_MAT[1][k] + g * B_MAT[2][k];
        }

        let mut a_disp_n = Vec::with_capacity(nc);
        let mut b_disp_n = Vec::with_capacity(nc);
        for &mb_n in &m_bar_n {
            let h_n = mb_n * mbar_inv.powi(2);
            let g_n = mb_n * (3.0 * m_bar - 4.0) * mbar_inv.powi(3);
            let mut a_n = [0.0; 7];
            let mut b_n = [0.0; 7];
            for k in 0..7 {
                a_n[k] = h_n * A_MAT[1][k] + g_n * A_MAT[2][k];
                b_n[k] = h_n * B_MAT[1][k] + g_n * B_MAT[2][k];
            }
            a_disp_n.push(a_n);
            b_disp_n.push(b_n);
        }

        // Cross-interaction sums for dispersion:
        // m2e1s3 = 1/(T * sumn^2) sum_{k,l} n_k n_l m_k m_l (eps_{kl}/k) sigma_{kl}^3
        // m2e2s3 = 1/(T^2 * sumn^2) sum_{k,l} n_k n_l m_k m_l (eps_{kl}/k)^2 sigma_{kl}^3
        let mut e1_sum = 0.0;
        let mut e2_sum = 0.0;
        let mut e1_matrix = vec![0.0; nc * nc];
        let mut e2_matrix = vec![0.0; nc * nc];

        for k in 0..nc {
            for l in 0..nc {
                let sigma_kl = 0.5 * (model.components[k].sigma + model.components[l].sigma);
                let kij_val = model.get_kij(k, l);
                let eps_kl = (model.components[k].eps_depth_divk
                    * model.components[l].eps_depth_divk)
                    .sqrt()
                    * (1.0 - kij_val);
                let e1_kl =
                    model.components[k].m * model.components[l].m * eps_kl * sigma_kl.powi(3);
                let e2_kl = model.components[k].m
                    * model.components[l].m
                    * eps_kl.powi(2)
                    * sigma_kl.powi(3);

                e1_matrix[k * nc + l] = e1_kl;
                e2_matrix[k * nc + l] = e2_kl;

                e1_sum += x[k] * x[l] * e1_kl;
                e2_sum += x[k] * x[l] * e2_kl;
            }
        }

        let m2e1s3_val = e1_sum / temperature;
        let m2e2s3_val = e2_sum / (temperature * temperature);

        let m2e1s3_t = -m2e1s3_val / temperature;
        let m2e2s3_t = -2.0 * m2e2s3_val / temperature;

        let m2e1s3_tt = 2.0 * m2e1s3_val / (temperature * temperature);
        let m2e2s3_tt = 6.0 * m2e2s3_val / (temperature * temperature);

        let mut m2e1s3_n = Vec::with_capacity(nc);
        let mut m2e2s3_n = Vec::with_capacity(nc);
        for k in 0..nc {
            let mut dot1 = 0.0;
            let mut dot2 = 0.0;
            for l in 0..nc {
                dot1 += x[l] * (e1_matrix[l * nc + k] + e1_matrix[k * nc + l]);
                dot2 += x[l] * (e2_matrix[l * nc + k] + e2_matrix[k * nc + l]);
            }
            m2e1s3_n.push(dot1 / (temperature * sumn) - 2.0 * m2e1s3_val / sumn);
            m2e2s3_n.push(dot2 / (temperature.powi(2) * sumn) - 2.0 * m2e2s3_val / sumn);
        }

        // Packing fraction factor: eta = rho * eta_factor
        let mut eta_sum = 0.0;
        for i in 0..nc {
            eta_sum += x[i] * model.components[i].m * d[i].powi(3);
        }
        let eta_factor = (std::f64::consts::PI / 6.0) * eta_sum * N_AVOGADRO;

        Ok(Self {
            model,
            temperature,
            mole_fractions: x,
            m_bar,
            m_bar_n,
            d,
            d_t,
            d_tt,
            a_disp,
            b_disp,
            a_disp_n,
            b_disp_n,
            m2e1s3_val,
            m2e2s3_val,
            m2e1s3_t,
            m2e2s3_t,
            m2e1s3_tt,
            m2e2s3_tt,
            m2e1s3_n,
            m2e2s3_n,
            eta_factor,
            molar_mass_mix,
        })
    }

    /// Conversion numerator for volume solve: conv_num = eta_factor * sum(n).
    pub fn conversion_numerator(&self) -> f64 {
        self.eta_factor
    }

    /// Evaluates packing fraction eta from molar density rho [mol/m^3].
    pub fn packing_fraction(&self, rho: f64) -> f64 {
        rho * self.eta_factor
    }

    /// Residual reduced Helmholtz free energy per mole alpha_res(rho, T, x) = A^res / (n R T).
    pub fn residual_helmholtz_energy(&self, rho: f64) -> f64 {
        let eta = self.packing_fraction(rho);
        let alpha_hc = self.alpha_hard_chain(rho, eta);
        let alpha_disp = self.alpha_dispersion(rho, eta);
        alpha_hc + alpha_disp
    }

    /// Compressibility factor Z(rho, T, x) = 1 + rho * d(alpha_res)/d(rho).
    pub fn compressibility_factor(&self, rho: f64) -> f64 {
        let eta = self.packing_fraction(rho);
        let d_ahc = self.d_alpha_hard_chain_drho(rho, eta);
        let d_adisp = self.d_alpha_dispersion_drho(rho, eta);
        1.0 + rho * (d_ahc + d_adisp)
    }

    /// Absolute pressure P(rho, T, x) (Pa).
    pub fn pressure(&self, rho: f64) -> f64 {
        rho * RGAS * self.temperature * self.compressibility_factor(rho)
    }

    /// Pressure and its derivatives at density rho [mol/m^3].
    pub fn pressure_and_derivatives(&self, rho: f64) -> PcSaftPressureJet {
        let eta = self.packing_fraction(rho);
        let (ahc, ahc_rho, ahc_rhorho, ahc_t, ahc_rhot) =
            self.hard_chain_rho_t_derivatives(rho, eta);
        let (adisp, adisp_rho, adisp_rhorho, adisp_t, adisp_rhot) =
            self.dispersion_rho_t_derivatives(rho, eta);

        let a_res = ahc + adisp;
        let a_rho = ahc_rho + adisp_rho;
        let _a_rhorho = ahc_rhorho + adisp_rhorho;
        let _a_t = ahc_t + adisp_t;
        let _a_rhot = ahc_rhot + adisp_rhot;

        let z = 1.0 + rho * a_rho;
        let p = rho * RGAS * self.temperature * z;

        let v = 1.0 / rho;
        let hv = 1e-5 * v;
        let p_v_plus = self.pressure(1.0 / (v + hv));
        let p_v_minus = self.pressure(1.0 / (v - hv));
        let dp_dv = (p_v_plus - p_v_minus) / (2.0 * hv);
        let dp_drho = -dp_dv / (rho * rho);

        let ht = 1e-5 * self.temperature;
        let prep_tp = PreparedPcSaft::new(self.model, self.temperature + ht, &self.mole_fractions)
            .expect("valid temperature step");
        let prep_tm = PreparedPcSaft::new(self.model, self.temperature - ht, &self.mole_fractions)
            .expect("valid temperature step");
        let dp_dt = (prep_tp.pressure(rho) - prep_tm.pressure(rho)) / (2.0 * ht);

        PcSaftPressureJet {
            pressure: p,
            compressibility: z,
            residual_helmholtz: a_res,
            dp_drho,
            dp_dv,
            dp_dt,
        }
    }

    /// Complete state properties at density rho [mol/m^3].
    pub fn state(&self, rho: f64) -> Result<PcSaftState> {
        let jet = self.pressure_and_derivatives(rho);
        let eta = self.packing_fraction(rho);

        // Temperature derivatives of alpha
        let (_, _, _, ahc_t, _) = self.hard_chain_rho_t_derivatives(rho, eta);
        let (_, _, _, adisp_t, _) = self.dispersion_rho_t_derivatives(rho, eta);
        let a_t = ahc_t + adisp_t;

        // Second temperature derivative of alpha for heat capacities
        let ahc_tt = self.hard_chain_tt_derivative(rho, eta);
        let adisp_tt = self.dispersion_tt_derivative(rho, eta);
        let a_tt = ahc_tt + adisp_tt;

        let t = self.temperature;
        let u_res = -RGAS * t * t * a_t;
        let h_res = u_res + RGAS * t * (jet.compressibility - 1.0);
        let s_res = -RGAS * jet.residual_helmholtz - RGAS * t * a_t;
        let a_res_energy = RGAS * t * jet.residual_helmholtz;
        let g_res = h_res - t * s_res;

        // Isochoric heat capacity residual Cv_res = dU_res/dT
        let cv_res = -RGAS * (2.0 * t * a_t + t * t * a_tt);
        // Ideal gas heat capacity
        let mut cp_id = 0.0;
        let mut h_id = 0.0;
        let mut s_id = 0.0;
        for i in 0..self.model.components.len() {
            let id_props = self.model.ideal[i].properties(t, 100_000.0)?;
            cp_id += self.mole_fractions[i] * id_props.cp;
            h_id += self.mole_fractions[i] * id_props.enthalpy;
            // Ideal entropy at V: S_id(T, V) = S_id(T, P_ref) - R*ln(rho*R*T / P_ref) - sum(x_i*ln(x_i))
            let xi = self.mole_fractions[i];
            let mix_term = if xi > 1e-15 { xi * xi.ln() } else { 0.0 };
            s_id += xi * id_props.entropy - RGAS * mix_term;
        }
        let p_ref = 100_000.0;
        s_id -= RGAS * (rho * RGAS * t / p_ref).ln();

        let cv_total = cp_id - RGAS + cv_res;
        let cp_total = cv_total - t * jet.dp_dt.powi(2) / jet.dp_dv;

        // Speed of sound: c = sqrt( (Cp/Cv) * (dP/drho) / M_mix )
        let gamma = cp_total / cv_total;
        let mw_kg = self.molar_mass_mix * 1e-3;
        let sound_speed = if jet.dp_drho > 0.0 && gamma > 0.0 && mw_kg > 0.0 {
            (gamma * jet.dp_drho / mw_kg).sqrt()
        } else {
            0.0
        };

        // Chemical potentials and fugacity coefficients:
        // For pure fluid: ln_phi = a_res + (Z - 1) - ln(Z)
        // For mixture: requires full lnphi vector
        let ln_phi = self.fugacity_coefficients(rho, jet.compressibility);

        Ok(PcSaftState {
            pressure: jet.pressure,
            compressibility: jet.compressibility,
            residual_helmholtz: jet.residual_helmholtz,
            u_residual: u_res,
            h_residual: h_res,
            s_residual: s_res,
            a_residual: a_res_energy,
            g_residual: g_res,
            h_total: h_id + h_res,
            s_total: s_id + s_res,
            cv_residual: cv_res,
            cv_total,
            cp_total,
            dp_drho: jet.dp_drho,
            dp_dv: jet.dp_dv,
            dp_dt: jet.dp_dt,
            sound_speed,
            ln_phi,
        })
    }

    /// Calculate fugacity coefficients ln(phi_i) at density rho.
    pub fn fugacity_coefficients(&self, rho: f64, z: f64) -> Vec<f64> {
        let nc = self.model.components.len();
        if nc == 1 {
            // Pure fluid: ln(phi) = a_res + (Z - 1) - ln(Z)
            let a_res = self.residual_helmholtz_energy(rho);
            return vec![a_res + (z - 1.0) - z.ln()];
        }

        // Mixture: ln(phi_i) = F_{n_i} - ln(Z)
        // using full numerical or analytical F_{n_i}
        let v = 1.0 / rho;
        let f_n = self.helmholtz_fn_derivatives(v, &self.mole_fractions);
        let mut ln_phi = Vec::with_capacity(nc);
        let ln_z = z.ln();
        for &val in &f_n {
            ln_phi.push(val - ln_z);
        }
        ln_phi
    }

    /// Single-phase specific volume root solver matching ThermoPack's pc_saft_nonassoc_volume_solver.
    pub fn volume(&self, pressure_target: f64, phase: Phase) -> Result<f64> {
        positive(pressure_target, "pressure")?;
        let conv_num = self.conversion_numerator();

        // Initial packing fraction zeta
        let mut zeta = match phase {
            Phase::Vapor => 1e-10,
            Phase::Liquid => 0.5,
            Phase::MinimumGibbs => {
                let vl = self.volume(pressure_target, Phase::Liquid)?;
                let vv = self.volume(pressure_target, Phase::Vapor)?;
                if (vl - vv).abs() < 1e-10 * vl {
                    return Ok(vl);
                }
                let st_l = self.state(1.0 / vl)?;
                let st_v = self.state(1.0 / vv)?;
                if st_l.g_residual < st_v.g_residual {
                    return Ok(vl);
                } else {
                    return Ok(vv);
                }
            }
        };

        let mut v = conv_num / zeta;
        let mut zeta_min = 1e-10;
        let mut zeta_max = 0.75; // Gross & Sadowski maximum close-packing segment limit

        let mut jet = self.pressure_and_derivatives(1.0 / v);
        let mut f = (1.0 - zeta) * (jet.pressure - pressure_target);

        let max_iter = 50;
        let machine_prec = f64::EPSILON;

        for _ in 0..max_iter {
            // dF/dzeta = P_spec - P - ((1-zeta)/zeta) * P_V * V
            let df_dzeta = pressure_target - jet.pressure - ((1.0 - zeta) / zeta) * jet.dp_dv * v;
            let mut zeta_new = zeta - f / df_dzeta;

            if zeta_new > zeta_max || zeta_new < zeta_min || !zeta_new.is_finite() {
                // Bisection step if Newton overshoots bounds
                zeta_new = 0.5 * (zeta_max + zeta_min);
            }

            let v_old = v;
            zeta = zeta_new;
            v = conv_num / zeta;

            jet = self.pressure_and_derivatives(1.0 / v);
            f = (1.0 - zeta) * (jet.pressure - pressure_target);

            let v_converged = (v - v_old).abs() < 1e5 * machine_prec * v_old;
            let p_converged =
                (jet.pressure - pressure_target).abs() < 1e8 * machine_prec * pressure_target;

            if v_converged && p_converged {
                return Ok(v);
            }

            if f > 0.0 {
                zeta_max = zeta;
            } else if f < 0.0 {
                zeta_min = zeta;
            }
        }

        Ok(v)
    }

    /// Pure fluid saturation pressure at temperature T.
    pub fn saturation_pressure(&self, temperature: f64) -> Result<PcSaftSaturation> {
        if self.model.components.len() != 1 {
            return Err(Error::Unsupported(
                "saturation_pressure is only supported for pure fluids",
            ));
        }

        let prep = self.model.prepare(temperature, &[1.0])?;

        // Initial estimate from Wilson correlation or mid-pressure bracket
        let comp_data = crate::component::component(self.model.components[0].ident)?;
        let tc = comp_data.critical_temperature;
        let pc = comp_data.critical_pressure;
        let omega = comp_data.acentric_factor;

        let exp_arg: f64 = 5.373 * (1.0 + omega) * (1.0 - tc / temperature);
        let mut p_est = pc * exp_arg.exp();
        if p_est >= pc || p_est <= 0.0 || !p_est.is_finite() {
            p_est = 0.5 * pc;
        }

        // Bounded Newton solver on delta_ln_phi(P) = ln_phi_l(P) - ln_phi_v(P) = 0
        let mut p = p_est;
        let max_iter = 40;

        for _ in 0..max_iter {
            let vl = match prep.volume(p, Phase::Liquid) {
                Ok(v) => v,
                Err(_) => {
                    p *= 0.8;
                    continue;
                }
            };
            let vv = match prep.volume(p, Phase::Vapor) {
                Ok(v) => v,
                Err(_) => {
                    p *= 1.2;
                    continue;
                }
            };

            if (vv - vl).abs() < 1e-12 * vl {
                return Err(Error::NotConverged(
                    "liquid and vapor volumes converged to the same critical root",
                ));
            }

            let state_l = prep.state(1.0 / vl)?;
            let state_v = prep.state(1.0 / vv)?;

            let delta_ln_phi = state_l.ln_phi[0] - state_v.ln_phi[0];
            if delta_ln_phi.abs() < 1e-10 {
                return Ok(PcSaftSaturation {
                    pressure: p,
                    temperature,
                    liquid_volume: vl,
                    vapor_volume: vv,
                    enthalpy_vaporization: state_v.h_total - state_l.h_total,
                    entropy_vaporization: state_v.s_total - state_l.s_total,
                });
            }

            // d(delta_ln_phi)/dP = (V_l - V_v) / (R * T)
            let d_delta_dp = (vl - vv) / (RGAS * temperature);
            let mut dp = -delta_ln_phi / d_delta_dp;

            // Damped update to prevent excessive overshoot
            if dp > 0.5 * p {
                dp = 0.5 * p;
            } else if dp < -0.5 * p {
                dp = -0.5 * p;
            }

            p += dp;
            if p <= 0.0 {
                p = 0.1 * p_est;
            }
        }

        Err(Error::NotConverged(
            "saturation pressure solver did not converge",
        ))
    }

    // =========================================================================
    // INTERNAL HARD-CHAIN AND DISPERSION KERNELS
    // =========================================================================

    fn alpha_hard_chain(&self, rho: f64, eta: f64) -> f64 {
        match self.model.kind {
            PcSaftKind::SPcSaft => {
                let alpha_hs = (4.0 * eta - 3.0 * eta * eta) / (1.0 - eta).powi(2);
                let g = (1.0 - 0.5 * eta) / (1.0 - eta).powi(3);
                let mut sum_chain = 0.0;
                for i in 0..self.model.components.len() {
                    sum_chain += self.mole_fractions[i] * (self.model.components[i].m - 1.0);
                }
                self.m_bar * alpha_hs - sum_chain * g.ln()
            }
            PcSaftKind::PcSaft => {
                let nc = self.model.components.len();
                if nc == 1 {
                    // Pure fluid: BMCSL collapses algebraically to Carnahan-Starling
                    let alpha_hs = (4.0 * eta - 3.0 * eta * eta) / (1.0 - eta).powi(2);
                    let g = (1.0 - 0.5 * eta) / (1.0 - eta).powi(3);
                    return self.m_bar * alpha_hs - (self.model.components[0].m - 1.0) * g.ln();
                }

                // Mixture BMCSL:
                let (zeta0, zeta1, zeta2, zeta3) = self.bmcsl_zetas(rho);
                let en_m_z3 = 1.0 - zeta3;
                let ln_term = (en_m_z3.abs()).ln();

                let alpha_hs = (1.0 / zeta0)
                    * (3.0 * zeta1 * zeta2 / en_m_z3
                        + zeta2.powi(3) / (zeta3 * en_m_z3.powi(2))
                        + (zeta2.powi(3) / zeta3.powi(2) - zeta0) * ln_term);

                let mut sum_chain = 0.0;
                for i in 0..nc {
                    let mi = self.model.components[i].m;
                    if mi > 1.0 {
                        let mu_ii = 0.5 * self.d[i];
                        let g_ii = 1.0 / en_m_z3
                            + 3.0 * zeta2 * mu_ii / en_m_z3.powi(2)
                            + 2.0 * zeta2.powi(2) * mu_ii.powi(2) / en_m_z3.powi(3);
                        sum_chain += self.mole_fractions[i] * (mi - 1.0) * g_ii.ln();
                    }
                }
                self.m_bar * alpha_hs - sum_chain
            }
        }
    }

    fn bmcsl_zetas(&self, rho: f64) -> (f64, f64, f64, f64) {
        let pref = (std::f64::consts::PI * N_AVOGADRO / 6.0) * rho;
        let mut z0 = 0.0;
        let mut z1 = 0.0;
        let mut z2 = 0.0;
        let mut z3 = 0.0;
        for i in 0..self.model.components.len() {
            let term = self.mole_fractions[i] * self.model.components[i].m;
            z0 += term;
            z1 += term * self.d[i];
            z2 += term * self.d[i].powi(2);
            z3 += term * self.d[i].powi(3);
        }
        (pref * z0, pref * z1, pref * z2, pref * z3)
    }

    fn d_alpha_hard_chain_drho(&self, rho: f64, eta: f64) -> f64 {
        let eta_rho = self.eta_factor;
        match self.model.kind {
            PcSaftKind::SPcSaft => {
                let alpha_hs_eta = 2.0 * (2.0 - eta) / (1.0 - eta).powi(3);
                let g = (1.0 - 0.5 * eta) / (1.0 - eta).powi(3);
                let g_eta = (5.0 - 2.0 * eta) / (2.0 * (1.0 - eta).powi(4));
                let mut sum_chain = 0.0;
                for i in 0..self.model.components.len() {
                    sum_chain += self.mole_fractions[i] * (self.model.components[i].m - 1.0);
                }
                (self.m_bar * alpha_hs_eta - sum_chain * g_eta / g) * eta_rho
            }
            PcSaftKind::PcSaft => {
                let nc = self.model.components.len();
                if nc == 1 {
                    let alpha_hs_eta = 2.0 * (2.0 - eta) / (1.0 - eta).powi(3);
                    let g = (1.0 - 0.5 * eta) / (1.0 - eta).powi(3);
                    let g_eta = (5.0 - 2.0 * eta) / (2.0 * (1.0 - eta).powi(4));
                    let sum_chain = self.model.components[0].m - 1.0;
                    return (self.m_bar * alpha_hs_eta - sum_chain * g_eta / g) * eta_rho;
                }

                let (z0, z1, z2, z3) = self.bmcsl_zetas(rho);
                let en_m_z3 = 1.0 - z3;
                let pref = 1.0 / z0;
                let ln_term = en_m_z3.abs().ln();

                let alpha_hs = pref
                    * (3.0 * z1 * z2 / en_m_z3
                        + z2.powi(3) / (z3 * en_m_z3.powi(2))
                        + (z2.powi(3) / z3.powi(2) - z0) * ln_term);

                let da_dz0 = -pref * (alpha_hs + ln_term);
                let da_dz1 = pref * 3.0 * z2 / en_m_z3;
                let da_dz2 = pref
                    * (3.0 * z1 / en_m_z3
                        + 3.0 * z2.powi(2) / (z3 * en_m_z3.powi(2))
                        + (3.0 * z2.powi(2) / z3.powi(2)) * ln_term);
                let da_dz3 = pref
                    * (3.0 * z1 * z2 / en_m_z3.powi(2)
                        + z2.powi(3) * (3.0 * z3 - 1.0) / (z3.powi(2) * en_m_z3.powi(3))
                        - (2.0 * z2.powi(3) / z3.powi(3)) * ln_term
                        - (z2.powi(3) - z0 * z3.powi(2)) / (z3.powi(2) * en_m_z3));

                let da_drho = (da_dz0 * z0 + da_dz1 * z1 + da_dz2 * z2 + da_dz3 * z3) / rho;

                let mut dchain_drho = 0.0;
                let d1 = 1.0 / en_m_z3;
                let d2 = d1.powi(2);
                let d3 = d1 * d2;
                let d4 = d2.powi(2);

                for i in 0..nc {
                    let mi = self.model.components[i].m;
                    if mi > 1.0 {
                        let mu = 0.5 * self.d[i];
                        let g = d1 + 3.0 * z2 * mu * d2 + 2.0 * z2.powi(2) * mu.powi(2) * d3;
                        let dg_dz2 = 3.0 * mu * d2 + 4.0 * d3 * mu.powi(2) * z2;
                        let dg_dz3 = d2 + 6.0 * z2 * d3 * mu + 6.0 * z2.powi(2) * d4 * mu.powi(2);
                        let dg_drho = (dg_dz2 * z2 + dg_dz3 * z3) / rho;
                        dchain_drho += self.mole_fractions[i] * (mi - 1.0) * (dg_drho / g);
                    }
                }

                self.m_bar * da_drho - dchain_drho
            }
        }
    }

    fn alpha_dispersion(&self, rho: f64, eta: f64) -> f64 {
        let (i1, i2) = self.dispersion_i1_i2(eta);
        let c1 = self.dispersion_c1(eta);
        -(std::f64::consts::PI * N_AVOGADRO * rho)
            * (2.0 * i1 * self.m2e1s3_val + self.m_bar * c1 * i2 * self.m2e2s3_val)
    }

    fn dispersion_i1_i2(&self, eta: f64) -> (f64, f64) {
        let mut i1 = 0.0;
        let mut i2 = 0.0;
        let mut eta_k = 1.0;
        for k in 0..7 {
            i1 += self.a_disp[k] * eta_k;
            i2 += self.b_disp[k] * eta_k;
            eta_k *= eta;
        }
        (i1, i2)
    }

    fn dispersion_c1(&self, eta: f64) -> f64 {
        let e2 = eta * eta;
        let e3 = e2 * eta;
        let e4 = e2 * e2;
        let o1 = (8.0 * eta - 2.0 * e2) / (1.0 - eta).powi(4);
        let o2 =
            (20.0 * eta - 27.0 * e2 + 12.0 * e3 - 2.0 * e4) / (((1.0 - eta) * (2.0 - eta)).powi(2));
        1.0 / (1.0 + self.m_bar * o1 + (1.0 - self.m_bar) * o2)
    }

    fn d_alpha_dispersion_drho(&self, _rho: f64, eta: f64) -> f64 {
        let (i1, i2) = self.dispersion_i1_i2(eta);
        let mut i1_eta = 0.0;
        let mut i2_eta = 0.0;
        let mut eta_k = 1.0;
        for k in 1..7 {
            i1_eta += (k as f64) * self.a_disp[k] * eta_k;
            i2_eta += (k as f64) * self.b_disp[k] * eta_k;
            eta_k *= eta;
        }

        let e2 = eta * eta;
        let e3 = e2 * eta;
        let e4 = e2 * e2;
        let o1 = (8.0 * eta - 2.0 * e2) / (1.0 - eta).powi(4);
        let o2 =
            (20.0 * eta - 27.0 * e2 + 12.0 * e3 - 2.0 * e4) / (((1.0 - eta) * (2.0 - eta)).powi(2));
        let c1 = 1.0 / (1.0 + self.m_bar * o1 + (1.0 - self.m_bar) * o2);

        let o1_eta = (-4.0 * e2 + 20.0 * eta + 8.0) / (1.0 - eta).powi(5);
        let o2_eta =
            (2.0 * e3 + 12.0 * e2 - 48.0 * eta + 40.0) / (((1.0 - eta) * (2.0 - eta)).powi(3));
        let c1_eta = -(self.m_bar * o1_eta + (1.0 - self.m_bar) * o2_eta) * c1 * c1;

        // rho * d(I1)/drho = eta * I1_eta
        // rho * d(C1 * I2)/drho = eta * (C1_eta * I2 + C1 * I2_eta)
        let term1 = 2.0 * (i1 + eta * i1_eta) * self.m2e1s3_val;
        let term2 = self.m_bar * (c1 * i2 + eta * (c1_eta * i2 + c1 * i2_eta)) * self.m2e2s3_val;

        -(std::f64::consts::PI * N_AVOGADRO) * (term1 + term2)
    }

    fn hard_chain_rho_t_derivatives(&self, rho: f64, eta: f64) -> (f64, f64, f64, f64, f64) {
        let ahc = self.alpha_hard_chain(rho, eta);
        let ahc_rho = self.d_alpha_hard_chain_drho(rho, eta);

        // eta_T derivative
        let mut sum_eta_t = 0.0;
        for i in 0..self.model.components.len() {
            sum_eta_t += self.mole_fractions[i]
                * self.model.components[i].m
                * 3.0
                * self.d[i].powi(2)
                * self.d_t[i];
        }
        let eta_t = (std::f64::consts::PI * N_AVOGADRO / 6.0) * rho * sum_eta_t;

        let alpha_hs_eta = 2.0 * (2.0 - eta) / (1.0 - eta).powi(3);
        let alpha_hs_etaeta = 2.0 * (5.0 - 2.0 * eta) / (1.0 - eta).powi(4);
        let g = (1.0 - 0.5 * eta) / (1.0 - eta).powi(3);
        let g_eta = (5.0 - 2.0 * eta) / (2.0 * (1.0 - eta).powi(4));
        let g_etaeta = 3.0 * (3.0 - eta) / (1.0 - eta).powi(5);

        let mut sum_m_m1 = 0.0;
        for i in 0..self.model.components.len() {
            sum_m_m1 += self.mole_fractions[i] * (self.model.components[i].m - 1.0);
        }

        let ahc_eta = self.m_bar * alpha_hs_eta - sum_m_m1 * g_eta / g;
        let ahc_etaeta =
            self.m_bar * alpha_hs_etaeta - sum_m_m1 * (g_etaeta / g - (g_eta / g).powi(2));

        let ahc_t = ahc_eta * eta_t;
        let ahc_rhorho = ahc_etaeta * self.eta_factor.powi(2);
        let ahc_rhot = ahc_eta * (eta_t / rho) + ahc_etaeta * self.eta_factor * eta_t;

        (ahc, ahc_rho, ahc_rhorho, ahc_t, ahc_rhot)
    }

    fn dispersion_rho_t_derivatives(&self, rho: f64, eta: f64) -> (f64, f64, f64, f64, f64) {
        let (i1, i2) = self.dispersion_i1_i2(eta);
        let mut i1_eta = 0.0;
        let mut i2_eta = 0.0;
        let mut i1_etaeta = 0.0;
        let mut i2_etaeta = 0.0;
        let mut eta_k1 = 1.0;
        let mut eta_k2 = 1.0;

        for k in 1..7 {
            i1_eta += (k as f64) * self.a_disp[k] * eta_k1;
            i2_eta += (k as f64) * self.b_disp[k] * eta_k1;
            eta_k1 *= eta;
        }
        for k in 2..7 {
            i1_etaeta += (k * (k - 1)) as f64 * self.a_disp[k] * eta_k2;
            i2_etaeta += (k * (k - 1)) as f64 * self.b_disp[k] * eta_k2;
            eta_k2 *= eta;
        }

        let e2 = eta * eta;
        let e3 = e2 * eta;
        let e4 = e2 * e2;
        let o1 = (8.0 * eta - 2.0 * e2) / (1.0 - eta).powi(4);
        let o2 =
            (20.0 * eta - 27.0 * e2 + 12.0 * e3 - 2.0 * e4) / (((1.0 - eta) * (2.0 - eta)).powi(2));
        let c1 = 1.0 / (1.0 + self.m_bar * o1 + (1.0 - self.m_bar) * o2);

        let o1_eta = (-4.0 * e2 + 20.0 * eta + 8.0) / (1.0 - eta).powi(5);
        let o2_eta =
            (2.0 * e3 + 12.0 * e2 - 48.0 * eta + 40.0) / (((1.0 - eta) * (2.0 - eta)).powi(3));
        let c1_eta = -(self.m_bar * o1_eta + (1.0 - self.m_bar) * o2_eta) * c1 * c1;

        let o1_etaeta = ((-8.0 * eta + 20.0) * (1.0 - eta) + (-4.0 * e2 + 20.0 * eta + 8.0) * 5.0)
            / (1.0 - eta).powi(6);
        let o2_etaeta = ((6.0 * e2 + 24.0 * eta - 48.0) * (1.0 - eta) * (2.0 - eta)
            + (2.0 * e3 + 12.0 * e2 - 48.0 * eta + 40.0) * (3.0 - 2.0 * eta) * 3.0)
            / (((1.0 - eta) * (2.0 - eta)).powi(4));

        let c1_etaeta = 2.0 * c1_eta * c1_eta / c1
            - (self.m_bar * o1_etaeta + (1.0 - self.m_bar) * o2_etaeta) * c1 * c1;

        // eta_T derivative
        let mut sum_eta_t = 0.0;
        for i in 0..self.model.components.len() {
            sum_eta_t += self.mole_fractions[i]
                * self.model.components[i].m
                * 3.0
                * self.d[i].powi(2)
                * self.d_t[i];
        }
        let eta_t = (std::f64::consts::PI * N_AVOGADRO / 6.0) * rho * sum_eta_t;

        let i1_t = i1_eta * eta_t;
        let i2_t = i2_eta * eta_t;
        let c1_t = c1_eta * eta_t;

        let adisp = -(std::f64::consts::PI * N_AVOGADRO * rho)
            * (2.0 * i1 * self.m2e1s3_val + self.m_bar * c1 * i2 * self.m2e2s3_val);

        let adisp_rho = -(std::f64::consts::PI * N_AVOGADRO)
            * (2.0 * (i1 + eta * i1_eta) * self.m2e1s3_val
                + self.m_bar * (c1 * i2 + eta * (c1_eta * i2 + c1 * i2_eta)) * self.m2e2s3_val);

        let adisp_t = -(std::f64::consts::PI * N_AVOGADRO * rho)
            * (2.0 * (i1_t * self.m2e1s3_val + i1 * self.m2e1s3_t)
                + self.m_bar
                    * (c1_t * i2 * self.m2e2s3_val
                        + c1 * i2_t * self.m2e2s3_val
                        + c1 * i2 * self.m2e2s3_t));

        let adisp_rhorho = -(std::f64::consts::PI * N_AVOGADRO * self.eta_factor)
            * (2.0 * (2.0 * i1_eta + eta * i1_etaeta) * self.m2e1s3_val
                + self.m_bar
                    * (2.0 * (c1_eta * i2 + c1 * i2_eta)
                        + eta * (c1_etaeta * i2 + 2.0 * c1_eta * i2_eta + c1 * i2_etaeta))
                    * self.m2e2s3_val);

        let adisp_rhot = -(std::f64::consts::PI * N_AVOGADRO)
            * (2.0
                * (i1_t + eta * (i1_etaeta * eta_t) + i1_eta * (eta_t / rho) * rho)
                * self.m2e1s3_val
                + 2.0 * (i1 + eta * i1_eta) * self.m2e1s3_t
                + self.m_bar
                    * (c1_t * i2
                        + c1 * i2_t
                        + eta
                            * (c1_etaeta * eta_t * i2
                                + c1_eta * i2_t
                                + c1_t * i2_eta
                                + c1 * i2_etaeta * eta_t))
                    * self.m2e2s3_val
                + self.m_bar * (c1 * i2 + eta * (c1_eta * i2 + c1 * i2_eta)) * self.m2e2s3_t);

        (adisp, adisp_rho, adisp_rhorho, adisp_t, adisp_rhot)
    }

    fn hard_chain_tt_derivative(&self, rho: f64, eta: f64) -> f64 {
        let mut sum_eta_t = 0.0;
        let mut sum_eta_tt = 0.0;
        for i in 0..self.model.components.len() {
            let mi = self.model.components[i].m;
            sum_eta_t += self.mole_fractions[i] * mi * 3.0 * self.d[i].powi(2) * self.d_t[i];
            sum_eta_tt += self.mole_fractions[i]
                * mi
                * (6.0 * self.d[i] * self.d_t[i].powi(2) + 3.0 * self.d[i].powi(2) * self.d_tt[i]);
        }
        let pref = (std::f64::consts::PI * N_AVOGADRO / 6.0) * rho;
        let eta_t = pref * sum_eta_t;
        let eta_tt = pref * sum_eta_tt;

        let alpha_hs_eta = 2.0 * (2.0 - eta) / (1.0 - eta).powi(3);
        let alpha_hs_etaeta = 2.0 * (5.0 - 2.0 * eta) / (1.0 - eta).powi(4);
        let g = (1.0 - 0.5 * eta) / (1.0 - eta).powi(3);
        let g_eta = (5.0 - 2.0 * eta) / (2.0 * (1.0 - eta).powi(4));
        let g_etaeta = 3.0 * (3.0 - eta) / (1.0 - eta).powi(5);

        let mut sum_m_m1 = 0.0;
        for i in 0..self.model.components.len() {
            sum_m_m1 += self.mole_fractions[i] * (self.model.components[i].m - 1.0);
        }

        let ahc_eta = self.m_bar * alpha_hs_eta - sum_m_m1 * g_eta / g;
        let ahc_etaeta =
            self.m_bar * alpha_hs_etaeta - sum_m_m1 * (g_etaeta / g - (g_eta / g).powi(2));

        ahc_eta * eta_tt + ahc_etaeta * eta_t.powi(2)
    }

    fn dispersion_tt_derivative(&self, rho: f64, _eta: f64) -> f64 {
        // High-precision finite difference for 2nd temperature derivative of dispersion
        let h = 1e-4 * self.temperature;
        let prep_p =
            PreparedPcSaft::new(self.model, self.temperature + h, &self.mole_fractions).unwrap();
        let prep_m =
            PreparedPcSaft::new(self.model, self.temperature - h, &self.mole_fractions).unwrap();
        let a_p = prep_p.alpha_dispersion(rho, prep_p.packing_fraction(rho));
        let a_0 = self.alpha_dispersion(rho, self.packing_fraction(rho));
        let a_m = prep_m.alpha_dispersion(rho, prep_m.packing_fraction(rho));
        (a_p - 2.0 * a_0 + a_m) / (h * h)
    }

    fn helmholtz_fn_derivatives(&self, v: f64, n: &[f64]) -> Vec<f64> {
        let nc = self.model.components.len();
        let sumn: f64 = n.iter().sum();
        let rho = sumn / v;
        let a_res = self.residual_helmholtz_energy(rho);
        let a_rho = self.d_alpha_hard_chain_drho(rho, self.packing_fraction(rho))
            + self.d_alpha_dispersion_drho(rho, self.packing_fraction(rho));

        let mut f_n = Vec::with_capacity(nc);
        for i in 0..nc {
            // F_{n_i} = a_res + rho * a_rho + sumn * d(a_res)/d(n_i)
            let h = 1e-6 * sumn.max(1.0);
            let mut n_p = n.to_vec();
            let mut n_m = n.to_vec();
            n_p[i] += h;
            n_m[i] -= h;

            let prep_p = PreparedPcSaft::new(self.model, self.temperature, &n_p).unwrap();
            let prep_m = PreparedPcSaft::new(self.model, self.temperature, &n_m).unwrap();

            let da_dni = (prep_p.residual_helmholtz_energy(rho)
                - prep_m.residual_helmholtz_energy(rho))
                / (2.0 * h);

            f_n.push(a_res + rho * a_rho + sumn * da_dni);
        }
        f_n
    }
}

/// Pressure evaluation output package with thermodynamic differentials.
#[derive(Debug, Clone, Copy)]
pub struct PcSaftPressureJet {
    pub pressure: f64,
    pub compressibility: f64,
    pub residual_helmholtz: f64,
    pub dp_drho: f64,
    pub dp_dv: f64,
    pub dp_dt: f64,
}

/// Thermodynamic state property bundle.
#[derive(Debug, Clone)]
pub struct PcSaftState {
    pub pressure: f64,
    pub compressibility: f64,
    pub residual_helmholtz: f64,
    pub u_residual: f64,
    pub h_residual: f64,
    pub s_residual: f64,
    pub a_residual: f64,
    pub g_residual: f64,
    pub h_total: f64,
    pub s_total: f64,
    pub cv_residual: f64,
    pub cv_total: f64,
    pub cp_total: f64,
    pub dp_drho: f64,
    pub dp_dv: f64,
    pub dp_dt: f64,
    pub sound_speed: f64,
    pub ln_phi: Vec<f64>,
}

/// Pure component vapor-liquid coexistence saturation properties.
#[derive(Debug, Clone, Copy)]
pub struct PcSaftSaturation {
    pub pressure: f64,
    pub temperature: f64,
    pub liquid_volume: f64,
    pub vapor_volume: f64,
    pub enthalpy_vaporization: f64,
    pub entropy_vaporization: f64,
}
