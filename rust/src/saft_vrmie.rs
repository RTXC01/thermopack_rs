// SAFT-VR Mie and SAFT-VRQ Mie Equation of State.
// Reference: Lafitte et al. (2013) doi:10.1063/1.4819786 and Aasen et al. (2019) doi:10.1063/1.5136079.
// Copyright thermotools and Contributors. Apache-2.0.

mod association;
mod kernel;
use self::association::MieRdf;
use crate::association::{AssociationParams, DeltaCombRule, assoc_helmholtz_value};
use crate::association::{
    assoc_helmholtz_all, assoc_helmholtz_first, assoc_helmholtz_pressure_jet,
};
/// Extensive reduced Helmholtz energy (mol), all TVN derivatives and site fractions.
pub type SaftVrMieHelmholtz = crate::association::AssociationHelmholtz;
use crate::component::component;
use crate::cubic::Phase;
use crate::error::{Error, Result};
use crate::ideal::IdealGas;
use crate::saft_vrmie_data;

// Universal physical constants matching ThermoPack (CODATA 2018 / 2019 SI)
const N_A: f64 = 6.02214076e23;
const K_B: f64 = 1.380649e-23;
const RGAS: f64 = N_A * K_B; // 8.31446261815324 J/(mol*K)
const H_PLANCK: f64 = 6.62607015e-34; // J*s
const PI: f64 = std::f64::consts::PI;
const MACHINE_PREC: f64 = 2.220446049250313e-16;

// Gauss-Kronrod 21 weights and positive nodes for Barker-Henderson diameter
#[allow(clippy::excessive_precision)]
const XK_21: [f64; 11] = [
    0.995657163025808080735527280689003,
    0.973906528517171720077964012084452,
    0.930157491355708226001207180059508,
    0.865063366688984510732096688423493,
    0.780817726586416897063717578345042,
    0.679409568299024406234327365114874,
    0.562757134668604683339000099272694,
    0.433395394129247190799265943165784,
    0.294392862701460198131126603103866,
    0.148874338981631210884826001129720,
    0.0,
];

#[allow(clippy::excessive_precision)]
const WK_21: [f64; 11] = [
    0.011694638867371874278064396062192,
    0.032558162307964727478818972459390,
    0.054755896574351996031381300244580,
    0.075039674810919952767043140916190,
    0.093125454583697605535065465083366,
    0.109387158802297641899210590325805,
    0.123491976262065851077958109831074,
    0.134709217311473325928054001771707,
    0.142775938577060080797094273138717,
    0.147739104901338491374841515972068,
    0.149445554002916905664936468389821,
];

// Effective packing fraction polynomial coefficients for Sutherland potential (row: 0..3, col: 0..3)
const LAM_COEFF: [[f64; 4]; 4] = [
    [0.81096, 1.0205, -1.9057, 1.0885],
    [1.7888, -19.341, 22.845, -6.1962],
    [-37.578, 151.26, -228.14, 106.98],
    [92.284, -463.50, 973.92, -677.64],
];

// Phi matrix for Lafitte alpha functions f_1 .. f_6 (row: 0..6, col: 0..5)
const PHI_MAT: [[f64; 6]; 7] = [
    [7.5365557, -359.44, 1550.9, -1.19932, -1911.28, 9236.9],
    [-37.60463, 1825.6, -5070.1, 9.063632, 21390.175, -129430.0],
    [71.745953, -3168.0, 6534.6, -17.9482, -51320.7, 357230.0],
    [-46.83552, 1884.2, -3288.7, 11.34027, 37064.54, -315530.0],
    [-2.467982, -0.82376, -2.7171, 20.52142, 1103.742, 1390.2],
    [-0.50272, -3.1935, 2.0883, -56.6377, -3264.61, -4518.2],
    [8.0956883, 3.709, 0.0, 40.53683, 2556.181, 4241.6],
];

/// Component parameters for SAFT-VR Mie and SAFT-VRQ Mie.
#[derive(Debug, Clone, Copy)]
pub struct SaftVrMieComponent {
    pub ident: &'static str,
    pub m: f64,
    pub sigma: f64,
    pub eps_depth_divk: f64,
    pub lambda_a: f64,
    pub lambda_r: f64,
    pub mass: f64,
    pub fh_order: i32,
    pub reference: &'static str,
    pub molar_mass: f64,
}

/// Thermodynamic Helmholtz energy terms for SAFT-VR Mie.
#[derive(Debug, Clone, Copy)]
pub struct SaftVrMieHelmholtzTerms {
    pub a_hs: f64,
    pub a1: f64,
    pub a2: f64,
    pub a3: f64,
    pub a_chain: f64,
    pub a_res: f64,
}

/// Pressure jet bundle with first volume, density, and temperature derivatives.
#[derive(Debug, Clone, Copy)]
pub struct SaftVrMiePressureJet {
    pub pressure: f64,
    pub compressibility: f64,
    pub residual_helmholtz: f64,
    pub dp_drho: f64,
    pub dp_dv: f64,
    pub dp_dt: f64,
}

/// Complete thermodynamic state properties bundle at (T, rho).
#[derive(Debug, Clone)]
pub struct SaftVrMieState {
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

/// Pure-fluid vapor-liquid saturation equilibrium result.
#[derive(Debug, Clone, Copy)]
pub struct SaftVrMieSaturation {
    pub pressure: f64,
    pub temperature: f64,
    pub liquid_volume: f64,
    pub vapor_volume: f64,
    pub enthalpy_vaporization: f64,
    pub entropy_vaporization: f64,
}

/// SAFT-VR Mie and SAFT-VRQ Mie model instance.
#[derive(Debug, Clone)]
pub struct SaftVrMie {
    components: Vec<SaftVrMieComponent>,
    association: AssociationParams,
    integral_coefficients: Vec<[f64; 66]>,
    kij: Vec<f64>,
    lij: Vec<f64>,
    gamma_ij: Vec<f64>,
    lambda_a_ij: Vec<f64>,
    lambda_r_ij: Vec<f64>,
    sigma_ij: Vec<f64>,
    eps_divk_ij: Vec<f64>,
    c_ij: Vec<f64>,
    alpha_ij: Vec<f64>,
    f_alpha_ij: Vec<[f64; 6]>,
    d_feynhibbs_param_ij: Vec<f64>,
    ideal: Vec<IdealGas>,
    exact_binary_dhs: bool,
}

impl SaftVrMie {
    /// Create a new SAFT-VR Mie model instance with default parameters.
    pub fn new(comps: &[&str]) -> Result<Self> {
        Self::new_with_options(comps, "Default")
    }

    /// Create a new SAFT-VR Mie model instance with a specific parameter reference.
    pub fn new_with_options(comps: &[&str], parameter_reference: &str) -> Result<Self> {
        if comps.is_empty() {
            return Err(Error::InvalidInput("components list must not be empty"));
        }

        let nc = comps.len();
        let mut components = Vec::with_capacity(nc);
        let mut ideal = Vec::with_capacity(nc);

        let mut schemes = Vec::with_capacity(nc);
        let mut eps_assoc = Vec::with_capacity(nc);
        let mut kab = Vec::with_capacity(nc);
        for &ident in comps {
            let record = saft_vrmie_data::lookup_pure(ident, Some(parameter_reference))
                .ok_or_else(|| {
                    Error::UnknownComponent(format!("SAFT-VR Mie record for {ident}"))
                })?;
            schemes.push(record.scheme);
            eps_assoc.push(record.eps_assoc_divk * crate::R);
            kab.push(record.kab);
            let comp_info = component(ident)?;
            let molar_mass = comp_info.molar_mass; // g/mol

            components.push(SaftVrMieComponent {
                ident: record.ident,
                m: record.m,
                sigma: record.sigma,
                eps_depth_divk: record.eps_depth_divk,
                lambda_a: record.lambda_a,
                lambda_r: record.lambda_r,
                mass: record.mass,
                fh_order: record.fh_order,
                reference: record.reference,
                molar_mass,
            });

            ideal.push(IdealGas::new(comp_info)?);
        }

        let association = AssociationParams::for_saft_vr_mie(
            &schemes,
            &eps_assoc,
            &kab,
            DeltaCombRule::Standard,
        )?;
        let integral_coefficients = if association.has_sites() {
            vec![[0.0; 66]; nc * nc]
        } else {
            Vec::new()
        };
        let mut s = Self {
            association,
            integral_coefficients,
            components,
            kij: vec![0.0; nc * nc],
            lij: vec![0.0; nc * nc],
            gamma_ij: vec![0.0; nc * nc],
            lambda_a_ij: vec![0.0; nc * nc],
            lambda_r_ij: vec![0.0; nc * nc],
            sigma_ij: vec![0.0; nc * nc],
            eps_divk_ij: vec![0.0; nc * nc],
            c_ij: vec![0.0; nc * nc],
            alpha_ij: vec![0.0; nc * nc],
            f_alpha_ij: vec![[0.0; 6]; nc * nc],
            d_feynhibbs_param_ij: vec![0.0; nc * nc],
            ideal,
            exact_binary_dhs: false,
        };

        // Initialize interaction parameters from database
        for i in 0..nc {
            for j in 0..nc {
                let idx = i * nc + j;
                if i != j {
                    s.kij[idx] = saft_vrmie_data::lookup_kij(
                        s.components[i].ident,
                        s.components[j].ident,
                        Some(parameter_reference),
                    );
                    s.lij[idx] = saft_vrmie_data::lookup_lij(
                        s.components[i].ident,
                        s.components[j].ident,
                        Some(parameter_reference),
                    );
                }
            }
        }

        s.recompute_cross_parameters();
        Ok(s)
    }

    /// Extensive residual F=Ares/(RT) and its ten TVN blocks at raw mole amounts.
    pub fn residual_helmholtz_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<SaftVrMieHelmholtz> {
        crate::error::positive(volume, "total volume")?;
        let state = self.prepare(temperature, moles)?;
        let total: f64 = moles.iter().sum();
        let rho = total / volume;
        let mut out = assoc_helmholtz_all(
            &self.association,
            &MieRdf::new(&state, rho),
            temperature,
            volume,
            moles,
        )?;
        for i in 0..moles.len() + 2 {
            for j in i..moles.len() + 2 {
                let f = state.kernel(rho, Some(i), Some(j))?.total;
                if i == j {
                    match i {
                        0 => {
                            out.f += total * f.f0;
                            out.f_t += total * f.f1;
                        }
                        1 => out.f_v += f.f1,
                        _ => out.f_n[i - 2] += f.f1,
                    }
                }
                match (i, j) {
                    (0, 0) => out.f_tt += total * f.f12,
                    (0, 1) => out.f_tv += f.f12,
                    (0, _) => out.f_tn[j - 2] += f.f12,
                    (1, 1) => out.f_vv += f.f12 / total,
                    (1, _) => out.f_vn[j - 2] += f.f12 / total,
                    _ => {
                        out.f_nn[(i - 2) * moles.len() + j - 2] += f.f12 / total;
                        if i != j {
                            out.f_nn[(j - 2) * moles.len() + i - 2] += f.f12 / total;
                        }
                    }
                }
            }
        }
        Ok(out)
    }

    /// Number of components in the mixture.
    #[inline]
    pub fn num_components(&self) -> usize {
        self.components.len()
    }

    /// Access component parameters.
    pub fn components(&self) -> &[SaftVrMieComponent] {
        &self.components
    }

    /// Set binary interaction parameter kij.
    pub fn set_kij(&mut self, i: usize, j: usize, val: f64) {
        let nc = self.num_components();
        self.kij[i * nc + j] = val;
        self.kij[j * nc + i] = val;
        self.recompute_cross_parameters();
    }

    /// Get binary interaction parameter kij.
    pub fn get_kij(&self, i: usize, j: usize) -> f64 {
        let nc = self.num_components();
        self.kij[i * nc + j]
    }

    /// Set binary interaction parameter lij.
    pub fn set_lij(&mut self, i: usize, j: usize, val: f64) {
        let nc = self.num_components();
        self.lij[i * nc + j] = val;
        self.lij[j * nc + i] = val;
        self.recompute_cross_parameters();
    }

    /// Get binary interaction parameter lij.
    pub fn get_lij(&self, i: usize, j: usize) -> f64 {
        let nc = self.num_components();
        self.lij[i * nc + j]
    }

    /// Recompute all cross-interaction parameters.
    fn recompute_cross_parameters(&mut self) {
        let nc = self.num_components();
        for i in 0..nc {
            for j in 0..nc {
                let idx = i * nc + j;
                let c_i = &self.components[i];
                let c_j = &self.components[j];

                // Attractive and repulsive exponents
                let la_ij = 3.0 + ((c_i.lambda_a - 3.0) * (c_j.lambda_a - 3.0)).sqrt();
                let lr_ij = 3.0
                    + (1.0 - self.gamma_ij[idx])
                        * ((c_i.lambda_r - 3.0) * (c_j.lambda_r - 3.0)).sqrt();
                self.lambda_a_ij[idx] = la_ij;
                self.lambda_r_ij[idx] = lr_ij;
                if self.association.has_sites() {
                    self.integral_coefficients[idx] = saft_vrmie_data::ASSOCIATION_B
                        .map(|row| row.into_iter().rev().fold(0.0, |acc, b| acc * lr_ij + b));
                }

                // Segment diameter
                let sig_ij = 0.5 * (c_i.sigma + c_j.sigma) * (1.0 - self.lij[idx]);
                self.sigma_ij[idx] = sig_ij;

                // Dispersion well depth (Lafitte combining rule)
                let si3 = c_i.sigma * c_i.sigma * c_i.sigma;
                let sj3 = c_j.sigma * c_j.sigma * c_j.sigma;
                let eps_ij = (1.0 - self.kij[idx])
                    * (si3 * sj3).sqrt()
                    * (c_i.eps_depth_divk * c_j.eps_depth_divk).sqrt()
                    / (sig_ij * sig_ij * sig_ij);
                self.eps_divk_ij[idx] = eps_ij;

                // Mie prefactor C_ij
                let exp_factor = la_ij / (lr_ij - la_ij);
                let c_fac = (lr_ij / (lr_ij - la_ij)) * (lr_ij / la_ij).powf(exp_factor);
                self.c_ij[idx] = c_fac;

                // van der Waals attractive constant alpha_ij
                let alpha = c_fac * (1.0 / (la_ij - 3.0) - 1.0 / (lr_ij - 3.0));
                self.alpha_ij[idx] = alpha;

                // f_1 .. f_6 polynomials of alpha
                let f_alpha = alpha_coefficients(alpha);
                self.f_alpha_ij[idx] = f_alpha;

                // Feynman-Hibbs quantum parameter
                let mass_product = c_i.mass * c_j.mass;
                if mass_product > 0.0 {
                    let reduced_mass = mass_product / (c_i.mass + c_j.mass);
                    self.d_feynhibbs_param_ij[idx] =
                        H_PLANCK * H_PLANCK / (96.0 * reduced_mass * K_B * PI * PI);
                } else {
                    self.d_feynhibbs_param_ij[idx] = 0.0;
                }
            }
        }
    }

    /// Whether cross hard-sphere diameters use exact Gauss-Kronrod quadrature on cross Mie potential.
    pub fn exact_binary_dhs(&self) -> bool {
        self.exact_binary_dhs
    }

    /// Set whether cross hard-sphere diameters use exact Gauss-Kronrod quadrature.
    pub fn set_exact_binary_dhs(&mut self, val: bool) {
        self.exact_binary_dhs = val;
    }

    /// Prepare state for evaluation at temperature T and mole fractions x.
    pub fn prepare(
        &self,
        temperature: f64,
        mole_fractions: &[f64],
    ) -> Result<PreparedSaftVrMie<'_>> {
        if temperature <= 0.0 || !temperature.is_finite() {
            return Err(Error::InvalidInput("temperature must be strictly positive"));
        }
        let nc = self.num_components();
        if mole_fractions.len() != nc {
            return Err(Error::InvalidInput("mole fractions length mismatch"));
        }
        if mole_fractions.iter().any(|&n| n < 0.0) {
            return Err(Error::InvalidInput("nonnegative mole amounts required"));
        }
        let sum_x: f64 = mole_fractions.iter().sum();
        if sum_x <= 0.0 || !sum_x.is_finite() {
            return Err(Error::InvalidInput(
                "mole fractions sum must be strictly positive",
            ));
        }
        let x: Vec<f64> = mole_fractions.iter().map(|&v| v / sum_x).collect();

        PreparedSaftVrMie::new(self, temperature, &x)
    }
}

/// State prepared at specific temperature T and mole fractions x.
#[derive(Debug, Clone)]
pub struct PreparedSaftVrMie<'a> {
    model: &'a SaftVrMie,
    temperature: f64,
    mole_fractions: Vec<f64>,
    d: Vec<f64>,    // nc * nc
    d_t: Vec<f64>,  // nc * nc
    d_tt: Vec<f64>, // nc * nc
    m_bar: f64,
    molar_mass_mix: f64,
    conv_eta: f64,  // eta = conv_eta * rho
    conv_zeta: f64, // zeta_bar = conv_zeta * rho
    xs: Vec<f64>,   // segment mole fractions xs_i = x_i * m_i / m_bar
}

impl<'a> PreparedSaftVrMie<'a> {
    fn new(model: &'a SaftVrMie, temperature: f64, x: &[f64]) -> Result<Self> {
        let nc = model.num_components();
        let mut d = vec![0.0; nc * nc];
        let mut d_t = vec![0.0; nc * nc];
        let mut d_tt = vec![0.0; nc * nc];

        // 1. Solve pure Barker-Henderson hard-sphere diameters
        for i in 0..nc {
            let idx = i * nc + i;
            let comp = &model.components[i];
            let c_fac = model.c_ij[idx];
            let d_param = model.d_feynhibbs_param_ij[idx];
            let (di, di_t, di_tt) = calc_hardsphere_diameter_pure(
                temperature,
                comp.sigma,
                comp.eps_depth_divk,
                comp.lambda_a,
                comp.lambda_r,
                c_fac,
                d_param,
                comp.fh_order,
            )?;
            d[idx] = di;
            d_t[idx] = di_t;
            d_tt[idx] = di_tt;
        }

        // 2. Mixture diameters
        for i in 0..nc {
            for j in 0..nc {
                if i != j {
                    let idx = i * nc + j;
                    let ii = i * nc + i;
                    let jj = j * nc + j;
                    d[idx] = 0.5 * (d[ii] + d[jj]);
                    d_t[idx] = 0.5 * (d_t[ii] + d_t[jj]);
                    d_tt[idx] = 0.5 * (d_tt[ii] + d_tt[jj]);
                }
            }
        }

        // 3. Mixture mean segment number and molar mass
        let mut m_bar = 0.0;
        let mut molar_mass_mix = 0.0;
        for (xi, comp) in x.iter().zip(&model.components) {
            m_bar += xi * comp.m;
            molar_mass_mix += xi * comp.molar_mass;
        }

        let xs: Vec<f64> = x
            .iter()
            .zip(&model.components)
            .map(|(xi, comp)| xi * comp.m / m_bar)
            .collect();

        // 4. Precompute packing fraction factors
        let mut sum_d3 = 0.0;
        let mut sum_sig3 = 0.0;
        for i in 0..nc {
            for j in 0..nc {
                let idx = i * nc + j;
                let dij = d[idx];
                let sij = model.sigma_ij[idx];
                sum_d3 += xs[i] * xs[j] * (dij * dij * dij);
                sum_sig3 += xs[i] * xs[j] * (sij * sij * sij);
            }
        }

        let conv_eta = (PI * N_A / 6.0) * m_bar * sum_d3;
        let conv_zeta = (PI * N_A / 6.0) * m_bar * sum_sig3;

        Ok(Self {
            model,
            temperature,
            mole_fractions: x.to_vec(),
            d,
            d_t,
            d_tt,
            m_bar,
            molar_mass_mix,
            conv_eta,
            conv_zeta,
            xs,
        })
    }

    /// Hard-sphere packing fraction eta = (pi * N_A / 6) * rho * m_bar * sum(xs_i xs_j d_ij^3).
    #[inline]
    pub fn packing_fraction(&self, rho: f64) -> f64 {
        self.conv_eta * rho
    }

    /// Dispersion packing fraction zeta_bar = (pi * N_A / 6) * rho * m_bar * sum(xs_i xs_j sigma_ij^3).
    #[inline]
    pub fn dispersion_packing_fraction(&self, rho: f64) -> f64 {
        self.conv_zeta * rho
    }

    /// Hard-sphere diameters (m).
    pub fn d(&self) -> &[f64] {
        &self.d
    }

    /// Temperature derivatives of hard-sphere diameters [m/K].
    pub fn d_t(&self) -> &[f64] {
        &self.d_t
    }

    /// Second temperature derivatives of hard-sphere diameters [m/K^2].
    pub fn d_tt(&self) -> &[f64] {
        &self.d_tt
    }

    /// Helmholtz free energy terms at density rho [mol/m^3].
    pub fn helmholtz_terms(&self, rho: f64) -> Result<SaftVrMieHelmholtzTerms> {
        if rho == 0.0 {
            return Ok(SaftVrMieHelmholtzTerms {
                a_hs: 0.,
                a1: 0.,
                a2: 0.,
                a3: 0.,
                a_chain: 0.,
                a_res: 0.,
            });
        }
        let terms = self.kernel(rho, None, None)?;
        let t = self.temperature;
        let a_assoc = if self.model.association.has_sites() && rho > 0.0 {
            assoc_helmholtz_value(
                &self.model.association,
                &MieRdf::new(self, rho),
                t,
                1.0 / rho,
                &self.mole_fractions,
            )?
            .0
        } else {
            0.0
        };
        Ok(SaftVrMieHelmholtzTerms {
            a_hs: terms.hs.f0,
            a1: terms.a1.f0,
            a2: terms.a2.f0,
            a3: terms.a3.f0,
            a_chain: terms.chain.f0,
            a_res: terms.total.f0 + a_assoc,
        })
    }

    /// Pressure and derivatives at density rho [mol/m^3].
    pub fn pressure_and_derivatives(&self, rho: f64) -> Result<SaftVrMiePressureJet> {
        let v = 1.0 / rho;
        let t = self.temperature;

        let vv = self.kernel(rho, Some(1), Some(1))?.total;
        let tv = self.kernel(rho, Some(0), Some(1))?.total;
        let mut f = vv.f0;
        let mut f_v = vv.f1;
        let mut f_vv = vv.f12;
        let mut f_tv = tv.f12;
        if self.model.association.has_sites() {
            let a = assoc_helmholtz_pressure_jet(
                &self.model.association,
                &MieRdf::new(self, rho),
                t,
                v,
                &self.mole_fractions,
            )?;
            f += a.f;
            f_v += a.f_v;
            f_vv += a.f_vv;
            f_tv += a.f_tv;
        }
        let p = RGAS * t * (rho - f_v);
        let z = 1. - f_v / rho;
        let dp_dv = -RGAS * t * (rho * rho + f_vv);
        let dp_drho = -dp_dv / (rho * rho);

        // dP/dT = rho * R - R * F_V - R * T * F_TV
        let dp_dt = rho * RGAS - RGAS * f_v - RGAS * t * f_tv;

        Ok(SaftVrMiePressureJet {
            pressure: p,
            compressibility: z,
            residual_helmholtz: f,
            dp_drho,
            dp_dv,
            dp_dt,
        })
    }

    /// Complete thermodynamic state properties at density rho [mol/m^3].
    pub fn state(&self, rho: f64) -> Result<SaftVrMieState> {
        let jet = self.pressure_and_derivatives(rho)?;
        let t = self.temperature;
        let tt = self.kernel(rho, Some(0), Some(0))?.total;
        let mut f_t = tt.f1;
        let mut f_tt = tt.f12;
        if self.model.association.has_sites() {
            let a = assoc_helmholtz_pressure_jet(
                &self.model.association,
                &MieRdf::new(self, rho),
                t,
                1. / rho,
                &self.mole_fractions,
            )?;
            f_t += a.f_t;
            f_tt += a.f_tt;
        }
        let u_res = -RGAS * t * t * f_t;
        let h_res = u_res + RGAS * t * (jet.compressibility - 1.0);
        let s_res = -RGAS * (t * f_t + jet.residual_helmholtz) + RGAS * (jet.compressibility).ln();
        let a_res_energy = RGAS * t * jet.residual_helmholtz;
        let g_res = h_res - t * s_res;

        let cv_res = -RGAS * (2.0 * t * f_t + t * t * f_tt);

        // Ideal gas contribution
        let mut cp_id = 0.0;
        let mut h_id = 0.0;
        let mut s_id = 0.0;
        let p_ref = 100_000.0;
        for i in 0..self.model.components.len() {
            let id_props = self.model.ideal[i].properties(t, p_ref)?;
            let xi = self.mole_fractions[i];
            cp_id += xi * id_props.cp;
            h_id += xi * id_props.enthalpy;
            let mix_term = if xi > 1e-15 { xi * xi.ln() } else { 0.0 };
            s_id += xi * id_props.entropy - RGAS * mix_term;
        }
        s_id -= RGAS * (rho * RGAS * t / p_ref).ln();

        let cv_total = cp_id - RGAS + cv_res;
        let cp_total = cv_total - t * (jet.dp_dt * jet.dp_dt) / jet.dp_dv;

        let gamma = cp_total / cv_total;
        let mw_kg = self.molar_mass_mix * 1e-3;
        let sound_speed = if jet.dp_drho > 0.0 && gamma > 0.0 && mw_kg > 0.0 {
            (gamma * jet.dp_drho / mw_kg).sqrt()
        } else {
            0.0
        };

        // Fugacity coefficients
        let ln_phi = self.fugacity_coefficients(rho, jet.compressibility)?;

        Ok(SaftVrMieState {
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

    /// Fugacity coefficients ln(phi_i) at density rho [mol/m^3].
    pub fn fugacity_coefficients(&self, rho: f64, z: f64) -> Result<Vec<f64>> {
        let nc = self.model.components.len();
        if nc == 1 {
            let terms = self.helmholtz_terms(rho)?;
            let ln_phi = terms.a_res + (z - 1.0) - z.ln();
            return Ok(vec![ln_phi]);
        }

        let mut ln_phi = (0..nc)
            .map(|i| Ok(self.kernel(rho, Some(2 + i), None)?.total.f1 - z.ln()))
            .collect::<Result<Vec<_>>>()?;
        if self.model.association.has_sites() {
            let a = assoc_helmholtz_first(
                &self.model.association,
                &MieRdf::new(self, rho),
                self.temperature,
                1. / rho,
                &self.mole_fractions,
            )?;
            for (value, assoc) in ln_phi.iter_mut().zip(a.f_n) {
                *value += assoc;
            }
        }
        Ok(ln_phi)
    }

    /// Single-phase volume root for a specified pressure and phase.
    pub fn volume(&self, pressure_target: f64, phase: Phase) -> Result<f64> {
        if pressure_target <= 0.0 || !pressure_target.is_finite() {
            return Err(Error::InvalidInput(
                "pressure target must be strictly positive",
            ));
        }

        let conv_eta = self.conv_eta;
        let mut zeta_max = 0.55;
        let mut zeta_min = 1e-8;

        let mut zeta = match phase {
            Phase::Vapor => {
                let z_ideal = conv_eta * pressure_target / (RGAS * self.temperature);
                z_ideal.clamp(zeta_min, 0.05)
            }
            Phase::Liquid => 0.45,
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

        let mut v = conv_eta / zeta;
        let mut jet = self.pressure_and_derivatives(1.0 / v)?;
        let mut f = (1.0 - zeta) * (jet.pressure - pressure_target);

        let max_iter = 60;
        for _ in 0..max_iter {
            let df_dzeta = pressure_target - jet.pressure - ((1.0 - zeta) / zeta) * jet.dp_dv * v;
            let mut zeta_new = zeta - f / df_dzeta;

            if zeta_new > zeta_max || zeta_new < zeta_min || !zeta_new.is_finite() {
                zeta_new = 0.5 * (zeta_max + zeta_min);
            }

            let v_old = v;
            zeta = zeta_new;
            v = conv_eta / zeta;

            jet = self.pressure_and_derivatives(1.0 / v)?;
            f = (1.0 - zeta) * (jet.pressure - pressure_target);

            let v_converged = (v - v_old).abs() < 1e5 * MACHINE_PREC * v_old;
            let p_converged =
                (jet.pressure - pressure_target).abs() < 1e8 * MACHINE_PREC * pressure_target;

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

    /// Pure-fluid saturation pressure at temperature T.
    pub fn saturation_pressure(&self, temperature: f64) -> Result<SaftVrMieSaturation> {
        if self.model.num_components() != 1 {
            return Err(Error::Unsupported(
                "saturation_pressure is only supported for pure fluids",
            ));
        }

        let prep = self.model.prepare(temperature, &[1.0])?;
        let comp_info = component(self.model.components[0].ident)?;
        let tc = comp_info.critical_temperature;
        let pc = comp_info.critical_pressure;
        let omega = comp_info.acentric_factor;

        let exp_arg: f64 = 5.373 * (1.0 + omega) * (1.0 - tc / temperature);
        let mut p_est = pc * exp_arg.exp();
        if p_est >= pc || p_est <= 0.0 || !p_est.is_finite() {
            p_est = 0.5 * pc;
        }

        let mut p = p_est;
        let max_iter = 50;

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
                return Ok(SaftVrMieSaturation {
                    pressure: p,
                    temperature,
                    liquid_volume: vl,
                    vapor_volume: vv,
                    enthalpy_vaporization: state_v.h_total - state_l.h_total,
                    entropy_vaporization: state_v.s_total - state_l.s_total,
                });
            }

            let d_delta_dp = (vl - vv) / (RGAS * temperature);
            let mut dp = -delta_ln_phi / d_delta_dp;
            if dp.abs() > 0.5 * p {
                dp = 0.5 * p * dp.signum();
            }

            p += dp;
            if p <= 0.0 {
                p = 0.5 * (p - dp);
            }
        }

        Err(Error::NotConverged(
            "saturation solver exceeded maximum iterations",
        ))
    }
}

// ----------------------------------------------------------------------------
// Internal numerical & potential helper routines
// ----------------------------------------------------------------------------

#[derive(Debug, Clone, Copy)]
struct MiePot {
    u_divk: f64,
    u_divk_t: f64,
    u_divk_tt: f64,
    u_divk_r: f64,
    u_divk_rr: f64,
    u_divk_tr: f64,
}

#[allow(clippy::too_many_arguments)]
fn calc_mie_potential_full(
    r: f64,
    temperature: f64,
    sigma: f64,
    eps_depth_divk: f64,
    lambda_a: f64,
    lambda_r: f64,
    c_fac: f64,
    d_param: f64,
    fh_order: i32,
) -> MiePot {
    let n = lambda_r;
    let m = lambda_a;
    let en_r = 1.0 / r;
    let r1 = sigma * en_r;

    let r_n = r1.powf(n);
    let r_m = r1.powf(m);
    let dr_n = -n * r_n * en_r;
    let dr_m = -m * r_m * en_r;
    let d2r_n = -(n + 1.0) * dr_n * en_r;
    let d2r_m = -(m + 1.0) * dr_m * en_r;

    let mut u = r_n - r_m;
    let mut u_t = 0.0;
    let mut u_tt = 0.0;
    let mut u_r = dr_n - dr_m;
    let mut u_rr = d2r_n - d2r_m;
    let mut u_tr = 0.0;

    if fh_order > 0 && d_param > 0.0 && temperature > 0.0 {
        let tinv = 1.0 / temperature;
        let d = d_param * tinv;
        let d_t = -d_param * tinv * tinv;
        let d_tt = 2.0 * d_param * (tinv * tinv * tinv);

        let q1_n = n * (n - 1.0);
        let q1_m = m * (m - 1.0);
        let en_r2 = en_r * en_r;

        let prod_q1 = en_r2 * (q1_n * r_n - q1_m * r_m);
        let prod_q1_r = -en_r2 * ((n + 2.0) * q1_n * r_n - (m + 2.0) * q1_m * r_m) * en_r;
        let prod_q1_rr = en_r2
            * ((n + 3.0) * (n + 2.0) * q1_n * r_n - (m + 3.0) * (m + 2.0) * q1_m * r_m)
            * en_r2;

        u += d * prod_q1;
        u_t += d_t * prod_q1;
        u_tt += d_tt * prod_q1;
        u_r += d * prod_q1_r;
        u_rr += d * prod_q1_rr;
        u_tr += d_t * prod_q1_r;

        if fh_order > 1 {
            let d2 = d * d;
            let d2_t = 2.0 * d * d_t;
            let d2_tt = 2.0 * (d * d_tt + d_t * d_t);

            let q2_n = 0.5 * n * (n - 1.0) * (n + 1.0) * (n + 2.0);
            let q2_m = 0.5 * m * (m - 1.0) * (m + 1.0) * (m + 2.0);
            let en_r4 = en_r2 * en_r2;

            let prod_q2 = en_r4 * (q2_n * r_n - q2_m * r_m);
            let prod_q2_r = -en_r4 * ((n + 4.0) * q2_n * r_n - (m + 4.0) * q2_m * r_m) * en_r;
            let prod_q2_rr = en_r4
                * ((n + 5.0) * (n + 4.0) * q2_n * r_n - (m + 5.0) * (m + 4.0) * q2_m * r_m)
                * en_r2;

            u += d2 * prod_q2;
            u_t += d2_t * prod_q2;
            u_tt += d2_tt * prod_q2;
            u_r += d2 * prod_q2_r;
            u_rr += d2 * prod_q2_rr;
            u_tr += d2_t * prod_q2_r;
        }
    }

    let pref = c_fac * eps_depth_divk;
    MiePot {
        u_divk: u * pref,
        u_divk_t: u_t * pref,
        u_divk_tt: u_tt * pref,
        u_divk_r: u_r * pref,
        u_divk_rr: u_rr * pref,
        u_divk_tr: u_tr * pref,
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn calc_hardsphere_diameter_pure(
    temperature: f64,
    sigma: f64,
    eps_depth_divk: f64,
    lambda_a: f64,
    lambda_r: f64,
    c_fac: f64,
    d_param: f64,
    fh_order: i32,
) -> Result<(f64, f64, f64)> {
    // 1. Solve effective sigma (where U(r) = 0)
    let (sigma_eff, sig_t, sig_tt) = if fh_order == 0 {
        (sigma, 0.0, 0.0)
    } else {
        let mut r = sigma;
        for _ in 0..25 {
            let pot = calc_mie_potential_full(
                r,
                temperature,
                sigma,
                eps_depth_divk,
                lambda_a,
                lambda_r,
                c_fac,
                d_param,
                fh_order,
            );
            let f = pot.u_divk;
            let df = pot.u_divk_r;
            let d2f = pot.u_divk_rr;
            let step = -f / df * (1.0 + f * d2f / (2.0 * df * df));
            r += step;
            if step.abs() < 1e-15 * sigma {
                break;
            }
        }
        let pot = calc_mie_potential_full(
            r,
            temperature,
            sigma,
            eps_depth_divk,
            lambda_a,
            lambda_r,
            c_fac,
            d_param,
            fh_order,
        );
        let s_t = -pot.u_divk_t / pot.u_divk_r;
        let s_tt =
            -(pot.u_divk_tt + pot.u_divk_rr * s_t * s_t + 2.0 * pot.u_divk_tr * s_t) / pot.u_divk_r;
        (r, s_t, s_tt)
    };

    // 2. Solve r0 (where exp(-beta u(r0)) = machine_prec)
    let ln_eps = MACHINE_PREC.ln();
    let mut rs = 0.70;
    for _ in 0..25 {
        let r = rs * sigma_eff;
        let pot = calc_mie_potential_full(
            r,
            temperature,
            sigma,
            eps_depth_divk,
            lambda_a,
            lambda_r,
            c_fac,
            d_param,
            fh_order,
        );
        let f = pot.u_divk / temperature + ln_eps;
        let df = sigma_eff * pot.u_divk_r / temperature;
        let d2f = sigma_eff * sigma_eff * pot.u_divk_rr / temperature;
        let step = -f / df * (1.0 + f * d2f / (2.0 * df * df));
        rs += step;
        if step.abs() < 1e-15 {
            break;
        }
    }

    let r0 = rs * sigma_eff;
    let pot0 = calc_mie_potential_full(
        r0,
        temperature,
        sigma,
        eps_depth_divk,
        lambda_a,
        lambda_r,
        c_fac,
        d_param,
        fh_order,
    );
    let r0_t = -(pot0.u_divk_t + ln_eps) / pot0.u_divk_r;
    let r0_tt = -(pot0.u_divk_tt + 2.0 * pot0.u_divk_tr * r0_t + pot0.u_divk_rr * r0_t * r0_t)
        / pot0.u_divk_r;

    // 3. Kronrod 21 Quadrature integration on [r0, sigma_eff]
    let mut d_int = 0.0;
    let mut d_t_int = 0.0;
    let mut d_tt_int = 0.0;

    let sigma_j2 = 0.5 * (sigma_eff - r0);
    let sigma_j2_t = 0.5 * (sig_t - r0_t);
    let sigma_j2_tt = 0.5 * (sig_tt - r0_tt);

    // Symmetric 21 points
    let mut x_nodes = [0.0; 21];
    let mut w_nodes = [0.0; 21];
    for i in 0..10 {
        x_nodes[i] = -XK_21[i];
        x_nodes[11 + i] = XK_21[9 - i];
        w_nodes[i] = WK_21[i];
        w_nodes[11 + i] = WK_21[9 - i];
    }
    x_nodes[10] = 0.0;
    w_nodes[10] = WK_21[10];

    let t2 = temperature * temperature;
    let t3 = t2 * temperature;

    for i in 0..21 {
        let xk = x_nodes[i];
        let wk = w_nodes[i];

        let r = sigma_j2 * xk + sigma_j2 + r0;
        let r_t = sigma_j2_t * xk + sigma_j2_t + r0_t;
        let r_tt = sigma_j2_tt * xk + sigma_j2_tt + r0_tt;

        let pot = calc_mie_potential_full(
            r,
            temperature,
            sigma,
            eps_depth_divk,
            lambda_a,
            lambda_r,
            c_fac,
            d_param,
            fh_order,
        );

        let pref = sigma_j2 * wk;
        let pref_t = sigma_j2_t * wk;
        let pref_tt = sigma_j2_tt * wk;

        let minus_exp = -(-pot.u_divk / temperature).exp();
        d_int += pref * (1.0 + minus_exp);

        let dt3 = pot.u_divk / t2 - pot.u_divk_t / temperature - pot.u_divk_r * r_t / temperature;
        let d_t1 = pref;
        let d_t2 = minus_exp;
        let d_t4 = pref_t;
        let d_t5 = 1.0 + minus_exp;
        d_t_int += d_t1 * d_t2 * dt3 + d_t4 * d_t5;

        let d_tt1 = pref_t;
        let d_tt2 = minus_exp * dt3;
        let d_tt3 = pot.u_divk_t / t2 + pot.u_divk_r * r_t / t2
            - 2.0 * pot.u_divk / t3
            - pot.u_divk_tt / temperature
            - pot.u_divk_tr * r_t / temperature
            + pot.u_divk_t / t2
            - pot.u_divk_tr * r_t / temperature
            - pot.u_divk_rr * (r_t * r_t) / temperature
            - pot.u_divk_r * r_tt / temperature
            + pot.u_divk_r * r_t / t2;
        let d_tt4 = pref_tt;
        let d_tt5 = d_tt2;
        d_tt_int += d_tt1 * d_t2 * dt3
            + d_t1 * d_tt2 * dt3
            + d_t1 * d_t2 * d_tt3
            + d_tt4 * d_t5
            + d_t4 * d_tt5;
    }

    Ok((r0 + d_int, r0_t + d_t_int, r0_tt + d_tt_int))
}

/// The six coefficients of the polynomials of the van der Waals attractive
/// constant, the routine `calcFunAlpha` of the source.
pub(crate) fn alpha_coefficients(alpha: f64) -> [f64; 6] {
    let mut f_alpha = [0.0; 6];
    let alpha_pow = [1.0, alpha, alpha * alpha, alpha * alpha * alpha];
    for k in 0..6 {
        let mut num = 0.0;
        for p in 0..4 {
            num += PHI_MAT[p][k] * alpha_pow[p];
        }
        let mut den = 1.0;
        for p in 1..4 {
            den += PHI_MAT[p + 3][k] * alpha_pow[p];
        }
        f_alpha[k] = num / den;
    }
    f_alpha
}
