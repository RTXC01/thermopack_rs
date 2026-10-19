//! Native Rust implementation of CPA (Cubic-Plus-Association).
//!
//! Combining cubic equations of state (SRK, PR) with Wertheim first-order
//! thermodynamic perturbation theory for associating fluids.
//!
//! # References
//! - Kontogeorgis et al. (1996, 1999, 2006, 2010): The Cubic-Plus-Association (CPA) Equation of State.
//! - Wertheim (1984, 1986): Fluids with Highly Directional Attractive Forces.
//! - Michelsen & Hendriks (2001): Physical Properties from Association Models.

use crate::component::component;
use crate::cpa_data::{AssocScheme, lookup_cpa_binary, lookup_cpa_pure};
use crate::cubic::Phase;
use crate::error::{Error, Result};
use crate::ideal::IdealGas;

/// Gas constant R [J / (mol * K)] matching ThermoPack Fortran exact constant.
pub const RGAS: f64 = 8.31446261815324;

/// Underlying cubic equation of state type for CPA.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CpaCubicKind {
    /// CPA-SRK (Soave-Redlich-Kwong reference, u=1, w=0).
    Srk,
    /// CPA-PR (Peng-Robinson reference, u=2, w=-1).
    Pr,
}

/// Single associating site descriptor.
#[derive(Debug, Clone, Copy)]
pub struct AssocSite {
    /// Component index owning this site.
    pub comp_idx: usize,
    /// Local site index within the component (1-based: 1, 2, 3, 4).
    pub site_idx: usize,
    /// Polarity: -1 = negative/lone-pair, +1 = positive/proton, 0 = bipolar.
    pub polarity: i8,
}

/// Component model parameters for CPA.
#[derive(Debug, Clone)]
pub struct CpaComponent {
    /// Component name.
    pub name: String,
    /// Energy parameter a0 [Pa * m^6 / mol^2].
    pub a0: f64,
    /// Covolume parameter b [m^3 / mol].
    pub b: f64,
    /// Critical temperature Tc [K] (used in Soave alpha).
    pub tc: f64,
    /// Classic alpha correlation parameter c1.
    pub c1: f64,
    /// Association energy epsilon [J / mol].
    pub eps: f64,
    /// Association volume beta [-].
    pub beta: f64,
    /// Association scheme.
    pub scheme: AssocScheme,
    /// Whether simplified hard-sphere RDF is used.
    pub simplified_rdf: bool,
    /// Molar mass [g / mol].
    pub molar_mass: f64,
}

/// Thermodynamic jet for pressure and its derivatives.
#[derive(Debug, Clone, Copy)]
pub struct CpaPressureJet {
    /// Pressure P [Pa].
    pub pressure: f64,
    /// Compressibility factor Z = P / (rho * R * T).
    pub compressibility: f64,
    /// Dimensionless residual Helmholtz energy per mole a_res = A_res / (R * T).
    pub residual_helmholtz: f64,
    /// Derivative dP/drho [Pa * m^3 / mol].
    pub dp_drho: f64,
    /// Derivative dP/dV [Pa * mol / m^3].
    pub dp_dv: f64,
    /// Derivative dP/dT [Pa / K].
    pub dp_dt: f64,
}

/// Full thermodynamic state properties at specified (T, rho).
#[derive(Debug, Clone)]
pub struct CpaState {
    /// Pressure P [Pa].
    pub pressure: f64,
    /// Compressibility factor Z [-].
    pub compressibility: f64,
    /// Dimensionless residual Helmholtz energy per mole.
    pub residual_helmholtz: f64,
    /// Residual internal energy U_res [J / mol].
    pub u_residual: f64,
    /// Residual enthalpy H_res [J / mol].
    pub h_residual: f64,
    /// Residual entropy S_res [J / (mol * K)].
    pub s_residual: f64,
    /// Residual Helmholtz energy A_res [J / mol].
    pub a_residual: f64,
    /// Residual Gibbs energy G_res [J / mol].
    pub g_residual: f64,
    /// Total enthalpy H_total = H_ideal + H_res [J / mol].
    pub h_total: f64,
    /// Total entropy S_total = S_ideal + S_res [J / (mol * K)].
    pub s_total: f64,
    /// Residual isochoric heat capacity Cv_res [J / (mol * K)].
    pub cv_residual: f64,
    /// Total isochoric heat capacity Cv_total [J / (mol * K)].
    pub cv_total: f64,
    /// Total isobaric heat capacity Cp_total [J / (mol * K)].
    pub cp_total: f64,
    /// Speed of sound [m / s].
    pub sound_speed: f64,
    /// Natural logarithm of fugacity coefficients ln(phi_i).
    pub ln_phi: Vec<f64>,
}

/// Saturation state for pure fluids.
#[derive(Debug, Clone, Copy)]
pub struct CpaSaturation {
    /// Saturation temperature T [K].
    pub temperature: f64,
    /// Saturation pressure P_sat [Pa].
    pub pressure: f64,
    /// Liquid molar volume [m^3 / mol].
    pub liquid_volume: f64,
    /// Vapor molar volume [m^3 / mol].
    pub vapor_volume: f64,
}

/// CPA Equation of State model.
#[derive(Debug, Clone)]
pub struct Cpa {
    pub cubic_kind: CpaCubicKind,
    pub components: Vec<CpaComponent>,
    pub kij_a: Vec<f64>,
    pub kij_eps: Vec<f64>,
    pub kij_beta: Vec<f64>,
    pub eps_comb_arithmetic: Vec<bool>,
    pub beta_comb_arithmetic: Vec<bool>,
    pub sites: Vec<AssocSite>,
    pub ideal: Vec<IdealGas>,
    pub simplified_rdf: bool,
}

impl Cpa {
    /// Create a new CPA model for the specified component names.
    ///
    /// Defaults to CPA-SRK and Default parameter sets.
    pub fn new(component_names: &[&str]) -> Result<Self> {
        Self::with_options(component_names, CpaCubicKind::Srk, None)
    }

    /// Create a new CPA model with specified cubic reference and parameter reference.
    pub fn with_options(
        component_names: &[&str],
        cubic_kind: CpaCubicKind,
        preferred_ref: Option<&str>,
    ) -> Result<Self> {
        if component_names.is_empty() {
            return Err(Error::InvalidInput("at least one component required"));
        }

        let nc = component_names.len();
        let mut components = Vec::with_capacity(nc);
        let mut ideal = Vec::with_capacity(nc);
        let mut has_self_assoc = false;
        let mut all_simplified = true;

        for &name in component_names {
            let rec = lookup_cpa_pure(name, preferred_ref).ok_or_else(|| {
                Error::UnknownComponent(format!("CPA parameters not found for {name}"))
            })?;

            let comp_info = component(name)?;
            let id_gas = IdealGas::new(comp_info)?;
            ideal.push(id_gas);

            let tc = comp_info.critical_temperature;

            if rec.scheme != AssocScheme::NoAssoc {
                has_self_assoc = true;
            }
            if !rec.simplified_rdf {
                all_simplified = false;
            }

            let (a0, b, c1) = if rec.scheme == AssocScheme::NoAssoc {
                let (oa, ob) = match cubic_kind {
                    CpaCubicKind::Srk => (0.4274802335655353, 0.0866403499649577),
                    CpaCubicKind::Pr => (0.4572355289213822, 0.0777960739038885),
                };
                let a0 = oa
                    * RGAS
                    * RGAS
                    * comp_info.critical_temperature
                    * comp_info.critical_temperature
                    / comp_info.critical_pressure;
                let b = ob * RGAS * comp_info.critical_temperature / comp_info.critical_pressure;
                let omega = comp_info.acentric_factor;
                let c1 = match cubic_kind {
                    CpaCubicKind::Srk => 0.480 + 1.574 * omega - 0.176 * omega * omega,
                    CpaCubicKind::Pr => 0.37464 + 1.54226 * omega - 0.26992 * omega * omega,
                };
                (a0, b, c1)
            } else {
                (rec.a0 * 1e-6, rec.b * 1e-3, rec.c1)
            };

            components.push(CpaComponent {
                name: name.to_ascii_uppercase(),
                a0,
                b,
                tc,
                c1,
                eps: if rec.scheme == AssocScheme::NoAssoc {
                    0.0
                } else {
                    rec.eps
                },
                beta: if rec.scheme == AssocScheme::NoAssoc {
                    0.0
                } else {
                    rec.beta
                },
                scheme: rec.scheme,
                simplified_rdf: rec.simplified_rdf,
                molar_mass: comp_info.molar_mass,
            });
        }

        if !has_self_assoc {
            return Err(Error::InvalidInput(
                "at least one CPA component must self-associate",
            ));
        }

        // Build binary interaction parameter matrices (nc x nc)
        let mut kij_a = vec![0.0; nc * nc];
        let mut kij_eps = vec![0.0; nc * nc];
        let mut kij_beta = vec![0.0; nc * nc];
        let mut eps_comb_arithmetic = vec![true; nc * nc];
        let mut beta_comb_arithmetic = vec![false; nc * nc];

        for i in 0..nc {
            for j in 0..nc {
                let idx = i * nc + j;
                if i == j {
                    continue;
                }
                if let Some(b_rec) = lookup_cpa_binary(&components[i].name, &components[j].name) {
                    kij_a[idx] = b_rec.kij_a;
                    kij_eps[idx] = b_rec.kij_eps;
                    kij_beta[idx] = b_rec.kij_beta;
                    eps_comb_arithmetic[idx] = b_rec.eps_comb_arithmetic;
                    beta_comb_arithmetic[idx] = b_rec.beta_comb_arithmetic;
                }
            }
        }

        // Bookkeep association sites
        let mut sites = Vec::new();
        for (i, comp) in components.iter().enumerate() {
            let n_sites = comp.scheme.num_sites();
            for s in 1..=n_sites {
                let pol = site_polarity(s, comp.scheme);
                sites.push(AssocSite {
                    comp_idx: i,
                    site_idx: s,
                    polarity: pol,
                });
            }
        }

        Ok(Self {
            cubic_kind,
            components,
            kij_a,
            kij_eps,
            kij_beta,
            eps_comb_arithmetic,
            beta_comb_arithmetic,
            sites,
            ideal,
            simplified_rdf: all_simplified,
        })
    }

    /// Number of components.
    #[inline]
    pub fn num_components(&self) -> usize {
        self.components.len()
    }

    /// Total number of active association sites across all components.
    #[inline]
    pub fn num_sites(&self) -> usize {
        self.sites.len()
    }

    /// Get binary interaction parameter kij for pair (i, j).
    #[inline]
    pub fn get_kij(&self, i: usize, j: usize) -> f64 {
        self.kij_a[i * self.num_components() + j]
    }

    /// Set binary interaction parameter kij for pair (i, j).
    pub fn set_kij(&mut self, i: usize, j: usize, val: f64) {
        let nc = self.num_components();
        self.kij_a[i * nc + j] = val;
        self.kij_a[j * nc + i] = val;
    }

    /// Prepare thermodynamic state at specified temperature and mole fractions.
    pub fn prepare(&self, temperature: f64, mole_fractions: &[f64]) -> Result<PreparedCpa<'_>> {
        PreparedCpa::new(self, temperature, mole_fractions)
    }
}

/// Evaluates site polarity according to ThermoPack `AssocSchemeUtils::polarity`.
fn site_polarity(site: usize, scheme: AssocScheme) -> i8 {
    match scheme {
        AssocScheme::Scheme1 => -1,
        AssocScheme::Scheme2B => {
            if site == 1 {
                -1
            } else {
                1
            }
        }
        AssocScheme::Scheme2C => {
            if site == 1 {
                0
            } else {
                -1
            }
        }
        AssocScheme::Scheme3B => {
            if site <= 2 {
                1
            } else {
                -1
            }
        }
        AssocScheme::Scheme4C => {
            if site <= 2 {
                -1
            } else {
                1
            }
        }
        AssocScheme::Scheme1ea => 1,
        _ => 0,
    }
}

/// Checks if two sites on the same component interact internally.
fn site_interaction_internal(site1: usize, site2: usize, scheme: AssocScheme) -> bool {
    match scheme {
        AssocScheme::Scheme1 => true,
        AssocScheme::Scheme2A => true,
        AssocScheme::Scheme2B => site1 != site2,
        AssocScheme::Scheme2C => site1 == 1 || site2 == 1,
        AssocScheme::Scheme3A => true,
        AssocScheme::Scheme3B => (site1 == 3) != (site2 == 3),
        AssocScheme::Scheme4A => true,
        AssocScheme::Scheme4B => (site1 == 4) != (site2 == 4),
        AssocScheme::Scheme4C => (site1 <= 2 && site2 >= 3) || (site2 <= 2 && site1 >= 3),
        _ => false,
    }
}

/// Checks if two sites on different components cross-interact.
fn cross_site_interaction(
    site1: usize,
    site2: usize,
    scheme1: AssocScheme,
    scheme2: AssocScheme,
) -> bool {
    if scheme1 == AssocScheme::NoAssoc || scheme2 == AssocScheme::NoAssoc {
        return false;
    }
    let p1 = site_polarity(site1, scheme1);
    let p2 = site_polarity(site2, scheme2);
    (p1 as i32) * (p2 as i32) <= 0
}

/// Prepared CPA state for allocation-free evaluations at fixed T and composition.
pub struct PreparedCpa<'a> {
    pub model: &'a Cpa,
    pub temperature: f64,
    pub mole_fractions: Vec<f64>,
    /// Pure component alpha parameters alpha_i(T)
    pub alpha_i: Vec<f64>,
    /// Pure component temperature-dependent energy parameters a_i(T) [Pa * m^6 / mol^2]
    pub a_i: Vec<f64>,
    /// Mixture covolume b_mix [m^3 / mol]
    pub b_mix: f64,
    /// Mixture energy parameter a_mix(T) [Pa * m^6 / mol^2]
    pub a_mix: f64,
    /// Temperature derivative da_mix/dT
    pub da_mix_dt: f64,
    /// Second temperature derivative d2a_mix/dT2
    pub d2a_mix_dt2: f64,
    /// Molar mass of mixture [g / mol]
    pub molar_mass_mix: f64,
    /// Precomputed site interaction matrix: whether site k interacts with site l
    pub interacts: Vec<bool>,
    /// Association energy matrix eps_kl [J / mol]
    pub eps_kl: Vec<f64>,
    /// Association volume matrix beta_kl [-]
    pub beta_kl: Vec<f64>,
    /// Cross-covolume matrix covol_kl [m^3 / mol]
    pub covol_kl: Vec<f64>,
    /// Boltzmann factor matrix exp(eps_kl / (R * T)) - 1.0
    pub boltzmann_factor: Vec<f64>,
}

impl<'a> PreparedCpa<'a> {
    /// Create and prepare a CPA state at temperature T [K] and mole fractions x.
    pub fn new(model: &'a Cpa, temperature: f64, mole_fractions: &[f64]) -> Result<Self> {
        if temperature <= 0.0 || !temperature.is_finite() {
            return Err(Error::InvalidInput("temperature must be strictly positive"));
        }
        let nc = model.num_components();
        if mole_fractions.len() != nc {
            return Err(Error::InvalidInput(
                "mole fractions length must match number of components",
            ));
        }

        let sum_x: f64 = mole_fractions.iter().sum();
        if sum_x <= 0.0 || !sum_x.is_finite() {
            return Err(Error::InvalidInput(
                "mole fractions sum must be strictly positive",
            ));
        }
        let x: Vec<f64> = mole_fractions.iter().map(|&xi| xi / sum_x).collect();

        // 1. Evaluate pure component alpha and a_i(T)
        let mut alpha_i = Vec::with_capacity(nc);
        let mut a_i = Vec::with_capacity(nc);
        let mut da_i_dt = Vec::with_capacity(nc);
        let mut d2a_i_dt2 = Vec::with_capacity(nc);

        for comp in &model.components {
            let tr = temperature / comp.tc;
            let sqrt_tr = tr.sqrt();
            let s = comp.c1;
            let xx = s * s;
            let yy = -2.0 * s * (s + 1.0);

            let alpha = 1.0 + xx * (tr - 1.0) + yy * (sqrt_tr - 1.0);
            let dalpha_dt = xx / comp.tc + 0.5 * yy / (temperature * comp.tc).sqrt();
            let d2alpha_dt2 = -yy / (4.0 * temperature * (temperature * comp.tc).sqrt());

            alpha_i.push(alpha);
            a_i.push(comp.a0 * alpha);
            da_i_dt.push(comp.a0 * dalpha_dt);
            d2a_i_dt2.push(comp.a0 * d2alpha_dt2);
        }

        // 2. Mixture cubic parameters
        let mut b_mix = 0.0;
        let mut molar_mass_mix = 0.0;
        for (i, &xi) in x.iter().enumerate() {
            b_mix += xi * model.components[i].b;
            molar_mass_mix += xi * model.components[i].molar_mass;
        }

        let mut a_mix = 0.0;
        let mut da_mix_dt = 0.0;
        let mut d2a_mix_dt2 = 0.0;

        for i in 0..nc {
            for j in 0..nc {
                let idx = i * nc + j;
                let kij = model.kij_a[idx];
                let a_ij = (a_i[i] * a_i[j]).sqrt() * (1.0 - kij);
                let da_ij_dt = 0.5 * (da_i_dt[i] * a_i[j] + a_i[i] * da_i_dt[j])
                    / (a_i[i] * a_i[j]).sqrt()
                    * (1.0 - kij);

                let term1 = da_i_dt[i] * a_i[j] + a_i[i] * da_i_dt[j];
                let dterm1 =
                    d2a_i_dt2[i] * a_i[j] + 2.0 * da_i_dt[i] * da_i_dt[j] + a_i[i] * d2a_i_dt2[j];
                let denom = (a_i[i] * a_i[j]).sqrt();
                let d2a_ij_dt2 = (0.5 * dterm1 / denom
                    - 0.25 * term1 * term1 / (denom * a_i[i] * a_i[j]))
                    * (1.0 - kij);

                a_mix += x[i] * x[j] * a_ij;
                da_mix_dt += x[i] * x[j] * da_ij_dt;
                d2a_mix_dt2 += x[i] * x[j] * d2a_ij_dt2;
            }
        }

        // 3. Precompute site-pair association constants
        let ns = model.sites.len();
        let mut interacts = vec![false; ns * ns];
        let mut eps_kl = vec![0.0; ns * ns];
        let mut beta_kl = vec![0.0; ns * ns];
        let mut covol_kl = vec![0.0; ns * ns];
        let mut boltzmann_factor = vec![0.0; ns * ns];

        for k in 0..ns {
            for l in 0..ns {
                let idx = k * ns + l;
                let sk = &model.sites[k];
                let sl = &model.sites[l];
                let ic = sk.comp_idx;
                let jc = sl.comp_idx;

                let pair_interacts = if ic == jc {
                    site_interaction_internal(sk.site_idx, sl.site_idx, model.components[ic].scheme)
                } else {
                    cross_site_interaction(
                        sk.site_idx,
                        sl.site_idx,
                        model.components[ic].scheme,
                        model.components[jc].scheme,
                    )
                };

                if pair_interacts {
                    interacts[idx] = true;
                    let (eps_pair, beta_pair) = if ic == jc {
                        (model.components[ic].eps, model.components[ic].beta)
                    } else {
                        let c_idx = ic * nc + jc;
                        let kij_e = model.kij_eps[c_idx];
                        let kij_b = model.kij_beta[c_idx];

                        let eps_p = if model.eps_comb_arithmetic[c_idx] {
                            0.5 * (model.components[ic].eps + model.components[jc].eps)
                        } else {
                            (model.components[ic].eps * model.components[jc].eps).sqrt()
                        } * (1.0 - kij_e);

                        let beta_p = if model.beta_comb_arithmetic[c_idx] {
                            0.5 * (model.components[ic].beta + model.components[jc].beta)
                        } else {
                            (model.components[ic].beta * model.components[jc].beta).sqrt()
                        } * (1.0 - kij_b);

                        (eps_p, beta_p)
                    };

                    let covol = 0.5 * (model.components[ic].b + model.components[jc].b);
                    eps_kl[idx] = eps_pair;
                    beta_kl[idx] = beta_pair;
                    covol_kl[idx] = covol;

                    let arg = eps_pair / (RGAS * temperature);
                    boltzmann_factor[idx] = arg.exp() - 1.0;
                }
            }
        }

        Ok(Self {
            model,
            temperature,
            mole_fractions: x,
            alpha_i,
            a_i,
            b_mix,
            a_mix,
            da_mix_dt,
            d2a_mix_dt2,
            molar_mass_mix,
            interacts,
            eps_kl,
            beta_kl,
            covol_kl,
            boltzmann_factor,
        })
    }

    /// Hard-sphere contact radial distribution function g(eta) and its derivatives.
    #[inline]
    pub fn g_rdf(&self, rho: f64) -> (f64, f64, f64) {
        let eta = 0.25 * self.b_mix * rho;
        if self.model.simplified_rdf {
            let den = 1.0 - 1.9 * eta;
            let d1 = 1.0 / den;
            let g = d1;
            let dg_deta = 1.9 * d1 * d1;
            let d2g_deta2 = 2.0 * 1.9 * 1.9 * d1.powi(3);
            (g, dg_deta, d2g_deta2)
        } else {
            let om_eta = 1.0 - eta;
            let d1 = 1.0 / om_eta;
            let g = (1.0 - 0.5 * eta) * d1.powi(3);
            let dg_deta = (2.5 - eta) * d1.powi(4);
            let d2g_deta2 = (9.0 - 3.0 * eta) * d1.powi(5);
            (g, dg_deta, d2g_deta2)
        }
    }

    /// Solve for non-bonded site fractions X_k at density rho [mol / m^3].
    pub fn solve_site_fractions(&self, rho: f64) -> Result<Vec<f64>> {
        let ns = self.model.sites.len();
        if ns == 0 {
            return Ok(Vec::new());
        }

        let (g, _, _) = self.g_rdf(rho);

        // Fast analytical branch for single-component 2B and 4C
        if self.model.components.len() == 1 {
            let scheme = self.model.components[0].scheme;
            if scheme == AssocScheme::Scheme2B {
                // Alcohol 2B: 1 interacting pair between site 1 and site 2
                let idx = 1; // site 0 and site 1
                let delta = g * self.beta_kl[idx] * self.covol_kl[idx] * self.boltzmann_factor[idx];
                let rho_del = rho * delta;
                let x_val = if rho_del < 1e-12 {
                    1.0 - rho_del + 2.0 * rho_del * rho_del
                } else {
                    ((1.0 + 4.0 * rho_del).sqrt() - 1.0) / (2.0 * rho_del)
                };
                return Ok(vec![x_val, x_val]);
            } else if scheme == AssocScheme::Scheme4C {
                // Water 4C: 2 interacting bonds per site
                let idx = 2; // site 0 interacts with site 2, 3
                let delta = g * self.beta_kl[idx] * self.covol_kl[idx] * self.boltzmann_factor[idx];
                let rho_del = rho * delta;
                let x_val = if rho_del < 1e-12 {
                    1.0 - 2.0 * rho_del + 8.0 * rho_del * rho_del
                } else {
                    ((1.0 + 8.0 * rho_del).sqrt() - 1.0) / (4.0 * rho_del)
                };
                return Ok(vec![x_val, x_val, x_val, x_val]);
            }
        }

        // General multi-site Newton solver
        // F_k(X) = X_k * (1 + rho * sum_l x_l * Delta_kl * X_l) - 1 = 0
        let mut x_k = vec![0.5; ns];
        let mut delta_mat = vec![0.0; ns * ns];
        for k in 0..ns {
            for l in 0..ns {
                let idx = k * ns + l;
                if self.interacts[idx] {
                    delta_mat[idx] =
                        g * self.beta_kl[idx] * self.covol_kl[idx] * self.boltzmann_factor[idx];
                }
            }
        }

        let max_iter = 30;
        let tol = 1e-12;

        for _ in 0..max_iter {
            // Compute residual vector r and Jacobian matrix J
            let mut r = vec![0.0; ns];
            let mut j_mat = vec![0.0; ns * ns];

            for k in 0..ns {
                let mut sum_term = 0.0;
                for l in 0..ns {
                    let comp_l = self.model.sites[l].comp_idx;
                    sum_term += self.mole_fractions[comp_l] * delta_mat[k * ns + l] * x_k[l];
                }
                r[k] = x_k[k] * (1.0 + rho * sum_term) - 1.0;

                for l in 0..ns {
                    let comp_l = self.model.sites[l].comp_idx;
                    let mut deriv =
                        rho * x_k[k] * self.mole_fractions[comp_l] * delta_mat[k * ns + l];
                    if k == l {
                        deriv += 1.0 + rho * sum_term;
                    }
                    j_mat[k * ns + l] = deriv;
                }
            }

            // Check residual convergence
            let max_r = r.iter().fold(0.0_f64, |acc, &val| acc.max(val.abs()));
            if max_r < tol {
                break;
            }

            // Solve J * dx = -r via Gaussian elimination
            let dx = solve_linear_system(&j_mat, &r, ns)?;
            for k in 0..ns {
                x_k[k] = (x_k[k] - dx[k]).clamp(1e-12, 1.0);
            }
        }

        Ok(x_k)
    }

    /// Residual Helmholtz energy terms at molar density rho [mol / m^3].
    pub fn helmholtz_terms(&self, rho: f64) -> Result<(f64, f64, f64)> {
        let b = self.b_mix;
        let brho = b * rho;
        if brho >= 1.0 {
            return Err(Error::InvalidInput("density exceeds maximum packing"));
        }

        let t = self.temperature;
        let rt = RGAS * t;
        let a = self.a_mix;

        // Cubic part
        let a_cubic = match self.model.cubic_kind {
            CpaCubicKind::Srk => -(1.0 - brho).ln() - (a / (b * rt)) * (1.0 + brho).ln(),
            CpaCubicKind::Pr => {
                let s2 = std::f64::consts::SQRT_2;
                -(1.0 - brho).ln()
                    - (a / (2.0 * s2 * b * rt))
                        * ((1.0 + (1.0 + s2) * brho) / (1.0 + (1.0 - s2) * brho)).ln()
            }
        };

        // Association part
        let mut a_assoc = 0.0;
        let ns = self.model.sites.len();
        if ns > 0 {
            let x_k = self.solve_site_fractions(rho)?;
            for (k, &xk) in x_k.iter().enumerate() {
                let comp_k = self.model.sites[k].comp_idx;
                let x_frac = xk.max(1e-30);
                a_assoc += self.mole_fractions[comp_k] * (x_frac.ln() - 0.5 * x_frac + 0.5);
            }
        }

        let a_res = a_cubic + a_assoc;
        Ok((a_cubic, a_assoc, a_res))
    }

    /// Pressure and derivatives at density rho [mol / m^3].
    pub fn pressure_and_derivatives(&self, rho: f64) -> Result<CpaPressureJet> {
        let v = 1.0 / rho;
        let t = self.temperature;

        let hv = 1e-5 * v;
        let rho_p = 1.0 / (v + hv);
        let rho_m = 1.0 / (v - hv);

        let (_, _, a_res) = self.helmholtz_terms(rho)?;
        let (_, _, a_res_p) = self.helmholtz_terms(rho_p)?;
        let (_, _, a_res_m) = self.helmholtz_terms(rho_m)?;

        let f_v = (a_res_p - a_res_m) / (2.0 * hv);
        let f_vv = (a_res_p - 2.0 * a_res + a_res_m) / (hv * hv);

        let p = rho * RGAS * t - RGAS * t * f_v;
        let z = p / (rho * RGAS * t);

        let dp_dv = -rho * rho * RGAS * t - RGAS * t * f_vv;
        let dp_drho = -dp_dv / (rho * rho);

        // Temperature derivative
        let ht = 1e-4 * t;
        let prep_tp = PreparedCpa::new(self.model, t + ht, &self.mole_fractions)?;
        let prep_tm = PreparedCpa::new(self.model, t - ht, &self.mole_fractions)?;

        let (_, _, a_tp_p) = prep_tp.helmholtz_terms(rho_p)?;
        let (_, _, a_tp_m) = prep_tp.helmholtz_terms(rho_m)?;
        let (_, _, a_tm_p) = prep_tm.helmholtz_terms(rho_p)?;
        let (_, _, a_tm_m) = prep_tm.helmholtz_terms(rho_m)?;

        let f_v_tp = (a_tp_p - a_tp_m) / (2.0 * hv);
        let f_v_tm = (a_tm_p - a_tm_m) / (2.0 * hv);
        let f_tv = (f_v_tp - f_v_tm) / (2.0 * ht);

        let dp_dt = rho * RGAS - RGAS * f_v - RGAS * t * f_tv;

        Ok(CpaPressureJet {
            pressure: p,
            compressibility: z,
            residual_helmholtz: a_res,
            dp_drho,
            dp_dv,
            dp_dt,
        })
    }

    /// Fugacity coefficients ln(phi_i) at density rho [mol / m^3].
    pub fn fugacity_coefficients(&self, rho: f64, z: f64) -> Result<Vec<f64>> {
        let nc = self.model.num_components();
        if nc == 1 {
            let (_, _, a_res) = self.helmholtz_terms(rho)?;
            let ln_phi = a_res + (z - 1.0) - z.ln();
            return Ok(vec![ln_phi]);
        }

        let mut ln_phi = Vec::with_capacity(nc);
        let dn = 1e-6;
        let v = 1.0 / rho;

        for i in 0..nc {
            let mut n_p = self.mole_fractions.clone();
            let mut n_m = self.mole_fractions.clone();
            n_p[i] += dn;
            n_m[i] -= dn;

            let sumn_p: f64 = n_p.iter().sum();
            let sumn_m: f64 = n_m.iter().sum();

            let x_p: Vec<f64> = n_p.iter().map(|&n| n / sumn_p).collect();
            let x_m: Vec<f64> = n_m.iter().map(|&n| n / sumn_m).collect();

            let prep_p = PreparedCpa::new(self.model, self.temperature, &x_p)?;
            let prep_m = PreparedCpa::new(self.model, self.temperature, &x_m)?;

            let rho_p = sumn_p / v;
            let rho_m = sumn_m / v;

            let (_, _, a_res_p) = prep_p.helmholtz_terms(rho_p)?;
            let (_, _, a_res_m) = prep_m.helmholtz_terms(rho_m)?;

            let f_tot_p = sumn_p * a_res_p;
            let f_tot_m = sumn_m * a_res_m;

            let d_f_dn = (f_tot_p - f_tot_m) / (2.0 * dn);
            ln_phi.push(d_f_dn - z.ln());
        }

        Ok(ln_phi)
    }

    /// Complete thermodynamic state properties at density rho [mol / m^3].
    pub fn state(&self, rho: f64) -> Result<CpaState> {
        let jet = self.pressure_and_derivatives(rho)?;
        let t = self.temperature;
        let ht = 1e-4 * t;

        let (_, _, a_res) = self.helmholtz_terms(rho)?;
        let prep_tp = PreparedCpa::new(self.model, t + ht, &self.mole_fractions)?;
        let prep_tm = PreparedCpa::new(self.model, t - ht, &self.mole_fractions)?;
        let (_, _, a_tp) = prep_tp.helmholtz_terms(rho)?;
        let (_, _, a_tm) = prep_tm.helmholtz_terms(rho)?;

        let f_t = (a_tp - a_tm) / (2.0 * ht);
        let f_tt = (a_tp - 2.0 * a_res + a_tm) / (ht * ht);

        let u_res = -RGAS * t * t * f_t;
        let h_res = u_res + RGAS * t * (jet.compressibility - 1.0);
        let s_res = -RGAS * (t * f_t + a_res) + RGAS * jet.compressibility.ln();
        let a_res_energy = RGAS * t * a_res;
        let g_res = h_res - t * s_res;

        let cv_res = -RGAS * (2.0 * t * f_t + t * t * f_tt);

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
        let cp_total = cv_total - t * jet.dp_dt.powi(2) / jet.dp_dv;
        let gamma = cp_total / cv_total;
        let mw_kg = self.molar_mass_mix * 1e-3;
        let sound_speed = if jet.dp_drho > 0.0 && gamma > 0.0 && mw_kg > 0.0 {
            (gamma * jet.dp_drho / mw_kg).sqrt()
        } else {
            0.0
        };

        let ln_phi = self.fugacity_coefficients(rho, jet.compressibility)?;

        Ok(CpaState {
            pressure: jet.pressure,
            compressibility: jet.compressibility,
            residual_helmholtz: a_res,
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
            sound_speed,
            ln_phi,
        })
    }

    /// Single-phase volume root for a specified pressure and phase.
    pub fn volume(&self, pressure_target: f64, phase: Phase) -> Result<f64> {
        if pressure_target <= 0.0 || !pressure_target.is_finite() {
            return Err(Error::InvalidInput(
                "pressure target must be strictly positive",
            ));
        }

        let b = self.b_mix;
        let v_min = b * 1.001;

        match phase {
            Phase::Vapor => {
                let mut v = RGAS * self.temperature / pressure_target;
                v = v.max(v_min * 2.0);

                for _ in 0..60 {
                    let rho = 1.0 / v;
                    let jet = self.pressure_and_derivatives(rho)?;
                    let diff = jet.pressure - pressure_target;
                    if diff.abs() < 1e-6 * pressure_target || diff.abs() < 0.05 {
                        return Ok(v);
                    }
                    let dv = -diff / jet.dp_dv;
                    if dv.abs() < 1e-12 * v {
                        return Ok(v);
                    }
                    let mut step = dv.clamp(-0.5 * v, 0.5 * v);
                    if v + step <= v_min {
                        step = 0.5 * (v_min - v);
                    }
                    v += step;
                }
                Ok(v)
            }
            Phase::Liquid => {
                let mut v = b * 1.05;

                for _ in 0..80 {
                    let rho = 1.0 / v;
                    let jet = self.pressure_and_derivatives(rho)?;
                    let diff = jet.pressure - pressure_target;
                    if diff.abs() < 1e-6 * pressure_target || diff.abs() < 0.05 {
                        return Ok(v);
                    }
                    let dv = -diff / jet.dp_dv;
                    if dv.abs() < 1e-12 * v {
                        return Ok(v);
                    }
                    let mut step = dv.clamp(-0.2 * (v - v_min), 0.5 * v);
                    if v + step <= v_min {
                        step = 0.5 * (v_min - v);
                    }
                    v += step;
                }
                Ok(v)
            }
            Phase::MinimumGibbs => {
                let vl = self.volume(pressure_target, Phase::Liquid)?;
                let vv = self.volume(pressure_target, Phase::Vapor)?;
                if (vl - vv).abs() < 1e-6 * vl {
                    return Ok(vl);
                }
                let st_l = self.state(1.0 / vl)?;
                let st_v = self.state(1.0 / vv)?;
                if st_l.g_residual <= st_v.g_residual {
                    Ok(vl)
                } else {
                    Ok(vv)
                }
            }
        }
    }

    /// Pure fluid saturation pressure and coexistence volumes at temperature T [K].
    pub fn saturation_pressure(&self, temperature: f64) -> Result<CpaSaturation> {
        if self.model.num_components() != 1 {
            return Err(Error::InvalidInput(
                "saturation calculation requires pure fluid",
            ));
        }

        // Initial guess from known correlation or 1 bar
        let mut p = 101_325.0;

        for _ in 0..50 {
            let vl = self.volume(p, Phase::Liquid)?;
            let vv = self.volume(p, Phase::Vapor)?;

            if (vl - vv).abs() < 1e-6 * vl {
                return Err(Error::NotConverged(
                    "liquid and vapor roots converged to identical state",
                ));
            }

            let jet_l = self.pressure_and_derivatives(1.0 / vl)?;
            let jet_v = self.pressure_and_derivatives(1.0 / vv)?;

            let phi_l = self.fugacity_coefficients(1.0 / vl, jet_l.compressibility)?[0];
            let phi_v = self.fugacity_coefficients(1.0 / vv, jet_v.compressibility)?[0];

            let delta_phi = phi_l - phi_v;
            if delta_phi.abs() < 1e-7 {
                return Ok(CpaSaturation {
                    temperature,
                    pressure: p,
                    liquid_volume: vl,
                    vapor_volume: vv,
                });
            }

            let d_delta_dp = (vl - vv) / (RGAS * temperature);
            let mut dp = -delta_phi / d_delta_dp;
            dp = dp.clamp(-0.3 * p, 0.5 * p);
            p += dp;
            if p <= 1e-2 {
                p = 1e-2;
            }
        }

        Err(Error::NotConverged(
            "saturation pressure failed to converge",
        ))
    }
}

/// Solves linear system A * x = b via Gaussian elimination with partial pivoting.
fn solve_linear_system(a: &[f64], b: &[f64], n: usize) -> Result<Vec<f64>> {
    let mut mat = a.to_vec();
    let mut rhs = b.to_vec();

    for i in 0..n {
        // Find pivot
        let mut max_val = mat[i * n + i].abs();
        let mut max_row = i;
        for k in (i + 1)..n {
            let val = mat[k * n + i].abs();
            if val > max_val {
                max_val = val;
                max_row = k;
            }
        }

        if max_val < 1e-18 {
            return Err(Error::NotConverged(
                "singular matrix in site fraction solver",
            ));
        }

        // Swap rows
        if max_row != i {
            for j in 0..n {
                mat.swap(i * n + j, max_row * n + j);
            }
            rhs.swap(i, max_row);
        }

        // Eliminate
        let pivot = mat[i * n + i];
        for k in (i + 1)..n {
            let factor = mat[k * n + i] / pivot;
            for j in i..n {
                mat[k * n + j] -= factor * mat[i * n + j];
            }
            rhs[k] -= factor * rhs[i];
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = rhs[i];
        for j in (i + 1)..n {
            sum -= mat[i * n + j] * x[j];
        }
        x[i] = sum / mat[i * n + i];
    }

    Ok(x)
}
