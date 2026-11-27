//! Native Rust implementation of CPA (Cubic-Plus-Association).
//!
//! Combines cubic equations of state (SRK, PR) and their mixing rules (`VdW`,
//! `HuronVidal1`, `HuronVidal2`, `WongSandler`, `HvWongSandler`) with the
//! shared Wertheim first-order thermodynamic perturbation theory association
//! engine (`crate::association`).
//!
//! # References
//! - Kontogeorgis et al. (1996, 1999, 2006, 2010): The Cubic-Plus-Association (CPA) Equation of State.
//! - Wertheim (1984, 1986): Fluids with Highly Directional Attractive Forces.
//! - Michelsen & Hendriks (2001): Physical Properties from Association Models.

use crate::alpha::AlphaCorrelation;
use crate::association::{
    self, AssocCombRule, AssocRdfFirst, AssocRdfSecond, AssocRdfValue, AssociationHelmholtz,
    AssociationParams, AssociationRdfProvider, DeltaCombRule,
};
use crate::component::component;
use crate::cpa_data::{
    AssocScheme, lookup_cpa_binary_for_eos, lookup_cpa_pure, lookup_cpa_pure_for_eos,
};
use crate::cubic::{Cubic, CubicHelmholtz, CubicKind, Phase, PreparedState};
use crate::error::{Error, Result, positive};
use crate::excess_gibbs::MixingRule;
use crate::ideal::IdealGas;
use crate::rdf::RdfDerivatives;
use crate::tuning::{HvPairParams, VolumeShiftParams, WsPairParams};
use crate::wong_sandler::WongSandlerData;

/// Gas constant R [J / (mol * K)] matching ThermoPack Fortran exact constant.
pub const RGAS: f64 = crate::R;

/// Underlying cubic equation of state type for CPA.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CpaCubicKind {
    /// CPA-SRK (Soave-Redlich-Kwong reference, u=1, w=0).
    Srk,
    /// CPA-PR (Peng-Robinson reference, u=2, w=-1).
    Pr,
}

impl CpaCubicKind {
    /// Corresponding [`CubicKind`].
    #[inline]
    pub fn cubic_kind(self) -> CubicKind {
        match self {
            Self::Srk => CubicKind::SoaveRedlichKwong,
            Self::Pr => CubicKind::PengRobinson,
        }
    }

    /// ThermoPack CPA EOS identifier (`"CPA-SRK"` or `"CPA-PR"`).
    #[inline]
    pub fn eos_name(self) -> &'static str {
        match self {
            Self::Srk => "CPA-SRK",
            Self::Pr => "CPA-PR",
        }
    }
}

/// Configuration options for [`Cpa::with_cpa_options`].
#[derive(Debug, Clone, Copy)]
pub struct CpaOptions<'a> {
    /// Underlying cubic reference equation (`Srk` or `Pr`).
    pub cubic_kind: CpaCubicKind,
    /// Cubic mixing rule (`VdW`, `HuronVidal1`, `HuronVidal2`, `WongSandler`, `HvWongSandler`).
    pub mixing_rule: MixingRule,
    /// Alpha correlation for non-associating components (defaults to `Classic`).
    pub alpha: Option<AlphaCorrelation>,
    /// Optional parameter reference string (e.g. `"DEFAULT"`, `"SINTEF"`, `"2B-SINTEF"`).
    pub parameter_reference: Option<&'a str>,
    /// Override simplified (`true`) vs original Carnahan-Starling (`false`) hard-sphere RDF.
    /// When `None`, follows the CPA pure-component database convention.
    pub simplified_rdf: Option<bool>,
    /// Combining rule for cross-association strength `Delta` (`Standard` or `Elliot`).
    pub delta_comb_rule: DeltaCombRule,
    /// Optional initial Wong-Sandler binary interaction parameters.
    pub wong_sandler: Option<&'a WongSandlerData>,
}

impl Default for CpaOptions<'_> {
    fn default() -> Self {
        Self {
            cubic_kind: CpaCubicKind::Srk,
            mixing_rule: MixingRule::VdW,
            alpha: None,
            parameter_reference: None,
            simplified_rdf: None,
            delta_comb_rule: DeltaCombRule::Standard,
            wong_sandler: None,
        }
    }
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
    /// Critical temperature Tc (K), used in Soave alpha.
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
    /// Pressure P (Pa).
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
    /// Pressure P (Pa).
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

/// Complete 10-block reduced residual Helmholtz energy `F = A^res / (R T)` and
/// all first and second `(T, V, n)` derivatives for CPA, plus non-bonded site
/// fractions `X_k`.
#[derive(Debug, Clone)]
pub struct CpaHelmholtz {
    /// Dimensionless reduced residual Helmholtz energy `F = A^res / (R T)`.
    pub f: f64,
    /// `dF / dT` [1/K].
    pub f_t: f64,
    /// `dF / dV` [1/m^3].
    pub f_v: f64,
    /// `dF / dn_i` [1/mol].
    pub f_n: Vec<f64>,
    /// `d^2F / dT^2` [1/K^2].
    pub f_tt: f64,
    /// `d^2F / (dT dV)` [1/(K m^3)].
    pub f_tv: f64,
    /// `d^2F / (dT dn_i)` [1/(K mol)].
    pub f_tn: Vec<f64>,
    /// `d^2F / dV^2` [1/m^6].
    pub f_vv: f64,
    /// `d^2F / (dV dn_i)` [1/(m^3 mol)].
    pub f_vn: Vec<f64>,
    /// `d^2F / (dn_i dn_j)` [1/mol^2], row-major `nc * nc`.
    pub f_nn: Vec<f64>,
    /// Non-bonded site fractions `X_k` (length `num_sites`).
    pub x_k: Vec<f64>,
}

/// Scalar thermodynamic quantity and its analytical `(T, P, n)` derivatives.
#[derive(Debug, Clone)]
pub struct CpaDerivatives {
    pub value: f64,
    pub dt: f64,
    pub dp: f64,
    pub dn: Vec<f64>,
}

/// Fugacity coefficients `ln(phi_i)` and their analytical `(T, P, n)` derivatives.
#[derive(Debug, Clone)]
pub struct CpaFugacityDerivatives {
    pub ln_phi: Vec<f64>,
    pub dt: Vec<f64>,
    pub dp: Vec<f64>,
    /// Row-major `nc * nc` matrix `d(ln phi_i) / dn_j`.
    pub dn: Vec<f64>,
}

/// Complete TP thermodynamic state with analytical `(T, P, n)` derivatives.
#[derive(Debug, Clone)]
pub struct CpaTpState {
    pub volume: f64,
    pub compressibility_factor: CpaDerivatives,
    pub specific_volume: CpaDerivatives,
    pub residual_enthalpy: CpaDerivatives,
    pub enthalpy: CpaDerivatives,
    pub residual_entropy: CpaDerivatives,
    pub entropy: CpaDerivatives,
    pub residual_gibbs: CpaDerivatives,
    pub ln_fugacity: CpaFugacityDerivatives,
    pub cp: f64,
    pub cv: f64,
    pub speed_of_sound: f64,
}

/// Saturation state for pure fluids.
#[derive(Debug, Clone, Copy)]
pub struct CpaSaturation {
    /// Saturation temperature T (K).
    pub temperature: f64,
    /// Saturation pressure P_sat (Pa).
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
    cubic: Cubic,
    assoc: AssociationParams,
    delta_comb_rule: DeltaCombRule,
}

impl Cpa {
    /// Create a new CPA-SRK model with van der Waals mixing and default parameters.
    pub fn new(component_names: &[&str]) -> Result<Self> {
        Self::with_options(component_names, CpaCubicKind::Srk, None)
    }

    /// Create a new CPA model with specified cubic reference and parameter reference.
    pub fn with_options(
        component_names: &[&str],
        cubic_kind: CpaCubicKind,
        preferred_ref: Option<&str>,
    ) -> Result<Self> {
        Self::with_cpa_options(
            component_names,
            CpaOptions {
                cubic_kind,
                parameter_reference: preferred_ref,
                ..Default::default()
            },
        )
    }

    /// Create a new CPA model with full cubic mixing, RDF, alpha, and association options.
    pub fn with_cpa_options(component_names: &[&str], options: CpaOptions<'_>) -> Result<Self> {
        if component_names.is_empty() {
            return Err(Error::InvalidInput("at least one component required"));
        }
        if !matches!(
            options.mixing_rule,
            MixingRule::VdW
                | MixingRule::HuronVidal1
                | MixingRule::HuronVidal2
                | MixingRule::WongSandler
                | MixingRule::HvWongSandler
        ) {
            return Err(Error::Unsupported(
                "selected mixing rule is not supported for CPA",
            ));
        }

        let nc = component_names.len();
        let eos_name = options.cubic_kind.eos_name();
        let param_ref = options.parameter_reference.unwrap_or("DEFAULT");
        let base_alpha = options.alpha.unwrap_or(AlphaCorrelation::Classic);

        let mut cubic = Cubic::new_cpa_base(
            options.cubic_kind.cubic_kind(),
            component_names,
            options.mixing_rule,
            base_alpha,
            param_ref,
            options.wong_sandler,
        )?;

        let mut components = Vec::with_capacity(nc);
        let mut ideal = Vec::with_capacity(nc);
        let mut has_self_assoc = false;
        let mut db_simplified_count = 0usize;
        let mut db_all_simplified = true;

        for (i, &name) in component_names.iter().enumerate() {
            let comp_info = component(name)?;
            let id_gas = IdealGas::new(comp_info)?;
            ideal.push(id_gas);

            // Follow Fortran `saft_interface.f90` + `cpa_parameters.f90`:
            // components absent from CPA pure DB default to `no_assoc` and standard cubic.
            // For `CPA-PR`, fallback to EOS-agnostic lookup preserves `Cpa::with_options(..., Pr, ...)`
            // compatibility when no `CPA-PR` pure rows exist in the upstream database.
            let rec = lookup_cpa_pure_for_eos(eos_name, name, options.parameter_reference)
                .or_else(|| lookup_cpa_pure(name, options.parameter_reference));

            let (scheme, simplified_flag, a0, b, c1, eps, beta) = if let Some(rec) = rec {
                db_simplified_count += 1;
                if !rec.simplified_rdf {
                    db_all_simplified = false;
                }
                if rec.scheme != AssocScheme::NoAssoc && rec.scheme != AssocScheme::Scheme1ea {
                    has_self_assoc = true;
                }
                if rec.scheme == AssocScheme::NoAssoc {
                    let a0 = cubic.pure_a0()[i];
                    let b = cubic.pure_bi()[i];
                    let c1 = cubic.alpha(i)?.parameters[0];
                    (rec.scheme, rec.simplified_rdf, a0, b, c1, 0.0, 0.0)
                } else {
                    let a0 = rec.a0 * 1.0e-6;
                    let b = rec.b * 1.0e-3;
                    let c1 = rec.c1;
                    cubic.set_pure_ab_override(i, a0, b)?;
                    cubic.set_alpha(i, AlphaCorrelation::ClassicFit, &[c1])?;
                    (rec.scheme, rec.simplified_rdf, a0, b, c1, rec.eps, rec.beta)
                }
            } else {
                let a0 = cubic.pure_a0()[i];
                let b = cubic.pure_bi()[i];
                let c1 = cubic.alpha(i)?.parameters[0];
                (AssocScheme::NoAssoc, true, a0, b, c1, 0.0, 0.0)
            };

            components.push(CpaComponent {
                name: comp_info.ident.to_string(),
                a0,
                b,
                tc: comp_info.critical_temperature,
                c1,
                eps,
                beta,
                scheme,
                simplified_rdf: simplified_flag,
                molar_mass: comp_info.molar_mass,
            });
        }

        if !has_self_assoc {
            return Err(Error::InvalidInput(
                "at least one CPA component must self-associate",
            ));
        }

        let simplified_rdf = options
            .simplified_rdf
            .unwrap_or(db_simplified_count == 0 || db_all_simplified);

        // Default combining rules from `cpa_parameters.f90:245-250`:
        // CPA-SRK: (AriComb, GeoComb), CPA-PR: (GeoComb, GeoComb).
        let default_eps_ari = options.cubic_kind == CpaCubicKind::Srk;
        let mut kij_a = vec![0.0; nc * nc];
        let mut kij_eps = vec![0.0; nc * nc];
        let mut kij_beta = vec![0.0; nc * nc];
        let mut eps_comb_arithmetic = vec![default_eps_ari; nc * nc];
        let mut beta_comb_arithmetic = vec![false; nc * nc];

        for i in 0..nc {
            for j in (i + 1)..nc {
                let ij = i * nc + j;
                let ji = j * nc + i;
                if let Some(b_rec) =
                    lookup_cpa_binary_for_eos(eos_name, &components[i].name, &components[j].name)
                {
                    kij_a[ij] = b_rec.kij_a;
                    kij_a[ji] = b_rec.kij_a;
                    kij_eps[ij] = b_rec.kij_eps;
                    kij_eps[ji] = b_rec.kij_eps;
                    kij_beta[ij] = b_rec.kij_beta;
                    kij_beta[ji] = b_rec.kij_beta;
                    eps_comb_arithmetic[ij] = b_rec.eps_comb_arithmetic;
                    eps_comb_arithmetic[ji] = b_rec.eps_comb_arithmetic;
                    beta_comb_arithmetic[ij] = b_rec.beta_comb_arithmetic;
                    beta_comb_arithmetic[ji] = b_rec.beta_comb_arithmetic;
                    if options.mixing_rule != MixingRule::HvWongSandler {
                        cubic.set_kij(i, j, b_rec.kij_a)?;
                    }
                }
            }
        }

        let sites = build_sites(&components);
        let assoc = build_cpa_assoc(
            &components,
            &eps_comb_arithmetic,
            &beta_comb_arithmetic,
            &kij_eps,
            &kij_beta,
            options.delta_comb_rule,
        )?;

        Ok(Self {
            cubic_kind: options.cubic_kind,
            components,
            kij_a,
            kij_eps,
            kij_beta,
            eps_comb_arithmetic,
            beta_comb_arithmetic,
            sites,
            ideal,
            simplified_rdf,
            cubic,
            assoc,
            delta_comb_rule: options.delta_comb_rule,
        })
    }

    fn rebuild_assoc(&mut self) -> Result<()> {
        self.sites = build_sites(&self.components);
        self.assoc = build_cpa_assoc(
            &self.components,
            &self.eps_comb_arithmetic,
            &self.beta_comb_arithmetic,
            &self.kij_eps,
            &self.kij_beta,
            self.delta_comb_rule,
        )?;
        Ok(())
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

    /// Underlying cubic equation of state.
    #[inline]
    pub fn cubic(&self) -> &Cubic {
        &self.cubic
    }

    /// Active cubic mixing rule.
    #[inline]
    pub fn mixing_rule(&self) -> MixingRule {
        self.cubic.mixing_rule
    }

    /// Active Wertheim association parameter structure.
    #[inline]
    pub fn association_params(&self) -> &AssociationParams {
        &self.assoc
    }

    /// Active `Delta` cross-association combining rule.
    #[inline]
    pub fn delta_comb_rule(&self) -> DeltaCombRule {
        self.delta_comb_rule
    }

    /// Switch between simplified (`true`) and original Carnahan-Starling (`false`)
    /// hard-sphere contact RDF.
    pub fn set_simplified_rdf(&mut self, simplified: bool) {
        self.simplified_rdf = simplified;
    }

    /// Set the cross-association `Delta` combining rule (`Standard` or `Elliot`).
    pub fn set_delta_comb_rule(&mut self, rule: DeltaCombRule) -> Result<()> {
        self.delta_comb_rule = rule;
        self.rebuild_assoc()
    }

    /// Get binary interaction parameter `kij` for pair `(i, j)`.
    #[inline]
    pub fn get_kij(&self, i: usize, j: usize) -> f64 {
        self.kij_a[i * self.num_components() + j]
    }

    /// Set symmetric binary interaction parameter `kij` for pair `(i, j)`.
    pub fn set_kij(&mut self, i: usize, j: usize, val: f64) {
        let nc = self.num_components();
        self.kij_a[i * nc + j] = val;
        self.kij_a[j * nc + i] = val;
        let _ = self.cubic.set_kij(i, j, val);
    }

    /// Checked symmetric `kij` setter returning an error on invalid indices or non-finite values.
    pub fn try_set_kij(&mut self, i: usize, j: usize, val: f64) -> Result<()> {
        self.cubic.set_kij(i, j, val)?;
        let nc = self.num_components();
        self.kij_a[i * nc + j] = val;
        self.kij_a[j * nc + i] = val;
        Ok(())
    }

    /// Set a single `(i, j)` entry of `kij` without modifying `(j, i)` (`thermopack_setkij`).
    pub fn set_kij_entry(&mut self, i: usize, j: usize, val: f64) -> Result<()> {
        self.cubic.set_kij_entry(i, j, val)?;
        let nc = self.num_components();
        self.kij_a[i * nc + j] = val;
        Ok(())
    }

    /// Get binary covolume interaction parameter `lij` for pair `(i, j)`.
    pub fn lij(&self, i: usize, j: usize) -> Result<f64> {
        self.cubic.lij(i, j)
    }

    /// Set symmetric binary covolume interaction parameter `lij` for pair `(i, j)`.
    pub fn set_lij(&mut self, i: usize, j: usize, val: f64) -> Result<()> {
        self.cubic.set_lij(i, j, val)
    }

    /// Set binary association interaction parameters `(kij_eps, kij_beta)` for pair `(i, j)`.
    pub fn set_kij_eps_beta(
        &mut self,
        i: usize,
        j: usize,
        kij_eps: f64,
        kij_beta: f64,
    ) -> Result<()> {
        let nc = self.num_components();
        if i >= nc || j >= nc || i == j {
            return Err(Error::InvalidInput("distinct component indices required"));
        }
        if !kij_eps.is_finite() || !kij_beta.is_finite() {
            return Err(Error::InvalidInput(
                "finite association interaction parameters required",
            ));
        }
        self.kij_eps[i * nc + j] = kij_eps;
        self.kij_eps[j * nc + i] = kij_eps;
        self.kij_beta[i * nc + j] = kij_beta;
        self.kij_beta[j * nc + i] = kij_beta;
        self.rebuild_assoc()
    }

    /// Set binary association combining rules (`arithmetic = true` / `geometric = false`)
    /// for `(eps, beta)` of pair `(i, j)`.
    pub fn set_assoc_comb_rules(
        &mut self,
        i: usize,
        j: usize,
        eps_arithmetic: bool,
        beta_arithmetic: bool,
    ) -> Result<()> {
        let nc = self.num_components();
        if i >= nc || j >= nc || i == j {
            return Err(Error::InvalidInput("distinct component indices required"));
        }
        self.eps_comb_arithmetic[i * nc + j] = eps_arithmetic;
        self.eps_comb_arithmetic[j * nc + i] = eps_arithmetic;
        self.beta_comb_arithmetic[i * nc + j] = beta_arithmetic;
        self.beta_comb_arithmetic[j * nc + i] = beta_arithmetic;
        self.rebuild_assoc()
    }

    /// Retrieve Huron-Vidal binary parameters for pair `(i, j)` (`thermopack_getHVparam`).
    pub fn hv_param(&self, i: usize, j: usize) -> Result<HvPairParams> {
        self.cubic.hv_param(i, j)
    }

    /// Tune Huron-Vidal binary parameters for pair `(i, j)` (`thermopack_setHVparam`).
    pub fn set_hv_param(&mut self, i: usize, j: usize, params: HvPairParams) -> Result<()> {
        self.cubic.set_hv_param(i, j, params)
    }

    /// Retrieve constant Wong-Sandler binary parameters for pair `(i, j)` (`thermopack_getWSparam`).
    pub fn ws_param(&self, i: usize, j: usize) -> Result<WsPairParams> {
        self.cubic.ws_param(i, j)
    }

    /// Tune constant Wong-Sandler binary parameters for pair `(i, j)` (`thermopack_setWSparam`).
    pub fn set_ws_param(&mut self, i: usize, j: usize, params: WsPairParams) -> Result<()> {
        self.cubic.set_ws_param(i, j, params)
    }

    /// Retrieve alpha parameters `[c1, c2, c3]` of component `i`.
    pub fn twu_param(&self, i: usize) -> Result<[f64; 3]> {
        self.cubic.twu_param(i)
    }

    /// Tune Twu alpha parameters for component `i`.
    pub fn set_twu_param(&mut self, i: usize, c1: f64, c2: f64, c3: f64) -> Result<()> {
        self.cubic.set_twu_param(i, c1, c2, c3)
    }

    /// Retrieve Mathias-Copeman alpha parameters `[c1, c2, c3]` of component `i`.
    pub fn mc_param(&self, i: usize) -> Result<[f64; 3]> {
        self.cubic.mc_param(i)
    }

    /// Tune Mathias-Copeman alpha parameters for component `i`.
    pub fn set_mc_param(&mut self, i: usize, c1: f64, c2: f64, c3: f64) -> Result<()> {
        self.cubic.set_mc_param(i, c1, c2, c3)
    }

    /// Set alpha correlation by name for component `i`.
    pub fn set_alpha_corr(&mut self, i: usize, name: &str, params: &[f64]) -> Result<()> {
        self.cubic.set_alpha_corr(i, name, params)
    }

    /// Set beta correlation by name for component `i`.
    pub fn set_beta_corr(&mut self, i: usize, name: &str, params: &[f64]) -> Result<()> {
        self.cubic.set_beta_corr(i, name, params)
    }

    /// Retrieve volume shift parameters for component `i`.
    pub fn volume_shift_parameters(&self, i: usize) -> Result<VolumeShiftParams> {
        self.cubic.volume_shift_parameters(i)
    }

    /// Set volume shift parameters for component `i`.
    pub fn set_volume_shift_parameters(
        &mut self,
        i: usize,
        params: VolumeShiftParams,
    ) -> Result<()> {
        self.cubic.set_volume_shift_parameters(i, params)
    }

    /// Prepare thermodynamic state at specified temperature and mole fractions.
    pub fn prepare(&self, temperature: f64, mole_fractions: &[f64]) -> Result<PreparedCpa<'_>> {
        PreparedCpa::new(self, temperature, mole_fractions)
    }

    /// Evaluate only the cubic (non-associating) 10-block reduced residual
    /// Helmholtz energy `F^cubic(T, V, n)` and derivatives.
    pub fn nonassoc_residual_helmholtz_tv(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<CubicHelmholtz> {
        self.cubic.residual_helmholtz_tv(temperature, volume, moles)
    }

    /// Evaluate only the Wertheim TPT1 association 10-block reduced residual
    /// Helmholtz energy `F^assoc(T, V, n)` and derivatives.
    pub fn assoc_residual_helmholtz_tv(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<AssociationHelmholtz> {
        let (_total, v_unit, composition) =
            validate_tvn(self.num_components(), temperature, volume, moles)?;
        if !self.assoc.has_sites() {
            let nc = self.num_components();
            return Ok(AssociationHelmholtz {
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
                x_k: Vec::new(),
            });
        }
        let rdf_provider = CpaRdfProvider::new(self, v_unit, &composition)?;
        association::assoc_helmholtz_all(&self.assoc, &rdf_provider, temperature, volume, moles)
    }

    /// Evaluate the complete CPA 10-block reduced residual Helmholtz energy
    /// `F(T, V, n) = F^cubic + F^assoc` and all first/second `(T, V, n)`
    /// derivatives in unit-mole coordinates with exact Euler scaling.
    pub fn residual_helmholtz_tv(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<CpaHelmholtz> {
        let (total, v_unit, composition) =
            validate_tvn(self.num_components(), temperature, volume, moles)?;
        let prepared = PreparedCpa::new(self, temperature, &composition)?;
        let mut h = prepared.residual_helmholtz_tv_unit(v_unit)?;
        if total != 1.0 {
            let inv_n = 1.0 / total;
            h.f *= total;
            h.f_t *= total;
            h.f_tt *= total;
            h.f_vv *= inv_n;
            for v in &mut h.f_vn {
                *v *= inv_n;
            }
            for v in &mut h.f_nn {
                *v *= inv_n;
            }
        }
        let all_finite = [h.f, h.f_t, h.f_v, h.f_tt, h.f_tv, h.f_vv]
            .iter()
            .all(|v| v.is_finite())
            && h.f_n.iter().all(|v| v.is_finite())
            && h.f_tn.iter().all(|v| v.is_finite())
            && h.f_vn.iter().all(|v| v.is_finite())
            && h.f_nn.iter().all(|v| v.is_finite())
            && h.x_k.iter().all(|v| v.is_finite());
        if !all_finite {
            return Err(Error::InvalidInput(
                "non-finite CPA residual Helmholtz derivative",
            ));
        }
        Ok(h)
    }

    /// Solve for non-bonded site fractions `X_k` at `(T, V, n)`.
    pub fn solve_assoc_x_k(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<Vec<f64>> {
        let (_, v_unit, composition) =
            validate_tvn(self.num_components(), temperature, volume, moles)?;
        if !self.assoc.has_sites() {
            return Ok(Vec::new());
        }
        let rdf_provider = CpaRdfProvider::new(self, v_unit, &composition)?;
        association::solve_for_x_k(
            &self.assoc,
            &rdf_provider,
            temperature,
            v_unit,
            &composition,
            None,
        )
    }

    /// Solve for molar volume [m^3/mol] at `(T, P, n)` for the requested `phase`.
    pub fn volume_tp(
        &self,
        temperature: f64,
        pressure: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<f64> {
        let prepared = self.prepare(temperature, moles)?;
        prepared.volume(pressure, phase)
    }

    /// Evaluate the complete TP thermodynamic state with exact analytical
    /// `(T, P, n)` derivatives at `(T, P, n)` for the requested `phase`.
    pub fn state_tp(
        &self,
        temperature: f64,
        pressure: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<CpaTpState> {
        positive(temperature, "temperature")?;
        positive(pressure, "pressure")?;
        let nc = self.num_components();
        if moles.len() != nc || moles.iter().any(|&v| !v.is_finite() || v < 0.0) {
            return Err(Error::InvalidInput(
                "finite non-negative mole numbers matching the model required",
            ));
        }
        let total: f64 = moles.iter().copied().sum();
        positive(total, "total moles")?;
        let x: Vec<f64> = moles.iter().map(|&v| v / total).collect();
        let prepared = PreparedCpa::new(self, temperature, &x)?;
        let v_molar = prepared.volume(pressure, phase)?;
        let v_total = v_molar * total;
        if !v_total.is_finite() || v_total <= 0.0 {
            return Err(Error::InvalidInput(
                "total volume overflowed or underflowed",
            ));
        }

        let h = prepared.residual_helmholtz_tv_unit(v_molar)?;
        let t = temperature;
        let v = v_molar;
        let inv_n = 1.0 / total;
        let p_eval = RGAS * t / v - RGAS * t * h.f_v;
        let p_v = -RGAS * t / (v * v) - RGAS * t * h.f_vv;
        let p_t = p_eval / t - RGAS * t * h.f_tv;
        let p_n: Vec<f64> = h
            .f_vn
            .iter()
            .map(|&fvn| RGAS * t / v - RGAS * t * fvn)
            .collect();

        if !p_v.is_finite() || p_v == 0.0 {
            return Err(Error::InvalidInput("singular pressure volume derivative"));
        }

        let v_t = -p_t / p_v;
        let v_p = 1.0 / p_v;
        let v_n_partial: Vec<f64> = p_n.iter().map(|&pni| -pni / p_v).collect();

        let z = p_eval * v / (RGAS * t);
        positive(z, "compressibility factor")?;
        let ln_z = z.ln();
        let z_t = (p_eval / (RGAS * t)) * (v_t - v / t);
        let z_p = (v + p_eval * v_p) / (RGAS * t);
        let z_dn: Vec<f64> = v_n_partial
            .iter()
            .map(|&vni| (p_eval / (RGAS * t)) * (vni - v) * inv_n)
            .collect();

        let h_res = -RGAS * t * t * h.f_t + p_eval * v - RGAS * t;
        let cv_res = -RGAS * (2.0 * t * h.f_t + t * t * h.f_tt);
        let cp_res = cv_res - t * p_t * p_t / p_v - RGAS;
        let h_res_t = cp_res;
        let h_res_p = v - t * v_t;
        let mut h_res_dn = Vec::with_capacity(nc);
        for i in 0..nc {
            h_res_dn.push(-RGAS * t * t * h.f_tn[i] + t * p_t * v_n_partial[i] - RGAS * t);
        }

        let s_res = -RGAS * (h.f + t * h.f_t - ln_z);
        let s_res_t = cp_res / t;
        let s_res_p = -v_t + RGAS / p_eval;
        let mut s_res_dn = Vec::with_capacity(nc);
        for i in 0..nc {
            s_res_dn.push(-RGAS * (h.f_n[i] + t * h.f_tn[i] - ln_z) + p_t * v_n_partial[i] - RGAS);
        }

        let g_res = h_res - t * s_res;
        let g_res_t = -s_res;
        let g_res_p = v - RGAS * t / p_eval;
        let g_res_dn: Vec<f64> = h.f_n.iter().map(|&fni| RGAS * t * (fni - ln_z)).collect();

        let mut h_id = 0.0;
        let mut cp_id = 0.0;
        let mut s_id = 0.0;
        let mut h_id_comp = Vec::with_capacity(nc);
        let mut s_id_comp = Vec::with_capacity(nc);
        for i in 0..nc {
            let id = self.ideal[i].properties(t, p_eval)?;
            h_id += x[i] * id.enthalpy;
            cp_id += x[i] * id.cp;
            if x[i] > 0.0 {
                s_id += x[i] * (id.entropy - RGAS * x[i].ln());
            }
            h_id_comp.push(id.enthalpy);
            s_id_comp.push(id.entropy - RGAS * x[i].max(1.0e-100).ln());
        }
        let mut h_tot_dn = Vec::with_capacity(nc);
        let mut s_tot_dn = Vec::with_capacity(nc);
        for i in 0..nc {
            h_tot_dn.push(h_res_dn[i] + h_id_comp[i]);
            s_tot_dn.push(s_res_dn[i] + s_id_comp[i]);
        }

        let mut ln_phi = Vec::with_capacity(nc);
        let mut dlnphi_dt = Vec::with_capacity(nc);
        let mut dlnphi_dp = Vec::with_capacity(nc);
        let mut dlnphi_dn = vec![0.0; nc * nc];
        for i in 0..nc {
            ln_phi.push(h.f_n[i] - ln_z);
            dlnphi_dt.push(h.f_tn[i] + 1.0 / t - v_n_partial[i] * p_t / (RGAS * t));
            dlnphi_dp.push(v_n_partial[i] / (RGAS * t) - 1.0 / p_eval);
            for j in 0..nc {
                let d_unit = h.f_nn[i * nc + j] + 1.0 + p_n[i] * p_n[j] / (RGAS * t * p_v);
                dlnphi_dn[i * nc + j] = d_unit * inv_n;
            }
        }

        let cv = cp_id - RGAS + cv_res;
        let cp = cp_id + cp_res;
        let dp_drho = -p_v * v * v;
        let mw_kg = prepared.molar_mass_mix * 1.0e-3;
        let sound_speed = if dp_drho > 0.0 && cp > 0.0 && cv > 0.0 && mw_kg > 0.0 {
            ((cp / cv) * dp_drho / mw_kg).sqrt()
        } else {
            0.0
        };

        let all_finite = [v, z, z_t, z_p, h_res, s_res, g_res, cp, cv, sound_speed]
            .iter()
            .all(|q| q.is_finite())
            && ln_phi.iter().all(|q| q.is_finite())
            && dlnphi_dt.iter().all(|q| q.is_finite())
            && dlnphi_dp.iter().all(|q| q.is_finite())
            && dlnphi_dn.iter().all(|q| q.is_finite());
        if !all_finite {
            return Err(Error::InvalidInput("non-finite CPA TP thermodynamic state"));
        }

        Ok(CpaTpState {
            volume: v,
            compressibility_factor: CpaDerivatives {
                value: z,
                dt: z_t,
                dp: z_p,
                dn: z_dn,
            },
            specific_volume: CpaDerivatives {
                value: v,
                dt: v_t,
                dp: v_p,
                dn: v_n_partial,
            },
            residual_enthalpy: CpaDerivatives {
                value: h_res,
                dt: h_res_t,
                dp: h_res_p,
                dn: h_res_dn,
            },
            enthalpy: CpaDerivatives {
                value: h_id + h_res,
                dt: cp,
                dp: h_res_p,
                dn: h_tot_dn,
            },
            residual_entropy: CpaDerivatives {
                value: s_res,
                dt: s_res_t,
                dp: s_res_p,
                dn: s_res_dn,
            },
            entropy: CpaDerivatives {
                value: s_id + s_res,
                dt: cp / t,
                dp: -v_t,
                dn: s_tot_dn,
            },
            residual_gibbs: CpaDerivatives {
                value: g_res,
                dt: g_res_t,
                dp: g_res_p,
                dn: g_res_dn,
            },
            ln_fugacity: CpaFugacityDerivatives {
                ln_phi,
                dt: dlnphi_dt,
                dp: dlnphi_dp,
                dn: dlnphi_dn,
            },
            cp,
            cv,
            speed_of_sound: sound_speed,
        })
    }

    /// Evaluate fugacity coefficients and their analytical `(T, P, n)` derivatives.
    pub fn ln_fugacity_derivatives_tp(
        &self,
        temperature: f64,
        pressure: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<CpaFugacityDerivatives> {
        Ok(self
            .state_tp(temperature, pressure, moles, phase)?
            .ln_fugacity)
    }
}

fn validate_tvn(
    nc: usize,
    temperature: f64,
    volume: f64,
    moles: &[f64],
) -> Result<(f64, f64, Vec<f64>)> {
    positive(temperature, "temperature")?;
    positive(volume, "volume")?;
    if moles.len() != nc {
        return Err(Error::InvalidInput(
            "the number of moles must match the components of the model",
        ));
    }
    if moles.iter().any(|&v| !v.is_finite() || v < 0.0) {
        return Err(Error::InvalidInput(
            "mole numbers must be finite and non-negative",
        ));
    }
    let total: f64 = moles.iter().copied().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(Error::InvalidInput(
            "the number of moles must be strictly positive",
        ));
    }
    let v_unit = volume / total;
    if !v_unit.is_finite() || v_unit <= 0.0 {
        return Err(Error::InvalidInput(
            "molar volume must be finite and strictly positive",
        ));
    }
    let composition = moles.iter().map(|&v| v / total).collect();
    Ok((total, v_unit, composition))
}

fn build_sites(components: &[CpaComponent]) -> Vec<AssocSite> {
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
    sites
}

fn build_cpa_assoc(
    components: &[CpaComponent],
    eps_comb_arithmetic: &[bool],
    beta_comb_arithmetic: &[bool],
    kij_eps: &[f64],
    kij_beta: &[f64],
    delta_comb_rule: DeltaCombRule,
) -> Result<AssociationParams> {
    let schemes: Vec<AssocScheme> = components.iter().map(|c| c.scheme).collect();
    let eps: Vec<f64> = components.iter().map(|c| c.eps).collect();
    let beta: Vec<f64> = components.iter().map(|c| c.beta).collect();
    let eps_rules: Vec<AssocCombRule> = eps_comb_arithmetic
        .iter()
        .map(|&ari| {
            if ari {
                AssocCombRule::AriComb
            } else {
                AssocCombRule::GeoComb
            }
        })
        .collect();
    let beta_rules: Vec<AssocCombRule> = beta_comb_arithmetic
        .iter()
        .map(|&ari| {
            if ari {
                AssocCombRule::AriComb
            } else {
                AssocCombRule::GeoComb
            }
        })
        .collect();
    AssociationParams::for_cpa(
        &schemes,
        &eps,
        &beta,
        &eps_rules,
        &beta_rules,
        kij_eps,
        kij_beta,
        delta_comb_rule,
    )
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
        AssocScheme::Scheme4B => {
            if site <= 3 {
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

/// Hard-sphere contact RDF provider for CPA wired to pure component covolumes
/// (`saft_rdf.f90:93-100` and `saft_globals.f90:51`).
struct CpaRdfProvider<'a> {
    model: &'a Cpa,
    v_unit: f64,
    b_assoc_mix: f64,
    eta: f64,
}

impl<'a> CpaRdfProvider<'a> {
    fn new(model: &'a Cpa, v_unit: f64, composition: &[f64]) -> Result<Self> {
        let mut b_assoc_mix = 0.0;
        for (comp, &xi) in model.components.iter().zip(composition.iter()) {
            b_assoc_mix += xi * comp.b;
        }
        let eta = 0.25 * b_assoc_mix / v_unit;
        if model.simplified_rdf {
            if !eta.is_finite() || 1.0 - 1.9 * eta <= 0.0 {
                return Err(Error::InvalidInput(
                    "density exceeds simplified CPA RDF packing limit",
                ));
            }
        } else if !eta.is_finite() || 1.0 - eta <= 0.0 {
            return Err(Error::InvalidInput(
                "density exceeds original CPA RDF packing limit",
            ));
        }
        Ok(Self {
            model,
            v_unit,
            b_assoc_mix,
            eta,
        })
    }

    #[inline]
    fn covol(&self, i: usize, j: usize) -> f64 {
        0.5 * (self.model.components[i].b + self.model.components[j].b)
    }

    #[inline]
    fn g_and_eta_derivs(&self) -> (f64, f64, f64) {
        let eta = self.eta;
        if self.model.simplified_rdf {
            let d1 = 1.0 / (1.0 - 1.9 * eta);
            let d2 = d1 * d1;
            (d1, 1.9 * d2, 2.0 * 1.9 * 1.9 * d2 * d1)
        } else {
            let d1 = 1.0 / (1.0 - eta);
            let d2 = d1 * d1;
            let d3 = d2 * d1;
            let d4 = d2 * d2;
            (
                (1.0 - 0.5 * eta) * d3,
                (2.5 - eta) * d4,
                (9.0 - 3.0 * eta) * d4 * d1,
            )
        }
    }
}

impl AssociationRdfProvider for CpaRdfProvider<'_> {
    fn pair_dependent(&self) -> bool {
        false
    }

    fn rdf_value(&self, i: usize, j: usize) -> Result<AssocRdfValue> {
        let (g, _, _) = self.g_and_eta_derivs();
        Ok(AssocRdfValue {
            g,
            covol: self.covol(i, j),
        })
    }

    fn rdf_first(&self, i: usize, j: usize) -> Result<AssocRdfFirst> {
        let nc = self.model.num_components();
        let (g, g_eta, _) = self.g_and_eta_derivs();
        let inv_v = 1.0 / self.v_unit;
        let eta_v = -0.25 * self.b_assoc_mix * inv_v * inv_v;
        let g_v = g_eta * eta_v;
        let mut g_n = Vec::with_capacity(nc);
        for comp in &self.model.components {
            g_n.push(g_eta * (0.25 * comp.b * inv_v));
        }
        Ok(AssocRdfFirst {
            g,
            g_t: 0.0,
            g_v,
            g_n,
            covol: self.covol(i, j),
            covol_t: 0.0,
        })
    }

    fn rdf_second(&self, i: usize, j: usize) -> Result<AssocRdfSecond> {
        let nc = self.model.num_components();
        let (g, g_eta, g_eta2) = self.g_and_eta_derivs();
        let inv_v = 1.0 / self.v_unit;
        let inv_v2 = inv_v * inv_v;
        let inv_v3 = inv_v2 * inv_v;
        let eta_v = -0.25 * self.b_assoc_mix * inv_v2;
        let eta_vv = 0.5 * self.b_assoc_mix * inv_v3;

        let mut eta_n = Vec::with_capacity(nc);
        let mut moles = Vec::with_capacity(nc);
        let mut volume_moles = Vec::with_capacity(nc);
        for comp in &self.model.components {
            let en = 0.25 * comp.b * inv_v;
            let evn = -0.25 * comp.b * inv_v2;
            eta_n.push(en);
            moles.push(g_eta * en);
            volume_moles.push(g_eta2 * eta_v * en + g_eta * evn);
        }
        let mut moles_second = vec![0.0; nc * nc];
        for r in 0..nc {
            for c in 0..nc {
                moles_second[r * nc + c] = g_eta2 * eta_n[r] * eta_n[c];
            }
        }
        Ok(AssocRdfSecond {
            rdf: RdfDerivatives {
                g,
                temperature: 0.0,
                volume: g_eta * eta_v,
                moles,
                temperature_second: 0.0,
                temperature_volume: 0.0,
                temperature_moles: vec![0.0; nc],
                volume_second: g_eta2 * eta_v * eta_v + g_eta * eta_vv,
                volume_moles,
                moles_second,
            },
            covol: self.covol(i, j),
            covol_t: 0.0,
            covol_tt: 0.0,
        })
    }
}

/// Prepared CPA state for fast analytical evaluations at fixed T and composition.
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
    /// Prepared underlying cubic state (VdW / HV1 / HV2 / WS / HVWS).
    cubic_state: PreparedState<'a>,
    /// Pure-component linear covolume average used in the association RDF.
    b_assoc_mix: f64,
}

impl<'a> PreparedCpa<'a> {
    /// Create and prepare a CPA state at temperature T (K) and mole fractions x.
    pub fn new(model: &'a Cpa, temperature: f64, mole_fractions: &[f64]) -> Result<Self> {
        positive(temperature, "temperature")?;
        let nc = model.num_components();
        if mole_fractions.len() != nc || mole_fractions.iter().any(|&v| !v.is_finite() || v < 0.0) {
            return Err(Error::InvalidInput(
                "finite non-negative mole fractions matching number of components required",
            ));
        }

        let sum_x: f64 = mole_fractions.iter().copied().sum();
        if sum_x <= 0.0 || !sum_x.is_finite() {
            return Err(Error::InvalidInput(
                "mole fractions sum must be strictly positive",
            ));
        }
        let x: Vec<f64> = mole_fractions.iter().map(|&xi| xi / sum_x).collect();

        let cubic_state = model.cubic.prepare(temperature, &x)?;
        let (a_mix, da_mix_dt, d2a_mix_dt2) = cubic_state.mixture_a_derivatives();
        let b_mix = cubic_state.eos_covolume();

        let mut alpha_i = Vec::with_capacity(nc);
        let mut a_i = Vec::with_capacity(nc);
        let mut b_assoc_mix = 0.0;
        let mut molar_mass_mix = 0.0;
        for (i, comp) in model.components.iter().enumerate() {
            let alpha_val = model.cubic.alpha(i)?.evaluate(temperature, comp.tc)?.value;
            alpha_i.push(alpha_val);
            a_i.push(model.cubic.pure_a0()[i] * alpha_val);
            b_assoc_mix += x[i] * comp.b;
            molar_mass_mix += x[i] * comp.molar_mass;
        }

        let ns = model.assoc.num_sites;
        let mut interacts = vec![false; ns * ns];
        let mut covol_kl = vec![0.0; ns * ns];
        let mut boltzmann_factor = vec![0.0; ns * ns];
        for k in 0..ns {
            let ic = model.assoc.site_to_comp[k];
            for l in 0..ns {
                let jc = model.assoc.site_to_comp[l];
                let idx = k * ns + l;
                if model.assoc.beta_kl[idx] > 0.0 {
                    interacts[idx] = true;
                    covol_kl[idx] = 0.5 * (model.components[ic].b + model.components[jc].b);
                    let arg = model.assoc.eps_kl[idx] / (RGAS * temperature);
                    boltzmann_factor[idx] = arg.exp_m1();
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
            eps_kl: model.assoc.eps_kl.clone(),
            beta_kl: model.assoc.beta_kl.clone(),
            covol_kl,
            boltzmann_factor,
            cubic_state,
            b_assoc_mix,
        })
    }

    #[inline]
    fn min_molar_volume(&self) -> f64 {
        let cubic_min = self.cubic_state.covolume();
        let rdf_min = if self.model.simplified_rdf {
            0.475 * self.b_assoc_mix
        } else {
            0.25 * self.b_assoc_mix
        };
        cubic_min.max(rdf_min)
    }

    /// Hard-sphere contact radial distribution function `g(eta)` and its
    /// first and second derivatives with respect to `eta = 0.25 * b_assoc * rho`.
    #[inline]
    pub fn g_rdf(&self, rho: f64) -> (f64, f64, f64) {
        let eta = 0.25 * self.b_assoc_mix * rho;
        if self.model.simplified_rdf {
            let den = 1.0 - 1.9 * eta;
            let d1 = 1.0 / den;
            let d2 = d1 * d1;
            (d1, 1.9 * d2, 2.0 * 1.9 * 1.9 * (d2 * d1))
        } else {
            let om_eta = 1.0 - eta;
            let d1 = 1.0 / om_eta;
            let d2 = d1 * d1;
            let d3 = d2 * d1;
            let d4 = d2 * d2;
            (
                (1.0 - 0.5 * eta) * d3,
                (2.5 - eta) * d4,
                (9.0 - 3.0 * eta) * crate::numerics::ipow(d1, 5),
            )
        }
    }

    /// Solve for non-bonded site fractions `X_k` at density `rho` [mol / m^3].
    pub fn solve_site_fractions(&self, rho: f64) -> Result<Vec<f64>> {
        positive(rho, "density")?;
        if !self.model.assoc.has_sites() {
            return Ok(Vec::new());
        }
        let v_unit = 1.0 / rho;
        let rdf_provider = CpaRdfProvider::new(self.model, v_unit, &self.mole_fractions)?;
        association::solve_for_x_k(
            &self.model.assoc,
            &rdf_provider,
            self.temperature,
            v_unit,
            &self.mole_fractions,
            None,
        )
    }

    /// Evaluate the 10-block reduced residual Helmholtz energy on a unit-mole
    /// basis (`N = 1`, molar volume `v_unit [m^3/mol]`).
    pub fn residual_helmholtz_tv_unit(&self, v_unit: f64) -> Result<CpaHelmholtz> {
        let cubic = self.cubic_state.residual_helmholtz_tv_unit(v_unit)?;
        if !self.model.assoc.has_sites() {
            return Ok(CpaHelmholtz {
                f: cubic.f,
                f_t: cubic.f_t,
                f_v: cubic.f_v,
                f_n: cubic.f_n,
                f_tt: cubic.f_tt,
                f_tv: cubic.f_tv,
                f_tn: cubic.f_tn,
                f_vv: cubic.f_vv,
                f_vn: cubic.f_vn,
                f_nn: cubic.f_nn,
                x_k: Vec::new(),
            });
        }
        let rdf_provider = CpaRdfProvider::new(self.model, v_unit, &self.mole_fractions)?;
        let assoc = association::assoc_helmholtz_all(
            &self.model.assoc,
            &rdf_provider,
            self.temperature,
            v_unit,
            &self.mole_fractions,
        )?;
        let nc = self.model.num_components();
        let mut f_n = cubic.f_n;
        let mut f_tn = cubic.f_tn;
        let mut f_vn = cubic.f_vn;
        let mut f_nn = cubic.f_nn;
        for i in 0..nc {
            f_n[i] += assoc.f_n[i];
            f_tn[i] += assoc.f_tn[i];
            f_vn[i] += assoc.f_vn[i];
        }
        for k in 0..(nc * nc) {
            f_nn[k] += assoc.f_nn[k];
        }
        Ok(CpaHelmholtz {
            f: cubic.f + assoc.f,
            f_t: cubic.f_t + assoc.f_t,
            f_v: cubic.f_v + assoc.f_v,
            f_n,
            f_tt: cubic.f_tt + assoc.f_tt,
            f_tv: cubic.f_tv + assoc.f_tv,
            f_tn,
            f_vv: cubic.f_vv + assoc.f_vv,
            f_vn,
            f_nn,
            x_k: assoc.x_k,
        })
    }

    /// Residual Helmholtz energy terms `(a_cubic, a_assoc, a_res)` at molar
    /// density `rho` [mol / m^3].
    pub fn helmholtz_terms(&self, rho: f64) -> Result<(f64, f64, f64)> {
        positive(rho, "density")?;
        let v_unit = 1.0 / rho;
        let cubic = self.cubic_state.residual_helmholtz_tv_unit(v_unit)?;
        let a_assoc = if self.model.assoc.has_sites() {
            let rdf_provider = CpaRdfProvider::new(self.model, v_unit, &self.mole_fractions)?;
            let (f_assoc, _) = association::assoc_helmholtz_value(
                &self.model.assoc,
                &rdf_provider,
                self.temperature,
                v_unit,
                &self.mole_fractions,
            )?;
            f_assoc
        } else {
            0.0
        };
        Ok((cubic.f, a_assoc, cubic.f + a_assoc))
    }

    /// Exact analytical pressure and derivatives at density `rho` [mol / m^3].
    pub fn pressure_and_derivatives(&self, rho: f64) -> Result<CpaPressureJet> {
        positive(rho, "density")?;
        let v_unit = 1.0 / rho;
        let t = self.temperature;
        let rt = RGAS * t;

        let (mut f, mut f_v, mut f_vv, mut f_tv) = self
            .cubic_state
            .residual_helmholtz_pressure_jet_unit(v_unit)?;

        if self.model.assoc.has_sites() {
            let rdf_provider = CpaRdfProvider::new(self.model, v_unit, &self.mole_fractions)?;
            let assoc = association::assoc_helmholtz_pressure_jet(
                &self.model.assoc,
                &rdf_provider,
                t,
                v_unit,
                &self.mole_fractions,
            )?;
            f += assoc.f;
            f_v += assoc.f_v;
            f_vv += assoc.f_vv;
            f_tv += assoc.f_tv;
        }

        let p = rho * rt - rt * f_v;
        let z = p / (rho * rt);
        let dp_dv = -rho * rho * rt - rt * f_vv;
        let dp_drho = -dp_dv / (rho * rho);
        let dp_dt = rho * RGAS - RGAS * f_v - rt * f_tv;

        Ok(CpaPressureJet {
            pressure: p,
            compressibility: z,
            residual_helmholtz: f,
            dp_drho,
            dp_dv,
            dp_dt,
        })
    }

    /// Exact analytical fugacity coefficients `ln(phi_i)` at density `rho` [mol / m^3].
    pub fn fugacity_coefficients(&self, rho: f64, z: f64) -> Result<Vec<f64>> {
        positive(rho, "density")?;
        positive(z, "compressibility factor")?;
        let v_unit = 1.0 / rho;
        let cubic = self.cubic_state.residual_helmholtz_tv_unit(v_unit)?;
        let mut f_n = cubic.f_n;
        if self.model.assoc.has_sites() {
            let rdf_provider = CpaRdfProvider::new(self.model, v_unit, &self.mole_fractions)?;
            let assoc = association::assoc_helmholtz_first(
                &self.model.assoc,
                &rdf_provider,
                self.temperature,
                v_unit,
                &self.mole_fractions,
            )?;
            for (dst, src) in f_n.iter_mut().zip(assoc.f_n.iter()) {
                *dst += *src;
            }
        }
        let ln_z = z.ln();
        Ok(f_n.into_iter().map(|fni| fni - ln_z).collect())
    }

    /// Complete analytical thermodynamic state properties at density `rho` [mol / m^3].
    pub fn state(&self, rho: f64) -> Result<CpaState> {
        positive(rho, "density")?;
        let v_unit = 1.0 / rho;
        let t = self.temperature;
        let rt = RGAS * t;
        let h = self.residual_helmholtz_tv_unit(v_unit)?;

        let p = rho * rt - rt * h.f_v;
        let z = p / (rho * rt);
        let dp_dv = -rho * rho * rt - rt * h.f_vv;
        let dp_drho = -dp_dv / (rho * rho);
        let dp_dt = rho * RGAS - RGAS * h.f_v - rt * h.f_tv;

        let u_res = -RGAS * t * t * h.f_t;
        let h_res = u_res + rt * (z - 1.0);
        let ln_z = if z > 0.0 { z.ln() } else { f64::NAN };
        let s_res = -RGAS * (t * h.f_t + h.f) + RGAS * ln_z;
        let a_res_energy = rt * h.f;
        let g_res = h_res - t * s_res;
        let cv_res = -RGAS * (2.0 * t * h.f_t + t * t * h.f_tt);

        let mut cp_id = 0.0;
        let mut h_id = 0.0;
        let mut s_id_tv = 0.0;
        let p_ref = 100_000.0;
        for i in 0..self.model.components.len() {
            let id_props = self.model.ideal[i].properties(t, p_ref)?;
            let xi = self.mole_fractions[i];
            cp_id += xi * id_props.cp;
            h_id += xi * id_props.enthalpy;
            let mix_term = if xi > 0.0 { xi * xi.ln() } else { 0.0 };
            s_id_tv += xi * id_props.entropy - RGAS * mix_term;
        }
        s_id_tv -= RGAS * (rho * rt / p_ref).ln();
        // Total entropy: S_total = S_ideal(T, V) + S_res_TV = S_ideal(T, P) + S_res_TP
        let s_total = s_id_tv - RGAS * (t * h.f_t + h.f);

        let cv_total = cp_id - RGAS + cv_res;
        let cp_total = cv_total - t * (dp_dt * dp_dt) / dp_dv;
        let gamma = cp_total / cv_total;
        let mw_kg = self.molar_mass_mix * 1.0e-3;
        let sound_speed = if dp_drho > 0.0 && gamma > 0.0 && mw_kg > 0.0 {
            (gamma * dp_drho / mw_kg).sqrt()
        } else {
            0.0
        };

        let ln_phi = h.f_n.iter().map(|&fni| fni - ln_z).collect();

        Ok(CpaState {
            pressure: p,
            compressibility: z,
            residual_helmholtz: h.f,
            u_residual: u_res,
            h_residual: h_res,
            s_residual: s_res,
            a_residual: a_res_energy,
            g_residual: g_res,
            h_total: h_id + h_res,
            s_total,
            cv_residual: cv_res,
            cv_total,
            cp_total,
            sound_speed,
            ln_phi,
        })
    }

    /// Single-phase volume root for a specified pressure and phase.
    pub fn volume(&self, pressure_target: f64, phase: Phase) -> Result<f64> {
        positive(pressure_target, "pressure target")?;
        let v_pole = self.min_molar_volume();
        let v_min = v_pole * (1.0 + 1.0e-7);

        let solve_from = |mut v: f64, is_liquid: bool| -> Result<f64> {
            let tol = 1.0e-12 * pressure_target + 1.0e-10;
            let mut jet = self.pressure_and_derivatives(1.0 / v)?;
            for _ in 0..100 {
                let diff = jet.pressure - pressure_target;
                if diff.abs() <= tol {
                    return Ok(v);
                }
                if !jet.dp_dv.is_finite() || jet.dp_dv == 0.0 {
                    break;
                }
                let dv = -diff / jet.dp_dv;
                if dv.abs() <= 1.0e-14 * v && diff.abs() <= 1.0e-9 * pressure_target + 1.0e-6 {
                    return Ok(v);
                }
                let mut step = if is_liquid {
                    dv.clamp(-0.35 * (v - v_pole), 0.5 * v)
                } else {
                    dv.clamp(-0.5 * v, 0.8 * v)
                };
                if v + step <= v_min {
                    step = 0.5 * (v_min - v);
                }
                let cand = (v + step).max(v_min);
                let cand_jet = self.pressure_and_derivatives(1.0 / cand)?;
                if (cand_jet.pressure - pressure_target).abs() > diff.abs()
                    && (jet.dp_dv > 0.0 || cand_jet.dp_dv > 0.0)
                {
                    // Push out of spinodal region towards the requested branch
                    v = if is_liquid {
                        (v_min + 0.5 * (v - v_min)).max(v_min)
                    } else {
                        (v * 1.5).max(v_min)
                    };
                    jet = self.pressure_and_derivatives(1.0 / v)?;
                } else {
                    v = cand;
                    jet = cand_jet;
                }
            }
            if (jet.pressure - pressure_target).abs() <= 1.0e-7 * pressure_target + 1.0e-2 {
                Ok(v)
            } else {
                Err(Error::NotConverged("CPA volume root"))
            }
        };

        match phase {
            Phase::Vapor => {
                let v0 = (RGAS * self.temperature / pressure_target).max(v_min * 2.0);
                solve_from(v0, false)
            }
            Phase::Liquid => {
                let v0 = (v_pole * 1.12).max(v_min * 1.01);
                solve_from(v0, true).or_else(|_| {
                    let v_vap = (RGAS * self.temperature / pressure_target).max(v_min * 2.0);
                    solve_from(v_vap, false)
                })
            }
            Phase::MinimumGibbs => {
                let vl = self.volume(pressure_target, Phase::Liquid)?;
                let vv = self.volume(pressure_target, Phase::Vapor)?;
                if (vl - vv).abs() < 1.0e-7 * vl {
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

    /// Pure fluid saturation pressure and coexistence volumes at temperature T (K).
    pub fn saturation_pressure(&self, temperature: f64) -> Result<CpaSaturation> {
        if self.model.num_components() != 1 {
            return Err(Error::InvalidInput(
                "saturation calculation requires pure fluid",
            ));
        }

        let mut p = 101_325.0;

        for _ in 0..60 {
            let vl = self.volume(p, Phase::Liquid)?;
            let vv = self.volume(p, Phase::Vapor)?;

            if (vl - vv).abs() < 1.0e-6 * vl {
                return Err(Error::NotConverged(
                    "liquid and vapor roots converged to identical state",
                ));
            }

            let jet_l = self.pressure_and_derivatives(1.0 / vl)?;
            let jet_v = self.pressure_and_derivatives(1.0 / vv)?;

            let phi_l = self.fugacity_coefficients(1.0 / vl, jet_l.compressibility)?[0];
            let phi_v = self.fugacity_coefficients(1.0 / vv, jet_v.compressibility)?[0];

            let delta_phi = phi_l - phi_v;
            if delta_phi.abs() < 1.0e-10 {
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
            if p <= 1.0e-2 {
                p = 1.0e-2;
            }
        }

        Err(Error::NotConverged(
            "saturation pressure failed to converge",
        ))
    }
}
