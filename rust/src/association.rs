//! Shared SAFT / CPA association engine ported from ThermoPack.
//!
//! Mirrors `src/assocschemeutils.f90`, `src/association_var.f90`, and
//! `src/saft_association.f90` in pure SI units (`K`, `m^3`, `mol`, `J`) with
//! Euler-homogeneous unit-mole evaluation so tiny total mole numbers never
//! suffer from `1/N^2` or `1/N^3` overflow or unscaled Newton tolerances.
//!
//! Designed for reuse across PC-SAFT / sPC-SAFT / PCP-SAFT / sPCP-SAFT,
//! SAFT-VR Mie, and CPA.

use crate::error::positive;
use crate::rdf::RdfDerivatives;
use crate::{Error, Result};

/// Universal gas constant [J/(mol K)] matching ThermoPack's `Rgas_default`.
const RGAS: f64 = crate::R;

/// Association scheme classification matching `src/assocschemeutils.f90`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssocScheme {
    /// No association sites (`no_assoc`).
    NoAssoc,
    /// Single self-associating site (`1`: AA).
    Scheme1,
    /// Two sites, all pairs interact (`2A`: AA, AB, BB).
    Scheme2A,
    /// Two oppositely charged sites (`2B`: AB).
    Scheme2B,
    /// Two sites, site 1 bipolar (`2C`: AA, AB).
    Scheme2C,
    /// Three sites, all pairs interact (`3A`).
    Scheme3A,
    /// Three sites: two positive, one negative (`3B`: AC, BC).
    Scheme3B,
    /// Four sites, all pairs interact (`4A`).
    Scheme4A,
    /// Four sites: three positive, one negative (`4B`: AD, BD, CD).
    Scheme4B,
    /// Four sites: two negative, two positive (`4C`: AC, AD, BC, BD).
    Scheme4C,
    /// Single electron-acceptor site with positive polarity (`1ea`, induced cross-association only).
    Scheme1ea,
}

impl AssocScheme {
    /// Number of association sites on a single molecule.
    pub const fn num_sites(self) -> usize {
        match self {
            AssocScheme::NoAssoc => 0,
            AssocScheme::Scheme1 | AssocScheme::Scheme1ea => 1,
            AssocScheme::Scheme2A | AssocScheme::Scheme2B | AssocScheme::Scheme2C => 2,
            AssocScheme::Scheme3A | AssocScheme::Scheme3B => 3,
            AssocScheme::Scheme4A | AssocScheme::Scheme4B | AssocScheme::Scheme4C => 4,
        }
    }

    /// String label matching `AssocSchemeUtils::get_assoc_string`.
    pub const fn as_str(self) -> &'static str {
        match self {
            AssocScheme::NoAssoc => "None",
            AssocScheme::Scheme1 => "1",
            AssocScheme::Scheme2A => "2A",
            AssocScheme::Scheme2B => "2B",
            AssocScheme::Scheme2C => "2C",
            AssocScheme::Scheme3A => "3A",
            AssocScheme::Scheme3B => "3B",
            AssocScheme::Scheme4A => "4A",
            AssocScheme::Scheme4B => "4B",
            AssocScheme::Scheme4C => "4C",
            AssocScheme::Scheme1ea => "1ea",
        }
    }

    /// Whether two 1-based sites `(site1, site2)` on like molecules interact.
    ///
    /// Mirrors `AssocSchemeUtils::site_interaction_internal`.
    pub const fn site_interaction_internal(self, site1: usize, site2: usize) -> bool {
        match self {
            AssocScheme::Scheme1
            | AssocScheme::Scheme2A
            | AssocScheme::Scheme3A
            | AssocScheme::Scheme4A => true,
            AssocScheme::Scheme2B => site1 != site2,
            AssocScheme::Scheme2C => site1 == 1 || site2 == 1,
            AssocScheme::Scheme3B => (site1 == 3) != (site2 == 3),
            AssocScheme::Scheme4B => (site1 == 4) != (site2 == 4),
            AssocScheme::Scheme4C => (site1 <= 2 && site2 >= 3) || (site2 <= 2 && site1 >= 3),
            AssocScheme::NoAssoc | AssocScheme::Scheme1ea => false,
        }
    }

    /// Polarity of a 1-based site (`-1`, `0`, or `+1`).
    ///
    /// Mirrors `AssocSchemeUtils::polarity`.
    pub const fn polarity(self, site: usize) -> i8 {
        match self {
            AssocScheme::Scheme1 => -1,
            AssocScheme::Scheme2B => match site {
                1 => -1,
                2 => 1,
                _ => 0,
            },
            AssocScheme::Scheme2C => match site {
                1 => 0,
                2 => -1,
                _ => 0,
            },
            AssocScheme::Scheme3B => {
                if site <= 2 {
                    1
                } else if site == 3 {
                    -1
                } else {
                    0
                }
            }
            AssocScheme::Scheme4C => {
                if site <= 2 {
                    -1
                } else if site >= 3 {
                    1
                } else {
                    0
                }
            }
            AssocScheme::Scheme1ea => 1,
            AssocScheme::NoAssoc
            | AssocScheme::Scheme2A
            | AssocScheme::Scheme3A
            | AssocScheme::Scheme4A
            | AssocScheme::Scheme4B => 0,
        }
    }

    /// Whether 1-based `site1` on a molecule with scheme `self` interacts with
    /// 1-based `site2` on an unlike molecule with scheme `other`.
    ///
    /// Mirrors `AssocSchemeUtils::cross_site_interaction`.
    pub const fn cross_site_interaction(
        self,
        site1: usize,
        other: AssocScheme,
        site2: usize,
    ) -> bool {
        if matches!(self, AssocScheme::NoAssoc) || matches!(other, AssocScheme::NoAssoc) {
            return false;
        }
        let pol1 = self.polarity(site1);
        let pol2 = other.polarity(site2);
        (pol1 as i16) * (pol2 as i16) <= 0
    }
}

/// Combining rule for binary association parameters (`eps`, `beta`).
///
/// Mirrors `geoComb`, `ariComb`, `defaultComb` in `src/assocschemeutils.f90`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AssocCombRule {
    /// Geometric mean (`sqrt(v1 * v2)`).
    GeoComb,
    /// Arithmetic mean (`0.5 * (v1 + v2)`).
    AriComb,
    /// EOS-specific default combining rule (e.g. PC-SAFT CR-1 / Wolbach-Sandler).
    #[default]
    DefaultComb,
}

impl AssocCombRule {
    /// Apply explicit arithmetic or geometric combining rule without intermediate
    /// overflow.
    pub fn apply(self, val1: f64, val2: f64) -> Result<f64> {
        if !val1.is_finite() || !val2.is_finite() {
            return Err(Error::InvalidInput(
                "non-finite parameter in association combining rule",
            ));
        }
        match self {
            AssocCombRule::AriComb => {
                let out = 0.5 * val1 + 0.5 * val2;
                if !out.is_finite() {
                    return Err(Error::InvalidInput(
                        "overflow in arithmetic association combining rule",
                    ));
                }
                Ok(out)
            }
            AssocCombRule::GeoComb => {
                if val1 < 0.0 || val2 < 0.0 {
                    return Err(Error::InvalidInput(
                        "negative parameter in geometric association combining rule",
                    ));
                }
                let out = val1.sqrt() * val2.sqrt();
                if !out.is_finite() {
                    return Err(Error::InvalidInput(
                        "overflow in geometric association combining rule",
                    ));
                }
                Ok(out)
            }
            AssocCombRule::DefaultComb => Err(Error::InvalidInput(
                "cannot apply DefaultComb directly without EOS-specific rule",
            )),
        }
    }
}

/// Cross-association strength combining rule (`STANDARD` vs `ELLIOT`).
///
/// Mirrors `association_var::STANDARD` and `association_var::ELLIOT`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DeltaCombRule {
    /// Standard combining of `eps_kl` and `beta_kl`.
    #[default]
    Standard,
    /// Elliot combining rule: `Delta_{kl} = sqrt(Delta_{k1,k2} * Delta_{l1,l2})`.
    Elliot,
}

/// Site-level association bookkeeping and interaction matrices.
///
/// Mirrors `association_var::association`.
#[derive(Debug, Clone)]
pub struct AssociationParams {
    pub num_components: usize,
    pub num_sites: usize,
    /// 0-based indices of components that carry at least one association site.
    pub comp_indices: Vec<usize>,
    /// Inclusive 0-based site range `(first, last)` per component, or `None`.
    pub comp_vs_sites: Vec<Option<(usize, usize)>>,
    /// 0-based component index owning each 0-based site.
    pub site_to_comp: Vec<usize>,
    /// Association volume matrix `beta_kl` (`num_sites * num_sites`, row-major).
    pub beta_kl: Vec<f64>,
    /// Association energy matrix `eps_kl` (`num_sites * num_sites`, row-major) `[J/mol]`.
    pub eps_kl: Vec<f64>,
    /// Cross-association `Delta` combining rule (`Standard` or `Elliot`).
    pub delta_comb_rule: DeltaCombRule,
}

impl AssociationParams {
    /// Construct empty (non-associating) parameters for `num_components`.
    pub fn empty(num_components: usize) -> Self {
        Self {
            num_components,
            num_sites: 0,
            comp_indices: Vec::new(),
            comp_vs_sites: vec![None; num_components],
            site_to_comp: Vec::new(),
            beta_kl: Vec::new(),
            eps_kl: Vec::new(),
            delta_comb_rule: DeltaCombRule::Standard,
        }
    }

    /// Whether the mixture has any association sites.
    #[inline]
    pub fn has_sites(&self) -> bool {
        self.num_sites > 0
    }

    #[allow(clippy::type_complexity)]
    fn build_bookkeeping(
        schemes: &[AssocScheme],
    ) -> (usize, Vec<usize>, Vec<Option<(usize, usize)>>, Vec<usize>) {
        let nc = schemes.len();
        let mut comp_indices = Vec::new();
        let mut comp_vs_sites = vec![None; nc];
        let mut site_to_comp = Vec::new();
        let mut num_sites = 0usize;

        for (ic, &scheme) in schemes.iter().enumerate() {
            let ns = scheme.num_sites();
            if ns > 0 {
                comp_indices.push(ic);
                let first = num_sites;
                let last = num_sites + ns - 1;
                comp_vs_sites[ic] = Some((first, last));
                for _ in 0..ns {
                    site_to_comp.push(ic);
                }
                num_sites += ns;
            }
        }
        (num_sites, comp_indices, comp_vs_sites, site_to_comp)
    }

    fn validate_params(&self) -> Result<()> {
        if self.eps_kl.iter().any(|&v| !v.is_finite() || v < 0.0) {
            return Err(Error::InvalidInput(
                "association energy parameters must be finite and non-negative",
            ));
        }
        if self.beta_kl.iter().any(|&v| !v.is_finite() || v < 0.0) {
            return Err(Error::InvalidInput(
                "association volume parameters must be finite and non-negative",
            ));
        }
        if self.delta_comb_rule == DeltaCombRule::Elliot {
            let nc = self.num_components;
            let ns = self.num_sites;
            for ic in 0..nc {
                let Some((k1, k2)) = self.comp_vs_sites[ic] else {
                    continue;
                };
                for jc in (ic + 1)..nc {
                    let Some((l1, l2)) = self.comp_vs_sites[jc] else {
                        continue;
                    };
                    if self.beta_kl[k1 * ns + k2] <= 0.0 || self.beta_kl[l1 * ns + l2] <= 0.0 {
                        return Err(Error::Unsupported(
                            "Elliot combining rule requires nonzero self-association beta",
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    /// Build association parameters for PC-SAFT / sPC-SAFT / PCP-SAFT / sPCP-SAFT.
    ///
    /// Mirrors `saft_interface::saft_setAssocParams` (`src/saft_interface.f90:273-336`).
    pub fn for_pc_saft(
        schemes: &[AssocScheme],
        eps: &[f64],
        beta: &[f64],
        sigma: &[f64],
        eps_comb_rules: &[AssocCombRule],
        beta_comb_rules: &[AssocCombRule],
        delta_comb_rule: DeltaCombRule,
    ) -> Result<Self> {
        let nc = schemes.len();
        if eps.len() != nc
            || beta.len() != nc
            || sigma.len() != nc
            || eps_comb_rules.len() != nc * nc
            || beta_comb_rules.len() != nc * nc
        {
            return Err(Error::InvalidInput(
                "association parameter dimension mismatch",
            ));
        }
        if sigma.iter().any(|&s| !s.is_finite() || s <= 0.0) {
            return Err(Error::InvalidInput(
                "segment diameters must be finite and strictly positive",
            ));
        }
        if eps.iter().any(|&v| !v.is_finite() || v < 0.0)
            || beta.iter().any(|&v| !v.is_finite() || v < 0.0)
        {
            return Err(Error::InvalidInput(
                "association parameters must be finite and non-negative",
            ));
        }

        let (num_sites, comp_indices, comp_vs_sites, site_to_comp) =
            Self::build_bookkeeping(schemes);
        let mut eps_kl = vec![0.0; num_sites * num_sites];
        let mut beta_kl = vec![0.0; num_sites * num_sites];

        for (k, &ic) in site_to_comp.iter().enumerate() {
            let (k_first, _) = comp_vs_sites[ic].expect("site belongs to component");
            let site_k = k - k_first + 1;
            for (l, &jc) in site_to_comp.iter().enumerate() {
                let (l_first, _) = comp_vs_sites[jc].expect("site belongs to component");
                let site_l = l - l_first + 1;
                let idx = k * num_sites + l;
                if ic == jc {
                    if schemes[ic].site_interaction_internal(site_k, site_l) {
                        eps_kl[idx] = eps[ic];
                        beta_kl[idx] = beta[ic];
                    }
                } else if schemes[ic].cross_site_interaction(site_k, schemes[jc], site_l) {
                    let pair = ic * nc + jc;
                    let eps_rule = eps_comb_rules[pair];
                    let beta_rule = beta_comb_rules[pair];
                    eps_kl[idx] = match eps_rule {
                        AssocCombRule::DefaultComb => 0.5 * eps[ic] + 0.5 * eps[jc],
                        other => other.apply(eps[ic], eps[jc])?,
                    };
                    beta_kl[idx] = match beta_rule {
                        AssocCombRule::DefaultComb => {
                            let ratio = (sigma[ic].sqrt() * sigma[jc].sqrt())
                                / (0.5 * sigma[ic] + 0.5 * sigma[jc]);
                            (beta[ic].sqrt() * beta[jc].sqrt()) * (ratio * ratio * ratio)
                        }
                        other => other.apply(beta[ic], beta[jc])?,
                    };
                }
            }
        }

        let out = Self {
            num_components: nc,
            num_sites,
            comp_indices,
            comp_vs_sites,
            site_to_comp,
            beta_kl,
            eps_kl,
            delta_comb_rule,
        };
        out.validate_params()?;
        Ok(out)
    }

    /// Build association parameters for SAFT-VR Mie (`saft_interface.f90:314-316`).
    pub fn for_saft_vr_mie(
        schemes: &[AssocScheme],
        eps: &[f64],
        beta: &[f64],
        delta_comb_rule: DeltaCombRule,
    ) -> Result<Self> {
        let nc = schemes.len();
        if eps.len() != nc || beta.len() != nc {
            return Err(Error::InvalidInput(
                "association parameter dimension mismatch",
            ));
        }
        if eps.iter().any(|&v| !v.is_finite() || v < 0.0)
            || beta.iter().any(|&v| !v.is_finite() || v < 0.0)
        {
            return Err(Error::InvalidInput(
                "association parameters must be finite and non-negative",
            ));
        }

        let (num_sites, comp_indices, comp_vs_sites, site_to_comp) =
            Self::build_bookkeeping(schemes);
        let mut eps_kl = vec![0.0; num_sites * num_sites];
        let mut beta_kl = vec![0.0; num_sites * num_sites];

        for (k, &ic) in site_to_comp.iter().enumerate() {
            let (k_first, _) = comp_vs_sites[ic].expect("site belongs to component");
            let site_k = k - k_first + 1;
            for (l, &jc) in site_to_comp.iter().enumerate() {
                let (l_first, _) = comp_vs_sites[jc].expect("site belongs to component");
                let site_l = l - l_first + 1;
                let idx = k * num_sites + l;
                if ic == jc {
                    if schemes[ic].site_interaction_internal(site_k, site_l) {
                        eps_kl[idx] = eps[ic];
                        beta_kl[idx] = beta[ic];
                    }
                } else if schemes[ic].cross_site_interaction(site_k, schemes[jc], site_l) {
                    eps_kl[idx] = eps[ic].sqrt() * eps[jc].sqrt();
                    let avg = 0.5 * beta[ic].cbrt() + 0.5 * beta[jc].cbrt();
                    beta_kl[idx] = avg * avg * avg;
                }
            }
        }

        let out = Self {
            num_components: nc,
            num_sites,
            comp_indices,
            comp_vs_sites,
            site_to_comp,
            beta_kl,
            eps_kl,
            delta_comb_rule,
        };
        out.validate_params()?;
        Ok(out)
    }

    /// Build association parameters for CPA (`saft_interface.f90:338-389`).
    #[allow(clippy::too_many_arguments)]
    pub fn for_cpa(
        schemes: &[AssocScheme],
        eps: &[f64],
        beta: &[f64],
        eps_comb_rules: &[AssocCombRule],
        beta_comb_rules: &[AssocCombRule],
        kij_eps: &[f64],
        kij_beta: &[f64],
        delta_comb_rule: DeltaCombRule,
    ) -> Result<Self> {
        let nc = schemes.len();
        if eps.len() != nc
            || beta.len() != nc
            || eps_comb_rules.len() != nc * nc
            || beta_comb_rules.len() != nc * nc
            || kij_eps.len() != nc * nc
            || kij_beta.len() != nc * nc
        {
            return Err(Error::InvalidInput(
                "association parameter dimension mismatch",
            ));
        }
        if eps.iter().any(|&v| !v.is_finite() || v < 0.0)
            || beta.iter().any(|&v| !v.is_finite() || v < 0.0)
            || kij_eps.iter().any(|&v| !v.is_finite())
            || kij_beta.iter().any(|&v| !v.is_finite())
        {
            return Err(Error::InvalidInput(
                "association parameters must be finite and non-negative",
            ));
        }

        let (num_sites, comp_indices, comp_vs_sites, site_to_comp) =
            Self::build_bookkeeping(schemes);
        let mut eps_kl = vec![0.0; num_sites * num_sites];
        let mut beta_kl = vec![0.0; num_sites * num_sites];

        for (k, &ic) in site_to_comp.iter().enumerate() {
            let (k_first, _) = comp_vs_sites[ic].expect("site belongs to component");
            let site_k = k - k_first + 1;
            for (l, &jc) in site_to_comp.iter().enumerate() {
                let (l_first, _) = comp_vs_sites[jc].expect("site belongs to component");
                let site_l = l - l_first + 1;
                let idx = k * num_sites + l;
                if ic == jc {
                    if schemes[ic].site_interaction_internal(site_k, site_l) {
                        eps_kl[idx] = eps[ic];
                        beta_kl[idx] = beta[ic];
                    }
                } else if schemes[ic].cross_site_interaction(site_k, schemes[jc], site_l) {
                    let pair = ic * nc + jc;
                    let eps_rule = match eps_comb_rules[pair] {
                        AssocCombRule::DefaultComb => AssocCombRule::AriComb,
                        r => r,
                    };
                    let beta_rule = match beta_comb_rules[pair] {
                        AssocCombRule::DefaultComb => AssocCombRule::GeoComb,
                        r => r,
                    };
                    eps_kl[idx] = eps_rule.apply(eps[ic], eps[jc])? * (1.0 - kij_eps[pair]);
                    beta_kl[idx] = beta_rule.apply(beta[ic], beta[jc])? * (1.0 - kij_beta[pair]);
                }
            }
        }

        let out = Self {
            num_components: nc,
            num_sites,
            comp_indices,
            comp_vs_sites,
            site_to_comp,
            beta_kl,
            eps_kl,
            delta_comb_rule,
        };
        out.validate_params()?;
        Ok(out)
    }
}

/// Scalar contact RDF and binary association covolume for pair `(i, j)`.
#[derive(Debug, Clone, Copy)]
pub struct AssocRdfValue {
    pub g: f64,
    pub covol: f64,
}

/// First-order contact RDF derivatives and binary association covolume for pair `(i, j)`.
#[derive(Debug, Clone)]
pub struct AssocRdfFirst {
    pub g: f64,
    pub g_t: f64,
    pub g_v: f64,
    pub g_n: Vec<f64>,
    pub covol: f64,
    pub covol_t: f64,
}

/// Second-order contact RDF derivatives and binary association covolume for pair `(i, j)`.
#[derive(Debug, Clone)]
pub struct AssocRdfSecond {
    pub rdf: RdfDerivatives,
    pub covol: f64,
    pub covol_t: f64,
    pub covol_tt: f64,
}

/// EOS-specific provider of contact radial distribution function `g_{ij}(T, V, n)`
/// and association covolume `covol_{ij}(T)` (`[m^3/mol]`) for the shared
/// association engine.
pub trait AssociationRdfProvider {
    /// Whether `g_{ij}` depends on component indices `(i, j)`.
    ///
    /// `true` for original PC-SAFT and SAFT-VR Mie; `false` for sPC-SAFT and CPA.
    fn pair_dependent(&self) -> bool;

    /// Value-only contact RDF and binary covolume at `(T, V, n)`.
    fn rdf_value(&self, i: usize, j: usize) -> Result<AssocRdfValue>;

    /// First-derivative contact RDF and binary covolume at `(T, V, n)`.
    fn rdf_first(&self, i: usize, j: usize) -> Result<AssocRdfFirst>;

    /// Full first- and second-derivative contact RDF and binary covolume at `(T, V, n)`.
    fn rdf_second(&self, i: usize, j: usize) -> Result<AssocRdfSecond>;
}

/// First-derivative association reduced residual Helmholtz bundle (`F = A_assoc / (R T)`).
#[derive(Debug, Clone)]
pub struct AssociationFirst {
    pub f: f64,
    pub f_t: f64,
    pub f_v: f64,
    pub f_n: Vec<f64>,
    pub x_k: Vec<f64>,
}

/// Temperature- and volume-derivative association bundle (for pressure / `Cp` / `Cv` / sound speed
/// without computing `n^2` mole Hessians).
#[derive(Debug, Clone)]
pub struct AssociationPressureJet {
    pub f: f64,
    pub f_t: f64,
    pub f_v: f64,
    pub f_tt: f64,
    pub f_tv: f64,
    pub f_vv: f64,
    pub x_k: Vec<f64>,
}

/// Complete 10-block association reduced residual Helmholtz energy (`F = A_assoc / (R T)`)
/// and solved non-bonded site fractions `X_k`.
#[derive(Debug, Clone)]
pub struct AssociationHelmholtz {
    pub f: f64,
    pub f_t: f64,
    pub f_v: f64,
    pub f_n: Vec<f64>,
    pub f_tt: f64,
    pub f_tv: f64,
    pub f_tn: Vec<f64>,
    pub f_vv: f64,
    pub f_vn: Vec<f64>,
    /// Row-major `nc * nc` second mole derivative matrix `F_nn[i * nc + j]`.
    pub f_nn: Vec<f64>,
    /// Non-bonded site fractions `X_k` in `(0, 1]`.
    pub x_k: Vec<f64>,
}

/// Direct derivatives of Michelsen's `Q(T, V, n, X)` holding `X` fixed.
#[derive(Debug, Clone)]
pub struct AssociationQDerivatives {
    pub q: f64,
    pub q_t: f64,
    pub q_v: f64,
    pub q_n: Vec<f64>,
    pub q_x: Vec<f64>,
    pub q_xx: Vec<f64>,
    pub q_xt: Vec<f64>,
    pub q_xv: Vec<f64>,
    pub q_xn: Vec<f64>,
    pub q_tt: f64,
    pub q_tv: f64,
    pub q_tn: Vec<f64>,
    pub q_vv: f64,
    pub q_vn: Vec<f64>,
    pub q_nn: Vec<f64>,
}

/// First derivatives of equilibrium site fractions `X_k(T, V, n)`.
#[derive(Debug, Clone)]
pub struct AssociationXDerivatives {
    pub x_t: Vec<f64>,
    pub x_v: Vec<f64>,
    /// Row-major `num_sites * nc`: `x_n[k * nc + i] = dX_k / dn_i`.
    pub x_n: Vec<f64>,
}

struct DeltaTables {
    delta: Vec<f64>,
    delta_t: Vec<f64>,
    delta_v: Vec<f64>,
    delta_n: Vec<f64>,
    delta_tt: Vec<f64>,
    delta_tv: Vec<f64>,
    delta_tn: Vec<f64>,
    delta_vv: Vec<f64>,
    delta_vn: Vec<f64>,
    delta_nn: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DeltaOrder {
    Value,
    First,
    Second,
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

fn calc_boltzmann_fac(assoc: &AssociationParams, t: f64) -> Result<Vec<f64>> {
    let ns = assoc.num_sites;
    let mut b = vec![1.0; ns * ns];
    let rt = RGAS * t;
    for k in 0..ns {
        for l in k..ns {
            let eps = assoc.eps_kl[k * ns + l];
            let val = if eps > 0.0 { (eps / rt).exp() } else { 1.0 };
            if !val.is_finite() {
                return Err(Error::InvalidInput(
                    "association Boltzmann factor overflowed to non-finite value",
                ));
            }
            b[k * ns + l] = val;
            b[l * ns + k] = val;
        }
    }
    Ok(b)
}

fn assemble_delta(
    assoc: &AssociationParams,
    provider: &dyn AssociationRdfProvider,
    t: f64,
    order: DeltaOrder,
) -> Result<DeltaTables> {
    let ns = assoc.num_sites;
    let nc = assoc.num_components;
    let boltz = calc_boltzmann_fac(assoc, t)?;
    let pair_dep = provider.pair_dependent();

    let mut delta = vec![0.0; ns * ns];
    let mut delta_t = if order != DeltaOrder::Value {
        vec![0.0; ns * ns]
    } else {
        Vec::new()
    };
    let mut delta_v = if order != DeltaOrder::Value {
        vec![0.0; ns * ns]
    } else {
        Vec::new()
    };
    let mut delta_n = if order != DeltaOrder::Value {
        vec![0.0; ns * ns * nc]
    } else {
        Vec::new()
    };
    let mut delta_tt = if order == DeltaOrder::Second {
        vec![0.0; ns * ns]
    } else {
        Vec::new()
    };
    let mut delta_tv = if order == DeltaOrder::Second {
        vec![0.0; ns * ns]
    } else {
        Vec::new()
    };
    let mut delta_vv = if order == DeltaOrder::Second {
        vec![0.0; ns * ns]
    } else {
        Vec::new()
    };
    let mut delta_tn = if order == DeltaOrder::Second {
        vec![0.0; ns * ns * nc]
    } else {
        Vec::new()
    };
    let mut delta_vn = if order == DeltaOrder::Second {
        vec![0.0; ns * ns * nc]
    } else {
        Vec::new()
    };
    let mut delta_nn = if order == DeltaOrder::Second {
        vec![0.0; ns * ns * nc * nc]
    } else {
        Vec::new()
    };

    // Cache single RDF if pair-independent (sPC-SAFT / CPA).
    let shared_val = if !pair_dep && order == DeltaOrder::Value {
        Some(provider.rdf_value(0, 0)?)
    } else {
        None
    };
    let shared_first = if !pair_dep && order == DeltaOrder::First {
        Some(provider.rdf_first(0, 0)?)
    } else {
        None
    };
    let shared_second = if !pair_dep && order == DeltaOrder::Second {
        Some(provider.rdf_second(0, 0)?)
    } else {
        None
    };

    let rt2 = RGAS * t * t;
    let rt3 = rt2 * t;
    let rt = RGAS * t;

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
            let expo = boltz[kl];
            let expo_m1 = if eps > 0.0 { (eps / rt).exp_m1() } else { 0.0 };

            match order {
                DeltaOrder::Value => {
                    let (g, covol) = if pair_dep {
                        let r = provider.rdf_value(ic, jc)?;
                        (r.g, r.covol)
                    } else {
                        let covol = provider.rdf_value(ic, jc)?.covol;
                        (shared_val.as_ref().unwrap().g, covol)
                    };
                    let d = g * beta * covol * expo_m1;
                    if !d.is_finite() {
                        return Err(Error::InvalidInput("non-finite association Delta"));
                    }
                    delta[kl] = d;
                    delta[lk] = d;
                }
                DeltaOrder::First => {
                    let owned;
                    let (r, covol, covol_t) = if pair_dep {
                        owned = provider.rdf_first(ic, jc)?;
                        (&owned, owned.covol, owned.covol_t)
                    } else {
                        let pair_r = provider.rdf_first(ic, jc)?;
                        (shared_first.as_ref().unwrap(), pair_r.covol, pair_r.covol_t)
                    };
                    let h = beta * covol * expo_m1;
                    let h_t = -eps * expo * covol * beta / rt2 + beta * covol_t * expo_m1;
                    let d = r.g * h;
                    let dt = r.g_t * h + r.g * h_t;
                    let dv = r.g_v * h;
                    if !d.is_finite() || !dt.is_finite() || !dv.is_finite() {
                        return Err(Error::InvalidInput("non-finite association Delta"));
                    }
                    delta[kl] = d;
                    delta[lk] = d;
                    delta_t[kl] = dt;
                    delta_t[lk] = dt;
                    delta_v[kl] = dv;
                    delta_v[lk] = dv;
                    for i in 0..nc {
                        let dn = r.g_n[i] * h;
                        if !dn.is_finite() {
                            return Err(Error::InvalidInput("non-finite association Delta_n"));
                        }
                        delta_n[kl * nc + i] = dn;
                        delta_n[lk * nc + i] = dn;
                    }
                }
                DeltaOrder::Second => {
                    let owned;
                    let (r, covol, covol_t, covol_tt) = if pair_dep {
                        owned = provider.rdf_second(ic, jc)?;
                        (&owned.rdf, owned.covol, owned.covol_t, owned.covol_tt)
                    } else {
                        let pair_r = provider.rdf_second(ic, jc)?;
                        (
                            &shared_second.as_ref().unwrap().rdf,
                            pair_r.covol,
                            pair_r.covol_t,
                            pair_r.covol_tt,
                        )
                    };
                    let h = beta * covol * expo_m1;
                    let h_t = -eps * expo * covol * beta / rt2 + beta * covol_t * expo_m1;
                    let h_tt = (2.0 + eps / rt) * eps / rt3 * expo * beta * covol
                        - 2.0 * eps * expo * covol_t * beta / rt2
                        + beta * covol_tt * expo_m1;

                    let d = r.g * h;
                    let dt = r.temperature * h + r.g * h_t;
                    let dv = r.volume * h;
                    let dtt = r.temperature_second * h + 2.0 * r.temperature * h_t + r.g * h_tt;
                    let dtv = r.temperature_volume * h + r.volume * h_t;
                    let dvv = r.volume_second * h;

                    if ![d, dt, dv, dtt, dtv, dvv].iter().all(|v| v.is_finite()) {
                        return Err(Error::InvalidInput("non-finite association Delta"));
                    }
                    delta[kl] = d;
                    delta[lk] = d;
                    delta_t[kl] = dt;
                    delta_t[lk] = dt;
                    delta_v[kl] = dv;
                    delta_v[lk] = dv;
                    delta_tt[kl] = dtt;
                    delta_tt[lk] = dtt;
                    delta_tv[kl] = dtv;
                    delta_tv[lk] = dtv;
                    delta_vv[kl] = dvv;
                    delta_vv[lk] = dvv;

                    for i in 0..nc {
                        let dn = r.moles[i] * h;
                        let dtn = r.temperature_moles[i] * h + r.moles[i] * h_t;
                        let dvn = r.volume_moles[i] * h;
                        delta_n[kl * nc + i] = dn;
                        delta_n[lk * nc + i] = dn;
                        delta_tn[kl * nc + i] = dtn;
                        delta_tn[lk * nc + i] = dtn;
                        delta_vn[kl * nc + i] = dvn;
                        delta_vn[lk * nc + i] = dvn;
                        for j in 0..nc {
                            let dnn = r.moles_second[i * nc + j] * h;
                            delta_nn[(kl * nc + i) * nc + j] = dnn;
                            delta_nn[(lk * nc + i) * nc + j] = dnn;
                        }
                    }
                }
            }
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
                if assoc.beta_kl[k1 * ns + k2] == 0.0 || assoc.beta_kl[l1 * ns + l2] == 0.0 {
                    return Err(Error::Unsupported(
                        "Elliot combining rule requires nonzero self-association beta",
                    ));
                }
                let ka = k1 * ns + k2;
                let lb = l1 * ns + l2;
                let da = delta[ka];
                let db = delta[lb];
                if da < 0.0 || db < 0.0 {
                    return Err(Error::InvalidInput("negative Delta in Elliot rule"));
                }
                for k in k1..=k2 {
                    for l in l1..=l2 {
                        let kl = k * ns + l;
                        let lk = l * ns + k;
                        if assoc.beta_kl[kl] <= 0.0 {
                            continue;
                        }
                        let d = (da * db).sqrt();
                        delta[kl] = d;
                        delta[lk] = d;
                        if d == 0.0 {
                            continue;
                        }
                        let inv_2d = 1.0 / (2.0 * d);
                        let inv_4d3 = 1.0 / (4.0 * d * d * d);

                        if order != DeltaOrder::Value {
                            let dt = (delta_t[ka] * db + da * delta_t[lb]) * inv_2d;
                            let dv = (delta_v[ka] * db + da * delta_v[lb]) * inv_2d;
                            delta_t[kl] = dt;
                            delta_t[lk] = dt;
                            delta_v[kl] = dv;
                            delta_v[lk] = dv;
                            for i in 0..nc {
                                let dn = (delta_n[ka * nc + i] * db + da * delta_n[lb * nc + i])
                                    * inv_2d;
                                delta_n[kl * nc + i] = dn;
                                delta_n[lk * nc + i] = dn;
                            }
                        }
                        if order == DeltaOrder::Second {
                            let st = delta_t[ka] * db + da * delta_t[lb];
                            let sv = delta_v[ka] * db + da * delta_v[lb];
                            let dtt = (delta_tt[ka] * db
                                + 2.0 * delta_t[ka] * delta_t[lb]
                                + da * delta_tt[lb])
                                * inv_2d
                                - st * st * inv_4d3;
                            let dtv = (delta_tv[ka] * db
                                + delta_t[ka] * delta_v[lb]
                                + delta_v[ka] * delta_t[lb]
                                + da * delta_tv[lb])
                                * inv_2d
                                - st * sv * inv_4d3;
                            let dvv = (delta_vv[ka] * db
                                + 2.0 * delta_v[ka] * delta_v[lb]
                                + da * delta_vv[lb])
                                * inv_2d
                                - sv * sv * inv_4d3;
                            delta_tt[kl] = dtt;
                            delta_tt[lk] = dtt;
                            delta_tv[kl] = dtv;
                            delta_tv[lk] = dtv;
                            delta_vv[kl] = dvv;
                            delta_vv[lk] = dvv;

                            for i in 0..nc {
                                let sn_i = delta_n[ka * nc + i] * db + da * delta_n[lb * nc + i];
                                let dtn = (delta_tn[ka * nc + i] * db
                                    + delta_t[ka] * delta_n[lb * nc + i]
                                    + delta_n[ka * nc + i] * delta_t[lb]
                                    + da * delta_tn[lb * nc + i])
                                    * inv_2d
                                    - st * sn_i * inv_4d3;
                                let dvn = (delta_vn[ka * nc + i] * db
                                    + delta_v[ka] * delta_n[lb * nc + i]
                                    + delta_n[ka * nc + i] * delta_v[lb]
                                    + da * delta_vn[lb * nc + i])
                                    * inv_2d
                                    - sv * sn_i * inv_4d3;
                                delta_tn[kl * nc + i] = dtn;
                                delta_tn[lk * nc + i] = dtn;
                                delta_vn[kl * nc + i] = dvn;
                                delta_vn[lk * nc + i] = dvn;
                                for j in 0..nc {
                                    let sn_j =
                                        delta_n[ka * nc + j] * db + da * delta_n[lb * nc + j];
                                    let dnn = (delta_nn[(ka * nc + i) * nc + j] * db
                                        + delta_n[ka * nc + i] * delta_n[lb * nc + j]
                                        + delta_n[ka * nc + j] * delta_n[lb * nc + i]
                                        + da * delta_nn[(lb * nc + i) * nc + j])
                                        * inv_2d
                                        - sn_i * sn_j * inv_4d3;
                                    delta_nn[(kl * nc + i) * nc + j] = dnn;
                                    delta_nn[(lk * nc + i) * nc + j] = dnn;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    Ok(DeltaTables {
        delta,
        delta_t,
        delta_v,
        delta_n,
        delta_tt,
        delta_tv,
        delta_tn,
        delta_vv,
        delta_vn,
        delta_nn,
    })
}

/// Solve a small dense linear system `A x = b` via Gaussian elimination with
/// partial pivoting (in-place on copies, zero `unsafe`).
fn solve_linear_system(n: usize, a_in: &[f64], b_in: &[f64]) -> Result<Vec<f64>> {
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut a = a_in.to_vec();
    let mut b = b_in.to_vec();

    for col in 0..n {
        let mut pivot_row = col;
        let mut pivot_abs = a[col * n + col].abs();
        for row in (col + 1)..n {
            let v = a[row * n + col].abs();
            if v > pivot_abs {
                pivot_abs = v;
                pivot_row = row;
            }
        }
        if !pivot_abs.is_finite() || pivot_abs <= 1e-300 {
            return Err(Error::InvalidInput(
                "singular linear system in association solver",
            ));
        }
        if pivot_row != col {
            for j in col..n {
                a.swap(col * n + j, pivot_row * n + j);
            }
            b.swap(col, pivot_row);
        }
        let diag = a[col * n + col];
        for row in (col + 1)..n {
            let factor = a[row * n + col] / diag;
            a[row * n + col] = 0.0;
            for j in (col + 1)..n {
                a[row * n + j] -= factor * a[col * n + j];
            }
            b[row] -= factor * b[col];
        }
    }

    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let mut sum = b[row];
        for j in (row + 1)..n {
            sum -= a[row * n + j] * x[j];
        }
        let val = sum / a[row * n + row];
        if !val.is_finite() {
            return Err(Error::InvalidInput(
                "non-finite solution in association linear solve",
            ));
        }
        x[row] = val;
    }
    Ok(x)
}

/// Solve `A X = B` for `nrhs` right-hand side columns stored as `rhs[col][row]`.
fn solve_linear_multi(n: usize, a_in: &[f64], rhs: &mut [Vec<f64>]) -> Result<()> {
    if n == 0 || rhs.is_empty() {
        return Ok(());
    }
    let mut a = a_in.to_vec();
    let nrhs = rhs.len();

    for col in 0..n {
        let mut pivot_row = col;
        let mut pivot_abs = a[col * n + col].abs();
        for row in (col + 1)..n {
            let v = a[row * n + col].abs();
            if v > pivot_abs {
                pivot_abs = v;
                pivot_row = row;
            }
        }
        if !pivot_abs.is_finite() || pivot_abs <= 1e-300 {
            return Err(Error::InvalidInput(
                "singular linear system in association derivatives",
            ));
        }
        if pivot_row != col {
            for j in col..n {
                a.swap(col * n + j, pivot_row * n + j);
            }
            for r in rhs.iter_mut() {
                r.swap(col, pivot_row);
            }
        }
        let diag = a[col * n + col];
        for row in (col + 1)..n {
            let factor = a[row * n + col] / diag;
            a[row * n + col] = 0.0;
            for j in (col + 1)..n {
                a[row * n + j] -= factor * a[col * n + j];
            }
            for r in rhs.iter_mut() {
                r[row] -= factor * r[col];
            }
        }
    }

    for row in (0..n).rev() {
        let diag = a[row * n + row];
        for r in rhs.iter_mut().take(nrhs) {
            let mut sum = r[row];
            for j in (row + 1)..n {
                sum -= a[row * n + j] * r[j];
            }
            let val = sum / diag;
            if !val.is_finite() {
                return Err(Error::InvalidInput(
                    "non-finite solution in association derivatives",
                ));
            }
            r[row] = val;
        }
    }
    Ok(())
}

/// Solve for equilibrium non-bonded site fractions `X_k` in `(0, 1]` at unit
/// total moles `(v_u = V / N, x = n / N)` given precomputed `Delta` matrix.
///
/// Mirrors `saft_association::solve_for_X_k` using Michelsen's modified Hessian
/// damped Newton method (80% reduction limit) with successive-substitution
/// fallback, plus exact infinite-dilution evaluation for zero-mole sites.
fn solve_x_k_from_delta(
    assoc: &AssociationParams,
    v_u: f64,
    x_mol: &[f64],
    delta: &[f64],
    tol: f64,
) -> Result<Vec<f64>> {
    let ns = assoc.num_sites;
    if ns == 0 {
        return Ok(Vec::new());
    }
    let mut m_u = vec![0.0; ns];
    let mut active_sites = Vec::with_capacity(ns);
    for (k, &ic) in assoc.site_to_comp.iter().enumerate() {
        let mk = x_mol[ic];
        m_u[k] = mk;
        if mk > 0.0 {
            active_sites.push(k);
        }
    }

    if active_sites.is_empty() {
        return Ok(vec![1.0; ns]);
    }

    // Check whether any active-active Delta entry is positive.
    // Even if `!any_delta` (e.g. a present non-self-associating `Scheme1ea`
    // component with an absent cross-associating species), active sites have
    // `X_k = 1.0`, while absent sites still need the infinite-dilution pass below.
    let mut any_delta = false;
    for &k in &active_sites {
        for &l in &active_sites {
            if delta[k * ns + l] > 0.0 {
                any_delta = true;
                break;
            }
        }
        if any_delta {
            break;
        }
    }

    let mut x_k = vec![1.0; ns];
    if any_delta {
        let na = active_sites.len();
        // Active scaled matrix: C_{ab} = (m_b / v_u) * Delta_{ab}
        // Site equation: 1 / X_a - 1 - sum_b C_{ab} X_b = 0.
        // Multiply row a by m_a gives symmetric Michelsen K_{ab} = m_a * C_{ab} = m_a m_b Delta_{ab} / v_u.
        let mut c_mat = vec![0.0; na * na];
        for (a, &ka) in active_sites.iter().enumerate() {
            for (b, &kb) in active_sites.iter().enumerate() {
                let val = m_u[kb] * delta[ka * ns + kb] / v_u;
                if !val.is_finite() || val < 0.0 {
                    return Err(Error::InvalidInput(
                        "invalid association strength matrix entry",
                    ));
                }
                c_mat[a * na + b] = val;
            }
        }

        let mut xa = vec![0.2; na];
        let mut resid = vec![0.0; na];
        let mut jac = vec![0.0; na * na];

        // Michelsen modified-Newton iterations on unweighted site residual
        // g_a = (1 / X_a - 1) - sum_b C_{ab} X_b = 0 so trace active components
        // (e.g. x_i = 1e-12) converge to full relative precision too!
        let mut converged = false;
        for _ in 0..40 {
            let mut max_rel_err = 0.0_f64;
            for a in 0..na {
                let mut sum_cx = 0.0;
                for b in 0..na {
                    sum_cx += c_mat[a * na + b] * xa[b];
                }
                let inv_x = 1.0 / xa[a];
                let g_a = (inv_x - 1.0) - sum_cx;
                resid[a] = g_a;
                let scale = inv_x.max(1.0);
                max_rel_err = max_rel_err.max(g_a.abs() / scale);

                for b in 0..na {
                    jac[a * na + b] = -c_mat[a * na + b];
                }
                // Michelsen modified diagonal: -(1 + sum_b C_{ab} X_b) / X_a
                jac[a * na + a] -= (1.0 + sum_cx) / xa[a];
            }

            if !max_rel_err.is_finite() {
                break;
            }
            if max_rel_err <= tol {
                converged = true;
                break;
            }

            // Solve jac * dx = -resid
            let rhs: Vec<f64> = resid.iter().map(|&v| -v).collect();
            let Ok(mut dx) = solve_linear_system(na, &jac, &rhs) else {
                break;
            };

            // Step limiting matching `saft_association::limit` (max 80% reduction, X <= 1).
            let mut scaling = 1.0_f64;
            for a in 0..na {
                if xa[a] + dx[a] < 0.2 * xa[a] {
                    scaling = scaling.min(-0.8 * xa[a] / dx[a]);
                }
                if xa[a] + dx[a] > 1.0 && dx[a].abs() > 1e-9 {
                    scaling = scaling.min((1.0 - xa[a]) / dx[a]);
                }
            }
            if scaling < 1.0 {
                for d in &mut dx {
                    *d *= scaling;
                }
            }

            let mut max_step = 0.0_f64;
            for a in 0..na {
                xa[a] = (xa[a] + dx[a]).clamp(1e-300, 1.0);
                max_step = max_step.max(dx[a].abs());
            }
            if max_step <= tol * 0.1 && max_rel_err <= tol * 10.0 {
                // Re-check residual once more on next iteration before accepting.
                continue;
            }
        }

        // Damped successive substitution fallback (`saft_association::succ_subs`) + Newton polish.
        if !converged {
            xa.fill(0.2);
            let omega = 0.2;
            for _ in 0..400 {
                let mut max_diff = 0.0_f64;
                for a in 0..na {
                    let mut sum_cx = 0.0;
                    for b in 0..na {
                        sum_cx += c_mat[a * na + b] * xa[b];
                    }
                    let x_new = (1.0 - omega) / (1.0 + sum_cx) + omega * xa[a];
                    max_diff = max_diff.max((x_new - xa[a]).abs());
                    xa[a] = x_new;
                }
                if max_diff < 1e6 * f64::EPSILON {
                    break;
                }
            }
            // Final Newton polish to tighten residual to `tol`.
            for _ in 0..20 {
                let mut max_rel_err = 0.0_f64;
                for a in 0..na {
                    let mut sum_cx = 0.0;
                    for b in 0..na {
                        sum_cx += c_mat[a * na + b] * xa[b];
                    }
                    let inv_x = 1.0 / xa[a];
                    let g_a = (inv_x - 1.0) - sum_cx;
                    resid[a] = g_a;
                    max_rel_err = max_rel_err.max(g_a.abs() / inv_x.max(1.0));
                    for b in 0..na {
                        jac[a * na + b] = -c_mat[a * na + b];
                    }
                    jac[a * na + a] -= (1.0 + sum_cx) / xa[a];
                }
                if max_rel_err.is_finite() && max_rel_err <= tol.max(1e-10) {
                    converged = true;
                    break;
                }
                let rhs: Vec<f64> = resid.iter().map(|&v| -v).collect();
                let Ok(dx) = solve_linear_system(na, &jac, &rhs) else {
                    break;
                };
                let mut scaling = 1.0_f64;
                for a in 0..na {
                    if xa[a] + dx[a] < 0.2 * xa[a] {
                        scaling = scaling.min(-0.8 * xa[a] / dx[a]);
                    }
                    if xa[a] + dx[a] > 1.0 && dx[a].abs() > 1e-9 {
                        scaling = scaling.min((1.0 - xa[a]) / dx[a]);
                    }
                }
                for a in 0..na {
                    xa[a] = (xa[a] + scaling * dx[a]).clamp(1e-300, 1.0);
                }
            }
        }

        if !converged {
            return Err(Error::NotConverged(
                "association site fraction solver did not converge",
            ));
        }

        for (a, &ka) in active_sites.iter().enumerate() {
            if !xa[a].is_finite() || xa[a] <= 0.0 || xa[a] > 1.0 + 1e-12 {
                return Err(Error::NotConverged(
                    "non-physical association site fraction",
                ));
            }
            x_k[ka] = xa[a].min(1.0);
        }
    }

    // Exact infinite-dilution limit for zero-mole associating sites:
    // X_k = 1 / (1 + sum_{l active} (m_u[l] / v_u) * Delta_{kl} * X_l).
    for k in 0..ns {
        if m_u[k] == 0.0 {
            let mut sum_cx = 0.0;
            for &l in &active_sites {
                sum_cx += (m_u[l] / v_u) * delta[k * ns + l] * x_k[l];
            }
            let val = 1.0 / (1.0 + sum_cx);
            if !val.is_finite() || val <= 0.0 {
                return Err(Error::InvalidInput(
                    "non-finite infinite-dilution site fraction",
                ));
            }
            x_k[k] = val;
        }
    }

    Ok(x_k)
}

/// Public entry point to solve for non-bonded site fractions `X_k` at `(T, V, n)`.
///
/// The `provider` must be constructed at unit total moles `(v_u = V / N, x = n / N)`
/// (since `g_{ij}` is 0-homogeneous in `(V, n)`).
pub fn solve_for_x_k(
    assoc: &AssociationParams,
    provider: &dyn AssociationRdfProvider,
    t: f64,
    v: f64,
    n: &[f64],
    tol: Option<f64>,
) -> Result<Vec<f64>> {
    let sumn = validate_state(t, v, n, assoc.num_components)?;
    if !assoc.has_sites() {
        return Ok(Vec::new());
    }
    let v_u = v / sumn;
    let x_mol: Vec<f64> = n.iter().map(|&ni| ni / sumn).collect();
    let dtab = assemble_delta(assoc, provider, t, DeltaOrder::Value)?;
    let tol = tol.unwrap_or(1e5 * f64::EPSILON);
    solve_x_k_from_delta(assoc, v_u, &x_mol, &dtab.delta, tol)
}

/// Compute value-only association reduced residual Helmholtz energy `F = A_assoc / (R T)`
/// and site fractions `X_k` without allocating derivative tensors.
///
/// `provider` must correspond to unit total moles `(v_u = V / N, x = n / N)`.
pub fn assoc_helmholtz_value(
    assoc: &AssociationParams,
    provider: &dyn AssociationRdfProvider,
    t: f64,
    v: f64,
    n: &[f64],
) -> Result<(f64, Vec<f64>)> {
    let sumn = validate_state(t, v, n, assoc.num_components)?;
    if !assoc.has_sites() {
        return Ok((0.0, Vec::new()));
    }
    let v_u = v / sumn;
    let x_mol: Vec<f64> = n.iter().map(|&ni| ni / sumn).collect();
    let dtab = assemble_delta(assoc, provider, t, DeltaOrder::Value)?;
    let x_k = solve_x_k_from_delta(assoc, v_u, &x_mol, &dtab.delta, 1e5 * f64::EPSILON)?;

    let mut q_u = 0.0;
    for (k, &ic) in assoc.site_to_comp.iter().enumerate() {
        let mk = x_mol[ic];
        if mk > 0.0 {
            q_u += mk * (x_k[k].ln() + 0.5 * (1.0 - x_k[k]));
        }
    }
    let f = sumn * q_u;
    if !f.is_finite() {
        return Err(Error::InvalidInput(
            "non-finite association Helmholtz energy",
        ));
    }
    Ok((f, x_k))
}

/// Compute association reduced residual Helmholtz energy and first `(T, V, n)` derivatives.
///
/// At equilibrium `Q_X = 0`, first derivatives require NO linear solve for `X_T, X_V, X_n`
/// and NO second-derivative `Delta` tensors.
/// `provider` must correspond to unit total moles `(v_u = V / N, x = n / N)`.
pub fn assoc_helmholtz_first(
    assoc: &AssociationParams,
    provider: &dyn AssociationRdfProvider,
    t: f64,
    v: f64,
    n: &[f64],
) -> Result<AssociationFirst> {
    let sumn = validate_state(t, v, n, assoc.num_components)?;
    let nc = assoc.num_components;
    if !assoc.has_sites() {
        return Ok(AssociationFirst {
            f: 0.0,
            f_t: 0.0,
            f_v: 0.0,
            f_n: vec![0.0; nc],
            x_k: Vec::new(),
        });
    }
    let v_u = v / sumn;
    let x_mol: Vec<f64> = n.iter().map(|&ni| ni / sumn).collect();
    let dtab = assemble_delta(assoc, provider, t, DeltaOrder::First)?;
    let x_k = solve_x_k_from_delta(assoc, v_u, &x_mol, &dtab.delta, 1e5 * f64::EPSILON)?;

    let ns = assoc.num_sites;
    let m_u: Vec<f64> = assoc.site_to_comp.iter().map(|&ic| x_mol[ic]).collect();
    let inv_v = 1.0 / v_u;

    let mut q_u = 0.0;
    for k in 0..ns {
        if m_u[k] > 0.0 {
            q_u += m_u[k] * (x_k[k].ln() + 0.5 * (1.0 - x_k[k]));
        }
    }

    let mut q_t_u = 0.0;
    let mut q_v_u = 0.0;
    let mut q_n_u = vec![0.0; nc];

    for k in 0..ns {
        for l in 0..ns {
            let kl = k * ns + l;
            let qty = 0.5 * m_u[k] * m_u[l] * x_k[k] * x_k[l] * inv_v;
            if qty == 0.0 {
                continue;
            }
            q_t_u -= qty * dtab.delta_t[kl];
            q_v_u += qty * (dtab.delta[kl] * inv_v - dtab.delta_v[kl]);
            for i in 0..nc {
                q_n_u[i] -= qty * dtab.delta_n[kl * nc + i];
            }
        }
    }
    for i in 0..nc {
        if let Some((k1, k2)) = assoc.comp_vs_sites[i] {
            for &xk in &x_k[k1..=k2] {
                q_n_u[i] += xk.ln();
            }
        }
    }

    let out = AssociationFirst {
        f: sumn * q_u,
        f_t: sumn * q_t_u,
        f_v: q_v_u,
        f_n: q_n_u,
        x_k,
    };
    if ![out.f, out.f_t, out.f_v].iter().all(|v| v.is_finite())
        || !out.f_n.iter().all(|v| v.is_finite())
    {
        return Err(Error::InvalidInput(
            "non-finite association first derivatives",
        ));
    }
    Ok(out)
}

/// Compute association `F, F_T, F_V, F_TT, F_TV, F_VV` (solving only for `X_T` and `X_V`,
/// avoiding `n^2` mole Hessian assembly).
///
/// `provider` must correspond to unit total moles `(v_u = V / N, x = n / N)`.
pub fn assoc_helmholtz_pressure_jet(
    assoc: &AssociationParams,
    provider: &dyn AssociationRdfProvider,
    t: f64,
    v: f64,
    n: &[f64],
) -> Result<AssociationPressureJet> {
    let sumn = validate_state(t, v, n, assoc.num_components)?;
    if !assoc.has_sites() {
        return Ok(AssociationPressureJet {
            f: 0.0,
            f_t: 0.0,
            f_v: 0.0,
            f_tt: 0.0,
            f_tv: 0.0,
            f_vv: 0.0,
            x_k: Vec::new(),
        });
    }
    let v_u = v / sumn;
    let x_mol: Vec<f64> = n.iter().map(|&ni| ni / sumn).collect();
    let dtab = assemble_delta(assoc, provider, t, DeltaOrder::Second)?;
    let x_k = solve_x_k_from_delta(assoc, v_u, &x_mol, &dtab.delta, 1e5 * f64::EPSILON)?;

    let ns = assoc.num_sites;
    let m_u: Vec<f64> = assoc.site_to_comp.iter().map(|&ic| x_mol[ic]).collect();
    let inv_v = 1.0 / v_u;
    let inv_v2 = inv_v * inv_v;

    let mut q_u = 0.0;
    for k in 0..ns {
        if m_u[k] > 0.0 {
            q_u += m_u[k] * (x_k[k].ln() + 0.5 * (1.0 - x_k[k]));
        }
    }

    let mut q_t_u = 0.0;
    let mut q_v_u = 0.0;
    let mut q_tt_u = 0.0;
    let mut q_tv_u = 0.0;
    let mut q_vv_u = 0.0;
    let mut q_xt = vec![0.0; ns];
    let mut q_xv = vec![0.0; ns];
    let mut q_xx = vec![0.0; ns * ns];

    for k in 0..ns {
        if m_u[k] == 0.0 {
            q_xx[k * ns + k] = 1.0;
            continue;
        }
        q_xx[k * ns + k] -= m_u[k] / (x_k[k] * x_k[k]);
        for l in 0..ns {
            let kl = k * ns + l;
            let k_mich = m_u[k] * m_u[l] * dtab.delta[kl] * inv_v;
            q_xx[kl] -= k_mich;

            let qty = 0.5 * m_u[k] * m_u[l] * x_k[k] * x_k[l] * inv_v;
            q_t_u -= qty * dtab.delta_t[kl];
            q_v_u += qty * (dtab.delta[kl] * inv_v - dtab.delta_v[kl]);
            q_tt_u -= qty * dtab.delta_tt[kl];
            q_tv_u += qty * (dtab.delta_t[kl] * inv_v - dtab.delta_tv[kl]);
            q_vv_u += qty
                * (2.0 * (-dtab.delta[kl] * inv_v + dtab.delta_v[kl]) * inv_v - dtab.delta_vv[kl]);

            q_xt[k] -= m_u[k] * m_u[l] * x_k[l] * dtab.delta_t[kl] * inv_v;
            q_xv[k] +=
                m_u[k] * m_u[l] * x_k[l] * (dtab.delta[kl] * inv_v2 - dtab.delta_v[kl] * inv_v);
        }
    }

    let mut rhs = vec![
        q_xt.iter().map(|&v| -v).collect::<Vec<_>>(),
        q_xv.iter().map(|&v| -v).collect::<Vec<_>>(),
    ];
    solve_linear_multi(ns, &q_xx, &mut rhs)?;
    let x_t = &rhs[0];
    let x_v = &rhs[1];

    let mut f_tt_u = q_tt_u;
    let mut f_tv_u = q_tv_u;
    let mut f_vv_u = q_vv_u;
    for k in 0..ns {
        f_tt_u += q_xt[k] * x_t[k];
        f_tv_u += q_xv[k] * x_t[k];
        f_vv_u += q_xv[k] * x_v[k];
    }

    let out = AssociationPressureJet {
        f: sumn * q_u,
        f_t: sumn * q_t_u,
        f_v: q_v_u,
        f_tt: sumn * f_tt_u,
        f_tv: f_tv_u,
        f_vv: f_vv_u / sumn,
        x_k,
    };
    if ![out.f, out.f_t, out.f_v, out.f_tt, out.f_tv, out.f_vv]
        .iter()
        .all(|v| v.is_finite())
    {
        return Err(Error::InvalidInput(
            "non-finite association pressure-jet derivatives",
        ));
    }
    Ok(out)
}

/// Evaluate `Q_derivatives_knowing_X` at `(T, V, n, X_k)` holding `X_k` fixed.
///
/// Mirrors `saft_association::Q_derivatives_knowing_X`.
/// `provider` must correspond to `(T, V, n)` (or unit-mole `(T, v_u, x)` when called internally).
pub fn q_derivatives_knowing_x(
    assoc: &AssociationParams,
    provider: &dyn AssociationRdfProvider,
    t: f64,
    v: f64,
    n: &[f64],
    x_k: &[f64],
    x_calculated: bool,
) -> Result<AssociationQDerivatives> {
    let _sumn = validate_state(t, v, n, assoc.num_components)?;
    let ns = assoc.num_sites;
    if x_k.len() != ns {
        return Err(Error::InvalidInput("X_k length mismatch"));
    }
    for &xk in x_k {
        if !xk.is_finite() || xk <= 0.0 {
            return Err(Error::InvalidInput("positive finite X_k required"));
        }
    }
    let dtab = assemble_delta(assoc, provider, t, DeltaOrder::Second)?;
    Ok(q_derivatives_from_delta(
        assoc,
        v,
        n,
        x_k,
        &dtab,
        x_calculated,
    ))
}

fn q_derivatives_from_delta(
    assoc: &AssociationParams,
    v: f64,
    n: &[f64],
    x_k: &[f64],
    dtab: &DeltaTables,
    x_calculated: bool,
) -> AssociationQDerivatives {
    let ns = assoc.num_sites;
    let nc = assoc.num_components;
    let inv_v = 1.0 / v;
    let inv_v2 = inv_v * inv_v;
    let m: Vec<f64> = assoc.site_to_comp.iter().map(|&ic| n[ic]).collect();
    let log_x: Vec<f64> = x_k.iter().map(|&x| x.ln()).collect();

    let mut qty = vec![0.0; ns * ns];
    let mut kmich = vec![0.0; ns * ns];
    for k in 0..ns {
        for l in k..ns {
            let kl = k * ns + l;
            let lk = l * ns + k;
            let qv = 0.5 * m[k] * m[l] * x_k[k] * x_k[l] * inv_v;
            qty[kl] = qv;
            qty[lk] = qv;
            let km = m[k] * m[l] * dtab.delta[kl] * inv_v;
            kmich[kl] = km;
            kmich[lk] = km;
        }
    }

    let mut q = 0.0;
    if x_calculated {
        for k in 0..ns {
            if m[k] > 0.0 {
                q += m[k] * (log_x[k] + 0.5 * (1.0 - x_k[k]));
            }
        }
    } else {
        for k in 0..ns {
            if m[k] > 0.0 {
                q += m[k] * (log_x[k] - x_k[k] + 1.0);
            }
            for l in 0..ns {
                q -= qty[k * ns + l] * dtab.delta[k * ns + l];
            }
        }
    }

    let mut q_t = 0.0;
    let mut q_v = 0.0;
    let mut q_tt = 0.0;
    let mut q_tv = 0.0;
    let mut q_vv = 0.0;
    let mut q_n = vec![0.0; nc];

    for k in 0..ns {
        for l in 0..ns {
            let kl = k * ns + l;
            let qv = qty[kl];
            q_t -= qv * dtab.delta_t[kl];
            q_v += qv * (dtab.delta[kl] * inv_v - dtab.delta_v[kl]);
            q_tt -= qv * dtab.delta_tt[kl];
            q_tv += qv * (dtab.delta_t[kl] * inv_v - dtab.delta_tv[kl]);
            q_vv += qv
                * (2.0 * (-dtab.delta[kl] * inv_v + dtab.delta_v[kl]) * inv_v - dtab.delta_vv[kl]);
            for i in 0..nc {
                q_n[i] -= qv * dtab.delta_n[kl * nc + i];
            }
        }
    }
    if x_calculated {
        for i in 0..nc {
            if let Some((k1, k2)) = assoc.comp_vs_sites[i] {
                for &lx in &log_x[k1..=k2] {
                    q_n[i] += lx;
                }
            }
        }
    } else {
        for i in 0..nc {
            if let Some((k1, k2)) = assoc.comp_vs_sites[i] {
                for k in k1..=k2 {
                    q_n[i] += log_x[k] - x_k[k] + 1.0;
                    for l in 0..ns {
                        q_n[i] -= m[l] * x_k[k] * x_k[l] * dtab.delta[k * ns + l] * inv_v;
                    }
                }
            }
        }
    }

    let mut q_x = vec![0.0; ns];
    if !x_calculated {
        for k in 0..ns {
            if m[k] == 0.0 {
                continue;
            }
            let mut val = m[k] * (1.0 / x_k[k] - 1.0);
            for l in 0..ns {
                val -= m[k] * m[l] * x_k[l] * dtab.delta[k * ns + l] * inv_v;
            }
            q_x[k] = val;
        }
    }

    let mut q_xx = vec![0.0; ns * ns];
    for k in 0..ns {
        for l in 0..ns {
            q_xx[k * ns + l] = -kmich[k * ns + l];
        }
        q_xx[k * ns + k] -= m[k] / (x_k[k] * x_k[k]);
    }

    let mut q_xt = vec![0.0; ns];
    let mut q_xv = vec![0.0; ns];
    for k in 0..ns {
        for l in 0..ns {
            let kl = k * ns + l;
            q_xt[k] -= m[k] * m[l] * x_k[l] * dtab.delta_t[kl] * inv_v;
            q_xv[k] += m[k] * m[l] * x_k[l] * (dtab.delta[kl] * inv_v2 - dtab.delta_v[kl] * inv_v);
        }
    }

    let mut q_xn = vec![0.0; ns * nc];
    for i in 0..nc {
        let sites_i = assoc.comp_vs_sites[i];
        for c_nu in 0..ns {
            let mut acc = 0.0;
            if !x_calculated
                && let Some((k1, k2)) = sites_i
                && c_nu >= k1
                && c_nu <= k2
            {
                acc += 1.0 / x_k[c_nu] - 1.0;
                for b_j in 0..ns {
                    acc -= m[b_j] * x_k[b_j] * dtab.delta[c_nu * ns + b_j] * inv_v;
                }
            }
            if let Some((k1, k2)) = sites_i {
                for a_k in k1..=k2 {
                    acc -= m[c_nu] * x_k[a_k] * dtab.delta[a_k * ns + c_nu] * inv_v;
                }
            }
            for b_j in 0..ns {
                acc -=
                    m[c_nu] * m[b_j] * x_k[b_j] * dtab.delta_n[(c_nu * ns + b_j) * nc + i] * inv_v;
            }
            q_xn[c_nu * nc + i] = acc;
        }
    }

    let mut q_tn = vec![0.0; nc];
    let mut q_vn = vec![0.0; nc];
    for i in 0..nc {
        let sites_i = assoc.comp_vs_sites[i];
        let mut acc_tn = 0.0;
        let mut acc_vn = 0.0;
        for k in 0..ns {
            let in_i = sites_i.is_some_and(|(k1, k2)| k >= k1 && k <= k2);
            if in_i {
                for l in 0..ns {
                    let kl = k * ns + l;
                    let fac = m[l] * x_k[k] * x_k[l] * inv_v;
                    acc_tn -= fac * dtab.delta_t[kl];
                    acc_vn += fac * (dtab.delta[kl] * inv_v - dtab.delta_v[kl]);
                }
            }
            for l in 0..ns {
                let kl = k * ns + l;
                let qv = qty[kl];
                acc_tn -= qv * dtab.delta_tn[kl * nc + i];
                acc_vn += qv * (dtab.delta_n[kl * nc + i] * inv_v - dtab.delta_vn[kl * nc + i]);
            }
        }
        q_tn[i] = acc_tn;
        q_vn[i] = acc_vn;
    }

    let mut q_nn = vec![0.0; nc * nc];
    for i in 0..nc {
        let sites_i = assoc.comp_vs_sites[i];
        for j in 0..nc {
            let sites_j = assoc.comp_vs_sites[j];
            let mut acc = 0.0;
            for k in 0..ns {
                let k_in_i = sites_i.is_some_and(|(k1, k2)| k >= k1 && k <= k2);
                let k_in_j = sites_j.is_some_and(|(l1, l2)| k >= l1 && k <= l2);
                if k_in_i && let Some((l1, l2)) = sites_j {
                    for l in l1..=l2 {
                        acc -= x_k[k] * x_k[l] * dtab.delta[k * ns + l] * inv_v;
                    }
                }
                if k_in_i {
                    for l in 0..ns {
                        acc -= m[l] * x_k[k] * x_k[l] * dtab.delta_n[(k * ns + l) * nc + j] * inv_v;
                    }
                }
                if k_in_j {
                    for l in 0..ns {
                        acc -= m[l] * x_k[k] * x_k[l] * dtab.delta_n[(k * ns + l) * nc + i] * inv_v;
                    }
                }
                for l in 0..ns {
                    let kl = k * ns + l;
                    acc -= qty[kl] * dtab.delta_nn[(kl * nc + i) * nc + j];
                }
            }
            q_nn[i * nc + j] = acc;
        }
    }

    AssociationQDerivatives {
        q,
        q_t,
        q_v,
        q_n,
        q_x,
        q_xx,
        q_xt,
        q_xv,
        q_xn,
        q_tt,
        q_tv,
        q_tn,
        q_vv,
        q_vn,
        q_nn,
    }
}

/// Compute first derivatives `(X_T, X_V, X_n)` of equilibrium site fractions `X_k`
/// at `(T, V, n)` using unit-mole scaling.
///
/// `provider` must correspond to unit total moles `(v_u = V / N, x = n / N)`.
pub fn x_derivatives_knowing_x(
    assoc: &AssociationParams,
    provider: &dyn AssociationRdfProvider,
    t: f64,
    v: f64,
    n: &[f64],
    x_k: &[f64],
) -> Result<AssociationXDerivatives> {
    let sumn = validate_state(t, v, n, assoc.num_components)?;
    let ns = assoc.num_sites;
    if ns == 0 {
        return Ok(AssociationXDerivatives {
            x_t: Vec::new(),
            x_v: Vec::new(),
            x_n: Vec::new(),
        });
    }
    if x_k.len() != ns {
        return Err(Error::InvalidInput("X_k length mismatch"));
    }
    let v_u = v / sumn;
    let x_mol: Vec<f64> = n.iter().map(|&ni| ni / sumn).collect();
    let dtab = assemble_delta(assoc, provider, t, DeltaOrder::First)?;
    let (x_t_u, x_v_u, x_n_u) = solve_x_first_derivatives_unit(assoc, v_u, &x_mol, x_k, &dtab)?;

    // X_k is 0-homogeneous in (V, n):
    // X_T(T, V, n) = X_T(T, v_u, x)
    // X_V(T, V, n) = X_V(T, v_u, x) / N
    // X_n(T, V, n) = X_n(T, v_u, x) / N
    let x_v: Vec<f64> = x_v_u.into_iter().map(|val| val / sumn).collect();
    let x_n: Vec<f64> = x_n_u.into_iter().map(|val| val / sumn).collect();
    if !x_t_u.iter().all(|v| v.is_finite())
        || !x_v.iter().all(|v| v.is_finite())
        || !x_n.iter().all(|v| v.is_finite())
    {
        return Err(Error::InvalidInput(
            "non-finite association site fraction derivatives",
        ));
    }
    Ok(AssociationXDerivatives {
        x_t: x_t_u,
        x_v,
        x_n,
    })
}

/// Solve for `(X_T, X_V, X_n)` in unit-mole coordinates `(v_u, x_mol)` using the
/// unweighted site residual `1/X_k - 1 = sum_l (m_l / v_u) Delta_{kl} X_l`,
/// which remains nonsingular even when some sites have `m_k = 0`.
fn solve_x_first_derivatives_unit(
    assoc: &AssociationParams,
    v_u: f64,
    x_mol: &[f64],
    x_k: &[f64],
    dtab: &DeltaTables,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    let ns = assoc.num_sites;
    let nc = assoc.num_components;
    let inv_v = 1.0 / v_u;
    let inv_v2 = inv_v * inv_v;
    let m: Vec<f64> = assoc.site_to_comp.iter().map(|&ic| x_mol[ic]).collect();

    // Linear operator A_{kl} = delta_{kl} / X_k^2 + (m_l / v_u) * Delta_{kl}.
    // Differentiating 1/X_k - 1 - sum_l (m_l / v_u) Delta_{kl} X_l = 0 gives:
    //   -sum_l A_{kl} dX_l - d(explicit)_k = 0  ==>  A * dX = -d(explicit).
    let mut a_mat = vec![0.0; ns * ns];
    let mut rhs_t = vec![0.0; ns];
    let mut rhs_v = vec![0.0; ns];
    let mut rhs_n: Vec<Vec<f64>> = vec![vec![0.0; ns]; nc];

    for k in 0..ns {
        a_mat[k * ns + k] += 1.0 / (x_k[k] * x_k[k]);
        let mut dt_acc = 0.0;
        let mut dv_acc = 0.0;
        for l in 0..ns {
            let kl = k * ns + l;
            a_mat[kl] += m[l] * dtab.delta[kl] * inv_v;
            dt_acc -= m[l] * x_k[l] * dtab.delta_t[kl] * inv_v;
            dv_acc += m[l] * x_k[l] * (dtab.delta[kl] * inv_v2 - dtab.delta_v[kl] * inv_v);
        }
        rhs_t[k] = dt_acc;
        rhs_v[k] = dv_acc;

        for i in 0..nc {
            let mut dn_acc = 0.0;
            if let Some((l1, l2)) = assoc.comp_vs_sites[i] {
                for l in l1..=l2 {
                    dn_acc -= x_k[l] * dtab.delta[k * ns + l] * inv_v;
                }
            }
            for l in 0..ns {
                dn_acc -= m[l] * x_k[l] * dtab.delta_n[(k * ns + l) * nc + i] * inv_v;
            }
            rhs_n[i][k] = dn_acc;
        }
    }

    let mut all_rhs = Vec::with_capacity(2 + nc);
    all_rhs.push(rhs_t);
    all_rhs.push(rhs_v);
    all_rhs.extend(rhs_n);
    solve_linear_multi(ns, &a_mat, &mut all_rhs)?;

    let x_t = all_rhs[0].clone();
    let x_v = all_rhs[1].clone();
    let mut x_n = vec![0.0; ns * nc];
    for i in 0..nc {
        for k in 0..ns {
            x_n[k * nc + i] = all_rhs[2 + i][k];
        }
    }
    Ok((x_t, x_v, x_n))
}

/// Compute all 10 reduced residual Helmholtz derivative blocks (`F..F_nn`) and
/// solved site fractions `X_k` for association (`saft_association::calcFder_assoc`).
///
/// Evaluated in unit-mole `(v_u = V / N, x = n / N)` coordinates and scaled via
/// exact Euler homogeneity so tiny total mole numbers never overflow or underflow.
/// `provider` must correspond to unit total moles `(v_u, x)`.
pub fn assoc_helmholtz_all(
    assoc: &AssociationParams,
    provider: &dyn AssociationRdfProvider,
    t: f64,
    v: f64,
    n: &[f64],
) -> Result<AssociationHelmholtz> {
    let sumn = validate_state(t, v, n, assoc.num_components)?;
    let nc = assoc.num_components;
    if !assoc.has_sites() {
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
    let v_u = v / sumn;
    let x_mol: Vec<f64> = n.iter().map(|&ni| ni / sumn).collect();
    let dtab = assemble_delta(assoc, provider, t, DeltaOrder::Second)?;
    let x_k = solve_x_k_from_delta(assoc, v_u, &x_mol, &dtab.delta, 1e5 * f64::EPSILON)?;

    let ns = assoc.num_sites;
    let q = q_derivatives_from_delta(assoc, v_u, &x_mol, &x_k, &dtab, true);
    let (x_t_u, x_v_u, x_n_u) = solve_x_first_derivatives_unit(assoc, v_u, &x_mol, &x_k, &dtab)?;

    // Second total derivatives via Michelsen's stationary-point formula:
    // F_{ab} = Q_{ab} + Q_{Xa}^T * X_b
    let mut f_tt_u = q.q_tt;
    let mut f_tv_u = q.q_tv;
    let mut f_vv_u = q.q_vv;
    for k in 0..ns {
        f_tt_u += q.q_xt[k] * x_t_u[k];
        f_tv_u += q.q_xv[k] * x_t_u[k];
        f_vv_u += q.q_xv[k] * x_v_u[k];
    }

    let mut f_tn_u = q.q_tn;
    let mut f_vn_u = q.q_vn;
    for i in 0..nc {
        let mut acc_tn = 0.0;
        let mut acc_vn = 0.0;
        for k in 0..ns {
            let qxn = q.q_xn[k * nc + i];
            acc_tn += qxn * x_t_u[k];
            acc_vn += qxn * x_v_u[k];
        }
        f_tn_u[i] += acc_tn;
        f_vn_u[i] += acc_vn;
    }

    let mut f_nn_u = q.q_nn;
    for i in 0..nc {
        for j in 0..nc {
            let mut acc = 0.0;
            for k in 0..ns {
                acc += q.q_xn[k * nc + i] * x_n_u[k * nc + j];
            }
            f_nn_u[i * nc + j] += acc;
        }
    }

    // Symmetrize minor floating-point asymmetry in F_nn.
    for i in 0..nc {
        for j in (i + 1)..nc {
            let avg = 0.5 * (f_nn_u[i * nc + j] + f_nn_u[j * nc + i]);
            f_nn_u[i * nc + j] = avg;
            f_nn_u[j * nc + i] = avg;
        }
    }

    // Exact Euler scaling from unit total moles to actual `sumn`.
    let out = AssociationHelmholtz {
        f: sumn * q.q,
        f_t: sumn * q.q_t,
        f_v: q.q_v,
        f_n: q.q_n,
        f_tt: sumn * f_tt_u,
        f_tv: f_tv_u,
        f_tn: f_tn_u,
        f_vv: f_vv_u / sumn,
        f_vn: f_vn_u.into_iter().map(|val| val / sumn).collect(),
        f_nn: f_nn_u.into_iter().map(|val| val / sumn).collect(),
        x_k,
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
            "non-finite association Helmholtz derivatives",
        ));
    }

    Ok(out)
}
