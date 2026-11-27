//! Excess Gibbs energy models and mixing rules (Huron-Vidal 1, Huron-Vidal 2, NRTL).
//! Ported from src/excess_gibbs.f90 and src/cbmix.f90.
//! Copyright thermotools and Contributors. Apache-2.0.

use crate::R as RGAS;
use crate::error::{Error, Result};
use crate::excess_gibbs_data::{GeCorrelation, lookup_excess_gibbs};

/// Mixing rule specification for cubic equations of state.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MixingRule {
    /// Classical van der Waals one-fluid mixing rules with binary interaction parameter kij.
    #[default]
    VdW,
    /// Huron-Vidal order 1 (HV1) with linear NRTL temperature dependence (tau = a/T + b).
    HuronVidal1,
    /// Huron-Vidal order 2 (HV2) with quadratic NRTL temperature dependence or Maribo-Mogensen form.
    HuronVidal2,
    /// Pure NRTL mixing rule.
    Nrtl,
    /// UNIFAC with the Flory-Huggins and Staverman-Guggenheim combinatorial terms.
    Unifac,
    /// UMR-PRU UNIFAC with the Staverman-Guggenheim combinatorial term only.
    Umr,
    /// VTPR UNIFAC residual-only group contribution mixing rule.
    Vtpr,
    /// The Wong-Sandler mixing rule of the cubic equation of state.
    WongSandler,
    /// Wong-Sandler mixing rule combined with the Huron-Vidal NRTL excess Gibbs formulation.
    HvWongSandler,
}

impl MixingRule {
    /// Resolve a mixing rule from its ThermoPack label or alias.
    pub fn from_name(name: &str) -> Result<Self> {
        match name.trim().to_ascii_uppercase().as_str() {
            "VDW" | "CLASSIC" => Ok(Self::VdW),
            "HV" | "HV1" | "HV0" | "HURON-VIDAL" | "HVCPA" => Ok(Self::HuronVidal1),
            "HV2" | "HURON-VIDAL2" | "HVCPA2" => Ok(Self::HuronVidal2),
            "NRTL" => Ok(Self::Nrtl),
            "UNIFAC" => Ok(Self::Unifac),
            "UMR" => Ok(Self::Umr),
            "VTPR" => Ok(Self::Vtpr),
            "WS" | "WONG-SANDLER" | "WONGSANDLER" | "WSCPA" => Ok(Self::WongSandler),
            "HVWONGSANDLER" | "HV-WS" | "HVWS" => Ok(Self::HvWongSandler),
            _ => Err(Error::InvalidInput("unknown cubic mixing rule")),
        }
    }

    /// Database code for the mixing rule.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::VdW => "VDW",
            Self::HuronVidal1 => "HV1",
            Self::HuronVidal2 => "HV2",
            Self::Nrtl => "NRTL",
            Self::Unifac => "UNIFAC",
            Self::Umr => "UMR",
            Self::Vtpr => "VTPR",
            Self::WongSandler => "WONG-SANDLER",
            Self::HvWongSandler => "HVWONGSANDLER",
        }
    }

    /// Whether this mixing rule requires excess Gibbs evaluations.
    #[inline]
    pub fn is_excess_gibbs(self) -> bool {
        matches!(
            self,
            Self::HuronVidal1
                | Self::HuronVidal2
                | Self::Nrtl
                | Self::Unifac
                | Self::Umr
                | Self::Vtpr
                | Self::WongSandler
                | Self::HvWongSandler
        )
    }

    /// Whether this mixing rule uses the UNIFAC group-contribution model.
    #[inline]
    pub fn is_unifac(self) -> bool {
        matches!(self, Self::Unifac | Self::Umr | Self::Vtpr)
    }

    /// Corresponding UNIFAC combinatorial variant, if this is a UNIFAC-family rule.
    #[inline]
    pub fn unifac_variant(self) -> Option<crate::unifac::UnifacVariant> {
        match self {
            Self::Unifac => Some(crate::unifac::UnifacVariant::Unifac),
            Self::Umr => Some(crate::unifac::UnifacVariant::Umr),
            Self::Vtpr => Some(crate::unifac::UnifacVariant::Vtpr),
            _ => None,
        }
    }

    /// Whether this mixing rule uses the Wong-Sandler mixing formulation.
    #[inline]
    pub fn is_wong_sandler(self) -> bool {
        matches!(self, Self::WongSandler | Self::HvWongSandler)
    }

    /// Whether this mixing rule uses covolume-weighted Huron-Vidal NRTL (`isHVmixModel`).
    #[inline]
    pub fn is_hv_model(self) -> bool {
        matches!(
            self,
            Self::HuronVidal1 | Self::HuronVidal2 | Self::HvWongSandler
        )
    }
}

/// The limit of the excess Gibbs model that a cubic equation uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeLimit {
    /// Infinite-pressure limit, used by the Huron-Vidal and NRTL models.
    Infinite,
    /// Zero-pressure limit, used by the UNIFAC mixing rule.
    Zero,
}

/// Limiting constant of the excess Gibbs model for a cubic equation.
///
/// The original `getInfinitLimitC` and `getZeroLimitC` select the constant
/// from the sub-equation of state. The UNIFAC mixing rule uses the
/// zero-pressure limit, every other excess Gibbs model uses the
/// infinite-pressure limit.
pub fn excess_gibbs_limit(kind: crate::cubic::CubicKind, limit: GeLimit) -> f64 {
    use crate::cubic::CubicKind;
    // The original `excess_gibbs::getInfinitLimitC` selects the constant of
    // the infinite pressure limit, `excess_gibbs::getZeroLimitC` the one of
    // the zero pressure limit.
    match limit {
        GeLimit::Infinite => match kind {
            CubicKind::VanDerWaals => 1.0,
            CubicKind::SoaveRedlichKwong | CubicKind::RedlichKwong => std::f64::consts::LN_2,
            CubicKind::PengRobinson => {
                let s2 = std::f64::consts::SQRT_2;
                -(s2 - 1.0).ln() / s2
            }
        },
        GeLimit::Zero => match kind {
            CubicKind::VanDerWaals => 0.55,
            CubicKind::SoaveRedlichKwong | CubicKind::RedlichKwong => 0.593,
            CubicKind::PengRobinson => 0.53087,
        },
    }
}

/// Binary interaction parameters for Excess Gibbs / Huron-Vidal.
#[derive(Debug, Clone, Copy)]
pub struct ExcessGibbsPair {
    pub kij: f64,
    pub has_record: bool,
    pub correlation: GeCorrelation,
    pub alpha_ij: f64,
    pub poly_ij: [f64; 3],
}

impl Default for ExcessGibbsPair {
    fn default() -> Self {
        Self {
            kij: 0.0,
            has_record: false,
            correlation: GeCorrelation::Polynomial,
            alpha_ij: 0.0,
            poly_ij: [0.0; 3],
        }
    }
}

/// Holds configured binary Excess Gibbs interaction data for a mixture.
#[derive(Debug, Clone)]
pub struct ExcessGibbsMixture {
    pub mrule: MixingRule,
    pub nc: usize,
    /// nc * nc matrix of pair parameters indexed as [i * nc + j] (interaction of i in j)
    pub pairs: Vec<ExcessGibbsPair>,
    /// Active UNIFAC parameter set, present for the UNIFAC mixing rule.
    pub unifac: Option<crate::unifac::Unifac>,
    /// The binary parameters of the Wong-Sandler mixing rule.
    pub wong_sandler: Option<crate::wong_sandler::WongSandlerData>,
}

impl ExcessGibbsMixture {
    /// Initialize excess Gibbs parameters for components using the database and optional reference.
    pub fn new(
        component_names: &[&str],
        mrule: MixingRule,
        eos_name: &str,
        preferred_ref: Option<&str>,
        kij_matrix: &[f64],
    ) -> Result<Self> {
        let nc = component_names.len();
        let mut pairs = vec![ExcessGibbsPair::default(); nc * nc];

        let unifac = if let Some(variant) = mrule.unifac_variant() {
            Some(crate::unifac::Unifac::new(component_names, variant)?)
        } else {
            None
        };

        if unifac.is_none() && !mrule.is_wong_sandler() {
            for i in 0..nc {
                for j in 0..nc {
                    let idx = i * nc + j;
                    if i == j {
                        continue;
                    }

                    let mut kij_val = if !kij_matrix.is_empty() {
                        kij_matrix[idx]
                    } else {
                        0.0
                    };

                    let name_i = component_names[i];
                    let name_j = component_names[j];

                    if let Some((rec, swapped)) =
                        lookup_excess_gibbs(eos_name, mrule.as_str(), name_i, name_j, preferred_ref)
                    {
                        kij_val = rec.kij;
                        let (alpha, poly) = if !swapped {
                            // (i, j) matches (uid1, uid2) in record
                            (rec.alpha_ij[0], rec.poly_ij)
                        } else {
                            // (i, j) matches (uid2, uid1) in record
                            (rec.alpha_ij[1], rec.poly_ji)
                        };

                        pairs[idx] = ExcessGibbsPair {
                            kij: kij_val,
                            has_record: true,
                            correlation: rec.correlation,
                            alpha_ij: alpha,
                            poly_ij: poly,
                        };
                    } else {
                        // Fallback to classic vdW with kij
                        pairs[idx] = ExcessGibbsPair {
                            kij: kij_val,
                            has_record: false,
                            correlation: GeCorrelation::Polynomial,
                            alpha_ij: 0.0,
                            poly_ij: [0.0; 3],
                        };
                    }
                }
            }
        } else if mrule == MixingRule::HvWongSandler && !kij_matrix.is_empty() {
            for i in 0..nc {
                for j in 0..nc {
                    if i != j {
                        pairs[i * nc + j].kij = kij_matrix[i * nc + j];
                    }
                }
            }
        }

        Ok(Self {
            mrule,
            nc,
            pairs,
            unifac,
            wong_sandler: None,
        })
    }
}

/// Output bundle of mixture cubic parameters and analytical derivatives from ExcessGibbsMix.
#[derive(Debug, Clone)]
pub struct ExcessGibbsCubicJet {
    pub a_mix: f64,
    pub da_dt: f64,
    pub d2a_dt2: f64,
    pub da_dni: Vec<f64>,
    pub d2a_dnid_t: Vec<f64>,
    pub d2a_dn_dn: Vec<f64>,
}

/// Excess Gibbs energy and derivatives consumed by the cubic mixing rule.
struct ExcessGibbsKernel {
    gibbs: f64,
    dg_dt: f64,
    d2g_dt2: f64,
    dg_dn: Vec<f64>,
    d2g_dn_dt: Vec<f64>,
    d2g_dn_dn: Vec<f64>,
}

/// The excess Gibbs energy of a mixture at the infinite pressure limit and
/// its derivatives.
///
/// The values are the ones of the original `excess_gibbs::gEinf`, which the
/// Wong-Sandler mixing rule consumes.
#[derive(Debug, Clone)]
pub struct ExcessGibbsJet {
    /// The excess Gibbs energy.
    pub gibbs: f64,
    /// The temperature derivative of the excess Gibbs energy.
    pub dg_dt: f64,
    /// The second temperature derivative of the excess Gibbs energy.
    pub d2g_dt2: f64,
    /// The derivatives with respect to the number of moles.
    pub dg_dn: Vec<f64>,
    /// The mixed temperature and mole derivatives.
    pub d2g_dn_dt: Vec<f64>,
    /// The second mole derivatives.
    pub d2g_dn_dn: Vec<f64>,
}

/// Computes the mixture cubic parameter a(T, x) and its exact derivatives via
/// Huron-Vidal, NRTL or UNIFAC.
/// Follows `ExcessGibbsMix` and `gExcess` from `src/excess_gibbs.f90`.
#[allow(clippy::too_many_arguments)]
pub fn excess_gibbs_cubic_mix(
    ge_mix: &ExcessGibbsMixture,
    kind: crate::cubic::CubicKind,
    temperature: f64,
    mole_fractions: &[f64],
    b_pure: &[f64],
    a_pure: &[f64],
    da_dt_pure: &[f64],
    d2a_dt2_pure: &[f64],
    b_mix: f64,
) -> Result<ExcessGibbsCubicJet> {
    excess_gibbs_cubic_mix_with_covolume(
        ge_mix,
        kind,
        temperature,
        mole_fractions,
        b_pure,
        a_pure,
        da_dt_pure,
        d2a_dt2_pure,
        b_mix,
        b_pure,
        &vec![0.; ge_mix.nc * ge_mix.nc],
    )
}

/// GE cubic coefficients with explicit first/second extensive covolume
/// derivatives. The ordinary linear-covolume entry point is retained above.
#[allow(clippy::too_many_arguments)]
pub fn excess_gibbs_cubic_mix_with_covolume(
    ge_mix: &ExcessGibbsMixture,
    kind: crate::cubic::CubicKind,
    temperature: f64,
    mole_fractions: &[f64],
    b_pure: &[f64],
    a_pure: &[f64],
    da_dt_pure: &[f64],
    d2a_dt2_pure: &[f64],
    b_mix: f64,
    b_partial: &[f64],
    b_second: &[f64],
) -> Result<ExcessGibbsCubicJet> {
    let nc = ge_mix.nc;
    if nc == 0 || mole_fractions.len() != nc || b_partial.len() != nc || b_second.len() != nc * nc {
        return Err(Error::InvalidInput(
            "component count mismatch in ExcessGibbs",
        ));
    }

    // The UNIFAC mixing rules (UNIFAC, UMR, VTPR) use the zero-pressure limit
    // of the cubic equation; every other excess Gibbs model uses the
    // infinite-pressure limit.
    let (kernel, lambda_ge) = if ge_mix.mrule.is_unifac() {
        let database = ge_mix.unifac.as_ref().ok_or(Error::InvalidInput(
            "UNIFAC mixing rule without a group database",
        ))?;
        let ge = database.excess_gibbs(temperature, mole_fractions)?;
        (
            ExcessGibbsKernel {
                gibbs: ge.gibbs,
                dg_dt: ge.dg_dt,
                d2g_dt2: ge.d2g_dt2,
                dg_dn: ge.dg_dn,
                d2g_dn_dt: ge.d2g_dn_dt,
                d2g_dn_dn: ge.d2g_dn_dn,
            },
            excess_gibbs_limit(kind, GeLimit::Zero),
        )
    } else {
        let lambda_ge = excess_gibbs_limit(kind, GeLimit::Infinite);
        (
            hv_nrtl_kernel(
                ge_mix,
                lambda_ge,
                temperature,
                mole_fractions,
                b_pure,
                a_pure,
                da_dt_pure,
                d2a_dt2_pure,
                true,
            )?,
            lambda_ge,
        )
    };

    let g_ex = kernel.gibbs;
    let dg_ex_dt = kernel.dg_dt;
    let d2g_ex_dt2 = kernel.d2g_dt2;
    let dg_ex_dni = kernel.dg_dn;
    let d2g_ex_dni_dt = kernel.d2g_dn_dt;
    let g_second = kernel.d2g_dn_dn;

    // Combine with cubic covolume B = b_mix to obtain a_mix and derivatives
    let mut sum_aa_b = 0.0;
    let mut sum_da_b = 0.0;
    let mut sum_d2a_b = 0.0;

    for i in 0..nc {
        let bi = b_pure[i];
        let zi = mole_fractions[i];
        sum_aa_b += zi * a_pure[i] / bi;
        sum_da_b += zi * da_dt_pure[i] / bi;
        sum_d2a_b += zi * d2a_dt2_pure[i] / bi;
    }

    let a_mix = b_mix * (sum_aa_b - g_ex / lambda_ge);
    let da_dt = b_mix * (sum_da_b - dg_ex_dt / lambda_ge);
    let d2a_dt2 = b_mix * (sum_d2a_b - d2g_ex_dt2 / lambda_ge);

    let mut da_dni = Vec::with_capacity(nc);
    let mut d2a_dnid_t = Vec::with_capacity(nc);

    for i in 0..nc {
        let bi = b_pure[i];
        // ai = (bi / B) * A + B * (a_i / b_i - dGex/dni / lambda_ge)
        let ai =
            (b_partial[i] / b_mix) * a_mix + b_mix * (a_pure[i] / bi - dg_ex_dni[i] / lambda_ge);
        let ait = (b_partial[i] / b_mix) * da_dt
            + b_mix * (da_dt_pure[i] / bi - d2g_ex_dni_dt[i] / lambda_ge);
        da_dni.push(ai);
        d2a_dnid_t.push(ait);
    }

    let mut d2a_dn_dn = vec![0.; nc * nc];
    for i in 0..nc {
        for j in 0..nc {
            d2a_dn_dn[i * nc + j] = (b_second[i * nc + j] * a_mix
                + b_partial[j] * da_dni[i]
                + b_partial[i] * da_dni[j])
                / b_mix
                - 2. * b_partial[i] * b_partial[j] * a_mix / (b_mix * b_mix)
                - b_mix * g_second[i * nc + j] / lambda_ge;
        }
    }
    Ok(ExcessGibbsCubicJet {
        a_mix,
        da_dt,
        d2a_dt2,
        da_dni,
        d2a_dnid_t,
        d2a_dn_dn,
    })
}

/// The excess Gibbs energy of a mixture at the infinite pressure limit and
/// its derivatives.
///
/// The original `excess_gibbs::gEinf` evaluates the excess Gibbs energy of
/// the Huron-Vidal, the NRTL and the Wong-Sandler models. The Wong-Sandler
/// mixing rule consumes the values.
#[allow(clippy::too_many_arguments)]
pub fn excess_gibbs_jet(
    ge_mix: &ExcessGibbsMixture,
    kind: crate::cubic::CubicKind,
    temperature: f64,
    mole_fractions: &[f64],
    b_pure: &[f64],
    a_pure: &[f64],
    da_dt_pure: &[f64],
    d2a_dt2_pure: &[f64],
) -> Result<ExcessGibbsJet> {
    let nc = ge_mix.nc;
    if nc == 0 || mole_fractions.len() != nc {
        return Err(Error::InvalidInput(
            "component count mismatch in ExcessGibbs",
        ));
    }
    if ge_mix.mrule.is_unifac() {
        return Err(Error::Unsupported(
            "the UNIFAC model does not deliver the excess Gibbs energy of the limit",
        ));
    }
    let kernel = hv_nrtl_kernel(
        ge_mix,
        excess_gibbs_limit(kind, GeLimit::Infinite),
        temperature,
        mole_fractions,
        b_pure,
        a_pure,
        da_dt_pure,
        d2a_dt2_pure,
        // The Wong-Sandler mixing rule consumes the second composition
        // derivatives of the excess Gibbs energy.
        true,
    )?;
    Ok(ExcessGibbsJet {
        gibbs: kernel.gibbs,
        dg_dt: kernel.dg_dt,
        d2g_dt2: kernel.d2g_dt2,
        dg_dn: kernel.dg_dn,
        d2g_dn_dt: kernel.d2g_dn_dt,
        d2g_dn_dn: kernel.d2g_dn_dn,
    })
}

/// Huron-Vidal and NRTL excess Gibbs energy with analytical derivatives.
/// Follows `gEinf` from `src/excess_gibbs.f90`.
#[allow(clippy::too_many_arguments)]
fn hv_nrtl_kernel(
    ge_mix: &ExcessGibbsMixture,
    lambda_ge: f64,
    temperature: f64,
    mole_fractions: &[f64],
    b_pure: &[f64],
    a_pure: &[f64],
    da_dt_pure: &[f64],
    d2a_dt2_pure: &[f64],
    second_composition: bool,
) -> Result<ExcessGibbsKernel> {
    let nc = ge_mix.nc;
    let t = temperature;

    // The five nc x nc matrices and the component vectors share two
    // allocations, so the preparation of a state does not allocate once per
    // matrix.
    let mut matrix_buffer: Vec<f64> = vec![0.0; 5 * nc * nc];
    let (rest, alpha) = matrix_buffer.split_at_mut(4 * nc * nc);
    let (rest, cij) = rest.split_at_mut(3 * nc * nc);
    let (rest, d2taudt2) = rest.split_at_mut(2 * nc * nc);
    let (rest, dtaudt) = rest.split_at_mut(nc * nc);
    let tau = rest;
    for value in cij.iter_mut() {
        *value = 1.0;
    }

    // 1. Calculate temperature-dependent tau_ji and C_ji
    if ge_mix.mrule == MixingRule::WongSandler {
        // The Wong-Sandler parameters come from the user, as in the original.
        // `getWSParam` fills `tau(r, c) = mixWS%f_tauij(r, c)` directly, and
        // `gEinf` sums `tau(j, i) * z_j * Cij(j, i)` for each component `i`.
        let ws = ge_mix
            .wong_sandler
            .as_ref()
            .ok_or(crate::Error::InvalidInput(
                "the Wong-Sandler mixing rule without binary parameters",
            ))?;
        for i in 0..nc {
            for j in 0..nc {
                if i == j {
                    continue;
                }
                let idx_ji = j * nc + i;
                let exponent = ws.alphaij[idx_ji];
                if exponent.abs() > 1.0e-20 {
                    let (value, derivative, second) = ws.tauij[idx_ji].evaluate(t);
                    tau[idx_ji] = value;
                    dtaudt[idx_ji] = derivative;
                    d2taudt2[idx_ji] = second;
                    alpha[idx_ji] = exponent;
                } else {
                    // The original reverts to a vanishing interaction.
                    tau[idx_ji] = 0.0;
                    dtaudt[idx_ji] = 0.0;
                    d2taudt2[idx_ji] = 0.0;
                    alpha[idx_ji] = 1.0;
                }
                cij[idx_ji] = (-alpha[idx_ji] * tau[idx_ji]).exp();
            }
        }
    } else {
        for i in 0..nc {
            let bi = b_pure[i];
            for j in 0..nc {
                if i == j {
                    continue;
                }
                let idx_ji = j * nc + i; // interaction of j in i
                let bj = b_pure[j];
                let pair = &ge_mix.pairs[idx_ji];

                if pair.has_record {
                    alpha[idx_ji] = pair.alpha_ij;
                    match pair.correlation {
                        GeCorrelation::Polynomial => {
                            let p = pair.poly_ij;
                            if ge_mix.mrule == MixingRule::HuronVidal1
                                || ge_mix.mrule == MixingRule::Nrtl
                            {
                                tau[idx_ji] = p[0] / t + p[1];
                                dtaudt[idx_ji] = -p[0] / (t * t);
                                d2taudt2[idx_ji] = 2.0 * p[0] / (t * t * t);
                            } else {
                                // HV2 2nd order polynomial: a/T + b + cT
                                tau[idx_ji] = p[0] / t + p[1] + p[2] * t;
                                dtaudt[idx_ji] = -p[0] / (t * t) + p[2];
                                d2taudt2[idx_ji] = 2.0 * p[0] / (t * t * t);
                            }
                        }
                        GeCorrelation::MariboMogensen => {
                            let p = pair.poly_ij;
                            let tc = p[2];
                            let t_ref = 298.15;
                            let d1 = 1.0 - t / tc;
                            let d2 = 1.0 - t_ref / tc;
                            let diff_term = d1 * d1 - d2 * d2;
                            let tau_val = (p[0] + p[1] * diff_term) / t;
                            let dtau_val = -tau_val / t - 2.0 * p[1] * (1.0 - t / tc) / (tc * t);
                            let d2tau_val = -2.0 * (dtau_val - p[1] / (tc * tc)) / t;

                            tau[idx_ji] = tau_val;
                            dtaudt[idx_ji] = dtau_val;
                            d2taudt2[idx_ji] = d2tau_val;
                        }
                    }
                    cij[idx_ji] = (-alpha[idx_ji] * tau[idx_ji]).exp();
                } else {
                    // Revert to classic vdW mixing rule matching
                    let kij = pair.kij;
                    let t1 = -2.0 * (bi * bj).sqrt() / (bi + bj);
                    let gii = -a_pure[i] * lambda_ge / bi;
                    // a_pure already contains alpha(T); only the covolumes and
                    // model constant belong in the temperature-independent factor.
                    let t2 = lambda_ge / bi;
                    let t3 = t2 * lambda_ge / bj;
                    let t4 = (t3 * a_pure[i] * a_pure[j]).sqrt();

                    tau[idx_ji] = (t1 * t4 * (1.0 - kij) - gii) / (RGAS * t);

                    let dt4_dt = if t4 > 1e-30 {
                        t3 * (a_pure[i] * da_dt_pure[j] + a_pure[j] * da_dt_pure[i]) / (2.0 * t4)
                    } else {
                        0.0
                    };

                    let d2t4_dt2 = if t4 > 1e-30 {
                        t3 * (a_pure[i] * d2a_dt2_pure[j]
                            + 2.0 * da_dt_pure[i] * da_dt_pure[j]
                            + a_pure[j] * d2a_dt2_pure[i])
                            / (2.0 * t4)
                            - dt4_dt * dt4_dt / t4
                    } else {
                        0.0
                    };

                    let dgii_dt = -t2 * da_dt_pure[i];
                    let d2gii_dt2 = -t2 * d2a_dt2_pure[i];

                    dtaudt[idx_ji] = ((1.0 - kij) * t1 * dt4_dt - dgii_dt) / (RGAS * t)
                        - ((1.0 - kij) * t1 * t4 - gii) / (RGAS * t * t);

                    d2taudt2[idx_ji] = ((1.0 - kij) * t1 * d2t4_dt2 - d2gii_dt2) / (RGAS * t)
                        - 2.0 * ((1.0 - kij) * t1 * dt4_dt - dgii_dt) / (RGAS * t * t)
                        + 2.0 * ((1.0 - kij) * t1 * t4 - gii) / (RGAS * t * t * t);

                    cij[idx_ji] = if ge_mix.mrule == MixingRule::Nrtl {
                        // Unlike fitted dimensionless C_ji, the source's
                        // no-record NRTL fallback stores b_j in L/mol. The
                        // diagonal stays one, so this unit cannot cancel.
                        1000.0 * bj
                    } else {
                        1.0
                    };
                }
            }
        }
    }

    // 2. Inner sums gamma1 .. gamma6 and fracin1 for each component
    let mut vector_buffer = vec![0.0; 9 * nc];
    let (rest, d2g_ex_dni_dt) = vector_buffer.split_at_mut(8 * nc);
    let (rest, dg_ex_dni) = rest.split_at_mut(7 * nc);
    let (rest, fracin1) = rest.split_at_mut(6 * nc);
    let (rest, gamma6) = rest.split_at_mut(5 * nc);
    let (rest, gamma5) = rest.split_at_mut(4 * nc);
    let (rest, gamma4) = rest.split_at_mut(3 * nc);
    let (rest, gamma3) = rest.split_at_mut(2 * nc);
    let (rest, gamma2) = rest.split_at_mut(nc);
    let gamma1 = rest;

    for i in 0..nc {
        for j in 0..nc {
            let idx_ji = j * nc + i;
            let bj = if matches!(ge_mix.mrule, MixingRule::Nrtl | MixingRule::WongSandler) {
                1.0
            } else {
                b_pure[j]
            };
            let zj = mole_fractions[j];
            let c_ji = cij[idx_ji];
            let tau_ji = tau[idx_ji];
            let a_ji = alpha[idx_ji];
            let dt_ji = dtaudt[idx_ji];
            let d2t_ji = d2taudt2[idx_ji];

            let term_base = bj * zj * c_ji;
            gamma1[i] += tau_ji * term_base;
            gamma2[i] += term_base;
            gamma3[i] += dt_ji * term_base * a_ji;
            gamma4[i] += dt_ji * term_base * (1.0 - tau_ji * a_ji);
            gamma5[i] += term_base * a_ji * (d2t_ji - a_ji * dt_ji * dt_ji);
            gamma6[i] += term_base
                * (d2t_ji * (1.0 - a_ji * tau_ji) - dt_ji * dt_ji * a_ji * (2.0 - tau_ji * a_ji));
        }

        if gamma2[i].abs() > 1e-30 {
            fracin1[i] = gamma1[i] * gamma3[i] / (gamma2[i] * gamma2[i]);
        }
    }

    // 3. Excess Gibbs energy and T derivatives
    let mut sum2i = 0.0;
    let mut sum3i = 0.0;
    let mut sum4i = 0.0;
    let mut sum5i = 0.0;
    let mut sum6i = 0.0;
    let mut sum7i = 0.0;
    let mut sum8i = 0.0;

    for i in 0..nc {
        let zi = mole_fractions[i];
        let g2 = gamma2[i];
        if g2.abs() < 1e-30 {
            continue;
        }
        sum2i += zi * gamma1[i] / g2;
        sum3i += zi * fracin1[i];
        sum4i += zi * gamma4[i] / g2;
        sum5i += zi * (fracin1[i] + gamma4[i] / g2);
        sum6i += zi * (2.0 * fracin1[i] * gamma3[i] / g2);
        sum7i += zi * (2.0 * gamma4[i] * gamma3[i] + gamma1[i] * gamma5[i]) / (g2 * g2);
        sum8i += zi * gamma6[i] / g2;
    }

    let g_ex = RGAS * t * sum2i;
    let dg_ex_dt = RGAS * (sum2i + t * (sum3i + sum4i));
    let d2g_ex_dt2 = RGAS * (2.0 * sum5i + t * (sum6i + sum7i + sum8i));

    // 4. Composition derivatives d(G_ex)/dni and d2(G_ex)/dni_dt
    for i in 0..nc {
        let mut sum1ji = 0.0;
        let mut sum2ji = 0.0;
        let mut sum3ji = 0.0;
        let mut sum4ji = 0.0;
        let mut sum5ji = 0.0;
        let mut sum6ji = 0.0;
        let mut sum7ji = 0.0;

        let bi = if matches!(ge_mix.mrule, MixingRule::Nrtl | MixingRule::WongSandler) {
            1.0
        } else {
            b_pure[i]
        };

        for j in 0..nc {
            let idx_ij = i * nc + j; // interaction of i in j
            let zj = mole_fractions[j];
            let g2j = gamma2[j];
            if g2j.abs() < 1e-30 {
                continue;
            }
            let c_ij = cij[idx_ij];
            let tau_ij = tau[idx_ij];
            let a_ij = alpha[idx_ij];
            let dt_ij = dtaudt[idx_ij];

            sum1ji += zj * bi * c_ij * tau_ij / g2j;
            sum2ji += zj * bi * c_ij * gamma1[j] / (g2j * g2j);
            sum3ji += zj * bi * c_ij * dt_ij * (1.0 - a_ij * tau_ij) / g2j;
            sum4ji += zj * bi * c_ij * tau_ij * gamma3[j] / (g2j * g2j);
            sum5ji += zj * bi * c_ij * dt_ij * a_ij * gamma1[j] / (g2j * g2j);
            sum6ji += zj * bi * c_ij * gamma4[j] / (g2j * g2j);
            sum7ji += zj * bi * c_ij * fracin1[j] / g2j;
        }

        let term_self = if gamma2[i].abs() > 1e-30 {
            gamma1[i] / gamma2[i]
        } else {
            0.0
        };
        dg_ex_dni[i] = RGAS * t * (term_self + sum1ji - sum2ji);

        let term_self_dt = if gamma2[i].abs() > 1e-30 {
            fracin1[i] + gamma4[i] / gamma2[i]
        } else {
            0.0
        };
        d2g_ex_dni_dt[i] = (1.0 / t) * dg_ex_dni[i]
            + RGAS * t * (term_self_dt + sum3ji + sum4ji + sum5ji - sum6ji - 2.0 * sum7ji);
    }

    // 5. Second composition derivatives of the extensive excess Gibbs energy.
    let mut d2g_ex_dnid_nj = Vec::new();
    if second_composition {
        d2g_ex_dnid_nj = vec![0.0; nc * nc];
        for i in 0..nc {
            let bi = if matches!(ge_mix.mrule, MixingRule::Nrtl | MixingRule::WongSandler) {
                1.0
            } else {
                b_pure[i]
            };
            for j in 0..nc {
                let bj = if matches!(ge_mix.mrule, MixingRule::Nrtl | MixingRule::WongSandler) {
                    1.0
                } else {
                    b_pure[j]
                };
                let mut sum_i = 0.0;
                let mut sum_j = 0.0;
                let mut sum_ij = 0.0;
                for ii in 0..nc {
                    let gamma2_ii = gamma2[ii] * gamma2[ii];
                    let c_ii_i = cij[i * nc + ii];
                    let c_ii_j = cij[j * nc + ii];
                    let z_ii = mole_fractions[ii];
                    sum_i += z_ii * (bj * c_ii_j) * tau[i * nc + ii] * bi * c_ii_i / gamma2_ii;
                    sum_j += z_ii * (bi * c_ii_i) * tau[j * nc + ii] * bj * c_ii_j / gamma2_ii;
                    sum_ij +=
                        z_ii * bi * bj * c_ii_i * c_ii_j * gamma1[ii] / (gamma2_ii * gamma2[ii]);
                }
                let g2i = gamma2[i] * gamma2[i];
                let g2j = gamma2[j] * gamma2[j];
                d2g_ex_dnid_nj[i * nc + j] = RGAS
                    * t
                    * (tau[j * nc + i] * bj * cij[j * nc + i] / gamma2[i]
                        - bj * cij[j * nc + i] * gamma1[i] / g2i
                        + tau[i * nc + j] * bi * cij[i * nc + j] / gamma2[j]
                        - sum_i
                        - bi * cij[i * nc + j] * gamma1[j] / g2j
                        - sum_j
                        + 2.0 * sum_ij);
            }
        }
    }

    Ok(ExcessGibbsKernel {
        gibbs: g_ex,
        dg_dt: dg_ex_dt,
        d2g_dt2: d2g_ex_dt2,
        dg_dn: dg_ex_dni.to_vec(),
        d2g_dn_dt: d2g_ex_dni_dt.to_vec(),
        d2g_dn_dn: d2g_ex_dnid_nj,
    })
}

/// Value-only Huron-Vidal / NRTL / Wong-Sandler excess Gibbs energy `g^E` (`[J/mol]`)
/// with zero matrix or derivative allocations.
fn hv_nrtl_value(
    ge_mix: &ExcessGibbsMixture,
    lambda_ge: f64,
    temperature: f64,
    mole_fractions: &[f64],
    b_pure: &[f64],
    a_pure: &[f64],
) -> Result<f64> {
    let nc = ge_mix.nc;
    let t = temperature;
    let ws = if ge_mix.mrule == MixingRule::WongSandler {
        Some(
            ge_mix
                .wong_sandler
                .as_ref()
                .ok_or(crate::Error::InvalidInput(
                    "the Wong-Sandler mixing rule without binary parameters",
                ))?,
        )
    } else {
        None
    };

    let mut sum2i = 0.0;
    for i in 0..nc {
        let bi = b_pure[i];
        let mut gamma1 = 0.0;
        let mut gamma2 = 0.0;
        for j in 0..nc {
            let (tau_ji, c_ji) = if i == j {
                (0.0, 1.0)
            } else if let Some(ws) = ws {
                let idx_ji = j * nc + i;
                let exponent = ws.alphaij[idx_ji];
                if exponent.abs() > 1.0e-20 {
                    let (value, _, _) = ws.tauij[idx_ji].evaluate(t);
                    (value, (-exponent * value).exp())
                } else {
                    (0.0, 1.0)
                }
            } else {
                let idx_ji = j * nc + i;
                let bj = b_pure[j];
                let pair = &ge_mix.pairs[idx_ji];
                if pair.has_record {
                    let tau_val = match pair.correlation {
                        GeCorrelation::Polynomial => {
                            let p = pair.poly_ij;
                            if ge_mix.mrule == MixingRule::HuronVidal1
                                || ge_mix.mrule == MixingRule::Nrtl
                            {
                                p[0] / t + p[1]
                            } else {
                                p[0] / t + p[1] + p[2] * t
                            }
                        }
                        GeCorrelation::MariboMogensen => {
                            let p = pair.poly_ij;
                            let tc = p[2];
                            let t_ref = 298.15;
                            let d1 = 1.0 - t / tc;
                            let d2 = 1.0 - t_ref / tc;
                            (p[0] + p[1] * (d1 * d1 - d2 * d2)) / t
                        }
                    };
                    (tau_val, (-pair.alpha_ij * tau_val).exp())
                } else {
                    let kij = pair.kij;
                    let t1 = -2.0 * (bi * bj).sqrt() / (bi + bj);
                    let gii = -a_pure[i] * lambda_ge / bi;
                    let t2 = lambda_ge / bi;
                    let t3 = t2 * lambda_ge / bj;
                    let t4 = (t3 * a_pure[i] * a_pure[j]).sqrt();
                    let tau_val = (t1 * t4 * (1.0 - kij) - gii) / (RGAS * t);
                    let c_val = if ge_mix.mrule == MixingRule::Nrtl {
                        1000.0 * bj
                    } else {
                        1.0
                    };
                    (tau_val, c_val)
                }
            };
            let weight_b = if matches!(ge_mix.mrule, MixingRule::Nrtl | MixingRule::WongSandler) {
                1.0
            } else {
                b_pure[j]
            };
            let term_base = weight_b * mole_fractions[j] * c_ji;
            gamma1 += tau_ji * term_base;
            gamma2 += term_base;
        }
        if !gamma1.is_finite() || !gamma2.is_finite() {
            return Err(Error::InvalidInput("temperature outside numerical range"));
        }
        if gamma2.abs() >= 1e-30 {
            sum2i += mole_fractions[i] * gamma1 / gamma2;
        }
    }
    let ge = RGAS * t * sum2i;
    if !ge.is_finite() {
        return Err(Error::InvalidInput("temperature outside numerical range"));
    }
    Ok(ge)
}

/// Value-only infinite-pressure excess Gibbs energy `g^E` (`[J/mol]`) for Wong-Sandler
/// without computing first or second derivative matrices.
pub(crate) fn excess_gibbs_value(
    ge_mix: &ExcessGibbsMixture,
    kind: crate::cubic::CubicKind,
    temperature: f64,
    mole_fractions: &[f64],
    b_pure: &[f64],
    a_pure: &[f64],
) -> Result<f64> {
    let nc = ge_mix.nc;
    if nc == 0 || mole_fractions.len() != nc {
        return Err(Error::InvalidInput(
            "component count mismatch in ExcessGibbs",
        ));
    }
    if ge_mix.mrule.is_unifac() {
        return Err(Error::Unsupported(
            "the UNIFAC model does not deliver the excess Gibbs energy of the limit",
        ));
    }
    hv_nrtl_value(
        ge_mix,
        excess_gibbs_limit(kind, GeLimit::Infinite),
        temperature,
        mole_fractions,
        b_pure,
        a_pure,
    )
}

/// Value-only cubic mixture attraction `a_mix` (`[Pa m^6 / mol^2]`) for excess-Gibbs
/// mixing rules (`HV1`, `HV2`, `NRTL`, `UNIFAC`) without computing derivative Hessians.
pub(crate) fn excess_gibbs_cubic_a_value(
    ge_mix: &ExcessGibbsMixture,
    kind: crate::cubic::CubicKind,
    temperature: f64,
    mole_fractions: &[f64],
    b_pure: &[f64],
    a_pure: &[f64],
    b_mix: f64,
) -> Result<f64> {
    let nc = ge_mix.nc;
    if nc == 0 || mole_fractions.len() != nc {
        return Err(Error::InvalidInput(
            "component count mismatch in ExcessGibbs",
        ));
    }
    let (g_ex, lambda_ge) = if ge_mix.mrule.is_unifac() {
        let database = ge_mix.unifac.as_ref().ok_or(Error::InvalidInput(
            "UNIFAC mixing rule without a group database",
        ))?;
        (
            database.excess_gibbs_value(temperature, mole_fractions)?,
            excess_gibbs_limit(kind, GeLimit::Zero),
        )
    } else {
        let lambda_ge = excess_gibbs_limit(kind, GeLimit::Infinite);
        (
            hv_nrtl_value(
                ge_mix,
                lambda_ge,
                temperature,
                mole_fractions,
                b_pure,
                a_pure,
            )?,
            lambda_ge,
        )
    };
    let mut sum_aa_b = 0.0;
    for i in 0..nc {
        sum_aa_b += mole_fractions[i] * a_pure[i] / b_pure[i];
    }
    let a = b_mix * (sum_aa_b - g_ex / lambda_ge);
    if !a.is_finite() {
        return Err(Error::InvalidInput("temperature outside numerical range"));
    }
    Ok(a)
}
