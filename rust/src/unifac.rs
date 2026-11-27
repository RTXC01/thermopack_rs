//! UNIFAC activity coefficient model with Flory-Huggins and Staverman-Guggenheim
//! combinatorial terms.
//! Ported from src/unifac.f90; the group database is generated into
//! `crate::unifac_data` from src/unifacdata.f90.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! The original model is evaluated with mole fractions, not molar amounts.
//! `ExcessGibbsMix` passes the composition vector straight through to
//! `Ge_UNIFAC_GH_SG`, so every quantity in this module is molar.

use crate::R as RGAS;
use crate::error::{Error, Result};
use crate::unifac_data::{UNIFAC_SUBGROUPS, lookup_unifac_component, lookup_unifac_interaction};

/// Combinatorial-term variant selected by the mixing-rule label.
///
/// The original `init_unifac` maps the label to two independent switches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnifacVariant {
    /// Label "UNIFAC": Flory-Huggins and Staverman-Guggenheim terms.
    Unifac,
    /// Label "UMR": Staverman-Guggenheim term only.
    Umr,
    /// Label "VTPR": no combinatorial term.
    Vtpr,
}

impl UnifacVariant {
    /// Resolve a source mixing-rule label. The original rejects any other
    /// label with a fatal error.
    pub fn from_label(label: &str) -> Result<Self> {
        match label {
            "UNIFAC" => Ok(Self::Unifac),
            "UMR" => Ok(Self::Umr),
            "VTPR" => Ok(Self::Vtpr),
            _ => Err(Error::InvalidInput("invalid unifac mixing rule")),
        }
    }

    /// Whether the Flory-Huggins combinatorial term is active.
    #[inline]
    pub fn flory_huggins(self) -> bool {
        matches!(self, Self::Unifac)
    }

    /// Whether the Staverman-Guggenheim combinatorial term is active.
    #[inline]
    pub fn staverman_guggenheim(self) -> bool {
        matches!(self, Self::Unifac | Self::Umr)
    }
}

/// Active UNIFAC parameter set for one mixture.
///
/// The original `init_unifac` restricts the model to the components present
/// in its own database and aborts otherwise; a missing component is an error
/// here as well.
#[derive(Debug, Clone)]
pub struct Unifac {
    /// Number of mixture components.
    pub nc: usize,
    /// Number of active sub-groups.
    pub ng: usize,
    /// Combinatorial-term variant.
    pub variant: UnifacVariant,
    /// Group counts of the active sub-groups, indexed as `i * ng + k`.
    vik: Vec<f64>,
    /// Attraction parameters of the active main groups, indexed as `j * ng + k`.
    ajk: Vec<f64>,
    /// Interaction parameters of the active main groups.
    bjk: Vec<f64>,
    /// Interaction parameters of the active main groups.
    cjk: Vec<f64>,
    /// Surface area parameters of the active sub-groups.
    qk: Vec<f64>,
    /// Surface area parameter of each component.
    qi: Vec<f64>,
    /// Volume parameter of each component.
    ri: Vec<f64>,
    /// Main group of each active sub-group.
    main_group_mapping: Vec<u8>,
    /// Main-group pairs without an interaction record. The source prints a
    /// message and uses zero parameters for these.
    missing_interactions: Vec<(u8, u8)>,
}

impl Unifac {
    /// Build the active parameter set for a component list.
    ///
    /// The component identifiers are the source `uid` values, for example
    /// "C2", "NC4" or "H2O".
    pub fn new(component_ids: &[&str], variant: UnifacVariant) -> Result<Self> {
        let nc = component_ids.len();
        if nc == 0 {
            return Err(Error::InvalidInput(
                "UNIFAC requires at least one component",
            ));
        }
        let mut records = Vec::with_capacity(nc);
        for id in component_ids {
            let record = lookup_unifac_component(id).ok_or_else(|| {
                Error::UnknownComponent(format!("UNIFAC component not found: {id}"))
            })?;
            records.push(record);
        }

        // The source counts how many sub-groups of each kind the components
        // use, then keeps the kinds that occur at least once.
        let mut active_groups: Vec<usize> = Vec::new();
        for kind in 0..UNIFAC_SUBGROUPS.len() {
            let count: u64 = records.iter().map(|r| r.v[kind] as u64).sum();
            if count > 0 {
                active_groups.push(kind);
            }
        }
        let ng = active_groups.len();

        let mut vik = vec![0.0; nc * ng];
        let mut qk = vec![0.0; ng];
        let mut ri = vec![0.0; nc];
        let mut qi = vec![0.0; nc];
        let mut main_group_mapping = vec![0u8; ng];
        for (k, &kind) in active_groups.iter().enumerate() {
            let parameter = &UNIFAC_SUBGROUPS[kind];
            qk[k] = parameter.qk;
            main_group_mapping[k] = parameter.main_group;
        }
        for i in 0..nc {
            for (k, &kind) in active_groups.iter().enumerate() {
                vik[i * ng + k] = f64::from(records[i].v[kind]);
            }
            for (k, &kind) in active_groups.iter().enumerate() {
                qi[i] += vik[i * ng + k] * qk[k];
                ri[i] += vik[i * ng + k] * UNIFAC_SUBGROUPS[kind].rk;
            }
        }

        let mut ajk = vec![0.0; ng * ng];
        let mut bjk = vec![0.0; ng * ng];
        let mut cjk = vec![0.0; ng * ng];
        let mut missing_interactions = Vec::new();
        for j in 0..ng {
            for k in 0..ng {
                let first = main_group_mapping[j];
                let second = main_group_mapping[k];
                match lookup_unifac_interaction(first, second) {
                    Some(record) => {
                        ajk[j * ng + k] = record.aij;
                        bjk[j * ng + k] = record.bij;
                        cjk[j * ng + k] = record.cij;
                    }
                    None => {
                        // The source reports the pair and uses zero parameters.
                        missing_interactions.push((first, second));
                    }
                }
            }
        }

        Ok(Self {
            nc,
            ng,
            variant,
            vik,
            ajk,
            bjk,
            cjk,
            qk,
            qi,
            ri,
            main_group_mapping,
            missing_interactions,
        })
    }

    /// Main-group pairs that the source database cannot resolve.
    pub fn missing_interactions(&self) -> &[(u8, u8)] {
        &self.missing_interactions
    }

    /// Main group of every active sub-group.
    pub fn main_group_mapping(&self) -> &[u8] {
        &self.main_group_mapping
    }

    /// Set main-group interaction parameters `(a_ij, b_ij, c_ij)` for all active
    /// sub-group pairs belonging to `(main_group_i, main_group_j)`, mirroring
    /// `unifac::setUNIFACgroupInteraction` without the Fortran `intent(out)`
    /// deallocation defect.
    pub fn set_group_interaction(
        &mut self,
        main_group_i: u8,
        main_group_j: u8,
        aij: f64,
        bij: f64,
        cij: f64,
    ) -> Result<()> {
        if !aij.is_finite() || !bij.is_finite() || !cij.is_finite() {
            return Err(Error::InvalidInput(
                "finite UNIFAC group interaction parameters required",
            ));
        }
        let ng = self.ng;
        for ii in 0..ng {
            if self.main_group_mapping[ii] != main_group_i {
                continue;
            }
            for jj in 0..ng {
                if self.main_group_mapping[jj] == main_group_j {
                    let idx = ii * ng + jj;
                    self.ajk[idx] = aij;
                    self.bjk[idx] = bij;
                    self.cjk[idx] = cij;
                }
            }
        }
        Ok(())
    }

    /// Retrieve main-group interaction parameters `(a_ij, b_ij, c_ij)` for an
    /// active pair of main groups, mirroring `unifac::getUNIFACgroupInteraction`.
    pub fn get_group_interaction(
        &self,
        main_group_i: u8,
        main_group_j: u8,
    ) -> Option<(f64, f64, f64)> {
        let ng = self.ng;
        for ii in 0..ng {
            if self.main_group_mapping[ii] != main_group_i {
                continue;
            }
            for jj in 0..ng {
                if self.main_group_mapping[jj] == main_group_j {
                    let idx = ii * ng + jj;
                    return Some((self.ajk[idx], self.bjk[idx], self.cjk[idx]));
                }
            }
        }
        None
    }

    /// Excess Gibbs energy and its analytical derivatives.
    ///
    /// All inputs and outputs are molar: the composition is a vector of mole
    /// fractions and the energies are per mole of mixture.
    ///
    /// The source-indexed loops below address several arrays of the same
    /// length in the same order; the iterator rewrite changes the summation
    /// order and therefore the rounding.
    #[allow(clippy::needless_range_loop)]
    pub fn excess_gibbs(
        &self,
        temperature: f64,
        mole_fractions: &[f64],
    ) -> Result<UnifacExcessGibbs> {
        let nc = self.nc;
        let ng = self.ng;
        if mole_fractions.len() != nc {
            return Err(Error::InvalidInput(
                "component count mismatch in UNIFAC excess Gibbs",
            ));
        }
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(Error::InvalidInput("UNIFAC temperature must be positive"));
        }
        let t = temperature;
        let n = mole_fractions;

        // Mixture energies and their temperature derivatives.
        let mut ejk = vec![0.0; ng * ng];
        let mut dejk_dt = vec![0.0; ng * ng];
        let mut d2ejk_dt2 = vec![0.0; ng * ng];
        for j in 0..ng {
            for k in 0..ng {
                let index = j * ng + k;
                let a = self.ajk[index];
                let b = self.bjk[index];
                let c = self.cjk[index];
                let e = (-(a / t + b + c * t)).exp();
                ejk[index] = e;
                let d = e * (a / (t * t) - c);
                dejk_dt[index] = d;
                d2ejk_dt2[index] = d * (a / (t * t) - c) - e * 2.0 * a / (t * t * t);
            }
        }

        // Some useful sums.
        let mut sum_vik_qk = vec![0.0; nc];
        for i in 0..nc {
            for k in 0..ng {
                sum_vik_qk[i] += self.vik[i * ng + k] * self.qk[k];
            }
        }
        let mut sum_nivik_qk = 0.0;
        for i in 0..nc {
            sum_nivik_qk += n[i] * sum_vik_qk[i];
        }

        // Component group fractions and the mixture group fractions.
        let mut theta_ij = vec![0.0; nc * ng];
        for i in 0..nc {
            for j in 0..ng {
                theta_ij[i * ng + j] = self.vik[i * ng + j] * self.qk[j] / sum_vik_qk[i];
            }
        }
        let mut sum_nlvlj = vec![0.0; ng];
        for j in 0..ng {
            let mut sum = 0.0;
            for i in 0..nc {
                sum += n[i] * self.vik[i * ng + j];
            }
            sum_nlvlj[j] = sum;
        }
        let mut theta_j = vec![0.0; ng];
        for j in 0..ng {
            theta_j[j] = sum_nlvlj[j] * self.qk[j] / sum_nivik_qk;
        }

        // Group and component group-activity sums.
        let mut lambda_k = vec![0.0; ng];
        for k in 0..ng {
            let mut sum = 0.0;
            for j in 0..ng {
                sum += theta_j[j] * ejk[j * ng + k];
            }
            lambda_k[k] = sum.ln();
        }
        let mut lambda_ik = vec![0.0; nc * ng];
        for i in 0..nc {
            for k in 0..ng {
                let mut sum = 0.0;
                for j in 0..ng {
                    sum += theta_ij[i * ng + j] * ejk[j * ng + k];
                }
                lambda_ik[i * ng + k] = sum.ln();
            }
        }

        // Excess Gibbs energy of the residual part.
        let mut sum_vq_lambda = vec![0.0; nc];
        for i in 0..nc {
            let mut sum = 0.0;
            for j in 0..ng {
                sum += self.vik[i * ng + j] * self.qk[j] * (lambda_k[j] - lambda_ik[i * ng + j]);
            }
            sum_vq_lambda[i] = sum;
        }
        let mut ge = 0.0;
        for i in 0..nc {
            ge -= n[i] * sum_vq_lambda[i];
        }

        // Differential of lambda_k with respect to n_i.
        let mut sum_vij_qj_ejk = vec![0.0; nc * ng];
        for i in 0..nc {
            for k in 0..ng {
                let mut sum = 0.0;
                for j in 0..ng {
                    sum += self.vik[i * ng + j] * self.qk[j] * ejk[j * ng + k];
                }
                sum_vij_qj_ejk[i * ng + k] = sum;
            }
        }
        let mut sum_nivij_qj_ejk = vec![0.0; ng];
        for k in 0..ng {
            let mut sum = 0.0;
            for i in 0..nc {
                sum += n[i] * sum_vij_qj_ejk[i * ng + k];
            }
            sum_nivij_qj_ejk[k] = sum;
        }
        let mut dlambda_k_dn = vec![0.0; nc * ng];
        for i in 0..nc {
            for k in 0..ng {
                dlambda_k_dn[i * ng + k] =
                    sum_vij_qj_ejk[i * ng + k] / sum_nivij_qj_ejk[k] - sum_vik_qk[i] / sum_nivik_qk;
            }
        }

        // dG/dn.
        let mut dg_dn = vec![0.0; nc];
        for i in 0..nc {
            let mut temp = 0.0;
            for j in 0..ng {
                temp += sum_nlvlj[j] * self.qk[j] * dlambda_k_dn[i * ng + j];
            }
            dg_dn[i] = -sum_vq_lambda[i] - temp;
        }

        // d2G/dn_i/dn_j.
        let mut d2g_dn_dn = vec![0.0; nc * nc];
        {
            let mut d2lambda_k_dn2 = vec![0.0; nc * nc * ng];
            for i in 0..nc {
                for j in 0..nc {
                    // The source expression is elementwise in the group index.
                    let second = sum_vik_qk[i] * sum_vik_qk[j] / (sum_nivik_qk * sum_nivik_qk);
                    for k in 0..ng {
                        let first = -sum_vij_qj_ejk[i * ng + k] * sum_vij_qj_ejk[j * ng + k]
                            / (sum_nivij_qj_ejk[k] * sum_nivij_qj_ejk[k]);
                        d2lambda_k_dn2[(i * nc + j) * ng + k] = first + second;
                    }
                }
            }
            for i in 0..nc {
                for j in 0..nc {
                    let mut sum_q_v_dlambda = 0.0;
                    for k in 0..ng {
                        sum_q_v_dlambda +=
                            self.qk[k] * self.vik[j * ng + k] * dlambda_k_dn[i * ng + k];
                    }
                    let mut sum_q_v_dlambda_swapped = 0.0;
                    for k in 0..ng {
                        sum_q_v_dlambda_swapped +=
                            self.qk[k] * self.vik[i * ng + k] * dlambda_k_dn[j * ng + k];
                    }
                    let mut temp = 0.0;
                    for k in 0..ng {
                        temp += sum_nlvlj[k] * d2lambda_k_dn2[(i * nc + j) * ng + k] * self.qk[k];
                    }
                    d2g_dn_dn[i * nc + j] = -(sum_q_v_dlambda + sum_q_v_dlambda_swapped) - temp;
                }
            }
        }

        // Temperature derivatives.
        let mut sum_vij_qj_dejk_dt = vec![0.0; nc * ng];
        let mut sum_vij_qj_d2ejk_dt2 = vec![0.0; nc * ng];
        for i in 0..nc {
            for k in 0..ng {
                let mut first = 0.0;
                let mut second = 0.0;
                for j in 0..ng {
                    first += self.vik[i * ng + j] * self.qk[j] * dejk_dt[j * ng + k];
                    second += self.vik[i * ng + j] * self.qk[j] * d2ejk_dt2[j * ng + k];
                }
                sum_vij_qj_dejk_dt[i * ng + k] = first;
                sum_vij_qj_d2ejk_dt2[i * ng + k] = second;
            }
        }
        let mut dlambda_ik_dt = vec![0.0; nc * ng];
        let mut d2lambda_ik_dt2 = vec![0.0; nc * ng];
        for i in 0..nc {
            for k in 0..ng {
                let first = sum_vij_qj_dejk_dt[i * ng + k] / sum_vij_qj_ejk[i * ng + k];
                dlambda_ik_dt[i * ng + k] = first;
                d2lambda_ik_dt2[i * ng + k] =
                    sum_vij_qj_d2ejk_dt2[i * ng + k] / sum_vij_qj_ejk[i * ng + k] - first * first;
            }
        }
        let mut sum_nivij_qj_dejk_dt = vec![0.0; ng];
        let mut sum_nivij_qj_d2ejk_dt2 = vec![0.0; ng];
        for k in 0..ng {
            let mut first = 0.0;
            let mut second = 0.0;
            for i in 0..nc {
                first += n[i] * sum_vij_qj_dejk_dt[i * ng + k];
                second += n[i] * sum_vij_qj_d2ejk_dt2[i * ng + k];
            }
            sum_nivij_qj_dejk_dt[k] = first;
            sum_nivij_qj_d2ejk_dt2[k] = second;
        }
        let mut dlambda_k_dt = vec![0.0; ng];
        let mut d2lambda_k_dt2 = vec![0.0; ng];
        for k in 0..ng {
            dlambda_k_dt[k] = sum_nivij_qj_dejk_dt[k] / sum_nivij_qj_ejk[k];
            d2lambda_k_dt2[k] =
                sum_nivij_qj_d2ejk_dt2[k] / sum_nivij_qj_ejk[k] - dlambda_k_dt[k] * dlambda_k_dt[k];
        }
        let mut sum_vij_qj_dldt_diff = vec![0.0; nc];
        for i in 0..nc {
            let mut sum = 0.0;
            for j in 0..ng {
                sum += self.vik[i * ng + j]
                    * self.qk[j]
                    * (dlambda_k_dt[j] - dlambda_ik_dt[i * ng + j]);
            }
            sum_vij_qj_dldt_diff[i] = sum;
        }

        // d2G/dn_i/dT.
        let mut d2g_dn_dt = vec![0.0; nc];
        {
            let mut d2lambda_k_dn_dt = vec![0.0; nc * ng];
            for i in 0..nc {
                for k in 0..ng {
                    d2lambda_k_dn_dt[i * ng + k] = sum_vij_qj_dejk_dt[i * ng + k]
                        / sum_nivij_qj_ejk[k]
                        - sum_vij_qj_ejk[i * ng + k] * sum_nivij_qj_dejk_dt[k]
                            / (sum_nivij_qj_ejk[k] * sum_nivij_qj_ejk[k]);
                }
            }
            for i in 0..nc {
                let mut temp = 0.0;
                for j in 0..ng {
                    temp += sum_nlvlj[j] * self.qk[j] * d2lambda_k_dn_dt[i * ng + j];
                }
                d2g_dn_dt[i] = -sum_vij_qj_dldt_diff[i] - temp;
            }
        }

        // dG/dT.
        let mut dg_dt = 0.0;
        for i in 0..nc {
            dg_dt -= n[i] * sum_vij_qj_dldt_diff[i];
        }

        // d2G/dT2.
        let mut d2g_dt2 = 0.0;
        for i in 0..nc {
            let mut sum = 0.0;
            for j in 0..ng {
                sum += self.vik[i * ng + j]
                    * self.qk[j]
                    * (d2lambda_k_dt2[j] - d2lambda_ik_dt2[i * ng + j]);
            }
            d2g_dt2 -= n[i] * sum;
        }

        // Final scaling of the residual part. The order follows the source:
        // the second temperature derivative uses the unscaled first one and
        // the mixed derivative uses the unscaled composition derivative.
        let unscaled_dg_dn = dg_dn.clone();
        let unscaled_d2g_dn_dt = d2g_dn_dt.clone();
        d2g_dt2 = RGAS * (2.0 * dg_dt + t * d2g_dt2);
        dg_dt = RGAS * (ge + dg_dt * t);
        for i in 0..nc {
            d2g_dn_dt[i] = RGAS * (unscaled_dg_dn[i] + unscaled_d2g_dn_dt[i] * t);
        }
        ge *= RGAS * t;
        dg_dn = dg_dn.iter().map(|value| value * RGAS * t).collect();
        d2g_dn_dn = d2g_dn_dn.iter().map(|value| value * RGAS * t).collect();

        // Combinatorial terms.
        if self.variant.flory_huggins() {
            let (ge_fh, dg_dn_fh, d2g_dn_dn_fh) = self.flory_huggins(n);
            ge += ge_fh * RGAS * t;
            dg_dt += ge_fh * RGAS;
            for i in 0..nc {
                dg_dn[i] += dg_dn_fh[i] * RGAS * t;
                d2g_dn_dt[i] += dg_dn_fh[i] * RGAS;
            }
            for i in 0..nc {
                for j in 0..nc {
                    d2g_dn_dn[i * nc + j] += d2g_dn_dn_fh[i * nc + j] * RGAS * t;
                }
            }
        }
        if self.variant.staverman_guggenheim() {
            let (ge_sg, dg_dn_sg, d2g_dn_dn_sg) = self.staverman_guggenheim(n);
            ge += ge_sg * RGAS * t;
            dg_dt += ge_sg * RGAS;
            for i in 0..nc {
                dg_dn[i] += dg_dn_sg[i] * RGAS * t;
                d2g_dn_dt[i] += dg_dn_sg[i] * RGAS;
            }
            for i in 0..nc {
                for j in 0..nc {
                    d2g_dn_dn[i * nc + j] += d2g_dn_dn_sg[i * nc + j] * RGAS * t;
                }
            }
        }

        let result = UnifacExcessGibbs {
            gibbs: ge,
            dg_dt,
            d2g_dt2,
            dg_dn,
            d2g_dn_dt,
            d2g_dn_dn,
        };
        if !result.is_finite() {
            return Err(Error::NotConverged("non-finite UNIFAC excess Gibbs energy"));
        }
        Ok(result)
    }

    /// Value-only molar excess Gibbs energy `g^E` (`[J/mol]`) without computing
    /// temperature or composition derivative matrices.
    #[allow(clippy::needless_range_loop)]
    pub(crate) fn excess_gibbs_value(
        &self,
        temperature: f64,
        mole_fractions: &[f64],
    ) -> Result<f64> {
        let nc = self.nc;
        let ng = self.ng;
        if mole_fractions.len() != nc {
            return Err(Error::InvalidInput(
                "component count mismatch in UNIFAC excess Gibbs",
            ));
        }
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(Error::InvalidInput("UNIFAC temperature must be positive"));
        }
        let t = temperature;
        let n = mole_fractions;

        let mut ejk = vec![0.0; ng * ng];
        for j in 0..ng {
            for k in 0..ng {
                let index = j * ng + k;
                let value = (-(self.ajk[index] / t + self.bjk[index] + self.cjk[index] * t)).exp();
                if !value.is_finite() {
                    return Err(Error::NotConverged("non-finite UNIFAC excess Gibbs energy"));
                }
                ejk[index] = value;
            }
        }

        let mut sum_vik_qk = vec![0.0; nc];
        let mut sum_nivik_qk = 0.0;
        for i in 0..nc {
            let mut s = 0.0;
            for k in 0..ng {
                s += self.vik[i * ng + k] * self.qk[k];
            }
            sum_vik_qk[i] = s;
            sum_nivik_qk += n[i] * s;
        }

        let mut theta_j = vec![0.0; ng];
        for j in 0..ng {
            let mut sum = 0.0;
            for i in 0..nc {
                sum += n[i] * self.vik[i * ng + j];
            }
            theta_j[j] = sum * self.qk[j] / sum_nivik_qk;
        }

        let mut lambda_k = vec![0.0; ng];
        for k in 0..ng {
            let mut sum = 0.0;
            for j in 0..ng {
                sum += theta_j[j] * ejk[j * ng + k];
            }
            lambda_k[k] = sum.ln();
        }

        let mut ge = 0.0;
        for i in 0..nc {
            let mut sum_vq_lambda = 0.0;
            for k in 0..ng {
                let mut sum = 0.0;
                for j in 0..ng {
                    let theta_ij = self.vik[i * ng + j] * self.qk[j] / sum_vik_qk[i];
                    sum += theta_ij * ejk[j * ng + k];
                }
                sum_vq_lambda += self.vik[i * ng + k] * self.qk[k] * (lambda_k[k] - sum.ln());
            }
            ge -= n[i] * sum_vq_lambda;
        }
        ge *= RGAS * t;

        if self.variant.flory_huggins() {
            let mut sumn = 0.0;
            let mut sumnr = 0.0;
            let mut sumlnr = 0.0;
            for i in 0..nc {
                sumn += n[i];
                sumnr += n[i] * self.ri[i];
                sumlnr += n[i] * self.ri[i].ln();
            }
            let ge_fh = sumlnr - sumn * sumnr.ln() + sumn * sumn.ln();
            ge += ge_fh * RGAS * t;
        }
        if self.variant.staverman_guggenheim() {
            const ZDIV2: f64 = 5.0;
            let mut sumnr = 0.0;
            let mut sumnq = 0.0;
            let mut sum = 0.0;
            for i in 0..nc {
                sumnr += n[i] * self.ri[i];
                sumnq += n[i] * self.qi[i];
                sum += n[i] * self.qi[i] * (self.qi[i] / self.ri[i]).ln();
            }
            let ge_sg = ZDIV2 * (sum - sumnq * sumnq.ln() + sumnq * sumnr.ln());
            ge += ge_sg * RGAS * t;
        }

        if !ge.is_finite() {
            return Err(Error::NotConverged("non-finite UNIFAC excess Gibbs energy"));
        }
        Ok(ge)
    }

    /// Flory-Huggins combinatorial term.
    ///
    /// The term is temperature independent, so it contributes no energy
    /// temperature derivative.
    #[allow(clippy::needless_range_loop)]
    fn flory_huggins(&self, n: &[f64]) -> (f64, Vec<f64>, Vec<f64>) {
        let nc = self.nc;
        let mut sumn = 0.0;
        for &value in n {
            sumn += value;
        }
        let mut sumnr = 0.0;
        for i in 0..nc {
            sumnr += n[i] * self.ri[i];
        }
        let mut lnr = vec![0.0; nc];
        let mut sumlnr = 0.0;
        for i in 0..nc {
            lnr[i] = self.ri[i].ln();
            sumlnr += n[i] * lnr[i];
        }
        let lnsumnr = sumnr.ln();
        let lnsumn = sumn.ln();
        let ge = sumlnr - sumn * lnsumnr + sumn * lnsumn;

        let mut dg_dn = vec![0.0; nc];
        for i in 0..nc {
            dg_dn[i] = lnr[i] - lnsumnr + lnsumn + 1.0 - sumn * self.ri[i] / sumnr;
        }
        let mut d2g_dn_dn = vec![0.0; nc * nc];
        for i in 0..nc {
            for j in i..nc {
                let value = -(self.ri[i] + self.ri[j]) / sumnr
                    + 1.0 / sumn
                    + sumn * self.ri[i] * self.ri[j] / (sumnr * sumnr);
                d2g_dn_dn[i * nc + j] = value;
                d2g_dn_dn[j * nc + i] = value;
            }
        }
        (ge, dg_dn, d2g_dn_dn)
    }

    /// Staverman-Guggenheim combinatorial term.
    ///
    /// The term is temperature independent, so it contributes no energy
    /// temperature derivative.
    #[allow(clippy::needless_range_loop)]
    fn staverman_guggenheim(&self, n: &[f64]) -> (f64, Vec<f64>, Vec<f64>) {
        const ZDIV2: f64 = 5.0;
        let nc = self.nc;
        let mut sumnr = 0.0;
        let mut sumnq = 0.0;
        for i in 0..nc {
            sumnr += n[i] * self.ri[i];
            sumnq += n[i] * self.qi[i];
        }
        let mut lnqor = vec![0.0; nc];
        let mut sum = 0.0;
        for i in 0..nc {
            lnqor[i] = (self.qi[i] / self.ri[i]).ln();
            sum += n[i] * self.qi[i] * lnqor[i];
        }
        let lnsumnr = sumnr.ln();
        let lnsumnq = sumnq.ln();
        let ge = ZDIV2 * (sum - sumnq * lnsumnq + sumnq * lnsumnr);

        let mut dg_dn = vec![0.0; nc];
        for i in 0..nc {
            dg_dn[i] = ZDIV2
                * self.qi[i]
                * (lnqor[i] - lnsumnq + lnsumnr - 1.0 + self.ri[i] * sumnq / (self.qi[i] * sumnr));
        }
        let mut d2g_dn_dn = vec![0.0; nc * nc];
        for i in 0..nc {
            for j in i..nc {
                let value = ZDIV2
                    * (-self.qi[i] * self.qi[j] / sumnq
                        + (self.qi[i] * self.ri[j] + self.qi[j] * self.ri[i]) / sumnr
                        - self.ri[i] * self.ri[j] * sumnq / (sumnr * sumnr));
                d2g_dn_dn[i * nc + j] = value;
                d2g_dn_dn[j * nc + i] = value;
            }
        }
        (ge, dg_dn, d2g_dn_dn)
    }
}

/// Excess Gibbs energy and its analytical derivatives for a UNIFAC mixture.
#[derive(Debug, Clone)]
pub struct UnifacExcessGibbs {
    /// Excess Gibbs energy [J/mol].
    pub gibbs: f64,
    /// First temperature derivative [J/(mol K)].
    pub dg_dt: f64,
    /// Second temperature derivative [J/(mol K2)].
    pub d2g_dt2: f64,
    /// Composition derivatives [J/mol], one per component.
    pub dg_dn: Vec<f64>,
    /// Mixed composition and temperature derivatives [J/(mol K)].
    pub d2g_dn_dt: Vec<f64>,
    /// Second composition derivatives [J/mol], indexed as `i * nc + j`.
    pub d2g_dn_dn: Vec<f64>,
}

impl UnifacExcessGibbs {
    /// Whether every reported quantity is finite.
    pub fn is_finite(&self) -> bool {
        self.gibbs.is_finite()
            && self.dg_dt.is_finite()
            && self.d2g_dt2.is_finite()
            && self.dg_dn.iter().all(|value| value.is_finite())
            && self.d2g_dn_dt.iter().all(|value| value.is_finite())
            && self.d2g_dn_dn.iter().all(|value| value.is_finite())
    }
}
