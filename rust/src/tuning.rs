//! Parameter tuning helpers mirroring `src/tuning.f90`.
//!
//! Provides structured parameter types and typed getters/setters on [`Cubic`]
//! for Huron-Vidal (`HV`), Wong-Sandler (`WS`), Twu and Mathias-Copeman alpha
//! correlations, named alpha/beta correlation setters, binary `kij`/`lij`,
//! volume-shift polynomials, UNIFAC main-group interactions, and critical
//! parameter redefinition.

use crate::cubic::Cubic;
use crate::error::{Error, Result, positive};
use crate::excess_gibbs_data::GeCorrelation;
use crate::translation::{TranslationKind, VolumeTranslation};
use crate::wong_sandler::{Fraction, WongSandlerData};
use crate::{AlphaCorrelation, Beta, BetaCorrelation, R};

/// Pairwise Huron-Vidal / NRTL excess Gibbs interaction parameters
/// (`thermopack_getHVparam` / `thermopack_setHVparam`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HvPairParams {
    /// Non-randomness factor `alpha_ij`.
    pub alpha_ij: f64,
    /// Non-randomness factor `alpha_ji`.
    pub alpha_ji: f64,
    /// Linear-inverse-T energy coefficient `aGE_ij` (`[K]`).
    pub a_ge_ij: f64,
    /// Linear-inverse-T energy coefficient `aGE_ji` (`[K]`).
    pub a_ge_ji: f64,
    /// Constant term `bGE_ij` (`[-]`).
    pub b_ge_ij: f64,
    /// Constant term `bGE_ji` (`[-]`).
    pub b_ge_ji: f64,
    /// Linear-T term `cGE_ij` (`[1/K]`, used by HV2 and HVWongSandler).
    pub c_ge_ij: f64,
    /// Linear-T term `cGE_ji` (`[1/K]`, used by HV2 and HVWongSandler).
    pub c_ge_ji: f64,
}

/// Pairwise constant Wong-Sandler parameters
/// (`thermopack_getWSparam` / `thermopack_setWSparam`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WsPairParams {
    /// NRTL non-randomness factor `alpha_ij`.
    pub alpha_ij: f64,
    /// NRTL non-randomness factor `alpha_ji`.
    pub alpha_ji: f64,
    /// Second-virial cross interaction `k_ij`.
    pub k_ij: f64,
    /// Second-virial cross interaction `k_ji`.
    pub k_ji: f64,
    /// Constant NRTL interaction energy parameter `tau_ij`.
    pub tau_ij: f64,
    /// Constant NRTL interaction energy parameter `tau_ji`.
    pub tau_ji: f64,
}

/// Volume translation polynomial parameters mirroring
/// `thermopack_get_volume_shift_parameters` / `thermopack_set_volume_shift_parameters`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VolumeShiftParams {
    /// Polynomial coefficients `[ciA, ciB, ciC, ciD, ciE, ciF]` in SI units
    /// (`m^3/mol`, `m^3/(mol K)`, ...).
    pub coefficients: [f64; 6],
    /// Polynomial degree/type (`Constant = 1`, `Linear = 2`, `Quadratic = 3`, `Quintic = 6`).
    pub ci_type: TranslationKind,
}

impl From<VolumeTranslation> for VolumeShiftParams {
    fn from(value: VolumeTranslation) -> Self {
        Self {
            coefficients: value.coefficients,
            ci_type: value.kind,
        }
    }
}

impl TryFrom<VolumeShiftParams> for VolumeTranslation {
    type Error = Error;

    fn try_from(value: VolumeShiftParams) -> Result<Self> {
        VolumeTranslation::new(value.ci_type, value.coefficients)
    }
}

impl Cubic {
    /// Retrieve Huron-Vidal / NRTL binary interaction parameters for pair `(i, j)`
    /// (`thermopack_getHVparam`).
    pub fn hv_param(&self, i: usize, j: usize) -> Result<HvPairParams> {
        self.validate_pair_indices(i, j)?;
        let ge = self
            .excess_gibbs
            .as_ref()
            .filter(|ge| !ge.mrule.is_unifac() && ge.mrule != crate::MixingRule::WongSandler)
            .ok_or(Error::Unsupported(
                "Huron-Vidal / NRTL parameters are not allocated for this mixing rule",
            ))?;
        let nc = self.components().len();
        let ij = &ge.pairs[i * nc + j];
        let ji = &ge.pairs[j * nc + i];
        Ok(HvPairParams {
            alpha_ij: ij.alpha_ij,
            alpha_ji: ji.alpha_ij,
            a_ge_ij: ij.poly_ij[0],
            a_ge_ji: ji.poly_ij[0],
            b_ge_ij: ij.poly_ij[1],
            b_ge_ji: ji.poly_ij[1],
            c_ge_ij: ij.poly_ij[2],
            c_ge_ji: ji.poly_ij[2],
        })
    }

    /// Tune Huron-Vidal / NRTL binary interaction parameters for pair `(i, j)`
    /// (`thermopack_setHVparam`). Marks both `(i, j)` and `(j, i)` as active
    /// polynomial correlations (`correlation = 1`).
    pub fn set_hv_param(&mut self, i: usize, j: usize, params: HvPairParams) -> Result<()> {
        self.validate_pair_indices(i, j)?;
        let vals = [
            params.alpha_ij,
            params.alpha_ji,
            params.a_ge_ij,
            params.a_ge_ji,
            params.b_ge_ij,
            params.b_ge_ji,
            params.c_ge_ij,
            params.c_ge_ji,
        ];
        if vals.iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput("finite HV parameters required"));
        }
        let nc = self.components().len();
        let ge = self
            .excess_gibbs
            .as_mut()
            .filter(|ge| !ge.mrule.is_unifac() && ge.mrule != crate::MixingRule::WongSandler)
            .ok_or(Error::Unsupported(
                "Huron-Vidal / NRTL parameters are not allocated for this mixing rule",
            ))?;
        let ij = &mut ge.pairs[i * nc + j];
        ij.alpha_ij = params.alpha_ij;
        ij.poly_ij = [params.a_ge_ij, params.b_ge_ij, params.c_ge_ij];
        ij.correlation = GeCorrelation::Polynomial;
        ij.has_record = true;

        let ji = &mut ge.pairs[j * nc + i];
        ji.alpha_ij = params.alpha_ji;
        ji.poly_ij = [params.a_ge_ji, params.b_ge_ji, params.c_ge_ji];
        ji.correlation = GeCorrelation::Polynomial;
        ji.has_record = true;
        Ok(())
    }

    /// Retrieve constant Wong-Sandler binary parameters for pair `(i, j)`
    /// (`thermopack_getWSparam`).
    pub fn ws_param(&self, i: usize, j: usize) -> Result<WsPairParams> {
        self.validate_pair_indices(i, j)?;
        if !self.mixing_rule.is_wong_sandler() {
            return Err(Error::Unsupported(
                "Wong-Sandler parameters are not allocated for this mixing rule",
            ));
        }
        let ws = self.wong_sandler.as_ref().ok_or(Error::Unsupported(
            "Wong-Sandler parameters are not initialized",
        ))?;
        let nc = self.components().len();
        let ij = i * nc + j;
        let ji = j * nc + i;
        Ok(WsPairParams {
            alpha_ij: ws.alphaij[ij],
            alpha_ji: ws.alphaij[ji],
            k_ij: ws.kij[ij].numerator[0],
            k_ji: ws.kij[ji].numerator[0],
            tau_ij: ws.tauij[ij].numerator[0],
            tau_ji: ws.tauij[ji].numerator[0],
        })
    }

    /// Tune constant Wong-Sandler binary parameters for pair `(i, j)`
    /// (`thermopack_setWSparam`).
    pub fn set_ws_param(&mut self, i: usize, j: usize, params: WsPairParams) -> Result<()> {
        self.validate_pair_indices(i, j)?;
        if !self.mixing_rule.is_wong_sandler() {
            return Err(Error::Unsupported(
                "Wong-Sandler parameters are not allocated for this mixing rule",
            ));
        }
        let vals = [
            params.alpha_ij,
            params.alpha_ji,
            params.k_ij,
            params.k_ji,
            params.tau_ij,
            params.tau_ji,
        ];
        if vals.iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput("finite WS parameters required"));
        }
        let nc = self.components().len();
        let ws = self
            .wong_sandler
            .get_or_insert_with(|| WongSandlerData::new(nc));
        ws.set_alphaij(i, j, params.alpha_ij)?;
        ws.set_alphaij(j, i, params.alpha_ji)?;
        ws.set_kij(i, j, Fraction::constant(params.k_ij))?;
        ws.set_kij(j, i, Fraction::constant(params.k_ji))?;
        ws.set_tauij(i, j, Fraction::constant(params.tau_ij))?;
        ws.set_tauij(j, i, Fraction::constant(params.tau_ji))?;
        let synced = ws.clone();
        if let Some(ge) = &mut self.excess_gibbs {
            ge.wong_sandler = Some(synced);
        }
        Ok(())
    }

    /// Retrieve the stored alpha parameters `[c1, c2, c3]` of component `i`
    /// (`thermopack_getTWUparam`).
    pub fn twu_param(&self, i: usize) -> Result<[f64; 3]> {
        Ok(self.alpha(i)?.parameters)
    }

    /// Set the Twu alpha correlation parameters `(c1, c2, c3)` for component `i`
    /// (`thermopack_setTWUparam`).
    pub fn set_twu_param(&mut self, i: usize, c1: f64, c2: f64, c3: f64) -> Result<()> {
        self.set_alpha(i, AlphaCorrelation::Twu, &[c1, c2, c3])
    }

    /// Retrieve the stored alpha parameters `[c1, c2, c3]` of component `i`
    /// (`thermopack_getMCparam`).
    pub fn mc_param(&self, i: usize) -> Result<[f64; 3]> {
        Ok(self.alpha(i)?.parameters)
    }

    /// Set the Mathias-Copeman alpha correlation parameters `(c1, c2, c3)` for
    /// component `i` (`thermopack_setMCparam`).
    pub fn set_mc_param(&mut self, i: usize, c1: f64, c2: f64, c3: f64) -> Result<()> {
        self.set_alpha(i, AlphaCorrelation::MathiasCopeman, &[c1, c2, c3])
    }

    /// Set the alpha correlation of component `i` by correlation name and parameter
    /// slice (`thermopack_set_alpha_corr`).
    pub fn set_alpha_corr(&mut self, i: usize, name: &str, params: &[f64]) -> Result<()> {
        let corr = AlphaCorrelation::from_name(name)?;
        self.set_alpha(i, corr, params)
    }

    /// Set the covolume beta correlation of component `i` by correlation name and
    /// parameter slice (`thermopack_set_beta_corr`).
    pub fn set_beta_corr(&mut self, i: usize, name: &str, params: &[f64]) -> Result<()> {
        let corr = BetaCorrelation::from_name(name)?;
        let beta_params = match corr {
            BetaCorrelation::Classic => {
                if !params.is_empty() && params.len() != 2 {
                    return Err(Error::InvalidInput(
                        "classic beta expects 0 or 2 parameters",
                    ));
                }
                [0.0, 0.0]
            }
            BetaCorrelation::Quantum => {
                if params.len() != 2 {
                    return Err(Error::InvalidInput(
                        "quantum beta requires exactly 2 parameters",
                    ));
                }
                [params[0], params[1]]
            }
        };
        let beta = Beta::new(corr, beta_params)?;
        self.set_beta(i, beta)
    }

    /// Retrieve volume translation parameters for component `i`
    /// (`thermopack_get_volume_shift_parameters`).
    pub fn volume_shift_parameters(&self, i: usize) -> Result<VolumeShiftParams> {
        Ok(VolumeShiftParams::from(self.volume_translation(i)?))
    }

    /// Set volume translation parameters for component `i` without altering the
    /// enabled flag (`thermopack_set_volume_shift_parameters`).
    pub fn set_volume_shift_parameters(
        &mut self,
        i: usize,
        params: VolumeShiftParams,
    ) -> Result<()> {
        let translation = VolumeTranslation::try_from(params)?;
        self.set_volume_translation(i, translation)
    }

    /// Retrieve UNIFAC main-group interaction parameters `(a_ij, b_ij, c_ij)` for
    /// an active main-group pair (`unifac::getUNIFACgroupInteraction`).
    pub fn unifac_group_interaction(
        &self,
        main_group_i: u8,
        main_group_j: u8,
    ) -> Result<(f64, f64, f64)> {
        let unifac = self
            .excess_gibbs
            .as_ref()
            .and_then(|ge| ge.unifac.as_ref())
            .ok_or(Error::Unsupported("model does not use UNIFAC mixing"))?;
        unifac
            .get_group_interaction(main_group_i, main_group_j)
            .ok_or(Error::InvalidInput(
                "UNIFAC main-group pair not active in mixture",
            ))
    }

    /// Set UNIFAC main-group interaction parameters `(a_ij, b_ij, c_ij)` for an
    /// active main-group pair (`unifac::setUNIFACgroupInteraction`).
    pub fn set_unifac_group_interaction(
        &mut self,
        main_group_i: u8,
        main_group_j: u8,
        aij: f64,
        bij: f64,
        cij: f64,
    ) -> Result<()> {
        let unifac = self
            .excess_gibbs
            .as_mut()
            .and_then(|ge| ge.unifac.as_mut())
            .ok_or(Error::Unsupported("model does not use UNIFAC mixing"))?;
        unifac.set_group_interaction(main_group_i, main_group_j, aij, bij, cij)
    }

    /// Redefine critical temperature, critical pressure, and acentric factor of
    /// component `i`, recomputing `a0[i]` and `bi[i]` (`cbselect::redefine_TcPcAcf_comp_cubic`).
    pub fn redefine_critical_parameters(
        &mut self,
        i: usize,
        tc: f64,
        pc: f64,
        acf: f64,
    ) -> Result<()> {
        positive(tc, "critical temperature")?;
        positive(pc, "critical pressure")?;
        if !acf.is_finite() {
            return Err(Error::InvalidInput("finite acentric factor required"));
        }
        let (omega_a, omega_b) = crate::cubic::omega_ab(self.kind());
        let rt_c = R * tc;
        let a0 = omega_a * rt_c * rt_c / pc;
        let bi = omega_b * R * tc / pc;
        self.set_component_criticals_internal(i, tc, pc, acf, a0, bi)
    }
}
