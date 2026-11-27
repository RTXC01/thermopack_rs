//! Complex cubic model initializers ported from `src/complexmodelinit.f90`.
//!
//! Provides `init_vtpr` (VTPR: PR + VTPR UNIFAC + Twu alpha + `b_exponent = 4/3`
//! plus Peneloux PR volume translation) and `init_umr` (UMR-PRU: PR + UMR UNIFAC
//! plus UMR alpha + `b_exponent = 2.0` + Peneloux PR volume translation).

use crate::cubic::{Cubic, CubicKind, CubicOptions};
use crate::error::Result;
use crate::excess_gibbs::MixingRule;
use crate::{Alpha, AlphaCorrelation};

/// Initialize the VTPR (Volume-Translated Peng-Robinson) group-contribution
/// model (`complexmodelinit::init_VTPR`) using `"DEFAULT"` parameter references.
///
/// Configures Peng-Robinson with [`MixingRule::Vtpr`], [`AlphaCorrelation::Twu`],
/// covolume power-mean exponent `b_exponent = 4.0 / 3.0`, and default Peneloux
/// PR volume translation.
pub fn init_vtpr(names: &[&str]) -> Result<Cubic> {
    init_vtpr_with_refs(names, "DEFAULT", "DEFAULT")
}

/// Initialize the VTPR model with explicit `kij_ref` and `alpha_ref` parameter
/// references (`complexmodelinit::init_VTPR`).
pub fn init_vtpr_with_refs(names: &[&str], _kij_ref: &str, alpha_ref: &str) -> Result<Cubic> {
    let mut model = Cubic::with_options(
        CubicKind::PengRobinson,
        names,
        CubicOptions {
            alpha: Some(AlphaCorrelation::Twu),
            parameter_reference: alpha_ref,
            volume_translation: false,
            mixing_rule: MixingRule::Vtpr,
            wong_sandler: None,
        },
    )?;
    model.set_b_exponent(4.0 / 3.0)?;
    model.enable_volume_translation("DEFAULT")?;
    Ok(model)
}

/// Initialize the UMR-PRU (Universal Mixing Rule Peng-Robinson UNIFAC) model
/// (`complexmodelinit::init_UMR`) using `"DEFAULT"` parameter references.
///
/// Configures Peng-Robinson with [`MixingRule::Umr`], covolume power-mean
/// exponent `b_exponent = 2.0`, default Peneloux PR volume translation, and
/// UMR alpha correlations:
/// - `H2O` uses Mathias-Copeman from the database (matching `cbalpha.f90`),
/// - non-water components use the UMR acentric polynomial
///   `m(w) = 0.384401 + 1.52276 w - 0.213808 w^2 + 0.034616 w^3 - 0.001976 w^4`
///   evaluated with the classic Soave functional form (`AlphaCorrelation::ClassicFit`),
///   avoiding the upstream `'ALPHA_UMR'` identifier mismatch and missing
///   `cbAlphaUMRIdx` evaluator branch in `src/complexmodelinit.f90` and `src/cbalpha.f90`.
pub fn init_umr(names: &[&str]) -> Result<Cubic> {
    init_umr_with_refs(names, "DEFAULT", "DEFAULT")
}

/// Initialize the UMR-PRU model with explicit `kij_ref` and `alpha_ref`
/// parameter references (`complexmodelinit::init_UMR`).
pub fn init_umr_with_refs(names: &[&str], _kij_ref: &str, alpha_ref: &str) -> Result<Cubic> {
    let mut model = Cubic::with_options(
        CubicKind::PengRobinson,
        names,
        CubicOptions {
            alpha: Some(AlphaCorrelation::Classic),
            parameter_reference: alpha_ref,
            volume_translation: false,
            mixing_rule: MixingRule::Umr,
            wong_sandler: None,
        },
    )?;
    let nc = model.components().len();
    for i in 0..nc {
        let comp = model.components()[i];
        let alpha = if comp.ident.eq_ignore_ascii_case("H2O") {
            Alpha::initialize(
                CubicKind::PengRobinson,
                &comp,
                AlphaCorrelation::Umr,
                alpha_ref,
            )?
        } else {
            let m = AlphaCorrelation::Umr.acentric_parameter(comp.acentric_factor)?;
            Alpha::new(AlphaCorrelation::ClassicFit, &[m])?
        };
        model.set_alpha_direct(i, alpha)?;
    }
    model.set_b_exponent(2.0)?;
    model.enable_volume_translation("DEFAULT")?;
    Ok(model)
}

impl Cubic {
    /// Initialize a VTPR (Volume-Translated Peng-Robinson) model with default
    /// references (`complexmodelinit::init_VTPR`).
    pub fn init_vtpr(names: &[&str]) -> Result<Self> {
        init_vtpr(names)
    }

    /// Initialize a VTPR model with explicit `kij_ref` and `alpha_ref`
    /// (`complexmodelinit::init_VTPR`).
    pub fn init_vtpr_with_refs(names: &[&str], kij_ref: &str, alpha_ref: &str) -> Result<Self> {
        init_vtpr_with_refs(names, kij_ref, alpha_ref)
    }

    /// Initialize a UMR-PRU model with default references
    /// (`complexmodelinit::init_UMR`).
    pub fn init_umr(names: &[&str]) -> Result<Self> {
        init_umr(names)
    }

    /// Initialize a UMR-PRU model with explicit `kij_ref` and `alpha_ref`
    /// (`complexmodelinit::init_UMR`).
    pub fn init_umr_with_refs(names: &[&str], kij_ref: &str, alpha_ref: &str) -> Result<Self> {
        init_umr_with_refs(names, kij_ref, alpha_ref)
    }
}
