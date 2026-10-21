//! Native thermodynamic models ported from the checked-in ThermoPack fork.
//!
//! All quantities use SI molar units. See `docs/rust/coverage.md` for the
//! explicit migration status; this crate is not yet a complete replacement.
//!
//! # The scope
//!
//! The crate holds the equations of state, the mixing rules, the phase
//! equilibrium and the property kernels of the fork:
//!
//! * the cubic equations of state with the van der Waals, the Huron-Vidal,
//!   the NRTL, the UNIFAC and the Wong-Sandler mixing rules, the pseudo
//!   components, the volume translation and the temperature dependent beta
//!   coefficient, together with the flash, the stability analysis, the
//!   saturation states, the isotherms, the isobars, the isentropes and the
//!   isenthalpies, the Joule-Thomson inversion and the spinodal,
//! * the cubic association model and the cubic plus association model,
//! * the PC-SAFT and the SAFT-VR Mie models with their parameter databases,
//!   the radial distribution function of the cubic association model, the
//!   hard-sphere model of Boublik with its packing fractions and its pair
//!   correlation function, and the quantum corrections of the hydrogen, the
//!   helium and the neon,
//! * the pure fluid equations of state of the Lee-Kesler tables, the
//!   multiparameter models, the MBWR model of the IAPWS, the PeTS model of
//!   the liquid hydrogen and the solid phases of the carbon dioxide and of the
//!   water,
//! * the utilities of the original: the Gauss and the Gauss-Kronrod rules,
//!   the hyperdual numbers, the chain rule and the ideal gas properties.
//!
//! # The parity
//!
//! Every kernel is ported from the original and qualified against the exact
//! source library. The reference probe `tools/rust_port/reference_probe.f90`
//! calls the original routines and reports their results, the validators of
//! `tools/rust_port` compare them with the native model and store the
//! qualified values in `tests/fixtures`. The Rust tests read those fixtures,
//! so the test suite fails without the reference library. The migration
//! status of every kernel, with the number of qualified cases and the largest
//! observed difference, is in `docs/rust/coverage.md`.
//!
//! # The performance
//!
//! The models prepare their state once and evaluate the properties of the
//! prepared state without an allocation. The preparation of a state is
//! dominated by the alpha correlation of the components and by the mixing
//! rule, and the property evaluation of the prepared state costs about
//! 20 ns for the pressure. The measured comparison with the source, the
//! implemented optimizations and the open targets are in
//! `docs/rust/performance.md`.
//!
//! ```
//! use thermopack_rs::{Cubic, CubicKind, Phase};
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let model = Cubic::new(CubicKind::PengRobinson, &["CO2", "C1"])?;
//! let state = model.prepare(300.0, &[0.8, 0.2])?;
//! let volume = state.volume(5e6, Phase::Vapor)?;
//! let mut ln_phi = [0.0; 2];
//! state.ln_fugacity_into(volume, &mut ln_phi)?;
//! assert!(volume > state.covolume());
//! # Ok(())
//! # }
//! ```
pub mod alpha;
pub mod beta;
pub mod bmcsl;
pub mod caloric;
pub mod chain_rule;
pub mod component;
pub mod cpa;
pub mod cpa_data;
pub mod critical;
pub mod cubic;
mod error;
pub mod excess_gibbs;
pub mod excess_gibbs_data;
pub mod flash;
pub mod fmt;
pub mod hardsphere_wca;
pub mod helmholtz;
mod hydrogen;
pub mod hyperdual;
pub mod ideal;
pub mod isolines;
pub mod joule_thompson;
pub mod lee_kesler;
pub mod ljs;
pub mod ljs_wca;
pub mod mbwr;
pub mod mixture_critical;
pub mod mixture_saturation;
pub mod mixture_spinodal;
pub mod multiparameter;
pub mod multiparameter_data;
pub mod mut_flash;
pub mod numerics;
pub mod pc_saft;
pub mod pc_saft_data;
pub mod pets;
pub mod quadrature;
pub mod quadrature_data;
pub mod rdf;
pub mod saft_vrmie;
pub mod saft_vrmie_data;
pub mod solid;
pub mod solid_co2;
pub mod solid_data;
pub mod solid_saturation;
pub mod solid_water;
pub mod spinodal;
pub mod spinodal_properties;
pub mod translation;
pub mod unifac;
pub mod unifac_data;
pub mod wong_sandler;

pub use alpha::{Alpha, AlphaCorrelation, AlphaValue};
pub use beta::{Beta, BetaCorrelation, BetaValue};
pub use caloric::{
    CaloricFlashOptions, EquilibriumProperties, EquilibriumResponse, VolumeFlashOptions,
};
pub use component::{Component, component, components};
pub use cpa::{
    Cpa, CpaComponent, CpaCubicKind, CpaPressureJet, CpaSaturation, CpaState, PreparedCpa,
};
pub use cpa_data::{AssocScheme, CpaBinaryRecord, CpaPureRecord};
pub use cubic::{
    Cubic, CubicKind, CubicOptions, Phase, PreparedState, Pressure, PseudoCritical, Residual,
    Saturation, TvProperties,
};
pub use error::{Error, Result};
pub use excess_gibbs::{ExcessGibbsMixture, MixingRule};
pub use excess_gibbs_data::{ExcessGibbsRecord, GeCorrelation};
pub use flash::{
    FlashOptions, FlashPhase, FlashResult, RachfordRiceIteration, Stability, StabilityTrial,
};
pub use helmholtz::{DensityHelmholtz, PureHelmholtz, PurePressure};
pub use ideal::{EntropyReference, IdealGas, IdealProperties};
pub use lee_kesler::{BwrCoeffs, BwrCoeffsDiff1, BwrCoeffsDiff2, LeeKesler, PreparedLeeKesler};
pub use mbwr::{
    Mbwr, MbwrDensityIteration, MbwrDensityParameters, MbwrDensityPath, MbwrEquilibrium,
    MbwrHelmholtzTerm, MbwrInitialization, MbwrKind, MbwrPressure, MbwrRecord, PreparedMbwr,
    PreparedMbwrDensity, mbwr_records,
};
pub use mbwr::{MbwrResidualTp, MbwrTpDerivatives, MbwrTpQuantity};
pub use mixture_saturation::{MixtureSaturation, SaturationKind};
pub use multiparameter::{
    HelmholtzJet, LennardJonesKind, LennardJonesParameters, Multiparameter,
    MultiparameterDensityIteration, MultiparameterState, PreparedMultiparameter,
};
pub use multiparameter_data::{
    AncillarySatDensity, IdealGasHelmholtz, MultiparameterRecord, lookup_multiparameter_record,
};
pub use pc_saft::{
    PcSaft, PcSaftComponent, PcSaftKind, PcSaftPressureJet, PcSaftSaturation, PcSaftState,
    PreparedPcSaft,
};
pub use pc_saft_data::{PcSaftKijRecord, PcSaftRecord};
pub use pets::{
    Pets, PetsContributions, PetsDiameter, PetsEquilibrium, PetsEquilibriumTv, PetsTvProperties,
    PetsVolumeIteration, PreparedPets, PureCriticalPoint,
};
pub use saft_vrmie::{
    PreparedSaftVrMie, SaftVrMie, SaftVrMieComponent, SaftVrMieHelmholtzTerms,
    SaftVrMiePressureJet, SaftVrMieSaturation, SaftVrMieState,
};
pub use saft_vrmie_data::{SaftVrMieKijRecord, SaftVrMieLijRecord, SaftVrMieRecord};
pub use solid::SolidLine;
pub use solid_co2::SolidCo2;
pub use solid_data::{
    SolidCorrelation, lookup_melting_correlation, lookup_sublimation_correlation,
};
pub use solid_water::{PreparedSolidWater, SolidWater, SolidWaterGibbs, SolidWaterState};
pub use translation::{TranslationKind, TranslationValue, VolumeTranslation};
pub use unifac::{Unifac, UnifacExcessGibbs, UnifacVariant};
pub use unifac_data::{
    UnifacComponent, UnifacInteraction, UnifacSubGroup, lookup_unifac_component,
    lookup_unifac_interaction,
};

/// Gas constant [J/(mol K)], using the exact expression from ThermoPack.
pub const R: f64 = 6.02214076e23 * 1.380649e-23;

/// The Avogadro constant in inverse mole.
pub const N_A: f64 = 6.02214076e23;
