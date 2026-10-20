//! Native thermodynamic models ported from the checked-in ThermoPack fork.
//!
//! All quantities use SI molar units. See `docs/rust/coverage.md` for the
//! explicit migration status; this crate is not yet a complete replacement.
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
pub mod caloric;
pub mod component;
pub mod cpa;
pub mod cpa_data;
pub mod cubic;
mod error;
pub mod excess_gibbs;
pub mod excess_gibbs_data;
pub mod flash;
pub mod fmt;
pub mod helmholtz;
mod hydrogen;
pub mod hyperdual;
pub mod ideal;
pub mod lee_kesler;
pub mod mbwr;
pub mod multiparameter;
pub mod multiparameter_data;
pub mod mut_flash;
pub mod numerics;
pub mod pc_saft;
pub mod pc_saft_data;
pub mod pets;
pub mod saft_vrmie;
pub mod saft_vrmie_data;
pub mod solid;
pub mod solid_co2;
pub mod solid_data;
pub mod solid_saturation;
pub mod solid_water;
pub mod translation;
pub mod unifac;
pub mod unifac_data;

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
