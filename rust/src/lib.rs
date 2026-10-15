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
pub mod component;
pub mod cubic;
mod error;
pub mod flash;
pub mod ideal;
pub mod numerics;

pub use component::{Component, component, components};
pub use cubic::{Cubic, CubicKind, Phase, PreparedState, Pressure, Residual, Saturation};
pub use error::{Error, Result};
pub use flash::{FlashOptions, FlashPhase, FlashResult, Stability, StabilityTrial};

/// Gas constant [J/(mol K)], using the exact expression from ThermoPack.
pub const R: f64 = 6.02214076e23 * 1.380649e-23;
