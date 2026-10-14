//! Native thermodynamic models ported from the checked-in ThermoPack fork.
//!
//! All quantities use SI molar units. See `docs/rust/coverage.md` for the
//! explicit migration status; this crate is not yet a complete replacement.
pub mod component;
pub mod cubic;
mod error;
pub mod ideal;
pub mod numerics;

pub use component::{Component, component, components};
pub use cubic::{Cubic, CubicKind, Phase, PreparedState, Pressure, Residual};
pub use error::{Error, Result};

/// Gas constant [J/(mol K)], using the exact expression from ThermoPack.
pub const R: f64 = 6.02214076e23 * 1.380649e-23;
