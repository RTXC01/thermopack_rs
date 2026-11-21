//! Reduced residual Helmholtz derivatives for pure-fluid models.
use crate::{Error, Result, error::positive};

/// Dimensionless residual Helmholtz energy per mole, alpha=Ares/(n R T),
/// differentiated in temperature (K) and molar density (mol/m3).
#[derive(Debug, Clone, Copy)]
pub struct DensityHelmholtz {
    pub value: f64,
    pub dt: f64,
    pub drho: f64,
    pub dtt: f64,
    pub dtrho: f64,
    pub drhorho: f64,
}

impl DensityHelmholtz {
    pub(crate) fn values(self) -> [f64; 6] {
        [
            self.value,
            self.dt,
            self.drho,
            self.dtt,
            self.dtrho,
            self.drhorho,
        ]
    }
    pub(crate) fn add(self, other: Self) -> Self {
        Self {
            value: self.value + other.value,
            dt: self.dt + other.dt,
            drho: self.drho + other.drho,
            dtt: self.dtt + other.dtt,
            dtrho: self.dtrho + other.dtrho,
            drhorho: self.drhorho + other.drhorho,
        }
    }
    pub(crate) fn extensive(self, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        let rho = moles / volume;
        let result = PureHelmholtz {
            value: moles * self.value,
            dt: moles * self.dt,
            dv: -rho * rho * self.drho,
            dn: self.value + rho * self.drho,
            dtt: moles * self.dtt,
            dtv: -rho * rho * self.dtrho,
            dtn: self.dt + rho * self.dtrho,
            dvv: (2. * rho * rho * self.drho + rho * rho * rho * self.drhorho) / volume,
            dvn: -(2. * rho * self.drho + rho * rho * self.drhorho) / volume,
            dnn: (2. * self.drho + rho * self.drhorho) / volume,
        };
        if result.values().iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput(
                "non-finite extensive Helmholtz derivatives",
            ));
        }
        Ok(result)
    }
}

/// Reduced extensive Helmholtz energy F=A/(R T), in mol. The returning
/// method specifies the total, ideal, or residual contribution.
/// Derivatives hold the remaining T (K), total V (m3), n (mol) fixed.
#[derive(Debug, Clone, Copy)]
pub struct PureHelmholtz {
    pub value: f64,
    pub dt: f64,
    pub dv: f64,
    pub dn: f64,
    pub dtt: f64,
    pub dtv: f64,
    pub dtn: f64,
    pub dvv: f64,
    pub dvn: f64,
    pub dnn: f64,
}
impl PureHelmholtz {
    pub(crate) fn values(self) -> [f64; 10] {
        [
            self.value, self.dt, self.dv, self.dn, self.dtt, self.dtv, self.dtn, self.dvv,
            self.dvn, self.dnn,
        ]
    }
}

/// Pure-fluid pressure and its T, total-V and mole-number derivatives.
#[derive(Debug, Clone, Copy)]
pub struct PurePressure {
    pub value: f64,
    pub dt: f64,
    pub dv: f64,
    pub dn: f64,
}
