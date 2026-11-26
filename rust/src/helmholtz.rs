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
    pub(crate) fn zero() -> Self {
        Self {
            value: 0.0,
            dt: 0.0,
            drho: 0.0,
            dtt: 0.0,
            dtrho: 0.0,
            drhorho: 0.0,
        }
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
    pub(crate) fn mul(self, other: Self) -> Self {
        Self {
            value: self.value * other.value,
            dt: self.dt * other.value + self.value * other.dt,
            drho: self.drho * other.value + self.value * other.drho,
            dtt: self.dtt * other.value + 2.0 * self.dt * other.dt + self.value * other.dtt,
            dtrho: self.dtrho * other.value
                + self.dt * other.drho
                + self.drho * other.dt
                + self.value * other.dtrho,
            drhorho: self.drhorho * other.value
                + 2.0 * self.drho * other.drho
                + self.value * other.drhorho,
        }
    }
    pub(crate) fn scale_inv_t_power(self, temperature: f64, power: i32) -> Self {
        let inv_t = 1.0 / temperature;
        let b0 = inv_t.powi(power);
        let k = f64::from(power);
        let b1 = -k * b0 * inv_t;
        let b2 = k * (k + 1.0) * b0 * inv_t * inv_t;
        Self {
            value: b0 * self.value,
            dt: b1 * self.value + b0 * self.dt,
            drho: b0 * self.drho,
            dtt: b2 * self.value + 2.0 * b1 * self.dt + b0 * self.dtt,
            dtrho: b1 * self.drho + b0 * self.dtrho,
            drhorho: b0 * self.drhorho,
        }
    }
    pub(crate) fn molar_tvn(self, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        let inv_v = 1.0 / volume;
        let inv_v2 = inv_v * inv_v;
        let rho = moles * inv_v;
        let result = PureHelmholtz {
            value: self.value,
            dt: self.dt,
            dv: -rho * inv_v * self.drho,
            dn: inv_v * self.drho,
            dtt: self.dtt,
            dtv: -rho * inv_v * self.dtrho,
            dtn: inv_v * self.dtrho,
            dvv: (2.0 * rho * self.drho + rho * rho * self.drhorho) * inv_v2,
            dvn: -(self.drho + rho * self.drhorho) * inv_v2,
            dnn: self.drhorho * inv_v2,
        };
        if result.values().iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput(
                "non-finite molar Helmholtz derivatives",
            ));
        }
        Ok(result)
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
    pub(crate) fn pressure(
        self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<PurePressure> {
        positive(temperature, "temperature")?;
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        let rt = crate::R * temperature;
        let value = -rt * self.dv + moles * rt / volume;
        let dt = -rt * self.dtv + value / temperature;
        let dv = -rt * self.dvv - moles * rt / (volume * volume);
        let dn = -rt * self.dvn + rt / volume;
        if ![value, dt, dv, dn].iter().all(|v| v.is_finite()) {
            return Err(Error::InvalidInput("non-finite pressure derivatives"));
        }
        Ok(PurePressure { value, dt, dv, dn })
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
