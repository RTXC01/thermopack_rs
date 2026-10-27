//! Owned pure-fluid MBWR total caloric properties at explicit TV states.
//!
//! Ported from `src/mbwr_additional.f90` and `src/eostv.f90` in the
//! ThermoPack source distribution. Copyright thermotools and Contributors.
//! Apache-2.0.
//!
//! ```
//! use thermopack_rs::{MbwrFluid,MbwrKind};
//! let fluid=MbwrFluid::new(MbwrKind::Mbwr32,"C1")?;
//! let prepared=fluid.prepare(300.0)?;
//! let state=prepared.caloric_tv(1e-3,1.0)?;
//! assert!((state.enthalpy.value-state.internal_energy.value
//!          -state.pressure.value*1e-3).abs()<1e-8);
//! # Ok::<(),thermopack_rs::Error>(())
//! ```
use crate::{
    Error, IdealGas, IdealProperties, Mbwr, MbwrKind, MbwrRecord, PreparedMbwr, PureHelmholtz,
    PurePressure, R, Result, error::positive,
};

mod mbwr_total_tp;
pub use mbwr_total_tp::{
    MbwrSourceTpProperties, MbwrSourceTpState, MbwrTpProperties, PreparedMbwrFluidTp,
    PreparedMbwrSourceTp,
};

/// A total TV property and its partial derivatives at fixed remaining T, V,
/// and n. Energy, entropy, and chemical-potential units depend on the property.
#[derive(Debug, Clone, Copy)]
pub struct MbwrTvQuantity {
    pub value: f64,
    pub dt: f64,
    pub dv: f64,
    pub dn: f64,
}

/// Total pure-fluid properties at explicit T, total volume, and mole number.
/// U, H, S, A, and G are extensive; pressure and chemical potential are
/// intensive. Derivatives are with respect to K, m³, and mol.
#[derive(Debug, Clone, Copy)]
pub struct MbwrTvProperties {
    pub pressure: PurePressure,
    pub helmholtz_energy: MbwrTvQuantity,
    pub internal_energy: MbwrTvQuantity,
    pub enthalpy: MbwrTvQuantity,
    pub entropy: MbwrTvQuantity,
    pub gibbs_energy: MbwrTvQuantity,
    pub chemical_potential: MbwrTvQuantity,
    /// (∂U/∂T)V,n [J/K].
    pub heat_capacity_v: f64,
    /// Constant-pressure heat capacity when the local pressure response is
    /// nonsingular [J/K].
    pub heat_capacity_p: Option<f64>,
    /// Residual ln(phi), defined only for positive pressure.
    pub ln_fugacity_coefficient: Option<f64>,
    /// Speed of sound [m/s] when the local squared response is finite and positive.
    pub speed_of_sound: Option<f64>,
    /// Total reduced Helmholtz jet F=A/(RT), differentiated in T,V,n.
    pub reduced_helmholtz: PureHelmholtz,
}

/// Owned pure MBWR coefficients and ideal-gas Cp metadata. Database-backed
/// construction keeps ThermoPack's raw zero H/S integration constants.
#[derive(Debug, Clone)]
pub struct MbwrFluid {
    eos: Mbwr,
    ideal: IdealGas,
    molar_mass: f64,
}

impl MbwrFluid {
    /// Select the first source fit with `Mbwr::new` alias lookup. Ideal metadata
    /// is then matched by that fit's exact canonical database identifier.
    pub fn new(kind: MbwrKind, identifier: &str) -> Result<Self> {
        Self::from_eos(Mbwr::new(kind, identifier)?)
    }

    /// Select a source fit and require a matching canonical component record.
    pub fn from_record(source_name: &str) -> Result<Self> {
        Self::from_eos(Mbwr::from_record(source_name)?)
    }

    /// Combine a coefficient model with explicitly supplied ideal metadata and
    /// molar mass [kg/mol]. Useful for legacy MBWR identifiers absent from the
    /// general component database.
    pub fn with_ideal_gas(eos: Mbwr, ideal: IdealGas, molar_mass_kg_per_mol: f64) -> Result<Self> {
        positive(molar_mass_kg_per_mol, "molar mass")?;
        Ok(Self {
            eos,
            ideal,
            molar_mass: molar_mass_kg_per_mol,
        })
    }

    fn from_eos(eos: Mbwr) -> Result<Self> {
        let record: &'static MbwrRecord = eos.record().ok_or(Error::InvalidInput(
            "database-backed MBWR record required for ideal metadata",
        ))?;
        let component = crate::components()
            .iter()
            .find(|component| component.ident.eq_ignore_ascii_case(record.identifier))
            .ok_or_else(|| Error::UnknownComponent(record.identifier.to_owned()))?;
        let ideal = IdealGas::with_integration_constants(component, 0.0, 0.0)?;
        Self::with_ideal_gas(eos, ideal, component.molar_mass / 1000.0)
    }

    pub fn mbwr(&self) -> &Mbwr {
        &self.eos
    }

    pub fn ideal_gas(&self) -> &IdealGas {
        &self.ideal
    }

    pub fn ideal_gas_mut(&mut self) -> &mut IdealGas {
        &mut self.ideal
    }

    pub fn molar_mass(&self) -> f64 {
        self.molar_mass
    }

    /// Prepare coefficients and ideal caloric values once at temperature T.
    pub fn prepare(&self, temperature: f64) -> Result<PreparedMbwrFluid<'_>> {
        let eos = self.eos.prepare(temperature)?;
        let ideal = self.ideal.properties(temperature, 1.0)?;
        Ok(PreparedMbwrFluid {
            eos,
            ideal,
            temperature,
            source_cp: ideal.cp,
            molar_mass: self.molar_mass,
        })
    }

    /// Total Helmholtz energy and its TVn gradient [J].
    pub fn helmholtz(&self, t: f64, v: f64, n: f64) -> Result<MbwrTvQuantity> {
        self.prepare(t)?.helmholtz(v, n)
    }

    /// Total reduced Helmholtz jet F=A/(RT) and its TVn derivatives.
    pub fn reduced_helmholtz(&self, t: f64, v: f64, n: f64) -> Result<PureHelmholtz> {
        self.prepare(t)?.reduced_helmholtz(v, n)
    }

    /// Total caloric properties and TVn gradients at supplied T, total V, n.
    pub fn caloric_tv(&self, t: f64, v: f64, n: f64) -> Result<MbwrTvProperties> {
        self.prepare(t)?.caloric_tv(v, n)
    }

    /// Source-order residual integration and source Cp derivative conventions.
    pub fn baseline_helmholtz(&self, t: f64, v: f64, n: f64) -> Result<MbwrTvQuantity> {
        self.prepare(t)?.baseline_helmholtz(v, n)
    }

    /// Source-order total reduced Helmholtz jet, including source Cp
    /// second-derivative conventions.
    pub fn baseline_reduced_helmholtz(&self, t: f64, v: f64, n: f64) -> Result<PureHelmholtz> {
        self.prepare(t)?.baseline_reduced_helmholtz(v, n)
    }

    /// Source-order residual integration and source Cp derivative conventions.
    pub fn baseline_caloric_tv(&self, t: f64, v: f64, n: f64) -> Result<MbwrTvProperties> {
        self.prepare(t)?.baseline_caloric_tv(v, n)
    }
}

/// Prepared total MBWR TV evaluator. The ideal H, S, Cp, and source Cp values
/// are cached once; each volume/mole evaluation adds one residual snapshot.
#[derive(Debug, Clone)]
pub struct PreparedMbwrFluid<'a> {
    eos: PreparedMbwr<'a>,
    ideal: IdealProperties,
    temperature: f64,
    source_cp: f64,
    molar_mass: f64,
}

impl PreparedMbwrFluid<'_> {
    pub fn helmholtz(&self, v: f64, n: f64) -> Result<MbwrTvQuantity> {
        let (_, f) = self.jets(v, n, false)?;
        helmholtz_quantity(self.temperature, f)
    }

    pub fn reduced_helmholtz(&self, v: f64, n: f64) -> Result<PureHelmholtz> {
        self.jets(v, n, false).map(|(_, total)| total)
    }

    pub fn caloric_tv(&self, v: f64, n: f64) -> Result<MbwrTvProperties> {
        self.evaluate(v, n, false)
    }

    pub fn baseline_helmholtz(&self, v: f64, n: f64) -> Result<MbwrTvQuantity> {
        let (_, f) = self.jets(v, n, true)?;
        helmholtz_quantity(self.temperature, f)
    }

    pub fn baseline_reduced_helmholtz(&self, v: f64, n: f64) -> Result<PureHelmholtz> {
        self.jets(v, n, true).map(|(_, total)| total)
    }

    pub fn baseline_caloric_tv(&self, v: f64, n: f64) -> Result<MbwrTvProperties> {
        self.evaluate(v, n, true)
    }

    fn evaluate(&self, v: f64, n: f64, baseline: bool) -> Result<MbwrTvProperties> {
        self.evaluate_with_residual(v, n, baseline)
            .map(|(properties, _)| properties)
    }

    fn evaluate_with_residual(
        &self,
        v: f64,
        n: f64,
        baseline: bool,
    ) -> Result<(MbwrTvProperties, PureHelmholtz)> {
        let (residual, f) = self.jets(v, n, baseline)?;
        let t = self.temperature;
        let rt = R * t;
        let pressure = PurePressure {
            value: -rt * f.dv,
            dt: -R * (f.dv + t * f.dtv),
            dv: -rt * f.dvv,
            dn: -rt * f.dvn,
        };
        let internal_energy = quantity(
            -rt * t * f.dt,
            -R * t * (2.0 * f.dt + t * f.dtt),
            -rt * t * f.dtv,
            -rt * t * f.dtn,
        );
        let entropy = quantity(
            -R * (f.value + t * f.dt),
            -R * (2.0 * f.dt + t * f.dtt),
            -R * (f.dv + t * f.dtv),
            -R * (f.dn + t * f.dtn),
        );
        let helmholtz_energy =
            quantity(rt * f.value, R * (f.value + t * f.dt), rt * f.dv, rt * f.dn);
        let pv = pressure.value * v;
        let enthalpy = quantity(
            internal_energy.value + pv,
            internal_energy.dt + v * pressure.dt,
            internal_energy.dv + pressure.value + v * pressure.dv,
            internal_energy.dn + v * pressure.dn,
        );
        let gibbs_energy = quantity(
            helmholtz_energy.value + pv,
            helmholtz_energy.dt + v * pressure.dt,
            helmholtz_energy.dv + pressure.value + v * pressure.dv,
            helmholtz_energy.dn + v * pressure.dn,
        );
        let chemical_potential =
            quantity(rt * f.dn, R * (f.dn + t * f.dtn), rt * f.dvn, rt * f.dnn);
        let heat_capacity_v = internal_energy.dt;
        let ideal_pv_scale = rt * n / (v * v);
        let heat_capacity_p = if pressure.dv.abs() > 64.0 * f64::EPSILON * ideal_pv_scale {
            let cp = heat_capacity_v - t * pressure.dt * pressure.dt / pressure.dv;
            cp.is_finite().then_some(cp)
        } else {
            None
        };
        let z = pressure.value * v / (n * rt);
        let ln_fugacity_coefficient = if pressure.value > 0.0 && z.is_finite() && z > 0.0 {
            // Preserve the dilute residual when rounding the total pressure
            // would make Z equal one. Near Z=0 use the positive pressure ratio
            // instead of subtracting nearly equal terms to recover Z.
            let delta_z = -v * residual.dv / n;
            let log_z = if !baseline && delta_z.abs() < 0.5 {
                delta_z.ln_1p()
            } else {
                z.ln()
            };
            let value = residual.dn - log_z;
            value.is_finite().then_some(value)
        } else {
            None
        };
        let speed_squared = -v * v / (n * self.molar_mass)
            * (pressure.dv - t * pressure.dt * pressure.dt / heat_capacity_v);
        let speed_of_sound = if speed_squared.is_finite() && speed_squared > 0.0 {
            Some(speed_squared.sqrt())
        } else {
            None
        };
        let properties = MbwrTvProperties {
            pressure,
            helmholtz_energy,
            internal_energy,
            enthalpy,
            entropy,
            gibbs_energy,
            chemical_potential,
            heat_capacity_v,
            heat_capacity_p,
            ln_fugacity_coefficient,
            speed_of_sound,
            reduced_helmholtz: f,
        };
        let values = [
            pressure.value,
            pressure.dt,
            pressure.dv,
            pressure.dn,
            helmholtz_energy.value,
            helmholtz_energy.dt,
            helmholtz_energy.dv,
            helmholtz_energy.dn,
            internal_energy.value,
            internal_energy.dt,
            internal_energy.dv,
            internal_energy.dn,
            enthalpy.value,
            enthalpy.dt,
            enthalpy.dv,
            enthalpy.dn,
            entropy.value,
            entropy.dt,
            entropy.dv,
            entropy.dn,
            gibbs_energy.value,
            gibbs_energy.dt,
            gibbs_energy.dv,
            gibbs_energy.dn,
            chemical_potential.value,
            chemical_potential.dt,
            chemical_potential.dv,
            chemical_potential.dn,
            heat_capacity_v,
        ];
        if values.iter().any(|value| !value.is_finite()) {
            return Err(Error::InvalidInput("non-finite MBWR TV caloric state"));
        }
        Ok((properties, residual))
    }

    fn jets(&self, v: f64, n: f64, baseline: bool) -> Result<(PureHelmholtz, PureHelmholtz)> {
        positive(v, "total volume")?;
        positive(n, "mole number")?;
        let residual = if baseline {
            self.eos.baseline_residual_helmholtz(v, n)?
        } else {
            self.eos.residual_helmholtz(v, n)?
        };
        let ideal = ideal_helmholtz(
            self.temperature,
            v,
            n,
            self.ideal,
            if baseline {
                self.source_cp
            } else {
                self.ideal.enthalpy_dt
            },
        )?;
        Ok((residual, add_helmholtz(residual, ideal)?))
    }
}

fn helmholtz_quantity(t: f64, f: PureHelmholtz) -> Result<MbwrTvQuantity> {
    let result = quantity(
        R * t * f.value,
        R * (f.value + t * f.dt),
        R * t * f.dv,
        R * t * f.dn,
    );
    if [result.value, result.dt, result.dv, result.dn]
        .iter()
        .any(|value| !value.is_finite())
    {
        return Err(Error::InvalidInput("non-finite MBWR TV Helmholtz energy"));
    }
    Ok(result)
}

fn ideal_helmholtz(
    t: f64,
    v: f64,
    n: f64,
    ideal: IdealProperties,
    cp_for_second_derivative: f64,
) -> Result<PureHelmholtz> {
    positive(t, "temperature")?;
    positive(v, "total volume")?;
    positive(n, "mole number")?;
    let rt = R * t;
    let log_concentration = (n * rt / v).ln();
    let f = n * (ideal.enthalpy / rt - ideal.entropy / R + log_concentration - 1.0);
    let ft = n * (-ideal.enthalpy / (R * t * t) + 1.0 / t);
    let ftt = n
        * (2.0 * ideal.enthalpy / (R * t * t * t)
            - cp_for_second_derivative / (R * t * t)
            - 1.0 / (t * t));
    let result = PureHelmholtz {
        value: f,
        dt: ft,
        dv: -n / v,
        dn: ideal.enthalpy / rt - ideal.entropy / R + log_concentration,
        dtt: ftt,
        dtv: 0.0,
        dtn: -ideal.enthalpy / (R * t * t) + 1.0 / t,
        dvv: n / (v * v),
        dvn: -1.0 / v,
        dnn: 1.0 / n,
    };
    if result.values().iter().any(|value| !value.is_finite()) {
        return Err(Error::InvalidInput("non-finite ideal TV Helmholtz state"));
    }
    Ok(result)
}

fn add_helmholtz(a: PureHelmholtz, b: PureHelmholtz) -> Result<PureHelmholtz> {
    let result = PureHelmholtz {
        value: a.value + b.value,
        dt: a.dt + b.dt,
        dv: a.dv + b.dv,
        dn: a.dn + b.dn,
        dtt: a.dtt + b.dtt,
        dtv: a.dtv + b.dtv,
        dtn: a.dtn + b.dtn,
        dvv: a.dvv + b.dvv,
        dvn: a.dvn + b.dvn,
        dnn: a.dnn + b.dnn,
    };
    if result.values().iter().any(|value| !value.is_finite()) {
        return Err(Error::InvalidInput("non-finite total TV Helmholtz state"));
    }
    Ok(result)
}

fn quantity(value: f64, dt: f64, dv: f64, dn: f64) -> MbwrTvQuantity {
    MbwrTvQuantity { value, dt, dv, dn }
}
