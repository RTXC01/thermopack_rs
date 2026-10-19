//! Original MBWR critical table/grid initialization, separate from solving
//! analytic critical inflection equations.
use super::*;
use crate::PureCriticalPoint;

const KR: f64 = 1000. * R;
#[derive(Debug, Clone, Copy)]
pub struct MbwrInitialization {
    pub parameters: MbwrDensityParameters,
    pub compressibility: f64,
    pub from_critical_table: bool,
    /// Source SRK initial-guess covolume in m3/mol.
    pub srk_covolume: f64,
    pub srk_alpha_coefficient: f64,
    /// Source SRK initial-guess attraction in Pa*m6/mol2.
    pub srk_attraction: f64,
}

impl Mbwr {
    /// Lazily evaluate the original default metadata path once. Both success
    /// and failure are cached. Coefficient-only construction remains available
    /// even for legacy identifiers missing general component metadata.
    pub fn source_initialization(&self) -> Result<MbwrInitialization> {
        self.source_initialization
            .get_or_init(|| self.calculate_source_initialization())
            .clone()
    }
    fn calculate_source_initialization(&self) -> Result<MbwrInitialization> {
        let record = self.record.ok_or(Error::InvalidInput(
            "custom MBWR initialization requires explicit metadata",
        ))?;
        let c = data::COMPONENT_INITIALIZATION
            .iter()
            .find(|(id, _)| *id == record.identifier)
            .ok_or_else(|| Error::UnknownComponent(record.identifier.to_owned()))?
            .1;
        let table = data::STORED_PROPERTIES.iter().find(|(kind, id, _)| {
            *kind == self.kind.coefficient_count() as u8 && *id == record.identifier
        });
        let (t, p, rho) = if let Some((_, _, a)) = table {
            (a[0], a[1], a[2])
        } else {
            // Preserve original backend-unit denominator/order; no SI round-trip.
            let rho = c[1] / (KR * c[2] * c[0]);
            self.critical_grid_backend(c[0], c[1], rho, c[3])?
        };
        let parameters = MbwrDensityParameters::new(t, p, rho * 1000., c[3])?;
        let result = MbwrInitialization {
            parameters,
            compressibility: p / (rho * KR * t),
            from_critical_table: table.is_some(),
            srk_covolume: 0.08664 * KR * t / p * 1e-3,
            srk_alpha_coefficient: 0.480 + 1.547 * c[3] - 0.176 * c[3].powi(2),
            srk_attraction: 0.42747 * (KR * t).powi(2) / p * 1e-6,
        };
        if [
            result.compressibility,
            result.srk_covolume,
            result.srk_alpha_coefficient,
            result.srk_attraction,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(Error::InvalidInput("non-finite source MBWR initialization"));
        }
        Ok(result)
    }
    /// Original 0.001 K / 1 mol/m3 grid search in its narrow Tc +/- 0.7 K,
    /// rhoc +/- 1500 mol/m3 window. Returns its coarse estimate, not an
    /// analytically refined critical point. Undefined source bounds are errors.
    pub fn baseline_critical_grid(
        &self,
        initial: MbwrDensityParameters,
    ) -> Result<PureCriticalPoint> {
        let (t, p, rho) = self.critical_grid_backend(
            initial.critical_temperature,
            initial.critical_pressure,
            initial.critical_density * 1e-3,
            initial.acentric_factor,
        )?;
        Ok(PureCriticalPoint {
            temperature: t,
            pressure: p,
            volume: 1e-3 / rho,
            compressibility: p / (rho * KR * t),
        })
    }
    fn critical_grid_backend(
        &self,
        tc: f64,
        pc: f64,
        rc: f64,
        acf: f64,
    ) -> Result<(f64, f64, f64)> {
        MbwrDensityParameters::new(tc, pc, rc * 1000., acf)?;
        let (tl, tu, rl, ru) = (tc - 0.7, tc + 0.7, rc - 1.5, rc + 1.5);
        if ![tl, tu, rl, ru].iter().all(|v| v.is_finite())
            || tl + 0.001 <= tl
            || tu - 0.001 >= tu
            || rl + 0.001 <= rl
            || ru - 0.001 >= ru
        {
            return Err(Error::InvalidInput(
                "unresolved MBWR critical grid increments",
            ));
        }
        let mut lower = None;
        let mut t = tl;
        loop {
            t += 0.001;
            let (coefficients, polynomial, _) = self.source_pressure_coefficients(t)?;
            let mut previous = 0.;
            let mut rho = rl;
            let mut decline = false;
            loop {
                rho += 0.001;
                let p = PreparedMbwrDensity::pressure_value_backend(&coefficients, polynomial, rho);
                if !p.is_finite() {
                    return Err(Error::NotConverged(
                        "non-finite source MBWR critical grid pressure",
                    ));
                }
                if p < previous {
                    lower = Some((t, p, rho));
                    decline = true;
                    break;
                }
                previous = p;
                if rho > ru {
                    break;
                }
            }
            if !decline {
                break;
            }
            if t > tu {
                return Err(Error::NotConverged(
                    "MBWR critical grid exceeds initial temperature window",
                ));
            }
        }
        let lower = lower.ok_or(Error::NotConverged(
            "source MBWR critical grid lower bound is uninitialized",
        ))?;
        t = tu;
        let upper = loop {
            t -= 0.001;
            let (coefficients, polynomial, _) = self.source_pressure_coefficients(t)?;
            let mut previous = 1e12;
            let mut rho = ru;
            let mut ascent = None;
            loop {
                rho -= 0.001;
                let p = PreparedMbwrDensity::pressure_value_backend(&coefficients, polynomial, rho);
                if !p.is_finite() {
                    return Err(Error::NotConverged(
                        "non-finite source MBWR critical grid pressure",
                    ));
                }
                if p > previous {
                    ascent = Some((t, p, rho));
                    break;
                }
                previous = p;
                if rho < rl {
                    break;
                }
            }
            if let Some(a) = ascent {
                break a;
            }
            if t < tl {
                return Err(Error::NotConverged(
                    "MBWR critical grid exceeds initial temperature window",
                ));
            }
        };
        let result = (
            (lower.0 + upper.0) / 2.,
            (lower.1 + upper.1) / 2.,
            (lower.2 + upper.2) / 2.,
        );
        if ![result.0, result.1, result.2]
            .iter()
            .all(|v| v.is_finite() && *v > 0.)
        {
            return Err(Error::NotConverged(
                "nonphysical source MBWR critical grid estimate",
            ));
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn default_initialization_is_lazy_and_cache_survives_cloning() {
        let model = Mbwr::new(MbwrKind::Mbwr32, "C1").unwrap();
        assert!(model.source_initialization.get().is_none());
        let a = model.source_initialization().unwrap();
        assert!(a.from_critical_table);
        assert!(model.source_initialization.get().is_some());
        let cloned = model.clone();
        assert!(cloned.source_initialization.get().is_some());
        let b = cloned.source_initialization().unwrap();
        assert_eq!(
            a.parameters.critical_temperature,
            b.parameters.critical_temperature
        );
        let missing = Mbwr::from_record("mbwrcx12").unwrap();
        assert!(missing.source_initialization.get().is_none());
        let err = missing.source_initialization().unwrap_err();
        assert!(missing.source_initialization.get().is_some());
        assert_eq!(missing.source_initialization().unwrap_err(), err);
        // The undefined default path must not prevent coefficient-only use.
        assert!(
            missing
                .pressure_density(300., 100.)
                .unwrap()
                .value
                .is_finite()
        );
    }
    #[test]
    fn grid_ignores_unused_pc_and_acf_and_rejects_undefined_bounds() {
        let model = Mbwr::new(MbwrKind::Mbwr32, "C1").unwrap();
        let data = model.source_density_parameters().unwrap();
        let a = model.baseline_critical_grid(data).unwrap();
        let b = model
            .baseline_critical_grid(MbwrDensityParameters {
                critical_pressure: 1e-300,
                acentric_factor: 1e200,
                ..data
            })
            .unwrap();
        assert_eq!(
            (a.temperature, a.pressure, a.volume),
            (b.temperature, b.pressure, b.volume)
        );
        let ideal = Mbwr::with_parameters(MbwrKind::Bender19, &[0.; 19], 0.1).unwrap();
        let err = ideal
            .baseline_critical_grid(MbwrDensityParameters::new(300., 1e6, 10000., 0.).unwrap())
            .unwrap_err();
        assert_eq!(
            err,
            Error::NotConverged("source MBWR critical grid lower bound is uninitialized")
        );
        assert!(
            model
                .baseline_critical_grid(MbwrDensityParameters {
                    critical_temperature: 1e30,
                    ..data
                })
                .is_err()
        );
        assert!(
            model
                .baseline_critical_grid(MbwrDensityParameters {
                    critical_density: f64::NAN,
                    ..data
                })
                .is_err()
        );
    }
}
