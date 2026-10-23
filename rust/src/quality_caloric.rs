//! Caloric constraints on the local fixed-phase-amount coexistence equations.
use super::*;
use crate::R;

/// Last specification in the fixed-quality coexistence corrector.
#[derive(Debug, Clone, Copy)]
pub enum QualitySpecification {
    Temperature(f64),
    Pressure(f64),
    LogK {
        component: usize,
        value: f64,
    },
    /// Molar entropy, J/(mol K).
    Entropy(f64),
    /// Molar enthalpy, J/mol.
    Enthalpy(f64),
}
impl From<EnvelopeConstraint> for QualitySpecification {
    fn from(c: EnvelopeConstraint) -> Self {
        match c {
            EnvelopeConstraint::Temperature(t) => Self::Temperature(t),
            EnvelopeConstraint::Pressure(p) => Self::Pressure(p),
            EnvelopeConstraint::LogK { component, value } => Self::LogK { component, value },
        }
    }
}
#[derive(Debug, Clone, Copy)]
pub struct QualityCaloricProperties {
    pub volume: f64,
    pub enthalpy: f64,
    pub entropy: f64,
    pub internal_energy: f64,
}
pub(super) struct PhaseCaloricJet {
    pub value: f64,
    pub dt: f64,
    pub dp: f64,
    pub partial: Vec<f64>,
}
pub(super) fn phase_jet(
    model: &Cubic,
    jet: &Jet<'_>,
    volume: f64,
    pressure: f64,
    entropy: bool,
    baseline: bool,
    derivatives: bool,
) -> Result<PhaseCaloricJet> {
    let state = jet.prepared_state()?;
    let t = state.temperature();
    let tv = state.caloric_tv(volume)?;
    let value = if entropy {
        tv.entropy
    } else {
        tv.internal_energy + pressure * volume
    };
    let mut result = PhaseCaloricJet {
        value,
        dt: 0.,
        dp: 0.,
        partial: Vec::new(),
    };
    if !derivatives {
        return Ok(result);
    }
    let p = state.pressure(volume)?;
    let vt = -p.dt / p.dv;
    let cp = if baseline {
        state.baseline_enthalpy_derivative_tp(volume)?
    } else {
        tv.heat_capacity_v - t * p.dt * p.dt / p.dv
    };
    result.dt = if entropy { cp / t } else { cp };
    result.dp = if entropy { -vt } else { volume - t * vt };
    result.partial = vec![0.; model.components().len()];
    for (i, &x) in state.composition().iter().enumerate() {
        if x == 0. {
            continue;
        }
        let ideal = model.ideal_gas(i)?.properties(t, pressure)?;
        result.partial[i] = if entropy {
            ideal.entropy - R * (x.ln() + jet.phi[i] + t * jet.dt[i])
        } else {
            ideal.enthalpy - R * t * t * jet.dt[i]
        };
    }
    if !value.is_finite()
        || !result.dt.is_finite()
        || !result.dp.is_finite()
        || result.partial.iter().any(|v| !v.is_finite())
    {
        return Err(Error::NotConverged("nonfinite quality caloric derivatives"));
    }
    Ok(result)
}
impl Cubic {
    /// Molar H/S/U/V at the reported phase volumes and pressure. Includes the
    /// active ideal-gas reference levels and the specified phase amounts.
    pub fn quality_properties(&self, state: &QualityState) -> Result<QualityCaloricProperties> {
        check_quality(state.vapor_fraction)?;
        positive(state.pressure, "quality state pressure")?;
        self.prepare(state.temperature, &state.feed)?;
        let left = self
            .prepare(state.temperature, &state.liquid)?
            .caloric_tv(state.liquid_volume)?;
        let right = self
            .prepare(state.temperature, &state.vapor)?
            .caloric_tv(state.vapor_volume)?;
        let b = state.vapor_fraction;
        let volume = (1. - b) * state.liquid_volume + b * state.vapor_volume;
        let internal_energy = (1. - b) * left.internal_energy + b * right.internal_energy;
        Ok(QualityCaloricProperties {
            volume,
            enthalpy: internal_energy + state.pressure * volume,
            entropy: (1. - b) * left.entropy + b * right.entropy,
            internal_energy,
        })
    }
    /// Warm correction at fixed feed and phase amount with any supported
    /// quality specification. Keeps the original feed without lever-rule drift.
    pub fn refine_quality_specification(
        &self,
        previous: &QualityState,
        specification: QualitySpecification,
    ) -> Result<QualityState> {
        self.prepare(previous.temperature, &previous.liquid)?;
        self.prepare(previous.temperature, &previous.vapor)?;
        self.quality_state_with_specification(
            &previous.feed,
            previous.vapor_fraction,
            specification,
            &previous.guess(),
        )
    }
    pub(super) fn pure_quality_caloric(
        &self,
        z: &[f64],
        beta: f64,
        target: f64,
        entropy: bool,
        guess: &QualityGuess,
    ) -> Result<QualityState> {
        if !target.is_finite() {
            return Err(Error::InvalidInput(
                "finite pure quality caloric target required",
            ));
        }
        let critical = crate::critical::critical_point(self, z)?;
        let lo = self.temperature_bounds().0;
        let hi = self
            .temperature_bounds()
            .1
            .min(critical.temperature * (1. - 1e-7));
        if hi <= lo {
            return Err(Error::InvalidInput(
                "pure caloric coexistence temperature interval",
            ));
        }
        let evaluate = |t: f64| -> Result<(QualityState, f64)> {
            let point = self.quality_pure(z, beta, EnvelopeConstraint::Temperature(t))?;
            let properties = self.quality_properties(&point)?;
            Ok((
                point,
                if entropy {
                    properties.entropy
                } else {
                    properties.enthalpy
                } - target,
            ))
        };
        let t = guess.temperature.clamp(lo, hi);
        let (initial, f0) = evaluate(t)?;
        let tolerance = 1e-10 * target.abs().max(1.);
        if f0.abs() < tolerance {
            return Ok(initial);
        }
        let mut step = (0.01 * t).max(0.01);
        let mut bracket = None;
        for _ in 0..40 {
            for next in [(t - step).max(lo), (t + step).min(hi)] {
                if next == t {
                    continue;
                }
                if let Ok((point, f)) = evaluate(next) {
                    if f.abs() < tolerance {
                        return Ok(point);
                    }
                    if f.signum() != f0.signum() {
                        bracket = Some(if next < t {
                            (next, t, f, f0)
                        } else {
                            (t, next, f0, f)
                        });
                        break;
                    }
                }
            }
            if bracket.is_some() {
                break;
            }
            if t - step <= lo && t + step >= hi {
                break;
            }
            step *= 1.6;
        }
        let (mut low, mut high, mut fl, _) =
            bracket.ok_or(Error::NotConverged("pure quality caloric bracket"))?;
        for iteration in 0..100 {
            let t = 0.5 * (low + high);
            let (mut point, f) = evaluate(t)?;
            if f.abs() < tolerance {
                point.iterations = iteration + 1;
                return Ok(point);
            }
            if f.signum() == fl.signum() {
                low = t;
                fl = f;
            } else {
                high = t;
            }
            if high - low < 8. * f64::EPSILON * t {
                break;
            }
        }
        Err(Error::NotConverged("pure quality caloric correction"))
    }
}
