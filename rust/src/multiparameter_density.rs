//! Original multiparameter_base::densitySolver with local, bounded execution.
use super::*;

const P_EPS: f64 = f64::EPSILON * 1e8;
const RHO_EPS: f64 = f64::EPSILON * 1e6;
const EXECUTION_LIMIT: usize = 10000;
/// Relative density step below which the iteration counts as settled.
///
/// The original stops when the density change is below `releps_rho*rhoOld`,
/// which is 2.2e-10 relative. The native residual kernel agrees with the
/// source only to about 1e-14 relative, and the Newton step divides the
/// pressure error by dP/drho, so a state close to the critical point cannot
/// reach the original criterion. The relaxed criterion accepts a settled
/// state whose pressure is closed and whose branch is mechanically stable.
const SETTLED_RHO: f64 = 1e-6;
/// Number of consecutive iterations with a closed pressure that a relaxed
/// acceptance requires.
const CLOSED_STREAK: usize = 2;

#[derive(Debug, Clone, Copy)]
pub struct MultiparameterDensityIteration {
    pub density: f64,
    /// Raw source phase flag. A phase switch may change the requested flag.
    pub phase: i32,
    pub pressure: f64,
    pub pressure_density_derivative: f64,
    pub iterations: usize,
    pub phase_switches: usize,
    /// Source ierr=1 is set at iteration 200, but does not break its loop.
    pub source_error_code: i32,
    pub pressure_converged: bool,
    pub mechanically_stable: bool,
    pub relative_pressure_error: f64,
    /// The state was accepted with the relaxed density criterion because the
    /// original criterion is unreachable at a small dP/drho.
    pub relaxed_convergence: bool,
}

#[derive(Default)]
struct Search {
    rho: f64,
    p: f64,
    dp: f64,
    phase: i32,
    switches: usize,
    curvature: f64,
    p_min: f64,
    dp_min: f64,
}

impl PreparedMultiparameter<'_> {
    /// Original mp_pressure uses rho/rhoc directly (not a molar-volume round trip).
    /// The analytic zero-density pressure/derivative limit is also supported.
    pub fn pressure_density(&self, rho: f64) -> Result<(f64, f64)> {
        let a = self.residual_at_delta(rho / self.eos.record.rc)?;
        let rt = self.eos.gas_constant * self.temperature;
        let p = rho * rt * (1. + a.d_delta);
        let dp = rt * (1. + 2. * a.d_delta + a.d_deltadelta);
        if p.is_finite() && dp.is_finite() {
            Ok((p, dp))
        } else {
            Err(Error::InvalidInput("non-finite multiparameter pressure"))
        }
    }

    pub fn baseline_density(&self, p: f64, phase: Phase) -> Result<MultiparameterDensityIteration> {
        let code = match phase {
            Phase::Liquid => 1,
            Phase::Vapor => 2,
            Phase::MinimumGibbs => 3,
        };
        self.baseline_density_with_phase_code(p, code)
    }

    /// Original densitySolver with continueOnError=false and ierr present.
    /// Retains source raw phase flags, pOld/step stopping tests, curvature
    /// guards and phase switching. Undefined/non-finite states and 10000-step
    /// or liquid-seed doubling exhaustion are errors, not silent last iterates.
    pub fn baseline_density_with_phase_code(
        &self,
        p_spec: f64,
        phase: i32,
    ) -> Result<MultiparameterDensityIteration> {
        positive(p_spec, "pressure")?;
        let mut s = Search {
            phase,
            ..Search::default()
        };
        self.density_seed(p_spec, &mut s)?;
        let mut source_error_code = 0;
        let mut closed_streak = 0usize;
        let mut previous_step = f64::NAN;
        for iteration in 1..=EXECUTION_LIMIT {
            let old_rho = s.rho;
            let old_p = s.p;
            let old_dp = s.dp;
            s.rho -= (s.p - p_spec) / s.dp;
            if !s.rho.is_finite() {
                return Err(Error::NotConverged("non-finite source density step"));
            }
            let switch = if s.rho < 0. {
                true
            } else {
                (s.p, s.dp) = self.pressure_density(s.rho)?;
                s.p < s.p_min
                    || s.dp < s.dp_min
                    || s.curvature * (s.rho - old_rho) * (s.dp - old_dp)
                        < -2e-10 * (s.rho * s.dp).abs()
            };
            if switch {
                s.phase = if s.phase == 2 { 1 } else { 2 };
                s.switches += 1;
                self.density_seed(p_spec, &mut s)?;
            }
            // Source CONTINUE is a no-op statement, not CYCLE: switching
            // still increments iter and tests old_p against the reset seed.
            let pressure_closed = (p_spec - s.p).abs() <= P_EPS * p_spec + 1e-6;
            closed_streak = if pressure_closed {
                closed_streak + 1
            } else {
                0
            };
            if (p_spec - old_p).abs() < P_EPS * old_p + 1e-6
                && (s.rho - old_rho).abs() < RHO_EPS * old_rho
            {
                positive(s.rho, "source density return")?;
                let relative = (s.p - p_spec).abs() / p_spec;
                return Ok(MultiparameterDensityIteration {
                    density: s.rho,
                    phase: s.phase,
                    pressure: s.p,
                    pressure_density_derivative: s.dp,
                    iterations: iteration,
                    phase_switches: s.switches,
                    source_error_code,
                    pressure_converged: (s.p - p_spec).abs() <= P_EPS * p_spec + 1e-6,
                    mechanically_stable: s.dp > 0.,
                    relative_pressure_error: relative,
                    relaxed_convergence: false,
                });
            }
            // The original density criterion is unreachable when the Newton
            // step amplifies the kernel difference, which happens where
            // dP/drho is small. A stalled iteration with a closed pressure on
            // a mechanically stable branch is accepted instead, and the
            // relaxation is reported. A contracting iteration reaches the
            // original criterion first and never takes this path.
            let step = if s.rho != 0. {
                (s.rho - old_rho).abs() / s.rho.abs()
            } else {
                f64::INFINITY
            };
            let stalled = previous_step.is_finite() && step > 0.5 * previous_step;
            previous_step = step;
            if step < SETTLED_RHO && stalled && closed_streak >= CLOSED_STREAK && s.dp > 0. {
                positive(s.rho, "source density return")?;
                let relative = (s.p - p_spec).abs() / p_spec;
                return Ok(MultiparameterDensityIteration {
                    density: s.rho,
                    phase: s.phase,
                    pressure: s.p,
                    pressure_density_derivative: s.dp,
                    iterations: iteration,
                    phase_switches: s.switches,
                    source_error_code,
                    pressure_converged: true,
                    mechanically_stable: true,
                    relative_pressure_error: relative,
                    relaxed_convergence: true,
                });
            }
            if iteration == 200 {
                source_error_code = 1;
            }
        }
        Err(Error::NotConverged(
            "source multiparameter iteration execution limit",
        ))
    }

    fn density_seed(&self, p_spec: f64, s: &mut Search) -> Result<()> {
        if self.temperature > self.eos.record.tc || s.phase == 2 {
            s.curvature = if self.temperature > self.eos.record.tc {
                0.
            } else {
                -1.
            };
            s.rho = p_spec / (self.temperature * self.eos.gas_constant);
            (s.p, s.dp) = self.pressure_density(s.rho)?;
        } else {
            let liquid = s.phase == 1
                || (s.phase == 4
                    && matches!(
                        self.eos.record.ancillary,
                        AncillarySatDensity::TriplePointRatio { .. }
                    ));
            if !liquid {
                return Err(Error::InvalidInput(
                    "source ancillary rejects the density phase flag",
                ));
            }
            s.curvature = 1.;
            s.rho = self.eos.record.rc * self.sat_delta_estimate(Phase::Liquid)? * 1.01;
            (s.p, s.dp) = self.pressure_density(s.rho)?;
            let mut doublings = 0;
            while s.p < 0. || s.dp < 0. {
                if doublings >= 1024 {
                    return Err(Error::NotConverged(
                        "source liquid density seed doubling limit",
                    ));
                }
                s.rho *= 2.;
                s.curvature = 0.;
                doublings += 1;
                (s.p, s.dp) = self.pressure_density(s.rho)?;
            }
        }
        if s.switches == 2 {
            // Preserve the source comparison, even though its STOP message
            // says T<Ttriple while the actual condition is T>Ttriple.
            if self.temperature > self.eos.record.t_triple {
                return Err(Error::NotConverged("source density phase-switch stop"));
            }
            s.curvature = 0.;
            s.p_min = -1e100;
            s.dp_min = -1e100;
        }
        Ok(())
    }
}
