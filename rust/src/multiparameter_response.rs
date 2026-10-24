//! Native equilibrium response and acoustic properties for pure multiparameter fluids.
use super::*;
use crate::caloric::EquilibriumResponse;
use crate::{FlashPhase, FlashResult, Saturation};

impl Multiparameter {
    /// Equilibrium TP response, including the phase amount response.
    /// A pure two-phase TP state is constrained to saturation and returns `None`.
    pub fn equilibrium_response(&self, flash: &FlashResult) -> Result<Option<EquilibriumResponse>> {
        validate_flash(flash)?;
        if flash.phase == FlashPhase::TwoPhase {
            validate_coexistence(self.saturation_pressure(flash.temperature)?, flash)?;
            return Ok(None);
        }
        let prep = self.prepare(flash.temperature)?;
        let v = homogeneous_volume(&prep, flash)?;
        let p = prep.pressure(v)?;
        let cp = prep.heat_capacity_p(v)?;
        if p.dv == 0. || !p.dv.is_finite() {
            return Ok(None);
        }
        let result = EquilibriumResponse {
            heat_capacity: cp,
            volume_dt: -p.dt / p.dv,
            volume_dp: 1. / p.dv,
        };
        Ok((cp > 0.
            && [result.heat_capacity, result.volume_dt, result.volume_dp]
                .iter()
                .all(|x| x.is_finite()))
        .then_some(result))
    }

    /// Equilibrium sound speed, including quality redistribution on coexistence.
    pub fn equilibrium_speed_of_sound(&self, flash: &FlashResult) -> Result<f64> {
        self.equilibrium_speed_of_sound_impl(flash)
    }

    /// Compatibility with the source TP sound wrapper's component-database mass.
    pub fn baseline_equilibrium_speed_of_sound(&self, flash: &FlashResult) -> Result<f64> {
        validate_flash(flash)?;
        let source_mass = self.source_tp_molar_mass()?;
        let prep = self.prepare(flash.temperature)?;
        if flash.phase != FlashPhase::TwoPhase {
            let v = homogeneous_volume(&prep, flash)?;
            return baseline_homogeneous_sound(&prep, flash.pressure, v, source_mass);
        }
        let sat = self.saturation_pressure(flash.temperature)?;
        validate_coexistence(sat, flash)?;
        // The source TP wrapper obtains its two phase states through the
        // phase-seeded isobaric density iterator, rather than the simultaneous
        // native coexistence solve. Its density conditioning matters near Tc.
        let vl = prep.volume(flash.pressure, Phase::Liquid)?;
        let vv = prep.volume(flash.pressure, Phase::Vapor)?;
        let l = baseline_phase(&prep, flash.pressure, vl)?;
        let v = baseline_phase(&prep, flash.pressure, vv)?;
        let beta = flash.vapor_fraction;
        let mix = |a: f64, b: f64| (1. - beta) * a + beta * b;
        let t = flash.temperature;
        let p = flash.pressure;
        let r = self.gas_constant();
        // Original state_functions::getStateFuncMatrix('SV') and source
        // symmetric indefinite solve. The pure composition block is zero.
        let mut matrix = [
            0.,
            t * (v.phi_dt - l.phi_dt),
            p * (v.phi_dp - l.phi_dp),
            0.,
            -t * mix(l.entropy_dt, v.entropy_dt) / r,
            p * mix(l.volume_dt, v.volume_dt) / r,
            0.,
            0.,
            p * p * mix(l.volume_dp, v.volume_dp) / (r * t),
        ];
        let mut rhs = [0., 0., p / (r * t)];
        let mut pivots = [0; 3];
        let mut work = [0.; 3];
        let solved = crate::numerics::solve_symmetric_in_place(
            &mut matrix,
            &mut rhs,
            &mut pivots,
            &mut work,
        );
        if matches!(solved, Ok(true)) && rhs[2].is_finite() && rhs[2] < 0. {
            let volume = mix(l.volume, v.volume);
            let c2 = -p * rhs[2] * volume * volume / source_mass;
            if c2.is_finite() && c2 > 0. {
                return Ok(c2.sqrt());
            }
        }
        // The original routine falls back to the liquid TP sound when the
        // matrix is singular, nonfinite, or gives a nonnegative dP/dV.
        baseline_homogeneous_sound(&prep, p, vl, source_mass)
    }

    fn source_tp_molar_mass(&self) -> Result<f64> {
        if self.gas_constant == 1. && self.lj.is_none() {
            return Ok(self.molar_mass());
        }
        let component = self.lj.map(|p| p.component).unwrap_or(self.record.ident);
        Ok(crate::component::component(component)?.molar_mass * 1e-3)
    }

    fn equilibrium_speed_of_sound_impl(&self, flash: &FlashResult) -> Result<f64> {
        validate_flash(flash)?;
        let t = flash.temperature;
        let mass = self.molar_mass();
        positive(mass, "molar mass")?;
        if flash.phase != FlashPhase::TwoPhase {
            let prep = self.prepare(t)?;
            let v = homogeneous_volume(&prep, flash)?;
            return prep.speed_of_sound(v);
        }
        let sat = self.saturation_pressure(t)?;
        validate_coexistence(sat, flash)?;
        let beta = flash.vapor_fraction;
        let prep = self.prepare(t)?;
        coexistence_sound(
            &prep,
            sat.dpdt,
            sat.liquid_volume,
            sat.vapor_volume,
            beta,
            mass,
        )
    }
}

struct BaselinePhase {
    volume: f64,
    volume_dt: f64,
    volume_dp: f64,
    entropy_dt: f64,
    entropy_dp: f64,
    phi_dt: f64,
    phi_dp: f64,
}

fn baseline_phase(prep: &PreparedMultiparameter<'_>, p: f64, v: f64) -> Result<BaselinePhase> {
    let a = prep.residual_jet(v)?;
    let ideal = prep.ideal_jet(v)?;
    let t = prep.temperature();
    let r = prep.eos.gas_constant;
    let rt = r * t;
    // Keep the source's dimensional derivative grouping. Near Tc it differs
    // numerically from simplifying these expressions with D=1+2a_d+a_dd.
    let vv = (a.d_deltadelta + 2. * a.d_delta) / v.powi(2);
    let tv = a.d_deltatau / (t * v);
    let pv = -rt * (vv + 1. / v.powi(2));
    let pt = p / t - rt * tv;
    if !pv.is_finite() || pv >= 0. {
        return Err(Error::NoPhysicalRoot);
    }
    let vn = -(rt * (vv * v + 1. / v)) / pv;
    let phi_dt = (-a.d_tau / t - tv * v) + (1. - vn * pt / r) / t;
    let phi_dp = vn / rt - 1. / p;
    let at = -ideal.d_tau / t - a.d_tau / t;
    let att =
        (ideal.d_tautau + 2. * ideal.d_tau) / t.powi(2) + (a.d_tautau + 2. * a.d_tau) / t.powi(2);
    let entropy_dt = (-pt / pv) * pt - r * (2. * at + t * att);
    let z = p / ((1. / v) * rt);
    let zt = -z * (1. / t - (-pt / pv) / v);
    let zp = z * (1. / p + 1. / (pv * v));
    let result = BaselinePhase {
        volume: z * rt / p,
        volume_dt: zt * rt / p + z * r / p,
        volume_dp: zp * rt / p - z * rt / p.powi(2),
        entropy_dt,
        entropy_dp: pt / pv,
        phi_dt,
        phi_dp,
    };
    if [
        result.volume,
        result.volume_dt,
        result.volume_dp,
        result.entropy_dt,
        result.phi_dt,
        result.phi_dp,
    ]
    .iter()
    .all(|x| x.is_finite())
    {
        Ok(result)
    } else {
        Err(Error::NoPhysicalRoot)
    }
}

fn baseline_homogeneous_sound(
    prep: &PreparedMultiparameter<'_>,
    p: f64,
    v: f64,
    mass: f64,
) -> Result<f64> {
    let q = baseline_phase(prep, p, v)?;
    let dtdp = -q.entropy_dp / q.entropy_dt;
    let vp = q.volume_dp + q.volume_dt * dtdp;
    if vp < 0. {
        Ok((-q.volume * q.volume / (mass * vp)).sqrt())
    } else {
        Ok(0.)
    }
}

fn coexistence_sound(
    prep: &PreparedMultiparameter<'_>,
    slope: f64,
    vl: f64,
    vv: f64,
    beta: f64,
    mass: f64,
) -> Result<f64> {
    let cv = (1. - beta) * phase_equilibrium_cv(prep, slope, vl)?
        + beta * phase_equilibrium_cv(prep, slope, vv)?;
    let volume = (1. - beta) * vl + beta * vv;
    // S_V|T = P_sat' and S_T|V = Cv_eq/T avoid cancellation of
    // divergent phase and quality derivatives near the critical point.
    let c2 = volume * volume * prep.temperature() * slope.powi(2) / (mass * cv);
    if c2.is_finite() && c2 > 0. {
        Ok(c2.sqrt())
    } else {
        Err(Error::NoPhysicalRoot)
    }
}

fn phase_equilibrium_cv(prep: &PreparedMultiparameter<'_>, slope: f64, volume: f64) -> Result<f64> {
    let p = prep.pressure(volume)?;
    if p.dv >= 0. || !p.dv.is_finite() {
        return Err(Error::NoPhysicalRoot);
    }
    let cv = prep.heat_capacity_v(volume)? - prep.temperature() * (p.dt - slope).powi(2) / p.dv;
    if cv.is_finite() && cv > 0. {
        Ok(cv)
    } else {
        Err(Error::NoPhysicalRoot)
    }
}

fn validate_flash(f: &FlashResult) -> Result<()> {
    positive(f.temperature, "temperature")?;
    positive(f.pressure, "pressure")?;
    if f.feed.len() != 1 || f.feed[0] != 1. || f.liquid.len() != 1 || f.vapor.len() != 1 {
        return Err(Error::InvalidInput("pure normalized flash required"));
    }
    if f.phase == FlashPhase::TwoPhase {
        if f.liquid[0] != 1. || f.vapor[0] != 1. {
            return Err(Error::InvalidInput("pure normalized flash required"));
        }
        let b = f.vapor_fraction;
        let l = f.liquid_fraction;
        if !b.is_finite()
            || !l.is_finite()
            || !(0. ..=1.).contains(&b)
            || !(0. ..=1.).contains(&l)
            || (b + l - 1.).abs() > 1e-10
        {
            return Err(Error::InvalidInput(
                "physical equilibrium phase fractions required",
            ));
        }
    }
    Ok(())
}

fn validate_coexistence(sat: Saturation, f: &FlashResult) -> Result<()> {
    if (sat.pressure - f.pressure).abs() > 1e-8 * f.pressure {
        return Err(Error::InvalidInput(
            "flash pressure is inconsistent with pure-fluid saturation",
        ));
    }
    Ok(())
}

fn homogeneous_volume(prep: &PreparedMultiparameter<'_>, f: &FlashResult) -> Result<f64> {
    match f.phase {
        FlashPhase::Vapor => prep.volume(f.pressure, Phase::Vapor),
        FlashPhase::Liquid => prep.volume(f.pressure, Phase::Liquid),
        FlashPhase::Single => prep
            .volume(f.pressure, Phase::Liquid)
            .or_else(|_| prep.volume(f.pressure, Phase::Vapor)),
        FlashPhase::TwoPhase => Err(Error::InvalidInput("homogeneous flash required")),
    }
}
