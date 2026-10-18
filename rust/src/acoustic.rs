//! Isentropic stiffness from the full fluid SV equilibrium matrix. Unlike a
//! TP Schur complement, this also handles pure two-phase coexistence.
use super::*;

impl Cubic {
    /// Single-phase sound speed from T (K), total volume (m3) and mole numbers.
    /// Does not establish stability. Valid TV states can have negative pressure.
    pub fn speed_of_sound_tv(&self, t: f64, volume: f64, moles: &[f64]) -> Result<f64> {
        self.sound_tv_impl(t, volume, moles, false)
    }

    /// Source TV sound convention with total volume and mole numbers.
    pub fn baseline_speed_of_sound_tv(&self, t: f64, volume: f64, moles: &[f64]) -> Result<f64> {
        self.sound_tv_impl(t, volume, moles, true)
    }

    fn sound_tv_impl(&self, t: f64, volume: f64, moles: &[f64], baseline: bool) -> Result<f64> {
        positive(volume, "total volume")?;
        if moles.len() != self.components().len() || moles.iter().any(|n| !n.is_finite() || *n < 0.)
        {
            return Err(Error::InvalidInput(
                "finite nonnegative mole numbers matching the model required",
            ));
        }
        let total: f64 = moles.iter().sum();
        positive(total, "total mole number")?;
        let z: Vec<_> = moles.iter().map(|n| n / total).collect();
        let state = self.prepare(t, &z)?;
        if baseline {
            state.baseline_speed_of_sound(volume / total)
        } else {
            state.speed_of_sound(volume / total)
        }
    }

    /// Equilibrium fluid sound speed [m/s], with thermal, mechanical and
    /// compositional phase redistribution. Assumes a converged equilibrium
    /// result; it does not establish stability. Supports pure coexistence.
    pub fn equilibrium_speed_of_sound(&self, flash: &FlashResult) -> Result<f64> {
        self.equilibrium_speed_of_sound_with_phases(flash, [Phase::Liquid, Phase::Vapor])
    }

    /// Fluid phase roots in [liquid-array, vapor-array] order. Supports
    /// supplied two-liquid/two-vapor states; does not establish equilibrium.
    pub fn equilibrium_speed_of_sound_with_phases(
        &self,
        flash: &FlashResult,
        roots: [Phase; 2],
    ) -> Result<f64> {
        self.equilibrium_sound_impl(flash, false, roots)
    }

    /// Source sound convention: reported H2/vdW Cp, zero for invalid
    /// single-phase stiffness, and feed-liquid fallback for a failed or
    /// nonnegative two-phase SV matrix response. This does not prove stability.
    pub fn baseline_equilibrium_speed_of_sound(&self, flash: &FlashResult) -> Result<f64> {
        self.baseline_equilibrium_speed_of_sound_with_phases(flash, [Phase::Liquid, Phase::Vapor])
    }

    /// Baseline acoustic convention with [liquid-array, vapor-array] fluid
    /// root overrides. Single-phase flags ignore the overrides as in the fork.
    pub fn baseline_equilibrium_speed_of_sound_with_phases(
        &self,
        flash: &FlashResult,
        roots: [Phase; 2],
    ) -> Result<f64> {
        self.equilibrium_sound_impl(flash, true, roots)
    }

    fn equilibrium_sound_impl(
        &self,
        flash: &FlashResult,
        baseline: bool,
        roots: [Phase; 2],
    ) -> Result<f64> {
        let t = flash.temperature;
        let p = flash.pressure;
        positive(p, "pressure")?;
        let feed = self.prepare(t, &flash.feed)?;
        let mass: f64 = flash
            .feed
            .iter()
            .zip(self.components())
            .map(|(x, c)| x * c.molar_mass * 1e-3)
            .sum();
        if flash.phase != FlashPhase::TwoPhase {
            let root = if flash.phase == FlashPhase::Vapor {
                Phase::Vapor
            } else {
                Phase::Liquid
            };
            let v = feed.volume(p, root)?;
            return single_sound(&feed, v, mass, baseline);
        }
        let bv = flash.vapor_fraction;
        let bl = flash.liquid_fraction;
        if !bv.is_finite() || !bl.is_finite() || bv < 0. || bl < 0. || (bv + bl - 1.).abs() > 1e-10
        {
            return Err(Error::InvalidInput(
                "physical equilibrium phase fractions required",
            ));
        }
        // Retain sound_velocity_2ph's first-vapor / second-liquid ordering.
        let first = acoustic_phase(self, t, p, &flash.vapor, roots[1], baseline)?;
        let second = acoustic_phase(self, t, p, &flash.liquid, roots[0], baseline)?;
        let n = flash.feed.len();
        let dim = n + 2;
        let mut matrix = vec![0.; dim * dim];
        let volume = bv * first.volume + bl * second.volume;
        let cp = bl * second.cp + bv * first.cp;
        let vt = bl * second.vt + bv * first.vt;
        let vp = bl * second.vp + bv * first.vp;
        matrix[n * dim + n] = -cp / R;
        matrix[(n + 1) * dim + n + 1] = p * p * vp / (R * t);
        matrix[n * dim + n + 1] = p * vt / R;
        matrix[(n + 1) * dim + n] = p * vt / R;
        for i in 0..n {
            if flash.feed[i] > 0. {
                let dt = t * (second.dt[i] - first.dt[i]);
                let dp = p * (second.dp[i] - first.dp[i]);
                matrix[i * dim + n] = dt;
                matrix[n * dim + i] = dt;
                matrix[i * dim + n + 1] = dp;
                matrix[(n + 1) * dim + i] = dp;
                for j in 0..n {
                    if flash.feed[j] > 0. && bv * bl >= f64::MIN_POSITIVE {
                        matrix[i * dim + j] = 1. / (bv * bl)
                            * (-1. + bl * first.dn[i * n + j] + bv * second.dn[i * n + j]);
                    }
                }
                let divisor = flash.vapor[i] * flash.liquid[i] * bv * bl;
                if divisor > f64::MIN_POSITIVE {
                    matrix[i * dim + i] += flash.feed[i] / divisor;
                }
            } else {
                matrix[i * dim + i] = 1.;
            }
        }
        let mut rhs = vec![0.; dim];
        rhs[n + 1] = p / (R * t);
        let mut pivots = vec![0; dim];
        let mut work = vec![0.; dim];
        let solved =
            numerics::solve_symmetric_in_place(&mut matrix, &mut rhs, &mut pivots, &mut work);
        if matches!(solved, Ok(true)) && rhs[n + 1].is_finite() && rhs[n + 1] < 0. {
            let squared = -p * rhs[n + 1] * volume * volume / mass;
            if squared.is_finite() && squared > 0. {
                return Ok(squared.sqrt());
            }
        }
        if baseline {
            let v = feed.volume(p, Phase::Liquid)?;
            single_sound(&feed, v, mass, true)
        } else {
            Err(Error::NotConverged("equilibrium sound SV response"))
        }
    }
}

fn single_sound(state: &PreparedState<'_>, volume: f64, mass: f64, baseline: bool) -> Result<f64> {
    if !baseline {
        return state.speed_of_sound(volume);
    }
    let cp = state.baseline_enthalpy_derivative_tp(volume)?;
    let pressure = state.pressure(volume)?;
    let vt = -pressure.dt / pressure.dv;
    let vp = 1. / pressure.dv;
    let derivative = vp + state.temperature() * vt * vt / cp;
    if derivative < 0. {
        Ok((-volume * volume / (mass * derivative)).sqrt())
    } else {
        Ok(0.)
    }
}

struct AcousticPhase {
    volume: f64,
    cp: f64,
    vt: f64,
    vp: f64,
    dt: Vec<f64>,
    dp: Vec<f64>,
    dn: Vec<f64>,
}

fn acoustic_phase(
    model: &Cubic,
    t: f64,
    p: f64,
    x: &[f64],
    phase: Phase,
    baseline: bool,
) -> Result<AcousticPhase> {
    let state = model.prepare(t, x)?;
    let volume = state.volume(p, phase)?;
    let cp = if baseline {
        state.baseline_enthalpy_derivative_tp(volume)?
    } else {
        state.heat_capacities(volume)?.0
    };
    let pressure = state.pressure(volume)?;
    let n = x.len();
    let mut dt = vec![0.; n];
    let mut dp = vec![0.; n];
    let mut dn = vec![0.; n * n];
    state.fugacity_derivatives_into(volume, &mut dt, &mut dp, &mut dn)?;
    Ok(AcousticPhase {
        volume,
        cp,
        vt: -pressure.dt / pressure.dv,
        vp: 1. / pressure.dv,
        dt,
        dp,
        dn,
    })
}
