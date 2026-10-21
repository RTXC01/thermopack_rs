//! Local Helmholtz spinodal points and volume continuation of their curves.
//! Equations from `src/spinodal.f90`; native solvers do not call Fortran.
use crate::{Cubic, Error, Phase, Result, error::positive, spinodal::SpinodalPoint};

fn eigenvalue(model: &Cubic, t: f64, v: f64, z: &[f64]) -> Result<f64> {
    Ok(model.prepare(t, z)?.mixture_stability(v)?.eigenvalue)
}

/// Refine the pressure and zero minimum stability eigenvalue from a local seed.
/// Positive pressure is required. This does not classify global phase stability
/// or select a liquid/vapor branch independently of the supplied T/V seed.
pub fn spinodal_point_with_guess(
    model: &Cubic,
    p: f64,
    z: &[f64],
    t: f64,
    v: f64,
) -> Result<SpinodalPoint> {
    positive(p, "spinodal pressure")?;
    positive(t, "spinodal temperature guess")?;
    positive(v, "spinodal volume guess")?;
    let mut q = [t.ln(), v.ln()];
    let evaluate = |q: [f64; 2]| -> Result<([f64; 2], crate::Pressure)> {
        let state = model.prepare(q[0].exp(), z)?;
        let volume = q[1].exp();
        let pressure = state.pressure(volume)?;
        Ok((
            [
                state.mixture_stability(volume)?.eigenvalue,
                (pressure.value - p) / p,
            ],
            pressure,
        ))
    };
    for _ in 0..100 {
        let (f, pressure) = evaluate(q)?;
        let norm = f[0].abs().max(f[1].abs());
        if norm < 2e-10 {
            return Ok(SpinodalPoint {
                temperature: q[0].exp(),
                volume: q[1].exp(),
                pressure: pressure.value,
            });
        }
        let step = 2e-5;
        let mut gradient = [0.0; 2];
        for j in 0..2 {
            let mut up = q;
            let mut down = q;
            up[j] += step;
            down[j] -= step;
            gradient[j] = (evaluate(up)?.0[0] - evaluate(down)?.0[0]) / (2.0 * step);
        }
        let a = gradient[0];
        let b = gradient[1];
        let c = q[0].exp() * pressure.dt / p;
        let d = q[1].exp() * pressure.dv / p;
        let determinant = a * d - b * c;
        if !determinant.is_finite() || determinant.abs() < 1e-18 {
            return Err(Error::NotConverged("spinodal pressure Jacobian"));
        }
        let mut delta = [
            (b * f[1] - d * f[0]) / determinant,
            (c * f[0] - a * f[1]) / determinant,
        ];
        let largest = delta[0].abs().max(delta[1].abs());
        if largest > 0.25 {
            for x in &mut delta {
                *x *= 0.25 / largest;
            }
        }
        let mut fraction = 1.0;
        let mut accepted = false;
        for _ in 0..24 {
            let trial = [q[0] + fraction * delta[0], q[1] + fraction * delta[1]];
            if let Ok((r, _)) = evaluate(trial)
                && r[0].abs().max(r[1].abs()) < norm
            {
                q = trial;
                accepted = true;
                break;
            }
            fraction *= 0.5;
        }
        if !accepted {
            return Err(Error::NotConverged("spinodal pressure backtracking"));
        }
    }
    Err(Error::NotConverged("spinodal pressure"))
}

pub(crate) fn mixture_spinodal_point(
    model: &Cubic,
    p: f64,
    phase: Phase,
    z: &[f64],
) -> Result<SpinodalPoint> {
    let mechanical = crate::spinodal::mechanical_spinodal_point(model, p, phase, z)?;
    let point = spinodal_point_with_guess(model, p, z, mechanical.temperature, mechanical.volume)?;
    Ok(point)
}

/// Reproduce the fork's mixture initialization and pressure-walk direction.
///
/// The original chooses a volume-step sign from the pressure error alone.
/// On a vapor branch this can walk through the turning region and return the
/// liquid-side limit. This explicit compatibility API preserves that selection;
/// [`crate::spinodal::spinodal_point`] selects the physical mechanical seed.
/// Iteration trajectories are not identical and work is bounded to 1000 steps.
pub fn baseline_mixture_spinodal_point(
    model: &Cubic,
    p: f64,
    phase: Phase,
    z: &[f64],
) -> Result<SpinodalPoint> {
    positive(p, "spinodal pressure")?;
    if phase == Phase::MinimumGibbs {
        return Err(Error::InvalidInput("a liquid or vapor branch is required"));
    }
    if z.iter().filter(|&&x| x > 0.0).count() == 1 {
        return crate::spinodal::spinodal_point(model, p, phase, z);
    }
    let critical = model.pseudo_critical(z)?;
    let saturation = model.pseudo_saturation_temperature(p, z)?;
    let t = saturation.temperature + if phase == Phase::Liquid { 0.1 } else { -0.1 };
    let v = crate::spinodal::limit_volume(model, t, critical.volume, phase == Phase::Liquid, z)?;
    let mut point = spinodal_temperature_at_volume(model, v, z, t)?;
    let sign = if point.pressure < p { 1.0 } else { -1.0 };
    for _ in 0..1000 {
        let v = point.volume * (sign * 0.05_f64).exp();
        let next = spinodal_temperature_at_volume(model, v, z, point.temperature)?;
        if (point.pressure - p) * (next.pressure - p) <= 0.0 {
            return spinodal_point_with_guess(model, p, z, next.temperature, next.volume);
        }
        point = next;
    }
    Err(Error::NotConverged("baseline spinodal pressure walk"))
}

/// Solve a local zero eigenvalue at fixed volume from a temperature seed.
pub fn spinodal_temperature_at_volume(
    model: &Cubic,
    v: f64,
    z: &[f64],
    t_guess: f64,
) -> Result<SpinodalPoint> {
    scalar_solve(model, z, t_guess, v, false)
}

/// Solve a local zero eigenvalue at fixed temperature from a volume seed.
/// Negative homogeneous pressures are valid and are returned as evaluated.
pub fn spinodal_volume_at_temperature(
    model: &Cubic,
    t: f64,
    z: &[f64],
    v_guess: f64,
) -> Result<SpinodalPoint> {
    scalar_solve(model, z, t, v_guess, true)
}

fn scalar_solve(
    model: &Cubic,
    z: &[f64],
    t: f64,
    v: f64,
    volume_unknown: bool,
) -> Result<SpinodalPoint> {
    positive(t, "spinodal temperature")?;
    positive(v, "spinodal volume")?;
    let mut q = if volume_unknown { v.ln() } else { t.ln() };
    let evaluate = |q: f64| {
        eigenvalue(
            model,
            if volume_unknown { t } else { q.exp() },
            if volume_unknown { q.exp() } else { v },
            z,
        )
    };
    for _ in 0..80 {
        let f = evaluate(q)?;
        if f.abs() < 2e-10 {
            let (t, v) = if volume_unknown {
                (t, q.exp())
            } else {
                (q.exp(), v)
            };
            return Ok(SpinodalPoint {
                temperature: t,
                volume: v,
                pressure: model.prepare(t, z)?.pressure(v)?.value,
            });
        }
        let h = 2e-5;
        let derivative = (evaluate(q + h)? - evaluate(q - h)?) / (2.0 * h);
        if derivative == 0.0 || !derivative.is_finite() {
            return Err(Error::NotConverged("spinodal scalar derivative"));
        }
        let step = (-f / derivative).clamp(-0.25, 0.25);
        let mut fraction = 1.0;
        let mut accepted = false;
        for _ in 0..24 {
            let trial = q + fraction * step;
            if let Ok(value) = evaluate(trial)
                && value.abs() < f.abs()
            {
                q = trial;
                accepted = true;
                break;
            }
            fraction *= 0.5;
        }
        if !accepted {
            return Err(Error::NotConverged("spinodal scalar backtracking"));
        }
    }
    Err(Error::NotConverged("spinodal scalar solve"))
}

/// Controls a local stability curve traced from the liquid branch toward vapor.
#[derive(Debug, Clone, Copy)]
pub struct SpinodalTraceOptions {
    pub initial_pressure: f64,
    pub initial_liquid_temperature: Option<f64>,
    pub minimum_vapor_temperature: Option<f64>,
    /// Positive increment in ln(molar volume). Default 0.02.
    pub log_volume_step: f64,
    pub maximum_points: usize,
}
impl Default for SpinodalTraceOptions {
    fn default() -> Self {
        Self {
            initial_pressure: 1e5,
            initial_liquid_temperature: None,
            minimum_vapor_temperature: None,
            log_volume_step: 0.02,
            maximum_points: 1000,
        }
    }
}
/// Why a bounded spinodal trace ended. A point limit is an explicit partial curve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpinodalTraceEnd {
    MinimumTemperature,
    InitialPressure,
    PointLimit,
}
#[derive(Debug, Clone)]
pub struct SpinodalCurve {
    pub points: Vec<SpinodalPoint>,
    pub end: SpinodalTraceEnd,
}

/// Trace a local spinodal at fixed composition by increasing volume, including
/// temperature turning points. Endpoints are refined to the requested minimum
/// temperature or return pressure. Sampling matches the source mixture volume
/// increment; pure curves use this same continuation rather than the source's
/// two 100-point temperature grids. No global phase-stability claim is made.
pub fn trace_spinodal(
    model: &Cubic,
    z: &[f64],
    options: SpinodalTraceOptions,
) -> Result<SpinodalCurve> {
    trace_spinodal_impl(model, z, options, false)
}

/// Mixture curve with the fork's pressure-reporting convention: its first
/// pressure is the input pressure even with an explicit initial temperature;
/// a temperature-limited last point retains the pre-refinement pressure.
/// Physical coordinates and bounded continuation still use the native solver.
/// The source's separate pure 200-point grid is not reproduced here.
pub fn baseline_trace_spinodal(
    model: &Cubic,
    z: &[f64],
    options: SpinodalTraceOptions,
) -> Result<SpinodalCurve> {
    if z.iter().filter(|&&x| x > 0.0).count() < 2 {
        return Err(Error::Unsupported("source pure spinodal sampling"));
    }
    trace_spinodal_impl(model, z, options, true)
}

fn trace_spinodal_impl(
    model: &Cubic,
    z: &[f64],
    options: SpinodalTraceOptions,
    baseline: bool,
) -> Result<SpinodalCurve> {
    positive(options.initial_pressure, "spinodal initial pressure")?;
    positive(options.log_volume_step, "spinodal volume step")?;
    if options.log_volume_step > 0.5 || options.maximum_points < 2 {
        return Err(Error::InvalidInput("spinodal trace bounds"));
    }
    let minimum = options.minimum_vapor_temperature.unwrap_or_else(|| {
        0.6 * z
            .iter()
            .zip(model.components())
            .map(|(z, c)| z * c.critical_temperature)
            .sum::<f64>()
    });
    positive(minimum, "spinodal minimum temperature")?;
    let mut point = if baseline {
        baseline_mixture_spinodal_point(model, options.initial_pressure, Phase::Liquid, z)?
    } else {
        crate::spinodal::spinodal_point(model, options.initial_pressure, Phase::Liquid, z)?
    };
    if let Some(t) = options.initial_liquid_temperature {
        point = spinodal_volume_at_temperature(model, t, z, point.volume)?;
    }
    let mut first = point;
    if baseline {
        first.pressure = options.initial_pressure;
    }
    let mut points = vec![first];
    for _ in 1..options.maximum_points {
        let v = (point.volume.ln() + options.log_volume_step).exp();
        let h = 2e-5;
        let gt = (eigenvalue(model, point.temperature * (1.0 + h), point.volume, z)?
            - eigenvalue(model, point.temperature * (1.0 - h), point.volume, z)?)
            / (2.0 * h * point.temperature);
        let gv = (eigenvalue(model, point.temperature, point.volume * (1.0 + h), z)?
            - eigenvalue(model, point.temperature, point.volume * (1.0 - h), z)?)
            / (2.0 * h * point.volume);
        let guess = point.temperature - gv * (v - point.volume) / gt;
        let next = spinodal_temperature_at_volume(model, v, z, guess)?;
        if next.temperature < minimum && next.temperature < point.temperature {
            let mut last = spinodal_volume_at_temperature(model, minimum, z, next.volume)?;
            if last.volume <= point.volume {
                return Err(Error::NotConverged("spinodal temperature endpoint branch"));
            }
            if baseline {
                last.pressure = next.pressure;
            }
            points.push(last);
            return Ok(SpinodalCurve {
                points,
                end: SpinodalTraceEnd::MinimumTemperature,
            });
        }
        if next.pressure < options.initial_pressure && next.pressure < point.pressure {
            let last = spinodal_point_with_guess(
                model,
                options.initial_pressure,
                z,
                next.temperature,
                next.volume,
            )?;
            if last.volume <= point.volume {
                return Err(Error::NotConverged("spinodal pressure endpoint branch"));
            }
            points.push(last);
            return Ok(SpinodalCurve {
                points,
                end: SpinodalTraceEnd::InitialPressure,
            });
        }
        points.push(next);
        point = next;
    }
    Ok(SpinodalCurve {
        points,
        end: SpinodalTraceEnd::PointLimit,
    })
}
