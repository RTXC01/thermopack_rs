//! Local pressure/entropy equations and homogeneous isentropic continuation.
//! Equations from `src/spinodal.f90::tv_meta_ps` and `map_meta_isentrope`.
//!
//! ```
//! use thermopack_rs::{Cubic, CubicKind, Phase};
//! use thermopack_rs::metastable::{IsentropeOptions, IsentropeEnd};
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let model = Cubic::new(CubicKind::PengRobinson, &["C1"])?;
//! let prepared = model.prepare(180.0, &[1.0])?;
//! let volume = prepared.volume(6e6, Phase::Liquid)?;
//! let entropy = prepared.caloric_tv(volume)?.entropy;
//! let seed = model.homogeneous_ps(6e6, entropy, &[1.0], 180.0, volume)?;
//! let curve = model.homogeneous_isentrope(&[1.0], entropy, seed,
//!     IsentropeOptions { points: 20, ..Default::default() })?;
//! assert_eq!(curve.end, IsentropeEnd::MinimumPressure);
//! # Ok(())
//! # }
//! ```
use crate::spinodal_properties::{
    SpinodalProperty, baseline_mixture_spinodal_intersection, baseline_pure_spinodal_intersection,
};
use crate::{Cubic, Error, FlashPhase, Phase, Result, error::positive};

/// A homogeneous PS solve or a point on an isentrope. No global phase flash is
/// implicit in this state. Negative homogeneous pressures are allowed.
#[derive(Debug, Clone, Copy)]
pub struct HomogeneousPsState {
    pub temperature: f64,
    pub volume: f64,
    /// Specified/reported pressure (Pa); actual EOS error is recorded below.
    pub pressure: f64,
    /// Entropy evaluated at the returned T/V (J/mol/K).
    pub entropy: f64,
    pub pressure_residual: f64,
    pub entropy_residual: f64,
    pub iterations: usize,
}

/// Termination of a homogeneous isentrope.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsentropeEnd {
    MinimumPressure,
    Spinodal,
    PointLimit,
}

#[derive(Debug, Clone)]
pub struct HomogeneousIsentrope {
    pub points: Vec<HomogeneousPsState>,
    pub end: IsentropeEnd,
}

/// Pressure grid and local-stability termination for a physical trace.
#[derive(Debug, Clone, Copy)]
pub struct IsentropeOptions {
    pub minimum_pressure: f64,
    pub points: usize,
    pub stop_at_spinodal: bool,
}
impl Default for IsentropeOptions {
    fn default() -> Self {
        Self {
            minimum_pressure: 1e5,
            points: 50,
            stop_at_spinodal: true,
        }
    }
}

fn state_at(
    model: &Cubic,
    z: &[f64],
    p: f64,
    s: f64,
    t: f64,
    v: f64,
    iterations: usize,
) -> Result<HomogeneousPsState> {
    let c = model.prepare(t, z)?.caloric_tv(v)?;
    Ok(HomogeneousPsState {
        temperature: t,
        volume: v,
        pressure: p,
        entropy: c.entropy,
        pressure_residual: c.pressure - p,
        entropy_residual: c.entropy - s,
        iterations,
    })
}

impl Cubic {
    /// Solve P(T,V)=p and S(T,V)=s from a positive local T/V seed. The analytic
    /// Jacobian uses S_T=Cv/T and S_V=P_T. This follows a local homogeneous
    /// root, including metastable or unstable roots, without a phase flash.
    pub fn homogeneous_ps(
        &self,
        p: f64,
        s: f64,
        z: &[f64],
        t_guess: f64,
        v_guess: f64,
    ) -> Result<HomogeneousPsState> {
        if !p.is_finite() || !s.is_finite() {
            return Err(Error::InvalidInput(
                "finite homogeneous PS specifications required",
            ));
        }
        positive(t_guess, "homogeneous PS temperature seed")?;
        positive(v_guess, "homogeneous PS volume seed")?;
        let (sp, ss) = (p.abs().max(1e5), s.abs().max(1.0));
        let mut q = [t_guess.ln(), v_guess.ln()];
        let evaluate = |q: [f64; 2]| -> Result<([f64; 2], [[f64; 2]; 2])> {
            let (t, v) = (q[0].exp(), q[1].exp());
            let state = self.prepare(t, z)?;
            let pressure = state.pressure(v)?;
            let c = state.caloric_tv(v)?;
            Ok((
                [(pressure.value - p) / sp, (c.entropy - s) / ss],
                [
                    [t * pressure.dt / sp, v * pressure.dv / sp],
                    [c.heat_capacity_v / ss, v * pressure.dt / ss],
                ],
            ))
        };
        for iteration in 0..100 {
            let (f, [[a, b], [c, d]]) = evaluate(q)?;
            let norm = f[0].abs().max(f[1].abs());
            if norm < 1e-10 {
                return state_at(self, z, p, s, q[0].exp(), q[1].exp(), iteration);
            }
            let determinant = a * d - b * c;
            if determinant == 0.0 || !determinant.is_finite() {
                return Err(Error::NotConverged("homogeneous PS Jacobian"));
            }
            let mut delta = [
                (b * f[1] - d * f[0]) / determinant,
                (c * f[0] - a * f[1]) / determinant,
            ];
            let largest = delta[0].abs().max(delta[1].abs());
            if largest > 0.4 {
                for x in &mut delta {
                    *x *= 0.4 / largest;
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
                return Err(Error::NotConverged("homogeneous PS backtracking"));
            }
        }
        Err(Error::NotConverged("homogeneous PS solve"))
    }

    /// Continue a pressure/entropy root from an explicit homogeneous seed.
    /// With stability termination enabled the seed must be locally stable;
    /// the first bracketed zero Helmholtz eigenvalue is refined as endpoint.
    pub fn homogeneous_isentrope(
        &self,
        z: &[f64],
        s: f64,
        seed: HomogeneousPsState,
        options: IsentropeOptions,
    ) -> Result<HomogeneousIsentrope> {
        if options.points < 2
            || !options.minimum_pressure.is_finite()
            || options.minimum_pressure >= seed.pressure
        {
            return Err(Error::InvalidInput(
                "decreasing isentrope pressure grid required",
            ));
        }
        let mut current =
            self.homogeneous_ps(seed.pressure, s, z, seed.temperature, seed.volume)?;
        let mut lambda = if options.stop_at_spinodal {
            stability(self, z, current)?
        } else {
            1.0
        };
        if lambda <= 0.0 {
            return Err(Error::InvalidInput(
                "locally stable isentrope seed required",
            ));
        }
        let initial = current.pressure;
        let step = (options.minimum_pressure - initial) / (options.points - 1) as f64;
        let mut points = vec![current];
        for i in 1..options.points {
            let target = if i + 1 == options.points {
                options.minimum_pressure
            } else {
                initial + step * i as f64
            };
            let mut reached = false;
            for _ in 0..128 {
                let mut pressure = target;
                let mut next = None;
                for _ in 0..24 {
                    let (t, v) = predict(self, z, current, pressure)?;
                    if let Ok(point) = self.homogeneous_ps(pressure, s, z, t, v)
                        && (point.temperature / current.temperature).ln().abs() <= 0.4
                        && (point.volume / current.volume).ln().abs() <= 0.4
                    {
                        next = Some(point);
                        break;
                    }
                    pressure = 0.5 * (pressure + current.pressure);
                }
                let Some(point) = next else {
                    return Err(Error::NotConverged("isentrope continuation step"));
                };
                let next_lambda = if options.stop_at_spinodal {
                    stability(self, z, point)?
                } else {
                    1.0
                };
                if lambda > 0.0 && next_lambda <= 0.0 {
                    let endpoint = spinodal_endpoint(self, z, s, current, point)?;
                    points.push(endpoint);
                    return Ok(HomogeneousIsentrope {
                        points,
                        end: IsentropeEnd::Spinodal,
                    });
                }
                current = point;
                lambda = next_lambda;
                if pressure == target {
                    points.push(current);
                    reached = true;
                    break;
                }
            }
            if !reached {
                return Err(Error::NotConverged("isentrope pressure subdivision"));
            }
        }
        Ok(HomogeneousIsentrope {
            points,
            end: IsentropeEnd::MinimumPressure,
        })
    }

    /// The fork's complete mapping convention: equilibrium PS initialization,
    /// majority-phase feed volume, and a source-reported spinodal endpoint.
    /// The first and last reported coordinates may not close the target entropy.
    /// Native interior correctors solve the homogeneous PS equations directly;
    /// original nonlinear trajectories and failure flags are not reproduced.
    pub fn baseline_metastable_isentrope(
        &self,
        z: &[f64],
        initial_pressure: f64,
        s: f64,
        minimum_pressure: f64,
        points: usize,
    ) -> Result<HomogeneousIsentrope> {
        positive(initial_pressure, "isentrope initial pressure")?;
        positive(minimum_pressure, "source isentrope minimum pressure")?;
        if points == 0 || minimum_pressure >= initial_pressure || !s.is_finite() {
            return Err(Error::InvalidInput("source isentrope specification"));
        }
        let endpoint = if z.iter().filter(|&&x| x > 0.0).count() == 1 {
            baseline_pure_spinodal_intersection(
                self,
                z,
                SpinodalProperty::Entropy,
                s,
                minimum_pressure,
                None,
                None,
            )?
        } else {
            baseline_mixture_spinodal_intersection(
                self,
                z,
                SpinodalProperty::Entropy,
                s,
                minimum_pressure,
                None,
            )?
        };
        let flash = self.ps_flash(initial_pressure, z, s)?;
        let phase = match flash.phase {
            FlashPhase::TwoPhase => {
                if flash.liquid_fraction > flash.vapor_fraction {
                    Phase::Liquid
                } else {
                    Phase::Vapor
                }
            }
            FlashPhase::Vapor => Phase::Vapor,
            _ => Phase::Liquid,
        };
        let v = self
            .prepare(flash.temperature, z)?
            .baseline_volume(initial_pressure, phase)?;
        let mut current = state_at(self, z, initial_pressure, s, flash.temperature, v, 0)?;
        let mut states = vec![current];
        if points == 1 {
            return Ok(HomogeneousIsentrope {
                points: states,
                end: IsentropeEnd::PointLimit,
            });
        }
        let end_pressure = endpoint.map_or(minimum_pressure, |p| p.pressure);
        let step = (end_pressure - initial_pressure) / (points - 1) as f64;
        for i in 1..points {
            if i + 1 == points
                && end_pressure != minimum_pressure
                && let Some(p) = endpoint
            {
                states.push(state_at(
                    self,
                    z,
                    p.pressure,
                    s,
                    p.temperature,
                    p.volume,
                    0,
                )?);
                break;
            }
            let pressure = current.pressure + step;
            let (t, v) = predict(self, z, current, pressure)?;
            current = self.homogeneous_ps(pressure, s, z, t, v)?;
            states.push(current);
        }
        Ok(HomogeneousIsentrope {
            points: states,
            end: if endpoint.is_some() && end_pressure != minimum_pressure {
                IsentropeEnd::Spinodal
            } else {
                IsentropeEnd::MinimumPressure
            },
        })
    }
}

fn predict(model: &Cubic, z: &[f64], old: HomogeneousPsState, p: f64) -> Result<(f64, f64)> {
    let state = model.prepare(old.temperature, z)?;
    let pressure = state.pressure(old.volume)?;
    let c = state.caloric_tv(old.volume)?;
    let determinant = pressure.dt * pressure.dt - pressure.dv * c.heat_capacity_v / old.temperature;
    let dt = pressure.dt / determinant / old.temperature * (p - old.pressure);
    let dv = -c.heat_capacity_v / old.temperature / determinant / old.volume * (p - old.pressure);
    let (t, v) = (old.temperature * dt.exp(), old.volume * dv.exp());
    if !t.is_finite() || !v.is_finite() || t <= 0.0 || v <= 0.0 {
        Ok((old.temperature, old.volume))
    } else {
        Ok((t, v))
    }
}

fn stability(model: &Cubic, z: &[f64], point: HomogeneousPsState) -> Result<f64> {
    Ok(model
        .prepare(point.temperature, z)?
        .mixture_stability(point.volume)?
        .eigenvalue)
}

fn spinodal_endpoint(
    model: &Cubic,
    z: &[f64],
    s: f64,
    mut high: HomogeneousPsState,
    mut low: HomogeneousPsState,
) -> Result<HomogeneousPsState> {
    for _ in 0..70 {
        let pressure = 0.5 * (high.pressure + low.pressure);
        let (t, v) = predict(model, z, high, pressure)?;
        let point = model.homogeneous_ps(pressure, s, z, t, v)?;
        let lambda = stability(model, z, point)?;
        if lambda.abs() < 1e-9 {
            return Ok(point);
        }
        if lambda > 0.0 {
            high = point;
        } else {
            low = point;
        }
    }
    Err(Error::NotConverged("isentrope spinodal endpoint"))
}
