//! Fork-order LV flash, retaining its early SS return and Gibbs-change stopping.
//! Kept separate from the chemical-potential-qualified native flash.
use super::*;

impl Cubic {
    /// Follow the fork's TP iteration order and stopping conventions.
    ///
    /// Successive substitution may return a stationary split with higher Gibbs
    /// energy than the feed. The Newton fallback stops on relative Gibbs change
    /// and can return a vanishing phase with unclosed chemical potentials.
    /// Inspect `fugacity_residual`; use `tp_flash` for the stricter physical
    /// acceptance criteria. Full raw source parity is still being qualified.
    pub fn baseline_tp_flash(&self, t: f64, p: f64, z: &[f64]) -> Result<FlashResult> {
        positive(p, "pressure")?;
        let prepared = self.prepare(t, z)?;
        let (roots, count) = prepared.volumes(p)?;
        let mut volume = roots[0];
        let mut feed_phase = FlashPhase::Single;
        if count > 1 {
            let l = prepared.residual(roots[0])?;
            let v = prepared.residual(roots[count - 1])?;
            if l.enthalpy - t * l.entropy < v.enthalpy - t * v.entropy {
                feed_phase = FlashPhase::Liquid;
            } else {
                feed_phase = FlashPhase::Vapor;
                volume = roots[count - 1];
            }
        }
        let n = z.len();
        let mut feed_phi = vec![0.; n];
        prepared.ln_fugacity_tp_into(p, volume, &mut feed_phi)?;
        let feed_g = gibbs(z, &feed_phi);
        let mut k = self.wilson_k(t, p)?;
        let mut rr = RachfordRiceIteration::new(n);
        rr.liquid_ln_fugacity.fill(1e-6);
        rr.vapor_ln_fugacity.fill(1e-6);
        rr.gibbs = 10. * feed_g;
        let mut dk = vec![0.; n];
        let mut dkm1 = vec![0.; n];
        let mut dkm2 = vec![0.; n];
        let mut stability: Option<Stability> = None;
        let mut stability_k = vec![0.; n];
        for iteration in 1..=10 {
            let previous_g = rr.gibbs;
            let k_old = k.clone();
            let liquid_old = rr.liquid_ln_fugacity.clone();
            let vapor_old = rr.vapor_ln_fugacity.clone();
            let dem = iteration % 5 == 0;
            if dem {
                extrapolate_k(&dk, &dkm1, &dkm2, &rr, &mut k);
            }
            self.rr_successive_substitution_iteration(t, p, z, &k, &mut rr, false)?;
            if rr.converged {
                if !rr.has_solution {
                    return Err(Error::NotConverged(
                        "source TP converged without an RR solution",
                    ));
                }
                return raw_split(self, t, p, z, rr.split, iteration);
            }
            if dem && rr.gibbs > previous_g {
                self.rr_successive_substitution_iteration(t, p, z, &k_old, &mut rr, false)?;
            }
            k.copy_from_slice(&rr.equilibrium_ratios);
            if rr.has_solution {
                dkm2.copy_from_slice(&dkm1);
                dkm1.copy_from_slice(&dk);
                for i in 0..n {
                    dk[i] = rr.liquid_ln_fugacity[i] - rr.vapor_ln_fugacity[i] - liquid_old[i]
                        + vapor_old[i];
                }
            }
            let smaller_g = rr.gibbs < feed_g;
            let valid_beta = rr.split.vapor_fraction > 0. && rr.split.vapor_fraction < 1.;
            let check = (iteration == 5 && (!smaller_g || !valid_beta)) || !rr.has_solution;
            if (check || (iteration == 10 && (!smaller_g || !valid_beta))) && stability.is_none() {
                let found = self.stability(t, p, z)?;
                if found.stable {
                    return Ok(single(t, p, z, feed_phase, volume, iteration));
                }
                let limit = -1000. * f64::EPSILON;
                let liquid_unstable = found.liquid_trial.tpd < limit;
                let vapor_unstable = found.vapor_trial.tpd < limit;
                for i in 0..n {
                    let liquid = if liquid_unstable {
                        found.liquid_trial.ln_fugacity[i]
                    } else {
                        feed_phi[i]
                    };
                    let vapor = if vapor_unstable {
                        found.vapor_trial.ln_fugacity[i]
                    } else {
                        feed_phi[i]
                    };
                    k[i] = safe_exp(liquid - vapor);
                }
                if liquid_unstable && vapor_unstable {
                    rr.liquid_ln_fugacity
                        .copy_from_slice(&found.liquid_trial.ln_fugacity);
                    rr.vapor_ln_fugacity
                        .copy_from_slice(&found.vapor_trial.ln_fugacity);
                } else {
                    rr.liquid_ln_fugacity.fill(1e-6);
                    rr.vapor_ln_fugacity.fill(1e-6);
                }
                stability_k.copy_from_slice(&k);
                stability = Some(found);
            } else if check {
                return source_newton(
                    self,
                    t,
                    p,
                    z,
                    rr.split,
                    &stability_k,
                    stability.as_ref(),
                    feed_g,
                    iteration,
                );
            }
            if iteration == 10 {
                return source_newton(
                    self,
                    t,
                    p,
                    z,
                    rr.split,
                    &stability_k,
                    stability.as_ref(),
                    feed_g,
                    iteration,
                );
            }
        }
        unreachable!("fixed ten-iteration SS loop always enters Newton")
    }
}

fn safe_exp(value: f64) -> f64 {
    value.min(f64::MAX.ln()).exp()
}

fn extrapolate_k(dk: &[f64], m1: &[f64], m2: &[f64], state: &RachfordRiceIteration, k: &mut [f64]) {
    let dot = |x: &[f64], y: &[f64]| x.iter().zip(y).map(|(a, b)| a * b).sum::<f64>();
    let (b01, b02, b11, b12, b22) = (
        dot(dk, m1),
        dot(dk, m2),
        dot(m1, m1),
        dot(m1, m2),
        dot(m2, m2),
    );
    let divisor = b11 * b22 - b12 * b12;
    if divisor.abs() > 1e-20 {
        let mu1 = (b02 * b12 - b01 * b22) / divisor;
        let mu2 = (b01 * b12 - b02 * b11) / divisor;
        for i in 0..k.len() {
            k[i] = safe_exp(
                state.liquid_ln_fugacity[i] - state.vapor_ln_fugacity[i]
                    + (dk[i] - m1[i] * mu2) / (1. - mu1 - mu2),
            );
        }
    }
}

fn single(
    t: f64,
    p: f64,
    z: &[f64],
    phase: FlashPhase,
    volume: f64,
    iterations: usize,
) -> FlashResult {
    let vapor = phase == FlashPhase::Vapor;
    FlashResult {
        temperature: t,
        pressure: p,
        feed: z.to_vec(),
        phase,
        liquid: if vapor { vec![0.; z.len()] } else { z.to_vec() },
        vapor: if vapor { z.to_vec() } else { vec![0.; z.len()] },
        vapor_fraction: match phase {
            FlashPhase::Liquid => 0.,
            FlashPhase::Vapor => 1.,
            _ => -1.,
        },
        liquid_fraction: match phase {
            FlashPhase::Liquid => 1.,
            FlashPhase::Vapor => 0.,
            _ => -1.,
        },
        liquid_volume: if vapor { None } else { Some(volume) },
        vapor_volume: if vapor { Some(volume) } else { None },
        iterations,
        fugacity_residual: 0.,
    }
}

fn raw_split(
    model: &Cubic,
    t: f64,
    p: f64,
    z: &[f64],
    split: numerics::RachfordRiceState,
    iterations: usize,
) -> Result<FlashResult> {
    let mut x = split.liquid.clone();
    let mut y = split.vapor.clone();
    normalize(&mut x)?;
    normalize(&mut y)?;
    let l = evaluate(model, t, p, &x, Phase::Liquid, false)?;
    let v = evaluate(model, t, p, &y, Phase::Vapor, false)?;
    let mut residual = 0.0_f64;
    for i in 0..z.len() {
        if z[i] > 0. {
            residual = residual.max((y[i].ln() + v.phi[i] - x[i].ln() - l.phi[i]).abs());
        }
    }
    Ok(FlashResult {
        temperature: t,
        pressure: p,
        feed: z.to_vec(),
        phase: FlashPhase::TwoPhase,
        liquid: split.liquid,
        vapor: split.vapor,
        vapor_fraction: split.vapor_fraction,
        liquid_fraction: split.liquid_fraction,
        liquid_volume: Some(l.volume),
        vapor_volume: Some(v.volume),
        iterations,
        fugacity_residual: residual,
    })
}

#[allow(clippy::too_many_arguments)]
fn source_newton(
    model: &Cubic,
    t: f64,
    p: f64,
    z: &[f64],
    mut split: numerics::RachfordRiceState,
    stability_k: &[f64],
    stability: Option<&Stability>,
    feed_g: f64,
    ss_iterations: usize,
) -> Result<FlashResult> {
    let n = z.len();
    let mut vapor = vec![0.; n];
    let mut liquid = vec![0.; n];
    if !(0. ..=1.).contains(&split.vapor_fraction) {
        let success = numerics::rachford_rice_extended(z, stability_k, &mut split, false)?;
        if !success || !(0. ..=1.).contains(&split.vapor_fraction) {
            let stability = stability.ok_or(Error::NotConverged(
                "source TP lacks a stability initializer",
            ))?;
            if split.vapor_fraction <= 0. {
                split.vapor_fraction = 0.001;
                for i in 0..n {
                    liquid[i] = ((1. - split.vapor_fraction)
                        * stability.liquid_trial.composition[i])
                        .min((1. - 1e-6) * z[i]);
                    vapor[i] = z[i] - liquid[i];
                }
            } else {
                split.vapor_fraction = 0.999;
                for i in 0..n {
                    vapor[i] = (split.vapor_fraction * stability.vapor_trial.composition[i])
                        .min((1. - 1e-6) * z[i]);
                    liquid[i] = z[i] - vapor[i];
                }
            }
        } else {
            for i in 0..n {
                vapor[i] = split.vapor_fraction * split.vapor[i];
                liquid[i] = z[i] - vapor[i];
            }
        }
        split.liquid_fraction = liquid.iter().sum();
        for i in 0..n {
            split.liquid[i] = liquid[i] / split.liquid_fraction;
            split.vapor[i] = vapor[i] / split.vapor_fraction;
        }
    } else {
        for i in 0..n {
            vapor[i] = split.vapor_fraction * split.vapor[i];
            liquid[i] = z[i] - vapor[i];
        }
    }
    // Preserve the smaller phase amount, as tp_solver::setV does.
    for i in 0..n {
        if vapor[i] >= liquid[i] {
            liquid[i] = split.liquid_fraction * split.liquid[i];
            vapor[i] = z[i] - liquid[i];
        } else {
            vapor[i] = split.vapor_fraction * split.vapor[i];
            liquid[i] = z[i] - vapor[i];
        }
    }
    let evaluate_amounts =
        |vapor: &[f64], liquid: &[f64], jacobian: bool| -> Result<(f64, Vec<f64>, Vec<f64>)> {
            let bv: f64 = vapor.iter().sum();
            let bl: f64 = liquid.iter().sum();
            let x: Vec<_> = liquid.iter().map(|v| v / bl).collect();
            let y: Vec<_> = vapor.iter().map(|v| v / bv).collect();
            let l = evaluate(model, t, p, &x, Phase::Liquid, jacobian)?;
            let v = evaluate(model, t, p, &y, Phase::Vapor, jacobian)?;
            let mut g = 0.;
            let mut gradient = vec![0.; n];
            let mut hessian = vec![0.; if jacobian { n * n } else { 0 }];
            for i in 0..n {
                if z[i] > 0. {
                    g += vapor[i] * (y[i].ln() + v.phi[i]) + liquid[i] * (x[i].ln() + l.phi[i]);
                    gradient[i] = y[i].ln() - x[i].ln() + v.phi[i] - l.phi[i];
                }
                if jacobian {
                    for j in 0..n {
                        if z[i] > 0. && z[j] > 0. {
                            hessian[i * n + j] = 1. / (bv * bl)
                                * (-1. + bv * l.dn[i * n + j] + bl * v.dn[i * n + j]);
                        }
                    }
                    hessian[i * n + i] += if z[i] > 0. {
                        z[i] / (x[i] * y[i] * bv * bl)
                    } else {
                        1.
                    };
                }
            }
            if !g.is_finite() {
                return Err(Error::NotConverged("non-finite source Gibbs objective"));
            }
            Ok((g, gradient, hessian))
        };
    let mut pivots = vec![0; n];
    let mut symmetric_work = vec![0.; n];
    for iteration in 1..=1000 {
        let (old_g, gradient, mut hessian) = evaluate_amounts(&vapor, &liquid, true)?;
        let mut step: Vec<_> = gradient.iter().map(|g| -g).collect();
        if n == 1 {
            step[0] /= if hessian[0] <= 0. {
                1000. * f64::EPSILON
            } else {
                hessian[0]
            };
        } else {
            let original = hessian.clone();
            if !numerics::solve_symmetric_in_place(
                &mut hessian,
                &mut step,
                &mut pivots,
                &mut symmetric_work,
            )? || step.iter().zip(&gradient).map(|(s, g)| s * g).sum::<f64>() >= 0.
            {
                step = source_modified_direction(&original, &gradient)?;
            }
        }
        let mut scale = 1.0_f64;
        for i in 0..n {
            if vapor[i] < liquid[i] {
                if vapor[i] + step[i] < 0. {
                    scale = scale.min(-vapor[i] / step[i]);
                } else if z[i] - vapor[i] - step[i] < 0. {
                    scale = scale.min((z[i] - vapor[i]) / step[i]);
                }
            } else if liquid[i] - step[i] < 0. {
                scale = scale.min(liquid[i] / step[i]);
            } else if z[i] - liquid[i] + step[i] < 0. {
                scale = scale.min((liquid[i] - z[i]) / step[i]);
            }
        }
        if scale < 1. {
            for s in &mut step {
                *s *= scale * (1. - 1e-12);
            }
        }
        let old_v = vapor.clone();
        let old_l = liquid.clone();
        let update = |alpha: f64, vapor: &mut [f64], liquid: &mut [f64]| {
            for i in 0..n {
                if old_v[i] >= old_l[i] {
                    liquid[i] = old_l[i] - alpha * step[i];
                    vapor[i] = z[i] - liquid[i];
                } else {
                    vapor[i] = old_v[i] + alpha * step[i];
                    liquid[i] = z[i] - vapor[i];
                }
            }
        };
        let mut alpha = 1.;
        update(alpha, &mut vapor, &mut liquid);
        let mut g = evaluate_amounts(&vapor, &liquid, false)?.0;
        let wolfe = 0.001 * step.iter().zip(&gradient).map(|(s, g)| s * g).sum::<f64>();
        for i in 1..=2 {
            if g >= old_g + alpha * wolfe {
                alpha *= 0.5_f64.powi(i);
                update(alpha, &mut vapor, &mut liquid);
                g = evaluate_amounts(&vapor, &liquid, false)?.0;
            } else {
                break;
            }
        }
        let tolerance = 300. * f64::EPSILON;
        let dummy = if old_g.abs() < 0.99 * tolerance {
            1.
        } else {
            0.
        };
        let error = ((g - old_g + dummy) / (old_g + dummy) - dummy).abs();
        if error < tolerance {
            if g - feed_g > 10. * f64::EPSILON {
                return Err(Error::NotConverged(
                    "source TP Newton Gibbs energy exceeds feed",
                ));
            }
            split.vapor_fraction = vapor.iter().sum();
            split.liquid_fraction = liquid.iter().sum();
            for i in 0..n {
                split.liquid[i] = liquid[i] / split.liquid_fraction;
                split.vapor[i] = vapor[i] / split.vapor_fraction;
            }
            return raw_split(model, t, p, z, split, ss_iterations + iteration);
        }
    }
    Err(Error::NotConverged(
        "source TP Newton Gibbs-change stopping",
    ))
}

// optimizer::modefiedCholeskyDecompNoPerm / solveModefiedCholesky.
fn source_modified_direction(a: &[f64], gradient: &[f64]) -> Result<Vec<f64>> {
    let n = gradient.len();
    let mut l = vec![0.; n * n];
    let mut c = vec![0.; n * n];
    let mut d = vec![0.; n];
    let mut gamma = 0.0_f64;
    let mut xi = 0.0_f64;
    for i in 0..n {
        l[i * n + i] = 1.;
        c[i * n + i] = a[i * n + i];
        gamma = gamma.max(a[i * n + i].abs());
        for j in 0..n {
            if i != j {
                xi = xi.max(a[i * n + j].abs());
            }
        }
    }
    let mu = 100. * f64::EPSILON;
    let delta = mu * (gamma + xi).max(1.);
    let beta = gamma.max(mu).max(xi / ((n * n) as f64 - 1.).sqrt()).sqrt();
    for j in 0..n {
        for s in 0..j {
            l[j * n + s] = c[j * n + s] / d[s];
        }
        let mut theta = 0.0_f64;
        for i in j + 1..n {
            c[i * n + j] = a[i * n + j] - (0..j).map(|s| l[j * n + s] * c[i * n + s]).sum::<f64>();
            theta = theta.max(c[i * n + j].abs());
        }
        d[j] = c[j * n + j].abs().max((theta / beta).powi(2)).max(delta);
        for i in j + 1..n {
            c[i * n + i] -= c[i * n + j].powi(2) / d[j];
        }
    }
    let mut step: Vec<_> = gradient.iter().map(|g| -g).collect();
    for i in 0..n {
        for j in 0..i {
            step[i] -= l[i * n + j] * step[j];
        }
    }
    for i in 0..n {
        step[i] /= d[i];
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            step[i] -= l[j * n + i] * step[j];
        }
    }
    if step.iter().any(|v| !v.is_finite()) {
        return Err(Error::NotConverged("source modified Cholesky direction"));
    }
    Ok(step)
}

#[cfg(test)]
mod tests {
    use super::source_modified_direction;

    #[test]
    fn source_ldl_preserves_positive_curvature_and_modifies_indefinite_curvature() {
        let direction = source_modified_direction(&[4., 1., 1., 3.], &[1., -2.]).unwrap();
        assert!((direction[0] + 5. / 11.).abs() < 1e-14);
        assert!((direction[1] - 9. / 11.).abs() < 1e-14);
        for a in [[1., 2., 2., -1.], [0.; 4]] {
            let gradient = [1., -2.];
            let direction = source_modified_direction(&a, &gradient).unwrap();
            assert!(
                direction
                    .iter()
                    .zip(gradient)
                    .map(|(d, g)| d * g)
                    .sum::<f64>()
                    < 0.
            );
        }
    }
}
