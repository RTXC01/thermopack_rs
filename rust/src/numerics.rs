//! Small numerical kernels with bounded work. Caller-owned RR buffers avoid
//! allocations during repeated solves.

#[path = "symmetric.rs"]
mod symmetric;
pub use symmetric::solve_symmetric_in_place;

/// Real roots of x^3 + a*x^2 + b*x + c, sorted in ascending order.
/// Unused entries are NaN. Returns the number of distinct roots.
pub fn cubic_roots(a: f64, b: f64, c: f64) -> ([f64; 3], usize) {
    let p = b - a * a / 3.0;
    let q = (2.0 * a * a * a / 27.0) - a * b / 3.0 + c;
    let half_q = q / 2.0;
    let third_p = p / 3.0;
    // Expanded polynomial discriminant avoids cancellation between O(1)
    // depressed-cubic terms when two physical roots are very small.
    let discriminant = -(a * a * b * b - 4.0 * b * b * b - 4.0 * a * a * a * c - 27.0 * c * c
        + 18.0 * a * b * c)
        / 108.0;
    let mut roots = [f64::NAN; 3];
    let mut count = if discriminant > 0.0 {
        // Compute the larger cube-root term first to avoid cancellation.
        let u = (-half_q - discriminant.sqrt().copysign(half_q)).cbrt();
        let v = if u != 0.0 { -third_p / u } else { 0.0 };
        roots[0] = u + v - a / 3.0;
        1
    } else if discriminant == 0.0 {
        let u = (-half_q).cbrt();
        roots[0] = 2.0 * u - a / 3.0;
        roots[1] = -u - a / 3.0;
        if u == 0.0 { 1 } else { 2 }
    } else if p < 0.0 {
        let radius = 2.0 * (-third_p).sqrt();
        let theta = (-half_q / (-third_p * third_p * third_p).sqrt())
            .clamp(-1.0, 1.0)
            .acos()
            / 3.0;
        for (i, root) in roots.iter_mut().enumerate() {
            *root = radius * (theta - 2.0 * std::f64::consts::PI * i as f64 / 3.0).cos() - a / 3.0;
        }
        roots.sort_by(f64::total_cmp);
        3
    } else {
        roots[0] = -a / 3.0;
        1
    };
    for root in &mut roots[..count] {
        for _ in 0..2 {
            let residual = ((*root + a) * *root + b) * *root + c;
            let derivative = (3.0 * *root + 2.0 * a) * *root + b;
            if derivative.abs() > 1e-14 {
                *root -= residual / derivative;
            }
        }
    }
    if count == 3 {
        let largest_index = (0..3)
            .max_by(|&i, &j| roots[i].abs().total_cmp(&roots[j].abs()))
            .unwrap();
        let largest = roots[largest_index];
        let smallest = roots.iter().map(|v| v.abs()).fold(f64::INFINITY, f64::min);
        if largest.abs() > 100.0 * smallest && largest != 0.0 {
            // Deflate the well-resolved root. Use coefficients to recover
            // the small-root sum instead of subtracting nearly equal roots.
            let product = -c / largest;
            let linear = (product - b) / largest;
            let discriminant = linear * linear - 4.0 * product;
            if discriminant >= 0.0 {
                let q = -0.5 * (linear + discriminant.sqrt().copysign(linear));
                if q != 0.0 {
                    roots = [largest, q, product / q];
                }
            }
        }
    }
    roots[..count].sort_by(f64::total_cmp);
    let mut distinct = 0;
    for i in 0..count {
        if distinct == 0
            || (roots[i] - roots[distinct - 1]).abs()
                > 32.0 * f64::EPSILON * roots[i].abs().max(roots[distinct - 1].abs())
        {
            roots[distinct] = roots[i];
            distinct += 1;
        }
    }
    count = distinct;
    roots[count..].fill(f64::NAN);
    (roots, count)
}

/// Rachford-Rice vapor fraction for a given set of K values.
/// A result at 0 or 1 indicates no interior root, not a stability verdict.
pub fn rachford_rice(z: &[f64], k: &[f64]) -> crate::Result<f64> {
    if z.is_empty()
        || z.len() != k.len()
        || z.iter().any(|x| !x.is_finite() || *x < 0.0)
        || k.iter().any(|x| !x.is_finite() || *x <= 0.0)
        || (z.iter().sum::<f64>() - 1.0).abs() > 1e-10
    {
        return Err(crate::Error::InvalidInput(
            "normalized z and positive K of equal length required",
        ));
    }
    let evaluate = |beta: f64| {
        let mut f = 0.0;
        let mut df = 0.0;
        for (&zi, &ki) in z.iter().zip(k) {
            if zi == 0. {
                continue;
            }
            let ratio = (ki - 1.0) / ((1.0 - beta) + beta * ki);
            f += zi * ratio;
            df -= zi * ratio * ratio;
        }
        (f, df)
    };
    if evaluate(0.0).0 <= 0.0 {
        return Ok(0.0);
    }
    if evaluate(1.0).0 >= 0.0 {
        return Ok(1.0);
    }
    let (mut low, mut high, mut beta) = (0.0, 1.0, 0.5);
    for _ in 0..100 {
        let (f, df) = evaluate(beta);
        if f.abs() < 1e-13 || high - low < 1e-14 {
            return Ok(beta);
        }
        if f > 0.0 {
            low = beta;
        } else {
            high = beta;
        }
        let newton = beta - f / df;
        beta = if newton > low && newton < high {
            newton
        } else {
            (low + high) * 0.5
        };
    }
    Err(crate::Error::NotConverged("Rachford-Rice"))
}

/// Workspace/results for the fork's extended Rachford-Rice solver.
/// Fractions may lie outside `[0,1]` during negative-flash initialization.
#[derive(Debug, Clone)]
pub struct RachfordRiceState {
    pub vapor_fraction: f64,
    pub liquid_fraction: f64,
    pub liquid: Vec<f64>,
    pub vapor: Vec<f64>,
}
impl RachfordRiceState {
    pub fn new(components: usize) -> Self {
        Self {
            vapor_fraction: 0.5,
            liquid_fraction: 0.5,
            liquid: vec![0.; components],
            vapor: vec![0.; components],
        }
    }
}

/// Port of tp_solver::rr_solve, including extended bounds and the liquid-
/// fraction formulation near a vapor endpoint. True means an RR result was
/// obtained, not thermodynamic stability or fractions inside `[0,1]`.
///
/// False preserves arrays not written by the source's failure path. Their
/// contents must not be interpreted as a new solution. `single_step` performs
/// the source's one-Newton-step mode using the supplied vapor fraction, and
/// leaves liquid_fraction unchanged as in that routine.
pub fn rachford_rice_extended(
    z: &[f64],
    k: &[f64],
    state: &mut RachfordRiceState,
    single_step: bool,
) -> crate::Result<bool> {
    if z.is_empty()
        || z.len() != k.len()
        || state.liquid.len() != z.len()
        || state.vapor.len() != z.len()
        || z.iter().any(|v| !v.is_finite() || *v < 0.)
        || k.iter().any(|v| !v.is_finite() || *v < 0.)
        || (z.iter().sum::<f64>() - 1.).abs() > 1e-10
        || (single_step && !state.vapor_fraction.is_finite())
    {
        return Err(crate::Error::InvalidInput(
            "normalized z, finite nonnegative K and matching RR buffers required",
        ));
    }
    let evaluate = |beta: f64, vapor: bool| {
        let (mut g, mut dg) = (0., 0.);
        for (&zi, &ki) in z.iter().zip(k) {
            let divisor = if vapor {
                1. - beta + beta * ki
            } else {
                beta + (1. - beta) * ki
            };
            if divisor == 0. {
                return None;
            }
            g += zi * (ki - 1.) / divisor;
            let ratio = (ki - 1.) / divisor;
            dg += zi * (ratio * ratio);
        }
        Some((g, if vapor { -dg } else { dg }))
    };
    if single_step {
        let beta = state.vapor_fraction;
        let (g, dg) = evaluate(beta, true).ok_or(crate::Error::NoPhysicalRoot)?;
        let beta = beta - g / dg;
        for (i, (&zi, &ki)) in z.iter().zip(k).enumerate() {
            let divisor = 1. - beta + beta * ki;
            state.liquid[i] = (1. - beta) * zi / divisor;
            state.vapor[i] = beta * zi * ki / divisor;
        }
        state.vapor_fraction = state.vapor.iter().sum();
        for (x, y) in state.liquid.iter_mut().zip(&mut state.vapor) {
            *x /= 1. - state.vapor_fraction;
            *y /= state.vapor_fraction;
        }
        if !state.vapor_fraction.is_finite()
            || state
                .liquid
                .iter()
                .chain(&state.vapor)
                .any(|v| !v.is_finite())
        {
            return Err(crate::Error::NotConverged(
                "non-finite single-step Rachford-Rice result",
            ));
        }
        return Ok(true);
    }
    let smallest = k.iter().copied().fold(f64::INFINITY, f64::min);
    let largest = k.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let (mut low, mut high) = if smallest < 1. && largest > 1. {
        if (smallest - 1.).abs() < 10. * f64::EPSILON && (largest - 1.).abs() < 10. * f64::EPSILON {
            return Ok(false);
        }
        (-1. / (largest - 1.), 1. / (1. - smallest))
    } else {
        let endpoint = |beta| {
            z.iter()
                .zip(k)
                .map(|(z, k)| z * (k - 1.) / (1. - beta + beta * k))
                .sum::<f64>()
        };
        if !(endpoint(0.) > 0. && endpoint(1.) < 0.) {
            return Ok(false);
        }
        (0., 1.)
    };
    for (&zi, &ki) in z.iter().zip(k) {
        if ki > 1. {
            low = low.max((ki * zi - 1.) / (ki - 1.));
        } else if ki < 1. {
            high = high.min((1. - zi) / (1. - ki));
        }
    }
    let mut vapor = true;
    if high > 0.5 {
        let g = z
            .iter()
            .zip(k)
            .map(|(z, k)| z * (k - 1.) / (0.5 + 0.5 * k))
            .sum::<f64>();
        if g < 0. {
            high = 0.5;
        } else {
            vapor = false;
            low = 1. - high;
            high = 0.5;
        }
    }
    let mut beta = 0.5 * (high + low);
    state.vapor_fraction = beta;
    let mut converged = false;
    for _ in 0..1500 {
        let Some((g, dg)) = evaluate(beta, vapor) else {
            return Ok(false);
        };
        if (vapor && g > 0.) || (!vapor && g < 0.) {
            low = beta;
        } else if (vapor && g < 0.) || (!vapor && g > 0.) {
            high = beta;
        }
        let old = beta;
        beta -= g / dg;
        if beta < low || beta > high {
            beta = 0.5 * (low + high);
        }
        state.vapor_fraction = beta;
        if g.abs() < 100. * f64::EPSILON && (beta - old).abs() < 100. * f64::EPSILON {
            converged = true;
            break;
        }
    }
    state.liquid_fraction = if vapor { 1. - beta } else { beta };
    if !vapor {
        state.vapor_fraction = 1. - beta;
    }
    if !converged {
        return Ok(false);
    }
    for (i, (&zi, &ki)) in z.iter().zip(k).enumerate() {
        let divisor = if vapor {
            1. - beta + beta * ki
        } else {
            beta + (1. - beta) * ki
        };
        state.liquid[i] = zi / divisor;
        state.vapor[i] = ki * zi / divisor;
    }
    Ok(true)
}

/// Solve a dense row-major linear system in place by partial pivoting.
/// The matrix is overwritten with its elimination factors; rhs holds x.
pub fn solve_linear_in_place(matrix: &mut [f64], rhs: &mut [f64]) -> crate::Result<()> {
    let n = rhs.len();
    if n == 0 || matrix.len() != n * n || matrix.iter().chain(rhs.iter()).any(|x| !x.is_finite()) {
        return Err(crate::Error::InvalidInput(
            "finite square linear system required",
        ));
    }
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&i, &j| matrix[i * n + k].abs().total_cmp(&matrix[j * n + k].abs()))
            .unwrap();
        if matrix[pivot * n + k].abs() < f64::MIN_POSITIVE {
            return Err(crate::Error::NotConverged("singular linear system"));
        }
        if pivot != k {
            for j in 0..n {
                matrix.swap(k * n + j, pivot * n + j);
            }
            rhs.swap(k, pivot);
        }
        for i in k + 1..n {
            let factor = matrix[i * n + k] / matrix[k * n + k];
            matrix[i * n + k] = factor;
            for j in k + 1..n {
                matrix[i * n + j] -= factor * matrix[k * n + j];
            }
            rhs[i] -= factor * rhs[k];
        }
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            rhs[i] -= matrix[i * n + j] * rhs[j];
        }
        rhs[i] /= matrix[i * n + i];
    }
    if rhs.iter().any(|x| !x.is_finite()) {
        Err(crate::Error::NotConverged("non-finite linear solution"))
    } else {
        Ok(())
    }
}

pub(crate) fn modified_newton_direction(hessian: &[f64], gradient: &[f64]) -> Vec<f64> {
    let n = gradient.len();
    // Equilibrate before testing positive definiteness. Trace phase amounts
    // can put valid Hessian diagonal entries twenty orders of magnitude apart.
    let scales: Vec<_> = (0..n)
        .map(|i| hessian[i * n + i].abs().max(1e-30).sqrt())
        .collect();
    let mut shift = 0.0;
    for _ in 0..16 {
        let mut l = vec![0.0; n * n];
        let mut positive = true;
        for i in 0..n {
            for j in 0..=i {
                let mut value = hessian[i * n + j] / (scales[i] * scales[j]);
                if i == j {
                    value += shift;
                }
                for k in 0..j {
                    value -= l[i * n + k] * l[j * n + k];
                }
                if i == j {
                    if value <= 1e-15 || !value.is_finite() {
                        positive = false;
                        break;
                    }
                    l[i * n + j] = value.sqrt();
                } else {
                    l[i * n + j] = value / l[j * n + j];
                }
            }
            if !positive {
                break;
            }
        }
        if positive {
            let mut step: Vec<_> = gradient.iter().zip(&scales).map(|(g, s)| -g / s).collect();
            for i in 0..n {
                for j in 0..i {
                    step[i] -= l[i * n + j] * step[j];
                }
                step[i] /= l[i * n + i];
            }
            for i in (0..n).rev() {
                for j in i + 1..n {
                    step[i] -= l[j * n + i] * step[j];
                }
                step[i] /= l[i * n + i];
            }
            for (value, scale) in step.iter_mut().zip(&scales) {
                *value /= scale;
            }
            return step;
        }
        shift = if shift == 0.0 { 1e-10 } else { 10.0 * shift };
    }
    gradient
        .iter()
        .zip(&scales)
        .map(|(g, s)| -g / (s * s))
        .collect()
}

#[cfg(test)]
mod tests {
    #[test]
    fn equilibrated_newton_step_resolves_trace_scale_separation() {
        let hessian = [1.0, 0.0, 0.0, 1e20];
        let gradient = [-1.0, -1.0];
        let step = super::modified_newton_direction(&hessian, &gradient);
        assert!((step[0] - 1.0).abs() < 1e-14);
        assert!((step[1] * 1e20 - 1.0).abs() < 1e-14);
    }
}
