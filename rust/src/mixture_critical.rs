//! Local mixture critical points from the Helmholtz stability matrix.
//!
//! Implements the conditions used by `critical::calcCriticalTV`: the smallest
//! eigenvalue of sqrt(z) A_nn/(RT) sqrt(z), and the third derivative along its
//! null direction, must vanish. These differ from a fixed-composition pressure
//! inflection. A local critical point need not be globally stable against other
//! phases; critical loci and critical end points require additional equations.

use crate::{Cubic, Error, PreparedState, Result, critical::CriticalPoint, error::positive};

/// Minimum local Helmholtz stability eigenpair at fixed T and total V.
#[derive(Debug, Clone)]
pub struct MixtureStability {
    pub eigenvalue: f64,
    /// Unit direction in sqrt(z)-scaled coordinates, with absent species zero.
    pub eigenvector: Vec<f64>,
}

impl PreparedState<'_> {
    /// Evaluate local mixture stability without the critical third derivative.
    /// A positive eigenvalue means local stability at fixed T/V; this does not
    /// establish global equilibrium against finite composition changes.
    pub fn mixture_stability(&self, volume: f64) -> Result<MixtureStability> {
        let z = self.composition();
        let n = z.len();
        let mut h = vec![0.0; n * n];
        self.log_fugacity_tv_into(volume, &mut vec![0.0; n], &mut h)?;
        let active: Vec<_> = z
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| (v > 0.0).then_some(i))
            .collect();
        let (eigenvalue, u) = scaled_eigenpair(z, &active, &h)?;
        let mut eigenvector = vec![0.0; n];
        for (i, &a) in active.iter().enumerate() {
            eigenvector[a] = u[i];
        }
        Ok(MixtureStability {
            eigenvalue,
            eigenvector,
        })
    }
}

fn scaled_eigenpair(z: &[f64], active: &[usize], hessian: &[f64]) -> Result<(f64, Vec<f64>)> {
    let n = z.len();
    let m = active.len();
    let mut matrix = vec![0.0; m * m];
    for (i, &a) in active.iter().enumerate() {
        for (j, &b) in active.iter().enumerate() {
            matrix[i * m + j] =
                z[a].sqrt() * z[b].sqrt() * 0.5 * (hessian[a * n + b] + hessian[b * n + a]);
        }
    }
    smallest_eigenpair(&mut matrix, m)
}

/// Local stability of a one-mole homogeneous state at fixed T and total V.
#[derive(Debug, Clone)]
pub struct MixtureCriticalConditions {
    /// Smallest eigenvalue of sqrt(z) A_nn/(RT) sqrt(z).
    pub eigenvalue: f64,
    /// Unit eigenvector in scaled composition coordinates. Absent species are
    /// exactly zero. Its largest absolute component is chosen positive.
    pub eigenvector: Vec<f64>,
    /// Third derivative of A/(RT) along d_i=sqrt(z_i)*eigenvector_i, with
    /// total volume fixed. Classic mixing uses a scalar hyperdual derivative;
    /// Wong-Sandler uses a five-point derivative of the analytic mole Hessian.
    /// This is not a pressure derivative at fixed composition.
    pub third_derivative: f64,
}

impl Cubic {
    /// Evaluate the local mixture critical conditions, including the active
    /// composition direction. Classic and Wong-Sandler mixing are supported.
    pub fn mixture_critical_conditions(
        &self,
        temperature: f64,
        volume: f64,
        composition: &[f64],
    ) -> Result<MixtureCriticalConditions> {
        let n = self.components().len();
        let state = self.prepare(temperature, composition)?;
        let mut log_fugacity = vec![0.0; n];
        let mut hessian = vec![0.0; n * n];
        state.log_fugacity_tv_into(volume, &mut log_fugacity, &mut hessian)?;
        let active: Vec<_> = composition
            .iter()
            .enumerate()
            .filter_map(|(i, &z)| (z > 0.0).then_some(i))
            .collect();
        let (eigenvalue, vector) = scaled_eigenpair(composition, &active, &hessian)?;
        let mut eigenvector = vec![0.0; n];
        let mut direction = vec![0.0; n];
        let mut step: f64 = 1.0e-3;
        for (i, &a) in active.iter().enumerate() {
            eigenvector[a] = vector[i];
            direction[a] = composition[a].sqrt() * vector[i];
            if direction[a] != 0.0 {
                step = step.min(0.1 * composition[a] / direction[a].abs());
            }
        }
        if let Some(third_derivative) = state.critical_third_directional(volume, &direction) {
            if !third_derivative.is_finite() {
                return Err(Error::NotConverged("non-finite critical third derivative"));
            }
            return Ok(MixtureCriticalConditions {
                eigenvalue,
                eigenvector,
                third_derivative,
            });
        }
        if step < 1.0e-12 {
            return Err(Error::NotConverged(
                "critical direction below numerical resolution",
            ));
        }
        let mut trial = vec![0.0; n];
        let mut contraction = |s: f64| -> Result<f64> {
            for i in 0..n {
                trial[i] = composition[i] + s * direction[i];
            }
            let amount: f64 = trial.iter().sum();
            for z in &mut trial {
                *z /= amount;
            }
            let shifted = self.prepare(temperature, &trial)?;
            shifted.log_fugacity_tv_into(volume / amount, &mut log_fugacity, &mut hessian)?;
            let mut value = 0.0;
            for &i in &active {
                for &j in &active {
                    value += direction[i] * hessian[i * n + j] * direction[j];
                }
            }
            // A_nn scales as 1/N; preparation uses a normalized one-mole state.
            Ok(value / amount)
        };
        let third_derivative = (contraction(-2.0 * step)? - 8.0 * contraction(-step)?
            + 8.0 * contraction(step)?
            - contraction(2.0 * step)?)
            / (12.0 * step);
        if !third_derivative.is_finite() {
            return Err(Error::NotConverged("non-finite critical third derivative"));
        }
        Ok(MixtureCriticalConditions {
            eigenvalue,
            eigenvector,
            third_derivative,
        })
    }
}

/// Refine a local mixture critical point from explicit positive T and molar V.
///
/// The composition must be normalized. A bounded Newton iteration in ln(T),
/// ln(V) uses a numerically differentiated Jacobian of the eigenvalue and third
/// derivative. Backtracking rejects invalid volumes and non-decreasing residuals.
/// The supplied guess selects a local basin; this does not enumerate a critical
/// locus or test global phase stability. The pure limit uses the same equations.
pub fn critical_point_with_guess(
    model: &Cubic,
    composition: &[f64],
    temperature: f64,
    volume: f64,
) -> Result<CriticalPoint> {
    positive(temperature, "critical temperature guess")?;
    positive(volume, "critical volume guess")?;
    let mut coordinates = [temperature.ln(), volume.ln()];
    let residual = |q: [f64; 2]| -> Result<[f64; 2]> {
        let c = model.mixture_critical_conditions(q[0].exp(), q[1].exp(), composition)?;
        Ok([c.eigenvalue, c.third_derivative])
    };
    for _ in 0..80 {
        let values = residual(coordinates)?;
        let norm = values[0].abs().max(values[1].abs());
        if norm < 2.0e-9 {
            let temperature = coordinates[0].exp();
            let volume = coordinates[1].exp();
            let pressure = model
                .prepare(temperature, composition)?
                .pressure(volume)?
                .value;
            positive(pressure, "critical pressure")?;
            return Ok(CriticalPoint {
                temperature,
                volume,
                pressure,
            });
        }
        let step = 2.0e-4;
        let mut jacobian = [[0.0; 2]; 2];
        for j in 0..2 {
            let mut upper = coordinates;
            let mut lower = coordinates;
            upper[j] += step;
            lower[j] -= step;
            let up = residual(upper)?;
            let down = residual(lower)?;
            for i in 0..2 {
                jacobian[i][j] = (up[i] - down[i]) / (2.0 * step);
            }
        }
        let [[a, b], [c, d]] = jacobian;
        let determinant = a * d - b * c;
        if !determinant.is_finite() || determinant.abs() < 1.0e-18 {
            return Err(Error::NotConverged("singular mixture critical Jacobian"));
        }
        let mut delta = [
            (b * values[1] - d * values[0]) / determinant,
            (c * values[0] - a * values[1]) / determinant,
        ];
        let largest = delta[0].abs().max(delta[1].abs());
        if largest > 0.3 {
            for d in &mut delta {
                *d *= 0.3 / largest;
            }
        }
        let mut accepted = false;
        let mut fraction = 1.0;
        for _ in 0..24 {
            let trial = [
                coordinates[0] + fraction * delta[0],
                coordinates[1] + fraction * delta[1],
            ];
            if let Ok(r) = residual(trial)
                && r[0].abs().max(r[1].abs()) < norm
            {
                coordinates = trial;
                accepted = true;
                break;
            }
            fraction *= 0.5;
        }
        if !accepted {
            return Err(Error::NotConverged("mixture critical backtracking"));
        }
    }
    Err(Error::NotConverged("mixture critical point"))
}

pub(crate) fn mixture_critical_point(model: &Cubic, composition: &[f64]) -> Result<CriticalPoint> {
    let guess = model.pseudo_critical(composition)?;
    let state = model.prepare(guess.temperature, composition)?;
    // The fork's default starts at four covolumes. Additional seeds are only
    // attempted after a failed local solve, without suppressing its equations.
    for (t, v) in [(1.0, 4.0), (0.9, 4.0), (1.1, 4.0), (1.0, 3.0), (1.0, 6.0)] {
        let volume = v * state.eos_covolume() - state.volume_translation().value;
        match critical_point_with_guess(model, composition, t * guess.temperature, volume) {
            Ok(point) => return Ok(point),
            Err(error @ Error::Unsupported(_)) => return Err(error),
            Err(_) => {}
        }
    }
    Err(Error::NotConverged("mixture critical point"))
}

// Symmetric Jacobi rotations. No foreign library or mutable global workspace.
// The active matrix excludes zero components instead of adding an arbitrary
// large diagonal as the source DSYEV caller does.
fn smallest_eigenpair(matrix: &mut [f64], n: usize) -> Result<(f64, Vec<f64>)> {
    if matrix.iter().any(|x| !x.is_finite()) {
        return Err(Error::InvalidInput("non-finite stability matrix"));
    }
    let mut vectors: Vec<f64> = vec![0.0; n * n];
    for i in 0..n {
        vectors[i * n + i] = 1.0;
    }
    for _ in 0..(100 * n * n) {
        let mut largest: f64 = 0.0;
        let (mut p, mut q) = (0, 0);
        for i in 0..n {
            for j in i + 1..n {
                if matrix[i * n + j].abs() > largest {
                    largest = matrix[i * n + j].abs();
                    p = i;
                    q = j;
                }
            }
        }
        let scale = (0..n).fold(1.0_f64, |v, i| v.max(matrix[i * n + i].abs()));
        if largest <= 2.0e-15 * scale {
            let index = (0..n)
                .min_by(|&i, &j| matrix[i * n + i].total_cmp(&matrix[j * n + j]))
                .ok_or(Error::InvalidInput("empty active composition"))?;
            let mut vector: Vec<_> = (0..n).map(|i| vectors[i * n + index]).collect();
            let pivot = vector
                .iter()
                .copied()
                .max_by(|a, b| a.abs().total_cmp(&b.abs()))
                .unwrap();
            if pivot < 0.0 {
                for v in &mut vector {
                    *v = -*v;
                }
            }
            return Ok((matrix[index * n + index], vector));
        }
        let app = matrix[p * n + p];
        let aqq = matrix[q * n + q];
        let apq = matrix[p * n + q];
        let tau = (aqq - app) / (2.0 * apq);
        let tangent = if tau >= 0.0 { 1.0 } else { -1.0 } / (tau.abs() + tau.hypot(1.0));
        let cosine = 1.0 / (1.0 + tangent * tangent).sqrt();
        let sine = tangent * cosine;
        matrix[p * n + p] = app - tangent * apq;
        matrix[q * n + q] = aqq + tangent * apq;
        matrix[p * n + q] = 0.0;
        matrix[q * n + p] = 0.0;
        for k in 0..n {
            if k != p && k != q {
                let kp = matrix[k * n + p];
                let kq = matrix[k * n + q];
                matrix[k * n + p] = cosine * kp - sine * kq;
                matrix[p * n + k] = matrix[k * n + p];
                matrix[k * n + q] = sine * kp + cosine * kq;
                matrix[q * n + k] = matrix[k * n + q];
            }
            let kp = vectors[k * n + p];
            let kq = vectors[k * n + q];
            vectors[k * n + p] = cosine * kp - sine * kq;
            vectors[k * n + q] = sine * kp + cosine * kq;
        }
    }
    Err(Error::NotConverged("symmetric stability eigenpair"))
}

#[cfg(test)]
mod tests {
    use super::smallest_eigenpair;
    #[test]
    fn symmetric_eigenpairs_include_degeneracy_and_large_active_sets() {
        for n in [1, 2, 5, 20] {
            // Rank-one perturbation of a diagonal identity: its one minimum
            // eigenpair is known without a second eigenvalue implementation.
            let v: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
            let norm = v.iter().map(|x| x * x).sum::<f64>();
            let mut a = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..n {
                    a[i * n + j] = -2.0 * v[i] * v[j] / norm;
                    if i == j {
                        a[i * n + j] += 1.0;
                    }
                }
            }
            let (lambda, u) = smallest_eigenpair(&mut a, n).unwrap();
            assert!((lambda + 1.0).abs() < 1e-13);
            for i in 0..n {
                assert!((u[i] - v[i] / norm.sqrt()).abs() < 1e-13);
            }
        }
        let (lambda, u) = smallest_eigenpair(&mut [2.0, 1.0, 1.0, 2.0], 2).unwrap();
        assert!((lambda - 1.0).abs() < 1e-14);
        assert!((u[0] + u[1]).abs() < 1e-14);
        let (lambda, _) = smallest_eigenpair(&mut [1.0, 0.0, 0.0, 1.0], 2).unwrap();
        assert_eq!(lambda, 1.0);
    }
}
