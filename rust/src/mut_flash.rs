//! mu-T flash: the component densities at a fixed temperature and chemical
//! potential.
//! Ported from `src/mut_solver.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! The original solves the logarithmic fugacity equation
//! `lnf(T, rho) = lnf_spec` with a damped Newton iteration and converts a
//! chemical potential specification into that equation at the initial guess.

use crate::R;
use crate::cubic::Cubic;
use crate::error::{Error, Result};

/// The result of a mu-T or lnf-T solve.
#[derive(Debug, Clone)]
pub struct MuTSolution {
    /// The solved component densities [mol/m3].
    pub density: Vec<f64>,
    /// The log fugacity coefficients at the solution.
    pub ln_fugacity: Vec<f64>,
    /// The number of Newton iterations.
    pub iterations: usize,
}

/// Maximum number of Newton iterations.
const MAX_ITERATIONS: usize = 200;
/// Relative convergence threshold of the fugacity residual.
const RELATIVE_TOLERANCE: f64 = 1.0e-20;
/// Absolute convergence threshold of the fugacity residual.
const ABSOLUTE_TOLERANCE: f64 = 1.0e-12;
/// Relative step of the finite-difference Jacobian.
const JACOBIAN_STEP: f64 = 1.0e-7;

/// The log fugacity coefficients and chemical potentials at a density.
pub fn log_fugacity_and_potential(
    model: &Cubic,
    temperature: f64,
    density: &[f64],
) -> Result<(Vec<f64>, Vec<f64>)> {
    let n = density.len();
    let total: f64 = density.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(Error::InvalidInput("the density guess must be positive"));
    }
    let composition: Vec<f64> = density.iter().map(|value| value / total).collect();
    if composition.iter().any(|value| *value <= 0.0) {
        return Err(Error::InvalidInput(
            "the density guess must have positive components",
        ));
    }
    let state = model.prepare(temperature, &composition)?;
    // The original works with a total volume of one cubic meter, so the molar
    // volume follows from the total number of moles.
    let molar_volume = 1.0 / total;
    let mut coefficient = vec![0.0; n];
    state.ln_fugacity_into(molar_volume, &mut coefficient)?;
    // The original reports the logarithm of the fugacity, not the fugacity
    // coefficient, and the chemical potential follows from the coefficient and
    // the ideal part.
    let pressure = state.pressure(molar_volume)?.value;
    let ln_pressure = pressure.ln();
    let mut ln_fugacity = vec![0.0; n];
    let mut chemical_potential = vec![0.0; n];
    let rt = R * temperature;
    for index in 0..n {
        ln_fugacity[index] = coefficient[index] + composition[index].ln() + ln_pressure;
        chemical_potential[index] = rt * (coefficient[index] - composition[index].ln() - 1.0);
    }
    Ok((ln_fugacity, chemical_potential))
}

/// Solve `lnf(T, rho) = ln_fugacity` for the component densities.
///
/// The Newton iteration uses a finite-difference Jacobian of the analytic log
/// fugacity coefficients, while the original uses the second composition
/// derivatives of the Helmholtz energy. Both converge to the same root of the
/// same residual.
pub fn solve_lnf_t(
    model: &Cubic,
    temperature: f64,
    ln_fugacity: &[f64],
    initial_density: &[f64],
) -> Result<MuTSolution> {
    let n = initial_density.len();
    if ln_fugacity.len() != n || n == 0 {
        return Err(Error::InvalidInput(
            "component count mismatch in the lnf solve",
        ));
    }
    if ln_fugacity.iter().any(|value| !value.is_finite()) {
        return Err(Error::InvalidInput(
            "the log fugacity target must be finite",
        ));
    }
    let mut density = initial_density.to_vec();
    if density
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(Error::InvalidInput("the density guess must be positive"));
    }
    let upper = 2.0 * density.iter().cloned().fold(0.0_f64, f64::max);

    let mut current = log_fugacity_and_potential(model, temperature, &density)?.0;
    let mut jacobian = vec![0.0; n * n];
    for iteration in 1..=MAX_ITERATIONS {
        let residual: Vec<f64> = current
            .iter()
            .zip(ln_fugacity)
            .map(|(value, target)| value - target)
            .collect();
        if residual
            .iter()
            .all(|value| value.abs() < ABSOLUTE_TOLERANCE + RELATIVE_TOLERANCE * value.abs())
        {
            return Ok(MuTSolution {
                density,
                ln_fugacity: current,
                iterations: iteration - 1,
            });
        }
        for index in 0..n {
            let mut step = JACOBIAN_STEP * density[index].abs().max(1.0);
            let mut shifted = None;
            for _ in 0..8 {
                let mut perturbed = density.clone();
                perturbed[index] += step;
                if let Ok((fugacity, _)) =
                    log_fugacity_and_potential(model, temperature, &perturbed)
                {
                    shifted = Some(fugacity);
                    break;
                }
                step *= 0.5;
            }
            let shifted = match shifted {
                Some(shifted) => shifted,
                None => {
                    // The column cannot be evaluated at a usable state, so the
                    // derivative stays at zero and the iteration may still
                    // converge.
                    for row in 0..n {
                        jacobian[row * n + index] = 0.0;
                    }
                    continue;
                }
            };
            for row in 0..n {
                jacobian[row * n + index] = (shifted[row] - current[row]) / step;
            }
        }
        let mut rhs = residual.clone();
        crate::numerics::solve_linear_in_place(&mut jacobian, &mut rhs)?;
        // The Newton step is the negative solution of J x = residual.
        let mut step: Vec<f64> = rhs.iter().map(|value| -value).collect();
        // Limit the step so that the densities stay positive and within the
        // original bounds.
        let largest = step
            .iter()
            .enumerate()
            .map(|(index, value)| value.abs() / density[index])
            .fold(0.0_f64, f64::max);
        if largest > 0.5 {
            let scale = 0.5 / largest;
            for value in step.iter_mut() {
                *value *= scale;
            }
        }
        // A step into a state without a positive pressure is not usable. The
        // original evaluates the formal log fugacity there, the native model
        // reports an error, so the step is halved until the state is usable.
        let previous = density.clone();
        let mut scale = 1.0;
        let mut attempts = 0;
        loop {
            for (index, value) in step.iter().enumerate() {
                let updated = previous[index] + scale * value;
                density[index] = updated.clamp(1.0e-12, upper);
            }
            match log_fugacity_and_potential(model, temperature, &density) {
                Ok((fugacity, _)) => {
                    current = fugacity;
                    break;
                }
                Err(_) if attempts < 12 => {
                    attempts += 1;
                    scale *= 0.5;
                }
                Err(error) => return Err(error),
            }
        }
    }
    Err(Error::NotConverged("the lnf-T iteration did not converge"))
}

/// Solve `mu(T, rho) = chemical_potential` for the component densities.
///
/// The original converts the chemical potential specification into a log
/// fugacity specification at the initial guess and then solves the log
/// fugacity equation.
pub fn solve_mu_t(
    model: &Cubic,
    temperature: f64,
    chemical_potential: &[f64],
    initial_density: &[f64],
) -> Result<MuTSolution> {
    let n = initial_density.len();
    if chemical_potential.len() != n || n == 0 {
        return Err(Error::InvalidInput(
            "component count mismatch in the mu solve",
        ));
    }
    if chemical_potential.iter().any(|value| !value.is_finite()) {
        return Err(Error::InvalidInput("the chemical potential must be finite"));
    }
    let (ln_fugacity, current_potential) =
        log_fugacity_and_potential(model, temperature, initial_density)?;
    let rt = R * temperature;
    let target: Vec<f64> = chemical_potential
        .iter()
        .zip(&current_potential)
        .zip(&ln_fugacity)
        .map(|((mu, current), fugacity)| (mu - current) / rt + fugacity)
        .collect();
    solve_lnf_t(model, temperature, &target, initial_density)
}
