//! Weeks-Chandler-Andersen (WCA) perturbation model for the splined
//! Lennard-Jones fluid (`LJS-WCA`).
//! Ported from `src/lj_splined.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.

// The coefficients of the source carry the precision of the source.
#![allow(clippy::excessive_precision)]

use crate::Phase;
use crate::chain_rule::PackingFractionDerivatives;
use crate::error::{Error, Result, positive};
use crate::hardsphere_wca::{
    CavityPreparedTemp, hard_sphere_diameter, packing_fraction_vdw_no_segments, pure_eta_per_rho,
    pure_hard_sphere_density, pure_hard_sphere_eta,
};
use crate::helmholtz::{DensityHelmholtz, PureHelmholtz, PurePressure};
use crate::ljs::{
    ALPHA_LJS, ARGON_EPS, ARGON_SIGMA, LjsState, LjsTpState, LjsTvState, LjsVolumeIteration,
    baseline_saft_volume_iteration, estimate_ljs_critical_parameters, evaluate_tp_state,
    evaluate_tv_state, ljs_pair_potential, ljs_potential_reduced, safe_tanh_derivs,
    select_physical_volume, solve_physical_volumes,
};

/// The reduced WCA potential split separation, `2^(1/6)`.
const R_SPLIT: f64 = 1.122_462_048_309_373;

/// The coefficients of the first order term of the source.
const A1_COEFFICIENTS: [f64; 12] = [
    -4.359342008565890e-01,
    -1.010684675331760,
    5.545568513310810e-01,
    -4.948459993882990e-01,
    -2.175255239422550e-02,
    2.784218910493600,
    -9.960315909054440e-01,
    7.619772900616760,
    -1.399502159877310e+01,
    8.957336641810310e-01,
    -5.449504201398070,
    8.809619552329060,
];

/// The coefficients of the second order term of the source.
const A2_COEFFICIENTS: [f64; 14] = [
    -8.549712604490140e-01,
    -3.309657322845890,
    2.150262062951070,
    1.819290452429280,
    4.881241094741120,
    -1.266600852072320e+01,
    2.720757103054650,
    -2.012625560669800e+01,
    3.826002401416580e+01,
    2.528308672696750,
    4.436429616364860,
    -2.714948088128380e+01,
    4.413607271085820e-02,
    9.921300491978420e-01,
];

/// The coefficients of the third order term of the source.
const A3_COEFFICIENTS: [f64; 13] = [
    -2.299701237101790e-01,
    -9.378822827801790e-01,
    2.598656446923390e-01,
    -2.508338793969640,
    4.339750952753030,
    -1.048143274408200e+01,
    4.115203492416100e+01,
    -3.548798441767570e+01,
    6.040169542924050,
    -9.971355447734860,
    1.004780979021150,
    7.443292068847610e-01,
    3.431138067711600e-01,
];

/// The coefficients of the fourth order term of the source.
const A4_COEFFICIENTS: [f64; 14] = [
    -3.887755601943300e-03,
    7.088380007122160e-01,
    -5.061184486003090,
    -5.601386794807250e-02,
    -1.855154519243070e-01,
    2.148962864405380e-01,
    6.920413400878820e-02,
    -6.212126555148940e-01,
    5.101710025798300,
    6.495708980657240,
    -2.494763026982940e+01,
    3.710757600727180e+01,
    3.683587731520040e-01,
    1.903500562425810,
];

/// The values of the reduced distance and their first two derivatives.
fn distance_basis(y: f64) -> [[f64; 3]; 3] {
    let mut basis = [[0.0; 3]; 3];
    basis[0][0] = 1.0;
    basis[1][0] = y;
    basis[2][0] = y * y;
    basis[1][1] = 1.0;
    basis[2][1] = 2.0 * y;
    basis[2][2] = 2.0;
    basis
}

/// The values of the reduced density and their first two derivatives.
fn density_basis(rho: f64) -> [[f64; 3]; 4] {
    let mut basis = [[0.0; 3]; 4];
    basis[0][0] = rho;
    basis[0][1] = 1.0;
    basis[1][2] = 2.0;
    for index in 1..4 {
        let order = (index + 1) as f64;
        basis[index][0] = basis[index - 1][0] * rho;
        basis[index][1] = order * basis[index - 1][0];
        if index > 1 {
            basis[index][2] = order * (order - 1.0) * basis[index - 2][0];
        }
    }
    basis
}

/// The first order term of the WCA perturbation family.
fn first_order(y: f64, rho: f64) -> (f64, f64, f64, f64, f64, f64) {
    let yv = distance_basis(y);
    let rhov = density_basis(rho);
    let mut p = [[0.0; 3]; 4];
    for order in 0..3 {
        p[0][order] = A1_COEFFICIENTS[0] * yv[0][order]
            + A1_COEFFICIENTS[1] * yv[1][order]
            + A1_COEFFICIENTS[2] * yv[2][order];
        p[1][order] = A1_COEFFICIENTS[3] * yv[0][order]
            + A1_COEFFICIENTS[4] * yv[1][order]
            + A1_COEFFICIENTS[5] * yv[2][order];
        p[2][order] = A1_COEFFICIENTS[6] * yv[0][order]
            + A1_COEFFICIENTS[7] * yv[1][order]
            + A1_COEFFICIENTS[8] * yv[2][order];
        p[3][order] = A1_COEFFICIENTS[9] * yv[0][order]
            + A1_COEFFICIENTS[10] * yv[1][order]
            + A1_COEFFICIENTS[11] * yv[2][order];
    }
    let mut a = 0.0;
    let mut a_y = 0.0;
    let mut a_r = 0.0;
    let mut a_yy = 0.0;
    let mut a_rr = 0.0;
    let mut a_ry = 0.0;
    for index in 0..4 {
        a += p[index][0] * rhov[index][0];
        a_y += p[index][1] * rhov[index][0];
        a_r += p[index][0] * rhov[index][1];
        a_yy += p[index][2] * rhov[index][0];
        a_rr += p[index][0] * rhov[index][2];
        a_ry += p[index][1] * rhov[index][1];
    }
    (a, a_r, a_y, a_rr, a_yy, a_ry)
}

/// The second order term of the WCA perturbation family.
fn second_order(y: f64, rho: f64) -> (f64, f64, f64, f64, f64, f64) {
    let c = &A2_COEFFICIENTS;
    let yv = distance_basis(y);
    let mut p = [[0.0; 3]; 4];
    for order in 0..3 {
        p[0][order] = c[0] * yv[0][order] + c[1] * yv[1][order] + c[2] * yv[2][order];
        p[1][order] = c[3] * yv[0][order] + c[4] * yv[1][order] + c[5] * yv[2][order];
        p[2][order] = c[6] * yv[0][order] + c[7] * yv[1][order] + c[8] * yv[2][order];
        p[3][order] = c[9] * yv[0][order] + c[10] * yv[1][order] + c[11] * yv[2][order];
    }
    let log_rho = rho.ln();
    let mut rho_p4 = [0.0; 4];
    rho_p4[0] = rho.powf(p[3][0]);
    rho_p4[1] = rho_p4[0] * p[3][1] * log_rho;
    rho_p4[2] = rho_p4[0] * p[3][2] * log_rho + rho_p4[0] * p[3][1] * p[3][1] * log_rho * log_rho;
    rho_p4[3] = p[3][1] * rho.powf(p[3][0] - 1.0) * (p[3][0] * log_rho + 1.0);
    let mut u = [[0.0; 3]; 3];
    u[0][0] = p[1][0] * (rho.powf(c[13]) + p[2][0] * rho_p4[0]);
    u[0][1] =
        p[1][0] * (c[13] * rho.powf(c[13] - 1.0) + p[2][0] * p[3][0] * rho.powf(p[3][0] - 1.0));
    u[0][2] = p[1][0]
        * (c[13] * (c[13] - 1.0) * rho.powf(c[13] - 2.0)
            + p[2][0] * p[3][0] * (p[3][0] - 1.0) * rho.powf(p[3][0] - 2.0));
    u[1][0] = p[1][1] * (rho.powf(c[13]) + p[2][0] * rho_p4[0])
        + p[1][0] * (p[2][1] * rho_p4[0] + p[2][0] * rho_p4[1]);
    u[2][0] = p[1][2] * (rho.powf(c[13]) + p[2][0] * rho_p4[0])
        + 2.0 * p[1][1] * (p[2][1] * rho_p4[0] + p[2][0] * rho_p4[1])
        + p[1][0] * (p[2][2] * rho_p4[0] + 2.0 * p[2][1] * rho_p4[1] + p[2][0] * rho_p4[2]);
    u[1][1] = p[1][1]
        * (c[13] * rho.powf(c[13] - 1.0) + p[2][0] * p[3][0] * rho.powf(p[3][0] - 1.0))
        + p[1][0] * (p[2][1] * p[3][0] * rho.powf(p[3][0] - 1.0) + p[2][0] * rho_p4[3]);
    let (th0, th1, th2) = safe_tanh_derivs(u[0][0]);
    let tanhu = [th0, th1, th2];
    let mut g = [[0.0; 3]; 3];
    g[0][0] = -rho * (p[0][0] + c[12]);
    g[0][1] = -p[0][0] - c[12];
    g[1][0] = -rho * p[0][1];
    g[2][0] = -rho * p[0][2];
    g[1][1] = -p[0][1];
    let a = p[0][0] * rho + g[0][0] * tanhu[0];
    let a_y = p[0][1] * rho + g[1][0] * tanhu[0] + g[0][0] * tanhu[1] * u[1][0];
    let a_r = p[0][0] + g[0][1] * tanhu[0] + g[0][0] * tanhu[1] * u[0][1];
    let a_yy = p[0][2] * rho
        + g[2][0] * tanhu[0]
        + 2.0 * g[1][0] * tanhu[1] * u[1][0]
        + g[0][0] * tanhu[1] * u[2][0]
        + g[0][0] * tanhu[2] * u[1][0] * u[1][0];
    let a_rr = 2.0 * g[0][1] * tanhu[1] * u[0][1]
        + g[0][0] * tanhu[1] * u[0][2]
        + g[0][0] * tanhu[2] * u[0][1] * u[0][1];
    let a_ry = p[0][1]
        + g[1][1] * tanhu[0]
        + g[1][0] * tanhu[1] * u[0][1]
        + g[0][1] * tanhu[1] * u[1][0]
        + g[0][0] * tanhu[2] * u[1][0] * u[0][1]
        + g[0][0] * tanhu[1] * u[1][1];
    (a, a_r, a_y, a_rr, a_yy, a_ry)
}

/// The third order term of the WCA perturbation family.
fn third_order(y: f64, rho: f64) -> (f64, f64, f64, f64, f64, f64) {
    let c = &A3_COEFFICIENTS;
    let yv = distance_basis(y);
    let mut p = [[0.0; 3]; 5];
    for order in 0..3 {
        p[0][order] = c[0] * yv[0][order] + c[1] * yv[1][order] + c[2] * yv[2][order];
        p[2][order] = c[3] * yv[0][order] + c[4] * yv[1][order];
        p[1][order] = c[5] * yv[0][order] + c[6] * yv[1][order] + c[7] * yv[2][order];
        p[3][order] = c[8] * yv[0][order] + c[9] * yv[1][order];
        p[4][order] = c[11] * yv[0][order] + c[12] * yv[1][order];
    }
    let log_rho = rho.ln();
    let mut rho_p = [[0.0; 4]; 2];
    for (index, entry) in rho_p.iter_mut().enumerate() {
        let source = index + 3;
        let base = rho.powf(p[source][0]);
        entry[0] = base;
        entry[1] = base * p[source][1] * log_rho;
        entry[2] =
            base * p[source][2] * log_rho + base * p[source][1] * p[source][1] * log_rho * log_rho;
        entry[3] = p[source][1] * rho.powf(p[source][0] - 1.0) * (p[source][0] * log_rho + 1.0);
    }
    let mut r = [[0.0; 3]; 3];
    for (index, entry) in rho_p.iter().enumerate() {
        let source = index + 3;
        let partner = source - 2;
        r[0][0] += p[partner][0] * entry[0];
        r[0][1] += p[partner][0] * p[source][0] * rho.powf(p[source][0] - 1.0);
        r[0][2] +=
            p[partner][0] * p[source][0] * (p[source][0] - 1.0) * rho.powf(p[source][0] - 2.0);
        r[1][0] += p[partner][1] * entry[0] + p[partner][0] * entry[1];
        r[2][0] +=
            2.0 * p[partner][1] * entry[1] + p[partner][2] * entry[0] + p[partner][0] * entry[2];
        r[1][1] += p[partner][1] * p[source][0] * rho.powf(p[source][0] - 1.0)
            + p[partner][0] * p[source][1] * rho.powf(p[source][0] - 1.0)
            + p[partner][0] * p[source][0] * entry[1] / rho;
    }
    let (th0, th1, th2) = safe_tanh_derivs(r[0][0]);
    let tanhr = [1.0 + c[10] * th0, c[10] * th1, c[10] * th2];
    let mut g = [[0.0; 3]; 3];
    g[0][0] = rho * p[0][0];
    g[0][1] = p[0][0];
    g[1][0] = rho * p[0][1];
    g[2][0] = rho * p[0][2];
    g[1][1] = p[0][1];
    let a = g[0][0] * tanhr[0];
    let a_y = g[1][0] * tanhr[0] + g[0][0] * tanhr[1] * r[1][0];
    let a_r = g[0][1] * tanhr[0] + g[0][0] * tanhr[1] * r[0][1];
    let a_yy = g[2][0] * tanhr[0]
        + 2.0 * g[1][0] * tanhr[1] * r[1][0]
        + g[0][0] * tanhr[1] * r[2][0]
        + g[0][0] * tanhr[2] * r[1][0] * r[1][0];
    let a_rr = 2.0 * g[0][1] * tanhr[1] * r[0][1]
        + g[0][0] * tanhr[1] * r[0][2]
        + g[0][0] * tanhr[2] * r[0][1] * r[0][1];
    let a_ry = g[1][1] * tanhr[0]
        + g[1][0] * tanhr[1] * r[0][1]
        + g[0][1] * tanhr[1] * r[1][0]
        + g[0][0] * tanhr[2] * r[1][0] * r[0][1]
        + g[0][0] * tanhr[1] * r[1][1];
    (a, a_r, a_y, a_rr, a_yy, a_ry)
}

/// The fourth order term of the WCA perturbation family.
fn fourth_order(y: f64, rho: f64) -> (f64, f64, f64, f64, f64, f64) {
    let c = &A4_COEFFICIENTS;
    let yv = distance_basis(y);
    let mut p = [[0.0; 3]; 4];
    for order in 0..3 {
        p[0][order] = c[0] * yv[0][order] + c[1] * yv[1][order] + c[2] * yv[2][order];
        p[1][order] = c[3] * yv[0][order] + c[4] * yv[1][order] + c[5] * yv[2][order];
        p[2][order] = c[6] * yv[0][order] + c[7] * yv[1][order] + c[8] * yv[2][order];
        p[3][order] = c[9] * yv[0][order] + c[10] * yv[1][order] + c[11] * yv[2][order];
    }
    let mut k = [[0.0; 3]; 3];
    k[0][0] = p[3][0] * rho.powf(c[13]);
    k[0][1] = c[13] * p[3][0] * rho.powf(c[13] - 1.0);
    k[0][2] = (c[13] - 1.0) * c[13] * p[3][0] * rho.powf(c[13] - 2.0);
    k[1][0] = p[3][1] * rho.powf(c[13]);
    k[2][0] = p[3][2] * rho.powf(c[13]);
    k[1][1] = c[13] * p[3][1] * rho.powf(c[13] - 1.0);
    let (th0, th1, th2) = safe_tanh_derivs(k[0][0]);
    let tanhk = [th0, th1, th2];
    let mut g = [[0.0; 2]; 3];
    g[0][0] = rho * p[2][0];
    g[0][1] = p[2][0];
    g[1][0] = rho * p[2][1];
    g[1][1] = p[2][1];
    g[2][0] = rho * p[2][2];
    let mut a = p[0][0] * rho.powf(c[12] + 1.0) + p[1][0] * rho;
    let mut a_y = p[0][1] * rho.powf(c[12] + 1.0) + p[1][1] * rho;
    let mut a_yy = p[0][2] * rho.powf(c[12] + 1.0) + p[1][2] * rho;
    let mut a_r = (c[12] + 1.0) * p[0][0] * rho.powf(c[12]) + p[1][0];
    let mut a_rr = c[12] * (c[12] + 1.0) * p[0][0] * rho.powf(c[12] - 1.0);
    let mut a_ry = (c[12] + 1.0) * p[0][1] * rho.powf(c[12]) + p[1][1];
    a += g[0][0] * tanhk[0];
    a_y += g[1][0] * tanhk[0] + g[0][0] * tanhk[1] * k[1][0];
    a_r += g[0][1] * tanhk[0] + g[0][0] * tanhk[1] * k[0][1];
    a_yy += g[2][0] * tanhk[0]
        + 2.0 * g[1][0] * tanhk[1] * k[1][0]
        + g[0][0] * tanhk[2] * k[1][0] * k[1][0]
        + g[0][0] * tanhk[1] * k[2][0];
    a_rr += 2.0 * g[0][1] * tanhk[1] * k[0][1]
        + g[0][0] * tanhk[1] * k[0][2]
        + g[0][0] * tanhk[2] * k[0][1] * k[0][1];
    a_ry += g[1][1] * tanhk[0]
        + g[1][0] * tanhk[1] * k[0][1]
        + g[0][1] * tanhk[1] * k[1][0]
        + g[0][0] * tanhk[2] * k[1][0] * k[0][1]
        + g[0][0] * tanhk[1] * k[1][1];
    (a, a_r, a_y, a_rr, a_yy, a_ry)
}

/// One perturbation term of the WCA family and its derivatives.
///
/// The index selects the term and follows the source, the first order term is
/// one and the fourth order term is four. The reduced temperature and the
/// reduced density are dimensionless. The method returns the term, its
/// derivative with respect to the reduced density, its derivative with respect
/// to the reduced distance and the second derivatives, in the order of the
/// source (`calc_LJs_WCA_ai_tr`).
pub fn perturbation_term(
    index: usize,
    reduced_temperature: f64,
    reduced_density: f64,
) -> Result<(f64, f64, f64, f64, f64, f64)> {
    if reduced_temperature <= 0.0 || !reduced_temperature.is_finite() {
        return Err(Error::InvalidInput(
            "the reduced temperature must be positive",
        ));
    }
    if reduced_density <= 0.0 || !reduced_density.is_finite() {
        return Err(Error::InvalidInput("the reduced density must be positive"));
    }
    let x0 = (2.0 / (1.0 + reduced_temperature.sqrt())).powf(1.0 / 6.0);
    let y = R_SPLIT - x0;
    let (value, a_r, a_y, a_rr, a_yy, a_ry) = match index {
        1 => {
            let (value, a_r, a_y, a_rr, a_yy, a_ry) = first_order(y, reduced_density);
            let prefactor = 2.0 * (1.0f64).atan();
            (
                value * prefactor,
                a_r * prefactor,
                a_y * prefactor,
                a_rr * prefactor,
                a_yy * prefactor,
                a_ry * prefactor,
            )
        }
        2 => second_order(y, reduced_density),
        3 => third_order(y, reduced_density),
        4 => fourth_order(y, reduced_density),
        _ => {
            return Err(Error::InvalidInput(
                "the perturbation term index must be one to four",
            ));
        }
    };
    Ok((value, a_r, a_y, a_rr, a_yy, a_ry))
}

/// One perturbation term of the WCA family with the prefactor of the source
/// and its derivatives (`calc_LJs_WCA_ai`).
///
/// The index selects the term (`1..=4`). `eps_divk` is in kelvin,
/// `reduced_distance` is `y = 2^(1/6) - d/sigma`, and `reduced_density` is
/// `rho* = N_A sigma^3 n/V`.
pub fn wca_term(
    index: usize,
    eps_divk: f64,
    reduced_distance: f64,
    reduced_density: f64,
) -> Result<(f64, f64, f64, f64, f64, f64)> {
    if reduced_density <= 0.0 || !reduced_density.is_finite() {
        return Err(Error::InvalidInput("the reduced density must be positive"));
    }
    let (value, value_r, value_y, value_rr, value_yy, value_ry) = match index {
        1 => first_order(reduced_distance, reduced_density),
        2 => second_order(reduced_distance, reduced_density),
        3 => third_order(reduced_distance, reduced_density),
        4 => fourth_order(reduced_distance, reduced_density),
        _ => {
            return Err(Error::InvalidInput(
                "the perturbation term index must be one to four",
            ));
        }
    };
    let prefactor = match index {
        1 => 2.0 * std::f64::consts::PI * eps_divk,
        2 => eps_divk * eps_divk,
        3 => eps_divk * eps_divk * eps_divk,
        _ => eps_divk * eps_divk * eps_divk * eps_divk,
    };
    Ok((
        prefactor * value,
        prefactor * value_r,
        prefactor * value_y,
        prefactor * value_rr,
        prefactor * value_yy,
        prefactor * value_ry,
    ))
}

/// Molar WCA perturbation term `a_i(T, rho)` (before dividing by `T^i`) and
/// its `(T, rho)` derivatives.
#[allow(clippy::too_many_arguments)]
pub fn wca_term_density(
    index: usize,
    sigma: f64,
    eps_divk: f64,
    diameter: f64,
    diameter_t: f64,
    diameter_tt: f64,
    rho: f64,
) -> Result<DensityHelmholtz> {
    positive(sigma, "segment diameter")?;
    positive(eps_divk, "well depth")?;
    positive(diameter, "hard-sphere diameter")?;
    positive(rho, "molar density")?;
    let scale = crate::N_A * sigma * sigma * sigma;
    let rho_star = scale * rho;
    let y = R_SPLIT - diameter / sigma;
    let y_t = -diameter_t / sigma;
    let y_tt = -diameter_tt / sigma;
    let (a, a_r, a_y, a_rr, a_yy, a_ry) = wca_term(index, eps_divk, y, rho_star)?;
    Ok(DensityHelmholtz {
        value: a,
        dt: a_y * y_t,
        drho: a_r * scale,
        dtt: a_yy * y_t * y_t + a_y * y_tt,
        dtrho: a_ry * y_t * scale,
        drhorho: a_rr * scale * scale,
    })
}

/// Extensive WCA perturbation term `n * a_i(T, V, n)` (before dividing by `T^i`)
/// and all 10 `(T, V, n)` derivatives.
///
/// Ported from `src/lj_splined.f90:calc_LJs_WCA_ai_TVN`.
#[allow(clippy::too_many_arguments)]
pub fn wca_term_tvn(
    index: usize,
    sigma: f64,
    eps_divk: f64,
    diameter: f64,
    diameter_t: f64,
    diameter_tt: f64,
    volume: f64,
    moles: f64,
) -> Result<PureHelmholtz> {
    positive(volume, "total volume")?;
    positive(moles, "mole number")?;
    let dh = wca_term_density(
        index,
        sigma,
        eps_divk,
        diameter,
        diameter_t,
        diameter_tt,
        moles / volume,
    )?;
    dh.extensive(volume, moles)
}

/// Term and hard-sphere reference controls for the `LJS-WCA` model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LjsWcaControl {
    /// Include the WCA cavity soft-repulsion integral (requires `enable_hard_sphere`).
    pub enable_cavity: bool,
    /// Include the hard-sphere reference (and cavity term if `enable_cavity` is true).
    pub enable_hard_sphere: bool,
    /// Include the first-order perturbation term `a1 / T`.
    pub enable_a1: bool,
    /// Include the second-order perturbation term `a2 / T^2`.
    pub enable_a2: bool,
    /// Include the third-order perturbation term `a3 / T^3`.
    pub enable_a3: bool,
    /// Include the fourth-order perturbation term `a4 / T^4`.
    pub enable_a4: bool,
    /// Use the Carnahan-Starling-Kolafa pure hard-sphere reference (`PURE_HS_CSK`)
    /// instead of Carnahan-Starling (`PURE_HS_CS`).
    pub kolafa_pure_reference: bool,
}

impl Default for LjsWcaControl {
    fn default() -> Self {
        Self {
            enable_cavity: true,
            enable_hard_sphere: true,
            enable_a1: true,
            enable_a2: true,
            enable_a3: true,
            enable_a4: true,
            kolafa_pure_reference: false,
        }
    }
}

/// The splined Lennard-Jones fluid model using Weeks-Chandler-Andersen (WCA)
/// perturbation theory (`LJS-WCA`).
///
/// Ported from `src/lj_splined.f90:ljs_wca_eos` and `calcFresLJs_WCA`.
#[derive(Debug, Clone, Copy)]
pub struct LjsWca {
    sigma: f64,
    eps: f64,
    control: LjsWcaControl,
}

impl Default for LjsWca {
    fn default() -> Self {
        Self::new()
    }
}

/// Temperature-prepared `LJS-WCA` state caching the WCA hard-sphere diameter,
/// packing-fraction factor, reduced distance, and cavity quadrature nodes.
#[derive(Debug, Clone)]
pub struct PreparedLjsWca {
    sigma: f64,
    eps: f64,
    control: LjsWcaControl,
    temperature: f64,
    diameter: f64,
    diameter_t: f64,
    diameter_tt: f64,
    eta_per_rho: f64,
    eta_per_rho_t: f64,
    eta_per_rho_tt: f64,
    rho_star_per_rho: f64,
    y: f64,
    y_t: f64,
    y_tt: f64,
    cavity: Option<CavityPreparedTemp>,
}

impl LjsWca {
    /// Create the default `LJS-WCA` model (Argon record: `sigma = 3.42e-10 m`, `eps/k = 124.0 K`).
    pub fn new() -> Self {
        Self {
            sigma: ARGON_SIGMA,
            eps: ARGON_EPS,
            control: LjsWcaControl::default(),
        }
    }

    /// Create an `LJS-WCA` model with custom `(sigma, eps/k)` parameters.
    pub fn with_parameters(sigma: f64, eps: f64) -> Self {
        Self {
            sigma,
            eps,
            control: LjsWcaControl::default(),
        }
    }

    /// Replace `(sigma, eps/k)` on the model (`ljs_wca_set_pure_params`).
    pub fn set_parameters(&mut self, sigma: f64, eps: f64) -> Result<()> {
        positive(sigma, "segment diameter")?;
        positive(eps, "well depth")?;
        self.sigma = sigma;
        self.eps = eps;
        Ok(())
    }

    /// Segment diameter `sigma` in metre.
    pub fn sigma(&self) -> f64 {
        self.sigma
    }

    /// Well depth `eps/k_B` in kelvin.
    pub fn eps(&self) -> f64 {
        self.eps
    }

    /// Current model control flags.
    pub fn control(&self) -> LjsWcaControl {
        self.control
    }

    /// Replace model control flags (`ljs_wca_model_control`).
    pub fn set_control(&mut self, control: LjsWcaControl) {
        self.control = control;
    }

    /// WCA hard-sphere diameter and its first and second temperature derivatives.
    pub fn diameter(&self, temperature: f64) -> Result<(f64, f64, f64)> {
        positive(self.sigma, "segment diameter")?;
        positive(self.eps, "well depth")?;
        positive(temperature, "temperature")?;
        Ok(hard_sphere_diameter(self.sigma, self.eps, temperature))
    }

    /// Hypothetical pure-fluid packing fraction `eta` at `(T, rho)` (`calc_ljx_WCA_zeta`).
    pub fn packing_fraction(&self, temperature: f64, density: f64) -> Result<f64> {
        positive(density, "molar density")?;
        let (d, _, _) = self.diameter(temperature)?;
        Ok(std::f64::consts::PI / 6.0 * crate::N_A * d * d * d * density)
    }

    /// Packing fraction and all `(T, V, n)` derivatives (`calcZetaX_vdW_no_segments`).
    pub fn packing_fraction_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<PackingFractionDerivatives> {
        let (d, d_t, d_tt) = self.diameter(temperature)?;
        packing_fraction_vdw_no_segments(volume, &[moles], &[d], &[d_t], &[d_tt])
    }

    /// Prepare temperature-only quantities (including the 122 cavity quadrature nodes
    /// when `enable_hard_sphere && enable_cavity`) for repeated evaluations at `temperature`.
    pub fn prepare(&self, temperature: f64) -> Result<PreparedLjsWca> {
        positive(self.sigma, "segment diameter")?;
        positive(self.eps, "well depth")?;
        positive(temperature, "temperature")?;
        let (d, d_t, d_tt) = hard_sphere_diameter(self.sigma, self.eps, temperature);
        let (eta_per_rho, eta_per_rho_t, eta_per_rho_tt) = pure_eta_per_rho(d, d_t, d_tt);
        let rho_star_per_rho = crate::N_A * self.sigma * self.sigma * self.sigma;
        let y = R_SPLIT - d / self.sigma;
        let y_t = -d_t / self.sigma;
        let y_tt = -d_tt / self.sigma;
        let cavity = if self.control.enable_hard_sphere && self.control.enable_cavity {
            Some(CavityPreparedTemp::new(
                self.sigma,
                self.eps,
                d,
                d_t,
                d_tt,
                temperature,
            )?)
        } else {
            None
        };
        Ok(PreparedLjsWca {
            sigma: self.sigma,
            eps: self.eps,
            control: self.control,
            temperature,
            diameter: d,
            diameter_t: d_t,
            diameter_tt: d_tt,
            eta_per_rho,
            eta_per_rho_t,
            eta_per_rho_tt,
            rho_star_per_rho,
            y,
            y_t,
            y_tt,
            cavity,
        })
    }

    /// Molar reduced residual Helmholtz energy `a_res(T, rho)` and its 6 `(T, rho)` derivatives.
    pub fn residual_density(&self, temperature: f64, rho: f64) -> Result<DensityHelmholtz> {
        self.prepare(temperature)?.residual_density(rho)
    }

    /// Extensive reduced residual Helmholtz energy `F = n * a_res` and all 10 `(T, V, n)`
    /// derivatives (`calcFresLJs_WCA`).
    pub fn residual_helmholtz_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<PureHelmholtz> {
        self.prepare(temperature)?
            .residual_helmholtz_tvn(volume, moles)
    }

    /// Molar reduced hard-sphere contribution and all 10 `(T, V, n)` derivatives
    /// (`calc_ljs_hard_sphere`).
    pub fn hard_sphere_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<PureHelmholtz> {
        self.prepare(temperature)?.hard_sphere_tvn(volume, moles)
    }

    /// Molar reduced WCA soft-repulsion (cavity) contribution and all 10 `(T, V, n)`
    /// derivatives (`calc_wca_soft_repulsion`).
    pub fn soft_repulsion_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<PureHelmholtz> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        let (d, d_t, d_tt) = self.diameter(temperature)?;
        crate::hardsphere_wca::cavity_integral_lj_fres(
            self.sigma,
            self.eps,
            d,
            d_t,
            d_tt,
            temperature,
            volume,
            moles,
            true,
        )
    }

    /// Molar reduced dispersion contribution (hard-sphere and cavity disabled) and all 10
    /// `(T, V, n)` derivatives (`calc_ljs_dispersion`).
    pub fn dispersion_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<PureHelmholtz> {
        self.prepare(temperature)?.dispersion_tvn(volume, moles)
    }

    /// Pressure (Pa) and compressibility factor `Z` at `(T, V, n)`.
    pub fn pressure(&self, temperature: f64, volume: f64, moles: f64) -> Result<(f64, f64)> {
        self.prepare(temperature)?.pressure(volume, moles)
    }

    /// Pressure and first `(T, V, n)` derivatives at `(T, V, n)`.
    pub fn pressure_tv(&self, temperature: f64, volume: f64, moles: f64) -> Result<PurePressure> {
        self.prepare(temperature)?.pressure_tv(volume, moles)
    }

    /// Full `(T, V, n)` thermodynamic state (caloric properties, speed of sound, fugacity).
    pub fn tv_state(&self, temperature: f64, volume: f64, moles: f64) -> Result<LjsTvState> {
        self.prepare(temperature)?.tv_state(volume, moles)
    }

    /// Legacy-compatible summary state at `(T, V, n)`.
    pub fn state(&self, temperature: f64, volume: f64, moles: f64) -> Result<LjsState> {
        Ok(self.tv_state(temperature, volume, moles)?.to_legacy_state())
    }

    /// Speed of sound (m/s) at `(T, V, n)`.
    pub fn speed_of_sound_tv(&self, temperature: f64, volume: f64, moles: f64) -> Result<f64> {
        Ok(self.tv_state(temperature, volume, moles)?.speed_of_sound)
    }

    /// All mechanically stable molar volumes (m^3/mol) at `(T, P)`, sorted ascending.
    pub fn volumes(&self, temperature: f64, pressure: f64) -> Result<Vec<f64>> {
        self.prepare(temperature)?.volumes(pressure)
    }

    /// Physically stable molar volume (m^3/mol) at `(T, P)` for `phase`.
    pub fn volume(&self, temperature: f64, pressure: f64, phase: Phase) -> Result<f64> {
        self.prepare(temperature)?.volume(pressure, phase)
    }

    /// Exact source SAFT non-association volume iteration (`pc_saft_nonassoc_volume_solver`).
    pub fn baseline_volume_iteration(
        &self,
        temperature: f64,
        pressure: f64,
        phase: Phase,
    ) -> Result<LjsVolumeIteration> {
        self.prepare(temperature)?
            .baseline_volume_iteration(pressure, phase)
    }

    /// Full `(T, P)` thermodynamic state and `(T, P)` derivatives for `phase`.
    pub fn state_tp(&self, temperature: f64, pressure: f64, phase: Phase) -> Result<LjsTpState> {
        self.prepare(temperature)?.state_tp(pressure, phase)
    }

    /// Estimated critical temperature (K) and critical molar volume (m^3/mol)
    /// (`saft_interface::estimate_critical_parameters`).
    pub fn estimate_critical_parameters(&self) -> Result<(f64, f64)> {
        estimate_ljs_critical_parameters(self.sigma, self.eps)
    }

    /// Fundamental Measure Theory compatibility flags `(is_fmt_consistent, na_enabled)`
    /// (`saft_interface::test_fmt_compatibility`).
    pub fn fmt_compatibility(&self) -> (bool, bool) {
        (true, false)
    }

    /// Dimensionless van der Waals attractive energy constant (`saft_interface::alpha`).
    pub fn alpha_vdw(&self) -> f64 {
        ALPHA_LJS
    }

    /// Reduced splined Lennard-Jones pair potential `u/eps` and reduced force `-du*/dr*`.
    pub fn reduced_potential(&self, separations: &[f64]) -> (Vec<f64>, Vec<f64>) {
        ljs_potential_reduced(separations)
    }

    /// Pair potential `u/k_B` (K) and radial force `-du/(k_B dr)` (K/m) with default cap `500.0`.
    pub fn pair_potential(&self, separations: &[f64]) -> Result<(Vec<f64>, Vec<f64>)> {
        ljs_pair_potential(self.sigma, self.eps, separations, 500.0)
    }

    /// Pair potential `u/k_B` (K) and radial force `-du/(k_B dr)` (K/m) with custom reduced cap.
    pub fn pair_potential_with_cap(
        &self,
        separations: &[f64],
        max_pot_val: f64,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        ljs_pair_potential(self.sigma, self.eps, separations, max_pot_val)
    }
}

impl PreparedLjsWca {
    /// Temperature (K) at which this state was prepared.
    pub fn temperature(&self) -> f64 {
        self.temperature
    }

    /// Cached WCA hard-sphere diameter and its first two temperature derivatives.
    pub fn diameter(&self) -> (f64, f64, f64) {
        (self.diameter, self.diameter_t, self.diameter_tt)
    }

    /// Cached molar packing-fraction factor `eta / rho` (m^3/mol).
    pub fn eta_per_rho(&self) -> f64 {
        self.eta_per_rho
    }

    /// Dispersion-only molar `(T, rho)` Helmholtz derivatives.
    fn dispersion_density(&self, rho: f64) -> Result<DensityHelmholtz> {
        positive(rho, "molar density")?;
        let scale = self.rho_star_per_rho;
        let rho_star = scale * rho;
        let mut total = DensityHelmholtz::zero();
        let enabled = [
            self.control.enable_a1,
            self.control.enable_a2,
            self.control.enable_a3,
            self.control.enable_a4,
        ];
        for (slot, &on) in enabled.iter().enumerate() {
            if !on {
                continue;
            }
            let index = slot + 1;
            let (a, a_r, a_y, a_rr, a_yy, a_ry) = wca_term(index, self.eps, self.y, rho_star)?;
            let term = DensityHelmholtz {
                value: a,
                dt: a_y * self.y_t,
                drho: a_r * scale,
                dtt: a_yy * self.y_t * self.y_t + a_y * self.y_tt,
                dtrho: a_ry * self.y_t * scale,
                drhorho: a_rr * scale * scale,
            };
            total = total.add(term.scale_inv_t_power(self.temperature, index as i32));
        }
        Ok(total)
    }

    /// Fast density-only `(a_res, da_res/drho, d^2 a_res/drho^2)` for volume solving.
    fn density_only(&self, rho: f64) -> Result<(f64, f64, f64)> {
        positive(rho, "molar density")?;
        let mut a = 0.0;
        let mut a_r = 0.0;
        let mut a_rr = 0.0;
        if self.control.enable_hard_sphere {
            let eta = self.eta_per_rho * rho;
            let (hs, hs_e, hs_ee) = pure_hard_sphere_eta(eta, self.control.kolafa_pure_reference)?;
            a += hs;
            a_r += hs_e * self.eta_per_rho;
            a_rr += hs_ee * self.eta_per_rho * self.eta_per_rho;
            if let Some(cav) = &self.cavity {
                let (c, c_e, c_ee) = cav.evaluate_eta_only(eta)?;
                a += c;
                a_r += c_e * self.eta_per_rho;
                a_rr += c_ee * self.eta_per_rho * self.eta_per_rho;
            }
        }
        let scale = self.rho_star_per_rho;
        let scale2 = scale * scale;
        let rho_star = scale * rho;
        let inv_t = 1.0 / self.temperature;
        let mut inv_t_pow = inv_t;
        let enabled = [
            self.control.enable_a1,
            self.control.enable_a2,
            self.control.enable_a3,
            self.control.enable_a4,
        ];
        for (slot, &on) in enabled.iter().enumerate() {
            if on {
                let (val, val_r, _, val_rr, _, _) = wca_term(slot + 1, self.eps, self.y, rho_star)?;
                a += inv_t_pow * val;
                a_r += inv_t_pow * val_r * scale;
                a_rr += inv_t_pow * val_rr * scale2;
            }
            inv_t_pow *= inv_t;
        }
        if !a.is_finite() || !a_r.is_finite() || !a_rr.is_finite() {
            return Err(Error::InvalidInput("non-finite LJS-WCA density evaluation"));
        }
        Ok((a, a_r, a_rr))
    }

    fn pressure_and_dp_drho(&self, rho: f64) -> Result<(f64, f64)> {
        let (_, a_r, a_rr) = self.density_only(rho)?;
        let rt = crate::R * self.temperature;
        let p = rt * (rho + rho * rho * a_r);
        let dp_drho = rt * (1.0 + 2.0 * rho * a_r + rho * rho * a_rr);
        Ok((p, dp_drho))
    }

    /// Molar reduced residual Helmholtz energy `a_res(T, rho)` and all 6 `(T, rho)` derivatives.
    pub fn residual_density(&self, rho: f64) -> Result<DensityHelmholtz> {
        positive(rho, "molar density")?;
        let mut total = DensityHelmholtz::zero();
        if self.control.enable_hard_sphere {
            let hs = pure_hard_sphere_density(
                self.eta_per_rho,
                self.eta_per_rho_t,
                self.eta_per_rho_tt,
                rho,
                self.control.kolafa_pure_reference,
            )?;
            total = total.add(hs);
            if let Some(cav) = &self.cavity {
                let c = cav.evaluate_density(
                    self.eta_per_rho,
                    self.eta_per_rho_t,
                    self.eta_per_rho_tt,
                    rho,
                )?;
                total = total.add(c);
            }
        }
        let disp = self.dispersion_density(rho)?;
        total = total.add(disp);
        Ok(total)
    }

    /// Extensive reduced residual Helmholtz energy `F = n * a_res` and all 10 `(T, V, n)` derivatives.
    pub fn residual_helmholtz_tvn(&self, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        self.residual_density(moles / volume)?
            .extensive(volume, moles)
    }

    /// Molar reduced hard-sphere contribution and all 10 `(T, V, n)` derivatives.
    pub fn hard_sphere_tvn(&self, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        let hs = pure_hard_sphere_density(
            self.eta_per_rho,
            self.eta_per_rho_t,
            self.eta_per_rho_tt,
            moles / volume,
            self.control.kolafa_pure_reference,
        )?;
        hs.molar_tvn(volume, moles)
    }

    /// Molar reduced WCA soft-repulsion contribution and all 10 `(T, V, n)` derivatives.
    pub fn soft_repulsion_tvn(&self, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        let rho = moles / volume;
        let dh = if let Some(cav) = &self.cavity {
            cav.evaluate_density(
                self.eta_per_rho,
                self.eta_per_rho_t,
                self.eta_per_rho_tt,
                rho,
            )?
        } else {
            let cav = CavityPreparedTemp::new(
                self.sigma,
                self.eps,
                self.diameter,
                self.diameter_t,
                self.diameter_tt,
                self.temperature,
            )?;
            cav.evaluate_density(
                self.eta_per_rho,
                self.eta_per_rho_t,
                self.eta_per_rho_tt,
                rho,
            )?
        };
        dh.molar_tvn(volume, moles)
    }

    /// Molar reduced dispersion contribution and all 10 `(T, V, n)` derivatives.
    pub fn dispersion_tvn(&self, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        self.dispersion_density(moles / volume)?
            .molar_tvn(volume, moles)
    }

    /// Pressure (Pa) and compressibility factor `Z` at `(V, n)`.
    pub fn pressure(&self, volume: f64, moles: f64) -> Result<(f64, f64)> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        let rho = moles / volume;
        let (p, _) = self.pressure_and_dp_drho(rho)?;
        let z = p / (rho * crate::R * self.temperature);
        Ok((p, z))
    }

    /// Pressure and first `(T, V, n)` derivatives at `(V, n)`.
    pub fn pressure_tv(&self, volume: f64, moles: f64) -> Result<PurePressure> {
        let f = self.residual_helmholtz_tvn(volume, moles)?;
        f.pressure(self.temperature, volume, moles)
    }

    /// Full `(T, V, n)` thermodynamic state at `(V, n)`.
    pub fn tv_state(&self, volume: f64, moles: f64) -> Result<LjsTvState> {
        let f = self.residual_helmholtz_tvn(volume, moles)?;
        evaluate_tv_state(self.temperature, volume, moles, f)
    }

    /// All mechanically stable molar volumes (m^3/mol) at `pressure`, sorted ascending.
    pub fn volumes(&self, pressure: f64) -> Result<Vec<f64>> {
        solve_physical_volumes(self.temperature, pressure, self.eta_per_rho, |rho| {
            self.pressure_and_dp_drho(rho)
        })
    }

    /// Physically stable molar volume (m^3/mol) at `pressure` for `phase`.
    pub fn volume(&self, pressure: f64, phase: Phase) -> Result<f64> {
        select_physical_volume(
            self.temperature,
            pressure,
            phase,
            self.eta_per_rho,
            |rho| self.pressure_and_dp_drho(rho),
            |rho| Ok(self.density_only(rho)?.0),
        )
    }

    /// Exact source SAFT non-association volume iteration (`pc_saft_nonassoc_volume_solver`).
    pub fn baseline_volume_iteration(
        &self,
        pressure: f64,
        phase: Phase,
    ) -> Result<LjsVolumeIteration> {
        baseline_saft_volume_iteration(self.temperature, pressure, phase, self.eta_per_rho, |rho| {
            self.pressure_and_dp_drho(rho)
        })
    }

    /// Full `(T, P)` thermodynamic state and `(T, P)` derivatives for `phase`.
    pub fn state_tp(&self, pressure: f64, phase: Phase) -> Result<LjsTpState> {
        let v = self.volume(pressure, phase)?;
        let f = self.residual_helmholtz_tvn(v, 1.0)?;
        evaluate_tp_state(self.temperature, pressure, v, f)
    }
}
