//! The perturbation terms of the White Bear family of the splined
//! Lennard-Jones model.
//! Ported from `src/lj_splined.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.

// The coefficients of the source carry the precision of the source.
#![allow(clippy::excessive_precision)]
//!
//! The White Bear family replaces the perturbation terms of the
//! Barker-Henderson theory with the soft repulsion of the White Bear model.
//! The source evaluates the four terms with polynomials in the reduced
//! distance of the hard-sphere centres and the reduced density, and the native
//! model evaluates the same terms.

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

/// The first order term of the White Bear family.
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

/// The second order term of the White Bear family.
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
    let e2u = (2.0 * u[0][0]).exp();
    let mut tanhu = [0.0; 3];
    tanhu[0] = (e2u - 1.0) / (e2u + 1.0);
    tanhu[1] = 1.0 - tanhu[0] * tanhu[0];
    tanhu[2] = -2.0 * tanhu[0] * tanhu[1];
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

/// The third order term of the White Bear family.
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
        // The source uses the polynomials four and five of its array, which
        // are the entries three and four of the native array.
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
    let e2r = (2.0 * r[0][0]).exp();
    let mut tanhr = [0.0; 3];
    tanhr[0] = (e2r - 1.0) / (e2r + 1.0);
    tanhr[1] = 1.0 - tanhr[0] * tanhr[0];
    tanhr[2] = -2.0 * tanhr[0] * tanhr[1];
    tanhr[0] = 1.0 + c[10] * tanhr[0];
    tanhr[1] *= c[10];
    tanhr[2] *= c[10];
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

/// The fourth order term of the White Bear family.
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
    let e2k = (2.0 * k[0][0]).exp();
    let mut tanhk = [0.0; 3];
    tanhk[0] = (e2k - 1.0) / (e2k + 1.0);
    tanhk[1] = 1.0 - tanhk[0] * tanhk[0];
    tanhk[2] = -2.0 * tanhk[0] * tanhk[1];
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

/// One perturbation term of the White Bear family and its derivatives.
///
/// The index selects the term and follows the source, the first order term is
/// one and the fourth order term is four. The reduced temperature and the
/// reduced density are dimensionless. The method returns the term, its
/// derivative with respect to the reduced density, its derivative with respect
/// to the reduced distance and the second derivatives, in the order of the
/// source.
pub fn perturbation_term(
    index: usize,
    reduced_temperature: f64,
    reduced_density: f64,
) -> Result<(f64, f64, f64, f64, f64, f64), crate::Error> {
    if reduced_temperature <= 0.0 {
        return Err(crate::Error::InvalidInput(
            "the reduced temperature must be positive",
        ));
    }
    if reduced_density <= 0.0 {
        return Err(crate::Error::InvalidInput(
            "the reduced density must be positive",
        ));
    }
    // The reduced distance of the hard-sphere centres of the source.
    let x0 = (2.0 / (1.0 + reduced_temperature.sqrt())).powf(1.0 / 6.0);
    let y = 2.0f64.powf(1.0 / 6.0) - x0;
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
            return Err(crate::Error::InvalidInput(
                "the perturbation term index must be one to four",
            ));
        }
    };
    Ok((value, a_r, a_y, a_rr, a_yy, a_ry))
}

/// One perturbation term of the White Bear family with the prefactor of the
/// source and its derivatives.
///
/// The index selects the term and follows the source, the first order term is
/// one and the fourth order term is four. The well depth is divided by the
/// Boltzmann constant, the reduced distance and the reduced density are
/// dimensionless. The method returns the term, its derivative with respect to
/// the reduced density, its derivative with respect to the reduced distance
/// and the second derivatives, in the order of the source.
pub fn wca_term(
    index: usize,
    eps_divk: f64,
    reduced_distance: f64,
    reduced_density: f64,
) -> Result<(f64, f64, f64, f64, f64, f64), crate::Error> {
    let (value, value_r, value_y, value_rr, value_yy, value_ry) = match index {
        1 => first_order(reduced_distance, reduced_density),
        2 => second_order(reduced_distance, reduced_density),
        3 => third_order(reduced_distance, reduced_density),
        4 => fourth_order(reduced_distance, reduced_density),
        _ => {
            return Err(crate::Error::InvalidInput(
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
