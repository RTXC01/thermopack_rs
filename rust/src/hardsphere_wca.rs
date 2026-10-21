//! The cavity terms of the White Bear hard-sphere reference.
//! Ported from `src/hardsphere_wca.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! The White Bear family adds the cavity of the hard-sphere reference and the
//! soft repulsion of the model. The source evaluates the cavity distribution
//! function, the Boltzmann factor of the Lennard-Jones potential and the
//! radius of the minimum of the cavity integral, and the native model
//! evaluates the same quantities.

/// The reduced separation of the source, the sixth root of two.
const R_SPLIT: f64 = 1.122_462_048_309_373;

/// The packing fraction of the source at which the cavity correlation changes.
const ETA_C: f64 = 0.52;

/// The cavity distribution function of the hard-sphere reference and its
/// derivatives.
#[derive(Debug, Clone, Copy)]
pub struct CavityDistribution {
    /// The value of the function.
    pub value: f64,
    /// The derivative with respect to the packing fraction.
    pub packing_fraction: f64,
    /// The second derivative with respect to the packing fraction.
    pub packing_fraction_second: f64,
    /// The derivative with respect to the temperature.
    pub temperature: f64,
    /// The second derivative with respect to the temperature.
    pub temperature_second: f64,
    /// The mixed derivative with respect to the packing fraction and the
    /// temperature.
    pub packing_fraction_temperature: f64,
}

/// The hard-sphere diameter of the White Bear reference and its derivatives.
///
/// The diameter is the segment diameter times the reduced diameter of the
/// source, and the derivatives are those with respect to the temperature.
pub fn hard_sphere_diameter(sigma: f64, eps_divk: f64, temperature: f64) -> (f64, f64, f64) {
    let t_star = temperature / eps_divk;
    let sqrt_t_star = t_star.sqrt();
    let denominator = 1.0 + sqrt_t_star;
    let d_star = 2.0f64.powf(1.0 / 6.0) * denominator.powf(-1.0 / 6.0);
    let d_star_t = -d_star / (12.0 * denominator * sqrt_t_star);
    let d_star_tt = -d_star_t * (7.0 / (12.0 * denominator * sqrt_t_star) + 1.0 / (2.0 * t_star));
    (
        sigma * d_star,
        sigma * d_star_t / eps_divk,
        sigma * d_star_tt / (eps_divk * eps_divk),
    )
}

/// The radius of the minimum of the cavity integral and its derivatives.
///
/// The radius is reduced by the segment diameter of the source, and the
/// derivatives are those with respect to the temperature.
pub fn minimum_cavity_radius(temperature: f64, eps_divk: f64) -> (f64, f64, f64) {
    let c = -f64::MIN_POSITIVE.ln() / eps_divk;
    let sqr = (c * temperature).sqrt();
    let x = 0.5 + 0.5 * sqr;
    let r = x.powf(-1.0 / 6.0);
    let r_t = -r * c / (24.0 * sqr * x);
    let r_tt = 7.0 * r_t * r_t / r - 0.5 * r_t / temperature;
    (r, r_t, r_tt)
}

/// The cavity correlation function of the hard spheres.
///
/// The packing fraction and the reduced separation are the inputs, and the
/// derivatives with respect to the separation and the temperature follow the
/// source.
pub fn cavity_distribution_function(
    eta: f64,
    r_star: f64,
    r_star_t: f64,
    r_star_tt: f64,
) -> CavityDistribution {
    let (value, packing_fraction, packing_fraction_second, r_first, r_second, r_mixed) = if eta
        <= ETA_C
    {
        let (a, b, c, a_e, b_e, c_e, a_ee, b_ee, c_ee) = cavity_parameters(eta);
        let log_value = a + b * r_star + c * r_star * r_star * r_star;
        let value = log_value.exp();
        let slope = a_e + b_e * r_star + c_e * r_star * r_star * r_star;
        let packing_fraction = value * slope;
        let packing_fraction_second =
            value * (slope * slope + a_ee + b_ee * r_star + c_ee * r_star * r_star * r_star);
        let r_first = value * (b + 3.0 * c * r_star * r_star);
        let r_second = value * ((b + 3.0 * c * r_star * r_star).powi(2) + 6.0 * c * r_star);
        let r_mixed =
            value * (slope * (b + 3.0 * c * r_star * r_star) + b_e + 3.0 * c_e * r_star * r_star);
        (
            value,
            packing_fraction,
            packing_fraction_second,
            r_first,
            r_second,
            r_mixed,
        )
    } else {
        let (a, b, c, a_e, b_e, c_e, a_ee, b_ee, c_ee) = cavity_parameters(ETA_C);
        let value = (a + b * r_star + c * r_star * r_star * r_star).exp();
        let slope = a_e + b_e * r_star + c_e * r_star * r_star * r_star;
        let packing_fraction = value * slope;
        let packing_fraction_second =
            value * (slope * slope + a_ee + b_ee * r_star + c_ee * r_star * r_star * r_star);
        let r_first = value * (b + 3.0 * c * r_star * r_star);
        let r_second = value * ((b + 3.0 * c * r_star * r_star).powi(2) + 6.0 * c * r_star);
        let r_mixed =
            value * (slope * (b + 3.0 * c * r_star * r_star) + b_e + 3.0 * c_e * r_star * r_star);
        let r_mixed_first = value
            * (2.0 * (b + 3.0 * c * r_star * r_star) * (b_e + 3.0 * c_e * r_star * r_star)
                + slope * (6.0 * c * r_star + (b + 3.0 * c * r_star * r_star).powi(2))
                + 6.0 * c_e * r_star);
        let r_mixed_second = packing_fraction_second * (b + 3.0 * c * r_star * r_star)
            + value
                * (2.0 * slope * (b_e + 3.0 * c_e * r_star * r_star)
                    + (b_ee + 3.0 * c_ee * r_star * r_star));
        let r_mixed_mixed = r_mixed_second * (b + 3.0 * c * r_star * r_star)
            + value
                * (6.0 * c * r_star * slope * slope
                    + 2.0
                        * (b + 3.0 * c * r_star * r_star)
                        * slope
                        * (b_e + 3.0 * c_e * r_star * r_star)
                    + 6.0 * c * r_star * (a_ee + b_ee * r_star + c_ee * r_star * r_star * r_star)
                    + (b + 3.0 * c * r_star * r_star) * (b_ee + 3.0 * c_ee * r_star * r_star)
                    + 2.0 * (b_e + 3.0 * c_e * r_star * r_star).powi(2)
                    + 2.0 * slope * 6.0 * c_e * r_star
                    + 6.0 * c_ee * r_star);
        let difference = eta - ETA_C;
        let value = value
            + packing_fraction * difference
            + 0.5 * packing_fraction_second * difference * difference;
        let packing_fraction = packing_fraction + packing_fraction_second * difference;
        let r_first =
            r_first + r_mixed * difference + 0.5 * r_mixed_second * difference * difference;
        let r_second =
            r_second + r_mixed_first * difference + 0.5 * r_mixed_mixed * difference * difference;
        let r_mixed = r_mixed + r_mixed_second * difference;
        (
            value,
            packing_fraction,
            packing_fraction_second,
            r_first,
            r_second,
            r_mixed,
        )
    };
    CavityDistribution {
        value,
        packing_fraction,
        packing_fraction_second,
        temperature: r_first * r_star_t,
        temperature_second: r_second * r_star_t * r_star_t + r_first * r_star_tt,
        packing_fraction_temperature: r_mixed * r_star_t,
    }
}

/// The Boltzmann factor of the Lennard-Jones potential and its derivatives.
///
/// The separation is reduced by the segment diameter, and the derivatives are
/// those with respect to the temperature.
pub fn lj_boltzmann_factor(
    r: f64,
    r_t: f64,
    r_tt: f64,
    temperature: f64,
    eps_divk: f64,
) -> (f64, f64, f64) {
    let (u, u_r, u_rr) = if r > R_SPLIT {
        (0.0, 0.0, 0.0)
    } else {
        let inverse = 1.0 / r;
        let inverse6 = inverse.powi(6);
        let inverse12 = inverse6 * inverse6;
        (
            4.0 * (inverse12 - inverse6) + 1.0,
            -4.0 * (12.0 * inverse12 - 6.0 * inverse6) * inverse,
            4.0 * (12.0 * 13.0 * inverse12 - 6.0 * 7.0 * inverse6) * inverse * inverse,
        )
    };
    let value = (-u * eps_divk / temperature).exp();
    let value_t =
        -value * (u_r * r_t * eps_divk / temperature - u * eps_divk / (temperature * temperature));
    let value_tt = value
        * ((u_r * r_t * eps_divk / temperature - u * eps_divk / (temperature * temperature))
            .powi(2)
            - (u_rr * r_t * r_t * eps_divk / temperature + u_r * r_tt * eps_divk / temperature
                - 2.0 * u_r * r_t * eps_divk / (temperature * temperature)
                + 2.0 * u * eps_divk / (temperature * temperature * temperature)));
    (value, value_t, value_tt)
}

/// The parameters of the cavity correlation function of the source.
fn cavity_parameters(eta: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    let denominator = (1.0 - eta).powi(3);
    let a = (3.0 - eta) / denominator - 3.0;
    let b = -3.0 * eta * (2.0 - eta) / denominator;
    let c = ((2.0 - eta) / denominator / 2.0).ln()
        - eta * (2.0 - 6.0 * eta + 3.0 * eta * eta) / denominator;
    let denominator = (1.0 - eta) * denominator;
    let a_e = -(2.0 * (eta - 4.0)) / denominator;
    let b_e = (3.0 * (eta * eta - 2.0 * eta - 2.0)) / denominator;
    let c_e =
        -((2.0 * (eta - 4.0) * eta + 7.0) * eta * eta + eta + 1.0) / ((eta - 2.0) * denominator);
    let denominator = (1.0 - eta) * denominator;
    let a_ee = -(6.0 * (eta - 5.0)) / denominator;
    let b_ee = (6.0 * (eta * eta - 2.0 * eta - 5.0)) / denominator;
    let c_ee = -(2.0 * eta.powi(5) - 10.0 * eta.powi(4) + 5.0 * eta.powi(3) + 31.0 * eta * eta
        - 29.0 * eta
        - 11.0)
        / ((eta - 2.0).powi(2) * denominator);
    (a, b, c, a_e, b_e, c_e, a_ee, b_ee, c_ee)
}
