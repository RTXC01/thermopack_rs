//! The pair correlation function of the additive hard-sphere model of
//! Boublik.
//!
//! The original `src/hardsphere_bmcsl.f90:calc_bmcsl_dgij_dzeta` evaluates the
//! reduced Helmholtz energy of a pair of components and its first and second
//! derivatives with respect to the packing fractions of the mixture. The
//! expression is the one of the paper of Boublik, "Hard-Sphere Equation of
//! State", J. Chem. Phys. 53, 471 (1970), as the original states.

/// The pair correlation function of the additive hard-sphere model and its
/// derivatives with respect to the packing fractions.
///
/// The derivatives are with respect to the second and the third packing
/// fraction and with respect to the reduced diameter of the pair, in the order
/// of the original `dg_dzeta` and `d2g_dzeta2`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BmcslPairCorrelation {
    /// The pair correlation function.
    pub value: f64,
    /// The first derivatives with respect to the packing fractions.
    pub dzeta: [f64; 3],
    /// The second derivatives with respect to the packing fractions.
    pub d2zeta: [[f64; 3]; 3],
}

/// The pair correlation function of the additive hard-sphere model.
///
/// `zeta_2` is the second packing fraction of the mixture, `zeta_3` the third
/// one and `mu` the reduced diameter of the pair, the ratio of the product of
/// the hard-sphere diameters and their sum.
pub fn boublik_pair_correlation(zeta_2: f64, zeta_3: f64, mu: f64) -> BmcslPairCorrelation {
    let div_diff_1 = 1.0 / (1.0 - zeta_3);
    let div_diff_2 = div_diff_1 * div_diff_1;
    let div_diff_3 = div_diff_1 * div_diff_2;
    let div_diff_4 = div_diff_2 * div_diff_2;
    let mu_2 = mu * mu;
    let zeta_2_2 = zeta_2 * zeta_2;

    let value = div_diff_1 + 3.0 * zeta_2 * div_diff_2 * mu + 2.0 * zeta_2_2 * mu_2 * div_diff_3;

    // The first derivatives with respect to (zeta_2, zeta_3, mu).
    let dzeta = [
        3.0 * mu * div_diff_2 + 4.0 * div_diff_3 * mu_2 * zeta_2,
        div_diff_2 + 6.0 * zeta_2 * div_diff_3 * mu + 6.0 * zeta_2_2 * div_diff_4 * mu_2,
        3.0 * zeta_2 * div_diff_2 + 4.0 * zeta_2_2 * div_diff_3 * mu,
    ];

    // The second derivatives with respect to (zeta_2, zeta_3, mu).
    let mut d2zeta = [[0.0; 3]; 3];
    d2zeta[0][0] = 4.0 * div_diff_3 * mu_2;
    d2zeta[1][1] = 2.0 * div_diff_3
        + 18.0 * zeta_2 * mu * div_diff_4
        + 24.0 * zeta_2_2 * mu_2 * div_diff_3 * div_diff_2;
    d2zeta[2][2] = 4.0 * zeta_2_2 * div_diff_3;
    d2zeta[2][0] = 3.0 * div_diff_2 + 8.0 * zeta_2 * mu * div_diff_3;
    d2zeta[0][2] = d2zeta[2][0];
    d2zeta[0][1] = 6.0 * div_diff_3 * mu + 12.0 * zeta_2 * div_diff_4 * mu_2;
    d2zeta[1][0] = d2zeta[0][1];
    d2zeta[2][1] = 6.0 * zeta_2 * div_diff_3 + 12.0 * zeta_2_2 * div_diff_4 * mu;
    d2zeta[1][2] = d2zeta[2][1];

    BmcslPairCorrelation {
        value,
        dzeta,
        d2zeta,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The derivatives with respect to the packing fractions are checked
    /// against a central difference of the expression.
    #[test]
    fn the_derivatives_match_the_difference_quotient() {
        let states = [(0.05, 0.1, 0.3), (0.12, 0.25, 0.4), (0.3, 0.4, 0.25)];
        for (zeta_2, zeta_3, mu) in states {
            let value = boublik_pair_correlation(zeta_2, zeta_3, mu);
            let step = 1.0e-6;
            for index in 0..3 {
                let mut up = [zeta_2, zeta_3, mu];
                let mut down = [zeta_2, zeta_3, mu];
                up[index] += step;
                down[index] -= step;
                let quotient = (boublik_pair_correlation(up[0], up[1], up[2]).value
                    - boublik_pair_correlation(down[0], down[1], down[2]).value)
                    / (2.0 * step);
                let difference = (quotient - value.dzeta[index]).abs();
                assert!(
                    difference < 1.0e-6,
                    "the first derivative {index} differs by {difference}"
                );
                for second in 0..3 {
                    let mut up = [zeta_2, zeta_3, mu];
                    let mut down = [zeta_2, zeta_3, mu];
                    up[index] += step;
                    down[index] -= step;
                    let quotient = (boublik_pair_correlation(up[0], up[1], up[2]).dzeta[second]
                        - boublik_pair_correlation(down[0], down[1], down[2]).dzeta[second])
                        / (2.0 * step);
                    let difference = (quotient - value.d2zeta[index][second]).abs();
                    assert!(
                        difference < 1.0e-4,
                        "the second derivative {index},{second} differs by {difference}"
                    );
                }
            }
        }
    }
}
