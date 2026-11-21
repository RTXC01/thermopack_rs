//! Lee-Kesler (1975) generalized corresponding states equation of state.
//!
//! Evaluates thermodynamic properties via Plöcker et al. (1982) critical mixing
//! rules and BWR generalized fluid equations for a simple fluid (w=0) and a
//! reference fluid (n-octane, w=0.3978).
//!
//! # Parity with ThermoPack
//! Parity against `src/leekesler.f90` and `.reference/thermopack.dll`:
//! - Plöcker critical temperature, volume, pressure, acentric factor mixing rules
//! - 12 BWR constants per fluid (simple and reference)
//! - Exact Newton-Raphson reduced volume solver with phase-consistency analysis
//! - Enthalpy and entropy departures, total caloric properties (with IdealGas)
//! - Component fugacity coefficients with analytical compositional derivatives

use crate::error::positive;
use crate::ideal::IdealGas;
use crate::{Component, Error, Phase, R, Result};

const ETA: f64 = 0.25;
const W_REF: f64 = 0.3978;
const VR_MINIMUM: f64 = 0.05;
const VR_MAXIMUM: f64 = 1.0e10;

// BWR constants: index 0 = simple fluid, index 1 = reference fluid
const B1: [f64; 2] = [0.1181193, 0.2026579];
const B2: [f64; 2] = [0.265728, 0.331511];
const B3: [f64; 2] = [0.154790, 0.027655];
const B4: [f64; 2] = [0.030323, 0.203488];
const C1: [f64; 2] = [0.0236744, 0.0313385];
const C2: [f64; 2] = [0.0186984, 0.0503618];
const C3: [f64; 2] = [0.0, 0.016901];
const C4: [f64; 2] = [0.042724, 0.041577];
const D1: [f64; 2] = [0.0000155488, 0.000048736];
const D2: [f64; 2] = [0.0000623689, 0.00000740336];
const BETA: [f64; 2] = [0.65392, 1.226];
const GAMMA: [f64; 2] = [0.060167, 0.03754];

#[derive(Debug, Clone, Copy)]
pub struct BwrCoeffs {
    pub b: f64,
    pub c: f64,
    pub d: f64,
    pub e: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct BwrCoeffsDiff1 {
    pub b_tr: f64,
    pub c_tr: f64,
    pub d_tr: f64,
    pub e_tr: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct BwrCoeffsDiff2 {
    pub b_trtr: f64,
    pub c_trtr: f64,
    pub d_trtr: f64,
    pub e_trtr: f64,
}

impl BwrCoeffs {
    #[inline]
    pub fn new(tr: f64, sr: usize) -> Self {
        let tr_inv = 1.0 / tr;
        let tr_inv2 = tr_inv * tr_inv;
        let tr_inv3 = tr_inv2 * tr_inv;
        Self {
            b: B1[sr] - B2[sr] * tr_inv - B3[sr] * tr_inv2 - B4[sr] * tr_inv3,
            c: C1[sr] - C2[sr] * tr_inv + C3[sr] * tr_inv3,
            d: D1[sr] + D2[sr] * tr_inv,
            e: C4[sr] * tr_inv3,
        }
    }
}

impl BwrCoeffsDiff1 {
    #[inline]
    pub fn new(tr: f64, sr: usize) -> Self {
        let tr_inv = 1.0 / tr;
        let tr_inv2 = tr_inv * tr_inv;
        let tr_inv3 = tr_inv2 * tr_inv;
        let tr_inv4 = tr_inv2 * tr_inv2;
        Self {
            b_tr: B2[sr] * tr_inv2 + 2.0 * B3[sr] * tr_inv3 + 3.0 * B4[sr] * tr_inv4,
            c_tr: C2[sr] * tr_inv2 - 3.0 * C3[sr] * tr_inv4,
            d_tr: -D2[sr] * tr_inv2,
            e_tr: -3.0 * C4[sr] * tr_inv4,
        }
    }
}

impl BwrCoeffsDiff2 {
    #[inline]
    pub fn new(tr: f64, sr: usize) -> Self {
        let tr_inv = 1.0 / tr;
        let tr_inv2 = tr_inv * tr_inv;
        let tr_inv3 = tr_inv2 * tr_inv;
        let tr_inv4 = tr_inv2 * tr_inv2;
        let tr_inv5 = tr_inv4 * tr_inv;
        Self {
            b_trtr: -2.0 * B2[sr] * tr_inv3 - 6.0 * B3[sr] * tr_inv4 - 12.0 * B4[sr] * tr_inv5,
            c_trtr: -2.0 * C2[sr] * tr_inv3 + 12.0 * C3[sr] * tr_inv5,
            d_trtr: 2.0 * D2[sr] * tr_inv3,
            e_trtr: 12.0 * C4[sr] * tr_inv5,
        }
    }
}

#[inline]
fn f_solver(moles: f64, vr: f64, b: f64, c: f64, d: f64, e: f64, sr: usize) -> f64 {
    let be = BETA[sr];
    let ga = GAMMA[sr];
    let ga_inv = 1.0 / ga;
    let vr_inv = 1.0 / vr;
    let vr_inv2 = vr_inv * vr_inv;
    let vr_inv5 = vr_inv2 * vr_inv2 * vr_inv;

    moles
        * (b * vr_inv + 0.5 * c * vr_inv2 + 0.2 * d * vr_inv5 + 0.5 * e * ga_inv * (be + 1.0)
            - e * (-ga * vr_inv2).exp() * (0.5 * ga_inv * (be + 1.0) + 0.5 * vr_inv2))
}

#[inline]
fn f_diff_tr(moles: f64, tr: f64, vr: f64, sr: usize) -> f64 {
    let d1 = BwrCoeffsDiff1::new(tr, sr);
    f_solver(moles, vr, d1.b_tr, d1.c_tr, d1.d_tr, d1.e_tr, sr)
}

#[inline]
fn f_diff2_tr(moles: f64, tr: f64, vr: f64, sr: usize) -> f64 {
    let d2 = BwrCoeffsDiff2::new(tr, sr);
    f_solver(moles, vr, d2.b_trtr, d2.c_trtr, d2.d_trtr, d2.e_trtr, sr)
}

#[inline]
fn f_diff_tr_vr(moles: f64, tr: f64, vr: f64, sr: usize) -> f64 {
    let d1 = BwrCoeffsDiff1::new(tr, sr);
    f_diff_vr(moles, vr, d1.b_tr, d1.c_tr, d1.d_tr, d1.e_tr, sr)
}

#[derive(Debug, Clone, Copy)]
struct FluidReducedState {
    z: f64,
    dz_dtr: f64,
    dz_dpr: f64,
    hr: f64,
    dhr_dtr: f64,
    dhr_dpr: f64,
    sr: f64,
    dsr_dtr: f64,
    dsr_dpr: f64,
    gr: f64,
    dgr_dtr: f64,
    dgr_dpr: f64,
    d2gr_dtr2: f64,
    d2gr_dtr_dpr: f64,
    d2gr_dpr2: f64,
}

#[inline]
fn evaluate_reduced_fluid(
    tr: f64,
    pr: f64,
    vr: f64,
    coeffs: &BwrCoeffs,
    sr: usize,
) -> Result<FluidReducedState> {
    if !tr.is_finite() || tr <= 0.0 || !pr.is_finite() || pr <= 0.0 || !vr.is_finite() || vr <= 0.0
    {
        return Err(Error::InvalidInput(
            "Lee-Kesler reduced state outside numerical range",
        ));
    }
    let f = f_solver(1.0, vr, coeffs.b, coeffs.c, coeffs.d, coeffs.e, sr);
    let f_tr = f_diff_tr(1.0, tr, vr, sr);
    let f_trtr = f_diff2_tr(1.0, tr, vr, sr);
    let f_vrvr = f_diff2_vr(1.0, vr, coeffs.b, coeffs.c, coeffs.d, coeffs.e, sr);
    let f_trvr = f_diff_tr_vr(1.0, tr, vr, sr);

    let z = pr * vr / tr;
    if !z.is_finite() || z <= 0.0 {
        return Err(Error::NotConverged(
            "Lee-Kesler reduced compressibility factor outside positive finite range",
        ));
    }
    let dp_dvr = -tr * (1.0 / (vr * vr) + f_vrvr);
    if !dp_dvr.is_finite() || dp_dvr == 0.0 {
        return Err(Error::NotConverged(
            "singular Lee-Kesler pressure-volume derivative",
        ));
    }
    let dp_dtr = pr / tr - tr * f_trvr;
    let dvr_dpr = 1.0 / dp_dvr;
    let dvr_dtr = -dp_dtr / dp_dvr;
    let dz_dpr = (vr + pr * dvr_dpr) / tr;
    let dz_dtr = (pr * dvr_dtr - z) / tr;

    let hr = tr * (-tr * f_tr + z - 1.0);
    let cvr = -tr * (2.0 * f_tr + tr * f_trtr);
    let dhr_dtr = cvr - 1.0 + tr * dp_dtr * dvr_dtr;
    let dhr_dpr = vr - tr * dvr_dtr;

    let sr_val = -(f + tr * f_tr) + z.ln();
    let dsr_dtr = dhr_dtr / tr;
    let dsr_dpr = 1.0 / pr - dvr_dtr;

    let gr = z - 1.0 - z.ln() + f;
    let tr2 = tr * tr;
    let pr2 = pr * pr;
    if !tr2.is_finite() || tr2 <= 0.0 || !pr2.is_finite() || pr2 <= 0.0 {
        return Err(Error::InvalidInput(
            "Lee-Kesler reduced state underflow/overflow",
        ));
    }
    let dgr_dtr = -hr / tr2;
    let dgr_dpr = (z - 1.0) / pr;
    let d2gr_dtr2 = -(dhr_dtr * tr - 2.0 * hr) / (tr2 * tr);
    let d2gr_dtr_dpr = dz_dtr / pr;
    let d2gr_dpr2 = (dz_dpr * pr - (z - 1.0)) / pr2;

    let state = FluidReducedState {
        z,
        dz_dtr,
        dz_dpr,
        hr,
        dhr_dtr,
        dhr_dpr,
        sr: sr_val,
        dsr_dtr,
        dsr_dpr,
        gr,
        dgr_dtr,
        dgr_dpr,
        d2gr_dtr2,
        d2gr_dtr_dpr,
        d2gr_dpr2,
    };
    if ![
        state.z,
        state.dz_dtr,
        state.dz_dpr,
        state.hr,
        state.dhr_dtr,
        state.dhr_dpr,
        state.sr,
        state.dsr_dtr,
        state.dsr_dpr,
        state.gr,
        state.dgr_dtr,
        state.dgr_dpr,
        state.d2gr_dtr2,
        state.d2gr_dtr_dpr,
        state.d2gr_dpr2,
    ]
    .iter()
    .all(|v| v.is_finite())
    {
        return Err(Error::InvalidInput(
            "non-finite Lee-Kesler reduced fluid derivatives",
        ));
    }
    Ok(state)
}

/// Scalar TP property value with analytical temperature, pressure, and
/// composition derivatives.
///
/// For intensive quantities (`Z`), `dn[i]` is the partial derivative
/// `dZ/dn_i`. For molar extensive quantities (`v`, `H`, `S`, `G`), `value`,
/// `dt`, and `dp` are molar (`[m^3/mol]`, `[J/mol]`, `[J/(mol K)]`) and
/// `dn[i]` is the partial molar derivative `d(N * value)/dn_i`.
#[derive(Debug, Clone)]
pub struct LeeKeslerDerivatives {
    pub value: f64,
    pub dt: f64,
    pub dp: f64,
    pub dn: Vec<f64>,
}

/// Log fugacity coefficients with analytical temperature, pressure, and
/// mole-number derivatives at constant `(T, P)`.
#[derive(Debug, Clone)]
pub struct LeeKeslerFugacityDerivatives {
    pub ln_phi: Vec<f64>,
    pub dt: Vec<f64>,
    pub dp: Vec<f64>,
    /// Row-major `nc * nc` matrix with entry `i * nc + j` equal to
    /// `d ln(phi_i) / d n_j`.
    pub dn: Vec<f64>,
}

/// Complete single-solve TP state for the Lee-Kesler equation of state.
#[derive(Debug, Clone)]
pub struct LeeKeslerState {
    pub compressibility_factor: LeeKeslerDerivatives,
    pub specific_volume: LeeKeslerDerivatives,
    pub residual_enthalpy: LeeKeslerDerivatives,
    pub enthalpy: LeeKeslerDerivatives,
    pub residual_entropy: LeeKeslerDerivatives,
    pub entropy: LeeKeslerDerivatives,
    pub residual_gibbs: LeeKeslerDerivatives,
    pub gibbs: LeeKeslerDerivatives,
    pub ln_fugacity: LeeKeslerFugacityDerivatives,
    pub heat_capacity_p: f64,
    pub heat_capacity_v: f64,
    pub speed_of_sound: Option<f64>,
}

#[inline]
fn f_diff_vr(moles: f64, vr: f64, b: f64, c: f64, d: f64, e: f64, sr: usize) -> f64 {
    let be = BETA[sr];
    let ga = GAMMA[sr];
    let vr_inv = 1.0 / vr;
    let vr_inv2 = vr_inv * vr_inv;
    let vr_inv3 = vr_inv * vr_inv2;
    let vr_inv6 = vr_inv3 * vr_inv3;

    -moles
        * (b * vr_inv2
            + c * vr_inv3
            + d * vr_inv6
            + e * vr_inv3 * (be + ga * vr_inv2) * (-ga * vr_inv2).exp())
}

#[inline]
fn f_diff2_vr(moles: f64, vr: f64, b: f64, c: f64, d: f64, e: f64, sr: usize) -> f64 {
    let be = BETA[sr];
    let ga = GAMMA[sr];
    let vr_inv = 1.0 / vr;
    let vr_inv2 = vr_inv * vr_inv;
    let vr_inv3 = vr_inv * vr_inv2;
    let vr_inv4 = vr_inv2 * vr_inv2;
    let vr_inv7 = vr_inv4 * vr_inv3;

    moles
        * (2.0 * b * vr_inv3
            + 3.0 * c * vr_inv4
            + 6.0 * d * vr_inv7
            + e * vr_inv4
                * (-ga * vr_inv2).exp()
                * (3.0 * be + ga * (5.0 - 2.0 * be) * vr_inv2 - 2.0 * ga * ga * vr_inv4))
}

#[inline]
fn f_diff3_vr(moles: f64, vr: f64, b: f64, c: f64, d: f64, e: f64, sr: usize) -> f64 {
    let be = BETA[sr];
    let ga = GAMMA[sr];
    let ga2 = ga * ga;
    let ga3 = ga * ga2;
    let vr_inv = 1.0 / vr;
    let vr_inv2 = vr_inv * vr_inv;
    let vr_inv4 = vr_inv2 * vr_inv2;
    let vr_inv5 = vr_inv4 * vr_inv;
    let vr_inv6 = vr_inv4 * vr_inv2;
    let vr_inv8 = vr_inv4 * vr_inv4;

    -moles
        * (6.0 * b * vr_inv4
            + 12.0 * c * vr_inv5
            + 42.0 * d * vr_inv8
            + e * vr_inv5
                * (-ga * vr_inv2).exp()
                * (12.0 * be
                    + ga * (30.0 - 18.0 * be) * vr_inv2
                    + ga2 * (4.0 * be - 26.0) * vr_inv4
                    + 4.0 * ga3 * vr_inv6))
}

#[inline]
fn fz_with_diff(pr: f64, tr: f64, z: f64, coeffs: &BwrCoeffs, sr: usize) -> (f64, f64, f64) {
    let vr = z * tr / pr;
    let fd_vr = f_diff_vr(1.0, vr, coeffs.b, coeffs.c, coeffs.d, coeffs.e, sr);
    let fd2_vr = f_diff2_vr(1.0, vr, coeffs.b, coeffs.c, coeffs.d, coeffs.e, sr);
    let fd3_vr = f_diff3_vr(1.0, vr, coeffs.b, coeffs.c, coeffs.d, coeffs.e, sr);

    let fz = z - (1.0 - vr * fd_vr);
    let tr_pr = tr / pr;
    let fzd = 1.0 + tr_pr * (fd_vr + vr * fd2_vr);
    let fzd2 = tr_pr * tr_pr * (2.0 * fd2_vr + vr * fd3_vr);
    (fz, fzd, fzd2)
}

#[inline]
fn fz(pr: f64, tr: f64, z: f64, coeffs: &BwrCoeffs, sr: usize) -> f64 {
    let vr = z * tr / pr;
    z - (1.0 - vr * f_diff_vr(1.0, vr, coeffs.b, coeffs.c, coeffs.d, coeffs.e, sr))
}

fn set_max_min_z(tr: f64, pr: f64, phase: usize) -> (f64, f64, f64) {
    let z_min = pr * VR_MINIMUM / tr;
    let z_max = pr * VR_MAXIMUM / tr;
    let z = if phase == 1 { z_min } else { z_max };
    (z, z_max, z_min)
}

fn z_initial(
    tr: f64,
    pr: f64,
    phase: usize,
    sr: usize,
    coeffs: &BwrCoeffs,
) -> (f64, f64, f64, bool, bool) {
    let tolerance = f64::EPSILON * 200.0;
    let max_iterations = 200;
    let safety_steps = 2;

    let mut z_min = pr * VR_MINIMUM / tr;
    let mut z_max = pr * VR_MAXIMUM / tr;
    let mut zd_min = z_min;
    let mut zd_max = z_max;

    let mut z = if phase == 1 { z_min } else { 1.0 };

    let mut f = fz(pr, tr, z, coeffs, sr);
    if f < 0.0 && phase == 2 {
        z -= 2.0 * f;
        f = fz(pr, tr, z, coeffs, sr);
    }
    let mut z_old = z * 2.0;
    let mut f_old = f;
    let mut f_diff = 1.0e100;
    let mut solved = false;
    let mut has_phase = false;

    for i in 1..=max_iterations {
        let f_diff_old = f_diff;
        let (f_val, fd, fd2) = fz_with_diff(pr, tr, z, coeffs, sr);
        f = f_val;
        f_diff = fd;

        if phase == 2 && f_diff.abs() > f_diff_old.abs() {
            z = 0.25 * z + 0.75 * z_old;
            let (f_val2, fd2_val, _) = fz_with_diff(pr, tr, z, coeffs, sr);
            f = f_val2;
            f_diff = fd2_val;
        }

        if f_diff >= 0.0 {
            if f == 0.0 {
                z_min = z;
                z_max = z;
            } else if f < 0.0 {
                z_min = z;
            } else {
                z_max = z;
            }
        } else if phase == 1 {
            z_max = z;
        } else {
            z_min = z;
        }

        if f * f_old <= 0.0 {
            z = if phase == 2 { z_max } else { z_min };
            has_phase = true;
            return (z, z_min, z_max, solved, has_phase);
        } else {
            f_old = f;
        }

        if phase == 1 {
            if f_diff < 0.0 {
                zd_max = z;
            } else {
                zd_min = z;
            }
        } else if f_diff < 0.0 {
            zd_min = z;
        } else {
            zd_max = z;
        }

        let dz = if phase == 2 && i <= safety_steps {
            let delta = 0.5 * f * fd2 / (f_diff * f_diff);
            -(f / f_diff) * (1.0 + delta)
        } else {
            (-f_diff / fd2).max(-z * 0.33)
        };

        let z_temp = z + dz;
        if ((z_temp - z) / z).abs() < tolerance || f_diff.abs() < tolerance {
            z = z_temp;
            if f * f_old > 0.0 {
                solved = true;
            }
            break;
        }

        z_old = z;
        if z_temp > zd_max {
            if z < zd_max {
                z = (z + zd_max) * 0.5;
            } else {
                z = (z + zd_min) * 0.5;
            }
        } else if z_temp < zd_min {
            if z > zd_min {
                z = (z + zd_min) * 0.5;
            } else {
                z = (z + zd_max) * 0.5;
            }
        } else {
            z = z_temp;
        }
        f_old = f;
    }

    (z, z_min, z_max, solved, has_phase)
}

#[allow(clippy::too_many_arguments)]
fn z_newt_raps(
    tr: f64,
    pr: f64,
    phase: usize,
    sr: usize,
    z_init: f64,
    mut z_max: f64,
    mut z_min: f64,
    coeffs: &BwrCoeffs,
) -> Result<f64> {
    let tolerance = f64::EPSILON * 200.0;
    let max_iterations = 200;

    if !z_init.is_finite() || z_init <= 0.0 || !z_min.is_finite() || !z_max.is_finite() {
        return Err(Error::NotConverged(
            "Lee-Kesler reduced volume initialization outside finite positive range",
        ));
    }

    let mut z = z_init;
    for _ in 1..=max_iterations {
        let (f_val, f_diff, f_diff2) = fz_with_diff(pr, tr, z, coeffs, sr);
        let mut f = f_val;
        let mut fd = f_diff;
        let mut fd2 = f_diff2;
        if !f.is_finite() || !fd.is_finite() || !fd2.is_finite() {
            return Err(Error::NotConverged(
                "Lee-Kesler BWR evaluation became non-finite",
            ));
        }
        if fd < 0.0 {
            z = if phase == 1 { z_max } else { z_min };
            let (f_val2, fd_val2, fd2_val2) = fz_with_diff(pr, tr, z, coeffs, sr);
            f = f_val2;
            fd = fd_val2;
            fd2 = fd2_val2;
            if !f.is_finite() || !fd.is_finite() || !fd2.is_finite() {
                return Err(Error::NotConverged(
                    "Lee-Kesler BWR evaluation became non-finite",
                ));
            }
        }
        if fd == 0.0 {
            return Err(Error::NotConverged(
                "singular Lee-Kesler BWR derivative during Newton solve",
            ));
        }
        let delta = 0.5 * f * fd2 / (fd * fd);
        let dz = -(f / fd) * (1.0 + delta);
        let z_temp = z + dz;
        if !z_temp.is_finite() {
            return Err(Error::NotConverged(
                "Lee-Kesler Newton step became non-finite",
            ));
        }

        if ((z_temp - z) / z).abs() < tolerance {
            z = z_temp;
            break;
        }

        let z_old = z;
        if z_temp > z_max {
            if z < z_max {
                z = (z + z_max) * 0.5;
            } else {
                z = (z + z_min) * 0.5;
            }
        } else if z_temp < z_min {
            if z > z_min {
                z = (z + z_min) * 0.5;
            } else {
                z = (z + z_max) * 0.5;
            }
        } else {
            z = z_temp;
        }

        let fold = f;
        let mut f_next = fz(pr, tr, z, coeffs, sr);
        if !f_next.is_finite() {
            return Err(Error::NotConverged(
                "Lee-Kesler BWR residual became non-finite",
            ));
        }
        if fold * f_next < 0.0 && f_next.abs() >= 0.9 * fold.abs() {
            z = 0.5 * (z_old + z);
            f_next = fz(pr, tr, z, coeffs, sr);
        }
        if f_next < 0.0 {
            z_min = z;
        } else {
            z_max = z;
        }
    }

    if !z.is_finite() || z <= 0.0 {
        return Err(Error::NotConverged(
            "Lee-Kesler compressibility root outside positive finite range",
        ));
    }
    let vr = z * tr / pr;
    if !vr.is_finite() || vr <= 0.0 {
        return Err(Error::NotConverged(
            "Lee-Kesler reduced volume outside positive finite range",
        ));
    }
    let residual = fz(pr, tr, z, coeffs, sr);
    if !residual.is_finite() || residual.abs() > 1.0e-7 * z.abs().max(1.0) {
        return Err(Error::NotConverged(
            "Lee-Kesler reduced volume solver did not converge",
        ));
    }
    Ok(vr)
}

fn calc_reduced_volume(
    tr: f64,
    pr: f64,
    phase: usize,
    coeffs_simp: &BwrCoeffs,
    coeffs_ref: &BwrCoeffs,
) -> Result<(f64, f64)> {
    if !tr.is_finite()
        || tr <= 0.0
        || !pr.is_finite()
        || pr <= 0.0
        || !(pr / tr).is_finite()
        || (pr / tr) <= 0.0
        || !(tr / pr).is_finite()
        || !(pr * pr).is_finite()
        || (pr * pr) <= 0.0
        || !(tr * tr).is_finite()
        || (tr * tr) <= 0.0
    {
        return Err(Error::InvalidInput(
            "Lee-Kesler reduced temperature or pressure outside numerical range",
        ));
    }
    let (mut z_simp, mut z_min_simp, mut z_max_simp, _, has_phase_simp) =
        z_initial(tr, pr, phase, 0, coeffs_simp);
    let (mut z_ref, mut z_min_ref, mut z_max_ref, _, has_phase_ref) =
        z_initial(tr, pr, phase, 1, coeffs_ref);

    let mut phase_simp = phase;
    let mut phase_ref = phase;
    let alt_phase = if phase == 1 { 2 } else { 1 };

    if has_phase_simp != has_phase_ref {
        if !has_phase_simp {
            phase_simp = alt_phase;
            let (zs, z_max_s, z_min_s) = set_max_min_z(tr, pr, phase_simp);
            z_simp = zs;
            z_max_simp = z_max_s;
            z_min_simp = z_min_s;

            let z_saved = z_ref;
            let z_max_saved = z_max_ref;
            let z_min_saved = z_min_ref;
            let (zr, z_min_r, z_max_r, _, has_phase_r) =
                z_initial(tr, pr, alt_phase, 1, coeffs_ref);
            if has_phase_r {
                phase_ref = alt_phase;
                z_ref = zr;
                z_min_ref = z_min_r;
                z_max_ref = z_max_r;
            } else {
                z_ref = z_saved;
                z_max_ref = z_max_saved;
                z_min_ref = z_min_saved;
            }
        } else {
            phase_ref = alt_phase;
            let (zr, z_max_r, z_min_r) = set_max_min_z(tr, pr, phase_ref);
            z_ref = zr;
            z_max_ref = z_max_r;
            z_min_ref = z_min_r;

            let z_saved = z_simp;
            let z_max_saved = z_max_simp;
            let z_min_saved = z_min_simp;
            let (zs, z_min_s, z_max_s, _, has_phase_s) =
                z_initial(tr, pr, alt_phase, 0, coeffs_simp);
            if has_phase_s {
                phase_simp = alt_phase;
                z_simp = zs;
                z_min_simp = z_min_s;
                z_max_simp = z_max_s;
            } else {
                z_simp = z_saved;
                z_max_simp = z_max_saved;
                z_min_simp = z_min_saved;
            }
        }
    } else if !has_phase_simp && !has_phase_ref {
        phase_simp = alt_phase;
        let (zs, z_max_s, z_min_s) = set_max_min_z(tr, pr, phase_simp);
        phase_ref = alt_phase;
        z_min_simp = z_min_s;
        z_max_simp = z_max_s;
        z_simp = zs;
        z_min_ref = z_min_simp;
        z_max_ref = z_max_simp;
        z_ref = z_simp;
    }

    let vr_simp = z_newt_raps(
        tr,
        pr,
        phase_simp,
        0,
        z_simp,
        z_max_simp,
        z_min_simp,
        coeffs_simp,
    )?;
    let vr_ref = z_newt_raps(
        tr, pr, phase_ref, 1, z_ref, z_max_ref, z_min_ref, coeffs_ref,
    )?;

    Ok((vr_simp, vr_ref))
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn ln_phi_m(
    t: f64,
    p: f64,
    tr: f64,
    pr: f64,
    vr: f64,
    moles: f64,
    coeffs: &BwrCoeffs,
    sr: usize,
) -> (f64, f64, f64) {
    let z_fac = vr * pr / tr;
    let f = f_solver(moles, vr, coeffs.b, coeffs.c, coeffs.d, coeffs.e, sr);
    let ln_phi = moles * (z_fac - 1.0 - z_fac.ln()) + f;
    let h_rt = -tr * f_diff_tr(moles, tr, vr, sr) + moles * (z_fac - 1.0);
    let dlnphi_dt = -h_rt / t;
    let dlnphi_dp = moles * (z_fac - 1.0) / p;
    (ln_phi, dlnphi_dt, dlnphi_dp)
}

#[inline]
fn f_diff_ni(
    moles: f64,
    tr: f64,
    vr: f64,
    coeffs: &BwrCoeffs,
    sr: usize,
    tr_i: f64,
    vr_i: f64,
) -> f64 {
    let f_n = f_solver(moles, vr, coeffs.b, coeffs.c, coeffs.d, coeffs.e, sr) / moles;
    let f_tr = f_diff_tr(moles, tr, vr, sr);
    let f_vr = f_diff_vr(moles, vr, coeffs.b, coeffs.c, coeffs.d, coeffs.e, sr);

    f_n + f_tr * tr_i + f_vr * vr_i
}

/// The Lee-Kesler Generalized Equation of State model.
#[derive(Debug, Clone)]
pub struct LeeKesler {
    components: Vec<Component>,
    kij: Vec<f64>,
    reference: String,
    ideals: Vec<IdealGas>,
}

impl LeeKesler {
    /// Create a new Lee-Kesler model with Default parameters.
    pub fn new(names: &[&str]) -> Result<Self> {
        Self::new_with_reference(names, "Default")
    }

    /// Create a new Lee-Kesler model with a specific parameter reference.
    pub fn new_with_reference(names: &[&str], reference: &str) -> Result<Self> {
        if names.is_empty() {
            return Err(Error::InvalidInput("at least one component required"));
        }
        let mut components = Vec::with_capacity(names.len());
        let mut ideals = Vec::with_capacity(names.len());
        for &name in names {
            let comp = *crate::component(name)?;
            ideals.push(IdealGas::new(&comp)?);
            components.push(comp);
        }

        let n = names.len();
        let mut kij = vec![1.0; n * n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    kij[i * n + j] = 1.0;
                } else {
                    let a = components[i].ident;
                    let b = components[j].ident;
                    let k =
                        crate::component::interaction_with_reference("LK", a, b, false, reference);
                    kij[i * n + j] = if k == 0.0 { 1.0 } else { k };
                }
            }
        }

        Ok(Self {
            components,
            kij,
            reference: reference.to_string(),
            ideals,
        })
    }

    pub fn components(&self) -> &[Component] {
        &self.components
    }

    pub fn reference(&self) -> &str {
        &self.reference
    }

    pub fn kij(&self, i: usize, j: usize) -> Result<f64> {
        let n = self.components.len();
        if i >= n || j >= n {
            return Err(Error::InvalidInput("component index out of bounds"));
        }
        Ok(self.kij[i * n + j])
    }

    pub fn set_kij(&mut self, i: usize, j: usize, value: f64) -> Result<()> {
        let n = self.components.len();
        if i >= n || j >= n {
            return Err(Error::InvalidInput("component index out of bounds"));
        }
        if !value.is_finite() || value <= 0.0 {
            return Err(Error::InvalidInput("finite positive kij required"));
        }
        self.kij[i * n + j] = value;
        self.kij[j * n + i] = value;
        Ok(())
    }

    /// Prepare a thermodynamic state for the specified mole fractions.
    ///
    /// Precomputes mixture critical values and their compositional derivatives,
    /// enabling fast, allocation-free evaluation across (T, P, Phase).
    pub fn prepare<'a>(&'a self, z: &[f64]) -> Result<PreparedLeeKesler<'a>> {
        let n = self.components.len();
        if z.len() != n {
            return Err(Error::InvalidInput("composition length mismatch"));
        }
        let total_moles: f64 = z.iter().copied().sum();
        if !total_moles.is_finite() || total_moles <= 0.0 {
            return Err(Error::InvalidInput("positive total moles required"));
        }
        for &zi in z {
            if !zi.is_finite() || zi < 0.0 {
                return Err(Error::InvalidInput("non-negative mole numbers required"));
            }
        }

        let moles = total_moles;
        let moles_inv = 1.0 / moles;

        // Individual pure critical volumes: vci = R * zci * tci / pci
        let mut vc = Vec::with_capacity(n);
        for comp in &self.components {
            vc.push(
                R * comp.critical_compressibility * comp.critical_temperature
                    / comp.critical_pressure,
            );
        }

        // Precompute pair critical volumes/temperatures and first/second
        // compositional derivatives in O(n^2) time.
        let mut vc_pair = vec![0.0; n * n];
        let mut vc_eta_tc_pair = vec![0.0; n * n];

        let mut vc_m = 0.0;
        let mut tc_m_sum = 0.0;
        let mut omega_m = 0.0;

        for j in 0..n {
            let tcj = self.components[j].critical_temperature;
            let vcj_cbrt = vc[j].cbrt();
            for k in 0..n {
                let tck = self.components[k].critical_temperature;
                let s_jk = vcj_cbrt + vc[k].cbrt();
                let vc_jk = 0.125 * (s_jk * s_jk * s_jk);
                let tc_jk = self.kij[j * n + k] * (tcj * tck).sqrt();
                let vc_eta_tc = vc_jk.powf(ETA) * tc_jk;

                let idx = j * n + k;
                vc_pair[idx] = vc_jk;
                vc_eta_tc_pair[idx] = vc_eta_tc;

                let n_prod = z[j] * z[k];
                vc_m += n_prod * vc_jk;
                tc_m_sum += n_prod * vc_eta_tc;
            }
            omega_m += z[j] * self.components[j].acentric_factor;
        }

        vc_m *= moles_inv * moles_inv;
        let vc_m_inv = 1.0 / vc_m;
        let vc_m_inv_eta = vc_m_inv.powf(ETA);

        let tc_m = tc_m_sum * moles_inv * moles_inv * vc_m_inv_eta;
        omega_m *= moles_inv;
        let zc_m = 0.2905 - 0.085 * omega_m;
        let pc_m = R * zc_m * tc_m * vc_m_inv;

        // First compositional derivatives
        let mut dvc_dni = Vec::with_capacity(n);
        let mut dtc_dni = Vec::with_capacity(n);
        let mut domega_dni = Vec::with_capacity(n);
        let mut dzc_dni = Vec::with_capacity(n);
        let mut dpc_dni = Vec::with_capacity(n);
        let mut comp_sum = Vec::with_capacity(n);

        for (i, comp) in self.components.iter().enumerate() {
            let mut v1 = 0.0;
            let mut csum = 0.0;
            for (l, &zl) in z.iter().enumerate() {
                let idx = i * n + l;
                v1 += zl * vc_pair[idx];
                csum += zl * vc_eta_tc_pair[idx];
            }

            let dvc = 2.0 * moles_inv * (moles_inv * v1 - vc_m);
            let t1 = -ETA * tc_m * vc_m_inv * dvc;
            let t2 = 2.0 * csum * vc_m_inv_eta * moles_inv * moles_inv;
            let t3 = -2.0 * tc_m * moles_inv;
            let dtc = t1 + t2 + t3;

            let domega = moles_inv * (comp.acentric_factor - omega_m);
            let dzc = -0.085 * domega;
            let dpc = pc_m * ((1.0 / zc_m) * dzc + (1.0 / tc_m) * dtc - vc_m_inv * dvc);

            dvc_dni.push(dvc);
            dtc_dni.push(dtc);
            domega_dni.push(domega);
            dzc_dni.push(dzc);
            dpc_dni.push(dpc);
            comp_sum.push(csum);
        }

        // Second compositional derivatives (symmetric n*n row-major matrices)
        let mut d2vc_dn = vec![0.0; n * n];
        let mut d2tc_dn = vec![0.0; n * n];
        let mut d2pc_dn = vec![0.0; n * n];
        let mut d2omega_dn = vec![0.0; n * n];
        let mut d2zc_dn = vec![0.0; n * n];

        let moles_inv2 = moles_inv * moles_inv;
        let moles_inv3 = moles_inv2 * moles_inv;
        let vc_m_inv2 = vc_m_inv * vc_m_inv;

        for i in 0..n {
            for j in i..n {
                let idx = i * n + j;
                let rev = j * n + i;

                let d2vc = 2.0 * moles_inv2 * (vc_pair[idx] - vc_m)
                    - 2.0 * moles_inv * (dvc_dni[i] + dvc_dni[j]);
                let d2w = -moles_inv * (domega_dni[i] + domega_dni[j]);
                let d2z = -0.085 * d2w;

                let d2tc = 2.0 * moles_inv2 * tc_m
                    - 2.0 * moles_inv * dtc_dni[j]
                    - ETA * vc_m_inv * dvc_dni[i] * dtc_dni[j]
                    + ETA * tc_m * vc_m_inv2 * dvc_dni[i] * dvc_dni[j]
                    - ETA * tc_m * vc_m_inv * d2vc
                    - (4.0 * moles_inv3 * vc_m_inv_eta
                        + 2.0 * ETA * dvc_dni[j] * moles_inv2 * vc_m_inv_eta * vc_m_inv)
                        * comp_sum[i]
                    + 2.0 * moles_inv2 * vc_m_inv_eta * vc_eta_tc_pair[idx];

                let d2pc = vc_m_inv
                    * (R * (d2z * tc_m
                        + dzc_dni[i] * dtc_dni[j]
                        + dzc_dni[j] * dtc_dni[i]
                        + zc_m * d2tc)
                        - d2vc * pc_m
                        - dvc_dni[i] * dpc_dni[j]
                        - dvc_dni[j] * dpc_dni[i]);

                d2vc_dn[idx] = d2vc;
                d2vc_dn[rev] = d2vc;
                d2omega_dn[idx] = d2w;
                d2omega_dn[rev] = d2w;
                d2zc_dn[idx] = d2z;
                d2zc_dn[rev] = d2z;
                d2tc_dn[idx] = d2tc;
                d2tc_dn[rev] = d2tc;
                d2pc_dn[idx] = d2pc;
                d2pc_dn[rev] = d2pc;
            }
        }

        Ok(PreparedLeeKesler {
            model: self,
            z: z.to_vec(),
            moles,
            tc_m,
            vc_m,
            pc_m,
            zc_m,
            omega_m,
            dtc_dni,
            dvc_dni,
            dpc_dni,
            dzc_dni,
            domega_dni,
            d2vc_dn,
            d2tc_dn,
            d2pc_dn,
            d2omega_dn,
            d2zc_dn,
        })
    }

    /// Direct evaluation of compressibility factor.
    pub fn compressibility_factor(&self, t: f64, p: f64, z: &[f64], phase: Phase) -> Result<f64> {
        self.prepare(z)?.compressibility_factor(t, p, phase)
    }

    /// Direct evaluation of compressibility factor with `(T, P, n)` derivatives.
    pub fn compressibility_factor_derivatives(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        phase: Phase,
    ) -> Result<LeeKeslerDerivatives> {
        self.prepare(z)?
            .compressibility_factor_derivatives(t, p, phase)
    }

    /// Direct evaluation of specific volume.
    pub fn specific_volume(&self, t: f64, p: f64, z: &[f64], phase: Phase) -> Result<f64> {
        self.prepare(z)?.specific_volume(t, p, phase)
    }

    /// Direct evaluation of specific volume with `(T, P, n)` derivatives.
    pub fn specific_volume_derivatives(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        phase: Phase,
    ) -> Result<LeeKeslerDerivatives> {
        self.prepare(z)?.specific_volume_derivatives(t, p, phase)
    }

    /// Direct evaluation of residual enthalpy with `(T, P, n)` derivatives.
    pub fn residual_enthalpy_derivatives(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        phase: Phase,
    ) -> Result<LeeKeslerDerivatives> {
        self.prepare(z)?.residual_enthalpy_derivatives(t, p, phase)
    }

    /// Direct evaluation of total enthalpy.
    pub fn enthalpy(&self, t: f64, p: f64, z: &[f64], phase: Phase) -> Result<f64> {
        self.prepare(z)?.enthalpy(t, p, phase)
    }

    /// Direct evaluation of total enthalpy with `(T, P, n)` derivatives.
    pub fn enthalpy_derivatives(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        phase: Phase,
    ) -> Result<LeeKeslerDerivatives> {
        self.prepare(z)?.enthalpy_derivatives(t, p, phase)
    }

    /// Direct evaluation of residual entropy with `(T, P, n)` derivatives.
    pub fn residual_entropy_derivatives(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        phase: Phase,
    ) -> Result<LeeKeslerDerivatives> {
        self.prepare(z)?.residual_entropy_derivatives(t, p, phase)
    }

    /// Direct evaluation of total entropy.
    pub fn entropy(&self, t: f64, p: f64, z: &[f64], phase: Phase) -> Result<f64> {
        self.prepare(z)?.entropy(t, p, phase)
    }

    /// Direct evaluation of total entropy with `(T, P, n)` derivatives.
    ///
    /// For absent or trace components (`x_i < 1e-100`), the ideal-gas mixing term
    /// in the partial molar entropy `dn[i]` uses ThermoPack's `max(x_i, 1e-100)`
    /// log regularization so `dn[i]` remains finite while the mixture value uses
    /// the exact `0 * ln(0) = 0` limit.
    pub fn entropy_derivatives(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        phase: Phase,
    ) -> Result<LeeKeslerDerivatives> {
        self.prepare(z)?.entropy_derivatives(t, p, phase)
    }

    /// Direct evaluation of residual Gibbs energy with `(T, P, n)` derivatives.
    pub fn residual_gibbs_derivatives(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        phase: Phase,
    ) -> Result<LeeKeslerDerivatives> {
        self.prepare(z)?.residual_gibbs_derivatives(t, p, phase)
    }

    /// Direct evaluation of total Gibbs energy with `(T, P, n)` derivatives.
    ///
    /// For absent or trace components (`x_i < 1e-100`), the ideal-gas mixing term
    /// in the chemical potential `dn[i]` uses ThermoPack's `max(x_i, 1e-100)`
    /// log regularization so `dn[i]` remains finite while the mixture value uses
    /// the exact `0 * ln(0) = 0` limit.
    pub fn gibbs_derivatives(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        phase: Phase,
    ) -> Result<LeeKeslerDerivatives> {
        self.prepare(z)?.gibbs_derivatives(t, p, phase)
    }

    /// Direct evaluation of log fugacity coefficients.
    pub fn ln_fugacity_coefficients(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        phase: Phase,
    ) -> Result<Vec<f64>> {
        self.prepare(z)?.ln_fugacity_coefficients(t, p, phase)
    }

    /// Direct evaluation of log fugacity coefficients with `(T, P, n)` derivatives.
    pub fn ln_fugacity_derivatives(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        phase: Phase,
    ) -> Result<LeeKeslerFugacityDerivatives> {
        self.prepare(z)?.ln_fugacity_derivatives(t, p, phase)
    }

    /// Direct evaluation of fugacity coefficients.
    pub fn fugacity_coefficients(
        &self,
        t: f64,
        p: f64,
        z: &[f64],
        phase: Phase,
    ) -> Result<Vec<f64>> {
        self.prepare(z)?.fugacity_coefficients(t, p, phase)
    }

    /// Direct evaluation of molar isobaric heat capacity `Cp` `[J/(mol K)]`.
    pub fn heat_capacity_p(&self, t: f64, p: f64, z: &[f64], phase: Phase) -> Result<f64> {
        self.prepare(z)?.heat_capacity_p(t, p, phase)
    }

    /// Direct evaluation of molar isochoric heat capacity `Cv` `[J/(mol K)]`.
    pub fn heat_capacity_v(&self, t: f64, p: f64, z: &[f64], phase: Phase) -> Result<f64> {
        self.prepare(z)?.heat_capacity_v(t, p, phase)
    }

    /// Direct evaluation of single-phase speed of sound `[m/s]`.
    pub fn speed_of_sound(&self, t: f64, p: f64, z: &[f64], phase: Phase) -> Result<f64> {
        self.prepare(z)?.speed_of_sound(t, p, phase)
    }

    /// Single-solve evaluation of all TP properties and analytical derivatives.
    pub fn state_tp(&self, t: f64, p: f64, z: &[f64], phase: Phase) -> Result<LeeKeslerState> {
        self.prepare(z)?.state_tp(t, p, phase)
    }
}

/// Prepared Lee-Kesler state with precomputed mixture critical properties.
#[derive(Debug, Clone)]
pub struct PreparedLeeKesler<'a> {
    model: &'a LeeKesler,
    z: Vec<f64>,
    moles: f64,
    tc_m: f64,
    vc_m: f64,
    pc_m: f64,
    zc_m: f64,
    omega_m: f64,
    dtc_dni: Vec<f64>,
    dvc_dni: Vec<f64>,
    dpc_dni: Vec<f64>,
    dzc_dni: Vec<f64>,
    domega_dni: Vec<f64>,
    d2vc_dn: Vec<f64>,
    d2tc_dn: Vec<f64>,
    d2pc_dn: Vec<f64>,
    d2omega_dn: Vec<f64>,
    d2zc_dn: Vec<f64>,
}

impl<'a> PreparedLeeKesler<'a> {
    pub fn tc_m(&self) -> f64 {
        self.tc_m
    }

    pub fn vc_m(&self) -> f64 {
        self.vc_m
    }

    pub fn pc_m(&self) -> f64 {
        self.pc_m
    }

    pub fn zc_m(&self) -> f64 {
        self.zc_m
    }

    pub fn omega_m(&self) -> f64 {
        self.omega_m
    }

    pub fn dtc_dni(&self) -> &[f64] {
        &self.dtc_dni
    }

    pub fn dvc_dni(&self) -> &[f64] {
        &self.dvc_dni
    }

    pub fn dpc_dni(&self) -> &[f64] {
        &self.dpc_dni
    }

    pub fn dzc_dni(&self) -> &[f64] {
        &self.dzc_dni
    }

    pub fn domega_dni(&self) -> &[f64] {
        &self.domega_dni
    }

    #[allow(clippy::too_many_arguments)]
    fn gibbs_dep_from_volumes(
        &self,
        t: f64,
        p: f64,
        tr: f64,
        pr: f64,
        vr_simp: f64,
        vr_ref: f64,
        coeffs_simp: &BwrCoeffs,
        coeffs_ref: &BwrCoeffs,
    ) -> Result<f64> {
        let (ln_phi_m_simp, _, _) = ln_phi_m(t, p, tr, pr, vr_simp, self.moles, coeffs_simp, 0);
        let (ln_phi_m_ref, _, _) = ln_phi_m(t, p, tr, pr, vr_ref, self.moles, coeffs_ref, 1);
        let g_dep =
            R * t * (ln_phi_m_simp + (self.omega_m / W_REF) * (ln_phi_m_ref - ln_phi_m_simp))
                / self.moles;
        if !g_dep.is_finite() {
            return Err(Error::InvalidInput(
                "non-finite Lee-Kesler Gibbs energy departure",
            ));
        }
        Ok(g_dep)
    }

    /// Solve `(tr, pr, coeffs_simp, coeffs_ref, vr_simp, vr_ref)` once for the
    /// requested phase. For [`Phase::MinimumGibbs`], liquid and vapor are each
    /// solved once and the winning `(vr_simp, vr_ref)` pair is reused directly
    /// without a third solve.
    fn solve_reduced_volumes(
        &self,
        t: f64,
        p: f64,
        phase: Phase,
    ) -> Result<(f64, f64, BwrCoeffs, BwrCoeffs, f64, f64)> {
        positive(t, "temperature")?;
        positive(p, "pressure")?;
        let tr = t / self.tc_m;
        let pr = p / self.pc_m;
        let coeffs_simp = BwrCoeffs::new(tr, 0);
        let coeffs_ref = BwrCoeffs::new(tr, 1);
        let (vr_simp, vr_ref) = match phase {
            Phase::Liquid => calc_reduced_volume(tr, pr, 1, &coeffs_simp, &coeffs_ref)?,
            Phase::Vapor => calc_reduced_volume(tr, pr, 2, &coeffs_simp, &coeffs_ref)?,
            Phase::MinimumGibbs => {
                let liq = calc_reduced_volume(tr, pr, 1, &coeffs_simp, &coeffs_ref).and_then(
                    |(vs, vr)| {
                        let g = self.gibbs_dep_from_volumes(
                            t,
                            p,
                            tr,
                            pr,
                            vs,
                            vr,
                            &coeffs_simp,
                            &coeffs_ref,
                        )?;
                        Ok((g, vs, vr))
                    },
                );
                let vap = calc_reduced_volume(tr, pr, 2, &coeffs_simp, &coeffs_ref).and_then(
                    |(vs, vr)| {
                        let g = self.gibbs_dep_from_volumes(
                            t,
                            p,
                            tr,
                            pr,
                            vs,
                            vr,
                            &coeffs_simp,
                            &coeffs_ref,
                        )?;
                        Ok((g, vs, vr))
                    },
                );
                match (liq, vap) {
                    (Ok((gl, vsl, vrl)), Ok((gv, vsv, vrv))) => {
                        if gl <= gv {
                            (vsl, vrl)
                        } else {
                            (vsv, vrv)
                        }
                    }
                    (Ok((_, vsl, vrl)), Err(_)) => (vsl, vrl),
                    (Err(_), Ok((_, vsv, vrv))) => (vsv, vrv),
                    (Err(e), Err(_)) => return Err(e),
                }
            }
        };
        Ok((tr, pr, coeffs_simp, coeffs_ref, vr_simp, vr_ref))
    }

    /// Compressibility factor Z.
    pub fn compressibility_factor(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        let (tr, pr, _, _, vr_simp, vr_ref) = self.solve_reduced_volumes(t, p, phase)?;
        let z_simp = pr * vr_simp / tr;
        let z_ref = pr * vr_ref / tr;
        let z = z_simp + (self.omega_m / W_REF) * (z_ref - z_simp);
        if !z.is_finite() {
            return Err(Error::InvalidInput(
                "non-finite Lee-Kesler compressibility factor",
            ));
        }
        Ok(z)
    }

    /// Specific volume [m^3/mol].
    pub fn specific_volume(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        let z = self.compressibility_factor(t, p, phase)?;
        let v = z * R * t / p;
        if !v.is_finite() {
            return Err(Error::InvalidInput("non-finite Lee-Kesler specific volume"));
        }
        Ok(v)
    }

    /// Residual enthalpy [J/mol].
    pub fn residual_enthalpy(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        let (tr, pr, _, _, vr_simp, vr_ref) = self.solve_reduced_volumes(t, p, phase)?;

        let z_simp = pr * vr_simp / tr;
        let z_ref = pr * vr_ref / tr;

        let h_simp =
            R * t * (-tr * f_diff_tr(self.moles, tr, vr_simp, 0) + self.moles * (z_simp - 1.0));
        let h_ref =
            R * t * (-tr * f_diff_tr(self.moles, tr, vr_ref, 1) + self.moles * (z_ref - 1.0));

        let h_res = (h_simp + (self.omega_m / W_REF) * (h_ref - h_simp)) / self.moles;
        if !h_res.is_finite() {
            return Err(Error::InvalidInput(
                "non-finite Lee-Kesler residual enthalpy",
            ));
        }
        Ok(h_res)
    }

    /// Total enthalpy [J/mol] (ideal gas mixture + residual).
    pub fn enthalpy(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        let h_res = self.residual_enthalpy(t, p, phase)?;
        let mut h_ideal = 0.0;
        let n = self.model.components.len();
        for i in 0..n {
            let props = self.model.ideals[i].properties(t, p)?;
            h_ideal += (self.z[i] / self.moles) * props.enthalpy;
        }
        Ok(h_ideal + h_res)
    }

    /// Residual entropy [J/(mol K)].
    pub fn residual_entropy(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        let (tr, pr, coeffs_simp, coeffs_ref, vr_simp, vr_ref) =
            self.solve_reduced_volumes(t, p, phase)?;

        let z_simp = pr * vr_simp / tr;
        let z_ref = pr * vr_ref / tr;

        let f_simp = f_solver(
            self.moles,
            vr_simp,
            coeffs_simp.b,
            coeffs_simp.c,
            coeffs_simp.d,
            coeffs_simp.e,
            0,
        );
        let s_vdep_simp = -R * (f_simp + tr * f_diff_tr(self.moles, tr, vr_simp, 0));
        let s_simp = s_vdep_simp + self.moles * R * z_simp.ln();

        let f_ref = f_solver(
            self.moles,
            vr_ref,
            coeffs_ref.b,
            coeffs_ref.c,
            coeffs_ref.d,
            coeffs_ref.e,
            1,
        );
        let s_vdep_ref = -R * (f_ref + tr * f_diff_tr(self.moles, tr, vr_ref, 1));
        let s_ref = s_vdep_ref + self.moles * R * z_ref.ln();

        let s_res = (s_simp + (self.omega_m / W_REF) * (s_ref - s_simp)) / self.moles;
        if !s_res.is_finite() {
            return Err(Error::InvalidInput(
                "non-finite Lee-Kesler residual entropy",
            ));
        }
        Ok(s_res)
    }

    /// Total entropy [J/(mol K)] (ideal gas mixture + residual).
    pub fn entropy(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        let s_res = self.residual_entropy(t, p, phase)?;
        let mut s_ideal = 0.0;
        let n = self.model.components.len();
        for i in 0..n {
            let zi = self.z[i] / self.moles;
            let props = self.model.ideals[i].properties(t, p)?;
            s_ideal += zi * props.entropy;
            if zi > 0.0 {
                s_ideal -= R * zi * zi.ln();
            }
        }
        Ok(s_ideal + s_res)
    }

    /// Natural logarithm of the component fugacity coefficients.
    pub fn ln_fugacity_coefficients(&self, t: f64, p: f64, phase: Phase) -> Result<Vec<f64>> {
        let (tr, pr, coeffs_simp, coeffs_ref, vr_simp, vr_ref) =
            self.solve_reduced_volumes(t, p, phase)?;

        let (ln_phi_m_simp, _, _) = ln_phi_m(t, p, tr, pr, vr_simp, self.moles, &coeffs_simp, 0);
        let (ln_phi_m_ref, _, _) = ln_phi_m(t, p, tr, pr, vr_ref, self.moles, &coeffs_ref, 1);

        let z_simp = pr * vr_simp / tr;
        let z_ref = pr * vr_ref / tr;
        let lg_z_simp = z_simp.ln();
        let lg_z_ref = z_ref.ln();

        let n = self.model.components.len();
        let mut ln_phi = Vec::with_capacity(n);

        let moles_inv = 1.0 / self.moles;
        let tc_m_inv = 1.0 / self.tc_m;
        let pc_m_inv = 1.0 / self.pc_m;

        for i in 0..n {
            let tr_i = -tr * tc_m_inv * self.dtc_dni[i];

            let vr_i_simp =
                vr_simp * (pc_m_inv * self.dpc_dni[i] - tc_m_inv * self.dtc_dni[i] - moles_inv);
            let f_diff_simp = f_diff_ni(self.moles, tr, vr_simp, &coeffs_simp, 0, tr_i, vr_i_simp);
            let ln_phi_simp_i = f_diff_simp - lg_z_simp;

            let vr_i_ref =
                vr_ref * (pc_m_inv * self.dpc_dni[i] - tc_m_inv * self.dtc_dni[i] - moles_inv);
            let f_diff_ref = f_diff_ni(self.moles, tr, vr_ref, &coeffs_ref, 1, tr_i, vr_i_ref);
            let ln_phi_ref_i = f_diff_ref - lg_z_ref;

            let val = ln_phi_simp_i
                + (self.omega_m / W_REF) * (ln_phi_ref_i - ln_phi_simp_i)
                + (self.domega_dni[i] / W_REF) * (ln_phi_m_ref - ln_phi_m_simp);
            if !val.is_finite() {
                return Err(Error::InvalidInput(
                    "non-finite Lee-Kesler fugacity coefficient",
                ));
            }
            ln_phi.push(val);
        }

        Ok(ln_phi)
    }

    /// Component fugacity coefficients phi_i = exp(ln_phi_i).
    pub fn fugacity_coefficients(&self, t: f64, p: f64, phase: Phase) -> Result<Vec<f64>> {
        let ln_phi = self.ln_fugacity_coefficients(t, p, phase)?;
        Ok(ln_phi.into_iter().map(f64::exp).collect())
    }

    /// Residual departure Gibbs energy G_dep [J/mol].
    pub fn gibbs_energy_departure(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        let (tr, pr, coeffs_simp, coeffs_ref, vr_simp, vr_ref) =
            self.solve_reduced_volumes(t, p, phase)?;
        self.gibbs_dep_from_volumes(t, p, tr, pr, vr_simp, vr_ref, &coeffs_simp, &coeffs_ref)
    }

    /// Precomputed second composition derivatives of mixture critical volume `vc_m` (row-major `nc * nc`).
    pub fn d2vc_dn(&self) -> &[f64] {
        &self.d2vc_dn
    }

    /// Precomputed second composition derivatives of mixture critical temperature `tc_m` (row-major `nc * nc`).
    pub fn d2tc_dn(&self) -> &[f64] {
        &self.d2tc_dn
    }

    /// Precomputed second composition derivatives of mixture critical pressure `pc_m` (row-major `nc * nc`).
    pub fn d2pc_dn(&self) -> &[f64] {
        &self.d2pc_dn
    }

    /// Precomputed second composition derivatives of mixture acentric factor `omega_m` (row-major `nc * nc`).
    pub fn d2omega_dn(&self) -> &[f64] {
        &self.d2omega_dn
    }

    /// Precomputed second composition derivatives of mixture critical compressibility `zc_m` (row-major `nc * nc`).
    pub fn d2zc_dn(&self) -> &[f64] {
        &self.d2zc_dn
    }

    fn solve_reduced_states(
        &self,
        t: f64,
        p: f64,
        phase: Phase,
    ) -> Result<(f64, f64, FluidReducedState, FluidReducedState)> {
        let (tr, pr, coeffs_simp, coeffs_ref, vr_simp, vr_ref) =
            self.solve_reduced_volumes(t, p, phase)?;

        let f0 = evaluate_reduced_fluid(tr, pr, vr_simp, &coeffs_simp, 0)?;
        let fr = evaluate_reduced_fluid(tr, pr, vr_ref, &coeffs_ref, 1)?;
        Ok((tr, pr, f0, fr))
    }

    #[inline]
    fn interp(&self, v0: f64, vr: f64) -> f64 {
        v0 + (self.omega_m / W_REF) * (vr - v0)
    }

    #[inline]
    fn d_omega(v0: f64, vr: f64) -> f64 {
        (vr - v0) / W_REF
    }

    #[inline]
    fn z_derivatives_from_reduced(
        &self,
        tr: f64,
        pr: f64,
        f0: &FluidReducedState,
        fr: &FluidReducedState,
    ) -> LeeKeslerDerivatives {
        let n = self.model.components.len();
        let tc_inv = 1.0 / self.tc_m;
        let pc_inv = 1.0 / self.pc_m;
        let z_val = self.interp(f0.z, fr.z);
        let dz_dtr = self.interp(f0.dz_dtr, fr.dz_dtr);
        let dz_dpr = self.interp(f0.dz_dpr, fr.dz_dpr);
        let dz_dw = Self::d_omega(f0.z, fr.z);
        let dz_dt = dz_dtr * tc_inv;
        let dz_dp = dz_dpr * pc_inv;
        let mut dz_dn = Vec::with_capacity(n);
        for i in 0..n {
            let dtr_i = -tr * tc_inv * self.dtc_dni[i];
            let dpr_i = -pr * pc_inv * self.dpc_dni[i];
            dz_dn.push(dz_dtr * dtr_i + dz_dpr * dpr_i + dz_dw * self.domega_dni[i]);
        }
        LeeKeslerDerivatives {
            value: z_val,
            dt: dz_dt,
            dp: dz_dp,
            dn: dz_dn,
        }
    }

    #[inline]
    fn v_derivatives_from_z(
        &self,
        t: f64,
        p: f64,
        z: &LeeKeslerDerivatives,
    ) -> LeeKeslerDerivatives {
        let rt_p = R * t / p;
        let v_val = z.value * rt_p;
        let dv_dt = (z.dt * t + z.value) * R / p;
        let dv_dp = (z.dp - z.value / p) * rt_p;
        let dv_dn =
            z.dn.iter()
                .map(|&dzi| (self.moles * dzi + z.value) * rt_p)
                .collect();
        LeeKeslerDerivatives {
            value: v_val,
            dt: dv_dt,
            dp: dv_dp,
            dn: dv_dn,
        }
    }

    #[inline]
    fn h_res_derivatives_from_reduced(
        &self,
        tr: f64,
        pr: f64,
        f0: &FluidReducedState,
        fr: &FluidReducedState,
    ) -> LeeKeslerDerivatives {
        let n = self.model.components.len();
        let tc_inv = 1.0 / self.tc_m;
        let pc_inv = 1.0 / self.pc_m;
        let hr = self.interp(f0.hr, fr.hr);
        let dhr_dtr = self.interp(f0.dhr_dtr, fr.dhr_dtr);
        let dhr_dpr = self.interp(f0.dhr_dpr, fr.dhr_dpr);
        let dhr_dw = Self::d_omega(f0.hr, fr.hr);

        let h_res_val = R * self.tc_m * hr;
        let dh_res_dt = R * dhr_dtr;
        let dh_res_dp = R * self.tc_m * pc_inv * dhr_dpr;
        let mut dh_res_dn = Vec::with_capacity(n);
        for i in 0..n {
            let dtr_i = -tr * tc_inv * self.dtc_dni[i];
            let dpr_i = -pr * pc_inv * self.dpc_dni[i];
            let dhr_i = dhr_dtr * dtr_i + dhr_dpr * dpr_i + dhr_dw * self.domega_dni[i];
            dh_res_dn.push(h_res_val + self.moles * R * (self.dtc_dni[i] * hr + self.tc_m * dhr_i));
        }
        LeeKeslerDerivatives {
            value: h_res_val,
            dt: dh_res_dt,
            dp: dh_res_dp,
            dn: dh_res_dn,
        }
    }

    #[inline]
    fn s_res_derivatives_from_reduced(
        &self,
        tr: f64,
        pr: f64,
        f0: &FluidReducedState,
        fr: &FluidReducedState,
    ) -> LeeKeslerDerivatives {
        let n = self.model.components.len();
        let tc_inv = 1.0 / self.tc_m;
        let pc_inv = 1.0 / self.pc_m;
        let sr_val = self.interp(f0.sr, fr.sr);
        let dsr_dtr = self.interp(f0.dsr_dtr, fr.dsr_dtr);
        let dsr_dpr = self.interp(f0.dsr_dpr, fr.dsr_dpr);
        let dsr_dw = Self::d_omega(f0.sr, fr.sr);

        let s_res_val = R * sr_val;
        let ds_res_dt = R * dsr_dtr * tc_inv;
        let ds_res_dp = R * dsr_dpr * pc_inv;
        let mut ds_res_dn = Vec::with_capacity(n);
        for i in 0..n {
            let dtr_i = -tr * tc_inv * self.dtc_dni[i];
            let dpr_i = -pr * pc_inv * self.dpc_dni[i];
            let dsr_i = dsr_dtr * dtr_i + dsr_dpr * dpr_i + dsr_dw * self.domega_dni[i];
            ds_res_dn.push(s_res_val + self.moles * R * dsr_i);
        }
        LeeKeslerDerivatives {
            value: s_res_val,
            dt: ds_res_dt,
            dp: ds_res_dp,
            dn: ds_res_dn,
        }
    }

    #[inline]
    fn g_res_derivatives_from_reduced(
        &self,
        t: f64,
        p: f64,
        tr: f64,
        pr: f64,
        f0: &FluidReducedState,
        fr: &FluidReducedState,
    ) -> LeeKeslerDerivatives {
        let n = self.model.components.len();
        let tc_inv = 1.0 / self.tc_m;
        let pc_inv = 1.0 / self.pc_m;
        let rt = R * t;
        let rt_p = rt / p;
        let z_val = self.interp(f0.z, fr.z);
        let s_res_val = R * self.interp(f0.sr, fr.sr);
        let gr = self.interp(f0.gr, fr.gr);
        let g_t = self.interp(f0.dgr_dtr, fr.dgr_dtr);
        let g_p = self.interp(f0.dgr_dpr, fr.dgr_dpr);
        let g_w = Self::d_omega(f0.gr, fr.gr);
        let mut dg_res_dn = Vec::with_capacity(n);
        for i in 0..n {
            let dtr_i = -tr * tc_inv * self.dtc_dni[i];
            let dpr_i = -pr * pc_inv * self.dpc_dni[i];
            let di_gr = g_t * dtr_i + g_p * dpr_i + g_w * self.domega_dni[i];
            dg_res_dn.push(rt * (gr + self.moles * di_gr));
        }
        LeeKeslerDerivatives {
            value: rt * gr,
            dt: -s_res_val,
            dp: z_val * rt_p - rt_p,
            dn: dg_res_dn,
        }
    }

    #[allow(clippy::too_many_arguments)]
    #[inline]
    fn stable_speed_of_sound(
        v_val: f64,
        dv_dt: f64,
        dv_dp: f64,
        ds_tot_dt: f64,
        ds_tot_dp: f64,
        cp: f64,
        cv: f64,
        mw_kg: f64,
    ) -> Option<f64> {
        if !dv_dp.is_finite()
            || dv_dp >= 0.0
            || !cp.is_finite()
            || cp <= 0.0
            || !cv.is_finite()
            || cv <= 0.0
            || !ds_tot_dt.is_finite()
            || ds_tot_dt <= 0.0
            || !mw_kg.is_finite()
            || mw_kg <= 0.0
        {
            return None;
        }
        let dtdp_s = -ds_tot_dp / ds_tot_dt;
        let dvdp_s = dv_dp + dv_dt * dtdp_s;
        if !dvdp_s.is_finite() || dvdp_s >= 0.0 {
            return None;
        }
        let w2 = -(v_val * v_val) / (mw_kg * dvdp_s);
        if w2.is_finite() && w2 > 0.0 {
            Some(w2.sqrt())
        } else {
            None
        }
    }

    /// Single-solve evaluation of all TP properties, analytical derivatives,
    /// heat capacities, and speed of sound.
    pub fn state_tp(&self, t: f64, p: f64, phase: Phase) -> Result<LeeKeslerState> {
        let (tr, pr, f0, fr) = self.solve_reduced_states(t, p, phase)?;
        let n = self.model.components.len();
        let tc_inv = 1.0 / self.tc_m;
        let pc_inv = 1.0 / self.pc_m;
        let moles = self.moles;

        let mut dtr_dn = Vec::with_capacity(n);
        let mut dpr_dn = Vec::with_capacity(n);
        for i in 0..n {
            dtr_dn.push(-tr * tc_inv * self.dtc_dni[i]);
            dpr_dn.push(-pr * pc_inv * self.dpc_dni[i]);
        }

        // 1-4. Z, v, residual H, residual S
        let compressibility_factor = self.z_derivatives_from_reduced(tr, pr, &f0, &fr);
        let specific_volume = self.v_derivatives_from_z(t, p, &compressibility_factor);
        let residual_enthalpy = self.h_res_derivatives_from_reduced(tr, pr, &f0, &fr);
        let residual_entropy = self.s_res_derivatives_from_reduced(tr, pr, &f0, &fr);

        let z_val = compressibility_factor.value;
        let v_val = specific_volume.value;
        let dv_dt = specific_volume.dt;
        let dv_dp = specific_volume.dp;
        let h_res_val = residual_enthalpy.value;
        let dh_res_dt = residual_enthalpy.dt;
        let dh_res_dp = residual_enthalpy.dp;
        let s_res_val = residual_entropy.value;
        let ds_res_dt = residual_entropy.dt;
        let ds_res_dp = residual_entropy.dp;

        // 5. Residual Gibbs energy g_res [J/mol] and ln(phi_i) with all derivatives
        let gr = self.interp(f0.gr, fr.gr);
        let g_t = self.interp(f0.dgr_dtr, fr.dgr_dtr);
        let g_p = self.interp(f0.dgr_dpr, fr.dgr_dpr);
        let g_w = Self::d_omega(f0.gr, fr.gr);
        let g_tt = self.interp(f0.d2gr_dtr2, fr.d2gr_dtr2);
        let g_tp = self.interp(f0.d2gr_dtr_dpr, fr.d2gr_dtr_dpr);
        let g_pp = self.interp(f0.d2gr_dpr2, fr.d2gr_dpr2);
        let g_tw = Self::d_omega(f0.dgr_dtr, fr.dgr_dtr);
        let g_pw = Self::d_omega(f0.dgr_dpr, fr.dgr_dpr);

        let rt = R * t;
        let rt_p = rt / p;
        let mut d_gr = Vec::with_capacity(n);
        let mut ln_phi = Vec::with_capacity(n);
        let mut dlnphi_dt = Vec::with_capacity(n);
        let mut dlnphi_dp = Vec::with_capacity(n);
        let mut dg_res_dn = Vec::with_capacity(n);

        for i in 0..n {
            let di_gr = g_t * dtr_dn[i] + g_p * dpr_dn[i] + g_w * self.domega_dni[i];
            d_gr.push(di_gr);
            let lnp = gr + moles * di_gr;
            ln_phi.push(lnp);
            dg_res_dn.push(rt * lnp);
            dlnphi_dt.push(-residual_enthalpy.dn[i] / (rt * t));
            dlnphi_dp.push((moles * compressibility_factor.dn[i] + z_val - 1.0) / p);
        }

        let mut dlnphi_dn = vec![0.0; n * n];
        for i in 0..n {
            for j in i..n {
                let idx = i * n + j;
                let d2tr = (2.0 * tr * self.dtc_dni[i] * self.dtc_dni[j] * tc_inv
                    - tr * self.d2tc_dn[idx])
                    * tc_inv;
                let d2pr = (2.0 * pr * self.dpc_dni[i] * self.dpc_dni[j] * pc_inv
                    - pr * self.d2pc_dn[idx])
                    * pc_inv;
                let d2w = self.d2omega_dn[idx];

                let dg_t_j = g_tt * dtr_dn[j] + g_tp * dpr_dn[j] + g_tw * self.domega_dni[j];
                let dg_p_j = g_tp * dtr_dn[j] + g_pp * dpr_dn[j] + g_pw * self.domega_dni[j];
                let dg_w_j = g_tw * dtr_dn[j] + g_pw * dpr_dn[j];

                let val = d_gr[i]
                    + d_gr[j]
                    + moles
                        * (dg_t_j * dtr_dn[i]
                            + dg_p_j * dpr_dn[i]
                            + dg_w_j * self.domega_dni[i]
                            + g_t * d2tr
                            + g_p * d2pr
                            + g_w * d2w);
                dlnphi_dn[idx] = val;
                dlnphi_dn[j * n + i] = val;
            }
        }

        let residual_gibbs = LeeKeslerDerivatives {
            value: rt * gr,
            dt: -s_res_val,
            dp: v_val - rt_p,
            dn: dg_res_dn,
        };
        let ln_fugacity = LeeKeslerFugacityDerivatives {
            ln_phi,
            dt: dlnphi_dt,
            dp: dlnphi_dp,
            dn: dlnphi_dn,
        };

        // 6. Ideal gas mixture contributions + total H, S, G, Cp, Cv, sound speed
        let mut h_id = 0.0;
        let mut dh_id_dt = 0.0;
        let mut s_id = 0.0;
        let mut mw_kg = 0.0;
        let mut dh_tot_dn = residual_enthalpy.dn.clone();
        let mut ds_tot_dn = residual_entropy.dn.clone();
        let mut dg_tot_dn = residual_gibbs.dn.clone();

        for i in 0..n {
            let xi = self.z[i] / moles;
            let props = self.model.ideals[i].properties(t, p)?;
            let s_bar_id = props.entropy - R * xi.max(1.0e-100).ln();

            h_id += xi * props.enthalpy;
            dh_id_dt += xi * props.enthalpy_dt;
            if xi > 0.0 {
                s_id += xi * (props.entropy - R * xi.ln());
            }
            mw_kg += xi * self.model.components[i].molar_mass * 1.0e-3;

            dh_tot_dn[i] += props.enthalpy;
            ds_tot_dn[i] += s_bar_id;
            dg_tot_dn[i] += props.enthalpy - t * s_bar_id;
        }

        let cp = dh_id_dt + dh_res_dt;
        let cv = cp + t * (dv_dt * dv_dt) / dv_dp;

        let enthalpy = LeeKeslerDerivatives {
            value: h_id + h_res_val,
            dt: cp,
            dp: dh_res_dp,
            dn: dh_tot_dn,
        };
        let ds_tot_dt = dh_id_dt / t + ds_res_dt;
        let ds_tot_dp = ds_res_dp - R / p;
        let entropy = LeeKeslerDerivatives {
            value: s_id + s_res_val,
            dt: ds_tot_dt,
            dp: ds_tot_dp,
            dn: ds_tot_dn,
        };
        let gibbs = LeeKeslerDerivatives {
            value: (h_id - t * s_id) + residual_gibbs.value,
            dt: -(s_id + s_res_val),
            dp: v_val,
            dn: dg_tot_dn,
        };

        let speed_of_sound =
            Self::stable_speed_of_sound(v_val, dv_dt, dv_dp, ds_tot_dt, ds_tot_dp, cp, cv, mw_kg);

        Ok(LeeKeslerState {
            compressibility_factor,
            specific_volume,
            residual_enthalpy,
            enthalpy,
            residual_entropy,
            entropy,
            residual_gibbs,
            gibbs,
            ln_fugacity,
            heat_capacity_p: cp,
            heat_capacity_v: cv,
            speed_of_sound,
        })
    }

    /// Compressibility factor `Z` with analytical `(T, P, n)` derivatives.
    pub fn compressibility_factor_derivatives(
        &self,
        t: f64,
        p: f64,
        phase: Phase,
    ) -> Result<LeeKeslerDerivatives> {
        let (tr, pr, f0, fr) = self.solve_reduced_states(t, p, phase)?;
        Ok(self.z_derivatives_from_reduced(tr, pr, &f0, &fr))
    }

    /// Molar specific volume `v` `[m^3/mol]` with analytical `(T, P, n)` derivatives.
    pub fn specific_volume_derivatives(
        &self,
        t: f64,
        p: f64,
        phase: Phase,
    ) -> Result<LeeKeslerDerivatives> {
        let (tr, pr, f0, fr) = self.solve_reduced_states(t, p, phase)?;
        let z = self.z_derivatives_from_reduced(tr, pr, &f0, &fr);
        Ok(self.v_derivatives_from_z(t, p, &z))
    }

    /// Molar residual enthalpy `h_res` `[J/mol]` with analytical `(T, P, n)` derivatives.
    pub fn residual_enthalpy_derivatives(
        &self,
        t: f64,
        p: f64,
        phase: Phase,
    ) -> Result<LeeKeslerDerivatives> {
        let (tr, pr, f0, fr) = self.solve_reduced_states(t, p, phase)?;
        Ok(self.h_res_derivatives_from_reduced(tr, pr, &f0, &fr))
    }

    /// Molar total enthalpy `h` `[J/mol]` with analytical `(T, P, n)` derivatives.
    pub fn enthalpy_derivatives(
        &self,
        t: f64,
        p: f64,
        phase: Phase,
    ) -> Result<LeeKeslerDerivatives> {
        let (tr, pr, f0, fr) = self.solve_reduced_states(t, p, phase)?;
        let mut h = self.h_res_derivatives_from_reduced(tr, pr, &f0, &fr);
        for i in 0..self.model.components.len() {
            let xi = self.z[i] / self.moles;
            let props = self.model.ideals[i].properties(t, p)?;
            h.value += xi * props.enthalpy;
            h.dt += xi * props.enthalpy_dt;
            h.dn[i] += props.enthalpy;
        }
        Ok(h)
    }

    /// Molar residual entropy `s_res` `[J/(mol K)]` with analytical `(T, P, n)` derivatives.
    pub fn residual_entropy_derivatives(
        &self,
        t: f64,
        p: f64,
        phase: Phase,
    ) -> Result<LeeKeslerDerivatives> {
        let (tr, pr, f0, fr) = self.solve_reduced_states(t, p, phase)?;
        Ok(self.s_res_derivatives_from_reduced(tr, pr, &f0, &fr))
    }

    /// Molar total entropy `s` `[J/(mol K)]` with analytical `(T, P, n)` derivatives.
    ///
    /// For absent or trace components (`x_i < 1e-100`), the ideal-gas mixing term
    /// in the partial molar entropy `dn[i]` uses ThermoPack's `max(x_i, 1e-100)`
    /// log regularization so `dn[i]` remains finite while the mixture value uses
    /// the exact `0 * ln(0) = 0` limit.
    pub fn entropy_derivatives(
        &self,
        t: f64,
        p: f64,
        phase: Phase,
    ) -> Result<LeeKeslerDerivatives> {
        let (tr, pr, f0, fr) = self.solve_reduced_states(t, p, phase)?;
        let mut s = self.s_res_derivatives_from_reduced(tr, pr, &f0, &fr);
        s.dp -= R / p;
        for i in 0..self.model.components.len() {
            let xi = self.z[i] / self.moles;
            let props = self.model.ideals[i].properties(t, p)?;
            if xi > 0.0 {
                s.value += xi * (props.entropy - R * xi.ln());
            }
            s.dt += xi * props.enthalpy_dt / t;
            s.dn[i] += props.entropy - R * xi.max(1.0e-100).ln();
        }
        Ok(s)
    }

    /// Molar residual Gibbs energy `g_res` `[J/mol]` with analytical `(T, P, n)` derivatives.
    pub fn residual_gibbs_derivatives(
        &self,
        t: f64,
        p: f64,
        phase: Phase,
    ) -> Result<LeeKeslerDerivatives> {
        let (tr, pr, f0, fr) = self.solve_reduced_states(t, p, phase)?;
        Ok(self.g_res_derivatives_from_reduced(t, p, tr, pr, &f0, &fr))
    }

    /// Molar total Gibbs energy `g` `[J/mol]` with analytical `(T, P, n)` derivatives.
    ///
    /// For absent or trace components (`x_i < 1e-100`), the ideal-gas mixing term
    /// in the chemical potential `dn[i]` uses ThermoPack's `max(x_i, 1e-100)`
    /// log regularization so `dn[i]` remains finite while the mixture value uses
    /// the exact `0 * ln(0) = 0` limit.
    pub fn gibbs_derivatives(&self, t: f64, p: f64, phase: Phase) -> Result<LeeKeslerDerivatives> {
        let (tr, pr, f0, fr) = self.solve_reduced_states(t, p, phase)?;
        let mut g = self.g_res_derivatives_from_reduced(t, p, tr, pr, &f0, &fr);
        g.dp += R * t / p;
        let mut h_id = 0.0;
        let mut s_id = 0.0;
        for i in 0..self.model.components.len() {
            let xi = self.z[i] / self.moles;
            let props = self.model.ideals[i].properties(t, p)?;
            let s_bar_id = props.entropy - R * xi.max(1.0e-100).ln();
            h_id += xi * props.enthalpy;
            if xi > 0.0 {
                s_id += xi * (props.entropy - R * xi.ln());
            }
            g.dn[i] += props.enthalpy - t * s_bar_id;
        }
        g.value += h_id - t * s_id;
        g.dt -= s_id;
        Ok(g)
    }

    /// Log fugacity coefficients with analytical `(T, P, n)` derivatives.
    pub fn ln_fugacity_derivatives(
        &self,
        t: f64,
        p: f64,
        phase: Phase,
    ) -> Result<LeeKeslerFugacityDerivatives> {
        let (tr, pr, f0, fr) = self.solve_reduced_states(t, p, phase)?;
        let n = self.model.components.len();
        let tc_inv = 1.0 / self.tc_m;
        let pc_inv = 1.0 / self.pc_m;
        let moles = self.moles;

        let z_val = self.interp(f0.z, fr.z);
        let dz_dtr = self.interp(f0.dz_dtr, fr.dz_dtr);
        let dz_dpr = self.interp(f0.dz_dpr, fr.dz_dpr);
        let dz_dw = Self::d_omega(f0.z, fr.z);

        let hr = self.interp(f0.hr, fr.hr);
        let dhr_dtr = self.interp(f0.dhr_dtr, fr.dhr_dtr);
        let dhr_dpr = self.interp(f0.dhr_dpr, fr.dhr_dpr);
        let dhr_dw = Self::d_omega(f0.hr, fr.hr);
        let h_res_val = R * self.tc_m * hr;

        let gr = self.interp(f0.gr, fr.gr);
        let g_t = self.interp(f0.dgr_dtr, fr.dgr_dtr);
        let g_p = self.interp(f0.dgr_dpr, fr.dgr_dpr);
        let g_w = Self::d_omega(f0.gr, fr.gr);
        let g_tt = self.interp(f0.d2gr_dtr2, fr.d2gr_dtr2);
        let g_tp = self.interp(f0.d2gr_dtr_dpr, fr.d2gr_dtr_dpr);
        let g_pp = self.interp(f0.d2gr_dpr2, fr.d2gr_dpr2);
        let g_tw = Self::d_omega(f0.dgr_dtr, fr.dgr_dtr);
        let g_pw = Self::d_omega(f0.dgr_dpr, fr.dgr_dpr);

        let rt2 = R * t * t;
        let mut dtr_dn = Vec::with_capacity(n);
        let mut dpr_dn = Vec::with_capacity(n);
        let mut d_gr = Vec::with_capacity(n);
        let mut ln_phi = Vec::with_capacity(n);
        let mut dt = Vec::with_capacity(n);
        let mut dp = Vec::with_capacity(n);

        for i in 0..n {
            let dtr_i = -tr * tc_inv * self.dtc_dni[i];
            let dpr_i = -pr * pc_inv * self.dpc_dni[i];
            dtr_dn.push(dtr_i);
            dpr_dn.push(dpr_i);

            let di_gr = g_t * dtr_i + g_p * dpr_i + g_w * self.domega_dni[i];
            d_gr.push(di_gr);
            ln_phi.push(gr + moles * di_gr);

            let dhr_i = dhr_dtr * dtr_i + dhr_dpr * dpr_i + dhr_dw * self.domega_dni[i];
            let dh_res_dni = h_res_val + moles * R * (self.dtc_dni[i] * hr + self.tc_m * dhr_i);
            dt.push(-dh_res_dni / rt2);

            let dz_i = dz_dtr * dtr_i + dz_dpr * dpr_i + dz_dw * self.domega_dni[i];
            dp.push((moles * dz_i + z_val - 1.0) / p);
        }

        let mut dn = vec![0.0; n * n];
        for i in 0..n {
            for j in i..n {
                let idx = i * n + j;
                let d2tr = (2.0 * tr * self.dtc_dni[i] * self.dtc_dni[j] * tc_inv
                    - tr * self.d2tc_dn[idx])
                    * tc_inv;
                let d2pr = (2.0 * pr * self.dpc_dni[i] * self.dpc_dni[j] * pc_inv
                    - pr * self.d2pc_dn[idx])
                    * pc_inv;
                let d2w = self.d2omega_dn[idx];

                let dg_t_j = g_tt * dtr_dn[j] + g_tp * dpr_dn[j] + g_tw * self.domega_dni[j];
                let dg_p_j = g_tp * dtr_dn[j] + g_pp * dpr_dn[j] + g_pw * self.domega_dni[j];
                let dg_w_j = g_tw * dtr_dn[j] + g_pw * dpr_dn[j];

                let val = d_gr[i]
                    + d_gr[j]
                    + moles
                        * (dg_t_j * dtr_dn[i]
                            + dg_p_j * dpr_dn[i]
                            + dg_w_j * self.domega_dni[i]
                            + g_t * d2tr
                            + g_p * d2pr
                            + g_w * d2w);
                dn[idx] = val;
                dn[j * n + i] = val;
            }
        }

        Ok(LeeKeslerFugacityDerivatives { ln_phi, dt, dp, dn })
    }

    /// Molar isobaric heat capacity `Cp` `[J/(mol K)]`.
    pub fn heat_capacity_p(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        let (_, _, f0, fr) = self.solve_reduced_states(t, p, phase)?;
        let dh_res_dt = R * self.interp(f0.dhr_dtr, fr.dhr_dtr);
        let mut dh_id_dt = 0.0;
        for i in 0..self.model.components.len() {
            dh_id_dt +=
                (self.z[i] / self.moles) * self.model.ideals[i].properties(t, p)?.enthalpy_dt;
        }
        Ok(dh_id_dt + dh_res_dt)
    }

    /// Molar isochoric heat capacity `Cv` `[J/(mol K)]`.
    pub fn heat_capacity_v(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        let (_, _, f0, fr) = self.solve_reduced_states(t, p, phase)?;
        let dh_res_dt = R * self.interp(f0.dhr_dtr, fr.dhr_dtr);
        let mut dh_id_dt = 0.0;
        for i in 0..self.model.components.len() {
            dh_id_dt +=
                (self.z[i] / self.moles) * self.model.ideals[i].properties(t, p)?.enthalpy_dt;
        }
        let cp = dh_id_dt + dh_res_dt;
        let z_val = self.interp(f0.z, fr.z);
        let dz_dt = self.interp(f0.dz_dtr, fr.dz_dtr) / self.tc_m;
        let dz_dp = self.interp(f0.dz_dpr, fr.dz_dpr) / self.pc_m;
        let rt_p = R * t / p;
        let dv_dt = (dz_dt * t + z_val) * R / p;
        let dv_dp = (dz_dp - z_val / p) * rt_p;
        Ok(cp + t * (dv_dt * dv_dt) / dv_dp)
    }

    /// Single-phase speed of sound `[m/s]` matching `src/speed_of_sound.f90`.
    ///
    /// Requires a mechanically and thermally stable single-phase state
    /// (`(dv/dP)_T < 0`, `Cp > 0`, `Cv > 0`); returns [`Error::NoPhysicalRoot`]
    /// on unstable or singular states.
    pub fn speed_of_sound(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        let (_, _, f0, fr) = self.solve_reduced_states(t, p, phase)?;
        let tc_inv = 1.0 / self.tc_m;
        let pc_inv = 1.0 / self.pc_m;
        let z_val = self.interp(f0.z, fr.z);
        let dz_dt = self.interp(f0.dz_dtr, fr.dz_dtr) * tc_inv;
        let dz_dp = self.interp(f0.dz_dpr, fr.dz_dpr) * pc_inv;
        let rt_p = R * t / p;
        let v_val = z_val * rt_p;
        let dv_dt = (dz_dt * t + z_val) * R / p;
        let dv_dp = (dz_dp - z_val / p) * rt_p;

        let dh_res_dt = R * self.interp(f0.dhr_dtr, fr.dhr_dtr);
        let ds_res_dt = R * self.interp(f0.dsr_dtr, fr.dsr_dtr) * tc_inv;
        let ds_res_dp = R * self.interp(f0.dsr_dpr, fr.dsr_dpr) * pc_inv;

        let mut dh_id_dt = 0.0;
        let mut mw_kg = 0.0;
        for i in 0..self.model.components.len() {
            let xi = self.z[i] / self.moles;
            dh_id_dt += xi * self.model.ideals[i].properties(t, p)?.enthalpy_dt;
            mw_kg += xi * self.model.components[i].molar_mass * 1.0e-3;
        }
        let cp = dh_id_dt + dh_res_dt;
        let cv = cp + t * (dv_dt * dv_dt) / dv_dp;
        let ds_tot_dt = dh_id_dt / t + ds_res_dt;
        let ds_tot_dp = ds_res_dp - R / p;
        Self::stable_speed_of_sound(v_val, dv_dt, dv_dp, ds_tot_dt, ds_tot_dp, cp, cv, mw_kg)
            .ok_or(Error::NoPhysicalRoot)
    }
}
