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

#[allow(dead_code)]
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

#[allow(dead_code)]
#[inline]
fn f_diff2_tr(moles: f64, tr: f64, vr: f64, sr: usize) -> f64 {
    let d2 = BwrCoeffsDiff2::new(tr, sr);
    f_solver(moles, vr, d2.b_trtr, d2.c_trtr, d2.d_trtr, d2.e_trtr, sr)
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
) -> f64 {
    let tolerance = f64::EPSILON * 200.0;
    let max_iterations = 200;

    let mut z = z_init;
    for _ in 1..=max_iterations {
        let (f_val, f_diff, f_diff2) = fz_with_diff(pr, tr, z, coeffs, sr);
        let mut f = f_val;
        let mut fd = f_diff;
        let mut fd2 = f_diff2;
        if fd < 0.0 {
            z = if phase == 1 { z_max } else { z_min };
            let (f_val2, fd_val2, fd2_val2) = fz_with_diff(pr, tr, z, coeffs, sr);
            f = f_val2;
            fd = fd_val2;
            fd2 = fd2_val2;
        }
        let delta = 0.5 * f * fd2 / (fd * fd);
        let dz = -(f / fd) * (1.0 + delta);
        let z_temp = z + dz;

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

    z * tr / pr
}

fn calc_reduced_volume(
    tr: f64,
    pr: f64,
    phase: usize,
    coeffs_simp: &BwrCoeffs,
    coeffs_ref: &BwrCoeffs,
) -> (f64, f64) {
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
    );
    let vr_ref = z_newt_raps(
        tr, pr, phase_ref, 1, z_ref, z_max_ref, z_min_ref, coeffs_ref,
    );

    (vr_simp, vr_ref)
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

        // Pair critical volume and temperature
        let mut vc_m = 0.0;
        let mut tc_m_sum = 0.0;
        let mut omega_m = 0.0;

        for j in 0..n {
            let tcj = self.components[j].critical_temperature;
            let vcj = vc[j];
            for k in 0..n {
                let tck = self.components[k].critical_temperature;
                let vck = vc[k];
                let vc_jk = 0.125 * (vcj.cbrt() + vck.cbrt()).powi(3);
                let tc_jk = self.kij[j * n + k] * (tcj * tck).sqrt();

                let n_prod = z[j] * z[k];
                vc_m += n_prod * vc_jk;
                tc_m_sum += n_prod * tc_jk * vc_jk.powf(ETA);
            }
            omega_m += z[j] * self.components[j].acentric_factor;
        }

        vc_m *= moles_inv * moles_inv;
        let vc_m_inv = 1.0 / vc_m;

        let tc_m = tc_m_sum * moles_inv * moles_inv * vc_m_inv.powf(ETA);
        omega_m *= moles_inv;
        let zc_m = 0.2905 - 0.085 * omega_m;
        let pc_m = R * zc_m * tc_m * vc_m_inv;

        // Compositional derivatives
        let mut dvc_dni = Vec::with_capacity(n);
        let mut dtc_dni = Vec::with_capacity(n);
        let mut domega_dni = Vec::with_capacity(n);
        let mut dzc_dni = Vec::with_capacity(n);
        let mut dpc_dni = Vec::with_capacity(n);

        for i in 0..n {
            let tci = self.components[i].critical_temperature;
            let vci = vc[i];

            let mut v1 = 0.0;
            let mut comp_sum = 0.0;
            for l in 0..n {
                let tcl = self.components[l].critical_temperature;
                let vcl = vc[l];
                let vc_il = 0.125 * (vci.cbrt() + vcl.cbrt()).powi(3);
                let tc_il = self.kij[i * n + l] * (tci * tcl).sqrt();

                v1 += z[l] * vc_il;
                comp_sum += z[l] * vc_il.powf(ETA) * tc_il;
            }

            let dvc = 2.0 * moles_inv * (moles_inv * v1 - vc_m);
            let t1 = -ETA * tc_m * vc_m_inv * dvc;
            let t2 = 2.0 * comp_sum * vc_m_inv.powf(ETA) * moles_inv * moles_inv;
            let t3 = -2.0 * tc_m * moles_inv;
            let dtc = t1 + t2 + t3;

            let domega = moles_inv * (self.components[i].acentric_factor - omega_m);
            let dzc = -0.085 * domega;
            let dpc = pc_m * ((1.0 / zc_m) * dzc + (1.0 / tc_m) * dtc - vc_m_inv * dvc);

            dvc_dni.push(dvc);
            dtc_dni.push(dtc);
            domega_dni.push(domega);
            dzc_dni.push(dzc);
            dpc_dni.push(dpc);
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
        })
    }

    /// Direct evaluation of compressibility factor.
    pub fn compressibility_factor(&self, t: f64, p: f64, z: &[f64], phase: Phase) -> Result<f64> {
        self.prepare(z)?.compressibility_factor(t, p, phase)
    }

    /// Direct evaluation of specific volume.
    pub fn specific_volume(&self, t: f64, p: f64, z: &[f64], phase: Phase) -> Result<f64> {
        self.prepare(z)?.specific_volume(t, p, phase)
    }

    /// Direct evaluation of total enthalpy.
    pub fn enthalpy(&self, t: f64, p: f64, z: &[f64], phase: Phase) -> Result<f64> {
        self.prepare(z)?.enthalpy(t, p, phase)
    }

    /// Direct evaluation of total entropy.
    pub fn entropy(&self, t: f64, p: f64, z: &[f64], phase: Phase) -> Result<f64> {
        self.prepare(z)?.entropy(t, p, phase)
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

    fn resolve_phase_int(&self, t: f64, p: f64, phase: Phase) -> Result<usize> {
        match phase {
            Phase::Liquid => Ok(1),
            Phase::Vapor => Ok(2),
            Phase::MinimumGibbs => {
                let g_liq = self.gibbs_energy_departure_for_int(t, p, 1);
                let g_vap = self.gibbs_energy_departure_for_int(t, p, 2);
                match (g_liq, g_vap) {
                    (Ok(gl), Ok(gv)) => Ok(if gl <= gv { 1 } else { 2 }),
                    (Ok(_), Err(_)) => Ok(1),
                    (Err(_), Ok(_)) => Ok(2),
                    (Err(e), Err(_)) => Err(e),
                }
            }
        }
    }

    fn gibbs_energy_departure_for_int(&self, t: f64, p: f64, phase_int: usize) -> Result<f64> {
        let tr = t / self.tc_m;
        let pr = p / self.pc_m;
        let coeffs_simp = BwrCoeffs::new(tr, 0);
        let coeffs_ref = BwrCoeffs::new(tr, 1);
        let (vr_simp, vr_ref) = calc_reduced_volume(tr, pr, phase_int, &coeffs_simp, &coeffs_ref);
        let (ln_phi_m_simp, _, _) = ln_phi_m(t, p, tr, pr, vr_simp, self.moles, &coeffs_simp, 0);
        let (ln_phi_m_ref, _, _) = ln_phi_m(t, p, tr, pr, vr_ref, self.moles, &coeffs_ref, 1);
        let g_dep =
            R * t * (ln_phi_m_simp + (self.omega_m / W_REF) * (ln_phi_m_ref - ln_phi_m_simp));
        Ok(g_dep / self.moles)
    }

    /// Compressibility factor Z.
    pub fn compressibility_factor(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        positive(t, "temperature")?;
        positive(p, "pressure")?;

        let phase_int = self.resolve_phase_int(t, p, phase)?;
        let tr = t / self.tc_m;
        let pr = p / self.pc_m;

        let coeffs_simp = BwrCoeffs::new(tr, 0);
        let coeffs_ref = BwrCoeffs::new(tr, 1);
        let (vr_simp, vr_ref) = calc_reduced_volume(tr, pr, phase_int, &coeffs_simp, &coeffs_ref);

        let z_simp = pr * vr_simp / tr;
        let z_ref = pr * vr_ref / tr;
        Ok(z_simp + (self.omega_m / W_REF) * (z_ref - z_simp))
    }

    /// Specific volume [m^3/mol].
    pub fn specific_volume(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        let z = self.compressibility_factor(t, p, phase)?;
        Ok(z * R * t / p)
    }

    /// Residual enthalpy [J/mol].
    pub fn residual_enthalpy(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        positive(t, "temperature")?;
        positive(p, "pressure")?;

        let phase_int = self.resolve_phase_int(t, p, phase)?;
        let tr = t / self.tc_m;
        let pr = p / self.pc_m;

        let coeffs_simp = BwrCoeffs::new(tr, 0);
        let coeffs_ref = BwrCoeffs::new(tr, 1);
        let (vr_simp, vr_ref) = calc_reduced_volume(tr, pr, phase_int, &coeffs_simp, &coeffs_ref);

        let z_simp = pr * vr_simp / tr;
        let z_ref = pr * vr_ref / tr;

        let h_simp =
            R * t * (-tr * f_diff_tr(self.moles, tr, vr_simp, 0) + self.moles * (z_simp - 1.0));
        let h_ref =
            R * t * (-tr * f_diff_tr(self.moles, tr, vr_ref, 1) + self.moles * (z_ref - 1.0));

        let h_res = h_simp + (self.omega_m / W_REF) * (h_ref - h_simp);
        Ok(h_res / self.moles)
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
        positive(t, "temperature")?;
        positive(p, "pressure")?;

        let phase_int = self.resolve_phase_int(t, p, phase)?;
        let tr = t / self.tc_m;
        let pr = p / self.pc_m;

        let coeffs_simp = BwrCoeffs::new(tr, 0);
        let coeffs_ref = BwrCoeffs::new(tr, 1);
        let (vr_simp, vr_ref) = calc_reduced_volume(tr, pr, phase_int, &coeffs_simp, &coeffs_ref);

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

        let s_res = s_simp + (self.omega_m / W_REF) * (s_ref - s_simp);
        Ok(s_res / self.moles)
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
            if zi > 1.0e-30 {
                s_ideal -= R * zi * zi.ln();
            }
        }
        Ok(s_ideal + s_res)
    }

    /// Natural logarithm of the component fugacity coefficients.
    pub fn ln_fugacity_coefficients(&self, t: f64, p: f64, phase: Phase) -> Result<Vec<f64>> {
        positive(t, "temperature")?;
        positive(p, "pressure")?;

        let phase_int = self.resolve_phase_int(t, p, phase)?;
        let tr = t / self.tc_m;
        let pr = p / self.pc_m;

        let coeffs_simp = BwrCoeffs::new(tr, 0);
        let coeffs_ref = BwrCoeffs::new(tr, 1);
        let (vr_simp, vr_ref) = calc_reduced_volume(tr, pr, phase_int, &coeffs_simp, &coeffs_ref);

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
        positive(t, "temperature")?;
        positive(p, "pressure")?;

        let phase_int = self.resolve_phase_int(t, p, phase)?;
        self.gibbs_energy_departure_for_int(t, p, phase_int)
    }
}
