//! Weeks-Chandler-Andersen (WCA) hard-sphere reference and cavity perturbation.
//! Ported from `src/hardsphere_wca.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! The WCA perturbation family splits the Lennard-Jones potential at the minimum
//! `r_split = 2^(1/6)` and evaluates a temperature-dependent effective
//! hard-sphere diameter, a Carnahan-Starling (or Carnahan-Starling-Kolafa)
//! hard-sphere reference, and a 61-point Gauss-Kronrod cavity soft-repulsion
//! integral together with all analytical first and second `(T, V, n)`
//! derivatives.

use crate::chain_rule::PackingFractionDerivatives;
use crate::error::{Error, Result, positive};
use crate::helmholtz::{DensityHelmholtz, PureHelmholtz};
use crate::quadrature::Quadrature;
use std::f64::consts::PI;
use std::sync::OnceLock;

/// The reduced WCA potential split separation, `2^(1/6)`.
const R_SPLIT: f64 = 1.122_462_048_309_373;

/// The packing fraction at which the cavity correlation switches to quadratic
/// Taylor extrapolation in `eta`.
const ETA_C: f64 = 0.52;

/// Number of Gauss-Kronrod quadrature nodes used by `calc_cavity_integral_LJ_Fres`.
const N_QUAD: usize = 61;

static KRONROD_61: OnceLock<([f64; N_QUAD], [f64; N_QUAD])> = OnceLock::new();

fn kronrod_61() -> &'static ([f64; N_QUAD], [f64; N_QUAD]) {
    KRONROD_61.get_or_init(|| {
        let pos = Quadrature::Kronrod61.positions();
        let wt = Quadrature::Kronrod61.weights();
        let mut x = [0.0; N_QUAD];
        let mut w = [0.0; N_QUAD];
        x.copy_from_slice(&pos);
        w.copy_from_slice(&wt);
        (x, w)
    })
}

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

/// Molar reduced cavity Helmholtz contribution `F_n(eta, T) = -12*eta*(I1 + I2)/d*^3`
/// and its partial derivatives with respect to packing fraction `eta` and
/// temperature `T`.
#[derive(Debug, Clone, Copy)]
pub struct CavityEtaDerivatives {
    /// Molar reduced cavity Helmholtz energy `F_n`.
    pub value: f64,
    /// First derivative `dF_n / d(eta)` at fixed `T`.
    pub eta: f64,
    /// Second derivative `d^2 F_n / d(eta)^2` at fixed `T`.
    pub eta_eta: f64,
    /// Explicit temperature derivative `dF_n / dT` at fixed `eta`.
    pub temperature: f64,
    /// Explicit second temperature derivative `d^2 F_n / dT^2` at fixed `eta`.
    pub temperature_second: f64,
    /// Mixed derivative `d^2 F_n / (d(eta) dT)`.
    pub eta_temperature: f64,
}

/// Precomputed temperature-only quadrature node data for the 122 nodes
/// (61 for `I1` on `[r_min, d*]` plus 61 for `I2` on `[d*, 2^(1/6)]`).
#[derive(Debug, Clone)]
pub(crate) struct CavityPreparedTemp {
    d_star: f64,
    d_star_t: f64,
    d_star_tt: f64,
    r_hs: [f64; 2 * N_QUAD],
    r_hs_t: [f64; 2 * N_QUAD],
    r_hs_tt: [f64; 2 * N_QUAD],
    weight: [f64; 2 * N_QUAD],
    weight_t: [f64; 2 * N_QUAD],
    weight_tt: [f64; 2 * N_QUAD],
}

impl CavityPreparedTemp {
    pub(crate) fn new(
        sigma: f64,
        eps_divk: f64,
        diameter: f64,
        diameter_t: f64,
        diameter_tt: f64,
        temperature: f64,
    ) -> Result<Self> {
        positive(sigma, "segment diameter")?;
        positive(eps_divk, "well depth")?;
        positive(diameter, "hard-sphere diameter")?;
        positive(temperature, "temperature")?;
        if !diameter_t.is_finite() || !diameter_tt.is_finite() {
            return Err(Error::InvalidInput(
                "finite hard-sphere diameter derivatives required",
            ));
        }

        let (x_vec, w_vec) = kronrod_61();
        let d_star = diameter / sigma;
        let d_star_t = diameter_t / sigma;
        let d_star_tt = diameter_tt / sigma;
        let (r_min0, r_min0_t, r_min0_tt) = minimum_cavity_radius(temperature, eps_divk);

        let mut prepared = Self {
            d_star,
            d_star_t,
            d_star_tt,
            r_hs: [0.0; 2 * N_QUAD],
            r_hs_t: [0.0; 2 * N_QUAD],
            r_hs_tt: [0.0; 2 * N_QUAD],
            weight: [0.0; 2 * N_QUAD],
            weight_t: [0.0; 2 * N_QUAD],
            weight_tt: [0.0; 2 * N_QUAD],
        };

        let inv_d = 1.0 / d_star;
        let inv_d2 = inv_d * inv_d;
        let inv_d3 = inv_d2 * inv_d;

        // Interval 1: [r_min, d_star], integrand factor e_lj_wca
        let r_diff_1 = 0.5 * (d_star - r_min0);
        let r_diff_1_t = 0.5 * (d_star_t - r_min0_t);
        let r_diff_1_tt = 0.5 * (d_star_tt - r_min0_tt);
        for i in 0..N_QUAD {
            let x = x_vec[i];
            let w = w_vec[i];
            let r = r_diff_1 * x + r_diff_1 + r_min0;
            let r_t = r_diff_1_t * x + r_diff_1_t + r_min0_t;
            let r_tt = r_diff_1_tt * x + r_diff_1_tt + r_min0_tt;
            let pref = r_diff_1 * w;
            let pref_t = r_diff_1_t * w;
            let pref_tt = r_diff_1_tt * w;
            let (e, e_t, e_tt) = lj_boltzmann_factor(r, r_t, r_tt, temperature, eps_divk);
            let r2 = r * r;
            let two_r_rt = 2.0 * r * r_t;
            prepared.r_hs[i] = r * inv_d;
            prepared.r_hs_t[i] = r_t * inv_d - d_star_t * r * inv_d2;
            prepared.r_hs_tt[i] =
                r_tt * inv_d - 2.0 * d_star_t * r_t * inv_d2 - d_star_tt * r * inv_d2
                    + 2.0 * d_star_t * d_star_t * r * inv_d3;
            prepared.weight[i] = pref * e * r2;
            prepared.weight_t[i] = pref_t * e * r2 + pref * e_t * r2 + pref * e * two_r_rt;
            prepared.weight_tt[i] = pref_tt * e * r2
                + 2.0 * pref_t * e_t * r2
                + 2.0 * pref_t * e * two_r_rt
                + pref * e_tt * r2
                + 2.0 * pref * e_t * two_r_rt
                + 2.0 * pref * e * (r_t * r_t + r * r_tt);
        }

        // Interval 2: [d_star, r_split], integrand factor (e_lj_wca - 1)
        let r_diff_2 = 0.5 * (R_SPLIT - d_star);
        let r_diff_2_t = -0.5 * d_star_t;
        let r_diff_2_tt = -0.5 * d_star_tt;
        for i in 0..N_QUAD {
            let x = x_vec[i];
            let w = w_vec[i];
            let r = r_diff_2 * x + r_diff_2 + d_star;
            let r_t = r_diff_2_t * x + r_diff_2_t + d_star_t;
            let r_tt = r_diff_2_tt * x + r_diff_2_tt + d_star_tt;
            let pref = r_diff_2 * w;
            let pref_t = r_diff_2_t * w;
            let pref_tt = r_diff_2_tt * w;
            let (e_raw, e_t, e_tt) = lj_boltzmann_factor(r, r_t, r_tt, temperature, eps_divk);
            let e = e_raw - 1.0;
            let r2 = r * r;
            let two_r_rt = 2.0 * r * r_t;
            let slot = N_QUAD + i;
            prepared.r_hs[slot] = r * inv_d;
            prepared.r_hs_t[slot] = r_t * inv_d - d_star_t * r * inv_d2;
            prepared.r_hs_tt[slot] =
                r_tt * inv_d - 2.0 * d_star_t * r_t * inv_d2 - d_star_tt * r * inv_d2
                    + 2.0 * d_star_t * d_star_t * r * inv_d3;
            prepared.weight[slot] = pref * e * r2;
            prepared.weight_t[slot] = pref_t * e * r2 + pref * e_t * r2 + pref * e * two_r_rt;
            prepared.weight_tt[slot] = pref_tt * e * r2
                + 2.0 * pref_t * e_t * r2
                + 2.0 * pref_t * e * two_r_rt
                + pref * e_tt * r2
                + 2.0 * pref * e_t * two_r_rt
                + 2.0 * pref * e * (r_t * r_t + r * r_tt);
        }

        Ok(prepared)
    }

    /// Evaluate full `(eta, T)` derivatives at packing fraction `eta`.
    pub(crate) fn evaluate_eta_t(&self, eta: f64) -> Result<CavityEtaDerivatives> {
        if !eta.is_finite() || eta < 0.0 {
            return Err(Error::InvalidInput(
                "nonnegative finite packing fraction required for cavity integral",
            ));
        }
        let params = if eta <= ETA_C {
            cavity_parameters(eta)
        } else {
            cavity_parameters(ETA_C)
        };
        let (a, b, c, a_e, b_e, c_e, a_ee, b_ee, c_ee) = params;
        let diff = eta - ETA_C;
        let extrapolated = eta > ETA_C;

        let mut i_val = 0.0;
        let mut i_e = 0.0;
        let mut i_ee = 0.0;
        let mut i_t = 0.0;
        let mut i_tt = 0.0;
        let mut i_et = 0.0;

        for idx in 0..(2 * N_QUAD) {
            let rs = self.r_hs[idx];
            let rs_t = self.r_hs_t[idx];
            let rs_tt = self.r_hs_tt[idx];
            let w = self.weight[idx];
            let w_t = self.weight_t[idx];
            let w_tt = self.weight_tt[idx];
            let rs2 = rs * rs;
            let rs3 = rs2 * rs;

            let y0 = (a + b * rs + c * rs3).exp();
            let slope = a_e + b_e * rs + c_e * rs3;
            let y_e0 = y0 * slope;
            let slope_ee = a_ee + b_ee * rs + c_ee * rs3;
            let y_ee0 = y0 * (slope * slope + slope_ee);
            let dr = b + 3.0 * c * rs2;
            let dre = b_e + 3.0 * c_e * rs2;
            let y_r0 = y0 * dr;
            let y_rr0 = y0 * (dr * dr + 6.0 * c * rs);
            let y_er0 = y0 * (slope * dr + dre);

            let (y, y_e, y_ee, y_r, y_rr, y_er) = if !extrapolated {
                (y0, y_e0, y_ee0, y_r0, y_rr0, y_er0)
            } else {
                let dree = b_ee + 3.0 * c_ee * rs2;
                let y_err0 =
                    y0 * (2.0 * dr * dre + slope * (6.0 * c * rs + dr * dr) + 6.0 * c_e * rs);
                let y_eer0 = y_ee0 * dr + y0 * (2.0 * slope * dre + dree);
                let y_eerr0 = y_eer0 * dr
                    + y0 * (6.0 * c * rs * slope * slope
                        + 2.0 * dr * slope * dre
                        + 6.0 * c * rs * slope_ee
                        + dr * dree
                        + 2.0 * dre * dre
                        + 12.0 * slope * c_e * rs
                        + 6.0 * c_ee * rs);
                let half_d2 = 0.5 * diff * diff;
                (
                    y0 + y_e0 * diff + y_ee0 * half_d2,
                    y_e0 + y_ee0 * diff,
                    y_ee0,
                    y_r0 + y_er0 * diff + y_eer0 * half_d2,
                    y_rr0 + y_err0 * diff + y_eerr0 * half_d2,
                    y_er0 + y_eer0 * diff,
                )
            };

            let y_t = y_r * rs_t;
            let y_tt = y_rr * rs_t * rs_t + y_r * rs_tt;
            let y_et = y_er * rs_t;

            i_val += w * y;
            i_e += w * y_e;
            i_ee += w * y_ee;
            i_t += w_t * y + w * y_t;
            i_et += w_t * y_e + w * y_et;
            i_tt += w_tt * y + 2.0 * w_t * y_t + w * y_tt;
        }

        let d3 = self.d_star * self.d_star * self.d_star;
        let inv_d3 = 1.0 / d3;
        let inv_d4 = inv_d3 / self.d_star;
        let inv_d5 = inv_d4 / self.d_star;

        let value = -12.0 * eta * i_val * inv_d3;
        let fn_e = -12.0 * (i_val + eta * i_e) * inv_d3;
        let fn_ee = -12.0 * (2.0 * i_e + eta * i_ee) * inv_d3;
        let fn_et = -12.0 * (i_t + eta * i_et) * inv_d3
            + 36.0 * (i_val + eta * i_e) * self.d_star_t * inv_d4;
        let fn_xt = -12.0 * eta * (i_t - 3.0 * i_val * self.d_star_t / self.d_star) * inv_d3;
        let fn_xtt = -12.0 * eta * i_tt * inv_d3
            + 72.0 * eta * i_t * self.d_star_t * inv_d4
            + 36.0 * eta * i_val * self.d_star_tt * inv_d4
            - 144.0 * eta * i_val * self.d_star_t * self.d_star_t * inv_d5;

        let out = CavityEtaDerivatives {
            value,
            eta: fn_e,
            eta_eta: fn_ee,
            temperature: fn_xt,
            temperature_second: fn_xtt,
            eta_temperature: fn_et,
        };
        if ![
            out.value,
            out.eta,
            out.eta_eta,
            out.temperature,
            out.temperature_second,
            out.eta_temperature,
        ]
        .iter()
        .all(|v| v.is_finite())
        {
            return Err(Error::InvalidInput("non-finite cavity integral evaluation"));
        }
        Ok(out)
    }

    /// Fast density-only `(value, d/deta, d^2/deta^2)` evaluation for pressure and
    /// volume root finding (skips all temperature derivatives).
    pub(crate) fn evaluate_eta_only(&self, eta: f64) -> Result<(f64, f64, f64)> {
        if !eta.is_finite() || eta < 0.0 {
            return Err(Error::InvalidInput(
                "nonnegative finite packing fraction required for cavity integral",
            ));
        }
        let (a, b, c, a_e, b_e, c_e, a_ee, b_ee, c_ee) = if eta <= ETA_C {
            cavity_parameters(eta)
        } else {
            cavity_parameters(ETA_C)
        };
        let diff = eta - ETA_C;
        let extrapolated = eta > ETA_C;
        let half_d2 = 0.5 * diff * diff;

        let mut i_val = 0.0;
        let mut i_e = 0.0;
        let mut i_ee = 0.0;

        for idx in 0..(2 * N_QUAD) {
            let rs = self.r_hs[idx];
            let rs3 = rs * rs * rs;
            let w = self.weight[idx];

            let y0 = (a + b * rs + c * rs3).exp();
            let slope = a_e + b_e * rs + c_e * rs3;
            let y_e0 = y0 * slope;
            let y_ee0 = y0 * (slope * slope + a_ee + b_ee * rs + c_ee * rs3);

            let (y, y_e, y_ee) = if !extrapolated {
                (y0, y_e0, y_ee0)
            } else {
                (
                    y0 + y_e0 * diff + y_ee0 * half_d2,
                    y_e0 + y_ee0 * diff,
                    y_ee0,
                )
            };

            i_val += w * y;
            i_e += w * y_e;
            i_ee += w * y_ee;
        }

        let inv_d3 = 1.0 / (self.d_star * self.d_star * self.d_star);
        let value = -12.0 * eta * i_val * inv_d3;
        let fn_e = -12.0 * (i_val + eta * i_e) * inv_d3;
        let fn_ee = -12.0 * (2.0 * i_e + eta * i_ee) * inv_d3;
        if !value.is_finite() || !fn_e.is_finite() || !fn_ee.is_finite() {
            return Err(Error::InvalidInput(
                "non-finite cavity integral density evaluation",
            ));
        }
        Ok((value, fn_e, fn_ee))
    }

    /// Evaluate the molar cavity Helmholtz energy and all 6 `(T, rho)` derivatives.
    pub(crate) fn evaluate_density(
        &self,
        eta_per_rho: f64,
        eta_per_rho_t: f64,
        eta_per_rho_tt: f64,
        rho: f64,
    ) -> Result<DensityHelmholtz> {
        if !rho.is_finite() || rho < 0.0 {
            return Err(Error::InvalidInput("nonnegative molar density required"));
        }
        let eta = eta_per_rho * rho;
        let eta_t = eta_per_rho_t * rho;
        let eta_tt = eta_per_rho_tt * rho;
        let cav = self.evaluate_eta_t(eta)?;
        Ok(DensityHelmholtz {
            value: cav.value,
            dt: cav.temperature + cav.eta * eta_t,
            drho: cav.eta * eta_per_rho,
            dtt: cav.temperature_second
                + 2.0 * cav.eta_temperature * eta_t
                + cav.eta_eta * eta_t * eta_t
                + cav.eta * eta_tt,
            dtrho: cav.eta_temperature * eta_per_rho
                + cav.eta * eta_per_rho_t
                + cav.eta_eta * eta_t * eta_per_rho,
            drhorho: cav.eta_eta * eta_per_rho * eta_per_rho,
        })
    }
}

/// The WCA hard-sphere diameter and its first and second temperature derivatives.
///
/// Ported from `src/hardsphere_wca.f90:calc_dhs_WCA`.
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

/// Pure-fluid molar hard-sphere packing fraction factor `eta / rho = (pi/6) N_A d^3`
/// and its first and second temperature derivatives.
pub(crate) fn pure_eta_per_rho(
    diameter: f64,
    diameter_t: f64,
    diameter_tt: f64,
) -> (f64, f64, f64) {
    let prefactor = PI / 6.0 * crate::N_A;
    let d2 = diameter * diameter;
    let eta1 = prefactor * d2 * diameter;
    let eta1_t = prefactor * 3.0 * d2 * diameter_t;
    let eta1_tt = prefactor * (6.0 * diameter * diameter_t * diameter_t + 3.0 * d2 * diameter_tt);
    (eta1, eta1_t, eta1_tt)
}

/// Pure-fluid hard-sphere reduced molar Helmholtz energy `a_hs(eta)` and its
/// first and second packing-fraction derivatives (`PURE_HS_CS` or `PURE_HS_CSK`).
///
/// Ported from `src/saftvrmie_hardsphere.f90:calc_pure_ahs_div_nRT_eta`.
pub fn pure_hard_sphere_eta(eta: f64, kolafa: bool) -> Result<(f64, f64, f64)> {
    if !eta.is_finite() || !(0.0..1.0).contains(&eta) {
        return Err(Error::InvalidInput(
            "hard-sphere packing fraction must be in [0, 1)",
        ));
    }
    let om = 1.0 - eta;
    let om2 = om * om;
    let om3 = om2 * om;
    let om4 = om2 * om2;
    if kolafa {
        let a = (2.5 / om2 + 10.0 / om + 2.0 * eta + 5.0 * om.ln()) / 3.0;
        let a_e = (12.0 - 6.0 * eta + eta * eta - 2.0 * eta * eta * eta) / (3.0 * om3);
        let a_ee = -(5.0 / 3.0) * (-6.0 + 2.0 * eta + eta * eta) / om4;
        Ok((a, a_e, a_ee))
    } else {
        let a = (4.0 * eta - 3.0 * eta * eta) / om2;
        let a_e = 2.0 * (2.0 - eta) / om3;
        let a_ee = (10.0 - 4.0 * eta) / om4;
        Ok((a, a_e, a_ee))
    }
}

/// Pure-fluid molar hard-sphere reduced Helmholtz energy and all 6 `(T, rho)`
/// derivatives.
pub fn pure_hard_sphere_density(
    eta_per_rho: f64,
    eta_per_rho_t: f64,
    eta_per_rho_tt: f64,
    rho: f64,
    kolafa: bool,
) -> Result<DensityHelmholtz> {
    if !rho.is_finite() || rho < 0.0 {
        return Err(Error::InvalidInput("nonnegative molar density required"));
    }
    let eta = eta_per_rho * rho;
    let eta_t = eta_per_rho_t * rho;
    let eta_tt = eta_per_rho_tt * rho;
    let (a, a_e, a_ee) = pure_hard_sphere_eta(eta, kolafa)?;
    Ok(DensityHelmholtz {
        value: a,
        dt: a_e * eta_t,
        drho: a_e * eta_per_rho,
        dtt: a_e * eta_tt + a_ee * eta_t * eta_t,
        dtrho: a_e * eta_per_rho_t + a_ee * eta_t * eta_per_rho,
        drhorho: a_ee * eta_per_rho * eta_per_rho,
    })
}

/// Hypothetical pure-fluid packing fraction `zeta_x` without segment numbers
/// and its full first, second, and volume third derivatives.
///
/// Ported from `src/hardsphere_wca.f90:calcZetaX_vdW_no_segments`.
/// `d`, `d_t`, and `d_tt` are symmetric `nc x nc` row-major diameter matrices.
pub fn packing_fraction_vdw_no_segments(
    volume: f64,
    moles: &[f64],
    d: &[f64],
    d_t: &[f64],
    d_tt: &[f64],
) -> Result<PackingFractionDerivatives> {
    positive(volume, "total volume")?;
    let nc = moles.len();
    if nc == 0 {
        return Err(Error::InvalidInput("at least one component required"));
    }
    if d.len() != nc * nc || d_t.len() != nc * nc || d_tt.len() != nc * nc {
        return Err(Error::InvalidInput(
            "diameter matrices must have length nc * nc",
        ));
    }
    let ns: f64 = moles.iter().sum();
    positive(ns, "total moles")?;

    let prefactor = PI * crate::N_A / (6.0 * volume * ns);
    let mut zx = 0.0;
    let mut zx_t = 0.0;
    let mut zx_tt = 0.0;
    let mut sum_d3 = vec![0.0; nc];
    let mut sum_d2_dt = vec![0.0; nc];

    for i in 0..nc {
        let mut s_d3 = 0.0;
        let mut s_d2_dt = 0.0;
        let mut s_tt = 0.0;
        for (j, &nj) in moles.iter().enumerate() {
            let idx = j * nc + i;
            let dij = d[idx];
            let dtij = d_t[idx];
            let dttij = d_tt[idx];
            let d2 = dij * dij;
            s_d3 += nj * d2 * dij;
            s_d2_dt += nj * d2 * dtij;
            s_tt += nj * (6.0 * dij * dtij * dtij + 3.0 * d2 * dttij);
        }
        sum_d3[i] = s_d3;
        sum_d2_dt[i] = s_d2_dt;
        zx += moles[i] * s_d3;
        zx_t += moles[i] * s_d2_dt;
        zx_tt += prefactor * moles[i] * s_tt;
    }
    zx *= prefactor;
    zx_t *= 3.0 * prefactor;

    let inv_v = 1.0 / volume;
    let inv_v2 = inv_v * inv_v;
    let inv_v3 = inv_v2 * inv_v;
    let zx_v = -zx * inv_v;
    let zx_vv = 2.0 * zx * inv_v2;
    let zx_tv = -zx_t * inv_v;
    let zx_vvv = -6.0 * zx * inv_v3;
    let zx_vvt = 2.0 * zx_t * inv_v2;
    let zx_vtt = -zx_tt * inv_v;

    let mut zx_n = vec![0.0; nc];
    let mut zx_tn = vec![0.0; nc];
    let mut zx_vn = vec![0.0; nc];
    let mut zx_vtn = vec![0.0; nc];
    for k in 0..nc {
        zx_n[k] = 2.0 * prefactor * sum_d3[k] - zx / ns;
        zx_tn[k] = 6.0 * prefactor * sum_d2_dt[k] - zx_t / ns;
        zx_vn[k] = -zx_n[k] * inv_v;
        zx_vtn[k] = -zx_tn[k] * inv_v;
    }

    let mut zx_nn = vec![0.0; nc * nc];
    let mut zx_vnn = vec![0.0; nc * nc];
    for k in 0..nc {
        for l in 0..nc {
            let dkl = d[k * nc + l];
            let val = 2.0 * prefactor * dkl * dkl * dkl - (zx_n[l] + zx_n[k]) / ns;
            zx_nn[k * nc + l] = val;
            zx_vnn[k * nc + l] = -val * inv_v;
        }
    }

    Ok(PackingFractionDerivatives {
        value: zx,
        temperature: zx_t,
        temperature_second: zx_tt,
        volume: zx_v,
        volume_second: zx_vv,
        volume_third: zx_vvv,
        temperature_volume: zx_tv,
        volume_temperature_second: zx_vvt,
        volume_temperature_third: zx_vtt,
        moles: zx_n,
        moles_second: zx_nn,
        temperature_moles: zx_tn,
        volume_moles: zx_vn,
        volume_moles_second: zx_vnn,
        volume_temperature_moles: zx_vtn,
    })
}

/// Evaluate the molar reduced cavity Helmholtz energy `F_n(eta, T)` and its
/// partial derivatives with respect to packing fraction `eta` and temperature `T`.
///
/// Uses the 61-point Gauss-Kronrod quadrature of `src/hardsphere_wca.f90:calc_cavity_integral_LJ_Fres`.
pub fn cavity_integral_eta_t(
    sigma: f64,
    eps_divk: f64,
    diameter: f64,
    diameter_t: f64,
    diameter_tt: f64,
    eta: f64,
    temperature: f64,
) -> Result<CavityEtaDerivatives> {
    let prepared = CavityPreparedTemp::new(
        sigma,
        eps_divk,
        diameter,
        diameter_t,
        diameter_tt,
        temperature,
    )?;
    prepared.evaluate_eta_t(eta)
}

/// Evaluate the molar reduced cavity Helmholtz energy and all 6 `(T, rho)`
/// derivatives for a pure fluid.
pub fn cavity_integral_density(
    sigma: f64,
    eps_divk: f64,
    diameter: f64,
    diameter_t: f64,
    diameter_tt: f64,
    temperature: f64,
    rho: f64,
) -> Result<DensityHelmholtz> {
    let prepared = CavityPreparedTemp::new(
        sigma,
        eps_divk,
        diameter,
        diameter_t,
        diameter_tt,
        temperature,
    )?;
    let (eta1, eta1_t, eta1_tt) = pure_eta_per_rho(diameter, diameter_t, diameter_tt);
    prepared.evaluate_density(eta1, eta1_t, eta1_tt, rho)
}

/// Evaluate the WCA cavity soft-repulsion Helmholtz contribution and all 10
/// `(T, V, n)` derivatives.
///
/// Ported from `src/hardsphere_wca.f90:calc_cavity_integral_LJ_Fres`.
/// When `disable_n_multiplication` is `false`, returns the extensive reduced
/// Helmholtz contribution `F = n * a_soft` (used by `calcFresLJs_WCA`,
/// `calcFresLJs_uv_theory`, and `calcFresLJ_uf_theory`). When `true`, returns
/// the molar reduced contribution `a_soft` and its `(T, V, n)` derivatives
/// (used by `calc_wca_soft_repulsion`).
#[allow(clippy::too_many_arguments)]
pub fn cavity_integral_lj_fres(
    sigma: f64,
    eps_divk: f64,
    diameter: f64,
    diameter_t: f64,
    diameter_tt: f64,
    temperature: f64,
    volume: f64,
    moles: f64,
    disable_n_multiplication: bool,
) -> Result<PureHelmholtz> {
    positive(volume, "total volume")?;
    positive(moles, "mole number")?;
    let dh = cavity_integral_density(
        sigma,
        eps_divk,
        diameter,
        diameter_t,
        diameter_tt,
        temperature,
        moles / volume,
    )?;
    if disable_n_multiplication {
        dh.molar_tvn(volume, moles)
    } else {
        dh.extensive(volume, moles)
    }
}

/// Evaluate the extensive WCA cavity Helmholtz contribution and all 10 `(T, V, n)`
/// derivatives using the WCA hard-sphere diameter.
///
/// Ported from `src/hardsphere_wca.f90:calc_cavity_integral_LJ_Fres_test_TVN`.
pub fn cavity_integral_lj_fres_tvn(
    sigma: f64,
    eps_divk: f64,
    temperature: f64,
    volume: f64,
    moles: f64,
) -> Result<PureHelmholtz> {
    positive(sigma, "segment diameter")?;
    positive(eps_divk, "well depth")?;
    positive(temperature, "temperature")?;
    let (d, d_t, d_tt) = hard_sphere_diameter(sigma, eps_divk, temperature);
    cavity_integral_lj_fres(
        sigma,
        eps_divk,
        d,
        d_t,
        d_tt,
        temperature,
        volume,
        moles,
        false,
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
        let dr = b + 3.0 * c * r_star * r_star;
        let r_first = value * dr;
        let r_second = value * (dr * dr + 6.0 * c * r_star);
        let r_mixed = value * (slope * dr + b_e + 3.0 * c_e * r_star * r_star);
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
        let dr = b + 3.0 * c * r_star * r_star;
        let dre = b_e + 3.0 * c_e * r_star * r_star;
        let r_first = value * dr;
        let r_second = value * (dr * dr + 6.0 * c * r_star);
        let r_mixed = value * (slope * dr + dre);
        let r_mixed_first =
            value * (2.0 * dr * dre + slope * (6.0 * c * r_star + dr * dr) + 6.0 * c_e * r_star);
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
                    + 2.0 * (dre * dre)
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
        let inverse6 = crate::numerics::ipow(inverse, 6);
        let inverse12 = inverse6 * inverse6;
        (
            4.0 * (inverse12 - inverse6) + 1.0,
            -4.0 * (12.0 * inverse12 - 6.0 * inverse6) * inverse,
            4.0 * (12.0 * 13.0 * inverse12 - 6.0 * 7.0 * inverse6) * inverse * inverse,
        )
    };
    let value = (-u * eps_divk / temperature).exp();
    let du_dt = u_r * r_t * eps_divk / temperature - u * eps_divk / (temperature * temperature);
    let value_t = -value * du_dt;
    let value_tt = value
        * (du_dt * du_dt
            - (u_rr * r_t * r_t * eps_divk / temperature + u_r * r_tt * eps_divk / temperature
                - 2.0 * u_r * r_t * eps_divk / (temperature * temperature)
                + 2.0 * u * eps_divk / (temperature * temperature * temperature)));
    (value, value_t, value_tt)
}

/// The parameters of the cavity correlation function of the source.
fn cavity_parameters(eta: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    let om = 1.0 - eta;
    let denominator = om * om * om;
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
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta2 * eta2;
    let em2 = eta - 2.0;
    let c_ee = -(2.0 * crate::numerics::ipow(eta, 5) - 10.0 * eta4 + 5.0 * eta3 + 31.0 * eta2
        - 29.0 * eta
        - 11.0)
        / ((em2 * em2) * denominator);
    (a, b, c, a_e, b_e, c_e, a_ee, b_ee, c_ee)
}
