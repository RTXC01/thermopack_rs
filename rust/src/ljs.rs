//! Lennard-Jones splined models: Barker-Henderson (`LJS-BH`), UV theory
//! (`LJS-UV`), and UF theory (`LJS-UF`, `LJ-UF`).
//! Ported from `src/lj_splined.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.

#![allow(clippy::excessive_precision)]

use crate::Phase;
use crate::chain_rule::PackingFractionDerivatives;
use crate::error::{Error, Result, positive};
use crate::hardsphere_wca::{
    CavityPreparedTemp, hard_sphere_diameter as wca_hard_sphere_diameter,
    packing_fraction_vdw_no_segments, pure_eta_per_rho, pure_hard_sphere_density,
    pure_hard_sphere_eta,
};
use crate::helmholtz::{DensityHelmholtz, PureHelmholtz, PurePressure};
use std::f64::consts::PI;

/// The polynomial coefficients of the source, in the order of the source
/// array. The source stores them in a four by fifteen array, and the first
/// polynomial of a row holds its first five coefficients, the second
/// polynomial the next five and the third the last five.
const POLYNOMIALS: [f64; 60] = [
    0.039813685,
    -0.106920241,
    1.291256810,
    0.000000000,
    0.469671161,
    -0.295578203,
    -2.763380451,
    21.299989640,
    -0.346120315,
    0.343721691,
    1.416116081,
    -11.284586250,
    -0.344778684,
    0.248954221,
    0.331595545,
    -0.258500011,
    -0.535278660,
    0.358663104,
    -0.248745451,
    0.996821621,
    0.000000000,
    0.000000000,
    0.007804235,
    0.000000000,
    -7.861555167,
    7.084518499,
    2.799119819,
    -540.982471100,
    10.084157010,
    -9.814532406,
    -4.455465501,
    599.182974200,
    0.518650855,
    0.172122707,
    1.855450051,
    -116.677210100,
    0.000000000,
    0.000000000,
    -0.350181778,
    3.804716369,
    0.000000000,
    0.000000000,
    0.000000000,
    0.000000000,
    32.718743790,
    -28.998566930,
    0.000000000,
    3131.407636000,
    -51.573804390,
    47.108198520,
    0.000000000,
    -4126.665342000,
    6.550964581,
    -7.263318791,
    0.000000000,
    1013.877746000,
    0.000000000,
    0.000000000,
    0.000000000,
    -70.079287230,
];

/// The segment diameter of the argon record of the source database in metre.
pub const ARGON_SIGMA: f64 = 3.42e-10;

/// The well depth of the argon record of the source database in kelvin.
pub const ARGON_EPS: f64 = 124.0;

/// The attractive exponent of the model.
const LAMBDA_A: f64 = 6.0;

/// The repulsive exponent of the model.
const LAMBDA_R: f64 = 12.0;

/// The prefactor of the potential of the model.
const C_FAC: f64 = 4.0;

/// The van der Waals attractive constant of the splined Lennard-Jones potential
/// (`src/lj_splined.f90:alpha_ljs`).
pub const ALPHA_LJS: f64 = 0.5357287350120993;

/// Overflow-safe `(tanh(x), sech^2(x), d^2(tanh)/dx^2)` preserving exact
/// `(e^(2x)-1)/(e^(2x)+1)` arithmetic on finite domains.
pub(crate) fn safe_tanh_derivs(x: f64) -> (f64, f64, f64) {
    let two_x = 2.0 * x;
    if two_x >= 700.0 {
        (1.0, 0.0, 0.0)
    } else if two_x <= -700.0 {
        (-1.0, 0.0, 0.0)
    } else {
        let e2x = two_x.exp();
        let t0 = (e2x - 1.0) / (e2x + 1.0);
        let t1 = 1.0 - t0 * t0;
        let t2 = -2.0 * t0 * t1;
        (t0, t1, t2)
    }
}

/// Overflow-safe `tanh(x(r))` and its first and second derivatives with respect
/// to `r`, matching Fortran's `(e2x - 1)/(e2x + 1)` expression on finite domains.
fn safe_tanh_chain(x: f64, x_r: f64, x_rr: f64) -> (f64, f64, f64) {
    let two_x = 2.0 * x;
    if two_x >= 700.0 {
        (1.0, 0.0, 0.0)
    } else if two_x <= -700.0 {
        (-1.0, 0.0, 0.0)
    } else {
        let e2x = two_x.exp();
        let e2x_p1 = e2x + 1.0;
        let tanh = (e2x - 1.0) / e2x_p1;
        let denom = e2x_p1 * e2x_p1;
        let tanh_r = 4.0 * e2x * x_r / denom;
        let tanh_rr = 2.0 * tanh_r * x_r + 4.0 * e2x * x_rr / denom
            - 16.0 * e2x * e2x * x_r * x_r / (denom * e2x_p1);
        (tanh, tanh_r, tanh_rr)
    }
}

/// Legacy summary state of the model at a temperature and a volume.
#[derive(Debug, Clone, Copy)]
pub struct LjsState {
    /// Pressure of the state in pascal.
    pub pressure: f64,
    /// Compressibility factor of the state.
    pub compressibility: f64,
    /// Molar reduced residual Helmholtz energy of the state.
    pub residual_helmholtz: f64,
    /// Molar residual internal energy in joule per mole.
    pub u_residual: f64,
    /// Molar residual enthalpy in joule per mole.
    pub h_residual: f64,
    /// Molar residual entropy in joule per mole and kelvin.
    pub s_residual: f64,
    /// Molar enthalpy of the state in joule per mole.
    pub h_total: f64,
    /// Molar entropy of the state in joule per mole and kelvin.
    pub s_total: f64,
}

/// Scalar property with its `(T, P)` partial derivatives.
#[derive(Debug, Clone, Copy)]
pub struct LjsPropertyDerivatives {
    /// Property value at `(T, P)`.
    pub value: f64,
    /// Partial derivative with respect to temperature at constant pressure.
    pub dt: f64,
    /// Partial derivative with respect to pressure at constant temperature.
    pub dp: f64,
}

/// Full `(T, V, n)` thermodynamic state for any LJS/LJ model.
#[derive(Debug, Clone, Copy)]
pub struct LjsTvState {
    /// Temperature in kelvin.
    pub temperature: f64,
    /// Total volume in cubic metre.
    pub volume: f64,
    /// Total mole number in mole.
    pub moles: f64,
    /// Pressure and `(T, V, n)` derivatives.
    pub pressure: PurePressure,
    /// Compressibility factor `Z = P V / (n R T)`.
    pub compressibility: f64,
    /// Extensive reduced residual Helmholtz energy and all 10 `(T, V, n)` derivatives.
    pub residual_helmholtz: PureHelmholtz,
    /// Molar reduced residual Helmholtz energy `a_res = F / n`.
    pub a_residual: f64,
    /// Molar residual internal energy in J/mol.
    pub u_residual: f64,
    /// Molar residual enthalpy in J/mol.
    pub h_residual: f64,
    /// Molar residual entropy `s(T, v) - s_id(T, P)` in J/(mol K).
    pub s_residual: f64,
    /// Molar residual isochoric heat capacity in J/(mol K).
    pub cv_residual: f64,
    /// Molar residual isobaric heat capacity in J/(mol K).
    pub cp_residual: f64,
    /// Total molar enthalpy in J/mol.
    pub h_total: f64,
    /// Total molar entropy in J/(mol K).
    pub s_total: f64,
    /// Total molar isochoric heat capacity in J/(mol K).
    pub cv_total: f64,
    /// Total molar isobaric heat capacity in J/(mol K).
    pub cp_total: f64,
    /// Logarithm of fugacity coefficient `ln(phi) = F_n - ln(Z)`.
    pub ln_fugacity_coefficient: f64,
    /// Speed of sound in m/s (`NaN` if mechanically unstable or `cp/cv <= 0`).
    pub speed_of_sound: f64,
}

impl LjsTvState {
    /// Convert to the legacy `LjsState` summary.
    pub fn to_legacy_state(&self) -> LjsState {
        LjsState {
            pressure: self.pressure.value,
            compressibility: self.compressibility,
            residual_helmholtz: self.a_residual,
            u_residual: self.u_residual,
            h_residual: self.h_residual,
            s_residual: self.s_residual,
            h_total: self.h_total,
            s_total: self.s_total,
        }
    }
}

/// Full `(T, P)` thermodynamic state and `(T, P)` derivatives for any LJS/LJ model.
#[derive(Debug, Clone, Copy)]
pub struct LjsTpState {
    /// Temperature in kelvin.
    pub temperature: f64,
    /// Pressure in pascal.
    pub pressure: f64,
    /// Molar volume in m^3/mol.
    pub molar_volume: f64,
    /// `(dv/dT)_P` in m^3/(mol K).
    pub dv_dt: f64,
    /// `(dv/dP)_T` in m^3/(mol Pa).
    pub dv_dp: f64,
    /// Compressibility factor `Z` and its `(T, P)` derivatives.
    pub compressibility: LjsPropertyDerivatives,
    /// Total molar enthalpy (J/mol) and its `(T, P)` derivatives.
    pub enthalpy: LjsPropertyDerivatives,
    /// Total molar entropy (J/(mol K)) and its `(T, P)` derivatives.
    pub entropy: LjsPropertyDerivatives,
    /// `ln(phi)` and its `(T, P)` derivatives.
    pub ln_fugacity_coefficient: LjsPropertyDerivatives,
    /// Total molar isochoric heat capacity in J/(mol K).
    pub cv_total: f64,
    /// Total molar isobaric heat capacity in J/(mol K).
    pub cp_total: f64,
    /// Speed of sound in m/s.
    pub speed_of_sound: f64,
    /// Underlying `(T, V, n = 1)` state.
    pub tv: LjsTvState,
}

/// Result of ThermoPack's exact source SAFT non-association volume iteration
/// (`src/saft_interface.f90:pc_saft_nonassoc_volume_solver`).
#[derive(Debug, Clone, Copy)]
pub struct LjsVolumeIteration {
    /// Final molar volume iterate in m^3/mol.
    pub volume: f64,
    /// Iteration count executed.
    pub iterations: usize,
    /// Whether `|V - Vold| < 1e5 * eps * Vold` was satisfied.
    pub volume_converged: bool,
    /// Whether `|P - P_spec| < 1e8 * eps * P_spec` was satisfied.
    pub pressure_converged: bool,
    /// Relative pressure residual `|P - P_spec| / P_spec`.
    pub relative_pressure_error: f64,
}

/// Argon ideal-gas model configured with constant `Cp_id = 2.5 * R`, matching
/// `init_LJs_bh` / `init_LJs_WCA`, together with molar mass in kg/mol.
pub(crate) fn ljs_ideal_gas() -> Result<(crate::ideal::IdealGas, f64)> {
    let component = crate::component::component("AR")?;
    let molar_mass_kg = component.molar_mass * 1.0e-3;
    let mut gas = crate::ideal::IdealGas::new(component)?;
    let mut coefficients = [0.0; 21];
    coefficients[0] = 2.5;
    gas.set_heat_capacity(8, &coefficients)?;
    Ok((gas, molar_mass_kg))
}

/// Compute the full `(T, V, n)` thermodynamic state from extensive `F(T, V, n)`.
pub(crate) fn evaluate_tv_state(
    temperature: f64,
    volume: f64,
    moles: f64,
    f: PureHelmholtz,
) -> Result<LjsTvState> {
    let pressure = f.pressure(temperature, volume, moles)?;
    let z = pressure.value * volume / (moles * crate::R * temperature);
    let a_res = f.value / moles;
    let f_t_molar = f.dt / moles;
    let f_tt_molar = f.dtt / moles;
    let u_res = -crate::R * temperature * temperature * f_t_molar;
    let h_res = u_res + crate::R * temperature * (z - 1.0);
    let s_res = if z > 0.0 {
        -crate::R * (temperature * f_t_molar + a_res) + crate::R * z.ln()
    } else {
        f64::NAN
    };
    let cv_res = -crate::R * temperature * (2.0 * f_t_molar + temperature * f_tt_molar);
    let cp_res =
        cv_res - temperature * pressure.dt * pressure.dt / (moles * pressure.dv) - crate::R;

    let (ideal, molar_mass_kg) = ljs_ideal_gas()?;
    let (h_id, s_id, cp_id) = if pressure.value > 0.0 {
        let props = ideal.properties(temperature, pressure.value)?;
        (props.enthalpy, props.entropy, props.cp)
    } else {
        let props = ideal.properties(temperature, 1.0e5)?;
        (props.enthalpy, f64::NAN, props.cp)
    };
    let cv_id = cp_id - crate::R;
    let h_total = h_id + h_res;
    let s_total = s_id + s_res;
    let cv_total = cv_id + cv_res;
    let cp_total = cp_id + cp_res;
    let ln_phi = if z > 0.0 { f.dn - z.ln() } else { f64::NAN };
    let c2 = -(volume * volume / (moles * molar_mass_kg)) * (cp_total / cv_total) * pressure.dv;
    let speed_of_sound = if c2 > 0.0 && c2.is_finite() {
        c2.sqrt()
    } else {
        f64::NAN
    };

    Ok(LjsTvState {
        temperature,
        volume,
        moles,
        pressure,
        compressibility: z,
        residual_helmholtz: f,
        a_residual: a_res,
        u_residual: u_res,
        h_residual: h_res,
        s_residual: s_res,
        cv_residual: cv_res,
        cp_residual: cp_res,
        h_total,
        s_total,
        cv_total,
        cp_total,
        ln_fugacity_coefficient: ln_phi,
        speed_of_sound,
    })
}

/// Compute the full `(T, P)` thermodynamic state and `(T, P)` derivatives at `n = 1 mol`.
pub(crate) fn evaluate_tp_state(
    temperature: f64,
    pressure: f64,
    molar_volume: f64,
    f: PureHelmholtz,
) -> Result<LjsTpState> {
    let tv = evaluate_tv_state(temperature, molar_volume, 1.0, f)?;
    let p_t = tv.pressure.dt;
    let p_v = tv.pressure.dv;
    let p_n = tv.pressure.dn;
    let dv_dt = -p_t / p_v;
    let dv_dp = 1.0 / p_v;
    let dv_dn = -p_n / p_v;
    let z = pressure * molar_volume / (crate::R * temperature);
    let dz_dt = z * (dv_dt / molar_volume - 1.0 / temperature);
    let dz_dp = z * (1.0 / pressure + dv_dp / molar_volume);
    let dh_dt = tv.cp_total;
    let dh_dp = molar_volume - temperature * dv_dt;
    let ds_dt = tv.cp_total / temperature;
    let ds_dp = -dv_dt;
    let dlnphi_dt = f.dtn + 1.0 / temperature - dv_dn * p_t / (crate::R * temperature);
    let dlnphi_dp = dv_dn / (crate::R * temperature) - 1.0 / pressure;

    Ok(LjsTpState {
        temperature,
        pressure,
        molar_volume,
        dv_dt,
        dv_dp,
        compressibility: LjsPropertyDerivatives {
            value: z,
            dt: dz_dt,
            dp: dz_dp,
        },
        enthalpy: LjsPropertyDerivatives {
            value: tv.h_total,
            dt: dh_dt,
            dp: dh_dp,
        },
        entropy: LjsPropertyDerivatives {
            value: tv.s_total,
            dt: ds_dt,
            dp: ds_dp,
        },
        ln_fugacity_coefficient: LjsPropertyDerivatives {
            value: tv.ln_fugacity_coefficient,
            dt: dlnphi_dt,
            dp: dlnphi_dp,
        },
        cv_total: tv.cv_total,
        cp_total: tv.cp_total,
        speed_of_sound: tv.speed_of_sound,
        tv,
    })
}

/// Refine a single bracket `[lo, hi]` (in packing fraction `eta`) where
/// `P(lo) <= p_spec <= P(hi)` using safeguarded Newton-Raphson.
fn refine_eta_root<F>(
    mut lo: f64,
    mut hi: f64,
    p_spec: f64,
    eta_per_rho: f64,
    p_and_dp_drho: &F,
) -> Result<Option<f64>>
where
    F: Fn(f64) -> Result<(f64, f64)>,
{
    let mut eta = 0.5 * (lo + hi);
    for _ in 0..80 {
        let rho = eta / eta_per_rho;
        let (p, dp_drho) = p_and_dp_drho(rho)?;
        let res = p - p_spec;
        if res.abs() <= 1.0e-13 * p_spec.max(1.0) || (hi - lo) <= 1.0e-14 * eta.max(1.0e-16) {
            return Ok(Some(eta));
        }
        if res < 0.0 {
            lo = eta;
        } else {
            hi = eta;
        }
        let dp_deta = dp_drho / eta_per_rho;
        let next = if dp_deta.is_finite() && dp_deta.abs() > 0.0 {
            let cand = eta - res / dp_deta;
            if cand > lo && cand < hi {
                cand
            } else {
                0.5 * (lo + hi)
            }
        } else {
            0.5 * (lo + hi)
        };
        eta = next;
    }
    let rho = eta / eta_per_rho;
    let (p, _) = p_and_dp_drho(rho)?;
    if (p - p_spec).abs() <= 1.0e-7 * p_spec.max(1.0) {
        Ok(Some(eta))
    } else {
        Ok(None)
    }
}

/// Locate a spinodal `dP/drho = 0` inside `[lo, hi]` where `dP/drho` changes sign.
fn bisect_spinodal<F>(
    mut lo: f64,
    mut hi: f64,
    lo_dp: f64,
    eta_per_rho: f64,
    p_and_dp_drho: &F,
) -> Result<(f64, f64)>
where
    F: Fn(f64) -> Result<(f64, f64)>,
{
    for _ in 0..48 {
        let mid = 0.5 * (lo + hi);
        let (_, dp) = p_and_dp_drho(mid / eta_per_rho)?;
        if dp == 0.0 || (hi - lo) <= 1.0e-13 * mid.max(1.0e-15) {
            lo = mid;
            hi = mid;
            break;
        }
        if (dp > 0.0) == (lo_dp > 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let eta = 0.5 * (lo + hi);
    let (p, _) = p_and_dp_drho(eta / eta_per_rho)?;
    Ok((eta, p))
}

/// Find all mechanically stable molar volumes (m^3/mol) at `(T, P)` on the
/// physical packing-fraction domain `eta in (0, 0.75]` (extended up to `0.90`
/// only when `P(0.75) < P_spec`). Returns molar volumes sorted ascending
/// (liquid first, vapor last).
pub(crate) fn solve_physical_volumes<F>(
    temperature: f64,
    pressure: f64,
    eta_per_rho: f64,
    p_and_dp_drho: F,
) -> Result<Vec<f64>>
where
    F: Fn(f64) -> Result<(f64, f64)>,
{
    positive(temperature, "temperature")?;
    positive(pressure, "pressure")?;
    positive(eta_per_rho, "packing fraction factor")?;

    let eta_ig = eta_per_rho * pressure / (crate::R * temperature);
    let eta_min = (0.01 * eta_ig).clamp(1.0e-16, 1.0e-10);
    let mut grid = Vec::with_capacity(110);
    // Log-spaced vapor grid from eta_min up to 1e-2
    let log_lo = eta_min.ln();
    let log_hi = (1.0e-2_f64).ln();
    let n_log = 20usize;
    for i in 0..n_log {
        let frac = i as f64 / n_log as f64;
        grid.push((log_lo + frac * (log_hi - log_lo)).exp());
    }
    // Include ideal-gas estimate directly when inside the physical range
    if eta_ig > eta_min && eta_ig < 0.75 {
        grid.push(eta_ig);
    }
    // Uniform dense-fluid grid from 1e-2 to 0.75
    let n_lin = 74usize;
    for i in 0..=n_lin {
        grid.push(0.01 + 0.74 * (i as f64 / n_lin as f64));
    }
    grid.sort_by(f64::total_cmp);
    grid.dedup_by(|a, b| (*a - *b).abs() <= 1.0e-15 * a.abs().max(b.abs()));

    let mut samples: Vec<(f64, f64, f64)> = Vec::with_capacity(grid.len() + 16);
    for &eta in &grid {
        let (p, dp) = p_and_dp_drho(eta / eta_per_rho)?;
        samples.push((eta, p, dp));
    }

    // If ultra-high pressure exceeds P(0.75), extend up to 0.90 (still below eta = 1 pole)
    if let Some(&(_, p_top, _)) = samples.last()
        && p_top < pressure
    {
        for step in 1..=15 {
            let eta = 0.75 + 0.01 * step as f64;
            let (p, dp) = p_and_dp_drho(eta / eta_per_rho)?;
            samples.push((eta, p, dp));
            if p >= pressure {
                break;
            }
        }
    }

    let mut root_etas = Vec::new();
    for pair in samples.windows(2) {
        let (e0, p0, dp0) = pair[0];
        let (e1, p1, dp1) = pair[1];
        let r0 = p0 - pressure;
        let r1 = p1 - pressure;

        if r0 == 0.0 && dp0 > 0.0 {
            root_etas.push(e0);
        }
        // Check for a spinodal inside (e0, e1) so narrow brackets near extrema are not missed.
        if (dp0 > 0.0 && dp1 < 0.0) || (dp0 < 0.0 && dp1 > 0.0) {
            let (es, ps) = bisect_spinodal(e0, e1, dp0, eta_per_rho, &p_and_dp_drho)?;
            let rs = ps - pressure;
            if dp0 > 0.0 && r0 <= 0.0 && rs >= 0.0 {
                if let Some(root) = refine_eta_root(e0, es, pressure, eta_per_rho, &p_and_dp_drho)?
                {
                    root_etas.push(root);
                }
            } else if dp1 > 0.0
                && rs <= 0.0
                && r1 >= 0.0
                && let Some(root) = refine_eta_root(es, e1, pressure, eta_per_rho, &p_and_dp_drho)?
            {
                root_etas.push(root);
            }
        } else if r0 < 0.0
            && r1 >= 0.0
            && (dp0 > 0.0 || dp1 > 0.0)
            && let Some(root) = refine_eta_root(e0, e1, pressure, eta_per_rho, &p_and_dp_drho)?
        {
            root_etas.push(root);
        }
    }

    // Deduplicate and verify mechanical stability (dP/drho > 0).
    root_etas.sort_by(f64::total_cmp);
    let mut volumes = Vec::new();
    let mut prev_eta: Option<f64> = None;
    for eta in root_etas {
        if let Some(prev) = prev_eta
            && (eta - prev).abs() <= 1.0e-9 * eta.max(1.0e-12)
        {
            continue;
        }
        let rho = eta / eta_per_rho;
        let (p, dp) = p_and_dp_drho(rho)?;
        if dp > 0.0 && (p - pressure).abs() <= 1.0e-6 * pressure.max(1.0) {
            volumes.push(eta_per_rho / eta);
            prev_eta = Some(eta);
        }
    }

    // Sort ascending molar volume (liquid = smallest v first, vapor = largest v last).
    volumes.sort_by(f64::total_cmp);
    if volumes.is_empty() {
        return Err(Error::NotConverged(
            "no mechanically stable physical volume root found",
        ));
    }
    Ok(volumes)
}

/// Select a physically stable molar volume (m^3/mol) according to `phase`.
pub(crate) fn select_physical_volume<F, G>(
    temperature: f64,
    pressure: f64,
    phase: Phase,
    eta_per_rho: f64,
    p_and_dp_drho: F,
    a_res_at_rho: G,
) -> Result<f64>
where
    F: Fn(f64) -> Result<(f64, f64)>,
    G: Fn(f64) -> Result<f64>,
{
    let roots = solve_physical_volumes(temperature, pressure, eta_per_rho, p_and_dp_drho)?;
    match phase {
        Phase::Liquid => Ok(roots[0]),
        Phase::Vapor => Ok(*roots.last().expect("non-empty roots")),
        Phase::MinimumGibbs => {
            if roots.len() == 1 {
                return Ok(roots[0]);
            }
            let mut best_v = roots[0];
            let mut best_g = f64::INFINITY;
            for v in roots {
                let rho = 1.0 / v;
                let a_res = a_res_at_rho(rho)?;
                let z = pressure * v / (crate::R * temperature);
                let g_reduced = a_res + (z - 1.0) - z.ln();
                if g_reduced < best_g {
                    best_g = g_reduced;
                    best_v = v;
                }
            }
            Ok(best_v)
        }
    }
}

/// Exact port of `src/saft_interface.f90:pc_saft_nonassoc_volume_solver` for
/// pure fluids (`n = 1 mol`).
pub(crate) fn baseline_saft_volume_iteration<F>(
    temperature: f64,
    pressure: f64,
    phase: Phase,
    eta_per_rho: f64,
    p_and_dp_drho: F,
) -> Result<LjsVolumeIteration>
where
    F: Fn(f64) -> Result<(f64, f64)>,
{
    positive(temperature, "temperature")?;
    positive(pressure, "pressure")?;
    positive(eta_per_rho, "packing fraction factor")?;

    let conv_num = eta_per_rho;
    let mut zeta = if matches!(phase, Phase::Vapor) {
        1.0e-10
    } else {
        0.5
    };
    let mut v = conv_num / zeta;
    let mut zeta_min = 1.0e-10;
    let mut zeta_max = 0.75;

    let (mut p, dp_drho) = p_and_dp_drho(1.0 / v)?;
    let mut p_v = -dp_drho / (v * v);
    let maxiter = 40usize;
    let mut iter = 0usize;
    let mut f = (1.0 - zeta) * (p - pressure);
    let mut v_converged;
    let mut p_converged;

    loop {
        let df_dzeta = pressure - p - ((1.0 - zeta) / zeta) * p_v * v;
        zeta -= f / df_dzeta;
        if !zeta.is_finite() || zeta > zeta_max || zeta < zeta_min {
            zeta = 0.5 * (zeta_max + zeta_min);
        }
        let v_old = v;
        v = conv_num / zeta;
        let (p_new, dp_drho_new) = p_and_dp_drho(1.0 / v)?;
        p = p_new;
        p_v = -dp_drho_new / (v * v);
        f = (1.0 - zeta) * (p - pressure);

        v_converged = (v - v_old).abs() < 1.0e5 * f64::EPSILON * v_old;
        p_converged = (p - pressure).abs() < 1.0e8 * f64::EPSILON * pressure;
        if v_converged && p_converged {
            break;
        }
        iter += 1;
        if v_converged || iter > maxiter {
            break;
        }
        if f > 0.0 {
            zeta_max = zeta;
        } else if f < 0.0 {
            zeta_min = zeta;
        }
    }

    let volume = conv_num / zeta;
    Ok(LjsVolumeIteration {
        volume,
        iterations: iter,
        volume_converged: v_converged,
        pressure_converged: p_converged,
        relative_pressure_error: (p - pressure).abs() / pressure,
    })
}

/// Estimated critical temperature (K) and critical molar volume (m^3/mol)
/// for splined Lennard-Jones models (`src/saft_interface.f90:estimate_critical_parameters`).
pub fn estimate_ljs_critical_parameters(sigma: f64, eps_divk: f64) -> Result<(f64, f64)> {
    positive(sigma, "segment diameter")?;
    positive(eps_divk, "well depth")?;
    let tc = 0.885 * eps_divk;
    let vc = sigma * sigma * sigma * crate::N_A / 0.333;
    Ok((tc, vc))
}

/// Reduced splined Lennard-Jones interaction potential `u/eps` and reduced
/// force `-du*/dr*` (`src/lj_splined.f90:ljs_potential_reduced`).
pub fn ljs_potential_reduced(separations: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let xs = (26.0_f64 / 7.0).powf(1.0 / 6.0);
    let xc = 67.0 / 48.0 * xs;
    let a = -24192.0 / 3211.0 / (xs * xs);
    let b = -387072.0 / 61009.0 / (xs * xs * xs);
    let mut potential = Vec::with_capacity(separations.len());
    let mut force = Vec::with_capacity(separations.len());
    for &separation in separations {
        if separation < xs {
            let inverse = 1.0 / crate::numerics::ipow(separation, 6);
            potential.push(4.0 * (inverse * inverse - inverse));
            force.push((inverse / separation) * (48.0 * inverse - 24.0));
        } else if separation < xc {
            let difference = separation - xc;
            potential.push(a * difference * difference + b * difference * difference * difference);
            force.push(-(2.0 * a * difference + 3.0 * b * difference * difference));
        } else {
            potential.push(0.0);
            force.push(0.0);
        }
    }
    (potential, force)
}

/// Pair potential `u/k_B` (K) and radial force `-du/(k_B dr)` (K/m) for the
/// splined Lennard-Jones potential (`src/saft_interface.f90:potential`).
pub fn ljs_pair_potential(
    sigma: f64,
    eps_divk: f64,
    separations: &[f64],
    max_pot_val: f64,
) -> Result<(Vec<f64>, Vec<f64>)> {
    positive(sigma, "segment diameter")?;
    positive(eps_divk, "well depth")?;
    if !max_pot_val.is_finite() || max_pot_val <= 0.0 {
        return Err(Error::InvalidInput("positive max potential cap required"));
    }
    let reduced: Vec<f64> = separations.iter().map(|&r| r / sigma).collect();
    let (mut pot, mut force) = ljs_potential_reduced(&reduced);
    let cap = max_pot_val * eps_divk;
    for i in 0..separations.len() {
        pot[i] *= eps_divk;
        force[i] *= eps_divk / sigma;
        if pot[i].is_nan() || pot[i] > cap {
            pot[i] = cap;
        }
    }
    Ok((pot, force))
}

/// The Lennard-Jones splined model with the Barker-Henderson perturbation
/// theory (`LJS-BH`).
#[derive(Debug, Clone, Copy)]
pub struct LjsBh {
    sigma: f64,
    eps: f64,
    control: LjsControl,
    kolafa_pure_reference: bool,
}

/// The terms of the `LJS-BH` model that the original lets the user switch on and off.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LjsControl {
    /// Whether the hard-sphere reference is included.
    pub enable_hard_sphere: bool,
    /// Whether the first order perturbation term is included.
    pub enable_a1: bool,
    /// Whether the second order perturbation term is included.
    pub enable_a2: bool,
    /// Whether the third order perturbation term is included.
    pub enable_a3: bool,
    /// Whether the correction of the second order term is included.
    pub enable_chi_correction: bool,
    /// Whether the Lafitte form of the third order term is used.
    pub use_lafitte_a3: bool,
}

impl Default for LjsControl {
    fn default() -> Self {
        Self {
            enable_hard_sphere: true,
            enable_a1: true,
            enable_a2: true,
            enable_a3: true,
            enable_chi_correction: true,
            use_lafitte_a3: false,
        }
    }
}

impl Default for LjsBh {
    fn default() -> Self {
        Self::new()
    }
}

/// Temperature-prepared `LJS-BH` state caching the Barker-Henderson diameter,
/// temperature derivatives, packing-fraction factors, and Lafitte coefficients.
#[derive(Debug, Clone, Copy)]
pub struct PreparedLjsBh {
    eps: f64,
    control: LjsControl,
    kolafa_pure_reference: bool,
    temperature: f64,
    diameter: f64,
    diameter_t: f64,
    diameter_tt: f64,
    tau: f64,
    tau_t: f64,
    tau_tt: f64,
    eta_per_rho: f64,
    eta_per_rho_t: f64,
    eta_per_rho_tt: f64,
    rho_star_per_rho: f64,
    zeta_bar_per_rho: f64,
    f_alpha: [f64; 6],
}

impl LjsBh {
    /// The model of the source database, the argon record.
    pub fn new() -> Self {
        Self {
            sigma: ARGON_SIGMA,
            eps: ARGON_EPS,
            control: LjsControl::default(),
            kolafa_pure_reference: false,
        }
    }

    /// The model with the segment diameter and the well depth of the user.
    pub fn with_parameters(sigma: f64, eps: f64) -> Self {
        Self {
            sigma,
            eps,
            control: LjsControl::default(),
            kolafa_pure_reference: false,
        }
    }

    /// Replace `(sigma, eps/k)` consistently across all internal derived quantities
    /// (`ljs_bh_set_sigma_eps`).
    pub fn set_parameters(&mut self, sigma: f64, eps: f64) -> Result<()> {
        positive(sigma, "segment diameter")?;
        positive(eps, "well depth")?;
        self.sigma = sigma;
        self.eps = eps;
        Ok(())
    }

    /// The terms of the model that the user switches.
    pub fn control(&self) -> LjsControl {
        self.control
    }

    /// Switches the terms of the model, as the original does.
    pub fn set_control(&mut self, control: LjsControl) {
        self.control = control;
    }

    /// Whether the Carnahan-Starling-Kolafa (`PURE_HS_CSK`) hard-sphere reference is used.
    pub fn kolafa_pure_reference(&self) -> bool {
        self.kolafa_pure_reference
    }

    /// Switch between Carnahan-Starling (`false`, default) and Carnahan-Starling-Kolafa (`true`).
    pub fn set_kolafa_pure_reference(&mut self, kolafa: bool) {
        self.kolafa_pure_reference = kolafa;
    }

    /// The segment diameter of the model in metre.
    pub fn sigma(&self) -> f64 {
        self.sigma
    }

    /// The well depth of the model divided by the Boltzmann constant in kelvin.
    pub fn eps(&self) -> f64 {
        self.eps
    }

    /// The Barker-Henderson hard-sphere diameter of the model and its first
    /// and second temperature derivatives in metre.
    pub fn diameter(&self, temperature: f64) -> Result<(f64, f64, f64)> {
        crate::saft_vrmie::calc_hardsphere_diameter_pure(
            temperature,
            self.sigma,
            self.eps,
            LAMBDA_A,
            LAMBDA_R,
            C_FAC,
            1.0,
            0,
        )
    }

    /// The Barker-Henderson hard-sphere diameter divided by the segment
    /// diameter at a reduced temperature.
    pub fn diameter_div_sigma(&self, reduced_temperature: f64) -> Result<f64> {
        let temperature = reduced_temperature * self.eps;
        let (diameter, _, _) = self.diameter(temperature)?;
        Ok(diameter / self.sigma)
    }

    /// Prepare temperature-only quantities for repeated density/volume evaluations.
    pub fn prepare(&self, temperature: f64) -> Result<PreparedLjsBh> {
        positive(self.sigma, "segment diameter")?;
        positive(self.eps, "well depth")?;
        positive(temperature, "temperature")?;
        let (d, d_t, d_tt) = self.diameter(temperature)?;
        let (eta_per_rho, eta_per_rho_t, eta_per_rho_tt) = pure_eta_per_rho(d, d_t, d_tt);
        let rho_star_per_rho = crate::N_A * self.sigma * self.sigma * self.sigma;
        let zeta_bar_per_rho = PI / 6.0 * rho_star_per_rho;
        let f_alpha = crate::saft_vrmie::alpha_coefficients(ALPHA_LJS);
        Ok(PreparedLjsBh {
            eps: self.eps,
            control: self.control,
            kolafa_pure_reference: self.kolafa_pure_reference,
            temperature,
            diameter: d,
            diameter_t: d_t,
            diameter_tt: d_tt,
            tau: 1.0 - d / self.sigma,
            tau_t: -d_t / self.sigma,
            tau_tt: -d_tt / self.sigma,
            eta_per_rho,
            eta_per_rho_t,
            eta_per_rho_tt,
            rho_star_per_rho,
            zeta_bar_per_rho,
            f_alpha,
        })
    }

    /// The reduced residual chemical potential `F_n = ln(phi) + ln(Z)` of the model
    /// (preserved for backward compatibility; see `LjsTvState::ln_fugacity_coefficient`
    /// for `ln(phi)`).
    pub fn fugacity_coefficient(&self, temperature: f64, volume: f64, moles: f64) -> Result<f64> {
        let f = self.residual_helmholtz_tvn(temperature, volume, moles)?;
        Ok(f.dn)
    }

    /// The packing fraction of the model at a temperature and a density.
    pub fn packing_fraction(&self, temperature: f64, density: f64) -> Result<f64> {
        positive(density, "the LJS density")?;
        let (diameter, _, _) = self.diameter(temperature)?;
        Ok(PI / 6.0 * crate::N_A * (diameter * diameter * diameter) * density)
    }

    /// Packing fraction and all `(T, V, n)` derivatives.
    pub fn packing_fraction_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<PackingFractionDerivatives> {
        let (d, d_t, d_tt) = self.diameter(temperature)?;
        packing_fraction_vdw_no_segments(volume, &[moles], &[d], &[d_t], &[d_tt])
    }

    /// The reduced interaction potential of the model and the reduced force.
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

    /// All mechanically stable molar volumes (m^3/mol) at `(T, P)`, sorted ascending.
    pub fn volumes(&self, temperature: f64, pressure: f64) -> Result<Vec<f64>> {
        self.prepare(temperature)?.volumes(pressure)
    }

    /// The physically stable molar volume (m^3/mol) of the model at `(T, P)`.
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

    /// Full `(T, V, n)` thermodynamic state (caloric properties, speed of sound, fugacity).
    pub fn tv_state(&self, temperature: f64, volume: f64, moles: f64) -> Result<LjsTvState> {
        self.prepare(temperature)?.tv_state(volume, moles)
    }

    /// Legacy summary state of the model at a temperature and a volume.
    pub fn state(&self, temperature: f64, volume: f64, moles: f64) -> Result<LjsState> {
        Ok(self.tv_state(temperature, volume, moles)?.to_legacy_state())
    }

    /// Speed of sound (m/s) at `(T, V, n)`.
    pub fn speed_of_sound_tv(&self, temperature: f64, volume: f64, moles: f64) -> Result<f64> {
        Ok(self.tv_state(temperature, volume, moles)?.speed_of_sound)
    }

    /// Full `(T, P)` thermodynamic state and `(T, P)` derivatives for `phase`.
    pub fn state_tp(&self, temperature: f64, pressure: f64, phase: Phase) -> Result<LjsTpState> {
        self.prepare(temperature)?.state_tp(pressure, phase)
    }

    /// The pressure (Pa) and compressibility factor of the model at `(T, V, n)`.
    pub fn pressure(&self, temperature: f64, volume: f64, moles: f64) -> Result<(f64, f64)> {
        self.prepare(temperature)?.pressure(volume, moles)
    }

    /// Pressure and first `(T, V, n)` derivatives at `(T, V, n)`.
    pub fn pressure_tv(&self, temperature: f64, volume: f64, moles: f64) -> Result<PurePressure> {
        self.prepare(temperature)?.pressure_tv(volume, moles)
    }

    /// Molar reduced residual Helmholtz energy `a_res(T, rho)` and all 6 `(T, rho)` derivatives.
    pub fn residual_density(&self, temperature: f64, rho: f64) -> Result<DensityHelmholtz> {
        self.prepare(temperature)?.residual_density(rho)
    }

    /// Extensive reduced residual Helmholtz energy `F = n * a_res` and all 10 `(T, V, n)`
    /// derivatives (`calcFresLJs_bh`).
    pub fn residual_helmholtz_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<PureHelmholtz> {
        self.prepare(temperature)?
            .residual_helmholtz_tvn(volume, moles)
    }

    /// Molar first-order term `a1` (before dividing by `T`) and all 10 `(T, V, n)` derivatives
    /// (`calc_LJs_a1`).
    pub fn a1_tvn(&self, temperature: f64, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        self.prepare(temperature)?.a1_tvn(volume, moles)
    }

    /// Molar second-order term `a2` (before dividing by `T^2`) and all 10 `(T, V, n)` derivatives
    /// (`calc_LJs_a2`).
    pub fn a2_tvn(&self, temperature: f64, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        self.prepare(temperature)?.a2_tvn(volume, moles)
    }

    /// Molar third-order term `a3` (before dividing by `T^3`) and all 10 `(T, V, n)` derivatives
    /// (`calc_LJs_a3`).
    pub fn a3_tvn(&self, temperature: f64, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        self.prepare(temperature)?.a3_tvn(volume, moles)
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

    /// Molar reduced dispersion contribution and all 10 `(T, V, n)` derivatives
    /// (`calc_ljs_dispersion`).
    pub fn dispersion_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<PureHelmholtz> {
        self.prepare(temperature)?.dispersion_tvn(volume, moles)
    }

    /// The molar reduced residual Helmholtz energy of the model and its
    /// temperature derivative.
    pub fn residual_helmholtz_energy_derivative(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<(f64, f64)> {
        positive(volume, "the LJS volume")?;
        positive(moles, "the LJS number of moles")?;
        let dh = self.residual_density(temperature, moles / volume)?;
        Ok((dh.value, dh.dt))
    }

    /// The molar reduced residual Helmholtz energy of the model.
    pub fn residual_helmholtz_energy(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<f64> {
        Ok(self
            .residual_helmholtz_energy_derivative(temperature, volume, moles)?
            .0)
    }

    /// The reduced perturbation terms of the model (`calc_ai_reduced_LJs_ex`).
    pub fn reduced_perturbation_terms(
        &self,
        reduced_temperature: f64,
        reduced_density: f64,
    ) -> Result<(f64, f64, f64)> {
        let temperature = reduced_temperature * self.eps;
        let (diameter, _, _) = self.diameter(temperature)?;
        let x0 = self.sigma / diameter;
        let tau = 1.0 - 1.0 / x0;
        let first = 2.0 * PI * polynomial(0, tau, reduced_density);
        let second = -PI * polynomial(1, tau, reduced_density);
        let correction = polynomial(3, tau, reduced_density);
        let eta = PI * reduced_density / (6.0 * x0 * x0 * x0);
        let compressibility = khs(eta);
        let second = second * (1.0 + correction) * compressibility;
        let third = polynomial(2, tau, reduced_density);
        Ok((first, second, third))
    }

    /// One polynomial of the spline of the excess energy and its partial derivatives.
    pub fn spline_polynomial(
        &self,
        index: usize,
        reduced_temperature: f64,
        reduced_density: f64,
    ) -> Result<(f64, f64, f64, f64, f64, f64)> {
        let temperature = reduced_temperature * self.eps;
        let (diameter, _, _) = self.diameter(temperature)?;
        let x0 = self.sigma / diameter;
        let tau = 1.0 - 1.0 / x0;
        Ok(polynomial_with_derivatives(
            index,
            tau,
            reduced_density,
            1.0,
        ))
    }

    /// Estimated critical temperature (K) and critical molar volume (m^3/mol).
    pub fn estimate_critical_parameters(&self) -> Result<(f64, f64)> {
        estimate_ljs_critical_parameters(self.sigma, self.eps)
    }

    /// Fundamental Measure Theory compatibility flags `(is_fmt_consistent, na_enabled)`.
    pub fn fmt_compatibility(&self) -> (bool, bool) {
        (true, false)
    }

    /// Dimensionless van der Waals attractive energy constant (`alpha_ljs`).
    pub fn alpha_vdw(&self) -> f64 {
        ALPHA_LJS
    }
}

impl PreparedLjsBh {
    /// Temperature (K) at which this state was prepared.
    pub fn temperature(&self) -> f64 {
        self.temperature
    }

    /// Cached Barker-Henderson diameter and its first two temperature derivatives.
    pub fn diameter(&self) -> (f64, f64, f64) {
        (self.diameter, self.diameter_t, self.diameter_tt)
    }

    /// Cached molar packing-fraction factor `eta / rho` (m^3/mol).
    pub fn eta_per_rho(&self) -> f64 {
        self.eta_per_rho
    }

    fn poly_density(&self, index: usize, prefactor: f64, rho: f64) -> DensityHelmholtz {
        let scale = self.rho_star_per_rho;
        let rho_star = scale * rho;
        let (val, val_r, val_x, val_rr, val_xx, val_rx) =
            polynomial_with_derivatives(index, self.tau, rho_star, prefactor);
        DensityHelmholtz {
            value: val,
            dt: val_x * self.tau_t,
            drho: val_r * scale,
            dtt: val_xx * self.tau_t * self.tau_t + val_x * self.tau_tt,
            dtrho: val_rx * self.tau_t * scale,
            drhorho: val_rr * scale * scale,
        }
    }

    fn khs_density(&self, rho: f64) -> DensityHelmholtz {
        let eta = self.eta_per_rho * rho;
        let eta_t = self.eta_per_rho_t * rho;
        let eta_tt = self.eta_per_rho_tt * rho;
        let (k, k_e, k_ee) = khs_derivatives(eta);
        DensityHelmholtz {
            value: k,
            dt: k_e * eta_t,
            drho: k_e * self.eta_per_rho,
            dtt: k_e * eta_tt + k_ee * eta_t * eta_t,
            dtrho: k_e * self.eta_per_rho_t + k_ee * eta_t * self.eta_per_rho,
            drhorho: k_ee * self.eta_per_rho * self.eta_per_rho,
        }
    }

    fn a1_density(&self, rho: f64) -> DensityHelmholtz {
        self.poly_density(0, 2.0 * PI * self.eps, rho)
    }

    fn a2_density(&self, rho: f64) -> DensityHelmholtz {
        let mut a2 = self.poly_density(1, -PI * self.eps * self.eps, rho);
        if self.control.enable_chi_correction {
            let mut one_plus_chi = self.poly_density(3, 1.0, rho);
            one_plus_chi.value += 1.0;
            a2 = a2.mul(one_plus_chi);
        }
        let k = self.khs_density(rho);
        a2.mul(k)
    }

    fn a3_density(&self, rho: f64) -> DensityHelmholtz {
        let eps3 = self.eps * self.eps * self.eps;
        if self.control.use_lafitte_a3 {
            let z_scale = self.zeta_bar_per_rho;
            let z = z_scale * rho;
            let f4 = self.f_alpha[3];
            let f5 = self.f_alpha[4];
            let f6 = self.f_alpha[5];
            let exp_term = (f5 * z + f6 * z * z).exp();
            let u_p = f5 + 2.0 * f6 * z;
            let u_pp = 2.0 * f6;
            let val = -eps3 * f4 * z * exp_term;
            let dz = -eps3 * f4 * exp_term * (1.0 + z * u_p);
            let dzz = -eps3 * f4 * exp_term * (2.0 * u_p + z * (u_p * u_p + u_pp));
            DensityHelmholtz {
                value: val,
                dt: 0.0,
                drho: dz * z_scale,
                dtt: 0.0,
                dtrho: 0.0,
                drhorho: dzz * z_scale * z_scale,
            }
        } else {
            self.poly_density(2, eps3, rho)
        }
    }

    fn dispersion_density(&self, rho: f64) -> DensityHelmholtz {
        let mut total = DensityHelmholtz::zero();
        if self.control.enable_a1 {
            total = total.add(self.a1_density(rho).scale_inv_t_power(self.temperature, 1));
        }
        if self.control.enable_a2 {
            total = total.add(self.a2_density(rho).scale_inv_t_power(self.temperature, 2));
        }
        if self.control.enable_a3 {
            total = total.add(self.a3_density(rho).scale_inv_t_power(self.temperature, 3));
        }
        total
    }

    fn density_only(&self, rho: f64) -> Result<(f64, f64, f64)> {
        positive(rho, "molar density")?;
        let mut a = 0.0;
        let mut a_r = 0.0;
        let mut a_rr = 0.0;
        let eta = self.eta_per_rho * rho;
        if self.control.enable_hard_sphere {
            let (hs, hs_e, hs_ee) = pure_hard_sphere_eta(eta, self.kolafa_pure_reference)?;
            a += hs;
            a_r += hs_e * self.eta_per_rho;
            a_rr += hs_ee * self.eta_per_rho * self.eta_per_rho;
        }
        let scale = self.rho_star_per_rho;
        let scale2 = scale * scale;
        let rho_star = scale * rho;
        let inv_t = 1.0 / self.temperature;
        if self.control.enable_a1 {
            let (v, vr, _, vrr, _, _) =
                polynomial_with_derivatives(0, self.tau, rho_star, 2.0 * PI * self.eps);
            a += inv_t * v;
            a_r += inv_t * vr * scale;
            a_rr += inv_t * vrr * scale2;
        }
        if self.control.enable_a2 {
            let (mut v, mut vr, _, mut vrr, _, _) =
                polynomial_with_derivatives(1, self.tau, rho_star, -PI * self.eps * self.eps);
            vr *= scale;
            vrr *= scale2;
            if self.control.enable_chi_correction {
                let (c, cr, _, crr, _, _) = polynomial_with_derivatives(3, self.tau, rho_star, 1.0);
                let g = 1.0 + c;
                let gr = cr * scale;
                let grr = crr * scale2;
                vrr = vrr * g + 2.0 * vr * gr + v * grr;
                vr = vr * g + v * gr;
                v *= g;
            }
            let (k, ke, kee) = khs_derivatives(eta);
            let kr = ke * self.eta_per_rho;
            let krr = kee * self.eta_per_rho * self.eta_per_rho;
            let inv_t2 = inv_t * inv_t;
            a += inv_t2 * (v * k);
            a_r += inv_t2 * (vr * k + v * kr);
            a_rr += inv_t2 * (vrr * k + 2.0 * vr * kr + v * krr);
        }
        if self.control.enable_a3 {
            let eps3 = self.eps * self.eps * self.eps;
            let inv_t3 = inv_t * inv_t * inv_t;
            if self.control.use_lafitte_a3 {
                let z_scale = self.zeta_bar_per_rho;
                let z = z_scale * rho;
                let f4 = self.f_alpha[3];
                let f5 = self.f_alpha[4];
                let f6 = self.f_alpha[5];
                let exp_term = (f5 * z + f6 * z * z).exp();
                let u_p = f5 + 2.0 * f6 * z;
                let u_pp = 2.0 * f6;
                a += inv_t3 * (-eps3 * f4 * z * exp_term);
                a_r += inv_t3 * (-eps3 * f4 * exp_term * (1.0 + z * u_p) * z_scale);
                a_rr += inv_t3
                    * (-eps3
                        * f4
                        * exp_term
                        * (2.0 * u_p + z * (u_p * u_p + u_pp))
                        * z_scale
                        * z_scale);
            } else {
                let (v, vr, _, vrr, _, _) =
                    polynomial_with_derivatives(2, self.tau, rho_star, eps3);
                a += inv_t3 * v;
                a_r += inv_t3 * vr * scale;
                a_rr += inv_t3 * vrr * scale2;
            }
        }
        if !a.is_finite() || !a_r.is_finite() || !a_rr.is_finite() {
            return Err(Error::InvalidInput("non-finite LJS-BH density evaluation"));
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
                self.kolafa_pure_reference,
            )?;
            total = total.add(hs);
        }
        total = total.add(self.dispersion_density(rho));
        Ok(total)
    }

    /// Extensive reduced residual Helmholtz energy `F = n * a_res` and all 10 `(T, V, n)` derivatives.
    pub fn residual_helmholtz_tvn(&self, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        self.residual_density(moles / volume)?
            .extensive(volume, moles)
    }

    /// Molar `a1` (before dividing by `T`) and all 10 `(T, V, n)` derivatives.
    pub fn a1_tvn(&self, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        self.a1_density(moles / volume).molar_tvn(volume, moles)
    }

    /// Molar `a2` (before dividing by `T^2`) and all 10 `(T, V, n)` derivatives.
    pub fn a2_tvn(&self, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        self.a2_density(moles / volume).molar_tvn(volume, moles)
    }

    /// Molar `a3` (before dividing by `T^3`) and all 10 `(T, V, n)` derivatives.
    pub fn a3_tvn(&self, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        self.a3_density(moles / volume).molar_tvn(volume, moles)
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
            self.kolafa_pure_reference,
        )?;
        hs.molar_tvn(volume, moles)
    }

    /// Molar reduced dispersion contribution and all 10 `(T, V, n)` derivatives.
    pub fn dispersion_tvn(&self, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        self.dispersion_density(moles / volume)
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

/// One polynomial of the spline of the excess energy.
fn polynomial(index: usize, tau: f64, rho: f64) -> f64 {
    polynomial_with_derivatives(index, tau, rho, 1.0).0
}

/// One polynomial of the spline of the excess energy with its partial
/// derivatives and the prefactor of the source.
fn polynomial_with_derivatives(
    index: usize,
    tau: f64,
    rho: f64,
    prefactor: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    let mut powers = [0.0; 5];
    powers[4] = rho;
    for position in (0..4).rev() {
        powers[position] = powers[position + 1] * rho;
    }
    let mut first = [0.0; 5];
    first[4] = 1.0;
    first[3] = 2.0 * rho;
    first[2] = 3.0 * powers[3];
    first[1] = 4.0 * powers[2];
    first[0] = 5.0 * powers[1];
    let mut second = [0.0; 5];
    second[4] = 0.0;
    second[3] = 2.0;
    second[2] = 6.0 * rho;
    second[1] = 12.0 * powers[3];
    second[0] = 20.0 * powers[2];
    let g = [
        [1.0, 0.0, 0.0],
        [tau, 1.0, 0.0],
        [tau * tau, 2.0 * tau, 2.0],
    ];
    let mut value = 0.0;
    let mut value_x = 0.0;
    let mut value_r = 0.0;
    let mut value_xx = 0.0;
    let mut value_rr = 0.0;
    let mut value_rx = 0.0;
    for (position, block) in [0usize, 5, 10].into_iter().enumerate() {
        let mut sum = 0.0;
        let mut sum_first = 0.0;
        let mut sum_second = 0.0;
        for column in 0..5 {
            let coefficient = POLYNOMIALS[(block + column) * 4 + index];
            sum += coefficient * powers[column];
            sum_first += coefficient * first[column];
            sum_second += coefficient * second[column];
        }
        value += sum * g[position][0];
        value_x += sum * g[position][1];
        value_r += sum_first * g[position][0];
        value_xx += sum * g[position][2];
        value_rr += sum_second * g[position][0];
        value_rx += sum_first * g[position][1];
    }
    (
        prefactor * value,
        prefactor * value_r,
        prefactor * value_x,
        prefactor * value_rr,
        prefactor * value_xx,
        prefactor * value_rx,
    )
}

/// Public Wilhelmsen spline polynomial `calc_ai_LJs_ex` with 1-based index `1..=4`.
pub fn ai_ljs_ex(
    index: usize,
    tau: f64,
    rho: f64,
    prefactor: f64,
) -> Result<(f64, f64, f64, f64, f64, f64)> {
    if !(1..=4).contains(&index) {
        return Err(Error::InvalidInput(
            "the Wilhelmsen polynomial index must be one to four",
        ));
    }
    Ok(polynomial_with_derivatives(index - 1, tau, rho, prefactor))
}

/// First-order BH perturbation term `calc_a1_LJs_ex(x0, rho, eps)`.
pub fn a1_ljs_ex(x0: f64, rho: f64, eps: f64) -> Result<f64> {
    positive(x0, "reduced distance x0")?;
    let tau = 1.0 - 1.0 / x0;
    Ok(polynomial_with_derivatives(0, tau, rho, 2.0 * PI * eps).0)
}

/// Second-order BH perturbation term `calc_a2_LJs_ex(x0, rho, eps)`.
pub fn a2_ljs_ex(x0: f64, rho: f64, eps: f64) -> Result<f64> {
    positive(x0, "reduced distance x0")?;
    let tau = 1.0 - 1.0 / x0;
    let a2 = polynomial_with_derivatives(1, tau, rho, -PI * eps * eps).0;
    let chi = polynomial_with_derivatives(3, tau, rho, 1.0).0;
    let eta = PI * rho / (6.0 * x0 * x0 * x0);
    Ok(a2 * (1.0 + chi) * khs(eta))
}

/// Third-order BH perturbation term `calc_a3_LJs_ex(x0, rho, eps)`.
pub fn a3_ljs_ex(x0: f64, rho: f64, eps: f64) -> Result<f64> {
    positive(x0, "reduced distance x0")?;
    let tau = 1.0 - 1.0 / x0;
    Ok(polynomial_with_derivatives(2, tau, rho, eps * eps * eps).0)
}

/// The isothermal compressibility factor of the hard-sphere reference, the
/// Carnahan-Starling expression of the original.
fn khs(eta: f64) -> f64 {
    let eta2 = eta * eta;
    let denominator = 1.0 + 4.0 * eta + 4.0 * eta2 - 4.0 * eta2 * eta + eta2 * eta2;
    let om = 1.0 - eta;
    let om2 = om * om;
    (om2 * om2) / denominator
}

/// Isothermal compressibility factor `Khs(eta)` and its first and second
/// packing-fraction derivatives.
fn khs_derivatives(eta: f64) -> (f64, f64, f64) {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let d = 1.0 + 4.0 * eta + 4.0 * eta2 - 4.0 * eta3 + eta2 * eta2;
    let d_e = 4.0 + 8.0 * eta - 12.0 * eta2 + 4.0 * eta3;
    let d_ee = 8.0 - 24.0 * eta + 12.0 * eta2;
    let om = 1.0 - eta;
    let om2 = om * om;
    let om3 = om2 * om;
    let n = om2 * om2;
    let n_ee = 12.0 * om2;
    let k = n / d;
    let k_e = 4.0 * (eta2 - 5.0 * eta - 2.0) * om3 / (d * d);
    let k_ee = (n_ee - 2.0 * k_e * d_e - k * d_ee) / d;
    (k, k_e, k_ee)
}

/// The u-fraction of the UV theory and its partial derivatives (`uv_phi`).
pub fn uv_phi(_eps_divk: f64, _temperature: f64, rho: f64) -> (f64, f64, f64, f64, f64, f64) {
    const C: [f64; 5] = [0.11072527, 2.5809656, 1.2333385, 2.3534453, 4.1241290];
    let x = C[0] * rho + C[1] * rho.powf(C[2]) + C[3] * rho.powf(C[4]);
    let x_r = C[0] + C[1] * C[2] * rho.powf(C[2] - 1.0) + C[3] * C[4] * rho.powf(C[4] - 1.0);
    let x_rr = C[1] * C[2] * (C[2] - 1.0) * rho.powf(C[2] - 2.0)
        + C[3] * C[4] * (C[4] - 1.0) * rho.powf(C[4] - 2.0);
    let (tanh, tanh_r, tanh_rr) = safe_tanh_chain(x, x_r, x_rr);
    (tanh, tanh_r, 0.0, tanh_rr, 0.0, 0.0)
}

/// The temperature dependent u-fraction of the UV theory and its partial
/// derivatives (`uv_phi_temp`).
pub fn uv_phi_temp(eps_divk: f64, temperature: f64, rho: f64) -> (f64, f64, f64, f64, f64, f64) {
    const C: [f64; 7] = [
        1.01236, 0.99164, 2.63226, 1.21219, 2.25005, 3.91411, 0.436593,
    ];
    let beta = eps_divk / temperature;
    let beta_t = -eps_divk / (temperature * temperature);
    let beta_tt = 2.0 * eps_divk / (temperature * temperature * temperature);
    let sqrt_c7 = C[6].sqrt();
    let denominator_sq = 1.0 + C[6] * beta * beta;
    let denominator = denominator_sq.sqrt();
    let sigma_beta = sqrt_c7 * beta / denominator;
    let sigma_beta_b =
        sqrt_c7 / denominator - sqrt_c7 * C[6] * beta * beta / (denominator * denominator_sq);
    let sigma_beta_bb = -3.0 * sqrt_c7 * C[6] * beta / (denominator * denominator_sq)
        + 3.0 * sqrt_c7 * C[6] * C[6] * (beta * beta * beta)
            / (denominator * denominator_sq * denominator_sq);
    let prefactor = C[0] + sigma_beta * (C[1] - C[0]);
    let x = C[2] * rho.powf(C[3]) + C[4] * rho.powf(C[5]);
    let x_r = C[2] * C[3] * rho.powf(C[3] - 1.0) + C[4] * C[5] * rho.powf(C[5] - 1.0);
    let x_rr = C[2] * C[3] * (C[3] - 1.0) * rho.powf(C[3] - 2.0)
        + C[4] * C[5] * (C[5] - 1.0) * rho.powf(C[5] - 2.0);
    let (tanh, tanh_r, tanh_rr) = safe_tanh_chain(x, x_r, x_rr);
    let value = prefactor * tanh;
    let value_b = sigma_beta_b * (C[1] - C[0]) * tanh;
    let value_bb = sigma_beta_bb * (C[1] - C[0]) * tanh;
    let value_r = prefactor * tanh_r;
    let value_rr = prefactor * tanh_rr;
    let value_rb = sigma_beta_b * (C[1] - C[0]) * tanh_r;
    (
        value,
        value_r,
        value_b * beta_t,
        value_rr,
        value_bb * beta_t * beta_t + value_b * beta_tt,
        value_rb * beta_t,
    )
}

/// The first order term of the UV theory and its partial derivatives (`uv_a1_u`).
pub fn uv_a1_u(
    sigma: f64,
    eps_divk: f64,
    temperature: f64,
    diameter: f64,
    diameter_t: f64,
    diameter_tt: f64,
    rho: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    let d = diameter / sigma;
    let d_t = diameter_t / sigma;
    let d_tt = diameter_tt / sigma;
    let (a1, a1_r, a1_d, a1_rr, a1_dd, a1_rd, _a1_rdd) = uv_a1_u_mult_t(eps_divk, d, rho);
    let a1_t = a1_d * d_t / temperature - a1 / (temperature * temperature);
    let a1_rt = a1_rd * d_t / temperature - a1_r / (temperature * temperature);
    let a1_tt = a1_dd * d_t * d_t / temperature + a1_d * d_tt / temperature
        - 2.0 * a1_d * d_t / (temperature * temperature)
        + 2.0 * a1 / (temperature * temperature * temperature);
    (
        a1 / temperature,
        a1_r / temperature,
        a1_t,
        a1_rr / temperature,
        a1_tt,
        a1_rt,
    )
}

/// The second virial contribution of the UV `a1` term (`uv_a1_b2`).
pub fn uv_a1_b2(
    sigma: f64,
    eps_divk: f64,
    temperature: f64,
    diameter: f64,
    diameter_t: f64,
    diameter_tt: f64,
) -> (f64, f64, f64) {
    let d = diameter / sigma;
    let d_t = diameter_t / sigma;
    let d_tt = diameter_tt / sigma;
    let (_a1, a1_r, _a1_d, _a1_rr, _a1_dd, a1_rd, a1_rdd) = uv_a1_u_mult_t(eps_divk, d, 0.0);
    let b2 = a1_r / temperature;
    let b2_t = a1_rd * d_t / temperature - a1_r / (temperature * temperature);
    let b2_tt = a1_rd * d_tt / temperature + a1_rdd * d_t * d_t / temperature
        - 2.0 * a1_rd * d_t / (temperature * temperature)
        + 2.0 * a1_r / (temperature * temperature * temperature);
    (b2, b2_t, b2_tt)
}

/// The UV `a1 * T` polynomial and its `(rho, d)` derivatives (`uv_a1_u_mult_T`).
pub fn uv_a1_u_mult_t(eps_divk: f64, d: f64, rho: f64) -> (f64, f64, f64, f64, f64, f64, f64) {
    const C: [f64; 8] = [
        -0.8513603, 0.2977649, -0.1287404, -0.3377645, 0.7343770, -0.9839062, -0.4592915, 0.7963438,
    ];
    let d3 = d * d * d;
    let a1_div_r = 2.0
        * PI
        * ((C[0] + C[1] * d3)
            + (C[2] + C[3] * d3) * rho
            + (C[4] + C[5] * d3) * rho * rho
            + (C[6] + C[7] * d3) * rho * rho * rho);
    let a1_div_r_r = 2.0
        * PI
        * ((C[2] + C[3] * d3)
            + 2.0 * (C[4] + C[5] * d3) * rho
            + 3.0 * (C[6] + C[7] * d3) * rho * rho);
    let a1 = a1_div_r * rho;
    let a1_r = a1_div_r + a1_div_r_r * rho;
    let a1_rr = 2.0 * a1_div_r_r
        + 2.0 * PI * rho * (2.0 * (C[4] + C[5] * d3) + 6.0 * (C[6] + C[7] * d3) * rho);
    let a1_d =
        6.0 * PI * rho * (C[1] + C[3] * rho + C[5] * rho * rho + C[7] * rho * rho * rho) * d * d;
    let a1_dd =
        12.0 * PI * rho * (C[1] + C[3] * rho + C[5] * rho * rho + C[7] * rho * rho * rho) * d;
    let a1_rd = 6.0
        * PI
        * d
        * d
        * ((C[1] + C[3] * rho + C[5] * rho * rho + C[7] * rho * rho * rho)
            + rho * (C[3] + 2.0 * C[5] * rho + 3.0 * C[7] * rho * rho));
    let a1_rdd = 12.0
        * PI
        * d
        * ((C[1] + C[3] * rho + C[5] * rho * rho + C[7] * rho * rho * rho)
            + rho * (C[3] + 2.0 * C[5] * rho + 3.0 * C[7] * rho * rho));
    (
        a1 * eps_divk,
        a1_r * eps_divk,
        a1_d * eps_divk,
        a1_rr * eps_divk,
        a1_dd * eps_divk,
        a1_rd * eps_divk,
        a1_rdd * eps_divk,
    )
}

/// The q-function of the UV theory and its temperature derivatives (`uv_q_of_T`).
pub fn uv_q_of_temperature(eps_divk: f64, temperature: f64) -> (f64, f64, f64) {
    const M: f64 = 12.0;
    const N: f64 = 6.0;
    const C_MATRIX: [[f64; 3]; 4] = [
        [1.92840364363978, 5.20120816141761e-01, 0.0],
        [
            4.43165896265079e-01,
            1.82526759234412e-01,
            1.29885156087242e-02,
        ],
        [0.0, 1.10319989659929e-02, 6.41039871789327e-03],
        [0.0, -7.97813995328348e-05, 1.85866741090323e-05],
    ];
    let dm = M - 7.0;
    let m_matrix = [1.0, dm, dm * dm, dm * dm * dm];
    let mut coefficients = [0.0; 4];
    coefficients[0] = (2.0 * PI * M / N).sqrt();
    for row in 0..3 {
        let mut sum = 0.0;
        for column in 0..4 {
            sum += C_MATRIX[column][row] * m_matrix[column];
        }
        coefficients[row + 1] = sum;
    }
    let reduced_temperature = temperature / eps_divk;
    let root = reduced_temperature.sqrt();
    let exponent = -1.0 / (2.0 * M);
    let f = 1.0
        + coefficients[0] * root
        + coefficients[1] * reduced_temperature
        + coefficients[2] * reduced_temperature * root
        + coefficients[3] * reduced_temperature * reduced_temperature;
    let f_t = 0.5 * coefficients[0] / root
        + coefficients[1]
        + 1.5 * coefficients[2] * root
        + 2.0 * coefficients[3] * reduced_temperature;
    let f_tt = -0.25 * coefficients[0] / (reduced_temperature * root)
        + 0.75 * coefficients[2] / root
        + 2.0 * coefficients[3];
    let q = f.powf(exponent);
    let q_t = exponent * f.powf(exponent - 1.0) * f_t;
    let q_tt = exponent * (exponent - 1.0) * f.powf(exponent - 2.0) * f_t * f_t
        + exponent * f.powf(exponent - 1.0) * f_tt;
    let scale = 2.0f64.powf(1.0 / 6.0);
    (
        scale * q,
        scale * q_t / eps_divk,
        scale * q_tt / (eps_divk * eps_divk),
    )
}

/// The hard-sphere second virial coefficient of the UV theory and its
/// temperature derivatives (`uv_hs_b2`).
pub fn uv_hard_sphere_b2(eps_divk: f64, temperature: f64) -> (f64, f64, f64) {
    let (q, q_t, q_tt) = uv_q_of_temperature(eps_divk, temperature);
    let prefactor = 2.0 * PI / 3.0;
    (
        prefactor * q * q * q,
        3.0 * prefactor * q * q * q_t,
        6.0 * prefactor * q * q_t * q_t + 3.0 * prefactor * q * q * q_tt,
    )
}

/// The second virial coefficient of the Lennard-Jones spline and its
/// temperature derivatives (`ljs_b2`).
pub fn ljs_b2(eps_divk: f64, temperature: f64) -> (f64, f64, f64) {
    const BI: [f64; 5] = [1.345, -1.336, -3.85, 1.295, -0.416];
    let reduced = eps_divk / temperature;
    let reduced_t = -eps_divk / (temperature * temperature);
    let reduced_tt = 2.0 * eps_divk / (temperature * temperature * temperature);
    let reduced2 = reduced * reduced;
    let reduced3 = reduced2 * reduced;
    let b2 = BI[0]
        + BI[1] * reduced
        + BI[2] * reduced2
        + BI[3] * reduced3
        + BI[4] * (reduced2 * reduced2);
    let b2_t = BI[1] * reduced_t
        + 2.0 * BI[2] * reduced * reduced_t
        + 3.0 * BI[3] * reduced2 * reduced_t
        + 4.0 * BI[4] * reduced3 * reduced_t;
    let b2_tt = BI[1] * reduced_tt
        + 2.0 * BI[2] * reduced_t * reduced_t
        + 2.0 * BI[2] * reduced * reduced_tt
        + 6.0 * BI[3] * reduced * reduced_t * reduced_t
        + 3.0 * BI[3] * reduced2 * reduced_tt
        + 12.0 * BI[4] * reduced2 * reduced_t * reduced_t
        + 4.0 * BI[4] * reduced3 * reduced_tt;
    (b2, b2_t, b2_tt)
}

/// The difference of the second virial coefficient of the UV theory and its
/// derivatives (`uv_delta_b2_overall`).
pub fn uv_delta_b2_overall(
    sigma: f64,
    eps_divk: f64,
    temperature: f64,
    diameter: f64,
    diameter_t: f64,
    diameter_tt: f64,
    rho: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    let (a1_b2, a1_b2_t, a1_b2_tt) = uv_a1_b2(
        sigma,
        eps_divk,
        temperature,
        diameter,
        diameter_t,
        diameter_tt,
    );
    let (ljs, ljs_t, ljs_tt) = ljs_b2(eps_divk, temperature);
    let (hs, hs_t, hs_tt) = uv_hard_sphere_b2(eps_divk, temperature);
    let value = ljs - hs - a1_b2;
    let value_t = ljs_t - hs_t - a1_b2_t;
    let value_tt = ljs_tt - hs_tt - a1_b2_tt;
    (
        value * rho,
        value,
        value_t * rho,
        0.0,
        value_tt * rho,
        value_t,
    )
}

/// The first order term of the UV theory and its derivatives (`calc_uv_WCA`).
#[allow(clippy::too_many_arguments)]
pub fn calc_uv_wca(
    sigma: f64,
    eps_divk: f64,
    rho: f64,
    temperature: f64,
    diameter: f64,
    diameter_t: f64,
    diameter_tt: f64,
    enable_virial_term: bool,
    use_temperature_dependent_u_fraction: bool,
) -> (f64, f64, f64, f64, f64, f64) {
    let (b2, b2_r, b2_t, b2_rr, b2_tt, b2_rt) = if enable_virial_term {
        uv_delta_b2_overall(
            sigma,
            eps_divk,
            temperature,
            diameter,
            diameter_t,
            diameter_tt,
            rho,
        )
    } else {
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    };
    let (u, u_r, u_t, u_rr, u_tt, u_rt) = if enable_virial_term {
        if use_temperature_dependent_u_fraction {
            uv_phi_temp(eps_divk, temperature, rho)
        } else {
            uv_phi(eps_divk, temperature, rho)
        }
    } else {
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    };
    let (a1, a1_r, a1_t, a1_rr, a1_tt, a1_rt) = uv_a1_u(
        sigma,
        eps_divk,
        temperature,
        diameter,
        diameter_t,
        diameter_tt,
        rho,
    );
    let value = a1 + (1.0 - u) * b2;
    let value_r = a1_r + (1.0 - u) * b2_r - u_r * b2;
    let value_rr = a1_rr + (1.0 - u) * b2_rr - 2.0 * u_r * b2_r - u_rr * b2;
    let value_t = a1_t + (1.0 - u) * b2_t - u_t * b2;
    let value_tt = a1_tt + (1.0 - u) * b2_tt - 2.0 * u_t * b2_t - u_tt * b2;
    let value_rt = a1_rt + (1.0 - u) * b2_rt - u_rt * b2 - u_t * b2_r - u_r * b2_t;
    (value, value_r, value_t, value_rr, value_tt, value_rt)
}

/// Extensive UV-theory attractive contribution `F1 = n * a1_uv` and all 10 `(T, V, n)`
/// derivatives (`src/lj_splined.f90:calc_uv_WCA_TVN`).
#[allow(clippy::too_many_arguments)]
pub fn calc_uv_wca_tvn(
    sigma: f64,
    eps_divk: f64,
    temperature: f64,
    volume: f64,
    moles: f64,
    diameter: f64,
    diameter_t: f64,
    diameter_tt: f64,
    enable_virial_term: bool,
    use_temperature_dependent_u_fraction: bool,
) -> Result<PureHelmholtz> {
    positive(sigma, "segment diameter")?;
    positive(eps_divk, "well depth")?;
    positive(temperature, "temperature")?;
    positive(volume, "total volume")?;
    positive(moles, "mole number")?;
    let scale = crate::N_A * sigma * sigma * sigma;
    let rho_star = scale * moles / volume;
    let (a, a_r, a_t, a_rr, a_tt, a_rt) = calc_uv_wca(
        sigma,
        eps_divk,
        rho_star,
        temperature,
        diameter,
        diameter_t,
        diameter_tt,
        enable_virial_term,
        use_temperature_dependent_u_fraction,
    );
    let dh = DensityHelmholtz {
        value: a,
        dt: a_t,
        drho: a_r * scale,
        dtt: a_tt,
        dtrho: a_rt * scale,
        drhorho: a_rr * scale * scale,
    };
    dh.extensive(volume, moles)
}

/// One of the perturbation terms of the Barker-Henderson theory (`calc_ais_LJs_ex`).
pub fn ais_ljs_ex(index: usize, reduced_distance: f64, reduced_density: f64) -> Result<f64> {
    positive(reduced_distance, "the reduced distance")?;
    positive(reduced_density, "the reduced density")?;
    if !(1..=4).contains(&index) {
        return Err(Error::InvalidInput("the term index must be one to four"));
    }
    let prefactor = if index <= 2 {
        1.0 / reduced_density
    } else {
        1.0
    };
    let tau = 1.0 - 1.0 / reduced_distance;
    let (value, _, _, _, _, _) =
        polynomial_with_derivatives(index - 1, tau, reduced_density, prefactor);
    Ok(value)
}

/// The correction polynomial of the second order term (`calc_chi_LJs_ex`).
pub fn chi_ljs_ex(reduced_distance: f64, reduced_density: f64) -> Result<f64> {
    positive(reduced_distance, "the reduced distance")?;
    positive(reduced_density, "the reduced density")?;
    let tau = 1.0 - 1.0 / reduced_distance;
    let (value, _, _, _, _, _) = polynomial_with_derivatives(3, tau, reduced_density, 1.0);
    Ok(value)
}

// ============================================================================
// UF-theory kernels (`calc_uf_WCA` and `calc_uf_WCA_TVN`)
// ============================================================================

/// U-expansion parameter matrix `(5 x 4)` for `LJS-UF` (row-major by `i = 0..5, col = 0..4`).
const UF_U_LJS: [[f64; 4]; 5] = [
    [
        -4.19307460e-01,
        2.00550080e+00,
        -2.13577726e+00,
        4.88149448e-01,
    ],
    [
        2.63021608e+01,
        -1.04365573e-01,
        -2.66579307e+01,
        -6.41433360e-03,
    ],
    [
        6.14703510e+01,
        5.73140940e-01,
        -6.23488657e+01,
        7.42870835e-03,
    ],
    [-2.82732512, -1.61216849, 4.93092625, 2.23270821e-01],
    [
        6.09001685e-01,
        5.06954476e+01,
        -5.14171759e+01,
        4.95777156e-01,
    ],
];

/// U-expansion parameter matrix `(5 x 4)` for `LJ-UF`.
const UF_U_LJ: [[f64; 4]; 5] = [
    [-0.37447216, -0.00081593, -0.53140336, 0.06813345],
    [-0.01631296, 0.0173549, -0.48321813, -0.10811097],
    [-0.11968918, 0.06900995, -0.09297811, -0.05479239],
    [0.42113581, -0.11298284, 0.06698083, 0.13611539],
    [-0.13660913, 0.02263133, 0.00933081, 0.21970435],
];

/// Mayer-f-expansion parameter matrix `(5 x 5)` for `LJS-UF` (`param[i][j]` matching Fortran `param(i, j)`).
const UF_MF_LJS: [[f64; 5]; 5] = [
    [
        -1.38956236e-03,
        4.39757272e-02,
        -1.77080263e-01,
        2.46524435e-01,
        -1.14440939e-01,
    ],
    [
        7.21649732e-01,
        -6.53666725e-01,
        3.60375455e+00,
        -5.16905418e+00,
        2.39307542e+00,
    ],
    [
        -4.46006380e-01,
        6.21534008e+00,
        -2.00895809e+01,
        2.80637697e+01,
        -1.34547589e+01,
    ],
    [
        1.02494323e+00,
        -1.04680007e+01,
        3.64523410e+01,
        -5.00655126e+01,
        2.34538041e+01,
    ],
    [
        -3.98081890e-01,
        4.57504276e+00,
        -1.57757874e+01,
        2.15679397e+01,
        -1.00079543e+01,
    ],
];

/// Mayer-f-expansion parameter matrix `(5 x 5)` for `LJ-UF`.
const UF_MF_LJ: [[f64; 5]; 5] = [
    [0.00326607, -0.00631528, -0.00134063, 0.00803155, -0.0041289],
    [0.98177233, 0.37428833, -0.0261357, -0.23353006, 0.11235306],
    [0.09812966, 0.4494362, 0.18812166, -0.0407614, -0.191682],
    [0.1091711, -0.10567363, -0.01328055, 0.04651761, 0.01508],
    [0.00246771, 0.06797539, -0.0037673, 0.08330469, -0.06683759],
];

fn uf_u_a1(param: &[[f64; 4]; 5], rho: f64, t_star: f64) -> (f64, f64, f64, f64, f64, f64) {
    let sqrt_t = t_star.sqrt();
    let mut ai = [0.0; 5];
    let mut ai_t = [0.0; 5];
    let mut ai_tt = [0.0; 5];
    for i in 0..5 {
        let p0 = param[i][0];
        let p1 = param[i][1];
        let p2 = param[i][2];
        let p3 = param[i][3];
        ai[i] = p0 + p1 * sqrt_t + p2 * t_star.powf(p3);
        ai_t[i] = 0.5 * p1 / sqrt_t + p2 * p3 * t_star.powf(p3 - 1.0);
        ai_tt[i] = -0.25 * p1 / (sqrt_t * t_star) + p2 * p3 * (p3 - 1.0) * t_star.powf(p3 - 2.0);
    }
    let mut rhoa = [0.0; 5];
    let mut rhoa_r = [0.0; 5];
    let mut rhoa_rr = [0.0; 5];
    rhoa[0] = rho;
    rhoa_r[0] = 1.0;
    rhoa_rr[0] = 0.0;
    for i in 1..5 {
        let order = (i + 1) as f64;
        rhoa[i] = rhoa[i - 1] * rho;
        rhoa_r[i] = order * rhoa[i - 1];
        rhoa_rr[i] = order * rhoa_r[i - 1];
    }
    let (mut a, mut a_r, mut a_t, mut a_rr, mut a_tt, mut a_rt) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for i in 0..5 {
        a += ai[i] * rhoa[i];
        a_r += ai[i] * rhoa_r[i];
        a_t += ai_t[i] * rhoa[i];
        a_rr += ai[i] * rhoa_rr[i];
        a_tt += ai_tt[i] * rhoa[i];
        a_rt += ai_t[i] * rhoa_r[i];
    }
    let two_pi = 2.0 * PI;
    (
        two_pi * a,
        two_pi * a_r,
        two_pi * a_t,
        two_pi * a_rr,
        two_pi * a_tt,
        two_pi * a_rt,
    )
}

fn uf_mf_a1(param: &[[f64; 5]; 5], rho: f64, t_star: f64) -> (f64, f64, f64, f64, f64, f64) {
    let mut rhoa = [0.0; 5];
    let mut rhoa_r = [0.0; 5];
    let mut rhoa_rr = [0.0; 5];
    rhoa[0] = rho;
    rhoa_r[0] = 1.0;
    rhoa_rr[0] = 0.0;
    for i in 1..5 {
        let order = (i + 1) as f64;
        rhoa[i] = rhoa[i - 1] * rho;
        rhoa_r[i] = order * rhoa[i - 1];
        rhoa_rr[i] = order * rhoa_r[i - 1];
    }
    let b = 1.0 / t_star;
    let mut beta = [0.0; 5];
    let mut beta_t = [0.0; 5];
    let mut beta_tt = [0.0; 5];
    beta[0] = t_star;
    beta[1] = 1.0;
    beta_t[0] = 1.0;
    for i in 2..5 {
        let k = (i - 1) as f64;
        beta[i] = beta[i - 1] * b;
        beta_t[i] = -k * beta[i] * b;
        beta_tt[i] = -(k + 1.0) * beta_t[i] * b;
    }
    let (mut a, mut a_r, mut a_t, mut a_rr, mut a_tt, mut a_rt) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for i in 0..5 {
        for j in 0..5 {
            let pij = param[i][j];
            a += pij * beta[i] * rhoa[j];
            a_r += pij * beta[i] * rhoa_r[j];
            a_rr += pij * beta[i] * rhoa_rr[j];
            a_rt += pij * beta_t[i] * rhoa_r[j];
            a_t += pij * beta_t[i] * rhoa[j];
            a_tt += pij * beta_tt[i] * rhoa[j];
        }
    }
    let neg_two_pi = -2.0 * PI;
    (
        neg_two_pi * a,
        neg_two_pi * a_r,
        neg_two_pi * a_t,
        neg_two_pi * a_rr,
        neg_two_pi * a_tt,
        neg_two_pi * a_rt,
    )
}

fn uf_phi_2param(param: [f64; 2], rho: f64) -> (f64, f64, f64) {
    let a = param[0];
    let b = param[1];
    let x = a * rho + b * rho * rho;
    let dxdr = a + 2.0 * b * rho;
    let d2xdr2 = 2.0 * b;
    safe_tanh_chain(x, dxdr, d2xdr2)
}

/// First-order UF-theory WCA perturbation term `a1 / eps` and its `(rho*, T*)`
/// derivatives (`src/lj_splined.f90:calc_uf_WCA`).
///
/// `is_lj = false` selects the splined Lennard-Jones potential (`LJS-UF`);
/// `is_lj = true` selects the full Lennard-Jones potential (`LJ-UF`).
pub fn calc_uf_wca(rho: f64, t_star: f64, is_lj: bool) -> (f64, f64, f64, f64, f64, f64) {
    let ((phi_u, phi_u_r, phi_u_rr), amf, au) = if is_lj {
        (
            uf_phi_2param([2.215, 2.033], rho),
            uf_mf_a1(&UF_MF_LJ, rho, t_star),
            uf_u_a1(&UF_U_LJ, rho, t_star),
        )
    } else {
        (
            uf_phi_2param([2.47, 1.4], rho),
            uf_mf_a1(&UF_MF_LJS, rho, t_star),
            uf_u_a1(&UF_U_LJS, rho, t_star),
        )
    };
    let one_m = 1.0 - phi_u;
    let a = phi_u * au.0 + one_m * amf.0;
    let a_r = phi_u * au.1 + one_m * amf.1 + phi_u_r * (au.0 - amf.0);
    let a_rr = phi_u_r * au.1 + phi_u * au.3 - phi_u_r * amf.1
        + one_m * amf.3
        + phi_u_r * (au.1 - amf.1)
        + phi_u_rr * (au.0 - amf.0);
    let a_t = phi_u * au.2 + one_m * amf.2;
    let a_tt = phi_u * au.4 + one_m * amf.4;
    let a_rt = phi_u * au.5 + one_m * amf.5 + phi_u_r * (au.2 - amf.2);
    (a, a_r, a_t, a_rr, a_tt, a_rt)
}

/// Extensive UF-theory attractive term `n * eps * a1_uf` (before dividing by `T`)
/// and all 10 `(T, V, n)` derivatives (`src/lj_splined.f90:calc_uf_WCA_TVN`).
pub fn calc_uf_wca_tvn(
    sigma: f64,
    eps_divk: f64,
    is_lj: bool,
    temperature: f64,
    volume: f64,
    moles: f64,
) -> Result<PureHelmholtz> {
    positive(sigma, "segment diameter")?;
    positive(eps_divk, "well depth")?;
    positive(temperature, "temperature")?;
    positive(volume, "total volume")?;
    positive(moles, "mole number")?;
    let scale = crate::N_A * sigma * sigma * sigma;
    let rho_star = scale * moles / volume;
    let t_star = temperature / eps_divk;
    let (a, a_r, a_t, a_rr, a_tt, a_rt) = calc_uf_wca(rho_star, t_star, is_lj);
    let dh = DensityHelmholtz {
        value: eps_divk * a,
        dt: a_t,
        drho: eps_divk * a_r * scale,
        dtt: a_tt / eps_divk,
        dtrho: a_rt * scale,
        drhorho: eps_divk * a_rr * scale * scale,
    };
    dh.extensive(volume, moles)
}

// ============================================================================
// Full UV-theory model (`LjsUv` / `PreparedLjsUv`, `calcFresLJs_uv_theory`)
// ============================================================================

/// Model controls for the splined Lennard-Jones UV-theory equation of state (`LJS-UV`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LjsUvControl {
    /// Include the WCA hard-sphere reference (and cavity term if `enable_cavity` is true).
    pub enable_hard_sphere: bool,
    /// Include the WCA cavity soft-repulsion integral (requires `enable_hard_sphere`).
    pub enable_cavity: bool,
    /// Include the UV-theory attractive contribution.
    pub enable_a1: bool,
    /// Include the low-density second-virial correction term `(1 - u) * Delta B2 * rho*`.
    pub enable_virial_term: bool,
    /// Use the temperature-dependent u-fraction (`uv_phi_temp`) instead of `uv_phi`.
    pub use_temperature_dependent_u_fraction: bool,
    /// Use the Carnahan-Starling-Kolafa pure hard-sphere reference (`PURE_HS_CSK`).
    pub kolafa_pure_reference: bool,
}

impl Default for LjsUvControl {
    fn default() -> Self {
        Self {
            enable_hard_sphere: true,
            enable_cavity: true,
            enable_a1: true,
            enable_virial_term: true,
            use_temperature_dependent_u_fraction: false,
            kolafa_pure_reference: false,
        }
    }
}

/// The splined Lennard-Jones UV-theory equation of state (`LJS-UV`).
///
/// Ported from `src/lj_splined.f90:calcFresLJs_uv_theory`.
#[derive(Debug, Clone, Copy)]
pub struct LjsUv {
    sigma: f64,
    eps: f64,
    control: LjsUvControl,
}

impl Default for LjsUv {
    fn default() -> Self {
        Self::new()
    }
}

/// Temperature-prepared `LJS-UV` state caching WCA hard-sphere diameter,
/// temperature-only second-virial difference, and cavity quadrature nodes.
#[derive(Debug, Clone)]
pub struct PreparedLjsUv {
    sigma: f64,
    eps: f64,
    control: LjsUvControl,
    temperature: f64,
    diameter: f64,
    diameter_t: f64,
    diameter_tt: f64,
    eta_per_rho: f64,
    eta_per_rho_t: f64,
    eta_per_rho_tt: f64,
    rho_star_per_rho: f64,
    delta_b2: f64,
    delta_b2_t: f64,
    delta_b2_tt: f64,
    cavity: Option<CavityPreparedTemp>,
}

impl LjsUv {
    /// Create the default `LJS-UV` model (Argon record).
    pub fn new() -> Self {
        Self {
            sigma: ARGON_SIGMA,
            eps: ARGON_EPS,
            control: LjsUvControl::default(),
        }
    }

    /// Create an `LJS-UV` model with custom `(sigma, eps/k)` parameters.
    pub fn with_parameters(sigma: f64, eps: f64) -> Self {
        Self {
            sigma,
            eps,
            control: LjsUvControl::default(),
        }
    }

    /// Replace `(sigma, eps/k)` on the model.
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
    pub fn control(&self) -> LjsUvControl {
        self.control
    }

    /// Replace model control flags (`ljs_uv_model_control` / `ljs_wca_model_control`).
    pub fn set_control(&mut self, control: LjsUvControl) {
        self.control = control;
    }

    /// WCA hard-sphere diameter and its first and second temperature derivatives.
    pub fn diameter(&self, temperature: f64) -> Result<(f64, f64, f64)> {
        positive(self.sigma, "segment diameter")?;
        positive(self.eps, "well depth")?;
        positive(temperature, "temperature")?;
        Ok(wca_hard_sphere_diameter(self.sigma, self.eps, temperature))
    }

    /// Hypothetical pure-fluid packing fraction `eta` at `(T, rho)`.
    pub fn packing_fraction(&self, temperature: f64, density: f64) -> Result<f64> {
        positive(density, "molar density")?;
        let (d, _, _) = self.diameter(temperature)?;
        Ok(PI / 6.0 * crate::N_A * d * d * d * density)
    }

    /// Packing fraction and all `(T, V, n)` derivatives.
    pub fn packing_fraction_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<PackingFractionDerivatives> {
        let (d, d_t, d_tt) = self.diameter(temperature)?;
        packing_fraction_vdw_no_segments(volume, &[moles], &[d], &[d_t], &[d_tt])
    }

    /// Prepare temperature-only quantities for repeated evaluations at `temperature`.
    pub fn prepare(&self, temperature: f64) -> Result<PreparedLjsUv> {
        positive(self.sigma, "segment diameter")?;
        positive(self.eps, "well depth")?;
        positive(temperature, "temperature")?;
        let (d, d_t, d_tt) = wca_hard_sphere_diameter(self.sigma, self.eps, temperature);
        let (eta_per_rho, eta_per_rho_t, eta_per_rho_tt) = pure_eta_per_rho(d, d_t, d_tt);
        let rho_star_per_rho = crate::N_A * self.sigma * self.sigma * self.sigma;
        let (delta_b2, delta_b2_t, delta_b2_tt) = if self.control.enable_virial_term {
            let (a1_b2, a1_b2_t, a1_b2_tt) =
                uv_a1_b2(self.sigma, self.eps, temperature, d, d_t, d_tt);
            let (ljs, ljs_t, ljs_tt) = ljs_b2(self.eps, temperature);
            let (hs, hs_t, hs_tt) = uv_hard_sphere_b2(self.eps, temperature);
            (
                ljs - hs - a1_b2,
                ljs_t - hs_t - a1_b2_t,
                ljs_tt - hs_tt - a1_b2_tt,
            )
        } else {
            (0.0, 0.0, 0.0)
        };
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
        Ok(PreparedLjsUv {
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
            delta_b2,
            delta_b2_t,
            delta_b2_tt,
            cavity,
        })
    }

    /// Molar reduced residual Helmholtz energy `a_res(T, rho)` and all 6 `(T, rho)` derivatives.
    pub fn residual_density(&self, temperature: f64, rho: f64) -> Result<DensityHelmholtz> {
        self.prepare(temperature)?.residual_density(rho)
    }

    /// Extensive reduced residual Helmholtz energy `F = n * a_res` and all 10 `(T, V, n)`
    /// derivatives (`calcFresLJs_uv_theory`).
    pub fn residual_helmholtz_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<PureHelmholtz> {
        self.prepare(temperature)?
            .residual_helmholtz_tvn(volume, moles)
    }

    /// Molar reduced hard-sphere contribution and all 10 `(T, V, n)` derivatives.
    pub fn hard_sphere_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<PureHelmholtz> {
        self.prepare(temperature)?.hard_sphere_tvn(volume, moles)
    }

    /// Molar reduced WCA soft-repulsion contribution and all 10 `(T, V, n)` derivatives.
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

    /// Molar reduced dispersion contribution and all 10 `(T, V, n)` derivatives.
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

    /// Full `(T, V, n)` thermodynamic state.
    pub fn tv_state(&self, temperature: f64, volume: f64, moles: f64) -> Result<LjsTvState> {
        self.prepare(temperature)?.tv_state(volume, moles)
    }

    /// Legacy summary state at `(T, V, n)`.
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

    /// Exact source SAFT non-association volume iteration.
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

    /// Estimated critical temperature (K) and critical molar volume (m^3/mol).
    pub fn estimate_critical_parameters(&self) -> Result<(f64, f64)> {
        estimate_ljs_critical_parameters(self.sigma, self.eps)
    }

    /// Fundamental Measure Theory compatibility flags `(is_fmt_consistent, na_enabled)`.
    pub fn fmt_compatibility(&self) -> (bool, bool) {
        (true, false)
    }

    /// Dimensionless van der Waals attractive energy constant (`alpha_ljs`).
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

impl PreparedLjsUv {
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

    fn dispersion_density(&self, rho: f64) -> Result<DensityHelmholtz> {
        positive(rho, "molar density")?;
        if !self.control.enable_a1 {
            return Ok(DensityHelmholtz::zero());
        }
        let scale = self.rho_star_per_rho;
        let rho_star = scale * rho;
        let (a1, a1_r, a1_t, a1_rr, a1_tt, a1_rt) = uv_a1_u(
            self.sigma,
            self.eps,
            self.temperature,
            self.diameter,
            self.diameter_t,
            self.diameter_tt,
            rho_star,
        );
        let (val, val_r, val_t, val_rr, val_tt, val_rt) = if self.control.enable_virial_term {
            let b2 = self.delta_b2 * rho_star;
            let b2_r = self.delta_b2;
            let b2_t = self.delta_b2_t * rho_star;
            let b2_tt = self.delta_b2_tt * rho_star;
            let b2_rt = self.delta_b2_t;
            let (u, u_r, u_t, u_rr, u_tt, u_rt) =
                if self.control.use_temperature_dependent_u_fraction {
                    uv_phi_temp(self.eps, self.temperature, rho_star)
                } else {
                    uv_phi(self.eps, self.temperature, rho_star)
                };
            (
                a1 + (1.0 - u) * b2,
                a1_r + (1.0 - u) * b2_r - u_r * b2,
                a1_t + (1.0 - u) * b2_t - u_t * b2,
                a1_rr - 2.0 * u_r * b2_r - u_rr * b2,
                a1_tt + (1.0 - u) * b2_tt - 2.0 * u_t * b2_t - u_tt * b2,
                a1_rt + (1.0 - u) * b2_rt - u_rt * b2 - u_t * b2_r - u_r * b2_t,
            )
        } else {
            (a1, a1_r, a1_t, a1_rr, a1_tt, a1_rt)
        };
        Ok(DensityHelmholtz {
            value: val,
            dt: val_t,
            drho: val_r * scale,
            dtt: val_tt,
            dtrho: val_rt * scale,
            drhorho: val_rr * scale * scale,
        })
    }

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
        if self.control.enable_a1 {
            let scale = self.rho_star_per_rho;
            let rho_star = scale * rho;
            let d_star = self.diameter / self.sigma;
            let (a1_t, a1_rt, _, a1_rrt, _, _, _) = uv_a1_u_mult_t(self.eps, d_star, rho_star);
            let inv_t = 1.0 / self.temperature;
            let mut val = a1_t * inv_t;
            let mut val_r = a1_rt * inv_t;
            let mut val_rr = a1_rrt * inv_t;
            if self.control.enable_virial_term {
                let b2 = self.delta_b2 * rho_star;
                let b2_r = self.delta_b2;
                let (u, u_r, _, u_rr, _, _) = if self.control.use_temperature_dependent_u_fraction {
                    uv_phi_temp(self.eps, self.temperature, rho_star)
                } else {
                    uv_phi(self.eps, self.temperature, rho_star)
                };
                val += (1.0 - u) * b2;
                val_r += (1.0 - u) * b2_r - u_r * b2;
                val_rr += -2.0 * u_r * b2_r - u_rr * b2;
            }
            a += val;
            a_r += val_r * scale;
            a_rr += val_rr * scale * scale;
        }
        if !a.is_finite() || !a_r.is_finite() || !a_rr.is_finite() {
            return Err(Error::InvalidInput("non-finite LJS-UV density evaluation"));
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
        total = total.add(self.dispersion_density(rho)?);
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

    /// Full `(T, V, n)` thermodynamic state.
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

    /// Exact source SAFT non-association volume iteration.
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

// ============================================================================
// Full UF-theory model (`LjsUf` / `PreparedLjsUf`, `calcFresLJ_uf_theory`)
// ============================================================================

/// Underlying pair potential for UF-theory (`LJS-UF` vs `LJ-UF`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LjUfPotential {
    /// Splined Lennard-Jones potential (`LJS-UF`, `lj_potential = false`).
    LennardJonesSplined,
    /// Full 12-6 Lennard-Jones potential (`LJ-UF`, `lj_potential = true`).
    LennardJones,
}

/// Model controls for the UF-theory equation of state (`LJS-UF` / `LJ-UF`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LjsUfControl {
    /// Include the WCA hard-sphere reference (and cavity term if `enable_cavity` is true).
    pub enable_hard_sphere: bool,
    /// Include the WCA cavity soft-repulsion integral (requires `enable_hard_sphere`).
    pub enable_cavity: bool,
    /// Include the UF-theory attractive contribution.
    pub enable_a1: bool,
    /// Use the Carnahan-Starling-Kolafa pure hard-sphere reference (`PURE_HS_CSK`).
    pub kolafa_pure_reference: bool,
}

impl Default for LjsUfControl {
    fn default() -> Self {
        Self {
            enable_hard_sphere: true,
            enable_cavity: true,
            enable_a1: true,
            kolafa_pure_reference: false,
        }
    }
}

/// The UF-theory equation of state for splined (`LJS-UF`) or full (`LJ-UF`)
/// Lennard-Jones fluids.
///
/// Ported from `src/lj_splined.f90:calcFresLJ_uf_theory`.
/// Note: upstream `ljx_ux_eos_constructor` reverses the arguments to Fortran
/// `INDEX`, preventing `init_ljs("UF")` and `init_lj("UF")` from selecting
/// `calcFresLJ_uf_theory`; this native model implements `calcFresLJ_uf_theory`
/// directly and cleanly for both potentials.
#[derive(Debug, Clone, Copy)]
pub struct LjsUf {
    potential: LjUfPotential,
    sigma: f64,
    eps: f64,
    control: LjsUfControl,
}

/// Temperature-prepared UF-theory state.
#[derive(Debug, Clone)]
pub struct PreparedLjsUf {
    potential: LjUfPotential,
    sigma: f64,
    eps: f64,
    control: LjsUfControl,
    temperature: f64,
    diameter: f64,
    diameter_t: f64,
    diameter_tt: f64,
    eta_per_rho: f64,
    eta_per_rho_t: f64,
    eta_per_rho_tt: f64,
    rho_star_per_rho: f64,
    cavity: Option<CavityPreparedTemp>,
}

impl LjsUf {
    /// Create a default `LJS-UF` or `LJ-UF` model (Argon record).
    pub fn new(potential: LjUfPotential) -> Self {
        Self {
            potential,
            sigma: ARGON_SIGMA,
            eps: ARGON_EPS,
            control: LjsUfControl::default(),
        }
    }

    /// Create an `LJS-UF` or `LJ-UF` model with custom `(sigma, eps/k)` parameters.
    pub fn with_parameters(potential: LjUfPotential, sigma: f64, eps: f64) -> Self {
        Self {
            potential,
            sigma,
            eps,
            control: LjsUfControl::default(),
        }
    }

    /// Underlying pair potential (`LJS-UF` vs `LJ-UF`).
    pub fn potential(&self) -> LjUfPotential {
        self.potential
    }

    /// Replace `(sigma, eps/k)` on the model.
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
    pub fn control(&self) -> LjsUfControl {
        self.control
    }

    /// Replace model control flags.
    pub fn set_control(&mut self, control: LjsUfControl) {
        self.control = control;
    }

    /// WCA hard-sphere diameter and its first and second temperature derivatives.
    pub fn diameter(&self, temperature: f64) -> Result<(f64, f64, f64)> {
        positive(self.sigma, "segment diameter")?;
        positive(self.eps, "well depth")?;
        positive(temperature, "temperature")?;
        Ok(wca_hard_sphere_diameter(self.sigma, self.eps, temperature))
    }

    /// Hypothetical pure-fluid packing fraction `eta` at `(T, rho)`.
    pub fn packing_fraction(&self, temperature: f64, density: f64) -> Result<f64> {
        positive(density, "molar density")?;
        let (d, _, _) = self.diameter(temperature)?;
        Ok(PI / 6.0 * crate::N_A * d * d * d * density)
    }

    /// Packing fraction and all `(T, V, n)` derivatives.
    pub fn packing_fraction_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<PackingFractionDerivatives> {
        let (d, d_t, d_tt) = self.diameter(temperature)?;
        packing_fraction_vdw_no_segments(volume, &[moles], &[d], &[d_t], &[d_tt])
    }

    /// Prepare temperature-only quantities for repeated evaluations at `temperature`.
    pub fn prepare(&self, temperature: f64) -> Result<PreparedLjsUf> {
        positive(self.sigma, "segment diameter")?;
        positive(self.eps, "well depth")?;
        positive(temperature, "temperature")?;
        let (d, d_t, d_tt) = wca_hard_sphere_diameter(self.sigma, self.eps, temperature);
        let (eta_per_rho, eta_per_rho_t, eta_per_rho_tt) = pure_eta_per_rho(d, d_t, d_tt);
        let rho_star_per_rho = crate::N_A * self.sigma * self.sigma * self.sigma;
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
        Ok(PreparedLjsUf {
            potential: self.potential,
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
            cavity,
        })
    }

    /// Molar reduced residual Helmholtz energy `a_res(T, rho)` and all 6 `(T, rho)` derivatives.
    pub fn residual_density(&self, temperature: f64, rho: f64) -> Result<DensityHelmholtz> {
        self.prepare(temperature)?.residual_density(rho)
    }

    /// Extensive reduced residual Helmholtz energy `F = n * a_res` and all 10 `(T, V, n)`
    /// derivatives (`calcFresLJ_uf_theory`).
    pub fn residual_helmholtz_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<PureHelmholtz> {
        self.prepare(temperature)?
            .residual_helmholtz_tvn(volume, moles)
    }

    /// Molar reduced hard-sphere contribution and all 10 `(T, V, n)` derivatives.
    pub fn hard_sphere_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: f64,
    ) -> Result<PureHelmholtz> {
        self.prepare(temperature)?.hard_sphere_tvn(volume, moles)
    }

    /// Molar reduced WCA soft-repulsion contribution and all 10 `(T, V, n)` derivatives.
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

    /// Molar reduced dispersion contribution and all 10 `(T, V, n)` derivatives.
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

    /// Full `(T, V, n)` thermodynamic state.
    pub fn tv_state(&self, temperature: f64, volume: f64, moles: f64) -> Result<LjsTvState> {
        self.prepare(temperature)?.tv_state(volume, moles)
    }

    /// Legacy summary state at `(T, V, n)`.
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

    /// Exact source SAFT non-association volume iteration.
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

    /// Estimated critical temperature (K) and critical molar volume (m^3/mol).
    pub fn estimate_critical_parameters(&self) -> Result<(f64, f64)> {
        estimate_ljs_critical_parameters(self.sigma, self.eps)
    }

    /// Fundamental Measure Theory compatibility flags `(is_fmt_consistent, na_enabled)`.
    pub fn fmt_compatibility(&self) -> (bool, bool) {
        (true, false)
    }

    /// Dimensionless van der Waals attractive energy constant (`alpha_ljs`).
    pub fn alpha_vdw(&self) -> f64 {
        ALPHA_LJS
    }
}

impl PreparedLjsUf {
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

    fn is_lj(&self) -> bool {
        matches!(self.potential, LjUfPotential::LennardJones)
    }

    fn dispersion_density(&self, rho: f64) -> Result<DensityHelmholtz> {
        positive(rho, "molar density")?;
        if !self.control.enable_a1 {
            return Ok(DensityHelmholtz::zero());
        }
        let scale = self.rho_star_per_rho;
        let rho_star = scale * rho;
        let t_star = self.temperature / self.eps;
        let (a, a_r, a_t, a_rr, a_tt, a_rt) = calc_uf_wca(rho_star, t_star, self.is_lj());
        let f1 = DensityHelmholtz {
            value: self.eps * a,
            dt: a_t,
            drho: self.eps * a_r * scale,
            dtt: a_tt / self.eps,
            dtrho: a_rt * scale,
            drhorho: self.eps * a_rr * scale * scale,
        };
        Ok(f1.scale_inv_t_power(self.temperature, 1))
    }

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
        if self.control.enable_a1 {
            let scale = self.rho_star_per_rho;
            let rho_star = scale * rho;
            let t_star = self.temperature / self.eps;
            let (v, vr, _, vrr, _, _) = calc_uf_wca(rho_star, t_star, self.is_lj());
            let factor = self.eps / self.temperature;
            a += factor * v;
            a_r += factor * vr * scale;
            a_rr += factor * vrr * scale * scale;
        }
        if !a.is_finite() || !a_r.is_finite() || !a_rr.is_finite() {
            return Err(Error::InvalidInput("non-finite UF density evaluation"));
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
        total = total.add(self.dispersion_density(rho)?);
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

    /// Full `(T, V, n)` thermodynamic state.
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

    /// Exact source SAFT non-association volume iteration.
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
