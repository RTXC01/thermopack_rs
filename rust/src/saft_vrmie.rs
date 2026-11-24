// SAFT-VR Mie and SAFT-VRQ Mie Equation of State.
// Reference: Lafitte et al. (2013) doi:10.1063/1.4819786 and Aasen et al. (2019) doi:10.1063/1.5136079.
// Copyright thermotools and Contributors. Apache-2.0.

mod association;
mod kernel;
use self::association::MieRdf;
use self::kernel::EvalMode;
use crate::association::{
    AssocScheme, AssociationParams, DeltaCombRule, assoc_helmholtz_all, assoc_helmholtz_first,
    assoc_helmholtz_pressure_jet, assoc_helmholtz_value, solve_for_x_k,
};
/// Extensive reduced Helmholtz energy (mol), all TVN derivatives and site fractions.
pub type SaftVrMieHelmholtz = crate::association::AssociationHelmholtz;
use crate::component::component;
use crate::cubic::Phase;
use crate::error::{Error, Result, positive};
use crate::ideal::IdealGas;
use crate::saft_vrmie_data;

// Universal physical constants matching ThermoPack (CODATA 2018 / 2019 SI)
const N_A: f64 = 6.02214076e23;
const K_B: f64 = 1.380649e-23;
const RGAS: f64 = N_A * K_B; // 8.31446261815324 J/(mol*K)
const H_PLANCK: f64 = 6.626069311e-34; // J*s, thermopack_constants::h_const
const PI: f64 = std::f64::consts::PI;
const MACHINE_PREC: f64 = 2.220446049250313e-16;
const DEFAULT_MAX_REDUCED_POTENTIAL: f64 = 500.0;

// Gauss-Kronrod 21 weights and positive nodes for Barker-Henderson diameter
#[allow(clippy::excessive_precision)]
const XK_21: [f64; 11] = [
    0.995657163025808080735527280689003,
    0.973906528517171720077964012084452,
    0.930157491355708226001207180059508,
    0.865063366688984510732096688423493,
    0.780817726586416897063717578345042,
    0.679409568299024406234327365114874,
    0.562757134668604683339000099272694,
    0.433395394129247190799265943165784,
    0.294392862701460198131126603103866,
    0.148874338981631210884826001129720,
    0.0,
];

#[allow(clippy::excessive_precision)]
const WK_21: [f64; 11] = [
    0.011694638867371874278064396062192,
    0.032558162307964727478818972459390,
    0.054755896574351996031381300244580,
    0.075039674810919952767043140916190,
    0.093125454583697605535065465083366,
    0.109387158802297641899210590325805,
    0.123491976262065851077958109831074,
    0.134709217311473325928054001771707,
    0.142775938577060080797094273138717,
    0.147739104901338491374841515972068,
    0.149445554002916905664936468389821,
];

#[allow(clippy::excessive_precision)]
const X_NODES_21: [f64; 21] = [
    -XK_21[0], -XK_21[1], -XK_21[2], -XK_21[3], -XK_21[4], -XK_21[5], -XK_21[6], -XK_21[7],
    -XK_21[8], -XK_21[9], XK_21[10], XK_21[9], XK_21[8], XK_21[7], XK_21[6], XK_21[5], XK_21[4],
    XK_21[3], XK_21[2], XK_21[1], XK_21[0],
];

#[allow(clippy::excessive_precision)]
const W_NODES_21: [f64; 21] = [
    WK_21[0], WK_21[1], WK_21[2], WK_21[3], WK_21[4], WK_21[5], WK_21[6], WK_21[7], WK_21[8],
    WK_21[9], WK_21[10], WK_21[9], WK_21[8], WK_21[7], WK_21[6], WK_21[5], WK_21[4], WK_21[3],
    WK_21[2], WK_21[1], WK_21[0],
];

/// Evaluate `exp(-z) - 1 + z` without catastrophic cancellation for small `|z|`.
#[inline]
fn expm1_minus(z: f64) -> f64 {
    if z.abs() < 1.0e-2 {
        let z2 = z * z;
        z2 * (0.5
            - z * (1.0 / 6.0
                - z * (1.0 / 24.0
                    - z * (1.0 / 120.0 - z * (1.0 / 720.0 - z * (1.0 / 5040.0 - z / 40320.0))))))
    } else {
        (-z).exp_m1() + z
    }
}

/// Evaluate the reduced Mie excess above the well minimum, `u*(r_min * (1 + tau)) + 1 >= 0`,
/// without catastrophic cancellation near `tau = 0`.
#[inline]
fn mie_well_excess_tau(tau: f64, la: f64, lr: f64) -> f64 {
    let ln_t = tau.ln_1p();
    (la * expm1_minus(lr * ln_t) - lr * expm1_minus(la * ln_t)) / (lr - la)
}

/// Evaluate `d[u*(r_min - d)]/dd > 0` for `d > 0` without cancellation near `d = 0`.
#[inline]
fn mie_well_left_deriv(d: f64, r_min: f64, la: f64, lr: f64) -> f64 {
    let tau = -d / r_min;
    let ln_t = tau.ln_1p();
    let r = r_min - d;
    let pref = (la * lr) / ((lr - la) * r);
    pref * (-la * ln_t).exp() * (-(lr - la) * ln_t).exp_m1()
}

/// Evaluate `u*(1 + delta)` for `delta >= 0` without cancellation near `r = 1`.
#[inline]
fn mie_shifted_unit_value(delta: f64, c: f64, la: f64, lr: f64) -> f64 {
    let s = delta.ln_1p();
    c * (-la * s).exp() * (-(lr - la) * s).exp_m1()
}

/// Evaluate `du*(1 + delta)/d(delta)` for `delta >= 0`.
#[inline]
fn mie_shifted_unit_slope(delta: f64, c: f64, la: f64, lr: f64) -> f64 {
    let s = delta.ln_1p();
    (c / (1.0 + delta)) * (-lr * (-lr * s).exp() + la * (-la * s).exp())
}

/// Fast positive-exponent power avoiding `libm::pow` for exact small integers
/// and sharing a single `ln(base)` across `(base^n, base^m)` when both are non-integer.
#[inline(always)]
fn pow_pos_pair(base: f64, n: f64, m: f64) -> (f64, f64) {
    let n_int = (1.0..=64.0).contains(&n) && n == (n as i32) as f64;
    let m_int = (1.0..=64.0).contains(&m) && m == (m as i32) as f64;
    let eval_int = |b: f64, e: f64| -> f64 {
        if e == 6.0 {
            let b2 = b * b;
            b2 * b2 * b2
        } else {
            b.powi(e as i32)
        }
    };
    match (n_int, m_int) {
        (true, true) => (eval_int(base, n), eval_int(base, m)),
        (true, false) => (eval_int(base, n), (m * base.ln()).exp()),
        (false, true) => ((n * base.ln()).exp(), eval_int(base, m)),
        (false, false) => {
            let ln_b = base.ln();
            ((n * ln_b).exp(), (m * ln_b).exp())
        }
    }
}

// Effective packing fraction polynomial coefficients for Sutherland potential (row: 0..3, col: 0..3)
const LAM_COEFF: [[f64; 4]; 4] = [
    [0.81096, 1.0205, -1.9057, 1.0885],
    [1.7888, -19.341, 22.845, -6.1962],
    [-37.578, 151.26, -228.14, 106.98],
    [92.284, -463.50, 973.92, -677.64],
];

// Phi matrix for Lafitte alpha functions f_1 .. f_6 (row: 0..6, col: 0..5)
const PHI_MAT: [[f64; 6]; 7] = [
    [7.5365557, -359.44, 1550.9, -1.19932, -1911.28, 9236.9],
    [-37.60463, 1825.6, -5070.1, 9.063632, 21390.175, -129430.0],
    [71.745953, -3168.0, 6534.6, -17.9482, -51320.7, 357230.0],
    [-46.83552, 1884.2, -3288.7, 11.34027, 37064.54, -315530.0],
    [-2.467982, -0.82376, -2.7171, 20.52142, 1103.742, 1390.2],
    [-0.50272, -3.1935, 2.0883, -56.6377, -3264.61, -4518.2],
    [8.0956883, 3.709, 0.0, 40.53683, 2556.181, 4241.6],
];

/// Component parameters for SAFT-VR Mie and SAFT-VRQ Mie.
#[derive(Debug, Clone, Copy)]
pub struct SaftVrMieComponent {
    pub ident: &'static str,
    pub m: f64,
    pub sigma: f64,
    pub eps_depth_divk: f64,
    pub lambda_a: f64,
    pub lambda_r: f64,
    pub mass: f64,
    pub fh_order: i32,
    pub reference: &'static str,
    pub molar_mass: f64,
}

/// Thermodynamic Helmholtz energy terms for SAFT-VR Mie.
#[derive(Debug, Clone, Copy)]
pub struct SaftVrMieHelmholtzTerms {
    pub a_hs: f64,
    pub a1: f64,
    pub a2: f64,
    pub a3: f64,
    pub a_chain: f64,
    pub a_res: f64,
}

/// Pressure jet bundle with first volume, density, and temperature derivatives.
#[derive(Debug, Clone, Copy)]
pub struct SaftVrMiePressureJet {
    pub pressure: f64,
    pub compressibility: f64,
    pub residual_helmholtz: f64,
    pub dp_drho: f64,
    pub dp_dv: f64,
    pub dp_dt: f64,
}

/// Complete thermodynamic state properties bundle at (T, rho).
///
/// `s_residual` and `g_residual` use the TP-referenced residual basis
/// (`s_res = -R*(T*F_T + F - ln Z)`), matching [`SaftVrMieTpState`], while
/// `s_total` is the true total molar entropy at `(T, rho)` (agreeing with
/// [`SaftVrMie::state_tp`] at the same state and satisfying `Delta H = T Delta S`
/// at coexistence).
#[derive(Debug, Clone)]
pub struct SaftVrMieState {
    pub pressure: f64,
    pub compressibility: f64,
    pub residual_helmholtz: f64,
    pub u_residual: f64,
    pub h_residual: f64,
    /// Residual molar entropy `[J/(mol K)]` on the TP residual basis.
    pub s_residual: f64,
    pub a_residual: f64,
    /// Residual molar Gibbs energy `[J/mol]` on the TP residual basis.
    pub g_residual: f64,
    pub h_total: f64,
    /// Total molar entropy `[J/(mol K)]` at `(T, rho)`.
    pub s_total: f64,
    pub cv_residual: f64,
    pub cv_total: f64,
    pub cp_total: f64,
    pub dp_drho: f64,
    pub dp_dv: f64,
    pub dp_dt: f64,
    pub sound_speed: f64,
    pub ln_phi: Vec<f64>,
}

/// Pure-fluid vapor-liquid saturation equilibrium result.
#[derive(Debug, Clone, Copy)]
pub struct SaftVrMieSaturation {
    pub pressure: f64,
    pub temperature: f64,
    pub liquid_volume: f64,
    pub vapor_volume: f64,
    pub enthalpy_vaporization: f64,
    pub entropy_vaporization: f64,
}

/// Scalar property with exact first `(T, P, n)` derivatives.
///
/// Following ThermoPack's Fortran TP conventions, the mole-derivative field
/// `dn` holds one of two bases depending on the property on [`SaftVrMieTpState`]:
/// - For molar properties of extensive quantities (`specific_volume`,
///   `residual_enthalpy`, `enthalpy`, `residual_entropy`, `entropy`, and
///   `residual_gibbs`), `dn[i]` is the **partial molar property**
///   `d(N * q) / dn_i` at constant `(T, P, n_{j != i})` (homogeneous of degree
///   `0` in total moles `N`, satisfying `sum_i x_i * dn[i] = value`). Call
///   [`SaftVrMieDerivatives::intensive_dn`] to obtain the raw-mole derivative
///   `dq / dn_i = (dn[i] - value) / N` of the intensive molar value.
/// - For `compressibility_factor`, `dn[i]` is directly the **intensive**
///   raw-mole derivative `dZ / dn_i` at constant `(T, P, n_{j != i})`
///   (homogeneous of degree `-1` in `N`, matching Fortran `eos::zfac`).
#[derive(Debug, Clone)]
pub struct SaftVrMieDerivatives {
    /// Property value (intensive molar basis).
    pub value: f64,
    /// Temperature derivative at constant `(P, n)`.
    pub dt: f64,
    /// Pressure derivative at constant `(T, n)`.
    pub dp: f64,
    /// Mole derivatives at constant `(T, P, n_{j != i})`: partial molar
    /// `d(N * q) / dn_i` for molar volume/enthalpy/entropy/Gibbs fields, or
    /// intensive `dZ / dn_i` for `compressibility_factor`.
    pub dn: Vec<f64>,
}

impl SaftVrMieDerivatives {
    /// Convert partial molar derivatives `d(N * q) / dn_i` stored in `self.dn`
    /// into raw-mole derivatives of the intensive molar value:
    /// `dq / dn_i = (d(N * q) / dn_i - q) / N`.
    ///
    /// Use this on partial-molar fields of [`SaftVrMieTpState`]
    /// (`specific_volume`, `residual_enthalpy`, `enthalpy`, `residual_entropy`,
    /// `entropy`, `residual_gibbs`) when differentiating the intensive molar
    /// quantity `q(T, P, n)` directly with respect to raw mole numbers `n_i`.
    pub fn intensive_dn(&self, total_moles: f64) -> Result<Vec<f64>> {
        positive(total_moles, "total moles")?;
        let inv_n = 1.0 / total_moles;
        Ok(self
            .dn
            .iter()
            .map(|&q_bar| (q_bar - self.value) * inv_n)
            .collect())
    }
}

/// Fugacity coefficients `ln(phi_i)` and exact first `(T, P, n)` derivatives.
#[derive(Debug, Clone)]
pub struct SaftVrMieFugacityDerivatives {
    /// Natural logarithm of fugacity coefficients `ln(phi_i)`.
    pub ln_phi: Vec<f64>,
    /// Temperature derivatives `d(ln phi_i)/dT` at constant `(P, n)`.
    pub dt: Vec<f64>,
    /// Pressure derivatives `d(ln phi_i)/dP` at constant `(T, n)`.
    pub dp: Vec<f64>,
    /// Row-major `nc x nc` intensive raw-mole derivatives `d(ln phi_i)/dn_j`
    /// at constant `(T, P)` (homogeneous of degree `-1` in total moles).
    pub dn: Vec<f64>,
}

/// Complete TP thermodynamic state bundle with analytical `(T, P, n)` derivatives.
///
/// Mole-derivative basis contract (matching ThermoPack Fortran `eos`):
/// - `compressibility_factor.dn` and `ln_fugacity.dn` are intensive raw-mole
///   derivatives (`dZ/dn_i`, `d(ln phi_i)/dn_j`, homogeneous of degree `-1`).
/// - `specific_volume.dn`, `residual_enthalpy.dn`, `enthalpy.dn`,
///   `residual_entropy.dn`, `entropy.dn`, and `residual_gibbs.dn` are partial
///   molar quantities `d(N * q) / dn_i` (homogeneous of degree `0`). Use
///   [`SaftVrMieDerivatives::intensive_dn`] to convert any partial molar vector
///   to raw-mole derivatives `dq / dn_i` of the intensive molar value.
#[derive(Debug, Clone)]
pub struct SaftVrMieTpState {
    /// Molar volume `[m^3/mol]`.
    pub volume: f64,
    /// Compressibility factor `Z`, `(T, P)` derivatives, and intensive raw-mole `dZ/dn_i`.
    pub compressibility_factor: SaftVrMieDerivatives,
    /// Molar volume `v`, `(T, P)` derivatives, and partial molar volumes `d(N*v)/dn_i`.
    pub specific_volume: SaftVrMieDerivatives,
    /// Residual molar enthalpy `[J/mol]`, `(T, P)` derivatives, and partial molar `d(N*h_res)/dn_i`.
    pub residual_enthalpy: SaftVrMieDerivatives,
    /// Total molar enthalpy `[J/mol]`, `(T, P)` derivatives, and partial molar `d(N*h)/dn_i`.
    pub enthalpy: SaftVrMieDerivatives,
    /// Residual molar entropy `[J/(mol K)]` (TP basis), `(T, P)` derivatives, and partial molar `d(N*s_res)/dn_i`.
    pub residual_entropy: SaftVrMieDerivatives,
    /// Total molar entropy `[J/(mol K)]`, `(T, P)` derivatives, and partial molar `d(N*s)/dn_i`.
    pub entropy: SaftVrMieDerivatives,
    /// Residual molar Gibbs energy `[J/mol]` (TP basis), `(T, P)` derivatives, and partial molar `d(N*g_res)/dn_i`.
    pub residual_gibbs: SaftVrMieDerivatives,
    /// Fugacity coefficients `ln(phi_i)` and intensive `(T, P, n)` derivatives.
    pub ln_fugacity: SaftVrMieFugacityDerivatives,
    /// Isobaric heat capacity `C_p` `[J/(mol K)]`.
    pub cp: f64,
    /// Isochoric heat capacity `C_v` `[J/(mol K)]`.
    pub cv: f64,
    /// Single-phase speed of sound `[m/s]`.
    pub speed_of_sound: f64,
}

/// Reference presets from `saftvrmie_options::set_hs_reference`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum HardSphereReference {
    /// Default Lafitte additive hard-sphere reference (`LAFITTE_HS_REF = 1`).
    #[default]
    Lafitte,
    /// Single-component effective hard-sphere reference (`SINGLE_COMP_HS_REF = 2`).
    SingleComponent,
    /// Additive reference with Leonard non-additive correction (`ADDITIVE_HS_REF = 3`).
    Additive,
    /// Non-additive Santos hard-sphere reference (`NON_ADD_HS_REF = 4`).
    NonAdditive,
    /// Additive reference with exact cross diameters and extra term (`ADDITIVE_EXACT_HS_REF = 5`).
    AdditiveExact,
}

/// Toggle switches for individual SAFT-VR Mie Helmholtz contributions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MieTerms {
    /// Hard-sphere reference contribution (`enable_hs`).
    pub hard_sphere: bool,
    /// First-order dispersion contribution (`enable_A1`).
    pub first: bool,
    /// Second-order dispersion contribution (`enable_A2`).
    pub second: bool,
    /// Third-order dispersion contribution (`enable_A3`).
    pub third: bool,
    /// Chain formation contribution (`enable_chain`).
    pub chain: bool,
}
impl Default for MieTerms {
    fn default() -> Self {
        Self {
            hard_sphere: true,
            first: true,
            second: true,
            third: true,
            chain: true,
        }
    }
}

/// Packing-fraction mixing rule (`saftvrmie_options::zeta_mixing_rule`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackingRule {
    /// Lafitte quadratic packing rule (`ZETA_LAFITTE = 1`).
    Lafitte,
    /// Linear volume-fraction mixing rule (`ZETA_LINEAR = 2`).
    Linear,
    /// Leonard cubed mean-diameter rule (`ZETA_LEONARD = 3`).
    Leonard,
}

/// Second-order quantum correction specification (`saftvrmie_options::quantum_correction_spec`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QuantumVariant {
    /// Standard Feynman-Hibbs prefactor 1/2 (`quantum_correction_spec = 0`).
    FeynmanHibbs,
    /// Jaen-Kahn prefactor 9/10 (`quantum_correction_spec = 1`).
    JaenKahn,
}
impl QuantumVariant {
    fn second_factor(self) -> f64 {
        match self {
            Self::FeynmanHibbs => 0.5,
            Self::JaenKahn => 0.9,
        }
    }
}

/// Independent submodel selections from `saftvrmie_options`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MieOptions {
    /// Packing-fraction mixing rule.
    pub packing: PackingRule,
    /// Leonard extra hard-sphere non-additivity correction (`enable_hs_extra`).
    pub extra_hard_sphere: bool,
    /// Exact effective cross-potential parameters (`exact_crosspot_eff`).
    pub exact_cross_potential: bool,
    /// Geometric epsilon combining rule (`!use_epsrule_Lafitte`).
    pub geometric_epsilon: bool,
    /// BMCSL isothermal compressibility (`Khs_EoS = KHS_EOS_BMCSL`).
    pub bmcsl_compressibility: bool,
    /// Kolafa-Labik-Malijevsky pure hard-sphere EoS (`pure_hs_EoS = PURE_HS_CSK`).
    pub kolafa_pure_reference: bool,
    /// Pair-sigma prefactor in third-order dispersion (`a3_model = A3_SIJ_PREFAC`).
    pub pair_prefactor_a3: bool,
    /// Second-order quantum correction prefactor variant.
    pub quantum_variant: QuantumVariant,
    /// Include quantum corrections in second-order dispersion (`quantum_correct_A2`).
    pub quantum_correct_a2: bool,
}
impl Default for MieOptions {
    fn default() -> Self {
        Self {
            packing: PackingRule::Lafitte,
            extra_hard_sphere: false,
            exact_cross_potential: true,
            geometric_epsilon: false,
            bmcsl_compressibility: false,
            kolafa_pure_reference: false,
            pair_prefactor_a3: false,
            quantum_variant: QuantumVariant::FeynmanHibbs,
            quantum_correct_a2: true,
        }
    }
}

/// SAFT-VR Mie and SAFT-VRQ Mie model instance.
#[derive(Debug, Clone)]
pub struct SaftVrMie {
    components: Vec<SaftVrMieComponent>,
    assoc_schemes: Vec<AssocScheme>,
    assoc_eps: Vec<f64>,
    assoc_beta: Vec<f64>,
    association: AssociationParams,
    integral_coefficients: Vec<[f64; 66]>,
    kij: Vec<f64>,
    lij: Vec<f64>,
    gamma_ij: Vec<f64>,
    lambda_a_ij: Vec<f64>,
    lambda_r_ij: Vec<f64>,
    sigma_ij: Vec<f64>,
    eps_divk_ij: Vec<f64>,
    c_ij: Vec<f64>,
    alpha_ij: Vec<f64>,
    f_alpha_ij: Vec<[f64; 6]>,
    d_feynhibbs_param_ij: Vec<f64>,
    ideal: Vec<IdealGas>,
    exact_binary_dhs: bool,
    reference: HardSphereReference,
    terms: MieTerms,
    options: MieOptions,
    cutoff: Option<(f64, bool)>,
}

impl SaftVrMie {
    /// Create a new SAFT-VR Mie model instance with default parameters.
    pub fn new(comps: &[&str]) -> Result<Self> {
        Self::new_with_options(comps, "Default")
    }

    /// Create a new SAFT-VR Mie model instance with a specific parameter reference.
    pub fn new_with_options(comps: &[&str], parameter_reference: &str) -> Result<Self> {
        if comps.is_empty() {
            return Err(Error::InvalidInput("components list must not be empty"));
        }

        let nc = comps.len();
        let mut components = Vec::with_capacity(nc);
        let mut ideal = Vec::with_capacity(nc);

        let mut schemes = Vec::with_capacity(nc);
        let mut eps_assoc = Vec::with_capacity(nc);
        let mut kab = Vec::with_capacity(nc);
        for &ident in comps {
            let record = saft_vrmie_data::lookup_pure(ident, Some(parameter_reference))
                .ok_or_else(|| {
                    Error::UnknownComponent(format!("SAFT-VR Mie record for {ident}"))
                })?;
            schemes.push(record.scheme);
            eps_assoc.push(record.eps_assoc_divk * crate::R);
            kab.push(record.kab);
            let comp_info = component(ident)?;
            let molar_mass = comp_info.molar_mass; // g/mol

            components.push(SaftVrMieComponent {
                ident: record.ident,
                m: record.m,
                sigma: record.sigma,
                eps_depth_divk: record.eps_depth_divk,
                lambda_a: record.lambda_a,
                lambda_r: record.lambda_r,
                mass: record.mass,
                fh_order: record.fh_order,
                reference: record.reference,
                molar_mass,
            });

            ideal.push(IdealGas::new(comp_info)?);
        }

        let association = AssociationParams::for_saft_vr_mie(
            &schemes,
            &eps_assoc,
            &kab,
            DeltaCombRule::Standard,
        )?;
        let integral_coefficients = if association.has_sites() {
            vec![[0.0; 66]; nc * nc]
        } else {
            Vec::new()
        };
        let mut s = Self {
            assoc_schemes: schemes,
            assoc_eps: eps_assoc,
            assoc_beta: kab,
            association,
            integral_coefficients,
            components,
            kij: vec![0.0; nc * nc],
            lij: vec![0.0; nc * nc],
            gamma_ij: vec![0.0; nc * nc],
            lambda_a_ij: vec![0.0; nc * nc],
            lambda_r_ij: vec![0.0; nc * nc],
            sigma_ij: vec![0.0; nc * nc],
            eps_divk_ij: vec![0.0; nc * nc],
            c_ij: vec![0.0; nc * nc],
            alpha_ij: vec![0.0; nc * nc],
            f_alpha_ij: vec![[0.0; 6]; nc * nc],
            d_feynhibbs_param_ij: vec![0.0; nc * nc],
            ideal,
            exact_binary_dhs: false,
            reference: HardSphereReference::Lafitte,
            terms: MieTerms::default(),
            options: MieOptions::default(),
            cutoff: None,
        };

        // Initialize interaction parameters from database
        for i in 0..nc {
            for j in 0..nc {
                let idx = i * nc + j;
                if i != j {
                    s.kij[idx] = saft_vrmie_data::lookup_kij(
                        s.components[i].ident,
                        s.components[j].ident,
                        Some(parameter_reference),
                    );
                    s.lij[idx] = saft_vrmie_data::lookup_lij(
                        s.components[i].ident,
                        s.components[j].ident,
                        Some(parameter_reference),
                    );
                }
            }
        }

        s.recompute_cross_parameters();
        Ok(s)
    }

    /// Initialize the Aasen FH0, FH1 or FH2 model with the Santos non-additive reference.
    pub fn new_quantum(comps: &[&str], order: u8, reference: Option<&str>) -> Result<Self> {
        Self::new_quantum_with_hs_ref(comps, order, false, reference)
    }

    /// Initialize SAFT-VRQ Mie with explicit `additive_hs_ref` selection matching
    /// `eoslibinit::init_quantum_saftvrmie`.
    pub fn new_quantum_with_hs_ref(
        comps: &[&str],
        order: u8,
        additive_hs_ref: bool,
        reference: Option<&str>,
    ) -> Result<Self> {
        let default_ref = match order {
            0 => "AASEN2019-FH0",
            1 => "AASEN2019-FH1",
            2 => "AASEN2019-FH2",
            _ => return Err(Error::InvalidInput("FH order must be 0, 1 or 2")),
        };
        let mut model = Self::new_with_options(
            comps,
            reference
                .filter(|r| !r.eq_ignore_ascii_case("default"))
                .unwrap_or(default_ref),
        )?;
        if model
            .components
            .iter()
            .any(|c| c.fh_order != i32::from(order) || (order > 0 && (c.mass <= 0. || c.m != 1.)))
        {
            return Err(Error::InvalidInput(
                "component parameter set incompatible with quantum order",
            ));
        }
        model.set_hard_sphere_reference(if additive_hs_ref {
            HardSphereReference::AdditiveExact
        } else {
            HardSphereReference::NonAdditive
        });
        Ok(model)
    }

    /// Active hard-sphere reference preset.
    pub fn hard_sphere_reference(&self) -> HardSphereReference {
        self.reference
    }

    /// Select a hard-sphere reference preset matching `saftvrmie_options::set_hs_reference`.
    ///
    /// Note that `HardSphereReference::NonAdditive` resets `kolafa_pure_reference`
    /// to `false` (`PURE_HS_CS`) and leaves the active `options.packing` rule
    /// unchanged, matching Fortran `set_hs_reference(NON_ADD_HS_REF)`.
    pub fn set_hard_sphere_reference(&mut self, reference: HardSphereReference) {
        self.reference = reference;
        match reference {
            HardSphereReference::Lafitte => {
                self.options.packing = PackingRule::Lafitte;
                self.exact_binary_dhs = false;
                self.options.extra_hard_sphere = false;
            }
            HardSphereReference::SingleComponent => {
                self.options.packing = PackingRule::Leonard;
                self.exact_binary_dhs = true;
                self.options.extra_hard_sphere = false;
            }
            HardSphereReference::Additive => {
                self.options.packing = PackingRule::Lafitte;
                self.exact_binary_dhs = false;
                self.options.extra_hard_sphere = true;
            }
            HardSphereReference::NonAdditive => {
                self.options.kolafa_pure_reference = false;
                self.exact_binary_dhs = true;
                self.options.extra_hard_sphere = false;
            }
            HardSphereReference::AdditiveExact => {
                self.options.packing = PackingRule::Lafitte;
                self.exact_binary_dhs = true;
                self.options.extra_hard_sphere = true;
            }
        }
    }

    /// Active independent submodel options.
    pub fn options(&self) -> MieOptions {
        self.options
    }

    /// Update independent submodel options, validating model consistency.
    pub fn set_options(&mut self, options: MieOptions) -> Result<()> {
        if options.bmcsl_compressibility && self.components.iter().any(|c| c.m != 1.) {
            return Err(Error::InvalidInput(
                "BMCSL compressibility requires spherical molecules",
            ));
        }
        if self.reference == HardSphereReference::SingleComponent
            && options.packing != PackingRule::Leonard
        {
            return Err(Error::InvalidInput(
                "single component reference requires Leonard packing",
            ));
        }
        self.options = options;
        self.recompute_cross_parameters();
        Ok(())
    }

    /// Active Helmholtz term toggles.
    pub fn terms(&self) -> MieTerms {
        self.terms
    }

    /// Enable or disable individual SAFT-VR Mie Helmholtz terms.
    pub fn set_terms(&mut self, terms: MieTerms) {
        self.terms = terms;
    }

    /// Active reduced potential cutoff `(r_cut, shift_enabled)`, or `None` for full Mie range.
    pub fn cutoff(&self) -> Option<(f64, bool)> {
        self.cutoff
    }

    /// Reduced cutoff distance in sigma units. `None` restores the infinite-range potential.
    pub fn set_cutoff(&mut self, cutoff: Option<(f64, bool)>) -> Result<()> {
        if let Some((r, _)) = cutoff
            && (self.num_components() != 1 || !r.is_finite() || r <= 0.)
        {
            return Err(Error::InvalidInput(
                "truncation requires a pure fluid and positive finite cutoff",
            ));
        }
        self.cutoff = cutoff;
        Ok(())
    }

    /// Whether the model has any active associating sites.
    pub fn has_association(&self) -> bool {
        self.association.has_sites()
    }

    /// Override association parameters `(scheme, eps [J/mol], beta [m^3])` for one component.
    pub fn set_assoc_params(
        &mut self,
        index: usize,
        scheme: AssocScheme,
        eps: f64,
        beta: f64,
    ) -> Result<()> {
        if index >= self.num_components() {
            return Err(Error::InvalidInput("component index out of bounds"));
        }
        let mut schemes = self.assoc_schemes.clone();
        let mut eps_vec = self.assoc_eps.clone();
        let mut beta_vec = self.assoc_beta.clone();
        schemes[index] = scheme;
        eps_vec[index] = eps;
        beta_vec[index] = beta;
        let association = AssociationParams::for_saft_vr_mie(
            &schemes,
            &eps_vec,
            &beta_vec,
            self.association.delta_comb_rule,
        )?;
        self.assoc_schemes = schemes;
        self.assoc_eps = eps_vec;
        self.assoc_beta = beta_vec;
        self.association = association;
        let nc = self.num_components();
        if self.association.has_sites() {
            if self.integral_coefficients.len() != nc * nc {
                self.integral_coefficients = vec![[0.0; 66]; nc * nc];
            }
            self.recompute_cross_parameters();
        } else {
            self.integral_coefficients.clear();
        }
        Ok(())
    }

    /// Select the cross-association Delta combining rule (`Standard` or `Elliot`).
    pub fn set_delta_comb_rule(&mut self, rule: DeltaCombRule) -> Result<()> {
        self.association = AssociationParams::for_saft_vr_mie(
            &self.assoc_schemes,
            &self.assoc_eps,
            &self.assoc_beta,
            rule,
        )?;
        Ok(())
    }

    /// Solve for equilibrium non-bonded association site fractions `X_k` at `(T, V, n)`.
    pub fn solve_assoc_x_k(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<Vec<f64>> {
        positive(volume, "total volume")?;
        let state = self.prepare(temperature, moles)?;
        if !self.association.has_sites() {
            return Ok(Vec::new());
        }
        let total: f64 = moles.iter().sum();
        let v_unit = volume / total;
        if !v_unit.is_finite() || v_unit <= 0.0 {
            return Err(Error::InvalidInput(
                "molar volume must be finite and strictly positive",
            ));
        }
        solve_for_x_k(
            &self.association,
            &MieRdf::new(&state, 1.0 / v_unit),
            temperature,
            v_unit,
            &state.mole_fractions,
            None,
        )
    }

    /// `(is_fmt_consistent, nonadditive_extra_enabled)` matching `test_fmt_compatibility`.
    pub fn fmt_compatibility(&self) -> (bool, bool) {
        if self.num_components() == 1 {
            return (true, false);
        }
        let additive = matches!(
            self.reference,
            HardSphereReference::Lafitte
                | HardSphereReference::Additive
                | HardSphereReference::AdditiveExact
        ) && self.options.packing == PackingRule::Lafitte;
        let extra = self.options.extra_hard_sphere;
        (additive && (!extra || self.exact_binary_dhs), extra)
    }

    /// Extensive residual F=Ares/(RT) and its ten TVN blocks at raw mole amounts.
    pub fn residual_helmholtz_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<SaftVrMieHelmholtz> {
        self.derivatives(temperature, volume, moles, EvalMode::Total)
    }

    /// Reduced molar dispersion energy and its TVN derivatives, as used by DFT.
    pub fn dispersion_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<SaftVrMieHelmholtz> {
        self.derivatives(temperature, volume, moles, EvalMode::Dispersion)
    }

    /// Reduced molar hard-sphere energy and its TVN derivatives, as used by DFT.
    pub fn hard_sphere_tvn(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<SaftVrMieHelmholtz> {
        self.derivatives(temperature, volume, moles, EvalMode::HardSphere)
    }

    fn derivatives_unit(
        &self,
        state: &PreparedSaftVrMie<'_>,
        v_unit: f64,
        mode: EvalMode,
    ) -> Result<SaftVrMieHelmholtz> {
        let nc = self.num_components();
        let rho = 1.0 / v_unit;
        let empty = AssociationParams::empty(nc);
        let model = if mode == EvalMode::Total {
            &self.association
        } else {
            &empty
        };
        let mut out = assoc_helmholtz_all(
            model,
            &MieRdf::new(state, rho),
            state.temperature,
            v_unit,
            &state.mole_fractions,
        )?;
        for i in 0..nc + 2 {
            for j in i..nc + 2 {
                let terms = state.eval_second(rho, i, j, mode)?;
                let f = match mode {
                    EvalMode::Total => terms.total,
                    EvalMode::Dispersion => terms.dispersion,
                    EvalMode::HardSphere => terms.hard_sphere,
                };
                if i == j {
                    match i {
                        0 => {
                            out.f += f.f0;
                            out.f_t += f.f1;
                        }
                        1 => out.f_v += f.f1,
                        _ => out.f_n[i - 2] += f.f1,
                    }
                }
                match (i, j) {
                    (0, 0) => out.f_tt += f.f12,
                    (0, 1) => out.f_tv += f.f12,
                    (0, _) => out.f_tn[j - 2] += f.f12,
                    (1, 1) => out.f_vv += f.f12,
                    (1, _) => out.f_vn[j - 2] += f.f12,
                    _ => {
                        out.f_nn[(i - 2) * nc + j - 2] += f.f12;
                        if i != j {
                            out.f_nn[(j - 2) * nc + i - 2] += f.f12;
                        }
                    }
                }
            }
        }
        Ok(out)
    }

    fn derivatives(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
        mode: EvalMode,
    ) -> Result<SaftVrMieHelmholtz> {
        positive(volume, "total volume")?;
        let state = self.prepare(temperature, moles)?;
        let total: f64 = moles.iter().sum();
        let v_unit = volume / total;
        if !v_unit.is_finite() || v_unit <= 0.0 {
            return Err(Error::InvalidInput(
                "molar volume must be finite and strictly positive",
            ));
        }
        let mut out = self.derivatives_unit(&state, v_unit, mode)?;
        if mode == EvalMode::Total {
            out.f *= total;
            out.f_t *= total;
            out.f_tt *= total;
            out.f_vv /= total;
            for val in &mut out.f_vn {
                *val /= total;
            }
            for val in &mut out.f_nn {
                *val /= total;
            }
        } else {
            // Molar quantities are homogeneous of degree 0 in (V, n).
            out.f_v /= total;
            out.f_tv /= total;
            out.f_vv = (out.f_vv / total) / total;
            for val in out.f_n.iter_mut().chain(&mut out.f_tn) {
                *val /= total;
            }
            for val in out.f_vn.iter_mut().chain(&mut out.f_nn) {
                *val = (*val / total) / total;
            }
        }
        if ![out.f, out.f_t, out.f_v, out.f_tt, out.f_tv, out.f_vv]
            .iter()
            .all(|v| v.is_finite())
            || !out.f_n.iter().all(|v| v.is_finite())
            || !out.f_tn.iter().all(|v| v.is_finite())
            || !out.f_vn.iter().all(|v| v.is_finite())
            || !out.f_nn.iter().all(|v| v.is_finite())
        {
            return Err(Error::InvalidInput(
                "non-finite SAFT-VR Mie Helmholtz derivative",
            ));
        }
        Ok(out)
    }

    /// Solve for the molar volume `[m^3/mol]` at `(T, P, n)` for the requested `phase`.
    pub fn volume_tp(
        &self,
        temperature: f64,
        pressure: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<f64> {
        let prepared = self.prepare(temperature, moles)?;
        prepared.volume(pressure, phase)
    }

    /// Evaluate the complete TP thermodynamic state with exact first `(T, P, n)`
    /// derivatives from all ten TVN Helmholtz blocks at `(T, P, n)` for `phase`.
    pub fn state_tp(
        &self,
        temperature: f64,
        pressure: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<SaftVrMieTpState> {
        positive(temperature, "temperature")?;
        positive(pressure, "pressure")?;
        let prepared = self.prepare(temperature, moles)?;
        let total: f64 = moles.iter().sum();
        let v_molar = prepared.volume(pressure, phase)?;
        let v_total = v_molar * total;
        if !v_total.is_finite() || v_total <= 0.0 {
            return Err(Error::InvalidInput(
                "total volume overflowed or underflowed",
            ));
        }

        let nc = self.num_components();
        let h = self.derivatives_unit(&prepared, v_molar, EvalMode::Total)?;
        let t = temperature;
        let v = v_molar;
        let inv_n = 1.0 / total;
        let p_eval = RGAS * t / v - RGAS * t * h.f_v;
        let p_v = -RGAS * t / (v * v) - RGAS * t * h.f_vv;
        let p_t = p_eval / t - RGAS * t * h.f_tv;
        let p_n: Vec<f64> = h
            .f_vn
            .iter()
            .map(|&fvn| RGAS * t / v - RGAS * t * fvn)
            .collect();

        let v_t = -p_t / p_v;
        let v_p = 1.0 / p_v;
        let v_n_partial: Vec<f64> = p_n.iter().map(|&pni| -pni / p_v).collect();

        let z = p_eval * v / (RGAS * t);
        let ln_z = z.ln();
        let z_t = (p_eval / (RGAS * t)) * (v_t - v / t);
        let z_p = (v + p_eval * v_p) / (RGAS * t);
        let z_dn: Vec<f64> = v_n_partial
            .iter()
            .map(|&vni| (p_eval / (RGAS * t)) * (vni - v) * inv_n)
            .collect();

        let h_res = -RGAS * t * t * h.f_t + p_eval * v - RGAS * t;
        let cv_res = -RGAS * (2.0 * t * h.f_t + t * t * h.f_tt);
        let cp_res = cv_res - t * p_t * p_t / p_v - RGAS;
        let h_res_t = cp_res;
        let h_res_p = v - t * v_t;
        let mut h_res_dn = Vec::with_capacity(nc);
        for i in 0..nc {
            h_res_dn.push(-RGAS * t * t * h.f_tn[i] + t * p_t * v_n_partial[i] - RGAS * t);
        }

        let s_res = -RGAS * (h.f + t * h.f_t - ln_z);
        let s_res_t = cp_res / t;
        let s_res_p = -v_t + RGAS / p_eval;
        let mut s_res_dn = Vec::with_capacity(nc);
        for i in 0..nc {
            s_res_dn.push(-RGAS * (h.f_n[i] + t * h.f_tn[i] - ln_z) + p_t * v_n_partial[i] - RGAS);
        }

        let g_res = h_res - t * s_res;
        let g_res_t = -s_res;
        let g_res_p = v - RGAS * t / p_eval;
        let g_res_dn: Vec<f64> = h.f_n.iter().map(|&fni| RGAS * t * (fni - ln_z)).collect();

        let mut h_id = 0.0;
        let mut cp_id = 0.0;
        let mut s_id = 0.0;
        let mut h_id_comp = Vec::with_capacity(nc);
        let mut s_id_comp = Vec::with_capacity(nc);
        for i in 0..nc {
            let id = self.ideal[i].properties(t, p_eval)?;
            let xi = prepared.mole_fractions[i];
            h_id += xi * id.enthalpy;
            cp_id += xi * id.cp;
            if xi > 0.0 {
                s_id += xi * (id.entropy - RGAS * xi.ln());
            }
            h_id_comp.push(id.enthalpy);
            s_id_comp.push(id.entropy - RGAS * xi.max(1.0e-100).ln());
        }
        let h_tot_dn: Vec<f64> = h_res_dn
            .iter()
            .zip(&h_id_comp)
            .map(|(r, id)| r + id)
            .collect();
        let s_tot_dn: Vec<f64> = s_res_dn
            .iter()
            .zip(&s_id_comp)
            .map(|(r, id)| r + id)
            .collect();

        let mut ln_phi = Vec::with_capacity(nc);
        let mut dlnphi_dt = Vec::with_capacity(nc);
        let mut dlnphi_dp = Vec::with_capacity(nc);
        let mut dlnphi_dn = vec![0.0; nc * nc];
        for i in 0..nc {
            ln_phi.push(h.f_n[i] - ln_z);
            dlnphi_dt.push(h.f_tn[i] + 1.0 / t - v_n_partial[i] * p_t / (RGAS * t));
            dlnphi_dp.push(v_n_partial[i] / (RGAS * t) - 1.0 / p_eval);
            for j in 0..nc {
                let d_unit = h.f_nn[i * nc + j] + 1.0 + p_n[i] * p_n[j] / (RGAS * t * p_v);
                dlnphi_dn[i * nc + j] = d_unit * inv_n;
            }
        }

        let cv = cp_id - RGAS + cv_res;
        let cp = cp_id + cp_res;
        let dp_drho = -p_v * v * v;
        let mw_kg = prepared.molar_mass_mix * 1.0e-3;
        let sound_speed = if dp_drho > 0.0 && cp > 0.0 && cv > 0.0 && mw_kg > 0.0 {
            ((cp / cv) * dp_drho / mw_kg).sqrt()
        } else {
            0.0
        };

        let all_finite = [
            v,
            z,
            z_t,
            z_p,
            v_t,
            v_p,
            h_res,
            h_res_t,
            h_res_p,
            s_res,
            s_res_t,
            s_res_p,
            g_res,
            g_res_t,
            g_res_p,
            cp,
            cv,
            sound_speed,
        ]
        .iter()
        .all(|q| q.is_finite())
            && z_dn.iter().all(|q| q.is_finite())
            && v_n_partial.iter().all(|q| q.is_finite())
            && h_res_dn.iter().all(|q| q.is_finite())
            && s_res_dn.iter().all(|q| q.is_finite())
            && g_res_dn.iter().all(|q| q.is_finite())
            && ln_phi.iter().all(|q| q.is_finite())
            && dlnphi_dt.iter().all(|q| q.is_finite())
            && dlnphi_dp.iter().all(|q| q.is_finite())
            && dlnphi_dn.iter().all(|q| q.is_finite());
        if !all_finite {
            return Err(Error::InvalidInput(
                "non-finite SAFT-VR Mie TP thermodynamic state",
            ));
        }

        Ok(SaftVrMieTpState {
            volume: v,
            compressibility_factor: SaftVrMieDerivatives {
                value: z,
                dt: z_t,
                dp: z_p,
                dn: z_dn,
            },
            specific_volume: SaftVrMieDerivatives {
                value: v,
                dt: v_t,
                dp: v_p,
                dn: v_n_partial,
            },
            residual_enthalpy: SaftVrMieDerivatives {
                value: h_res,
                dt: h_res_t,
                dp: h_res_p,
                dn: h_res_dn,
            },
            enthalpy: SaftVrMieDerivatives {
                value: h_id + h_res,
                dt: cp,
                dp: h_res_p,
                dn: h_tot_dn,
            },
            residual_entropy: SaftVrMieDerivatives {
                value: s_res,
                dt: s_res_t,
                dp: s_res_p,
                dn: s_res_dn,
            },
            entropy: SaftVrMieDerivatives {
                value: s_id + s_res,
                dt: cp / t,
                dp: -v_t,
                dn: s_tot_dn,
            },
            residual_gibbs: SaftVrMieDerivatives {
                value: g_res,
                dt: g_res_t,
                dp: g_res_p,
                dn: g_res_dn,
            },
            ln_fugacity: SaftVrMieFugacityDerivatives {
                ln_phi,
                dt: dlnphi_dt,
                dp: dlnphi_dp,
                dn: dlnphi_dn,
            },
            cp,
            cv,
            speed_of_sound: sound_speed,
        })
    }

    /// Evaluate fugacity coefficients and their exact `(T, P, n)` derivatives at `(T, P, n)`.
    pub fn ln_fugacity_derivatives_tp(
        &self,
        temperature: f64,
        pressure: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<SaftVrMieFugacityDerivatives> {
        Ok(self
            .state_tp(temperature, pressure, moles, phase)?
            .ln_fugacity)
    }

    /// Pure-component Barker-Henderson hard-sphere diameters `(d, d_T, d_TT)` at `T`,
    /// matching `saft_interface::calc_hard_sphere_diameter`.
    pub fn hard_sphere_diameters(
        &self,
        temperature: f64,
    ) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
        positive(temperature, "temperature")?;
        let nc = self.num_components();
        let mut d = Vec::with_capacity(nc);
        let mut d_t = Vec::with_capacity(nc);
        let mut d_tt = Vec::with_capacity(nc);
        for i in 0..nc {
            let pair = self.potential(i, i);
            let se = pair.effective_sigma(temperature);
            let (di, dit, ditt) = pair.diameter(temperature, se)?;
            d.push(di);
            d_t.push(dit);
            d_tt.push(ditt);
        }
        Ok((d, d_t, d_tt))
    }

    /// Pair hard-sphere diameter `(d_ij, d_T_ij, d_TT_ij)` at `T`,
    /// matching `saft_interface::calc_hard_sphere_diameter_ij`.
    pub fn hard_sphere_diameter_ij(
        &self,
        i: usize,
        j: usize,
        temperature: f64,
    ) -> Result<(f64, f64, f64)> {
        positive(temperature, "temperature")?;
        let nc = self.num_components();
        if i >= nc || j >= nc {
            return Err(Error::InvalidInput("component index out of bounds"));
        }
        let dummy_x = vec![1.0 / nc as f64; nc];
        let prep = PreparedSaftVrMie::new(self, temperature, dummy_x)?;
        let k = i * nc + j;
        Ok((prep.d[k], prep.d_t[k], prep.d_tt[k]))
    }

    /// Effective quantum/classical collision diameter `(sigma_eff, dT, dTT)` for pair `(i, j)`.
    pub fn effective_sigma_ij(
        &self,
        i: usize,
        j: usize,
        temperature: f64,
    ) -> Result<(f64, f64, f64)> {
        positive(temperature, "temperature")?;
        let nc = self.num_components();
        if i >= nc || j >= nc {
            return Err(Error::InvalidInput("component index out of bounds"));
        }
        let dummy_x = vec![1.0 / nc as f64; nc];
        let prep = PreparedSaftVrMie::new(self, temperature, dummy_x)?;
        let [v, vt, vtt] = prep.effective[i * nc + j].0;
        Ok((v, vt, vtt))
    }

    /// Effective quantum/classical well depth `(eps_eff/kB, dT, dTT)` for pair `(i, j)`.
    pub fn effective_epsilon_ij(
        &self,
        i: usize,
        j: usize,
        temperature: f64,
    ) -> Result<(f64, f64, f64)> {
        positive(temperature, "temperature")?;
        let nc = self.num_components();
        if i >= nc || j >= nc {
            return Err(Error::InvalidInput("component index out of bounds"));
        }
        let dummy_x = vec![1.0 / nc as f64; nc];
        let prep = PreparedSaftVrMie::new(self, temperature, dummy_x)?;
        let [v, vt, vtt] = prep.effective[i * nc + j].1;
        Ok((v, vt, vtt))
    }

    /// Dimensionless van der Waals attractive energy matrix `alpha_ij(T)` (`nc x nc`),
    /// matching `saft_interface::alpha` (including truncated/shifted classical potentials).
    ///
    /// For shifted cutoffs, cutoffs at, below, or unresolvably close to `r_min`
    /// (`r_cut - r_min <= 1e-12 * r_min`) return `Error::NotConverged`. Near `r_min`
    /// (`r_cut <= 2 * r_min`), the shifted well depth `eps_eff = u*(r_cut) + 1` and
    /// inner zero offset `r_min - sigma_eff` are solved relative to `r_min`, using
    /// 21-point Gauss-Kronrod quadrature on shallow wells (`eps_eff < 0.25`). For
    /// deep wells and large finite cutoffs (`r_cut > 2 * r_min`), `u_cut = -u*(r_cut)`
    /// is evaluated directly, the inner zero offset `sigma_eff - 1` is solved near
    /// `sigma = 1`, and the tail uses `r_cut^(3 - lambda)` to avoid `1 - eps_eff`
    /// cancellation and `r_cut^3` overflow.
    pub fn alpha_vdw(&self, temperature: f64) -> Result<Vec<f64>> {
        positive(temperature, "temperature")?;
        let nc = self.num_components();
        if let Some((r_cut, shift)) = self.cutoff {
            if self.components.iter().any(|c| c.fh_order > 0) {
                return Err(Error::Unsupported(
                    "truncated alpha_vdw is not supported for quantum potentials",
                ));
            }
            let mut out = vec![0.0; nc * nc];
            for i in 0..nc {
                for j in i..nc {
                    let idx = i * nc + j;
                    let la = self.lambda_a_ij[idx];
                    let lr = self.lambda_r_ij[idx];
                    let c = self.c_ij[idx];
                    let a = if shift {
                        let r_min = (lr / la).powf(1.0 / (lr - la));
                        let delta_c = r_cut - r_min;
                        if !r_min.is_finite() || !delta_c.is_finite() || delta_c <= 1.0e-12 * r_min
                        {
                            return Err(Error::NotConverged(
                                "cannot resolve shifted potential well for cutoff radius",
                            ));
                        }
                        if r_cut <= 2.0 * r_min {
                            let eps_eff = mie_well_excess_tau(delta_c / r_min, la, lr);
                            if !eps_eff.is_finite() || eps_eff <= 0.0 || eps_eff >= 1.0 {
                                return Err(Error::NotConverged(
                                    "cannot resolve shifted potential well depth for cutoff radius",
                                ));
                            }
                            // Solve g(d) = (u*(r_min - d) + 1) - eps_eff = 0 for the positive
                            // offset d = r_min - sigma_eff in (0, r_min - 1), avoiding catastrophic
                            // subtraction of rs - r_min for shallow wells.
                            let d_max = r_min - 1.0;
                            let mut lo = 0.0_f64;
                            let mut hi = d_max;
                            let mut d_s = delta_c.min(0.5 * d_max);
                            let mut converged = false;
                            for _ in 0..60 {
                                let f = mie_well_excess_tau(-d_s / r_min, la, lr) - eps_eff;
                                if f.abs() <= 1.0e-13 * eps_eff {
                                    converged = true;
                                    break;
                                }
                                if f < 0.0 {
                                    lo = d_s;
                                } else {
                                    hi = d_s;
                                }
                                let df = mie_well_left_deriv(d_s, r_min, la, lr);
                                let mut next = if df.is_finite() && df > 0.0 {
                                    d_s - f / df
                                } else {
                                    0.5 * (lo + hi)
                                };
                                if !next.is_finite() || next <= lo || next >= hi {
                                    next = 0.5 * (lo + hi);
                                }
                                let step = (next - d_s).abs();
                                d_s = next;
                                if step <= 1.0e-14 * d_s || (hi - lo) <= 1.0e-14 * d_s {
                                    let f_chk = mie_well_excess_tau(-d_s / r_min, la, lr) - eps_eff;
                                    converged =
                                        f_chk.is_finite() && f_chk.abs() <= 1.0e-11 * eps_eff;
                                    break;
                                }
                            }
                            if !converged || !d_s.is_finite() || d_s <= 0.0 || d_s >= d_max {
                                return Err(Error::NotConverged(
                                    "cannot solve shifted potential zero for cutoff radius",
                                ));
                            }
                            let sigma_eff = r_min - d_s;
                            if eps_eff < 0.25 {
                                // Integrate the strictly positive integrand
                                // (1 - (u*(r_min + x) + 1) / eps_eff) * (r_min + x)^2 over
                                // x in [-d_s, 0] and [0, delta_c] with 21-point Gauss-Kronrod.
                                let quad_x = |x_lo: f64, x_hi: f64| -> f64 {
                                    let mid = 0.5 * (x_lo + x_hi);
                                    let half = 0.5 * (x_hi - x_lo);
                                    let mut acc = 0.0;
                                    for k in 0..21 {
                                        let x = mid + half * X_NODES_21[k];
                                        let r = r_min + x;
                                        let excess = mie_well_excess_tau(x / r_min, la, lr);
                                        let depth = (1.0 - excess / eps_eff).max(0.0);
                                        acc += W_NODES_21[k] * depth * r * r;
                                    }
                                    half * acc
                                };
                                (quad_x(-d_s, 0.0) + quad_x(0.0, delta_c)) / sigma_eff.powi(3)
                            } else {
                                let u_cut = c * (r_cut.powf(-la) - r_cut.powf(-lr));
                                let rc_a = r_cut.powf(3.0 - la);
                                let rc_r = r_cut.powf(3.0 - lr);
                                let sigma3 = sigma_eff.powi(3);
                                let num = c / (sigma_eff.powf(la - 3.0) * (la - 3.0))
                                    - c / (sigma_eff.powf(lr - 3.0) * (lr - 3.0))
                                    + c * (lr / (3.0 * (lr - 3.0))) * rc_r
                                    - c * (la / (3.0 * (la - 3.0))) * rc_a
                                    + u_cut * sigma3 / 3.0;
                                num / (eps_eff * sigma3)
                            }
                        } else {
                            // Deep-well / large-cutoff path (r_cut > 2 * r_min): evaluate
                            // u_cut = -u*(r_cut) directly, solve u*(1 + ds) + u_cut = 0 for
                            // ds = sigma_eff - 1 >= 0, and evaluate the cutoff tail via
                            // r_cut^(3 - lambda) without 1 - eps_eff cancellation or r_cut^3 overflow.
                            let u_cut = c * (r_cut.powf(-la) - r_cut.powf(-lr));
                            if !u_cut.is_finite() || u_cut < 0.0 || u_cut >= 1.0 {
                                return Err(Error::NotConverged(
                                    "cannot resolve shifted potential well depth for cutoff radius",
                                ));
                            }
                            let eps_eff = 1.0 - u_cut;
                            let init_slope = c * (lr - la);
                            let d_max = r_min - 1.0;
                            let delta_s = if u_cut <= 1.0e-30 * init_slope {
                                0.0
                            } else {
                                let mut lo = 0.0_f64;
                                let mut hi = d_max;
                                let mut ds = (u_cut / init_slope).min(0.5 * d_max);
                                if !ds.is_finite() || ds <= 0.0 {
                                    ds = 0.5 * d_max;
                                }
                                let mut converged = false;
                                for _ in 0..60 {
                                    let f = mie_shifted_unit_value(ds, c, la, lr) + u_cut;
                                    if f.abs() <= 1.0e-13 * u_cut {
                                        converged = true;
                                        break;
                                    }
                                    if f > 0.0 {
                                        lo = ds;
                                    } else {
                                        hi = ds;
                                    }
                                    let df = mie_shifted_unit_slope(ds, c, la, lr);
                                    let mut next = if df.is_finite() && df < 0.0 {
                                        ds - f / df
                                    } else {
                                        0.5 * (lo + hi)
                                    };
                                    if !next.is_finite() || next <= lo || next >= hi {
                                        next = 0.5 * (lo + hi);
                                    }
                                    let step = (next - ds).abs();
                                    ds = next;
                                    if step <= 1.0e-14 * ds || (hi - lo) <= 1.0e-14 * ds {
                                        let f_chk = mie_shifted_unit_value(ds, c, la, lr) + u_cut;
                                        converged =
                                            f_chk.is_finite() && f_chk.abs() <= 1.0e-11 * u_cut;
                                        break;
                                    }
                                }
                                if !converged || !ds.is_finite() || ds < 0.0 || ds >= d_max {
                                    return Err(Error::NotConverged(
                                        "cannot solve shifted potential zero for cutoff radius",
                                    ));
                                }
                                ds
                            };
                            let sigma_eff = 1.0 + delta_s;
                            let rc_a = r_cut.powf(3.0 - la);
                            let rc_r = r_cut.powf(3.0 - lr);
                            let sigma3 = sigma_eff.powi(3);
                            let num = c / (sigma_eff.powf(la - 3.0) * (la - 3.0))
                                - c / (sigma_eff.powf(lr - 3.0) * (lr - 3.0))
                                + c * (lr / (3.0 * (lr - 3.0))) * rc_r
                                - c * (la / (3.0 * (la - 3.0))) * rc_a
                                + u_cut * sigma3 / 3.0;
                            num / (eps_eff * sigma3)
                        }
                    } else {
                        c / (la - 3.0) - c / (lr - 3.0) + c * r_cut.powf(3.0 - lr) / (lr - 3.0)
                            - c * r_cut.powf(3.0 - la) / (la - 3.0)
                    };
                    if !a.is_finite() || (shift && a <= 0.0) {
                        return Err(Error::NotConverged(
                            "non-finite or non-positive truncated van der Waals alpha parameter",
                        ));
                    }
                    out[idx] = a;
                    out[j * nc + i] = a;
                }
            }
            return Ok(out);
        }

        if self.components.iter().all(|c| c.fh_order == 0) {
            return Ok(self.alpha_ij.clone());
        }

        let dummy_x = vec![1.0 / nc as f64; nc];
        let prep = PreparedSaftVrMie::new(self, temperature, dummy_x)?;
        let mut out = vec![0.0; nc * nc];
        for i in 0..nc {
            for j in i..nc {
                let idx = i * nc + j;
                let la = self.lambda_a_ij[idx];
                let lr = self.lambda_r_ij[idx];
                let sigma = self.sigma_ij[idx];
                let sigma_eff = prep.effective[idx].0[0];
                let eps_eff = prep.effective[idx].1[0];
                let order = self.components[i].fh_order as usize;
                let dfh = self.d_feynhibbs_param_ij[idx] / (sigma * sigma * temperature);
                let second_fac = self.options.quantum_variant.second_factor();
                let qf = |lam: f64| {
                    let q1 = lam * (lam - 1.0) * dfh;
                    [1.0, q1, second_fac * (lam + 1.0) * (lam + 2.0) * dfh * q1]
                };
                let qa = qf(la);
                let qr = qf(lr);
                let ratio = sigma / sigma_eff;
                let mut sum = 0.0;
                for k in 0..=order {
                    let a = la + 2.0 * k as f64;
                    let r = lr + 2.0 * k as f64;
                    sum += qa[k] * ratio.powf(a) / (a - 3.0) - qr[k] * ratio.powf(r) / (r - 3.0);
                }
                let val = sum * self.c_ij[idx] * self.eps_divk_ij[idx] / eps_eff;
                if !val.is_finite() {
                    return Err(Error::NotConverged(
                        "non-finite quantum van der Waals alpha parameter",
                    ));
                }
                out[idx] = val;
                out[j * nc + i] = val;
            }
        }
        Ok(out)
    }

    /// Pair interaction potential `u(r)/k_B [K]` and force `-du/dr / k_B [K/m]`
    /// at separation `r [m]` and temperature `T [K]`, matching ThermoPack
    /// `saft_interface::potential` (including Feynman-Hibbs, truncation/shift,
    /// and the default reduced core cap of `500.0 * eps_ij / k_B`).
    ///
    /// Matching Fortran `saft_interface::potential`, the potential value is capped
    /// at `500.0 * eps_divk_ij` (and `NaN` is replaced by that cap), while the
    /// returned force `-du/dr / k_B` remains the uncapped analytical derivative
    /// of the Mie/quantum-Mie potential (or `0.0` for `r > r_cut`). Use
    /// [`SaftVrMie::pair_potential_with_cap`] with `None` for the raw uncapped
    /// potential or `Some(cap)` for a custom reduced cap (`max_pot_val_in`).
    pub fn pair_potential(
        &self,
        i: usize,
        j: usize,
        r: f64,
        temperature: f64,
    ) -> Result<(f64, f64)> {
        self.pair_potential_with_cap(i, j, r, temperature, Some(DEFAULT_MAX_REDUCED_POTENTIAL))
    }

    /// Pair interaction potential `u(r)/k_B [K]` and uncapped analytical force
    /// `-du/dr / k_B [K/m]` with an explicit reduced potential cap
    /// `max_reduced_potential` (in multiples of `eps_ij / k_B`, matching
    /// `saft_interface::potential`'s `max_pot_val_in`), or `None` for the raw
    /// uncapped potential.
    pub fn pair_potential_with_cap(
        &self,
        i: usize,
        j: usize,
        r: f64,
        temperature: f64,
        max_reduced_potential: Option<f64>,
    ) -> Result<(f64, f64)> {
        positive(r, "separation")?;
        positive(temperature, "temperature")?;
        if let Some(cap) = max_reduced_potential
            && !cap.is_finite()
        {
            return Err(Error::InvalidInput("reduced potential cap must be finite"));
        }
        let nc = self.num_components();
        if i >= nc || j >= nc {
            return Err(Error::InvalidInput("component index out of bounds"));
        }
        let pair = self.potential(i, j);
        let (mut u, force) = if let Some((rc, shift)) = self.cutoff {
            let r_cut = rc * pair.sigma;
            if r > r_cut {
                (0.0, 0.0)
            } else {
                let pot = pair.evaluate(r, temperature);
                let mut u = pot.u_divk;
                if shift {
                    u -= pair.evaluate(r_cut, temperature).u_divk;
                }
                (u, -pot.u_divk_r)
            }
        } else {
            let pot = pair.evaluate(r, temperature);
            (pot.u_divk, -pot.u_divk_r)
        };
        if let Some(max_pot_val) = max_reduced_potential {
            let cap = max_pot_val * pair.eps_depth_divk;
            if u.is_nan() || u > cap {
                u = cap;
            }
        }
        Ok((u, force))
    }

    /// Set the segment mass `[kg]` used for quantum corrections on component `index`.
    pub fn set_mass(&mut self, index: usize, mass: f64) -> Result<()> {
        if index >= self.num_components() {
            return Err(Error::InvalidInput("component index out of bounds"));
        }
        if !mass.is_finite() || mass < 0.0 {
            return Err(Error::InvalidInput("mass must be finite and non-negative"));
        }
        self.components[index].mass = mass;
        self.recompute_cross_parameters();
        Ok(())
    }

    /// Quantum de Boer parameter `Lambda = h / (sigma * sqrt(m * eps))` for component `index`.
    pub fn de_boer_parameter(&self, index: usize) -> Result<f64> {
        let comp = self
            .components
            .get(index)
            .ok_or(Error::InvalidInput("component index out of bounds"))?;
        if comp.mass <= 0.0 {
            return Err(Error::InvalidInput(
                "component has non-positive quantum mass",
            ));
        }
        let eps = comp.eps_depth_divk * K_B;
        Ok(H_PLANCK / (comp.sigma * (comp.mass * eps).sqrt()))
    }

    /// Adjust component `index` quantum mass to match the specified de Boer parameter `Lambda`.
    pub fn set_de_boer_parameter(&mut self, index: usize, lambda: f64) -> Result<()> {
        if index >= self.num_components() {
            return Err(Error::InvalidInput("component index out of bounds"));
        }
        if !lambda.is_finite() || lambda < MACHINE_PREC {
            return Err(Error::InvalidInput(
                "de Boer parameter must be strictly positive",
            ));
        }
        let comp = &self.components[index];
        let eps = comp.eps_depth_divk * K_B;
        let sigma = comp.sigma;
        let mass = (H_PLANCK * H_PLANCK) / (sigma * sigma * lambda * lambda * eps);
        self.components[index].mass = mass;
        self.recompute_cross_parameters();
        Ok(())
    }

    /// Thermal de Broglie wavelength `h / sqrt(2 pi m k_B T) [m]` for component `index`.
    pub fn de_broglie_wavelength(&self, index: usize, temperature: f64) -> Result<f64> {
        positive(temperature, "temperature")?;
        let comp = self
            .components
            .get(index)
            .ok_or(Error::InvalidInput("component index out of bounds"))?;
        let mass = if comp.mass > 0.0 {
            comp.mass
        } else {
            1.0e-3 * comp.molar_mass / N_A
        };
        Ok(H_PLANCK / (2.0 * PI * mass * K_B * temperature).sqrt())
    }

    /// Site ownership, association energies (J/mol) and bonding volumes (m^3).
    pub fn association_parameters(&self) -> &AssociationParams {
        &self.association
    }

    /// The Dufal association integral for a pair of associating components
    /// and its derivatives with respect to T, total V and raw mole amounts.
    pub fn association_rdf_tvn(
        &self,
        i: usize,
        j: usize,
        t: f64,
        v: f64,
        n: &[f64],
    ) -> Result<crate::rdf::RdfDerivatives> {
        use crate::association::AssociationRdfProvider;
        for index in [i, j] {
            if self
                .association
                .comp_vs_sites
                .get(index)
                .and_then(|s| *s)
                .is_none()
            {
                return Err(Error::InvalidInput(
                    "RDF requires an associating component pair",
                ));
            }
        }
        positive(v, "total volume")?;
        let p = self.prepare(t, n)?;
        let sum: f64 = n.iter().sum();
        let v_unit = v / sum;
        if !v_unit.is_finite() || v_unit <= 0.0 {
            return Err(Error::InvalidInput(
                "molar volume must be finite and strictly positive",
            ));
        }
        let mut r = MieRdf::new(&p, 1.0 / v_unit).rdf_second(i, j)?.rdf;
        r.volume /= sum;
        r.temperature_volume /= sum;
        r.volume_second = (r.volume_second / sum) / sum;
        for val in r.moles.iter_mut().chain(&mut r.temperature_moles) {
            *val /= sum;
        }
        for val in r.volume_moles.iter_mut().chain(&mut r.moles_second) {
            *val = (*val / sum) / sum;
        }
        let all_finite = [
            r.g,
            r.temperature,
            r.volume,
            r.temperature_second,
            r.temperature_volume,
            r.volume_second,
        ]
        .iter()
        .all(|q| q.is_finite())
            && r.moles.iter().all(|q| q.is_finite())
            && r.temperature_moles.iter().all(|q| q.is_finite())
            && r.volume_moles.iter().all(|q| q.is_finite())
            && r.moles_second.iter().all(|q| q.is_finite());
        if !all_finite {
            return Err(Error::InvalidInput(
                "non-finite SAFT-VR Mie association RDF derivative",
            ));
        }
        Ok(r)
    }

    /// Number of components in the mixture.
    #[inline]
    pub fn num_components(&self) -> usize {
        self.components.len()
    }

    /// Access component parameters.
    pub fn components(&self) -> &[SaftVrMieComponent] {
        &self.components
    }

    #[inline(always)]
    fn potential(&self, i: usize, j: usize) -> PairPotential {
        let k = i * self.num_components() + j;
        PairPotential {
            sigma: self.sigma_ij[k],
            eps_depth_divk: self.eps_divk_ij[k],
            lambda_a: self.lambda_a_ij[k],
            lambda_r: self.lambda_r_ij[k],
            c_fac: self.c_ij[k],
            d_param: self.d_feynhibbs_param_ij[k],
            fh_order: self.components[i].fh_order,
            second_factor: self.options.quantum_variant.second_factor(),
        }
    }

    /// Repulsive-exponent binary interaction in the Mie combining rule.
    pub fn set_gammaij(&mut self, i: usize, j: usize, value: f64) {
        let nc = self.num_components();
        self.gamma_ij[i * nc + j] = value;
        self.gamma_ij[j * nc + i] = value;
        self.recompute_cross_parameters();
    }

    /// Get repulsive-exponent binary interaction parameter `gamma_ij`.
    pub fn get_gammaij(&self, i: usize, j: usize) -> f64 {
        let nc = self.num_components();
        self.gamma_ij[i * nc + j]
    }

    /// Set binary interaction parameter kij.
    pub fn set_kij(&mut self, i: usize, j: usize, val: f64) {
        let nc = self.num_components();
        self.kij[i * nc + j] = val;
        self.kij[j * nc + i] = val;
        self.recompute_cross_parameters();
    }

    /// Get binary interaction parameter kij.
    pub fn get_kij(&self, i: usize, j: usize) -> f64 {
        let nc = self.num_components();
        self.kij[i * nc + j]
    }

    /// Set binary interaction parameter lij.
    pub fn set_lij(&mut self, i: usize, j: usize, val: f64) {
        let nc = self.num_components();
        self.lij[i * nc + j] = val;
        self.lij[j * nc + i] = val;
        self.recompute_cross_parameters();
    }

    /// Get binary interaction parameter lij.
    pub fn get_lij(&self, i: usize, j: usize) -> f64 {
        let nc = self.num_components();
        self.lij[i * nc + j]
    }

    /// Recompute all cross-interaction parameters.
    fn recompute_cross_parameters(&mut self) {
        let nc = self.num_components();
        for i in 0..nc {
            for j in 0..nc {
                let idx = i * nc + j;
                let c_i = &self.components[i];
                let c_j = &self.components[j];

                // Attractive and repulsive exponents
                let la_ij = 3.0 + ((c_i.lambda_a - 3.0) * (c_j.lambda_a - 3.0)).sqrt();
                let lr_ij = 3.0
                    + (1.0 - self.gamma_ij[idx])
                        * ((c_i.lambda_r - 3.0) * (c_j.lambda_r - 3.0)).sqrt();
                self.lambda_a_ij[idx] = la_ij;
                self.lambda_r_ij[idx] = lr_ij;
                if self.association.has_sites() {
                    self.integral_coefficients[idx] = saft_vrmie_data::ASSOCIATION_B
                        .map(|row| row.into_iter().rev().fold(0.0, |acc, b| acc * lr_ij + b));
                }

                // Segment diameter
                let sig_ij = 0.5 * (c_i.sigma + c_j.sigma) * (1.0 - self.lij[idx]);
                self.sigma_ij[idx] = sig_ij;

                // Dispersion well depth (Lafitte combining rule)
                let si3 = c_i.sigma * c_i.sigma * c_i.sigma;
                let sj3 = c_j.sigma * c_j.sigma * c_j.sigma;
                let eps_ij = (1.0 - self.kij[idx])
                    * (si3 * sj3).sqrt()
                    * (c_i.eps_depth_divk * c_j.eps_depth_divk).sqrt()
                    / (sig_ij * sig_ij * sig_ij);
                self.eps_divk_ij[idx] = if self.options.geometric_epsilon {
                    (1. - self.kij[idx]) * (c_i.eps_depth_divk * c_j.eps_depth_divk).sqrt()
                } else {
                    eps_ij
                };

                // Mie prefactor C_ij
                let exp_factor = la_ij / (lr_ij - la_ij);
                let c_fac = (lr_ij / (lr_ij - la_ij)) * (lr_ij / la_ij).powf(exp_factor);
                self.c_ij[idx] = c_fac;

                // van der Waals attractive constant alpha_ij
                let alpha = c_fac * (1.0 / (la_ij - 3.0) - 1.0 / (lr_ij - 3.0));
                self.alpha_ij[idx] = alpha;

                // f_1 .. f_6 polynomials of alpha
                let f_alpha = alpha_coefficients(alpha);
                self.f_alpha_ij[idx] = f_alpha;

                // Feynman-Hibbs quantum parameter
                let mass_product = c_i.mass * c_j.mass;
                if mass_product > 0.0 {
                    let reduced_mass = mass_product / (c_i.mass + c_j.mass);
                    self.d_feynhibbs_param_ij[idx] =
                        H_PLANCK * H_PLANCK / (96.0 * reduced_mass * K_B * PI * PI);
                } else {
                    self.d_feynhibbs_param_ij[idx] = 0.0;
                }
            }
        }
    }

    /// Whether cross hard-sphere diameters use exact Gauss-Kronrod quadrature on cross Mie potential.
    pub fn exact_binary_dhs(&self) -> bool {
        self.exact_binary_dhs
    }

    /// Set whether cross hard-sphere diameters use exact Gauss-Kronrod quadrature.
    pub fn set_exact_binary_dhs(&mut self, val: bool) {
        self.exact_binary_dhs = val;
    }

    /// Prepare state for evaluation at temperature T and mole fractions x.
    pub fn prepare(
        &self,
        temperature: f64,
        mole_fractions: &[f64],
    ) -> Result<PreparedSaftVrMie<'_>> {
        if temperature <= 0.0 || !temperature.is_finite() {
            return Err(Error::InvalidInput("temperature must be strictly positive"));
        }
        if self.reference == HardSphereReference::SingleComponent && !self.exact_binary_dhs {
            return Err(Error::InvalidInput(
                "single component reference requires exact pair diameters",
            ));
        }
        let nc = self.num_components();
        if mole_fractions.len() != nc {
            return Err(Error::InvalidInput("mole fractions length mismatch"));
        }
        if mole_fractions.iter().any(|&n| n < 0.0) {
            return Err(Error::InvalidInput("nonnegative mole amounts required"));
        }
        let sum_x: f64 = mole_fractions.iter().sum();
        if sum_x <= 0.0 || !sum_x.is_finite() {
            return Err(Error::InvalidInput(
                "mole fractions sum must be strictly positive",
            ));
        }
        let x: Vec<f64> = mole_fractions.iter().map(|&v| v / sum_x).collect();

        PreparedSaftVrMie::new(self, temperature, x)
    }
}

/// State prepared at specific temperature T and mole fractions x.
#[derive(Debug, Clone)]
pub struct PreparedSaftVrMie<'a> {
    model: &'a SaftVrMie,
    temperature: f64,
    mole_fractions: Vec<f64>,
    exact_d: Vec<[f64; 3]>,
    effective: Vec<([f64; 3], [f64; 3])>,
    d: Vec<f64>,    // nc * nc
    d_t: Vec<f64>,  // nc * nc
    d_tt: Vec<f64>, // nc * nc
    m_bar: f64,
    molar_mass_mix: f64,
    conv_eta: f64,  // eta = conv_eta * rho
    conv_zeta: f64, // zeta_bar = conv_zeta * rho
    xs: Vec<f64>,   // segment mole fractions xs_i = x_i * m_i / m_bar
}

impl<'a> PreparedSaftVrMie<'a> {
    fn new(model: &'a SaftVrMie, temperature: f64, x: Vec<f64>) -> Result<Self> {
        let nc = model.num_components();
        let mut d = vec![0.0; nc * nc];
        let mut d_t = vec![0.0; nc * nc];
        let mut d_tt = vec![0.0; nc * nc];

        let extra = model.options.extra_hard_sphere;
        let need_exact_d = extra && nc > 1;
        let mut exact_d = if need_exact_d {
            vec![[0.; 3]; nc * nc]
        } else {
            Vec::new()
        };
        let mut effective = vec![([0.; 3], [0.; 3]); nc * nc];
        for i in 0..nc {
            for j in i..nc {
                let k = i * nc + j;
                let eps = model.eps_divk_ij[k];
                let la = model.lambda_a_ij[k];
                let lr = model.lambda_r_ij[k];
                let pair = model.potential(i, j);
                let fh = pair.fh_order;
                if fh != model.components[j].fh_order {
                    return Err(Error::InvalidInput("inconsistent FH orders in mixture"));
                }
                if i != j && fh > 0 && !model.options.exact_cross_potential {
                    continue;
                }
                let se = pair.effective_sigma(temperature);
                let ee = if fh == 0 {
                    [eps, 0., 0.]
                } else {
                    let potential = |r| pair.evaluate(r, temperature);
                    let mut r = se.0 * (lr / la).powf(1. / (lr - la));
                    for _ in 0..35 {
                        let u = potential(r);
                        let step = u.u_divk_r / u.u_divk_rr;
                        r -= step;
                        if step.abs() < 1e-14 * pair.sigma {
                            break;
                        }
                    }
                    let u = potential(r);
                    [
                        -u.u_divk,
                        -u.u_divk_t,
                        -u.u_divk_tt + u.u_divk_tr.powi(2) / u.u_divk_rr,
                    ]
                };
                effective[k] = ([se.0, se.1, se.2], ee);
                effective[j * nc + i] = effective[k];
            }
        }
        for i in 0..nc {
            for j in i..nc {
                let k = i * nc + j;
                if i != j
                    && model.components[i].fh_order > 0
                    && !model.options.exact_cross_potential
                {
                    for a in 0..3 {
                        effective[k].0[a] =
                            0.5 * (effective[i * nc + i].0[a] + effective[j * nc + j].0[a]);
                        effective[k].1[a] =
                            0.5 * (effective[i * nc + i].1[a] + effective[j * nc + j].1[a]);
                    }
                    effective[j * nc + i] = effective[k];
                }
                if i == j || model.exact_binary_dhs || need_exact_d {
                    let [se, st, stt] = effective[k].0;
                    let (v, vt, vtt) =
                        model.potential(i, j).diameter(temperature, (se, st, stt))?;
                    if need_exact_d {
                        exact_d[k] = [v, vt, vtt];
                        exact_d[j * nc + i] = exact_d[k];
                    }
                    if i == j || model.exact_binary_dhs {
                        d[k] = v;
                        d_t[k] = vt;
                        d_tt[k] = vtt;
                        let sym = j * nc + i;
                        d[sym] = v;
                        d_t[sym] = vt;
                        d_tt[sym] = vtt;
                    }
                }
            }
        }
        if !model.exact_binary_dhs {
            for i in 0..nc {
                let ii = i * nc + i;
                for j in i + 1..nc {
                    let jj = j * nc + j;
                    let k = i * nc + j;
                    let sym = j * nc + i;
                    d[k] = 0.5 * (d[ii] + d[jj]);
                    d_t[k] = 0.5 * (d_t[ii] + d_t[jj]);
                    d_tt[k] = 0.5 * (d_tt[ii] + d_tt[jj]);
                    d[sym] = d[k];
                    d_t[sym] = d_t[k];
                    d_tt[sym] = d_tt[k];
                }
            }
        }

        // 3. Mixture mean segment number and molar mass
        let mut m_bar = 0.0;
        let mut molar_mass_mix = 0.0;
        for (xi, comp) in x.iter().zip(&model.components) {
            m_bar += xi * comp.m;
            molar_mass_mix += xi * comp.molar_mass;
        }

        let xs: Vec<f64> = x
            .iter()
            .zip(&model.components)
            .map(|(xi, comp)| xi * comp.m / m_bar)
            .collect();

        // 4. Precompute packing fraction factors
        let mut sum_d3 = 0.0;
        let mut sum_sig3 = 0.0;
        for i in 0..nc {
            for j in 0..nc {
                let idx = i * nc + j;
                let dij = d[idx];
                let sij = effective[idx].0[0];
                match model.options.packing {
                    PackingRule::Lafitte => {
                        sum_d3 += xs[i] * xs[j] * dij.powi(3);
                        sum_sig3 += xs[i] * xs[j] * sij.powi(3);
                    }
                    PackingRule::Linear if i == j => {
                        sum_d3 += xs[i] * dij.powi(3);
                        sum_sig3 += xs[i] * sij.powi(3);
                    }
                    PackingRule::Leonard => {
                        sum_d3 += xs[i] * xs[j] * dij;
                        sum_sig3 += xs[i] * xs[j] * sij;
                    }
                    _ => {}
                }
            }
        }

        if model.options.packing == PackingRule::Leonard {
            sum_d3 = sum_d3.powi(3);
            sum_sig3 = sum_sig3.powi(3);
        }
        let conv_eta = (PI * N_A / 6.0) * m_bar * sum_d3;
        let conv_zeta = (PI * N_A / 6.0) * m_bar * sum_sig3;

        Ok(Self {
            model,
            temperature,
            mole_fractions: x,
            d,
            d_t,
            d_tt,
            exact_d,
            effective,
            m_bar,
            molar_mass_mix,
            conv_eta,
            conv_zeta,
            xs,
        })
    }

    /// Hard-sphere packing fraction eta = (pi * N_A / 6) * rho * m_bar * sum(xs_i xs_j d_ij^3).
    #[inline]
    pub fn packing_fraction(&self, rho: f64) -> f64 {
        self.conv_eta * rho
    }

    /// Dispersion packing fraction zeta_bar = (pi * N_A / 6) * rho * m_bar * sum(xs_i xs_j sigma_ij^3).
    #[inline]
    pub fn dispersion_packing_fraction(&self, rho: f64) -> f64 {
        self.conv_zeta * rho
    }

    /// Hard-sphere diameters (m).
    pub fn d(&self) -> &[f64] {
        &self.d
    }

    /// Temperature derivatives of hard-sphere diameters [m/K].
    pub fn d_t(&self) -> &[f64] {
        &self.d_t
    }

    /// Second temperature derivatives of hard-sphere diameters [m/K^2].
    pub fn d_tt(&self) -> &[f64] {
        &self.d_tt
    }

    /// Helmholtz free energy terms at density rho [mol/m^3].
    pub fn helmholtz_terms(&self, rho: f64) -> Result<SaftVrMieHelmholtzTerms> {
        if rho == 0.0 {
            return Ok(SaftVrMieHelmholtzTerms {
                a_hs: 0.,
                a1: 0.,
                a2: 0.,
                a3: 0.,
                a_chain: 0.,
                a_res: 0.,
            });
        }
        let terms = self.eval_scalar(rho)?;
        let t = self.temperature;
        let a_assoc = if self.model.association.has_sites() && rho > 0.0 {
            assoc_helmholtz_value(
                &self.model.association,
                &MieRdf::new(self, rho),
                t,
                1.0 / rho,
                &self.mole_fractions,
            )?
            .0
        } else {
            0.0
        };
        Ok(SaftVrMieHelmholtzTerms {
            a_hs: terms.hs,
            a1: terms.a1,
            a2: terms.a2,
            a3: terms.a3,
            a_chain: terms.chain,
            a_res: terms.total + a_assoc,
        })
    }

    /// Pressure and derivatives at density rho [mol/m^3].
    pub fn pressure_and_derivatives(&self, rho: f64) -> Result<SaftVrMiePressureJet> {
        let v = 1.0 / rho;
        let t = self.temperature;

        let vv = self.eval_second(rho, 1, 1, EvalMode::Total)?.total;
        let tv = self.eval_second(rho, 0, 1, EvalMode::Total)?.total;
        let mut f = vv.f0;
        let mut f_v = vv.f1;
        let mut f_vv = vv.f12;
        let mut f_tv = tv.f12;
        if self.model.association.has_sites() {
            let a = assoc_helmholtz_pressure_jet(
                &self.model.association,
                &MieRdf::new(self, rho),
                t,
                v,
                &self.mole_fractions,
            )?;
            f += a.f;
            f_v += a.f_v;
            f_vv += a.f_vv;
            f_tv += a.f_tv;
        }
        let p = RGAS * t * (rho - f_v);
        let z = 1. - f_v / rho;
        let dp_dv = -RGAS * t * (rho * rho + f_vv);
        let dp_drho = -dp_dv / (rho * rho);

        // dP/dT = rho * R - R * F_V - R * T * F_TV
        let dp_dt = rho * RGAS - RGAS * f_v - RGAS * t * f_tv;

        Ok(SaftVrMiePressureJet {
            pressure: p,
            compressibility: z,
            residual_helmholtz: f,
            dp_drho,
            dp_dv,
            dp_dt,
        })
    }

    /// Complete thermodynamic state properties at density rho [mol/m^3].
    pub fn state(&self, rho: f64) -> Result<SaftVrMieState> {
        let t = self.temperature;
        let v = 1.0 / rho;
        let nc = self.model.components.len();

        let vv = self.eval_second(rho, 1, 1, EvalMode::Total)?.total;
        let tv = self.eval_second(rho, 0, 1, EvalMode::Total)?.total;
        let tt = self.eval_second(rho, 0, 0, EvalMode::Total)?.total;
        let mut f = vv.f0;
        let mut f_v = vv.f1;
        let mut f_vv = vv.f12;
        let mut f_tv = tv.f12;
        let mut f_t = tt.f1;
        let mut f_tt = tt.f12;

        let mut assoc_fn: Option<Vec<f64>> = None;
        if self.model.association.has_sites() {
            let rdf = MieRdf::new(self, rho);
            if nc == 1 {
                let a = assoc_helmholtz_pressure_jet(
                    &self.model.association,
                    &rdf,
                    t,
                    v,
                    &self.mole_fractions,
                )?;
                f += a.f;
                f_t += a.f_t;
                f_v += a.f_v;
                f_tt += a.f_tt;
                f_tv += a.f_tv;
                f_vv += a.f_vv;
            } else {
                let a =
                    assoc_helmholtz_all(&self.model.association, &rdf, t, v, &self.mole_fractions)?;
                f += a.f;
                f_t += a.f_t;
                f_v += a.f_v;
                f_tt += a.f_tt;
                f_tv += a.f_tv;
                f_vv += a.f_vv;
                assoc_fn = Some(a.f_n);
            }
        }

        let p = RGAS * t * (rho - f_v);
        let z = 1.0 - f_v / rho;
        let dp_dv = -RGAS * t * (rho * rho + f_vv);
        let dp_drho = -dp_dv / (rho * rho);
        let dp_dt = rho * RGAS - RGAS * f_v - RGAS * t * f_tv;

        let u_res = -RGAS * t * t * f_t;
        let h_res = u_res + RGAS * t * (z - 1.0);
        let s_res_tv = -RGAS * (t * f_t + f);
        let s_res = s_res_tv + RGAS * z.ln();
        let a_res_energy = RGAS * t * f;
        let g_res = h_res - t * s_res;

        let cv_res = -RGAS * (2.0 * t * f_t + t * t * f_tt);

        // Ideal gas contribution evaluated at TV ideal pressure rho * R * T
        let mut cp_id = 0.0;
        let mut h_id = 0.0;
        let mut s_id_tv = 0.0;
        let p_ref = 100_000.0;
        for i in 0..nc {
            let id_props = self.model.ideal[i].properties(t, p_ref)?;
            let xi = self.mole_fractions[i];
            cp_id += xi * id_props.cp;
            h_id += xi * id_props.enthalpy;
            let mix_term = if xi > 0.0 { xi * xi.ln() } else { 0.0 };
            s_id_tv += xi * id_props.entropy - RGAS * mix_term;
        }
        s_id_tv -= RGAS * (rho * RGAS * t / p_ref).ln();

        let cv_total = cp_id - RGAS + cv_res;
        let cp_total = cv_total - t * (dp_dt * dp_dt) / dp_dv;

        let gamma = cp_total / cv_total;
        let mw_kg = self.molar_mass_mix * 1e-3;
        let sound_speed = if dp_drho > 0.0 && gamma > 0.0 && mw_kg > 0.0 {
            (gamma * dp_drho / mw_kg).sqrt()
        } else {
            0.0
        };

        // Fugacity coefficients
        let ln_z = z.ln();
        let ln_phi = if nc == 1 {
            vec![f + (z - 1.0) - ln_z]
        } else {
            let mut vals = (0..nc)
                .map(|i| Ok(self.eval_first(rho, 2 + i)?.total.f1 - ln_z))
                .collect::<Result<Vec<_>>>()?;
            if let Some(afn) = assoc_fn {
                for (value, assoc) in vals.iter_mut().zip(afn) {
                    *value += assoc;
                }
            }
            vals
        };

        Ok(SaftVrMieState {
            pressure: p,
            compressibility: z,
            residual_helmholtz: f,
            u_residual: u_res,
            h_residual: h_res,
            s_residual: s_res,
            a_residual: a_res_energy,
            g_residual: g_res,
            h_total: h_id + h_res,
            s_total: s_id_tv + s_res_tv,
            cv_residual: cv_res,
            cv_total,
            cp_total,
            dp_drho,
            dp_dv,
            dp_dt,
            sound_speed,
            ln_phi,
        })
    }

    /// Fugacity coefficients ln(phi_i) at density rho [mol/m^3].
    pub fn fugacity_coefficients(&self, rho: f64, z: f64) -> Result<Vec<f64>> {
        let nc = self.model.components.len();
        if nc == 1 {
            let terms = self.helmholtz_terms(rho)?;
            let ln_phi = terms.a_res + (z - 1.0) - z.ln();
            return Ok(vec![ln_phi]);
        }

        let ln_z = z.ln();
        let mut ln_phi = (0..nc)
            .map(|i| Ok(self.eval_first(rho, 2 + i)?.total.f1 - ln_z))
            .collect::<Result<Vec<_>>>()?;
        if self.model.association.has_sites() {
            let a = assoc_helmholtz_first(
                &self.model.association,
                &MieRdf::new(self, rho),
                self.temperature,
                1. / rho,
                &self.mole_fractions,
            )?;
            for (value, assoc) in ln_phi.iter_mut().zip(a.f_n) {
                *value += assoc;
            }
        }
        Ok(ln_phi)
    }

    /// Dimensionless Gibbs energy difference `(G(vl) - G(vv)) / (R*T)` at a
    /// common imposed pressure `pressure`:
    /// `(F(T, vl) - F(T, vv)) + P*(vl - vv)/(R*T) - ln(vl / vv)`.
    ///
    /// Because each root satisfies `d[F + P*v/(R*T) - ln v]/dv = 0` at `P(v) = P`,
    /// small pressure residuals on stiff dense liquid roots enter only at second
    /// order rather than leaking first-order into `ln(Z)`.
    fn imposed_pressure_gibbs_diff(&self, pressure: f64, vl: f64, vv: f64) -> Result<f64> {
        let fl = self.helmholtz_terms(1.0 / vl)?.a_res;
        let fv = self.helmholtz_terms(1.0 / vv)?.a_res;
        Ok((fl - fv) + pressure * (vl - vv) / (RGAS * self.temperature) - (vl / vv).ln())
    }

    /// Single-phase volume root for a specified pressure and phase.
    pub fn volume(&self, pressure_target: f64, phase: Phase) -> Result<f64> {
        if pressure_target <= 0.0 || !pressure_target.is_finite() {
            return Err(Error::InvalidInput(
                "pressure target must be strictly positive",
            ));
        }

        let conv_eta = self.conv_eta;
        if phase == Phase::MinimumGibbs {
            let vl = self.volume(pressure_target, Phase::Liquid)?;
            let vv = self.volume(pressure_target, Phase::Vapor)?;
            if (vl - vv).abs() < 1e-10 * vl {
                return Ok(vl);
            }
            let dg = self.imposed_pressure_gibbs_diff(pressure_target, vl, vv)?;
            return Ok(if dg < 0.0 { vl } else { vv });
        }

        let zeta_ideal = conv_eta * pressure_target / (RGAS * self.temperature);
        if !zeta_ideal.is_finite() || zeta_ideal <= 0.0 {
            return Err(Error::NotConverged(
                "SAFT-VR Mie volume target underflowed representable packing range",
            ));
        }
        let mut zeta_min = (1.0e-10_f64).min(0.1 * zeta_ideal).max(1.0e-140);
        let mut zeta_max = 0.75;

        let mut zeta = match phase {
            Phase::Vapor => zeta_ideal.clamp(zeta_min, 0.05),
            Phase::Liquid => 0.45,
            Phase::MinimumGibbs => unreachable!(),
        };

        let mut v = conv_eta / zeta;
        if !v.is_finite() || v <= 0.0 || !(1.0 / v).is_finite() || (1.0 / v) <= 0.0 {
            return Err(Error::NotConverged(
                "SAFT-VR Mie volume target exceeded representable molar volume",
            ));
        }
        let mut jet = self.pressure_and_derivatives(1.0 / v)?;
        let mut f = (1.0 - zeta) * (jet.pressure - pressure_target);

        let max_iter = 80;
        for _ in 0..max_iter {
            let df_dzeta = pressure_target - jet.pressure - ((1.0 - zeta) / zeta) * jet.dp_dv * v;
            let mut zeta_new = if df_dzeta.is_finite() && df_dzeta != 0.0 {
                zeta - f / df_dzeta
            } else {
                0.5 * (zeta_max + zeta_min)
            };

            if zeta_new > zeta_max || zeta_new < zeta_min || !zeta_new.is_finite() {
                zeta_new = 0.5 * (zeta_max + zeta_min);
            }

            let v_old = v;
            zeta = zeta_new;
            v = conv_eta / zeta;
            if !v.is_finite() || v <= 0.0 || !(1.0 / v).is_finite() || (1.0 / v) <= 0.0 {
                return Err(Error::NotConverged(
                    "SAFT-VR Mie volume solver left representable molar volume range",
                ));
            }

            jet = self.pressure_and_derivatives(1.0 / v)?;
            f = (1.0 - zeta) * (jet.pressure - pressure_target);

            let v_converged = (v - v_old).abs() <= 1e2 * MACHINE_PREC * v_old;
            let p_tol = (1e5 * MACHINE_PREC * pressure_target)
                .max(1e2 * MACHINE_PREC * (jet.dp_dv * v).abs())
                .min(1.0e-6 * pressure_target);
            let p_converged =
                jet.pressure.is_finite() && (jet.pressure - pressure_target).abs() <= p_tol;

            if v_converged && p_converged {
                return Ok(v);
            }

            if f > 0.0 {
                zeta_max = zeta;
            } else if f < 0.0 {
                zeta_min = zeta;
            }
        }

        if v.is_finite()
            && v > 0.0
            && jet.pressure.is_finite()
            && (jet.pressure - pressure_target).abs() <= 1e-6 * pressure_target
        {
            Ok(v)
        } else {
            Err(Error::NotConverged(
                "SAFT-VR Mie volume solver did not converge",
            ))
        }
    }

    /// Solve `P(conv_eta / zeta) = pressure_target` strictly within a monotone
    /// branch bracket `[zeta_lo, zeta_hi]`.
    fn volume_in_zeta_bracket(
        &self,
        pressure_target: f64,
        mut zeta_min: f64,
        mut zeta_max: f64,
        zeta_init: f64,
    ) -> Result<f64> {
        if !pressure_target.is_finite()
            || pressure_target <= 0.0
            || !zeta_min.is_finite()
            || !zeta_max.is_finite()
            || zeta_min <= 0.0
            || zeta_min >= zeta_max
            || zeta_max >= 1.0
        {
            return Err(Error::NotConverged(
                "SAFT-VR Mie volume solver received invalid packing bracket",
            ));
        }
        let conv_eta = self.conv_eta;
        let mut zeta = zeta_init.clamp(zeta_min, zeta_max);
        let mut v = conv_eta / zeta;
        if !v.is_finite() || v <= 0.0 || !(1.0 / v).is_finite() || (1.0 / v) <= 0.0 {
            return Err(Error::NotConverged(
                "SAFT-VR Mie volume target exceeded representable molar volume",
            ));
        }
        let mut jet = self.pressure_and_derivatives(1.0 / v)?;
        let mut f = (1.0 - zeta) * (jet.pressure - pressure_target);

        let max_iter = 80;
        for _ in 0..max_iter {
            let df_dzeta = pressure_target - jet.pressure - ((1.0 - zeta) / zeta) * jet.dp_dv * v;
            let mut zeta_new = if df_dzeta.is_finite() && df_dzeta != 0.0 {
                zeta - f / df_dzeta
            } else {
                0.5 * (zeta_max + zeta_min)
            };

            if zeta_new > zeta_max || zeta_new < zeta_min || !zeta_new.is_finite() {
                zeta_new = 0.5 * (zeta_max + zeta_min);
            }

            let v_old = v;
            zeta = zeta_new;
            v = conv_eta / zeta;
            if !v.is_finite() || v <= 0.0 || !(1.0 / v).is_finite() || (1.0 / v) <= 0.0 {
                return Err(Error::NotConverged(
                    "SAFT-VR Mie volume solver left representable molar volume range",
                ));
            }

            jet = self.pressure_and_derivatives(1.0 / v)?;
            f = (1.0 - zeta) * (jet.pressure - pressure_target);

            let v_converged = (v - v_old).abs() <= 1e2 * MACHINE_PREC * v_old;
            let p_tol = (1e5 * MACHINE_PREC * pressure_target)
                .max(1e2 * MACHINE_PREC * (jet.dp_dv * v).abs())
                .min(1.0e-6 * pressure_target);
            let p_converged =
                jet.pressure.is_finite() && (jet.pressure - pressure_target).abs() <= p_tol;

            if v_converged && p_converged {
                return Ok(v);
            }

            if f > 0.0 {
                zeta_max = zeta;
            } else if f < 0.0 {
                zeta_min = zeta;
            }
        }

        if v.is_finite()
            && v > 0.0
            && jet.pressure.is_finite()
            && (jet.pressure - pressure_target).abs() <= 1e-6 * pressure_target
        {
            Ok(v)
        } else {
            Err(Error::NotConverged(
                "SAFT-VR Mie volume solver did not converge",
            ))
        }
    }

    /// Locate the mechanical spinodals on this pure-fluid isotherm by searching
    /// for `dP/drho < 0` in the physical fluid packing range (`eta <= 0.56`,
    /// staying below unphysical high-packing association/dispersion polynomial poles).
    ///
    /// Returns `Ok(None)` when `dP/drho >= 0` everywhere (supercritical state,
    /// where vapor-liquid coexistence does not exist).
    fn locate_spinodals(&self) -> Result<Option<IsothermSpinodals>> {
        let conv_eta = self.conv_eta;
        let dp_drho_at = |eta: f64| -> Result<f64> {
            Ok(self.pressure_and_derivatives(eta / conv_eta)?.dp_drho)
        };

        // Scan across the physical fluid packing range to spot or bracket the
        // fluid-phase minimum of dP/drho without stepping into unphysical high-eta poles.
        const N_GRID: usize = 40;
        let eta_lo = 0.005_f64;
        let eta_hi = 0.56_f64;
        let mut grid_eta = [0.0_f64; N_GRID];
        let mut grid_dp = [0.0_f64; N_GRID];
        let mut min_idx = 0usize;
        let mut min_dp = f64::INFINITY;
        for k in 0..N_GRID {
            let eta = eta_lo + (eta_hi - eta_lo) * (k as f64) / ((N_GRID - 1) as f64);
            let dp = dp_drho_at(eta)?;
            grid_eta[k] = eta;
            grid_dp[k] = dp;
            if dp < min_dp {
                min_dp = dp;
                min_idx = k;
            }
        }

        let mut eta_mid = if min_dp < 0.0 {
            Some(grid_eta[min_idx])
        } else {
            None
        };

        // If no coarse point was negative, minimize dP/drho around min_idx via
        // golden-section search to catch narrow near-critical unstable loops.
        if eta_mid.is_none() {
            let mut a = if min_idx > 0 {
                grid_eta[min_idx - 1]
            } else {
                1.0e-4
            };
            let mut b = if min_idx + 1 < N_GRID {
                grid_eta[min_idx + 1]
            } else {
                0.58
            };
            let inv_phi = 0.5 * (5.0_f64.sqrt() - 1.0);
            let mut c = b - inv_phi * (b - a);
            let mut d = a + inv_phi * (b - a);
            let mut fc = dp_drho_at(c)?;
            let mut fd = dp_drho_at(d)?;
            for _ in 0..55 {
                if fc < 0.0 {
                    eta_mid = Some(c);
                    break;
                }
                if fd < 0.0 {
                    eta_mid = Some(d);
                    break;
                }
                if (b - a) <= 1.0e-10 * b {
                    break;
                }
                if fc < fd {
                    b = d;
                    d = c;
                    fd = fc;
                    c = b - inv_phi * (b - a);
                    fc = dp_drho_at(c)?;
                } else {
                    a = c;
                    c = d;
                    fc = fd;
                    d = a + inv_phi * (b - a);
                    fd = dp_drho_at(d)?;
                }
            }
        }

        let Some(eta_unstable) = eta_mid else {
            return Ok(None);
        };

        // Vapor spinodal (local P maximum): dP/drho crosses + -> - on (0, eta_unstable).
        let mut lo_v = 1.0e-8_f64.min(0.1 * eta_unstable);
        while dp_drho_at(lo_v)? <= 0.0 && lo_v > 1.0e-30 {
            lo_v *= 0.1;
        }
        if dp_drho_at(lo_v)? <= 0.0 {
            return Ok(None);
        }
        let mut hi_v = eta_unstable;
        for _ in 0..64 {
            if (hi_v - lo_v) <= 1.0e-12 * hi_v {
                break;
            }
            let mid = 0.5 * (lo_v + hi_v);
            if dp_drho_at(mid)? > 0.0 {
                lo_v = mid;
            } else {
                hi_v = mid;
            }
        }
        let eta_spin_v = lo_v;
        let p_max = self
            .pressure_and_derivatives(eta_spin_v / conv_eta)?
            .pressure;

        // Liquid spinodal (local P minimum): find the first positive-dP/drho grid point
        // strictly above eta_unstable, bounded strictly by ETA_LIQUID_CAP so the liquid
        // bracket never inverts or steps into unphysical high-packing poles.
        const ETA_LIQUID_CAP: f64 = 0.72;
        let mut lo_l = eta_unstable;
        let mut hi_l = None;
        for k in 0..N_GRID {
            if grid_eta[k] > eta_unstable && grid_dp[k] > 0.0 {
                hi_l = Some(grid_eta[k]);
                break;
            }
        }
        let mut hi_l = match hi_l {
            Some(eta) => eta,
            None => {
                if eta_unstable >= ETA_LIQUID_CAP {
                    return Ok(None);
                }
                let mut probe = (eta_unstable + 0.02).min(ETA_LIQUID_CAP);
                while probe < ETA_LIQUID_CAP && dp_drho_at(probe)? <= 0.0 {
                    probe = (probe + 0.02).min(ETA_LIQUID_CAP);
                }
                if dp_drho_at(probe)? <= 0.0 {
                    return Ok(None);
                }
                probe
            }
        };
        for _ in 0..64 {
            if (hi_l - lo_l) <= 1.0e-12 * hi_l {
                break;
            }
            let mid = 0.5 * (lo_l + hi_l);
            if dp_drho_at(mid)? < 0.0 {
                lo_l = mid;
            } else {
                hi_l = mid;
            }
        }
        let eta_spin_l = hi_l;
        if eta_spin_l >= ETA_LIQUID_CAP {
            return Ok(None);
        }
        let p_min = self
            .pressure_and_derivatives(eta_spin_l / conv_eta)?
            .pressure;

        // Choose a safe upper liquid-branch packing bound strictly above eta_spin_l
        // and at or below ETA_LIQUID_CAP where dP/drho > 0.
        let mut eta_l_max = eta_spin_l;
        for k in 0..N_GRID {
            if grid_eta[k] > eta_spin_l {
                if grid_dp[k] <= 0.0 {
                    break;
                }
                eta_l_max = grid_eta[k];
            }
        }
        if eta_l_max <= eta_spin_l {
            eta_l_max = (eta_spin_l + 0.05).min(ETA_LIQUID_CAP);
        }

        if !eta_spin_v.is_finite()
            || !eta_spin_l.is_finite()
            || !eta_l_max.is_finite()
            || eta_spin_v <= 0.0
            || eta_spin_v >= eta_spin_l
            || eta_spin_l >= eta_l_max
            || eta_l_max > ETA_LIQUID_CAP
            || !p_max.is_finite()
            || !p_min.is_finite()
            || p_max <= 0.0
            || p_max <= p_min
        {
            return Ok(None);
        }
        Ok(Some(IsothermSpinodals {
            eta_v: eta_spin_v,
            p_max,
            eta_l: eta_spin_l,
            p_min,
            eta_l_max,
        }))
    }

    /// Pure-fluid saturation pressure at temperature T.
    pub fn saturation_pressure(&self, temperature: f64) -> Result<SaftVrMieSaturation> {
        if self.model.num_components() != 1 {
            return Err(Error::Unsupported(
                "saturation_pressure is only supported for pure fluids",
            ));
        }

        let owned_prep;
        let prep = if (self.temperature - temperature).abs() == 0.0 {
            self
        } else {
            owned_prep = self.model.prepare(temperature, &[1.0])?;
            &owned_prep
        };
        let comp_info = component(self.model.components[0].ident)?;
        let tc = comp_info.critical_temperature;
        let pc = comp_info.critical_pressure;
        let omega = comp_info.acentric_factor;

        let exp_arg: f64 = 5.373 * (1.0 + omega) * (1.0 - tc / temperature);
        let mut p = pc * exp_arg.exp();
        if p >= pc || p <= 0.0 || !p.is_finite() {
            p = 0.5 * pc;
        }

        let mut p_lo: Option<f64> = None;
        let mut p_hi: Option<f64> = None;
        let mut spinodals: Option<IsothermSpinodals> = None;
        let max_iter = 60;

        for _ in 0..max_iter {
            let (vl, vv) = if let Some(sp) = spinodals {
                let zeta_ideal = prep.conv_eta * p / (RGAS * temperature);
                let zeta_min_v = (1.0e-10_f64)
                    .min(0.1 * zeta_ideal)
                    .max(1.0e-140)
                    .min(0.5 * sp.eta_v);
                let vl = prep.volume_in_zeta_bracket(
                    p,
                    sp.eta_l,
                    sp.eta_l_max,
                    0.5 * (sp.eta_l + sp.eta_l_max),
                )?;
                let zeta_v_init = if zeta_min_v < 0.9 * sp.eta_v {
                    zeta_ideal.clamp(zeta_min_v, 0.9 * sp.eta_v)
                } else {
                    0.5 * (zeta_min_v + sp.eta_v)
                };
                let vv = prep.volume_in_zeta_bracket(p, zeta_min_v, sp.eta_v, zeta_v_init)?;
                (vl, vv)
            } else {
                let vl_res = prep.volume(p, Phase::Liquid);
                let vv_res = prep.volume(p, Phase::Vapor);
                match (vl_res, vv_res) {
                    (Ok(vl), Ok(vv)) if vv > vl && (vv - vl) > 1e-10 * vl => (vl, vv),
                    _ => {
                        // Switch to exact mechanical spinodal bounds so near-critical
                        // coexistence is preserved and supercritical states are rejected.
                        let Some(sp) = prep.locate_spinodals()? else {
                            return Err(Error::NotConverged(
                                "supercritical temperature: no vapor-liquid coexistence",
                            ));
                        };
                        spinodals = Some(sp);
                        if sp.p_min > 0.0 {
                            p_lo = Some(p_lo.map_or(sp.p_min, |lo| lo.max(sp.p_min)));
                        }
                        p_hi = Some(p_hi.map_or(sp.p_max, |hi| hi.min(sp.p_max)));
                        p = if sp.p_min > 0.0 {
                            0.5 * (sp.p_min + sp.p_max)
                        } else if p <= 0.0 || p >= sp.p_max {
                            0.5 * sp.p_max
                        } else {
                            p
                        };
                        continue;
                    }
                }
            };

            let dg = prep.imposed_pressure_gibbs_diff(p, vl, vv)?;
            if dg > 0.0 {
                p_lo = Some(p_lo.map_or(p, |lo| lo.max(p)));
            } else if dg < 0.0 {
                p_hi = Some(p_hi.map_or(p, |hi| hi.min(p)));
            }

            let d_dg_dlnp = p * (vl - vv) / (RGAS * temperature);
            let d_ln_p = -dg / d_dg_dlnp;
            let bracket_tight =
                matches!((p_lo, p_hi), (Some(lo), Some(hi)) if (hi - lo) <= 1e-13 * hi);
            if (dg.abs() < 1e-13 && d_ln_p.abs() < 1e-12) || (bracket_tight && dg.abs() < 1e-10) {
                let state_l = prep.state(1.0 / vl)?;
                let state_v = prep.state(1.0 / vv)?;
                let p_closed = state_l.pressure.is_finite()
                    && state_v.pressure.is_finite()
                    && (state_l.pressure - p).abs() <= 1.0e-6 * p
                    && (state_v.pressure - p).abs() <= 1.0e-6 * p;
                if vl < vv && state_l.dp_drho > 0.0 && state_v.dp_drho > 0.0 && p_closed {
                    return Ok(SaftVrMieSaturation {
                        pressure: p,
                        temperature,
                        liquid_volume: vl,
                        vapor_volume: vv,
                        enthalpy_vaporization: state_v.h_total - state_l.h_total,
                        entropy_vaporization: state_v.s_total - state_l.s_total,
                    });
                }
                return Err(Error::NotConverged(
                    "saturation roots failed mechanical stability or pressure closure check",
                ));
            }

            let mut p_next = if d_ln_p.is_finite() {
                p * d_ln_p.clamp(-2.0, 2.0).exp()
            } else {
                f64::NAN
            };
            if !p_next.is_finite()
                || p_next <= 0.0
                || p_lo.is_some_and(|lo| p_next <= lo)
                || p_hi.is_some_and(|hi| p_next >= hi)
            {
                p_next = match (p_lo, p_hi) {
                    (Some(lo), Some(hi)) if hi <= 10.0 * lo => 0.5 * (lo + hi),
                    (Some(lo), Some(hi)) => (lo * hi).sqrt(),
                    _ => p * d_ln_p.clamp(-1.0, 1.0).exp(),
                };
            }
            p = p_next;
        }

        Err(Error::NotConverged(
            "saturation solver exceeded maximum iterations",
        ))
    }
}

// ----------------------------------------------------------------------------
// Internal numerical & potential helper routines
// ----------------------------------------------------------------------------

#[derive(Debug, Clone, Copy)]
struct IsothermSpinodals {
    eta_v: f64,
    p_max: f64,
    eta_l: f64,
    p_min: f64,
    eta_l_max: f64,
}

#[derive(Debug, Clone, Copy)]
struct MiePot {
    u_divk: f64,
    u_divk_t: f64,
    u_divk_tt: f64,
    u_divk_r: f64,
    u_divk_rr: f64,
    u_divk_tr: f64,
}

#[derive(Clone, Copy)]
struct PairPotential {
    sigma: f64,
    eps_depth_divk: f64,
    lambda_a: f64,
    lambda_r: f64,
    c_fac: f64,
    d_param: f64,
    fh_order: i32,
    second_factor: f64,
}
impl PairPotential {
    #[inline(always)]
    fn evaluate_classical(self, r: f64) -> (f64, f64, f64) {
        let n = self.lambda_r;
        let m = self.lambda_a;
        let en_r = 1.0 / r;
        let r1 = self.sigma * en_r;
        let (r_n, r_m) = pow_pos_pair(r1, n, m);
        let dr_n = -n * r_n * en_r;
        let dr_m = -m * r_m * en_r;
        let d2r_n = -(n + 1.0) * dr_n * en_r;
        let d2r_m = -(m + 1.0) * dr_m * en_r;
        let pref = self.c_fac * self.eps_depth_divk;
        (
            (r_n - r_m) * pref,
            (dr_n - dr_m) * pref,
            (d2r_n - d2r_m) * pref,
        )
    }

    #[inline(always)]
    fn evaluate(self, r: f64, temperature: f64) -> MiePot {
        let Self {
            sigma,
            eps_depth_divk,
            lambda_a,
            lambda_r,
            c_fac,
            d_param,
            fh_order,
            second_factor,
        } = self;
        let n = lambda_r;
        let m = lambda_a;
        let en_r = 1.0 / r;
        let r1 = sigma * en_r;

        let (r_n, r_m) = pow_pos_pair(r1, n, m);
        let dr_n = -n * r_n * en_r;
        let dr_m = -m * r_m * en_r;
        let d2r_n = -(n + 1.0) * dr_n * en_r;
        let d2r_m = -(m + 1.0) * dr_m * en_r;

        let mut u = r_n - r_m;
        let mut u_t = 0.0;
        let mut u_tt = 0.0;
        let mut u_r = dr_n - dr_m;
        let mut u_rr = d2r_n - d2r_m;
        let mut u_tr = 0.0;

        if fh_order > 0 && d_param > 0.0 && temperature > 0.0 {
            let tinv = 1.0 / temperature;
            let d = d_param * tinv;
            let d_t = -d_param * tinv * tinv;
            let d_tt = 2.0 * d_param * (tinv * tinv * tinv);

            let q1_n = n * (n - 1.0);
            let q1_m = m * (m - 1.0);
            let en_r2 = en_r * en_r;

            let prod_q1 = en_r2 * (q1_n * r_n - q1_m * r_m);
            let prod_q1_r = -en_r2 * ((n + 2.0) * q1_n * r_n - (m + 2.0) * q1_m * r_m) * en_r;
            let prod_q1_rr = en_r2
                * ((n + 3.0) * (n + 2.0) * q1_n * r_n - (m + 3.0) * (m + 2.0) * q1_m * r_m)
                * en_r2;

            u += d * prod_q1;
            u_t += d_t * prod_q1;
            u_tt += d_tt * prod_q1;
            u_r += d * prod_q1_r;
            u_rr += d * prod_q1_rr;
            u_tr += d_t * prod_q1_r;

            if fh_order > 1 {
                let d2 = d * d;
                let d2_t = 2.0 * d * d_t;
                let d2_tt = 2.0 * (d * d_tt + d_t * d_t);

                let q2_n = second_factor * n * (n - 1.0) * (n + 1.0) * (n + 2.0);
                let q2_m = second_factor * m * (m - 1.0) * (m + 1.0) * (m + 2.0);
                let en_r4 = en_r2 * en_r2;

                let prod_q2 = en_r4 * (q2_n * r_n - q2_m * r_m);
                let prod_q2_r = -en_r4 * ((n + 4.0) * q2_n * r_n - (m + 4.0) * q2_m * r_m) * en_r;
                let prod_q2_rr = en_r4
                    * ((n + 5.0) * (n + 4.0) * q2_n * r_n - (m + 5.0) * (m + 4.0) * q2_m * r_m)
                    * en_r2;

                u += d2 * prod_q2;
                u_t += d2_t * prod_q2;
                u_tt += d2_tt * prod_q2;
                u_r += d2 * prod_q2_r;
                u_rr += d2 * prod_q2_rr;
                u_tr += d2_t * prod_q2_r;
            }
        }

        let pref = c_fac * eps_depth_divk;
        MiePot {
            u_divk: u * pref,
            u_divk_t: u_t * pref,
            u_divk_tt: u_tt * pref,
            u_divk_r: u_r * pref,
            u_divk_rr: u_rr * pref,
            u_divk_tr: u_tr * pref,
        }
    }

    #[inline]
    fn effective_sigma(self, temperature: f64) -> (f64, f64, f64) {
        let sigma = self.sigma;
        let fh_order = self.fh_order;
        // 1. Solve effective sigma (where U(r) = 0)
        if fh_order == 0 {
            (sigma, 0.0, 0.0)
        } else {
            let mut r = sigma;
            for _ in 0..25 {
                let pot = self.evaluate(r, temperature);
                let f = pot.u_divk;
                let df = pot.u_divk_r;
                let d2f = pot.u_divk_rr;
                let step = -f / df * (1.0 + f * d2f / (2.0 * df * df));
                r += step;
                if step.abs() < 1e-15 * sigma {
                    break;
                }
            }
            let pot = self.evaluate(r, temperature);
            let s_t = -pot.u_divk_t / pot.u_divk_r;
            let s_tt = -(pot.u_divk_tt + pot.u_divk_rr * s_t * s_t + 2.0 * pot.u_divk_tr * s_t)
                / pot.u_divk_r;
            (r, s_t, s_tt)
        }
    }

    // Keep diameter self-contained so LLVM always optimizes its inner 21-point
    // Gauss-Kronrod loop independently of unrelated callers in this module.
    #[inline(never)]
    fn diameter(
        self,
        temperature: f64,
        (sigma_eff, sig_t, sig_tt): (f64, f64, f64),
    ) -> Result<(f64, f64, f64)> {
        let ln_eps = MACHINE_PREC.ln();
        let inv_t = 1.0 / temperature;
        let inv_t2 = inv_t * inv_t;
        let inv_t3 = inv_t2 * inv_t;

        if self.fh_order == 0 {
            // Classical Mie fast path: u_t = u_tt = u_tr = 0 and sig_t = sig_tt = 0.
            let mut rs = 0.70;
            for _ in 0..25 {
                let r = rs * sigma_eff;
                let (u, u_r, u_rr) = self.evaluate_classical(r);
                let f = u * inv_t + ln_eps;
                let df = sigma_eff * u_r * inv_t;
                let d2f = sigma_eff * sigma_eff * u_rr * inv_t;
                let step = -f / df * (1.0 + f * d2f / (2.0 * df * df));
                rs += step;
                if step.abs() < 1e-15 {
                    break;
                }
            }

            let r0 = rs * sigma_eff;
            let (_, u0_r, u0_rr) = self.evaluate_classical(r0);
            let r0_t = -ln_eps / u0_r;
            let r0_tt = -(u0_rr * r0_t * r0_t) / u0_r;

            let mut d_int = 0.0;
            let mut d_t_int = 0.0;
            let mut d_tt_int = 0.0;

            let sigma_j2 = 0.5 * (sigma_eff - r0);
            let sigma_j2_t = -0.5 * r0_t;
            let sigma_j2_tt = -0.5 * r0_tt;

            for i in 0..21 {
                let xk = X_NODES_21[i];
                let wk = W_NODES_21[i];
                let one_minus_xk = 1.0 - xk;

                let r = sigma_j2 * xk + sigma_j2 + r0;
                let r_t = 0.5 * one_minus_xk * r0_t;
                let r_tt = 0.5 * one_minus_xk * r0_tt;

                let (u, u_r, u_rr) = self.evaluate_classical(r);

                let pref = sigma_j2 * wk;
                let pref_t = sigma_j2_t * wk;
                let pref_tt = sigma_j2_tt * wk;

                let minus_exp = -(-u * inv_t).exp();
                let one_plus_mexp = 1.0 + minus_exp;
                d_int += pref * one_plus_mexp;

                let dt3 = u * inv_t2 - u_r * r_t * inv_t;
                d_t_int += pref * minus_exp * dt3 + pref_t * one_plus_mexp;

                let d_tt2 = minus_exp * dt3;
                let d_tt3 = 2.0 * u_r * r_t * inv_t2
                    - 2.0 * u * inv_t3
                    - u_rr * (r_t * r_t) * inv_t
                    - u_r * r_tt * inv_t;
                d_tt_int += pref_t * minus_exp * dt3
                    + pref * d_tt2 * dt3
                    + pref * minus_exp * d_tt3
                    + pref_tt * one_plus_mexp
                    + pref_t * d_tt2;
            }

            return Ok((r0 + d_int, r0_t + d_t_int, r0_tt + d_tt_int));
        }

        // Quantum Feynman-Hibbs path: solve r0 (where exp(-beta u(r0)) = machine_prec).
        let mut rs = 0.70;
        for _ in 0..25 {
            let r = rs * sigma_eff;
            let pot = self.evaluate(r, temperature);
            let f = pot.u_divk * inv_t + ln_eps;
            let df = sigma_eff * pot.u_divk_r * inv_t;
            let d2f = sigma_eff * sigma_eff * pot.u_divk_rr * inv_t;
            let step = -f / df * (1.0 + f * d2f / (2.0 * df * df));
            rs += step;
            if step.abs() < 1e-15 {
                break;
            }
        }

        let r0 = rs * sigma_eff;
        let pot0 = self.evaluate(r0, temperature);
        let r0_t = -(pot0.u_divk_t + ln_eps) / pot0.u_divk_r;
        let r0_tt = -(pot0.u_divk_tt + 2.0 * pot0.u_divk_tr * r0_t + pot0.u_divk_rr * r0_t * r0_t)
            / pot0.u_divk_r;

        // Kronrod 21 quadrature integration on [r0, sigma_eff].
        let mut d_int = 0.0;
        let mut d_t_int = 0.0;
        let mut d_tt_int = 0.0;

        let sigma_j2 = 0.5 * (sigma_eff - r0);
        let sigma_j2_t = 0.5 * (sig_t - r0_t);
        let sigma_j2_tt = 0.5 * (sig_tt - r0_tt);

        for i in 0..21 {
            let xk = X_NODES_21[i];
            let wk = W_NODES_21[i];

            let r = sigma_j2 * xk + sigma_j2 + r0;
            let r_t = sigma_j2_t * xk + sigma_j2_t + r0_t;
            let r_tt = sigma_j2_tt * xk + sigma_j2_tt + r0_tt;

            let pot = self.evaluate(r, temperature);

            let pref = sigma_j2 * wk;
            let pref_t = sigma_j2_t * wk;
            let pref_tt = sigma_j2_tt * wk;

            let minus_exp = -(-pot.u_divk * inv_t).exp();
            d_int += pref * (1.0 + minus_exp);

            let dt3 = pot.u_divk * inv_t2 - pot.u_divk_t * inv_t - pot.u_divk_r * r_t * inv_t;
            let d_t1 = pref;
            let d_t2 = minus_exp;
            let d_t4 = pref_t;
            let d_t5 = 1.0 + minus_exp;
            d_t_int += d_t1 * d_t2 * dt3 + d_t4 * d_t5;

            let d_tt1 = pref_t;
            let d_tt2 = minus_exp * dt3;
            let d_tt3 = 2.0 * pot.u_divk_t * inv_t2 + 2.0 * pot.u_divk_r * r_t * inv_t2
                - 2.0 * pot.u_divk * inv_t3
                - pot.u_divk_tt * inv_t
                - 2.0 * pot.u_divk_tr * r_t * inv_t
                - pot.u_divk_rr * (r_t * r_t) * inv_t
                - pot.u_divk_r * r_tt * inv_t;
            let d_tt4 = pref_tt;
            let d_tt5 = d_tt2;
            d_tt_int += d_tt1 * d_t2 * dt3
                + d_t1 * d_tt2 * dt3
                + d_t1 * d_t2 * d_tt3
                + d_tt4 * d_t5
                + d_t4 * d_tt5;
        }

        Ok((r0 + d_int, r0_t + d_t_int, r0_tt + d_tt_int))
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn calc_hardsphere_diameter_pure(
    temperature: f64,
    sigma: f64,
    eps_depth_divk: f64,
    lambda_a: f64,
    lambda_r: f64,
    c_fac: f64,
    d_param: f64,
    fh_order: i32,
) -> Result<(f64, f64, f64)> {
    let pair = PairPotential {
        sigma,
        eps_depth_divk,
        lambda_a,
        lambda_r,
        c_fac,
        d_param,
        fh_order,
        second_factor: 0.5,
    };
    pair.diameter(temperature, pair.effective_sigma(temperature))
}

/// The six coefficients of the polynomials of the van der Waals attractive
/// constant, the routine `calcFunAlpha` of the source.
pub(crate) fn alpha_coefficients(alpha: f64) -> [f64; 6] {
    let mut f_alpha = [0.0; 6];
    let alpha_pow = [1.0, alpha, alpha * alpha, alpha * alpha * alpha];
    for k in 0..6 {
        let mut num = 0.0;
        for p in 0..4 {
            num += PHI_MAT[p][k] * alpha_pow[p];
        }
        let mut den = 1.0;
        for p in 1..4 {
            den += PHI_MAT[p + 3][k] * alpha_pow[p];
        }
        f_alpha[k] = num / den;
    }
    f_alpha
}
