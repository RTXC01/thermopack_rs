//! Extended Corresponding States Principle (ExtCSP) equation of state.
//!
//! Pure-SI native port of `src/extcsp.f90` (Michelsen & Mollerup, ch. 2 / pp. 115-120),
//! combining a two-parameter cubic shape-factor EOS (`SRK` or `PR` with any supported
//! cubic mixing rule: `VdW`, `HuronVidal1`, `HuronVidal2`, `Nrtl`, `Unifac`,
//! `WongSandler`) with a pure-fluid reference EOS:
//! - NIST multiparameter EOS (`Multiparameter`, e.g. `"C3"` propane)
//! - Pure-fluid cubic EOS (`SRK` or `PR`)
//! - Pure-fluid MBWR-19 (Bender) or MBWR-32 EOS (`Mbwr`)
//!
//! Unlike the upstream Fortran `src/extcsp.f90`, this implementation operates
//! exclusively in base SI units (`K`, `Pa`, `m^3`, `mol`, `J`) and fixes four
//! upstream defects documented in `docs/rust/baseline-defects.md`:
//! 1. Cubic reference EOS initialization (`refEosType` left unset in `csp_init`).
//! 2. MBWR reference EOS derivative overwrite in `calcRefEqDiff`.
//! 3. `MBWR19` selection prefix check in `extcsp_eos_allocate_and_init`.
//! 4. Uninitialized `F` seed, wrong `N`/`sumn` symbol in `F_function_TWU`, and
//!    hardcoded classic derivatives in `a0_diff` for Twu alpha correlations.

use crate::error::positive;
use crate::{
    Alpha, AlphaCorrelation, Component, Cubic, CubicKind, CubicOptions, Error, IdealGas, Mbwr,
    MbwrDensityParameters, MbwrKind, MixingRule, Multiparameter, Phase, R, Result,
};

/// Reference equation of state specification for [`ExtCspOptions`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExtCspReferenceKind {
    /// NIST multiparameter reference equation (standard choice: `"C3"`).
    NistMeos,
    /// Pure-fluid two-parameter cubic reference equation (`SRK` or `PR`).
    Cubic {
        kind: CubicKind,
        alpha: AlphaCorrelation,
    },
    /// Pure-fluid MBWR-19 (Bender) or MBWR-32 reference equation.
    Mbwr(MbwrKind),
}

/// Configuration options for [`ExtCsp::with_options`].
#[derive(Debug, Clone, Copy)]
pub struct ExtCspOptions<'a> {
    /// Two-parameter cubic equation used for shape-factor calculation (`SRK` or `PR`).
    pub shape_eos: CubicKind,
    /// Mixing rule applied to the mixture shape-factor cubic EOS.
    pub shape_mixing: MixingRule,
    /// Alpha correlation used for the shape-factor cubic EOS and reference shape fluid.
    pub shape_alpha: AlphaCorrelation,
    /// Reference fluid identifier (default `"C3"`).
    pub reference_component: &'a str,
    /// Reference equation of state backend (default [`ExtCspReferenceKind::NistMeos`]).
    pub reference_eos: ExtCspReferenceKind,
    /// Parameter set reference passed to component/interaction/alpha databases.
    pub parameter_reference: &'a str,
    /// Optional binary parameters when `shape_mixing` is [`MixingRule::WongSandler`].
    pub wong_sandler: Option<&'a crate::wong_sandler::WongSandlerData>,
}

impl Default for ExtCspOptions<'_> {
    fn default() -> Self {
        Self {
            shape_eos: CubicKind::SoaveRedlichKwong,
            shape_mixing: MixingRule::VdW,
            shape_alpha: AlphaCorrelation::Classic,
            reference_component: "C3",
            reference_eos: ExtCspReferenceKind::NistMeos,
            parameter_reference: "DEFAULT",
            wong_sandler: None,
        }
    }
}

#[derive(Debug, Clone)]
enum ReferenceBackend {
    Multiparameter(Box<Multiparameter>),
    Cubic(Box<Cubic>),
    Mbwr {
        model: Box<Mbwr>,
        density_parameters: MbwrDensityParameters,
    },
}

/// Reduced residual Helmholtz energy `F = A_res / (R T)` and all first and
/// second derivatives with respect to `(T, V, n)` at an extensive TV state.
#[derive(Debug, Clone)]
pub struct ExtCspHelmholtz {
    pub f: f64,
    pub f_t: f64,
    pub f_v: f64,
    pub f_n: Vec<f64>,
    pub f_tt: f64,
    pub f_vv: f64,
    pub f_tv: f64,
    pub f_tn: Vec<f64>,
    pub f_vn: Vec<f64>,
    /// Row-major `nc * nc` symmetric matrix with entry `i * nc + j` equal to
    /// `d^2 F / (dn_i dn_j)`.
    pub f_nn: Vec<f64>,
}

/// Pressure `[Pa]` and its exact first derivatives at an extensive TV state.
#[derive(Debug, Clone)]
pub struct ExtCspPressure {
    pub value: f64,
    pub dt: f64,
    pub dv: f64,
    pub dn: Vec<f64>,
}

/// Extensive caloric state at fixed `(T, V, n)`.
#[derive(Debug, Clone)]
pub struct ExtCspCaloricTv {
    pub pressure: f64,
    pub internal_energy: f64,
    pub enthalpy: f64,
    pub entropy: f64,
    pub helmholtz_energy: f64,
    pub gibbs_energy: f64,
    pub heat_capacity_v: f64,
    pub heat_capacity_p: f64,
}

/// Scalar TP property value with analytical `(T, P, n)` derivatives.
///
/// For intensive quantities (`Z`), `dn[i] = dZ / dn_i`. For extensive
/// properties (`V`, `H`, `S`, `G`), `value`, `dt`, `dp` are extensive for the
/// supplied mole vector `n` (equal to molar when `sum(n) == 1.0`), and
/// `dn[i]` is the partial molar derivative `d(value) / dn_i`.
#[derive(Debug, Clone)]
pub struct ExtCspDerivatives {
    pub value: f64,
    pub dt: f64,
    pub dp: f64,
    pub dn: Vec<f64>,
}

/// Log fugacity coefficients with analytical `(T, P, n)` derivatives at
/// constant `(T, P)`.
#[derive(Debug, Clone)]
pub struct ExtCspFugacityDerivatives {
    pub ln_phi: Vec<f64>,
    pub dt: Vec<f64>,
    pub dp: Vec<f64>,
    /// Row-major `nc * nc` symmetric matrix with entry `i * nc + j` equal to
    /// `d ln(phi_i) / dn_j`.
    pub dn: Vec<f64>,
}

/// Extended Corresponding States Principle equation of state.
#[derive(Debug, Clone)]
pub struct ExtCsp {
    shape_eos: Cubic,
    ref_component: Component,
    ref_alpha: Alpha,
    ref_classic_m0: f64,
    reference: ReferenceBackend,
    tc0: f64,
    pc0: f64,
    ac0: f64,
    bc0: f64,
}

impl ExtCsp {
    /// Create a default Extended CSP model (`SRK`, `VdW`, `Classic`, `"C3"`, `NIST_MEOS`).
    pub fn new(components: &[&str]) -> Result<Self> {
        Self::with_options(components, ExtCspOptions::default())
    }

    /// Create an Extended CSP model with explicit shape-factor and reference EOS options.
    pub fn with_options(components: &[&str], options: ExtCspOptions<'_>) -> Result<Self> {
        if !matches!(
            options.shape_eos,
            CubicKind::SoaveRedlichKwong | CubicKind::PengRobinson
        ) {
            return Err(Error::InvalidInput(
                "ExtCSP shape EOS must be SoaveRedlichKwong or PengRobinson",
            ));
        }
        let shape_eos = Cubic::with_options(
            options.shape_eos,
            components,
            CubicOptions {
                alpha: Some(options.shape_alpha),
                mixing_rule: options.shape_mixing,
                parameter_reference: options.parameter_reference,
                volume_translation: false,
                wong_sandler: options.wong_sandler,
            },
        )?;
        let shape_eos_ref = Cubic::with_options(
            options.shape_eos,
            &[options.reference_component],
            CubicOptions {
                alpha: Some(options.shape_alpha),
                mixing_rule: MixingRule::VdW,
                parameter_reference: options.parameter_reference,
                volume_translation: false,
                ..Default::default()
            },
        )?;
        let ref_component = shape_eos_ref.components()[0];
        let ref_alpha = *shape_eos_ref.alpha(0)?;
        let classic_corr = match options.shape_eos {
            CubicKind::SoaveRedlichKwong => AlphaCorrelation::Soave,
            CubicKind::PengRobinson => AlphaCorrelation::PengRobinson,
            _ => AlphaCorrelation::Classic,
        };
        let ref_classic_m0 = classic_corr.acentric_parameter(ref_component.acentric_factor)?;

        let (reference, tc0, pc0) = match options.reference_eos {
            ExtCspReferenceKind::NistMeos => {
                let meos = Multiparameter::new(options.reference_component)?;
                let tc0 = meos.tc();
                let pc0 = meos.pc();
                (ReferenceBackend::Multiparameter(Box::new(meos)), tc0, pc0)
            }
            ExtCspReferenceKind::Cubic { kind, alpha } => {
                if !matches!(kind, CubicKind::SoaveRedlichKwong | CubicKind::PengRobinson) {
                    return Err(Error::InvalidInput(
                        "ExtCSP cubic reference EOS must be SoaveRedlichKwong or PengRobinson",
                    ));
                }
                let cb = Cubic::with_options(
                    kind,
                    &[options.reference_component],
                    CubicOptions {
                        alpha: Some(alpha),
                        mixing_rule: MixingRule::VdW,
                        parameter_reference: options.parameter_reference,
                        volume_translation: false,
                        ..Default::default()
                    },
                )?;
                let tc0 = cb.components()[0].critical_temperature;
                let pc0 = cb.components()[0].critical_pressure;
                (ReferenceBackend::Cubic(Box::new(cb)), tc0, pc0)
            }
            ExtCspReferenceKind::Mbwr(kind) => {
                let mbwr = Mbwr::new(kind, options.reference_component)?;
                let density_parameters = mbwr.source_density_parameters()?;
                let tc0 = density_parameters.critical_temperature;
                let pc0 = density_parameters.critical_pressure;
                (
                    ReferenceBackend::Mbwr {
                        model: Box::new(mbwr),
                        density_parameters,
                    },
                    tc0,
                    pc0,
                )
            }
        };

        positive(tc0, "reference critical temperature")?;
        positive(pc0, "reference critical pressure")?;
        let (omega_a, omega_b) = crate::cubic::omega_ab(options.shape_eos);
        let rt0 = R * tc0;
        let bc0 = omega_b * rt0 / pc0;
        let ac0 = omega_a * (rt0 * rt0) / pc0;

        Ok(Self {
            shape_eos,
            ref_component,
            ref_alpha,
            ref_classic_m0,
            reference,
            tc0,
            pc0,
            ac0,
            bc0,
        })
    }

    /// Components of the mixture.
    pub fn components(&self) -> &[Component] {
        self.shape_eos.components()
    }

    /// Reference component metadata.
    pub fn reference_component(&self) -> Component {
        self.ref_component
    }

    /// Reference EOS critical temperature `[K]`.
    pub fn reference_critical_temperature(&self) -> f64 {
        self.tc0
    }

    /// Reference EOS critical pressure `[Pa]`.
    pub fn reference_critical_pressure(&self) -> f64 {
        self.pc0
    }

    /// Reference shape attraction constant `ac0` `[Pa m^6 / mol^2]`.
    pub fn reference_ac0(&self) -> f64 {
        self.ac0
    }

    /// Reference shape covolume constant `bc0` `[m^3 / mol]`.
    pub fn reference_bc0(&self) -> f64 {
        self.bc0
    }

    /// Borrow the underlying cubic shape-factor EOS.
    pub fn shape_eos(&self) -> &Cubic {
        &self.shape_eos
    }

    /// Mutably borrow the underlying cubic shape-factor EOS.
    pub fn shape_eos_mut(&mut self) -> &mut Cubic {
        &mut self.shape_eos
    }

    /// Get binary attraction interaction parameter `kij`.
    pub fn kij(&self, i: usize, j: usize) -> Result<f64> {
        self.shape_eos.kij(i, j)
    }

    /// Set binary attraction interaction parameter `kij`.
    pub fn set_kij(&mut self, i: usize, j: usize, value: f64) -> Result<()> {
        self.shape_eos.set_kij(i, j, value)
    }

    /// Set binary covolume interaction parameter `lij`.
    pub fn set_lij(&mut self, i: usize, j: usize, value: f64) -> Result<()> {
        self.shape_eos.set_lij(i, j, value)
    }

    /// Active ideal-gas parameters for component `index`.
    pub fn ideal_gas(&self, index: usize) -> Result<&IdealGas> {
        self.shape_eos.ideal_gas(index)
    }

    /// Mutably access ideal-gas parameters for component `index`.
    pub fn ideal_gas_mut(&mut self, index: usize) -> Result<&mut IdealGas> {
        self.shape_eos.ideal_gas_mut(index)
    }

    fn validate_moles(&self, moles: &[f64]) -> Result<(f64, Vec<f64>)> {
        if moles.len() != self.shape_eos.components().len() {
            return Err(Error::InvalidInput("composition length mismatch"));
        }
        let mut sum = 0.0;
        for &ni in moles {
            if !ni.is_finite() || ni < 0.0 {
                return Err(Error::InvalidInput(
                    "non-negative finite mole numbers required",
                ));
            }
            sum += ni;
        }
        positive(sum, "total moles")?;
        let x: Vec<f64> = moles.iter().map(|&ni| ni / sum).collect();
        Ok((sum, x))
    }

    fn is_classic_m0_alpha(&self) -> bool {
        matches!(
            self.ref_alpha.correlation,
            AlphaCorrelation::Classic
                | AlphaCorrelation::Soave
                | AlphaCorrelation::PengRobinson
                | AlphaCorrelation::GraboskiDaubert
                | AlphaCorrelation::PengRobinson78
                | AlphaCorrelation::ClassicFit
        )
    }

    /// Solve `F * alpha_ref(N * T / F) = D / (H * ac0)` for the temperature shape factor `F` `[mol]`.
    fn solve_temperature_shape_factor(&self, h: f64, d: f64, t: f64, moles: f64) -> Result<f64> {
        positive(h, "volume shape factor H")?;
        positive(d, "mixture attraction parameter D")?;
        let target = d / (h * self.ac0);
        if self.is_classic_m0_alpha() {
            let m0 = self.ref_alpha.parameters[0];
            let root = (m0 * (t * moles / self.tc0).sqrt() + target.sqrt()) / (1.0 + m0);
            let f = root * root;
            positive(f, "temperature shape factor F")?;
            return Ok(f);
        }

        // General alpha (including Twu): solve g(F) = F * alpha0(N*T/F) - target = 0
        // starting from the classic-alpha closed-form seed.
        let m0 = self.ref_classic_m0;
        let root_guess = (m0 * (t * moles / self.tc0).sqrt() + target.sqrt()) / (1.0 + m0);
        let mut f = (root_guess * root_guess).max(1.0e-12 * moles);
        let tol = 1.0e-14 * target.max(1.0);
        let eval_g = |f_val: f64| -> Result<(f64, f64)> {
            let t0 = moles * t / f_val;
            let av = self.ref_alpha.evaluate(t0, self.tc0)?;
            let g = f_val * av.value - target;
            let dg = av.value - t0 * av.dt;
            Ok((g, dg))
        };

        for _ in 0..80 {
            let (g, dg) = eval_g(f)?;
            if !g.is_finite() || !dg.is_finite() || dg == 0.0 {
                break;
            }
            if g.abs() <= tol {
                return Ok(f);
            }
            let step = g / dg;
            let mut next = f - step;
            let mut damping = 1.0;
            while (!next.is_finite() || next <= 0.0) && damping > 1.0e-8 {
                damping *= 0.5;
                next = f - damping * step;
            }
            if next <= 0.0 || !next.is_finite() {
                break;
            }
            let (g_next, _) = eval_g(next)?;
            if g_next.abs() > g.abs() {
                // Backtrack in log space to keep F positive.
                let mut alpha = 0.5;
                let mut accepted = false;
                for _ in 0..24 {
                    let cand = f - alpha * step;
                    if cand > 0.0
                        && cand.is_finite()
                        && let Ok((gc, _)) = eval_g(cand)
                        && gc.abs() < g.abs()
                    {
                        next = cand;
                        accepted = true;
                        break;
                    }
                    alpha *= 0.5;
                }
                if !accepted {
                    break;
                }
            }
            let rel = (next - f).abs() / next.max(1.0e-30);
            f = next;
            if rel <= 1.0e-14 {
                let (gf, _) = eval_g(f)?;
                if gf.abs() <= 1.0e-11 * target.max(1.0) {
                    return Ok(f);
                }
            }
        }

        // Bracketed bisection fallback using monotonicity of F * alpha0(N*T/F).
        let mut lo = f * 0.5;
        let mut hi = f * 2.0;
        for _ in 0..80 {
            if eval_g(lo).is_ok_and(|(g, _)| g <= 0.0) {
                break;
            }
            lo *= 0.5;
        }
        for _ in 0..80 {
            if eval_g(hi).is_ok_and(|(g, _)| g >= 0.0) {
                break;
            }
            hi *= 2.0;
        }
        let (glo, _) = eval_g(lo)?;
        let (ghi, _) = eval_g(hi)?;
        if glo > 0.0 || ghi < 0.0 {
            return Err(Error::NotConverged(
                "ExtCSP temperature shape factor F solver failed to bracket root",
            ));
        }
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            let (gm, _) = eval_g(mid)?;
            if gm.abs() <= tol || (hi - lo) <= 1.0e-14 * mid {
                return Ok(mid);
            }
            if gm < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(0.5 * (lo + hi))
    }

    /// Evaluate reference reduced residual Helmholtz `alpr(T0, v0)` and its
    /// first and second `(T0, v0)` derivatives `[alpr, alpr_t, alpr_v, alpr_tt, alpr_tv, alpr_vv]`.
    fn reference_alpr_jet(&self, t0: f64, v0: f64) -> Result<[f64; 6]> {
        match &self.reference {
            ReferenceBackend::Multiparameter(meos) => {
                let f = meos.prepare(t0)?.residual_helmholtz(v0, 1.0)?;
                Ok([f.value, f.dt, f.dv, f.dtt, f.dtv, f.dvv])
            }
            ReferenceBackend::Cubic(cb) => cb.prepare(t0, &[1.0])?.pure_residual_helmholtz_tv(v0),
            ReferenceBackend::Mbwr { model, .. } => {
                let f = model.prepare(t0)?.residual_helmholtz(v0, 1.0)?;
                Ok([f.value, f.dt, f.dv, f.dtt, f.dtv, f.dvv])
            }
        }
    }

    /// Solve reference EOS molar volume `v0` `[m^3/mol]` at `(t0, p0, phase)`.
    fn reference_volume(&self, t0: f64, p0: f64, phase: Phase) -> Result<f64> {
        match &self.reference {
            ReferenceBackend::Multiparameter(meos) => {
                let prep = meos.prepare(t0)?;
                match phase {
                    Phase::Liquid | Phase::Vapor => {
                        let iter = prep.baseline_density(p0, phase)?;
                        if iter.pressure_converged
                            && iter.mechanically_stable
                            && iter.density > 0.0
                            && iter.density.is_finite()
                        {
                            Ok(1.0 / iter.density)
                        } else {
                            Err(Error::NotConverged(
                                "ExtCSP multiparameter reference density did not converge",
                            ))
                        }
                    }
                    Phase::MinimumGibbs => {
                        let vl = prep.volume(p0, Phase::Liquid);
                        let vv = prep.volume(p0, Phase::Vapor);
                        match (vl, vv) {
                            (Ok(v1), Ok(v2)) => {
                                let g1 = prep.ln_fugacity(v1)?;
                                let g2 = prep.ln_fugacity(v2)?;
                                Ok(if g1 <= g2 { v1 } else { v2 })
                            }
                            (Ok(v1), Err(_)) => Ok(v1),
                            (Err(_), Ok(v2)) => Ok(v2),
                            (Err(e), Err(_)) => Err(e),
                        }
                    }
                }
            }
            ReferenceBackend::Cubic(cb) => cb.prepare(t0, &[1.0])?.volume(p0, phase),
            ReferenceBackend::Mbwr {
                model,
                density_parameters,
            } => {
                let prep_d = model.prepare_baseline_density(t0, *density_parameters)?;
                let solve_phase = |ph: Phase| -> Result<f64> {
                    let iter = prep_d.density(p0, ph)?;
                    if iter.pressure_converged
                        && iter.mechanically_stable
                        && iter.density > 0.0
                        && iter.density.is_finite()
                    {
                        Ok(1.0 / iter.density)
                    } else {
                        Err(Error::NotConverged(
                            "ExtCSP MBWR reference density did not converge",
                        ))
                    }
                };
                match phase {
                    Phase::Liquid | Phase::Vapor => solve_phase(phase),
                    Phase::MinimumGibbs => {
                        let vl = solve_phase(Phase::Liquid);
                        let vv = solve_phase(Phase::Vapor);
                        match (vl, vv) {
                            (Ok(v1), Ok(v2)) => {
                                let prep = model.prepare(t0)?;
                                let f1 = prep.residual_helmholtz(v1, 1.0)?;
                                let f2 = prep.residual_helmholtz(v2, 1.0)?;
                                let z1 = p0 * v1 / (R * t0);
                                let z2 = p0 * v2 / (R * t0);
                                let g1 = f1.value + z1 - 1.0 - z1.ln();
                                let g2 = f2.value + z2 - 1.0 - z2.ln();
                                Ok(if g1 <= g2 { v1 } else { v2 })
                            }
                            (Ok(v1), Err(_)) => Ok(v1),
                            (Err(_), Ok(v2)) => Ok(v2),
                            (Err(e), Err(_)) => Err(e),
                        }
                    }
                }
            }
        }
    }

    /// Compute shape factors `(H, F)` (`[mol]`) at `(T, n)`.
    pub fn shape_factors(&self, t: f64, moles: &[f64]) -> Result<(f64, f64)> {
        positive(t, "temperature")?;
        let (n_tot, x) = self.validate_moles(moles)?;
        let jet = self
            .shape_eos
            .prepare(t, &x)?
            .mixture_parameter_jet(n_tot)?;
        let h = jet.b / self.bc0;
        let f = self.solve_temperature_shape_factor(h, jet.d, t, n_tot)?;
        Ok((h, f))
    }

    /// Reduced residual Helmholtz energy `F = A_res / (R T)` and all 10 first
    /// and second derivative blocks with respect to `(T, V, n)`.
    ///
    /// Inputs are temperature `t` `[K]`, total volume `volume` `[m^3]`, and
    /// mole vector `moles` `[mol]`.
    pub fn residual_helmholtz_tv(
        &self,
        t: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<ExtCspHelmholtz> {
        positive(t, "temperature")?;
        positive(volume, "volume")?;
        let (n_tot, x) = self.validate_moles(moles)?;
        let nc = moles.len();

        let jet = self
            .shape_eos
            .prepare(t, &x)?
            .mixture_parameter_jet(n_tot)?;
        let inv_bc0 = 1.0 / self.bc0;

        let h = jet.b * inv_bc0;
        let h_t = jet.b_t * inv_bc0;
        let h_tt = jet.b_tt * inv_bc0;
        let h_i: Vec<f64> = jet.b_i.iter().map(|&v| v * inv_bc0).collect();
        let h_it: Vec<f64> = jet.b_it.iter().map(|&v| v * inv_bc0).collect();
        let h_ij: Vec<f64> = jet.b_ij.iter().map(|&v| v * inv_bc0).collect();

        let f_shape = self.solve_temperature_shape_factor(h, jet.d, t, n_tot)?;
        let t0 = n_tot * t / f_shape;
        let v0 = volume / h;

        let alpha_val = self.ref_alpha.evaluate(t0, self.tc0)?;
        let a0 = self.ac0 * alpha_val.value;
        let a0_t0 = self.ac0 * alpha_val.dt;
        let a0_t0t0 = self.ac0 * alpha_val.dtt;

        let g_f = a0 - t0 * a0_t0;
        if !g_f.is_finite() || g_f == 0.0 {
            return Err(Error::InvalidInput(
                "singular ExtCSP temperature shape factor",
            ));
        }

        // Q = D / H and its (T, n) derivatives.
        let inv_h = 1.0 / h;
        let q = jet.d * inv_h;
        let q_t = (jet.d_t - q * h_t) * inv_h;
        let q_tt = (jet.d_tt - 2.0 * q_t * h_t - q * h_tt) * inv_h;

        let mut q_i = vec![0.0; nc];
        let mut q_it = vec![0.0; nc];
        for i in 0..nc {
            q_i[i] = (jet.d_i[i] - q * h_i[i]) * inv_h;
            q_it[i] = (jet.d_it[i] - q_t * h_i[i] - q_i[i] * h_t - q * h_it[i]) * inv_h;
        }

        // First derivatives of F and t0 = (N * T) / F.
        let inv_f = 1.0 / f_shape;
        let inv_gf = 1.0 / g_f;

        let f_t = (q_t - a0_t0 * n_tot) * inv_gf;
        let t0_t = (n_tot - t0 * f_t) * inv_f;
        let f_tt = (q_tt - f_shape * a0_t0t0 * t0_t * t0_t) * inv_gf;
        let t0_tt = (-f_tt * t0 - 2.0 * f_t * t0_t) * inv_f;

        let mut f_i = vec![0.0; nc];
        let mut t0_i = vec![0.0; nc];
        let mut f_it = vec![0.0; nc];
        let mut t0_it = vec![0.0; nc];
        for i in 0..nc {
            f_i[i] = (q_i[i] - a0_t0 * t) * inv_gf;
            t0_i[i] = (t - t0 * f_i[i]) * inv_f;
            f_it[i] = (q_it[i] - a0_t0 - f_shape * a0_t0t0 * t0_i[i] * t0_t) * inv_gf;
            t0_it[i] = (1.0 - f_it[i] * t0 - f_i[i] * t0_t - f_t * t0_i[i]) * inv_f;
        }

        // v0 = V / H and all its (T, V, n) derivatives.
        let v0_v = inv_h;
        let v0_t = -v0 * h_t * inv_h;
        let v0_tv = -h_t * inv_h * inv_h;
        let v0_tt = (2.0 * h_t * h_t * inv_h - h_tt) * v0 * inv_h;

        let mut v0_i = vec![0.0; nc];
        let mut v0_vi = vec![0.0; nc];
        let mut v0_ti = vec![0.0; nc];
        for i in 0..nc {
            v0_i[i] = -v0 * h_i[i] * inv_h;
            v0_vi[i] = -h_i[i] * inv_h * inv_h;
            v0_ti[i] = (2.0 * h_t * h_i[i] * inv_h - h_it[i]) * v0 * inv_h;
        }

        // Reference EOS reduced residual Helmholtz jet at (t0, v0).
        let [a, a_t, a_v, a_tt, a_tv, a_vv] = self.reference_alpr_jet(t0, v0)?;

        // Combine F_res(T, V, n) = N * alpr(t0(T, n), v0(T, V, n)).
        let da_t = a_t * t0_t + a_v * v0_t;
        let da_v = a_v * v0_v;
        let d2a_tt = a_tt * t0_t * t0_t
            + 2.0 * a_tv * t0_t * v0_t
            + a_vv * v0_t * v0_t
            + a_t * t0_tt
            + a_v * v0_tt;
        let d2a_vv = a_vv * v0_v * v0_v;
        let d2a_tv = (a_tv * t0_t + a_vv * v0_t) * v0_v + a_v * v0_tv;

        let mut da_i = vec![0.0; nc];
        let mut f_n = vec![0.0; nc];
        let mut f_tn = vec![0.0; nc];
        let mut f_vn = vec![0.0; nc];
        for i in 0..nc {
            let dai = a_t * t0_i[i] + a_v * v0_i[i];
            da_i[i] = dai;
            f_n[i] = a + n_tot * dai;
            let d2a_ti = a_tt * t0_t * t0_i[i]
                + a_tv * (t0_t * v0_i[i] + t0_i[i] * v0_t)
                + a_vv * v0_t * v0_i[i]
                + a_t * t0_it[i]
                + a_v * v0_ti[i];
            f_tn[i] = da_t + n_tot * d2a_ti;
            let d2a_vi = (a_tv * t0_i[i] + a_vv * v0_i[i]) * v0_v + a_v * v0_vi[i];
            f_vn[i] = da_v + n_tot * d2a_vi;
        }

        let mut f_nn = vec![0.0; nc * nc];
        for i in 0..nc {
            for j in i..nc {
                let idx = i * nc + j;
                let q_ij =
                    (jet.d_ij[idx] - q_i[i] * h_i[j] - q_i[j] * h_i[i] - q * h_ij[idx]) * inv_h;
                let f_ij = (q_ij - f_shape * a0_t0t0 * t0_i[i] * t0_i[j]) * inv_gf;
                let t0_ij = (-f_ij * t0 - f_i[i] * t0_i[j] - f_i[j] * t0_i[i]) * inv_f;
                let v0_ij = (2.0 * h_i[i] * h_i[j] * inv_h - h_ij[idx]) * v0 * inv_h;

                let d2a_ij = a_tt * t0_i[i] * t0_i[j]
                    + a_tv * (t0_i[i] * v0_i[j] + t0_i[j] * v0_i[i])
                    + a_vv * v0_i[i] * v0_i[j]
                    + a_t * t0_ij
                    + a_v * v0_ij;
                let val = da_i[i] + da_i[j] + n_tot * d2a_ij;
                f_nn[idx] = val;
                f_nn[j * nc + i] = val;
            }
        }

        let result = ExtCspHelmholtz {
            f: n_tot * a,
            f_t: n_tot * da_t,
            f_v: n_tot * da_v,
            f_n,
            f_tt: n_tot * d2a_tt,
            f_vv: n_tot * d2a_vv,
            f_tv: n_tot * d2a_tv,
            f_tn,
            f_vn,
            f_nn,
        };
        if !result.f.is_finite()
            || !result.f_t.is_finite()
            || !result.f_v.is_finite()
            || !result.f_tt.is_finite()
            || !result.f_vv.is_finite()
            || !result.f_tv.is_finite()
            || result.f_n.iter().any(|v| !v.is_finite())
            || result.f_tn.iter().any(|v| !v.is_finite())
            || result.f_vn.iter().any(|v| !v.is_finite())
            || result.f_nn.iter().any(|v| !v.is_finite())
        {
            return Err(Error::InvalidInput(
                "non-finite ExtCSP Helmholtz derivatives",
            ));
        }
        Ok(result)
    }

    /// Pressure `[Pa]` and exact `(T, V, n)` derivatives at `(T, V, n)`.
    pub fn pressure_tv(&self, t: f64, volume: f64, moles: &[f64]) -> Result<ExtCspPressure> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let f = self.residual_helmholtz_tv(t, volume, moles)?;
        let rt = R * t;
        let value = rt * (n_tot / volume - f.f_v);
        let dt = value / t - rt * f.f_tv;
        let dv = -rt * (f.f_vv + n_tot / (volume * volume));
        let dn: Vec<f64> = f
            .f_vn
            .iter()
            .map(|&fvn| rt * (1.0 / volume - fvn))
            .collect();
        Ok(ExtCspPressure { value, dt, dv, dn })
    }

    /// Extensive total volume `V` `[m^3]` solving `P(T, V, n) = p` for the
    /// requested phase. Equal to molar volume `[m^3/mol]` when `sum(n) == 1.0`.
    pub fn volume_tp(&self, t: f64, p: f64, moles: &[f64], phase: Phase) -> Result<f64> {
        positive(t, "temperature")?;
        positive(p, "pressure")?;
        let (n_tot, _) = self.validate_moles(moles)?;
        let (h, f_shape) = self.shape_factors(t, moles)?;
        let t0 = n_tot * t / f_shape;
        let p0 = h * p / f_shape;
        let v0 = self.reference_volume(t0, p0, phase)?;
        let v = h * v0;
        positive(v, "ExtCSP volume")?;
        Ok(v)
    }

    /// Molar specific volume `[m^3/mol]` at `(T, P, n, phase)`.
    pub fn specific_volume(&self, t: f64, p: f64, moles: &[f64], phase: Phase) -> Result<f64> {
        let (n_tot, _) = self.validate_moles(moles)?;
        Ok(self.volume_tp(t, p, moles, phase)? / n_tot)
    }

    /// Extensive total volume `V` `[m^3]` with exact `(T, P, n)` derivatives.
    /// When `sum(n) == 1.0`, `value`, `dt`, `dp` are molar and `dn` are partial
    /// molar volumes.
    pub fn specific_volume_derivatives(
        &self,
        t: f64,
        p: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<ExtCspDerivatives> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let v = self.volume_tp(t, p, moles, phase)?;
        let f = self.residual_helmholtz_tv(t, v, moles)?;
        let rt = R * t;
        let dp_dv = -rt * (f.f_vv + n_tot / (v * v));
        let dp_dt = p / t - rt * f.f_tv;
        let dv_dt = -dp_dt / dp_dv;
        let dv_dp = 1.0 / dp_dv;
        let dv_dn: Vec<f64> = f
            .f_vn
            .iter()
            .map(|&fvn| -rt * (1.0 / v - fvn) / dp_dv)
            .collect();
        Ok(ExtCspDerivatives {
            value: v,
            dt: dv_dt,
            dp: dv_dp,
            dn: dv_dn,
        })
    }

    /// Compressibility factor `Z = P V / (N R T)` at `(T, P, n, phase)`.
    pub fn compressibility_factor(
        &self,
        t: f64,
        p: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<f64> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let v = self.volume_tp(t, p, moles, phase)?;
        Ok(p * v / (n_tot * R * t))
    }

    /// Compressibility factor `Z` with exact `(T, P, n)` derivatives.
    pub fn compressibility_factor_derivatives(
        &self,
        t: f64,
        p: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<ExtCspDerivatives> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let vd = self.specific_volume_derivatives(t, p, moles, phase)?;
        let v = vd.value;
        let z = p * v / (n_tot * R * t);
        let dz_dt = -z * (1.0 / t - vd.dt / v);
        let dz_dp = z * (1.0 / p + vd.dp / v);
        let dz_dn: Vec<f64> = vd
            .dn
            .iter()
            .map(|&dvdni| -z * (1.0 / n_tot - dvdni / v))
            .collect();
        Ok(ExtCspDerivatives {
            value: z,
            dt: dz_dt,
            dp: dz_dp,
            dn: dz_dn,
        })
    }

    /// Natural logarithm of component fugacity coefficients `ln(phi_i)` at `(T, P, n, phase)`.
    pub fn ln_fugacity_tp(&self, t: f64, p: f64, moles: &[f64], phase: Phase) -> Result<Vec<f64>> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let v = self.volume_tp(t, p, moles, phase)?;
        let f = self.residual_helmholtz_tv(t, v, moles)?;
        let z = p * v / (n_tot * R * t);
        positive(z, "compressibility factor for fugacity")?;
        let ln_z = z.ln();
        Ok(f.f_n.iter().map(|&fn_i| fn_i - ln_z).collect())
    }

    /// Log fugacity coefficients with exact `(T, P, n)` derivatives at `(T, P, n, phase)`.
    pub fn ln_fugacity_derivatives_tp(
        &self,
        t: f64,
        p: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<ExtCspFugacityDerivatives> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let nc = moles.len();
        let v = self.volume_tp(t, p, moles, phase)?;
        let f = self.residual_helmholtz_tv(t, v, moles)?;
        let rt = R * t;
        let z = p * v / (n_tot * rt);
        positive(z, "compressibility factor for fugacity")?;
        let ln_z = z.ln();

        let dp_dv = -rt * (f.f_vv + n_tot / (v * v));
        let dp_dt = p / t - rt * f.f_tv;
        let dp_dn: Vec<f64> = f.f_vn.iter().map(|&fvn| rt * (1.0 / v - fvn)).collect();
        let dv_dn: Vec<f64> = dp_dn.iter().map(|&dpni| -dpni / dp_dv).collect();

        let ln_phi: Vec<f64> = f.f_n.iter().map(|&fn_i| fn_i - ln_z).collect();
        let dt: Vec<f64> = (0..nc)
            .map(|i| f.f_tn[i] + (1.0 - dv_dn[i] * dp_dt / R) / t)
            .collect();
        let dp: Vec<f64> = (0..nc).map(|i| dv_dn[i] / rt - 1.0 / p).collect();
        let mut dn = vec![0.0; nc * nc];
        let inv_n = 1.0 / n_tot;
        for i in 0..nc {
            for j in 0..nc {
                dn[i * nc + j] = f.f_nn[i * nc + j] + inv_n - dv_dn[j] * dp_dn[i] / rt;
            }
        }
        Ok(ExtCspFugacityDerivatives { ln_phi, dt, dp, dn })
    }

    /// Extensive residual enthalpy `H_res` `[J]` at `(T, P, n, phase)`.
    pub fn residual_enthalpy(&self, t: f64, p: f64, moles: &[f64], phase: Phase) -> Result<f64> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let v = self.volume_tp(t, p, moles, phase)?;
        let f = self.residual_helmholtz_tv(t, v, moles)?;
        Ok(-R * t * t * f.f_t + p * v - n_tot * R * t)
    }

    /// Extensive residual enthalpy `H_res` `[J]` with exact `(T, P, n)` derivatives.
    pub fn residual_enthalpy_derivatives(
        &self,
        t: f64,
        p: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<ExtCspDerivatives> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let v = self.volume_tp(t, p, moles, phase)?;
        let f = self.residual_helmholtz_tv(t, v, moles)?;
        let rt = R * t;
        let dp_dv = -rt * (f.f_vv + n_tot / (v * v));
        let dp_dt = p / t - rt * f.f_tv;
        let dv_dt = -dp_dt / dp_dv;
        let value = -rt * t * f.f_t + p * v - n_tot * rt;
        let dt = t * (dv_dt * dp_dt - R * (2.0 * f.f_t + t * f.f_tt + n_tot / t));
        let dp = v - t * dv_dt;
        let dn: Vec<f64> = f
            .f_vn
            .iter()
            .zip(f.f_tn.iter())
            .map(|(&fvn, &ftn)| {
                let dv_dni = -rt * (1.0 / v - fvn) / dp_dv;
                t * (dv_dni * dp_dt - R * (t * ftn + 1.0))
            })
            .collect();
        Ok(ExtCspDerivatives { value, dt, dp, dn })
    }

    /// Extensive total enthalpy `H` `[J]` at `(T, P, n, phase)`.
    pub fn enthalpy(&self, t: f64, p: f64, moles: &[f64], phase: Phase) -> Result<f64> {
        let h_res = self.residual_enthalpy(t, p, moles, phase)?;
        let mut h_id = 0.0;
        for (i, &ni) in moles.iter().enumerate() {
            if ni == 0.0 {
                continue;
            }
            h_id += ni * self.shape_eos.ideal_gas(i)?.properties(t, p)?.enthalpy;
        }
        Ok(h_id + h_res)
    }

    /// Extensive total enthalpy `H` `[J]` with exact `(T, P, n)` derivatives.
    pub fn enthalpy_derivatives(
        &self,
        t: f64,
        p: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<ExtCspDerivatives> {
        let mut res = self.residual_enthalpy_derivatives(t, p, moles, phase)?;
        for (i, &ni) in moles.iter().enumerate() {
            let ideal = self.shape_eos.ideal_gas(i)?.properties(t, p)?;
            res.value += ni * ideal.enthalpy;
            res.dt += ni * ideal.enthalpy_dt;
            res.dn[i] += ideal.enthalpy;
        }
        Ok(res)
    }

    /// Extensive residual entropy `S_res` `[J/K]` at `(T, P, n, phase)`.
    pub fn residual_entropy(&self, t: f64, p: f64, moles: &[f64], phase: Phase) -> Result<f64> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let v = self.volume_tp(t, p, moles, phase)?;
        let f = self.residual_helmholtz_tv(t, v, moles)?;
        let z = p * v / (n_tot * R * t);
        positive(z, "compressibility factor for entropy")?;
        Ok(R * (-f.f - t * f.f_t + n_tot * z.ln()))
    }

    /// Extensive residual entropy `S_res` `[J/K]` with exact `(T, P, n)` derivatives.
    pub fn residual_entropy_derivatives(
        &self,
        t: f64,
        p: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<ExtCspDerivatives> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let v = self.volume_tp(t, p, moles, phase)?;
        let f = self.residual_helmholtz_tv(t, v, moles)?;
        let rt = R * t;
        let z = p * v / (n_tot * rt);
        positive(z, "compressibility factor for entropy")?;
        let ln_z = z.ln();

        let dp_dv = -rt * (f.f_vv + n_tot / (v * v));
        let dp_dt = p / t - rt * f.f_tv;
        let dv_dt = -dp_dt / dp_dv;
        let value = R * (-f.f - t * f.f_t + n_tot * ln_z);
        let dt = dv_dt * dp_dt - R * (2.0 * f.f_t + t * f.f_tt + n_tot / t);
        let dp = n_tot * R / p - dv_dt;
        let dn: Vec<f64> = (0..moles.len())
            .map(|i| {
                let dv_dni = -rt * (1.0 / v - f.f_vn[i]) / dp_dv;
                dv_dni * dp_dt - R * (f.f_n[i] + t * f.f_tn[i] + 1.0 - ln_z)
            })
            .collect();
        Ok(ExtCspDerivatives { value, dt, dp, dn })
    }

    /// Extensive total entropy `S` `[J/K]` at `(T, P, n, phase)`.
    pub fn entropy(&self, t: f64, p: f64, moles: &[f64], phase: Phase) -> Result<f64> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let s_res = self.residual_entropy(t, p, moles, phase)?;
        let mut s_id = 0.0;
        for (i, &ni) in moles.iter().enumerate() {
            if ni == 0.0 {
                continue;
            }
            let xi = ni / n_tot;
            let ideal = self.shape_eos.ideal_gas(i)?.properties(t, p)?;
            s_id += ni * (ideal.entropy - R * xi.ln());
        }
        Ok(s_id + s_res)
    }

    /// Extensive total entropy `S` `[J/K]` with exact `(T, P, n)` derivatives.
    pub fn entropy_derivatives(
        &self,
        t: f64,
        p: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<ExtCspDerivatives> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let mut res = self.residual_entropy_derivatives(t, p, moles, phase)?;
        res.dp -= n_tot * R / p;
        for (i, &ni) in moles.iter().enumerate() {
            let xi = ni / n_tot;
            let ideal = self.shape_eos.ideal_gas(i)?.properties(t, p)?;
            if ni > 0.0 {
                res.value += ni * (ideal.entropy - R * xi.ln());
            }
            res.dt += ni * ideal.enthalpy_dt / t;
            res.dn[i] += ideal.entropy - R * xi.max(1.0e-100).ln();
        }
        Ok(res)
    }

    /// Extensive residual Gibbs energy `G_res` `[J]` at `(T, P, n, phase)`.
    pub fn residual_gibbs(&self, t: f64, p: f64, moles: &[f64], phase: Phase) -> Result<f64> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let v = self.volume_tp(t, p, moles, phase)?;
        let f = self.residual_helmholtz_tv(t, v, moles)?;
        let z = p * v / (n_tot * R * t);
        positive(z, "compressibility factor for Gibbs energy")?;
        Ok(R * t * (f.f + n_tot * (z - 1.0 - z.ln())))
    }

    /// Extensive residual Gibbs energy `G_res` `[J]` with exact `(T, P, n)` derivatives.
    pub fn residual_gibbs_derivatives(
        &self,
        t: f64,
        p: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<ExtCspDerivatives> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let v = self.volume_tp(t, p, moles, phase)?;
        let f = self.residual_helmholtz_tv(t, v, moles)?;
        let rt = R * t;
        let z = p * v / (n_tot * rt);
        positive(z, "compressibility factor for Gibbs energy")?;
        let ln_z = z.ln();
        let value = rt * (f.f + n_tot * (z - 1.0 - ln_z));
        let dt = -R * (-f.f - t * f.f_t + n_tot * ln_z);
        let dp = v - n_tot * rt / p;
        let dn: Vec<f64> = f.f_n.iter().map(|&fn_i| rt * (fn_i - ln_z)).collect();
        Ok(ExtCspDerivatives { value, dt, dp, dn })
    }

    /// Extensive total Gibbs energy `G` `[J]` at `(T, P, n, phase)`.
    pub fn gibbs(&self, t: f64, p: f64, moles: &[f64], phase: Phase) -> Result<f64> {
        let h = self.enthalpy(t, p, moles, phase)?;
        let s = self.entropy(t, p, moles, phase)?;
        Ok(h - t * s)
    }

    /// Extensive total Gibbs energy `G` `[J]` with exact `(T, P, n)` derivatives.
    pub fn gibbs_derivatives(
        &self,
        t: f64,
        p: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<ExtCspDerivatives> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let mut res = self.residual_gibbs_derivatives(t, p, moles, phase)?;
        let mut s_id = 0.0;
        let mut h_id = 0.0;
        for (i, &ni) in moles.iter().enumerate() {
            let xi = ni / n_tot;
            let ideal = self.shape_eos.ideal_gas(i)?.properties(t, p)?;
            let s_bar_i = ideal.entropy - R * xi.max(1.0e-100).ln();
            if ni > 0.0 {
                h_id += ni * ideal.enthalpy;
                s_id += ni * (ideal.entropy - R * xi.ln());
            }
            res.dn[i] += ideal.enthalpy - t * s_bar_i;
        }
        res.value += h_id - t * s_id;
        res.dt -= s_id;
        res.dp += n_tot * R * t / p;
        Ok(res)
    }

    /// Extensive caloric state `(P, U, H, S, A, G, Cv, Cp)` at fixed `(T, V, n)`.
    pub fn caloric_tv(&self, t: f64, volume: f64, moles: &[f64]) -> Result<ExtCspCaloricTv> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let f = self.residual_helmholtz_tv(t, volume, moles)?;
        let rt = R * t;
        let pressure = rt * (n_tot / volume - f.f_v);
        let dp_dv = -rt * (f.f_vv + n_tot / (volume * volume));
        let dp_dt = pressure / t - rt * f.f_tv;

        let u_res = -rt * t * f.f_t;
        let s_res = -R * (f.f + t * f.f_t);
        let a_res = rt * f.f;
        let cv_res = -rt * (2.0 * f.f_t + t * f.f_tt);

        let p_ref = 1.0e5;
        let mut u_id = 0.0;
        let mut s_id = 0.0;
        let mut cv_id = 0.0;
        for (i, &ni) in moles.iter().enumerate() {
            if ni == 0.0 {
                continue;
            }
            let ideal = self.shape_eos.ideal_gas(i)?.properties(t, p_ref)?;
            u_id += ni * (ideal.enthalpy - rt);
            s_id += ni * (ideal.entropy - R * (ni * rt / (volume * p_ref)).ln());
            cv_id += ni * (ideal.enthalpy_dt - R);
        }
        let a_id = u_id - t * s_id;

        let internal_energy = u_id + u_res;
        let entropy = s_id + s_res;
        let helmholtz_energy = a_id + a_res;
        let enthalpy = internal_energy + pressure * volume;
        let gibbs_energy = helmholtz_energy + pressure * volume;
        let heat_capacity_v = cv_id + cv_res;
        let heat_capacity_p = heat_capacity_v - t * (dp_dt * dp_dt) / dp_dv;

        Ok(ExtCspCaloricTv {
            pressure,
            internal_energy,
            enthalpy,
            entropy,
            helmholtz_energy,
            gibbs_energy,
            heat_capacity_v,
            heat_capacity_p,
        })
    }

    /// Extensive heat capacities `(Cp, Cv)` `[J/K]` at `(T, P, n, phase)`.
    pub fn heat_capacities_tp(
        &self,
        t: f64,
        p: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<(f64, f64)> {
        let v = self.volume_tp(t, p, moles, phase)?;
        let cal = self.caloric_tv(t, v, moles)?;
        Ok((cal.heat_capacity_p, cal.heat_capacity_v))
    }

    /// Single-phase speed of sound `[m/s]` at `(T, V, n)`.
    pub fn speed_of_sound_tv(&self, t: f64, volume: f64, moles: &[f64]) -> Result<f64> {
        let (n_tot, _) = self.validate_moles(moles)?;
        let f = self.residual_helmholtz_tv(t, volume, moles)?;
        let rt = R * t;
        let pressure = rt * (n_tot / volume - f.f_v);
        let dp_dv = -rt * (f.f_vv + n_tot / (volume * volume));
        let dp_dt = pressure / t - rt * f.f_tv;
        let cv_res = -rt * (2.0 * f.f_t + t * f.f_tt);

        let mut cv_id = 0.0;
        let mut mass = 0.0;
        for (i, &ni) in moles.iter().enumerate() {
            if ni == 0.0 {
                continue;
            }
            let ideal = self.shape_eos.ideal_gas(i)?.properties(t, 1.0e5)?;
            cv_id += ni * (ideal.enthalpy_dt - R);
            mass += ni * self.shape_eos.components()[i].molar_mass * 1.0e-3;
        }
        let cv = cv_id + cv_res;
        let dp_dv_s = dp_dv - t * (dp_dt * dp_dt) / cv;
        let w2 = -volume * volume * dp_dv_s / mass;
        if !w2.is_finite() || w2 <= 0.0 {
            return Err(Error::NoPhysicalRoot);
        }
        Ok(w2.sqrt())
    }

    /// Single-phase speed of sound `[m/s]` at `(T, P, n, phase)`.
    pub fn speed_of_sound_tp(&self, t: f64, p: f64, moles: &[f64], phase: Phase) -> Result<f64> {
        let v = self.volume_tp(t, p, moles, phase)?;
        self.speed_of_sound_tv(t, v, moles)
    }
}
