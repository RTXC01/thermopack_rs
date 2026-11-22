//! PC-SAFT, sPC-SAFT, PCP-SAFT, and sPCP-SAFT equations of state ported from
//! ThermoPack.
//!
//! Implements:
//! - Original PC-SAFT (Gross & Sadowski 2001, doi: 10.1021/ie0003887)
//! - Simplified sPC-SAFT (von Solms et al. 2003, doi: 10.1021/ie020753p)
//! - Polar PCP-SAFT / sPCP-SAFT with Gross (2005) quadrupole-quadrupole,
//!   Gross & Vrabec (2006) dipole-dipole, and Vrabec & Gross (2008)
//!   dipole-quadrupole contributions
//! - Wertheim TPT1 association (Michelsen & Hendriks 2001) via the shared
//!   [`crate::association`] module with exact first and second `(T, V, n)`
//!   derivatives
//! - Full 10-block `(T, V, n)` residual Helmholtz derivatives evaluated in
//!   unit-mole `(v = V/N, x = n/N)` space with exact Euler homogeneity scaling
//!   so all evaluations remain total-mole-scaling invariant across extreme
//!   mole-number scales.

use crate::association::{
    self, AssocCombRule, AssocRdfFirst, AssocRdfSecond, AssocRdfValue, AssocScheme,
    AssociationHelmholtz, AssociationParams, AssociationRdfProvider, DeltaCombRule,
};
use crate::bmcsl::{self, HsDiameters, PackingFractions};
use crate::error::positive;
use crate::ideal::IdealGas;
use crate::multipol::{MultipolHelmholtz, MultipolParams};
use crate::pc_saft_data::{select_pc_saft_kij, select_pc_saft_kij_record, select_pc_saft_record};
use crate::rdf::RdfDerivatives;
use crate::{Error, Phase, Result};

/// Avogadro's constant [1/mol] matching ThermoPack's `thermopack_constants.f90`.
pub const N_AVOGADRO: f64 = 6.022_140_76e23;
/// Universal gas constant [J/(mol K)] matching ThermoPack's `thermopack_constants.f90`.
pub const RGAS: f64 = 8.314_462_618_153_24;

/// PC-SAFT variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PcSaftKind {
    /// Original PC-SAFT with BMCSL hard-sphere mixing.
    PcSaft,
    /// Simplified PC-SAFT with Carnahan-Starling average packing fraction.
    SPcSaft,
    /// Original PCP-SAFT with BMCSL hard-sphere mixing and multipolar terms.
    PcpSaft,
    /// Simplified sPCP-SAFT with Carnahan-Starling mixing and multipolar terms.
    SPcpSaft,
}

impl PcSaftKind {
    /// Returns `true` for the simplified variants (`SPcSaft`, `SPcpSaft`).
    pub const fn is_simplified(self) -> bool {
        matches!(self, Self::SPcSaft | Self::SPcpSaft)
    }

    /// Returns `true` for the polar variants (`PcpSaft`, `SPcpSaft`).
    pub const fn is_polar(self) -> bool {
        matches!(self, Self::PcpSaft | Self::SPcpSaft)
    }

    /// ThermoPack EOS label string.
    pub const fn label(self) -> &'static str {
        match self {
            Self::PcSaft => "PC-SAFT",
            Self::SPcSaft => "sPC-SAFT",
            Self::PcpSaft => "PCP-SAFT",
            Self::SPcpSaft => "sPCP-SAFT",
        }
    }
}

/// Construction options for [`PcSaft::with_options`].
#[derive(Debug, Clone, Copy)]
pub struct PcSaftOptions<'a> {
    /// PC-SAFT variant (`PcSaft`, `SPcSaft`, `PcpSaft`, `SPcpSaft`).
    pub kind: PcSaftKind,
    /// Parameter-database reference string (e.g. `"Default"`, `"Gross2005/Gross2006"`).
    pub parameter_reference: &'a str,
    /// Whether to include the Wertheim TPT1 association term (default `true`).
    pub enable_association: bool,
    /// Optional explicit override for polar terms (`None` follows `kind.is_polar()`).
    pub enable_polar: Option<bool>,
    /// Combining rule for the cross-association strength $\Delta_{kl}$ (default `Standard`).
    pub delta_comb_rule: DeltaCombRule,
}

impl Default for PcSaftOptions<'static> {
    fn default() -> Self {
        Self {
            kind: PcSaftKind::PcSaft,
            parameter_reference: "Default",
            enable_association: true,
            enable_polar: None,
            delta_comb_rule: DeltaCombRule::Standard,
        }
    }
}

/// Universal dispersion polynomial coefficients from Gross & Sadowski 2001.
/// Dimension: `a_mat(row: 0..2, col: 0..6)`.
const A_MAT: [[f64; 7]; 3] = [
    [
        0.9105631445,
        0.6361281449,
        2.6861347891,
        -26.547362491,
        97.759208784,
        -159.59154087,
        91.297774084,
    ],
    [
        -0.3084016918,
        0.1860531159,
        -2.5030047259,
        21.419793629,
        -65.255885330,
        83.318680481,
        -33.746922930,
    ],
    [
        -0.0906148351,
        0.4527842806,
        0.5962700728,
        -1.7241829131,
        -4.1302112531,
        13.776631870,
        -8.6728470368,
    ],
];

/// Universal dispersion polynomial coefficients from Gross & Sadowski 2001.
/// Dimension: `b_mat(row: 0..2, col: 0..6)`.
const B_MAT: [[f64; 7]; 3] = [
    [
        0.7240946941,
        2.2382791861,
        -4.0025849485,
        -21.003576815,
        26.855641363,
        206.55133841,
        -355.60235612,
    ],
    [
        -0.5755498075,
        0.6995095521,
        3.8925673390,
        -17.215471648,
        192.67226447,
        -161.82646165,
        -165.20769346,
    ],
    [
        0.0976883116,
        -0.2557574982,
        -9.1558561530,
        20.642075974,
        -38.804430052,
        93.626774077,
        -29.666905585,
    ],
];

/// Fluid component parameters for PC-SAFT / sPC-SAFT / PCP-SAFT / sPCP-SAFT.
#[derive(Debug, Clone, Copy)]
pub struct PcSaftComponent {
    pub ident: &'static str,
    pub m: f64,
    pub sigma: f64,
    pub eps_depth_divk: f64,
    /// Association energy $\epsilon^{\text{assoc}}/k_B$ `[K]` (0 when non-associating).
    pub eps: f64,
    /// Dimensionless association volume $\beta^{\text{assoc}}$ (0 when non-associating).
    pub beta: f64,
    /// Association site scheme.
    pub assoc_scheme: AssocScheme,
    /// Dipole moment $\mu$ `[Debye]`.
    pub mu: f64,
    /// Quadrupole moment $Q$ `[Debye·Å]`.
    pub q: f64,
    pub reference: &'static str,
    pub molar_mass: f64,
}

/// Complete 10-block reduced residual Helmholtz energy $F = A^{\text{res}}/(RT)$
/// `[mol]` and all first and second `(T, V, n)` derivatives.
#[derive(Debug, Clone)]
pub struct PcSaftHelmholtz {
    pub f: f64,
    pub f_t: f64,
    pub f_v: f64,
    pub f_n: Vec<f64>,
    pub f_tt: f64,
    pub f_tv: f64,
    pub f_tn: Vec<f64>,
    pub f_vv: f64,
    pub f_vn: Vec<f64>,
    /// Row-major `nc * nc` second mole derivatives $F_{n_i n_j}$.
    pub f_nn: Vec<f64>,
    /// Converged non-bonded site fractions $X_k$ (empty when non-associating).
    pub x_k: Vec<f64>,
}

impl PcSaftHelmholtz {
    fn zeros(nc: usize) -> Self {
        Self {
            f: 0.0,
            f_t: 0.0,
            f_v: 0.0,
            f_n: vec![0.0; nc],
            f_tt: 0.0,
            f_tv: 0.0,
            f_tn: vec![0.0; nc],
            f_vv: 0.0,
            f_vn: vec![0.0; nc],
            f_nn: vec![0.0; nc * nc],
            x_k: Vec::new(),
        }
    }

    fn add_assign_parts(&mut self, other: &Self) {
        self.f += other.f;
        self.f_t += other.f_t;
        self.f_v += other.f_v;
        self.f_tt += other.f_tt;
        self.f_tv += other.f_tv;
        self.f_vv += other.f_vv;
        for (a, b) in self.f_n.iter_mut().zip(&other.f_n) {
            *a += *b;
        }
        for (a, b) in self.f_tn.iter_mut().zip(&other.f_tn) {
            *a += *b;
        }
        for (a, b) in self.f_vn.iter_mut().zip(&other.f_vn) {
            *a += *b;
        }
        for (a, b) in self.f_nn.iter_mut().zip(&other.f_nn) {
            *a += *b;
        }
    }

    fn add_assign_assoc(&mut self, assoc: &AssociationHelmholtz) {
        self.f += assoc.f;
        self.f_t += assoc.f_t;
        self.f_v += assoc.f_v;
        self.f_tt += assoc.f_tt;
        self.f_tv += assoc.f_tv;
        self.f_vv += assoc.f_vv;
        for (a, b) in self.f_n.iter_mut().zip(&assoc.f_n) {
            *a += *b;
        }
        for (a, b) in self.f_tn.iter_mut().zip(&assoc.f_tn) {
            *a += *b;
        }
        for (a, b) in self.f_vn.iter_mut().zip(&assoc.f_vn) {
            *a += *b;
        }
        for (a, b) in self.f_nn.iter_mut().zip(&assoc.f_nn) {
            *a += *b;
        }
        self.x_k = assoc.x_k.clone();
    }

    fn add_assign_multipol(&mut self, mpol: &MultipolHelmholtz) {
        self.f += mpol.f;
        self.f_t += mpol.f_t;
        self.f_v += mpol.f_v;
        self.f_tt += mpol.f_tt;
        self.f_tv += mpol.f_tv;
        self.f_vv += mpol.f_vv;
        for (a, b) in self.f_n.iter_mut().zip(&mpol.f_n) {
            *a += *b;
        }
        for (a, b) in self.f_tn.iter_mut().zip(&mpol.f_tn) {
            *a += *b;
        }
        for (a, b) in self.f_vn.iter_mut().zip(&mpol.f_vn) {
            *a += *b;
        }
        for (a, b) in self.f_nn.iter_mut().zip(&mpol.f_nn) {
            *a += *b;
        }
    }

    fn scale_from_unit(&mut self, total: f64) -> Result<()> {
        let inv_n = 1.0 / total;
        self.f *= total;
        self.f_t *= total;
        self.f_tt *= total;
        self.f_vv *= inv_n;
        for v in &mut self.f_vn {
            *v *= inv_n;
        }
        for v in &mut self.f_nn {
            *v *= inv_n;
        }
        let ok = self.f.is_finite()
            && self.f_t.is_finite()
            && self.f_v.is_finite()
            && self.f_tt.is_finite()
            && self.f_tv.is_finite()
            && self.f_vv.is_finite()
            && self.f_n.iter().all(|v| v.is_finite())
            && self.f_tn.iter().all(|v| v.is_finite())
            && self.f_vn.iter().all(|v| v.is_finite())
            && self.f_nn.iter().all(|v| v.is_finite())
            && self.x_k.iter().all(|v| v.is_finite());
        if !ok {
            return Err(Error::InvalidInput(
                "non-finite PC-SAFT Helmholtz energy or derivative",
            ));
        }
        Ok(())
    }
}

/// Scalar TP thermodynamic property with first `(T, P, n)` derivatives.
#[derive(Debug, Clone)]
pub struct PcSaftDerivatives {
    pub value: f64,
    pub dt: f64,
    pub dp: f64,
    pub dn: Vec<f64>,
}

/// Fugacity coefficients `ln(phi_i)` with first `(T, P, n)` derivatives.
#[derive(Debug, Clone)]
pub struct PcSaftFugacityDerivatives {
    pub ln_phi: Vec<f64>,
    pub dt: Vec<f64>,
    pub dp: Vec<f64>,
    /// Row-major `nc * nc` derivatives $\partial \ln \phi_i / \partial n_j$.
    pub dn: Vec<f64>,
}

/// Complete TP thermodynamic state with exact first derivatives.
#[derive(Debug, Clone)]
pub struct PcSaftTpState {
    pub volume: f64,
    pub compressibility_factor: PcSaftDerivatives,
    pub specific_volume: PcSaftDerivatives,
    pub residual_enthalpy: PcSaftDerivatives,
    pub enthalpy: PcSaftDerivatives,
    pub residual_entropy: PcSaftDerivatives,
    pub entropy: PcSaftDerivatives,
    pub residual_gibbs: PcSaftDerivatives,
    pub ln_fugacity: PcSaftFugacityDerivatives,
    pub cp: f64,
    pub cv: f64,
    pub speed_of_sound: f64,
}

/// PC-SAFT / sPC-SAFT / PCP-SAFT / sPCP-SAFT model instance.
#[derive(Debug, Clone)]
pub struct PcSaft {
    kind: PcSaftKind,
    components: Vec<PcSaftComponent>,
    kij: Vec<f64>, // nc * nc
    eps_comb_rules: Vec<AssocCombRule>,
    beta_comb_rules: Vec<AssocCombRule>,
    enable_association: bool,
    delta_comb_rule: DeltaCombRule,
    assoc: AssociationParams,
    multipol: Option<MultipolParams>,
    ideal: Vec<IdealGas>,
}

impl PcSaft {
    /// Create a new original PC-SAFT model with default parameters (including association).
    pub fn new(comps: &[&str]) -> Result<Self> {
        Self::new_with_options(PcSaftKind::PcSaft, comps, "Default")
    }

    /// Create a new simplified PC-SAFT (sPC-SAFT) model with default parameters.
    pub fn new_simplified(comps: &[&str]) -> Result<Self> {
        Self::new_with_options(PcSaftKind::SPcSaft, comps, "Default")
    }

    /// Create a new original polar PCP-SAFT model (using ThermoPack's
    /// `"Gross2005/Gross2006"` default reference override).
    pub fn new_polar(comps: &[&str]) -> Result<Self> {
        Self::new_with_options(PcSaftKind::PcpSaft, comps, "Default")
    }

    /// Create a new simplified polar sPCP-SAFT model (using ThermoPack's
    /// `"Gross2005/Gross2006"` default reference override).
    pub fn new_simplified_polar(comps: &[&str]) -> Result<Self> {
        Self::new_with_options(PcSaftKind::SPcpSaft, comps, "Default")
    }

    /// Create a new PC-SAFT model with specified variant and parameter reference.
    pub fn new_with_options(
        kind: PcSaftKind,
        comps: &[&str],
        parameter_reference: &str,
    ) -> Result<Self> {
        Self::with_options(
            comps,
            PcSaftOptions {
                kind,
                parameter_reference,
                enable_association: true,
                enable_polar: None,
                delta_comb_rule: DeltaCombRule::Standard,
            },
        )
    }

    /// Create a non-associating, non-polar PC-SAFT or sPC-SAFT model matching
    /// ThermoPack's `pc_saft_nonassoc` (`F_PC_SAFT_TVn` / `F_sPC_SAFT_TVn`) kernels.
    pub fn new_nonassoc(
        kind: PcSaftKind,
        comps: &[&str],
        parameter_reference: &str,
    ) -> Result<Self> {
        Self::with_options(
            comps,
            PcSaftOptions {
                kind,
                parameter_reference,
                enable_association: false,
                enable_polar: Some(false),
                delta_comb_rule: DeltaCombRule::Standard,
            },
        )
    }

    /// Create a PC-SAFT model with explicit [`PcSaftOptions`].
    pub fn with_options(comps: &[&str], options: PcSaftOptions<'_>) -> Result<Self> {
        if comps.is_empty() {
            return Err(Error::InvalidInput("components list cannot be empty"));
        }
        let nc = comps.len();
        let enable_polar = options.enable_polar.unwrap_or(options.kind.is_polar());

        // Match Fortran `eoslibinit::init_pcsaft` lines 983-985:
        // when polar is enabled and parameter_reference is "DEFAULT", override
        // to "Gross2005/Gross2006" so polar database entries are preferred.
        let effective_ref =
            if enable_polar && options.parameter_reference.eq_ignore_ascii_case("DEFAULT") {
                "Gross2005/Gross2006"
            } else {
                options.parameter_reference
            };

        let mut components = Vec::with_capacity(nc);
        let mut ideal = Vec::with_capacity(nc);
        const SMALL: f64 = 1.0e-12;

        for &ident in comps {
            let record = select_pc_saft_record(ident, effective_ref).ok_or_else(|| {
                Error::UnknownComponent(format!("PC-SAFT parameter not found for fluid: {ident}"))
            })?;
            let comp_data = crate::component::component(ident)?;
            // Match `pc_saft_parameters.f90:191-192`: fall back to `compdatadb`
            // dipole/quadrupole moments when the PC-SAFT record leaves them zero.
            let mu = if record.mu.abs() < SMALL {
                comp_data.mu_dipole
            } else {
                record.mu
            };
            let q = if record.q.abs() < SMALL {
                comp_data.q_quadrupole
            } else {
                record.q
            };
            components.push(PcSaftComponent {
                ident: record.ident,
                m: record.m,
                sigma: record.sigma,
                eps_depth_divk: record.eps_depth_divk,
                eps: record.eps,
                beta: record.beta,
                assoc_scheme: record.assoc_scheme,
                mu,
                q,
                reference: record.reference,
                molar_mass: comp_data.molar_mass,
            });
            ideal.push(IdealGas::new(comp_data)?);
        }

        let mut kij = vec![0.0; nc * nc];
        let mut eps_comb_rules = vec![AssocCombRule::DefaultComb; nc * nc];
        let mut beta_comb_rules = vec![AssocCombRule::DefaultComb; nc * nc];
        for i in 0..nc {
            for j in 0..nc {
                if i != j {
                    kij[i * nc + j] =
                        select_pc_saft_kij(components[i].ident, components[j].ident, effective_ref);
                    if let Some(rec) = select_pc_saft_kij_record(
                        components[i].ident,
                        components[j].ident,
                        effective_ref,
                    ) {
                        eps_comb_rules[i * nc + j] = rec.eps_comb_rule;
                        beta_comb_rules[i * nc + j] = rec.beta_comb_rule;
                    }
                }
            }
        }

        let assoc = build_association_params(
            &components,
            &eps_comb_rules,
            &beta_comb_rules,
            options.enable_association,
            options.delta_comb_rule,
        )?;
        let multipol = if enable_polar {
            Some(build_multipol_params(&components, &kij, true, true, true)?)
        } else {
            None
        };

        Ok(Self {
            kind: options.kind,
            components,
            kij,
            eps_comb_rules,
            beta_comb_rules,
            enable_association: options.enable_association,
            delta_comb_rule: options.delta_comb_rule,
            assoc,
            multipol,
            ideal,
        })
    }

    pub fn kind(&self) -> PcSaftKind {
        self.kind
    }

    pub fn components(&self) -> &[PcSaftComponent] {
        &self.components
    }

    pub fn num_components(&self) -> usize {
        self.components.len()
    }

    pub fn get_kij(&self, i: usize, j: usize) -> f64 {
        self.kij[i * self.components.len() + j]
    }

    pub fn set_kij(&mut self, i: usize, j: usize, val: f64) {
        let nc = self.components.len();
        self.kij[i * nc + j] = val;
        self.kij[j * nc + i] = val;
        if let Some(mp) = &self.multipol {
            let (dd, qq, dq) = (mp.enable_dd, mp.enable_qq, mp.enable_dq);
            if let Ok(rebuilt) = build_multipol_params(&self.components, &self.kij, dd, qq, dq) {
                self.multipol = Some(rebuilt);
            }
        }
    }

    /// Returns the active association parameter table.
    pub fn association(&self) -> &AssociationParams {
        &self.assoc
    }

    /// Returns `true` if this model has at least one active association site.
    pub fn has_association(&self) -> bool {
        self.assoc.has_sites()
    }

    /// Returns the multipolar parameter table when polar contributions are enabled.
    pub fn multipol(&self) -> Option<&MultipolParams> {
        self.multipol.as_ref()
    }

    /// Returns `true` if this model has any active multipolar (DD/QQ/DQ) interaction.
    pub fn has_polar(&self) -> bool {
        self.multipol.as_ref().is_some_and(|m| m.is_active())
    }

    /// Updates the pure-component association parameters for component `index`
    /// and rebuilds the association parameter table.
    pub fn set_assoc_params(
        &mut self,
        index: usize,
        eps: f64,
        beta: f64,
        scheme: AssocScheme,
    ) -> Result<()> {
        if index >= self.components.len() {
            return Err(Error::InvalidInput("component index out of bounds"));
        }
        if !eps.is_finite() || eps < 0.0 || !beta.is_finite() || beta < 0.0 {
            return Err(Error::InvalidInput(
                "association parameters must be finite and non-negative",
            ));
        }
        self.components[index].eps = eps;
        self.components[index].beta = beta;
        self.components[index].assoc_scheme = scheme;
        self.enable_association = true;
        self.assoc = build_association_params(
            &self.components,
            &self.eps_comb_rules,
            &self.beta_comb_rules,
            self.enable_association,
            self.delta_comb_rule,
        )?;
        Ok(())
    }

    /// Sets the cross-association combining rules for pair `(i, j)` and rebuilds
    /// the association parameter table.
    pub fn set_assoc_comb_rules(
        &mut self,
        i: usize,
        j: usize,
        eps_rule: AssocCombRule,
        beta_rule: AssocCombRule,
    ) -> Result<()> {
        let nc = self.components.len();
        if i >= nc || j >= nc {
            return Err(Error::InvalidInput("component index out of bounds"));
        }
        self.eps_comb_rules[i * nc + j] = eps_rule;
        self.eps_comb_rules[j * nc + i] = eps_rule;
        self.beta_comb_rules[i * nc + j] = beta_rule;
        self.beta_comb_rules[j * nc + i] = beta_rule;
        self.assoc = build_association_params(
            &self.components,
            &self.eps_comb_rules,
            &self.beta_comb_rules,
            self.enable_association,
            self.delta_comb_rule,
        )?;
        Ok(())
    }

    /// Sets the cross-association strength combining rule (`Cr1` or `Elliot`) and
    /// rebuilds the association parameter table.
    pub fn set_delta_comb_rule(&mut self, rule: DeltaCombRule) -> Result<()> {
        self.delta_comb_rule = rule;
        self.assoc = build_association_params(
            &self.components,
            &self.eps_comb_rules,
            &self.beta_comb_rules,
            self.enable_association,
            self.delta_comb_rule,
        )?;
        Ok(())
    }

    /// Updates dipole (`mu` `[D]`) and quadrupole (`q` `[D·Å]`) moments for component
    /// `index` and rebuilds the multipolar parameters if polar mode is active.
    pub fn set_polar_moments(&mut self, index: usize, mu: f64, q: f64) -> Result<()> {
        if index >= self.components.len() {
            return Err(Error::InvalidInput("component index out of bounds"));
        }
        if !mu.is_finite() || !q.is_finite() {
            return Err(Error::InvalidInput("polar moments must be finite"));
        }
        self.components[index].mu = mu;
        self.components[index].q = q;
        if let Some(mp) = &self.multipol {
            let (dd, qq, dq) = (mp.enable_dd, mp.enable_qq, mp.enable_dq);
            self.multipol = Some(build_multipol_params(
                &self.components,
                &self.kij,
                dd,
                qq,
                dq,
            )?);
        }
        Ok(())
    }

    /// Toggles individual multipolar sub-terms (dipole-dipole, quadrupole-quadrupole,
    /// dipole-quadrupole).
    pub fn set_multipol_control(&mut self, dd: bool, qq: bool, dq: bool) -> Result<()> {
        self.multipol = Some(build_multipol_params(
            &self.components,
            &self.kij,
            dd,
            qq,
            dq,
        )?);
        Ok(())
    }

    /// Prepare state for thermodynamic evaluations at temperature `T` and composition `x`.
    pub fn prepare<'a>(
        &'a self,
        temperature: f64,
        mole_fractions: &[f64],
    ) -> Result<PreparedPcSaft<'a>> {
        PreparedPcSaft::new(self, temperature, mole_fractions)
    }

    /// Evaluates the complete 10-block reduced residual Helmholtz energy
    /// $F(T, V, n) = A^{\text{res}}/(RT)$ `[mol]` and all first and second
    /// `(T, V, n)` derivatives, summing hard-chain, dispersion, association,
    /// and multipolar contributions.
    pub fn residual_helmholtz_tv(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<PcSaftHelmholtz> {
        let (total, v_unit, composition) =
            validate_tvn(self.num_components(), temperature, volume, moles)?;
        let prepared = PreparedPcSaft::new(self, temperature, &composition)?;
        let mut out = prepared.nonassoc_helmholtz_unit(v_unit)?;

        if self.assoc.has_sites() {
            let rdf_provider = PcSaftRdfProvider::new(&prepared, v_unit);
            let assoc = association::assoc_helmholtz_all(
                &self.assoc,
                &rdf_provider,
                temperature,
                v_unit,
                &composition,
            )?;
            out.add_assign_assoc(&assoc);
        }

        if let Some(mp) = &self.multipol
            && mp.is_active()
        {
            let mpol = mp.fres_all(temperature, v_unit, &composition)?;
            out.add_assign_multipol(&mpol);
        }

        out.scale_from_unit(total)?;
        Ok(out)
    }

    /// Evaluates only the non-associating, non-polar (hard-chain + dispersion)
    /// 10-block reduced residual Helmholtz energy $F(T, V, n)$ and derivatives.
    pub fn nonassoc_residual_helmholtz_tv(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<PcSaftHelmholtz> {
        let (total, v_unit, composition) =
            validate_tvn(self.num_components(), temperature, volume, moles)?;
        let prepared = PreparedPcSaft::new(self, temperature, &composition)?;
        let mut out = prepared.nonassoc_helmholtz_unit(v_unit)?;
        out.scale_from_unit(total)?;
        Ok(out)
    }

    /// Evaluates only the Wertheim TPT1 association 10-block reduced residual
    /// Helmholtz energy $F^{\text{assoc}}(T, V, n)$ and derivatives.
    pub fn assoc_residual_helmholtz_tv(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<AssociationHelmholtz> {
        let (_total, v_unit, composition) =
            validate_tvn(self.num_components(), temperature, volume, moles)?;
        if !self.assoc.has_sites() {
            let nc = self.num_components();
            return Ok(AssociationHelmholtz {
                f: 0.0,
                f_t: 0.0,
                f_v: 0.0,
                f_n: vec![0.0; nc],
                f_tt: 0.0,
                f_tv: 0.0,
                f_tn: vec![0.0; nc],
                f_vv: 0.0,
                f_vn: vec![0.0; nc],
                f_nn: vec![0.0; nc * nc],
                x_k: Vec::new(),
            });
        }
        let prepared = PreparedPcSaft::new(self, temperature, &composition)?;
        let rdf_provider = PcSaftRdfProvider::new(&prepared, v_unit);
        association::assoc_helmholtz_all(&self.assoc, &rdf_provider, temperature, volume, moles)
    }

    /// Evaluates only the Gross/Vrabec multipolar 10-block reduced residual
    /// Helmholtz energy $F^{\text{polar}}(T, V, n)$ and derivatives.
    pub fn multipol_residual_helmholtz_tv(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<MultipolHelmholtz> {
        let _ = validate_tvn(self.num_components(), temperature, volume, moles)?;
        match &self.multipol {
            Some(mp) if mp.is_active() => mp.fres_all(temperature, volume, moles),
            _ => {
                let nc = self.num_components();
                Ok(MultipolHelmholtz {
                    f: 0.0,
                    f_t: 0.0,
                    f_v: 0.0,
                    f_n: vec![0.0; nc],
                    f_tt: 0.0,
                    f_tv: 0.0,
                    f_tn: vec![0.0; nc],
                    f_vv: 0.0,
                    f_vn: vec![0.0; nc],
                    f_nn: vec![0.0; nc * nc],
                })
            }
        }
    }

    /// Solves for the non-bonded association site fractions $X_k$ at `(T, V, n)`.
    pub fn solve_assoc_x_k(
        &self,
        temperature: f64,
        volume: f64,
        moles: &[f64],
    ) -> Result<Vec<f64>> {
        let (_, v_unit, composition) =
            validate_tvn(self.num_components(), temperature, volume, moles)?;
        if !self.assoc.has_sites() {
            return Ok(Vec::new());
        }
        let prepared = PreparedPcSaft::new(self, temperature, &composition)?;
        let rdf_provider = PcSaftRdfProvider::new(&prepared, v_unit);
        association::solve_for_x_k(
            &self.assoc,
            &rdf_provider,
            temperature,
            v_unit,
            &composition,
            None,
        )
    }

    /// Solves for the molar volume [m^3/mol] at `(T, P, n)` for the requested `phase`.
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

    /// Evaluates the complete TP state with exact first `(T, P, n)` derivatives
    /// at `(T, P, n)` for the requested `phase`.
    pub fn state_tp(
        &self,
        temperature: f64,
        pressure: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<PcSaftTpState> {
        positive(temperature, "temperature")?;
        positive(pressure, "pressure")?;
        let nc = self.num_components();
        if moles.len() != nc {
            return Err(Error::InvalidInput(
                "moles length does not match components",
            ));
        }
        if moles.iter().any(|&v| !v.is_finite() || v < 0.0) {
            return Err(Error::InvalidInput(
                "mole numbers must be finite and non-negative",
            ));
        }
        let total: f64 = moles.iter().sum();
        if !total.is_finite() || total <= 0.0 {
            return Err(Error::InvalidInput(
                "total moles must be finite and strictly positive",
            ));
        }
        let x: Vec<f64> = moles.iter().map(|&v| v / total).collect();
        let prepared = PreparedPcSaft::new(self, temperature, &x)?;
        let v_molar = prepared.volume(pressure, phase)?;
        let v_total = v_molar * total;
        if !v_total.is_finite() || v_total <= 0.0 {
            return Err(Error::InvalidInput(
                "total volume overflowed or underflowed",
            ));
        }

        // Compute unit-mole 10-block Helmholtz energy at `(T, v_molar, x)` for
        // scale-invariant intensive derivatives, then apply exact `1/total`
        // scaling to extensive mole derivatives.
        let mut h = prepared.nonassoc_helmholtz_unit(v_molar)?;
        if self.assoc.has_sites() {
            let rdf_provider = PcSaftRdfProvider::new(&prepared, v_molar);
            let assoc = association::assoc_helmholtz_all(
                &self.assoc,
                &rdf_provider,
                temperature,
                v_molar,
                &x,
            )?;
            h.add_assign_assoc(&assoc);
        }
        if let Some(mp) = &self.multipol
            && mp.is_active()
        {
            let mpol = mp.fres_all(temperature, v_molar, &x)?;
            h.add_assign_multipol(&mpol);
        }

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

        // Residual enthalpy per mole and exact partial molar derivatives
        // matching ThermoPack's `enthalpy(..., residual=True, dhdn=True)`.
        let h_res = -RGAS * t * t * h.f_t + p_eval * v - RGAS * t;
        let cv_res = -RGAS * (2.0 * t * h.f_t + t * t * h.f_tt);
        let cp_res = cv_res - t * p_t * p_t / p_v - RGAS;
        let h_res_t = cp_res;
        let h_res_p = v - t * v_t;
        let mut h_res_dn = Vec::with_capacity(nc);
        for i in 0..nc {
            let h_i = -RGAS * t * t * h.f_tn[i] + t * p_t * v_n_partial[i] - RGAS * t;
            h_res_dn.push(h_i);
        }

        // Residual entropy (TP basis) per mole and partial molar derivatives
        // matching ThermoPack's `entropy(..., residual=True, dsdn=True)`.
        let s_res = -RGAS * (h.f + t * h.f_t - ln_z);
        let s_res_t = cp_res / t;
        let s_res_p = -v_t + RGAS / p_eval;
        let mut s_res_dn = Vec::with_capacity(nc);
        for i in 0..nc {
            let s_i = -RGAS * (h.f_n[i] + t * h.f_tn[i] - ln_z) + p_t * v_n_partial[i] - RGAS;
            s_res_dn.push(s_i);
        }

        // Residual Gibbs energy (TP basis) per mole and partial molar derivatives.
        let g_res = h_res - t * s_res;
        let g_res_t = -s_res;
        let g_res_p = v - RGAS * t / p_eval;
        let g_res_dn: Vec<f64> = h.f_n.iter().map(|&fni| RGAS * t * (fni - ln_z)).collect();

        // Ideal gas mixture properties at `(T, P, x)`.
        let mut h_id = 0.0;
        let mut cp_id = 0.0;
        let mut s_id = 0.0;
        let mut h_id_comp = Vec::with_capacity(nc);
        let mut s_id_comp = Vec::with_capacity(nc);
        for i in 0..nc {
            let id = self.ideal[i].properties(t, p_eval)?;
            h_id += x[i] * id.enthalpy;
            cp_id += x[i] * id.cp;
            let xi_ln = if x[i] > 0.0 { x[i].ln() } else { 0.0 };
            s_id += x[i] * (id.entropy - RGAS * xi_ln);
            h_id_comp.push(id.enthalpy);
            s_id_comp.push(id.entropy - RGAS * xi_ln);
        }
        let mut h_tot_dn = Vec::with_capacity(nc);
        let mut s_tot_dn = Vec::with_capacity(nc);
        for i in 0..nc {
            h_tot_dn.push(h_res_dn[i] + h_id_comp[i]);
            s_tot_dn.push(s_res_dn[i] + s_id_comp[i]);
        }

        // Fugacity coefficients and `(T, P, n)` derivatives.
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

        let all_finite = v.is_finite()
            && z.is_finite()
            && z_t.is_finite()
            && z_p.is_finite()
            && h_res.is_finite()
            && s_res.is_finite()
            && g_res.is_finite()
            && cp.is_finite()
            && cv.is_finite()
            && sound_speed.is_finite()
            && ln_phi.iter().all(|q| q.is_finite())
            && dlnphi_dt.iter().all(|q| q.is_finite())
            && dlnphi_dp.iter().all(|q| q.is_finite())
            && dlnphi_dn.iter().all(|q| q.is_finite());
        if !all_finite {
            return Err(Error::InvalidInput(
                "non-finite PC-SAFT TP thermodynamic state",
            ));
        }

        Ok(PcSaftTpState {
            volume: v,
            compressibility_factor: PcSaftDerivatives {
                value: z,
                dt: z_t,
                dp: z_p,
                dn: z_dn,
            },
            specific_volume: PcSaftDerivatives {
                value: v,
                dt: v_t,
                dp: v_p,
                dn: v_n_partial,
            },
            residual_enthalpy: PcSaftDerivatives {
                value: h_res,
                dt: h_res_t,
                dp: h_res_p,
                dn: h_res_dn,
            },
            enthalpy: PcSaftDerivatives {
                value: h_id + h_res,
                dt: cp,
                dp: h_res_p,
                dn: h_tot_dn,
            },
            residual_entropy: PcSaftDerivatives {
                value: s_res,
                dt: s_res_t,
                dp: s_res_p,
                dn: s_res_dn,
            },
            entropy: PcSaftDerivatives {
                value: s_id + s_res,
                dt: cp / t,
                dp: -v_t,
                dn: s_tot_dn,
            },
            residual_gibbs: PcSaftDerivatives {
                value: g_res,
                dt: g_res_t,
                dp: g_res_p,
                dn: g_res_dn,
            },
            ln_fugacity: PcSaftFugacityDerivatives {
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

    /// Evaluates fugacity coefficients and their `(T, P, n)` derivatives at `(T, P, n)`.
    pub fn ln_fugacity_derivatives_tp(
        &self,
        temperature: f64,
        pressure: f64,
        moles: &[f64],
        phase: Phase,
    ) -> Result<PcSaftFugacityDerivatives> {
        Ok(self
            .state_tp(temperature, pressure, moles, phase)?
            .ln_fugacity)
    }
}

fn validate_tvn(
    nc: usize,
    temperature: f64,
    volume: f64,
    moles: &[f64],
) -> Result<(f64, f64, Vec<f64>)> {
    positive(temperature, "temperature")?;
    positive(volume, "volume")?;
    if moles.len() != nc {
        return Err(Error::InvalidInput(
            "the number of moles must match the components of the model",
        ));
    }
    if moles.iter().any(|&v| !v.is_finite() || v < 0.0) {
        return Err(Error::InvalidInput(
            "mole numbers must be finite and non-negative",
        ));
    }
    let total: f64 = moles.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(Error::InvalidInput(
            "the number of moles must be strictly positive",
        ));
    }
    let v_unit = volume / total;
    if !v_unit.is_finite() || v_unit <= 0.0 {
        return Err(Error::InvalidInput(
            "molar volume must be finite and strictly positive",
        ));
    }
    let composition = moles.iter().map(|&v| v / total).collect();
    Ok((total, v_unit, composition))
}

fn build_association_params(
    components: &[PcSaftComponent],
    eps_comb_rules: &[AssocCombRule],
    beta_comb_rules: &[AssocCombRule],
    enable_association: bool,
    delta_comb_rule: DeltaCombRule,
) -> Result<AssociationParams> {
    let nc = components.len();
    if !enable_association {
        return Ok(AssociationParams::empty(nc));
    }
    let schemes: Vec<AssocScheme> = components.iter().map(|c| c.assoc_scheme).collect();
    let eps: Vec<f64> = components.iter().map(|c| c.eps).collect();
    let beta: Vec<f64> = components.iter().map(|c| c.beta).collect();
    let sigma: Vec<f64> = components.iter().map(|c| c.sigma).collect();
    AssociationParams::for_pc_saft(
        &schemes,
        &eps,
        &beta,
        &sigma,
        eps_comb_rules,
        beta_comb_rules,
        delta_comb_rule,
    )
}

fn build_multipol_params(
    components: &[PcSaftComponent],
    kij: &[f64],
    enable_dd: bool,
    enable_qq: bool,
    enable_dq: bool,
) -> Result<MultipolParams> {
    let nc = components.len();
    let m: Vec<f64> = components.iter().map(|c| c.m).collect();
    let mu: Vec<f64> = components.iter().map(|c| c.mu).collect();
    let q: Vec<f64> = components.iter().map(|c| c.q).collect();
    let mut sigma_ij = vec![0.0; nc * nc];
    let mut eps_divk_ij = vec![0.0; nc * nc];
    for i in 0..nc {
        for j in 0..nc {
            sigma_ij[i * nc + j] = 0.5 * (components[i].sigma + components[j].sigma);
            eps_divk_ij[i * nc + j] = (components[i].eps_depth_divk * components[j].eps_depth_divk)
                .sqrt()
                * (1.0 - kij[i * nc + j]);
        }
    }
    let mut mp = MultipolParams::new(&mu, &q, &m, &sigma_ij, &eps_divk_ij)?;
    mp.enable_dd = enable_dd;
    mp.enable_qq = enable_qq;
    mp.enable_dq = enable_dq;
    Ok(mp)
}

/// RDF provider wired to a [`PreparedPcSaft`] at fixed `(T, v_unit, x)`.
///
/// For original `PcSaft` / `PcpSaft`, precomputes the BMCSL `PackingFractions`
/// once so pair correlations `g_ij` are evaluated without recomputing packing
/// moments across all site pairs.
struct PcSaftRdfProvider<'a, 'b> {
    prepared: &'a PreparedPcSaft<'b>,
    v_unit: f64,
    diameters: HsDiameters,
    packing: Option<PackingFractions>,
}

impl<'a, 'b> PcSaftRdfProvider<'a, 'b> {
    fn new(prepared: &'a PreparedPcSaft<'b>, v_unit: f64) -> Self {
        let diameters = HsDiameters {
            d: prepared.d.clone(),
            d_t: prepared.d_t.clone(),
            d_tt: prepared.d_tt.clone(),
        };
        let packing = if prepared.model.kind.is_simplified() {
            None
        } else {
            let segments: Vec<f64> = prepared.model.components.iter().map(|c| c.m).collect();
            Some(bmcsl::packing_fractions(
                v_unit,
                &prepared.mole_fractions,
                &diameters,
                &segments,
            ))
        };
        Self {
            prepared,
            v_unit,
            diameters,
            packing,
        }
    }

    #[inline]
    fn covol(&self, i: usize, j: usize) -> f64 {
        let sig = 0.5
            * (self.prepared.model.components[i].sigma + self.prepared.model.components[j].sigma);
        N_AVOGADRO * sig * sig * sig
    }
}

impl AssociationRdfProvider for PcSaftRdfProvider<'_, '_> {
    fn pair_dependent(&self) -> bool {
        !self.prepared.model.kind.is_simplified()
    }

    fn rdf_value(&self, i: usize, j: usize) -> Result<AssocRdfValue> {
        let covol = self.covol(i, j);
        let g = if let Some(packing) = &self.packing {
            let mu = self.diameters.reduced_diameter(i, j);
            bmcsl::boublik_pair_correlation(
                packing.moments.zet[2],
                packing.moments.zet[3],
                mu.value,
            )
            .value
        } else {
            let eta = self.prepared.packing_fraction(1.0 / self.v_unit);
            let om = 1.0 - eta;
            (1.0 - 0.5 * eta) / (om * om * om)
        };
        Ok(AssocRdfValue { g, covol })
    }

    fn rdf_first(&self, i: usize, j: usize) -> Result<AssocRdfFirst> {
        let second = self.rdf_second(i, j)?;
        Ok(AssocRdfFirst {
            g: second.rdf.g,
            g_t: second.rdf.temperature,
            g_v: second.rdf.volume,
            g_n: second.rdf.moles,
            covol: second.covol,
            covol_t: second.covol_t,
        })
    }

    fn rdf_second(&self, i: usize, j: usize) -> Result<AssocRdfSecond> {
        let covol = self.covol(i, j);
        let rdf = if let Some(packing) = &self.packing {
            let count = self.prepared.model.num_components();
            let correlation = bmcsl::pair_correlation(packing, &self.diameters, i, j);
            let mut moles_second = Vec::with_capacity(count * count);
            for row in &correlation.d_nn[..count] {
                moles_second.extend_from_slice(&row[..count]);
            }
            RdfDerivatives {
                g: correlation.value,
                temperature: correlation.d_t,
                volume: correlation.d_v,
                moles: correlation.d_n,
                temperature_second: correlation.d_tt,
                temperature_volume: correlation.d_tv,
                temperature_moles: correlation.d_tn,
                volume_second: correlation.d_vv,
                volume_moles: correlation.d_vn,
                moles_second,
            }
        } else {
            crate::rdf::spcsaft_rdf_unit(self.prepared, self.v_unit)
        };
        Ok(AssocRdfSecond {
            rdf,
            covol,
            covol_t: 0.0,
            covol_tt: 0.0,
        })
    }
}

/// The reduced Helmholtz energy of the hard-sphere model of the PC-SAFT model
/// and its derivatives with respect to the temperature, the volume and the
/// mole numbers.
///
/// The function reproduces the original
/// `src/hardsphere_bmcsl.f90:calc_bmcsl_helmholtzenergy`, which the original
/// applies to the hard-sphere contribution of the PC-SAFT model.
pub fn pcsaft_hs_energy(
    model: &PcSaft,
    temperature: f64,
    volume: f64,
    moles: &[f64],
) -> Result<bmcsl::HardSphereEnergyDerivatives> {
    let (_, _, composition) = validate_tvn(model.num_components(), temperature, volume, moles)?;
    let prepared = PreparedPcSaft::new(model, temperature, &composition)?;
    let diameters = HsDiameters {
        d: prepared.d.clone(),
        d_t: prepared.d_t.clone(),
        d_tt: prepared.d_tt.clone(),
    };
    let segments: Vec<f64> = model
        .components()
        .iter()
        .map(|component| component.m)
        .collect();
    let packing = bmcsl::packing_fractions(volume, moles, &diameters, &segments);
    Ok(bmcsl::hard_sphere_energy(&packing))
}

/// The derivatives of the hard-chain contribution of the PC-SAFT model in the
/// volume representation of the original.
#[derive(Debug, Clone, Copy)]
struct HardChainDerivatives {
    d_t: f64,
    d_tt: f64,
    d_v: f64,
    d_vv: f64,
    d_tv: f64,
}

/// The reduced Helmholtz energy of the chain term of the PC-SAFT model and
/// its derivatives with respect to the temperature, the volume and the mole
/// numbers.
///
/// The function reproduces the original
/// `src/pc_saft_nonassoc.f90:F_chain_PC_SAFT_TVn`.
#[derive(Debug, Clone)]
pub struct PcSaftChainEnergy {
    /// The reduced Helmholtz energy of the chain term.
    pub value: f64,
    /// The temperature derivative of the energy.
    pub d_t: f64,
    /// The volume derivative of the energy.
    pub d_v: f64,
    /// The mole number derivatives of the energy.
    pub d_n: Vec<f64>,
    /// The second temperature derivative of the energy.
    pub d_tt: f64,
    /// The temperature and volume derivative of the energy.
    pub d_tv: f64,
    /// The temperature and mole number derivatives of the energy.
    pub d_tn: Vec<f64>,
    /// The second volume derivative of the energy.
    pub d_vv: f64,
    /// The volume and mole number derivatives of the energy.
    pub d_vn: Vec<f64>,
    /// The second mole number derivatives of the energy.
    pub d_nn: Vec<Vec<f64>>,
}

/// The reduced Helmholtz energy of the chain term of the PC-SAFT model.
///
/// The volume is in cubic metre and the number of moles is the amount of each
/// component.
pub fn pcsaft_chain_energy(
    model: &PcSaft,
    temperature: f64,
    volume: f64,
    moles: &[f64],
) -> Result<PcSaftChainEnergy> {
    let (_, _, composition) = validate_tvn(model.num_components(), temperature, volume, moles)?;
    let count = model.num_components();
    let prepared = PreparedPcSaft::new(model, temperature, &composition)?;
    let diameters = HsDiameters {
        d: prepared.d.clone(),
        d_t: prepared.d_t.clone(),
        d_tt: prepared.d_tt.clone(),
    };
    let segments: Vec<f64> = model
        .components()
        .iter()
        .map(|component| component.m)
        .collect();
    let packing = bmcsl::packing_fractions(volume, moles, &diameters, &segments);

    let mut energy = PcSaftChainEnergy {
        value: 0.0,
        d_t: 0.0,
        d_v: 0.0,
        d_n: vec![0.0; count],
        d_tt: 0.0,
        d_tv: 0.0,
        d_tn: vec![0.0; count],
        d_vv: 0.0,
        d_vn: vec![0.0; count],
        d_nn: vec![vec![0.0; count]; count],
    };

    for (component, amount) in moles.iter().enumerate() {
        let seg = model.components()[component].m;
        if seg <= 1.0 {
            continue;
        }
        let correlation = bmcsl::pair_correlation(&packing, &diameters, component, component);
        let lng = correlation.logarithm();
        let weight = amount * (seg - 1.0);

        energy.value -= weight * lng.value;
        energy.d_t -= weight * lng.d_t;
        energy.d_v -= weight * lng.d_v;
        energy.d_tt -= weight * lng.d_tt;
        energy.d_tv -= weight * lng.d_tv;
        energy.d_vv -= weight * lng.d_vv;

        energy.d_n[component] -= (seg - 1.0) * lng.value;
        energy.d_tn[component] -= (seg - 1.0) * lng.d_t;
        energy.d_vn[component] -= (seg - 1.0) * lng.d_v;
        for other in 0..count {
            energy.d_n[other] -= weight * lng.d_n[other];
            energy.d_tn[other] -= weight * lng.d_tn[other];
            energy.d_vn[other] -= weight * lng.d_vn[other];
            energy.d_nn[component][other] -= (seg - 1.0) * lng.d_n[other];
            energy.d_nn[other][component] -= (seg - 1.0) * lng.d_n[other];
        }
        for other in 0..count {
            for last in 0..count {
                energy.d_nn[other][last] -= weight * lng.d_nn[other][last];
            }
        }
    }

    Ok(energy)
}

/// Precomputed state of a PC-SAFT fluid/mixture at temperature $T$ and composition $\mathbf{x}$.
/// Construct with [`PcSaft::prepare`] or [`PreparedPcSaft::new`].
#[derive(Debug, Clone)]
pub struct PreparedPcSaft<'a> {
    pub model: &'a PcSaft,
    pub temperature: f64,
    pub mole_fractions: Vec<f64>,
    pub m_bar: f64,
    composition_sum: f64,
    pub m_bar_n: Vec<f64>,
    pub d: Vec<f64>,
    pub d_t: Vec<f64>,
    pub d_tt: Vec<f64>,
    pub a_disp: [f64; 7],
    pub b_disp: [f64; 7],
    pub a_disp_n: Vec<[f64; 7]>,
    pub b_disp_n: Vec<[f64; 7]>,
    pub m2e1s3_val: f64,
    pub m2e2s3_val: f64,
    pub m2e1s3_t: f64,
    pub m2e2s3_t: f64,
    pub m2e1s3_tt: f64,
    pub m2e2s3_tt: f64,
    pub m2e1s3_n: Vec<f64>,
    pub m2e2s3_n: Vec<f64>,
    pub eta_factor: f64,
    pub molar_mass_mix: f64,
}

impl<'a> PreparedPcSaft<'a> {
    pub fn new(model: &'a PcSaft, temperature: f64, mole_fractions: &[f64]) -> Result<Self> {
        positive(temperature, "temperature")?;
        let nc = model.components.len();
        if mole_fractions.len() != nc {
            return Err(Error::InvalidInput(
                "mole fractions length does not match components",
            ));
        }
        if mole_fractions.iter().any(|&v| !v.is_finite() || v < 0.0) {
            return Err(Error::InvalidInput(
                "mole fractions must be finite and non-negative",
            ));
        }
        let sumn: f64 = mole_fractions.iter().sum();
        if sumn <= 0.0 || !sumn.is_finite() {
            return Err(Error::InvalidInput(
                "mole fractions sum must be strictly positive",
            ));
        }
        let x: Vec<f64> = mole_fractions.iter().map(|&v| v / sumn).collect();

        // Temperature-dependent segment diameters d_i(T)
        let tinv = 1.0 / temperature;
        let tinv2 = tinv * tinv;
        let tinv3 = tinv2 * tinv;
        let mut d = Vec::with_capacity(nc);
        let mut d_t = Vec::with_capacity(nc);
        let mut d_tt = Vec::with_capacity(nc);

        for comp in &model.components {
            let eps = comp.eps_depth_divk;
            let sig_expo = comp.sigma * (-3.0 * eps * tinv).exp();
            let di = comp.sigma - 0.12 * sig_expo;
            let di_t = -0.36 * eps * tinv2 * sig_expo;
            let di_tt = sig_expo * (0.72 * eps - tinv * 1.08 * eps * eps) * tinv3;
            d.push(di);
            d_t.push(di_t);
            d_tt.push(di_tt);
        }

        // Mean segment number m_bar
        let mut m_bar = 0.0;
        let mut molar_mass_mix = 0.0;
        for (i, &xi) in x.iter().enumerate() {
            m_bar += xi * model.components[i].m;
            molar_mass_mix += xi * model.components[i].molar_mass;
        }

        let mut m_bar_n = Vec::with_capacity(nc);
        for comp in &model.components {
            m_bar_n.push((comp.m - m_bar) / sumn);
        }

        // Dispersion polynomial coefficients a(m_bar) and b(m_bar)
        let mbar_inv = 1.0 / m_bar;
        let mbar_inv2 = mbar_inv * mbar_inv;
        let mbar_inv3 = mbar_inv2 * mbar_inv;
        let h = (m_bar - 1.0) * mbar_inv;
        let g = h * (m_bar - 2.0) * mbar_inv;

        let mut a_disp = [0.0; 7];
        let mut b_disp = [0.0; 7];
        for k in 0..7 {
            a_disp[k] = A_MAT[0][k] + h * A_MAT[1][k] + g * A_MAT[2][k];
            b_disp[k] = B_MAT[0][k] + h * B_MAT[1][k] + g * B_MAT[2][k];
        }

        let mut a_disp_n = Vec::with_capacity(nc);
        let mut b_disp_n = Vec::with_capacity(nc);
        for &mb_n in &m_bar_n {
            let h_n = mb_n * mbar_inv2;
            let g_n = mb_n * (3.0 * m_bar - 4.0) * mbar_inv3;
            let mut a_n = [0.0; 7];
            let mut b_n = [0.0; 7];
            for k in 0..7 {
                a_n[k] = h_n * A_MAT[1][k] + g_n * A_MAT[2][k];
                b_n[k] = h_n * B_MAT[1][k] + g_n * B_MAT[2][k];
            }
            a_disp_n.push(a_n);
            b_disp_n.push(b_n);
        }

        // Cross-interaction sums for dispersion:
        // m2e1s3 = 1/(T * sumn^2) sum_{k,l} n_k n_l m_k m_l (eps_{kl}/k) sigma_{kl}^3
        // m2e2s3 = 1/(T^2 * sumn^2) sum_{k,l} n_k n_l m_k m_l (eps_{kl}/k)^2 sigma_{kl}^3
        let mut e1_sum = 0.0;
        let mut e2_sum = 0.0;
        let mut e1_matrix = vec![0.0; nc * nc];
        let mut e2_matrix = vec![0.0; nc * nc];

        for k in 0..nc {
            for l in 0..nc {
                let sigma_kl = 0.5 * (model.components[k].sigma + model.components[l].sigma);
                let sigma_kl3 = sigma_kl * sigma_kl * sigma_kl;
                let kij_val = model.get_kij(k, l);
                let eps_kl = (model.components[k].eps_depth_divk
                    * model.components[l].eps_depth_divk)
                    .sqrt()
                    * (1.0 - kij_val);
                let mk_ml = model.components[k].m * model.components[l].m;
                let e1_kl = mk_ml * eps_kl * sigma_kl3;
                let e2_kl = mk_ml * eps_kl * eps_kl * sigma_kl3;

                e1_matrix[k * nc + l] = e1_kl;
                e2_matrix[k * nc + l] = e2_kl;

                e1_sum += x[k] * x[l] * e1_kl;
                e2_sum += x[k] * x[l] * e2_kl;
            }
        }

        let t2 = temperature * temperature;
        let m2e1s3_val = e1_sum / temperature;
        let m2e2s3_val = e2_sum / t2;

        let m2e1s3_t = -m2e1s3_val / temperature;
        let m2e2s3_t = -2.0 * m2e2s3_val / temperature;

        let m2e1s3_tt = 2.0 * m2e1s3_val / t2;
        let m2e2s3_tt = 6.0 * m2e2s3_val / t2;

        let mut m2e1s3_n = Vec::with_capacity(nc);
        let mut m2e2s3_n = Vec::with_capacity(nc);
        for k in 0..nc {
            let mut dot1 = 0.0;
            let mut dot2 = 0.0;
            for l in 0..nc {
                dot1 += x[l] * (e1_matrix[l * nc + k] + e1_matrix[k * nc + l]);
                dot2 += x[l] * (e2_matrix[l * nc + k] + e2_matrix[k * nc + l]);
            }
            m2e1s3_n.push(dot1 / (temperature * sumn) - 2.0 * m2e1s3_val / sumn);
            m2e2s3_n.push(dot2 / (t2 * sumn) - 2.0 * m2e2s3_val / sumn);
        }

        // Packing fraction factor: eta = rho * eta_factor
        let mut eta_sum = 0.0;
        for i in 0..nc {
            let di = d[i];
            eta_sum += x[i] * model.components[i].m * di * di * di;
        }
        let eta_factor = (std::f64::consts::PI / 6.0) * eta_sum * N_AVOGADRO;

        Ok(Self {
            model,
            temperature,
            mole_fractions: x,
            m_bar,
            composition_sum: sumn,
            m_bar_n,
            d,
            d_t,
            d_tt,
            a_disp,
            b_disp,
            a_disp_n,
            b_disp_n,
            m2e1s3_val,
            m2e2s3_val,
            m2e1s3_t,
            m2e2s3_t,
            m2e1s3_tt,
            m2e2s3_tt,
            m2e1s3_n,
            m2e2s3_n,
            eta_factor,
            molar_mass_mix,
        })
    }

    /// Conversion numerator for volume solve: `conv_num = eta_factor` (unit-mole basis).
    pub fn conversion_numerator(&self) -> f64 {
        self.eta_factor
    }

    /// Evaluates packing fraction `eta` from molar density `rho` [mol/m^3].
    pub fn packing_fraction(&self, rho: f64) -> f64 {
        rho * self.eta_factor
    }

    /// Residual reduced Helmholtz free energy per mole `alpha_res(rho, T, x) = A^res / (n R T)`,
    /// including hard-chain, dispersion, association (when active), and multipolar (when active).
    pub fn residual_helmholtz_energy(&self, rho: f64) -> f64 {
        let eta = self.packing_fraction(rho);
        let mut alpha = self.alpha_hard_chain(rho, eta) + self.alpha_dispersion(rho, eta);
        let v_unit = 1.0 / rho;
        if self.model.assoc.has_sites() {
            let rdf = PcSaftRdfProvider::new(self, v_unit);
            if let Ok((f_assoc, _)) = association::assoc_helmholtz_value(
                &self.model.assoc,
                &rdf,
                self.temperature,
                v_unit,
                &self.mole_fractions,
            ) {
                alpha += f_assoc;
            } else {
                return f64::NAN;
            }
        }
        if let Some(mp) = &self.model.multipol
            && mp.is_active()
        {
            if let Ok(f_mp) = mp.fres_value(self.temperature, v_unit, &self.mole_fractions) {
                alpha += f_mp;
            } else {
                return f64::NAN;
            }
        }
        alpha
    }

    /// Compressibility factor `Z(rho, T, x) = 1 + rho * d(alpha_res)/d(rho)`.
    pub fn compressibility_factor(&self, rho: f64) -> f64 {
        let eta = self.packing_fraction(rho);
        let mut d_a_drho =
            self.d_alpha_hard_chain_drho(rho, eta) + self.d_alpha_dispersion_drho(rho, eta);
        let v_unit = 1.0 / rho;
        if self.model.assoc.has_sites() {
            let rdf = PcSaftRdfProvider::new(self, v_unit);
            if let Ok(jet) = association::assoc_helmholtz_pressure_jet(
                &self.model.assoc,
                &rdf,
                self.temperature,
                v_unit,
                &self.mole_fractions,
            ) {
                d_a_drho += -jet.f_v / (rho * rho);
            } else {
                return f64::NAN;
            }
        }
        if let Some(mp) = &self.model.multipol
            && mp.is_active()
        {
            if let Ok(jet) = mp.fres_pressure_jet(self.temperature, v_unit, &self.mole_fractions) {
                d_a_drho += -jet.f_v / (rho * rho);
            } else {
                return f64::NAN;
            }
        }
        1.0 + rho * d_a_drho
    }

    /// The second density derivative and the mixed density and temperature
    /// derivative of the reduced residual Helmholtz energy `(alpha_rhorho, alpha_rhot)`.
    pub fn residual_helmholtz_density_derivatives(&self, rho: f64) -> (f64, f64) {
        let eta = self.packing_fraction(rho);
        let (ahc_rhorho, ahc_rhot) = self.hard_chain_second_derivatives(rho, eta);
        let (_, _, adisp_rhorho, _, adisp_rhot) = self.dispersion_rho_t_derivatives(rho, eta);
        let mut a_rhorho = ahc_rhorho + adisp_rhorho;
        let mut a_rhot = ahc_rhot + adisp_rhot;
        let v_unit = 1.0 / rho;
        let rho2 = rho * rho;
        let rho3 = rho2 * rho;
        let rho4 = rho2 * rho2;

        if self.model.assoc.has_sites() {
            let rdf = PcSaftRdfProvider::new(self, v_unit);
            if let Ok(jet) = association::assoc_helmholtz_pressure_jet(
                &self.model.assoc,
                &rdf,
                self.temperature,
                v_unit,
                &self.mole_fractions,
            ) {
                a_rhorho += 2.0 * jet.f_v / rho3 + jet.f_vv / rho4;
                a_rhot += -jet.f_tv / rho2;
            } else {
                return (f64::NAN, f64::NAN);
            }
        }
        if let Some(mp) = &self.model.multipol
            && mp.is_active()
        {
            if let Ok(jet) = mp.fres_pressure_jet(self.temperature, v_unit, &self.mole_fractions) {
                a_rhorho += 2.0 * jet.f_v / rho3 + jet.f_vv / rho4;
                a_rhot += -jet.f_tv / rho2;
            } else {
                return (f64::NAN, f64::NAN);
            }
        }
        (a_rhorho, a_rhot)
    }

    /// Absolute pressure `P(rho, T, x)` `[Pa]`.
    pub fn pressure(&self, rho: f64) -> f64 {
        rho * RGAS * self.temperature * self.compressibility_factor(rho)
    }

    /// Pressure and its derivatives at density `rho` [mol/m^3].
    pub fn pressure_and_derivatives(&self, rho: f64) -> PcSaftPressureJet {
        let eta = self.packing_fraction(rho);
        let (ahc, ahc_rho) = self.hard_chain_energy(rho, eta);
        let (adisp, adisp_rho, adisp_rhorho, _, adisp_rhot) =
            self.dispersion_rho_t_derivatives(rho, eta);
        let (ahc_rhorho, ahc_rhot) = self.hard_chain_second_derivatives(rho, eta);

        let mut a_res = ahc + adisp;
        let mut a_rho = ahc_rho + adisp_rho;
        let mut a_rhorho = ahc_rhorho + adisp_rhorho;
        let mut a_rhot = ahc_rhot + adisp_rhot;

        let v_unit = 1.0 / rho;
        let rho2 = rho * rho;
        let rho3 = rho2 * rho;
        let rho4 = rho2 * rho2;

        if self.model.assoc.has_sites() {
            let rdf = PcSaftRdfProvider::new(self, v_unit);
            if let Ok(jet) = association::assoc_helmholtz_pressure_jet(
                &self.model.assoc,
                &rdf,
                self.temperature,
                v_unit,
                &self.mole_fractions,
            ) {
                a_res += jet.f;
                a_rho += -jet.f_v / rho2;
                a_rhorho += 2.0 * jet.f_v / rho3 + jet.f_vv / rho4;
                a_rhot += -jet.f_tv / rho2;
            } else {
                return PcSaftPressureJet {
                    pressure: f64::NAN,
                    compressibility: f64::NAN,
                    residual_helmholtz: f64::NAN,
                    dp_drho: f64::NAN,
                    dp_dv: f64::NAN,
                    dp_dt: f64::NAN,
                };
            }
        }
        if let Some(mp) = &self.model.multipol
            && mp.is_active()
        {
            if let Ok(jet) = mp.fres_pressure_jet(self.temperature, v_unit, &self.mole_fractions) {
                a_res += jet.f;
                a_rho += -jet.f_v / rho2;
                a_rhorho += 2.0 * jet.f_v / rho3 + jet.f_vv / rho4;
                a_rhot += -jet.f_tv / rho2;
            } else {
                return PcSaftPressureJet {
                    pressure: f64::NAN,
                    compressibility: f64::NAN,
                    residual_helmholtz: f64::NAN,
                    dp_drho: f64::NAN,
                    dp_dv: f64::NAN,
                    dp_dt: f64::NAN,
                };
            }
        }

        let z = 1.0 + rho * a_rho;
        let p = rho * RGAS * self.temperature * z;
        let dp_drho = RGAS * self.temperature * (1.0 + 2.0 * rho * a_rho + rho2 * a_rhorho);
        let dp_dv = -dp_drho * rho2;
        let dp_dt = rho * RGAS * (z + self.temperature * rho * a_rhot);

        PcSaftPressureJet {
            pressure: p,
            compressibility: z,
            residual_helmholtz: a_res,
            dp_drho,
            dp_dv,
            dp_dt,
        }
    }

    /// Complete state properties at density `rho` [mol/m^3].
    pub fn state(&self, rho: f64) -> Result<PcSaftState> {
        positive(rho, "density")?;
        let v_unit = 1.0 / rho;
        let eta = self.packing_fraction(rho);

        let (ahc, ahc_rho) = self.hard_chain_energy(rho, eta);
        let (adisp, adisp_rho, adisp_rhorho, adisp_t, adisp_rhot) =
            self.dispersion_rho_t_derivatives(rho, eta);
        let (ahc_rhorho, ahc_rhot) = self.hard_chain_second_derivatives(rho, eta);
        let (ahc_t, ahc_tt) = self.hard_chain_t_derivatives(rho, eta);
        let adisp_tt = self.dispersion_tt_derivative(rho, eta);

        let mut a_res = ahc + adisp;
        let mut a_rho = ahc_rho + adisp_rho;
        let mut a_rhorho = ahc_rhorho + adisp_rhorho;
        let mut a_rhot = ahc_rhot + adisp_rhot;
        let mut a_t = ahc_t + adisp_t;
        let mut a_tt = ahc_tt + adisp_tt;

        let rho2 = rho * rho;
        let rho3 = rho2 * rho;
        let rho4 = rho2 * rho2;

        let mut assoc_fn: Option<Vec<f64>> = None;
        if self.model.assoc.has_sites() {
            let rdf = PcSaftRdfProvider::new(self, v_unit);
            let assoc = association::assoc_helmholtz_all(
                &self.model.assoc,
                &rdf,
                self.temperature,
                v_unit,
                &self.mole_fractions,
            )?;
            a_res += assoc.f;
            a_rho += -assoc.f_v / rho2;
            a_rhorho += 2.0 * assoc.f_v / rho3 + assoc.f_vv / rho4;
            a_rhot += -assoc.f_tv / rho2;
            a_t += assoc.f_t;
            a_tt += assoc.f_tt;
            assoc_fn = Some(assoc.f_n);
        }

        let mut mpol_fn: Option<Vec<f64>> = None;
        if let Some(mp) = &self.model.multipol
            && mp.is_active()
        {
            let mpol = mp.fres_all(self.temperature, v_unit, &self.mole_fractions)?;
            a_res += mpol.f;
            a_rho += -mpol.f_v / rho2;
            a_rhorho += 2.0 * mpol.f_v / rho3 + mpol.f_vv / rho4;
            a_rhot += -mpol.f_tv / rho2;
            a_t += mpol.f_t;
            a_tt += mpol.f_tt;
            mpol_fn = Some(mpol.f_n);
        }

        let t = self.temperature;
        let z = 1.0 + rho * a_rho;
        let p = rho * RGAS * t * z;
        let dp_drho = RGAS * t * (1.0 + 2.0 * rho * a_rho + rho2 * a_rhorho);
        let dp_dv = -dp_drho * rho2;
        let dp_dt = rho * RGAS * (z + t * rho * a_rhot);

        let u_res = -RGAS * t * t * a_t;
        let h_res = u_res + RGAS * t * (z - 1.0);
        let s_res = -RGAS * a_res - RGAS * t * a_t;
        let a_res_energy = RGAS * t * a_res;
        let g_res = h_res - t * s_res;

        let cv_res = -RGAS * (2.0 * t * a_t + t * t * a_tt);
        let mut cp_id = 0.0;
        let mut h_id = 0.0;
        let mut s_id = 0.0;
        for i in 0..self.model.components.len() {
            let id_props = self.model.ideal[i].properties(t, 100_000.0)?;
            cp_id += self.mole_fractions[i] * id_props.cp;
            h_id += self.mole_fractions[i] * id_props.enthalpy;
            let xi = self.mole_fractions[i];
            let mix_term = if xi > 1e-15 { xi * xi.ln() } else { 0.0 };
            s_id += xi * id_props.entropy - RGAS * mix_term;
        }
        let p_ref = 100_000.0;
        s_id -= RGAS * (rho * RGAS * t / p_ref).ln();

        let cv_total = cp_id - RGAS + cv_res;
        let cp_total = cv_total - t * (dp_dt * dp_dt) / dp_dv;

        let gamma = cp_total / cv_total;
        let mw_kg = self.molar_mass_mix * 1e-3;
        let sound_speed = if dp_drho > 0.0 && gamma > 0.0 && mw_kg > 0.0 {
            (gamma * dp_drho / mw_kg).sqrt()
        } else {
            0.0
        };

        let ln_phi = if self.model.components.len() == 1 {
            vec![a_res + (z - 1.0) - z.ln()]
        } else {
            let mut fn_tot = self.nonassoc_helmholtz_fn_derivatives(rho);
            if let Some(afn) = assoc_fn {
                for (a, b) in fn_tot.iter_mut().zip(afn) {
                    *a += b;
                }
            }
            if let Some(mfn) = mpol_fn {
                for (a, b) in fn_tot.iter_mut().zip(mfn) {
                    *a += b;
                }
            }
            let ln_z = z.ln();
            for v in &mut fn_tot {
                *v -= ln_z;
            }
            fn_tot
        };

        Ok(PcSaftState {
            pressure: p,
            compressibility: z,
            residual_helmholtz: a_res,
            u_residual: u_res,
            h_residual: h_res,
            s_residual: s_res,
            a_residual: a_res_energy,
            g_residual: g_res,
            h_total: h_id + h_res,
            s_total: s_id + s_res,
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

    /// Calculate fugacity coefficients `ln(phi_i)` at density `rho` and compressibility `z`.
    pub fn fugacity_coefficients(&self, rho: f64, z: f64) -> Vec<f64> {
        let nc = self.model.components.len();
        if nc == 1 {
            let a_res = self.residual_helmholtz_energy(rho);
            return vec![a_res + (z - 1.0) - z.ln()];
        }

        let mut ln_phi = self.helmholtz_fn_derivatives(rho);
        let ln_z = z.ln();
        for value in &mut ln_phi {
            *value -= ln_z;
        }
        ln_phi
    }

    /// Single-phase specific volume root solver matching ThermoPack's
    /// `saft_master_volume_solver` (`saft_volume_solver` / `pc_saft_nonassoc_volume_solver`).
    pub fn volume(&self, pressure_target: f64, phase: Phase) -> Result<f64> {
        positive(pressure_target, "pressure")?;
        let conv_num = self.conversion_numerator();

        if phase == Phase::MinimumGibbs {
            let vl = self.volume(pressure_target, Phase::Liquid)?;
            let vv = self.volume(pressure_target, Phase::Vapor)?;
            if (vl - vv).abs() < 1e-10 * vl {
                return Ok(vl);
            }
            let jet_l = self.pressure_and_derivatives(1.0 / vl);
            let jet_v = self.pressure_and_derivatives(1.0 / vv);
            let g_tp_l =
                jet_l.residual_helmholtz + jet_l.compressibility - 1.0 - jet_l.compressibility.ln();
            let g_tp_v =
                jet_v.residual_helmholtz + jet_v.compressibility - 1.0 - jet_v.compressibility.ln();
            return Ok(if g_tp_l < g_tp_v { vl } else { vv });
        }

        let zeta_ideal = pressure_target * conv_num / (RGAS * self.temperature);
        let mut zeta_min = (1.0e-10_f64).min(0.1 * zeta_ideal).max(1.0e-30);
        let mut zeta_max = 0.75; // Gross & Sadowski maximum close-packing segment limit

        let mut zeta = match phase {
            Phase::Vapor => (1.0e-10_f64).min(zeta_ideal.max(1.0e-30)),
            Phase::Liquid => 0.5,
            Phase::MinimumGibbs => unreachable!(),
        };

        let mut v = conv_num / zeta;
        let mut jet = self.pressure_and_derivatives(1.0 / v);
        let mut f = (1.0 - zeta) * (jet.pressure - pressure_target);

        let max_iter = 80;
        let machine_prec = f64::EPSILON;

        for _ in 0..max_iter {
            // dF/dzeta = P_spec - P - ((1 - zeta) / zeta) * P_V * V
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
            v = conv_num / zeta;

            jet = self.pressure_and_derivatives(1.0 / v);
            f = (1.0 - zeta) * (jet.pressure - pressure_target);

            let v_converged = (v - v_old).abs() < 1e5 * machine_prec * v_old;
            let p_converged =
                (jet.pressure - pressure_target).abs() < 1e8 * machine_prec * pressure_target;

            if v_converged && p_converged && v.is_finite() {
                return Ok(v);
            }

            if f > 0.0 {
                zeta_max = zeta;
            } else if f < 0.0 {
                zeta_min = zeta;
            }
        }

        if v.is_finite()
            && (jet.pressure - pressure_target).abs() <= 1e-6 * pressure_target.max(1.0)
        {
            Ok(v)
        } else {
            Err(Error::NotConverged(
                "PC-SAFT volume solver did not converge",
            ))
        }
    }

    /// Pure fluid saturation pressure at temperature `T`.
    pub fn saturation_pressure(&self, temperature: f64) -> Result<PcSaftSaturation> {
        if self.model.components.len() != 1 {
            return Err(Error::Unsupported(
                "saturation_pressure is only supported for pure fluids",
            ));
        }

        let prep = self.model.prepare(temperature, &[1.0])?;

        let comp_data = crate::component::component(self.model.components[0].ident)?;
        let tc = comp_data.critical_temperature;
        let pc = comp_data.critical_pressure;
        let omega = comp_data.acentric_factor;

        let exp_arg: f64 = 5.373 * (1.0 + omega) * (1.0 - tc / temperature);
        let mut p_est = pc * exp_arg.exp();
        if p_est >= pc || p_est <= 0.0 || !p_est.is_finite() {
            p_est = 0.5 * pc;
        }

        let mut p = p_est;
        let max_iter = 60;

        for _ in 0..max_iter {
            let vl = match prep.volume(p, Phase::Liquid) {
                Ok(v) => v,
                Err(_) => {
                    p *= 0.8;
                    continue;
                }
            };
            let vv = match prep.volume(p, Phase::Vapor) {
                Ok(v) => v,
                Err(_) => {
                    p *= 1.2;
                    continue;
                }
            };

            if (vv - vl).abs() < 1e-12 * vl {
                return Err(Error::NotConverged(
                    "liquid and vapor volumes converged to the same critical root",
                ));
            }

            let state_l = prep.state(1.0 / vl)?;
            let state_v = prep.state(1.0 / vv)?;

            let delta_ln_phi = state_l.ln_phi[0] - state_v.ln_phi[0];
            if delta_ln_phi.abs() < 1e-10 {
                return Ok(PcSaftSaturation {
                    pressure: p,
                    temperature,
                    liquid_volume: vl,
                    vapor_volume: vv,
                    enthalpy_vaporization: state_v.h_total - state_l.h_total,
                    entropy_vaporization: state_v.s_total - state_l.s_total,
                });
            }

            let d_delta_dp = (vl - vv) / (RGAS * temperature);
            let mut dp = -delta_ln_phi / d_delta_dp;

            if dp > 0.5 * p {
                dp = 0.5 * p;
            } else if dp < -0.5 * p {
                dp = -0.5 * p;
            }

            p += dp;
            if p <= 0.0 {
                p = 0.1 * p_est;
            }
        }

        Err(Error::NotConverged(
            "saturation pressure solver did not converge",
        ))
    }

    // =========================================================================
    // FULL 10-BLOCK UNIT-MOLE (T, V, n) NON-ASSOCIATING HELMHOLTZ KERNEL
    // =========================================================================

    fn nonassoc_helmholtz_unit(&self, v_unit: f64) -> Result<PcSaftHelmholtz> {
        let mut out = if self.model.kind.is_simplified() {
            self.simplified_hard_chain_helmholtz_unit(v_unit)
        } else {
            self.original_hard_chain_helmholtz_unit(v_unit)
        };
        let disp = self.dispersion_helmholtz_unit(v_unit);
        out.add_assign_parts(&disp);
        Ok(out)
    }

    fn original_hard_chain_helmholtz_unit(&self, v_unit: f64) -> PcSaftHelmholtz {
        let nc = self.model.components.len();
        let x = &self.mole_fractions;
        let diameters = HsDiameters {
            d: self.d.clone(),
            d_t: self.d_t.clone(),
            d_tt: self.d_tt.clone(),
        };
        let segments: Vec<f64> = self.model.components.iter().map(|c| c.m).collect();
        let packing = bmcsl::packing_fractions(v_unit, x, &diameters, &segments);
        let hs = bmcsl::hard_sphere_energy(&packing);

        // F_hs = M * a_hs where M = sum_i x_i m_i = m_bar (at unit total moles).
        let m_tot = self.m_bar;
        let mut out = PcSaftHelmholtz::zeros(nc);
        out.f = m_tot * hs.value;
        out.f_t = m_tot * hs.d_t;
        out.f_v = m_tot * hs.d_v;
        out.f_tt = m_tot * hs.d_tt;
        out.f_tv = m_tot * hs.d_tv;
        out.f_vv = m_tot * hs.d_vv;
        for i in 0..nc {
            let mi = segments[i];
            out.f_n[i] = m_tot * hs.d_n[i] + mi * hs.value;
            out.f_tn[i] = m_tot * hs.d_tn[i] + mi * hs.d_t;
            out.f_vn[i] = m_tot * hs.d_vn[i] + mi * hs.d_v;
            for j in 0..nc {
                let mj = segments[j];
                out.f_nn[i * nc + j] = m_tot * hs.d_nn[i][j] + mi * hs.d_n[j] + mj * hs.d_n[i];
            }
        }

        // F_chain = -sum_{i: m_i > 1} x_i * (m_i - 1) * ln(g_ii)
        for i in 0..nc {
            let seg = segments[i];
            if seg <= 1.0 {
                continue;
            }
            let corr = bmcsl::pair_correlation(&packing, &diameters, i, i);
            let lng = corr.logarithm();
            let seg_m1 = seg - 1.0;
            let weight = x[i] * seg_m1;

            out.f -= weight * lng.value;
            out.f_t -= weight * lng.d_t;
            out.f_v -= weight * lng.d_v;
            out.f_tt -= weight * lng.d_tt;
            out.f_tv -= weight * lng.d_tv;
            out.f_vv -= weight * lng.d_vv;

            out.f_n[i] -= seg_m1 * lng.value;
            out.f_tn[i] -= seg_m1 * lng.d_t;
            out.f_vn[i] -= seg_m1 * lng.d_v;
            for j in 0..nc {
                out.f_n[j] -= weight * lng.d_n[j];
                out.f_tn[j] -= weight * lng.d_tn[j];
                out.f_vn[j] -= weight * lng.d_vn[j];
                out.f_nn[i * nc + j] -= seg_m1 * lng.d_n[j];
                out.f_nn[j * nc + i] -= seg_m1 * lng.d_n[j];
            }
            for j in 0..nc {
                for k in 0..nc {
                    out.f_nn[j * nc + k] -= weight * lng.d_nn[j][k];
                }
            }
        }

        out
    }

    fn simplified_hard_chain_helmholtz_unit(&self, v_unit: f64) -> PcSaftHelmholtz {
        let nc = self.model.components.len();
        let x = &self.mole_fractions;
        let inv_v = 1.0 / v_unit;
        let pref = (std::f64::consts::PI * N_AVOGADRO / 6.0) * inv_v;

        let mut eta = 0.0;
        let mut eta_t = 0.0;
        let mut eta_tt = 0.0;
        let mut eta_x = Vec::with_capacity(nc);
        let mut eta_tx = Vec::with_capacity(nc);
        let mut eta_vx = Vec::with_capacity(nc);
        let mut c_chain = 0.0;

        for i in 0..nc {
            let m = self.model.components[i].m;
            let d = self.d[i];
            let dt = self.d_t[i];
            let dtt = self.d_tt[i];
            let d2 = d * d;
            let d3 = d2 * d;
            let ex = pref * m * d3;
            let etx = pref * m * (3.0 * d2 * dt);
            let ettx = pref * m * (6.0 * d * dt * dt + 3.0 * d2 * dtt);
            eta += x[i] * ex;
            eta_t += x[i] * etx;
            eta_tt += x[i] * ettx;
            c_chain += x[i] * (m - 1.0);
            eta_x.push(ex);
            eta_tx.push(etx);
            eta_vx.push(-ex * inv_v);
        }

        let eta_v = -eta * inv_v;
        let eta_vv = 2.0 * eta * inv_v * inv_v;
        let eta_tv = -eta_t * inv_v;

        let om = 1.0 - eta;
        let om2 = om * om;
        let om3 = om2 * om;
        let om4 = om2 * om2;
        let om5 = om4 * om;

        let a_hs = (4.0 * eta - 3.0 * eta * eta) / om2;
        let a_hs_1 = 2.0 * (2.0 - eta) / om3;
        let a_hs_2 = 2.0 * (5.0 - 2.0 * eta) / om4;

        let g = (1.0 - 0.5 * eta) / om3;
        let g_1 = (5.0 - 2.0 * eta) / (2.0 * om4);
        let g_2 = 3.0 * (3.0 - eta) / om5;
        let ln_g = g.ln();
        let ln_g_1 = g_1 / g;
        let ln_g_2 = g_2 / g - ln_g_1 * ln_g_1;

        let m_tot = self.m_bar;
        let w1 = m_tot * a_hs_1 - c_chain * ln_g_1;
        let w2 = m_tot * a_hs_2 - c_chain * ln_g_2;

        let mut out = PcSaftHelmholtz::zeros(nc);
        out.f = m_tot * a_hs - c_chain * ln_g;
        out.f_t = w1 * eta_t;
        out.f_v = w1 * eta_v;
        out.f_tt = w2 * eta_t * eta_t + w1 * eta_tt;
        out.f_tv = w2 * eta_t * eta_v + w1 * eta_tv;
        out.f_vv = w2 * eta_v * eta_v + w1 * eta_vv;

        let mut coeff_1 = Vec::with_capacity(nc);
        for i in 0..nc {
            let mi = self.model.components[i].m;
            let c1i = mi * a_hs_1 - (mi - 1.0) * ln_g_1;
            coeff_1.push(c1i);
            out.f_n[i] = mi * a_hs - (mi - 1.0) * ln_g + w1 * eta_x[i];
            out.f_tn[i] = c1i * eta_t + w2 * eta_t * eta_x[i] + w1 * eta_tx[i];
            out.f_vn[i] = c1i * eta_v + w2 * eta_v * eta_x[i] + w1 * eta_vx[i];
        }
        for i in 0..nc {
            for j in 0..nc {
                out.f_nn[i * nc + j] =
                    coeff_1[i] * eta_x[j] + coeff_1[j] * eta_x[i] + w2 * eta_x[i] * eta_x[j];
            }
        }
        out
    }

    fn dispersion_helmholtz_unit(&self, v_unit: f64) -> PcSaftHelmholtz {
        let nc = self.model.components.len();
        let x = &self.mole_fractions;
        let t = self.temperature;
        let inv_t = 1.0 / t;
        let inv_t2 = inv_t * inv_t;
        let inv_v = 1.0 / v_unit;
        let pref_eta = (std::f64::consts::PI * N_AVOGADRO / 6.0) * inv_v;

        let mut eta = 0.0;
        let mut eta_t = 0.0;
        let mut eta_tt = 0.0;
        let mut eta_x = Vec::with_capacity(nc);
        let mut eta_tx = Vec::with_capacity(nc);
        let mut eta_vx = Vec::with_capacity(nc);

        for i in 0..nc {
            let m = self.model.components[i].m;
            let d = self.d[i];
            let dt = self.d_t[i];
            let dtt = self.d_tt[i];
            let d2 = d * d;
            let d3 = d2 * d;
            let ex = pref_eta * m * d3;
            let etx = pref_eta * m * (3.0 * d2 * dt);
            let ettx = pref_eta * m * (6.0 * d * dt * dt + 3.0 * d2 * dtt);
            eta += x[i] * ex;
            eta_t += x[i] * etx;
            eta_tt += x[i] * ettx;
            eta_x.push(ex);
            eta_tx.push(etx);
            eta_vx.push(-ex * inv_v);
        }

        let eta_v = -eta * inv_v;
        let eta_vv = 2.0 * eta * inv_v * inv_v;
        let eta_tv = -eta_t * inv_v;

        // m_bar and its unit-mole derivatives:
        // m_bar_x[i] = m_i - m_bar, m_bar_xx[i,j] = -(m_bar_x[i] + m_bar_x[j]).
        let mbar = self.m_bar;
        let mbar_x: Vec<f64> = self.model.components.iter().map(|c| c.m - mbar).collect();

        // Scalar derivatives of a_k(m_bar), b_k(m_bar) up to 2nd order.
        let mi = 1.0 / mbar;
        let mi2 = mi * mi;
        let mi3 = mi2 * mi;
        let mi4 = mi2 * mi2;
        let h_m = (mbar - 1.0) * mi;
        let h_m1 = mi2;
        let h_m2 = -2.0 * mi3;
        let g_m = h_m * (mbar - 2.0) * mi;
        let g_m1 = (3.0 * mbar - 4.0) * mi3;
        let g_m2 = (-6.0 * mbar + 12.0) * mi4;

        let mut ak = [0.0; 7];
        let mut ak_m = [0.0; 7];
        let mut ak_mm = [0.0; 7];
        let mut bk = [0.0; 7];
        let mut bk_m = [0.0; 7];
        let mut bk_mm = [0.0; 7];
        for k in 0..7 {
            ak[k] = A_MAT[0][k] + h_m * A_MAT[1][k] + g_m * A_MAT[2][k];
            ak_m[k] = h_m1 * A_MAT[1][k] + g_m1 * A_MAT[2][k];
            ak_mm[k] = h_m2 * A_MAT[1][k] + g_m2 * A_MAT[2][k];
            bk[k] = B_MAT[0][k] + h_m * B_MAT[1][k] + g_m * B_MAT[2][k];
            bk_m[k] = h_m1 * B_MAT[1][k] + g_m1 * B_MAT[2][k];
            bk_mm[k] = h_m2 * B_MAT[1][k] + g_m2 * B_MAT[2][k];
        }

        let mut i1 = 0.0;
        let mut i1_e = 0.0;
        let mut i1_ee = 0.0;
        let mut i1_m = 0.0;
        let mut i1_em = 0.0;
        let mut i1_mm = 0.0;
        let mut i2 = 0.0;
        let mut i2_e = 0.0;
        let mut i2_ee = 0.0;
        let mut i2_m = 0.0;
        let mut i2_em = 0.0;
        let mut i2_mm = 0.0;

        let mut powers = [1.0; 7];
        for k in 1..7 {
            powers[k] = powers[k - 1] * eta;
        }
        for k in 0..7 {
            let pk = powers[k];
            i1 += ak[k] * pk;
            i1_m += ak_m[k] * pk;
            i1_mm += ak_mm[k] * pk;
            i2 += bk[k] * pk;
            i2_m += bk_m[k] * pk;
            i2_mm += bk_mm[k] * pk;
            if k >= 1 {
                let kf = k as f64;
                let pk1 = powers[k - 1];
                i1_e += kf * ak[k] * pk1;
                i1_em += kf * ak_m[k] * pk1;
                i2_e += kf * bk[k] * pk1;
                i2_em += kf * bk_m[k] * pk1;
            }
            if k >= 2 {
                let k2f = (k * (k - 1)) as f64;
                let pk2 = powers[k - 2];
                i1_ee += k2f * ak[k] * pk2;
                i2_ee += k2f * bk[k] * pk2;
            }
        }

        // C1(eta, mbar) = 1 / (1 + mbar * o1(eta) + (1 - mbar) * o2(eta))
        let e2 = eta * eta;
        let e3 = e2 * eta;
        let e4 = e2 * e2;
        let om = 1.0 - eta;
        let om2 = om * om;
        let om4 = om2 * om2;
        let om5 = om4 * om;
        let om6 = om4 * om2;
        let q = om * (2.0 - eta);
        let q2 = q * q;
        let q3 = q2 * q;
        let q4 = q2 * q2;

        let o1 = (8.0 * eta - 2.0 * e2) / om4;
        let o1_e = (-4.0 * e2 + 20.0 * eta + 8.0) / om5;
        let o1_ee = ((-8.0 * eta + 20.0) * om + (-4.0 * e2 + 20.0 * eta + 8.0) * 5.0) / om6;

        let o2 = (20.0 * eta - 27.0 * e2 + 12.0 * e3 - 2.0 * e4) / q2;
        let o2_e = (2.0 * e3 + 12.0 * e2 - 48.0 * eta + 40.0) / q3;
        let o2_ee = ((6.0 * e2 + 24.0 * eta - 48.0) * q
            + (2.0 * e3 + 12.0 * e2 - 48.0 * eta + 40.0) * (3.0 - 2.0 * eta) * 3.0)
            / q4;

        let denom = 1.0 + mbar * o1 + (1.0 - mbar) * o2;
        let d_e = mbar * o1_e + (1.0 - mbar) * o2_e;
        let d_ee = mbar * o1_ee + (1.0 - mbar) * o2_ee;
        let d_m = o1 - o2;
        let d_em = o1_e - o2_e;

        let c1 = 1.0 / denom;
        let c1_2 = c1 * c1;
        let c1_3 = c1_2 * c1;
        let c1_e = -c1_2 * d_e;
        let c1_m = -c1_2 * d_m;
        let c1_ee = 2.0 * c1_3 * d_e * d_e - c1_2 * d_ee;
        let c1_em = 2.0 * c1_3 * d_e * d_m - c1_2 * d_em;
        let c1_mm = 2.0 * c1_3 * d_m * d_m;

        // G1(eta, mbar) = 2 * I1
        let g1 = [
            2.0 * i1,
            2.0 * i1_e,
            2.0 * i1_m,
            2.0 * i1_ee,
            2.0 * i1_em,
            2.0 * i1_mm,
        ];

        // P = C1 * I2, G2(eta, mbar) = mbar * P
        let p = c1 * i2;
        let p_e = c1_e * i2 + c1 * i2_e;
        let p_m = c1_m * i2 + c1 * i2_m;
        let p_ee = c1_ee * i2 + 2.0 * c1_e * i2_e + c1 * i2_ee;
        let p_em = c1_em * i2 + c1_e * i2_m + c1_m * i2_e + c1 * i2_em;
        let p_mm = c1_mm * i2 + 2.0 * c1_m * i2_m + c1 * i2_mm;
        let g2 = [
            mbar * p,
            mbar * p_e,
            p + mbar * p_m,
            mbar * p_ee,
            p_e + mbar * p_em,
            2.0 * p_m + mbar * p_mm,
        ];

        // Symmetric cross-interaction matrices E1, E2 and S1(T, x), S2(T, x)
        let mut e1_mat = vec![0.0; nc * nc];
        let mut e2_mat = vec![0.0; nc * nc];
        let mut q1 = 0.0;
        let mut q2_sum = 0.0;
        for i in 0..nc {
            for j in 0..nc {
                let sig = 0.5 * (self.model.components[i].sigma + self.model.components[j].sigma);
                let sig3 = sig * sig * sig;
                let eps = (self.model.components[i].eps_depth_divk
                    * self.model.components[j].eps_depth_divk)
                    .sqrt()
                    * (1.0 - self.model.get_kij(i, j));
                let mi_mj = self.model.components[i].m * self.model.components[j].m;
                let e1_ij = mi_mj * eps * sig3;
                let e2_ij = mi_mj * eps * eps * sig3;
                e1_mat[i * nc + j] = e1_ij;
                e2_mat[i * nc + j] = e2_ij;
                q1 += x[i] * x[j] * e1_ij;
                q2_sum += x[i] * x[j] * e2_ij;
            }
        }

        let s1 = q1 * inv_t;
        let s1_t = -s1 * inv_t;
        let s1_tt = 2.0 * s1 * inv_t2;
        let s2 = q2_sum * inv_t2;
        let s2_t = -2.0 * s2 * inv_t;
        let s2_tt = 6.0 * s2 * inv_t2;

        let mut s1_x = vec![0.0; nc];
        let mut s1_tx = vec![0.0; nc];
        let mut s2_x = vec![0.0; nc];
        let mut s2_tx = vec![0.0; nc];
        for i in 0..nc {
            let mut dot1 = 0.0;
            let mut dot2 = 0.0;
            for j in 0..nc {
                dot1 += x[j] * e1_mat[i * nc + j];
                dot2 += x[j] * e2_mat[i * nc + j];
            }
            s1_x[i] = 2.0 * dot1 * inv_t;
            s1_tx[i] = -s1_x[i] * inv_t;
            s2_x[i] = 2.0 * dot2 * inv_t2;
            s2_tx[i] = -2.0 * s2_x[i] * inv_t;
        }

        // Lift G1(eta, mbar) and G2(eta, mbar) into (T, v, x) derivatives,
        // accumulate H = G1 * S1 + G2 * S2.
        let mut h = PcSaftHelmholtz::zeros(nc);
        for (g_jet, s_val, s_t, s_tt, s_x, s_tx, e_mat, s_scale) in [
            (g1, s1, s1_t, s1_tt, &s1_x, &s1_tx, &e1_mat, 2.0 * inv_t),
            (g2, s2, s2_t, s2_tt, &s2_x, &s2_tx, &e2_mat, 2.0 * inv_t2),
        ] {
            let [gv, ge, gm, gee, gem, gmm] = g_jet;
            let gt = ge * eta_t;
            let gv_d = ge * eta_v;
            let gtt = gee * eta_t * eta_t + ge * eta_tt;
            let gtv = gee * eta_t * eta_v + ge * eta_tv;
            let gvv = gee * eta_v * eta_v + ge * eta_vv;

            h.f += gv * s_val;
            h.f_t += gt * s_val + gv * s_t;
            h.f_v += gv_d * s_val;
            h.f_tt += gtt * s_val + 2.0 * gt * s_t + gv * s_tt;
            h.f_tv += gtv * s_val + gv_d * s_t;
            h.f_vv += gvv * s_val;

            let mut gx = Vec::with_capacity(nc);
            for i in 0..nc {
                let gxi = ge * eta_x[i] + gm * mbar_x[i];
                let gtxi = (gee * eta_x[i] + gem * mbar_x[i]) * eta_t + ge * eta_tx[i];
                let gvxi = (gee * eta_x[i] + gem * mbar_x[i]) * eta_v + ge * eta_vx[i];
                gx.push(gxi);

                h.f_n[i] += gxi * s_val + gv * s_x[i];
                h.f_tn[i] += gtxi * s_val + gt * s_x[i] + gxi * s_t + gv * s_tx[i];
                h.f_vn[i] += gvxi * s_val + gv_d * s_x[i];
            }
            for i in 0..nc {
                for j in 0..nc {
                    let mbar_ij = -(mbar_x[i] + mbar_x[j]);
                    let gxij = gee * eta_x[i] * eta_x[j]
                        + gem * (eta_x[i] * mbar_x[j] + eta_x[j] * mbar_x[i])
                        + gmm * mbar_x[i] * mbar_x[j]
                        + gm * mbar_ij;
                    let s_ij = s_scale * e_mat[i * nc + j];
                    h.f_nn[i * nc + j] +=
                        gxij * s_val + gx[i] * s_x[j] + gx[j] * s_x[i] + gv * s_ij;
                }
            }
        }

        // F_disp = -(pi * N_A / v) * H(T, v, x)
        let factor = -std::f64::consts::PI * N_AVOGADRO * inv_v;
        let mut out = PcSaftHelmholtz::zeros(nc);
        out.f = factor * h.f;
        out.f_t = factor * h.f_t;
        out.f_v = factor * (h.f_v - h.f * inv_v);
        out.f_tt = factor * h.f_tt;
        out.f_tv = factor * (h.f_tv - h.f_t * inv_v);
        out.f_vv = factor * (h.f_vv - 2.0 * h.f_v * inv_v + 2.0 * h.f * inv_v * inv_v);
        for i in 0..nc {
            out.f_n[i] = factor * h.f_n[i];
            out.f_tn[i] = factor * h.f_tn[i];
            out.f_vn[i] = factor * (h.f_vn[i] - h.f_n[i] * inv_v);
            for j in 0..nc {
                out.f_nn[i * nc + j] = factor * h.f_nn[i * nc + j];
            }
        }
        out
    }

    // =========================================================================
    // INTERNAL HARD-CHAIN AND DISPERSION KERNELS (FAST DENSITY PATH)
    // =========================================================================

    fn alpha_hard_chain(&self, rho: f64, eta: f64) -> f64 {
        if self.model.kind.is_simplified() {
            let om = 1.0 - eta;
            let om2 = om * om;
            let om3 = om2 * om;
            let alpha_hs = (4.0 * eta - 3.0 * eta * eta) / om2;
            let g = (1.0 - 0.5 * eta) / om3;
            let mut sum_chain = 0.0;
            for i in 0..self.model.components.len() {
                sum_chain += self.mole_fractions[i] * (self.model.components[i].m - 1.0);
            }
            self.m_bar * alpha_hs - sum_chain * g.ln()
        } else {
            let nc = self.model.components.len();
            if nc == 1 {
                let om = 1.0 - eta;
                let om2 = om * om;
                let om3 = om2 * om;
                let alpha_hs = (4.0 * eta - 3.0 * eta * eta) / om2;
                let g = (1.0 - 0.5 * eta) / om3;
                let m0 = self.model.components[0].m;
                let chain = if m0 > 1.0 { (m0 - 1.0) * g.ln() } else { 0.0 };
                return self.m_bar * alpha_hs - chain;
            }

            let (zeta0, zeta1, zeta2, zeta3) = self.bmcsl_zetas(rho);
            let en_m_z3 = 1.0 - zeta3;
            let en_m_z3_2 = en_m_z3 * en_m_z3;
            let en_m_z3_3 = en_m_z3_2 * en_m_z3;
            let z2_2 = zeta2 * zeta2;
            let z2_3 = z2_2 * zeta2;
            let z3_2 = zeta3 * zeta3;
            let ln_term = en_m_z3.abs().ln();

            let alpha_hs = (1.0 / zeta0)
                * (3.0 * zeta1 * zeta2 / en_m_z3
                    + z2_3 / (zeta3 * en_m_z3_2)
                    + (z2_3 / z3_2 - zeta0) * ln_term);

            let mut sum_chain = 0.0;
            for i in 0..nc {
                let mi = self.model.components[i].m;
                if mi > 1.0 {
                    let mu_ii = 0.5 * self.d[i];
                    let g_ii = 1.0 / en_m_z3
                        + 3.0 * zeta2 * mu_ii / en_m_z3_2
                        + 2.0 * z2_2 * mu_ii * mu_ii / en_m_z3_3;
                    sum_chain += self.mole_fractions[i] * (mi - 1.0) * g_ii.ln();
                }
            }
            self.m_bar * alpha_hs - sum_chain
        }
    }

    fn bmcsl_zetas(&self, rho: f64) -> (f64, f64, f64, f64) {
        let pref = (std::f64::consts::PI * N_AVOGADRO / 6.0) * rho;
        let mut z0 = 0.0;
        let mut z1 = 0.0;
        let mut z2 = 0.0;
        let mut z3 = 0.0;
        for i in 0..self.model.components.len() {
            let term = self.mole_fractions[i] * self.model.components[i].m;
            let di = self.d[i];
            let di2 = di * di;
            z0 += term;
            z1 += term * di;
            z2 += term * di2;
            z3 += term * di2 * di;
        }
        (pref * z0, pref * z1, pref * z2, pref * z3)
    }

    fn d_alpha_hard_chain_drho(&self, rho: f64, eta: f64) -> f64 {
        let eta_rho = self.eta_factor;
        if self.model.kind.is_simplified() {
            let om = 1.0 - eta;
            let om2 = om * om;
            let om3 = om2 * om;
            let om4 = om2 * om2;
            let alpha_hs_eta = 2.0 * (2.0 - eta) / om3;
            let g = (1.0 - 0.5 * eta) / om3;
            let g_eta = (5.0 - 2.0 * eta) / (2.0 * om4);
            let mut sum_chain = 0.0;
            for i in 0..self.model.components.len() {
                sum_chain += self.mole_fractions[i] * (self.model.components[i].m - 1.0);
            }
            (self.m_bar * alpha_hs_eta - sum_chain * g_eta / g) * eta_rho
        } else {
            let nc = self.model.components.len();
            if nc == 1 {
                let om = 1.0 - eta;
                let om2 = om * om;
                let om3 = om2 * om;
                let om4 = om2 * om2;
                let alpha_hs_eta = 2.0 * (2.0 - eta) / om3;
                let g = (1.0 - 0.5 * eta) / om3;
                let g_eta = (5.0 - 2.0 * eta) / (2.0 * om4);
                let m0 = self.model.components[0].m;
                let sum_chain = if m0 > 1.0 { m0 - 1.0 } else { 0.0 };
                return (self.m_bar * alpha_hs_eta - sum_chain * g_eta / g) * eta_rho;
            }

            let (z0, z1, z2, z3) = self.bmcsl_zetas(rho);
            let en_m_z3 = 1.0 - z3;
            let en_m_z3_2 = en_m_z3 * en_m_z3;
            let en_m_z3_3 = en_m_z3_2 * en_m_z3;
            let z2_2 = z2 * z2;
            let z2_3 = z2_2 * z2;
            let z3_2 = z3 * z3;
            let z3_3 = z3_2 * z3;
            let pref = 1.0 / z0;
            let ln_term = en_m_z3.abs().ln();

            let alpha_hs = pref
                * (3.0 * z1 * z2 / en_m_z3
                    + z2_3 / (z3 * en_m_z3_2)
                    + (z2_3 / z3_2 - z0) * ln_term);

            let da_dz0 = -pref * (alpha_hs + ln_term);
            let da_dz1 = pref * 3.0 * z2 / en_m_z3;
            let da_dz2 = pref
                * (3.0 * z1 / en_m_z3
                    + 3.0 * z2_2 / (z3 * en_m_z3_2)
                    + (3.0 * z2_2 / z3_2) * ln_term);
            let da_dz3 = pref
                * (3.0 * z1 * z2 / en_m_z3_2 + z2_3 * (3.0 * z3 - 1.0) / (z3_2 * en_m_z3_3)
                    - (2.0 * z2_3 / z3_3) * ln_term
                    - (z2_3 - z0 * z3_2) / (z3_2 * en_m_z3));

            let da_drho = (da_dz0 * z0 + da_dz1 * z1 + da_dz2 * z2 + da_dz3 * z3) / rho;

            let mut dchain_drho = 0.0;
            let d1 = 1.0 / en_m_z3;
            let d2 = d1 * d1;
            let d3 = d1 * d2;
            let d4 = d2 * d2;

            for i in 0..nc {
                let mi = self.model.components[i].m;
                if mi > 1.0 {
                    let mu = 0.5 * self.d[i];
                    let mu2 = mu * mu;
                    let g = d1 + 3.0 * z2 * mu * d2 + 2.0 * z2_2 * mu2 * d3;
                    let dg_dz2 = 3.0 * mu * d2 + 4.0 * d3 * mu2 * z2;
                    let dg_dz3 = d2 + 6.0 * z2 * d3 * mu + 6.0 * z2_2 * d4 * mu2;
                    let dg_drho = (dg_dz2 * z2 + dg_dz3 * z3) / rho;
                    dchain_drho += self.mole_fractions[i] * (mi - 1.0) * (dg_drho / g);
                }
            }

            self.m_bar * da_drho - dchain_drho
        }
    }

    fn alpha_dispersion(&self, rho: f64, eta: f64) -> f64 {
        let (i1, i2) = self.dispersion_i1_i2(eta);
        let c1 = self.dispersion_c1(eta);
        -(std::f64::consts::PI * N_AVOGADRO * rho)
            * (2.0 * i1 * self.m2e1s3_val + self.m_bar * c1 * i2 * self.m2e2s3_val)
    }

    fn dispersion_i1_i2(&self, eta: f64) -> (f64, f64) {
        let mut i1 = 0.0;
        let mut i2 = 0.0;
        let mut eta_k = 1.0;
        for k in 0..7 {
            i1 += self.a_disp[k] * eta_k;
            i2 += self.b_disp[k] * eta_k;
            eta_k *= eta;
        }
        (i1, i2)
    }

    fn dispersion_c1(&self, eta: f64) -> f64 {
        let e2 = eta * eta;
        let e3 = e2 * eta;
        let e4 = e2 * e2;
        let om = 1.0 - eta;
        let om2 = om * om;
        let om4 = om2 * om2;
        let q = om * (2.0 - eta);
        let o1 = (8.0 * eta - 2.0 * e2) / om4;
        let o2 = (20.0 * eta - 27.0 * e2 + 12.0 * e3 - 2.0 * e4) / (q * q);
        1.0 / (1.0 + self.m_bar * o1 + (1.0 - self.m_bar) * o2)
    }

    fn dispersion_c1_eta(&self, eta: f64) -> f64 {
        let e2 = eta * eta;
        let e3 = e2 * eta;
        let om = 1.0 - eta;
        let om2 = om * om;
        let om5 = om2 * om2 * om;
        let q = om * (2.0 - eta);
        let q3 = q * q * q;
        let c1 = self.dispersion_c1(eta);
        let o1_eta = (-4.0 * e2 + 20.0 * eta + 8.0) / om5;
        let o2_eta = (2.0 * e3 + 12.0 * e2 - 48.0 * eta + 40.0) / q3;
        -(self.m_bar * o1_eta + (1.0 - self.m_bar) * o2_eta) * c1 * c1
    }

    fn dispersion_c1_etaeta(&self, eta: f64) -> f64 {
        let e2 = eta * eta;
        let e3 = e2 * eta;
        let om = 1.0 - eta;
        let om2 = om * om;
        let om6 = om2 * om2 * om2;
        let q = om * (2.0 - eta);
        let q2 = q * q;
        let q4 = q2 * q2;
        let c1 = self.dispersion_c1(eta);
        let c1_eta = self.dispersion_c1_eta(eta);
        let o1_etaeta = ((-8.0 * eta + 20.0) * om + (-4.0 * e2 + 20.0 * eta + 8.0) * 5.0) / om6;
        let o2_etaeta = ((6.0 * e2 + 24.0 * eta - 48.0) * q
            + (2.0 * e3 + 12.0 * e2 - 48.0 * eta + 40.0) * (3.0 - 2.0 * eta) * 3.0)
            / q4;
        2.0 * c1_eta * c1_eta / c1
            - (self.m_bar * o1_etaeta + (1.0 - self.m_bar) * o2_etaeta) * c1 * c1
    }

    fn d_alpha_dispersion_drho(&self, _rho: f64, eta: f64) -> f64 {
        let (i1, i2) = self.dispersion_i1_i2(eta);
        let mut i1_eta = 0.0;
        let mut i2_eta = 0.0;
        let mut eta_k = 1.0;
        for k in 1..7 {
            i1_eta += (k as f64) * self.a_disp[k] * eta_k;
            i2_eta += (k as f64) * self.b_disp[k] * eta_k;
            eta_k *= eta;
        }

        let e2 = eta * eta;
        let e3 = e2 * eta;
        let e4 = e2 * e2;
        let om = 1.0 - eta;
        let om2 = om * om;
        let om4 = om2 * om2;
        let om5 = om4 * om;
        let q = om * (2.0 - eta);
        let q2 = q * q;
        let q3 = q2 * q;

        let o1 = (8.0 * eta - 2.0 * e2) / om4;
        let o2 = (20.0 * eta - 27.0 * e2 + 12.0 * e3 - 2.0 * e4) / q2;
        let c1 = 1.0 / (1.0 + self.m_bar * o1 + (1.0 - self.m_bar) * o2);

        let o1_eta = (-4.0 * e2 + 20.0 * eta + 8.0) / om5;
        let o2_eta = (2.0 * e3 + 12.0 * e2 - 48.0 * eta + 40.0) / q3;
        let c1_eta = -(self.m_bar * o1_eta + (1.0 - self.m_bar) * o2_eta) * c1 * c1;

        let term1 = 2.0 * (i1 + eta * i1_eta) * self.m2e1s3_val;
        let term2 = self.m_bar * (c1 * i2 + eta * (c1_eta * i2 + c1 * i2_eta)) * self.m2e2s3_val;

        -(std::f64::consts::PI * N_AVOGADRO) * (term1 + term2)
    }

    fn hard_chain_energy(&self, rho: f64, eta: f64) -> (f64, f64) {
        (
            self.alpha_hard_chain(rho, eta),
            self.d_alpha_hard_chain_drho(rho, eta),
        )
    }

    fn hard_chain_t_derivatives(&self, rho: f64, eta: f64) -> (f64, f64) {
        if self.model.kind.is_simplified() || self.model.components.len() == 1 {
            let (ahc_eta, ahc_etaeta) = self.pure_or_simplified_hard_chain_eta_derivatives(eta);
            let eta_t = self.packing_fraction_t_derivative(rho);
            let eta_tt = self.packing_fraction_tt_derivative(rho);
            (
                ahc_eta * eta_t,
                ahc_eta * eta_tt + ahc_etaeta * eta_t * eta_t,
            )
        } else {
            let derivatives = self.hard_chain_bmcsl_derivatives(rho);
            (derivatives.d_t, derivatives.d_tt)
        }
    }

    fn hard_chain_second_derivatives(&self, rho: f64, eta: f64) -> (f64, f64) {
        if self.model.kind.is_simplified() || self.model.components.len() == 1 {
            let (ahc_eta, ahc_etaeta) = self.pure_or_simplified_hard_chain_eta_derivatives(eta);
            let eta_t = self.packing_fraction_t_derivative(rho);
            (
                ahc_etaeta * self.eta_factor * self.eta_factor,
                ahc_eta * (eta_t / rho) + ahc_etaeta * self.eta_factor * eta_t,
            )
        } else {
            let derivatives = self.hard_chain_bmcsl_derivatives(rho);
            let rho2 = rho * rho;
            (
                2.0 * derivatives.d_v / (rho2 * rho) + derivatives.d_vv / (rho2 * rho2),
                -derivatives.d_tv / rho2,
            )
        }
    }

    fn packing_fraction_tt_derivative(&self, rho: f64) -> f64 {
        let mut sum_eta_tt = 0.0;
        for i in 0..self.model.components.len() {
            let segments = self.model.components[i].m;
            let di = self.d[i];
            let dti = self.d_t[i];
            sum_eta_tt += self.mole_fractions[i]
                * segments
                * (6.0 * di * dti * dti + 3.0 * di * di * self.d_tt[i]);
        }
        (std::f64::consts::PI * N_AVOGADRO / 6.0) * rho * sum_eta_tt
    }

    fn packing_fraction_t_derivative(&self, rho: f64) -> f64 {
        let mut sum_eta_t = 0.0;
        for i in 0..self.model.components.len() {
            let di = self.d[i];
            sum_eta_t +=
                self.mole_fractions[i] * self.model.components[i].m * 3.0 * di * di * self.d_t[i];
        }
        (std::f64::consts::PI * N_AVOGADRO / 6.0) * rho * sum_eta_t
    }

    fn pure_or_simplified_hard_chain_eta_derivatives(&self, eta: f64) -> (f64, f64) {
        let om = 1.0 - eta;
        let om2 = om * om;
        let om3 = om2 * om;
        let om4 = om2 * om2;
        let om5 = om4 * om;

        let alpha_hs_eta = 2.0 * (2.0 - eta) / om3;
        let alpha_hs_etaeta = 2.0 * (5.0 - 2.0 * eta) / om4;
        let g = (1.0 - 0.5 * eta) / om3;
        let g_eta = (5.0 - 2.0 * eta) / (2.0 * om4);
        let g_etaeta = 3.0 * (3.0 - eta) / om5;

        let mut sum_m_m1 = 0.0;
        for i in 0..self.model.components.len() {
            let mi = self.model.components[i].m;
            if self.model.kind.is_simplified() || mi > 1.0 {
                sum_m_m1 += self.mole_fractions[i] * (mi - 1.0);
            }
        }

        let r = g_eta / g;
        let ahc_eta = self.m_bar * alpha_hs_eta - sum_m_m1 * r;
        let ahc_etaeta = self.m_bar * alpha_hs_etaeta - sum_m_m1 * (g_etaeta / g - r * r);
        (ahc_eta, ahc_etaeta)
    }

    fn hard_chain_bmcsl_derivatives(&self, rho: f64) -> HardChainDerivatives {
        let volume = 1.0 / rho;
        let moles = &self.mole_fractions;
        let diameters = HsDiameters {
            d: self.d.clone(),
            d_t: self.d_t.clone(),
            d_tt: self.d_tt.clone(),
        };
        let segments: Vec<f64> = self
            .model
            .components
            .iter()
            .map(|component| component.m)
            .collect();
        let moments = bmcsl::packing_moments(volume, moles, &diameters, &segments);
        let (hard_sphere_t, hard_sphere_tt) = bmcsl::hard_sphere_energy_t_derivatives(&moments);
        let (hard_sphere_v, hard_sphere_vv, hard_sphere_tv) =
            bmcsl::hard_sphere_energy_v_derivatives(&moments);

        let mut derivatives = HardChainDerivatives {
            d_t: self.m_bar * hard_sphere_t,
            d_tt: self.m_bar * hard_sphere_tt,
            d_v: self.m_bar * hard_sphere_v,
            d_vv: self.m_bar * hard_sphere_vv,
            d_tv: self.m_bar * hard_sphere_tv,
        };
        for component in 0..self.model.components.len() {
            let seg = self.model.components[component].m;
            if seg <= 1.0 {
                continue;
            }
            let temperature =
                bmcsl::pair_correlation_t_derivatives(&moments, &diameters, component, component);
            let volume =
                bmcsl::pair_correlation_v_derivatives(&moments, &diameters, component, component);
            let weight = self.mole_fractions[component] * (seg - 1.0);
            derivatives.d_t -= weight * temperature.log_d_t;
            derivatives.d_tt -= weight * temperature.log_d_tt;
            derivatives.d_v -= weight * volume.log_d_v;
            derivatives.d_vv -= weight * volume.log_d_vv;
            derivatives.d_tv -= weight * volume.log_d_tv;
        }
        derivatives
    }

    fn dispersion_rho_t_derivatives(&self, rho: f64, eta: f64) -> (f64, f64, f64, f64, f64) {
        let (i1, i2) = self.dispersion_i1_i2(eta);
        let mut i1_eta = 0.0;
        let mut i2_eta = 0.0;
        let mut i1_etaeta = 0.0;
        let mut i2_etaeta = 0.0;
        let mut eta_k1 = 1.0;
        let mut eta_k2 = 1.0;

        for k in 1..7 {
            i1_eta += (k as f64) * self.a_disp[k] * eta_k1;
            i2_eta += (k as f64) * self.b_disp[k] * eta_k1;
            eta_k1 *= eta;
        }
        for k in 2..7 {
            i1_etaeta += (k * (k - 1)) as f64 * self.a_disp[k] * eta_k2;
            i2_etaeta += (k * (k - 1)) as f64 * self.b_disp[k] * eta_k2;
            eta_k2 *= eta;
        }

        let c1 = self.dispersion_c1(eta);
        let c1_eta = self.dispersion_c1_eta(eta);
        let c1_etaeta = self.dispersion_c1_etaeta(eta);

        let eta_t = self.packing_fraction_t_derivative(rho);

        let i1_t = i1_eta * eta_t;
        let i2_t = i2_eta * eta_t;
        let c1_t = c1_eta * eta_t;

        let adisp = -(std::f64::consts::PI * N_AVOGADRO * rho)
            * (2.0 * i1 * self.m2e1s3_val + self.m_bar * c1 * i2 * self.m2e2s3_val);

        let adisp_rho = -(std::f64::consts::PI * N_AVOGADRO)
            * (2.0 * (i1 + eta * i1_eta) * self.m2e1s3_val
                + self.m_bar * (c1 * i2 + eta * (c1_eta * i2 + c1 * i2_eta)) * self.m2e2s3_val);

        let adisp_t = -(std::f64::consts::PI * N_AVOGADRO * rho)
            * (2.0 * (i1_t * self.m2e1s3_val + i1 * self.m2e1s3_t)
                + self.m_bar
                    * (c1_t * i2 * self.m2e2s3_val
                        + c1 * i2_t * self.m2e2s3_val
                        + c1 * i2 * self.m2e2s3_t));

        let adisp_rhorho = -(std::f64::consts::PI * N_AVOGADRO * self.eta_factor)
            * (2.0 * (2.0 * i1_eta + eta * i1_etaeta) * self.m2e1s3_val
                + self.m_bar
                    * (2.0 * (c1_eta * i2 + c1 * i2_eta)
                        + eta * (c1_etaeta * i2 + 2.0 * c1_eta * i2_eta + c1 * i2_etaeta))
                    * self.m2e2s3_val);

        let f_eta = 2.0 * i1_eta * self.m2e1s3_val
            + self.m_bar * (c1_eta * i2 + c1 * i2_eta) * self.m2e2s3_val;
        let f_t = 2.0 * i1 * self.m2e1s3_t + self.m_bar * c1 * i2 * self.m2e2s3_t;
        let f_etaeta = 2.0 * i1_etaeta * self.m2e1s3_val
            + self.m_bar
                * (c1_etaeta * i2 + 2.0 * c1_eta * i2_eta + c1 * i2_etaeta)
                * self.m2e2s3_val;
        let f_eta_t =
            2.0 * i1_eta * self.m2e1s3_t + self.m_bar * (c1_eta * i2 + c1 * i2_eta) * self.m2e2s3_t;

        let adisp_rhot = -(std::f64::consts::PI * N_AVOGADRO)
            * (f_eta * eta_t
                + f_t
                + rho * self.eta_factor * (f_etaeta * eta_t + f_eta_t)
                + f_eta * eta_t);

        (adisp, adisp_rho, adisp_rhorho, adisp_t, adisp_rhot)
    }

    fn dispersion_tt_derivative(&self, rho: f64, eta: f64) -> f64 {
        let (i1, i2) = self.dispersion_i1_i2(eta);
        let mut i1_eta = 0.0;
        let mut i2_eta = 0.0;
        let mut eta_k1 = 1.0;
        for k in 1..7 {
            i1_eta += (k as f64) * self.a_disp[k] * eta_k1;
            i2_eta += (k as f64) * self.b_disp[k] * eta_k1;
            eta_k1 *= eta;
        }
        let mut i1_etaeta = 0.0;
        let mut i2_etaeta = 0.0;
        let mut eta_k2 = 1.0;
        for k in 2..7 {
            i1_etaeta += (k * (k - 1)) as f64 * self.a_disp[k] * eta_k2;
            i2_etaeta += (k * (k - 1)) as f64 * self.b_disp[k] * eta_k2;
            eta_k2 *= eta;
        }

        let c1 = self.dispersion_c1(eta);
        let c1_eta = self.dispersion_c1_eta(eta);
        let c1_etaeta = self.dispersion_c1_etaeta(eta);

        let eta_t = self.packing_fraction_t_derivative(rho);
        let eta_tt = self.packing_fraction_tt_derivative(rho);

        let f_eta = 2.0 * i1_eta * self.m2e1s3_val
            + self.m_bar * (c1_eta * i2 + c1 * i2_eta) * self.m2e2s3_val;
        let f_etaeta = 2.0 * i1_etaeta * self.m2e1s3_val
            + self.m_bar
                * (c1_etaeta * i2 + 2.0 * c1_eta * i2_eta + c1 * i2_etaeta)
                * self.m2e2s3_val;
        let f_eta_t =
            2.0 * i1_eta * self.m2e1s3_t + self.m_bar * (c1_eta * i2 + c1 * i2_eta) * self.m2e2s3_t;
        let f_tt = 2.0 * i1 * self.m2e1s3_tt + self.m_bar * c1 * i2 * self.m2e2s3_tt;

        let second = f_etaeta * eta_t * eta_t + f_eta * eta_tt + 2.0 * f_eta_t * eta_t + f_tt;
        -(std::f64::consts::PI * N_AVOGADRO * rho) * second
    }

    /// First mole-number derivatives of the extensive residual Helmholtz energy
    /// (hard-chain + dispersion + association + multipolar).
    fn helmholtz_fn_derivatives(&self, rho: f64) -> Vec<f64> {
        let mut out = self.nonassoc_helmholtz_fn_derivatives(rho);
        let v_unit = 1.0 / rho;
        if self.model.assoc.has_sites() {
            let rdf = PcSaftRdfProvider::new(self, v_unit);
            if let Ok(first) = association::assoc_helmholtz_first(
                &self.model.assoc,
                &rdf,
                self.temperature,
                v_unit,
                &self.mole_fractions,
            ) {
                for (a, b) in out.iter_mut().zip(first.f_n) {
                    *a += b;
                }
            } else {
                out.fill(f64::NAN);
                return out;
            }
        }
        if let Some(mp) = &self.model.multipol
            && mp.is_active()
        {
            if let Ok(mpol_fn) = mp.fres_fn(self.temperature, v_unit, &self.mole_fractions) {
                for (a, b) in out.iter_mut().zip(mpol_fn) {
                    *a += b;
                }
            } else {
                out.fill(f64::NAN);
                return out;
            }
        }
        out
    }

    fn nonassoc_helmholtz_fn_derivatives(&self, rho: f64) -> Vec<f64> {
        let v = 1.0 / rho;
        let n = &self.mole_fractions;
        let eta = self.packing_fraction(rho);
        let nc = n.len();
        let (i1, i2) = self.dispersion_i1_i2(eta);
        let c1 = self.dispersion_c1(eta);
        let e2 = eta * eta;
        let om = 1.0 - eta;
        let om2 = om * om;
        let om3 = om2 * om;
        let om4 = om2 * om2;
        let q_disp = om * (2.0 - eta);
        let o1 = (8.0 * eta - 2.0 * e2) / om4;
        let o2 = (20.0 * eta - 27.0 * e2 + 12.0 * e2 * eta - 2.0 * e2 * e2) / (q_disp * q_disp);
        let mut powers = [1.0; 7];
        for k in 1..7 {
            powers[k] = powers[k - 1] * eta;
        }
        let mut i1_eta = 0.0;
        let mut i2_eta = 0.0;
        for k in 1..7 {
            i1_eta += k as f64 * self.a_disp[k] * powers[k - 1];
            i2_eta += k as f64 * self.b_disp[k] * powers[k - 1];
        }
        let c1_eta = self.dispersion_c1_eta(eta);
        let (z0, z1, z2, z3) = self.bmcsl_zetas(rho);
        let q = 1.0 - z3;
        let q2 = q * q;
        let q3 = q2 * q;
        let q4 = q2 * q2;
        let z2_2 = z2 * z2;
        let simplified = self.model.kind.is_simplified();
        let hs = if !simplified {
            Some(bmcsl::boublik_helmholtz_energy([z0, z1, z2, z3]))
        } else {
            None
        };
        let moment_prefactor = std::f64::consts::PI * N_AVOGADRO / (6.0 * v);
        let alpha_hs = (4.0 * eta - 3.0 * eta * eta) / om2;
        let g = (1.0 - 0.5 * eta) / om3;
        let ahc_eta = if simplified {
            self.pure_or_simplified_hard_chain_eta_derivatives(eta).0
        } else {
            0.0
        };
        let mut chain_log = vec![0.0; nc];
        let mut chain_z2 = 0.0;
        let mut chain_z3 = 0.0;
        if !simplified {
            for j in 0..nc {
                let weight = self.model.components[j].m - 1.0;
                if weight <= 0.0 {
                    continue;
                }
                let mu = 0.5 * self.d[j];
                let mu2 = mu * mu;
                let gij = 1.0 / q + 3.0 * z2 * mu / q2 + 2.0 * z2_2 * mu2 / q3;
                let dg2 = 3.0 * mu / q2 + 4.0 * z2 * mu2 / q3;
                let dg3 = 1.0 / q2 + 6.0 * z2 * mu / q3 + 6.0 * z2_2 * mu2 / q4;
                chain_log[j] = weight * gij.ln();
                chain_z2 += n[j] * weight * dg2 / gij;
                chain_z3 += n[j] * weight * dg3 / gij;
            }
        }
        let disp_value = self.alpha_dispersion(rho, eta);
        let mut out = chain_log;
        for (i, value) in out.iter_mut().enumerate() {
            let log_chain = *value;
            let mi = self.model.components[i].m;
            let di = self.d[i];
            let di2 = di * di;
            let eta_n = moment_prefactor * mi * di2 * di;
            let hard = if simplified {
                mi * alpha_hs - (mi - 1.0) * g.ln() + ahc_eta * eta_n
            } else {
                let hs = hs.as_ref().expect("BMCSL derivatives for PC-SAFT");
                let dz = [
                    moment_prefactor * mi,
                    moment_prefactor * mi * di,
                    moment_prefactor * mi * di2,
                    eta_n,
                ];
                mi * hs.value
                    + self.m_bar * hs.dzeta.iter().zip(dz).map(|(a, b)| a * b).sum::<f64>()
                    - log_chain
                    - chain_z2 * dz[2]
                    - chain_z3 * dz[3]
            };
            let dm = self.m_bar_n[i] * self.composition_sum;
            let mut di1 = i1_eta * eta_n;
            let mut di2 = i2_eta * eta_n;
            for (k, &power) in powers.iter().enumerate() {
                di1 += self.a_disp_n[i][k] * self.composition_sum * power;
                di2 += self.b_disp_n[i][k] * self.composition_sum * power;
            }
            let dc1 = c1_eta * eta_n - c1 * c1 * (o1 - o2) * dm;
            let de1 = self.m2e1s3_n[i] * self.composition_sum;
            let de2 = self.m2e2s3_n[i] * self.composition_sum;
            let disp = 2.0 * disp_value
                - std::f64::consts::PI
                    * N_AVOGADRO
                    * rho
                    * (2.0 * (di1 * self.m2e1s3_val + i1 * de1)
                        + (dm * c1 * i2 + self.m_bar * dc1 * i2 + self.m_bar * c1 * di2)
                            * self.m2e2s3_val
                        + self.m_bar * c1 * i2 * de2);
            *value = hard + disp;
        }
        out
    }
}

/// Pressure evaluation output package with thermodynamic differentials.
#[derive(Debug, Clone, Copy)]
pub struct PcSaftPressureJet {
    pub pressure: f64,
    pub compressibility: f64,
    pub residual_helmholtz: f64,
    pub dp_drho: f64,
    pub dp_dv: f64,
    pub dp_dt: f64,
}

/// Thermodynamic state property bundle.
#[derive(Debug, Clone)]
pub struct PcSaftState {
    pub pressure: f64,
    pub compressibility: f64,
    pub residual_helmholtz: f64,
    pub u_residual: f64,
    pub h_residual: f64,
    pub s_residual: f64,
    pub a_residual: f64,
    pub g_residual: f64,
    pub h_total: f64,
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

/// Pure component vapor-liquid coexistence saturation properties.
#[derive(Debug, Clone, Copy)]
pub struct PcSaftSaturation {
    pub pressure: f64,
    pub temperature: f64,
    pub liquid_volume: f64,
    pub vapor_volume: f64,
    pub enthalpy_vaporization: f64,
    pub entropy_vaporization: f64,
}

#[cfg(test)]
mod mole_derivative_tests {
    use super::*;

    #[test]
    fn chain_segment_threshold_matches_energy_derivative() {
        for kind in [PcSaftKind::PcSaft, PcSaftKind::SPcSaft] {
            let mut model = PcSaft::new_nonassoc(kind, &["CO2", "C1", "C3"], "Default").unwrap();
            model.components[0].m = 0.75;
            model.components[1].m = 1.0;
            let n = [0.2, 0.3, 0.5];
            let rho = 3000.0;
            let volume = 1.0 / rho;
            let prepared = model.prepare(320.0, &n).unwrap();
            let analytic = prepared.helmholtz_fn_derivatives(rho);
            let energy = |amounts: [f64; 3]| {
                let total: f64 = amounts.iter().sum();
                total
                    * model
                        .prepare(320.0, &amounts)
                        .unwrap()
                        .residual_helmholtz_energy(total / volume)
            };
            for i in 0..3 {
                let h = 1e-5;
                let mut plus = n;
                plus[i] += h;
                let mut minus = n;
                minus[i] -= h;
                let numeric = (energy(plus) - energy(minus)) / (2.0 * h);
                assert!(
                    (analytic[i] - numeric).abs() < 2e-7 * (1.0 + numeric.abs()),
                    "{kind:?} component {i}: {} {numeric}",
                    analytic[i]
                );
            }
        }
    }
}
