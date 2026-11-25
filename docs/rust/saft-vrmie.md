# SAFT-VR Mie and SAFT-VRQ Mie

This module implements native Rust versions of the SAFT-VR Mie equation of state
(Lafitte et al., *J. Chem. Phys.* 139, 154504, 2013), Wertheim TPT1 association
with the Dufal et al. (2015) Mie radial distribution integral, and the quantum
extension SAFT-VRQ Mie (Aasen et al., 2019, 2020) for quantum fluids
($\text{H}_2, \text{He}, \text{Ne}$).

## Theory and Formulation

The extensive reduced residual Helmholtz energy is evaluated as:

$$F^{\text{res}}(T, V, \mathbf{n}) = \frac{A^{\text{res}}}{R T} = n_{\text{tot}} \left( \bar{m} a^{\text{mono}} + a^{\text{chain}} + a^{\text{assoc}} \right)$$

where $\bar{m} = \sum_i x_i m_i$ is the average segment number per mole,
$a^{\text{mono}} = a^{\text{hs}} + a^{\text{hs,extra}} + \beta a_1 + \beta^2 a_2 + \beta^3 a_3$
is the monomer contribution per segment (with optional non-additive extra
hard-sphere term), $a^{\text{chain}}$ is the chain formation contribution per
mole, and $a^{\text{assoc}}$ is the Wertheim TPT1 association contribution per
mole.

### Core Submodels and Controls

- **Hard-sphere references (`HardSphereReference`)**: `Lafitte` (default),
  `SingleComponent`, `Additive`, `NonAdditive` (preserves active `PackingRule`,
  matching `set_hs_reference(NON_ADD_HS_REF)`), and `AdditiveExact`, plus
  explicit `extra_hard_sphere` and `kolafa_pure_reference` controls.
- **Term toggles (`MieTerms`)**: independent switches for `hard_sphere`,
  `first` ($a_1$), `second` ($a_2$), `third` ($a_3$), and `chain`.
- **Submodel options (`MieOptions`)**: `PackingRule::{Lafitte, Linear}`,
  `bmcsl_compressibility` ($K_{\text{hs}}$ Carnahan-Starling vs BMCSL),
  `pair_prefactor_a3`, `geometric_epsilon`, `exact_cross_potential`,
  `QuantumVariant::{FeynmanHibbs, JaenKahn}`, `quantum_correct_a2`, and
  `extra_hard_sphere`.
- **Potential truncation and shift**: pure-fluid `set_cutoff(Some((r_cut, shift)))`
  for truncated and truncated-shifted Mie/quantum-Mie dispersion and pair
  potentials (shifted `alpha_vdw` rejects cutoffs at, below, or within
  $10^{-12} r_{\min}$ of $r_{\min}$; near $r_{\min}$ ($r_{\text{cut}} \le 2 r_{\min}$),
  it solves for the inner shifted zero in offset coordinates
  $d_s = r_{\min} - \sigma_{\text{eff}}$ scaled to
  $\varepsilon_{\text{eff}} = u^*(r_{\text{cut}}) + 1$ and uses
  cancellation-free 21-point Gauss-Kronrod quadrature on each side of $r_{\min}$
  for shallow wells $\varepsilon_{\text{eff}} < 0.25$; for deep wells and large
  finite cutoffs ($r_{\text{cut}} > 2 r_{\min}$), it evaluates
  $u_{\text{cut}} = -u^*(r_{\text{cut}})$ directly, solves for
  $\delta_s = \sigma_{\text{eff}} - 1 \ge 0$, and writes the tail via
  $r_{\text{cut}}^{3 - \lambda}$ to avoid $1 - \varepsilon_{\text{eff}}$
  cancellation and $r_{\text{cut}}^3$ overflow).
- **Association (`rust/src/association.rs`, `rust/src/saft_vrmie/association.rs`)**:
  reuses the shared Wertheim TPT1 site solver (`AssocScheme`, `Standard` and
  `Elliot` $\Delta$ combining rules, exact infinite-dilution zero-mole site
  fractions, and fallback-only machine-precision Newton polish so routine
  SAFT-VR Mie, PC-SAFT, and CPA calls incur zero extra residual/Jacobian work)
  paired with the Dufal et al. (2015) polynomial Mie association integral
  (`association_rdf_tvn` and `rdf::saftvrmie_rdf`). `PreparedSaftVrMie`
  precomputes the 11 temperature-dependent Dufal outer coefficients per
  associating pair in Veltkamp-Dekker double-double precision, then evaluates
  the 11-step density polynomial via a fast standard `f64` Horner path with a
  running cancellation monitor (`abs_sum > 1024 * |val|`) that falls back to
  11-step Dekker double-double arithmetic only on cancellation-sensitive dense
  states. Controls include `has_association`, `set_assoc_params`,
  `set_delta_comb_rule`, and `solve_assoc_x_k`.
- **Quantum extensions**: `new_quantum` and `new_quantum_with_hs_ref`
  (FH0/FH1/FH2, `FeynmanHibbs` and `JaenKahn` variants, `H_PLANCK = 6.626069311e-34`),
  plus `set_mass`, `de_boer_parameter`, `set_de_boer_parameter`, and
  `de_broglie_wavelength`.
- **Model-side FMT/DFT hooks**: `fmt_compatibility`, `dispersion_tvn`,
  `hard_sphere_tvn`, `hard_sphere_diameters`, `hard_sphere_diameter_ij`,
  `effective_sigma_ij`, `effective_epsilon_ij`, `alpha_vdw`, `pair_potential`
  (matching `saft_interface::potential` with default `500.0 * eps_ij/k_B` core
  cap and uncapped analytical force), and `pair_potential_with_cap` (custom cap
  or `None` for raw uncapped potential).

## Database and Parameters

Generated deterministically from `src/saftvrmie_datadb.f90` and
`src/saftvrmie_association.f90` by `tools/rust_port/generate_saftvrmie_data.py`
into `rust/src/saft_vrmie_data.rs`:

- **77 pure fluid records** ($m, \sigma, \epsilon/k_B, \lambda_r, \lambda_a$,
  Feynman-Hibbs order and molecular mass, plus association scheme, $\epsilon_{\text{assoc}}$,
  and $\beta_{\text{assoc}}$ for all 8 associating fluids: 4 $\times$ `4C`,
  3 $\times$ `3B`, 1 $\times$ `4B`).
- **18 binary $k_{ij}$ records** and **12 binary $l_{ij}$ records** (with
  `set_kij`/`get_kij`, `set_lij`/`get_lij`, `set_gammaij`/`get_gammaij`).
- **66 $\times$ 7 Dufal `ASSOCIATION_B`** polynomial coefficient table.

## Extensive TVN Derivatives and Normal TP/TV Paths

- `residual_helmholtz_tvn(T, V, n)` computes $F^{\text{res}}$ and all ten first
  and second derivative blocks ($F, F_T, F_V, F_n, F_{TT}, F_{TV}, F_{Tn},
  F_{VV}, F_{Vn}, F_{nn}$) plus solved site fractions $X_k$ in unit-mole
  coordinates with exact Euler homogeneity scaling (verified invariant down to
  $10^{-140}$ mol, returning explicit errors on underflow/overflow).
- `dispersion_tvn(T, V, n)` and `hard_sphere_tvn(T, V, n)` return reduced molar
  Helmholtz energies and all ten TVN derivative blocks for DFT coupling.
- `volume_tp`, `state_tp`, and `ln_fugacity_derivatives_tp` consume all ten
  unit-mole TVN blocks to provide normal TP properties ($v, Z, H^{\text{res}},
  H, S^{\text{res}}, S, G^{\text{res}}, \ln\boldsymbol{\phi}, C_p, C_v$, speed
  of sound) together with exact analytical $(T, P, \mathbf{n})$ first
  derivatives. Following ThermoPack Fortran `eos` conventions, `dn` on
  `compressibility_factor` and `ln_fugacity` stores intensive raw-mole
  derivatives ($\partial q / \partial n_i$, degree $-1$ in $N$), while `dn` on
  `specific_volume`, `residual_enthalpy`, `enthalpy`, `residual_entropy`,
  `entropy`, and `residual_gibbs` stores partial molar quantities
  $\bar{q}_i = \partial(N q)/\partial n_i$ (degree $0$ in $N$). Call
  `SaftVrMieDerivatives::intensive_dn(total_moles)` to convert any partial
  molar vector into $\partial q / \partial n_i = (\bar{q}_i - q)/N$.
- `PreparedSaftVrMie` caches temperature-dependent Barker-Henderson diameters
  (21-point Gauss-Kronrod quadrature with self-contained classical/quantum
  diameter paths and shared `ln`/integer-exponent power evaluation) and
  effective quantum parameters for fast repeated TV evaluations
  (`helmholtz_terms` via zero-overhead `f64`, analytic `pressure_and_derivatives`
  and `state` via stack-buffered `Dual1`/`Dual2` without re-preparing
  temperature caches, `volume` with a strict $\le 10^{-6} P_{\text{target}}$
  relative pressure closure cap on every returned root, noise-guarded bracketing
  with spacing-bounded representable-`f64` ULP recovery (`<= 14` evaluations in
  deeply sub-ULP rejected dense regimes, `<= 50` otherwise, stopping as soon as
  the $10^{-6}$ cap is met), and common-pressure `MinimumGibbs` that checks for
  a clear vapor win against the bounded liquid candidate before paying for dense
  ULP recovery; and pure `saturation_pressure` using the stabilized
  common-pressure Gibbs difference $(F_l - F_v) + P(v_l - v_v)/(RT) - \ln(v_l/v_v)$
  with guarded mechanical-spinodal branch bracketing ($\eta \le 0.72$), verified
  liquid and vapor pressure closure ($\le 10^{-6}$ relative), ordinary Gibbs
  convergence $|\Delta G/(RT)| < 10^{-13}$ (or $< 10^{-10}$ on a tight bracket),
  a documented cold-liquid pressure-snap branch (when 1 ULP of liquid volume
  exceeds $10^{-6} P_{\text{sat}}$ yet $|P_l - p_{\text{eq}}| \le 5\times 10^{-6} p_{\text{eq}}$,
  reporting $P_l$ with a re-solved vapor root so both pressures close within
  $10^{-6}$ while $|\Delta G/(RT)| \le 5\times 10^{-6}$), and distinct errors
  for supercritical isotherms versus liquid packing-cap exhaustion). For fixed
  ultra-low liquid pressure targets where the physical 1-ULP pressure spacing
  $|\partial P/\partial v|\,\mathrm{ulp}(v_l)$ exceeds $2\times 10^{-6} P_{\text{target}}$
  (e.g., $\text{NH}_3$ at $120\text{ K}$ where $\Delta P_{\text{ulp}} \approx 1.18\times 10^{-6}\text{ Pa}$
  versus $\pm 3.05\times 10^{-7}\text{ Pa}$ at $0.305\text{ Pa}$ and
  $\pm 2.76\times 10^{-7}\text{ Pa}$ at $0.275800421086915\text{ Pa}$),
  `volume` returns `Ok` when a representable stable root lies inside the $10^{-6}$
  window and `Error::NotConverged` when the straddling `f64` volumes contain no
  closing root. `SaftVrMieState::s_total` combines ideal TV entropy with TV
  residual entropy so it agrees with `state_tp.entropy.value` and satisfies
  $\Delta H = T\Delta S$ at coexistence.

## Exact-Source Qualification

`tools/rust_port/validate_saftvrmie.py` compiles an isolated exact-source
Fortran probe (`tools/rust_port/saftvrmie_probe.f90`) linked against
`.reference` and qualifies **67 live cases** (**988 scalar comparisons** plus
**3 explicitly marked `UNSUPPORTED_SOURCE` slots** for classical
`de_Broglie_wavelength`), persisted in `tests/fixtures/saftvrmie_tvn.csv` and
verified in `tests/saft_vrmie.rs` with unit-aware per-slot tolerances,
non-finite rejection, tamper-detection tests, and an independently pinned
67-case roster:

- **TVN residual, dispersion, hard-sphere, and Dufal RDF** (43 cases covering
  classical, associating `NH3`/`H2O,MEOH`/`ETOH,C1`, Elliot rule, non-unit mole
  scaling, all 5 hard-sphere references, all 5 term masks, cutoff/shift,
  FH0/FH1/FH2, `additive-exact-quantum`, `linear-then-nonadd`, and all 9
  independent submodel options).
- **Normal TP states** (8 cases covering `CO2,C1` liquid/vapor, `NH3`,
  `H2O,MEOH` standard and Elliot, `ETOH,C1`, `H2` FH1, `H2,NE` FH2 across
  $Z, H^{\text{res}}, H, S^{\text{res}}, S, \ln\boldsymbol{\phi}$, all
  $(T, P, n)$ derivatives, and speed of sound).
- **Model-side FMT/DFT hooks** (5 cases covering classical, truncated-shifted
  at $r_{\text{cut}} = 3.5$ and $r_{\text{cut}} = 1000$,
  and quantum exact/arithmetic-cross diameters, $d_T$, $d_{TT}$, $d_{ij}$,
  $d_{ij,T}$, effective $\sigma_{ij}/\epsilon_{ij}$, $\alpha_{ij}$, and de
  Broglie wavelength with unit-resolved tolerances down to $10^{-29}$).
- **Pure-fluid saturation** (6 cases covering cold, boiling, and near-critical
  associating `H2O` at 300 K, 325 K, 373.15 K, 679.0 K, 679.03 K and `CO2` at
  260 K against Fortran `safe_bubP` and `specificVolume`, exact-pressure root
  checks to `1e-8`, plus $\Delta H = T\Delta S$; `tests/saft_vrmie.rs` also
  verifies `H2O` at 679.02 K and 679.05 K and supercritical rejection at
  680 K/700 K and `CO2` 320 K).
- **Pair potential and force** (5 cases covering `CO2` capped core
  $r = 0.5\sigma$ at $500\epsilon/k_B$ with uncapped force, attractive well,
  truncated-shifted inside/beyond cutoff, and FH2 `H2,NE` cross potential
  against Fortran `saft_interface::potential`).

Five verified upstream Fortran defects (extra hard-sphere $a_{VV}$,
`calc_saftvrmie_hard_sphere` $F_{TV}$ sign and unguarded optional `F`,
`test_fmt_compatibility` unassigned output, classical `de_Broglie_wavelength`
zero-mass check, and zero-mole $X_k = 1$) are documented in
[baseline-defects.md](baseline-defects.md) and independently verified by finite
differences and thermodynamic identities in `tests/saft_vrmie.rs`.

## Prepared-State Release Performance

Measured on Windows GNU Rust 1.100 (`--release`, thin LTO, `codegen-units=1`,
stack-buffered for $n_c \le 8$, symmetric pair loop, single association solve
per state). Matched 5-round interleaved `prepare` comparisons (`20,000`
iterations per round, excluding constructor time) comparing `a2ca7d9`
(out-of-line `PairPotential::evaluate`) against the self-contained
`PairPotential::diameter` + `pow_pos_pair` implementation:

| Workload | `a2ca7d9` median (range) | Current median (range) |
|---|---|---|
| `CO2/C1` `prepare` ($T=280\text{ K}, x=[0.4, 0.6]$) | 11.90 µs (11.30–12.42 µs) | 9.77 µs (9.56–9.94 µs) |
| `C1` `prepare` ($T=160\text{ K}, x=[1.0]$) | 4.79 µs (4.55–5.02 µs) | 4.08 µs (4.01–4.29 µs) |
| `H2` FH1 `prepare` ($T=40\text{ K}, x=[1.0]$) | 2.98 µs (2.70–3.29 µs) | 2.53 µs (2.37–2.69 µs) |

Prepared-state evaluation latencies (same release profile):

| Workload | Operation | Latency |
|---|---|---|
| `CO2/C1` ($T=280\text{ K}, x=[0.4, 0.6], \rho=5000$) | `helmholtz_terms` | 3.54 µs |
| `CO2/C1` | `pressure_and_derivatives` (analytic) | 12.2 µs |
| `CO2/C1` | `state` (analytic) | 26.8–30.4 µs |
| `CO2/C1` | `residual_helmholtz_tvn` (all 10 blocks) | 71.5 µs |
| `CO2/C1` | `state_tp` (root + all 10 blocks) | 146.1 µs |
| `H2` FH1 ($T=40\text{ K}, x=[1.0], \rho=4000$) | `helmholtz_terms` | 0.81 µs |
| `H2` FH1 | `pressure_and_derivatives` (analytic) | 4.03 µs |
| `H2` FH1 | `state` (analytic) | 6.48 µs |
| `H2` FH1 | `residual_helmholtz_tvn` (all 10 blocks) | 14.4 µs |
| `H2O/MEOH` (7 sites, $T=380\text{ K}, x=[0.4, 0.6]$) | `state` (analytic + association) | 28.6–41.6 µs |

Compared to the `57f0649` baseline where `PreparedSaftVrMie::pressure_and_derivatives`
and `state` evaluated central finite differences by calling `PreparedSaftVrMie::new`
twice and `helmholtz_terms` seven times, the analytic `Dual1`/`Dual2` kernel
eliminates all repeated temperature preparations and speeds up prepared
pressure/state evaluations while removing finite-difference error.

## Remaining Limitations

- **Mixture phase equilibria**: mixture TP/PH/PS/UV/SV flashes, phase
  envelopes, spinodals, and critical-point solvers for SAFT-VR Mie are not yet
  wired (only pure `saturation_pressure` and single-phase TV/TP paths are
  implemented).
- **`KHS_EOS_SANTOS`**: unimplemented in upstream Fortran
  (`src/saftvrmie_dispersion.f90:3209` calls `stoperror`); omitted here.
- **Truncated `alpha_vdw` for quantum potentials**: returns
  `Error::Unsupported`, matching upstream `calc_alpha_ts` scope.
- **Full DFT convolution integration** (linking FMT weighted densities with
  spatial dispersion convolutions) and **Python/C/C++ FFI bindings** remain
  pending.
