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
  `SingleComponent`, `Additive`, `NonAdditive`, and `AdditiveExact`, plus
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
  potentials.
- **Association (`rust/src/association.rs`, `rust/src/saft_vrmie/association.rs`)**:
  reuses the shared Wertheim TPT1 site solver (`AssocScheme`, `Standard` and
  `Elliot` $\Delta$ combining rules, exact infinite-dilution zero-mole site
  fractions) paired with the Dufal et al. (2015) polynomial Mie association
  integral (`association_rdf_tvn` and `rdf::saftvrmie_rdf`). Controls include
  `has_association`, `set_assoc_params`, `set_delta_comb_rule`, and
  `solve_assoc_x_k`.
- **Quantum extensions**: `new_quantum` and `new_quantum_with_hs_ref`
  (FH0/FH1/FH2, `FeynmanHibbs` and `JaenKahn` variants, `H_PLANCK = 6.626069311e-34`),
  plus `set_mass`, `de_boer_parameter`, `set_de_boer_parameter`, and
  `de_broglie_wavelength`.
- **Model-side FMT/DFT hooks**: `fmt_compatibility`, `dispersion_tvn`,
  `hard_sphere_tvn`, `hard_sphere_diameters`, `hard_sphere_diameter_ij`,
  `effective_sigma_ij`, `effective_epsilon_ij`, `alpha_vdw`, and
  `pair_potential`.

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
  derivatives.
- `PreparedSaftVrMie` caches temperature-dependent Barker-Henderson diameters
  (21-point Gauss-Kronrod quadrature) and effective quantum parameters for fast
  repeated TV evaluations (`helmholtz_terms` via zero-overhead `f64`, analytic
  `pressure_and_derivatives` and `state` via stack-buffered `Dual1`/`Dual2`
  without re-preparing temperature caches, `volume` with `Liquid`/`Vapor`/`MinimumGibbs`,
  and pure `saturation_pressure`).

## Exact-Source Qualification

`tools/rust_port/validate_saftvrmie.py` compiles an isolated exact-source
Fortran probe (`tools/rust_port/saftvrmie_probe.f90`) linked against
`.reference` and qualifies **54 live cases** (**938 scalar comparisons**),
persisted in `tests/fixtures/saftvrmie_tvn.csv` and verified in
`tests/saft_vrmie.rs`:

- **TVN residual, dispersion, hard-sphere, and Dufal RDF** (42 cases covering
  classical, associating `NH3`/`H2O,MEOH`/`ETOH,C1`, Elliot rule, non-unit mole
  scaling, all 5 hard-sphere references, all 5 term masks, cutoff/shift,
  FH0/FH1/FH2, `additive-exact-quantum`, and all 9 independent submodel
  options): worst scaled error **`2.57e-10`** (dominated by the central
  finite-difference check replacing Fortran's defective extra hard-sphere
  $F_{VV}$).
- **Normal TP states** (8 cases covering `CO2,C1` liquid/vapor, `NH3`,
  `H2O,MEOH` standard and Elliot, `ETOH,C1`, `H2` FH1, `H2,NE` FH2 across
  $Z, H^{\text{res}}, H, S^{\text{res}}, S, \ln\boldsymbol{\phi}$, all
  $(T, P, n)$ derivatives, and speed of sound): worst scaled error **`7.21e-12`**.
- **Model-side FMT/DFT hooks** (4 cases covering classical, truncated-shifted,
  and quantum exact/arithmetic-cross diameters, effective $\sigma_{ij}/\epsilon_{ij}$,
  $\alpha_{ij}$, and de Broglie wavelength): worst scaled error **`4.05e-16`**.

Five verified upstream Fortran defects (extra hard-sphere $a_{VV}$,
`calc_saftvrmie_hard_sphere` $F_{TV}$ sign and unguarded optional `F`,
`test_fmt_compatibility` unassigned output, classical `de_Broglie_wavelength`
zero-mass check, and zero-mole $X_k = 1$) are documented in
[baseline-defects.md](baseline-defects.md) and independently verified by finite
differences and thermodynamic identities in `tests/saft_vrmie.rs`.

## Prepared-State Release Performance

Measured on the optimized generic `Scalar`/`Dual1`/`Dual2` kernel (`--release`,
thin LTO, stack-buffered for $n_c \le 8$, symmetric pair loop, single
association solve per state):

| Workload | Operation | Latency |
|---|---|---|
| `CO2/C1` ($T=280\text{ K}, x=[0.4, 0.6], \rho=5000$) | `prepare` | 10.6 µs |
| `CO2/C1` | `helmholtz_terms` | 3.54 µs |
| `CO2/C1` | `pressure_and_derivatives` (analytic) | 12.2 µs |
| `CO2/C1` | `state` (analytic) | 26.8 µs |
| `CO2/C1` | `residual_helmholtz_tvn` (all 10 blocks) | 71.5 µs |
| `CO2/C1` | `state_tp` (root + all 10 blocks) | 146.1 µs |
| `H2` FH1 ($T=40\text{ K}, x=[1.0], \rho=4000$) | `prepare` | 2.66 µs |
| `H2` FH1 | `helmholtz_terms` | 0.81 µs |
| `H2` FH1 | `pressure_and_derivatives` (analytic) | 4.03 µs |
| `H2` FH1 | `state` (analytic) | 6.48 µs |
| `H2` FH1 | `residual_helmholtz_tvn` (all 10 blocks) | 14.4 µs |
| `H2O/MEOH` (7 sites, $T=380\text{ K}, x=[0.4, 0.6]$) | `state` (analytic + association) | 35.1 µs |

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
