# PC-SAFT and sPC-SAFT

This module implements native Rust versions of the Perturbed-Chain Statistical Associating Fluid Theory (PC-SAFT) of Gross & Sadowski (2001) and its simplified variant sPC-SAFT of von Solms et al. (2003).

## Theory and Formulation

The residual Helmholtz energy per mole $\alpha^{\text{res}} = A^{\text{res}} / (n R T)$ is composed of hard-chain and dispersion contributions:

$$\alpha^{\text{res}} = \alpha^{\text{hc}} + \alpha^{\text{disp}}$$

### 1. Hard-Chain Contribution ($\alpha^{\text{hc}}$)

Segments have temperature-dependent diameters $d_i(T)$ according to the Chen-Kreglewski expression:

$$d_i(T) = \sigma_i \left[1 - 0.12 \exp\left(-3 \frac{\epsilon_i / k_B}{T}\right)\right]$$

For simplified PC-SAFT (sPC-SAFT):
An average packing fraction $\eta = \rho \frac{\pi}{6} N_A \sum_i x_i m_i d_i^3$ is used with Carnahan-Starling for both the hard-sphere reference and the contact radial distribution function $g(\eta)$:

$$\alpha^{\text{hs}} = \frac{4\eta - 3\eta^2}{(1-\eta)^2}, \quad g = \frac{1 - \eta/2}{(1-\eta)^3}$$

$$\alpha^{\text{hc}} = \bar{m} \alpha^{\text{hs}} - \sum_i x_i (m_i - 1) \ln g$$

where $\bar{m} = \sum_i x_i m_i$.

For original PC-SAFT:
The hard-sphere mixture is evaluated with the Boublik-Mansoori-Carnahan-Starling-Leland (BMCSL) moments:

$$\zeta_l = \rho \frac{\pi N_A}{6} \sum_i x_i m_i d_i^l, \quad l \in \{0, 1, 2, 3\}$$

$$\alpha^{\text{hs}}_{\text{BMCSL}} = \frac{1}{\zeta_0} \left[ \frac{3 \zeta_1 \zeta_2}{1-\zeta_3} + \frac{\zeta_2^3}{\zeta_3 (1-\zeta_3)^2} + \left(\frac{\zeta_2^3}{\zeta_3^2} - \zeta_0\right) \ln(1-\zeta_3) \right]$$

The radial distribution function at contact for segment $i$ on molecule $i$ is evaluated with $\mu_i = d_i / 2$:

$$g_{ii} = \frac{1}{1-\zeta_3} + 3 \zeta_2 \mu_i \frac{1}{(1-\zeta_3)^2} + 2 \zeta_2^2 \mu_i^2 \frac{1}{(1-\zeta_3)^3}$$

$$\alpha^{\text{hc}} = \bar{m} \alpha^{\text{hs}}_{\text{BMCSL}} - \sum_i x_i (m_i - 1) \ln g_{ii}$$

> [!NOTE]
> For pure fluids ($nc = 1$), $\zeta_l = \zeta_0 d^l$ and $\zeta_2 \mu = \frac{1}{2} \zeta_3 = \frac{1}{2} \eta$. The $\ln(1-\zeta_3)$ term vanishes identically ($\zeta_2^3 / \zeta_3^2 - \zeta_0 = 0$), and $g_{11}$ reduces identically to $(1 - \eta/2)/(1-\eta)^3$. Original PC-SAFT and sPC-SAFT are mathematically identical for all pure components to full floating-point precision ($0.0$ difference verified against the exact reference DLL).

### 2. Dispersion Contribution ($\alpha^{\text{disp}}$)

Both original PC-SAFT and sPC-SAFT share the identical dispersion formulation:

$$\alpha^{\text{disp}} = -\pi \rho N_A \left[ 2 I_1 m_2 \epsilon_1 \sigma^3 + \bar{m} C_1 I_2 m_2 \epsilon_2 \sigma^3 \right]$$

where:
- Universal dispersion polynomials $I_1(\eta) = \sum_{k=0}^6 a_k(\bar{m}) \eta^k$ and $I_2(\eta) = \sum_{k=0}^6 b_k(\bar{m}) \eta^k$.
- The coefficients $a_k(\bar{m})$ and $b_k(\bar{m})$ interpolate universal $7 \times 3$ matrices $a_{mat}$ and $b_{mat}$.
- The compressibility factor $C_1(\eta, \bar{m}) = \left[1 + \bar{m} \frac{8\eta - 2\eta^2}{(1-\eta)^4} + (1-\bar{m}) \frac{20\eta - 27\eta^2 + 12\eta^3 - 2\eta^4}{((1-\eta)(2-\eta))^2}\right]^{-1}$.
- Lorentz-Berthelot combining rules:
  $$\sigma_{ij} = \frac{\sigma_i + \sigma_j}{2}, \quad \frac{\epsilon_{ij}}{k_B} = \sqrt{\frac{\epsilon_i}{k_B} \frac{\epsilon_j}{k_B}} (1 - k_{ij})$$
- Dispersion sums:
  $$m_2 \epsilon_1 \sigma^3 = \frac{1}{T} \sum_{i,j} x_i x_j m_i m_j \left(\frac{\epsilon_{ij}}{k_B}\right) \sigma_{ij}^3$$
  $$m_2 \epsilon_2 \sigma^3 = \frac{1}{T^2} \sum_{i,j} x_i x_j m_i m_j \left(\frac{\epsilon_{ij}}{k_B}\right)^2 \sigma_{ij}^3$$

## Database and Parameters

Extracted from `src/pc_saft_datadb.f90` via `tools/rust_port/generate_pcsaft_data.py`:
- 77 pure component records ($m$, $\sigma$, $\epsilon / k_B$).
- 34 binary interaction $k_{ij}$ records (e.g. $\text{CO}_2 + \text{C}_1$, $\text{CO}_2 + \text{C}_3$, $\text{H}_2\text{S} + \text{alkanes}$, $\text{NH}_3 + \text{H}_2\text{O}$).
- Reference resolution supports `Default` fallback or specific citations (e.g. `Gross_Sadowski2001`, `Tang_Gross2010`, `Aasen2024`).

## Performance and State Properties

`PreparedPcSaft` caches temperature and composition functions. Scalar pressure
and residual-energy queries allocate no memory; state and fugacity results own
output vectors. Construct prepared states with `prepare` or `new`:
- `residual_helmholtz_energy(rho)`: $\alpha^{\text{res}}$
- `pressure(rho)`: $P = \rho R T (1 + \rho \frac{\partial \alpha^{\text{res}}}{\partial \rho})$
- `pressure_and_derivatives(rho)`: $P$, $Z$, $\frac{\partial P}{\partial \rho}$, $\frac{\partial P}{\partial V}$, $\frac{\partial P}{\partial T}$
- `state(rho)`: complete state bundle ($P, Z, U^{\text{res}}, H^{\text{res}}, S^{\text{res}}, A^{\text{res}}, G^{\text{res}}, H^{\text{total}}, S^{\text{total}}, C_v, C_p, w, \ln \boldsymbol{\phi}$)
- `volume(P, phase)`: single-phase bounded $(1-\zeta)$ Newton solver from `pc_saft_nonassoc_volume_solver`, with packing-fraction bounds $[10^{-10}, 0.75]$.
- `saturation_pressure(T)`: pure fluid vapor-liquid coexistence saturation pressure with Clapeyron initialization and damped Newton root finding.

## The radial distribution function

The original evaluates the radial distribution function of the PC-SAFT model
in src/pc_saft_nonassoc.f90:g_pc_saft_TVn, which calls the pair correlation
function of the Boublik hard-sphere model in
src/hardsphere_bmcsl.f90:calc_bmcsl_gij. The port needs three pieces of
that module:

- the hard-sphere diameter of every component with its temperature
  derivatives, src/hardsphere_bmcsl.f90:calc_dhs,
- the packing fraction and the segment fractions of the mixture with their
  derivatives, calc_bmcsl_zeta_and_derivatives, which the native model
  already evaluates for the hard-chain contribution of the association term
  in PreparedPcSaft::bmcsl_zetas,
- the reduced Helmholtz energy of the pair and its first and second
  derivatives with respect to the packing fractions,
  calc_bmcsl_dgij_dzeta, with the chain rule of the original to the
  temperature, the volume and the mole numbers.

All three pieces are ported, in the module rust/src/bmcsl.rs, and the native
function rdf::pcsaft_rdf reproduces the original g_pc_saft_TVn in the volume
representation of the original.

The module of the original publishes the pair correlation function
calc_bmcsl_gij and the packing fractions
calc_bmcsl_zeta_and_derivatives, but keeps the reduced pair Helmholtz
energy calc_bmcsl_dgij_dzeta and its derivatives with respect to the
packing fractions private. The reference probe therefore calls the published
routines of the original, and the qualification compares the pair and its
temperature, volume and mole number derivatives.

The qualification is tools/rust_port/validate_pcsaft_rdf.py. The reference
probe tools/rust_port/reference_probe.f90 calls calc_dhs,
calc_bmcsl_zeta_and_derivatives and calc_bmcsl_gij of the original for the
active model of the PC-SAFT model.

- 6 cases: the components, the temperature, the volume, the number of moles
  and the components of the pair.
- The compared quantities are the function and its nine derivatives.
- All 6 cases match with a relative 1e-12, and the observed worst relative
  difference is 4.1e-16.

The committed fixture is tests/fixtures/pcsaft_rdf.csv; the test is
tests/pcsaft_rdf.rs.

## The pair correlation function of the hard-sphere model

The module rust/src/bmcsl.rs holds the reduced pair Helmholtz energy of the
additive hard-sphere model of Boublik and its first and second derivatives
with respect to the packing fractions, in the function
bmcsl::boublik_pair_correlation. The expression and its derivatives follow
src/hardsphere_bmcsl.f90:calc_bmcsl_dgij_dzeta.

The derivatives of the kernel are verified against a central difference of the
expression, because the original keeps the routine private and delivers the
pair and its temperature, volume and mole number derivatives only. The chain
rule of the original is verified against the exact-source library, as
described above.

The module holds the chain rule as well. The function
bmcsl::packing_fractions evaluates the moments zeta_0 to zeta_3 of the number
density of the mixture and their derivatives with respect to the temperature,
the volume and the mole numbers, in the order of the original
src/hardsphere_bmcsl.f90:calc_bmcsl_zeta_and_derivatives. The function
bmcsl::pair_correlation applies the chain rule to the pair of components and
returns the pair correlation function with its first and second derivatives,
in the order of the original src/hardsphere_bmcsl.f90:calc_bmcsl_gij. The
reduced diameter of the pair and its temperature derivatives are the ones of
bmcsl::HsDiameters::reduced_diameter.

The chain rule is verified against central differences of the expression, with
relative steps of the temperature, the volume and the mole numbers. The second
derivatives of the volume and the mole numbers need a larger step, because the
packing fractions scale with the inverse of the volume.

## The reduced Helmholtz energy of the hard-sphere model

The module rust/src/bmcsl.rs holds the reduced Helmholtz energy of the
additive hard-sphere model of Boublik and its first and second derivatives
with respect to the packing fractions, in the function
bmcsl::boublik_helmholtz_energy. The expression and its derivatives follow
src/hardsphere_bmcsl.f90:calc_bmcsl_dalpha_dzeta, which the original keeps
private.

The function bmcsl::hard_sphere_energy applies the chain rule of the original
src/hardsphere_bmcsl.f90:calc_bmcsl_helmholtzenergy to the temperature, the
volume and the mole numbers. The native function pc_saft::pcsaft_hs_energy
delivers the energy and its derivatives for the model of the PC-SAFT model.

The qualification is tools/rust_port/validate_pcsaft_hs.py. The reference
probe tools/rust_port/reference_probe.f90 calls calc_dhs,
calc_bmcsl_zeta_and_derivatives and calc_bmcsl_helmholtzenergy of the original
for the active model of the PC-SAFT model.

- 6 cases: the components, the temperature, the volume and the number of
  moles.
- The compared quantities are the energy and its nine derivatives.
- All 6 cases match with a relative 1e-12, and the observed worst relative
  difference is 7.7e-15.

The committed fixture is tests/fixtures/pcsaft_hs.csv; the test is
tests/pcsaft_hs.rs.

## The chain term of the PC-SAFT model

The chain term of the PC-SAFT model is the sum over the components with more
than one segment of the logarithm of the pair correlation function of the
component with itself. The original evaluates it in
src/pc_saft_nonassoc.f90:F_chain_PC_SAFT_TVn, which calls the published
routine calc_bmcsl_lngij of the hard-sphere module.

The native function pc_saft::pcsaft_chain_energy reproduces the original. The
logarithm of the pair correlation function and its derivatives follow the
original src/hardsphere_bmcsl.f90:calc_bmcsl_lngij and are the function
bmcsl::PairCorrelation::logarithm of the native module.

The qualification is tools/rust_port/validate_pcsaft_chain.py. The reference
probe tools/rust_port/reference_probe.f90 calls calc_dhs,
calc_bmcsl_zeta_and_derivatives and F_chain_PC_SAFT_TVn of the original for
the active model of the PC-SAFT model.

- 6 cases: the components, the temperature, the volume and the number of
  moles.
- The compared quantities are the energy and its nine derivatives.
- All 6 cases match with a relative 1e-12, and the observed worst relative
  difference is 3.2e-16.

The committed fixture is tests/fixtures/pcsaft_chain.csv; the test is
tests/pcsaft_chain.rs.

The native model of the PC-SAFT model still evaluates the temperature
derivative of the hard-chain contribution with the expressions of the
simplified PC-SAFT model, which are exact for a pure fluid only. The port of
the chain term and of the hard-sphere energy is the first piece of the work
that replaces them for a mixture.

## The state of the PC-SAFT model

The original evaluates the reduced residual Helmholtz energy of the PC-SAFT
model and its derivatives in src/pc_saft_nonassoc.f90:F_PC_SAFT_TVn. The
native model works in the density representation, so the qualification
compares the quantities after the change of variables:

- the reduced residual Helmholtz energy of the mixture is F/sumn,
- the derivative with respect to the density is -F_V/rho^2,
- the derivative with respect to the temperature is F_T/sumn.

The qualification is tools/rust_port/validate_pcsaft_state.py. The reference
probe tools/rust_port/reference_probe.f90 calls F_PC_SAFT_TVn of the original
for the active model of the PC-SAFT model.

- 6 cases: the components, the temperature, the volume and the number of
  moles.
- The compared quantities are the reduced residual Helmholtz energy of the
  mixture, its density derivative, its temperature derivative, its second
  density derivative and its mixed density and temperature derivative.
- All 6 cases match with a relative 1e-12, and the observed worst relative
  difference is 3.5e-13.

The committed fixture is tests/fixtures/pcsaft_state.csv; the test is
tests/pcsaft_state.rs.

The qualification found a deviation of the temperature derivative of a
mixture of about one per mille, because the native model evaluated the
hard-chain contribution with the expressions of the simplified PC-SAFT model.
The model now applies the chain rule of the hard-sphere model of Boublik for
a mixture, and the qualification matches again. The pure fluid keeps the
expressions of the simplified model, which are exact for a pure fluid, and
which are cheaper.

The chain rule needs the packing fractions of the mixture and the pair
correlation function of every component with itself. The evaluation of the
temperature derivatives is free of an allocation, through the functions
bmcsl::packing_moments, bmcsl::hard_sphere_energy_t_derivatives and
bmcsl::pair_correlation_t_derivatives of the native module.

The pressure derivatives of the model are analytic. The density and the
temperature derivatives of the hard-chain contribution follow the chain rule
of the hard-sphere model of Boublik for a mixture, and the packing fraction
of the simplified model otherwise. The mixed density and temperature
derivative of the dispersion term follows the chain rule of the packing
fraction and of the mean segment numbers, and needed a correction: the
expression of the native model missed the term of the derivative of the
coefficient of the dispersion term with respect to the packing fraction. The
derivative of the pressure with respect to the density follows from the
compressibility factor, and the derivative with respect to the volume from
the derivative with respect to the density.

The [analytic first mole derivative qualification](pc-saft-mole-derivatives.md)
covers mixture fugacity in both variants, source failures, missing ideal-data
construction and the remaining derivative scope.

## Association, Multipolar PCP-SAFT, and Complete 10-Block Derivatives

The PC-SAFT family now covers all four ThermoPack variants (`PcSaft`, `SPcSaft`,
`PcpSaft`, `SPcpSaft`):

- **Complete parameter database** (`tools/rust_port/generate_pcsaft_data.py` ->
  `rust/src/pc_saft_data.rs`): extracts all 77 pure records including `eos`,
  `eps`, `beta`, `assoc_scheme`, `mu`, `Q` (with `compdatadb.f90`
  `mu_dipole`/`q_quadrupole` fallback matching `pc_saft_parameters.f90:191-192`
  and `"Gross2005/Gross2006"` polar default reference override matching
  `eoslibinit.f90:983-985`) and all 34 binary records including `eps_comb_rule`
  and `beta_comb_rule` (`DefaultComb`, `AriComb`, `GeoComb`).
- **Shared Wertheim TPT1 association engine** (`rust/src/association.rs`):
  reusable across PC-SAFT, SAFT-VR Mie, and CPA; supports all 11 association
  schemes (`NoAssoc` through `Scheme4C` plus `Scheme1ea`), CR-1 (`Standard`)
  and `Elliot` combining rules, unweighted Michelsen Newton site solver with
  successive-substitution fallback, exact infinite-dilution `X_k` for zero-mole
  species, and full 10-block `(T, V, n)` Helmholtz derivatives (`F..F_nn`).
- **Multipolar perturbation module** (`rust/src/multipol.rs`): implements
  Gross (2005) quadrupole-quadrupole, Gross & Vrabec (2006) dipole-dipole, and
  Vrabec & Gross (2008) dipole-quadrupole Padé contributions with exact
  10-block `(T, V, n)` derivatives via `crate::hyperdual::HyperDual` and
  fast-skip when no polar terms are active.
- **Unit-mole Euler scaling**: `residual_helmholtz_tv`, `assoc_residual_helmholtz_tv`,
  `multipol_residual_helmholtz_tv`, and `state_tp` evaluate in unit-mole
  `(v_u = V/N, x = n/N)` coordinates and apply exact Euler homogeneity scaling,
  preserving scaling invariance down to `N = 1e-140` mol while returning an
  explicit error on genuine floating-point overflow (`N = 1e-310` mol).
- **Exact-source qualification** (`tools/rust_port/validate_pcsaft.py` ->
  `tests/fixtures/pcsaft_full.csv`, `tests/pc_saft_assoc_polar.rs`):
  35 exact-source cases across all four variants (non-associating, associating
  pure and binary/ternary mixtures including `NH3+H2O` `ariComb`/`geoComb`,
  polar `CO2`/`ACETONE`/`H2S`, and polar+associating states) compare all 10
  `(T, V, n)` derivative blocks from `query_pcsaft_full_fres`, solved `X_k`,
  and complete TP state derivatives (`Z, v, H_res, H, S_res, S, ln_phi` with
  `(T, P, n)` derivatives and speed of sound), matching with worst unit-aware
  scaled error `3.55e-10`. Each `PCSAFT_FULL` fixture row evaluates the TV
  residual Helmholtz kernel on raw extensive `moles` and the TP columns on
  normalized unit-mole `x = moles / sum(moles)` matching ThermoPack's TP
  oracle; direct raw-amount non-unit-mole `state_tp` scaling (including
  `N = 1e-140` mol) is tested directly in `tests/pc_saft_assoc_polar.rs`.
  Simplified `rdf::spcsaft_rdf` is qualified in 3
  cases against `query_spcsaft_rdf` with worst scaled error `1.16e-14`.
  Use `PcSaft::new_nonassoc` when explicitly evaluating the non-associating
  `F_PC_SAFT_TVn` / `F_sPC_SAFT_TVn` kernels on associating fluids.
