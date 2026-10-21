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

`PreparedPcSaft` provides allocation-free, cache-friendly evaluations:
- `residual_helmholtz_energy(rho)`: $\alpha^{\text{res}}$
- `pressure(rho)`: $P = \rho R T (1 + \rho \frac{\partial \alpha^{\text{res}}}{\partial \rho})$
- `pressure_and_derivatives(rho)`: $P$, $Z$, $\frac{\partial P}{\partial \rho}$, $\frac{\partial P}{\partial V}$, $\frac{\partial P}{\partial T}$
- `state(rho)`: complete state bundle ($P, Z, U^{\text{res}}, H^{\text{res}}, S^{\text{res}}, A^{\text{res}}, G^{\text{res}}, H^{\text{total}}, S^{\text{total}}, C_v, C_p, w, \ln \boldsymbol{\phi}$)
- `volume(P, phase)`: single-phase bounded $(1-\zeta)$ Newton solver from `pc_saft_nonassoc_volume_solver`, guaranteed to converge within physical limits $[10^{-10}, 0.75]$.
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
