# SAFT-VR Mie and SAFT-VRQ Mie

This module implements native Rust versions of the SAFT-VR Mie equation of state (Lafitte et al., J. Chem. Phys. 139, 154504, 2013) and its quantum extension SAFT-VRQ Mie (Aasen et al., 2019, 2020) for quantum fluids ($\text{H}_2, \text{He}, \text{Ne}$).

## Theory and Formulation

The residual Helmholtz energy per mole of mixture is given by:

$$a^{\text{res}} = \frac{A^{\text{res}}}{n R T} = \bar{m} a^{\text{mono}} + a^{\text{chain}}$$

where $\bar{m} = \sum_i x_i m_i$ is the average segment number per mole, $a^{\text{mono}}$ is the monomer residual Helmholtz energy per segment, and $a^{\text{chain}}$ is the chain formation contribution per mole of mixture.

### 1. Monomer Contribution ($a^{\text{mono}}$)

The monomer fluid consists of spherical segments interacting via the Mie $(\lambda_r, \lambda_a)$ pair potential:

$$\phi_{\text{Mie}}(r) = \mathcal{C} \epsilon \left[ \left(\frac{\sigma}{r}\right)^{\lambda_r} - \left(\frac{\sigma}{r}\right)^{\lambda_a} \right], \quad \mathcal{C} = \left(\frac{\lambda_r}{\lambda_r - \lambda_a}\right) \left(\frac{\lambda_r}{\lambda_a}\right)^{\frac{\lambda_a}{\lambda_r - \lambda_a}}$$

The monomer residual energy is expanded via high-temperature perturbation theory up to third order:

$$a^{\text{mono}} = a^{\text{hs}} + \beta a_1 + \beta^2 a_2 + \beta^3 a_3$$

where $\beta = 1 / (k_B T)$ or $1 / T$ when energy parameters are expressed in Kelvin.

#### Hard-Sphere Reference ($a^{\text{hs}}$)
Evaluated using the BMCSL hard-sphere mixture moments:

$$\zeta_l = \rho \frac{\pi N_A}{6} \sum_i x_i m_i d_{ii}^l, \quad l \in \{0, 1, 2, 3\}$$

$$a^{\text{hs}} = \frac{1}{\zeta_0} \left[ \frac{3 \zeta_1 \zeta_2}{1 - \zeta_3} + \frac{\zeta_2^3}{\zeta_3 (1 - \zeta_3)^2} + \left(\frac{\zeta_2^3}{\zeta_3^2} - \zeta_0\right) \ln(1 - \zeta_3) \right]$$

Temperature-dependent segment diameters $d_{ij}(T)$ are calculated via Barker-Henderson numerical quadrature:

$$d_{ij}(T) = \int_0^{\sigma_{ij}} \left[1 - \exp\left(-\frac{\phi_{ij}(r)}{k_B T}\right)\right] dr$$

computed allocation-free via 32-point Gauss-Legendre quadrature. For quantum fluids (SAFT-VRQ Mie), Feynman-Hibbs quantum effective potentials account for non-classical segment swelling and dispersion modification.

#### Dispersion Terms ($a_1, a_2, a_3$)
- **First order ($a_1$)**: Sutherland reference integral $a_{1,s}(\lambda)$ combined with Barker-Henderson core correction $B(\lambda, x_0)$ with $x_0 = \sigma / d$:
  $$\tilde{a}_1 = \mathcal{C} \left[ x_0^{\lambda_a} (a_{s, a} + B_a) - x_0^{\lambda_r} (a_{s, r} + B_r) \right]$$
  Mixture cross-interactions are scaled by the pair packing fraction $\eta_{ij} = \frac{\pi N_A \bar{m} \rho d_{ij}^3}{6}$.
- **Second order ($a_2$)**: Macroscopic Carnahan-Starling compressibility $K_{\text{hs}}(\eta)$ and fluctuation factor $\chi$:
  $$a_2 = \frac{1}{2} K_{\text{hs}} (1 + \chi) \epsilon \mathcal{C}^2 \eta_{ij} \tilde{a}_{22}$$
- **Third order ($a_3$)**: Empirical exponential formulation parameterized by $\alpha$:
  $$a_3 = -\epsilon^3 f_4(\alpha) \bar{\zeta} \exp(f_5(\alpha) \bar{\zeta} + f_6(\alpha) \bar{\zeta}^2)$$

### 2. Chain Formation ($a^{\text{chain}}$)

Covalent bonding of segments into chains gives:

$$a^{\text{chain}} = -\sum_i x_i (m_i - 1) \ln g_i(\sigma)$$

The contact radial distribution function $g(\sigma)$ follows the Lafitte expansion:

$$g(\sigma) = g_d \exp\left(\frac{\beta \epsilon g_1 + \beta^2 \epsilon^2 g_2}{g_d}\right)$$

where:
- $g_d = \exp(k_0 + k_1 x_0 + k_2 x_0^2 + k_3 x_0^3)$ from Boublik-Carnahan-Starling hard-sphere contact RDF.
- $g_1 = g_{11} + g_{12}$:
  $$g_{11} = \frac{1}{4\epsilon} (\tilde{a}_1 + \eta \tilde{a}_{1, \eta})$$
  $$g_{12} = \frac{\mathcal{C}}{12\epsilon} \left[ -\lambda_a x_0^{\lambda_a} (a_{s, a} + B_a) + \lambda_r x_0^{\lambda_r} (a_{s, r} + B_r) \right]$$
- $g_2 = (1 + \gamma_c)(g_{21} + g_{22})$:
  $$g_{21} = \frac{\mathcal{C}^2}{8\epsilon} \left[ K_{\text{hs}}(\eta) (\tilde{a}_{22} + \eta \tilde{a}_{22, \eta}) + \eta \frac{\partial K_{\text{hs}}}{\partial \eta} \tilde{a}_{22} \right]$$
  $$g_{22} = K_{\text{hs}} \frac{\mathcal{C}^2}{12\epsilon} \left[ -\lambda_a x_0^{2\lambda_a} (a_{s, 2a} + B_{2a}) + (\lambda_a + \lambda_r) x_0^{\lambda_a + \lambda_r} (a_{s, ar} + B_{ar}) - \lambda_r x_0^{2\lambda_r} (a_{s, 2r} + B_{2r}) \right]$$

> [!NOTE]
> The analytical derivation of $g_{21}$ incorporates both the direct derivative $\tilde{a}_{22, \eta}$ and the derivative of the Carnahan-Starling compressibility factor $\frac{\partial K_{\text{hs}}}{\partial \eta}$. This achieves machine-epsilon agreement ($5.55 \times 10^{-17}$) against the ThermoPack Fortran reference kernel (`calcG21_if`).

## Database and Parameters

Extracted from `src/saftvrmie_datadb.f90` and `saftvrqmie_datadb.f90`:
- **104 pure fluid records**: $m, \sigma, \epsilon/k_B, \lambda_r, \lambda_a$ and Feynman-Hibbs mass/de Broglie parameters for quantum fluids.
- Binary interaction parameter support for $k_{ij}$ and $l_{ij}$ (combining rules for energy and diameter).

## Performance and State Properties

`PreparedSaftVrMie` provides allocation-free, cache-friendly evaluations:
- `residual_helmholtz_energy(rho)`: $a^{\text{res}}$
- `pressure_and_derivatives(rho)`: $P, Z, \frac{\partial P}{\partial \rho}, \frac{\partial P}{\partial V}, \frac{\partial P}{\partial T}$
- `state(rho)`: complete thermodynamic state bundle ($P, Z, U^{\text{res}}, H^{\text{res}}, S^{\text{res}}, A^{\text{res}}, G^{\text{res}}, H^{\text{total}}, S^{\text{total}}, C_v, C_p, w, \ln \boldsymbol{\phi}$)
- `volume(P, phase)`: robust bounded root solver for liquid, vapor, and minimum Gibbs energy phases.
- `saturation_pressure(T)`: vapor-liquid coexistence with equal chemical potentials and Clapeyron initialization.
- Exact differential fugacity coefficients:
  $$\ln \phi_i = \left(\frac{\partial F^{\text{res}}}{\partial n_i}\right)_{T, V} - \ln Z$$
  where $F^{\text{res}} = \frac{n_{\text{tot}} A^{\text{res}}}{R T}$, agreeing with Fortran reference DLL down to $10^{-10}$ relative difference.
