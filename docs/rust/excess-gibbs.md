# Excess Gibbs / Huron-Vidal Mixing Rules

This module implements native Rust versions of advanced $g^E$-based cubic mixing rules: **Huron-Vidal 1 (HV1)**, **Huron-Vidal 2 (HV2)**, and **NRTL**, integrated directly into the `Cubic` equation of state framework.

## Theory and Formulation

Conventional van der Waals one-fluid quadratic mixing rules ($a_{mix} = \sum_i \sum_j x_i x_j a_{ij}$) are known to fail for strongly non-ideal, polar, and asymmetric mixtures (e.g., water + hydrocarbons, alcohols, glycols). The Huron-Vidal mixing rule equates the excess Helmholtz free energy of the equation of state at infinite pressure with the excess Gibbs free energy from a liquid activity coefficient model ($g_\infty^E$):

$$a_{mix}(T, \mathbf{x}) = b_{mix} \left( \sum_i x_i \frac{a_i(T)}{b_i} - \frac{g_\infty^E(T, \mathbf{x})}{C^*} \right)$$

where:
- $b_{mix} = \sum_i x_i b_i$ is the linear mixture covolume.
- $a_i(T)$ and $b_i$ are the pure-component energy and covolume parameters.
- $C^*$ is an EOS-dependent integration constant derived from the infinite-pressure limit of the Helmholtz free energy:
  - **SRK**: $C^* = \ln 2 \approx 0.6931471805599453$
  - **PR**: $C^* = \frac{1}{\sqrt{2}} \ln(1 + \sqrt{2}) \approx 0.6232252401402305$

### 1. Excess Gibbs Free Energy Models ($g_\infty^E$)

In ThermoPack, $g_\infty^E$ is evaluated using the **NRTL** (Non-Random Two-Liquid) local composition formulation:

$$g_\infty^E(T, \mathbf{x}) = R T \sum_i x_i \frac{\sum_j x_j G_{ji} \tau_{ji}}{\sum_k x_k G_{ki}} = R T \sum_i x_i \frac{N_i}{D_i}$$

where:
- $\tau_{ji}$ is the dimensionless interaction energy parameter:
  - **HV1 / NRTL (constant / linear)**: $\tau_{ji} = \frac{a_{ji}^{GE}}{T} + b_{ji}^{GE}$
  - **HV2 (quadratic polynomial)**: $\tau_{ji} = \frac{a_{ji}^{GE}}{T} + b_{ji}^{GE} + c_{ji}^{GE} T$
  - **Maribo-Mogensen correlation**: $\tau_{ji} = \frac{a_{ji}^{GE}}{T} \left( \frac{T - b_{ji}^{GE}}{c_{ji}^{GE}} \right)$
- $G_{ji} = \exp(-\alpha_{ji} \tau_{ji})$ is the local composition Boltzmann factor, with non-randomness parameter $\alpha_{ji} = \alpha_{ij}$.
- $D_i = \sum_k x_k G_{ki}$ and $N_i = \sum_j x_j G_{ji} \tau_{ji}$.

### 2. Analytical Temperature and Composition Derivatives

To preserve the zero-allocation performance and exact thermodynamic consistency of `PreparedState`, all derivatives are evaluated analytically:

1. **Temperature Derivatives**:
   $$\frac{\partial \tau_{ji}}{\partial T} = -\frac{a_{ji}^{GE}}{T^2} + c_{ji}^{GE}, \quad \frac{\partial^2 \tau_{ji}}{\partial T^2} = \frac{2 a_{ji}^{GE}}{T^3}$$
   $$\frac{\partial G_{ji}}{\partial T} = -\alpha_{ji} G_{ji} \frac{\partial \tau_{ji}}{\partial T}$$
   $$\frac{\partial g_\infty^E}{\partial T} = \frac{g_\infty^E}{T} + R T \sum_i x_i \left[ \frac{\partial N_i / \partial T}{D_i} - \frac{N_i \partial D_i / \partial T}{D_i^2} \right]$$

2. **Mole Number Derivatives**:
   Because $n g_\infty^E$ is homogeneous of degree 1 in mole numbers $\mathbf{n}$:
   $$\frac{\partial (n g_\infty^E)}{\partial n_k} = \mu_k^{E,\infty} = R T \left[ \frac{N_k}{D_k} + \sum_i x_i \frac{G_{ki}(\tau_{ki} D_i - N_i)}{D_i^2} \right]$$

3. **Cubic Mixture Derivatives**:
   $$\frac{\partial a_{mix}}{\partial T} = b_{mix} \left( \sum_i x_i \frac{\partial a_i / \partial T}{b_i} - \frac{\partial g_\infty^E / \partial T}{C^*} \right)$$
   $$\frac{\partial^2 a_{mix}}{\partial T^2} = b_{mix} \left( \sum_i x_i \frac{\partial^2 a_i / \partial T^2}{b_i} - \frac{\partial^2 g_\infty^E / \partial T^2}{C^*} \right)$$
   $$\frac{\partial a_{mix}}{\partial n_k} = b_k \left( \sum_i x_i \frac{a_i}{b_i} - \frac{g_\infty^E}{C^*} \right) + b_{mix} \left( \frac{a_k}{b_k} - \frac{\partial (n g_\infty^E) / \partial n_k}{C^*} \right)$$

### 4. The Second Composition Derivatives

The Huron-Vidal and the NRTL models deliver the second composition
derivatives of the excess Gibbs energy,
$$\frac{\partial^2 g_\infty^E}{\partial n_i \partial n_j},$$
which the Wong-Sandler mixing rule consumes in the sum `D` of its appendix
A9. They are evaluated analytically in the same kernel and are exposed
together with the remaining derivatives in the `ExcessGibbsJet` type and in
the function `excess_gibbs::excess_gibbs_jet`.

## Architecture and Database

- `rust/src/excess_gibbs.rs`: Core analytical $g_\infty^E$ equations, exact first and second derivatives, and cubic parameter combination.
- `rust/src/excess_gibbs_data.rs`: All 19 source-derived binary interaction parameters extracted from `src/mixdatadb.f90` via `tools/rust_port/generate_excess_gibbs_data.py`.
- `rust/src/cubic.rs`: Integration into `CubicKind::PengRobinson`, `CubicKind::SoaveRedlichKwong`, and `PreparedState`.
- `rust/src/wong_sandler.rs`: The Wong-Sandler mixing rule, which consumes the excess Gibbs energy of the infinite pressure limit; see the [Wong-Sandler notes](wong-sandler.md).
- `tests/excess_gibbs.rs`: Automated verification against `.reference/thermopack.dll` across SRK, PR, HV1, HV2, and NRTL.
- `tests/wong_sandler.rs`: Automated verification of the mixture coefficients of the Wong-Sandler rule against `.reference/thermopack.dll`.

## Performance Characteristics

- **Allocation-Free Hot Path**: When calling `eos.prepare(T, z)`, all composition- and temperature-dependent mixing quantities ($a, a_T, a_{TT}, a_{n_i}, a_{n_i T}, a_{n_i n_j}$) are precomputed into contiguous stack/slice storage. Pressure scans and volume roots execute with zero heap allocations.
- **Precision**: Validated against ThermoPack reference DLL to within machine precision ($< 10^{-10}$ relative error on pressure and fugacity coefficients).
