# NIST Multiparameter Equations of State

This module implements native Rust versions of reference-quality NIST multiparameter fundamental equations of state formulated in the reduced Helmholtz free energy $\alpha(\tau, \delta) = \alpha^0(\tau, \delta) + \alpha^{\text{res}}(\tau, \delta)$.

## Implemented Fluids

The implementation embeds exact coefficients for 5 primary fluids from `src/multiparameter_*.f90`:

| Fluid | Formula | Source Reference | $T_c$ (K) | $P_c$ (MPa) | $\rho_c$ ($\text{mol}/\text{m}^3$) |
|---|---|---|---|---|---|
| **Propane** | $\text{C}_3\text{H}_8$ | Lemmon, McLinden & Wagner (2009) | 369.89 | 4.2512 | 5000.0 |
| **R134a** | $\text{C}_2\text{H}_2\text{F}_4$ | Tillner-Roth & Baehr (1994) | 374.21 | 4.05928 | 4978.83 |
| **Normal Hydrogen** | $\text{n-H}_2$ | Leachman et al. (2009) | 33.145 | 1.2964 | 15508.0 |
| **Ortho Hydrogen** | $\text{o-H}_2$ | Leachman et al. (2009) | 33.22 | 1.31065 | 15445.0 |
| **Para Hydrogen** | $\text{p-H}_2$ | Leachman et al. (2009) | 32.938 | 1.2858 | 15538.0 |

## Mathematical Formulation

The thermodynamic state is parameterized by dimensionless variables:
$$\tau = \frac{T_c}{T}, \quad \delta = \frac{\rho}{\rho_c} = \frac{1}{v \rho_c}$$

### 1. Residual Helmholtz Energy ($\alpha^{\text{res}}$)

The residual Helmholtz energy is represented as a sum of three analytical term types:

1. **Polynomial terms**:
   $$\alpha_{\text{pol}}^{\text{res}} = \sum_i n_i \tau^{t_i} \delta^{d_i}$$
2. **Single-exponential terms**:
   $$\alpha_{\text{exp}}^{\text{res}} = \sum_i n_i \tau^{t_i} \delta^{d_i} \exp\left(-\delta^{l_i}\right)$$
3. **Double-exponential (Gaussian bell) terms**:
   $$\alpha_{\text{expexp}}^{\text{res}} = \sum_i n_i \tau^{t_i} \delta^{d_i} \exp\left(-\eta_i (\delta - \epsilon_i)^2 - \beta_i (\tau - \gamma_i)^2\right)$$

All first- and second-order derivatives with respect to $\tau$ and $\delta$ are computed simultaneously using loop-fused analytical formulas:
- $\alpha_\delta^{\text{res}} = \delta \frac{\partial \alpha^{\text{res}}}{\partial \delta}$
- $\alpha_{\delta\delta}^{\text{res}} = \delta^2 \frac{\partial^2 \alpha^{\text{res}}}{\partial \delta^2}$
- $\alpha_\tau^{\text{res}} = \tau \frac{\partial \alpha^{\text{res}}}{\partial \tau}$
- $\alpha_{\tau\tau}^{\text{res}} = \tau^2 \frac{\partial^2 \alpha^{\text{res}}}{\partial \tau^2}$
- $\alpha_{\delta\tau}^{\text{res}} = \delta \tau \frac{\partial^2 \alpha^{\text{res}}}{\partial \delta \partial \tau}$

### 2. Ideal Gas Helmholtz Energy ($\alpha^0$)

The ideal gas contribution is modeled according to the molecular structure:

- **Einstein Summation (Propane, Hydrogen isotopes)**:
  $$\alpha^0 = \ln \delta + c_\tau \ln \tau + a_1 + a_2 \tau + \sum_i v_i \ln\left(1 - e^{-b_i \tau}\right)$$
  where $c_\tau = 3.0$ for propane (nonlinear) and $c_\tau = 1.5$ for hydrogen (diatomic).
- **Power Law (R134a)**:
  $$\alpha^0 = \ln \delta + a_1 + a_2 \tau + a_3 \ln \tau + a_4 \tau^{-0.5} + a_5 \tau^{-0.75}$$

### 3. Thermodynamic State Functions

- **Pressure**:
  $$P(\rho, T) = \rho R T \left(1 + \alpha_\delta^{\text{res}}\right)$$
  $$\frac{\partial P}{\partial \rho} = R T \left(1 + 2 \alpha_\delta^{\text{res}} + \alpha_{\delta\delta}^{\text{res}}\right), \quad \frac{\partial P}{\partial v} = -\rho^2 \frac{\partial P}{\partial \rho}$$
  $$\frac{\partial P}{\partial T} = \rho R \left(1 + \alpha_\delta^{\text{res}} - \alpha_{\delta\tau}^{\text{res}}\right)$$
- **Caloric Properties**:
  $$U = R T \tau \left(\alpha_\tau^0 + \alpha_\tau^{\text{res}}\right)$$
  $$S = R \left[\tau \left(\alpha_\tau^0 + \alpha_\tau^{\text{res}}\right) - \left(\alpha^0 + \alpha^{\text{res}}\right)\right]$$
  $$H = U + P v$$
  $$C_v = -R \left(\alpha_{\tau\tau}^0 + \alpha_{\tau\tau}^{\text{res}}\right)$$
  $$C_p = C_v + R \frac{\left(1 + \alpha_\delta^{\text{res}} - \alpha_{\delta\tau}^{\text{res}}\right)^2}{1 + 2 \alpha_\delta^{\text{res}} + \alpha_{\delta\delta}^{\text{res}}}$$
- **Speed of Sound**:
  $$w = \sqrt{\frac{R T}{M} \left[1 + 2 \alpha_\delta^{\text{res}} + \alpha_{\delta\delta}^{\text{res}} - \frac{\left(1 + \alpha_\delta^{\text{res}} - \alpha_{\delta\tau}^{\text{res}}\right)^2}{\alpha_{\tau\tau}^0 + \alpha_{\tau\tau}^{\text{res}}}\right]}$$
- **Fugacity Coefficient**:
  $$\ln \phi = \alpha^{\text{res}} + Z - 1 - \ln Z \quad \text{where } Z = 1 + \alpha_\delta^{\text{res}}$$

## Performance and Precision

- **Zero Allocations**: `eos.prepare(T)` returns `PreparedMultiparameter`, which caches all temperature-dependent exponential powers and factors. Subsequent calls to `pressure`, `enthalpy`, `entropy`, `speed_of_sound`, and `ln_fugacity` perform zero memory allocations.
- **Differential Parity**: Verified against `.reference/thermopack.dll` across subcritical liquid, dense gas, and supercritical isotherms to within machine precision ($< 10^{-10}$ relative error on pressure, enthalpy, and entropy).
