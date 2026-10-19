# Lee-Kesler Generalized Corresponding States Equation of State

This module implements a native Rust version of the Lee-Kesler (1975) generalized corresponding states equation of state and the Plöcker et al. (1982) critical property mixing rules.

## Theoretical Background

The Lee-Kesler equation of state predicts the compressibility factor and residual properties of non-polar and slightly polar fluids using a three-parameter corresponding states principle based on critical temperature $T_c$, critical pressure $P_c$, and acentric factor $\omega$:

$$Z = Z^{(0)}(T_r, P_r) + \frac{\omega}{\omega^{(r)}} \left[ Z^{(r)}(T_r, P_r) - Z^{(0)}(T_r, P_r) \right]$$

where:
- $Z^{(0)}$ is the compressibility factor of a simple fluid ($\omega = 0$).
- $Z^{(r)}$ is the compressibility factor of a reference fluid (n-octane with $\omega^{(r)} = 0.3978$).
- $T_r = T / T_{cM}$ is the reduced temperature.
- $P_r = P / P_{cM}$ is the reduced pressure.

Both simple and reference fluids are represented by a modified 12-constant Benedict-Webb-Rubin (BWR) equation of state:

$$Z = 1 + \frac{B}{v_r} + \frac{C}{v_r^2} + \frac{D}{v_r^5} + \frac{c_4}{T_r^3 v_r^2} \left(\beta + \frac{\gamma}{v_r^2}\right) \exp\left(-\frac{\gamma}{v_r^2}\right)$$

where $v_r = \frac{P_{cM} v}{R T_{cM}}$ is the modified reduced volume, and:

$$B = b_1 - \frac{b_2}{T_r} - \frac{b_3}{T_r^2} - \frac{b_4}{T_r^3}$$
$$C = c_1 - \frac{c_2}{T_r} + \frac{c_3}{T_r^3}$$
$$D = d_1 + \frac{d_2}{T_r}$$

## Plöcker Critical Mixing Rules (1982)

For mixtures, pseudocritical properties $T_{cM}, v_{cM}, P_{cM}, \omega_M$ and $Z_{cM}$ are determined via Plöcker's mixing rules with parameter $\eta = 0.25$:

$$v_{cM} = \frac{1}{n^2} \sum_j \sum_k n_j n_k v_{c,jk}$$

$$T_{cM} = \frac{1}{v_{cM}^\eta n^2} \sum_j \sum_k n_j n_k v_{c,jk}^\eta T_{c,jk}$$

$$\omega_M = \frac{1}{n} \sum_j n_j \omega_j$$

$$Z_{cM} = 0.2905 - 0.085 \omega_M$$

$$P_{cM} = R \frac{Z_{cM} T_{cM}}{v_{cM}}$$

Binary interaction parameters modify the cross critical temperature:

$$T_{c,jk} = k_{jk} \sqrt{T_{c,j} T_{c,k}}, \quad v_{c,jk} = \frac{1}{8} \left( v_{c,j}^{1/3} + v_{c,k}^{1/3} \right)^3$$

where $v_{c,i} = R Z_{c,i} T_{c,i} / P_{c,i}$ is the pure component critical volume.

### Compositional Derivatives

For fugacity coefficients $\ln \phi_i = \left(\frac{\partial F}{\partial n_i}\right)_{T, V} - \ln Z$, analytical compositional derivatives of all pseudocritical mixing properties are evaluated:

$$\frac{\partial v_{cM}}{\partial n_i} = \frac{2}{n} \left( \frac{1}{n} \sum_l n_l v_{c,il} - v_{cM} \right)$$

$$\frac{\partial T_{cM}}{\partial n_i} = -\eta \frac{T_{cM}}{v_{cM}} \frac{\partial v_{cM}}{\partial n_i} + \frac{2}{n^2 v_{cM}^\eta} \sum_l n_l v_{c,il}^\eta T_{c,il} - \frac{2 T_{cM}}{n}$$

$$\frac{\partial \omega_M}{\partial n_i} = \frac{\omega_i - \omega_M}{n}$$

$$\frac{\partial Z_{cM}}{\partial n_i} = -0.085 \frac{\partial \omega_M}{\partial n_i}$$

$$\frac{\partial P_{cM}}{\partial n_i} = P_{cM} \left( \frac{1}{Z_{cM}} \frac{\partial Z_{cM}}{\partial n_i} + \frac{1}{T_{cM}} \frac{\partial T_{cM}}{\partial n_i} - \frac{1}{v_{cM}} \frac{\partial v_{cM}}{\partial n_i} \right)$$

## Numerical Implementation & Phase Resolution

- **Newton-Raphson Solver**: Matches ThermoPack's `rootSelection = 3` logic in `src/leekesler.f90`:
  1. For liquid, searches starting at $z_{\text{min}} = P_r v_{r,\text{min}} / T_r$.
  2. For vapor, starts at $z = 1.0$ and steps toward roots.
  3. Uses second-order Householder/Halley acceleration steps: $\delta = \frac{1}{2} f f'' / (f')^2$, $\Delta z = -\frac{f}{f'} (1 + \delta)$ with bounds bracketing $[z_{\text{min}}, z_{\text{max}}]$.
  4. Phase fallback and Gibbs minimization: if a single-phase volume is non-existent, gracefully falls back to the existing phase; when `Phase::MinimumGibbs` is requested, evaluates departure Gibbs energy $G^{\text{dep}} = H^{\text{res}} - T S^{\text{res}}$ to identify the stable root.
- **Allocation-Free Evaluation**: `model.prepare(z)` precomputes all mixture critical properties ($T_{cM}, v_{cM}, P_{cM}, \omega_M, Z_{cM}$) and their $n$-element compositional derivatives once. Subsequent property calls (`compressibility_factor`, `specific_volume`, `enthalpy`, `entropy`, `ln_fugacity_coefficients`) require zero heap allocations.

## Validation

The test suite in `tests/lee_kesler.rs` is automatically generated against the compiled reference library `.reference/thermopack.dll`:
- Tested across pure C1 and C3, binary C1+C2 (default $k_{ij}=1.0$), binary C1+C2 with Plöcker (1982) calibrated interaction parameters ($k_{12}=1.052$), binary C1+C3, and CO2+C1.
- State conditions include compressed liquids, subcritical vapors, supercritical fluids, and high-pressure dense gases.
- Compressibility factors $Z$, specific volumes $v$, total enthalpies $H$, total entropies $S$, and individual component fugacity coefficients $\ln \phi_i$ agree with the Fortran reference to within $10^{-10}$ relative tolerance (machine precision).
