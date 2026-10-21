# CPA (Cubic-Plus-Association)

This module implements a native Rust version of the Cubic-Plus-Association (CPA) equation of state of Kontogeorgis et al. (1996, 2006), combining a cubic physical equation of state (SRK or PR) with Wertheim's first-order thermodynamic perturbation theory (TPT1) for hydrogen bonding and multi-site association.

## Theory and Formulation

The residual Helmholtz free energy per mole $\alpha^{\text{res}} = A^{\text{res}} / (n R T)$ is partitioned into a physical cubic contribution and a chemical association contribution:

$$\alpha^{\text{res}} = \alpha^{\text{cubic}} + \alpha^{\text{assoc}}$$

### 1. Physical Cubic Contribution ($\alpha^{\text{cubic}}$)

The physical part is modeled with either the Soave-Redlich-Kwong (SRK) or Peng-Robinson (PR) cubic equation of state:

For SRK:
$$\alpha^{\text{cubic}} = -\ln(1 - b \rho) - \frac{a(T)}{b R T} \ln(1 + b \rho)$$

For PR:
$$\alpha^{\text{cubic}} = -\ln(1 - b \rho) - \frac{a(T)}{2\sqrt{2} b R T} \ln\left[\frac{1 + (1+\sqrt{2}) b \rho}{1 + (1-\sqrt{2}) b \rho}\right]$$

where $\rho$ is the molar density ($\text{mol}/\text{m}^3$), $b$ is the covolume ($\text{m}^3/\text{mol}$), and $a(T)$ is the temperature-dependent energy parameter:

$$a_i(T) = a_{0,i} \left[1 + c_{1,i} \left(1 - \sqrt{\frac{T}{T_{c,i}}}\right)\right]^2$$

In mixtures, standard van der Waals one-fluid mixing rules apply:
$$b = \sum_i x_i b_i$$
$$a(T) = \sum_i \sum_j x_i x_j \sqrt{a_i(T) a_j(T)} (1 - k_{ij, a})$$

> [!IMPORTANT]
> **Non-Associating Components in CPA Mixtures**:
> When a non-associating component (e.g. $\text{CO}_2$, $\text{CH}_4$, $\text{N}_2$) is present in a CPA mixture, ThermoPack strictly does **not** use the pure-fluid CPA parameters (which were fitted for pure fluid or with association). Instead, ThermoPack uses the **standard cubic parameters** derived from the fluid's critical properties:
> $$b_i = \Omega_b \frac{R T_{c,i}}{P_{c,i}}, \quad a_{0,i} = \Omega_a \frac{R^2 T_{c,i}^2}{P_{c,i}}, \quad c_{1,i} = m(\omega_i)$$
> For SRK: $\Omega_a \approx 0.427480$, $\Omega_b \approx 0.086640$, $m = 0.480 + 1.574 \omega - 0.176 \omega^2$.
> This behavior is replicated identically in `Cpa::with_options`, yielding exact $10^{-11}$ machine-precision agreement with ThermoPack Fortran.

### 2. Wertheim Association Contribution ($\alpha^{\text{assoc}}$)

According to Wertheim's multi-site association theory:

$$\alpha^{\text{assoc}} = \sum_{i} x_i \sum_{A \in \Gamma_i} \left(\ln X_i^A - \frac{X_i^A}{2} + \frac{1}{2}\right)$$

where $\Gamma_i$ is the set of association sites on molecule $i$, and $X_i^A$ is the fraction of site $A$ on molecule $i$ not bonded to any other site.

The site fractions satisfy the coupled non-linear Wertheim mass-action equations:

$$X_i^A = \frac{1}{1 + \rho \sum_j x_j \sum_{B \in \Gamma_j} X_j^B \Delta^{A_i B_j}}$$

The association strength between site $A$ on molecule $i$ and site $B$ on molecule $j$ is given by:

$$\Delta^{A_i B_j} = g(\eta) \left[\exp\left(\frac{\epsilon^{A_i B_j}}{R T}\right) - 1\right] \beta^{A_i B_j} b_{ij}$$

where:
- $\eta = \frac{b \rho}{4}$ is the reduced packing fraction.
- $g(\eta)$ is the radial distribution function at contact:
  - **Simplified CPA**: $g(\eta) = \frac{1}{1 - 1.9\eta}$
  - **Standard CPA**: $g(\eta) = \frac{1 - \eta/2}{(1 - \eta)^3}$ (Carnahan-Starling)
- $b_{ij} = \frac{b_i + b_j}{2}$ is the cross-covolume.
- $\epsilon^{A_i B_j}$ is the association energy parameter ($\text{J}/\text{mol}$).
- $\beta^{A_i B_j}$ is the dimensionless association volume parameter.

### 3. Cross-Association Combining Rules

For cross-interactions between site $A$ on component $i$ and site $B$ on component $j$ ($i \ne j$):
- Association energy:
  $$\epsilon_{ij} = \frac{\epsilon_i + \epsilon_j}{2} (1 - k_{ij, \epsilon}) \quad \text{or} \quad \sqrt{\epsilon_i \epsilon_j} (1 - k_{ij, \epsilon})$$
- Association volume:
  $$\beta_{ij} = \sqrt{\beta_i \beta_j} (1 - k_{ij, \beta}) \quad \text{or} \quad \frac{\beta_i + \beta_j}{2} (1 - k_{ij, \beta})$$
- **Site Polarity**: Association is allowed if and only if sites have complementary polarities ($p_A p_B \le 0$).

### 4. Fast Closed-Form Analytical Branches

For single-component associating fluids, `PreparedCpa` evaluates closed-form analytical solutions:
- **Scheme 2B** (Alcohols, Glycols: 1 positive and 1 negative site):
  $$X = \frac{\sqrt{1 + 4 \rho \Delta} - 1}{2 \rho \Delta}$$
- **Scheme 4C** (Water: 2 positive and 2 negative sites):
  $$X = \frac{\sqrt{1 + 8 \rho \Delta} - 1}{4 \rho \Delta}$$

For arbitrary multi-component, multi-site mixtures, an allocation-free Newton-Raphson solver with analytical Jacobian and Gaussian elimination resolves all site fractions in 3–6 iterations with quadratic convergence.

## Database and Parameter Sets

Extracted directly from `src/compdatadb.f90` and `src/mixdatadb.f90` via `tools/rust_port/generate_cpa_data.py`:
- **25 pure component CPA records**:
  - Water ($\text{H}_2\text{O}$): 4C scheme.
  - Alcohols ($\text{MEOH}$, $\text{ETOH}$, $\text{1PROPANOL}$, $\text{2PROPANOL}$, $\text{1BUTANOL}$, etc.): 2B scheme.
  - Glycols ($\text{MEG}$, $\text{DEG}$, $\text{TEG}$): 4C scheme.
  - Acids ($\text{ACTACID}$): 1 site dimerizing scheme.
  - Solvents and non-associating fluids ($\text{CO}_2$, $\text{CH}_4$, $\text{C}_2\text{H}_6$, etc.).
- **12 binary interaction records**:
  - $k_{ij, a}$, $k_{ij, \epsilon}$, $k_{ij, \beta}$ parameters and combining rule flags for key systems such as $\text{H}_2\text{O} + \text{CO}_2$, $\text{H}_2\text{O} + \text{ETOH}$, $\text{H}_2\text{O} + \text{MEOH}$, $\text{H}_2\text{O} + \text{MEG}$, $\text{H}_2\text{O} + \text{TEG}$.

## Verification and Test Coverage

Differential qualification against `.reference/thermopack.dll` in `tests/cpa.rs`:
- `test_cpa_pure_h2o_states`: Pressure, compressibility, residual Helmholtz energy, and fugacity coefficients across vapor, dense vapor, and supercritical states ($< 10^{-11}$ error).
- `test_cpa_pure_etoh_states`: Ethanol (2B) state properties across vapor, dense vapor, and gas phases.
- `test_cpa_pure_meoh_states`: Methanol (2B) state properties across multiple densities.
- `test_cpa_binary_h2o_co2`: Water + Carbon Dioxide binary mixture (associating + non-associating fluid) with $k_{ij}$ and induced interaction rules ($< 10^{-10}$ error on pressure and $\ln \phi$).
- `test_cpa_binary_h2o_etoh`: Water + Ethanol cross-associating binary mixture with cross $k_{ij}$ interaction rules.
- `test_cpa_pure_h2o_volume_and_saturation`: Water single-phase vapor/liquid volume roots at 1 atm and vapor-liquid saturation pressure ($P_{\text{sat}} \approx 1.002195 \times 10^5$ Pa at 373.15 K).

## The hard-sphere contact radial distribution function

The original `src/saft_rdf.f90` routes a request for the radial distribution
function of a pair of components to the model of the active thermodynamic
model. The branch of the cubic association model is the hard-sphere contact
function of the mixture, which needs no temperature dependence:

$$g = \frac{1}{1 - 1.9\,\eta} \quad \text{(simplified)}, \qquad g = \frac{1 - \eta/2}{(1 - \eta)^3} \quad \text{(original)},$$

with the packing fraction $\eta = \sum_i n_i b_i / (4 V)$ of the mixture. The
native module `rdf` provides the function and its volume and mole number
derivatives as `rdf::cpa_rdf`; the temperature derivatives vanish, as in the
original.

`tools/rust_port/validate_cpa_rdf.py` compares the native function against the
exact-source library. The reference probe calls the original `master_saft_rdf`
for the active model of the cubic association model.

- 5 cases: the component, the temperature, the volume and the number of moles.
- The compared quantities are the function and its nine derivatives.
- All 5 cases match with a relative 1e-12, and the observed worst relative
  difference is 4.6e-16.

The committed fixture is `tests/fixtures/cpa_rdf.csv`; the test is
`tests/cpa_rdf.rs`. The branches of the remaining models of the original
dispatcher are pending.
