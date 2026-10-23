# Pure multiparameter equations of state

`Multiparameter` implements the five NIST fits (C3, R134A, N-H2, O-H2 and
P-H2), the full Lennard-Jones fit and the truncated-shifted LJ fit from
`src/multiparameter_*.f90`. Both LJ fits have source parameter sets for Ne,
Ar, Kr and Xe. The public reduced `init_LJTS` initializer is also available
before noble-gas SI scaling.

[Pure coexistence](meos-saturation.md) adds saturation pressure, cold inverse
temperature, both phase volumes and analytic Clapeyron slopes for all 14
variants. Its direct public-source comparisons and independently checked
recovery states have separate fixtures and explicit critical/branch limits.

```rust
use thermopack_rs::{Multiparameter, LennardJonesKind};
let eos = Multiparameter::lennard_jones("AR", LennardJonesKind::Full)?;
let isotherm = eos.prepare(150.0)?;
let state = isotherm.state(1e-3)?;
assert!(state.pressure.value > 0.0);
# Ok::<(), thermopack_rs::Error>(())
```

`new` retains the NIST fluid aliases. `lennard_jones` accepts a noble-gas
identifier and `Full` or `TruncatedShifted`. Each model owns its metadata and
scaling; coefficient arrays remain static and shared. Consequently `record()`
now borrows model-owned metadata instead of returning a static reference.
Copying a model copies its metadata without allocating coefficient storage.
Unsupported LJ components are native errors; the source constructor prints a
message and continues with unset parameters, which is not a defined result.

## Helmholtz kernels and units

The polynomial, single-exponential and Gaussian residual terms use
`tau=Tc/T` and `delta=rho/rhoc`. `HelmholtzJet` contains alpha and scaled
derivatives tau*alpha_tau, tau^2*alpha_tautau, delta*alpha_delta,
delta^2*alpha_deltadelta and tau*delta*alpha_taudelta. The ideal kernel uses
Einstein terms, R134A power terms, or the monatomic LJ specialization.
The LJ ideal coefficients are ordered as in the source expression, including
LJTS entropy alignment at T*=0.8 and p*=0.001 before SI scaling.

The constructor scales Tc/triple T with epsilon/k and rhoc with
1/(N_A*sigma^3), calculates Pc from the EOS at the stored Tc/rhoc, and
calculates triple densities from the ancillary expressions. LJ limits, sigma,
epsilon/k and the original mass field are exposed separately. The source
limits are metadata, not enforced phase-stability or validity guarantees.
Noble-gas models use SI K, Pa, mol/m3, m3/mol and molar energies. The reduced
LJTS constructor keeps R=1 and M=1: it uses the original reduced input/output
units and leaves the original zero Pc/triple-density metadata unchanged.

`prepare(T)` uses fixed-size inline caches for temperature powers, Gaussian
exponents/derivatives and the ideal temperature jet. Preparation and state
queries allocate no heap storage. Each ideal oscillator shares one negative
exponential across its energy and first/second derivatives. Capacity checks
cover every compiled fit;
adding a larger fit requires updating those capacities. Residual terms retain
source family sums and the combined Gaussian exponent. These choices matter
for Cp near criticality and avoid the large amplification of ulp-level
rounding that occurred with a single accumulated sum and split exponentials.

`state(v)` shares one residual jet and the cached ideal jet across pressure
and its T/v derivatives, U, H, S, A, G, Cv, Cp, fugacity and sound. Scalar
methods remain available. `residual_helmholtz(V,n)` and `ideal_helmholtz(V,n)`
return reduced extensive F=A/(RT) and all ten first/second T/V/n quantities.
Total V and positive mole amount are explicit; the internal molar volume is
V/n. F has units mol, while alpha and the scaled jets are dimensionless.

A state can be an unstable or extrapolated homogeneous state. Negative
pressure/Cp are retained; this API does not choose a phase. Non-finite
required properties are errors. Optional fugacity is absent when Z<=0, and
optional sound is absent when its squared value is negative or non-finite.
LJ ancillary densities return one above Tc as the source does. C3/R134A
supercritical ancillary formulas are non-finite in the source and produce
native errors. Hydrogen ancillary densities retain their triple-point ratios.
Only liquid/vapor ancillary requests are supported.

## Source sound convention

`constructor_LJ` stores `compdb%mw` (g/mol) in `meos%molarMass`, whose sound
routine expects kg/mol. The direct source sound result is therefore smaller
than physical SI sound by sqrt(1000). `baseline_speed_of_sound` and the bulk
state's corresponding field preserve that direct source result. Normal
`speed_of_sound`, its state field and `molar_mass()` use SI kg/mol.
The qualification explicitly compares the baseline sound; its physical
counterpart is checked independently from the pressure/Cp/Cv identity.
This convention concerns the direct meos method, not a claim about every
higher-level sound wrapper.

## Density iteration

The original `multiparameter_base::densitySolver` solves for the density at a
specified pressure and phase flag. The native `baseline_density` and
`baseline_density_with_phase_code` reproduce its structure: the phase
switching on an unstable branch, the raw source phase flags, the curvature
guard, the pressure and density stopping tests, the source `ierr` at
iteration 200 and the two-switch rule that turns off the robustness measures
above the triple point. The execution is bounded at 10000 steps and 1024
liquid seed doublings, and an undefined or non-finite state is an error
rather than a last iterate.

The original density criterion needs the density change to fall below
`releps_rho*rhoOld`, which is 2.2e-10 relative, while the pressure criterion
is met. The native residual kernel agrees with the source only to about 1e-14
relative, and the Newton step divides that pressure difference by dP/drho, so
a state close to the critical point cannot reach the original density
criterion. The native model therefore accepts a stalled iteration, with a
closed pressure on a mechanically stable branch, and reports it through
`relaxed_convergence`. A contracting iteration reaches the original criterion
first and never takes that path.

The qualification covers 1694 attempted states over the fourteen model
configurations, with reduced temperatures from 0.6 to 3 and reduced pressures
from 1e-4 to 10, plus the phase flags that the original does not support and
selected below-triple extrapolations.

* 1669 states return a density and match the exact-source library. The
  density tolerance adds the conditioning term `1e-14*|p|/|dP/drho|`, which
  is the uncertainty that the kernel agreement allows at a state with a small
  dP/drho.
* 25 states terminate in the source because its ancillary rejects the phase
  flag, and the native model reports the same error.
* No state is censored, and the relaxed criterion is used in 11 of the 1669
  states, all of them at the reduced critical point.

The committed fixture is `tests/fixtures/meos_density.csv`; the test is
`tests/meos_density.rs`.

## Qualification and numerical limits

The source matrix contains 5082 homogeneous states over fourteen model
configurations: five NIST fluids, eight scaled LJ fits and reduced LJTS.
Eleven T/Tc values from 0.45 to 8, eleven rho/rhoc values from 1e-8 to 3.2,
and mole amounts 0.01, 1 and 3.7 include critical, dilute, dense, unstable and
extrapolated states. All 46 selected quantities match with relative 2e-10
and quantity-specific absolute floors in metadata. All fourteen constructor
metadata vectors match their qualified common fields and LJ-specific limits,
scaling and term counts. The remaining zero entries in NIST metadata vectors
do not qualify source bounds/term counts.

There are no failed reference processes. The fixtures retain 465 non-finite
source fugacities, 231 non-finite source sound values and 330 non-finite
returns for each C3/R134A supercritical ancillary branch as empty fields.
Native optional values/errors are checked explicitly instead of discarding
those states. Independent tests verify TV and extensive T/V/n derivatives,
Maxwell/caloric and sound identities, extensivity, scalar/bulk agreement,
reference scaling, reduced units, bounds and error behavior.

The source snapshot runs in process-main context, as the oracle does.
On this Windows host, CRT exp rounding differs by an ulp on a newly created
thread. The ill-conditioned C3 critical Cp demonstrates the amplification:
main -841538410.1318147 versus worker -841538399.2540895 J/(mol K), while
underlying kernels differ only at floating-point rounding scale. Physical
finite-difference tests run in normal test threads. Bitwise identity across
thread contexts, compilers or all extrapolations is not promised. Hosted CI
has not been observed here.

These data replace the previous broad documentation claims based on a small
1e-4 smoke test. The legacy tests remain as smoke tests. Volume solving has
only the earlier C3 closure checks; complete root selection, raw source
iteration trajectories, mixture/ideal-mixture behavior, alternate reference
states and multiparameter flashes/envelopes remain pending.

## Native state measurements

The [benchmark record](benchmark-meos.json) compares 24-output homogeneous
state loops for six profiles at T=1.4 Tc and 100 densities rho/rhoc from 0.05
to 1.535. Each sample evaluates 100000 states; five samples follow a warm-up.
All 24 fields match at 600 source and 1200 native preflight states and enter
timed checksums. The source calls its original jet, pressure, TV Helmholtz,
Cv/Cp and sound methods. Native state evaluation shares the jets. Source
NIST models retain their existing internal temperature caches; Rust fresh
mode builds preparation per point, while prepared mode reuses it. Model
construction, Python and FFI are outside both timers. The measured sound
field follows the direct source mass convention.

| Profile | Source ns/state | Fresh Rust ns/state | Prepared Rust ns/state | Source / fresh | Source / prepared |
|---|---:|---:|---:|---:|---:|
| C3 | 1731.14 | 3370.43 | 1162.23 | 0.51 | 1.49 |
| R134A | 6299.14 | 2264.37 | 1360.66 | 2.78 | 4.63 |
| N-H2 | 1527.93 | 2662.62 | 755.37 | 0.57 | 2.02 |
| LJ argon | 10612.67 | 3693.23 | 1384.09 | 2.87 | 7.67 |
| LJTS argon | 10013.24 | 3361.48 | 1348.94 | 2.98 | 7.42 |
| Reduced LJTS | 9821.10 | 3347.51 | 1315.71 | 2.93 | 7.46 |

Prepared native state queries are faster in all six measured profiles.
Fresh preparation remains slower than the internally cached source for C3
and N-H2; these results do not support a claim of universal speedup.
The record covers this Windows GNU host, these isotherms and the complete
24-output state workload. It excludes phase solving, the extensive derivative
APIs, singular responses, changing-temperature trajectories and the entire
library. Source/harness/Cargo fingerprints protect measurement inputs.

```
python tools/rust_port/generate_lj_meos_data.py
python tools/rust_port/validate_meos.py --compiler-directory /path/to/compiler/bin
cargo test --test meos --test meos_source --test multiparameter
python tools/rust_port/benchmark_meos.py --compiler /path/to/gfortran
```
