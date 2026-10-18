# PeTS Helmholtz and TV implementation

`Pets::new()` constructs this fork's pure-argon PeTS model: segment diameter
3.42e-10 m and epsilon/k=136 K, with the argon ideal-gas reference state.
`with_parameters` and `set_parameters` accept positive finite diameter (m) and
well depth divided by Boltzmann's constant (K). Each model owns its parameters
and ideal reference state. Changing one model leaves clones independent.

```rust
use thermopack_rs::Pets;
let model = Pets::new()?;
let prepared = model.prepare(300.0)?;
let f = prepared.residual_helmholtz(1e-3, 1.0)?;
let state = prepared.caloric_tv(1e-3, 1.0)?;
assert!(state.pressure.value > 0.0);
# Ok::<(), thermopack_rs::Error>(())
```

TV calls accept total volume V (m3) and mole number n (mol), rather than a
normalized composition or molar volume. `PureHelmholtz` returns the reduced
extensive residual F=Ares/(R*T), in mol, and all nine first/second T,V,n
derivatives. `DensityHelmholtz` describes dimensionless molar alpha and its
T/rho derivatives. The shared conversion to extensive derivatives follows
the chain rule with rho=n/V; independent Euler and mixed-partial tests cover
every derivative.

`caloric_tv` returns pressure with T/V/n derivatives, total U/H/A (J), S (J/K),
Cv/Cp (J/K), chemical potential (J/mol) and the pure log fugacity coefficient.
Fugacity is `None` at nonpositive pressure; Cp is `None` at singular/non-finite
isobaric response. These calls evaluate homogeneous states, including unstable
loops, without asserting phase stability. The domain requires T,V,n>0 and
packing fraction below one. Non-finite values produce errors. `ideal_gas_mut`
exposes the existing owned Cp/reference setters.

The model follows `pets.f90`:

```
d(T) = sigma * [1 - 0.127112544 * exp(-3.052785558 * epsilon_k / T)]
eta = (pi/6) * d^3 * N_A * rho
alpha_hs = (4*eta - 3*eta^2)/(1-eta)^2
alpha_disp = -pi*N_A*sigma^3*rho * [2*I1(eta)*q + C1(eta)*I2(eta)*q^2]
q = epsilon_k/T
C1 = 1 / [1 + (8*eta-2*eta^2)/(1-eta)^4]
```

The seventh-coefficient I1/I2 polynomials and all analytic derivatives come
from the fork. One evaluation shares diameter, packing fraction and polynomial
terms across hard-sphere and dispersion contributions. Simultaneous Horner
evaluation produces each polynomial and its first two derivatives.
`prepare(T)` caches temperature terms for repeated density sweeps. The kernel
and prepared calls allocate no scratch memory and make no Fortran calls.

`virial_coefficients(T)` evaluates analytic zero-density limits B and C, with
P/(RT)=rho+B*rho^2+C*rho^3+.... The native zero-density density derivatives are
finite analytic limits; the source kernel's expressions are singular at
exactly rho=0. `baseline_virial_coefficients` reproduces the fork's n=1 mol,
V=1e7 m3 estimator and preserves its density-polynomial and derivative
operation order. Its C subtracts large nearly equal terms: at T/epsilon_k=0.02,
ordinary algebraic reassociation changed the returned value by rounding bins
of order 1e-10–1e-9 m6/mol2. The compatibility path retains the source order;
the main density kernel uses the faster shared evaluation. Tests independently
bound finite-volume cancellation and check the analytic limit.

`speed_of_sound_tv` evaluates the physical isentropic pressure response and
returns an error for invalid stiffness. `baseline_speed_of_sound_tv` preserves
the source zero return for nonnegative/NaN stiffness. Constant-entropy pressure
perturbations independently verify physical sound; tests also verify extensive
scaling, caloric derivatives, constant-pressure Cp, Maxwell relations, chemical
potential and ideal-reference shifts.

`potential(r)` returns the truncated/shifted pair potential divided by k (K)
and force divided by k (K/m), with cutoff 2.5 sigma. `potential_into` writes
caller-owned buffers and reproduces the source's first-cutoff termination:
after encountering an out-of-cutoff distance, later values remain zero even
for unsorted input. `lj_truncation_correction(T,V,n)` returns the reduced
extensive full-LJ minus truncated/shifted correction and its ten quantities
using the source's two correction contributions.

The Python PeTS wrapper mistakenly resolves its setter to the Fortran getter:
its cached Python parameters can change without changing the active EOS.
Rust setters and reference probes use the actual `pets_set_pure_params`
routine's semantics. Python binding compatibility, including that wrapper
behavior, remains separate from the native numerical kernel.

Qualification compares 588 full state cases, each with 56 outputs, and 16
vector-potential cases against the exact-source build. Four parameter pairs,
seven reduced temperatures 0.02–10, seven packing fractions 1e-10–0.85 and
n=0.3/1/3.5 cover default/custom models, dilute/dense states and extensive
scaling. All 604 attempted source calls are finite and match. There are 120
nonpositive-pressure states and 72 zero-sound outputs. Quantity-specific
absolute tolerances and 5e-10 relative tolerance are stored in metadata;
virial C uses 1e-11 m6/mol2 absolute tolerance without widening it for the
cancellation-sensitive cases. The original source remains unmodified.

This qualifies the PeTS Helmholtz/TV, potential, correction and virial paths.
Pressure roots and pure critical/saturation calculations now have separate
[equilibrium qualification](pets-equilibrium.md), including source last-iterate
compatibility and strict pressure closure. [Pure TP/PH/PS flashes](pets-flash.md)
also have separate qualification, as do [UV/SV flashes](pets-isochoric.md).
Other inherited SAFT
utilities, higher derivatives and language bindings remain pending. This does
not establish a complete PeTS interface replacement or a complete ThermoPack
replacement.

Reproduce after building the exact-source library and probe:

```
python tools/rust_port/validate_pets.py --compiler-directory /path/to/compiler/bin
cargo test --test pets
python tools/rust_port/benchmark_pets.py --compiler /path/to/gfortran
```

The native-loop [benchmark report](benchmark-pets.json) records the full
kernel implementation committed as `f653b85`. It measures the
ten-quantity Helmholtz kernel at 300 K, one mole and 1000 total volumes from
5e-5 to 1.049e-3 m3. All 1000 state outputs are compared before timing, then
one million evaluations per sample compare all ten accumulated checksums.
Five samples follow a discarded warm-up. Both loops are native; construction,
data generation, comparisons and I/O are outside the timers. Reports retain
source/Cargo hashes, compiler flags and raw samples. Source changes during a
measurement cause the driver to reject the run. Timings qualify this kernel
workload; higher-level PeTS performance remains unmeasured.

| Full Helmholtz operation | ns/state | Speedup over source |
|---|---:|---:|
| Source Fortran | 4481 | 1.00x |
| Rust including temperature preparation | 149 | 30.07x |
| Rust using prepared temperature terms | 89.8 | 49.90x |

The local comparison uses the recorded gfortran `-O3` reference build and
Rust release/thin-LTO profile. Shared analytic terms and fixed-size scalar
derivatives avoid the source's repeated packing evaluations and dynamic
component-sized scratch arrays. The measured operation returns F plus every
first/second T,V,n derivative; it is not a scalar-energy-only comparison.
