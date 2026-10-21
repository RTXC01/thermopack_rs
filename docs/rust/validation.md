# Reproducing validation

The native crate is dependency-free. Rust 1.88 or newer is required.

```
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test --all-targets
cargo test --release --tests
cargo test --doc
cargo doc --no-deps
python tools/rust_port/verify_fixtures.py
```

`cargo test` checks committed reference fixtures without a Fortran runtime.
`tests/fixtures/*.metadata.json` records source hashes, compiler, flags,
the pinned LAPACK commit, the reference library hash and tolerances.
Fixture SHA-256 values cover exact CSV bytes with LF line endings, enforced
by `.gitattributes` and the generators on every platform. The integrity checker
also verifies row counts; native CI runs it without building Fortran.

To regenerate fixtures from this fork, install gfortran (the local qualification
uses 14.3.0; the source warns about version 15), Python and NumPy:

```
git submodule update --init --depth 1 external/lapack
python tools/rust_port/build_reference.py --compiler /path/to/gfortran
python tools/rust_port/build_probe.py --compiler /path/to/gfortran
python tools/rust_port/validate_cubic.py
python tools/rust_port/validate_saturation.py
python tools/rust_port/validate_derivatives.py
python tools/rust_port/validate_flash.py
python tools/rust_port/validate_initialization.py
python tools/rust_port/validate_caloric.py
python tools/rust_port/validate_isochoric.py
python tools/rust_port/validate_pets.py
python tools/rust_port/validate_pets_equilibrium.py
python tools/rust_port/validate_pets_flash.py
python tools/rust_port/validate_pets_isochoric.py
python tools/rust_port/validate_mbwr.py
python tools/rust_port/validate_mbwr_properties.py
python tools/rust_port/benchmark_isochoric.py --compiler /path/to/gfortran --preflight-only
python tools/rust_port/benchmark_isochoric.py --compiler /path/to/gfortran --tcpr --preflight-only
python tools/rust_port/validate_ideal.py
python tools/rust_port/validate_ideal_configuration.py
python tools/rust_port/validate_alpha.py
python tools/rust_port/validate_alpha_properties.py
python tools/rust_port/validate_translation.py
python tools/rust_port/validate_beta.py
python tools/rust_port/validate_rr.py
python tools/rust_port/validate_rr_iteration.py
python tools/rust_port/validate_symmetric.py
python tools/rust_port/validate_response.py
python tools/rust_port/validate_sound.py
python tools/rust_port/validate_translated_flash.py
python tools/rust_port/generate_unifac_data.py
python tools/rust_port/generate_solid_data.py
python tools/rust_port/validate_unifac.py
python tools/rust_port/validate_unifac_cubic.py
python tools/rust_port/validate_solid_lines.py
python tools/rust_port/validate_solid_co2.py
python tools/rust_port/validate_solid_saturation.py
python tools/rust_port/validate_pseudo_components.py
python tools/rust_port/validate_isolines.py
python tools/rust_port/validate_spinodal.py
python tools/rust_port/validate_joule_thompson.py
python tools/rust_port/validate_virial.py
python tools/rust_port/validate_virial_matrix.py
python tools/rust_port/validate_ljs.py
python tools/rust_port/validate_ljs_spline.py
python tools/rust_port/validate_ljs_helmholtz.py
python tools/rust_port/validate_ljs_state.py
python tools/rust_port/validate_ljs_fugacity.py
python tools/rust_port/validate_ljs_caloric.py
python tools/rust_port/validate_ljs_potential.py
python tools/rust_port/validate_ljs_wca_term.py
python tools/rust_port/validate_wca_term.py
python tools/rust_port/validate_white_bear.py
python tools/rust_port/validate_uv_phi.py
python tools/rust_port/validate_uv_terms.py
python tools/rust_port/validate_quadrature.py
python tools/rust_port/validate_critical.py
python tools/rust_port/validate_mut_solve.py
python tools/rust_port/validate_hyperdual.py
```

On Windows add `--compiler-directory /path/to/gfortran/bin` to the validation
commands so Python can resolve runtime DLLs. The local isolated compiler is
under `.toolchains/gcc14/mingw64/bin` from the
[WinLibs standalone distribution](https://winlibs.com/). No global installation
or PATH modification is required; set PATH for the current shell to find cargo.
Build artifacts and toolchain archives are ignored by Git. The reference loader
injects platform settings in its Python process; it does not edit the wrapper.

The builder uses the upstream Makefile dependency order and compiles only the
pinned LAPACK/BLAS routines required by ThermoPack. A fingerprint of all source
hashes, compiler, flags and LAPACK revision invalidates stale cached objects.

## Qualified cases

* 2088 fluid acoustic attempts: 817 equilibrium and 1258 TV results match;
  13 source TP terminations are recorded. Pure qualities include zero/one
  endpoints, and root overrides exercise supplied liquid/liquid and vapor/vapor
  evaluations. TV cases include 29 nonpositive-pressure states, 10 with positive
  source sound speed, and total-mole scaling. Independent constant-entropy
  finite differences check physical mixture/pure acoustic responses. See
  [sound conventions and measurements](sound.md).
* 441 equilibrium TP-response attempts: 435 finite results match phase, reported
  Cp, dV/dT and dV/dP; six source terminations are recorded. Independent native
  finite differences in 72 classical/translated states check physical Cp,
  volume T/P derivatives, internal-energy T/P derivatives and Maxwell symmetry.
  See [equilibrium response](equilibrium-response.md).
* 226 upper symmetric indefinite systems match pinned LAPACK factorization,
  signed pivots, solutions and singularity status, with no source failures.
  Native tests independently verify componentwise backward error and untouched
  lower-triangle entries. See [symmetric solver](symmetric-solver.md).
* 704 extended/one-step Rachford–Rice attempts: 698 finite source results/statuses
  match; six non-finite source updates are recorded. See
  [extended RR conventions](rachford-rice.md). This does not qualify full flash
  convergence or source output fields left undefined on a failed RR solve.
* 2170 LV RR substitution attempts: 2110 finite source results/statuses match
  both flags, fractions, phase compositions, K updates, log fugacity and Gibbs
  energy. Sixty non-finite source one-step calls are recorded. All four
  convergence/RR-status combinations are covered. Full TP restart/Newton and
  optional non-vapor phaseY paths remain pending.
* `baseline_tp_flash` matches raw ordering, both fractions and composition
  arrays for all 2173 classical/translated TP snapshots and 718/721 quantum
  snapshots. `validate_baseline_tp.py` rechecks every frozen result and currently
  returns failure for three vanishing-vapor composition differences. See
  [fork-order TP qualification](phase-equilibrium.md#fork-order-tp-compatibility).
* 1320 cubic states: PR, SRK, SRK with RK alpha, vdW; six pure fluids and five
  mixtures, 200/300/500 K, 1 kPa through 20 MPa, liquid/vapor requests, zero
  mole-fraction components and nonzero binary interactions.
* Quantities: volume, pressure, analytic dP/dT and dP/dv, residual/total H/S,
  log fugacity coefficients and the baseline enthalpy temperature derivative.
* 48 pure saturation states: methane, CO2 and propane; four cubic models;
  reduced temperatures 0.65, 0.8, 0.95 and 0.999; coexistence pressure, phase
  volumes and Clapeyron slope.
* Independent finite-difference and physical-identity tests, low-pressure
  cubic roots, input errors, Rachford-Rice closure and concurrent model use.
* 1320 fugacity T/P/mole-number derivative cases, independently checked with
  finite differences, Hessian symmetry, Gibbs-Duhem and nonzero lij.
* 1046 attempted TP flashes: 1043 numerical source results, including the
  upstream eight-component near-critical regressions. One physical split
  requires exchanging phase labels; three source calls terminate their process.
  See [phase equilibrium qualification](phase-equilibrium.md).
* 64 initialization cases: 28 fixed-composition pseudo-critical states and
  36 pure saturation-temperature inversions across all four cubic models.
* 912 attempted PH/PS flashes: 697 closing source results match temperature,
  phase code, both fractions and both composition arrays. The other 215
  source calls abort, raise calculation errors or return states that fail
  independent caloric closure. They remain recorded failures, not parity
  successes. The grid covers eight feeds, four temperatures, three pressures
  and pure latent-heat qualities. Metadata contains every failed request.
* Equilibrium Cp includes phase redistribution and agrees with independent
  finite differences of equilibrium H and S. Native PH/PS results also close
  their specified H or S independently of the reference comparisons.
* 2616 attempted ideal-property cases covering every default fluid, eight
  temperatures (including the 30 K hydrogen branch), and three pressures:
  2592 finite source results across 108 fluids match H, S, Cp and the reported
  entropy derivatives. BUTANAL's 24 attempts terminate the reference with its
  invalid default Cp correlation 100. Native H/S temperature derivatives are
  independently checked for all 108 supported fluids, including the explicit
  type-10 hydrogen factor-1000 discrepancy.
* 720 configured ideal-property cases: all 12 correlations, all four hydrogen
  isomers for type 10, six ordered setter sequences, four temperatures and two
  pressures. They compare properties, reference getters, correlation IDs and
  all 21 coefficients; correlation IDs/coefficients require exact equality.
  Native tests check parameter isolation, invalid setter atomicity and PH/PS
  closure with customized correlations and references.
* 4312 alpha initializer/setter attempts: 4232 numerical results and 80 recorded
  source terminations in unimplemented paths. Another 3240 EOS cases qualify
  alpha variants through complete properties and fugacity derivatives. See
  [alpha scope and source conventions](alpha.md); this isolates the alpha
  initializer and does not by itself qualify translated/quantum constructors.
* 1977 volume-translation parameter cases and 3120 full EOS/property cases
  qualify polynomial ci correlations, the Rackett fallback, full tcPR
  construction, reference selection, setters and fugacity derivatives.
  See [volume translation](volume-translation.md) for API and numerical scope.
* 3960 translated flash/saturation attempts: 1130 TP, 1293 PH, 1326 PS and
  72 pure saturation results match; 86 source terminations and 53 nonclosing
  PH results are recorded separately. Native H/S closure and two-phase
  material balances are independently checked in the committed fixture tests.
* 4036 beta/critical-data parameter and 4512 complete property cases qualify
  explicit beta settings and full QuantumCubic construction, including analytic
  derivatives, owned critical/mass overrides and nonzero kij/lij. All attempted
  source calls are numerical and match. See [beta.md](beta.md). The separate
  `validate_quantum_flash.py` diagnostic currently reports unresolved numerical
  and branch differences and is not a passing qualification gate.

Each comparison uses `abs(error) <= atol + rtol * max(abs(reference), abs(rust))`.
The [MBWR qualification](mbwr.md) adds exact raw/integrated coefficients for
65 stored fits and 6760 states with 26 pressure/Helmholtz outputs. All source
calls are finite. Source-order methods match the full matrix; the separate
stable integral has 318 recorded summation differences at extreme dilute,
cryogenic extrapolations. Independent tests verify finite differences,
Helmholtz pressure identities, Euler/scaling and analytic zero-density limits.
The [MBWR residual-state qualification](mbwr-properties.md) adds 5983 finite
state-function results with 18 outputs, retaining 1817 nonpositive-pressure
source results. Original entropy response signs and raw mole-response residues
are compared explicitly; the source-snapshot executable uses process-main
context to match the oracle's Windows CRT math behavior. Physical finite-
difference checks run on normal test workers.
The [MBWR equilibrium qualification](mbwr-equilibrium.md) adds 65 metadata rows,
49 independently solved original-pressure critical points, and 3406 closed/stable
source density iterates that belong to native bounded root sets. Fifty unqualified
source returns and 60 requested-phase selection differences are retained. This
qualifies root membership, not full source phase-selection or trajectory fidelity.
The [MBWR source density qualification](mbwr-density.md) separately compares
8839 finite raw iterates and exact phase flags across all 65 fits, retaining 438
unclosed source returns in the comparison and recording one source NaN/native
error. Explicit initialization, including synthetic metadata for 17 legacy fits,
and original default-initializer limitations are documented.
The [MBWR default-initialization qualification](mbwr-initialization.md) adds
32 finite initializations and 27 explicit source critical-grid results, checking
native errors for 16 original terminations and 17 undefined missing-component
lookups. All 65 default initialization statuses are retained in a separate
fixture; the source bounds defect is explicitly distinguished from valid output.
The [PeTS flash qualification](pets-flash.md) adds 923 finite/closing TP/PH/PS
states from 1035 source attempts, retaining 112 failed/nonclosing calls.
Saturated endpoints use an explicit physical-state comparison convention;
metadata preserves the original phase flags and arrays for every finite result.
The [PeTS UV/SV qualification](pets-isochoric.md) adds 399 finite/closing states
from 666 source calls; 260 source terminations and seven unqualified returns
are retained. Native tests also verify analytic equilibrium Cv against fixed-V
differences of both U and S, homogeneous/latent recovery and caloric setters.
The [PeTS qualification](pets.md) adds 588 full states (56 quantities each)
and 16 vector-potential cases. All attempted source calls are finite. Independent
tests verify extensive and mixed derivatives, caloric/mechanical identities,
constant-entropy sound, constant-pressure Cp and analytic virial limits.
The [UV/SV and direct TV qualification](isochoric.md) adds 409 finite source
flash results and 141 TV property results. Metadata retains 137 failed source
calls. Native tests independently check U/S/V closure, component balance,
chemical equilibrium and caloric derivatives; single-phase comparison
canonicalization and remaining raw output limits are explicit.
The [UNIFAC qualification](unifac.md) adds 1540 unit-level excess Gibbs
comparisons against the exact-source library: 20 pure components, all 190
binary pairs and five ternary mixtures at three temperatures, with the
`UNIFAC`, `UMR` and `VTPR` labels. Every case matches within a relative
tolerance of 2e-10, covering the excess Gibbs energy, both temperature
derivatives, the composition derivatives, the mixed derivative and the second
composition derivatives. The source reports energies per kilomole, so the
native per-mole values are scaled by 1000 for the comparison. A second
qualification compares the mixing rule inside the cubic equation of state:
270 of 288 PR and SRK states match with worst scaled errors of 1.5e-11 for the
pressure and 2.4e-12 for the fugacity coefficients; the 18 remaining cases
terminate in the source because its UNIFAC component database spells toluene
`TOLUENE` while the equation of state uses `TOLU`, and the native model
reports the same missing component.
The [solid line qualification](solid.md) adds 156 melting and sublimation
line comparisons against the exact-source library: 23 melting and 16
sublimation components at four temperature limits with seven points per line.
152 cases match, including the temperature grid and every pressure, with a
relative tolerance of 1e-10 and an absolute tolerance of 1e-4 Pa that cover
the cancellation of nearly equal terms in some melting curves; the 4
remaining cases are the hydrogen sulfide sublimation lines, where the source
stops because the correlation label has no fourth character and the native
model reports the same error.
The [mu-T qualification](mut-flash.md) adds 32 log fugacity solves over four
mixtures and four temperatures, with the specification taken at the true state
and the solve started from a different guess. All 32 states agree with the
exact-source library to a worst relative difference of 7.4e-13. The original's
chemical potential variant derives its specification from the initial guess
and is therefore a fixed point, which the native model reproduces and the
test verifies.
The [solid carbon dioxide qualification](solid.md#solid-carbon-dioxide) adds
7 states with 13 quantities each: the Gibbs energy of solid carbon dioxide
with its first and second temperature and pressure derivatives, the mixed
derivative, the specific volume, the entropy, the enthalpy, the internal
energy, the Helmholtz energy, the heat capacity and the speed of sound. All
7 states match the exact-source library within a relative 1e-11 after the
source reference shift.
The [solid saturation qualification](solid.md#solid-saturation-line) adds 6
single component cases of the original `solidPointOnEnvelope`: three
temperature specifications and three pressure specifications on the
solid-fluid line of carbon dioxide. All 6 cases match, the solved pressure
within a relative 1e-6 and the solved temperature within an absolute 1e-8 K
and a relative 1e-9, with an observed worst relative difference of 3.3e-8.
The reference levels of the solid model follow the liquid at the triple point
of the original, 216.592 K and 513602.16 Pa.
The [pseudo component qualification](pseudo-components.md) adds 30 cases of the
original `init_cubic_pseudo`: a single pseudo component, two pseudo components
and a pseudo component with a real fluid, in the Peng-Robinson and the
Soave-Redlich-Kwong equations of state. The molar volume, the enthalpy, the
entropy and the compressibility factor match with a relative 1e-10, and the
observed worst relative difference is 7.9e-15.
The [iso-line qualification](isolines.md) adds 13 cases of the original
`isotherm`, `isobar`, `isenthalp` and `isentrope`: isotherms and isobars of
two pure fluids and of two mixtures in the Peng-Robinson and the
Soave-Redlich-Kwong equations of state, and the isenthalp and the isentrope
of methane and of the methane-propane mixture. The temperature, the pressure,
the molar volume, the molar entropy and the molar enthalpy of every point
match with a relative 1e-6, and the observed worst relative difference is
1.4e-9. The not-a-number point of the original, which a grid of one point
produces, is reproduced.
The [stability limit qualification](spinodal.md) adds 12 cases of the original
`initial_stab_limit_point`: the liquid and the vapor branch of the meta
stability limit of methane and of carbon dioxide in the Peng-Robinson
equation of state and the two branches of methane in the Soave-Redlich-Kwong
equation of state. The temperature and the molar volume of the point match
with an absolute 1e-6 K and an absolute 1e-10 cubic metre per mole, and the
observed worst relative difference is 1.4e-9, which is the tolerance of the
source solver.
The [Joule-Thomson qualification](spinodal.md#the-joule-thomson-inversion)
adds 78 cases of the original `map_jt_inversion`: the inversion curve of
nitrogen in the Peng-Robinson and the Soave-Redlich-Kwong equations of state.
The temperature and the pressure of the point match with an absolute 1e-6 K
and an absolute 1e-4 Pa, and the observed worst relative difference is 3.3e-9.
The [virial coefficient qualification](equilibrium-response.md#virial-coefficients)
adds 10 cases of the original `virial_coefficients`: methane, carbon dioxide
and nitrogen, two mixtures, in the Peng-Robinson and the Soave-Redlich-Kwong
equations of state. The second virial coefficient matches with a relative 1e-4
and an absolute 1e-8 cubic metre per mole, and the observed worst relative
difference is 1.8e-4, the truncation error of the finite volume estimate that
the native model and the original share. The third virial coefficient of the
original carries the cancellation of that estimate and is not reproduced.
The [virial matrix qualification](equilibrium-response.md#the-matrix-of-the-second-virial-coefficients)
adds 6 cases of the original `secondVirialCoeffMatrix`: methane, the
methane-propane mixture and the carbon dioxide-methane-nitrogen mixture in
the Peng-Robinson equation of state and the methane-propane mixture in the
Soave-Redlich-Kwong equation of state. The matrix of the pure and the cross
second virial coefficients matches with a relative 1e-3, and the observed
worst relative difference is 1.5e-4.
The [Lennard-Jones splined qualification](ljs.md) adds 8 cases of the original
`calc_ai_reduced_LJs_ex` and `ljs_bh_get_bh_diameter_div_sigma`: the three
reduced perturbation terms and the reduced Barker-Henderson hard-sphere
diameter of the splined Lennard-Jones model. The four quantities match with a
relative 1e-10, and the observed worst relative difference is 1.5e-12.
The [Lennard-Jones spline qualification](ljs.md) adds 8 cases of the original
`calc_ai_LJs_ex`: the four polynomials of the spline of the excess energy and
their partial derivatives. The polynomial, the first and the second
derivative with respect to the reduced temperature, the first and the second
derivative with respect to the reduced density and the mixed derivative match
with a relative 1e-10, and the observed worst relative difference is 3.1e-14.
The [Lennard-Jones Helmholtz energy qualification](ljs.md) adds 12 cases of the
original `calcFresLJs_bh`: the reduced residual Helmholtz energy of the
splined Lennard-Jones model over the dilute and the dense fluid, and the same
states with the Lafitte form of the third order term. The energy and its
temperature derivative match with a relative 1e-9, and the observed worst
relative difference is 4.4e-14.
The [Lennard-Jones state qualification](ljs.md) adds 8 cases of the original
state of the splined Lennard-Jones model: the pressure of the state and the
compressibility factor, with the analytic volume derivative of the reduced
residual Helmholtz energy. The pressure matches with a relative 1e-6 and the
compressibility factor with a relative 1e-6, and the observed worst relative
difference is 1.9e-13.
The [Lennard-Jones chemical potential qualification](ljs.md) adds 8 cases of
the original `calcFresLJs_bh`: the reduced residual chemical potential of the
single component, which the original reports as the derivative of the reduced
energy with respect to the number of moles. The quantity matches with a
relative 1e-6, and the observed worst relative difference is 1.1e-15.
The [Lennard-Jones caloric qualification](ljs.md) adds 8 cases of the original
state of the splined Lennard-Jones model: the molar enthalpy and the molar
entropy of the state. The two quantities match with a relative 1e-8, and the
observed worst relative difference is 4.2e-15.
The [Lennard-Jones potential qualification](ljs.md) adds 4 cases of the
original `ljs_potential_reduced`: the reduced interaction potential and the
reduced force over the attractive branch, the minimum and the smoothing of the
cutoff. The two quantities match with a relative 1e-12, and the observed worst
relative difference is 2.8e-15.
The [White Bear term qualification](ljs.md) adds 8 cases of the original
`calc_LJs_WCA_ai_tr`: the four perturbation terms of the White Bear family
over the reduced temperature and the reduced density. The terms match with a
relative 1e-12, and the observed worst relative difference is zero.
The [White Bear cavity term qualification](ljs.md) adds 5 cases of the original
`calc_dhs_WCA`, `minimum_r_cavity_int_LJ`, `cavity_distribution_function` and
`lj_boltzmann_factor`: the hard-sphere diameter of the reference, the radius of
the minimum of the cavity integral, the cavity correlation function of the hard
spheres and the Boltzmann factor of the Lennard-Jones potential, on both sides
of the packing fraction where the source changes the correlation. The fifteen
values match with a relative 1e-11, and the observed worst relative difference
is 1.6e-14.
The [UV theory fraction qualification](ljs.md) adds 5 cases of the original
`uv_phi`: the u-fraction of the UV theory and its five partial derivatives.
The quantities match with a relative 1e-12, and the observed worst relative
difference is 4.1e-15.
The [UV theory term qualification](ljs.md) adds 5 cases of the original
`uv_a1_u` and `uv_a1_b2`: the first and the second order terms of the UV
theory and their partial derivatives. The nine quantities match with a relative
1e-11, and the observed worst relative difference is 2.5e-16.
The [UV theory virial qualification](ljs.md) adds 4 cases of the original
`uv_q_of_T`, `uv_hs_b2` and `ljs_b2`: the q-function of the UV theory, the
hard-sphere second virial coefficient and the second virial coefficient of the
Lennard-Jones spline with their temperature derivatives. The nine quantities
match with a relative 1e-11, and the observed worst relative difference is
3.7e-16.
The [White Bear term qualification](ljs.md) adds 5 cases of the original
`calc_LJs_WCA_ai`: the same perturbation terms with the prefactor of the
source and their five partial derivatives. The quantities match with a relative
1e-12, and the observed worst relative difference is 1.0e-15.
The [quadrature qualification](quadrature.md) adds 14 cases of the original
`get_quadrature_positions`, `get_quadrature_weights` and
`integrate_with_differential`: the seven Gauss rules and the seven Gauss-Kronrod
rules. The nodes, the weights, the integral and the differential match with a
relative 1e-13, and the observed worst relative difference is 3.3e-16.
The [critical point qualification](spinodal.md#critical-points) adds 7 cases of
the original `calcCriticalTV`: methane, carbon dioxide, propane and water in
the Peng-Robinson equation of state, nitrogen and methane in the
Soave-Redlich-Kwong equation of state and a mixture with a zero mole fraction
of the second component. The temperature, the molar volume and the pressure of
the critical point match with a relative 1e-6, and the observed worst relative
difference is 3.8e-7, which is the tolerance of the source solver.
The [multiparameter density qualification](multiparameter.md#density-iteration)
adds 1694 attempted density solves over the fourteen model configurations: 1669
states match the exact-source library, including the density, the raw phase
flag, the source `ierr` and the final pressure and its density derivative. The
density tolerance adds the conditioning term `1e-14*|p|/|dP/drho|`, because
the Newton step divides the kernel's 1e-14 relative pressure difference by a
small dP/drho. The 25 remaining states terminate in the source because its
ancillary rejects the phase flag, and the native model reports the same error.

Tolerances are quantity-specific and recorded in fixture metadata and tests.
They are not a promise of bitwise identity across compilers or all states.
The vdW baseline derivative defect is preserved by an explicit compatibility
method; see [baseline-defects.md](baseline-defects.md).

The [ice Ih qualification](solid-water.md) adds 780 finite source states with
twelve outputs, six reference scenarios, all first/second Gibbs derivatives,
and an explicit baseline for the source Cp mass-conversion defect. Physical Cp
and derivative identities have independent finite-difference checks. Its
native-loop benchmark preflights 100 source and 200 fresh/prepared native states.

The [FMT qualification](fmt.md) adds 1288 finite grid returns across three
functionals and all seven combinations of optional outputs, with up to 43
quantities per return. All 56 non-finite source returns are retained as errors.
The analytic path has independent gradient/Hessian and bulk-limit tests;
the source baseline retains operand order and reciprocal derivative floors.

The [pure multiparameter qualification](multiparameter.md) adds 5082 states
with 46 quantities and fourteen metadata vectors for five NIST fluids, eight
LJ noble-gas fits and reduced LJTS. It includes extensive T/V/n derivatives,
critical states, three mole amounts, non-finite source outputs and an explicit
baseline for the source's LJ sound mass units. Source snapshots run in the
oracle's process-main context; physical identity tests also run in normal
test threads. This supersedes the earlier loose NIST smoke-test evidence.

CI contains native Windows/Linux/macOS checks and a Linux exact-source build
that regenerates comparisons. It has been configured locally; hosted CI has
not run yet. Original pFUnit/Python example suites and models outside the
coverage ledger's qualified subset remain pending.

The [hyperdual qualification](hyperdual.md) adds 6 states covering the first,
second, third and zero differentiation order, with the twelve elementary
functions of the calculus and their eight components each. All 6 states agree
with the exact-source library and the largest observed relative difference is
zero. The third order derivative of the arc tangent reproduces the double
power that the original writes.
