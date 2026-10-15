# Reproducing validation

The native crate is dependency-free. Rust 1.88 or newer is required.

```
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test --all-targets
cargo test --release --tests
cargo test --doc
cargo doc --no-deps
```

`cargo test` checks committed reference fixtures without a Fortran runtime.
`tests/fixtures/*.metadata.json` records source hashes, compiler, flags,
the pinned LAPACK commit, the reference library hash and tolerances.

To regenerate fixtures from this fork, install gfortran (the local qualification
uses 14.3.0; the source warns about version 15), Python and NumPy:

```
git submodule update --init --depth 1 external/lapack
python tools/rust_port/build_reference.py --compiler /path/to/gfortran
python tools/rust_port/validate_cubic.py
python tools/rust_port/validate_saturation.py
python tools/rust_port/validate_derivatives.py
python tools/rust_port/validate_flash.py
python tools/rust_port/validate_initialization.py
python tools/rust_port/validate_caloric.py
python tools/rust_port/validate_ideal.py
python tools/rust_port/validate_ideal_configuration.py
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

Each comparison uses `abs(error) <= atol + rtol * max(abs(reference), abs(rust))`.
Tolerances are quantity-specific and recorded in fixture metadata and tests.
They are not a promise of bitwise identity across compilers or all states.
The vdW baseline derivative defect is preserved by an explicit compatibility
method; see [baseline-defects.md](baseline-defects.md).

CI contains native Windows/Linux/macOS checks and a Linux exact-source build
that regenerates comparisons. It has been configured locally; hosted CI has
not run yet. Original pFUnit/Python example suites and models outside the
coverage ledger's qualified subset remain pending.
