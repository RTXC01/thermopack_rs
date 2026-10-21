# The Wong-Sandler mixing rule

The Wong-Sandler rule replaces the binary attraction and the binary covolume of
a cubic equation of state with its own coefficients. It is the only mixing
rule of the original that does not keep the `(1 - kij)` structure of the
binary attraction, so it needs its own port instead of a parameter of the
classic one.

The original is `src/wong_sandler.f90`:

* `WongSandlerMix` evaluates the sums `Q` and `D` of the appendix A8 and A9 of
  the paper of Wong and Sandler, "A theoretically correct mixing rule for
  cubic equations of state", AIChE Journal 38 (1992) 671-680, and derives
  the mixture coefficients from them,
  $$b = \frac{Q}{n - D}, \qquad a = \frac{Q \, D \, R T}{n - D},$$
* `fidel_alpha` solves `alpha = 1 / (2 + g12 g21)` for the second order term,
* the binary parameters `mixWS%f_kij`, `mixWS%f_tauij` and `mixWS%alphaij`
  are fractions of the temperature.

The native port lives in `rust/src/wong_sandler.rs` with the types
`Fraction`, `WongSandlerData` and `WongSandlerCoefficients` and the function
`wong_sandler_mix`. The rule is reached through the cubic model:

```rust
let mut data = WongSandlerData::new(2);
data.set_kij(0, 1, Fraction::constant(0.05))?;
data.set_tauij(0, 1, Fraction::constant(0.001))?;
data.set_alphaij(0, 1, 0.3)?;
let options = CubicOptions {
    mixing_rule: MixingRule::WongSandler,
    wong_sandler: Some(&data),
    ..Default::default()
};
let model = Cubic::with_options(CubicKind::PengRobinson, &["CO2", "NC4"], options)?;
let coefficients = model.wong_sandler_coefficients(350.0, &[0.3, 0.7])?;
```

The model returns the mixture attraction and covolume with their temperature
derivatives and their mole number derivatives, which is the interface the
original `WongSandlerMix` offers through `cbeos%a`, `cbeos%b` and their
derivatives.

## Conventions of the original

The port keeps the conventions of the original, because they differ from the
conventions of the rest of the native library.

* The original keeps the attraction of the cubic equation of state in
  pascal litre squared per mole squared and the covolume in litres per mole,
  the native library uses pascal cubic metre per mole squared and cubic metre
  per mole. The values of the original therefore carry a factor of `1e-6` for
  the attraction and of `1e-3` for the covolume.
* The original uses `kRgas = 8314.4621` J per kilomole and kelvin, the native
  library uses `R = 8.3144621` J per mole and kelvin. The excess Gibbs energy
  of the original is therefore a thousand times larger, while the
  dimensionless combinations of the rule are unchanged.
* The limiting constant of the rule is the one of the infinite pressure limit,
  `excess_gibbs::getInfinitLimitC`: `1.0` for the van der Waals equation,
  `ln 2` for the Soave-Redlich-Kwong equation and `-ln(sqrt(2) - 1) / sqrt(2)`
  for the Peng-Robinson equation. The zero pressure limit
  `excess_gibbs::getZeroLimitC` belongs to the UNIFAC mixing rule and is not
  used here.
* The original writes the parameters of the second order term into the
  matrix element `(i, j)` and reads them transposed in `excess_gibbs::gEinf`.
  The native port keeps this convention, so the parameters of a pair have to
  be set in both directions.
* The original evaluates the alpha function of the components during a state
  evaluation. The mixing rule needs the alpha function, therefore the original
  calls `cbCalcAlphaTerm` before the rule. The native library evaluates the
  alpha correlation directly and does not need the extra call.

## Validation

`tools/rust_port/validate_wong_sandler.py` compares the fourteen mixture
coefficients of the native port against the exact source library. The
reference probe `tools/rust_port/reference_probe.f90` calls the original
`cbCalcAlphaTerm`, `excess_gibbs::gEinf` and `WongSandlerMix` and reports the
mixture coefficients of the original.

The eight cases cover the Peng-Robinson, the Soave-Redlich-Kwong and the van
der Waals equation with the carbon dioxide and the normal butane or the
propane, with and without binary interaction parameters, with and without the
parameters of the second order term. The agreement holds for the attraction
and the covolume, their temperature derivatives and their mole number
derivatives, with a worst relative difference of `1.26e-12` against a
tolerance of `1e-10`.

`tools/rust_port/validate_wong_sandler_state.py` compares the state of a
mixture with the rule against the exact source library. The reference probe
evaluates the pressure with its volume and temperature derivatives, the
enthalpy, the entropy, the fugacity coefficients and their temperature,
pressure and mole number derivatives at a molar volume.

The six cases cover the three equations of state with liquid-like and
gas-like molar volumes, with and without binary interaction parameters and
with and without the parameters of the second order term. All ninety-six
values agree, with a worst relative difference of `5.0e-13` against a
tolerance of `1e-9`.

The fixtures `tests/fixtures/wong_sandler.csv` and
`tests/fixtures/wong_sandler_state.csv` and the tests
`tests/wong_sandler.rs` and `tests/wong_sandler_state.rs` hold the qualified
values.

## Scope

The port delivers the coefficients of the rule and the state of a mixture with
the rule. The prepared state of the cubic equation of state takes the
attraction, the covolume and their derivatives from the rule, and the second
mole derivatives of the rule replace the binary coefficients of the classic
mixing rules in the second mole derivatives of the fugacity.

The original offers the rule also for the cubic part of the cubic
association model, `cbMixWSCPA`, and for the combination of the
Huron-Vidal and the Wong-Sandler rules, `cbMixHVWongSandler`. The three
variants of the original share the routine `WongSandlerMix` and differ in the
excess Gibbs model of the infinite pressure limit, which the rule consumes:
the parameters of the rule for `cbMixWongSandler` and `cbMixWSCPA`, the
Huron-Vidal model for `cbMixHVWongSandler`. The native library offers the
rule for the cubic equation of state with the parameters of the rule only.

