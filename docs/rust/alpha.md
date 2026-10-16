# Cubic alpha correlations

`Cubic::with_alpha(kind, components, requested, reference)` ports the alpha
initializer `cbAlpha::tpInitAlphaCorr`. It selects alpha parameters only;
binary interactions, ideal data, geometry and volume translation keep their
existing defaults. `Cubic::set_alpha(i, correlation, coefficients)` applies an
explicit correlation to one component (zero-based index). `alpha(i)` exposes
the requested/effective correlation, coefficients, fallback status and database
reference. Models own these settings independently.

`Alpha::evaluate(T,Tc)` returns alpha, d(alpha)/dT and d2(alpha)/dT2. These feed
the existing mixing, pressure, residual, fugacity, saturation and flash kernels.
Pure caloric flashes preserve a configured alpha when forming a pure-component
saturation model. Unsupported ideal data are retained as explicit errors and
do not block pressure/alpha calculations; caloric requests still return the
error rather than inventing a heat capacity.

## Source behavior

| Request | Initializer behavior | Explicit coefficients |
|---|---|---|
| CLASSIC | Maps to PR, SOAVE or VDW; computes the acentric parameter | One fitted classic coefficient |
| TWU | Reads the compiled TWU table; falls back to classic if absent/zero | Three Twu coefficients |
| MC | Also reads the TWU table in this fork, then uses the MC formula | Three MC coefficients; higher terms drop at T >= Tc |
| GERG | Falls back to classic | Source evaluation is unimplemented; Rust returns an error |
| CLASSICFIT, GD, SOAVE, PR, PT | Leaves the single coefficient zero | One classic coefficient |
| UMR | Water reads the MC table; other fluids store UMR parameters | UMR evaluation is unimplemented in the source, including explicit water UMR |
| RK, VDW | Selects the corresponding function | One stored coefficient, unused by evaluation |
| SW | Leaves its coefficient zero; uses temperature-dependent kappa | One kappa0 coefficient, with a branch at Tc |
| PR78 | Computes the acentric fit, with a branch at omega=0.491 | One classic coefficient |

`CubicKind::RedlichKwong` still denotes SRK geometry/interactions with RK alpha
by default. Selecting an alpha named PT or SW does not implement the PT or SW
EOS geometry. Those geometries remain pending.

`AlphaCorrelation::acentric_parameter(omega)` exposes the source acentric fit
separately. This distinction matters: an explicit initializer name such as PR
does not populate that fit, whereas CLASSIC for a PR model does. Raw source
selection behavior is preserved rather than replaced with textbook defaults.

The generator embeds all 192 TWU and 65 MC records from `compdatadb.f90`.
Selection preserves source order: the earliest exact non-DEFAULT reference
token wins, preferring the caller's earlier token; otherwise the last default
record wins. Case and slash-separated reference priorities are tested.

## Qualification scope

The Fortran validation probe links the existing exact-source library, sharing
its active model and globals. It invokes `tpInitAlphaCorr` directly and reads
`cbCalcAlphaTerm` results and effective parameters. Its source, compiler,
library and binary hashes are stored with fixtures.

* 4312 initializer/setter attempts: 4232 finite numerical source results match
  effective method, active coefficients, alpha and both derivatives. The 80
  recorded source terminations cover 56 non-water UMR initialization requests
  and 24 explicit UMR/GERG setters. Rust returns errors for all of them.
* 3240 EOS cases: all four existing geometries, 15 alpha configurations, three
  pure/mixture feeds, 250/400/700 K, 1 kPa/1 MPa/20 MPa and liquid/vapor roots.
  All match volume, pressure and derivatives, residual/total H/S, baseline Cp,
  log fugacity and its T/P/mole-number derivatives.
* Independent alpha/pressure/H finite differences, nonzero kij/lij mixing,
  configured pure latent-heat quality and mixture PH/PS round trips pass.

This alpha initializer is distinct from the public `eoslibinit::init_cubic`
constructor. That constructor forces TWU and volume translation when its
reference contains tcPR; QuantumCubic also changes critical data and beta.
Complete tcPR construction is qualified separately in the
[volume translation comparisons](volume-translation.md). QuantumCubic
construction/properties are qualified in [beta.md](beta.md); its expanded
cryogenic flash grid and language-interface compatibility remain pending.
The alpha probe deliberately
isolates the initializer; its results alone do not qualify translated/quantum
models. All original default-model fixtures continue to run alongside these
additional cases.
