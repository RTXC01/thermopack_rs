# Covolume beta factors and QuantumCubic

`Beta::new(BetaCorrelation::Quantum, [xi, Tadj])` evaluates the source's
temperature-dependent covolume factor, including its first two derivatives:

```
beta(T) = [(1 + xi/(T + Tadj)) / (1 + xi/(Tc + Tadj))]^3
bi(T)   = bi(Tc) * beta(T)
```

Classic beta is one, with zero derivatives. `Cubic::set_beta(i, beta)` changes
one component, using a zero-based index; `beta(i)` returns its active settings.
Parameters are owned, finite values. A finite but singular correlation returns
an evaluation error. The prepared state caches beta-dependent covolumes and
their derivatives only when a non-classic correlation is active.

The initializer's quantum parameters are:

| Component | xi | Tadj (K) |
|---|---:|---:|
| HE | 1.4912 | 3.2634 |
| NE | 0.4673 | 2.4634 |
| H2 | 3.0696 | 12.682 |
| D2 | 1.6501 | 7.309 |
| P-H2 | 3.0710 | 12.688 |
| O-H2 | 3.0702 | 12.684 |
| Other components | 0 | 0 |

## Full constructor

`Cubic::with_options(..., CubicOptions { parameter_reference: "QuantumCubic",
..Default::default() })` forces TWU alpha, quantum beta and volume translation,
as in `eoslibinit::init_cubic`. It also applies the source's Tc/Pc/molecular-mass
overrides before calculating ideal reference offsets, alpha and EOS constants.
N-H2 and E-H2 retain the source's zero beta parameters and unchanged critical
data. Exact effective critical data and masses are included in reference tests.
The source's default-Cp lookup behavior is preserved.
`Cubic::quantum(names)` matches the separate Python `qcubic` helper: it selects
PR with the combined `QuantumCubic/tcPR/DEFAULT` reference priority and owns
a 2 to 999 K caloric interval. `temperature_bounds()` and
`set_temperature_bounds(min,max)` inspect/change model-local bounds. Native
`ph_flash`/`ps_flash` use those bounds; `caloric_options()` provides a starting
point for explicit per-call options. General cubic constructors still default
to 80 to 999 K. The expanded cryogenic diagnostic explicitly uses 0.5 to 999 K
in both implementations.
The public reference is limited to the source constructor's 100-byte buffer;
keywords beyond that buffer do not force a quantum model. Boundary cases are
qualified against Fortran. The isolated alpha initializer keeps its separate
reference-selection scope.

Each model now owns its component records, so `Cubic::components()` returns
`&[Component]` with effective critical data. The global `component()` and
`components()` functions still expose the immutable compiled database. This
changes the earlier experimental `&[&'static Component]` model accessor type;
ordinary field access is unchanged. Configured pure caloric submodels copy
effective critical data, beta, alpha, ideal parameters and translation settings.

Constructor overrides are implemented for the existing PR/SRK/RK geometries
with classical mixing. Other geometries/mixing rules remain in the coverage ledger.
The source rejects unsupported volume-shift geometry/fallback combinations;
native construction returns an explicit error for those combinations. Explicit
beta setters also work with the existing van der Waals geometry.

## Derivative propagation

Let `I(b,v)` be the cubic attraction integral and `A_b = R*T/(v-b) - a*I_b`.
In addition to the constant-b expressions, the thermal corrections are:

* `P_T += [R*T/(v-b)^2 + a*D_b/D^2] * bT`, with the cubic denominator D.
* `Hres -= T*A_b*bT` and `Sres -= A_b*bT`.
* `Cvres -= T*[2*(R/(v-b)-aT*I_b)*bT +
  (R*T/(v-b)^2-a*I_bb)*bT^2 + A_b*bTT]`.
* `Cpres = Cvres - T*P_T^2/P_v - R`, using the corrected pressure derivative.

Fugacity temperature derivatives include bT and component covolume-partial
derivatives. Mole-number derivatives use the current bi(T), including nonzero
lij mixing. Temperature-dependent translation is applied after these terms,
including its cross terms with the corrected pressure derivative.

## Qualification and remaining solver work

* 4036 parameter cases match active methods, coefficients, beta/betaT/betaTT,
  critical data and molecular masses. They cover every compiled fluid under
  default/QuantumCubic PR/SRK/RK constructors and explicit setters under all
  four existing geometries.
* 4512 property cases match liquid/vapor volumes, pressure derivatives,
  residual/total H/S, baseline Cp and fugacity T/P/mole-number derivatives.
  Cases include 3 to 700 K, 1 kPa to 10 MPa, cryogenic mixtures, hydrogen
  isomers, zero fractions, custom beta/ci, nonzero kij/lij and all four existing
  geometries. Every attempted parameter/property call returned a numerical
  result; all passed the recorded tolerances.
  Combined-reference PR cases call the actual Python `qcubic` helper as the
  oracle; single-reference/public constructor cases remain separately encoded.
* Independent finite differences check beta, P_T, Cp, Cv and fugacity T/P/n
  derivatives with simultaneous variable covolume and translation.
* 36 pure saturation states match pressure, both physical volumes and the
  Clapeyron slope: six quantum fluids, three geometries and two reduced
  temperatures. The saturation qualifier is independent of caloric flash
  temperature/phase-branch selection.

The expanded cryogenic TP/PH/PS grid remains unqualified. Some native calls
fail to converge, and some closing caloric specifications select a different
temperature branch. The source also returns higher-Gibbs stationary splits
and vanishing-phase outputs without chemical-potential closure in documented
cases. These observations are pending compatibility work, not successful
flash comparisons. The reference data and failing validator are retained for
continued solver development; the property qualification above does not imply
that all QuantumCubic flashes are interchangeable with the fork.
The retained single-reference grid contains 2161 finite/closing source returns
and 305 recorded source failures. Current native comparisons leave 82 cases
unresolved: 5 TP, 37 PH and 40 PS. The 36 pure saturation cases all match.

`validate_quantum_flash.py` records successful and failed native comparisons
separately, with source fingerprints, and returns nonzero while differences
remain. Its metadata also records source returns at a different temperature
than the state used to generate H/S; caloric closure alone cannot distinguish
multiple closing roots. [quantum_flash_audit.json](../../tests/fixtures/quantum_flash_audit.json)
contains independently recalculated source Gibbs and chemical-potential checks.
Use `python tools/rust_port/validate_quantum_flash.py --fixtures-only` to recheck
the frozen source returns without Fortran. Every current difference is retained
in [quantum_native_differences.json](../../tests/fixtures/quantum_native_differences.json).
