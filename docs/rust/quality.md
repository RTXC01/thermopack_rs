# Fixed-quality cubic states

The native cubic API now solves local two-phase states with a specified molar
vapor fraction beta. The source specification is `saturation::sat_newton`,
including the log-ratio coordinate used by `saturation_curve` continuation.
General fixed-quality envelope tracing remains pending; this checkpoint
provides its point corrector and warm continuation primitives.

## API

```rust
use thermopack_rs::{Cubic, CubicKind};
use thermopack_rs::envelope::EnvelopeConstraint;
# fn main() -> Result<(), thermopack_rs::Error> {
let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"])?;
let state = model.quality_temperature(1e6, &[0.5, 0.5], 0.4)?;
let next = model.refine_quality_state(
    &state, EnvelopeConstraint::Temperature(state.temperature + 0.1)
)?;
assert_eq!(state.feed, next.feed);
# Ok(())
# }
```

`quality_temperature` fixes pressure in Pa; `quality_pressure` fixes
temperature in K. Both require a normalized feed and beta in [0,1]. Zero feed
species stay exactly absent from both phases. Returned volumes are m³/mol,
and `dpdt` is Pa/K along the branch at fixed feed and beta. The original feed
is stored in `QualityState`, avoiding cumulative lever-rule drift when a
previous state becomes the next guess.

`quality_state_with_guess` accepts a `QualityGuess` containing T, P and y/x
ratios. Its constraint can fix T, P or one active component's log K, allowing
local continuation around pressure/temperature turns. Positive ratios are
required only for active species. `refine_quality_state` uses the stored feed,
quality and phase ratios. These are local root calculations: global VLLE
stability or selection among multiple separated branches is not established.

Cold initialization starts from the nearest bubble/dew boundary, with a
pseudo-saturation seed to preserve the source's branch choice for quantum
mixtures. It first attempts the requested quality directly, then tries sixteen
quality increments from the boundary if needed. Exact pure feeds use the
existing pseudo-saturation T/P solve and report the specified quality; log-K
constraints are undefined there. Those pure helpers retain the existing
saturation temperature-bound convention.

## Equations and derivatives

For log-ratio coordinate u_i and fixed beta, define

```text
K_i = exp(u_i)
X_i = z_i / (1 - beta + beta K_i)
Y_i = K_i X_i
F_i = u_i + ln(phi_V,i) - ln(phi_L,i)
F_RR = sum(Y_i - X_i)
```

The last equation fixes a selected log K, log T or log P. Both phase
compositions respond at interior quality. Fugacity evaluation normalizes X
and Y independently away from material-balance closure. The physical analytic
Jacobian includes their separate sums in the chain rule. Independent finite
differences exercise it away from the RR root, at beta=0/1 and interior beta,
with UNIFAC, HV and translated models.

The safeguarded physical Newton solve uses a maximum residual of 1e-10,
bounded log steps and backtracking. It rejects a trivial K=1 solution.
Phase preparation reuses the existing model-local envelope jets and their
analytic fugacity derivatives. Scratch vectors and trial buffers are reused
within each correction. No Fortran is called by the Rust solver.

## Source conventions and remaining differences

`baseline_quality_state_with_guess` preserves raw cubic root selection and
the fork's unscaled Newton steps, up to ten line-search halvings, post-update
L2 tolerance of 1e-8 and 500-iteration limit. The fork normalizes phase amounts
inside `eos::thermo` but returns composition derivatives at unit total. Its
quality Jacobian does not rescale those derivatives by the original phase
sums. The baseline path preserves that convention; the physical path uses
the normalized chain rule.

This distinction matters near coincident phases. Eight methane/water cases
at 350 K and beta=0.1/0.4 have very small phase-composition separation. The
physical cold initializer returns an error there while the source accepts a
near-trivial state. These eight calls remain explicit initialization gaps,
not qualified physical cold results.

In the local release build, 576 of 577 finite source-baseline returns pass
every strict coordinate comparison. PR NRTL methane/water at beta=0.4 has a
raw trace-composition difference of about 9.0e-8, exceeding the fixed ordinary
tolerance. It is retained in fixture metadata and is **not counted as strict
parity**. Development evaluation order also exposes a small composition miss
in the beta=0.1 NRTL cases. Regression tests retain bounded observations for
the eight near-coincident states; all other 569 states must pass the ordinary
tolerances. These observations establish neither an exact nonlinear
trajectory nor a rigorous conditioning allowance for those raw differences.

Full source failure/status identity, critical/coincident-phase regularization,
general fixed-quality curves, entropy/enthalpy constraints, non-vapor second
phases and other EOS families remain pending.

## Qualification

`validate_quality.py` performs 588 unique source attempts: 510 fixed T/P
initializations and 78 warm log-K constraints. Of these, 577 return finite,
closing source states. Eleven source failures are retained: six failed bubble
initializations and five solver bound statuses. No calls were censored, no
finite returns were discarded, and no returned volume-pressure warnings were
found.

The release baseline has 576 strict matches and the explicit raw difference
above. Physical cold initialization matches 491 of 499 finite fixed T/P
returns; the other eight are recorded errors. The matrix includes PR, SRK,
RK and VDW geometries, binary/ternary/zero feeds, beta=0/1 and interior values,
translated cubics, UNIFAC, HV1/HV2/NRTL and quantum H2/Ne and H2/D2.

Ordinary comparisons use relative 2e-7 and absolute floors 2e-6 K, 0.1 Pa,
1e-12 m³/mol and 2e-8 for mole fractions. Slopes use relative 2e-6 and
absolute 1e-3 Pa/K. No general tolerance is enlarged for the raw gaps.
Independent tests close phase pressure, fugacity and material balance; compare
fixed-quality slopes against warm finite differences; check beta=0/1 against
bubble/dew states; and exercise pure dispatch, invalid inputs and exact zeros.

The validator fails on unexpected mismatches. Its two narrow expected-gap
predicates preserve the recorded NRTL composition misses and the eight cold
errors, and always save their raw differences. A passing validator therefore
means no unexpected regression, not complete source parity.

## Performance

Local Windows release medians per correction were:

| Workload | Preflight states | Source | Rust | Speedup |
|---|---:|---:|---:|---:|
| PR methane/ethane | 17 | 116.37 µs | 15.04 µs | 7.74x |
| SRK methane/ethane | 18 | 117.02 µs | 16.14 µs | 7.25x |
| Translated PR methane/propane | 5 | 125.74 µs | 15.44 µs | 8.14x |
| UNIFAC PR methane/propane | 6 | 351.01 µs | 81.86 µs | 4.29x |

`benchmark_quality.py` compares native executables through source
`sat_newton`/`newton_extrapolate` and Rust
`baseline_quality_state_with_guess`. It uses identical source boundary guesses,
models and inputs. Preflight compares every returned T/P, quality, both
volumes, both compositions and slope on matching interior-quality workloads.

Construction, parsing, seed generation and Python/FFI are excluded. Both
executables perform the full corrector, phase reporting and analytic slope
calculation. Timed samples cycle 1000 corrections per profile and compare all
ten separate output sums. Five measured samples follow one warm-up. The
[report](quality-benchmark.json) retains source identity, inputs, samples and
source/harness/Cargo fingerprints. This does not measure cold initialization,
complete envelope tracing or whole-library performance.
