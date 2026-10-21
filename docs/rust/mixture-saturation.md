# Cubic mixture bubble and dew states

`Cubic::bubble_pressure(T,x)`, `dew_pressure(T,y)`,
`bubble_temperature(P,x)` and `dew_temperature(P,y)` solve separated VLE
branches with the specified phase composition fixed. Each returns
`MixtureSaturation`: T in K, P in Pa, both normalized compositions, both
molar volumes in m3/mol, a fugacity residual and the branch derivative dP/dT.
Input phase fractions must be nonnegative and sum to one within 1e-10;
accepted inputs are normalized and exactly absent components stay absent.

```rust
use thermopack_rs::{Cubic,CubicKind};
let eos = Cubic::new(CubicKind::PengRobinson, &["C1","C2"])?;
let bubble = eos.bubble_pressure(180.0, &[0.5,0.5])?;
let next = eos.refine_saturation_pressure(180.1, &bubble)?;
assert!(next.vapor[0] > next.liquid[0]);
# Ok::<(), thermopack_rs::Error>(())
```

The existing pure and pseudo-pure methods remain separate. A fixed phase
with one active component uses that pure limit, even in a larger component
list. No phase fraction or two-phase feed split is assigned at saturation.

## Equations and safeguards

The unknowns are the logarithms of active incipient phase amounts and the
unknown pressure or temperature. The equations enforce component fugacity
equality and unit incipient amount. Analytic fugacity derivatives in T/P/n
form the dense Jacobian. Inactive variables are removed instead of being
replaced by small positive amounts. Logarithms keep trial amounts positive.
Newton steps are limited and backtracked; accepting a residual also requires
a small predicted coordinate correction, which matters in flat equations.

Wilson estimates start ordinary solves. Pseudo-pure initialization and
alternative scalar seeds are fallbacks. Difficult temperature solves can
first follow a pressure boundary using its analytic dP/dT, then impose the
exact requested pressure with the coupled Newton solve. Fixed-temperature
queries reuse the preparation of the specified phase. Matrix, derivative and
residual buffers are reused within a solve, and line-search trials omit the
Jacobian calculation.

Automatic branch selection checks the trial-composition-weighted derivative
of the tangent-plane distance in the specification direction. It selects an
onset-oriented branch rather than an arbitrary retrograde root. This local
test does not establish global stability against all other phases.
The four `_with_guess` methods accept a scalar initial estimate; a scalar
alone does not uniquely identify a retrograde branch. The two
`refine_saturation_pressure(T, previous)` and
`refine_saturation_temperature(P, previous)` methods also use the previous
incipient composition and can follow another separated branch. The parent
composition and bubble/dew kind are taken from the previous state.

The returned dP/dT comes from differentiating the full saturation equations,
including the changing incipient composition. A mixture's slope generally
cannot be obtained by applying the pure-fluid Clapeyron expression directly.
`iterations` counts the successful final coupled solve; it excludes failed
initial guesses and preliminary continuation work.

Inverted/coincident phase volumes are rejected. When both relative volume
separation and maximum composition separation are below 1e-4, the result is
reported as requiring critical-point treatment. Tiny log-fugacity residuals
in that region can correspond to very different boundary coordinates.
There is no claim of critical-point convergence, complete envelopes, VLLE,
LLE, all retrograde branches or raw source iteration/status fidelity here.

## Exact-fork qualification

The matrix uses the original safe bubble/dew wrappers in isolated workers.
It attempts 1008 calculations across PR, SRK, RK and van der Waals, translated
PR/SRK, and UNIFAC PR/SRK. Binary and ternary hydrocarbon/CO2 mixtures include
asymmetric, zero-component and embedded pure compositions. Pressure solves
are followed by the corresponding source temperature inversion when a
qualified pressure exists.

All 959 separated boundary returns match T/P, both volumes and both phase
compositions at relative 2e-7 with absolute floors of 1e-6 K, 1e-2 Pa,
1e-12 m3/mol and 1e-9 composition. The largest observed relative differences
are 2.2e-9 in T, 2.0e-8 in P, 2.6e-8 in volume and 2.0e-8 in composition.
The metadata retains all 41 source failures and eight unqualified returns
(seven fugacity-closure failures and one coalescing result). There are no
deadline-censored calls. Of the 959 boundary returns, 21 have pressure-root
warnings from the separate source `specific_volume` routine. Those original
volumes are compared as returned and match within the same volume tolerances;
their warnings are retained instead of discarding valid T/P/composition
results. The other 938 also pass the source volume pressure-closure check.
Every qualified native result independently closes pressure and fugacity
and has negative dP/dv in each phase. Native results for 15 of the 49
failed/unqualified source attempts are recorded without a parity claim.

One retained source inversion illustrates why separation is checked: the
van der Waals C1/C2 dew-temperature solve returns almost identical phases
near 212.76 K despite a separated 220 K boundary at that pressure.
Its volume separation is about 1.1e-6 and its composition separation about
1.7e-7. Such a return does not identify a numerically resolved VLE boundary.

Native tests compare all qualified fixture states and independently check
pressure, fugacity and normalization closure; log-variable Jacobians by
finite differences; fixed-composition slopes; inverse specifications; pure
limits; absence of components; high-pressure initialization; both pressure
branches of a retrograde example and warm branch continuation. Other cubic
mixing rules, quantum models and global multiphase/critical behavior are not
qualified by this matrix. Hosted CI has not been observed here.

## Native solve measurements

The [benchmark record](benchmark-mixture-saturation.json) compares native
executable loops for PR/SRK C1/C2 at fixed parent composition 0.5/0.5.
Each of the four methods uses 100 points: T=160..199.6 K or P=0.1..3.07 MPa.
There are 1000 independent solves per sample, five samples after a warm-up.
All eight returned state fields match in 800 source and 800 native preflight
states; separate timed checksums include each field. Both implementations
initialize each solve automatically. Model construction, Python and FFI are
outside the timers. The source safe wrappers are followed by two volume calls;
native results already include those volumes and additionally compute dP/dT.

| Model/method | Source microseconds/solve | Rust microseconds/solve | Source / Rust |
|---|---:|---:|---:|
| PR BP | 200.72 | 18.50 | 10.85 |
| PR DP | 169.50 | 16.76 | 10.11 |
| PR BT | 171.80 | 55.43 | 3.10 |
| PR DT | 173.49 | 56.04 | 3.10 |
| SRK BP | 187.51 | 16.80 | 11.16 |
| SRK DP | 159.12 | 16.66 | 9.55 |
| SRK BT | 158.77 | 51.87 | 3.06 |
| SRK DT | 158.81 | 54.33 | 2.92 |

BP/DP denote bubble/dew pressure; BT/DT denote temperature. These measurements
cover this separated binary workload on the recorded Windows GNU host.
Translated/activity models, critical/retrograde fallbacks, warm refinement
and whole-library performance are outside this timing claim. Source, native,
harness and Cargo fingerprints are retained with the samples.

```
python tools/rust_port/validate_mixture_saturation.py --compiler-directory /path/to/compiler/bin
cargo test --test mixture_saturation
python tools/rust_port/benchmark_mixture_saturation.py --compiler /path/to/gfortran
```
