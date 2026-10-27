# Total MBWR TP properties and density integration

`MbwrFluid` now connects the [total TV caloric layer](mbwr-caloric.md) to
local TP derivatives, bounded stable pressure branches and the fork's original
density iteration. These are pure homogeneous states. Coexistence and flashes
remain pending.

```rust
use thermopack_rs::{MbwrFluid, MbwrKind, Phase};
let fluid = MbwrFluid::new(MbwrKind::Mbwr32, "C1")?;
let prepared = fluid.prepare_tp(150., 30000.)?;
let vapor = prepared.caloric_tp(1e6, 1., Phase::Vapor)?;
assert!((vapor.pressure - 1e6).abs() < 1e-3);
# Ok::<(), thermopack_rs::Error>(())
```

## Physical TP evaluation

`caloric_tvp(T,V,n)` evaluates a supplied total volume and maps the TV gradient
to TP with V_T=-P_T/P_V, V_P=1/P_V and V_n=-P_n/P_V. Positive EOS pressure is
required. The result contains total V, U/H/S/A/G, chemical potential, Z and
ln(phi), with optional T/P/n derivative records, plus Cv, Cp and local sound.
Energies are J, entropy and extensive heat capacities J/K, chemical potential
J/mol, pressure Pa, volume m3 and molar density mol/m3. T/P/n partials keep the
other two coordinates fixed. There is no mixture-composition interpretation.

`prepare_tp(T,maximum_density)` prepares thermal coefficients, ideal values and
the bounded pressure branches once. Its `caloric_tp(P,n,phase)` selects the
lowest/highest stable density for vapor/liquid or the minimum residual Gibbs
root for `Phase::MinimumGibbs`. Total volume is n times the selected molar
volume. The convenience method takes T/P/n/phase/maximum_density. The upper
density is an explicit search boundary, not an EOS validity limit. Minimum
Gibbs selection covers the roots in that domain, not other phases or solids.
Returned pressure is the evaluated EOS pressure; the existing bounded solver
controls closure against the requested pressure.

At roundoff-unresolved P_V (64 machine epsilons times nRT/V^2), V/U/H/S/A/Z
derivatives are absent while finite values remain available. Exact pure-fluid
identities preserve finite first G, chemical-potential and ln(phi) derivatives
at this limit. Fugacity derivatives use the shared residual jet, avoiding the
cancellation in total H minus ideal H and in Z minus one at extreme dilution.
Each state shares one residual evaluation across its caloric quantities.

`baseline_caloric_tvp` maps the source-order TV Hessian in the same way, retaining
the original Cp derivative convention. Explicit type-10 ideal H2 metadata has
the previously documented factor-1000 H/S-primitive discrepancy; physical
derivatives follow the returned values, while this source TVP view follows the
unchanged Fortran Hessian. Raw zero ideal integration constants, reference/Cp
setters and missing-metadata errors are inherited from the TV layer.

## Classical source TP evaluation

`baseline_caloric_tp_at_volume(T,P,V,n)` evaluates total H/S and residual Z and
ln(phi) at an explicitly imposed P. P need not close the EOS. It preserves
the actual `single_component.f90` behavior: total H/S values include ideal
contributions and T derivatives add source Cp, but the fork incorrectly writes
the ideal P/n derivative additions into local buffers instead of the outputs.
Consequently the returned S_P and H_n/S_n omit their ideal contributions.
The separate sign defect in `MBWR_Sres` is also retained. These derivative
conventions are not thermodynamic identities.

`prepare_baseline_tp(T,MbwrDensityParameters)` caches the source density
coefficients and ideal values. Its `caloric_tp(P,n,phase)` returns
`MbwrSourceTpState { iteration, properties }`. The raw density diagnostics
include the reported phase, path, EOS pressure, closure and mechanical slope;
`properties` is a separate `Result` so a failed caloric assembly cannot discard
a finite source iterate. Convenience methods accept explicit initialization or
use the original cached default initialization, including its known errors.
The source iteration accepts liquid/vapor flags. Use the bounded physical API
for minimum-Gibbs selection.

Non-finite native derivative triplets are represented as absent records.
The existing raw `Mbwr::baseline_residual_tp` API retains its original behavior. All new total-state APIs use actual n and
extensive H/S. The fork's public `eos` TP wrapper normalizes composition and
returns molar properties; direct `single_component` calls qualify the extensive
interface here.

## Source evidence and unresolved derivative blocks

The committed matrix retains **9,752 requests and 25,174 native outputs**:

| Cohort | Requests | Evidence |
| --- | ---: | --- |
| Positive-pressure TVP | 4,959 | Every positive-pressure state from the earlier 65-fit TV matrix; independent original TVP H/S/ln(phi) calls plus complete original TV jets |
| Nonpositive EOS pressure | 1,035 | Both TVP APIs reject; the explicitly imposed positive-P classical values are checked separately |
| Classical source density/TP | 3,434 | All 65 fits, liquid/vapor flags, near-critical and supercritical requests, scaled mole amounts and reference/Cp edits; actual original `enthalpy_single`/`entropy_single` dispatch |
| Bounded TP integration | 324 | Methane, propane and R134a, both endpoints and minimum Gibbs, two mole scales and six reduced temperatures; native-proposed volumes independently checked with original EOS/TVP routines |

The source density cohort includes 113 unclosed returns, 1,175 reported phase
changes and 884 requests with explicitly supplied argon ideal/synthetic
initialization data for missing legacy metadata. They remain in the matrix.
They qualify algebra and source behavior, not missing-fluid physical accuracy.

**93 requests retain unqualified H/Z/ln(phi) derivative blocks.** They are
selected solely from the original density results: positive finite rho,
|P_rho| <= 64 machine epsilons times RT, and failed pressure closure. The
classification does not depend on native mismatches or record names. Large
implicit derivatives at these unclosed extrema vary with compilation context.
All original derivative values remain in the golden vectors and audit.
Native derivative triplets can be finite or absent; they are not counted as
qualified numerical parity. Density/phase outputs, all H/S/Z/ln(phi) values
and every S derivative still face the strict comparison.

For example, legacy fit `mbwrcx21` with explicit AR ideal data and
(Tc,Pc,rhoc,acf)=(300 K,5 MPa,10000 mol/m3,0.1), T=299.97 K, P=4 MPa, n=1
and vapor input returns rho=4197.516024649706 mol/m3. Its original EOS pressure
is about -0.856 MPa and P_rho about 1.68e-11 Pa/(mol/m3). The fork returns
huge finite H/Z derivatives; the recorded native response has absent derivative
records. Other build contexts can produce different finite triplets.

The archive records all 93 cases in `unqualified_source_derivatives`, with
actual vector differences separately retained in `known_mismatches`.
**9,659 requests are fully qualified**; the other 93 have the declared twelve
unqualified coordinates. Coincidental agreement in a snapshot does not make
those derivative blocks parity successes. No request is omitted and no
numeric tolerance was widened.

Numeric comparisons use 2e-7 relative plus 2e-10 times natural ideal scales;
optional flags are exact and sound has a 2e-8 m/s floor. Raw source density uses
the existing 2e-9 comparison and exact phase flags. Source/objects/modules,
compiler, helper, generator dependencies, executables, schedule and native
inputs are hash-bound. The gzip archive retains complete process output and
all vectors; the CSV replays independently of Fortran. Timeouts are censored,
retained and fail qualification without automatic retries. Offline negative
tests reject omitted requests, paired gold/native tampering, phase changes and
false conditioning classifications.

```text
python tools/rust_port/qualify_mbwr_total_tp.py --compiler /path/to/gfortran --source-build-dir /path/to/reference
python tools/rust_port/audit_mbwr_total_tp.py
python -m unittest discover -s tools/rust_port -p test_mbwr_total_tp_audit.py
cargo test --release --test mbwr_total_tp --test mbwr_total_tp_source
```

No speed claim is made for this layer before native-loop measurements.
MBWR coexistence/flashes, broader CSP integration, existing language bindings
and the unresolved source derivative blocks remain within the full migration goal.
