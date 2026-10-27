# Total MBWR caloric state and derivatives

`MbwrFluid` owns a coefficient model, ideal-gas data and molar mass. It adds
total homogeneous U/H/S/A/G, chemical potential, pressure, heat capacities and
local sound response to the [residual kernels](mbwr.md). Pressure and each
potential have analytic first T/total-V/n derivatives. `reduced_helmholtz` returns all ten
value/first/second derivatives of extensive F=A/(RT), including mixed and mole
derivatives. These evaluations solve no density or phase-equilibrium problem.

```rust
use thermopack_rs::{MbwrFluid, MbwrKind};
let fluid = MbwrFluid::new(MbwrKind::Mbwr32, "C1")?;
let prepared = fluid.prepare(300.0)?;
let q = prepared.caloric_tv(1e-3, 1.0)?;
assert!((q.enthalpy.value-q.internal_energy.value-q.pressure.value*1e-3).abs()<1e-8);
# Ok::<(), thermopack_rs::Error>(())
```

`new` retains the coefficient lookup's aliases and first-fit selection;
`from_record` selects an exact fit. **48 of the 65 stored fits** have matching
canonical general component/Cp data. Their constructors use that data without
running the original critical-grid initialization. The other seventeen lack
default ideal metadata and return an error. `with_ideal_gas` accepts explicit
owned ideal data and molar mass in kg/mol for those fits and custom models.
It does not invent missing fluid data. Clones own independent mutable ideal
parameters; borrowed prepared states cache the ideal values and MBWR thermal
coefficients once per temperature. Each state shares one residual snapshot
and one ideal snapshot across its properties and gradients.

## Units and source conventions

T is K, V is total m3 and n is mol. U/H/A/G are J, S and Cv/Cp are J/K,
pressure is Pa, and chemical potential is J/mol. `MbwrTvQuantity` holds a
value and derivatives with respect to T, V and n at fixed remaining inputs.
The chemical potential is intensive; the five caloric potentials are
extensive. Positive finite T/V/n are required. Finite negative-pressure or
unstable mathematical EOS states remain observable. Their evaluation does not
establish fluid stability, fit validity or physical accuracy in extrapolation.

The original MBWR constructor copies **raw zero ideal H/S integration
constants**. It does not perform the 298.15 K reference normalization used
by other model constructors. For example its C1 ideal H at 298.15 K is
10025.726002805888 J/mol, and S at 1 bar is 112.14282031044416 J/(mol K).
Native constructors preserve this behavior. `IdealGas::with_integration_constants`
also exposes explicit finite constants. Existing formation-enthalpy and
standard-entropy setters perform the original reference normalization; a Cp
setter preserves integration constants, so changing the correlation can change
H/S at the reference temperature. Setter ordering is qualified explicitly.

`caloric_tv` and `reduced_helmholtz` use the stable integral backend and actual
derivatives of the returned ideal values. `baseline_caloric_tv` and
`baseline_reduced_helmholtz` use original residual summation and Cp conventions.
For the source's type-10 ideal-H2 correlation, H/S primitives carry a factor
1000 absent from its Cp routine. The physical view uses the actual H derivative
in F_TT; the baseline view reproduces the original inconsistent F_TT/Cv.
This affects explicitly supplied type-10 metadata; the 48 default MBWR fits
have other Cp types. These TV derivatives also differ from the separately
preserved sign error in the source [TP residual entropy](mbwr-properties.md).

Cp is optional at a singular or roundoff-unresolved P_V; Cv and finite
potentials/gradients remain available. The unresolved threshold is 64 machine
epsilons times nRT/V². `ln_fugacity_coefficient` is optional for nonpositive
pressure. Sound is optional unless its local squared response is finite and
positive. The source's zero return for an invalid sound response maps to
`None`; its original raw return is retained in the audit. A real local sound
response alone does not prove thermodynamic stability.

The physical fugacity calculation uses `log1p(Z-1)` from the residual jet in
the dilute region, retaining the first virial term even when total pressure
rounding makes Z equal one. The source-compatible view retains its original
pressure-ratio expression. A separate low-density regression checks this limit.

## Qualification

**5,994 explicit TV states** pass both views, or 11,988 native state comparisons.
The coefficient matrix covers all 65 fits at six temperatures (75..900 K),
five densities (0.01..15000 mol/m3) and three mole amounts. Missing-data fits
use explicit argon ideal metadata for algebraic checks, not a claim that the
missing fluid has argon's properties. Additional cases exercise reference/Cp
setter ordering and explicit hydrogen-isomer metadata. The total includes
72 setter cases, 24 type-10 cases and 1,035 nonpositive-pressure states; none
are discarded. Forty-five coordinates include the complete total Helmholtz
jet, pressure and all six thermal value/gradient blocks, Cv and optional responses.

A standalone Fortran program uses unchanged `eostv` dispatch, `MBWR_Fres`
and `Fideal` routines from the exact reference objects. It only allocates a
selected coefficient record and fresh component data, bypassing unrelated
critical initialization. It resets component storage for every request to
prevent previous setter constants from leaking into later cases. Thirty
callable public constructors produce matching full TV snapshots; sixteen
known constructor terminations are retained separately. Alternate fits cannot
all be selected by the public first-match constructor. This kernel-backed
thermal qualification does not claim public critical initialization succeeds.

Physical expectations are derived from the **original** Helmholtz Hessian,
with the explicit type-10 Cp correction, never from native-generated gold.
Independent finite differences check all six potentials' T/V/n gradients for
all 48 default-capable fits, the total Hessian, scaling and thermodynamic
identities. Comparison tolerance is 2e-7 relative plus 2e-10 times each
quantity's natural ideal scale; optional flags are exact and sound has a
2e-8 m/s floor. No tolerance was widened to hide differences.

The CSV holds both independent expected views. The deterministic gzip audit
archive retains all source/native vectors, complete process output, constructor
outcomes, original object/module/source hashes, compiler identity, executable
hashes and native input hashes. A small metadata file binds its SHA-256 and
the CSV. Offline audit and negative tests reject omitted cases, changed source
vectors, and substituting source Cp for the physical correction.

```text
python tools/rust_port/qualify_mbwr_caloric.py --compiler /path/to/gfortran --source-build-dir /path/to/reference
python tools/rust_port/audit_mbwr_caloric.py
python -m unittest discover -s tools/rust_port -p test_mbwr_caloric_audit.py
cargo test --release --test mbwr_caloric --test mbwr_caloric_source
```

Source reuse validates the coefficient objects, compiler, helper, generator,
constructor-status fixture and exact schedule. Timeouts are retained as
censored observations and fail qualification without retries. Historical
source-executable identity remains attached to a reused collection.
The original Python internal-energy wrapper mispackages its mole derivative
as a volume-derivative pointer; the probe calls the Fortran routine directly.

The [total TP layer](mbwr-total-tp.md) now connects these calorics to local TP
derivatives and bounded/original density selection. MBWR coexistence/flashes,
broad CSP integration and existing language interfaces remain pending.
Performance has not yet been measured for this layer. It evaluates total
homogeneous caloric states and does not replace the complete MBWR equilibrium API.
