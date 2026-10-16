# Volume translation and reference-aware cubic construction

`Cubic::with_options(kind, names, CubicOptions { .. })` selects source
parameter references for alpha, ci and classical binary interactions.
As in `eoslibinit::init_cubic`, a reference containing tcPR forces TWU alpha
and enables Peneloux translation even if the caller's flags request otherwise.
QuantumCubic changes beta and critical data as well; its constructor and
properties are qualified separately in [beta.md](beta.md). Expanded cryogenic
flash qualification remains pending.

`enable_volume_translation(reference)` selects the compiled ci record or
the source Rackett-factor fallback. `init_peneloux_volume_translation(ref)`
preserves the Python wrapper's behavior: its source drops the reference
argument and initializes DEFAULT parameters. These two APIs are deliberately
distinct. PR and SRK/RK support the source Rackett fallback; the source rejects
other geometries when it needs that fallback.

`set_volume_translation(i, parameters)` sets six coefficients and a constant,
linear, quadratic or quintic correlation. Setting coefficients does not enable
translation. `disable_volume_translation` retains the coefficients;
`set_volume_translation_enabled` controls their use. Settings are owned by the
model and cannot mutate another model or the embedded database. Indices are
zero based. All coefficients must be finite.

## Numerical convention

The physical volume v and unshifted EOS volume ve obey `ve = v + c(T,z)`, with
`c = sum(z_i*c_i)`. A prepared state caches c, cT, cTT and each component's
polynomial values/derivatives. The polynomial is evaluated with Horner's method.
Roots are solved/polished in ve and returned as physical volumes. Nonpositive
physical roots are rejected explicitly.

At fixed T,P, the corrections to the unshifted properties are:

* `G = G0 - P*c`
* `H = H0 + P*(T*cT - c)`
* `S = S0 + P*cT`
* `Cp = Cp0 + P*T*cTT`
* `ln(phi_i) = ln(phi_i0) - P*c_i/(R*T)`

At fixed physical volume, `P_T = P0_T + P0_v*cT` and `P_v = P0_v`.
The heat-capacity correction is
`Cv = Cv0 + 2*T*P0_T*cT + T*P0_v*cT^2 + P*T*cTT`.
Temperature/pressure fugacity derivatives include the corresponding ci terms;
the mole-number derivatives at fixed T,P are unchanged because the shift is
linear in component amounts.

Pure/fixed-composition coexistence pressures are unchanged by translation.
Returned phase volumes and caloric properties are shifted. Configured alpha
and ci coefficients are preserved when a pure caloric flash constructs its
internal saturation model. The source pseudo-critical helper continues to
report its unshifted estimated volume, matching `TP_CalcPseudo`; it is not a
complete mixture critical-point solver.

## Data and qualification

Numerical constants come from compiled `compdatadb.f90`/`mixdatadb.f90`, including
all 191 ci records, all 111 ideal Cp records and non-default classical PR/SRK
interactions. Reference priority follows the original exact-token/default
selection rules. Component metadata also embeds the compiled Rackett factor.

Cp selection has a separate source defect: `init_cubic` forwards a reference
of length 100 to `getCpDBindex`, whose database labels have length 40.
`isRef` calls Fortran `INDEX` without trimming; a longer substring cannot
match. The constructor therefore always takes the last `Default` Cp record,
including for empty references. Native construction preserves this behavior.
The two blank-label alternative Cp records for CO2/NC16 are embedded but are
not selected by this constructor. Empty/slash references are explicitly tested;
blank references do not select blank-label alpha/ci entries.

The linked Fortran probe reads active ci parameters, c/cT/cTT, correlation and
enabled/temperature-dependent flags. Its provenance is recorded with fixtures.

* 1977 parameter cases qualify all database records, Rackett fallback,
  polynomial setters, ignored wrapper references and disabled-state retention.
* 3120 property cases qualify full PR/SRK/RK constructors, tcPR behavior,
  constant through quintic shifts, pure/mixture/zero-fraction feeds, liquid/vapor
  volumes, pressure derivatives, residual/total H/S, baseline Cp and fugacity
  T/P/mole-number derivatives. Every attempted call returned a numerical result
  and passed the recorded tolerances.
* Independent finite differences verify P_T, Cp, Cv and fugacity derivatives;
  saturation/latent-heat quality and mixture PH/PS round trips pass.
* 3960 translated flash/saturation attempts return 3821 qualified numerical
  results: 1130 TP, 1293 PH, 1326 PS and 72 pure saturation states. Exact phase
  flags, both fractions/composition arrays, temperature recovery, physical
  saturation volumes and Clapeyron slopes match. Native caloric results also
  close H or S independently; two-phase material balances are checked.

The remaining 139 source attempts comprise 86 process terminations (22 TP,
22 PH and 42 PS) and 53 returned PH states that fail independent caloric
closure. Metadata records each request and result. These are unresolved cases,
not successful numerical comparisons; source-failure compatibility remains
pending. Qualification covers the recorded grid and does not establish all
solver convergence, phase-label or constructor combinations. See the coverage
ledger for the remaining model families, beta functions and interfaces.
