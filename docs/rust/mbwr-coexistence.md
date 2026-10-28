# Pure MBWR coexistence

`Mbwr::saturation_pressure(T, maximum_density)` solves vapor/liquid coexistence
in SI units. `saturation_pressure_bounded` also accepts an explicit pressure
window. `prepare_equilibrium(T, maximum_density)` shares prepared pressure
branches between density and coexistence calculations. `MbwrFluid` exposes the
same operations; ideal reference offsets cancel between phases.

`saturation_temperature(P, minimum_temperature, maximum_temperature,
maximum_density)` inverts the coexistence curve within the caller's bracket.
Both bracket endpoints must have resolvable coexistence states and bracket the
requested pressure. This method assumes the selected coexistence branch is
continuous and monotone within that bracket; arbitrary nonmonotone custom EOS
curves are not qualified. The returned pressure is the computed pressure, not
an overwritten input.

```rust
use thermopack_rs::Mbwr;
let model = Mbwr::from_record("mbwrcx1")?;
let sat = model.saturation_pressure(150., 45_000.)?;
assert!(sat.liquid_volume < sat.vapor_volume);
# Ok::<(), thermopack_rs::Error>(())
```

The solver reuses analytically isolated mechanically stable pressure branches.
It solves residual chemical-potential equality in log pressure with safeguarded
Newton steps and bisection, then checks pressure closure, positive pressure
slopes and Gibbs energy against all stable roots found in the explicit density
domain. The vapor branch is paired with each higher stable branch. Distinct
simultaneous global candidates cause an ambiguity error. The domain bounds
search; it does not establish the fitted EOS's physical validity at high density.
This is vapor/liquid equilibrium, not a general liquid/liquid or multiphase solver.

The analytic Clapeyron slope uses residual temperature derivatives; ideal
enthalpy cancels. Phase separation scales the chemical-potential tolerance to
avoid accepting coalesced roots near criticality. Unresolved endpoints return
an error. At machine-resolved density roots, severe pressure-term cancellation
can use a bound formed from the absolute original pressure terms. Ordinary
roots retain their existing fast path. This does not widen the external source
comparison tolerances. No speed improvement is claimed without measurements.

## Independent evidence

The committed fixture retains **1,008 requests** across all 65 stored fits:
390 physical-domain probes, 276 actual public source constructor/solver calls,
and 342 internally configured original-coefficient solver calls. These overlap
by operation but have distinct provenance. Known critical fits use independently
qualified critical seeds; missing metadata is explicitly supplied and does not
create a claim of default constructor support.

* **404 original solver matches** compare temperature, pressure, two volumes and
  Clapeyron slope against the original Fortran results. Relative tolerance is
  2e-6, with absolute floors of 1e-7 K, 1e-4 Pa, 1e-11 and 1e-10 m3/mol,
  and 1e-4 Pa/K respectively.
* **920 native returns** pass separate original pressure, residual Helmholtz,
  fugacity, mechanical-stability, chemical-potential and Clapeyron checks.
  Original kernels also check proposed stable roots for lower Gibbs energy.
  These are existence/equation certificates, not additional solver-coordinate
  matches. There are 1,012 distinct original-kernel invocations.
* **88 native error observations** remain unqualified. They are retained, not
  treated as correct error gold. All 245 subcritical physical probes with
  independently known critical data return; 49 supercritical and 26 probes
  using missing-critical metadata account for the physical cohort's 75 errors.

All **618 original solver attempts** remain in the archive: 404 qualified
matches, 96 initialization failures, 96 nonconverged calls, 14 coincident or
swapped-phase returns and eight in-domain metastable returns. The fourteen
collapsed results can have successful source status and finite T/P; a derived
Clapeyron slope then becomes undefined. Eight otherwise closed source states
have an independently demonstrated lower-G stable root in the selected density
domain. Their branch-selection behavior is recorded separately. Full compatibility
with the source solver's raw branch-selection semantics remains pending.
No requests were censored or silently retried in this record.

The CSV `original_solver` column contains numeric gold only for the 404 matches.
`native_observed` stores observations and is never used as numeric gold. The
standalone Rust replay checks all requests, requires the 920 certified inputs
to return, compares only independent gold, and allows either an error or a
physics-checked return for the 88 unqualified observations. Offline audit
recomputes classifications from raw process output and binds kernel requests,
root queries, cohorts, source/build identities, fixture hashes and counts.
The compressed audit evidence retains failed source calls as well as successes.

## Reproduction and remaining scope

Run `cargo test --release --test mbwr_coexistence_source` and
`python tools/rust_port/audit_mbwr_coexistence.py` for offline qualification.
The adversarial audit tests are selected by `test_mbwr_coexistence_audit.py`.
To recollect independent source evidence, first build the reference following
[the migration guide](README.md), then run:

```text
python tools/rust_port/qualify_mbwr_coexistence.py --compiler PATH_TO_GFORTRAN --source-build-dir .reference
```

The original reference is linked without modifying its source. Requests have
an eight-second process limit with termination/reaping; censored observations
are retained and fail qualification. Cached reuse verifies frozen request and
source/build identities. A new native implementation requires a new native
qualification record; historical hashes are never refreshed by hand.

The public source constructor refines active MBWR critical metadata through
`init_fallback_and_redefine_criticals` and `redefine_critical_parameters` before
saturation. It also updates the fallback cubic and acentric factor. The internal
oracle invokes those original routines with qualified seeds; it does not merely
substitute the database's coarse critical constants. See [baseline findings](baseline-defects.md).

Raw source saturation iteration/phase-selection compatibility, saturation curve
spacing and critical endpoints, flash integration, broader stability, CSP,
existing language interfaces and performance qualification remain open.
