# Pure multiparameter TP, PH and PS flashes

`Multiparameter` provides `tp_flash(T,P)`, `ph_flash(P,H)`, `ps_flash(P,S)`,
`saturated_flash(T,vapor_fraction)`, and `equilibrium_properties(&flash)`.
The caloric methods also accept `CaloricFlashOptions` through their
`_with_options` variants. They cover the five NIST fluids, eight physical LJ
variants and reduced LJTS described in [the kernel notes](multiparameter.md).
Inputs use the model's units: SI for physical models, reduced units for LJTS
with R=1. H/S use the original **DEFAULT** caloric reference convention.
Alternative reference shifts remain pending.

```rust
use thermopack_rs::Multiparameter;
let propane = Multiparameter::new("C3")?;
let split = propane.saturated_flash(280.0, 0.3)?;
let q = propane.equilibrium_properties(&split)?;
let recovered = propane.ph_flash(split.pressure, q.enthalpy)?;
assert!((recovered.temperature - 280.0).abs() < 1e-5);
# Ok::<(), thermopack_rs::Error>(())
```

## Phase selection and solver bounds

TP compares dimensionless Gibbs energies of the mechanically stable,
pressure-closing liquid- and vapor-seeded source density iterations. It
removes the temperature-only ideal contribution before comparing roots.
It deduplicates coincident candidates and returns `Single` when only one
distinct candidate is available. This is a fluid branch selection policy;
it does not search all the extrapolated EOS's interior stable roots or solids.
Use [bounded pressure branches](meos-equilibrium.md) for that separate root
enumeration capability.

Pure coexistence T/P do not determine quality. At those points the source's
general Gibbs/stability solver may select a different endpoint than the Rust
Gibbs comparison, including differences between floating-point builds.
The qualification ledger retains those choices separately. Specified-quality
states are supplied by `saturated_flash`; source TP endpoint trajectory parity
is not claimed at coexistence.

PH/PS first invert the coexistence curve below the **stored** model critical
pressure, matching `getCriticalParam` after the source MEOS constructor
redefines component critical data. The local refined critical solver is not
the gate. Reduced LJTS has no stored positive Pc, so its gate evaluates the
EOS at its stored Tc/rhoc, equivalent to scaling the physical LJTS fit.

Targets between the liquid and vapor endpoint H or S use the lever rule.
Like the source pure solvers, this two-phase return bypasses homogeneous
temperature bounds. Homogeneous branches use a local expanding bracket,
Cp-based Newton updates and bisection. A stable log-volume solve anchored to
the corresponding coexistence endpoint avoids the source density iterator's
branch switching near an endpoint. Failed local continuation falls back to
the phase-seeded density iteration.

When temperature-only iteration reaches its precision limit, a bounded
coupled Newton correction solves pressure and H/S in log T/log rho. Its
analytic Jacobian is checked by finite differences. Accepted states require
positive pressure and mechanical stability; the caloric tolerance remains
H/(R T) or S/R, using the model's R. The coupled pressure closure uses the
scale |P| + R T rho and a roundoff guard, including reduced units.

Physical default homogeneous bounds are the source's 80..999 K. Explicit
options are useful for cryogenic hydrogen. The native reduced-LJTS defaults
are 0.3 Tc..3 Tc; this is a native convenience, since this fork has no public
reduced-model flash wrapper. Explicit option bounds always use model units.
Out-of-bound homogeneous targets return an error. Scalar iteration counts
exclude the nested coexistence/density work; latent returns have count zero.

## Properties and qualification

`equilibrium_properties` returns molar V/H/S/U and optional homogeneous Cp.
It recomputes phase volumes for this model rather than trusting cached
`FlashResult` volumes. Each caloric state evaluates ideal/residual jets once;
singular/nonpositive Cp becomes `None` without discarding finite H/S/U.
Two-phase Cp is `None`. Source-style `Single` fractions are -1/-1; homogeneous
PH/PS compositions are both [1], as in the source pure solvers.

Reproduce the audit with:

```text
python tools/rust_port/validate_meos_flash.py --allow-known-gaps
python tools/rust_port/verify_fixtures.py
cargo test --test meos_flash
```

Build the exact-source library and probe first, as described in
[the migration guide](README.md). On Windows pass the compiler runtime
directory through `--compiler-directory`. Source caches bind the library,
worker, generator and seed fixture hashes; `--reuse-source` rejects stale
caches. Each source call has a four-second deadline. Timeouts are censored
observations, killed/reaped without retrying the same case.

The grid includes subcritical/supercritical TP seeds, guessed and unguessed
PH/PS, specified qualities including both endpoints, default bounds and
independent source TV-kernel targets at 0.999 and 0.9999 Tc. Those direct
targets let the actual PH/PS solver run even when source bubble-pressure
seed preparation fails. Reduced LJTS uses the identical physical NE LJTS
fit with explicit T/P/V/H/S scales; no original reduced flash wrapper is
claimed.

The main metadata retains every request and raw source phase flag,
composition, fraction, warning and terminal error. Seed preparation failures
are distinguished from executed flash failures. Only finite, mechanically
stable, pressure-closing and caloric-closing source returns qualify for
coordinate comparison. Interior two-phase canonicalization uses a 1e-8
fraction threshold; raw amounts remain in the ledger. Matched coordinates,
coexistence TP endpoint choices and source-kernel-qualified native recoveries
have separate fixtures. Remaining near-endpoint coordinate differences are
explicit gaps, rather than silently widened tolerances.

Every finite native result in this audit is independently evaluated through
the original pressure/caloric TV kernels. Recovery rows establish native
physical/specification closure where an original solver fails or returns a
nonclosing state; they do not establish original solver coordinate parity.
The source-kernel normalized caloric and specification checks use 2e-10.

The committed audit contains **2,100 source attempts**: 1,708 coordinate
matches, 36 explicit near-endpoint coordinate gaps and 46 underdetermined
coexistence TP choices among 1,790 qualified source returns. The remaining
source outcomes are 210 saturation-seed failures, 44 PH errors, two terminal
PS density-solver errors and 54 nonclosing returns. There were no censored
calls. All 100 executed failed/nonclosing caloric specifications have native
recoveries qualified through the source TV kernels. Altogether **1,890 finite
native flash states** passed those independent kernel checks; there are no
native solver errors on the qualified source grid.

The 36 coordinate gaps concern near-endpoint phase/fraction representation:
R134A at 0.999 Tc, full LJ Ar/Kr/Xe at 0.9999 Tc and all five LJTS variants at
0.999 Tc, with both caloric specifications and guess modes. They pass the
bulk volume, temperature, H/S and source-kernel closure checks, but do not
meet the strict phase/fraction/trace-phase-volume coordinate criteria.
`meos_flash_known_gaps.json` declares exact grid keys and source outcomes.
`--allow-known-gaps` also enforces per-field/bulk/fraction bounds; it rejects
new case keys, native errors, coincident phases, unlisted source outcomes
and censored calls. A successful audit in this mode is explicitly a partial
surface qualification, not an assertion that every original call succeeded
or that all outputs are interchangeable.

No TP/PH/PS speedup measurement is claimed in this change. Allocation-free
scalar trial states and shared jets remove avoidable work, but end-to-end
timing against the source remains a separate task. Whole-fork replacement,
all-root/global-stability policy, other MEOS flashes, response derivatives,
phase envelopes and existing language interfaces remain incomplete.
