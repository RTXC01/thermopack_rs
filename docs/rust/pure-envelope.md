# Pure-fluid dispatch through the envelope API

`Cubic::phase_envelope` and `baseline_phase_envelope` now handle pure feeds,
including a selected component embedded in a multicomponent model. The
mixture continuation path remains available for mixtures. The common return
type retains full phase compositions, both volumes, slopes and an explicit
termination reason.

## Physical and source behavior

| Behavior | `phase_envelope` | `baseline_phase_envelope` |
|---|---|---|
| Pure-feed classification | Exactly one positive mole fraction | Fork's `max(z) == sum(z)` |
| Grid | Linear in specified T or P; requested pressure spacing 20 kPa | Same source grid |
| Point cap of at least two | Resamples the complete pure interval, then applies bounds | Resamples to the database critical point |
| One-point cap | Initial-state prefix with `PointLimit` | Only the database critical endpoint |
| Critical endpoint | EOS pressure inflection, `CriticalPoint` | Database T/P and raw roots, `SourceCriticalParameters` |
| Maximum pressure | Closes the boundary, `MaximumPressure` | Ignored by the source pure branch |
| Model maximum temperature | Closes the boundary, `MaximumTemperature` | Ignored by the source pure grid |
| Minimum-temperature override | Validates the initial state against the clamped lower bound | Ignored by the source pure grid |
| Maximum mixture log step | Unused for pure grids | Ignored by the source pure branch |
| Initial bubble/dew kind | Returned pure points use bubble labels | Fork always returns beta=0 and K=1 |

The independent [pure-curve API](pure-curves.md) additionally offers
logarithmic grids. The general source envelope API does not select that
sampling mode. Its maximum point count can enlarge the actual grid spacing.

Physical tracing selects the first encountered pressure or upper-temperature
boundary. An initial state outside its bounds is an error. Accepted prefixes
are retained on later failures. Minimum temperature does not terminate an
ascending pure curve; it constrains the start. Source initialization can
still fail under its global solver bounds even though the pure grid ignores
the envelope's optional temperature override.

For example, the source pure-water branch continues to its database critical
pressure even when the envelope maximum pressure is 15 MPa. The physical API
stops at the requested 15 MPa boundary. A point cap of one likewise represents
different source and physical conventions, as shown above.

## Critical starts and source initialization

A physical start within 32 machine epsilons of the EOS critical T or P returns
one coincident-phase critical node. The dedicated physical pure-curve API has
the same critical-start convention; it still requires a point budget of at
least two for ordinary grid construction.

The source can accept a coarse metastable `PureSat` estimate at exactly Tc.
It then sees a zero fugacity difference when the ordinary roots coalesce.
Consequently its first point can be a single-phase state below Pc before the
database endpoint is appended. For SRK methane, the source starts at
4,598,811.687 Pa and appends 4,598,837 Pa. PR CO2 has an analogous difference
of about 41.75 Pa. Substituting Pc directly misses those source arrays.

The native baseline initializer reproduces the source's metastable cubic
extremum selection, coarse energy tolerance and line search. It preserves
the source ordering of the pseudo-critical calculation in L/mol, including
the grouping of squared factors in its omega iteration. A one-ulp change in
the pseudo-critical temperature can otherwise select a different initializer
branch. This source-order helper is separate from the established SI kernels.
Physical root discovery is unchanged.

Critical-start qualification here covers C1/CO2 with the four classical models.
Broader alpha/GE/quantum critical starts and complete raw failure/iteration
identity still require qualification. Native explicit outcomes for failed
source calls are not counted as matching successful source results.

## Rounding-level trace components

For feeds such as `[1.0, 1e-18]`, the source's sum rounds to the largest mole
fraction, so it takes the pure branch. Its temperature/pressure grid follows
the dominant component while K remains one for every component. The native
baseline preserves the normalized original feed in both phase arrays and
uses that feed for the reported raw volumes and fugacity residual.

Seven such source cases are retained in this matrix. Their trace-component
fugacity differences can be large: the largest observed maximum residual is
about 9.48. These source-classified curves are not established multicomponent
VLE equilibria. Physical dispatch counts every positive species and keeps
these feeds on the mixture path; convergence at that degeneracy is separately
qualified. No positive trace species is silently discarded to force a pure
physical result.

At one VDW methane/ethane critical node, the derived full fugacity residual is
0 in the source and approximately 1.72e-6 natively. The allowed critical-root
rounding changes the trace-component diagnostic even though both returned
volumes satisfy their declared source-polynomial bounds. This difference is
retained separately in metadata. The original envelope routine does not return
that residual, and ordinary raw diagnostic parity is not claimed there.
Each native residual is independently recomputed from its reported volumes.
Pure-dispatch `iterations` is zero because the shared scalar grid does not
retain per-node iteration counts. Slopes follow the selected pure component,
including when baseline classification approximates a rounding-level trace feed.

## Qualification

`validate_pure_envelope.py` attempts 180 source calls: 169 return finite
arrays and 11 terminate in the source. There are no censored or nonfinite
returns in this matrix. The successful set comprises **153 general-envelope
calls and 16 dedicated pure-curve edge calls**, with **13,917 stored states**.

Every returned count, all T/P and both phase volumes pass the documented
tolerances, including the source critical-root conditioning convention. All
general-envelope phase labels and normalized composition arrays match exactly.
The cases cover four classical geometries, translated models, selected
components with UNIFAC/HV1/HV2/NRTL, quantum models, critical T/P starts,
one-point grids, ignored controls and rounding-level trace feeds.

There are 175 nodes at database-critical T/P. All 350 phase-volume comparisons
check the source cubic polynomial backward error and independently recompute
the displacement bound; 348 exceed ordinary volume tolerance and need that
explicit allowance. All other coordinates retain relative 2e-7 with absolute
floors 2e-6 K, 0.1 Pa and 1e-12 m³/mol. Coarse critical-temperature initializer
states below Pc are ordinary comparisons, with no critical-volume allowance.
Their coordinates pass the same fixed tolerances.

Independent tests cover physical boundary ordering/closure, pressure
inflections, full phase arrays, one-point semantics and trace classification.
Mixture envelope and existing pure-curve fixtures remain regression checks.
The source failure records and the critical diagnostic difference are part of
the saved evidence, not successful full-library parity claims.

## Performance

The local Windows release medians per complete curve were:

| Workload | Points | Source | Rust | Speedup |
|---|---:|---:|---:|---:|
| PR methane, pressure start | 100 | 5.519 ms | 0.975 ms | 5.66x |
| SRK CO2, temperature start | 100 | 5.317 ms | 0.840 ms | 6.33x |
| Translated PR water, pressure start | 100 | 6.221 ms | 1.055 ms | 5.90x |
| Quantum H2, pressure start | 64 | 4.347 ms | 0.948 ms | 4.58x |

`benchmark_pure_curve.py --envelope` times native executable loops through the
general envelope APIs. PR methane, SRK CO2, translated PR water and quantum H2
produce identical 100/100/100/64-point source/native grids. All 364 states are
compared before timing. Source beta/K and native labels/compositions are
checked; eight critical phase volumes use the same explicit conditioning.
Each timed sample compares counts and separate T/P/volume sums.

Cold initialization, complete continuation, the source critical endpoint,
both volumes and native phase/residual reporting are included. Model
construction and Python/FFI are excluded. Five measured samples follow one
warm-up, with twenty complete curves per sample. The
[report](pure-envelope-benchmark.json) records every sample, source identity
and native/harness/Cargo hashes. These are scoped workloads, not whole-library
performance or proof of complete raw iteration fidelity.

General fixed quality, global VLLE/solid stability, other EOS families and
broader raw failure behavior remain in the [coverage ledger](coverage.md).
