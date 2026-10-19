# MBWR default initialization and source critical grid

`Mbwr::source_initialization()` follows the fork's `readDbParameters` metadata
path: look up general component data, select the stored MBWR critical table when
available, otherwise run `MBWR_criticalParameters`, then calculate the original
critical compressibility and SRK density-seed constants. It returns an owned
`MbwrInitialization`. `source_density_parameters()` returns its parameter block
for the [source density solver](mbwr-density.md).

```rust
use thermopack_rs::{Mbwr, MbwrKind, Phase};
let model = Mbwr::new(MbwrKind::Mbwr32, "IC4")?;
let initialization = model.source_initialization()?;
assert!(!initialization.from_critical_table);
let solver = model.prepare_baseline_density(300.0, initialization.parameters)?;
let state = solver.density(1e5, Phase::Vapor)?;
# Ok::<(), thermopack_rs::Error>(())
```

Initialization is lazy: constructing coefficients does not trigger a critical
search. A model-local `OnceLock` caches either the successful result or the
error, including across cloning. Models share no global active EOS. The cache
is valid because model coefficients/record metadata are immutable; explicit
custom density parameters remain independently supported.

The parameter block uses K, Pa and mol/m3. `compressibility` is dimensionless;
`srk_covolume` uses m3/mol, `srk_alpha_coefficient` is dimensionless, and
`srk_attraction` uses Pa*m6/mol2. `from_critical_table` distinguishes stored
metadata from the coarse search. `source_critical_point()` continues to expose
only the explicit table; `critical_point()` separately refines the physical
inflection equations. The initializer does not silently substitute the latter.

## Coarse grid semantics

`baseline_critical_grid(initial)` ports the original critical search explicitly.
Its temperature window is initial Tc +/- 0.7 K and density window is initial
rhoc +/- 1500 mol/m3. It uses repeated 0.001 K and 1 mol/m3 increments, not
index-derived coordinates, preserving the original floating-point accumulation.
The lower scan starts with P_old=0 and the upper scan with P_old=1e12 Pa. It
retains the original first-decline/first-ascent rules and averages the last lower
and first upper estimates. The supplied initial Pc and acentric factor are not
used by the pressure-grid equations; tests verify this independence.

The implementation prepares source-order pressure coefficients once per
temperature and evaluates only pressure in the inner loop. It reuses the native
integer-power cache, omitting both unused derivative calculations and unrelated
SRK preparation. No timing speedup for initialization has been qualified yet.

The source can exit the first lower-temperature scan without ever assigning
its lower-bound variables. Rust reports an explicit error for that undefined
result. It also rejects non-finite pressures, unresolvable floating-point grid
increments, invalid domains and nonphysical final estimates. An out-of-window
search returns an error corresponding to the source termination; it does not
expand the source's search window or replace its result with another solver.

## Original missing-component defect

`getCompDBindex` returns zero when no general component record exists.
`readDbParameters` then accesses `compdb(compDbIdx)%acf` without checking that
index. Seventeen legacy MBWR identifiers therefore have an out-of-bounds default
initializer, despite having usable pressure coefficients. The initial isolated
audit observed an arbitrary finite F2 initialization and an unbounded search
for R113. Those undefined behaviors are not reproducible numerical targets.

The validation probe calls the original lookup first and reports this condition
before executing the undefined access. Rust reports a missing-component error;
its coefficient-only kernels, physical critical solve with explicit guesses,
and density solve with supplied metadata remain available. No arbitrary zero
acentric factor or substitute critical data is installed as a default.

## Qualification

`validate_mbwr_initialization.py` covers all 65 records:

* All 32 finite original default initializations match eight outputs: Tc, Pc,
  rhoc, Zc, acentric factor and the three SRK constants. Eighteen use a critical
  table and fourteen use the coarse grid.
* All 16 original Fortran terminations have native error results. Their raw
  exit codes and source output remain in fixture metadata.
* All 17 undefined missing-component initializers are diagnosed by the probe
  before `compdb(0)` and have native errors. They are not numerical successes
  or original source termination counts.
* All 27 direct original critical-grid searches match, including stored-fit
  seeds and T/rho offsets for both methane families. Their output is the source
  coarse estimate, not the analytically refined critical point.

Relative tolerance is 2e-10, with quantity-specific absolute tolerances in the
fixtures. The standalone snapshot executable runs in the oracle's process-main
floating context. Normal tests verify lazy initialization, cached success/error,
cloning, coefficient-only availability after failure, unused-input independence,
uninitialized source bounds and stagnant grid increments. Broader host/runtime
rounding behavior is not a bitwise portability promise.

Recorded source terminations keep their raw exit code and backtrace frame count,
but the per-frame load addresses are replaced by a placeholder. Those addresses
depend on address-space layout randomization, so retaining them would make the
fixture metadata differ on every regeneration. Regenerating twice on the same
baseline now produces byte-identical fixtures.

The optional `meta_extrem` density branch, `MBWR_stableDensity`, full caloric,
coexistence/flash and CSP integration remain pending within the full migration.

```
python tools/rust_port/validate_mbwr_initialization.py --compiler-directory /path/to/compiler/bin
cargo test --test mbwr_source_initialization
```
