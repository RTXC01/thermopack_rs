# MBWR density measurements

The optimized port is 1.48–2.59x faster than the initial Rust implementation on
the measured methane workload. Some MBWR19 cases still trail Fortran. Both the
initial and optimized reports are retained; these measurements do not establish
a whole-library speed claim.

The [initial report](benchmark-mbwr-density.json) compares native executable
loops, excluding Python/FFI and model initialization. Both implementations
return density, phase flag, pressure and pressure slope. There are 1200
per-state preflight comparisons per implementation. Every state closes pressure
and has nonnegative P_rho; raw values and exact phase flags agree. Timed
checksums include all four outputs.

Workloads are C1 with MBWR19/32, liquid/vapor seeds at 150 K and a vapor seed
at 300 K, 100 pressures from 0.5 to 0.995 MPa, and 20000 calls per sample.
Five samples follow one warm-up. Fresh mode includes temperature preparation
per query; prepared mode reuses it. Original model/critical/SRK initialization
is excluded in both; native fresh preparation also validates supplied metadata
and recalculates its SRK constants. The report records source, harness and
Cargo fingerprints and rejects modifications during measurement.

Initial median microseconds per state:

| Family / state | Fortran fresh | Rust fresh | Fortran prepared | Rust prepared |
|---|---:|---:|---:|---:|
| MBWR19 liquid | 1.017 | 1.691 | 0.851 | 1.473 |
| MBWR19 vapor | 0.608 | 0.968 | 0.435 | 0.948 |
| MBWR19 supercritical | 0.525 | 0.779 | 0.342 | 0.695 |
| MBWR32 liquid | 1.504 | 2.518 | 1.176 | 2.317 |
| MBWR32 vapor | 0.895 | 1.541 | 0.701 | 1.375 |
| MBWR32 supercritical | 0.778 | 1.189 | 0.573 | 1.144 |

These are host/workload-specific Windows GNU-toolchain results, not a general
library speed claim. Critical searches, bounded branch discovery, fallback-heavy
states, other fluids and other model families are outside this workload.
The initial port computes pressure's second density derivative even when the
Newton caller needs only pressure and its first derivative; this is an
optimization target to test while retaining the raw source comparisons.

## Optimized implementation

The [optimized report](benchmark-mbwr-density-optimized.json) measures the
checked-in combination of three changes: skip second density derivatives when
the caller needs only P/P_rho, share the original positive integer-power
multiplication trees rather than calling runtime `powi` repeatedly, and compile
separate fixed-size pressure loops for MBWR19/32. Polynomial/Gaussian summation
order, source phase rules and stopping conditions are retained. No fast-math or
precision reduction is used. All 8839 finite raw source outputs and the one
non-finite-source/native-error case continue to match at unchanged tolerances;
native path counts remain unchanged.

Final median microseconds per state:

| Family / state | Fortran fresh | Rust fresh | Fortran prepared | Rust prepared |
|---|---:|---:|---:|---:|
| MBWR19 liquid | 1.012 | 1.051 | 0.897 | 0.993 |
| MBWR19 vapor | 0.620 | 0.514 | 0.450 | 0.470 |
| MBWR19 supercritical | 0.532 | 0.462 | 0.343 | 0.404 |
| MBWR32 liquid | 1.445 | 1.186 | 1.194 | 1.100 |
| MBWR32 vapor | 0.940 | 0.595 | 0.788 | 0.555 |
| MBWR32 supercritical | 0.908 | 0.501 | 0.688 | 0.440 |

The same-run MBWR32 Fortran/Rust ratios are 1.22–1.81 for fresh and 1.09–1.56
for prepared calls. MBWR19 fresh vapor/supercritical calls improve, but its
fresh liquid and all three prepared workloads still trail Fortran. Those remain
performance gaps. The 1.48–2.59x comparison against initial Rust uses medians
from separate runs; host scheduling/clock variation limits fine-grained claims.
In particular, use the raw samples rather than inferring small compiler wins
from rounded table entries.

The [derivative-only experiment](benchmark-mbwr-density-derivatives.json),
[cached-power experiment](benchmark-mbwr-density-powers.json) and
[iterator experiment](benchmark-mbwr-density-iterator.json) retain intermediate
measurements and exact source fingerprints. Skipping derivatives alone regressed
the initial code on this host. Assembly inspection showed runtime integer-power
helper calls; the final power cache removes those calls. The iterator rewrite
of the source-indexed polynomial loop also regressed timings, so the final
fixed-index loop has a narrowly scoped Clippy allowance. These observations
guide this implementation, not a general rule for Rust iterators or compilers.

The final benchmark again passes 1200 per-state comparisons per implementation,
pressure closure/stability, exact phase flags, four-output timed checksums and
all source/harness/Cargo fingerprint guards. CI runs correctness preflights;
timing ratios are not hardware-independent pass/fail gates.

```
python tools/rust_port/benchmark_mbwr_density.py --compiler /path/to/gfortran
python tools/rust_port/benchmark_mbwr_density.py --compiler /path/to/gfortran --preflight-only
```
