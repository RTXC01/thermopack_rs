# MBWR density measurements

The initial faithful density port is slower than the original Fortran on the
measured methane workload. This baseline is retained before optimization;
switching languages alone does not establish a speed improvement.

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

```
python tools/rust_port/benchmark_mbwr_density.py --compiler /path/to/gfortran
python tools/rust_port/benchmark_mbwr_density.py --compiler /path/to/gfortran --preflight-only
```
