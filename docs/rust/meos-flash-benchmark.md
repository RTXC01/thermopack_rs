# Pure MEOS flash performance qualification

`benchmark_meos_flash.py` compares original Fortran and native Rust loops on a
shared source-qualified manifest. The six profiles are NIST propane/R134A,
full LJ krypton/argon, LJTS krypton and reduced LJTS. Six separate lanes cover
TP, cold PH, guessed PH, cold PS, guessed PS, and prepared scalar TV H/S.
These measurements do not cover all fourteen models or the whole library.
UV/SV and response/sound are separate, currently unmeasured capabilities.

## Run and measurement scope

Build the exact reference library and probe as described in the migration
README. Pass the source build directory when sharing an existing local build:

```text
python tools/rust_port/benchmark_meos_flash.py --compiler /path/to/gfortran --source-build-dir /path/to/reference --preflight-only
python tools/rust_port/benchmark_meos_flash.py --compiler /path/to/gfortran --source-build-dir /path/to/reference --measure --repeats 5 --output docs/rust/meos-flash-baseline.json
python -m unittest discover -s tools/rust_port -p test_benchmark_meos_flash.py
```

Default invocation performs preflight without collecting times. `--measure`
is explicit and must run on the integrated tree with other builds/qualification
jobs idle. Freeze all Rust, Cargo, harness, source object, metadata and input
files throughout measurement. The driver checks these fingerprints after
preflight and again after timing. A report path must be fresh; historical
reports are never overwritten by the driver.

Both workers initialize each model once per process. Parsing, construction,
Python/FFI and process startup are excluded from their timers. Each group
cycles 1,000 fresh calls through the same ordered rows. One warm-up precedes
five measured samples by default, and source/native process order alternates.
Raw samples and all eight aggregate output checksums are retained. Expected
checksums use independently preflighted outputs and exact per-row visit
counts, including incomplete final cycles. Summation tolerances account for
absolute term magnitudes so cancellation does not hide a changed workload.

TP/PH/PS time solving plus equilibrium H/S evaluation. Cold calls use source
T=0 and Rust `None`; guessed calls share identical explicit guesses. The
source sets the same PH/PS temperature bounds per operation. Its model
activation is included. The scalar TV lane evaluates H/S at both fixture
phase volumes, making four scalar calls. Rust temperature preparation is
outside its timer. The source has no corresponding prepared-state object;
its direct scalar calls are timed. That lane does not set flash bounds and
is not the public `equilibrium_properties` API.

Reduced source results use the original physical neon LJTS model and explicit
T/P/V/H/S scales. The driver verifies the committed neon epsilon/k, sigma and
R against the worker's conversion constants. Source PH/PS status -1 is accepted
only after finite, physical and specification checks, matching the public
wrappers; statuses are retained. The reviewed cohort has no warning statuses.

## Correctness and retained scope

The reviewed preflight includes **1,484 rows**, all eligible for timing:
742 flash requests and 742 scalar property requests. The flash cohort includes
**220 endpoint-seeded caloric requests** and **180 near-critical caloric rows**.
There are no per-input mismatch exclusions. All source/native outputs match
within the explicit per-quantity tolerances, and **2,968 independent checks**
verify their actual T/V outputs against the original TV kernels: H/S,
mechanical/thermal stability, pressure and coexistence chemical potential.
Original kernel states, complete raw preflight vectors and statuses are retained.

Endpoint membership is derived from original request metadata (specified
quality 0 or 1), rather than equal volume slots: ordinary single-phase
returns also have equal slots. At an endpoint, a tiny source/native quality or phase-flag
difference must not compare the occupied phase against the other phase's
volume slot. The preflight compares weighted bulk volume and canonical phase
fields only for endpoint inputs, requires the raw split to be within 2e-6 of
an endpoint, and retains every original pair. Material splits are never
collapsed. One-phase quality is undefined and is compared as canonical zero.
Pressure must preserve the input, source caloric specifications close within
2e-7 in R T or R, and native specifications within 3e-10. Every raw result is
checked before canonicalization.

Twenty otherwise-selected pure coexistence TP requests are excluded
because T/P alone does not specify quality. Thirty-four source-failure recovery
rows are outside the primary successful fixture. The separate 36-case endpoint
coordinate-gap audit is fingerprinted; its gap-only rows are not in the primary
fixture. Forty coarse grid-cell intersections with that audit are counted,
using P/(R Tc rhoc), and are not filtered. This is a common successful fluid
branch workload, not a global-stability or complete trajectory qualification.

The report includes all numeric rows and raw outputs, compiler/platform/flags,
Git head, every Rust/Cargo/harness/input fingerprint, actual source and native
executables, DLL/probe identities, and complete source object/module paths and hashes.
Source object names retain their full paths to prevent basename collisions.

## Recorded Windows GNU baseline

The [immutable baseline record](meos-flash-baseline.json) measures the integrated
`5ae8d7a` tree on this Windows GNU host, with one warm-up and five measured
samples. Its checksums and complete frozen-input checks passed. The standalone
report audit verifies recorded cohorts, visit counts, raw sums, medians and
provenance structure without rebuilding Fortran:

```text
python tools/rust_port/compare_meos_flash_benchmarks.py
python -m unittest discover -s tools/rust_port -p test_compare_meos_flash_benchmarks.py
```

The following ratios are **original median / Rust median**; values above one
mean lower observed Rust duration. They describe only this common successful
cohort, these outputs and this host.

| Model | TP | PH cold | PH guessed | PS cold | PS guessed | Scalar TV H/S |
|---|---:|---:|---:|---:|---:|---:|

| C3 | 1.45 | 0.97 | 0.74 | 0.98 | 0.58 | 0.93 |
| R134A | 2.63 | 1.62 | 1.34 | 1.78 | 1.35 | 0.98 |
| Full LJ KR | 5.07 | 1.45 | 1.48 | 1.49 | 1.50 | 1.61 |
| Full LJ AR | 5.26 | 1.48 | 1.53 | 1.56 | 1.54 | 1.62 |
| LJTS KR | 5.10 | 2.91 | 2.86 | 2.83 | 2.75 | 1.67 |
| Reduced LJTS | 5.16 | 2.86 | 2.80 | 2.88 | 2.83 | 1.67 |

Rust has lower observed medians in 30 of 36 groups. TP is 1.45..5.26 times
faster here, and most caloric groups improve, but propane guessed PH/PS regress:
216.08 versus 159.39 microseconds for PH, and 280.53 versus 162.98 microseconds
for PS (about 36% and 72% longer). Those are optimization targets, not omitted
results. The small cold-propane and R134A scalar differences lie within the
observed run-to-run spread; they are not strong evidence of a repeatable win
or loss. Some scalar runs have substantial variability, and all raw samples,
including the warm-up, remain available. No whole-library speedup is inferred.

Within each lane, original and Rust workers use identical inputs. Different
lanes can have different successful-case subsets, so cold/guessed medians are
not a controlled experiment of guess quality by themselves. The record retains
the exact per-group counts and original source request context.
