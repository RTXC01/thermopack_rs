# Native cubic and flash benchmarks

Fluid acoustic workloads have separate current [measurements and scope](sound.md#native-loop-measurements).

Local measurements for the recorded source versions. See [benchmark.json](benchmark.json) for
machine/compiler details, reference hashes and individual samples. Full tcPR
construction has its own [report](benchmark-tcpr.json). These
numbers qualify the documented workloads, not the complete library.

| Native operation | Fortran ns/state | Rust ns/state | Speedup |
|---|---:|---:|---:|
| Pressure + T/V derivatives, including mixing | 922 | 422 | 2.18x |
| Vapor volume, including mixing | 1045 | 614 | 1.70x |
| Pressure + T/V derivatives, prepared state | 922 | 13.9 | 66.19x |
| Vapor volume, prepared state | 1045 | 209 | 5.00x |
| TP flash, including stability and split solving | 105167 | 39495 | 2.66x |
| PH flash, including temperature recovery | 480804 | 297131 | 1.62x |
| PS flash, including temperature recovery | 2382249 | 281709 | 8.46x |
| Ideal H/S/Cp, uncached reference constants | 163 | 474 | 0.34x |
| Ideal H/S/Cp, cached reference constants | 163 | 240 | 0.68x |

PR, CO2/methane/nitrogen at mole fractions 0.7/0.2/0.1, T=300 K.
Pressure scans use 1000 molar volumes from 1e-4 to 1.099e-3 m3/mol;
volume scans use 1000 pressures from 1 MPa to 1.999 MPa. Each sample performs
one million states. The TP workload uses PR methane/ethane/propane at
0.5/0.3/0.2, P=1 MPa, with 100 temperatures from 180 to 259.2 K. Each
TP sample performs 10000 complete flashes, including allocation, stability,
phase decisions and split iterations; phase-fraction/phase-code checksums
must agree. PH/PS use H/S targets generated at those same 100 states outside
the timer. Each caloric sample performs 1000 complete calls with default
initial guesses, including temperature recovery, stability, split solving and
allocation. Before timing, both implementations must recover every target
temperature within 1e-5 K; the source must also close H/(RT) or S/R within
1e-8. The fork's PS wrapper permits status -1, so that status is allowed only
when these independent checks pass. Temperature/phase-fraction/phase-code
checksums must agree. Reported values are medians of five samples after a
discarded warm-up. All timed loops are native code. Checksums of pressure,
derivatives and volume must agree before timings are accepted.

The ideal-gas workload uses pure propane's default DIPPR correlation at 1000
temperatures from 150 to 849.3 K, P=0.1 MPa, one million states per sample.
Both paths return H, S and Cp, with their checksums compared to native Fortran.
The cached Rust path is the one used by `Cubic` model calculations. It is
1.98x faster than standalone `Component::ideal`, but the Fortran ideal kernel
is still about 1.47x faster; this remains a measured optimization gap.

Rust uses release thin LTO and one codegen unit. Fortran uses the recorded
gfortran `-O3` flags. Neither uses fast-math or reduced precision. The reference
uses its public numerical interface, which recalculates state parameters;
the prepared Rust path intentionally amortizes temperature/composition mixing.
The rows that include mixing are the fairer comparison for changing states.
The Fortran baseline has not been exhaustively retuned for this workload.

Implemented optimizations: cache temperature/composition parameters, traverse
only nonzero kij pairs after factoring the ideal attractive sum, evaluate
volume/pressure kernels without allocation, calculate roots analytically with
polishing and stable deflation, embed database constants rather than parse
them during calculations. Caloric flashes limit preliminary homogeneous
iterations across a phase-transition jump and reuse ideal/residual evaluations
for H/S/Cp. Owned ideal parameters cache reference integration constants;
DIPPR shares two exponential factors across Cp/H/S. The alpha evaluator is
inlined into parameter preparation to reduce the cost of per-component
correlation dispatch. The first table uses default classic alpha; neither table
qualifies every alpha correlation. Preliminary local PH/PS
times were approximately 864/871 microseconds; the current default report is 297/282
microseconds with the same tolerances and reference fixtures. Reusing
preparation does not reduce precision.

The tcPR workload uses exactly the same state grids, but both implementations
initialize `parameter_reference="tcPR"`. This forces TWU alpha and volume
translation, including their temperature derivatives. Ideal C3 measurements
still use the default correlation in both runs. Initialization is outside the
timed state loops; the rows with mixing include per-call parameter preparation.
The harness validates numeric row names explicitly because tcPR initialization
can print informational lines containing commas. Every required timing row,
finite checksum and caloric preflight must still pass.

| Full tcPR operation | Fortran ns/state | Rust ns/state | Speedup |
|---|---:|---:|---:|
| Pressure + T/V derivatives, including mixing | 1450 | 956 | 1.52x |
| Vapor volume, including mixing | 1492 | 1162 | 1.28x |
| TP flash | 128192 | 56493 | 2.27x |
| PH flash | 621815 | 407865 | 1.52x |
| PS flash | 3055259 | 404745 | 7.55x |

TWU's two powers share `ln(T/Tc)` before exponentiation, saving one logarithm
per component evaluation. Separate multiplicative factors retain the source's
evaluation order. Ratios outside the normal floating-point range use the
original power functions, including their zero-exponent semantics. Alpha and
full translated property/flash reference fixtures pass unchanged. The
[pre-optimization tcPR measurements](benchmark-tcpr-before-shared-log.json)
and [post-optimization measurements](benchmark-tcpr-after-shared-log.json)
retain the source hashes and five raw samples for comparison. These archived
reports isolate that earlier optimization from subsequent model changes.
On that workload, sharing the logarithm reduced the observed median Rust
pressure call including mixing from 1272 to 1026 ns (19%) and the volume call
including mixing from 1450 to 1232 ns (15%). These are local before/after
measurements; the individual samples remain available in both reports.

Translated preparation now borrows active coefficient records directly, avoiding
copies of the public getter's coefficients and provenance metadata for every
component. The [earlier coefficient-copying measurements](benchmark-tcpr-before-borrowed-translation.json)
retain the comparison samples. Inline pressure helpers also allow root-polishing
calls to omit unused temperature-derivative terms. Neither change reduces
precision or changes the public parameter getter. Borrowing the coefficients
reduced the measured translated pressure call including preparation from
1389 to 956 ns (31%) and the volume call from
1573 to 1162 ns (26%).

Reproduce with:

```
python tools/rust_port/benchmark.py --compiler /path/to/gfortran --repeats 5
python tools/rust_port/benchmark.py --compiler /path/to/gfortran --repeats 5 --tcpr
```

Next optimization targets are reusable preparation buffers, separately
requested derivative paths and profiling the flash/SAFT hot loops once those
implementations pass parity. No SIMD or thread-level speedups are claimed yet.
