# Native cubic and flash benchmarks

Measured locally on 2026-10-15. See [benchmark.json](benchmark.json) for
machine/compiler details, reference hashes and individual samples. These
numbers qualify the documented workloads, not the complete library.

| Native operation | Fortran ns/state | Rust ns/state | Speedup |
|---|---:|---:|---:|
| Pressure + T/V derivatives, including mixing | 978 | 429 | 2.28x |
| Vapor volume, including mixing | 1101 | 614 | 1.79x |
| Pressure + T/V derivatives, prepared state | 978 | 8.5 | 115x |
| Vapor volume, prepared state | 1101 | 197 | 5.58x |
| TP flash, including stability and split solving | 107113 | 39066 | 2.74x |
| PH flash, including temperature recovery | 512543 | 329913 | 1.55x |
| PS flash, including temperature recovery | 2500192 | 321718 | 7.77x |
| Ideal H/S/Cp, uncached reference constants | 181 | 479 | 0.38x |
| Ideal H/S/Cp, cached reference constants | 181 | 232 | 0.78x |

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
2.07x faster than standalone `Component::ideal`, but the Fortran ideal kernel
is still about 1.28x faster; this remains a measured optimization gap.

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
correlation dispatch. These timings use default classic alpha; they do not
qualify the performance of every new correlation. Preliminary local PH/PS
times were approximately 864/871 microseconds; the current report is 330/322
microseconds with the same tolerances and reference fixtures. Reusing
preparation does not reduce precision.

Reproduce with:

```
python tools/rust_port/benchmark.py --compiler /path/to/gfortran --repeats 5
```

Next optimization targets are reusable preparation buffers, separately
requested derivative paths and profiling the flash/SAFT hot loops once those
implementations pass parity. No SIMD or thread-level speedups are claimed yet.
