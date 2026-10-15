# Native cubic and flash benchmarks

Measured locally on 2026-10-15. See [benchmark.json](benchmark.json) for
machine/compiler details, reference hashes and individual samples. These
numbers qualify one workload, not the complete library.

| Native operation | Fortran ns/state | Rust ns/state | Speedup |
|---|---:|---:|---:|
| Pressure + T/V derivatives, including mixing | 933 | 395 | 2.36x |
| Vapor volume, including mixing | 1059 | 588 | 1.80x |
| Pressure + T/V derivatives, prepared state | 933 | 8.8 | 106x |
| Vapor volume, prepared state | 1059 | 204 | 5.20x |
| TP flash, including stability and split solving | 102289 | 39696 | 2.58x |
| PH flash, including temperature recovery | 482025 | 357215 | 1.35x |
| PS flash, including temperature recovery | 2333910 | 348610 | 6.69x |

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
for H/S/Cp. This reduced the preliminary local PH/PS times from approximately
864/871 microseconds to the reported 357/349 microseconds, with the same
tolerances and reference fixtures. Reusing preparation does not reduce precision.

Reproduce with:

```
python tools/rust_port/benchmark.py --compiler /path/to/gfortran --repeats 5
```

Next optimization targets are reusable preparation buffers, separately
requested derivative paths and profiling the flash/SAFT hot loops once those
implementations pass parity. No SIMD or thread-level speedups are claimed yet.
