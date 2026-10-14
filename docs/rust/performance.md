# First cubic benchmark

Measured locally on 2026-10-14. See [benchmark.json](benchmark.json) for
machine/compiler details, reference hashes and individual samples. These
numbers qualify one workload, not the complete library.

| Native operation | Fortran ns/state | Rust ns/state | Speedup |
|---|---:|---:|---:|
| Pressure + T/V derivatives, including mixing | 915 | 384 | 2.38x |
| Vapor volume, including mixing | 1046 | 592 | 1.77x |
| Pressure + T/V derivatives, prepared state | 915 | 8.4 | 109x |
| Vapor volume, prepared state | 1046 | 197 | 5.31x |

PR, CO2/methane/nitrogen at mole fractions 0.7/0.2/0.1, T=300 K.
Pressure scans use 1000 molar volumes from 1e-4 to 1.099e-3 m3/mol;
volume scans use 1000 pressures from 1 MPa to 1.999 MPa. Each sample performs
one million states. Reported values are medians of five samples after a
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
them during calculations. Reusing preparation does not reduce precision.

Reproduce with:

```
python tools/rust_port/benchmark.py --compiler /path/to/gfortran --repeats 5
```

Next optimization targets are reusable preparation buffers, separately
requested derivative paths and profiling the flash/SAFT hot loops once those
implementations pass parity. No SIMD or thread-level speedups are claimed yet.
