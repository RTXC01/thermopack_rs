# Residual TP multiparameter ideal-mixture kernels

`MultiparameterIdealMixture` ports the internal rule in
`src/multiparameter_idealmix.f90`. It supports C3, R134A and normal/ortho/para
hydrogen, including repeated constituents and pure limits. These are
composition-weighted pure-fluid states at **common temperature and pressure**;
they are not a common-volume Helmholtz mixture or an interaction model.

```rust
use thermopack_rs::{MultiparameterIdealMixture, Phase};
let model = MultiparameterIdealMixture::new(&["C3", "R134A"])?;
let q = model.residual_tp(450.0, 1e5, &[0.3, 0.7], Phase::Vapor)?;
assert!(q.compressibility_factor > 0.0);
# Ok::<(), thermopack_rs::Error>(())
```

Inputs are SI and mole fractions must be finite, nonnegative and sum to one
within 1e-12. Every constituent is evaluated, including zero-weight entries,
as in the source. Only liquid/vapor seeded homogeneous branches are supported.
Each root must close pressure and have negative P_V. The original comments
restrict the rule to homogeneous regions; global stability is not established.

The result contains Z, residual molar H/S/G, their T/P derivatives, pure
component volumes and fugacity information. Linear composition coefficients
describe `sum(x_i q_i)`. They are **not** constrained mole-number derivatives;
a simplex perturbation adding dx to i and subtracting dx from j has derivative
`q_i - q_j`. Mixture entropy of mixing and total caloric properties are outside
this API. One density solve and one residual jet per constituent supply the
shared properties and derivatives; performance is not measured yet.

`ln_fugacity_coefficients` explicitly reports natural logs. Representable
exponentials are separately available as `fugacity_coefficients`; the source
mixture routine itself returns phi. Its matrix of composition derivatives of
ln(phi) is zero. Physical `residual_gibbs_energy_dp` is
`sum(x_i v_i) - R T/P`. The separately named `baseline_gibbs_dp` reproduces
the source bug that returns Gres in this derivative slot.

## Exact-source qualification boundary

This fork's public `NIST_MEOS_MIX` constructor terminates because its allocator
forwards that unsupported label to the pure-fluid allocation helper. The
qualification program records that failure first, then separately constructs
a genuine `meos_idealmix` instance and allocates its constituents with the
supported `NIST_MEOS` label. It calls the unchanged original routines, linked
from the exact reference objects. No original sources or DLL are modified.
This qualifies the internal residual kernels, not the broken public wrapper.

The audit includes **734 requests**, all matching, and **1,784 independent
original TV-kernel checks** of native constituent volumes and residual H/S/G.
Five pure limits, all ten binary pairs, ternary and five-component cases,
zero fractions, both requested phases, common supercritical points, cryogenic
hydrogen and a compressed C3/R134A state are included. Source roots must have
positive Cv, negative P_V and pressure closure. All requests and raw outputs
are retained; no source failures are filtered out of the kernel cohort.
The 600 K R134A cases exceed its stored 455 K fit limit; they check numerical
extrapolation parity and do not validate the empirical fit's accuracy there.

Tolerance is 2e-8 relative plus explicit per-coordinate absolute floors in
the fixture/metadata. Analytic physical T/P derivatives and constrained
composition changes are also checked by independent finite differences.
The audit retains object/module hashes, original source hashes, compiler
identity, helper and native executable hashes, request schedule, native input
hashes, all original states and the public constructor failure.

```text
python tools/rust_port/qualify_meos_ideal_mixture.py --compiler /path/to/gfortran --source-build-dir /path/to/reference
python -m unittest discover -s tools/rust_port -p test_meos_ideal_mixture_audit.py
cargo test --release --test meos_ideal_mixture --test meos_ideal_mixture_source
```

The committed fixture and metadata can be replayed without Fortran. A source
cache may be reused only when its original objects, compiler, helper,
generator and complete request schedule still match. A cached collection
retains the hash of the executable actually used for that collection,
regardless of whether relinking produces a different executable hash.

Total H/S/G, mixture TV kernels, complete mole derivatives, reference setters,
stability and phase solvers, public language interfaces and performance
qualification remain pending. The [source defect ledger](baseline-defects.md)
separates demonstrated bugs from older observations not reproduced by this
class-safe probe.
