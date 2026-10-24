# Pure multiparameter response and sound

`Multiparameter::equilibrium_response(&flash)` returns homogeneous Cp,
dV/dT at fixed pressure, and dV/dP at fixed temperature. Pure coexistence
constrains pressure to the saturation curve, so independent TP response is
undefined: the method validates the supplied pressure and returns `None`,
including at quality endpoints. Singular or nonpositive Cp gives no response.

`equilibrium_speed_of_sound(&flash)` uses the fit's physical molar mass and
includes quality redistribution along an isentrope. Homogeneous sound uses
the Helmholtz derivative identity directly. For coexistence, with saturation
slope p', it evaluates

```text
Cv_eq,i = Cv_i - T (P_T,i - p')² / P_V,i
Cv_eq   = (1-beta) Cv_eq,L + beta Cv_eq,V
w²      = V_mix² T p'² / (M Cv_eq)
```

Here P_V is negative on each mechanically stable phase. These expressions
avoid subtracting large phase and quality derivatives near the critical point.
Endpoints retain the coexistence-side response; their sound can differ from
the homogeneous phase's sound. The phase flag therefore matters.

All methods require a normalized pure feed, positive finite T/P and, for a
split, pure phase compositions and physical fractions summing to one. They
recompute phase volumes rather than trusting cached volumes in `FlashResult`.
Homogeneous volume follows the supplied phase's source density seed, with
liquid then vapor fallback for `Single`. This does not prove global stability.
Physical coexistence sound uses native coexistence at T after requiring the
supplied P to agree within 1e-8 relative.

## Original TP wrapper compatibility

`baseline_equilibrium_speed_of_sound` reproduces the original public TP
wrapper's **component-database** mass, dimensional derivative grouping and
pure SV matrix solve, including its liquid-sound fallback. Coexistence phase
volumes use the original phase-seeded density iteration at the supplied T/P.
Preserving that pressure matters near the critical point.

The fit and database masses differ slightly for propane, R134A and hydrogen.
For example, the selected hydrogen database records use 2.0159 g/mol.
The public TP wrapper converts database g/mol to kg/mol correctly for LJ
models too. A separate, direct internal MEOS TV sound routine has the LJ
constructor mass-unit defect described in [the kernel notes](multiparameter.md);
it is a different compatibility target. Reduced LJTS uses R=M=1 in native
reduced units. Physical models use SI.

## Qualification

The committed `multiparameter_response_source.csv` contains **854 strict
matches** across all fourteen variants: 224 homogeneous phase controls and
630 coexistence cases at nine T/Tc values from 0.45 through 0.9999 and qualities
0, 0.1, 0.5, 0.9 and 1. All 854 source calls returned qualified states; there
are no accepted errors, censored calls or coordinate differences in this grid.
Response columns come from the original public enthalpy/volume wrappers;
compatibility sound comes from the original public TP sound wrapper.
Reduced expectations scale the physical **neon LJTS** wrapper using its
committed R, epsilon/k, sigma and database mass.

Each native output is also checked using original TV kernels at the actual
native phase volumes: pressure, mechanical/thermal stability, coexistence
chemical potential, response identities and fit-mass physical sound. Those
854 checks are retained in the metadata, together with all raw wrapper
results, warnings, requests, scales and source/native fingerprints. Physical
sound expectations use these original kernels rather than native gold values.

Relative tolerances are 2e-7 for Cp/V_T/V_P, 2e-8 for physical sound and
pressure, and 2e-6 for source-compatible TP sound. Absolute floors are
quantity-specific (notably 1e-12 for V_T and 1e-16 for V_P); reduced pressure
floors scale with units. Rust replay needs no Fortran runtime. Independent
isentropic central differences cover all variants, five reduced temperatures
and three interior qualities. Additional tests cover endpoint sound, mass
conventions and malformed inputs.

```text
python tools/rust_port/qualify_multiparameter_response.py --compiler-directory /path/to/compiler/bin
python tools/rust_port/qualify_multiparameter_response.py --reuse-source --compiler-directory /path/to/compiler/bin
python -m unittest discover -s tools/rust_port -p test_meos_response_audit.py
cargo test --test multiparameter_response
```

Source reuse verifies the DLL/probe, worker, generator and numeric-input
fingerprints and the complete request schedule. Every native error,
unqualified source return or comparison failure fails qualification. A source
call has a four-second execution deadline; a timeout is retained as censored,
the worker is killed and reaped, and that case is never silently retried.

This is a bounded pure-fluid qualification. Broad singular/fallback behavior,
mixture response, global stability, alternative reference shifts and performance
measurements remain pending. No speedup claim is made for these methods.
