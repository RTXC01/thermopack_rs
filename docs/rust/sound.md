# Fluid sound speed

`Cubic::equilibrium_speed_of_sound(&flash)` computes acoustic speed in m/s
under thermal, mechanical and chemical phase equilibrium. It solves the full
SV matrix in phase amounts, log temperature and log pressure. This avoids the
singular independent-TP response of pure coexistence. It uses analytic volume,
caloric and log-fugacity derivatives, with a native symmetric-indefinite solve.
The supplied state is assumed to be converged; sound speed does not establish
global stability.

`baseline_equilibrium_speed_of_sound` preserves the fork's reported H2/vdW Cp
and its fluid fallback rules. A failed or nonnegative two-phase SV response
uses feed-liquid single-phase sound. Invalid single-phase isentropic stiffness
returns zero. The physical API uses actual caloric derivatives and returns an
error for an invalid stiffness or failed SV response. Both have `..._with_phases`
variants accepting roots in `[liquid_array_root, vapor_array_root]` order;
single-phase flags ignore overrides. Supplied liquid/liquid or vapor/vapor
compositions still need separate equilibrium/stability qualification.

`Cubic::speed_of_sound_tv(T,V,n)` takes total volume in m3 and mole numbers;
`baseline_speed_of_sound_tv` retains source Cv and zero-return conventions.
Prepared-state methods take molar volume. The TV path evaluates Cv directly
without invoking positive-pressure TP residual properties. Negative pressure
is therefore allowed for a valid EOS volume; positive acoustic stiffness alone
does not prove a stable phase. Physical single-phase stiffness is
`P_V|S = P_V|T - T*P_T^2/Cv`, and `c^2 = -v^2*P_V|S / molar_mass_kg`.
The TP and TV paths share untranslated Cv/beta thermal terms and include
temperature-dependent volume translation. Mole-number/volume scaling leaves
sound speed unchanged.

The reference grid attempts 2088 calls: all 817 finite equilibrium and 1258 TV
results match. Thirteen source TP terminations are recorded separately.
The grid covers PR/SRK/RK/vdW, default/tcPR/QuantumCubic constructors, custom
quantum beta and polynomial translation, nonzero kij/lij, pure qualities
0/0.01/0.5/0.99/1, zero feeds, fluid-root overrides and extensive scaling.
Twenty-nine TV states have nonpositive pressure; ten return positive sound.
Absolute/relative tolerances are 1e-5 m/s and 2e-7; phase flags are exact.
Independent constant-entropy pressure perturbations check physical sound for
mixtures and pure coexistence, and agree with TP-response identities where
those are defined. Solid acoustics, other EOS families and exhaustive critical
or failed-source trajectories remain pending.

## Native-loop measurements

Local medians of five measured samples after a discarded warm-up, with matching
source/Rust checksums in every repetition. Both loops execute native code;
Python/ctypes and printing are outside timers. See [default report](benchmark-sound.json)
and [tcPR report](benchmark-sound-tcpr.json) for compiler, source fingerprints,
raw samples and checksums. Neither build uses fast-math or reduced precision.

| Operation | Fortran ns/state | Rust ns/state | Speedup |
|---|---:|---:|---:|
| Default PR TV, including preparation | 2650 | 903 | 2.93x |
| Default PR TV, prepared T/composition | 2650 | 414 | 6.40x |
| Default PR equilibrium, precomputed split | 20253 | 4443 | 4.56x |
| Full tcPR TV, including preparation | 3782 | 2049 | 1.85x |
| Full tcPR TV, prepared T/composition | 3782 | 419 | 9.03x |
| Full tcPR equilibrium, precomputed split | 25645 | 7594 | 3.38x |

TV uses CO2/C1/N2 at 0.7/0.2/0.1 and 300 K, scanning 1000 molar volumes from
1e-4 to 1.099e-3 m3/mol, 500000 evaluations per sample. The prepared path
amortizes temperature/composition mixing; the preparation-inclusive row is the
fairer comparison when those change. Equilibrium uses C1/C2/C3 at 0.5/0.3/0.2,
1 MPa and 100 precomputed TP splits from 180 to 259.2 K, 10000 sound calls per
sample. Phase preparation, volume roots, derivatives and, for two-phase states,
the SV solve remain inside both timers; TP flash initialization is excluded from both. These are
documented workloads, not a speed guarantee for every state or EOS.

Reproduce with `python tools/rust_port/benchmark_sound.py --compiler /path/to/gfortran`
and add `--tcpr` for full tcPR. The matching Rust entry point is
`cargo bench --bench sound`, with `-- --tcpr` for that constructor.
