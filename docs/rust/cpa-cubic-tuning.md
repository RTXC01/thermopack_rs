# CPA second derivatives & mixing, UMR/VTPR, HV/WS, apparent composition, tuning, and complex model initializers

This document records the scope, qualification, verified upstream defects, and
release hot-path performance of the native Rust ports covering:

- Complete **CPA** second `(T, V, n)` Helmholtz derivatives, TP state
  derivatives, and cubic mixing (`VdW`, `HuronVidal1`, `HuronVidal2`,
  `WongSandler` / `cbMixWSCPA`, `HvWongSandler`),
- **UMR / VTPR** and **HV / Wong-Sandler** (`cbMixHVWongSandler`,
  `fidel_alpha`) cubic mixing rules and power-mean covolume mixing
  (`b_exponent` / `lowcase_bij` / `calc_bmix_cubic`),
- **`src/apparent_composition.f90`** (`rust/src/apparent_composition.rs`),
- **`src/tuning.f90`** (`rust/src/tuning.rs`, plus `Unifac` group-interaction
  getters/setters and critical-parameter redefinition),
- **`src/complexmodelinit.f90`** (`rust/src/complex_model.rs`: `init_vtpr`,
  `init_umr`).

## Architecture & reuse

- `rust/src/cpa.rs` composes `Cubic` with the shared Wertheim TPT1 association
  engine in `rust/src/association.rs` (`AssociationParams::for_cpa`,
  `AssociationRdfProvider`), replacing all prior finite-difference derivatives in
  `PreparedCpa` (`pressure_and_derivatives`, `fugacity_coefficients`, `state`,
  `volume`, `saturation_pressure`) with exact analytical derivatives and fixing
  the prior `PreparedCpa::state` `s_total` double-counting of `R ln Z`.
- `Cubic::residual_helmholtz_tv` and `PreparedState::residual_helmholtz_tv_unit`
  provide the complete 10-block reduced residual Helmholtz energy (`F..F_nn`) in
  unit-mole coordinates with exact Euler scaling across all cubic mixing rules,
  alpha correlations, quantum beta, power-mean `b_exponent`, and volume
  translation.
- `wong_sandler::hv_wong_sandler_mix` and `hv_wong_sandler_ab_values` implement
  the exact analytic equivalent of Fortran `calc_Q_HVNRTL` without hyperdual
  overhead, and fix pre-existing asymmetric `(tau_ij, alpha_ij)` transpose
  indexing in `rust/src/excess_gibbs.rs` (`hv_nrtl_kernel` and `hv_nrtl_value`).
- `apparent_composition::ApparentComposition` ports `apparent_constructor`,
  `apparent_to_real_mole_numbers`, `real_to_apparent_diff`,
  `real_to_apparent_differentials`, `TP_lnfug_apparent`, and `getModFugacity`
  while guarding zero stoichiometry entries so absent salts do not produce `NaN`
  for present salts.

## Live exact-source parity & independent qualification

Pinned live Fortran parity is implemented in
`tools/rust_port/cpa_cubic_tuning_probe.f90` and
`tools/rust_port/validate_cpa_cubic_tuning.py`, integrated into both the
`native` (`--self-test` under `python` and `python -O`) and `source-parity`
jobs of `.github/workflows/rust.yml`.

- **Committed fixture**: `tests/fixtures/cpa_cubic_tuning.csv` (verified by
  `tools/rust_port/verify_fixtures.py` and `tests/fixtures/cpa_cubic_tuning.metadata.json`).
- **57 live exact-source cases (1,246 scalar comparisons)** across 8 categories:
  - `CPA_TVN` (19 cases): pure `H2O` (4C, 2B `SINTEF`, original CS RDF), `MEOH`,
    `ETOH`, `NH3` (`3B-SINTEF`, `4B-SINTEF`), binary `H2O+CO2`/`CO2+H2O`,
    `H2O+ETOH`, `H2O+MEOH` (`Standard` and `Elliot`), `H2O+C1`, `H2O+C3`,
    ternary `H2O+MEOH+C1`, and CPA with `HV1`, `HV2`, and `WSCPA` (`1e-10` rel /
    `1e-11` abs on all 10 Helmholtz blocks and `X_k`).
  - `CPA_TP` (10 cases): liquid and vapor `v`, `Z`, `H_res`, `H`, `S_res`, `S`,
    `ln(phi)` with all analytical `(T, P, n)` derivatives and speed of sound,
    including `WSCPA` (`1e-8` rel / `1e-9` abs).
  - `CUBIC_TVN` (8 cases) & `CUBIC_TP` (4 cases): `UMR`, `VTPR`, `init_VTPR`,
    and tuned binary/ternary `HVWongSandler` (`1e-10` TVN, `1e-8` TP).
  - `FIDEL` (4 cases): `fidel_alpha` (`1e-12` rel).
  - `BMIX` (2 cases): `lowcase_bij` and `calc_bmix_cubic` at `s = 4/3` and
    `s = 2.0` (`1e-12` rel).
  - `APPARENT` (3 cases): `H2O+NaCl`, `H2O+NaCl+CaCl2` (common ion `Cl-`), and
    identity mapping across all mole, Helmholtz-differential, fugacity, and
    `mod_fugacity` outputs (`1e-10` rel).
  - `TUNING` (7 cases): `HV2`, `WS`, `TWU`/`MC`/`QUANTUM` beta, alpha names
    (`SW`/`PR78`), asymmetric `kij` + `lij`, quadratic volume shift, and active
    UNIFAC group-interaction getter (`1e-10` rel).
- **7 live upstream Fortran defect probes** verified directly against
  `.reference/thermopack.dll` (`init_umr`, `cpa_hv1_init`, `cpa_hv2_init`,
  `cpa_setkijandji`, `cpa_setlijandji`, `cpa_setkij_single`,
  `vtpr_bmix_defect`, `apparent_nan`, `scheme_1ea`). All 8 verified upstream
  defects are documented in [baseline-defects.md](baseline-defects.md).
- **Independent Rust verification (`tests/cpa_cubic_tuning.rs`)**:
  - Central finite-difference checks of all 10 TVN Helmholtz blocks and Maxwell
    symmetry across CPA (`VdW`, original RDF, `Elliot`, `HV2`, `WS`, `CPA-PR`
    `HvWongSandler`), `init_vtpr` (`s = 4/3`), `init_umr` (`s = 2.0`), and
    `QuantumCubic` + `b_exponent = 2.0`, plus `fugacity_derivatives_into`
    cross-checks against full TP states.
  - Euler homogeneity (`F = sum_i n_i F_n[i] + V F_V`), scaling invariance from
    `1e-140` to `1e12` mol, and exact zero-mole infinite-dilution `X_k` / `F_n`
    limits.
  - Multi-salt absent-ion `ApparentComposition` finite-difference and non-`NaN`
    checks, `set_unifac_group_interaction`, and `redefine_critical_parameters`.

## Hot-path release performance (`cargo bench --bench cpa`)

Measured with `cargo bench --bench cpa` (`--release` profile, Windows x86_64).
Adding an allocation-free scalar cubic pressure jet
(`PreparedState::residual_helmholtz_pressure_jet_unit`) and reusing Newton
iterates in `PreparedCpa::volume` reduced `PreparedCpa::pressure_and_derivatives`
latency by 1.54x (`4.59 us -> 2.98 us`) and `Cpa::state_tp` latency by 1.75x
(`88.2 us -> 50.3 us`) with bit-identical checksums:

| Benchmark | Median | Min | Max |
|---|---:|---:|---:|
| `cpa_h2o_pressure_and_derivs` | 2.98 us | 2.87 us | 3.60 us |
| `cpa_h2o_state` | 5.22 us | 4.91 us | 6.87 us |
| `cpa_bin_tvn_prepared` (10-block TVN) | 6.99 us | 6.85 us | 7.10 us |
| `cpa_bin_state_tp` (root + full TP derivs) | 50.29 us | 45.46 us | 63.82 us |
| `cpa_hv2_tvn_prepared` | 6.88 us | 6.73 us | 7.03 us |
| `cubic_hvws_prepare_tvn` (analytic Q) | 2.23 us | 1.42 us | 2.27 us |
| `cubic_vtpr_prepare_tvn` (`s = 4/3`) | 3.87 us | 3.69 us | 4.32 us |
| `cubic_umr_prepare_tvn` (`s = 2.0`) | 3.74 us | 3.50 us | 3.85 us |
| `apparent_tp_lnfug_3x5` | 0.28 us | 0.26 us | 0.29 us |

## Remaining gaps outside this PR

- Generic mixture flashes, phase envelopes, spinodals, and critical-point
  solvers for CPA remain with the equilibrium/flash roadmap entries.
- Electrolyte EOS models consuming `ApparentComposition` remain with their own
  roadmap entry.
