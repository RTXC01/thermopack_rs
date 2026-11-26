# Lennard-Jones splined models (`lj_splined` and `hardsphere_wca`)

The splined Lennard-Jones family in ThermoPack (`src/lj_splined.f90` and
`src/hardsphere_wca.f90`) provides pure-fluid perturbation equations of state
for the splined Lennard-Jones (and full Lennard-Jones UF) fluid:

* **`LJS-BH`** (`ljs::LjsBh`, `PreparedLjsBh`): Barker-Henderson perturbation
  theory with Wilhelmsen spline polynomials (`a1`, `a2` with `chi` correction
  and Carnahan-Starling `Khs`, and `a3` spline or Lafitte form).
* **`LJS-WCA`** (`ljs_wca::LjsWca`, `PreparedLjsWca`): Weeks-Chandler-Andersen
  (WCA) perturbation theory with WCA effective hard-sphere diameter, pure
  Carnahan-Starling (`PURE_HS_CS`) or Carnahan-Starling-Kolafa (`PURE_HS_CSK`)
  hard-sphere reference, 61-point Gauss-Kronrod cavity soft-repulsion integral,
  and four-order WCA perturbation polynomials (`a1`–`a4`).
* **`LJS-UV`** (`ljs::LjsUv`, `PreparedLjsUv`): van Westen UV perturbation
  theory combining the WCA hard-sphere reference and cavity soft repulsion with
  the UV attractive `a1` term and low-density second-virial blending
  (`uv_phi` or temperature-dependent `uv_phi_temp`).
* **`LJS-UF` and `LJ-UF`** (`ljs::LjsUf`, `PreparedLjsUf`,
  `LjUfPotential::{LennardJonesSplined, LennardJones}`): UF perturbation theory
  combining the WCA hard-sphere reference and cavity soft repulsion with the
  u-expansion / Mayer-f-expansion blended attractive term (`calc_uf_WCA`).

> Note on terminology: earlier notes referred to the WCA perturbation family as
> "White Bear". In ThermoPack `src/lj_splined.f90` and `src/hardsphere_wca.f90`,
> the perturbation models are **Weeks-Chandler-Andersen (WCA)** models using a
> Carnahan-Starling (`PURE_HS_CS`) or Carnahan-Starling-Kolafa (`PURE_HS_CSK`)
> pure hard-sphere reference. White Bear FMT weighted-density functionals live
> in `src/fmt.f90` (`rust/src/fmt.rs`); spatial DFT convolutions are external
> to ThermoPack.

## Source scope

| Original routine | Rust |
| --- | --- |
| `src/lj_splined.f90:init_LJs_bh`, `ljs_bh_set_sigma_eps`, `ljs_bh_get_pure_params` | `ljs::LjsBh::{new, with_parameters, set_parameters, sigma, eps}` |
| `src/lj_splined.f90:ljs_bh_model_control` | `ljs::LjsBh::{control, set_control, set_kolafa_pure_reference}`, `ljs::LjsControl` |
| `src/lj_splined.f90:calcFresLJs_bh` (all 10 TVN blocks) | `ljs::LjsBh::{residual_density, residual_helmholtz_tvn, residual_helmholtz_energy, residual_helmholtz_energy_derivative, pressure, pressure_tv, fugacity_coefficient, state, tv_state, speed_of_sound_tv}`, `ljs::PreparedLjsBh` |
| `src/lj_splined.f90:calc_LJs_a1`, `calc_LJs_a2`, `calc_LJs_a3` | `ljs::LjsBh::{a1_tvn, a2_tvn, a3_tvn}`, `ljs::PreparedLjsBh::{a1_tvn, a2_tvn, a3_tvn}` |
| `src/lj_splined.f90:calc_ai_LJs_ex`, `calc_a1_LJs_ex`, `calc_a2_LJs_ex`, `calc_a3_LJs_ex`, `calc_ais_LJs_ex`, `calc_chi_LJs_ex`, `calc_ai_reduced_LJs_ex` | `ljs::{ai_ljs_ex, a1_ljs_ex, a2_ljs_ex, a3_ljs_ex, ais_ljs_ex, chi_ljs_ex}`, `ljs::LjsBh::{spline_polynomial, reduced_perturbation_terms}` |
| `src/lj_splined.f90:ljs_bh_get_bh_diameter_div_sigma`, `calc_LJs_bh_zeta` | `ljs::LjsBh::{diameter, diameter_div_sigma, packing_fraction, packing_fraction_tvn}` |
| `src/lj_splined.f90:init_LJs_WCA`, `ljs_wca_model_control`, `ljs_wca_set_pure_params`, `calc_ljx_WCA_zeta` | `ljs_wca::LjsWca::{new, with_parameters, set_parameters, control, set_control, diameter, packing_fraction, packing_fraction_tvn}`, `ljs_wca::LjsWcaControl` |
| `src/lj_splined.f90:calcFresLJs_WCA`, `calc_LJs_WCA_ai`, `calc_LJs_WCA_ai_tr`, `calc_LJs_WCA_ai_TVN` | `ljs_wca::{perturbation_term, wca_term, wca_term_density, wca_term_tvn}`, `ljs_wca::LjsWca`, `ljs_wca::PreparedLjsWca` |
| `src/lj_splined.f90:calcFresLJs_uv_theory`, `ljs_uv_model_control`, `calc_uv_WCA`, `calc_uv_WCA_TVN`, `uv_phi`, `uv_phi_temp`, `uv_a1_u`, `uv_a1_u_mult_T`, `uv_a1_b2`, `uv_q_of_T`, `uv_hs_b2`, `ljs_b2`, `uv_delta_b2_overall` | `ljs::LjsUv`, `ljs::LjsUvControl`, `ljs::PreparedLjsUv`, `ljs::{calc_uv_wca, calc_uv_wca_tvn, uv_phi, uv_phi_temp, uv_a1_u, uv_a1_u_mult_t, uv_a1_b2, uv_q_of_temperature, uv_hard_sphere_b2, ljs_b2, uv_delta_b2_overall}` |
| `src/lj_splined.f90:calcFresLJ_uf_theory`, `calc_uf_WCA`, `calc_uf_WCA_TVN` | `ljs::LjsUf`, `ljs::LjUfPotential`, `ljs::LjsUfControl`, `ljs::PreparedLjsUf`, `ljs::{calc_uf_wca, calc_uf_wca_tvn}` |
| `src/lj_splined.f90:calc_ljs_dispersion`, `calc_ljs_hard_sphere`, `calc_wca_soft_repulsion` | `dispersion_tvn`, `hard_sphere_tvn`, `soft_repulsion_tvn` on `LjsBh`, `LjsWca`, `LjsUv`, `LjsUf` (and prepared counterparts) |
| `src/lj_splined.f90:ljs_potential_reduced`, `alpha_ljs`; `src/saft_interface.f90:{potential, alpha, estimate_critical_parameters, test_fmt_compatibility}` | `ljs::{ALPHA_LJS, ljs_potential_reduced, ljs_pair_potential, estimate_ljs_critical_parameters}` and model methods `reduced_potential`, `pair_potential`, `pair_potential_with_cap`, `alpha_vdw`, `estimate_critical_parameters`, `fmt_compatibility` |
| `src/hardsphere_wca.f90:calc_dhs_WCA`, `minimum_r_cavity_int_LJ`, `cavity_distribution_function`, `lj_boltzmann_factor`, `calcZetaX_vdW_no_segments`, `calc_cavity_integral_LJ_Fres`, `calc_cavity_integral_LJ_Fres_test_TVN` | `hardsphere_wca::{hard_sphere_diameter, minimum_cavity_radius, cavity_distribution_function, lj_boltzmann_factor, packing_fraction_vdw_no_segments, cavity_integral_eta_t, cavity_integral_density, cavity_integral_lj_fres, cavity_integral_lj_fres_tvn, pure_hard_sphere_eta, pure_hard_sphere_density}` |
| `src/saft_interface.f90:pc_saft_nonassoc_volume_solver`, `specificVolume`, `zfac`, `enthalpy`, `entropy`, `thermo`, `sound_velocity_2ph` | `volume`, `volumes`, `baseline_volume_iteration`, `state_tp`, `tv_state`, `speed_of_sound_tv` on all four model types |

## Volume solvers and state functions

All four models (`LjsBh`, `LjsWca`, `LjsUv`, `LjsUf`) share two explicit
volume-solve interfaces:

* **Physical stable solver** (`volume(T, P, phase)` and `volumes(T, P)`): scans
  the physical packing-fraction domain from an unclipped dilute ideal-gas lower
  bracket `eta_lo <= 0.5 * eta_ig` up to `0.75` (extended up to `0.90` only when
  ultra-high pressure exceeds `P(0.75)`, always staying below the
  Carnahan-Starling pole at `eta = 1`), refines both adjacent and sub-cell
  near-critical spinodal pairs (`dP/drho = 0`), decomposes the domain into
  strictly monotone mechanically stable branches, polishes every mechanically
  stable root (`dP/dV < 0`) to pressure-relative tolerance with safeguarded
  Newton-Raphson, and selects `Phase::Liquid` (smallest molar volume),
  `Phase::Vapor` (largest molar volume), or `Phase::MinimumGibbs` (minimum
  `a_res + Z - 1 - ln(Z)`). This fixes the earlier `LjsBh::volume` behavior
  that scanned past the `eta = 1` pole for liquid states.
* **Exact source compatibility iteration** (`baseline_volume_iteration`):
  reproduces `src/saft_interface.f90:pc_saft_nonassoc_volume_solver`
  (`zeta_min = 1e-10`, `zeta_max = 0.75`, vapor seed `1e-10`, liquid seed
  `0.5`, 40-iteration cap) and reports both source convergence flags explicitly.

Full `(T, V, n)` states (`tv_state` -> `LjsTvState`) and `(T, P)` states
(`state_tp` -> `LjsTpState`) use ThermoPack's Argon ideal-gas record with
constant `Cp_id = 2.5 * R`, reporting `Z`, `h`, `s`, `ln(phi)` with analytical
`(T, P)` derivatives, `Cv`, `Cp`, and speed of sound (which returns `NaN` on
mechanically unstable states where `dP/dV >= 0` or `Cp/Cv <= 0`). Legacy
`state()` preserves the positive-pressure error (`Error::InvalidInput("pressure")`)
at nonpositive pressures while `tv_state()` exposes negative-pressure states
explicitly; `LjsBh::fugacity_coefficient` (which returns `F_n = ln(phi) + ln(Z)`)
is preserved unchanged for backward compatibility.

## Validation

### Complete live exact-source suite (`tools/rust_port/validate_ljs_complete.py`)

`tools/rust_port/validate_ljs_complete.py` compiles `tools/rust_port/ljs_probe.f90`
against `.reference/thermopack.dll` (GNU Fortran 14.3.0, `-fdefault-real-8 -fdefault-double-8`,
one fresh process per case) and qualifies **72 cases** (857 scalar comparisons)
committed to `tests/fixtures/ljs_complete.csv` (`tests/ljs_complete.rs`):

* Standalone `calc_cavity_integral_LJ_Fres_test_TVN` (4 cases, including
  `eta > 0.52` quadratic Taylor extrapolation) and binary
  `calcZetaX_vdW_no_segments` (1 case, 25 derivative outputs).
* Standalone UF kernels `calc_uf_WCA` and `calc_uf_WCA_TVN` for `LJS-UF` and
  `LJ-UF` (4 cases).
* `LJS-BH` (17 cases): full 10-block `Fres`, `use_Lafitte_a3`,
  `enable_chi_correction = false`, individual term switches, `PURE_HS_CSK`
  (Kolafa), custom `(sigma, eps)`, `calc_ljs_dispersion`, `calc_ljs_hard_sphere`,
  liquid/vapor/Lafitte TP states (20 outputs each: `v`, `Z`, `H_res`, `H`,
  `S_res`, `S`, `ln(phi)` with `(T, P)` derivatives, and speed of sound), and
  model hooks (`d`, `d_T`, `d_TT`, `alpha`, `Tc`, `vc`, FMT compatibility,
  capped and uncapped pair potential/force).
* `LJS-WCA` (18 cases): full 10-block `Fres`, `eta > 0.52`, `enable_cavity = false`,
  `enable_hs = false`, `PURE_HS_CSK` (Kolafa), individual `a1`–`a4` terms,
  custom `(sigma, eps)`, `calc_ljs_dispersion`, `calc_ljs_hard_sphere`,
  `calc_wca_soft_repulsion`, liquid/vapor/supercritical TP states, and hooks.
* `LJS-UV` (15 cases): full 10-block `Fres`, `eta > 0.52`,
  `use_temperature_dependent_u_fraction = true`, `enable_virial_term = false`,
  `enable_cavity = false`, `PURE_HS_CSK` (Kolafa), custom `(sigma, eps)`,
  all three DFT hooks, liquid/vapor/temp-u TP states, and hooks.
* `LJS-UF` and `LJ-UF` (13 cases): full 10-block `Fres`, `enable_cavity = false`,
  `PURE_HS_CSK` (Kolafa), DFT hooks, and liquid/vapor TP states.

All 72 cases pass within `atol = 1e-12, rtol = 1e-10` on all TVN/DFT/kernel/hook
outputs and `atol = 1e-9, rtol = 1e-8` on all 20-output TP states. In addition,
`tests/ljs_complete.rs` verifies:

* Central finite differences of all 10 TVN derivatives and Euler homogeneity
  across all five models and cavity branches (`eta < 0.52` and `eta > 0.52`).
* Exact molar additivity `F/n = a_hs + a_soft + a_disp` (WCA/UV/UF) and
  `F/n = a_hs + a_disp` (BH).
* Central finite differences of all 10 `(T, P)` derivatives and Maxwell /
  `Cp - Cv` identities across all five models.
* Stable `P(T, v(T, P)) = P` round trips across liquid, vapor, supercritical,
  low-pressure, and `MinimumGibbs` states for all five models, plus exact
  agreement of `LjsBh::volume` and `LjsBh::baseline_volume_iteration` with all
  8 liquid and vapor volumes in `tests/fixtures/ljs_state.csv`.

### Existing LJS / WCA / UV / cavity fixtures (preserved)

All previously qualified fixtures remain active and passing without regression:

* `tests/fixtures/ljs.csv` (8 cases, worst rel `1.5e-12`)
* `tests/fixtures/ljs_spline.csv` (8 cases, worst rel `3.1e-14`)
* `tests/fixtures/ljs_helmholtz.csv` (12 cases, worst rel `4.4e-14`)
* `tests/fixtures/ljs_state.csv` (8 cases, worst rel `1.9e-13`; liquid and vapor
  volumes now also directly matched by `LjsBh::volume`)
* `tests/fixtures/ljs_fugacity.csv` (8 cases, worst rel `1.1e-15`)
* `tests/fixtures/ljs_caloric.csv` (8 cases, worst rel `4.2e-15`)
* `tests/fixtures/ljs_potential.csv` (4 cases, worst rel `2.8e-15`)
* `tests/fixtures/ljs_wca_term.csv` (8 cases, worst rel `0`) and
  `tests/fixtures/wca_term.csv` (5 cases, worst rel `1.0e-15`)
* `tests/fixtures/white_bear.csv` (5 cases, worst rel `1.6e-14`)
* `tests/fixtures/uv_phi.csv` (5 cases, worst rel `4.1e-15`),
  `tests/fixtures/uv_terms.csv` (5 cases, worst rel `2.5e-16`),
  `tests/fixtures/uv_virials.csv` (4 cases, worst rel `3.7e-16`), and
  `tests/fixtures/uv_virial_difference.csv` (5 cases, worst rel `1.5e-15`)
* `tests/fixtures/ais_ljs.csv` (5 cases, worst rel `0`)
* `tests/fixtures/chain_rule.csv` (1 case, worst rel `0`)

## Release timings (`cargo bench --bench ljs`)

Measured on Windows x86_64 MSVC Rust release profile (thin LTO,
`codegen-units = 1`, 5 rounds, median / min / max). Baseline was measured at
commit `cb6efa1a3de0ea731baa1c6f025620f522ade402` on identical workloads and
verified by checksum:

| Workload | Baseline `cb6efa1` median (ns) | New median (ns) | Min / Max (ns) | Checksum |
| --- | ---: | ---: | ---: | --- |
| `bh_helmholtz_deriv` | 1540.5 | 1596.2 | 1561.6 / 1653.2 | `-1.922154e3` |
| `bh_pressure` | 1517.5 | 1527.9 | 1524.8 / 1554.2 | `7.098104e4` |
| `bh_state` | 3179.8 | 1726.3 | 1709.9 / 1799.5 | `2.232673e6` |
| `bh_volume_vapor` | 1,338,514.0 | 14,492.0 | 14,460.0 / 15,095.0 | `1.390804e-1` |
| `wca_term_a2` | 123.5 | 164.2 | 162.8 / 171.7 | `-6.929716e7` |
| `cavity_dist_fn` | 31.6 | 40.8 | 37.3 / 41.5 | `1.734144e6` |
| `bh_prepared_tvn` (new) | — | 127.4 | 124.9 / 148.1 | `-5.924331e11` |
| `bh_volume_liquid` (new) | — | 19,682.2 | 19,540.7 / 20,744.6 | `6.608616e-2` |
| `bh_baseline_volume` (new) | — | 2,485.1 | 2,460.6 / 2,615.7 | `6.608616e-2` |
| `cavity_integral_tvn` (new) | — | 3,582.7 | 3,541.3 / 4,055.0 | `6.133363e10` |
| `wca_tvn_fresh` (new) | — | 4,191.3 | 4,112.1 / 4,517.3 | `-3.562757e11` |
| `wca_tvn_prepared` (new) | — | 1,980.0 | 1,954.2 / 2,197.5 | `-7.125514e11` |
| `wca_volume_liquid_prep` (new) | — | 221,902.0 | 220,373.0 / 225,908.0 | `6.533929e-3` |
| `uv_tvn_prepared` (new) | — | 1,640.1 | 1,637.4 / 1,688.7 | `-6.402049e11` |
| `uv_volume_vapor_prep` (new) | — | 94,940.5 | 94,172.5 / 125,087.5 | `1.691944e-1` |
| `uf_tvn_prepared` (new) | — | 2,258.4 | 1,816.9 / 2,279.1 | `-6.764593e11` |

`bh_state` is 1.84x faster (single `prepare(T)` call), and `bh_volume_vapor` is
92.4x faster while simultaneously fixing the `eta > 1` liquid-root pole defect.
Preparing temperature-only WCA cavity quadrature nodes (`wca_tvn_prepared`)
cuts full 10-block TVN evaluation time by 2.12x compared to fresh calls.

## Known upstream defects and scope boundaries

* **Upstream Fortran defects** (independently confirmed live and documented in
  [baseline-defects.md](baseline-defects.md)):
  * `ljx_ux_eos_constructor` reverses `INDEX` arguments so `init_lj("UF")`
    silently runs UV theory instead of UF theory.
  * `assocIndices_bookkeeping` omits `eosLJS_UF`, causing `init_ljs("UF")` to
    abort with `"At least one CPA-component must self-associate."`.
  * `init_lj("UV")` (`"LJ-UV"`) is absent from `src/eosdata.f90`.
  * `ljs_bh_set_pure_params` leaves `sigma_eff` and `sigma_ij_cube` stale.
  * `calc_cavity_integral_LJ_Fres_test` returns zero for `F_e`, `F_ee`, `F_Te`.
  * `lj_boltzmann_factor` leaves `u_r`/`u_rr` uninitialized for `r > 2^(1/6)`.
  * `calc_ljs_dispersion` / `calc_ljs_hard_sphere` assume `F_T`/`F_V`/`F_n` are
    present whenever second mixed/mole derivatives are requested.
* **Genuine scope boundaries**:
  * All `src/lj_splined.f90` models are pure-fluid (`nc = 1`) equations of state
    in ThermoPack (`packing_fraction_vdw_no_segments` supports general `nc >= 1`
    as in `src/hardsphere_wca.f90`).
  * Spatial DFT convolutions (integrating FMT weighted densities and dispersion
    kernels over spatial grids) do not exist in `src/lj_splined.f90`,
    `src/hardsphere_wca.f90`, or ThermoPack's Fortran `src/` tree; all
    model-side DFT/FMT hooks (`dispersion_tvn`, `hard_sphere_tvn`,
    `soft_repulsion_tvn`, `diameter`, `alpha_vdw`, `pair_potential`,
    `fmt_compatibility`, `estimate_critical_parameters`) are implemented and
    validated.
  * Generic saturation/envelope/flash solvers for LJS models are not part of
    `src/lj_splined.f90` or `src/hardsphere_wca.f90` and remain with their
    respective equilibrium milestones.
