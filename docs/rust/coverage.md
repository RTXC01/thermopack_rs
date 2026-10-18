# Coverage ledger

Status definitions: **pending** = no native implementation; **native** = Rust
implemented, awaiting reference qualification; **validated** = tested against
the exact source baseline for documented cases. Validation of a subset does
not qualify an entire family.

| Capability | Source | Status |
|---|---|---|
| Extended Rachford–Rice and one-step initialization | tp_solver::rr_solve | validated: 698 finite/status comparisons, six non-finite source calls recorded; 2–20 components and negative-flash roots; full flash integration pending |
| LV RR successive-substitution iteration | tp_solver::rr_successive_substitution_iteration | validated: 2110 finite results/statuses including all four convergence/RR flag combinations; 60 non-finite source calls recorded; complete TP iteration/restart path and optional non-vapor phaseY pending |
| Component/interaction database | fluids/, binaries/, compdatadb.f90, mixdatadb.f90 | native: 109 fluids with compiled numerical constants; all 111 Cp and 191 ci records; default/non-default classical PR/SRK interactions; JSON aliases; other model parameter sets pending |
| PR, SRK, RK, van der Waals; classic alpha/vdW mixing | cbmix, cbalpha, cubic, cbhelm | validated: 1320 state cases; RK = SRK with RK alpha in the reference |
| Cubic alpha correlations and parameters | cbalpha, compdatadb | validated alpha initializer/setters: 4232 numerical cases, 80 source terminations recorded; 3240 EOS/derivative cases; all 257 TWU/MC records embedded; explicit UMR/GERG evaluation returns source-unimplemented errors |
| Volume translation and reference-aware cubic construction | volume_shift, eoslibinit, cbhelm | validated: 1977 parameter/derivative, 3120 full property and 3821 flash/saturation cases; constant/linear/quadratic/quintic ci, Rackett fallback and complete tcPR construction; 139 source-failed flash attempts recorded |
| Cubic covolume beta and QuantumCubic construction | cbbeta, eoslibinit, cbmix, cbhelm | validated: 4036 parameter/critical-data, 4512 property/derivative and 36 pure saturation cases; explicit classic/quantum setters, full constructor and qcubic helper with combined references/model-local bounds; cryogenic flash qualification remains pending |
| Other cubic geometry, mixing, pseudo fluids | cubic_eos, cbselect, excess_gibbs, wong_sandler | pending |
| Ideal gas, H2, caloric reference states | ideal, idealh2 | native correlations 1-12; 2592 default-property comparisons across 108 fluids and 720 configured/reference-state cases; 24 BUTANAL source failures recorded; cubic constructor's unconditional default-Cp fallback preserved |
| T-V / T-P properties and derivatives, chemical potentials, sound | eos, eostv, single_phase, speed_of_sound | validated cubic P, dP/dT, dP/dv, H, S, log fugacity and its T/P/n derivatives, baseline dH/dT; other derivatives pending; physical sound native, source sound compatibility pending |
| Stable TP/PH/PS/UV/SV/mu-T flashes, grid flashes | stability, tp_solver, ph_solver, ps_solver, uv_solver, sv_solver, mut_solver | validated cubic TP: 1043 physical splits with one phase permutation and three reference aborts; PH/PS: 697 closing source results from 912 attempts; 215 reference failures recorded; raw ordering, failed-source compatibility and UV/SV/mu-T/grid flashes pending |
| Pure/mixture saturation and phase envelopes, VLLE | saturation*, puresaturation, binaryplot, multi_phase_envelope_tv | validated pure cubic saturation P, phase volumes, Clapeyron slope in 48 cases and saturation T in 36 cases; fixed-composition pseudo saturation native; mixture saturation/envelopes/VLLE pending |
| Critical points, spinodals, isolines, Joule-Thomson | critical, spinodal, isolines, joule_thompson_inversion | fixed-composition cubic pseudo-critical helper validated in 28 cases; true mixture critical points and remaining curves pending |
| CPA, association, polar interactions | cpa_parameters, saft_association, association_var, multipol | pending |
| PC-SAFT (including polar/simplified variants) | pc_saft* | pending |
| SAFT-VR Mie / VRQ Mie | saftvrmie* | pending |
| Multiparameter and ideal mixtures | multiparameter* | pending |
| MBWR, CSP, Lee-Kesler | mbwr*, extcsp, leekesler | pending |
| PETS, LJ spline BH/WCA/UV | pets, lj_splined, hardsphere_wca | pending |
| Solids, melting/sublimation, hydrates | solideos, *_gibbs, solid*, hydrate* | pending |
| DFT/FMT, quadratures, hyperdual utilities | fmt, saft_rdf, quadratures, hyperdual* | pending |
| Python/C/C++ compatibility | addon/pycThermopack, include/, addon/cppThermopack | pending |

The exhaustive machine-readable surface is in inventory.json. All rows must
be qualified before describing this project as a complete replacement.
