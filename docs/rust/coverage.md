# Coverage ledger

Status definitions: **pending** = no native implementation; **native** = Rust
implemented, awaiting reference qualification; **validated** = tested against
the exact source baseline for documented cases. Validation of a subset does
not qualify an entire family.

| Capability | Source | Status |
|---|---|---|
| Extended Rachford–Rice and one-step initialization | tp_solver::rr_solve | validated: 698 finite/status comparisons, six non-finite source calls recorded; 2–20 components and negative-flash roots; integrated into baseline TP compatibility |
| LV RR successive-substitution iteration | tp_solver::rr_successive_substitution_iteration | validated: 2110 finite results/statuses including all four convergence/RR flag combinations; 60 non-finite source calls recorded; integrated into baseline TP; full trajectory fidelity and optional non-vapor phaseY pending |
| Symmetric indefinite linear solve | pinned LAPACK DSYTF2, DSYTRS2, DSYCONV | validated: 226 factor/pivot/solution/status cases, 1–32 variables, no source failures; allocation-free with caller buffers; unblocked upper triangle only; large blocked LAPACK paths pending |
| Equilibrium TP response | state_functions::dhdt_twoPhase, dvdt_twoPhase, dvdp_twoPhase | validated baseline Cp/dVdT/dVdP in 435 finite source results; six source terminations recorded; native physical response and U/S derivatives independently checked by finite differences; pure latent TP response is undefined; general critical/fallback trajectories pending |
| Fixed-volume caloric properties and UV/SV flashes | eostv, uv_solver, sv_solver, state_functions | validated cubic subset: 141 TV P/U/S and 409 finite UV/SV outputs, plus 400 default/tcPR benchmark preflight comparisons; 137 failed source calls recorded; pure latent/zero feeds and physical derivatives; non-finite caloric residuals rejected; single-phase output canonicalized explicitly; full raw output/trajectory fidelity, broad quantum/critical and other model families pending |
| Fluid sound speed | speed_of_sound::sound_velocity_2ph, speed_of_sound_TV | validated baseline: 817 equilibrium and 1258 TV results; pure coexistence including endpoints, zero feeds, fluid-root overrides, quantum beta/translation and extensive scaling; 13 source TP failures recorded; physical derivatives checked along constant-entropy paths; solid acoustics and broad critical/multiphase trajectories pending |
| Fork-order TP compatibility | tp_solver::twoPhaseTPflash, mod_newton_search; optimizer | native baseline_tp_flash: raw phase ordering/fractions/compositions match all 1043 classical and 1130 translated snapshots, plus 718/721 quantum snapshots; three ghost-vapor composition differences retained; reuses native stability search and unblocked symmetric solve, so complete source path fidelity remains pending |
| Component/interaction database | fluids/, binaries/, compdatadb.f90, mixdatadb.f90 | native: 109 fluids with compiled numerical constants; all 111 Cp and 191 ci records; default/non-default classical PR/SRK interactions; JSON aliases; other model parameter sets pending |
| PR, SRK, RK, van der Waals; classic alpha/vdW mixing | cbmix, cbalpha, cubic, cbhelm | validated: 1320 state cases; RK = SRK with RK alpha in the reference |
| Cubic alpha correlations and parameters | cbalpha, compdatadb | validated alpha initializer/setters: 4232 numerical cases, 80 source terminations recorded; 3240 EOS/derivative cases; all 257 TWU/MC records embedded; explicit UMR/GERG evaluation returns source-unimplemented errors |
| Volume translation and reference-aware cubic construction | volume_shift, eoslibinit, cbhelm | validated: 1977 parameter/derivative, 3120 full property and 3821 flash/saturation cases; constant/linear/quadratic/quintic ci, Rackett fallback and complete tcPR construction; 139 source-failed flash attempts recorded |
| Cubic covolume beta and QuantumCubic construction | cbbeta, eoslibinit, cbmix, cbhelm | validated: 4036 parameter/critical-data, 4512 property/derivative and 36 pure saturation cases; explicit classic/quantum setters, full constructor and qcubic helper with combined references/model-local bounds; cryogenic flash qualification remains pending |
| Other cubic geometry, mixing, pseudo fluids | cubic_eos, cbselect, excess_gibbs, wong_sandler | pending |
| Ideal gas, H2, caloric reference states | ideal, idealh2 | native correlations 1-12; 2592 default-property comparisons across 108 fluids and 720 configured/reference-state cases; 24 BUTANAL source failures recorded; cubic constructor's unconditional default-Cp fallback preserved |
| T-V / T-P properties and derivatives, chemical potentials, sound | eos, eostv, single_phase, speed_of_sound | validated cubic P, dP/dT, dP/dv, H, S, log fugacity and its T/P/n derivatives, baseline dH/dT; fluid sound qualification above; other derivatives and solid acoustics pending |
| Stable TP/PH/PS/UV/SV/mu-T flashes, grid flashes | stability, tp_solver, ph_solver, ps_solver, uv_solver, sv_solver, mut_solver | validated cubic TP: 1043 physical splits with one phase permutation and three reference aborts; PH/PS: 697 closing source results from 912 attempts; 215 reference failures recorded; UV/SV scope above; full raw ordering, failed-source compatibility and mu-T/grid flashes pending |
| Pure/mixture saturation and phase envelopes, VLLE | saturation*, puresaturation, binaryplot, multi_phase_envelope_tv | validated pure cubic saturation P, phase volumes, Clapeyron slope in 48 cases and saturation T in 36 cases; fixed-composition pseudo saturation native; mixture saturation/envelopes/VLLE pending |
| Critical points, spinodals, isolines, Joule-Thomson | critical, spinodal, isolines, joule_thompson_inversion | fixed-composition cubic pseudo-critical helper validated in 28 cases; true mixture critical points and remaining curves pending |
| CPA, association, polar interactions | cpa_parameters, saft_association, association_var, multipol | pending |
| PC-SAFT (including polar/simplified variants) | pc_saft* | pending |
| SAFT-VR Mie / VRQ Mie | saftvrmie* | pending |
| Multiparameter and ideal mixtures | multiparameter* | pending |
| MBWR, CSP, Lee-Kesler | mbwr*, extcsp, leekesler | pending |
| PeTS Helmholtz and homogeneous TV | pets, eostv, speed_of_sound | validated: 588 states with 56 outputs and 16 vector-potential cases, no failed source calls; full first/second pure derivatives, TV caloric/mechanical/fugacity/sound properties, LJ correction, source-order and analytic virials, default/custom parameters; higher derivatives and other inherited SAFT interfaces pending |
| PeTS pure pressure roots and coexistence | saft_interface, critical, saturation | validated subset: all 294 finite raw last iterates, 237 independently pressure-closing stable roots, three critical points and 36 saturation/inverse states; 57 unclosed/unstable source root returns retained; default/custom parameter scaling, strict/source iteration APIs and analytic critical-density derivatives; broad degeneracy and complete source trajectory/flag fidelity pending |
| PeTS pure TP/PH/PS and equilibrium properties | tp_solver, ph_solver, ps_solver, state_functions | validated subset: 923 finite/closing source states from 1035 attempts (126 TP, 395 PH, 402 PS); 112 failed/nonclosing source calls recorded; explicit endpoint canonicalization with all raw source phase flags/arrays preserved; native temperature/specification and latent-quality closure; broad degeneracy and raw phase-flag/trajectory fidelity pending |
| PeTS equilibrium TV and UV/SV | uv_solver, sv_solver, state_functions | validated subset: 399 finite/closing source states from 666 attempts (224 UV, 175 SV); 260 source terminations and seven unqualified returns retained; pure scalar isochore solve with analytic equilibrium Cv, homogeneous/latent U/S/V closure, endpoint canonicalization and all raw source flags retained; broad degeneracy and raw phase-flag/trajectory fidelity pending |
| LJ spline BH/WCA/UV | lj_splined, hardsphere_wca | pending |
| Solids, melting/sublimation, hydrates | solideos, *_gibbs, solid*, hydrate* | pending |
| DFT/FMT, quadratures, hyperdual utilities | fmt, saft_rdf, quadratures, hyperdual* | pending |
| Python/C/C++ compatibility | addon/pycThermopack, include/, addon/cppThermopack | pending |

The exhaustive machine-readable surface is in inventory.json. All rows must
be qualified before describing this project as a complete replacement.
