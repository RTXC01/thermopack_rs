# Coverage ledger

Status definitions: **pending** = no native implementation; **native** = Rust
implemented, awaiting reference qualification; **validated** = tested against
the exact source baseline for documented cases. Validation of a subset does
not qualify an entire family.

| Capability | Source | Status |
|---|---|---|
| Component/interaction database | fluids/, binaries/, compdatadb.f90, mixdatadb.f90 | native: 109 fluids with compiled numerical constants/default Cp and PR/SRK interactions; JSON aliases; non-default sets pending |
| PR, SRK, RK, van der Waals; classic alpha/vdW mixing | cbmix, cbalpha, cubic, cbhelm | validated: 1320 state cases; RK = SRK with RK alpha in the reference |
| Other cubic, alpha, beta, mixing, translation, pseudo fluids | cubic_eos, cbselect, cbbeta, excess_gibbs, wong_sandler, volume_shift | pending |
| Ideal gas, H2, caloric reference states | ideal, idealh2 | native correlations 1-12; 2592 default-property comparisons across 108 fluids and 720 configured/reference-state cases; 24 BUTANAL source failures recorded; non-default database reference selection pending |
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
