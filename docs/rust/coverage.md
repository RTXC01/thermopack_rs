# Coverage ledger

Status definitions: **pending** = no native implementation; **native** = Rust
implemented, awaiting reference qualification; **validated** = tested against
the exact source baseline for documented cases. Validation of a subset does
not qualify an entire family.

| Capability | Source | Status |
|---|---|---|
| Component/interaction database | fluids/, binaries/, compdatadb.f90, mixdatadb.f90 | pending |
| PR, SRK, RK, van der Waals; classic alpha/vdW mixing | cbmix, cbalpha, cubic, cbhelm | pending |
| Other cubic, alpha, beta, mixing, translation, pseudo fluids | cubic_eos, cbselect, cbbeta, excess_gibbs, wong_sandler, volume_shift | pending |
| Ideal gas, H2, caloric reference states | ideal, idealh2 | pending |
| T-V / T-P properties and derivatives, chemical potentials, sound | eos, eostv, single_phase, speed_of_sound | pending |
| Stable TP/PH/PS/UV/SV/mu-T flashes, grid flashes | stability, tp_solver, ph_solver, ps_solver, uv_solver, sv_solver, mut_solver | pending |
| Pure/mixture saturation and phase envelopes, VLLE | saturation*, puresaturation, binaryplot, multi_phase_envelope_tv | pending |
| Critical points, spinodals, isolines, Joule-Thomson | critical, spinodal, isolines, joule_thompson_inversion | pending |
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
