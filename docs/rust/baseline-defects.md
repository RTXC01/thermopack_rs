# Baseline discrepancies

Exact baseline: d68c794. These are observed differences, not hypothetical
limitations. Compatibility and mathematical consistency are tested separately.

## Alpha initialization and constructor overrides

The alpha initializer reads TWU records for both TWU and MC requests. Several
explicit correlation names leave their fitted coefficient zero. Non-water UMR
and explicitly set GERG/UMR terminate on evaluation because the source lacks
those evaluator branches. Native alpha selection preserves the effective
method/parameters and returns explicit errors for unimplemented evaluation.

Separately, the public cubic constructor overrides alpha to TWU and enables
translation for tcPR references. QuantumCubic also changes beta and critical
data. Full tcPR construction is now qualified separately from the isolated
alpha kernel; QuantumCubic construction/properties are qualified separately
while its cryogenic flash grid remains pending. See [beta.md](beta.md), [alpha.md](alpha.md) and
[volume-translation.md](volume-translation.md) for scope and source probes.

The Python `init_peneloux_volume_translation(reference)` wrapper drops its
reference argument before calling Fortran, selecting DEFAULT parameters.
The public cubic constructor does pass its reference. Native APIs preserve
both behaviors through distinct entry points, with explicit documentation.

## Cryogenic TP returns are not always stable equilibrium states

The exact fork returns a PR H2/NE split at 3 K, 1 kPa, z=0.6/0.4 with
`Gsplit/(RT)=-22.8596711116`, while the homogeneous liquid feed has
`Gfeed/(RT)=-26.0314067678`. Chemical-potential residuals are zero, but the
split has higher Gibbs energy by 3.1717356562 RT. For z=0/0.5/0.5 in HE/H2/NE,
the corresponding increase is 0.5135988542 RT. The SS path in
`tp_solver::twoPhaseTPflash` can return immediately on convergence before its
later Gibbs comparison is reached.

At 80 K, 0.1 MPa, z=0.5/0.3/0.2 in C1/CO2/N2 under QuantumCubic construction,
the source also returns TwoPhase with a vapor fraction below 1e-12 and a
vanishing-phase composition that does not close chemical potentials. The PR
case has a maximum dimensionless residual of 14.5247. SRK/RK exhibit the same
kind of return. This is a raw phase-output convention/defect, not a valid VLE
state merely because total H/S or material balance closes.

All five observations are independently recalculated by
`tools/rust_port/audit_quantum_flash.py` and stored with the exact source
fingerprints in [quantum_flash_audit.json](../../tests/fixtures/quantum_flash_audit.json).
The raw cases remain in the reference data and native failure counts.
Complete source-convention compatibility for these returns remains pending;
the validated QuantumCubic property kernel does not qualify its entire flash
surface. Nonmonotonic cryogenic H/S curves also admit multiple closing
temperatures, so equality to the state used to generate H/S is insufficient
to establish equality to the source solver's returned branch.

## Compiled fluid constants differ from JSON

`src/compdatadb.f90` rounds molecular masses when it is generated. For example,
R1234yf is 114.0416 g/mol in the compiled fork and 114.041593 in JSON. This
affects its mass-based ideal caloric correlation and exceeds the chosen parity
tolerances. DME and four hydrogen isomers also have rounded compiled masses.
The Rust generator reads numerical fluid constants and default Cp records from
the compiled database, retaining JSON names/aliases. It records hashes of both
sources. The selected default Cp record still follows last-default-wins order.

The public cubic constructor passes a length-100 reference to the Cp database's
length-40 labels. Its case-sensitive `INDEX(label, reference)` cannot match,
so `getCpDBindex` always uses the last `Default` record. Empty/slash references
are qualified for CO2/NC16, which also have unlabelled alternate Cp records;
the constructor cannot select those records. Native construction preserves
this default fallback rather than interpreting the reference as Cp selection.

## Hydrogen ideal H/S and Cp use different scaling

For correlation 10, `src/idealh2.f90` multiplies its H and S primitives by 1000
but not Cp. Consequently numerical dH/dT is 1000 times the reported ideal Cp.
The routine uses its own fitting gas constant, 8.3144772 J/(mol K). N-H2 and
E-H2 select this correlation by default; P-H2 and O-H2 select correlation 12.
The source also omits its last two oscillator terms at T <= 30 K.

`Component::ideal` preserves source H, S and Cp, including reference states.
Its `enthalpy_dt` field explicitly reports the derivative of those H values.
Native `heat_capacities` and caloric solves use that derivative;
`baseline_enthalpy_derivative_tp` preserves the source Cp convention. This
maintains the distinction between reference parity and value consistency.
The default ideal-property grid qualifies 2592 cases for 108 fluids. BUTANAL's
default correlation is 100, which terminates the reference and returns an
explicit unsupported-correlation error in Rust; all 24 attempts are recorded.

## Ideal reference getters and setter order

`set_ideal_cp_correlation` changes coefficients without recalculating the H/S
integration constants. The reference-energy setters do recalculate them,
except when both declared references are zero or a sentinel exceeds 1e15.
Therefore setter order affects returned values. The native owned `IdealGas`
preserves these rules and caches the integration constants. The API rejects
oversized arrays rather than applying the Python wrapper's silent truncation.

`get_standard_entropy` subtracts `R*ln(1.01325)` when the declared reference
pressure is 1ATM. Converting actual entropy to 1bar would add this term.
`IdealGas::standard_entropy` preserves the source getter; use
`properties(298.15, 1e5).entropy` for the actual value at 1bar. All correlations
and six setter sequences are qualified in 720 source comparisons, including
nonzero raw API44 enthalpy constants and the DIPPR entropy integration constant.

## Duplicate interaction defaults

For PR CO2/H2O, both JSON and compiled `src/mixdatadb.f90` contain duplicate
default kij entries, 0.074 and -0.065. `exact_substring_match` deliberately
skips DEFAULT, so `get_binary_db_idx` selects the **last** default, -0.065.
Selecting the first matching JSON record instead gives 0.074.
Rust interaction constants therefore come from the compiled Fortran database,
in its selection order, with last-default-wins lookup. Updating JSON without regenerating Fortran changes
what a source-based port would otherwise select. Neither source is modified.

## van der Waals enthalpy temperature derivative

`src/cubic.f90`, `cbCalcDerivatives_svol`, vdW branch assigns
`fftt = -2*DD/n1**2`; differentiating its own `fft = DD/(T*n1)` instead gives
`-2*DD/(T**2*n1)`, since attraction is temperature independent.

At methane, T=200 K and P=1 MPa, the fork's `enthalpy(..., dhdt=True)`
returns approximately 37754.1077 J/(mol K). Differentiating the enthalpy
values at fixed pressure gives approximately 35.0868 J/(mol K).

`PreparedState::baseline_enthalpy_derivative_tp` preserves the fork's result.
`PreparedState::heat_capacities` returns physical Cp/Cv and is tested against
enthalpy/internal-energy finite differences. The compatibility correction is
`2*a*T/(1000*v**2) - 2*a/(T*v)` in SI molar units. The factor 1000 comes from
the source's conversion of cubic volume to L/mol. Native sound speed uses
physical heat capacities. Additional upstream derivative interfaces involving
the erroneous F_TT remain pending in the coverage ledger.

## TP flash reference aborts and phase ordering

At T=200 K, P=1 MPa, z(CO2,H2O)=(0.9,0.1), PR, SRK and SRK with RK
alpha terminate the calling process: the modified Newton solution has higher
Gibbs energy than the feed. The complete outputs are in
`tests/fixtures/tp_flash.metadata.json`. These calls are attempted and recorded,
but cannot provide numerical expected values. Rust returns a result or an
explicit error without terminating its caller.

At SRK, T=300 K, P=5 MPa and the same composition, the baseline reports
betaV approximately 0.0990607025 and x(CO2) approximately 0.998871541;
Rust reports betaV approximately 0.900939298 and x(CO2) approximately
0.000781095. Swapping phases gives the same fractions and compositions.
Both phase compositions have a unique volume root, so root selection cannot
enforce the labels. Qualification records this permutation explicitly.
Exact raw phase-order compatibility is still pending.

## PH/PS failures and non-closing successful returns

The 912-case caloric grid records 14 reference process terminations, 106
wrapper calculation failures and 95 results that return successfully but do
not close the specified H or S. Of the latter, 82 are vdW and 13 are PR/SRK/RK.
All failed requests, source output where available, and non-closing residuals
are retained in `tests/fixtures/caloric_flash.metadata.json`.

The vdW derivative defect above also affects caloric solvers; it is not the
explanation for every failure in this grid. The oracle computes the requested
property independently at the returned state and requires residuals below
1e-8 in H/(RT) or S/R before accepting a numerical result. It recognizes only
the wrapper's known calculation exceptions/warning; unexpected Python errors
fail validation. PS's permitted status -1 is also subject to closure checks.

The 697 closing results qualify raw output parity for those cases only.
Native caloric solvers use physical Cp and require closure. Reproducing the
fork's non-closing outputs remains an explicit compatibility gap, rather
than silently treating them as a correct physical solution.

## MBWR integrated-term cancellation

The original MBWR Helmholtz expansion subtracts large density-independent
polynomial/exponential terms. At extreme dilute extrapolations from 0.5 to
10 K, rounding in that subtraction differs from a stable analytic integral
in 318 of the 6760 qualified kernel states. The stable path independently
satisfies derivative and pressure identities there. The differences and both
sets of values are retained in the MBWR comparison report.

Native `baseline_residual_density`, `baseline_residual_helmholtz` and
`baseline_pressure_tv` preserve the original integrated coefficients and
summation order. All 6760 source-order states match without loosening the
stable-path tolerances. See [MBWR formulation and qualification](mbwr.md) for
the exact coefficient checks, zero-density extension and remaining interfaces.

## MBWR entropy pressure-response signs

`MBWR_Sres` computes `P_V=-R*(T*F_VV-n*T/V^2)` and
`P_n=-R*(T*F_Vn+T/V)`. These differ from the physical pressure derivatives
used in `MBWR_Hres`, `MBWR_zfac` and `MBWR_lnphi`. In the ideal-coefficient
limit, its residual entropy value is zero but its T/P responses become
`-2*n*R/T` and `2*V/T`, instead of zero.

`Mbwr::baseline_residual_tp` preserves those source expressions.
`residual_tp` uses physical pressure response and verifies S_T=H_T/T against
independent finite differences. See [MBWR state functions](mbwr-properties.md)
for scope, derivative singularities and source snapshot execution context.

## MBWR default initialization bounds

`readDbParameters` uses `compdb(getCompDBindex(id))` without checking for a
zero lookup result. Seventeen legacy MBWR identifiers lack that general data;
their default initialization accesses out of bounds. The reference probe
diagnoses this before execution. Rust returns an error while retaining the
coefficient-only kernels and explicit-metadata density solver.

The coarse critical search can also reach its result without assigning a lower
bound. Native code reports that undefined path explicitly. The qualified matrix
contains 32 finite default initializations, 16 source terminations and 17
undefined lookups, plus 27 direct critical-grid comparisons. See
[MBWR initialization](mbwr-initialization.md) for the original loop rules,
model-local cache and the distinction from analytic critical refinement.

## Ice Ih heat-capacity mass conversion

`sh2o_d2gdT2` multiplies its mass-specific result by the H2O molar mass.
`sh2o_heat_capacity` then calculates `-T*sh2o_d2gdT2` and multiplies by the
molar mass again. Its value is therefore 0.01801528 times physical molar Cp.
Native `SolidWater::baseline_heat_capacity` preserves that source value;
`heat_capacity` uses -T*G_TT and matches finite differences of enthalpy.
All other public source values and reference shifts are compared directly.
See [ice Ih qualification](solid-water.md).

## FMT hyperdual reciprocal derivatives

`PowerHyperDualReal` uses the actual value in x**(-1), but replaces |x|<1e-15
by a signed floor when calculating the derivative powers. The resulting
derivatives differ from mathematical differentiation of the FMT expression
near tiny coefficient denominators. The native `baseline_energy_density`
preserves the floor and source operation order; `energy_density` computes
analytic derivatives without the reciprocal floor. Both retain the source
max(1-n3,1e-12) branch and raw, cancellation-sensitive White Bear formulas.
See [FMT qualification](fmt.md) for singular states and numerical scope.

## Hyperdual arc tangent third derivative

`atanHyperDual` writes `deriv**3**v1%f1` for the last term of its third
order derivative, which is a double power: the exponent of the derivative is
three raised to the first derivative of the argument. The native
`hyperdual::HyperDual::atan` preserves that expression, and the
[hyperdual qualification](hyperdual.md) confirms that all eight components of
all twelve elementary functions agree with the exact-source library.

## Lennard-Jones meos sound mass units

`constructor_LJ` copies `compdb%mw` in g/mol into `meos%molarMass`, a field
documented and used as kg/mol by the direct `speed_of_sound` method. Its
sound is smaller by sqrt(1000). The Rust `baseline_speed_of_sound` preserves
that source value; normal sound and mass access use SI units. The source
matrix checks the baseline and an independent thermodynamic identity checks
physical sound. See [pure multiparameter qualification](multiparameter.md).

## MBWR ideal references and Python TV packing

`init_multiparameter` copies component data with zero H/S integration constants
and does not call `set_reference_energies`. Its MBWR calorics therefore use raw
Cp integrals rather than the formation/standard references applied by other
constructors. Reference setters normalize both constants; Cp setters preserve
them. The [total MBWR layer](mbwr-caloric.md) qualifies this ordering explicitly.

The Python `internal_energy_tv` wrapper appends `np.array(dedv_c)` when the
mole derivative is requested, instead of returning the `dedn_c` array. On the
qualified NumPy runtime this raises a pointer-buffer-format error. The original
Fortran routine supplies all three correct TVn derivatives; the caloric audit
calls that routine directly and retains the complete derivative block.

## Multiparameter ideal mixture is unreachable and inconsistent

`src/multiparameter_idealmix.f90` mixes the multiparameter equations of state
with an ideal mixing rule, but the fork cannot construct that model. The
`NIST_MEOS_MIX` label reaches `single_eos_allocate_and_init`, which accepts
that label for a multicomponent model, and that routine then calls
`single_eos_alloc` with the same label. `single_eos_alloc` only handles
`NIST_MEOS`, `LJ_MEOS` and `LJTS_MEOS` and stops with "Wrong input for
multiparameter EoS." for `NIST_MEOS_MIX`. The label therefore always
terminates the caller, and the mixture model is unreachable through the
public API.

The internal `..._Gres` routine assigns weighted residual Gibbs energy to
its pressure derivative. The public non-residual `TP_CalcGibbs` mixture path
also uses an uninitialized `sumne` while adding ideal contributions.
`..._fugacity` returns coefficients `phi`, although its caller names the
output `lnfug`. These conventions require explicit treatment in a port.

A later class-safe probe of the unchanged routines does **not** reproduce
the earlier reported shifted-array observations. Nor are H/S pressure
derivatives undefined in that probe: the local derivatives are passed
unconditionally to the pure component routines, which assign them before
they are used. The previous observations are not treated as a source-kernel
contract. The new 734-case audit retains complete vectors, root checks and
the failed public construction call.

The [native residual mixture kernel](meos-ideal-mixture.md) is qualified
against those internal routines with a probe that constructs the correct
dynamic type using supported pure labels. This bypasses the allocation
defect only in the validation helper; it does not modify the original
library or establish public mixture-wrapper parity. Total properties and
mixture phase solvers remain pending.
## Iso-line grid of a single point is not a number

src/isolines.f90:linspace divides by the number of intervals of the grid,
so a grid of one point holds a not-a-number value instead of the limit of
the grid. The pure fluid branch of isotherm and isobar splits the grid at
the saturation state, and the number of points of the two branches follows
the original formula

    n_vap = ceiling(n (P_sat - P_min) / (P_max - P_min))

reduced by one when it exceeds half of the grid. That count can be one, and
the original then reports a not-a-number pressure with not-a-number volume,
entropy and enthalpy for that point. The isotherm of methane at 150 K from
1e5 Pa to 1e7 Pa with ten points is such a case, because the saturation
pressure of methane at 150 K is 1047350.03 Pa and the vapor branch holds a
single point.

The native model reproduces the not-a-number value, because the parity of the
grid depends on it, and the qualification compares the not-a-number values of
both models.
See [iso-line qualification](isolines.md).

## The pressure of a molar volume of the cubic equation of state

src/cubic.f90:cbCalcPressure takes the molar volume as the argument  and
computes the pressure with p = cbPress(cbeos,T,V), where V is neither a
local variable of the routine nor an argument of it. The module does not
switch off the implicit typing of the routine, so V is an uninitialized
implicit real. The pressure of the routine therefore does not follow from the
supplied volume, and the derivatives dpdv, dpdt and d2pdv2 follow the
mixture parameters of the recalculated state, not the state of the argument.

The state of a mixture with the Wong-Sandler mixing rule is therefore not
evaluated with cbCalcPressure in the reference probe. The probe uses the
generic routine eostv::pressure, which evaluates the pressure of the state
from the volume and the mole numbers, as the fugacity and the caloric
properties of the same state do. The native library evaluates the pressure
from the volume of the prepared state and is qualified against the generic
route.
See [Wong-Sandler notes](wong-sandler.md).

## Mixture spinodal pressure walk and curve endpoint pressures

`src/spinodal.f90::nestedLoopVT` selects the volume-step sign from the pressure
error alone, even on the vapor branch where pressure has the opposite slope.
All 54 vapor requests in the mixture spinodal grid return the liquid-side
boundary. `baseline_mixture_spinodal_point` follows this source seed/walk
selection; the physical API starts from the requested mechanical branch.
Both roots satisfy local stability and pressure equations, so closure alone
cannot establish branch-selection parity.

`map_stability_limit` replaces T/V at a temperature-limited last point but
retains the pressure evaluated before refinement, and assigns the input P0
to the first point even with an explicit starting temperature. All 22 complete
qualified mixture curves contain a stale endpoint pressure. Raw source P and
evaluated EOS P are both retained in the fixtures. The compatibility curve
preserves reporting; the physical curve evaluates the final pressure. See
[mixture spinodal qualification](mixture-spinodal.md) for scope and limits.

## Spinodal property intersection returned volume

`propertyFunctionWrapperSpinodal` predicts T from the scalar trial volume,
corrects its local V at fixed T to evaluate the property, and stores only T
back into the parameter vector. `bracketSolveForPropertySpinodal` returns the
uncorrected trial as Vs. Thus a successful source result can fail both the
spinodal and property conditions at its reported coordinates. The new pure
compatibility search preserves this convention; all successful raw coordinates
are retained and compared. In the qualified grid 145 returns have a minimum
Helmholtz eigenvalue exceeding 1e-6 in magnitude. Physical searches refine
temperature at fixed volume and close both constraints. See
[spinodal property notes](spinodal-properties.md) for near-grid ties and the
separate critical multiple-root rounding allowance.

## Excess Gibbs fallback interaction initialization

`cbselect::SelectBinary` requests mixing label `Classic` for the attraction
interactions used by HV1/HV2/NRTL fallback pairs. This fork's `kijdb` stores
those classical records with label `vdW`, so the request returns zero even
when the normal vdW constructor finds a nonzero interaction. The GE path
also leaves covolume interactions at zero. For PR methane/water, for example,
the vdW constructor selects kij=0.5, while HV1 selects zero.

The Rust GE constructor preserves these defaults. Explicit `set_kij` updates
both the public parameter and fallback evaluation. The fallback uses the
geometric mean of the already temperature-dependent attractions; the previous
Rust formula counted those attractions twice. Four exact-source pressure and
entropy states and independent classic/HV fallback identity checks now cover
this previously unqualified path, including temperature derivatives.

## Metastable-isentrope first and terminal coordinates

`map_meta_isentrope` starts with an equilibrium PS flash, then evaluates the
majority phase's volume root at the feed composition. In a two-phase mixture
this reported homogeneous state need not have the target entropy. Its assigned
spinodal endpoint also inherits the uncorrected property-search volume above.
In the current 97 complete qualified maps, 25 first/last coordinates have
pressure or entropy errors above the documented closure threshold. Fixtures
retain those errors. The compatibility mapper preserves coordinates while
recording actual residuals; physical continuation corrects its explicit seed
and every endpoint. See [metastable notes](metastable.md).

## Pure spinodal extrema selection and pressure coordinates

`locate_spinodal_prop_min_max_pure_fluid` keeps the last branch turning
segment, or its first point when the branch is monotone. It compares only
these two branch candidates, rather than including every endpoint. Its final
`Ps_min = pressure(Ts_max,vs_min,z)` mixes maximum temperature and minimum
volume. All 462 finite returns in the current grid exhibit a minimum-pressure
coordinate mismatch. The compatibility API preserves this formula.

For pressure, both candidates are intended to be at P0. Strict ordering by
their residuals makes the labels depend on iteration roundoff. Of 83 such
source returns, 50 have reversed native primary labels. The native result
marks this ambiguity; qualification checks both complete assignments while
retaining the original ordering and derived alternative pressures. This is
candidate/reporting parity, not full raw ordered-result parity. See
[pure extrema notes](spinodal-extrema.md).

## Envelope pressure-cap stagnation and formulation labels

`envelope_plot_epc` clamps log(P) to log(Pmax), but checks the exponentiated
pressure against Pmax using a strict comparison. With Pmax=2 MPa on this host,
rounding leaves the pressure just below the cap. Both PR and SRK methane/ethane
cases fill the 1000-point buffer with repeated terminal states. Baseline
envelope tracing preserves the arrays; the physical API recognizes the
logarithmic bound and reports `MaximumPressure`.

The source can switch its bubble/dew formulation once both compositions have
one EOS volume root, before their density order changes. Baseline phase labels
and volumes retain this convention. Physical tracing orders reported phases
by molar volume and changes the fixed-phase label consistently. This reporting
choice is separate from global stability. See [envelope notes](envelope.md).
