# Baseline discrepancies

Exact baseline: d68c794. These are observed differences, not hypothetical
limitations. Compatibility and mathematical consistency are tested separately.

## Compiled fluid constants differ from JSON

`src/compdatadb.f90` rounds molecular masses when it is generated. For example,
R1234yf is 114.0416 g/mol in the compiled fork and 114.041593 in JSON. This
affects its mass-based ideal caloric correlation and exceeds the chosen parity
tolerances. DME and four hydrogen isomers also have rounded compiled masses.
The Rust generator reads numerical fluid constants and default Cp records from
the compiled database, retaining JSON names/aliases. It records hashes of both
sources. The selected default Cp record still follows last-default-wins order.

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
