# Baseline discrepancies

Exact baseline: d68c794. These are observed differences, not hypothetical
limitations. Compatibility and mathematical consistency are tested separately.

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
