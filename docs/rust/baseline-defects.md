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
