# Numerical contract

Use the fork's gas constant `6.02214076e23 * 1.380649e-23` J/(mol K),
not a rounded textbook value. Temperatures are K, pressures Pa, volumes
m3/mol, energies J/mol and compositions mole fractions.

Conventional cubic pressure is `RT/(v-b) - a(T)/[(v+d1*b)(v+d2*b)]`.
For PR, `d1=1+sqrt(2)`, `d2=1-sqrt(2)`; for SRK/RK, `d1=1`, `d2=0`;
for vdW both are zero. Critical omega coefficients are calculated with
the same iteration as `cbCalcOmegaZc` rather than rounded 0.45724/0.07780.
Classic alpha coefficients come from `get_alpha_param_acentric`.

Database selection must match the source's reference-selection rules. Aliases,
critical data overrides, non-default reference sets and all kij/lij variants
must eventually be covered. No silent approximations for unsupported models.
For DEFAULT selection, the source skips exact-substring matching and retains
the last default entry. Interaction constants use the compiled database's
rounded values, rather than potentially higher-precision JSON values.
Numerical pure-fluid constants and default Cp records likewise use the compiled
`compdatadb.f90`, while JSON supplies aliases. Type-10 hydrogen ideal values
have a source scaling discrepancy: `IdealProperties::cp` preserves the source
correlation and `enthalpy_dt` differentiates its returned H. Native caloric
iterations use the latter.

Analytic derivatives are compared to central finite differences independently
of parity tests. Cubic roots are filtered by admissible volume and polished
against the pressure equation. Stable phase selection minimizes residual Gibbs
energy among admissible roots; this alone does not establish mixture stability.
Flashes require tangent-plane stability tests before being production-qualified.
See [baseline-defects.md](baseline-defects.md) for the vdW derivative
compatibility method; physical heat capacities remain consistent with values.

The reference is a local build from the baseline source. Fixtures record source
hashes, compiler, flags, case inputs, quantities and tolerances. Relative errors
are insufficient for quantities near zero. Phase-boundary comparisons must
record phase identity, closure residuals and branch selection.
