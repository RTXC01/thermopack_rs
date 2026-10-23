# Activity-based cubic fugacity derivatives

The prepared cubic model now retains the full attraction Hessian for UNIFAC,
HV1/HV2 and NRTL mixing. `PreparedState::fugacity_derivatives_into` consumes
this Hessian when evaluating derivatives with respect to mole numbers at
fixed temperature and pressure. Caller-owned output buffers remain supported.

The previous code supplied the classic geometric-mean attraction Hessian
even when the mixture attraction came from an excess Gibbs model. That
inconsistent Jacobian could slow or prevent nonlinear continuation, although
pressure and fugacity values were correct. The existing UNIFAC evaluator
already returned its analytic extensive Gibbs Hessian; the cubic mixer now
uses it. HV/NRTL Hessian terms also had transposed interaction indices in
their inner sums and two quotient terms. These indices now follow
`excess_gibbs::gEinf`.

For extensive attraction A, covolume B and excess Gibbs energy G, the mixer
uses the source identity

```text
A_ij = (B_ij A + B_j A_i + B_i A_j)/B
       - 2 B_i B_j A/B² - B G_ij/C
```

Here C is the cubic excess-Gibbs limiting constant. The full-covolume entry
point also accepts first and second covolume derivatives; the original
linear-covolume helper remains available. Prepared states store the resulting
A_ij matrix once. The classical and Wong-Sandler coefficient paths retain
their separate Hessians. Mixed GE/quantum beta constructions and arbitrary
nondefault covolume interactions are not qualified by the present cases.

## Qualification

`tools/rust_port/validate_ge_derivatives.py` compares all 378 finite source
returns against the library compiled from this fork. Every T, P and mole
derivative of every log fugacity coefficient matches. There were no source
failures, censored calls, nonfinite returns or source volume-pressure warnings.

The cases cover 250, 300 and 400 K; 0.1, 1 and 10 MPa; liquid and vapor root
selection; PR/SRK UNIFAC binary and ternary mixtures, including a zero species
and translated models; SRK/PR HV1/HV2 methane/water; PR HV1 water/ethylene
glycol with Hemmingsen2011 data; and PR NRTL hydrogen sulfide/propane with
Dicko2012 data. Default missing-pair HV and NRTL fallbacks are included,
alongside a PR H2S/propane/methane model with fitted and missing NRTL pairs.

The no-record NRTL branch stores its off-diagonal `C_ji = b_j` in the
source's L/mol units. Native covolumes use m³/mol, so this branch explicitly
multiplies by 1000. Its diagonal is dimensionless one, preventing that scale
from cancelling between numerator and denominator. Fitted NRTL pairs and the
HV fallback retain their own conventions. The added 126 source derivative
states and independent mole perturbations exercise this correction.

Tolerances are relative 2e-7 with absolute floors 1e-9 K⁻¹, 1e-11 Pa⁻¹ and
1e-8 mol⁻¹. Worst scaled errors, using max(1, |source|, |native|), are
7.36e-15 for T, 1.20e-17 for P and 2.28e-13 for mole derivatives. Fixture
metadata records exact source/library identity, coefficients provenance,
call outcomes and the fixture hash.

Independent tests perturb each mole number and renormalize the phase while
holding T/P fixed. They compare the resulting fugacity differences with every
matrix entry, and also check symmetry and the Gibbs-Duhem row sums. Existing
GE values, Wong-Sandler coefficients/states and mixture saturation tests
remain separate regressions. This qualification does not establish every
activity-model critical Hessian or full phase-envelope algorithm parity.
