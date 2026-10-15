# Phase stability and caloric flashes

The native solver follows the mathematical framework in `stability.f90` and
`tp_solver.f90`, with immutable Rust model parameters and explicit errors.

## Derivatives

The extensive dimensionless residual Helmholtz function at fixed T,V is
`F = -N*ln(1-B/V) - D*I(B,V)/(R*T)`, where `D=sum(n_i*n_j*a_ij)` and
`B=sum(n_i*n_j*b_ij)/N`. The implementation supports both linear covolume
mixing and nonzero symmetric lij parameters. Differentiating this expression
gives the chemical-potential Hessian. At fixed T,P:

`d ln(phi_i)/d n_j = F_ij + 1/N + P_i*P_j/(R*T*P_V)`.

The T/P derivatives follow the same implicit volume elimination. These
derivatives are checked against finite differences, symmetry and Gibbs-Duhem,
and against 1320 exact-source cases including absent components. The
attraction-integral derivatives use a series in B/V for dilute gases, avoiding
cancellation between almost equal terms.

`ln_fugacity_tp_into` accepts the imposed pressure and solved volume. This
avoids reconstructing low pressure by subtracting large attractive and
repulsive contributions in a dense liquid. That cancellation otherwise adds
noise to tangent-plane line searches close to a trivial solution.

## Stability and split solving

`Cubic::stability(T,P,z)` selects the minimum-Gibbs feed volume root and
performs liquid-like and vapor-like Wilson trials. It minimizes the modified
tangent-plane objective in `alpha_i=2*sqrt(W_i)` coordinates. Three successive
substitutions precede modified Newton iterations with an analytic composition
Hessian, a positive-definite search matrix and a backtracking line search.
Diagonal equilibration before Cholesky factorization prevents a trace species'
large curvature from imposing excessive damping on ordinary components.
As in the fork, a sufficiently negative trial can terminate early: it proves
instability without completing minimization. The returned trial composition
is normalized, while `tpd` refers to the unnormalized modified objective.

`Cubic::tp_flash(T,P,z)` uses negative stability trials to initialize K values,
then successive substitution/Rachford-Rice. Difficult splits minimize total
Gibbs energy in vapor mole-number coordinates, subject to `0 < V_i < z_i`.
The analytic Hessian includes both phases' ideal composition terms and
fugacity derivatives. A converged split must close component balances, satisfy
equal chemical potentials and not exceed the feed Gibbs energy.

Zero-feed components are excluded from minimization and remain exactly zero.
Nonconvergence returns `Error::NotConverged`; it never terminates the process.
`FlashOptions` controls iteration limits and fugacity convergence tolerance.
The implementation does not yet cover full three-phase/global stability,
electrolytes, UV/SV/mu-T flashes or phase envelopes.

## Result conventions and qualification

`FlashPhase` preserves the fork's integer codes: two-phase=0, liquid=1,
vapor=2, and single/undefined type=4. When only one feed volume root exists,
the source returns `Single`, betaV=betaL=-1, x=z and y=0. The Rust result
preserves that convention. Fractions are physical only for the other codes.

The fixture grid attempts 1046 cases across four cubic models, pure fluids,
binary/ternary/heavier mixtures and zero fractions. It includes the exact
eight-component near-critical points from `unittests/test_tp_flash.pf`.

* 1043 source calls return numerical results. Phase flags match exactly.
* 1042 also match the raw phase ordering, fractions and compositions within
  the recorded absolute/relative tolerances (2e-7 each).
* One SRK CO2/water case has exchanged phase labels, but the same physical
  split. Exchange is allowed in the physical comparison only after independently
  checking that both compositions have unique volume roots. Its raw ordering
  compatibility remains unresolved; it is not described as exact label parity.
* Three CO2/water calls terminate the reference process. These failures are
  recorded separately and are not counted as numerical parity successes.

See [baseline-defects.md](baseline-defects.md) and fixture metadata. The oracle
runs in an isolated persistent child process. A Fortran STOP produces a
recorded failure; only the following independent case starts a fresh worker.
Python exceptions fail validation rather than being classified as Fortran
failures. Native input validation, material/fugacity closure and zero/trace
component regressions are also tested independently of fixtures.

## PH/PS and equilibrium caloric properties

`ph_flash(P,z,H)` and `ps_flash(P,z,S)` recover temperature and the phase split,
using SI molar enthalpy and entropy with the model's reference states. Options
control a temperature guess, bounds (80 to 999 K by default), iteration limits,
caloric residual tolerance and the inner TP solver. The native solver uses
physical derivatives, including for vdW; see the baseline defect notes.

The solver first tries a cheap homogeneous caloric solve, accepting it only
after a TP stability check and caloric closure. Otherwise it brackets the
equilibrium target and uses safeguarded Newton steps with nested TP flashes.
The preliminary homogeneous search is limited to 16 root iterations: its
minimum-Gibbs branch can jump at coexistence and leave no homogeneous root.
Shared ideal integrals and residual derivatives are evaluated once per phase.
A failed TP trial may be retried toward a valid bracket endpoint. Every
accepted value comes from an actual successful evaluation; an unresolved
interval returns an explicit error. This is not a complete multiphase solver.

`equilibrium_properties` returns V, H, S, U and optional equilibrium Cp.
For a split, Cp includes phase amount and composition redistribution:
`Cp = Cp_frozen + R*T^2 * delta^T * Hessian^-1 * delta`, where delta is the
difference of phase log-fugacity temperature derivatives and the Hessian is
in vapor mole-number coordinates. Independent H/S finite differences test
this response. Pure two-phase latent heat has no finite isobaric Cp, so its
heat-capacity field is `None`.

Pure caloric flashes invert saturation pressure and calculate quality from
latent H or S, including feeds with only one nonzero component. The
`pseudo_critical` helper finds the cubic EOS critical point at fixed
composition; it is not a true mixture critical point. `pseudo_saturation_*`
also holds both phase compositions fixed and is an initialization estimate,
not a mixture bubble/dew calculation. Qualification includes 28 pseudo-critical
cases and 36 pure saturation-temperature inversions.

The fork's PH and PS output conventions differ. Pure PH/PS fills both phase
composition arrays with z. For mixtures, PS retains TP conventions. PH fills
both arrays with zero for single-phase results; below pseudo-Pc it classifies
the phase using pseudo-Vc. Above pseudo-Pc its Single result uses betaV=-1,
betaL=0. These are compatibility sentinels, not physical fractions.
`equilibrium_properties` uses the feed for single-phase results.

Of 912 attempted source PH/PS cases, 697 return states that independently
close H or S within 1e-8 in H/(RT) or S/R. All 697 match native temperature,
raw phase codes, fractions and composition arrays within fixture tolerances.
The remaining 215 are recorded source failures or non-closing results; they
do not establish numerical parity. Exact compatibility with those outputs
and broad global phase stability remain pending.
