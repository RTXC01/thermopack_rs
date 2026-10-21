# Stability limits

The stability limit of a fluid is the meta-stability limit, the curve where
the isotherm has an extremum and the first volume derivative of the pressure
vanishes. The source maps that limit and reports the point of it at a given
pressure.

## Source scope

| Original routine | Rust |
| --- | --- |
| `src/spinodal.f90:initial_stab_limit_point` | `spinodal::spinodal_point` |
| `src/spinodal.f90:map_stability_limit` | not yet ported |

The original solves the point of the limit with a nested loop that walks the
limit line from an initial point, and it refines the crossing with a
bracketing solver on the same condition. The native model solves the two
conditions of the point directly, the pressure of the state and the vanishing
first volume derivative of the pressure, with a bisection on the volume and a
bisection on the temperature.

The original works with the stability matrix of the mixture, the Hessian of
the Helmholtz energy with respect to the mole numbers, and the lowest
eigenvalue of that matrix. For a single component the original uses the
first volume derivative of the pressure instead, in `stabFunSingleV`, and the
native model follows that condition. The native model therefore reports the
same point as the Python interface of the fork, which exposes the original
routine as `spinodal_point`.

## The two branches

The liquid branch of the limit has the smaller molar volume and the
isotherm minimum, and the vapor branch has the larger molar volume and the
isotherm maximum. The pressure of both branches rises with the temperature
and reaches the critical pressure at the critical temperature. The pressure
of the liquid branch falls without bound towards the zero temperature, and
the pressure of the vapor branch falls to the zero pressure.

## Validation

`tools/rust_port/validate_spinodal.py` compares the native model against the
exact-source library. The oracle calls the Python interface of the fork,
which reaches the original routine.

* 12 cases: the liquid and the vapor branch of methane and of carbon dioxide
  in the Peng-Robinson equation of state and the two branches of methane in
  the Soave-Redlich-Kwong equation of state.
* The compared quantities are the temperature and the molar volume of the
  point. All 12 cases match with an absolute 1e-6 K on the temperature and an
  absolute 1e-10 cubic metre per mole on the volume, together with a relative
  1e-8. The observed worst relative difference is 1.4e-9, which is the
  tolerance of the source solver.

The committed fixture is `tests/fixtures/spinodal.csv`; the test is
`tests/spinodal.rs`.

## Critical points

The critical point of a pure fluid is the state where the isotherm has an
inflection, so the first and the second volume derivative of the pressure
vanish. The original solves the two conditions with the stability matrix of
the mixture, in `calcCriticalTV`, and for a single component the two
conditions of that matrix are the volume derivatives of the pressure.

The native model, in `critical`, solves the two conditions with a Newton
iteration in the temperature and the molar volume, starting from the pseudo
critical state of the model. The native pure Jacobian is a centred difference.
The source mixture solver differentiates its eigenvector and uses perturbed
composition states; its Jacobian is not the same pure pressure formulation.
True local mixture critical points are now available through the same public
function; see [mixture critical points](mixture-critical.md) for their separate
Helmholtz equations and qualification.

`tools/rust_port/validate_critical.py` compares the native model against the
exact-source library. The oracle calls the Python interface of the fork,
which reaches the original routine.

* 7 cases: methane, carbon dioxide, propane and water in the Peng-Robinson
  equation of state, nitrogen and methane in the Soave-Redlich-Kwong equation
  of state, and a mixture with a zero mole fraction of the second component.
* The compared quantities are the temperature, the molar volume and the
  pressure of the critical point. All 7 cases match with a relative 1e-6,
  and the observed worst relative difference is 3.8e-7, which is the
  tolerance of the source solver.

The committed fixture is `tests/fixtures/critical.csv`; the test is
`tests/critical.rs`.

## The Joule-Thomson inversion

The inversion curve of the source, in `src/joule_thompson_inversion.f90`, is the
curve of the states where the Joule-Thomson coefficient vanishes. The original
solves the condition

$$
T \left(\frac{\partial P}{\partial T}\right)_v
+ v \left(\frac{\partial P}{\partial v}\right)_T = 0
$$

for the temperature at a given molar volume, and it maps the curve with that
solution. The native model, in `joule_thompson`, solves the same condition
with a Newton iteration in the logarithm of the temperature, and it uses the
same centred difference for the derivative of the condition that the original
uses when it extrapolates along the curve.

`tools/rust_port/validate_joule_thompson.py` compares the native model
against the exact-source library. The oracle maps the inversion curve of the
original through the Python interface of the fork, and the native model
solves the condition at every molar volume of that curve with a perturbed
initial guess.

* 78 cases: the inversion curve of nitrogen in the Peng-Robinson and the
  Soave-Redlich-Kwong equations of state.
* The compared quantities are the temperature and the pressure of the point.
  All 78 cases match with an absolute 1e-6 K on the temperature and an
  absolute 1e-4 Pa on the pressure, together with a relative 1e-7. The
  observed worst relative difference is 3.3e-9.

The first point of the curve of the original is the zero pressure solution,
which the original computes from the virial coefficients and not from the
inversion condition, so the comparison starts at the second point. The
original also fails to map the curve of methane and of carbon dioxide, and
the nitrogen curve is the reference.

The committed fixture is `tests/fixtures/joule_thompson.csv`; the test is
`tests/joule_thompson.rs`.

## Remaining scope

* The mapping of the whole limit line, the source `map_stability_limit`, is
  not yet ported, and the native model reports a single point instead of a
  curve.
* The limit of a mixture needs the stability matrix of the mixture with its
  eigenvalues, and the native model reports the single component case only.
* Critical loci, critical end points and a global phase-stability classification
  of mixture critical points remain pending; the local point solve is native.
* The curve of the Joule-Thomson inversion is not mapped by the native model,
  and it reports a single point instead.
