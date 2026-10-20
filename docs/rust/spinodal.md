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

## Remaining scope

* The mapping of the whole limit line, the source `map_stability_limit`, is
  not yet ported, and the native model reports a single point instead of a
  curve.
* The limit of a mixture needs the stability matrix of the mixture with its
  eigenvalues, and the native model reports the single component case only.
* The Joule-Thomson inversion of the source, which builds on the limit line,
  is not yet ported.
