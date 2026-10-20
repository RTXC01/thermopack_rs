# Iso-lines

An iso-line is a curve in the temperature-pressure plane that holds one
property of the fluid constant. The source maps an isotherm, an isobar, an
isenthalp and an isentrope.

## Source scope

| Original routine | Rust |
| --- | --- |
| `src/isolines.f90:isotherm` | `isolines::isotherm` |
| `src/isolines.f90:isobar` | `isolines::isobar` |
| `src/isolines.f90:isenthalp` | not yet ported |
| `src/isolines.f90:isentrope` | not yet ported |

The isenthalp and the isentrope of the original solve a pressure-enthalpy and
a pressure-entropy flash at every point of the grid. The native model has no
pressure-enthalpy and no pressure-entropy flash, so those two lines are not
yet ported.

## Method

The original splits the grid of a pure fluid at the saturation state, so that
every point is evaluated in the phase that is stable there, and it runs a
temperature-pressure flash at every point of the grid of a mixture, with the
properties of the two phases weighted by the phase fractions. The native
model follows the same construction.

The number of points of the two branches of a pure fluid follows the original:

    n_vap = ceiling(n (P_sat - P_min) / (P_max - P_min))

and the count is reduced by one when it exceeds half of the grid. The
saturation state of the original is the bubble point of the pure fluid, and
the native model uses its fixed composition saturation solve.

## The single point grid defect

The original `linspace` divides by the number of intervals, so a grid of one
point holds a not-a-number value instead of the limit of the grid. The split
of a pure fluid grid at the saturation state can leave a single point in one
of the two branches, and the original then reports a not-a-number pressure
with not-a-number properties. The native model reproduces that value, and the
qualification compares the not-a-number values of both.

## Validation

`tools/rust_port/validate_isolines.py` compares the native model against the
exact-source library. The oracle calls the Python interface of the fork,
which reaches the original routines.

* 9 cases: isotherms and isobars of methane and of carbon dioxide in the
  Peng-Robinson and the Soave-Redlich-Kwong equations of state, and isotherms
  and an isobar of the methane-propane and the methane-carbon dioxide
  mixtures.
* The compared quantities are the pressure or the temperature of the point,
  the molar volume, the molar entropy and the molar enthalpy. All 9 cases
  match with a relative tolerance of 1e-6, and the observed worst relative
  difference is 8.9e-10. The difference is the saturation state of the
  original and of the native model, which agree to a relative 1e-8.

The committed fixture is `tests/fixtures/isolines.csv`; the test is
`tests/isolines.rs`.

## Remaining scope

* The isenthalp and the isentrope need the pressure-enthalpy and the
  pressure-entropy flash.
* The saturation-line crossing of the original, which inserts the saturation
  point into a mixture grid, is inactive in the original, because the number
  of points of the grid is never reduced there.
* The iso-lines of the multiparameter, the PeTS and the other models are not
  yet ported.
