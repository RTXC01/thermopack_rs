# Pseudo components

A pseudo component is a component whose critical data are supplied by the
user instead of the fluid database. Gas processing applications use them for
the heavy fractions of a reservoir fluid and for the pseudo components of a
distillation column.

## Source scope

| Original routine | Rust |
| --- | --- |
| `src/eoslibinit.f90:init_cubic_pseudo` | `cubic::Cubic::new_pseudo` |

The source routine replaces the critical temperature, the critical pressure,
the acentric factor and the molar mass of every component named `PSEUDO` and
then selects the cubic equation of state for the resulting component list. The
lists are indexed by the position of the component, and the values of the
other components are ignored.

The native constructor `Cubic::new_pseudo` follows the same convention. The
component named `PSEUDO` takes the critical temperature, the critical
pressure, the acentric factor and the molar mass of its list entry, and the
remaining components keep the values of the source database. The critical
temperature and the critical pressure set the parameters of the cubic
equation, and the acentric factor enters the classic alpha correlation of the
equation of state. The molar mass is the property of the ideal gas part.

The source applies the pseudo values to the component list of the preceding
initialization, because the Python interface calls the cubic initialization
first and the pseudo initialization second. The native constructor takes the
component list and the pseudo values in one call, so the result does not
depend on a previous initialization.

## Validation

`tools/rust_port/validate_pseudo_components.py` compares the native model
against the exact-source library. The oracle builds the same model through the
Python interface of the fork with `init_pseudo` and evaluates the properties
of the mixture.

* 30 cases: a single pseudo component, two pseudo components, a pseudo
  component with methane and a pseudo component with carbon dioxide, in the
  Peng-Robinson and the Soave-Redlich-Kwong equations of state, in the liquid
  and the vapor phase.
* The quantities are the molar volume, the enthalpy, the entropy and the
  compressibility factor. All 30 cases match with a relative tolerance of
  1e-10 and quantity-specific absolute tolerances. The observed worst
  relative difference is 7.9e-15.

The committed fixture is `tests/fixtures/pseudo_components.csv`; the test is
`tests/pseudo_components.rs`.

## Remaining scope

* The Wong-Sandler mixing rule and the other cubic geometry variants of
  `src/cubic_eos.f90` are not yet ported.
* The pseudo component is not connected to the pseudo critical data of the
  quantum cubic construction, and the critical data of the source can be
  redefined after the initialization with `redefine_critical_parameters`.
