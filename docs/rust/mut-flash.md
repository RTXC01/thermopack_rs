# mu-T flash

The mu-T flash solves for the component densities at a fixed temperature and a
fixed specification. The native model implements the log fugacity variant of
the original solver and the chemical potential variant that the original
derives from it.

## Source scope

| Original routine | Rust |
| --- | --- |
| `src/mut_solver.f90:solve_lnf_t` | `mut_flash::solve_lnf_t` |
| `src/mut_solver.f90:solve_mu_t` | `mut_flash::solve_mu_t` |
| `src/mut_solver.f90:lnf_t_function` | `mut_flash::log_fugacity_and_potential` |
| `src/eosTV.f90:thermo_tv` | `mut_flash::log_fugacity_and_potential` |
| `src/eosTV.f90:chemical_potential_tv` | `mut_flash::log_fugacity_and_potential` |

## Model

The specification is the logarithm of the fugacity of every component, and the
solver finds the densities that reproduce it at the given temperature. The
original evaluates the specification with `thermo_tv`, which reports the
logarithm of the fugacity, and converts a chemical potential specification
into that logarithm at the initial guess.

The native model follows the same structure: the specification is evaluated at
the state, and a damped Newton iteration solves the residual. The iteration
uses a finite-difference Jacobian of the analytic log fugacity, while the
original uses the second composition derivatives of the Helmholtz energy, which
are a different approximation of the same derivative. Both converge to the same
root of the same residual.

The chemical potential of the mixture follows from the fugacity coefficients:
the ideal part is `-R T (ln(x) + 1)` and the residual part is `R T ln(phi)`, so
the fugacity coefficients carry the whole chemical potential.

## Source limitations reproduced

* **The chemical potential variant is a fixed point.** The original derives
  the log fugacity specification from the chemical potential at the initial
  guess, so the specification is already satisfied there and the solve returns
  the initial guess. The native model keeps that behaviour, and the test
  verifies it.
* **The log fugacity is not the fugacity coefficient.** The original reports
  the logarithm of the fugacity. The native model reports the same quantity
  and uses it in the specification.
* **The original solves in a formal region.** The original evaluates the log
  fugacity at any state, including the states with a negative pressure inside
  the two-phase region. The native model reports an error there, because its
  fugacity coefficient needs a positive pressure. The solver halves a step
  that leaves the usable region.

## Validation

`tools/rust_port/validate_mut_solve.py` compares the native solver against the
exact-source library through a reference probe that builds the source model,
evaluates the target specification and solves the problem from a different
guess.

* 32 states over four mixtures and four temperatures, with the specification
  taken at the true state and the solve started from a different guess.
* All 32 states agree within a relative tolerance of 1e-6; the worst observed
  relative difference is 7.4e-13.
* The committed fixture is `tests/fixtures/mut_solve.csv`; the test is
  `tests/mut_flash.rs`.

## Remaining scope

* The isotherm mapping from the saturation line to the spinodal needs the
  spinodal point location, which is not ported.
* The mixture envelopes of the multiparameter models are not ported.
