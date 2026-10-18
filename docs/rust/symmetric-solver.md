# Symmetric indefinite solve

`numerics::solve_symmetric_in_place(matrix, rhs, pivots, work)` ports the upper
unblocked Bunch–Kaufman factorization from the pinned LAPACK DSYTF2, followed
by DSYTRS2's conversion, permutations and triangular/diagonal solves. It uses
native Rust, without a BLAS/LAPACK runtime. The BSD-3-Clause copyright notices
and terms are retained in `rust/LAPACK-LICENSE`; the crate's combined license
declaration reflects this port alongside Apache-2.0 ThermoPack code.

The matrix is row-major and only its upper triangle is used. The lower triangle
is never read or changed. A call overwrites the upper triangle with the source
factorization and fills signed one-based pivots. Positive entries denote 1x1
blocks; paired negative entries denote 2x2 blocks. The RHS becomes the solution
on `Ok(true)`. On `Ok(false)`, a zero pivot was found: the RHS stays unchanged
and the factor/pivots describe the singular factorization. Invalid dimensions,
non-finite inputs in the used triangle/RHS and non-finite computed outputs
return explicit errors. The caller must treat an error as having no solution.

Provide exactly n pivot slots and at least n floating-point work slots for an
n-variable system. All buffers can be reused; the routine allocates no memory.
It does not estimate conditioning or perform iterative refinement. It remains
unblocked for all sizes, rather than reproducing LAPACK's workspace-dependent
large blocked path. Full blocked-factor trajectory parity is not claimed.

The linked reference probe calls the exact pinned DSYTF2/DSYTRS2 routines.
All 226 attempts return finite results/statuses and match: dimensions 1, 2, 3,
8, 20 and 32; positive definite, indefinite and diagonally scaled systems;
2x2 pivots, permutations, singular systems and arbitrary unused lower entries.
Status and pivots require exact equality. Factor/solution tolerances are 1e-12
absolute plus 2e-12 relative. Four systems are singular and retain their input
RHS. Native tests independently check componentwise backward error, including
scaled systems, rather than assuming small forward error in ill-conditioned
problems. See `tests/fixtures/symmetric.metadata.json` for provenance.

The fork-order TP Newton path uses this solver and falls back to the source's
modified LDL direction for a singular system or non-descent direction. The
change preserves every previously matching raw TP comparison, but does not
resolve its three ghost-vapor composition differences. This is a verified
numerical foundation; no speedup is claimed without measurement.
