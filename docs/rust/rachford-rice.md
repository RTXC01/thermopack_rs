# Extended Rachford–Rice

`numerics::rachford_rice_extended(z, K, workspace, single_step)` ports
`tp_solver::rr_solve`. It finds roots outside [0,1] when the extended material
balance permits a negative-flash initializer. This differs from the existing
`rachford_rice`, which clips to a physical endpoint if there is no interior root.
Neither function establishes phase stability or fugacity equilibrium.

`RachfordRiceState::new(n)` allocates reusable composition buffers. The solve
itself allocates no memory. The source's composition-dependent bounds, choice
of vapor/liquid fraction, 1500-iteration limit and epsilon-scaled convergence
tests are retained. The liquid fraction is stored separately: it can remain
positive when subtracting the vapor fraction from one rounds to zero.

The one-step mode uses the supplied vapor fraction, computes a Newton update
and recomputes normalized phase amounts. It deliberately preserves the input
liquid-fraction field, matching the fork. Non-finite updates return an explicit
native error. A false return means there is no new RR result; the source's
output fields are undefined on that path and are not qualified. Native arrays
not written by the failure path retain their earlier contents and must not be
interpreted as a successful solution.

The direct linked Fortran probe qualifies 698 finite results/statuses from 704
attempts. These include deterministic random mixtures with 2, 3, 8 and 20
components; K values across 24 decades; explicit trace/zero fractions; near-one
K values; no-root cases and one-step updates. Six source one-step calls produce
non-finite values, recorded separately. Success flags require exact equality;
fractions/compositions use the stored absolute/relative tolerances. Independent
checks cover negative fractions, material balance, buffer reuse and a liquid
fraction below machine resolution at a vapor endpoint.

This kernel is groundwork for reproducing source flash initialization and
convergence paths. It has not yet replaced the optimized physical RR helper in
the native TP solver; the 82 cryogenic flash differences remain pending.
