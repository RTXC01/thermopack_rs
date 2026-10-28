# Analytic PC-SAFT mixture mole derivatives

Mixture fugacity in `PreparedPcSaft` now uses analytic first derivatives of the
extensive reduced residual Helmholtz energy, F = Ares/(RT), at fixed T and V.
Both original PC-SAFT and simplified sPC-SAFT use this path. Pure-fluid
fugacity retains its existing expression. The scope is the nonassociation
hard-chain/dispersion model; polar and association contributions remain open.

The previous implementation rebuilt two perturbed prepared states for each
component. The new implementation differentiates the BMCSL or simplified
hard-sphere terms, contact values, mean segment count, universal dispersion
coefficients and pair moments. Temperature/composition preparation retains
its existing O(nc squared) pair sums; the prepared gradient query takes O(nc)
work and one output vector allocation. The contact gradient sums are shared
between components, and the same vector stores contact logarithms, then mole
gradients, then fugacity coefficients. This is a complexity and allocation
statement. The [measured gradient workload](pc-saft-mole-performance.md)
records prepared and fresh paths separately against the previous Rust code
and the original Fortran kernel.

Construct prepared states with `PcSaft::prepare` or `PreparedPcSaft::new`.
A private stored composition total converts the existing public derivative
fields back to normalized one-mole coefficients, preserving preparation with
unnormalized positive totals. Adding this field changes direct struct-literal
construction in this in-progress Rust API. Existing derivative-field meanings
and model construction methods are preserved.

## Qualification

The frozen source corpus contains 346 isolated attempts: nineteen targeted
inputs in both variants, plus all 77 stored parameter records in binary
mixtures at dilute and dense packing conditions in both variants. The dense
volume is chosen near packing fraction 0.325 from source-derived coefficients.
The other component follows the requested reference's default lookup rules.
Targeted inputs include ternary/quaternary mixtures, zero-feed components,
nonuniform compositions and nondefault parameter references.

There are 344 finite successful original-kernel results. Native construction
and all component gradients match 340 rows within
`abs(native - source) < 2e-11 * (1 + abs(source))`. Four BUTANAL rows have valid
source kernel results but the native constructor requires unavailable default
ideal heat-capacity data. They remain explicit unsupported-construction gaps,
not parity successes or evidence that those original kernels fail. Two further
Default C2_1 source constructor failures remain in the complete raw archive;
nondefault C2_1 reference cases are compared successfully.

The helper selects the actual original PCSAFT_eos or sPCSAFT_eos class and calls
F_PC_SAFT_TVn or F_sPC_SAFT_TVn with only F_n requested. It links the unchanged
library compiled from fork baseline d68c794. This is an internal nonassociation
kernel comparison, including components whose public model might also add
association. It does not qualify full associating public dispatch.

The mixture identity ln(phi_i) = F_n,i - ln(Z) lets the Rust replay pass Z=1
as an adapter for testing F_n directly, including mathematical density states
where physical pressure is nonpositive. This adapter does not describe the
physical fugacity at those states. Separate tests use actual EOS Z and check
finite differences of extensive energy, Euler's identity, composition scaling,
zero-feed limits and fixed-pressure Gibbs-Duhem behavior. Custom m below and
at one checks the original mixture chain's positive-(m-1) convention and the
simplified model's own convention.

## Reproduction and evidence limits

Run `cargo test --test pcsaft_mole_derivatives`,
`python tools/rust_port/capture_pcsaft_mole_derivatives.py --audit`, and
`python -m unittest discover -s tools/rust_port -p test_pcsaft_mole_audit.py`.
The CSV stores only independent original values. The JSONL archive retains
every request, stdout, stderr, return code and stored outcome, including failures.
The audit recomputes outcomes from raw process output, checks the complete
request matrix and provenance, and reconstructs the CSV. Censored processes
fail qualification; the collector kills/reaps each timed-out process and does
not silently retry. The call limit is 45 seconds.

Metadata binds the exact reference DLL, original source hashes, compiler and
helper binary/module/source identities, flags, CSV and raw archive hashes.
The offline audit verifies archived relationships and exact helper/worker
source bytes; it does not recompile or independently regenerate source outputs.
The source DLL and module files are read-only. Helper source bytes are preserved
by Git attributes, including their historical final newline convention.

To capture a new record, use the compiler matching the original reference
metadata and run the collector with `--source` pointing to the original build
checkout and `--compiler-directory` pointing to that compiler's bin directory.
A fresh record replaces the artifacts together; historical hashes are never
updated by hand to mask changed input or output. Expected source construction
failures remain recorded and are not treated as numerical parity successes.

Complete F_nn/F_Tn/F_Vn derivatives, polar and association contributions,
BUTANAL ideal-data construction, broad phase/flash integration and the existing
language interfaces remain pending. A first-gradient qualification does not
close the full PC-SAFT family or the overall Rust rewrite.
