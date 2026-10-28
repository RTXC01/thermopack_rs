# PC-SAFT first-gradient measurements

The [frozen record](pcsaft-mole-performance.json) compares the complete internal
nonassociation F_n vector for 32 mixture profiles. On this machine, all profile
medians favor the analytic Rust implementation. This is a bounded kernel
measurement, not a whole-library or full associating-model speed claim.

| Comparison | Minimum / median / maximum ratio across profiles |
|---|---:|
| Previous Rust / analytic Rust, prepared | 11.17 / 21.95 / 72.41 |
| Previous Rust / analytic Rust, fresh preparation | 4.33 / 8.38 / 20.36 |
| Original Fortran / analytic Rust, prepared | 15.92 / 22.64 / 41.54 |
| Original Fortran / analytic Rust, fresh preparation | 5.66 / 6.45 / 8.20 |

Each ratio uses that profile's seven-trial median. The middle column statistic
is the median of 32 profile ratios, not aggregate application throughput.
Across profiles, prepared analytic medians range from 265 to 964 ns per call;
fresh analytic medians range from 1,145 to 2,485 ns. Trials show substantial
variation: the largest sample standard deviation divided by its group median
is about 21% for prepared analytic Rust and 45% for original Fortran. The raw
minimum, maximum, standard deviation and all samples remain in the record.
No outliers were removed and no failed timed calls were omitted.

## Workload and measurement boundaries

Profiles have 2, 3, 4 or 8 components, both PC-SAFT variants, equal feeds or a
zero first component, T=350 K and densities 100 or 1,000 mol/m3. The JSON's
`dilute`/`dense` labels distinguish those two densities within this workload;
they do not establish dense-liquid or phase-stability coverage. The separate
[derivative qualification](pc-saft-mole-derivatives.md) includes the much denser
packing-fraction controls. Neither benchmark chooses a phase or solves a flash.

Model construction is outside every timed loop. Prepared Rust reuses a
prepared temperature/composition state. Fresh Rust constructs that state on
every evaluation while reusing the model. Original Fortran invokes its normal
F_PC_SAFT_TVn/F_sPC_SAFT_TVn with F_n requested and retains the model's ordinary
internal cache behavior. Its model initialization and Python/ctypes overhead
are outside the native Fortran loop. Thus prepared reuse is an explicit API
advantage; the fresh path records its preparation cost rather than hiding it.

The mixture identity `fugacity_coefficients(rho, 1) = F_n` is the Rust adapter.
Pure-fluid special cases are excluded. All vector components contribute to a
weighted checksum in both languages, with black-box barriers in the Rust
loop. This avoids timing only pressure while unused gradient work disappears.
Seven rotated/reversed implementation orders use 1,000 warmups and 20,000
iterations per timed call. There are 96 full-vector preflights and 1,120 valid
timed calls across 160 groups. The candidate is checked against original
values at 2e-11*(1+abs(source)); the historical finite-difference path has a
separate 2e-6 bound. Every timed vector and checksum is checked again.

The run used an AMD Ryzen 7 1700, Windows 10, Rust 1.99.0/LLVM 23.1.1 and
GNU Fortran 14.3.0. The original DLL is the unchanged d68c794 reference.
Rust uses the repository's release library profile and an identical standalone
harness compiled with -O and thin LTO; the Fortran helper uses -O3 and the
reference's double-precision conventions. Exact commands, compiler identities,
source hashes, Git differences, executable/library/module hashes and host
information are retained in the record. Baseline Git is 9bd3c68; the candidate
is that checkout plus the captured analytic-gradient diff.

## Reproduce and audit

`benchmark_pcsaft_mole_derivatives.py build --help` describes the coordinated
build. Build mode compiles both release libraries, compiles the identical
standalone Rust harness against each rlib, links the original-loop helper and
writes a manifest. The two checkouts and original reference build are explicit
inputs. Build completes before any timing starts.

Use measurement mode with `--manifest`, `--output`, `--samples 7`,
`--iterations 20000`, `--warmup 1000` and `--timeout 45`. `--preflight-only`
checks the vectors without collecting the seven-trial timing record. The
run rejects changed frozen inputs, failed/censored calls, incomplete groups,
nonfinite outputs and checksum mismatches. Process timeouts kill and reap the
child. Calls and failures are retained, with no silent retries.

Run `python tools/rust_port/audit_pcsaft_mole_benchmark.py` to audit the
committed historical record offline. Mutation tests are in
`test_pcsaft_mole_benchmark_audit.py`. The audit binds parsed values to raw
stdout and recomputes accuracy checks, checksums, sample identities/order,
statistics and ratios. Captured before/after hashes establish internal
consistency; the offline audit does not independently attest execution on the
historical hardware or rebuild the old binaries. Optional current-input checks
are separate because subsequent development can legitimately change the tree.

Future measurements still need other temperatures, liquid/critical conditions,
association/polar terms, full state and solver workloads, batch/caller-buffer
APIs, other machines and the existing language interfaces. These timings do
not establish complete interchangeability with ThermoPack.

The exact four measured harness inputs are preserved in the compressed
source snapshot beside the record, with hashes bound to its build manifest.
Current tool files have normalized line endings and trailing whitespace; this
formatting cleanup does not replace or refresh historical input hashes. Restore
the archived bytes in an isolated checkout when reproducing that exact input
record. The snapshot metadata also binds the unmodified timing record.
