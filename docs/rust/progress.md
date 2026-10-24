# Overall implementation progress

The current **engineering estimate is about 35%**, with uncertainty of roughly
ten percentage points. This estimates progress toward a complete native Rust
replacement of this fork, including all model families, public capabilities,
numerical qualification and measured performance. It is not the percentage of
functions already interchangeable with ThermoPack, or an estimate of elapsed
or remaining development time.

The [coverage ledger](coverage.md) and [source inventory](inventory.json) are
the evidence for the estimate. A tested subset does not complete its whole
model family. The following weights make the estimate's assumptions explicit;
the per-area percentages are engineering judgments, not automated test scores.

| Workstream | Weight | Estimated completion | Weighted contribution |
|---|---:|---:|---:|
| Native EOS kernels, models and parameter databases | 30% | 65% | 19.5 points |
| State properties and complete derivative capabilities | 20% | 55% | 11.0 points |
| Solvers, stability, phase diagrams and model integration | 25% | 15% | 3.75 points |
| Existing language interfaces and distribution compatibility | 10% | 0% | 0 points |
| Whole-surface parity/performance qualification and delivery | 15% | 5% | 0.75 points |
| **Total** | **100%** | | **35 points** |

Many cubic capabilities and several non-cubic state kernels have strong
reference evidence. Remaining work includes broader non-cubic solvers and
derivatives, association/polar and corresponding-state capabilities, complete
phase diagrams and global stability, remaining solid/hydrate and DFT behavior,
and the existing Python/C/C++ interfaces. Small qualified additions should not
be mistaken for closing these much larger gaps.

After each commit the chat reports the rounded overall estimate and the
concrete completed change. It may remain unchanged across several commits.
The pure MEOS TP/PH/PS addition qualifies a bounded surface across fourteen
models, with 1,890 independent source-kernel checks. Its 36 endpoint
coordinate gaps, incomplete response/phase-diagram capabilities and unmeasured
flash performance leave the overall estimate at about 35%.
Reassess the workstream estimates when a material gap closes, and record why;
do not increment the percentage merely because a commit or test was added.
The final 100% requires an exhaustive completion audit against the full user
objective and the checked-in fork. No such completion claim is made here.

Pure MEOS response and sound add 854 strict source comparisons across fourteen
variants. This closes a bounded response subset; broad singular cases, global
stability, mixture integration and performance measurements remain pending.
The overall estimate remains about 35%.

Pure MEOS fixed-volume equilibrium and UV/SV add 560 original-kernel native
closure checks, 240 qualified source returns and 320 source-failure recoveries.
The 59 endpoint phase-flag differences and missing global-stability, broad
derivative/integration and performance evidence keep the overall estimate
at about 35%.
