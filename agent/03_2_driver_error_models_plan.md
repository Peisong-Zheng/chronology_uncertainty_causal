# Common-driver chronology models: focused Z-only sensitivity

Authorized on 10 November 2026. Implement in the root script `03_2_check_driver_error_models.py`. This branch extends the completed 03_1 crossing while keeping X and Y accurately dated. It asks whether the prescribed Z chronology model changes the ability of conditional methods to exclude spurious X-Y links and retain a genuine direct X-to-Y link. It does not assess recovery of the full three-variable network.

## Fixed systems and comparison

Retain **Common driver only** and **Common driver + direct X to Y**. Their equations, gains, 20-year delay, Z/X/Y memories of 500/200/200 years, zero measurement noise, sampling positions, 30-kyr histories, 29-1 kyr analysis window and 1401-point analysis grid remain those of 03. Reuse the three generated histories with seeds 20261009, 20261301 and 20261302. No process history or sampling mask is regenerated.

For every realization and system, X/Y are the existing irregular accurate-age analysis inputs. Only the ages assigned to its sampled Z values change. Each Z age scheme is applied to all three realizations and both systems. Both original-uniform and irregular accurate-age references remain available, with the latter used for paired score changes and in-figure baseline decisions.

The tested assignments, ordered X/Y/Z, are `(accurate, accurate, C)`, `(accurate, accurate, L)`, `(accurate, accurate, E)` and `(accurate, accurate, J)`. The C-C background extension, an equal-RMS shape control, memory changes and nonlinear systems remain future decisions; none is part of this run.

## Chronology definitions and random pairing

Use the existing `make_age_uncertainty`, `simulate_observations` and `draw_age_ensemble` functions without changing their numerical algorithms. Preserve the main study's level-3 definitions:

| Type | Proposal standard deviation |
|---|---|
| C | Constant 40 years. |
| L | Linear increase from zero at 0 kyr BP to 40 years at 30 kyr BP. |
| E | Exponential increase over the same 0-30 kyr bounds, curvature 2, maximum 40 years. |
| J | Five smooth compact peaks, each with a 10-year maximum and a 2.8-kyr defined interval. |

J centers use the first five established fractions `[0.5, 0.2, 0.8, 0.35, 0.65]` of the 29-1 kyr analysis interval. Their non-overlapping defined intervals cover 50% of that interval. This is not a threshold-exceedance fraction or a step discontinuity. C/L/E have their specified old-end maximum at 30 kyr, rather than the 29-kyr analysis boundary. These are idealized uncertainty envelopes, not archive classifications.

For each of three scheme seeds, reuse the original Z observation/ensemble seed formulas, with record index 2: `SeedSequence([seed, 3103, 2])` and `SeedSequence([seed, 3104, 2])`. Retain the normalized-sample-order correlation length 1/60. Reusing these seeds across types pairs their underlying random draws; differing whole-curve rejection histories can change which proposals are accepted. Do not assume member m represents the same candidate draw across types when rejections differ.

Each type/scheme receives its own uncertain observed Z chronology. Its ten ordered members are drawn around that observed chronology, not around the true ages. Regenerate the C arrays locally as a verification only and require exact equality with the existing C arrays; then copy the saved C files and results. Generate nine new L/E/J schemes with native calls. Retain acceptance rates, attempted/accepted counts, order rejections and observed-center order violations. Do not sort, clip, repair or rescale accepted curves. Require every analysis input to cover the common grid without extrapolation; an invalid input must stop preparation rather than silently shrink its window or be replaced with a favorable draw.

## Interpretation of uncertainty magnitude

This is a comparison of the **complete prescribed uncertainty settings**, including their different amplitudes and distributions over age. It is not an isolated test of envelope shape. Equal peak SDs do not imply equal overall proposal RMS; the J peak amplitude also differs from C/L/E.

Report the proposal-envelope RMS over the common analysis interval, the observed-center RMS error, the RMS of each member's error relative to true ages, and the effective spread of the ten members. For an error series e on the 1401-point analysis grid, `RMS = sqrt(mean(e**2))`. Effective spread is the mean pointwise population SD across the ten age-error curves on that same grid; it is descriptive ensemble spread, not an uncertainty in the mean or a calibrated credible interval. Truth is used only for diagnostics. No observed center, member, or accepted input is selected or scaled using its realized error or causal result.

X/Y timing error remains zero. Relative X-Z and Y-Z age errors therefore equal the Z errors, whereas relative X-Y error is zero. Three different age schemes per type are retained separately: their temporal error patterns can matter in addition to their RMS values.

## Scale and exact reuse

- Uncertain inputs: 2 systems x 3 realizations x 3 schemes x 4 types x 10 members = **720**.
- Accurate-age references: 2 systems x 3 realizations x 2 references = **12**.
- Total cases: **732**.
- Full imports from 03_1: 180 C cases and 12 references = **192**.
- New L/E/J inputs: **540**, with three conditional method calls each = **1620 new fits**.
- Four pairwise methods, including CCM, reuse each realization/system's accurate-age evidence after checking actual parsed X/Y arrays and estimator seeds. No pairwise fit is repeated.
- Each case retains seven complete native archives and 26 directional rows: **5124 archives** and **19032 rows**. Of these archives, **3504 are exact copies** and **1620 are new conditional results**. Keep full native Z-related fields without adding a full-network recovery metric.

The conditional methods remain native CPU PCMCI+/GPDC, corrected LKIF and the existing statsmodels conditional GC wrapper. Use 32 outer processes and one numerical thread per worker, existing GPDC null caches, no nested pools and no GPU substitution. Method seeds remain fixed within realization as in 03. Consequently a row comparison retains its original estimator seed and is not a pure process variance component.

## Outputs and figures

Save inputs, the twelve full Z age packages, diagnostics, native results, lag tables, committed member edges and the three existing kinds of X-Y summaries under `data/processed/03_2_check_driver_error_models/`. Include `age_type` in case identities, result rows and all type-specific groupings. The complete 26-row member table remains the final commit. Source files from 03 and 03_1 are read-only.

Exact X-Y recovery means neither direction for Common driver only, and only X to Y for Common driver + direct X to Y, with no unresolved evidence. Divide by all requested members. Keep explicit incorrect arrows, missing true arrows, unresolved evidence, failed fits and pending cases. Native-score differences use every finite successful pair relative to the same realization's irregular accurate-age reference. No cross-method physical interpretation of score magnitudes is intended.

Three figures will be exported as PNG/PDF with lowercase panel letters and detailed English captions:

1. **00_error_models:** show the four proposal envelopes, actual Z error magnitudes and effective ensemble spread/acceptance. State the unequal overall uncertainty and the difference between the full 0-30 kyr envelope domain and the 29-1 kyr analysis interval.
2. **01_xy_recovery:** two systems by three conditional methods. Within each panel, rows are realizations and columns are C/L/E/J. Show exact correct/requested counts out of 30, unresolved counts, and each row's single accurate-age reference. Thirty means three schemes with ten members each, not thirty independent histories.
3. **02_scheme_detail:** retain the same six method/system panels and realization/type layout, but show three short bars per cell for the three schemes, with correct/10 and unresolved counts. This preserves interactions like the original GC R1/A2 anomaly instead of hiding them in a pooled cell.

Scheme colours and ordering remain consistent. Missing results are explicit and never displayed as observed zero success. Pairwise evidence is unchanged and retained in tables rather than repeated as additional result panels. Captions belong in `agent/figure_review/08_driver_error_models.md`.

## Validation and execution

Before production, record a preservation snapshot of pre-existing scientific sources and outputs. Check source scientific settings, packages, platform, complete results and source hashes; no old metadata is edited to bypass a guard. Verify all 732 prepared identities, native C reconstruction, all 192 exact imported inputs, new envelope and chronology generation, system/realization pairing, order/coverage, and unchanged X/Y as read by the actual default CSV parser. Validate failure retention, summary denominators, type grouping, input/source guards, partial commit behavior and same-run resume.

Compare one predetermined new input's full three conditional native outputs against direct calls to the unchanged 03 methods and Windows-spawned calls. This is a small implementation check, not selective scientific refitting. Run the 540 new cases only after this preflight passes. Freeze the calculation source before the run; presentation-only updates must not change calculation settings or outputs.

After production, audit every native archive and copied file, rebuild all three summaries exactly, check all plotted counts against saved directional rows, and inspect the actual PNG/PDF figures. Compare hashes, sizes and modification times of every pre-existing scientific file to the initial snapshot. Retain unfavorable results and unresolved decisions without selectively repeating them. Write an English result report explaining accurate-age exceptions, type and scheme differences, and the limited generality of three existing histories.

## Execution and recovery record

Preflight passed in 102.59 seconds at `data/benchmarks/driver_error_models_checks/20261110T133558073115Z/report.json`. The first production invocation began at 13:38:47 UTC on 10 November 2026 under `data/logs/driver_error_models_20261110T133847522941Z/`. The user accidentally closed the program at about 13:49 and explicitly requested continuation. A read-only checkpoint audit found 297 complete cases (192 full imports and 105 new L/E/J cases), compared with the last periodically written status count of 294. All 7722 committed directional rows had successful fit status. All 732 inputs and the source, platform, settings and computation guards remained valid; no fits were repeated in the audit.

Recovery started at 13:54:03 UTC under `data/logs/driver_error_models_resume_20261110T135403003994Z/`, loaded the 297 committed cases and launched only the 435 incomplete cases with the same 32-worker settings. The old invocation remains recorded as interrupted. An additional snapshot protects the 3267 files belonging to already committed cases (inputs, native archives, lag tables and member tables); their hashes, sizes and modification times must be checked again after recovery. Incomplete native output is not a complete case and may be replaced while finishing that case. The resumed status timer describes the recovery invocation, not total time across both invocations.

Recovery completed all 732 cases at 14:36:22 UTC. The resumed scientific script took 2535.030 seconds (42.25 minutes); the interrupted invocation had about 10.2 minutes of recorded monitor time. These are separate execution records, and incomplete work may have been repeated during recovery. The final outputs contain the planned 1620 new conditional archives and 3504 exact source copies. All 19032 directional rows have successful fit status. The artifact audit passed all native evidence and exact reconstruction of the three summary tables. See the [results report](03_2_driver_error_models_results.md) and the [figure captions](figure_review/08_driver_error_models.md).

Final artifact, figure and both preservation audits passed in the recovery directory. All 28595 pre-existing scientific files and the 3267 files belonging to the 297 pre-recovery commits retain their content hashes, sizes and modification times. The final figure check covered all 72 type cells, 216 scheme cells and 18 accurate reference entries without rerunning any causal fits.
