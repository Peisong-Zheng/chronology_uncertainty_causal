# Common-driver chronology models: X-Y results

Completed on 10 November 2026. The [plan](03_2_driver_error_models_plan.md), [root script](../03_2_check_driver_error_models.py) and [figure captions](figure_review/08_driver_error_models.md) define the comparison. X/Y retain accurate ages; only Z has uncertain ages. The outcome is the X-Y relation, with no assessment of full three-variable network recovery.

## What changes when the Z error model changes

In **Common driver only**, the C-type exceptions become less frequent under the prescribed L/E/J settings. For PCMCI+, realization R2 gives 4/30 correct X-Y decisions with C, compared with 24/30, 29/30 and 30/30 with L, E and J. For GC, realization R1 gives 18/30, 28/30, 30/30 and 30/30. LKIF correctly retains neither X-Y arrow for all four types and all three realizations.

In **Common driver + direct X to Y**, the response is less uniform. PCMCI+ and GC retain the true forward arrow in all 360 uncertain inputs, but PCMCI+ still has unresolved or extra reverse evidence. Its R2 counts are 18/30, 20/30, 13/30 and 14/30 for C/L/E/J; that realization's irregular accurate-age reference was already unresolved. LKIF detects the true forward arrow in R1 for 28/30 C inputs, but only 1/30 L inputs and none of the E/J inputs. Both accurate-age references already miss the forward link in this LKIF realization, so the C result cannot be taken as evidence that larger dating errors improve inference.

These comparisons concern complete error settings. Their overall magnitudes differ substantially: median realized Z RMS errors are 52.5, 29.9, 22.5 and 6.2 years for C/L/E/J. The results do not isolate envelope shape at a common error magnitude.

## Comparison and reading order

Each type has three Z age schemes, A1/A2/A3, generated with the same corresponding observation and ensemble seed formulas across types. Each scheme consists of an uncertain observed center and ten ordered members drawn around it. Every scheme is applied to every one of the three saved process realizations, R1/R2/R3, under both systems. No dynamics, sampling masks or X/Y chronology arrays were regenerated. Memories remain 500 years for Z and 200 years for X/Y, and each realization retains its existing estimator seed.

The study contains 720 uncertain inputs and 12 accurate references. All 180 C-type inputs and both accurate references for each realization/system were imported from 03_1. The nine L/E/J chronology packages are new. The same three histories are reused throughout; the crossing does not create nine independent climate histories.

Read [Figure 00: error envelopes and diagnostics](../figures/03_2_check_driver_error_models/00_error_models.png) first, then [Figure 01: X-Y recovery by type](../figures/03_2_check_driver_error_models/01_xy_recovery.png). Figure 01 contains all 72 realization-by-type cells across the six system/method panels. [Figure 02: scheme detail](../figures/03_2_check_driver_error_models/02_scheme_detail.png) resolves those cells into all 216 separate scheme counts. Each figure has a PDF companion.

In the following tables, every type entry gives **correct/30; U = unresolved count**. Thirty means three age schemes with ten requested members each for one fixed realization. It is not thirty independent histories. Exact recovery requires neither X/Y direction in Common driver only, or only X to Y in Common driver + direct X to Y, without unresolved evidence. Failed or pending cases would remain in the denominator; none occurred in the completed results. An unresolved state can still retain a directional arrow, so unresolved and directional-error counts are not generally additive.

## Common driver only

All conditional accurate-age references retain neither X-Y direction, under both original-uniform and irregular accurate sampling.

| Method | Realization | Accurate, correct/1 | C | L | E | J |
|---|---|---:|---|---|---|---|
| PCMCI+ | R1 | 1/1 | 29/30; U=0 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 |
| PCMCI+ | R2 | 1/1 | 4/30; U=16 | 24/30; U=2 | 29/30; U=0 | 30/30; U=0 |
| PCMCI+ | R3 | 1/1 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 |
| LKIF | R1 | 1/1 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 |
| LKIF | R2 | 1/1 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 |
| LKIF | R3 | 1/1 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 |
| GC | R1 | 1/1 | 18/30; U=0 | 28/30; U=0 | 30/30; U=0 | 30/30; U=0 |
| GC | R2 | 1/1 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 |
| GC | R3 | 1/1 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 |

The directional breakdown explains the remaining exceptions. Across the three realizations, PCMCI+ has five incorrect X-to-Y arrows, six incorrect Y-to-X arrows and sixteen unresolved states under C. Under L it has two incorrect arrows in each direction and two unresolved states, all in R2. Under E it has only one incorrect Y-to-X arrow, again in R2; under J there are no X-Y errors or unresolved states. GC's errors are entirely Y-to-X: twelve under C and two under L, all in R1. LKIF has no incorrect X-Y arrows under these settings.

The scheme decomposition matters. Entries below are correct counts for A1/A2/A3, each out of ten:

| Existing exception | C | L | E | J |
|---|---|---|---|---|
| PCMCI+, R2 | 3 / 0 / 1 | 8 / 10 / 6 | 10 / 10 / 9 | 10 / 10 / 10 |
| GC, R1 | 8 / 1 / 9 | 10 / 8 / 10 | 10 / 10 / 10 | 10 / 10 / 10 |

Thus the original C-type R2 sensitivity is not equally strong under every prescribed type. The previously prominent GC R1/A2 exception also weakens: its nine incorrect reverse arrows under C become two under L and none under E/J. The C R1/A2 reverse support remains 0.9 despite no true X-Y link. It is one example of strong within-scheme agreement on a wrong direction, not a population false-positive rate.

## Common driver plus direct X to Y

The accurate column uses each realization's single irregular accurate-age reference. The original-uniform reference is also saved and is relevant to the two exceptions discussed below.

| Method | Realization | Accurate, correct/1 | C | L | E | J |
|---|---|---|---|---|---|---|
| PCMCI+ | R1 | 1/1 | 25/30; U=3 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 |
| PCMCI+ | R2 | 0/1: unresolved | 18/30; U=10 | 20/30; U=10 | 13/30; U=16 | 14/30; U=16 |
| PCMCI+ | R3 | 1/1 | 17/30; U=10 | 22/30; U=8 | 21/30; U=9 | 29/30; U=1 |
| LKIF | R1 | 0/1: true arrow missed | 28/30; U=0 | 1/30; U=0 | 0/30; U=0 | 0/30; U=0 |
| LKIF | R2 | 1/1 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 |
| LKIF | R3 | 1/1 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 |
| GC | R1 | 1/1 | 28/30; U=0 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 |
| GC | R2 | 1/1 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 |
| GC | R3 | 1/1 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 | 30/30; U=0 |

PCMCI+ retains the true X-to-Y arrow in **360/360** uncertain inputs. Under C/L/E/J, its extra Y-to-X counts are 7/0/1/0, and its unresolved counts are 23/18/25/17. Its lower exact-recovery counts therefore reflect added reverse or unresolved evidence, not missed forward arrows. GC also retains X-to-Y in 360/360 inputs; its only incomplete X-Y decisions are the two C-type R1/A2 inputs with an added reverse arrow. All L/E/J GC inputs have the correct complete X-Y decision in this experiment.

Two baseline exceptions prevent a simple monotonic interpretation:

- **PCMCI+, R2:** the original-uniform input gives the correct X-Y state, but the irregular accurate-age input already has unresolved evidence while retaining X-to-Y. Under C/L/E/J, correct A1/A2/A3 counts are respectively 5/6/7, 6/6/8, 3/4/6 and 5/5/4 out of ten. Smaller overall Z errors do not remove this baseline limitation consistently.
- **LKIF, R1:** both accurate references miss the true X-to-Y arrow. Its C/L/E/J missed-forward counts are 2/29/30/30, with no reverse or unresolved evidence. Under L, the sole forward detection occurs in A2; under E and J, all three schemes retain neither X-Y direction. Changing Z timing changes the conditioning input. The C result's greater forward detection does not establish more reliable adjustment or a general benefit of dating uncertainty.

LKIF retains X-to-Y in every R2/R3 input under every type. These realization-specific exceptions would be obscured by reporting only a pooled method ranking.

## Prescribed scale, realized error and accepted members

The native envelopes are evaluated over the full retained 0-30 kyr BP range. C/L/E have a prescribed maximum of 40 years, with E curvature two. J has five fixed 2.8-kyr intervals, each with a prescribed 10-year smooth peak, spanning 50% of the 29-1 kyr analysis window by total interval width. The envelopes were not rescaled to equal RMS. All diagnostics below use the common 29-1 kyr grid.

| Type | Proposal-envelope RMS, yr | Observed-center RMS median, yr | Member RMS median [min, max], yr | Effective spread median, yr |
|---|---:|---:|---:|---:|
| C | 40.0 | 33.7 | 52.5 [44.3, 71.3] | 35.0 |
| L | 22.7 | 19.8 | 29.9 [21.8, 41.5] | 17.5 |
| E | 16.8 | 14.7 | 22.5 [14.7, 34.4] | 11.8 |
| J | 5.0 | 3.8 | 6.2 [3.3, 9.6] | 2.6 |

Observed-center and effective-spread medians summarize the three schemes. Member RMS statistics summarize their thirty saved curves. Effective spread is the mean over the analysis grid of the pointwise population SD across ten members, calculated separately for each scheme; it measures ensemble spread rather than error relative to truth. The displayed min-max intervals are descriptive ranges, not confidence intervals.

All twelve scheme/type packages accepted 10/10 proposals, giving **120 accepted curves from 120 attempts**, with no observed-center order inversions. In this particular run, no type-specific rejection disrupted the shared candidate sequence. The implementation still rejects whole unordered curves and does not assume matched accepted indices when rejection histories differ. Every input covers the common analysis grid without extrapolation.

The lower false-arrow counts under L/E/J in Common driver only are compatible with their lower error magnitudes, but this experiment cannot separate magnitude from temporal allocation of the error. The direct-link exceptions also show that lower chronology error does not guarantee a higher complete-recovery count for a particular estimator and realization. A future equal-proposal-RMS control would answer a different, more specific question and was not included here.

## Exact controls and retained evidence

All **10080** paired X-Y native-score changes against the corresponding irregular accurate-age reference are usable. The **5760 pairwise changes are exactly zero**, as required by unchanged X/Y values, times and estimator seeds. CCM remains a pairwise comparator and supplies no evidence about conditioning on Z. The remaining 4320 score changes come from the three conditional methods; complete directional values remain in [paired_score_changes.csv](../data/processed/03_2_check_driver_error_models/paired_score_changes.csv).

The 180 C inputs and twelve references reproduce all **2688 imported X-Y directional records** exactly, including **2520 C-type records**, when compared with the 03_1 scientific fields. All type/scheme counts are retained in [cell_summary.csv](../data/processed/03_2_check_driver_error_models/cell_summary.csv); [xy_case_states.csv](../data/processed/03_2_check_driver_error_models/xy_case_states.csv) preserves individual correct, missed, extra and unresolved decisions. Native archives retain full three-variable outputs, although this report does not score the complete network.

The final dataset has **732 complete cases, 5124 native archives and 19032 successful directional rows**, including 10248 X-Y directional rows. Of the native archives, **3504 are exact source copies and 1620 are final newly computed conditional results**. The latter count describes the completed dataset, not the cumulative number of method invocations across the interrupted and resumed runs.

## Execution history and scope

The [first invocation](../data/logs/driver_error_models_20261110T133847522941Z/) was externally interrupted after the user accidentally closed the program. Its last supervisor sample recorded 612.8 seconds, while inspection before recovery found 297 complete cases: 192 full source imports and 105 new inputs. The periodic status file had last reported 294; completeness was established from committed case artifacts.

The [resumed invocation](../data/logs/driver_error_models_resume_20261110T135403003994Z/) retained those 297 commits and completed the remaining 435 inputs. Its script runtime was **2535.0 seconds (42.25 minutes)**. This is the recovery-stage runtime, not a total experiment runtime. The first invocation's sampled elapsed time is not an exact interruption duration, and uncommitted cases may repeat work during recovery; the two records should not be reported as a precisely measured uninterrupted benchmark or exactly 1620 cumulative calls.

Both invocations used 32 outer workers with one numerical thread per worker. During recovery, the sampled minimum available physical memory was **5.05 GiB**, and peak sampled system commit was **60.43 GiB**. These are system-wide resource measurements, not memory attributed solely to this experiment. The production scientific settings and computation were unchanged across recovery.

The [preflight](../data/benchmarks/driver_error_models_checks/20261110T133558073115Z/report.json) and final [artifact audit](../data/logs/driver_error_models_resume_20261110T135403003994Z/artifact_audit.json) passed. The final audit independently reconstructed all 720 uncertain inputs and twelve chronology packages, checked all 5124 native archives and 284016 serialized arrays, verified the 3504 source-native copies, and exactly rebuilt all three summary tables. No causal fit was repeated during auditing.

The [preservation audit](../data/logs/driver_error_models_resume_20261110T135403003994Z/preservation_report.json) confirmed unchanged hashes, sizes and modification times for all **28595 pre-existing scientific files**. The separate [resume-preservation audit](../data/logs/driver_error_models_resume_20261110T135403003994Z/resume_preservation_report.json) confirmed that the **3267 files belonging to the 297 pre-resume committed cases** remained unchanged. The [figure review](../data/logs/driver_error_models_resume_20261110T135403003994Z/figure_review.json) also passed the numerical and visual checks of the final PNG/PDF figures.

This report's interpretation is limited to the existing three realizations, three schemes per type, fixed linear systems and unequal prescribed error scales. It does not estimate population false-positive rates, a universal method ranking or independent process-versus-chronology variance components. Row comparisons also retain the original realization-specific estimator seeds.
