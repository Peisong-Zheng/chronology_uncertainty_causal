# Memory sensitivity: results of the compact all-method experiment

Completed on 9 November 2026. The calculation is in [02_8_check_memory_sensitivity.py](../02_8_check_memory_sensitivity.py); the approved comparison is described in [the experiment plan](02_8_memory_sensitivity_plan.md).

## Main finding

The results support the concern that the excellent CCM result in 02_3 depends partly on the selected process-memory scales, especially the intrinsic memory of Y. Shortening X alone has a modest effect on complete CCM graph recovery. Shortening Y from 200 to 40 years causes a much larger loss, primarily because the lag screen can no longer resolve the reverse direction consistently. This is sensitivity of the current estimator settings to the generating dynamics and chronology perturbations; it does not establish a universal memory requirement for CCM.

The other methods respond differently. Granger causality continues to detect the true forward direction but increasingly adds a reverse arrow. PCMCI+ combines reverse arrows, unresolved evidence and forward losses. LKIF does not show a general decline in complete graph recovery as memory becomes shorter. A single detection or recovery fraction would conceal these differences.

## Comparison and baseline performance

| Label | X memory (years) | Y intrinsic memory (years) |
|---|---:|---:|
| Reference | 500 | 200 |
| Short X | 100 | 200 |
| Short Y | 500 | 40 |
| Short X and Y | 100 | 40 |

Each setting contains linear, nonlinear and independent systems, each with three process seeds. Every process realization has an original-uniform accurate-age reference, an irregular-sampling accurate-age reference and ten chronology members. The experiment therefore contains **432 input pairs**, each analyzed by PCMCI+, LKIF, CCM and GC. All conditions reuse the same sampling locations and the same saved level-4 C-C chronology arrays, with 80-year per-record proposal SD. The realized relative X/Y error is not bounded by 80 years; the ten member-specific relative-error RMS values range from 112.4 to 173.9 years.

The original and irregular accurate-age references give the same graph decisions throughout this experiment. Nearly all recover the prescribed graph in all three seeds. There are two exceptions, both in the Short Y condition:

- CCM recovers complete graphs in 2/3 seeds for each coupled form. In seed 20261301, the forward arrow is retained, but the reverse lag evidence is unresolved even with accurate ages. Its mean reverse optimum is +20 years; only 78% of library optima are in the future for the linear system and 84% for the nonlinear system, below the fixed 90% criterion. The remaining optima are at zero.
- LKIF recovers 2/3 nonlinear graphs. Seed 20261302 already contains a reverse arrow with accurate ages.

These cases remain in all summaries. In contrast, **Short X and Y has correct accurate-age graphs in all three seeds for every method and process form**. Its subsequent deterioration therefore supplies a particularly clear comparison of added chronology effects.

Even within the two Short Y seeds whose accurate CCM graphs are correct, only 9/20 chronology graphs remain correct for each coupled form. The baseline exception therefore cannot explain the whole Short Y decline. This breakdown explains the result; it does not replace the primary all-seed comparison or exclude the difficult seed.

## Complete graph recovery after chronology perturbation

Each entry is the number of correct graphs among **30 requested inputs**, formed by ten chronology members for each of three process seeds. A correct coupled graph contains only X to Y; a correct independent graph contains neither arrow. Neither may contain unresolved evidence. These pooled counts are descriptive, not 30 independent process histories. The saved tables retain results for each seed. At the user's request, the current main figures use pooled displays without seed-specific marker shapes.

| Process | Method | 500/200 | 100/200 | 500/40 | 100/40 |
|---|---|---:|---:|---:|---:|
| Linear | PCMCI+ | 9 | 10 | 4 | 7 |
| Linear | LKIF | 10 | 13 | 15 | 16 |
| Linear | CCM | 30 | 28 | 14 | 13 |
| Linear | GC | 24 | 14 | 11 | 5 |
| Nonlinear | PCMCI+ | 10 | 12 | 4 | 5 |
| Nonlinear | LKIF | 7 | 14 | 11 | 16 |
| Nonlinear | CCM | 30 | 28 | 13 | 11 |
| Nonlinear | GC | 24 | 14 | 10 | 5 |
| Independent | PCMCI+ | 30 | 28 | 30 | 30 |
| Independent | LKIF | 29 | 29 | 29 | 29 |
| Independent | CCM | 30 | 28 | 30 | 30 |
| Independent | GC | 23 | 25 | 28 | 29 |

For CCM, the seed-specific correct counts are:

| Process | 500/200 | 100/200 | 500/40 | 100/40 |
|---|---|---|---|---|
| Linear | 10, 10, 10 | 10, 9, 9 | 5, 5, 4 | 4, 4, 5 |
| Nonlinear | 10, 10, 10 | 10, 9, 9 | 5, 4, 4 | 5, 4, 2 |

Each triplet follows seeds 20261009, 20261301 and 20261302, with ten requested members per seed. The deterioration with short Y appears in all three histories; it is not produced by a single poor seed.

## Why CCM changes

All **240 coupled chronology inputs remain significant in both directions in the original Tp=0 CCM test**. The loss of complete graphs arises at the directional lag screen. This does not mean the original scores are unchanged: their magnitudes decrease, especially with both memories short, while remaining significant against the existing swap null.

Accurate-age lag optima move closer to zero when Y memory is shortened. Across-seed median forward optima are -140 years for the reference linear system and -160 years for the reference nonlinear system, compared with -60 years for both forms under Short Y or Short X and Y. The corresponding reverse median moves from +100 to +20 years. These reconstruction optima are distinct from the imposed 20-year equation delay. With perturbations larger than this separation, library optima more often straddle zero, making the direction unresolved.

The following counts separate retained arrows from unresolved evidence. **F/R/U** means forward retained / reverse retained / input with unresolved evidence, each out of 30. These categories overlap and must not be added.

| Process | Method | 500/200 F/R/U | 100/200 F/R/U | 500/40 F/R/U | 100/40 F/R/U |
|---|---|---|---|---|---|
| Linear | PCMCI+ | 23 / 11 / 9 | 17 / 7 / 6 | 29 / 20 / 13 | 20 / 17 / 3 |
| Linear | LKIF | 11 / 4 / 0 | 16 / 4 / 0 | 18 / 8 / 0 | 22 / 9 / 0 |
| Linear | CCM | 30 / 0 / 0 | 28 / 0 / 2 | 26 / 1 / 16 | 24 / 1 / 16 |
| Linear | GC | 30 / 6 / 0 | 30 / 16 / 0 | 30 / 19 / 0 | 30 / 25 / 0 |
| Nonlinear | PCMCI+ | 22 / 11 / 6 | 15 / 6 / 4 | 17 / 13 / 13 | 11 / 8 / 7 |
| Nonlinear | LKIF | 7 / 1 / 0 | 16 / 3 / 0 | 13 / 6 / 0 | 21 / 7 / 0 |
| Nonlinear | CCM | 30 / 0 / 0 | 29 / 0 / 2 | 23 / 0 / 17 | 21 / 1 / 19 |
| Nonlinear | GC | 30 / 6 / 0 | 30 / 16 / 0 | 30 / 20 / 0 | 30 / 25 / 0 |

Much of CCM's short-Y loss concerns unresolved reverse evidence, even when the forward arrow survives. With both memories short, the reverse direction is unresolved in 16/30 linear and 17/30 nonlinear cases; only one reverse arrow is retained in each form. Thus 13/30 and 11/30 correct graphs should not be described as detecting the true forward relation in only 13 or 11 cases.

The original forward CCM score also loses more under short memories. The median across the three seed-specific median changes, relative to each condition's own irregular accurate-age reference, is:

| Process | 500/200 | 100/200 | 500/40 | 100/40 |
|---|---:|---:|---:|---:|
| Linear: change in Tp=0 rho | -0.030 | -0.132 | -0.109 | -0.269 |
| Nonlinear: change in Tp=0 rho | -0.013 | -0.078 | -0.107 | -0.259 |

These are changes in the original CCM statistic, not differences in the maximum scanned correlation. This table retains the median of the three per-seed medians; the revised Figure 04 instead shows the pooled distribution and median of all thirty paired differences, so its summary need not equal this table exactly. Score loss and directional ambiguity are related diagnostics but are not interchangeable.

## What the other methods show

**GC:** All coupled chronology inputs retain X to Y. Its decline in complete graph recovery is entirely due to additional Y-to-X arrows: 6/30 in each reference form, increasing to 25/30 with both memories short. In the independent system, however, shorter memory improves correct no-arrow recovery from 23/30 to 29/30. The response to memory therefore depends on the underlying causal structure.

**PCMCI+:** Shorter X slightly increases complete graph recovery while reducing forward retention. For example, linear forward retention falls from 23/30 to 17/30, but reverse and unresolved outcomes also decline. With Short Y, linear forward retention is 29/30, yet only 4/30 complete graphs are correct because reverse arrows and unresolved evidence are common. For the nonlinear system with both memories short, forward retention is also reduced to 11/30. Complete graph recovery alone does not distinguish these failure modes.

**LKIF:** The reference chronology results mainly lose the forward arrow. Shorter memories generally increase forward retention and complete graph recovery in this experiment, although reverse arrows also become more frequent. This is not evidence that chronology perturbation benefits LKIF: all four coupled-memory conditions have negative median forward native-score changes relative to their own accurate-age references. Baseline native scores and the existing analytical decisions differ between the generating systems. The nonlinear Short Y baseline exception must also remain visible.

**Independent controls:** PCMCI+ retains false arrows in two Short X cases; its other memory groups are correct in all 30 cases. LKIF retains one reverse arrow per memory group. CCM retains no arrows in any independent chronology group; its two Short X incorrect graphs have unresolved evidence rather than retained false arrows. GC has seven, five, two and one incorrect independent graphs in memory-table order. These small conditional counts do not calibrate population false-positive rates.

## Interpretation limits and the next discriminating comparison

Changing intrinsic memory while keeping K and unforced background SD fixed requires changing the AR coefficients, per-step coupling and innovation SDs. It does not preserve the driven variance of Y. Median forced/background Y SD ratios across seeds are:

| Process | 500/200 | 100/200 | 500/40 | 100/40 |
|---|---:|---:|---:|---:|
| Linear | 1.939 | 1.333 | 2.310 | 1.991 |
| Nonlinear | 1.825 | 1.280 | 2.228 | 1.959 |

The high ratio under Short Y makes a simple reduction of driven amplitude an insufficient explanation for its CCM decline. Nevertheless, this experiment does not isolate memory from every change in spectrum, amplitude or effective estimator resolution. Y's observed autocorrelation also includes memory transmitted from X; the specified Y timescale describes its unforced autoregressive decay.

Short Y also has higher accurate-age forward CCM skill in the linear system: rho ranges from 0.748 to 0.858 across seeds, compared with 0.599 to 0.706 in the reference condition. The loss of directional robustness therefore occurs despite strong baseline cross-map skill.

All methods keep their original configurations. In particular, CCM uses the same 40-year embedding-history span, 200-year temporal exclusion, L=300 and +/-300-year scan. Their ratios to process memory change. Three process histories share one set of ten chronology perturbations, and only one uncertainty type pair and level are tested. The evidence therefore supports sensitivity of the current pipeline in this controlled branch, not a general ranking of optimally tuned algorithms.

The most focused next comparison would add the same design at a smaller chronology amplitude, concentrating on 500/200 versus 100/40. This would test whether the loss recedes when age error becomes small relative to the memory and lag-profile separation. It is a proposed follow-up, **not part of this completed run**; increasing member count alone would not answer that scale question.

## Execution, verification and saved evidence

The run completed all 432 input pairs and 3456 directional rows, with no recorded fit failures or pending cases. It used 32 outer Windows workers with one numerical thread each, native CPU GPDC and the validated CCM neighbor-reuse implementation. Calculation and figure export took **1458.999 seconds (24.317 minutes)**; minimum sampled available physical memory was **8.90 GiB**. No worker benchmark or GPU estimator migration was performed.

The preflight verified exact reproduction of all 36 retained first-seed reference inputs, paired innovations and clocks, failure/resume semantics, and complete native-output agreement for two full-setting cases under direct and spawned execution. Its [report](../data/benchmarks/memory_checks/20261109T114552073978Z/report.json) records the checks. Production evidence, source snapshot and resource measurements are under [the run directory](../data/logs/memory_sensitivity_20261109T114908564938Z/); see its `artifact_audit.json`, `preservation_report.json` and `figure_review.json` for final acceptance. The pre-existing main-study scientific files and results remain unchanged.

All 2880 member-minus-baseline score contrasts and 2160 cross-memory contrasts are usable. [Saved result tables](../data/processed/02_8_check_memory_sensitivity/) include per-seed summaries, complete directional rows, process and clock diagnostics, 1728 native archives and all CCM lag curves. The [run status](../data/processed/02_8_check_memory_sensitivity/run_status.json) records completion.

Four research figures retain panel labels, numerical annotations, legends and explanatory text, with PNG/PDF companions:

1. [Complete graph recovery](../figures/02_8_check_memory_sensitivity/01_graph_recovery.png).
2. [CCM raw and screened decision states](../figures/02_8_check_memory_sensitivity/02_ccm_raw_screened_retention.png).
3. [CCM lag profiles and peak-sign diagnostics](../figures/02_8_check_memory_sensitivity/03_ccm_lag_profiles.png).
4. [Paired native-score changes](../figures/02_8_check_memory_sensitivity/04_paired_native_score_changes.png).

The readability revision keeps all results and the same four figure groups. Figure 01 now links each accurate-age reference to its perturbed result and labels correct counts directly. Figure 02 separates retained, raw-nonsignificant, lag-rejected, unresolved and unavailable directional states using requested-input denominators. Figure 03 supplements the unchanged curves with saved per-input library-sign fractions. Figure 04 makes small independent-system changes visible through an explicitly expanded y-scale; it does not change score units.

[Detailed English captions](figure_review/05_memory_sensitivity.md) define all aggregation and interval conventions. To regenerate the branch's figures and summaries from the matching committed run, set `RUN_CAUSAL=False` in the root script; this retains the scientific configuration checks and does not fit causal estimators.
