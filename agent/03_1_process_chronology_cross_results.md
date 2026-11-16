# Process-realization and Z-chronology cross: X-Y results

Completed on 10 November 2026. See the [design](03_1_process_chronology_cross_plan.md), [script](../03_1_check_process_chronology_cross.py) and [figure captions](figure_review/07_process_chronology_cross.md). This report assesses the X-Y relation only.

## What the cross resolves

The original common-driver-only PCMCI+ failure at R2/A2 mostly follows **realization R2 across the three saved age schemes**. Correct X-Y counts for R2 are 3/10, 0/10 and 1/10 under A1, A2 and A3, respectively. R1 gives 10/10, 9/10 and 10/10; R3 gives 10/10 throughout. Thus, A2 is not a universally damaging chronology: the original diagonal paired a particularly susceptible realization with its worst of these three tested schemes.

This row pattern is method-specific. Conditional LKIF correctly retains neither X-Y arrow in all 90 common-driver-only inputs. GC instead shows a strong R1/A2 exception: nine of ten inputs acquire an incorrect Y-to-X arrow, compared with two under R1/A1 and one under R1/A3. The effect of a chronology assignment therefore depends on the saved history and method, rather than on a single universally adverse age scheme.

When the direct X-to-Y link is present, PCMCI+ retains that forward arrow in **all 90** crossed inputs. Its lower complete-recovery counts arise from extra reverse arrows or unresolved evidence. The R2 irregular accurate-age reference was already unresolved, so that row cannot be described as entirely chronology-induced deterioration.

## Fixed comparison and how to read the tables

R1, R2 and R3 are the existing process realizations with seeds 20261009, 20261301 and 20261302. A1, A2 and A3 are the existing complete Z age schemes associated with those seeds: each contains one uncertain observed chronology and ten ordered members drawn around it. Every scheme is applied to every realization, with the same cross in both systems. Memories remain 500 years for Z and 200 years for X/Y; all dynamics, sampling and estimator settings remain those of 03. Only Z timing changes within a row. X/Y values and ages are fixed, and each row retains its original estimator seed.

Each crossed cell has ten requested members. Exact recovery requires neither X-Y arrow in the common-driver-only system, or only X-to-Y in the direct-link system, with no unresolved evidence. Failed and pending cases remain in the denominator; **none occurred**. Directional errors and unresolved evidence are separately reported and may coexist. The accurate column is each realization's own single irregular-sampling accurate-age reference, not ten additional independent observations.

### Common driver only

Each crossed entry reads **Correct / Extra X-to-Y / Extra Y-to-X / Unresolved**, all counts out of ten. Correct and unresolved are mutually exclusive; the other counts should not generally be added. Every conditional accurate-age reference is correct for this system, under both original-uniform and irregular sampling.

| Method | Realization | Accurate, correct/1 | A1 | A2 | A3 |
|---|---|---:|---|---|---|
| PCMCI+ | R1 | 1/1 | 10 / 0 / 0 / 0 | 9 / 0 / 1 / 0 | 10 / 0 / 0 / 0 |
| PCMCI+ | R2 | 1/1 | 3 / 2 / 2 / 3 | 0 / 1 / 1 / 8 | 1 / 2 / 2 / 5 |
| PCMCI+ | R3 | 1/1 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 |
| LKIF | R1 | 1/1 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 |
| LKIF | R2 | 1/1 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 |
| LKIF | R3 | 1/1 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 |
| GC | R1 | 1/1 | 8 / 0 / 2 / 0 | 1 / 0 / 9 / 0 | 9 / 0 / 1 / 0 |
| GC | R2 | 1/1 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 |
| GC | R3 | 1/1 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 |

PCMCI+'s descriptive row totals are 29/30, 4/30 and 30/30 correct, whereas its column totals are 23/30, 19/30 and 21/30. Sixteen of the thirty R2 inputs have unresolved evidence. This locates the original anomaly within the tested cross but does not identify its detailed statistical mechanism or estimate a process-versus-clock variance fraction. GC's R1/A2 reverse-arrow support is 0.9 despite the absence of any true X-Y link: high within-scheme support is not sufficient evidence of causal correctness.

### Common driver plus direct X-to-Y

Each crossed entry reads **Correct / Missed X-to-Y / Extra Y-to-X / Unresolved**, out of ten. A missed arrow is a usable decision that does not retain the true forward arrow; it is not a failed computation. Counts of missed/extra arrows and unresolved evidence may overlap.

| Method | Realization | Accurate, correct/1 | A1 | A2 | A3 |
|---|---|---|---|---|---|
| PCMCI+ | R1 | 1/1 | 7 / 0 / 0 / 3 | 8 / 0 / 2 / 0 | 10 / 0 / 0 / 0 |
| PCMCI+ | R2 | 0/1: unresolved; X-to-Y retained | 5 / 0 / 1 / 4 | 6 / 0 / 1 / 3 | 7 / 0 / 0 / 3 |
| PCMCI+ | R3 | 1/1 | 3 / 0 / 1 / 6 | 5 / 0 / 2 / 3 | 9 / 0 / 0 / 1 |
| LKIF | R1 | 0/1: X-to-Y missed | 10 / 0 / 0 / 0 | 9 / 1 / 0 / 0 | 9 / 1 / 0 / 0 |
| LKIF | R2 | 1/1 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 |
| LKIF | R3 | 1/1 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 |
| GC | R1 | 1/1 | 10 / 0 / 0 / 0 | 8 / 0 / 2 / 0 | 10 / 0 / 0 / 0 |
| GC | R2 | 1/1 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 |
| GC | R3 | 1/1 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 | 10 / 0 / 0 / 0 |

PCMCI+'s row totals are 25/30, 18/30 and 17/30 correct; its column totals are 15/30, 19/30 and 26/30. A3 gives the highest correct count in each row for this system. These are descriptive comparisons of the three particular schemes, not evidence that their generating uncertainty model has different parameters.

Two accurate-age exceptions matter. PCMCI+ recovers the direct-link R2 relation on the original uniform record but already has unresolved evidence after accurate-age sampling/interpolation. LKIF misses X-to-Y in R1 under both accurate references. Its subsequent 28/30 forward detections across uncertain Z schemes do not establish that age error improves causal inference: changing Z timing changes the adjustment itself. GC retains the forward arrow in all 90 crossed inputs; its only incomplete X-Y decisions are the two added reverse arrows at R1/A2.

## Chronology and exact controls

All schemes use the same C-type 40-year proposal SD. Their realized Z errors relative to truth differ because each includes an uncertain observed center and its conditional members. Median [minimum, maximum] RMS errors over the analysis window are:

| Scheme | Realized Z RMS error, years |
|---|---:|
| A1 | 49.5 [45.9, 55.4] |
| A2 | 57.3 [46.4, 71.3] |
| A3 | 52.8 [44.3, 62.0] |

These are ranges across ten saved curves, not confidence intervals. X/Y clock errors are zero, so relative X-Z and Y-Z RMS errors equal the listed Z errors. A2 has the largest median here, but the schemes also differ in the temporal pattern of their errors; this cross does not isolate error magnitude as a causal factor.

All **2520** paired X-Y native-score changes against each realization's own irregular accurate-age reference are usable. The **1440 pairwise changes are exactly zero**, as required by unchanged X/Y inputs and direct reuse of their native results. This includes CCM, which is not a conditional estimator in this comparison. The remaining 1080 changes are from the three conditional methods; complete directional scores remain in [paired_score_changes.csv](../data/processed/03_1_check_process_chronology_cross/paired_score_changes.csv).

The twelve accurate references and sixty original diagonal cases were imported. Their **1008 X-Y directional records**, including **840 diagonal records**, exactly match the original 03 scientific fields. Thus, the off-diagonal results extend the existing comparison without changing its baselines or original diagonal.

## Interpretation and execution record

This cross shows that the original PCMCI+ common-driver-only anomaly is not confined to its original Z chronology assignment, while GC has a different realization-scheme interaction. It does not establish a universal method ranking, a population error rate, or independent variance components. There are three existing process histories and three existing Z schemes, not nine independent histories. Members reuse the same observed center within a scheme; schemes are reused across rows and paired systems. Row comparisons also retain their original estimator seeds. No conclusion about the complete three-variable network is drawn here.

The calculation completed **192/192 inputs**, with **1344 native archives: 360 new conditional fits and 984 exact archive copies**. All **2688 saved X-Y directional rows** have status `ok`. Runtime was **799.4 seconds (13.32 minutes)** with 32 workers; minimum sampled available physical memory was **11.26 GiB**, and peak sampled system commit was **26.83 GiB**.

The [preflight report](../data/benchmarks/process_chronology_cross_checks/20261110T123801313973Z/report.json) passed independent reconstruction of all 180 crossed Z inputs, byte-identical import of all 72 reference/diagonal inputs, complete clock preservation, same-input direct-versus-Windows-spawn native equality, failure retention and resume guards. [Run logs, source snapshot, resources and final audit files](../data/logs/process_chronology_cross_20261110T124145323892Z/) are stored separately from the [completed scientific outputs](../data/processed/03_1_check_process_chronology_cross/).

The final [artifact audit](../data/logs/process_chronology_cross_20261110T124145323892Z/artifact_audit.json) passed all 192 committed cases, 4992 full directional rows, 1344 native archives, exact source copies and reconstruction of the three summary tables. The [preservation audit](../data/logs/process_chronology_cross_20261110T124145323892Z/preservation_report.json) confirmed unchanged sizes, modification times and SHA256 hashes for all 26,464 pre-existing scientific files (1,111,881,237 bytes). The [figure review](../data/logs/process_chronology_cross_20261110T124145323892Z/figure_review.json) passed all 54 crossed cells, 18 accurate-age labels, and visual checks of both PNG and PDF versions. No scientific fits were repeated during these audits.

The two figures are [the crossed design and actual age assignments](../figures/03_1_check_process_chronology_cross/00_cross_design.png) and [all six X-Y recovery matrices](../figures/03_1_check_process_chronology_cross/01_xy_recovery_matrix.png), with PDF companions and [English captions](figure_review/07_process_chronology_cross.md).
