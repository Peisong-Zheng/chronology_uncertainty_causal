# Uncertainty and sampling figure review

Rolled back on 8 November 2026 at the user's request. All eleven figures restore their original colors, line widths, annotations, numerical cell labels, titles, footers and panel arrangements from the pre-style archive, while retaining lowercase panel labels and corrected scientific provenance. The captions are retained and reconciled to the restored figures. The 02_3 matrices remain 3 by 4 panels and the 02_4 recovery view remains one 20 by 12 matrix. Redrawing uses saved results only; no causal fit, estimator setting or numerical result is changed.

## Scope and provenance

- The four figures in `figures/02_3_check_uncertainty_escalation/` use the completed Windows L=300 run: 1,920 input pairs, including 480 level-2 pairs imported from 02, with ten requested chronology members per process/type-pair/level.
- The seven figures in `figures/02_4_check_sampling_sensitivity/` use the retained 2026-11-04 run. Its saved `parameters.json` specifies `CCM_LIBRARY_FRACTION=0.8` and `fixed_gap_paired_sampling_v1`, not the L=300 settings now present in the script. The corresponding effective library length is 1,119 states for these inputs. A cosmetic redraw must use these saved results and report their original settings. Changing to L=300 requires a separately identified scientific run; it is not part of this figure revision.
- Ten sampling-mask groups in 02_4 are not ten chronology members and are not ten independent climate histories. Each process form uses one realization. Do not pool these two types of repetitions or interpret either as calibrated accuracy across process realizations.
- Original research-sized canvases and Matplotlib typography are restored. Lowercase panel labels remain. No figure or dense panel is split, and useful information remains in the graphic as well as being explained in its caption. The historical review observations below describe the pre-label figures; they are not instructions to remove their annotations again.

## 02_3 uncertainty escalation

### 00_uncertainty_envelopes

**Files and source:** `figures/02_3_check_uncertainty_escalation/00_uncertainty_envelopes.png` and `.pdf`; `02_3_check_uncertainty_escalation.py`, `plot_uncertainty_envelopes()`.

**Before revision.** The four panels lack letters. The long two-line supertitle, four legends and footer repeat methods. J uses line widths of 1.1–5 points to expose coincident curves; these visually imply different prominence and still make nested membership hard to follow. The zero-to-30 kyr domain is valid for the envelope but differs from the 29–1 kyr analysis window, which is not identified in the panels.

**Restored presentation.** Keep the original 2 by 2 C/L/E/J arrangement, all four legends, numerical SD and J coverage labels, two-line title and explanatory footer. Restore the original level colors and the nested J curves with their original stroke widths; the widest curves are drawn first so coincident levels remain identifiable. Retain panel labels a-d. No interval-membership strip replaces the curves.

**Current policy:** keep the original helpful in-figure explanations and numerical annotations; the caption complements them rather than replacing them.

**Caption.** Chronology-uncertainty envelopes used in the four-level experiment: (a) constant SD (C), (b) an SD increasing linearly with age (L), (c) an exponential increase with fixed curvature 2 (E), and (d) localized smooth uncertainty intervals (J). Levels 1–4 have C/L/E peak proposal SDs of 10, 20, 40 and 80 years, respectively; J retains a 10-year peak SD while increasing the number of 2.8-kyr intervals from 1 to 3, 5 and 7. In (d), colored curves show the nested subsets of localized intervals at each level; wider strokes are drawn first to expose coincident curves. The defined interval widths cover 10%, 30%, 50% and 70% of the 29-1 kyr analysis window. Colors identify levels, and the full envelope is shown over the 30-kyr retained process. These are proposal SDs before conditioning whole chronology curves on stratigraphic order, rather than hard error bounds or guaranteed realized spreads. The types are not matched in total uncertainty, and increasing J coverage intentionally increases total uncertainty.

### 01_x_to_y_support

**Files and source:** `figures/02_3_check_uncertainty_escalation/01_x_to_y_support.png` and `.pdf`; `02_3_check_uncertainty_escalation.py`, `plot_results()`, `q_x_to_y` branch.

**Before revision.** The 3×4 layout is scientifically useful but all 768 cells have text, often with a star and a tiny triangle. Long headers and four footer lines dominate the margins. Process, method and truth are repeated in every title. Two decimal places exceed the useful precision of ten-member support. Small green unresolved triangles are difficult to recognize on dark blue cells. There are no panel letters.

**Restored presentation.** Restore the original 3 by 4 process/method grid, all numerical q values, threshold stars, level tick labels and repeated ordered-pair labels. Retain the original blue 0-1 color scale, green unresolved triangles and orange missing/unusable squares, explicit truth subtitles, completion status, level definitions and denominator explanations. Keep row-major panel labels a-l. The existing settings header now explicitly reports the saved CCM library length (L=300 for this run).

**Current policy:** keep the original helpful in-figure explanations and numerical annotations; the caption complements them rather than replacing them.

**Caption.** Retention of the X → Y arrow under four levels of chronology uncertainty. Panels (a–d), (e–h) and (i–l) show the linear, nonlinear and independent process realizations, respectively; columns show PCMCI+, LKIF, CCM and Granger causality. Rows within each panel identify the ordered uncertainty-type assignment to X and Y, and columns identify the common level assigned to both records. Color gives q, the fraction retaining X → Y among chronology members with usable decisions in both directions; asterisks mark q ≥ 0.9 with at least 90% of the ten requested members usable. Cell numerals give q, green triangles identify groups containing unresolved evidence, and orange squares identify missing or unusable members. Small fractions give usable/requested counts for incomplete groups; a dash identifies a pending calculation, whereas NA denotes no usable decision. The coupled models have X → Y only, with a 20-year delay, whereas the independent model has no cross-variable coupling; CCM uses L=300 and the paired 100-library lag screen. These fractions describe sensitivity within one process realization per form, not a calibrated probability of causation or accuracy across independent histories.

### 02_y_to_x_support

**Files and source:** `figures/02_3_check_uncertainty_escalation/02_y_to_x_support.png` and `.pdf`; `02_3_check_uncertainty_escalation.py`, `plot_results()`, `q_y_to_x` branch.

**Before revision.** The same typography and annotation crowding affects the reverse-direction figure. Its interpretation particularly depends on distinguishing a low retained-arrow fraction from unresolved evidence, because neither is evidence for a reverse causal relationship. Retaining the same visual grammar as the forward figure is more important than maximizing numerical text.

**Restored presentation.** Restore exactly the same visual grammar as the forward-support figure: numerical support and threshold stars, full process/method/truth titles, type-pair and level ticks, original blue scale, unresolved/missing markers and explanatory footers. Retain row-major labels a-l and the saved L=300 header. Larger values still mean more retained reverse arrows.

**Current policy:** keep the original helpful in-figure explanations and numerical annotations; the caption complements them rather than replacing them.

**Caption.** Retention of the reverse-direction arrow Y → X under chronology uncertainty. Panels (a–d), (e–h) and (i–l) correspond to the linear, nonlinear and independent process realizations, with PCMCI+, LKIF, CCM and Granger causality arranged by column. Each cell represents one ordered X/Y uncertainty-type pair at a common level, and its color gives the retained-arrow fraction q among members with usable decisions in both directions. Numerical cell labels give q, and asterisks mark q ≥ 0.9 when at least nine of ten requested members are usable; green triangles and orange squares retain the separate indications of unresolved evidence and missing or unusable members, respectively; small fractions give usable/requested counts for incomplete groups. A dash indicates a pending calculation and NA indicates no usable decision. Y → X is absent from all three generating models, but these within-record fractions are not false-positive rates estimated across independent process realizations. The common color scale and chronology settings match the forward-direction figure, and CCM uses L=300 with 100 paired libraries for lag screening; a non-retained arrow must not be equated with a resolved absence when an unresolved marker is present.

### 03_correct_member_graph

**Files and source:** `figures/02_3_check_uncertainty_escalation/03_correct_member_graph.png` and `.pdf`; `02_3_check_uncertainty_escalation.py`, `plot_results()`, `fraction_correct_requested` branch.

**Before revision.** The repeated `n/10` labels make the figure resemble a numerical table and compete with unresolved markers. Its denominator differs from the two support figures and must remain explicit. The current colorbar conveys this correctly, but the many footnotes obscure rather than reinforce the distinction.

**Restored presentation.** Restore all correct/requested fractions in the original 3 by 4 blue heatmaps, with green unresolved triangles, orange missing/unusable squares, truth subtitles, level definitions, completion status and denominator explanation. Keep row-major labels a-l and the saved L=300 header. No support-threshold marker is added to this different quantity.

**Current policy:** keep the original helpful in-figure explanations and numerical annotations; the caption complements them rather than replacing them.

**Caption.** Complete recovery of the generating graph across chronology-uncertainty types and levels. Panels (a–d), (e–h) and (i–l) show the linear, nonlinear and independent process realizations; columns show PCMCI+, LKIF, CCM and Granger causality. Color is the number of members recovering the complete graph divided by the ten requested members for each ordered type pair and level. Recovery requires X → Y alone in the coupled systems and neither cross-variable arrow in the independent system, with no unresolved evidence; unusable, unresolved and missing members remain in the denominator and do not count as recovered. Cell labels give correct/requested counts. Green triangles and orange squares indicate groups containing unresolved and missing or unusable members, respectively, while dashes denote pending calculations. CCM uses L=300 with lag screening, and these recovery fractions characterize one realization per process form rather than accuracy across independent histories.

## 02_4 sampling sensitivity

### 01_sampling_age_density

**Files and source:** `figures/02_4_check_sampling_sensitivity/01_sampling_age_density.png` and `.pdf`; `02_4_check_sampling_sensitivity.py`, first figure in `plot_sampling_designs()`.

**Before revision.** The plot uses implementation strings such as `cv025` and `none` as titles, repeats age-axis labels, lacks panel letters and has a sentence-length title. Individual masks are so faint that their number and interpretation are unclear. The large upper whitespace in the CV=0.25 row arises from sharing the full density range with CV=0.75; this supports honest amplitude comparison and should not simply be cropped independently.

**Restored presentation.** Restore six original panels with full CV/gradient titles, age labels, pale blue individual masks, the original heavier blue mean, dashed 80-year-density reference and fixed-count supertitle. Keep original shared ranges and add labels a-f; the caption explains the unlabelled reference line and the mask repetitions.

**Current policy:** keep the original helpful in-figure explanations and numerical annotations; the caption complements them rather than replacing them.

**Caption.** Age distribution of observations under matched irregular-sampling designs. Panels (a–c) use gap CV=0.25 and panels (d–f) use CV=0.75; columns show no age gradient, preferentially sparse sampling at the old end, and preferentially sparse sampling at the young end. Every mask contains 351 observations over 29–1 kyr BP, including both endpoints, with 175 gaps of each of 60 and 100 years for CV=0.25 or 20 and 140 years for CV=0.75. Thin curves show the realized densities in 4-kyr bins for the A and B masks of ten paired groups, the heavy curve is their mean, and the dashed line is the 12.5 observations kyr⁻¹ reference associated with an 80-year interval. The age-gradient designs rearrange the same gaps with paired weighted orderings; they do not prescribe a linear density trend in BP age. These masks describe sampling variation, with exact ages and zero measurement noise, rather than chronology uncertainty or independent process realizations.

### 02_sampling_reconstruction

**Files and source:** `figures/02_4_check_sampling_sensitivity/02_sampling_reconstruction.png` and `.pdf`; `02_4_check_sampling_sensitivity.py`, second figure in `plot_sampling_designs()`.

**Before revision.** The AA and AB X traces are intentionally identical, but the duplication is not explained. The x ticks are crowded at 0.25-kyr intervals with two decimal places. All traces use the same orange, which loses the project's X/Y convention. Repeated labels and the long title reduce space, and the first-panel legend overlaps the data region.

**Restored presentation.** Restore original gray native curves, orange interpolated curves and filled orange observation points, full X/Y and sampling titles, first-panel legend, age/proxy axis labels, original zoom/ticks and scientific supertitle. Retain the original per-panel y ranges and labels a-f. The intentional AA/AB identity of X remains visible and explained.

**Current policy:** keep the original helpful in-figure explanations and numerical annotations; the caption complements them rather than replacing them.

**Caption.** Effect of sampling-mask alignment and age-dependent density on reconstruction of a fixed linear coupled process. Panels (a–c) show X and panels (d–f) show Y over 24–22 kyr BP; columns compare no-gradient AA, no-gradient AB and old-end-sparse AA designs from mask group 0 at CV=0.75. Gray curves are the original 20-year records, colored points are the sampled observations, and colored lines are their linear interpolations to the common 20-year analysis grid. AA uses the same sampling times for X and Y, whereas AB retains mask A for X and applies mask B to Y; therefore X is identical in the first two columns. The old-end-sparse design retains the same gap multiset and total observation count while rearranging gap order. Ages are exact and measurement noise is zero, so the displayed differences arise from sampling and interpolation rather than chronology errors.

### 03_density_native_scores

**Files and source:** `figures/02_4_check_sampling_sensitivity/03_density_native_scores.png` and `.pdf`; `02_4_check_sampling_sensitivity.py`, first figure in `plot_sampling_results()`.

**Before revision.** Twelve unlabelled panels carry long method/score titles. Both directions use identical circular marks and solid lines, leaving color as the sole distinction. Method-dependent vertical ranges are appropriate, but the legend and y-axis labels do not make the within-method interpretation sufficiently explicit.

**Restored presentation.** Restore full method/native-score titles, the original blue/orange direction colors with circular points and solid lines, regular-gap ticks, grids, legend and common-grid supertitle. Keep the original within-column shared y scales and row-major labels a-l. Each regular-gap point is still one evaluation, without invented replication intervals.

**Current policy:** keep the original helpful in-figure explanations and numerical annotations; the caption complements them rather than replacing them.

**Caption.** Native causal-method scores under synchronous regular sampling at 20-, 40-, 80- and 160-year intervals. Panels (a–d), (e–h) and (i–l) show the linear, nonlinear and independent process realizations, respectively; columns show PCMCI+, LKIF, CCM and Granger causality. Blue circles and solid lines represent X → Y, and orange circles and solid lines represent Y → X. Scores are the maximum GPDC statistic, maximum absolute normalized information flow, original zero-shift CCM correlation, and log residual-variance ratio, respectively; their magnitudes are comparable within a method, not across methods. Each case uses exact sampling ages and is linearly interpolated onto the same 1,401-point, 20-year grid over 29–1 kyr BP; the 20-year control reproduces the original uniform input. Each point is a single evaluation of a fixed process realization, not an ensemble mean, and the retained CCM results use a library fraction of 0.8, corresponding to 1,119 states.

### 04_correct_graph_requested

**Files and source:** `figures/02_4_check_sampling_sensitivity/04_correct_graph_requested.png` and `.pdf`; `02_4_check_sampling_sensitivity.py`, second figure in `plot_sampling_results()`.

**Before revision.** A single 20 by 12 matrix has long code-like scenario labels and two-line process/method columns. The continuous block layout hides the three-process grouping. Every cell contains counts, sometimes augmented by unresolved counts. A panel letter is not meaningful for this single-axis matrix. The original matrix layout is retained as requested.

**Restored presentation.** Restore the original single 20 by 12 viridis matrix, including every correct/requested count, any U/F/P annotations, numbered full scenario names, repeated process/method column labels, title, colorbar and denominator/status footer. No artificial panel letter is added to a single-axis figure. The requested count stays visible in every cell.

**Current policy:** keep the original helpful in-figure explanations and numerical annotations; the caption complements them rather than replacing them.

**Caption.** Complete recovery of the generating graph under accurate-age sampling controls. Columns are grouped into linear, nonlinear and independent processes, each with PCMCI+, LKIF, CCM and Granger causality; row groups show regular-density controls, matched irregular/asynchronous masks and matched age-gradient masks. Color gives correct graphs divided by requested cases: one case for each regular interval and ten mask groups for each irregular or gradient assignment, as indicated by the denominators in the cell labels. Cell numerators report the correct-graph counts. Recovery requires X to Y alone in the coupled models and neither cross-variable arrow in the independent model, with no unresolved evidence; unresolved, unusable and pending cases remain in the denominator and do not count as recovered. U, F and P give counts of unresolved, failed/unusable and pending cases where present, while gray cells identify unreported or pending conditions rather than zero recovery. AA/BB use synchronous masks, AB/BA exchange the two record-specific masks, and gradient controls rearrange the same gap multiset while remaining synchronous. CCM uses the retained L=1,119 results; all comparisons reuse one process realization per form and quantify sampling sensitivity rather than accuracy across independent histories.

### 05_paired_native_score_contrasts

**Files and source:** `figures/02_4_check_sampling_sensitivity/05_paired_native_score_contrasts.png` and `.pdf`; `02_4_check_sampling_sensitivity.py`, third figure in `plot_sampling_results()`.

**Before revision.** This is the most direct paired comparison but its formulas appear as a long footer, and its title repeats the legend meaning. Six multiline category labels are crowded beneath each of the four bottom panels. Identical marks distinguish directions only by color. Native-score units vary by method, so a shared interpretation of vertical magnitudes would be misleading.

**Restored presentation.** Restore individual blue/orange mask-group points, diamond medians, original thick IQR bars, the zero reference, category/CV labels, full title and formula/interpretation footer. Keep within-method y scales, incomplete usable/requested labels and row-major labels a-l. No p-values or confidence intervals are introduced.

**Current policy:** keep the original helpful in-figure explanations and numerical annotations; the caption complements them rather than replacing them.

**Caption.** Paired changes in native scores caused by asynchronous sampling and age-dependent density. Panels (a–d), (e–h) and (i–l) show the linear, nonlinear and independent process realizations, with PCMCI+, LKIF, CCM and Granger causality arranged by column; blue and orange marks denote X → Y and Y → X. Small points represent ten paired mask groups, diamonds give the median, and vertical segments show the 25th–75th percentiles, which describe variation among masks rather than confidence intervals. For each group, the asynchronous contrast is (s_AB + s_BA − s_AA − s_BB)/2; each gradient contrast is the mean of the synchronous AA and BB scores for that gradient minus the corresponding no-gradient synchronous mean, using the same group and gap multiset. Contrasts are shown separately for gap CV=0.25 and 0.75, and the horizontal line marks no score change. Only contrasts with all required finite successful component scores contribute, with incomplete usable/requested counts retained when applicable; CCM uses the original zero-shift correlation and the retained L=1,119 configuration. The score differences have method-specific meanings and describe paired sampling changes within one history per process form.

### 06_ccm_raw_screened_retention

**Files and source:** `figures/02_4_check_sampling_sensitivity/06_ccm_raw_screened_retention.png` and `.pdf`; `02_4_check_sampling_sensitivity.py`, fourth figure in `plot_sampling_results()`.

**Before revision.** Forty bars per panel obscure the simple raw-versus-screened comparison, while the all-zero independent panels look empty. Scenario numbers 1–20 require another figure to decode. The title can be misread as a usable-case denominator, although the actual plotted denominator is all requested cases. There are no panel labels, and large upper margins carry repeated prose.

**Restored presentation.** Restore paired bars in all six original panels: open blue bars for raw zero-shift retention and filled blue bars for retention after screening. Restore numeric scenario references 1-20, titles, legend, requested-case axis label and eligibility footer. Keep labels a-f. Zero-height bars represent computed zeros; numerical recovery/status tables and the caption distinguish these from missing results.

**Current policy:** keep the original helpful in-figure explanations and numerical annotations; the caption complements them rather than replacing them.

**Scenario key for the restored numeric ticks:** the numbers refer to the identically numbered rows of `04_correct_graph_requested`.

| Rows | Sampling conditions |
|---|---|
| 1-4 | Synchronous regular gaps of 20, 40, 80 and 160 years |
| 5-8 | No-gradient CV=0.25 masks in AA, BB, AB, BA order |
| 9-12 | No-gradient CV=0.75 masks in AA, BB, AB, BA order |
| 13-16 | CV=0.25 old-sparse AA/BB, then young-sparse AA/BB |
| 17-20 | CV=0.75 old-sparse AA/BB, then young-sparse AA/BB |

**Caption.** CCM arrow retention before and after lag screening across the accurate-age sampling controls. Panels (a,b), (c,d) and (e,f) show the linear, nonlinear and independent processes; left panels represent X → Y and right panels represent Y → X. Open bars show raw zero-shift retention and filled bars show retention after the lag screen, with each numerator restricted to cases having usable decisions in both directions at the corresponding stage and each denominator including all requested cases. Numeric ticks 1-20 refer to the identically numbered recovery-matrix rows: regular gaps (1-4), no-gradient irregular masks at CV=0.25 (5-8) and 0.75 (9-12), and synchronous old-/young-end-sparse masks at CV=0.25 (13-16) and 0.75 (17-20). Regular-density conditions have one requested case and each irregular or gradient condition ten mask groups; zero-height bars denote computed zero retention for complete cases; pending cases are marked P. The retained calculations use L=1,119, 100 cut-and-swap null surrogates and 100 paired libraries over shifts from −300 to +300 years; the screen can retain an originally significant arrow only when its lag evidence is temporally compatible. Lag screening changes directional retention but does not replace the original zero-shift correlation or p value, and these fractions describe sampling variation within one realization per process form.



### 07_irregular_sampling_native_scores

**Files and source:** `figures/02_4_check_sampling_sensitivity/07_irregular_sampling_native_scores.png` and `.pdf`; `02_4_check_sampling_sensitivity.py`, `plot_sampling_individual_scores()`.

**Before revision.** Eight category labels repeatedly append CV text; twelve panels lack letters. Median/IQR marks use the same color and shape as raw points, and two dashed reference lines are described by a neutral legend entry that does not make their directional assignment explicit. A long subtitle/footer repeats mark definitions and scientific cautions.

**Restored presentation.** Restore all eight AA/BB/AB/BA-by-CV tick labels, original blue/orange points, diamond medians, thick IQRs, direction-colored dashed 80-year reference lines, full method-score titles, supertitle and scientific footer. Keep original within-method y limits and row-major labels a-l; all successful finite scores remain included.

**Current policy:** keep the original helpful in-figure explanations and numerical annotations; the caption complements them rather than replacing them.

**Caption.** Native scores for individual irregular-sampling pairings with fixed total point count and gap distribution. Panels (a–d), (e–h) and (i–l) show the linear, nonlinear and independent process realizations, with PCMCI+, LKIF, CCM and Granger causality arranged by column. Within each CV block, AA and BB apply synchronous sampling masks to X and Y, whereas AB and BA assign different masks and exchange their record assignment; all masks contain 351 observations over 29–1 kyr BP. Blue and orange marks show X → Y and Y → X, small points represent the ten mask groups, diamonds show medians and vertical segments show interquartile ranges; direction-colored dashed lines give the regular 80-year reference scores. Summaries include all finite successful native scores, regardless of directional significance, and the CCM quantity is the original zero-shift correlation, retained independently of the later lag-screen decision. The interquartile ranges describe mask variation rather than sampling confidence intervals or uncertainty across independent histories; CCM uses the saved L=1,119 configuration, and score magnitudes should be compared within each method.

## Redraw and validation

Run `.venv/Scripts/python.exe agent/redraw_sensitivity_figures.py` from the project root to regenerate these eleven figures from the existing CSV tables. The runner imports plotting functions only, does not call either experiment's main routine or summary writers, and does not execute estimators. It reads saved envelope constants for 02_3 and retains the historical sampling results without changing the current L=300 production setting or bypassing a resume guard. All existing panel counts and scientific comparisons remain intact. The restored plots retain all original scientific annotations and use explicit panel labels; captions describe those restored marks and preserve the distinction between the L=300 uncertainty run and L=1,119 sampling run.
