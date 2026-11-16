# Memory sensitivity figure captions

These figures belong to `02_8_check_memory_sensitivity.py`. Their reading order was revised at the user's request on 9 November 2026: complete graph recovery, CCM raw/screened decision states, CCM lag profiles and peak-sign diagnostics, and paired native-score changes. Seed-specific circle/square/triangle symbols were removed from the first, second and fourth figures. All three process histories remain in the calculations and saved tables; no causal estimates were recomputed. The subsequent readability revision preserves all panel groups and adds explicit decision states and existing peak-sign diagnostics. Direction colors are now consistent: blue for X to Y and orange for Y to X; the all-method recovery overview uses neutral baseline symbols and green perturbed results.

The experiment uses four memory combinations, three process forms and three process seeds. Each combination has an original-uniform accurate-age reference, an irregular-sampling accurate-age reference and ten chronology members, giving 432 input pairs analyzed by PCMCI+, LKIF, CCM and Granger causality. The chronology condition reuses the saved level-4 C-C clocks with an 80-year proposal SD; actual chronology errors are recorded separately and are not bounded by 80 years.

Memory labels are ordered pairs (X memory, Y intrinsic memory), in years: (500, 200), (100, 200), (500, 40) and (100, 40). They are not ratios. Coupled cases retain K=4 and a 20-year equation delay; per-step coefficients and innovations preserve the specified long-run gain and unforced background standard deviations. Sampling locations and chronology arrays are shared across conditions. Three process histories with ten shared chronology members do not constitute thirty independent histories.

The figures retain numerical counts, explanatory text and lowercase panel labels. Native method scores have different definitions and scales. Pooled proportions and score distributions are descriptive; no displayed interval is a confidence interval for population performance. Seed-specific outputs and baseline exceptions remain documented in `agent/02_8_memory_sensitivity_results.md` and the saved per-seed tables.

All four revised PNGs were visually reviewed against these captions. Renderer checks found no clipped text; the eight lag-profile panel groups retain their original curve and shading coordinates exactly. The figure-only audit is `data/benchmarks/memory_figure_checks/20261109T180042688693Z/report.json`.

## 01_graph_recovery

**Files:** `figures/02_8_check_memory_sensitivity/01_graph_recovery.png` and `.pdf`.

**Panel mapping:** Panels (a-d), (e-h) and (i-l) show linear, nonlinear and independent systems. Columns are PCMCI+, LKIF, CCM and Granger causality.

**Caption.** Complete causal graph recovery before and after chronology perturbation under four memory combinations. Open gray circles show original-uniform accurate-age references, dark filled circles show irregular-sampling accurate-age references, and green circles show chronology-perturbed inputs; thin gray connectors highlight each accurate-to-perturbed comparison. Fractions pool correct and requested counts across three histories, with three inputs per accurate-age condition and thirty per chronology condition. Green labels give chronology correct/requested counts, while tick annotations retain original (O) and irregular (I) accurate-age counts and unresolved/failed/pending chronology-case counts (U/F/P). Correct coupled graphs contain only X to Y and correct independent graphs contain neither arrow, with no unresolved evidence; failures and pending inputs remain in the denominator. These descriptive complete-graph fractions are not per-history directional support q.

The (500, 40) accurate-age exceptions remain visible: CCM is correct in two of three histories for each coupled form, and LKIF in two of three nonlinear histories. Reference, Short X, Short Y and Short X + Y identify the four fixed memory combinations; connectors do not imply a continuous parameter sweep.

## 02_ccm_raw_screened_retention

**Files:** `figures/02_8_check_memory_sensitivity/02_ccm_raw_screened_retention.png` and `.pdf`.

**Panel mapping:** Rows show linear, nonlinear and independent systems. Columns show the (500, 200), (100, 200), (500, 40) and (100, 40)-year memory pairs; panels run from (a) to (l).

**Caption.** CCM directional outcomes before and after the lag screen, displayed as mutually exclusive states for every requested input. Within each panel, paired bars compare accurate ages (A, three inputs) with chronology perturbations (C, thirty inputs) for raw X-to-Y, screened X-to-Y, raw Y-to-X and screened Y-to-X outcomes. Blue and orange indicate significant raw tests or retained screened arrows in the respective directions; light gray indicates an originally nonsignificant test, medium gray a significant direction rejected by the lag screen, hatched gold unresolved evidence after the screen, and hatched dark red an unavailable outcome. Segment heights divide counts by all requested inputs in that directional group, including unavailable outcomes, and segment labels retain the counts; these state fractions do not use the paired-valid denominator of q. Original evidence remains available if a subsequent lag scan fails, while the corresponding screened result is marked unavailable; a valid nonsignificant raw test is not reclassified as unresolved merely because its diagnostic lag is unstable. Case-level U/F/P notes count each input once across both directions and are not the sum of directional unresolved segments.

This display separates unresolved evidence from an explicit lag rejection, which the previous retained-only dots could not show. Neither a nonsignificant test nor a rejected candidate establishes absence as a general causal fact. No cross-history q threshold is applied to the bars. All scientific decisions, including the existing per-input 90% library-sign criterion, remain unchanged.

## 03_ccm_lag_profiles

**Files:** `figures/02_8_check_memory_sensitivity/03_ccm_lag_profiles.png` and `.pdf`.

**Panel mapping:** Eight panel groups retain four memory rows and two process columns. Rows are (500, 200), (100, 200), (500, 40) and (100, 40) years; columns are linear and nonlinear. Each lettered group contains its original lag profile, a separate diagnostic strip, and numerical notes.

**Caption.** CCM lag profiles and per-input library-peak-sign diagnostics under different process memories. Blue and orange identify X to Y and Y to X; dashed curves show irregular accurate-age references and solid curves show chronology-perturbed inputs. Curve aggregation is unchanged: each input averages its hundred paired-library correlations at every lag, chronology curves take the pointwise median within each history, and displayed lines take the pointwise median across three histories; shading is the range of those three chronology history-summary curves, not a confidence interval. The separate 0-1 strip plots each input's saved fraction of hundred library optima that are in the past for X to Y or in the future for Y to X, using open circles for three accurate-age inputs and filled circles for thirty chronology inputs; vertical jitter only separates points visually. The 0.9 strip guide applies to each input's own hundred libraries, not to pooled library results or a cross-history retention fraction; final decisions also require valid complete curves, the matching interior mean-curve peak and the original CCM significance test. Notes retain available/requested profiles and separately summarized peak-shift medians; the zero and minus-20-year profile guides identify temporal zero and the imposed equation delay, which need not coincide with a reconstruction optimum.

The two axes have different meanings: profile x is target shift in years, whereas strip x is a fraction of library peaks. A high X-to-Y past fraction supports temporal compatibility; a high Y-to-X future fraction supports rejecting that reverse candidate, subject to the other screen conditions. The gray minus-20-year line is not a required peak location. Native strength remains original Tp=0 rho; the peak-sign strip and lag-scan maxima do not replace it.

## 04_paired_native_score_changes

**Files:** `figures/02_8_check_memory_sensitivity/04_paired_native_score_changes.png` and `.pdf`.

**Panel mapping:** The process rows and method columns match Figure 01. The independent-system row explicitly uses an expanded vertical scale.

**Caption.** Change in native method score caused by chronology perturbation relative to each member's own memory- and history-matched irregular-sampling accurate-age reference. Blue and orange marks show X to Y and Y to X; small points are individual usable paired differences, and large points and bars show their pooled median and interquartile range. Tick annotations retain usable/requested counts; a difference requires both its member score and matched reference score to be finite and successful. Linear and nonlinear panels share vertical limits within each method, whereas the clearly labeled independent row expands its scale to reveal small changes near zero; units are unchanged and differ between methods. The thirty requested differences come from three histories and ten shared clock pairs, so their IQR combines within-history chronology variation and between-history variation and is not a confidence interval across thirty independent histories. CCM retains original Tp=0 rho, including successful raw scores with a later unresolved or failed lag screen, rather than using the maximum scanned rho.

Differences are calculated before pooling. The displayed median is the median of all usable paired differences, not the median of three history medians or a difference between pooled original scores. Negative differences indicate lower native scores, not an automatic increase in graph error; expanded independent panels should be read using their own numerical y-axis values.
