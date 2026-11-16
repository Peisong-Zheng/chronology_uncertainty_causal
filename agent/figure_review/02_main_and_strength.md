# Figure review for the causal pipeline and coupling sensitivity

Updated on 8 November 2026 after the user requested a visual rollback. The original research-figure layouts, colors, fonts, line widths, diagnostic annotations, titles and explanatory footers are restored from the pre-style archive. Lowercase panel labels and verified scientific-label corrections are retained. The detailed captions below remain available alongside the restored figure annotations.

Redraw entry point: `agent/redraw_main_figures.py`. It uses saved results only and explicitly disables estimator/generator calls. Current 02 results use CCM L=300; the retained 02_1 comparison uses its saved L=1119. No scientific fits are rerun.

## Provenance and scientific constraints

- The current 02 baseline and chronology outputs use CCM library length **300**, with 10 libraries for the original Tp=0 score, 100 cut-and-swap null records per direction and 100 paired random libraries for the 31-shift lag screen. The age comparison uses 10 members and uncertainty level 2.
- The existing 02_1 figure uses the **older 0.8 library fraction**. Its saved CCM diagnostics report **1119 states**, not 300. A style-only redraw must use these saved scores and state that setting. Running the currently edited source with recomputation enabled would silently create a different scientific comparison and is outside this plan.
- The main simulated cases have one process realization per form; the coupling branch has three process seeds. Chronology members are not independent climate histories.
- A direction can coexist with unresolved evidence. Direction support, complete resolved graph states, absent retained arrows and failed/missing results have different meanings and must remain distinct.
- Ages decrease from left to right. The known coupled-model delay is 20 years. A CCM curve's optimal cross-map lag is not an estimate that must equal this imposed delay.


## 02 / 00_model_truth

**Files and source:** `figures/02_check_causal_pipeline/00_model_truth.png` and `.pdf`; `02_check_causal_pipeline.py`, `plot_process_truth`.

**Restored figure:** Restored the original six-panel network/response figure, including delay, missing-link descriptions, per-step coefficients, response formulas and both explanatory footer lines. Lowercase labels a-f remain.

**Caption:** Known causal structure and deterministic response of the simulated bivariate systems. Panels (a–c) show the prescribed cross-variable links for the linear, nonlinear and independent cases, and (d–f) show the equilibrium mean response of Y when the standardized driver z is held constant. The coupled cases have X→Y only, a 20-year delay and long-run gain K=4, whereas the independent case has no cross-variable forcing. The response is K g(z), where g(z)=z for the linear case and g(z)=tanh(2z)/σ[tanh(2Z)] for the nonlinear case, with Z standard normal and z defined using the theoretical stationary standard deviation of X. These curves describe the held-driver mean, not a realized noisy trajectory; independent, zero-mean innovations remain present in the simulated records. K=c/(1−aY) relates the per-step coupling c to Y's memory and is not a common scale for causal-method scores.

## 02 / 01_baseline_original_uniform

**Files and source:** `figures/02_check_causal_pipeline/01_baseline_original_uniform.png` and `.pdf`; `plot_baseline_results`, called by `plot_checks`.

**Restored figure:** Restored the four method panels, blue/orange directions, numeric directional p values in every process tick label, full figure title and status legend/footer. Labels a-d remain; the title explicitly identifies CCM L=300.

**Caption:** Directional causal-method scores for the original uniformly sampled records. Panels (a–d) show PCMCI+, Liang–Kleeman information flow, convergent cross mapping (CCM) and Granger causality, respectively, for the linear, nonlinear and independent cases over 29–1 kyr BP on the native 20-year grid. Blue circles represent X→Y and orange squares Y→X; filled markers denote retained directions and open markers directions not retained, with distinct symbols reserved for unresolved evidence and failed lag checks. CCM points retain the original Tp=0 correlation averaged over 10 libraries of 300 states, even when the subsequent lag screen rejects the direction; its original significance test uses 100 cut-and-swap surrogates per direction. Each process form supplies one realization, and the coupled-model truth is X→Y with a 20-year delay. The vertical scales are method-specific and match the corresponding irregular-sampling baseline; their magnitudes are not interchangeable measures of physical coupling. Directional p values remain visible beneath each process category; CCM p0 refers to the original Tp=0 surrogate test.

## 02 / 02_baseline_irregular_accurate

**Files and source:** `figures/02_check_causal_pipeline/02_baseline_irregular_accurate.png` and `.pdf`; `plot_baseline_results`, called by `plot_checks`.

**Restored figure:** Restored the four method panels and full diagnostic labels, using the same score limits as the original-uniform baseline. Numeric p values, the accurate-age interpolation title and the status legend/footer are visible again. Labels a-d remain and the CCM library is correctly identified as L=300.

**Caption:** Directional causal-method scores after irregular sampling with accurate ages. Panels (a–d) show PCMCI+, Liang–Kleeman information flow, CCM and Granger causality for the same linear, nonlinear and independent process realizations used in the original-uniform baseline. Irregular observations have retained gaps of 20 or 40 years and are interpolated onto the common 20-year grid over 29–1 kyr BP, without age perturbations or measurement noise. Blue circles and orange squares identify X→Y and Y→X, respectively; filled and open markers indicate retained and non-retained directions, while unresolved evidence and failures remain separately identified. CCM uses libraries of 300 states and displays the original Tp=0 score rather than the largest lag-scanned correlation. Method-specific axis limits match the original-uniform comparison so that changes reflect sampling and interpolation rather than rescaling. Directional p values remain visible beneath each process category; CCM p0 refers to the original Tp=0 surrogate test.

## 02 / 03_irregular_sampling_effect

**Files and source:** `figures/02_check_causal_pipeline/03_irregular_sampling_effect.png` and `.pdf`; the reconstruction section of `plot_checks`.

**Restored figure:** Restored all six original record comparisons, orange observed points and reconstructions, full-window RMSE/SD and increment-variance diagnostics in the titles, and the explanation of the detail window. Labels a-f remain.

**Caption:** Effect of mild irregular sampling and interpolation on the simulated records. Panels (a–c) show X and (d–f) show Y for the linear, nonlinear and independent systems, respectively; the shared X driver is repeated to align the process comparisons. Black solid lines are the original uniform records, dashed lines are reconstructions on the 20-year grid, and orange circles are observations with accurate ages. The displayed window is 15.6–14.4 kyr BP, with older ages at left, whereas the reconstruction diagnostics use the entire 29–1 kyr BP analysis window. Full-window RMSE normalized by the original standard deviation is 3.4% for X and 2.6%, 3.0% and 6.8% for the three Y records; the corresponding increment-variance ratios are 0.965 for X and 0.960, 0.958 and 0.956 for Y. These changes isolate sampling and interpolation because ages are exact and measurement noise is zero.

## 02 / 04_age_ensemble

**Files and source:** `figures/02_check_causal_pipeline/04_age_ensemble.png` and `.pdf`; `plot_age_checks`, figure 04 section.

**Restored figure:** Restored the original eight chronology panels with blue intervals, gray member traces, black center errors, red dashed zero references, acceptance/member counts in every title, the level-specific title and explanatory footer. Labels a-h remain.

**Caption:** Accepted working chronology ensembles for the two records at uncertainty level 2. Panels (a,b), (c,d), (e,f) and (g,h) show the constant (C), linearly increasing (L), exponentially increasing (E) and localized-peak (J) age-uncertainty envelopes for X and Y, respectively. Thin gray curves illustrate the first five accepted chronologies, blue shaded bands span the pointwise 5th–95th percentiles of all 10 accepted members, and black curves show the observed chronology's error relative to true age; the red dashed zero line marks exact age. C/L/E have maximum proposal standard deviation 20 years, while J has three 2.8-kyr intervals with peak standard deviation 10 years and total defined-interval coverage of 30% of the analysis window. Members are proposed around the observed chronology, with correlated errors along sample order and independent clocks for X and Y, and whole curves violating stratigraphic order are rejected. All eight displayed clocks accepted 10 of 10 proposals in this run; these finite-member bands describe within-record chronology variation, not confidence intervals for a causal effect.

## 02 / 05_resampling

**Files and source:** `figures/02_check_causal_pipeline/05_resampling.png` and `.pdf`; `plot_age_checks`, figure 05 section.

**Restored figure:** Restored four full/detail panels with blue member traces, orange observed points, the accurate-age reference, full L-J settings, requested member count and true-delay title. Labels a-d remain.

**Caption:** Propagation of chronology variation into uniformly resampled records for the linear-system L–J example. Panels (a,b) show X with linearly increasing age uncertainty and (c,d) show Y with localized-peak uncertainty, over the full 29–1 kyr BP window and the fixed 15.6–14.4 kyr BP detail interval. Pale blue curves show the 10 chronology-member reconstructions, black curves the same irregular observations resampled using their accurate ages, and orange markers the observed proxy values at their uncertain observed ages. Every reconstruction covers the identical 1401-point, 20-year analysis grid without extrapolation; proxy values and dynamics are held fixed. At level 2 the X proposal standard deviation reaches 20 years, and Y has three 2.8-kyr intervals with 10-year peak standard deviation. The coupled-system delay remains 20 years, independently of the chronology-error amplitude.

## 02 / 06_direction_support

**Files and source:** `figures/02_check_causal_pipeline/06_direction_support.png` and `.pdf`; `plot_age_checks`, figure 06 section.

**Restored figure:** Restored the original 12-panel grouped bars, all ordered pair labels, per-panel accurate-age reference decisions, full titles, chronology settings and explanatory footer. High-unresolved flags and incomplete counts/hatching remain. Labels a-l remain.

**Caption:** Directional support across the 16 ordered chronology-type pairs at uncertainty level 2. Panels (a-d), (e-h) and (i-l) show the linear, nonlinear and independent systems, respectively, with PCMCI+, Liang-Kleeman information flow, CCM and Granger causality arranged from left to right. Blue and orange bars give X-to-Y and Y-to-X support; pair labels specify X's uncertainty type first and Y's second, so C-J and J-C are different comparisons. Each q is the fraction retaining that direction among members with usable decisions in both directions; the dashed line is the q=0.9 reporting threshold, which additionally requires at least 90% of the ten requested members to be usable. U marks high unresolved support and incomplete groups retain usable/requested counts and hatching; high q does not remove unresolved evidence or establish correctness across independent histories. C/L/E have maximum proposal SD 20 years, J has three 10-year peaks over intervals covering 30% of the window, and CCM uses L=300; the accurate-age reference retains only X-to-Y in coupled systems and neither direction in the independent system.

## 02 / 07_member_states

**Files and source:** `figures/02_check_causal_pipeline/07_member_states.png` and `.pdf`; `plot_age_checks`, figure 07 section; scientific state definitions in `toolbox/evaluation.py`, `summarize_ensemble`.

**Restored figure:** Restored the original 12-panel stacked bars, all six original colors, full unresolved-category legend, repeated process/method titles and the explanatory footer. Labels a-l remain; the requested-member denominator is unchanged.

**Caption:** Complete member outcomes under the 16 ordered level-2 chronology-type pairs. Panels (a-d), (e-h) and (i-l) show the linear, nonlinear and independent systems, respectively; columns show PCMCI+, Liang-Kleeman information flow, CCM and Granger causality. Every stacked bar partitions the ten requested chronology members into forward-only, reverse-only, both-direction, neither-direction, unresolved and failed/missing categories, identified by the original six-color legend. The four resolved arrow states require no unresolved evidence; a usable member with an unoriented/conflicting PCMCI+ contemporaneous link or inconclusive significant CCM lag screen belongs to the unresolved category even if an arrow is also retained. Failed or missing members stay in the denominator, while neither means a usable resolved result with no retained arrow. Pair labels give X's uncertainty type first; C/L/E have maximum proposal SD 20 years and J has three 10-year peaks with 30% defined-interval coverage. Only resolved forward-only outcomes match the coupled-system truth and only resolved neither outcomes match the independent truth; these are within-record fractions, not calibration across independent process histories.

## 02 / 08_age_diagnostics

**Files and source:** `figures/02_check_causal_pipeline/08_age_diagnostics.png` and `.pdf`; `plot_age_checks`, figure 08 section; saved `age_diagnostics.csv`.

**Restored figure:** Restored the original three diagnostics, original saved category order, detailed titles/y labels, in-panel legend and finite-member interpretation footer. Labels a-c remain.

**Caption:** Diagnostics of the accepted level-2 chronology ensembles. Panel (a) shows the fraction of whole-curve proposals accepted without violating age order; (b) shows the fraction of true sample ages within the pointwise 5th–95th percentile interval of the 10 accepted members, evaluated over all sample positions and separately over positions with nonzero input uncertainty; and (c) shows the median ratio of accepted-member standard deviation to proposal standard deviation over the nonzero-uncertainty positions. Categories retain the saved clock order shown on the x axis; each label identifies its record and uncertainty type. The dashed lines in (b) and (c) indicate nominal interval content 0.9 and an SD ratio of 1, respectively. All displayed clocks accepted 10 of 10 proposals, but interval inclusion and spread remain finite-member, within-record diagnostics. Correlated sample ages are not independent coverage trials, and zero-uncertainty sections can increase all-sample inclusion without improving calibration in uncertain sections.

## 02 / 09_ccm_baseline_lags

**Files and source:** `figures/02_check_causal_pipeline/09_ccm_baseline_lags.png` and `.pdf`; `plot_ccm_baseline_lags`.

**Restored figure:** Restored six lag panels with both solid blue/orange curves, library-SD bands, p0, screen outcomes and past/future fractions in the annotation boxes. Full titles and the interpretation footer are restored. Labels a-f remain; the title correctly states 100 matched libraries of 300 states.

**Caption:** CCM lag profiles for the two accurate-age baselines. Panels (a–c) show the original-uniform linear, nonlinear and independent records, and (d–f) show the corresponding irregular observations interpolated using their accurate ages. Blue solid curves represent reconstruction of X from Y's embedding (the X→Y test), and orange solid curves the reverse; lines are mean cross-map skill and bands are ±1 standard deviation over 100 paired random libraries of 300 states, not confidence intervals. The same library seed is used across the 31 shifts from −300 to +300 years on common anchors, with embedding dimension 3, backward delay of one 20-year step and temporal exclusion of 10 steps. The vertical zero line, pale −40 to 0 year embedding span and dotted optimal-mean-lag lines provide temporal references. Negative interior peaks can support temporal compatibility, whereas positive peaks reject it and zero, boundary, tied or unstable peaks remain inconclusive; the screening rule cannot upgrade an originally nonsignificant Tp=0 test. The original surrogate p value tests Tp=0 skill rather than the selected peak, and the optimal cross-map lag should not be equated with the model's imposed 20-year delay. Annotation boxes retain the original p0 values, lag-screen states and past/future fractions for both directions.

## 02 / 10_ccm_filter_support

**Files and source:** `figures/02_check_causal_pipeline/10_ccm_filter_support.png` and `.pdf`; `plot_ccm_filter_comparison`.

**Restored figure:** Restored the original six-panel grouped bar display: gray original support, blue/orange retained support and green unresolved fractions. Usable/requested counts are printed for every pair; full chronology settings and the paired-denominator footer are restored. Labels a-f remain.

**Caption:** Effect of the CCM lag screen on directional support at uncertainty level 2. Panels (a,b), (c,d) and (e,f) show the linear, nonlinear and independent systems, with X→Y on the left and Y→X on the right. Gray bars show the proportions of paired usable members significant in the original Tp=0 test; blue and orange bars show retained X?Y and Y?X support after lag screening, respectively. Green bars show significant responses with unresolved lag evidence, and each ordered X/Y pair is annotated with its usable/requested count. Both support proportions use the same member denominator, excluding lag failures from both quantities; all groups in this run have 10 paired usable members out of 10 requested. The dashed 0.9 line is a reporting threshold, not a calibrated causal probability, and the unresolved fraction is not interchangeable with a retained arrow. CCM uses libraries of 300 states, 100 original swap-null surrogates per direction and 100 paired libraries for the lag screen; chronology settings are C/L/E maximum proposal standard deviation 20 years and three J intervals with 10-year peaks and 30% defined-interval coverage.

## 02_1 / 01_strength_ranking

**Files and source:** `figures/02_1_check_coupling_strength/01_strength_ranking.png` and `.pdf`; `02_1_check_coupling_strength.py`, `plot_strength_scores`; saved `parameters.json`, `method_diagnostics.json` and `strength_results.csv`.

**Restored figure:** Restored all eight original panels with full method/process titles, native score labels, blue/orange paired-seed trajectories, the original seven-entry status legend and interpretation footer. Labels a-h remain. The displayed L=1119 is read from saved method diagnostics, independently of the current L=300 production setting.

**Caption:** Native causal-method scores under increasing prescribed coupling gain. Panels (a–d) show the linear systems and (e–h) the nonlinear systems, with PCMCI+, Liang–Kleeman information flow, CCM and Granger causality arranged from left to right. Each line connects K=1, 2 and 4 for one of three process seeds with paired stochastic forcing across gain settings; blue solid circles denote the true X→Y direction and orange dashed squares the reverse. Filled/open symbols distinguish retained directions from nonsignificant original tests, while separate symbols identify CCM lag incompatibility, unresolved evidence and lag-check failure without discarding a successful original score. The records are analyzed on the original 20-year grid over 29–1 kyr BP, with the imposed coupled-model delay fixed at 20 years. This retained coupling experiment uses the earlier CCM library fraction 0.8, corresponding to 1119 states, with 10 original libraries, 100 directional swap-null surrogates and 100 paired lag libraries. Scores have method-specific definitions and are not estimates of K; these three paired realizations test score ordering rather than calibrating effect sizes or population error rates.

## Rollback checks

- The 12 original figures and their panel counts are retained: 6, 4, 4, 6, 8, 4, 12, 12, 3, 6 and 6 for main figures 00-10; 8 for the strength figure.
- All nonplot calculation functions and scientific top-level settings match the pre-style archive. Changes relative to that archive are confined to panel labels and accurate CCM library-length display.
- The original source and figure archive remains at `data/archives/pre_figure_style_20261108T121553339241Z/figures_and_sources.zip`.
