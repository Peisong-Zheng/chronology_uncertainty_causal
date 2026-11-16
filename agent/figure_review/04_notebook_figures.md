# Notebook figure styles and captions

Rollback updated on 8 November 2026 at the user's request. This document covers 17 figure definitions. The seven existing figures from 02_5, 02_6 and 02_7 have been redrawn from saved arrays/tables and visually checked, with PNG and PDF outputs. The ten figures defined by 01 and 02_2 have updated plotting cells but remain unrendered: their active figure/data directories are absent. No estimator, simulation, chronology generation or benchmark was rerun. All original panel counts and layouts are retained, including the six-panel 02_7 numerical-agreement figure.

Cell numbers are zero-based JSON cell indices. Panel letters follow reading order; standalone single-axis figures omit a redundant letter. Captions for the seven rendered figures match the restored panels; the ten code-only captions still need a visual check when those figures are eventually generated. Root notebooks remain free of generated outputs.

The pre-style colors, fonts, line widths, titles, explanatory annotations and axis scales have been restored. Lowercase panel labels and the written captions are retained. Scientific corrections remain: library lengths keep their saved or notebook-specific provenance, GPDC residuals are labelled in standardized-input units, and exact-zero absolute-error panels have nonnegative axes. PNGs use the original 180 dpi; PDF companions remain available. No global figure-style override is applied.

## 01 Inspect simulated records

These six figures are an exploratory generation tutorial, not a rendering of the current main chronology ensemble. In particular, the notebook uses a separate process seed, a 30–0 kyr window, C/L/E scales of 50 years, J peaks of 100 years and a J interval-width sweep. The main study instead uses its saved four-level configurations and the 29–1 kyr analysis window. The tutorial's existing scientific prose and settings have not been revised in this style-only change; its historical framing should not be read as the current production status.

### 01.1 True processes

**Files and source:** `figures/01_inspect_simulated_records/01_true_processes.png` and `.pdf`; notebook cell 6, first figure. Retained layout: three vertically stacked axes; (a) linear, (b) nonlinear, (c) independent. Status: plotting code updated; not executed.

**Restored presentation:** Restored the original three-row layout, full 30?0 kyr BP range, thin X/Y traces, process headings and per-panel legends. Lowercase panel labels remain; this source-only figure has not been regenerated.

**Caption draft:** Simulated records used to illustrate the observation workflow. Panels (a–c) show the linear, nonlinear and independent systems, respectively, with X in blue and Y in orange. Each retained record spans 30 kyr at 20-year spacing after warm-up; X and uncoupled Y have intrinsic memory times of 500 and 200 years. The coupled systems contain X→Y only, with a 20-year equation delay and long-run gain K=4, whereas the independent system has no cross-variable coupling. The scenarios reuse process innovations for a paired comparison and contain no measurement or chronology error. Ages decrease from left to right, so physical time advances toward younger ages.

### 01.2 Coupling transforms

**Files and source:** `figures/01_inspect_simulated_records/02_coupling_transforms.png` and `.pdf`; cell 6, second figure. One axis; no standalone panel letter required.

**Restored presentation:** Restored the original blue linear and orange nonlinear curves, axis labels and two-item legend. The single-axis figure has no redundant panel letter.

**Caption draft:** Coupling transforms applied to the simulated driver. The blue curve shows the linear transform and the orange curve the saturating hyperbolic-tangent transform. Both transforms are normalized using their theoretical stationary standard deviations under the specified driver distribution, rather than the variance of an individual realized record. The plot isolates the functional response g(X); its slope and amplitude are not native causal-method scores, and the per-step coupling coefficient also includes the long-run gain and Y memory coefficient.

### 01.3 Irregular samples

**Files and source:** `figures/01_inspect_simulated_records/03_irregular_samples.png` and `.pdf`; cell 9, first figure. Two by two: (a) X full record, (b) X detail, (c) Y full record, (d) Y detail. Status: plotting code updated; not executed.

**Restored presentation:** Restored the original gray process traces, filled sample markers, explicit detail-age heading, two first-row legends and detail-window shading. All four panels retain lowercase letters.

**Caption draft:** Independent irregular sampling of the linear coupled system. Panels (a, b) show X and panels (c, d) show Y, with full records in the left column and the 15.6–14.4 kyr BP interval in the right column. Gray lines show the true uniform records, and colored symbols mark the retained observations at their exact ages; shading in the full-record panels locates the detailed interval. Sampling targets are 27.5 years apart before rounding, with 10% jitter, producing gaps of 20 or 40 years on the native grid. X and Y use independent sampling locations, and their last retained ages need not equal zero. No age or measurement error has been added at this stage.

### 01.4 Sampling intervals

**Files and source:** `figures/01_inspect_simulated_records/04_sampling_intervals.png` and `.pdf`; cell 9, second figure. One axis; no standalone letter. Status: plotting code updated; not executed.

**Restored presentation:** Restored the original interval axis in kyr, X/Y sample points, black solid target line and its full ?Target before grid rounding? label. No extra separation or averaging of coincident gaps is applied.

**Caption draft:** Realized sampling intervals as a function of age. Symbols give each consecutive X or Y gap at its midpoint age, and the black solid line shows the 27.5-year target before rounding to the 20-year native grid. Independent jittered draws produce the discrete 20- and 40-year retained gaps, so the realized mean spacing need not equal the pre-rounding target. The interval axis is expressed in kyr (0.020 and 0.040 kyr correspond to 20 and 40 years), while the horizontal axis is age in kyr BP. Coincident symbols represent overlapping sampling intervals rather than additional averaging.

### 01.5 Age-uncertainty envelopes

**Files and source:** `figures/01_inspect_simulated_records/05_age_uncertainty_envelopes.png` and `.pdf`; cell 11. Two vertical axes: (a) C/L/E/J envelopes, (b) J coverage variants. Status: plotting code updated; not executed.

**Restored presentation:** Restored the original color cycle, solid curves, kyr SD units, explicit J coverage title and full coverage legends. Both original panels retain labels; the envelope calculations are unchanged.

**Caption draft:** Exploratory age-error envelopes and the effect of widening localized uncertainty intervals. Panel (a) compares constant (C), linear (L), exponential (E) and localized (J) proposal standard deviations over 30–0 kyr BP. C has a constant SD of 50 years, L and E reach 50 years at the old end, and J has three smooth peaks with a maximum SD of 100 years and a combined defined-interval coverage of 30%. Panel (b) holds those peak heights and centers fixed while changing their total interval coverage to 10%, 30% and 50% by changing interval width. Coverage denotes the summed defined interval lengths divided by record duration; it is neither a probability nor a threshold-exceedance fraction. These illustrative envelopes are not RMS-matched and are distinct from the main study's four-level peak-count design.

### 01.6 Simulated observations

**Files and source:** `figures/01_inspect_simulated_records/06_simulated_observations.png` and `.pdf`; cell 14. Two by two: (a) X observation detail, (b) X age error, (c) Y observation detail, (d) Y age error. Status: plotting code updated; not executed.

**Restored presentation:** Restored the original true-age traces, observed-age symbols, full input-SD and age-error legends, detail headings and kyr error units. All four panels retain labels and the original independent legends.

**Caption draft:** Effects of assigning uncertain ages to fixed sampled values in the linear system. Panels (a, c) compare observations plotted at their assigned ages with the same values at their true ages in the 15.6–14.4 kyr BP detail window; panels (b, d) show the corresponding age-error realizations for X and Y over the full record. X uses the linear envelope with an old-end proposal SD of 50 years, and Y uses three localized peaks with a maximum SD of 100 years and 30% defined-interval coverage. The gray bands indicate the input ±1 SD envelopes, not confidence limits or hard bounds, and the zero lines mark exact ages. Errors are correlated in sample order within each record and independent between X and Y; values and sample identities remain fixed and measurement error is zero. These illustrative observation draws can violate age order and are distinct from the accepted, stratigraphically ordered chronology members used in the main analysis.

## 02_2 Historical CCM lag investigation

Four figures are defined, not five: cell 5 defines the `save_figure` helper and does not create a figure. The current source retains its historical library choices: L=1112 for the embedding comparison and L=1119 for the selected E=3, tau=-1 setting. No 02_2 outputs are currently present in the active data/figure tree. Retain these as parameter-investigation figures; do not relabel them L=300 or infer realized findings from unexecuted current cells. Numerical settings below describe the source definition, not a new execution.

### 02_2.1 Original records

**Files and source:** `figures/02_2_inspect_ccm_lags/01_original_records.png` and `.pdf`; cell 6. Three vertical axes: (a) linear, (b) nonlinear, (c) independent. Status: plotting code updated; not executed.

**Restored presentation:** Restored thin blue/orange X/Y traces, original coupling/no-coupling headings, proxy-value labels and one legend. The 29?1 kyr BP limits and all three panel labels remain.

**Caption draft:** Uniform records for the CCM lag-parameter investigation. Panels (a–c) show the linear, nonlinear and independent systems, with X in blue and Y in orange. Each series contains 1,401 points over 29–1 kyr BP at 20-year spacing, without sampling, chronology or measurement perturbations. The coupled systems have X→Y only with a 20-year equation delay and K=4; the independent system has no cross-variable coupling. Shared process innovations permit matched comparisons across forms, while the trajectories represent one process realization per form. Physical time advances from older to younger ages.

### 02_2.2 Embedding comparison

**Files and source:** `figures/02_2_inspect_ccm_lags/02_embedding_comparison.png` and `.pdf`; cell 9. Two by three: (a–c) linear and (d–f) nonlinear; columns are (E, tau)=(3,−5), (3,−1), (2,−1). Status: plotting code updated; not executed.

**Restored presentation:** Restored explicit E, tau and embedding-span headings, the explanatory supertitle, target shifts in kyr, original blue/orange directions and shared rho scale. All six labels, pointwise bands and zero references remain.

**Caption draft:** Sensitivity of CCM lag profiles to the embedding history. Panels (a–c) show the linear system and (d–f) the nonlinear system; columns use (E, tau)=(3,−5), (3,−1) and (2,−1), corresponding to 200-, 40- and 20-year embedding spans. Lines show the mean cross-map correlation and shading the pointwise 5th–95th percentiles across 20 paired random libraries of 1,112 states, with a 200-year temporal exclusion and target shifts from −500 to +500 years. For X→Y, X(t+Tp) is reconstructed from Y's history, so negative target shifts refer to past X; the reverse direction is defined analogously. The vertical reference marks Tp=0, and shared seeds and anchors make the lag comparisons paired. The shaded ranges describe library-subsampling variability, not confidence intervals or calibrated causal probabilities.

### 02_2.3 Library count and search window

**Files and source:** `figures/02_2_inspect_ccm_lags/03_count_window_sensitivity.png` and `.pdf`; cell 14. Three by two: linear (a, b), nonlinear (c, d), independent (e, f); X→Y left, Y→X right. Status: plotting code updated; not executed.

**Restored presentation:** Restored original count colors, full count/sign legend, process/direction headings and peak-sign-stability supertitle. The shared fraction scale and 0.9 reference are unchanged; all six panel labels remain.

**Caption draft:** Stability of CCM peak signs across random-library counts and lag-search bounds. Rows show the linear, nonlinear and independent systems, and columns show X→Y and Y→X, with panels (a–f) ordered by row. Solid circles and dashed squares give fractions with negative and positive interior peaks, respectively, for 20, 50 and 100 libraries of 1,119 states. Counts use nested prefixes of the same 100 draws, and narrower search windows use subsets of a common ±500-year scan without changing its anchors. The horizontal line marks the 0.9 peak-sign rule; boundary, zero, tied or otherwise unresolved peaks account for the remaining fraction. This diagnostic evaluates an engineering screening rule, and it cannot turn a nonsignificant original Tp=0 CCM result into a retained arrow.

### 02_2.4 Selected lag distributions

**Files and source:** `figures/02_2_inspect_ccm_lags/04_selected_lag_distributions.png` and `.pdf`; cell 17. Three by two: linear (a, b), nonlinear (c, d), independent (e, f); mean profiles left, per-library optima right. Status: plotting code updated; not executed.

**Restored presentation:** Restored the original six-panel layout, per-panel scales, target shifts in kyr, direction colors, legends and explanatory supertitle. Original-seed profiles and optimum distributions retain panel letters; no holdout calculation was run.

**Caption draft:** CCM lag curves and library-wise optimum locations for the selected historical settings. Rows show the original linear, nonlinear and independent records; panels (a, c, e) show mean rho with pointwise 5th–95th percentile ranges, and panels (b, d, f) show the fraction of the 100 requested libraries whose unique optimum occurs at each plotted target shift. The embedding uses E=3 and tau=−1, with L=1,119 states, a 200-year temporal exclusion and shifts from −300 to +300 years in 20-year increments. Blue solid and vermilion dashed traces represent X→Y and Y→X; the zero line separates reconstruction of past and future source values. Boundary optima remain unresolved, and ties or failed curves receive no unique optimum, so the displayed frequencies need not sum to one. These panels show one record per process form; the two additional validation realizations defined in the notebook are not displayed here, and the optimum need not equal the model's 20-year equation delay.

## 02_5 GPDC backend comparison

All three existing PNGs were visually reviewed. Saved `data/processed/02_5_compare_gpdc_backends/run_info.json` confirms six fixed tests, four backends, T=1391, a shared 10,000-sample native null distribution and an RTX 4080 CUDA run. These are implementation diagnostics rather than full PCMCI+ searches or an evaluation of directional graph recovery. The original full backend identifiers remain in the legends.

### 02_5.1 Regular inputs

**Files and source:** `figures/02_5_compare_gpdc_backends/01_regular_inputs.png` and `.pdf`; notebook cell 4. Three vertical axes: (a) linear, (b) nonlinear, (c) independent.

**Restored presentation:** Restored the original three-row layout, thin blue/orange traces, process names on the ordinates and common age axis. Lowercase labels and the added PDF companion remain.

**Caption:** Uniform inputs used to compare GPDC regression backends. Panels (a–c) show the linear, nonlinear and independent records, respectively, with X in blue and Y in orange. The saved accurate-age records contain 1,401 observations from 29 to 1 kyr BP at 20-year spacing, without chronology or measurement perturbations. The backend comparison uses the same in-memory inputs for each implementation and retains 1,391 aligned test samples after accounting for the specified lag support. These are three process histories, not chronology-ensemble members or repeated GPU trials.

### 02_5.2 Backend comparison

**Files and source:** `figures/02_5_compare_gpdc_backends/02_backend_comparison.png` and `.pdf`; cell 16, first figure. One by three: (a) GPDC statistic, (b) native empirical p, (c) maximum residual difference.

**Restored presentation:** Restored the original connected points, full backend identifiers, case-label rotation, grid and logarithmic error panel, including the exact-zero-floor explanation on its ordinate. All three panel labels remain.

**Caption:** Comparison of native GPDC, stock GPDCtorch and aligned fixed-kernel GPyTorch regressions on six conditional-independence tests. Panels (a–c) show the GPDC statistic, its native empirical p value and the maximum absolute residual difference from native GPDC, respectively; cases combine three process forms and both source–target orderings. Each test compares source(t−1) and target(t), conditional on source(t−2) and target(t−1), using 1,391 samples and the same 10,000-draw native null distribution. The aligned CPU and GPU fits use float64 arithmetic and the current fixed GP model, whereas stock GPDCtorch uses its trained float32 configuration; residual differences are measured after the native input-standardization sequence. The dotted line in (c) marks an absolute diagnostic scale of 10−8, and any exact zero is displayed at 10−17 solely to permit a logarithmic axis. A reported p value of zero means no null exceedance in the finite reference sample, and agreement on these fixed tests does not establish equivalence of complete PCMCI+ searches.

### 02_5.3 Residual alignment

**Files and source:** `figures/02_5_compare_gpdc_backends/03_residual_alignment.png` and `.pdf`; cell 16, second figure. Two vertical panels: (a) Y regression residuals, (b) aligned-minus-native residuals.

**Restored presentation:** Restored the original long conditional-test title, full backend legends, thin traces and raw residual-difference scale. Both panels retain labels; the corrected ?standardized-input units? wording remains because residuals themselves are not standardized.

**Caption:** Regression residuals for the fixed nonlinear X(t−1), Y(t) conditional test. Panel (a) compares Y residuals from native GPDC, stock GPDCtorch and the aligned CPU and GPU implementations over the displayed old-age segment; panel (b) shows the aligned residuals minus native GPDC. All regressions condition on X(t−2) and Y(t−1) and use the native input-standardization sequence; the two panels therefore have different vertical scales, with differences in (b) expressed in units of 10−14. The displayed ages extend from the first aligned test sample at 28.8 kyr BP to 24 kyr BP, with older ages on the left. The difference curves diagnose numerical agreement for this test and are not uncertainty intervals or evidence about the physical magnitude of a causal effect.

## 02_6 CCM library-length sensitivity

Both existing PNGs were visually reviewed and checked against saved `settings.json`: L=350, 700 and 1119; 100 libraries per length; 31 target shifts; 1,371 common anchors; E=3, tau=−1 and 200-year exclusion. These are historical diagnostic lengths, not the current production L=300. Both figures now have PNG and PDF exports from the saved results.

### 02_6.1 Rho by library length

**Files and source:** `figures/02_6_check_ccm_library_length/01_rho_by_library_length.png` and `.pdf`; cell 10. Three by two: linear (a, b), nonlinear (c, d), independent (e, f); X→Y left, Y→X right.

**Restored presentation:** Restored the original solid blue/orange/green curves, panel-specific rho scales, original legend and explicit median/5?95% supertitle. All six panel labels and the added PDF companion remain.

**Caption:** Effect of CCM library length on the distribution of lagged cross-map skill. Rows show the linear, nonlinear and independent records, and columns show X→Y and Y→X, with panels (a–f) ordered by row. For each library length L=350, 700 or 1119, lines give the median rho and shading the pointwise 5th–95th percentiles across 100 random libraries sampled without replacement within each draw. Shared seeds and 1,371 common anchors pair the shifts from −300 to +300 years; E=3, tau=−1 and the 200-year temporal exclusion are fixed. The vertical line marks Tp=0, and each panel retains its own rho scale, as stated in the figure title. These ranges describe library variability on one fixed history per process form, not chronology uncertainty, confidence intervals or a full CCM convergence test.

### 02_6.2 Peak-lag distributions

**Files and source:** `figures/02_6_check_ccm_library_length/02_peak_lag_distributions.png` and `.pdf`; cell 12. Same six-panel mapping as 02_6.1.

**Restored presentation:** Restored the original circles and solid connecting lines, common fraction scale, gray boundary bins, original legend and explanatory supertitle. All six panel labels and the PDF companion remain.

**Caption:** Distributions of library-wise optimal CCM target shifts for three library lengths. Rows show the linear, nonlinear and independent systems, while columns show X→Y and Y→X, with panels (a–f) ordered by row. Each point gives the number of unique maxima at that 20-year shift divided by all 100 requested libraries at L=350, 700 or 1119; connecting segments only guide the eye between discrete bins. Gray strips identify the ±300-year search-boundary bins, and the vertical line marks zero shift. Boundary optima remain unresolved, while tied or failed curves have no assigned unique optimum and can make the plotted frequencies sum to less than one. The distribution measures library-subsampling variability on fixed records; its width need not decrease monotonically with library length, and its peak is not an estimate of the 20-year equation delay.

## 02_7 CCM neighbor-query reuse

Both existing figures were redrawn from saved settings, arrays and timing tables, and visually checked. The notebook's standalone CCM functions and all non-plot scientific cells are unchanged; only its panel-label import and plotting cell differ from the pre-style source. The pre-style source remains archived. The saved experiment is a sequential CPU lag-scan comparison at L=300, not a GPU test or a complete production-speed benchmark. Both figures now have PNG and PDF exports.

### 02_7.1 Numerical agreement

**Files and source:** `figures/02_7_check_ccm_neighbor_reuse/01_numerical_agreement.png` and `.pdf`; cell 14, first figure. Retained layout: three by two. Rows are linear, nonlinear and independent; panels (a, c, e) show median curves and (b, d, f) maximum absolute differences across libraries.

**Restored presentation:** Restored the original direction colors, native/reused curves, long panel headings and explicit L=300 supertitle. The retained correction keeps absolute-difference axes nonnegative and computed zeros visible, without adding epsilon to the data. All six labels and the PDF remain.

**Caption:** Agreement of native CCM lag scans with exact neighbor-query reuse. Rows show the linear, nonlinear and independent records; panels (a, c, e) show median cross-map skill and panels (b, d, f) show the maximum absolute native–reuse difference across libraries at each shift. Blue and vermilion distinguish X→Y and Y→X; solid and dashed/symbol curves distinguish native and reused-query medians. The displayed first-pass scans use 100 paired libraries of 300 states, 31 shifts from −300 to +300 years, E=3, tau=−1 and a 200-year temporal exclusion. All 37,200 paired rho comparisons across three records and two passes were exactly equal, with matching seeds, anchors and 1,200 peak rows; the right-column curves are computed zeros. Vertical references mark Tp=0. This comparison validates the tested sequential CPU lag scans with PyEDM 2.5.7, rather than full CCM surrogate computations, GPU behavior or production throughput.

### 02_7.2 Lag-scan timings

**Files and source:** `figures/02_7_check_ccm_neighbor_reuse/02_lag_scan_timings.png` and `.pdf`; cell 14, second figure. Retained layout: one by three, with (a) linear, (b) nonlinear and (c) independent.

**Restored presentation:** Restored median bars, both original timing observations, numerical speedups in panel titles and the explanatory two-pass supertitle. All three labels and the added PDF companion remain; timings were not rerun.

**Caption:** Wall time of native and reused-query CCM lag scans. Panels (a–c) show the linear, nonlinear and independent inputs; bars show median times and black points show each of two timing passes, whose execution order was reversed on the second pass. Each scan computes 100 libraries, 31 shifts and both directions at L=300 using one CPU numerical thread; timing includes preparation, query-cache handling and cleanup but excludes file saving and plotting. Median paired speedups were 3.81, 4.00 and 4.12 for the linear, nonlinear and independent inputs, respectively. Both passes are shown because two measurements do not define a reliable confidence interval. These timings concern the lag scan alone and are not a measured speedup of full CCM, PCMCI+, a 32-worker experiment or a GPU implementation.

## Reproducing the seven cached figures

Run `.\.venv\Scripts\python.exe agent\redraw_notebook_figures.py` from the project root. The maintained entry point loads saved CSV/NPZ results and explicitly executes only 02_5 plotting cells 4/16, 02_6 cells 10/12 and 02_7 cell 14. It does not run a full notebook or import causal-estimator implementations. The numerical data remain unchanged.

01 and 02_2 have restored plotting definitions with panel labels and captions but were not executed, because generating their missing inputs would run scientific work. Their scientific parameters retain their tutorial/historical meanings. A later manuscript version can choose final dimensions and figure selection; this revision preserves every original research panel.
