# Figure style and caption catalogue

Updated on 10 November 2026. At the user's request, the broad visual restyling has been rolled back to the pre-style research figures. Original colors, fonts, line widths, titles, numerical annotations, legends, explanatory footers and panel layouts are restored. Lowercase panel letters, verified scientific-label corrections and the detailed caption files are retained. Later requested branch-specific additions are described below.

The catalogue now covers 45 exported figures and 10 additional figure definitions in the two unexecuted inspection notebooks, including four memory-sensitivity figures, six fixed-memory common-driver figures, two realization-by-Z-age-scheme figures and three Z-error-model figures. The common-driver figures retain the original four views and add two requested method comparisons on 10 November 2026. PNG/PDF versions of the same figure count once. Historical archives, temporary checks and third-party environment assets are outside the current figure set. Dense panels remain intact. Journal-specific simplification is deferred until manuscript preparation; the current figures should remain informative when viewed on their own.

## Complete catalogue

Each linked document retains the detailed review and English captions. The rollback notes identify the current representation; earlier simplification proposals are superseded. Captions complement the restored in-figure explanations. Captions for the ten unexecuted definitions remain code-based drafts.

| Source | Existing exported figures | Additional definitions without current exports | Detailed plans and captions |
|---|---:|---:|---|
| 02 causal pipeline | 11 | 0 | [Main pipeline](02_main_and_strength.md) |
| 02_1 coupling strength | 1 | 0 | [Coupling strength](02_main_and_strength.md) |
| 02_3 uncertainty escalation | 4 | 0 | [Uncertainty and sampling](03_uncertainty_and_sampling.md) |
| 02_4 sampling sensitivity | 7 | 0 | [Uncertainty and sampling](03_uncertainty_and_sampling.md) |
| 01 simulated-record inspection | 0 | 6 | [Notebook figures](04_notebook_figures.md) |
| 02_2 CCM parameter inspection | 0 | 4 | [Notebook figures](04_notebook_figures.md) |
| 02_5 GPDC comparison | 3 | 0 | [Notebook figures](04_notebook_figures.md) |
| 02_6 CCM library length | 2 | 0 | [Notebook figures](04_notebook_figures.md) |
| 02_7 CCM neighbor reuse | 2 | 0 | [Notebook figures](04_notebook_figures.md) |
| 02_8 memory sensitivity | 4 | 0 | [Memory sensitivity](05_memory_sensitivity.md) |
| 03 common driver | 6 | 0 | [Common driver](06_common_driver.md) |
| 03_1 realization by Z age scheme | 2 | 0 | [Crossed design and X-Y recovery](07_process_chronology_cross.md) |
| 03_2 Z error models | 3 | 0 | [Error models, recovery and scheme detail](08_driver_error_models.md) |
| Total | 45 | 10 | 55 figure definitions |

The original review covered all 30 existing PNGs alongside plotting code and saved settings. The 10 missing exports received a code-based review; their original plot content is restored with panel labels, without running their calculations. The original 55 figure files are retained in `data/archives/pre_figure_style_20261108T121553339241Z/figures_and_sources.zip`, together with the pre-edit source. The revision adds five previously missing PDF companions, giving 30 PNG/PDF pairs.

## Redraw from saved results

From the project root:

```powershell
.\.venv\Scripts\python.exe agent\redraw_figures.py
```

This invokes the three maintained `agent/redraw_*_figures.py` scripts to read existing CSV/NPZ/JSON results and call plotting code only. It does not launch the research scripts' calculation entry points or regenerate chronology members. The notebook runner executes only explicitly identified plotting cells with saved arrays. Missing required artifacts produce an error, not a new experiment. The 01 and 02_2 notebooks are deliberately excluded because their current saved artifacts are absent.

That launcher redraws the original 30 figures. For the four new 02_8 figures, set `RUN_CAUSAL=False` in `02_8_check_memory_sensitivity.py` and run that script with its matching saved scientific settings. It regenerates summaries and figures from committed results without causal fits. To preserve summary files as well, import the script and call `plot_results` with the saved summary tables and diagnostics directly, as used for the presentation revisions. All four actual PNGs were reviewed with their captions on 9 November; the four PDF companions passed structural checks. The original production review is `data/logs/memory_sensitivity_20261109T114908564938Z/figure_review.json`.

The two 03_1 figures show the realization-by-Z-age-scheme design and exact X-Y recovery matrices. Their 54 crossed cells and 18 accurate-age reference labels were checked directly against saved directional evidence. The actual PNGs and rendered PDF companions passed visual review on 10 November; see `data/logs/process_chronology_cross_20261110T124145323892Z/figure_review.json`. For this completed branch, `RUN_CAUSAL=False` reads the committed results without causal fitting; calling `plot_results` directly also leaves the saved scientific tables untouched. The original 03 figures and results remain unchanged.

The current reading order is graph recovery (01), CCM decision states (02), lag profiles and peak-sign diagnostics (03), and paired score changes (04). The displays retain all three histories without seed-specific shapes. The readability revision adds accurate-to-perturbed connectors, requested-input state bars with unresolved evidence, saved per-input library-sign diagnostics, and an explicitly expanded independent-system score scale. Numerical curves, scientific settings and all 3616 branch data files remain unchanged. The latest review is `data/benchmarks/memory_figure_checks/20261109T180042688693Z/report.json`; the preceding figures and source are retained in `data/archives/memory_figures_before_readability_20261109T180042688693Z/`. The earlier order-change archive and report with timestamp `20261109T174054354018Z` remain available.

## Journal requirements and project choices

The following sources were checked on 8 November 2026. Current AGU instructions take precedence over older print-production guidance. The concrete style settings below are project recommendations, not a claim that GRL mandates one palette or an exact Matplotlib font size.

| Source | Relevant guidance | Consequence for this project |
|---|---|---|
| [AGU text and graphics requirements](https://www.agu.org/publications/authors/journals/text-graphics-requirements) | Use lowercase letters for figure parts; combine all parts of a figure in one file; keep figure titles/captions out of the graphic; put readily transferable explanatory material in captions. Accessible color choices and redundant patterns/symbols are encouraged. PDF is an accepted format. | Add panel letters and simplify in-figure prose. Keep necessary axes, units, legends and short facet headings. Use line styles or symbols alongside color. |
| [AGU submission checklist](https://www.agu.org/publications/authors/journals/submission-checklists) | GRL has a 12-publication-unit limit; revised submissions require complete individual figure files. | Treat the current collection as a research atlas. Select and combine a small main-text set later; retain diagnostics in supporting information or the repository. The limit does not mean 12 figures are available in addition to the text. |
| [Wiley figure preparation](https://authors.wiley.com/author-resources/Journal-Authors/Prepare/manuscript-preparation-guidelines.html/figure-preparation.html) | Prefer PDF for line art. Recommend readable final-size text, figures approximately 80–180 mm wide, 600 dpi for rasterized line art and 300 dpi for images. | Use vector PDF for plots, with PNG previews. Choose 85 or 170 mm widths and inspect the actual output at those sizes. Raster resolution does not apply to vector paths. |
| [Older AGU graphics guide hosted by Wiley](https://onlinelibrary.wiley.com/page/journal/21699402/homepage/graphics.htm) | Lists Arial/Helvetica/Times/Symbol, font embedding, a 0.5-point minimum line, 85/170 mm column widths, and print-era GRL 6-point lettering. It also discusses CMYK and older submission systems. | Use it as supplementary technical context. Do not present the old 6-point GRL setting or CMYK workflow as a newly verified current requirement. Select larger readable lettering and embedded fonts; retain RGB research masters unless production requests a conversion. |

The current AGU page does not specify an obligatory font family, point size or color palette. Our choices below emphasize readability and consistent scientific meaning. Full submission requirements should be checked again when the manuscript is ready.

## Lessons from the previous GRL project

Read-only reference: `D:/VScode/Last_deglaciation`. The reviewed code includes `Figure_1_Ens.ipynb` (cell 9), `Figure2_Ens_v6.ipynb` (cell 33), and `Figure3_v12.ipynb` (cell 15); cell numbers here are zero-based. Their embedded figure outputs were also inspected. The folder contains several exploratory versions, so these are concrete style examples rather than an assertion about which notebook produced every final published panel.

Useful features to retain are bold lowercase panel letters, aligned time axes, a common age direction, stable colors within a figure, uncertainty shown as shaded intervals, and PDF export. In particular, the stacked Figure 2 and Figure 3 examples organize related records around a shared age axis rather than repeating a full explanatory title for every trace.

Do not copy every numerical setting. The reference notebooks use canvases such as 14 by 10 inches or 5 by 16 inches, 12–14-point panel letters and some 1.8–2.5-point strokes. Those values only make sense together with the final scaling. They also contain saturated colors, multiple y axes and substantial in-panel annotations; these are not necessary for the current simulation figures. Use the reference's panel discipline and scientific hierarchy, then set sizes for the current layouts explicitly.

## Current figure policy after rollback

- Preserve each figure's original research canvas, panel grid, colors, fonts, line widths and informative annotations. Do not simplify a figure merely because its caption contains the same information.
- Retain bold lowercase panel letters in reading order. Colorbars and secondary axes are not separate panels; the single sampling matrix needs no redundant letter.
- Keep verified corrections to scientific labels, units and saved-result provenance. In particular, current 02/02_3 results use CCM L=300, whereas retained 02_1/02_4 results use L=1119.
- Retain numerical cell values, counts, p-value labels, parameter summaries, explanatory titles and footers from the original research figures. Distinguish missing, unresolved and true-zero outcomes.
- Keep the detailed caption bank and PNG/PDF exports. No dense figure is split. A future manuscript pass can make a separately reviewed selection and layout.

`toolbox/plotting.py` now supplies only panel labels and PNG/PDF export. It does not override plotting styles or impose a common palette. Root scripts and notebook plotting cells retain their original composition. The standalone CCM estimator in 02_7 is unchanged.

The journal guidance and previous-project observations above remain useful references for a later publication pass. They do not authorize removing information from the current research figures.

## Provenance and scientific interpretation

| Figure family | What the present figures represent | Consequence for redrawing |
|---|---|---|
| 02 and 02_3 | Completed Windows L=300, N=10 pilot. The four-level result has 1920/1920 pairs, including 480 reused level-2 pairs. Completion audit and preservation checks passed. | Redraw from current saved results; no causal fits are needed for a visual revision. |
| 02_1 and 02_4 | Retained results from the earlier library-fraction 0.8 setting; current root code has since changed to L=300. | Read saved settings and use those values in captions. Do not relabel these as L=300, bypass resume guards, or launch a new experiment just to restyle them. |
| 01 and 02_2 | Historical illustrative/parameter-investigation definitions without current exported figures. Some settings intentionally differ from the main chronology study. | Treat captions as code-grounded drafts, explicitly identify exploratory scope, and decide separately whether these illustrations still belong in the final set. |
| 02_5 | Fixed conditional-test GPDC comparison, not a complete PCMCI+ GPU validation. | Timing captions must state the tested operation and avoid implying whole-pipeline acceleration. |
| 02_6 | Library-length sensitivity pilot at its recorded library sizes. | Preserve each tested library length and the paired-library interpretation. |
| 02_7 | Standalone exact neighbor-search reuse validation. | Only panel-label imports and plotting code change; its standalone scientific functions stay unchanged. Redraw from saved arrays and timings, retaining all six agreement panels. |
| 02_8 | Completed 432-pair all-method memory sensitivity, with three process seeds and ten shared level-4 C-C clocks per condition. | Preserve separate accurate-age baselines, per-seed results, requested denominators and full CCM lag-screen evidence. The four figures and captions retain the research appearance. |
| 03 | Completed 192-input linear common-driver experiment at one fixed memory setting, with three independent process/clock blocks. | Keep pairwise versus conditional evidence distinct, separate XY subgraph outcomes from full-graph recovery, and label CCM as pairwise. Captions retain counts, baseline exceptions and within-block chronology interpretation. |

In the main 02 study, chronology members describe sensitivity within one process realization per form. The sampling branch's ten mask groups are sampling repetitions, not ten independent climate histories. Its coupled systems have X→Y only with a 20-year imposed delay; the independent process has no cross-coupling. The 02_8 and 03 branches each retain three process histories, with their clock pairing specified separately. In 03, both structures have Z→X and Z→Y, and H1 additionally has X→Y. Age increases toward the past, so plots put older ages on the left and use `Age (kyr BP)` without reversing causal direction.

The four actual 03 figures passed visual and table-consistency review: eight PNG/PDF exports, 22 panel labels, no text outside the canvas, and all 1344 displayed XY case states consistent with saved results. The read-only review is recorded in `data/logs/common_driver_20261109T211815766544Z/figure_review.json`.

The 10 November 03 revision uses descriptive system names, clarifies the CCM reconstruction-score panels, and adds conditional-minus-pairwise exact XY recovery bars (04) and method recovery curves (05). The original directional/state plots remain. Exact recovery uses the existing correct_xy rule with all requested cases in the denominator; unresolved states stay visible separately. The pre-revision source, captions, report and figures are preserved in `data/archives/common_driver_figures_before_revision_20261110T105738358740Z/`. Review and data-preservation checks are recorded in `data/benchmarks/common_driver_figure_checks/20261110T105738358740Z/`.

Support q uses members with usable decisions in both directions. Correct-graph fractions use the requested member count and exclude unresolved evidence. Neither q≥0.9 nor the lag-sign stability threshold is a calibrated causal probability. Bands must be identified as actual quantiles across libraries/chronologies/records, not generically called confidence intervals. Native scores are comparable within a method, not as a common physical causal-strength scale.

## Review and later publication work

The three 03_2 figures passed numerical and visual review of all six PNG/PDF exports, 15 panel labels, 72 type cells, 216 scheme cells, 18 accurate reference entries and twelve chronology packages. Actual chronology diagnostics and requested denominators agree with saved arrays and directional results. The review is recorded in `data/logs/driver_error_models_resume_20261110T135403003994Z/figure_review.json`; the final English captions are in [08_driver_error_models.md](08_driver_error_models.md). No causal fits or production re-export were needed for this review.

1. Check all revised PNGs at their research viewing scale, panel order, margins, legends, units, colorbar limits and caption agreement. Verify PDF font embedding and preservation of scientific source files.
2. Keep source metadata authoritative when redrawing the older strength/sampling results. Existing calculation guards remain intact; no estimator calls are allowed as a side effect of redrawing.
3. Leave the ten unexecuted 01/02_2 figures identified as restored plotting code with panel labels but without new output. A drawing pass must not silently initiate their numerical experiments.
4. When manuscript assembly begins, select its figure subset and explicitly redesign at 85/170 mm, including grayscale/color-vision checks at final size. Revisit splitting only if requested at that stage. Keep the full research atlas available.

## Acceptance criteria for this drawing pass

- Every multi-panel output has unique, ordered lowercase labels and a matching caption; single-panel plots do not receive an unnecessary `a`.
- Text is readable on the retained research canvases; no large canvas is shrunk to journal width or split without authorization.
- Preserve original scales, category order and visual encodings; genuinely different method units remain separate.
- Main, supporting and implementation-check figures are clearly distinguished without deleting any scientific evidence.
- No figure silently changes library size, chronology count, process realization, sampling construction, lag convention, normalization or denominator.
- Keep the original title/footer prose and numerical annotations alongside the detailed captions.
- New exports and captions are visually checked together. Draft panel mappings are updated if the implementation changes the layout.

## TBK excitation diagnostic

The five subsequent TBK excitation figures have [English captions](10_tbk_excitation.md). All five PNGs and five independently rendered PDFs were reviewed against the completed 04_1 results on 13 November 2026. The review and final plotting-source snapshot are in `data/logs/tbk_excitation_20261113T121730984237Z/`; earlier figures and their numerical results were preserved. Retained directions, contemporaneous/lag-screen uncertainty, original CCM skill and GC residual diagnostics remain distinct.

## Earlier styling validation (superseded appearance)

The report in `data/benchmarks/figure_style_checks/20261108T123051789632Z/report.json` records the completed redraw check. All 30 original panel layouts are preserved; the single-matrix sampling figure remains a single panel without a redundant letter. All 30 rendered figures have ordered panel labels where appropriate and no visible text outside the canvas. Every figure now has a PNG and a PDF with embedded TrueType fonts and no Type 3 fonts. All 20,606 pre-existing files under `data/processed/` retain their content hashes and modification times, and the four scientific toolbox modules are byte-identical to the pre-style archive. All notebook code cells compile, with root notebook outputs kept empty. The 02_7 reference guard accepts the explicit presentation changes but rejects an injected scientific library-size change.

## Rollback validation

The pre-style archive supplies the restored plot content. All 30 existing figures are redrawn from saved data; no estimators or chronology generators are run. The rollback check is recorded in `data/benchmarks/figure_rollback_checks/20261108/report.json`. It checks panel labels, visible text bounds, unchanged scientific code/settings and preservation of all pre-existing processed data.
