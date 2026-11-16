# Fixed-memory common-driver experiment: results

**System names (10 November 2026).** The figures use **Common driver only** and **Common driver + direct X to Y**. Tables abbreviate these as **Common only** and **Common + direct**. They are two generating systems, not labels for a statistical null and alternative hypothesis. Historical `h0_common_driver` / `h1_direct` file keys are retained so saved inputs, native results and resume identities remain unchanged.

Completed on 9 November 2026. The calculation is in [03_check_common_driver.py](../03_check_common_driver.py), with the comparison specified in the [experiment plan](03_common_driver_plan.md).

## Main finding

Including an accurately dated common driver removes the misleading X-Y links in Common only for all three conditional methods and all three histories. Dating uncertainty can undermine that adjustment even when X and Y themselves remain unchanged. The clearest Z-only loss occurs for PCMCI+, which retains the correct X-Y state in 20/30 inputs, with eight unresolved cases and two cases containing incorrect directional arrows. This loss is concentrated in one process/clock block; LKIF remains correct in 30/30 and GC in 28/30. The result demonstrates a possible failure of adjustment under chronology uncertainty, rather than a uniform failure across methods or histories.

When a genuine X-to-Y link is present, the methods lose different information. Conditional PCMCI+ retains the forward arrow in every perturbed input but often adds reverse or unresolved evidence. Conditional LKIF mainly loses the forward arrow when X/Y chronologies are uncertain. Conditional GC retains the forward arrow throughout, with relatively few added reverse arrows. Complete X-Y recovery, directional detection and recovery of the entire three-node network therefore need separate interpretation.

## Fixed comparison and actual chronology errors

Common only contains Z-to-X and Z-to-Y; Common + direct adds X-to-Y. Intrinsic memories are fixed at 500 years for Z and 200 years for X and Y. Long-run gains are 2 for each common-driver branch and 0/1 for the X-to-Y branch, with a 20-year direct equation delay. Both structures share innovations, the synchronous sampling mask and chronology arrays within each block. There is no memory or nonlinear sweep in this experiment.

Three independent process/clock blocks use seeds 20261009, 20261301 and 20261302. Each structure/block has two accurate-age references and ten members for each of three uncertainty assignments, giving **192 inputs**. The thirty inputs per perturbed condition are ten conditional chronology members from each of three histories, not thirty independent histories. The original-uniform reference and the irregular-sampling accurate-age reference are kept separately.

All nine record/block clocks accepted ten of ten proposals, giving 90 accepted curves from 90 proposals; the observed centers had no order inversions in this run. Members were drawn around those uncertain observed centers. The C-type proposal SD is 40 years, but realized errors relative to truth are larger and are not bounded by 40 years. Over the analysis window, median individual-record RMS errors are 55.3 years for X, 54.7 for Y and 52.5 for Z.

| Assignment | Relative X-Y RMS, years | Relative X-Z RMS, years | Relative Y-Z RMS, years |
|---|---:|---:|---:|
| Z only | 0 | 52.5 [44.3, 71.3] | 52.5 [44.3, 71.3] |
| X/Y only | 76.5 [57.7, 93.3] | 55.3 [40.2, 68.3] | 54.7 [42.0, 70.0] |
| All three | 76.5 [57.7, 93.3] | 74.7 [55.2, 98.9] | 79.2 [60.4, 100.0] |

Entries are medians [minimum, maximum] across thirty member/block inputs. These are realized relative-clock diagnostics, not confidence intervals. The exact equality of X/Y inputs between the accurate and Z-only conditions, and between X/Y-only and all-three conditions, is preserved in the pairwise results.

## Accurate-age limitations

The full original process is VAR(1); interpolating sampled observations and omitting Z can change the effective statistical model. Accordingly, inaccurate pairwise results with accurate ages are baseline limitations, not chronology-induced errors.

- **Common only:** all conditional methods recover the correct X-Y state in all three histories, for both accurate-age references. Pairwise PCMCI+, LKIF and GC recover the expected no-X/Y-arrow state in only 1/3, 0/3 and 0/3 irregular accurate-age inputs, respectively.
- **Common + direct, LKIF:** conditional analysis misses X-to-Y in seed 20261009 for both accurate-age references, giving 2/3 correct X-Y states. The original-uniform adjusted directional p value is approximately 0.0598. This seed is retained throughout.
- **Common + direct, PCMCI+:** conditional recovery is 3/3 for original-uniform data but 2/3 for irregular accurate-age inputs. In seed 20261301, the forward arrow remains detected but unresolved evidence appears after sampling/interpolation. This limitation predates chronology perturbation.
- **Full-network GC:** Common only has 3/3 correct X-Y states but only 2/3 correct complete networks under either accurate reference. Seed 20261302 also detects X-to-Z, with original-uniform p approximately 0.0167.

Pairwise PCMCI+ also changes several individual states between the original and irregular accurate references, even where aggregate correct counts are unchanged. Saved case-level results retain those differences.

## Does conditioning improve the X-Y decision?

Each entry below is **pairwise correct count -> conditional correct count**, out of the requested count in the heading. A correct X-Y state has neither arrow in Common only and only X-to-Y in Common + direct, with no unresolved evidence. Failed or pending cases would remain in the denominator; none occurred here.

| Structure | Method | Accurate ages, n=3 | Z only, n=30 | X/Y only, n=30 | All three, n=30 |
|---|---|---:|---:|---:|---:|
| Common only | PCMCI+ | 1 -> 3 | 10 -> 20 | 6 -> 4 | 6 -> 5 |
| Common only | LKIF | 0 -> 3 | 0 -> 30 | 15 -> 29 | 15 -> 28 |
| Common only | GC | 0 -> 3 | 0 -> 28 | 0 -> 22 | 0 -> 26 |
| Common + direct | PCMCI+ | 0 -> 2 | 0 -> 22 | 3 -> 15 | 3 -> 8 |
| Common + direct | LKIF | 1 -> 2 | 10 -> 30 | 9 -> 5 | 9 -> 6 |
| Common + direct | GC | 0 -> 3 | 0 -> 30 | 14 -> 29 | 14 -> 28 |

The accurate-age column uses the irregular-sampling reference. Z-only pairwise counts simply repeat each unchanged baseline result ten times. Likewise, pairwise results are identical between X/Y-only and all-three assignments; they contain no information about Z's clock.

For Common only with only Z uncertain, conditional PCMCI+ correct counts are **10, 0, 10** across the three seed blocks. Eight of the ten incorrect cases in the middle block are unresolved, while the other two contain directional errors. The corresponding conditional GC counts are 8, 10, 10 and LKIF counts are 10, 10, 10. This heterogeneity limits generalization from the pooled proportions.

## What changes: forward, reverse or unresolved evidence?

Cells give **F/R/U**, the forward-retained, reverse-retained and unresolved-case counts. F, R and U can overlap and must not be added. Accurate-age cells have three requested inputs; perturbed cells have thirty. All failed and pending counts are zero.

| Structure | Analysis | Method | Accurate | Z only | X/Y only | All three |
|---|---|---|---:|---:|---:|---:|
| Common only | Pairwise | PCMCI+ | 1/2/0 | 10/20/0 | 9/15/10 | 9/15/10 |
| Common only | Conditional | PCMCI+ | 0/0/0 | 1/1/8 | 3/6/18 | 4/13/9 |
| Common only | Pairwise | LKIF | 2/2/0 | 20/20/0 | 9/9/0 | 9/9/0 |
| Common only | Conditional | LKIF | 0/0/0 | 0/0/0 | 1/0/0 | 1/1/0 |
| Common only | Pairwise | GC | 3/3/0 | 30/30/0 | 18/30/0 | 18/30/0 |
| Common only | Conditional | GC | 0/0/0 | 0/2/0 | 7/1/0 | 3/1/0 |
| Common + direct | Pairwise | PCMCI+ | 3/2/1 | 30/20/10 | 25/16/14 | 25/16/14 |
| Common + direct | Conditional | PCMCI+ | 3/0/1 | 30/1/7 | 30/10/9 | 30/11/15 |
| Common + direct | Pairwise | LKIF | 3/2/0 | 30/20/0 | 25/18/0 | 25/18/0 |
| Common + direct | Conditional | LKIF | 2/0/0 | 30/0/0 | 5/3/0 | 6/3/0 |
| Common + direct | Pairwise | GC | 3/3/0 | 30/30/0 | 30/16/0 | 30/16/0 |
| Common + direct | Conditional | GC | 3/0/0 | 30/0/0 | 30/1/0 | 30/2/0 |

For Common + direct, PCMCI+'s all-three correct count of 8/30 does **not** mean it retains the true forward direction in only eight cases: it retains that arrow in all thirty, while reverse and unresolved evidence prevent complete X-Y recovery. LKIF's loss under X/Y uncertainty includes forward non-detection in 25/30 or 24/30 cases. Its Z-only result rises to 30/30 forward retention despite its 2/3 accurate baseline; this does not establish a benefit of dating error, because changing the control record also changes how common-driver information is removed, and the full network is usually not recovered.

There is also one explicit example of high support for a wrong conditional direction. In Common only, seed 20261301 with all three clocks uncertain has PCMCI+ **q(Y-to-X)=0.9**, with 10/10 usable members and two cases containing unresolved evidence. The true X-Y graph contains neither arrow. This is a conditional example of support not guaranteeing correctness, not an estimated population frequency.

## The X-Y subgraph is not the complete network

These are correct **six-direction network** counts for the conditional methods, including the required Z-to-X and Z-to-Y links and absence of reverse links into Z.

| Structure | Method | Accurate, n=3 | Z only, n=30 | X/Y only, n=30 | All three, n=30 |
|---|---|---:|---:|---:|---:|
| Common only | PCMCI+ | 3 | 8 | 2 | 2 |
| Common only | LKIF | 3 | 4 | 6 | 1 |
| Common only | GC | 2 | 18 | 16 | 18 |
| Common + direct | PCMCI+ | 2 | 6 | 7 | 4 |
| Common + direct | LKIF | 2 | 4 | 0 | 0 |
| Common + direct | GC | 3 | 25 | 28 | 21 |

For example, Common only/Z-only LKIF correctly rejects both X-Y arrows in 30/30 inputs but recovers the full network in only 4/30. Its apparent success on the primary negative-control question therefore does not imply that all physical pathways were identified. The saved six-direction results retain the additional losses and arrows involving Z.

## Pairwise CCM is a separate comparator

The original Tp=0 CCM test is significant in both directions for **all 192 inputs**, including Common only. The lag screen changes those results as follows; cells again give F/R/U.

| Structure | Accurate, n=3 | Z only, n=30 | X/Y only, n=30 | All three, n=30 |
|---|---:|---:|---:|---:|
| Common only | 0/1/3 | 0/10/30 | 0/15/24 | 0/15/24 |
| Common + direct | 3/0/0 | 30/0/0 | 30/0/3 | 30/0/3 |

In Common only, every accurate-age input already has unresolved evidence and one retains the reverse arrow. Under X/Y uncertainty, fifteen reverse arrows are retained and twenty-four cases have unresolved evidence; no case recovers an unambiguous no-arrow X-Y state. These are limitations of the pairwise workflow in this common-driver system, present before chronology errors and modified by them.

In Common + direct, all thirty X/Y-perturbed inputs retain the forward arrow, but three have unresolved reverse evidence, leaving 27/30 complete X-Y states. Median original forward rho decreases from 0.793 in the irregular accurate reference to 0.765 with X/Y uncertainty; corresponding reverse medians are 0.823 and 0.791. These are Tp=0 scores, not scanned maximum correlations. Z-only results remain exactly unchanged because this CCM analysis does not use Z. Good Common + direct performance therefore does not demonstrate successful adjustment for a common driver.

## Scope and saved evidence

The conclusions are conditional on one linear structure family, one memory setting, one coupling choice, one sampling scheme, C-type 40-year proposals and the fixed estimator configurations. Three independent process/clock blocks provide a compact mechanistic comparison, while the ten members within each block measure chronology sensitivity. The design does not establish calibrated error rates, a universal method ranking or robustness to asymmetric/shorter response memory and nonlinear common driving. Per-block q and native scores remain available; scores have different meanings across methods, and adding Z can also change LKIF's normalization.

The run completed **192/192 cases**, **4992 successful directional rows** and **1344 native archives**, with 864 actual method calls and 480 exact-input archive reuses. All 4680 member-minus-accurate native-score contrasts are usable. Script runtime was **2046.1 seconds (34.10 minutes)** using 32 workers; minimum sampled available physical memory was **10.13 GiB** and peak sampled system commit was 26.69 GiB.

The [preflight report](../data/benchmarks/common_driver_checks/20261109T211251881262Z/report.json) passed input/innovation/clock pairing, native multivariate comparisons, Windows-spawn equality, failure retention and resume checks. The [run directory](../data/logs/common_driver_20261109T211815766544Z/) contains the source snapshot, logs and resource records; [current run status](../data/processed/03_check_common_driver/run_status.json) records completion. The final [artifact audit](../data/logs/common_driver_20261109T211815766544Z/artifact_audit.json) passed all 192 commits, 4992 directional rows, 1344 native archives, exact input regeneration, 480 byte-identical copied archives, all CCM lag/null/seed arrays, 3072 GC fit records and exact reconstruction of four summary tables. The [preservation audit](../data/logs/common_driver_20261109T211815766544Z/preservation_report.json) confirmed that 24,305 pre-existing files were unchanged. The [figure review](../data/logs/common_driver_20261109T211815766544Z/figure_review.json) passed all four PNG/PDF figure pairs and their 22 panel labels.

The first four figures follow the experiment's reading order:

1. [Design and input controls](../figures/03_check_common_driver/00_design_and_inputs.png).
2. [Effect of adding Z](../figures/03_check_common_driver/01_effect_of_conditioning.png).
3. [Complete X-Y decision states](../figures/03_check_common_driver/02_xy_decision_states.png).
4. [Pairwise CCM evidence](../figures/03_check_common_driver/03_pairwise_ccm_evidence.png).

PDF companions are stored alongside the PNGs. [English captions and reading notes](figure_review/06_common_driver.md) define the displayed counts and intervals; [saved result tables](../data/processed/03_check_common_driver/) retain per-block support, paired score changes, full networks and the complete native results.


## Revised comparisons and GC implementation (10 November 2026)

The original directional comparison and complete-state figure remain available. Two additional figures use the existing exact X-Y recovery definition, with no estimator calls or decision changes: [conditioning gain](../figures/03_check_common_driver/04_conditioning_gain.png) shows conditional minus pairwise recovery in percentage points, and [method comparison](../figures/03_check_common_driver/05_xy_recovery_comparison.png) shows one recovery curve per method, separating systems and analysis modes. The complete-state figure explains whether a low recovery fraction reflects missing, reversed, extra or unresolved evidence. The same strict success definition works when the true X-Y subgraph has no edges; F1 would require an arbitrary convention for its otherwise undefined perfect no-edge case.

For example, Common only/Z-only PCMCI+ improves from 10/30 pairwise to 20/30 conditional correct decisions, a +33.3 percentage-point adjustment gain. This remains worse than its 3/3 conditional accurate-age reference. A positive conditioning gain and deterioration from the accurate-age reference are therefore compatible.

The installed GC package is **statsmodels 0.14.5**, with fixed order one. Its [VAR interface supports multivariate conditional causality tests](https://www.statsmodels.org/stable/generated/statsmodels.tsa.vector_ar.var_model.VARResults.test_causality.html). Pairwise GC uses native `grangercausalitytests` and its `ssr_ftest`; the three-variable version fits native VAR and OLS models. For X to Y conditional on Z, it compares

$$Y_t=c+aY_{t-1}+bX_{t-1}+dZ_{t-1}+\epsilon_t$$

with the restriction b=0, retaining the histories of Y and Z on exactly the same response rows. The formal conditional p value comes from the native single-equation `OLS.f_test`, preserving the pairwise test's single-equation convention. With 1401 input rows, the pairwise test uses F(1,1397) and the conditional test uses F(1,1396). Native `VARResults.test_causality` is also called and saved as a diagnostic; its system F test instead uses denominator degrees of freedom 4188 here and is not the source of the formal p value. Both fits and tests are provided by statsmodels, rather than a custom regression estimator. See the [pairwise interface](https://www.statsmodels.org/stable/generated/statsmodels.tsa.stattools.grangercausalitytests.html) and [OLS F-test documentation](https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLSResults.f_test.html). The current online documentation may describe a newer release; the installed 0.14.5 source and saved fit fields were checked directly.

Figure 03c,d retain the continuous Tp=0 reconstruction scores, before the lag screen. Blue X-to-Y evidence means reconstructing X from Y's embedding; orange means reconstructing Y from X's embedding. Small dots are input-specific scores, large dots their medians, and whiskers span the middle 50%. Higher rho means more accurate reconstruction, not a larger direct causal effect. For example, Common + direct has an accurate-age reverse median rho of 0.823, yet none of its three reverse arrows survives the lag screen. The reconstruction score and the directional causal decision answer different questions.

The [figure-revision audit](../data/benchmarks/common_driver_figure_checks/20261110T105738358740Z/report.json) passed all six figures and twelve exports, all 1344 saved X-Y states, 24 conditioning-gain entries, 56 recovery entries and the CCM score summaries. All 2145 pre-existing branch data files and both causal toolbox modules retain their exact bytes, sizes and modification times; the scientific calculation digest is unchanged. All generation and fitting entry points were blocked during redraw and review. The preceding presentation is retained in the [revision archive](../data/archives/common_driver_figures_before_revision_20261110T105738358740Z/).
