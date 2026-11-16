# Common-driver chronology experiment: fixed-memory first version

**System names (10 November 2026).** The figures use **Common driver only** and **Common driver + direct X to Y**. Tables abbreviate these as **Common only** and **Common + direct**. They are two generating systems, not labels for a statistical null and alternative hypothesis. Historical `h0_common_driver` / `h1_direct` file keys are retained so saved inputs, native results and resume identities remain unchanged.

Authorized on 9 November 2026. This replaces the proposed four-memory matrix with one fixed linear system. The root script is `03_check_common_driver.py`; all new outputs have that stem. Existing 02 experiments and their saved outputs are preserved.

## Question and comparison

Can an observed common driver still remove misleading X-Y evidence when its chronology, the response chronologies, or all three chronologies are uncertain? Does adjustment retain a genuine X-to-Y link when that link is added?

Common only contains Z-to-X and Z-to-Y. Common + direct adds only X-to-Y. Within each process seed, the two structures share Z, X, all standardized innovations, sampling positions and clocks. Common only is a common-driver negative control for the X-Y link, not a system of independent variables.

## Fixed linear dynamics

All times below are in years. With a 20-year step,

$$a_i=\exp(-20/T_i),\quad (T_Z,T_X,T_Y)=(500,200,200).$$

$$Z_t=a_ZZ_{t-1}+\sqrt{1-a_Z^2}\epsilon_{Z,t},$$

$$X_t=a_XX_{t-1}+(1-a_X)2Z_{t-1}+\sqrt{1-a_X^2}\epsilon_{X,t},$$

$$Y_t=a_YY_{t-1}+(1-a_Y)[2Z_{t-1}+dX_{t-1}]+\sqrt{1-a_Y^2}\epsilon_{Y,t}.$$

The innovations are independent standard normals across variables and time; d=0 in Common only and d=1 in Common + direct. Values use fixed units in which each variable's unforced background SD is one. They are not normalized to their realized total SD. Direct effects have a 20-year equation delay. The Z-X-Y mediated path has two steps without adding another delay parameter. The stable triangular transition matrix has diagonal eigenvalues a_Z,a_X,a_Y. Save its theoretical stationary covariance, actual variability/autocorrelation and a conditional direct-signal diagnostic. Parameters are fixed before inspecting chronology outcomes; unfavorable baselines are retained.

Use 30 kyr of retained data, the existing 29-1 kyr analysis window and 20-year grid (1401 points). All conditions use a common 5000-year warm-up (250 steps). Memory, nonlinearity, coupling strength, measurement noise and sampling density are not swept. Measurement noise is zero. Equal response memories and common-driver gains make Common only exchangeable in X/Y while independent innovations prevent identical series.

## Inputs and clocks

Use three process seeds: 20261009, 20261301, 20261302. There are two accurate-age references per structure/seed: original uniform data and mildly irregular observations interpolated to the common grid. A single synchronous sampling mask is fixed across records and histories, with the existing 27.5-year target spacing and 10% jitter, rounded to native steps.

Each process seed identifies an independent chronology block. X, Y and Z receive independent uncertain observed chronologies, generated with the existing observation model. Around each uncertain observed chronology, draw ten ordered members with C-type 40-year proposal SD and normalized-sample-order correlation length 1/60. Use the same existing whole-curve rejection rule; never sort ages or change values. Reuse each block's clocks across Common only/Common + direct and the following assignments:

| Condition | X/Y clocks | Z clock |
|---|---|---|
| irregular_accurate | accurate | accurate |
| z_only | accurate | uncertain |
| xy_only | independently uncertain | accurate |
| all_uncertain | independently uncertain | independently uncertain |

The z_only condition leaves X/Y inputs exactly unchanged. xy_only leaves Z unchanged; its X/Y inputs match all_uncertain member by member. All accepted inputs must cover the same grid without extrapolation. Record proposal acceptance, observed-center errors, effective member spread, realized individual errors and all three relative-error RMS values. Proposal SD is not realized error or a hard bound. Chronology members are conditional sensitivity draws, not independent climate histories or a calibrated posterior.

Total inputs: 2 structures x 3 histories x (2 references + 3 assignments x 10 members) = **192**. The three independent process/clock blocks support exploratory conditional comparisons, not precise population false-positive rates.

## Analysis

On every input, compare pairwise X/Y PCMCI+, LKIF and GC with their full X/Y/Z conditional counterparts. Retain pairwise CCM as a separate comparator; it does not condition on Z. Pairwise fits use the existing toolbox unchanged. New multivariate wrappers are in `toolbox/multivariate.py`, preserving old two-variable implementation/results.

- PCMCI+: native CPU GPDC, fixed optimizer=None, alpha/pc_alpha 0.05, lag 0..5 and 10000 native null samples. Cache only the compatible native null distribution. Per-direction Holm families and native graph arrows remain required; zero-lag unresolved evidence belongs to its variable pair.
- LKIF: corrected author package, native multivariate fit, maximum lag 5 and existing analytical p values, normalized score and per-direction Holm correction.
- GC: fixed order 1. Full and source-excluded target regressions share response rows and retain all remaining variables' lags. Native statsmodels OLS equation F tests match the pairwise SSR-F null and denominator degrees of freedom; do not switch to the VAR system-level F denominator. Score is log(restricted SSR/full SSR). Save complete fitted arrays and joint VAR diagnostics.
- CCM: unchanged paired X/Y pipeline, E=3, tau=-1, exclusion=10, L=300, ten original libraries, 100 swap surrogates, 100 paired libraries across 31 shifts and existing 0.9 sign rule. Retain complete raw/null/lag arrays and unresolved states.

The full uniform process is VAR(1). Interpolation/chronology deformation can invalidate that specification. Marginalizing Z need not leave a bivariate VAR(1), so pairwise baseline errors cannot all be attributed to chronology. Native scores are comparable within a method and analysis, not interchangeable causal-strength units. Adding Z can change normalized LKIF scores; adjusted-minus-unadjusted score changes are descriptive rather than a pure effect estimate.

Each case has seven native method archives and 26 directional rows: four pairwise fits x 2, plus three conditional fits x 6. Exact unchanged X/Y inputs permit a small paired reuse: z_only copies the same block/structure's irregular-accurate pairwise evidence, and all_uncertain copies its matching xy_only pairwise evidence. Check input arrays and estimator seeds before copying, retain every native archive, and preserve scientific failures. Execute the 72 source cases before the 120 dependent cases. This requires 864 actual native method calls, while retaining all 1344 case-method archives. Full-network diagnostics use the six conditional edges; primary outcomes concern the X/Y pair. Preserve significant arrows and unresolved evidence separately. Failures/pending cases stay in requested denominators. Report per-history q only from jointly usable X/Y decisions, with counts; pooled detection/state fractions are not q.

## Figures and summaries

1. **00_design_and_inputs:** Common only/Common + direct graphs, fixed equations/settings and representative process/only-Z perturbation diagnostics. Select the first seed and member 0 in advance, with a fixed 17-13 kyr window; do not choose an attractive event after fitting.
2. **01_effect_of_conditioning:** Common only false X/Y arrows and Common + direct forward retention/reverse errors, comparing pairwise and conditional analyses within each chronology assignment and method. Direct labels show counts; unresolved cases remain visible.
3. **02_xy_decision_states:** complete X/Y states (neither, X-to-Y only, Y-to-X only, both, unresolved, unavailable) reveal what detection fractions hide. Display requested denominators and distinguish accurate references from thirty member inputs.
4. **03_pairwise_ccm_evidence:** original versus lag-screened CCM outcomes and Tp=0 scores, explicitly labeled as an unadjusted pairwise comparator.

Keep lowercase panel labels, units, helpful counts and readable legends. Use consistent direction colors and a clear distinction between analysis mode and chronology assignment. Export PNG/PDF. Detailed English captions are stored in `agent/figure_review/06_common_driver.md`. Save per-seed directional support, paired score changes, X/Y states, full conditional graph states and complete native evidence alongside member-level results.

## Execution and validation

Use 32 outer processes with one numerical thread each, native GPDC caching and existing exact CCM neighbor reuse; no worker benchmark or GPU migration. Case workers own separate files. The 26-row edge table is the final commit after all seven native archives; explicit scientific failures are retained. Bounded Windows atomic-replace retries and configuration/input/native-version checks protect resumption.

Before production, verify theoretical dynamics and paired innovations, all 192 prepared identities, fixed-grid coverage, clock pairing, the unchanged X/Y and Z controls, direct native multivariate outputs, Windows-spawn equality, failure retention and resume guards. Use `agent/verify_common_driver.py`. After production, audit all 1344 native archives and 4992 rows, reconstruct summaries, review all figures/captions and check preservation of old outputs. A completed execution with unresolved or unsuccessful scientific decisions is not rerun selectively.

Memory and nonlinear variants remain future sensitivity branches. The first implementation and its figures must keep the single fixed memory explicit.

## Completed execution

The 9 November 2026 Windows run completed all 192 inputs in 2046.07 seconds with 32 workers. Minimum sampled available physical memory was 10.13 GiB. The paired-input/native/Windows-spawn preflight passed at `data/benchmarks/common_driver_checks/20261109T211251881262Z/report.json`.

The production audit checked all 4992 directional rows, 1344 native archives, 480 exact archive copies, 192 input files against preflight, complete CCM arrays, 3072 GC fit records and exact reconstruction of all four summary tables. Every method call completed successfully; unresolved scientific decisions remain in the results. All 24,305 pre-existing scientific files retain their bytes, sizes and modification times. The four figures and captions passed visual and table-consistency review. Reports, source snapshot and resource logs are in `data/logs/common_driver_20261109T211815766544Z/`; the interpretation is in [03_common_driver_results.md](03_common_driver_results.md).


## Figure revision (10 November 2026)

The original four views are retained with descriptive system names. Two additional saved-result views show the conditional-minus-pairwise gain in exact X-Y recovery (percentage points) and method-wise exact recovery curves, with the two systems and analysis modes kept separate. These use the existing correct_xy rule and requested denominators, not a new estimator or significance threshold. Figure 03 now explicitly labels reconstruction direction, per-input rho, median and descriptive IQR. Captions and the results document explain the metric and GC interface.
