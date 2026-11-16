# Chronology uncertainty and causal inference: engineering guide

This is a scientific research project, not a general-purpose application. Keep the scientific calculation visible and the code simple. The user prefers explicit, readable Python even when it is longer.

## Current authority and handoff

- **Windows is the authoritative computation platform** from the 2026-10-26 migration onward: Intel i9-14900KF, 32 GB RAM, NVIDIA RTX4080.
- Generate new processes, observations and age ensembles on Windows. The user does **not** require identical Mac/Windows chronology curves. Do not spend time transferring all Mac arrays to enforce that identity.
- Optimization must still preserve the scientific result for the **same input on the same platform**. Validate the complete numerical outputs, not just significant/non-significant arrows.
- The old Mac experiment 05 was stopped at the user's request. Its data are retained locally under `data/legacy_macos/`, with figures under `figures/legacy_macos/`. Do not delete them or silently merge them into Windows results.
- The retained Mac run has **1,009 complete input pairs: 480 reused from 02 and 529 newly calculated**. Partial native files are also preserved. It is not a complete four-level result.
- The default is **10 generated and analyzed chronology members**, unless the user requests another count. Each process form still uses one independent process realization in the main experiments.
- Do not restart long experiments on the Mac. Short implementation checks are acceptable. Full runs and worker benchmarks belong on Windows.

## Files and responsibilities

| File | Purpose |
|---|---|
| `toolbox/process.py` | Generate the uniform process and known model truth. |
| `toolbox/observation.py` | Irregular sampling, uncertainty envelopes, observed records and ordered chronology ensembles. |
| `toolbox/causal.py` | Thin calls to native causal packages, uniform resampling and CCM lag screening with bounded native neighbor-query reuse; no plotting. |
| `toolbox/multivariate.py` | Three-variable native PCMCI+/GPDC, LKIF and conditional GC wrappers for 03; existing pairwise wrappers are unchanged. |
| `toolbox/evaluation.py` | Ensemble support, complete member states, native-score summaries and independent-realization evaluation. |
| `toolbox/tbk_model.py` | Near-verbatim annual four-box model core, with original noise and extended prescribed events for longer runs. |
| `toolbox/tbk_observation.py` | Fixed TBK sampling and independent ordered chronology layers around uncertain observed centers. |
| `01_inspect_simulated_records.ipynb` | Inspect generation and simulated observations. |
| `02_check_causal_pipeline.py` | Two accurate-age baselines and 16 ordered type pairs at uncertainty level 2. |
| `02_1_check_coupling_strength.py` | Original-uniform-record native-score ranking under K=1/2/4, with three process seeds. Former 03. |
| `02_2_inspect_ccm_lags.ipynb` | CCM embedding/library/lag-window parameter investigation. Former 04. |
| `02_3_check_uncertainty_escalation.py` | All four levels and 16 ordered pairs; imports completed level 2 from 02. Former 05. |
| `02_4_check_sampling_sensitivity.py` | Accurate-age density, irregularity, asynchronous sampling and age-gradient controls; 492 input pairs. |
| `02_5_compare_gpdc_backends.ipynb` | Six fixed conditional tests on saved original uniform records; native GPDC, stock GPDCtorch and aligned float64 GPyTorch on CPU/CUDA. No PCMCI+ search or production replacement. |
| `02_6_check_ccm_library_length.ipynb` | Compact lag-variability diagnostic on the three saved original uniform records: libraries of 350/700/1119 states, 100 paired draws and 31 shifts; no surrogate tests or production changes. |
| `02_7_check_ccm_neighbor_reuse.ipynb` | Notebook-local reference for exact native KDTree query reuse at L=300; scientific functions retained, with style-only plotting updates on 2026-11-08. |
| `02_8_check_memory_sensitivity.py` | Independent 432-pair memory-by-chronology comparison, all four methods, three paired process seeds and reused level-4 C-C clocks. |
| `03_check_common_driver.py` | Fixed-memory linear common-driver comparison, with and without X-to-Y, 192 inputs and pairwise versus conditional analyses. |
| `03_1_check_process_chronology_cross.py` | Cross three saved simulation realizations with three saved Z age schemes; X-Y outcomes only, 120 new conditional inputs and reused references/diagonals. |
| `03_2_check_driver_error_models.py` | Z-only C/L/E/J sensitivity with the same three-by-three crossing; 732 cases including 192 complete source imports and 540 new conditional inputs. |
| `04_0_inspect_tbk_records.ipynb` | Inspect four-box model, original process noise and graded chronology perturbations without causal fitting. |
| `04_check_tbk_chronology.py` | 123-case TBK endpoint chronology pilot, all four pairwise methods and complete native evidence. |
| `04_1_check_tbk_excitation.py` | Accurate-age TBK excitation diagnostic: three paired histories, three noise amplitudes, a known-forcing-response subtraction, and a separate baseline ParCorr comparison. |
| `SYSTEM_DESIGN.md` | Plain-English equations and derivations for the linear/nonlinear generating systems, memory, innovation scales and long-run gain. |
| `agent/run_experiment.py` | Portable sequential launcher: one numerical thread per worker, console plus one saved log per script. |
| `agent/benchmark_workers.py` | Windows 16/24/32-worker throughput comparison on 96 fixed prepared02 inputs. |
| `agent/verify_optimizations.py` | Portable same-input comparison of uncached/fresh native calls and optimized calls. |
| `agent/verify_ccm_neighbor_reuse.py` | Local full-output acceptance against the 2026-11-07 source freeze, including L=300/1119, failures and Windows spawn. Requires that archive and saved test inputs. |
| `agent/explore_cuda.py` | Optional CUDA hardware/matrix-operation probe; **not** a GPU causal estimator. |

Branch numbering expresses a relationship to 02, not a requirement to run every notebook before the main experiment. Production calculations do not import code or data from `tmp/`.

`agent/` contains maintained, versioned support scripts and requirements. Temporary agent work still belongs in `tmp/`. Scientific scripts and notebooks remain at the project root; `toolbox/`, `data/` and `figures/` retain their existing locations.

`Note/` contains local Chinese reading notes requested by the user and is excluded from Git. Write Chinese documentation only in this directory. All new or revised project documentation, notebook prose, code comments and figure text outside `Note/` must be in English. Do not translate or rewrite historical artifacts as part of an unrelated task.

## Scientific question and scope

The main question is how equal or different chronology-uncertainty types assigned to two records affect retained causal arrows and native method scores. The ensemble provides within-record support. One highly consistent ensemble does not establish causal correctness across independently generated climate histories.

The main chronology study deliberately keeps sampling, proxy values, dynamics and estimator settings fixed. Experiment 02_4 is a separate accurate-age sampling sensitivity branch. There is no main-study systematic age shift, feedback network, five-variable network, tipping experiment or measurement-noise sweep. The separately authorized 03 branch tests three-variable common driving at fixed memory. The subsequent 04 branch uses the existing four-box TBK physical model; real-data validation remains later work.

The 2026-11-10 research discussion is recorded in `agent/research_direction_and_tbk_notes.md`. The accepted question is whether chronology uncertainty changes climate-mechanism interpretation. Additional three-variable mistiming sensitivity remains deferred. The user subsequently authorized the TBK plan in `agent/04_tbk_chronology_plan.md`, with the latest duration revision retaining 30 kyr after a 9 kyr startup and sampling every 20 years on the common 29-1 kyr grid (1401 points). Both southern T1 and northern T4 receive independent C-C uncertainty at every nonzero level (20/50/100/200 years). Keep the four-box core, annual integration and intrinsic noise; no measurement noise or feedback into the externally prescribed transport is added.

The original notebook is preserved at `data/archives/tbk_original_20261110T214805311811Z/One_D_model_TBK.ipynb`. New files are `toolbox/tbk_model.py` (cell-12 core, southern diffusion 0.005), `toolbox/tbk_observation.py` (shared sampling and chronology), `04_0_inspect_tbk_records.ipynb` (inspection only), and `04_check_tbk_chronology.py` (123 input pilot). The model defaults remain the original 10 kyr for migration checks. The 39-kyr run extends only the existing event table from 20 to 26 intervals; it preserves the 1500-year cycle and update order. Every same-input migration array must match exactly.

TBK uses three independent process/chronology blocks, one accurate input and ten members at each of four nonzero levels per block. CPU GPDC (maximum lag five), LKIF (five), GC (one), and CCM (E=3, tau=-1, exclusion=10, L=300, original 10 libraries, 100 surrogates) retain existing wrappers. CCM scans 61 shifts -30..30 steps (+/-600 years), 100 paired libraries, with raw and screened results retained. Age-error correlation is 1000 years = 1/30 of the sampled support; members surround uncertain observed centers. Whole-curve order and coverage rejection is explicit. A two-layer member has unconditioned relative inter-record error SD 2*sigma; report achieved errors separately. Use 32 workers and single numerical threads. Interpret endpoint evidence and baseline-to-uncertain changes, not complete four-node graph recovery. Preserve previous 02/03 sources/results.

### Time and data conventions

- Age is in **kyr BP**, with 1950 CE as zero; positive values are older.
- Data rows advance forward in physical time: `time` increases and `age` decreases.
- `calendar_year = 1950 - 1000 * age`. Do not reverse the causal direction when plotting BP ages.
- Plots put older ages on the left and use `Age (kyr BP)`.
- Native process step and analysis step are 20 years (`0.020 kyr`). The retained process spans 30 kyr, 1501 points.
- Analysis uses the fixed 29–1 kyr window, 1401 points. Every accepted chronology member must cover the same grid. No member-specific cropping or extrapolation.
- The accurate-age references are: (1) original uniform data, cropped only; (2) irregular observations with their accurate ages, interpolated onto the same grid. There is no separate downsampled baseline.
- Main irregular sampling uses 27.5-year pre-rounding target gaps at both ends with 10% jitter; native retained gaps are 20/40 years. This mild sampling was chosen to separate chronology effects. The separate 02_4 branch tests stronger sampling changes, including age-dependent density, without changing this main-study setting.

### Dynamics and truth

- Coupled systems have X→Y only; the independent case has no cross-variable coupling.
- The true coupled-model delay is 20 years, one native step. It is not the uncertainty SD.
- X/Y memory times are 500/200 years. Warm-up is at least five times the longest memory, plus needed lag history.
- `ar = exp(-dt / memory)`; innovation SD is scaled by `sqrt(1-ar**2)` to maintain the specified unforced background SD.
- `Y[t] = ar_y * Y[t-1] + coupling * g(X[t-lag_steps]) + innovation`.
- Main long-run coupling gain K = 4; per-step coefficient is `(1-ar_y)*K`. Linear and nonlinear g are normalized using the theoretical driver distribution. K is not any estimator's native score.
- Measurement noise is zero in the main experiment. Do not silently alter noise, gain, sampling or history to clean up false arrows.

### Chronology model and ensemble

| Level | C/L/E peak proposal SD | J peak count | J defined-interval coverage |
|---|---:|---:|---:|
| 1 | 10 years | 1 | 10% |
| 2 | 20 years | 3 | 30% |
| 3 | 40 years | 5 | 50% |
| 4 | 80 years | 7 | 70% |

- C is constant SD, L increases linearly with age, E increases exponentially with fixed curvature 2. They are idealized envelopes, not rigid assignments to archive types.
- J uses smooth compact Gaussian-like peaks: peak SD 10 years, each interval 2.8 kyr wide. Centers are added in fractions `[0.5,0.2,0.8,0.35,0.65,0.05,0.95]` of the29–1 kyr window. Earlier peaks remain fixed. Coverage is total interval width/28 kyr, not a threshold-exceedance statistic.
- The four levels are not RMS-matched across types. Increasing J count intentionally increases total uncertainty.
- At each level use all 16 **ordered** X/Y pairs. C–J and J–C are different. X and Y use the same level index; do not introduce cross-level combinations.
- Whole chronology proposals have temporally correlated errors in sample order; X/Y clocks are independent. Reject a whole curve if it violates stratigraphic order. Never repair by sorting/clipping ages or moving proxy values.
- Chronology members are drawn around an observed, uncertain chronology, not around secretly supplied true ages. Truth is available for diagnostics only.
- Reuse each record/type/level clock across process forms and pair assignments to hold incidental variation fixed.
- A proposal SD is not a hard bound or automatically the SD after order conditioning. Retain acceptance, coverage and effective-spread diagnostics.
- `np.linalg.eigh` is used in covariance factorization. Same seed across different numerical platforms need not give the same individual curve. Windows results are authoritative; validate optimizations on identical local inputs instead.

## Methods and interpretation

Use the existing packages. Do not replace an algorithm with a hand-written estimator, a different conditional test, or an approximate GPU method without a separate scientific decision.

- **PCMCI+:** Tigramite GPDC, maximum lag 5 grid steps (100 years), alpha/pc_alpha 0.05, 10000 native null samples. The GP kernel optimizer is disabled (`optimizer=None`); conditional regressions still fit each input. Directional Holm families include cross-lags 0..5, and a retained direction additionally requires the native graph arrow. Preserve unresolved/conflicting contemporaneous links.
- **LKIF:** author's corrected package, maximum lag 5. Keep the existing analytical significance calculation and native normalization. Do not silently replace it with a different bootstrap scheme.
- **Granger:** statsmodels, fixed VAR order 1, with the existing wrapper's test/score. Do not switch to BIC or vary order by chronology member.
- **CCM:** PyEDM 2.5.7; E=3,tau=-1, temporal exclusion 10 steps. New production configuration uses **fixed library length 300**, as selected by the user on 2026-11-07, consistently for original scores, nulls and lag screening. Previously completed results used fraction 0.8 (1119 states here) and remain preserved; they cannot be resumed or imported as L=300 results. Historical parameter notebooks retain their own settings. Original Tp=0 skill averages 10 libraries; the directional swap-null test uses 100 surrogates. X→Y means reconstructing X from Y's embedding.
- **CCM lag screen:** both directions, 31 target shifts from −15 to +15 steps (±300 years), 100 paired random libraries. The same library seed is used across shifts, and all shifts use common anchors. Negative interior peaks support temporal compatibility; positive peaks reject the arrow; boundary/zero/tied/unstable peaks remain inconclusive. The 90% sign criterion is an engineering rule, not a calibrated causal probability.
- The screen cannot upgrade an originally nonsignificant CCM arrow. Keep raw Tp=0 p/rho, raw votes, screened votes and failure/unresolved flags. Native CCM strength remains Tp=0 rho, not the maximum scanned rho.
- Score magnitudes from PCMCI, LKIF, CCM and GC are not the same physical causal strength. Compare within-method score changes, spreads or rankings.

### Ensemble evaluation

- `q_x_to_y` / `q_y_to_x` use members with usable decisions in **both** directions.
- q≥0.9 is a reporting threshold, requiring at least 90% usable requested members. Always retain the numerical q and counts.
- Complete-graph recovery counts only X→Y for coupled processes, and neither arrow for independent processes, with no unresolved evidence. The four-level plot divides the correct count by the **requested** member count.
- Unresolved evidence can coexist with a retained arrow. Do not merge unresolved, missing and genuinely absent arrows.
- Native-score summaries use all finite successful scores, not only significant members.
- Increasing chronology members improves Monte Carlo stability for a record. It does not substitute for independent process realizations or calibrate false-positive rates.

## Performance implementation

### Sampling sensitivity branch (02_4)

- The detailed comparison contract is in `agent/02_4_sampling_sensitivity_plan.md`. Use the existing three process forms and one process realization per form, exact ages, zero measurement noise and unchanged estimators. Ten A/B mask groups are sampling repetitions, not independent process histories or chronology members.
- Density controls use synchronous regular 20/40/80/160-year gaps. The irregular controls have exactly 351 observations over 29–1 kyr, including both endpoints: 175 gaps each of 60/100 years (CV 0.25) or 20/140 years (CV 0.75). Compare AA, BB, AB, BA at each fixed gap multiset.
- Old-sparse and young-sparse controls relocate the same long/short intervals, using paired Gumbel draws and log-weight slope +2/-2 in interval order. Keep X/Y synchronous (AA, BB) for this second stage. Save actual density by BP-age bins; the weight is not a prescribed linear density in physical age.
- Total default case count is 12 density + 240 irregular/asynchronous + 240 age-gradient = 492. All inputs are interpolated to the same 1401-point grid, without extrapolation. The 20-year control exactly recovers the original uniform input.
- Save complete native arrays, nested results and GC scientific fit fields as per-method JSON gzip files; the eight-row per-case edge table is the final commit. Retain explicit failures. Regenerated inputs, scientific settings, packages, implementation revision and numerical platform must match before resume.
- Read committed edge tables with pandas `float_precision="round_trip"` when resuming. The default CSV parser can change the last floating-point digit when rebuilding aggregate tables. This precision setting applies to result loading; the validated estimator input-loading path remains unchanged.
- Native-score paired differences require all component scores to be finite and successful. Report usable/requested counts and descriptive medians/IQRs. Correct graph recovery divides by requested cases, excludes unresolved evidence and never silently removes failed/pending cases.
- 02_4 uses the existing 32-worker setting, single numerical threads, native GPDC null cache and CCM neighbor-query reuse. Observe memory pressure during its run; do not repeat the completed 24-worker benchmark.

### Memory sensitivity branch (02_8)

The authorized 2026-11-09 design is in `agent/02_8_memory_sensitivity_plan.md`, with completed results in `agent/02_8_memory_sensitivity_results.md`. Keep the main 02/02_3 parameters and saved results unchanged. The independent branch compares X/Y memories of 500/200, 100/200, 500/40 and 100/40 years, with K=4 and fixed unforced background SDs. Derive its AR coefficients, per-step coupling and innovation SDs from these settings. Use the same 125-step warm-up and paired standardized innovations for all conditions.

Use three process seeds (20261009, 20261301, 20261302), linear/nonlinear/independent forms, two accurate-age references and ten existing level-4 C-C chronology members: 432 input pairs, 1728 full-setting method calls and 3456 directional rows. The exact same sampling masks and X/Y chronology arrays are reused across memories and seeds; do not regenerate them around truth or scale their correlation length with process memory. All four existing methods run on every input, including native CPU GPDC and optimized CCM L=300. Keep 32 outer processes and single numerical threads; do not repeat the 24-worker benchmark.

Compare chronology scores to each memory/seed's own irregular accurate-age reference. Keep per-seed results in the saved tables rather than treating chronology members as independent histories. The user's 2026-11-09 presentation revision removes seed-specific marker shapes from the main figures: use descriptive pooled proportions and score distributions, retaining all three histories and the 3-history x 10-shared-clock interpretation. Save complete native arrays and GC fits, retain failed/unresolved outcomes, and examine CCM raw versus screened results and all lag curves. The eight-row case edge table is the final commit. Validate reference inputs, pairing, native-output equality and resume guards with `agent/verify_memory_sensitivity.py` before production. The memory run used 02_8 only; its figures and data have their own stem directories. The subsequent readability revision uses requested-directional-input state fractions in Figure 02 (not q), keeps raw evidence separate from lag rejection and unresolved evidence, and shows saved per-input 100-library sign fractions beneath Figure 03 profiles. Its 0.9 diagnostic guide is not a cross-history support rule. Figure 04 explicitly expands the independent-row scale without changing units; the main curves and all numerical outputs remain unchanged.

The 2026-11-09 production run completed 432/432 pairs and 3456 successful directional rows in 1458.999 seconds with 32 workers; minimum sampled available physical memory was 8.90 GiB. The complete artifact audit checked 1728 native archives, all CCM arrays, GC fits and exact reconstruction of all six summary tables. All 2880 paired score changes and 2160 cross-memory contrasts are usable. Preflight report: `data/benchmarks/memory_checks/20261109T114552073978Z/report.json`; production audit, preservation and figure review: `data/logs/memory_sensitivity_20261109T114908564938Z/`. All 20,616 recorded pre-existing scientific files/results retain their hashes, sizes and modification times. The four new PNG/PDF figures have reviewed English captions. Three seeds and a shared clock set support conditional sensitivity conclusions; preserve the accurate-age CCM and LKIF exceptions in the short-Y condition.

### Common-driver branch (03)

The 2026-11-09 authorized fixed-memory design is in `agent/03_common_driver_plan.md`. Use only linear dynamics in this first version: Z/X/Y intrinsic memories 500/200/200 years, unforced background SDs one, Z-to-X and Z-to-Y long-run gains two, and X-to-Y gain zero (H0) or one (H1). Direct edges have one 20-year step delay; no feedback. The same 250-step warm-up and standardized innovations pair H0/H1. Do not add the previously discussed memory matrix or nonlinear variants to this first run.

Three process/clock blocks (20261009, 20261301, 20261302), two structures, two accurate-age references and ten members in each of z_only/xy_only/all_uncertain yield 192 inputs. One synchronous sampling mask uses the existing 27.5-year pre-rounding target and 10% jitter. C-type proposal SD is 40 years with fixed sample-order correlation length 1/60. Each record/block has an independently generated uncertain observed chronology and members around that observed center; reuse them across structures/assignments. Record actual and relative error diagnostics. No extrapolation, sorting or silently recentered clocks.

Run existing pairwise X/Y PCMCI+, LKIF, GC and CCM, plus three-variable conditional PCMCI+, LKIF and GC: seven native archives and 26 directional rows per input. Keep CPU GPDC/optimizer=None/10000 null draws, maximum lag five, GC order one and full CCM L=300 settings. Conditional GC uses native OLS single-equation F tests to match the existing pairwise SSR-F denominator, not the VAR system-level F as its reported test. PCMCI unresolved evidence is pair-local. CCM remains explicitly unadjusted and is not ranked as a conditional direct-edge estimator.

Exact unchanged pairwise evidence is copied for z_only from its irregular-accurate case and for all_uncertain from its matching xy_only case, after actual X/Y arrays and estimator seeds match. Run the 72 source cases before the 120 dependent cases using 32 outer workers and single numerical threads. Preserve all original scores, failures and native archives; do not refit just to change an unfavorable result. This gives 864 native calls and 1344 saved archives. A 26-row case table is the final commit. Validation and read-only output auditing are in `agent/verify_common_driver.py`.

The first four figures show the design/inputs, adjustment effects, complete X/Y decision states and separately labeled pairwise CCM evidence. Their captions are in `agent/figure_review/06_common_driver.md`. Keep three-history/ten-member interpretation, counts, unresolved states and panel labels. Old 02 data and figures are preserved. The original 03 run used the launcher with 03 only; its completed `run_status.json` remains authoritative.

The first 03 run completed 192/192 inputs and 4992 successful directional rows on 2026-11-09 in 2046.07 seconds; minimum sampled available physical memory was 10.13 GiB. Preflight is at `data/benchmarks/common_driver_checks/20261109T211251881262Z/report.json`. The production audit, source snapshot, resource log, preservation report and four-figure review are in `data/logs/common_driver_20261109T211815766544Z/`. All 1344 native archives, 480 exact reuse copies, complete CCM arrays, 3072 GC fit records and four reconstructed summary tables passed. All 24,305 pre-existing scientific files are unchanged. See `agent/03_common_driver_results.md` for findings and accurate-age exceptions. These are three process/clock blocks at one memory/gain/chronology setting; neither pooled member fractions nor correct XY subgraphs establish general false-positive rates or full three-node graph recovery. Do not restart a completed run to replace unfavorable or unresolved evidence.

The 2026-11-10 presentation revision uses **Common driver only** and **Common driver + direct X to Y** as system names. Keep historical h0/h1 storage identifiers unchanged. The original four figure definitions remain, with CCM rho explicitly described as reconstruction skill. Added figures 04 and 05 compare conditional-minus-pairwise exact XY recovery (percentage points) and method-wise exact XY recovery curves. Success uses the existing correct_xy definition, includes no unresolved evidence, and divides by requested cases. Keep systems and pairwise/conditional modes separate; CCM belongs only in the pairwise panel. Curves join categorical chronology assignments, not an increasing uncertainty-severity axis. Revision archive: `data/archives/common_driver_figures_before_revision_20261110T105738358740Z/`; checks: `data/benchmarks/common_driver_figure_checks/20261110T105738358740Z/`. This is a saved-result presentation change, with no causal refitting or changes to native data.

### Simulation-realization by Z age-scheme diagnostic (03_1)

The user authorized this compact follow-up on 2026-11-10; its design is in `agent/03_1_process_chronology_cross_plan.md`. Call the factors **simulation realization** and **Z age scheme**, avoiding the ambiguous shorthand process/clock in figures. Retain both existing systems and all original scientific settings. Only Z has uncertain ages. For each system, cross the three saved realizations with all three saved Z chronology packages, each containing its uncertain observed center and ten ordered members. X/Y keep accurate ages. Do not generate new dynamics, clocks or sampling masks.

There are 180 crossed uncertain cases plus 12 accurate references = 192. Reuse 60 original diagonal cases and all 12 references in full; only 120 off-diagonal cases require three conditional fits each (360 calls). Pairwise PCMCI+, LKIF, GC and CCM evidence is exactly reusable because X/Y and their estimator seeds stay unchanged. Keep all seven native archives and 26 directional rows per case (1344 archives, 4992 rows); a complete edge table is the final commit. Import unfavorable or failed evidence unchanged. Check exact diagonal reconstructions and actual CSV-parsed X/Y inputs before reuse, as well as source settings, packages, platform and source hashes.

Primary outcomes concern the X-Y relation only. Do not expand this branch into full-network recovery, Z-edge performance, memory or nonlinearity sweeps. Native Z-related fields remain preserved as part of complete estimator outputs. Save exact XY states, requested-denominator recovery, within-cell usable counts/support and paired native-score changes against each realization's accurate reference. Three-by-three displays do not represent nine independent histories; row patterns also retain each realization's original estimator seed.

Use 32 outer workers and one numerical thread, native CPU GPDC and existing null caches. The original 03_1 launcher selected that branch only; source 03 and all 02 results remain read-only. Preflight and final audits use `agent/verify_process_chronology_cross.py`. The two figures explain the crossing and show six 3x3 XY recovery matrices, with correct/10, unresolved counts, diagonal reuse outlines and each row's accurate-age baseline. No causal calculation is needed for later figure-only updates.

The run completed 192/192 cases on 2026-11-10 in 799.423 seconds, including 360 new conditional calls. Minimum sampled available physical memory was 11.26 GiB. The preflight is `data/benchmarks/process_chronology_cross_checks/20261110T123801313973Z/report.json`; the source snapshot, production log and resource samples are in `data/logs/process_chronology_cross_20261110T124145323892Z/`. See `agent/03_1_process_chronology_cross_results.md` for the complete X-Y comparison. Retain these outputs and their original source 03 evidence; do not selectively rerun unfavorable members.

The final artifact audit passed all 1344 native archives (984 byte-identical copies and 360 new conditional archives), 4992 successful directional rows, complete CCM/GC evidence and exact reconstruction of the three summary tables. The two PNG/PDF figures and all 54 crossed cells plus 18 accurate references passed review. Hashes, sizes and modification times of all 26,464 pre-existing scientific files are unchanged; the artifact, figure and preservation reports are in the production log directory. Unresolved evidence remains explicit and is not a failed computation.

### Common-driver chronology models (03_2)

The user authorized the Z-only type comparison on 2026-11-10; see `agent/03_2_driver_error_models_plan.md`. Keep both linear systems and their 500/200/200-year memories, histories, masks, grid, methods and estimator seeds unchanged. X/Y retain accurate ages. For each of C/L/E/J, cross three realizations with three complete Z age schemes, each with ten members: 720 uncertain cases plus 12 accurate references = 732. Copy all 180 C cases and 12 references from completed 03_1. Only 540 new L/E/J inputs require conditional PCMCI+, LKIF and GC (1620 calls); all pairwise evidence is exactly reused after input/seed checks. Save 5124 complete native archives and 19032 directional rows; 3504 archives are source copies.

Use the existing level-3 envelopes: C/L/E 40-year maximum, 0-30 kyr envelope bounds, E curvature 2; J five smooth 10-year peaks with 2.8-kyr intervals at the established first five fractions of 29-1 kyr (50% defined-interval coverage). Preserve the original Z observation and ensemble seed formulas across types and correlation length 1/60. Draw around each type's uncertain observed center, reject whole unordered curves, and never repair or rescale ages. Reconstruct C exactly before copying it. Shared candidate seeds do not guarantee matched accepted members if type-specific rejection histories differ. Check coverage without extrapolation and retain acceptance/observed/member-error/effective-spread diagnostics.

This compares complete prescribed error settings, not pure shape at equal RMS. Do not add equal-RMS rescaling, uncertain X/Y backgrounds, memory/nonlinear sweeps or full-network recovery to this run. X-Y states and native-score changes stay type-specific, with requested denominators and explicit unresolved/failed outcomes. The three figures show envelopes and actual error, aggregate correct/30 matrices, and the separate three scheme counts per cell. The original 03_2 launcher selected that branch only, with 32 workers and one numerical thread. The current launcher subsequently moved to the authorized 04 TBK pilot. Preflight/final checks use `agent/verify_driver_error_models.py`; old 02/03/03_1 sources and outputs are read-only.

Preflight passed at `data/benchmarks/driver_error_models_checks/20261110T133558073115Z/report.json`. The initial production invocation in `data/logs/driver_error_models_20261110T133847522941Z/` was accidentally closed by the user. A read-only audit found 297 complete cases, 435 pending, and 7722 successful committed directional rows; input/source/parameter/platform checks passed. The user requested continuation. Recovery started at 13:54:03 UTC on 2026-11-10 in `data/logs/driver_error_models_resume_20261110T135403003994Z/`, loading all 297 complete cases without refitting. Check the current branch run_status and this recovery log, not the interrupted invocation, and do not start a duplicate. The final audit must additionally preserve hashes, sizes and mtimes of the 3267 files belonging to those pre-recovery committed cases. No scientific code or parameter changed for recovery.

Recovery completed all 732 cases at 14:36:22 UTC on 2026-11-10. Its 2535.030 seconds (42.25 minutes) describe the resumed invocation only; the initial monitor recorded about 10.2 minutes before interruption. All 19032 directional rows are successful, with 5124 complete native archives (3504 exact source copies and 1620 final new conditional results). Artifact, figure and both preservation audits passed in the recovery log directory. All 28595 pre-existing scientific files and all 3267 files belonging to the 297 pre-recovery commits retain their hashes, sizes and modification times. Three PNG/PDF figures have reviewed captions in `agent/figure_review/08_driver_error_models.md`; findings and baseline exceptions are in `agent/03_2_driver_error_models_results.md`. Minimum sampled available physical memory was 5.05 GiB and peak sampled system commit was 60.43 GiB, both whole-machine measurements. The completed run must not be restarted to replace unfavorable or unresolved evidence.

### TBK physical-model branch (04)

The accepted design is `agent/04_tbk_chronology_plan.md`; completed findings are in `agent/04_tbk_chronology_results.md`. The annual cell-12 core and process noise are preserved in `toolbox/tbk_model.py`. Default 10-kyr migration and extended 39-kyr arrays match the archived original kernel exactly. The 39-kyr experiment retains model times 9000..39000 and samples T1/T4 every 20 years, then uses the fixed 29-1 kyr, 1401-point grid. Synthetic BP ages do not identify historical DO events. Prescribed q has no temperature feedback; endpoint arrows are not full-network recovery or a calibrated AMOC-feedback test.

The run completed 123/123 inputs on 2026-11-11 in 1460.470 seconds (24.34 minutes), using 32 workers; minimum sampled available physical memory was 6.892 GiB. All 492 native archives, 984 directional rows and seven summary tables passed exact audit. There are no native or lag-scan failures, but 198 method-case states have unresolved evidence (123 PCMCI+ and 75 CCM). All 36,677 pre-existing scientific files retain hashes, sizes and modification times. Preflight: `data/benchmarks/tbk_chronology_checks/20261111T165732295477Z/report.json`; production source, logs and audits: `data/logs/tbk_chronology_20261111T171153343045Z/`. The inspection notebook and four production figures have separate captions in `agent/figure_review/09_tbk_chronology.md`.

On 2026-11-13 the user requested a more visible process-noise illustration. Only panel c of inspection Figure 00 uses a separately integrated SD=0.05 trajectory against SD=0. The prepared chronology inputs, panels a/b, remaining inspection figures and completed 04 production results retain SD=0.005. Keep the illustrative amplitude explicit and separate from production settings. The archived notebook's ten runs are nine SD=0.005 realizations and one noise-free reference, not ten noise levels. The Chinese local reading guide is `Note/04_tbk_figure_guide.md`.

Preserve accurate-age limitations: PCMCI+ retains both lagged arrows but also unresolved contemporaneous links; LKIF retains only North-to-South in all three accurate cases; raw CCM is strong and bidirectionally significant but lag-screened states are unresolved. GC retains both directions in all accurate cases, yet all 123 VAR(1) fits reject residual whiteness, including accurate ages. Do not call its subsequent bidirectional-retention curve an accuracy ranking. The 20-year uncertainty level adds LKIF bidirectional evidence in 30/30 inputs; this is not proof that dating errors improve causal correctness. Keep raw flags, unresolved states, requested denominators and each history's baseline. Do not rerun unfavorable results or tune histories after seeing them.

### TBK excitation diagnostic (04_1)

The authorized 13 November 2026 TBK follow-up is specified in `agent/04_1_tbk_excitation_plan.md`. It has twelve accurate-age inputs: the same three histories at process-noise SD 0.005/0.05/0.5 and each original-noise history minus its same-parameter noise-free response. Keep annual integration, physical coefficients, forcing, the 39-kyr duration/9-kyr startup, 20-year sampling and the 29-1 kyr analysis window fixed. Noise amplitudes share standardized innovations. The subtracted response is a model-assisted diagnostic, not real-proxy preprocessing. Reuse the three original accurate inputs and twelve full native archives only after exact input/settings/source checks. Run all four original methods on nine new inputs, plus native PCMCI+/ParCorr on the three original inputs as a separately labeled diagnostic; production GPDC remains unchanged. There are 51 native archives and 102 directional rows, with 39 new fits and 12 source copies. Use the 32-worker limit and single numerical threads. Preserve 04/04_0 sources and results, retain unresolved states separately from directional arrows, and do not interpret the figure's resolved-bidirectional fraction as method accuracy. Information-set, sampling and renewed chronology experiments remain deferred.

The 04_1 run completed 12/12 primary inputs and 3/3 baseline ParCorr controls on 2026-11-13 in 71.771 seconds with a 32-worker limit, existing GPDC caches and 12 copied native archives. All 51 native archives and 102 directional rows passed audit; no fit failed. Six scientific and seven figure tables reconstruct exactly, complete CCM arrays and GC fits are retained, and all 1,026 inventoried prior TBK files retain their bytes, sizes and modification times. Minimum sampled available physical memory was 10.004 GiB. Preflight and final-source acceptance are in `data/benchmarks/tbk_excitation_checks/20261113T121452115736Z/`; production, artifact/preservation audits and five-figure PNG/PDF review are in `data/logs/tbk_excitation_20261113T121730984237Z/`. See `agent/04_1_tbk_excitation_results.md` and `agent/figure_review/10_tbk_excitation.md`.

The fixed-memory diagnostic changes the dense baseline GPDC graph: tenfold noise retains only North-to-South with contemporaneous conflicts, hundredfold noise retains neither direction, and the known-forcing-subtracted response retains North-to-South in only one history. Baseline ParCorr retains only North-to-South in all three, with one contemporaneous conflict. CCM skill drops from about 0.999 to 0.004-0.098 after forcing subtraction, with no screened arrows; short-term autocorrelation remains substantial. These support forcing and conditional-test sensitivity, not a claim that more noise recovers feedback or that intrinsic memory is irrelevant. Proposed ParCorr checks on the three subtracted-response inputs and full-state/sampling comparisons have not been authorized or run.

### GPDC null cache

`prepare_gpdc_null(filename, sample_sizes=[T], sig_samples=..., seed=...)` calls the native generator once in the parent. Pass its returned path to `run_pcmciplus(..., null_dist_filename=...)`.

Current effective T is 1401−2×5=1391. Use one T and seed per file to preserve the native random sequence. Cache metadata checks T, seed, sample count, relevant packages and platform. Do not reuse an unrelated null distribution. PCMCI deep-copies its conditional test, so simply retaining an external GPDC object is not sufficient.

Every member still receives its own conditional GP regressions. Do not cache residuals or conditional-test results across different series.

### CCM neighbor-query reuse

`run_ccm` calls `filter_ccm_by_lag`, which uses `scan_ccm_reusing_neighbors` in `toolbox/causal.py`. One native object is retained per lag. The outer loop visits libraries, and each library is used for all lags before its query cache is cleared. Exact array bytes, shapes, dtypes and query k identify the cached SciPy results; native random sampling, exclusions, weights, predictions and correlations still execute unchanged. All 100×31×2 rho values remain available. The original Tp=0 and swap-null stage is unchanged.

The helper is validated for PyEDM 2.5.7 and temporarily replaces that module's KDTree factory, restoring it in `finally`. It supports sequential calls inside outer worker processes, not concurrent threads. Per-library/per-lag errors remain explicit, with failed curves left as NaN; a failed lag scan must not discard successful Tp=0/null evidence. The standalone helper exposes cache statistics, but those implementation counters are not mixed into the existing scientific result dictionaries. `02_7` retains the independent pre-migration scientific functions and settings. The 2026-11-08 figure revision changes only style imports and its final plotting cell; the verification script checks scientific cells against the freeze while allowing those explicit presentation changes.

Do not replace this with 31 calls at `sample=100`: that does not preserve the required individual paired library curves and random streams. Do not change embedding, nulls, library count or lag range as a performance shortcut.

### Multiprocessing and scale

- 02 and 02_3 use 32 outer workers on this 14900KF, selected from the measured 16/24/32-worker comparison below. Recheck throughput and memory pressure if the workload or hardware changes.
- Use one BLAS/OpenMP thread per worker. The formal launcher sets these variables before importing numerical packages.
- Keep an outer process pool; do not substitute threads. The current PyEDM sequential path uses mutable module-level worker data.
- Keep PyEDM internal parallelism off. Repeated nested pools can waste resources and memory.
- Keep worker functions at module scope and preserve `if __name__ == "__main__"` for Windows spawn.
- Benchmark 16/24/32 on the target PC with 96 identical inputs. Observe memory/commit and paging. Select the fastest configuration with acceptable memory pressure; do not pick 32 from the logical CPU count alone.
- Main case counts are derived from configuration:02 is 16×3×N;02_3 is 4×16×3×N, including imported level 2. N = 10 gives 480 and 1920 pairs respectively.

## Windows environment and launch procedure

Use 64-bit Python 3.11 and Git for Windows. Run commands from the project root. Do not copy `.venv` from Mac, and do not install unconstrained latest scientific packages.

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r agent\requirements.txt
```

`agent/requirements.txt` fixes the scientific versions. NumPy 1.26.4 is intentional for LKIF compatibility. LKIF is installed from Git commit `95c088ef8cc5e5dab2c6279158e98338ee51124a`; do not replace it with the known problematic PyPI wheel. The local Windows `.venv` was installed with 64-bit Python 3.11.7 on 2026-10-26; pinned requirements, the LKIF commit and `pip check` were verified. Notebook installation in the production `.venv` remains optional and unverified; the separate `.venv-gpdc` notebook environment is described below.

For notebooks, optionally install and register the project kernel:

```powershell
.\.venv\Scripts\python.exe -m pip install -r agent\requirements-notebook.txt
.\.venv\Scripts\python.exe -m ipykernel install --user --name chronology-causal --display-name "Python 3 (project .venv)"
```

Recommended first Windows workflow:

1. Run `.\.venv\Scripts\python.exe agent\verify_optimizations.py`. This is a few-case numerical check, not a full experiment.
2. In 02 set `PREPARE_ONLY=True`; leave `RUN_BASELINES=True` and `RESUME=True`. Set the launcher's `FOLLOWUP_EXPERIMENTS=[]` during preparation. Run the launcher below. It computes the small baseline set and prepares all chronology inputs, then stops before the ensemble fits.
3. Run `.\.venv\Scripts\python.exe agent\benchmark_workers.py`. Choose the worker count from measured throughput and memory behavior. The96-case benchmark is real causal computation and may take several minutes per worker setting.
4. Set `PREPARE_ONLY=False` in 02 and the selected `N_WORKERS` in both 02 and 02_3; keep matching chronology counts. With `EXPERIMENT="02_check_causal_pipeline.py"` and `FOLLOWUP_EXPERIMENTS=["02_3_check_uncertainty_escalation.py"]`, launch once to reuse the prepared baseline, complete level 2, then calculate the other levels automatically.
5. To launch only 02_3 after a completed local 02, set `EXPERIMENT="02_3_check_uncertainty_escalation.py"` and `FOLLOWUP_EXPERIMENTS=[]`. It imports the completed local level-2 results.

```powershell
.\.venv\Scripts\python.exe agent\run_experiment.py
```

The launcher uses its own interpreter, creates a timestamped log per script under `data/logs/`, and streams output to the terminal. The current launcher selects only 04_1 for the authorized TBK excitation diagnostic; 04, 03_2, 03, 03_1 and 02_8 are complete. The L=300 main chronology pilot is complete; no main-study calculation is restarted. Any nonzero exit stops the remaining sequence. Keep that terminal open for a manual launch and run only one instance per experiment output folder. Ctrl+C requests termination of the child process group; completed member files remain. If Windows cannot deliver the console signal, the launcher falls back to `taskkill /T /F` for the live child tree. Check `run_status.json`, not an old log selection, for the latest saved progress. On Windows the launcher temporarily prevents automatic system sleep across the whole sequence, then releases the request on exit. Display sleep remains allowed; no persistent power setting is changed. If this API fails, the launcher warns and continues.

The strength script and lag notebook are optional branches; they need not be rerun before every 02 experiment.

## Resume and fresh experiments

Before the 2026-11-07 migration, 28 current source/support files were frozen and SHA-256 verified in `data/archives/code_freeze_20261107T200542517660Z/source.zip`; its manifest records each file and the uncommitted Git state. This is a source archive, not a copy of generated results. The migrated 02/02_3 revision is `native_null_cache_and_ccm_neighbor_reuse_v2`; 02_4 uses `fixed_gap_paired_sampling_ccm_neighbors_v2`. Preserve existing data and figure folders before a fresh L=300 run, then run 02 before 02_3. Old metadata is deliberately rejected. The optional 02_1 still has `RECOMPUTE=True`, so preserve its old output folders before launching it as well. No new production ensemble was launched during migration; N remains 10 and workers remain 32.

The user authorized the fresh pilot after migration. On 2026-11-07, the previous 02/02_3 data and figure folders (17,517 files) were moved to `data/archives/pre_ccm300_pilot_20261107T202152109Z/`, retaining their relative paths; every file hash was verified before and after moving. 02_1, 02_4 and legacy macOS outputs were left in place. The new sequence started at 20:24:52 UTC with L=300, N=10, 32 workers and CPU GPDC. Its original supervisor status, launcher log, source snapshot and five-second resource log are in `data/logs/pilot_ccm300_20261107T202452271853Z/`. That invocation stopped at 21:46:58 UTC after a WinError 5 replacing 02_3's run_status.json. Actual committed results were intact: 02 had 480/480 pairs; 02_3 had 1440/1920 across levels 1/2/3, with all 480 level-4 pairs still pending.

Recovery started at 23:13:49 UTC on 2026-11-07 under `data/logs/pilot_ccm300_resume_20261107T231349346354Z/`. Check its pilot_status.json and launcher_console.log for current progress; do not start a duplicate. The only scientific-script change adds bounded, five-second retries for Windows access/sharing errors at the three atomic output replacements in 02_3. Real file-lock tests reproduced WinError 5, recovered after release and preserved both files on timeout. AST comparison confirms unchanged scientific calculations/settings and toolbox source. Resume keeps the original configuration/input guards and all committed fits, including any recorded failures. The supervisor schedules a complete artifact audit and SHA256/mtime checks of 12,960 pre-existing per-member files after completion. Both invocations and the verification report `tmp/windows_replace_check_20261107T231314106838Z/report.json` remain available. These are background launches; completion must be confirmed from current status/audit files.

- `RESUME=True` resumes complete per-member outputs only when the scientific configuration, packages, implementation revision, numerical environment and regenerated inputs match.
- A completed member can include an explicitly recorded failed fit. Keep it as evidence; do not automatically rerun only unfavorable or failed scientific results.
- A member without its complete committed result table is incomplete and may be recalculated. Partial native files are not valid complete votes.
- Worker count and run-control flags can change without redefining the scientific experiment. Changing chronology count, dynamics, chronology settings, methods, package versions or platform requires a fresh run.
- To increase ensemble size, preserve/rename the old output folders, set matching counts in 02 and 02_3, and run 02 first. Do not edit saved metadata to force an old result through a guard.
- For 02_3 figures only, set `RUN_CAUSAL=False` with matching saved settings. This does not recompute estimates. Pending levels must remain explicitly partial.
- This is a simple same-run resume mechanism, not a general cache search across historical experiments.

## Files to synchronize

Git should contain scripts, source notebooks without embedded generated output, toolbox, the maintained `agent/` support scripts and requirements, this guide and the small directory README files. This repository does not require a transferred Mac dataset to start the Windows pipeline.

Keep `.venv/`,`.venv-gpu/`, `tmp/`, generated data, figures, caches, logs and historical archives out of Git. Ignoring a file is not deleting it. If a cloud folder is used instead of Git, configure its exclusions separately; `.gitignore` does not control a cloud-sync application.

The original repository had no tracked files or commits at migration time. Check `git status` before creating the first curated commit. Do not upload credentials, local environment directories or the legacy experiment archive by accident.

All temporary agent work belongs in `tmp/` locally. The Windows agent must be able to understand/run the project from this guide and checked-in source; no personal skill installation or `/Users/...` path is a runtime requirement. Readable comments and notebook prose should be English; Chinese documentation belongs only in the Git-ignored `Note/` directory. Discussion with the user can use their chosen language.

Use plain functions, explicit parameters and readable intermediate calculations. Do not introduce a class/configuration framework merely to remove a few repeated lines. Give every configurable scientific parameter a comment explaining its meaning and units. Notebook sections should have a title and one or two plain English sentences. Keep plotting in root analysis scripts/notebooks, and agent scratch calculations under `tmp/`. Future real-climate inputs will need an explicitly chosen storage/ignore policy; the current broad data exclusion is for this synthetic migration.

## Optional CUDA exploration

Production GPDC/sklearn GPR, PyEDM KDTree, LKIF and statsmodels currently run on CPU. Installing CuPy does not switch these packages to CUDA.

If a compatible NVIDIA driver and CUDA Toolkit 12.x are installed, use a separate environment. The pip wheel requires the CUDA dependencies described in the [CuPy 13.6 installation guide](https://docs.cupy.dev/en/v13.6.0/install.html); installing the Python wheel alone is not a complete CUDA setup.

```powershell
py -3.11 -m venv .venv-gpu
.\.venv-gpu\Scripts\python.exe -m pip install -r agent\requirements-gpu.txt
.\.venv-gpu\Scripts\python.exe agent\explore_cuda.py
```

The optional requirements pin CuPy 13.6 to remain compatible with the existing NumPy 1.x line. The probe checks CUDA availability, compares a float64 RBF-covariance/Cholesky solve and records timings including transfers. **Its speed ratio is not an end-to-end PCMCI/CCM speedup.** It does not implement a new causal estimator.

Keep this experiment independent of CPU production. A real GPU backend would need a separate validated implementation from an existing package, full numerical/decision comparisons and a measured batch-throughput benefit. Do not replace methods or lower precision just to use the GPU.

The small `02_5_compare_gpdc_backends.ipynb` uses `.venv-gpdc` and `agent/requirements-gpdc.txt` instead of CuPy. It reuses the three original-uniform 02 records and matching T=1391 GPDC null cache. Six fixed conditional tests compare current GPDC, stock GPDCtorch, and a notebook-local aligned GPyTorch regression on CPU/CUDA. It does not run PCMCI+ search or change production. Select `.venv-gpdc/Scripts/python.exe` in the IDE and run from the project root. The root notebook has no generated outputs; its executed copy, tables and full regression arrays are under `data/processed/02_5_compare_gpdc_backends/`.

## Validation already performed and limits

The Mac implementation check used two fixed 1401-row inputs (linear and independent). Cached/uncached PCMCI arrays and graphs were exactly equal; fresh/reused-object CCM matched the full 100×31×2 rho array, nulls, seeds and arrow results exactly. Cache mismatches were rejected.

Observed single-process times were approximately 18.7→4.9 seconds and 17.7→3.8 seconds for PCMCI+, with a one-time 16.8-second null preparation. CCM changed 49.2→47.0 seconds and 51.2→49.9 seconds. These checks suggest about 1.3× combined improvement, **not** a 4× improvement of the whole workflow.

The earlier setup, resume, failure-retention and metadata-rejection checks used small mocked workflows. A separate real-estimator integration check used two spawned processes and reduced test-only settings: both scripts saved eight successful directional rows. Launcher logging, exit codes and temporary sleep-control lifecycle were checked locally or with mocks. Source and notebook-cell syntax checks passed; notebook execution was not repeated.

Windows acceptance on 2026-10-26 passed the complete two-case optimization comparison, with exact numerical equality, in `data/benchmarks/optimization_checks/20261026T214439840630Z/report.json`. Both accurate-age baselines produced all 48 directional rows with status `ok`, and all 480 level-2 chronology inputs were ready.

The same 96 prepared input pairs gave these Windows worker measurements:

| Workers | Wall time (minutes) | Input pairs/hour |
|---|---:|---:|
| 16 | 9.65 | 596.96 |
| 24 | 7.90 | 728.78 |
| 32 | 7.50 | 768.20 |
| 32, repeat | 7.56 | 762.12 |

All 88,182 fingerprinted components across the 96 pairs matched exactly between worker counts and the completed 32-worker repeat. Minimum available physical memory at 32 workers was about 11.9 GiB, with no sustained paging pressure observed. The user stopped the second 24-worker measurement; it is incomplete and is not additional throughput evidence. These measurements selected 32 workers for production.

The production sequence completed with `PREPARE_ONLY=False` and 32 workers: the saved status files report 480/480 pairs for 02 and 1920/1920 for 02_3, including 480 reused level-2 pairs, with no recorded run errors. Read each experiment's current `run_status.json` for subsequent progress. Sequential launch, signal-failure fallback and sleep-control lifecycle passed mocked checks. A real benchmark process tree was terminated with `taskkill`; native Windows Ctrl+C handling still needs acceptance. Execution of the original notebooks and the optional CuPy probe remains unverified; the isolated Torch notebook check below is separate.

On 2026-11-07, the isolated GPDC notebook completed all 24 case/backend combinations on the RTX 4080 (Torch 2.5.1+cu124, GPyTorch 1.13). Aligned float64, zero-mean, unit RBF/noise, untrained dense GP fits matched the current GPDC regression vectors to at most about 2.4e-14, with exactly equal scores and p values for all six fixed tests. Default trained float32 GPDCtorch changed the fits and scores. Single-pass median time per fixed test was about 0.28 s for current GPDC and 0.10 s for aligned CUDA; these are not full PCMCI+ or 32-worker throughput measurements. The environment passed `pip check`, the source notebook passed clean-kernel execution, and the figures were visually checked. Full production migration remains unvalidated.

Retain the original scientific outputs if a future optimization fails the numerical comparison. Investigate discrepancies rather than changing significance thresholds or dropping difficult cases.

The 2026-11-07 toolbox migration passed `agent/verify_ccm_neighbor_reuse.py` against the pre-edit source freeze. Five complete CCM cases (three uniform forms and one nonlinear level-4 input at L=300, plus one independent level-4 input at L=1119) matched every DataFrame, array and scalar exactly, including original scores/nulls and all lag curves. Each comparison covered five DataFrames, eight arrays and 130 scalar leaves. Two Windows-spawned workers reproduced the same complete results. Single-library/single-lag failure retention, original evidence preservation, oversized libraries, repeated states, bounded cache use and interruption restoration passed. The report is `data/benchmarks/ccm_migration_checks/20261107T201205780934Z/report.json`; 02_7 was byte-identical to the freeze at migration acceptance; the later style-only revision preserves its scientific cells. Same-input full-CCM times were about 26.5--35.4 s before and 14.1--18.1 s after migration, single-pass sequential measurements rather than a 32-worker benchmark. Old configuration/source guards were checked and reject the retained L=1119 runs. GPDC/GPU migration remains pending. The subsequent L=300, N=10 pilot completed successfully, as recorded below.

Sampling-branch preflight on 2026-11-04 passed: 492 unique inputs, fixed gap multisets/endpoints, shared masks, and exact recovery of the existing 02 original inputs by the 20-year control. Two full-setting cases (linear CV 0.75 AB; nonlinear CV 0.25 old-sparse AA) matched exactly between direct execution and two spawned workers across all eight complete native payloads and directional rows, excluding runtime only. The check included the full CCM 100×31×2 lag arrays, 100×2 nulls and GC regression fields. Sixteen isolated checks also passed for failure retention, committed outputs, resume/configuration guards and serialization errors. Reports are under `data/benchmarks/sampling_checks/20261104/`. The 32-worker production run completed all 492/492 pairs in 3739.757 seconds (62.329 minutes), with all 3936 directional rows successful and all 1440 paired score contrasts usable. Minimum available physical memory was 9.67 GiB. The full artifact audit passed, including all 1968 native archives. After fixing last-digit CSV result-loading drift with round-trip parsing, an actual resume preserved the bytes and modification times of all 2460 committed member/native files. The final comparison is `artifact_audit_20261104T215527447629Z.json` in the same reports directory. The earlier audit retains the full production runtime; the current run_status records the latest resume invocation.

## Figure styling and saved-result redraws

On 2026-11-08, the user rejected the broad figure restyling because useful information was missing. Restore and retain the original research colors, fonts, line widths, titles, numerical annotations, legends and explanatory footers, with lowercase panel labels and verified scientific-label corrections. The 30 exported figures retain their original layouts; the 10 definitions in 01/02_2 are updated without execution. `toolbox/plotting.py` supplies only panel labels and PNG/PDF export, without global style overrides. Do not split dense figures or simplify away helpful information without a later user request. Captions and source/provenance notes remain in `agent/figure_review/`.

Run `.\.venv\Scripts\python.exe agent\redraw_figures.py` from the project root to redraw the 30 saved figures without starting causal fits or generating chronologies. The three component `agent/redraw_*_figures.py` scripts load committed artifacts and call plotting functions or explicitly selected notebook plot cells. Keep old 02_1/02_4 L=1119 provenance; their appearance does not adopt the current L=300 estimator default. Previous figures and source are archived in `data/archives/pre_figure_style_20261108T121553339241Z/`. All 20,606 scientific data files retained their hashes and modification times. The earlier style check is retained as historical evidence; the current rollback check is `data/benchmarks/figure_rollback_checks/20261108/report.json`.

The L=300 chronology pilot has completed: 02 has 480/480 pairs and 02_3 has 1920/1920 including reused level 2. The recovery supervisor reports complete, its full artifact audit passed, and its 12,960 pre-recovery member files were unchanged. Do not interpret the older failed invocation as current progress.
