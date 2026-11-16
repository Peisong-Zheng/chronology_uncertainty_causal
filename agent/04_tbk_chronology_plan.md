# TBK chronology experiment plan

Authorized implementation specification, 10 November 2026. The user accepted the proposal and requested consideration of record-length alignment. The accepted revision retains 30 kyr after a 9 kyr startup, sampled at 20 years, matching the earlier experiments. Implementation, validation and the 123-input pilot are complete. Full output and preservation audits passed; findings are reported in [04 results](04_tbk_chronology_results.md). The [discussion record](research_direction_and_tbk_notes.md) preserves the scientific motivation and preceding decisions.

## Scientific question and scope

**Can chronology uncertainty change the apparent direction of influence between the southern and northern records of a climate model with reciprocal temperature interactions?** A particular concern is apparent one-way influence in a system that contains feedback. This is an investigation of causal evidence and its physical interpretation, not a claim that every pairwise method estimates the full four-box equation graph.

The user requests retention of the original four physical TBK boxes, process noise and dynamics. Analyze the southernmost T1 and northernmost T4 records. Both records have independent chronology uncertainty in every nonzero level; the earlier southern-only/northern-only design is superseded. Accurate dates remain an essential analytical reference, not a realistic-observation scenario. Hold process noise fixed and add no measurement noise. Additional three-variable mistiming experiments remain deferred.

The first version varies chronology amplitude only. It does not vary envelope type, cross-record error correlation, process-noise strength, coupling, forcing periodicity or memory. It does not fit a four-node causal graph or introduce temperature feedback into the prescribed transport forcing.

## Proposed files and migration contract

| Destination | Responsibility |
|---|---|
| `data/archives/tbk_original_<UTC>/One_D_model_TBK.ipynb` | Complete original notebook, including all cells and outputs, plus its checksum and a short source-variant note. Verify the archival copy before removing the exact root source file. |
| `toolbox/tbk_model.py` | Model integration only: original numerical equations, update order, forcing schedule, noise generation and boundary operations; return states, forcing, fluxes, noise and parameters. No dates, causal fits or figures. |
| `toolbox/tbk_observation.py` | Small shared preparation module for extracting the two records, fixed sampling, synthetic age coordinates and ordered chronology ensembles. Reuse `toolbox/observation.py`; do not implement a second chronology kernel. |
| `04_0_inspect_tbk_records.ipynb` | Short, clean-kernel notebook calling the model and chronology preparation, with visible editable settings and diagnostic figures. No causal fitting by default. |
| `04_check_tbk_chronology.py` | Root experiment: prepare inputs, call the existing four causal methods, save complete outputs, summarize and plot. Parameters remain visible near the top. |
| `agent/verify_tbk_migration.py` | Maintained same-input migration and experiment-preflight checks; no new general framework. |

The separate chronology module is proposed so the inspection notebook and production script use the same preparation functions. Chronology choices and examples remain visible in notebook cells; they do not enter the original model module. This interprets the user's requested separation in a way that avoids maintaining two copies of notebook-local chronology functions.

Archive first, then migrate. Update old notebook links when the archive move actually occurs. Production must not import archived code or `tmp/`. Maintain original authorship and provenance in the extracted model.

The core should be moved with minimal changes: keep the annual step, loop ordering, arithmetic ordering, ghost-box storage and replacement, gain=0, prescribed forcing update at the end of the step, and original random-array shape. Add a function boundary, explicit parameters, reproducible random input and returned arrays. Remove plotting from the extracted calculation. Do not vectorize/rederive the integration or correct indexing during migration. Preserve the original temperature increment interpretation of the noise, including the absence of sqrt(dt) scaling.

Some legacy diagnostic arrays, including target and global/area means, leave their final entry unfilled. Preserve those arrays exactly during migration. Any complete diagnostic recomputed for plotting must be separately named and derived explicitly from the saved states; do not silently repair legacy values in the model output.

For reproducibility, a local legacy NumPy RandomState or supplied pre-generated noise can preserve the legacy random stream without modifying global random state. Validate against the archived cell using identical draws, including ghost-box draws. Explicit seeding does not establish that an unseeded historical output has been reproduced; the acceptance target is the same input on this Windows platform.

## Select and document the original model variant

Two versions coexist in the notebook:

| Zero-based cells | Southern diffusion coefficient | Slowest temperature decay mode, years |
|---|---:|---:|
| 3 and 40 | 0.0005 | approximately 1500 |
| 8 and 12 | 0.005 | approximately 154 |

**Recommended primary source: cell 12, with 0.005.** It generated the saved `Temps` arrays used in most subsequent analyses. Cell 40 later changes a separate `Temp` array, so selecting the last cell would not reproduce those `Temps` analyses. This recommendation is based on provenance, not on favorable causal results. Preserve both original variants in the archive and make the existing diffusion parameter explicit; do not run a two-version sensitivity matrix in the first experiment.

The decay times above were calculated from the existing linear update matrix, without executing a model trajectory. For the recommended version, the four temperature modes are approximately 2.8, 9.2, 99.4 and 153.9 years. These are relaxation modes, not observed-record ACF times or prescribed cross-variable delays. The externally imposed complete switching cycle is approximately 1500 years, consisting of 1000 off and 500 on.

The global mean restoring term creates reciprocal mathematical dependence between all temperature boxes, including the endpoints. Intermediate boxes and prescribed forcing remain unobserved in the two-record analysis. A full-state annual equation edge, a coarse-resolution predictive relation and a total physical response are different targets. No single northern-to-southern arrow is designated as the unique truth.

## Model duration, observations and reference data

Keep model defaults corresponding to the original 10000-year integration for migration checks. For the scientific experiment, **use a 39000-year integration**, discard the first 9000 years for analysis, and retain coordinates from model time 9000 to 39000 inclusive: 30001 annual points over 30000 years. Preserve the existing `time=1..nsteps` convention and verify exact indices; do not silently shift forcing or output labels by a year.

This is an explicit experiment-duration setting, not a change to the model equations. Nine thousand years is six full forcing cycles and more than five of the longest decay times in either retained original variant. It reduces initial-condition transients without subtracting a fitted trend or modifying the forcing. The original 20 prescribed event intervals end at 30000 years. Extend that event table to 26 intervals using exactly the original 1000-off/500-on cycle, keeping the first 20 intervals and update order unchanged. This necessary duration extension is verified separately from the exact default 10-kyr migration. The discarded interval is not used as causal data.

Use the full retained 30 kyr as the observation support, with **synchronous observations every 20 years**, sampled from the annual solution without averaging or added value error. Define synthetic age as `(39000 - model_time) / 1000`, giving 30 to 0 kyr BP for this artificial record. The origin is a bookkeeping choice, not attribution to actual historical DO events. Rows remain in forward physical-time order, with age decreasing.

Analyze a fixed **29 to 1 kyr BP** grid at **20-year spacing**, giving **1401 points**. The one-kyr edge margins provide room for uncertain chronologies; they do not guarantee coverage, which must be checked. At zero uncertainty, observations and the analysis grid coincide exactly, so only one accurate-age estimator input per history is needed. The annual solution is a model diagnostic, not a second differently resolved benchmark mixed into the main denominator.

The 20-year observation interval and record length align with the other experiments. Alignment is useful for keeping sample count and the fixed CCM library fraction comparable; it is not a physical requirement. Annual integration remains unchanged, while 20-year observations do not resolve the fastest 2.8-year relaxation mode. The annual model graph is therefore not used as the unique statistical graph truth at the observed resolution.

## Chronology levels and their meaning

Use **C-C** constant proposal envelopes and independent southern/northern errors. The same nonzero level applies to both records. Proposed levels are:

| Level | Per-record observed-center proposal SD, years | Per-record member-around-center proposal SD, years | Idealized inter-record member error SD, years |
|---|---:|---:|---:|
| Accurate | 0 | 0 | 0 |
| 1 | 20 | 20 | 40 |
| 2 | 50 | 50 | 100 |
| 3 | 100 | 100 | 200 |
| 4 | 200 | 200 | 400 |

Four nonzero levels span modest timing uncertainty through a several-hundred-year relative mismatch, compared with approximately centennial temperature relaxation and a 500-year forcing-on duration. They are an idealized scale experiment, not calibrated uncertainty bounds for a specific real archive. Do not choose or change levels according to where a causal algorithm begins to fail.

The two layers are explicit. For record i and true age a[k],

$$
\widehat a_{i,k}=a_k+\sigma z^{\mathrm{obs}}_{i,k},\qquad
a^{(m)}_{i,k}=\widehat a_{i,k}+\sigma z^{(m)}_{i,k}.
$$

Before any ordering/coverage conditioning, independent layers give SD sqrt(2)*sigma for a single member relative to truth, and **2*sigma** for the difference between two independently uncertain members. The observed centers' relative SD is sqrt(2)*sigma. These are pointwise ensemble calculations, not the realized temporal RMS of any particular record, confidence limits or an assurance about accepted distributions.

Use a squared-exponential age-error correlation with a proposed **1000-year physical length scale**, held constant across levels. On the full 30-kyr synchronous sample support this corresponds to normalized sample-order length 1/30 in the existing implementation. It describes smoothly varying chronology distortion, not independent sample jitter, and is independent of the model's memory. If the support or sampling changes later, recompute and document this mapping rather than blindly reusing 0.1 or the old 1/60. This physical correlation length is a fixed first-version choice, not an inferred archive property.

### Generation, pairing and rejection

Each independent physical history has its own independent southern and northern chronology package. Within a history, reuse the same standardized observed-center candidate stream across levels, scaled only by sigma. Use separate independent member candidate streams, also shared across levels. Keep southern/northern and observed/member streams independent of each other and of process noise. Preserve candidate identifiers as well as accepted member numbers.

`simulate_observations` does not itself reject an unordered observed center. The TBK preparation layer must explicitly reject and redraw whole invalid centers. `draw_age_ensemble` performs member-order rejection; coverage on the fixed 29-1 kyr grid requires a further explicit check. To retain candidate IDs and avoid repeated identical eigendecompositions, the new preparation module calls the existing `_age_error_factor` once and uses the exact existing Gaussian multiplication in an explicit loop; its proposal values are checked against the native functions. Keep order and coverage rejection counts distinct, with a fixed finite attempt budget and explicit incomplete status. Never sort, clip, recenter or rescale an accepted curve to obtain a desired error, and never extrapolate proxy values.

Each level is conditioned on its own acceptance rules. Do not select valid curves at the highest level and simply scale them down: that would also condition lower levels on high-level validity. Matching accepted member numbers do not guarantee matching candidates if rejection histories differ. The primary paired contrast is always each case against its own physical history's accurate-age reference; an amplitude-only cross-level member contrast requires matching both observed-center and member candidate IDs for both records. A matched member perturbation with a different accepted center does not meet that requirement.

Save nominal sigma, center RMS relative to truth, member RMS relative to truth, member spread about the center, inter-record relative-error RMS, maximum absolute errors, local age spacing/stretch, acceptance counts and common-window coverage. Plot observed-center and member diagnostics separately. Do not label nominal sigma as the achieved chronology RMS.

No temperature values, forcing times or model integration steps are modified by chronology generation. No new process history is generated for an age member. Joint/shared chronology errors and other envelope shapes remain later questions.

## Causal methods and fixed settings

Run the existing pairwise PCMCI+, LKIF, GC and CCM wrappers on T1/T4, with X=South and Y=North. Keep alpha=0.05, native score definitions, numerical precision and the current scientific significance rules. No conditional three-record analysis, new estimator, GPU replacement or four-box graph search is proposed.

Because this branch uses a 20-year observation step, parameter values must be stated in steps and years:

| Method | Proposed main settings |
|---|---|
| PCMCI+ | Native CPU GPDC, optimizer=None, 10000 null samples, alpha/pc_alpha=0.05, maximum lag 5 steps=100 years. Existing directional Holm and native-arrow rules; preserve unresolved contemporaneous evidence. |
| LKIF | Corrected native implementation, maximum lag 5 steps; retain its native lag-index convention, analytical significance, normalization and all native arrays. The nominal model history spans 100 years; native lag labels are not physical propagation-time estimates. |
| Granger | Existing statsmodels fixed order 1, no BIC or member-specific selection. Its single-step history is 20 years, matching the earlier experiments. Preserve SSR-F, native score and VAR stability/residual-whiteness diagnostics. |
| CCM original test | PyEDM 2.5.7, E=3, tau=-1 step (20-year delay; 40-year embedding span), exclusion 10 steps=200 years, L=300, mean of 10 libraries, 100 cut-and-swap surrogates. |
| CCM lag screen | 61 target shifts from -600 to +600 years in 20-year increments: integer steps `-30,-29,...,29,30`. One hundred paired libraries on common anchors, existing 90% sign rule and all failure/unresolved categories. |

The CCM embedding span and exclusion retain their earlier physical durations. The wider lag scan is proposed in advance because the highest chronology level has a relative mismatch scale of hundreds of years. The +/-600-year interval remains below half the 1500-year forcing cycle, limiting exposure to adjacent-cycle replicas but not removing periodic ambiguity. It is not a hard bound on dating error, nor a guarantee that all peaks are interior. Retain boundary outcomes rather than extending a particular member's scan after seeing it. This is an explicitly different branch configuration from the earlier 31-shift +/-300-year scan; old results cannot be imported as these settings.

With 1401 points and these settings, Tp=0 has 1399 embedded states. The lag scan has 1341 common anchors, so L=300 is approximately **22.4%** of that pool. This is the fixed library length previously chosen by the user, now again close to one quarter. Verify actual native anchors in preflight for the wider scan.

The 100-year PCMCI/LKIF history and one-step GC are specified analysis models, not assertions that all TBK response memory lies within those spans. Preserve baseline fit limitations. Do not increase histories or change conditional tests per member to recover desired arrows. A future fixed-history or observation-resolution sensitivity would be reported separately if baseline diagnostics show it matters.

This is a comparison of specified analysis pipelines with explicitly unequal history models, not a ranking of optimally tuned algorithms or a test holding every method's temporal history equal. Changes caused by chronology are paired within each fixed pipeline.

The negative-lag screen can be applied to reciprocal systems: [Ye et al. (2015)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4592974/) explicitly studied bidirectional examples. However, neither their examples nor our existing checks calibrate the project's 90% rule for this periodically forced TBK system. Keep raw Tp=0 evidence and screened evidence separate. Neither a positive/zero/boundary peak nor a nonsignificant reverse test proves that physical feedback is absent. Cut-and-swap significance is relative to that surrogate construction; it does not test the absence of all shared periodic forcing. Do not silently replace the surrogate method to improve results.

## Primary comparisons and metrics

The primary analysis is within history: the same sampled values and estimator settings under accurate and uncertain dates. Both records are uncertain in all nonzero inputs. Use requested case counts as denominators and retain all baseline histories, including unfavorable or unresolved accurate-age outputs.

### Directional evidence states

For a compact mutually exclusive display, classify each method/case in this order:

1. Failed or unavailable computation.
2. Unresolved evidence in either direction, if computation is available.
3. Both directions retained, with no unresolved evidence.
4. South-to-North only, with no unresolved evidence.
5. North-to-South only, with no unresolved evidence.
6. Neither direction retained, with no unresolved evidence.

Keep the original two directional retention flags alongside this state because unresolved evidence can coexist with a retained arrow. Report their marginal fractions separately; an unresolved classification must not erase a successful directional detection.

The main method-comparison curve is the **fraction retaining both directions without unresolved evidence**, supplemented by the two one-way fractions, unresolved fraction and failures. Call this bidirectional evidence retention, not full-network accuracy or a calibrated feedback-detection rate. Do not use F1 as the primary metric: this experiment has neither a complete observable graph nor a suitable set of true absent endpoint directions for such an interpretation.

### Paired changes and ensemble support

- Save each accurate-age state and its transitions to all perturbed states. Distinguish a newly one-way result from a history already one-way with accurate ages. A focused count of both-to-one transitions must state its baseline-bidirectional denominator; it cannot replace the all-history main results.
- For each direction, save delta score = perturbed score minus that history's accurate-age score, using all finite successful scores, not just significant members. Plot PCMCI, LKIF, GC and CCM separately using their own score definitions and panel scales; these scores are dimensionless but not directly interchangeable. Do not divide by near-zero baseline scores or combine native scores into a physical influence scale.
- Retain within-history q_S_to_N and q_N_to_S, joint usable/requested counts and unresolved counts using the existing ensemble conventions. q>=0.9 is a reporting rule; with ten members it describes this finite ensemble and is not a causal probability. Do not treat 30 members from three histories as 30 independent histories.
- For CCM, retain original rho/p, raw directional significance, screened retention, lag rejection, unresolved categories, all 100x61x2 rho values and library peak-sign fractions. The original Tp=0 rho remains the strength statistic; do not substitute the maximum scanned rho.
- Show nominal levels and actual inter-record relative RMS together. Associations between realized error and evidence changes are descriptive at this pilot size, not universal error thresholds.

Feedback is a property of the model equations and its response to perturbations. Statistical evidence from two endpoint records also contains common-forcing and hidden-state effects. An accurate-age baseline is a reference for measuring added chronology sensitivity, not a substitute for physical truth. If no method reports two directions at accurate dates, the experiment still diagnoses a baseline limitation and subsequent evidence changes; do not redesign the physical system or omit those histories to manufacture a both-to-one example.

## Pilot size and interpretation

Propose three independent process-noise histories using recorded seeds, each with independent southern/northern chronology centers and ten members per nonzero level. This is three independent history/chronology blocks, not a crossed variance-decomposition experiment. Add one deterministic trajectory for model illustration and migration checks only; do not include it in stochastic outcome proportions or fit it as an extra main case.

$$
N_{\mathrm{cases}}=3\,[1+4\times10]=123.
$$

The main experiment comprises 492 full-setting method calls, 492 native archives and 984 directional rows. Each nonzero level has 30 requested cases; the accurate-age level has three. Label these unequal denominators clearly. Chronology members are paired within each history and do not inflate the independent-history count. Three blocks are a compact pilot; broader error-rate claims require subsequent independent-history/center replication, not merely more members.

Use 32 outer worker processes, single numerical threads, PyEDM internal parallelism off, native GPDC null caching and validated CCM neighbor reuse. For T=1401 and PCMCI maximum lag 5, the expected native effective null sample size is 1391; verify against Tigramite before generating the new cache. Observe memory use, but do not repeat the completed 24-worker benchmark. Estimate runtime from the full-setting preflight rather than promise an unmeasured duration.

## Notebook and figure design

The inspection notebook should be short and callable from a clean kernel:

1. Show the explicit model version, parameters, noise, seed, integration length and retained window; call the toolbox model.
2. Plot the transport forcing and all four accurately dated temperatures, highlighting the two endpoints. Include a fixed one-cycle zoom. Clearly distinguish deterministic illustration from a noisy realization and mark discarded startup/support margins.
3. Call the separate chronology preparation. Show southern and northern age offsets for the preselected first history/member across all levels, plus actual error summaries across all accepted members. A member ID is illustrative, not an independent climate realization.
4. Plot the same endpoint values on accurate and perturbed age axes at all five levels, with consistent colors, axes, fixed zoom interval and panel labels. Make level, nominal SD and realized relative RMS visible; do not choose examples after inspecting causal outcomes.
5. Display checks for zero-error reconstruction, unchanged sample values, ordered ages and fixed-grid coverage. No causal fits are required to inspect or change the plotted error level.

Four experiment figures are sufficient initially:

| Figure | Purpose |
|---|---|
| 00 Model and records | Explain four-box interactions, shared prescribed forcing, selected endpoints, response shapes and chronology distortion. |
| 01 Evidence versus uncertainty | Four method panels with requested-denominator state proportions; an aligned method-comparison view of bidirectional evidence retention. Accurate-age counts are explicit. |
| 02 Within-method score changes | Both directions' paired score changes versus error level, with descriptive distributions and zero-change guides; no shared strength scale. |
| 03 CCM evidence and timing | Separate raw significance from lag decisions, show two directional lag profiles across levels and boundary/zero/sign-instability diagnostics. |

Keep baseline-to-perturbed transitions in a readable summary table and annotate consequential transitions in Figure 01. Save all individual states and profiles. Preserve research-style informative labels, panel letters and PNG/PDF companions; write English captions under `agent/figure_review/` after verifying the actual figures. Do not call member/library spreads confidence intervals.

## Acceptance checks and staged implementation

1. **Archive and migrate only.** Verify the original notebook's byte-identical archive before moving it from the root. Compare full model arrays from the extracted core and the archived source with identical random inputs: all physical/ghost temperatures, q/target, diffusion/transport/restoring fluxes, noise and global/area means. Include noisy and noise-free configurations and both retained original diffusion parameters. The default 10-kyr calculation must match exactly, not merely resemble the old plots.
2. **Prepare and inspect.** Verify the 39-kyr integration indexing and prescribed events without changing their update order. Check the 30-kyr retained interval, 1401-point analysis grid, exact zero-error input, fixed values, chronology streams, separate center/member conditioning and coverage. Inspect model and chronology figures before causal production.
3. **Full-setting preflight.** Run a selected accurate-age input and a highest-level input through all four methods; compare complete native outputs between direct calls and Windows-spawn execution, and native CCM queries versus cached reuse for the new tau/lag grid. Check resume/configuration guards, failure/unresolved retention, new null metadata and result serialization. Inspect GC stability/whiteness, PCMCI unresolved evidence and CCM boundary/periodicity diagnostics. Success means correct computation and interpretable diagnostics, not detection of both directions.
4. **Freeze this branch's specification.** Record the model variant, retained window, chronology generation, physical/step settings and seeds before the 123-case run. Any consequential revision motivated by an implementation problem must be documented before production; do not revise levels or methods after seeing unfavorable chronology results.
5. **Run and audit.** Save prepared records and complete native outputs under `data/processed/04_check_tbk_chronology/`, with figures under its matching stem. Use eight directional rows as the final per-case commit, preserve failures and resume only matching complete inputs/settings. Audit reconstruction of summaries and preservation of previous 02/03 results.

The user has authorized these implementation and pilot steps. Relative to the initial discussion draft, the observation interval is now 20 years and retained duration is 30 kyr. The model source, startup, uncertainty levels and wider CCM scan remain explicit branch choices. No previous 02/03 results are recalculated or overwritten.
