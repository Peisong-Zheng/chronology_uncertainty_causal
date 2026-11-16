# Project support tools

This directory contains maintained, versioned support code and pinned dependencies. Temporary agent work belongs in `tmp/`. Scientific scripts and notebooks remain at the project root, and the existing `toolbox/`, `data/` and `figures/` locations are unchanged.

The [research direction and TBK discussion record](research_direction_and_tbk_notes.md) tracks the 10 November mechanism-focused framing, the decision to retain four TBK boxes and use their southern/northern endpoints, the verified process-noise implementation, and proposed chronology controls. Further three-variable mistiming sensitivity is deferred; the TBK branch has been authorized for implementation and a compact pilot.

The [04 TBK chronology plan](04_tbk_chronology_plan.md) replaces one-sided errors with graded uncertainty in both endpoint records. It proposes archival and minimal model extraction, separate shared chronology preparation, an inspection notebook, a 123-input pilot and explicit evidence-state metrics. The accepted duration revision uses 39 kyr annual integration, 9 kyr startup, 30 kyr retained support and the common 29-1 kyr/20-year grid. Exact migration, native-output and Windows-spawn preflight have passed. The 123-input pilot is complete; full native-output and preservation audits passed.

- `redraw_figures.py`: redraw all 30 existing figures from saved results through three plot-only component scripts, without running estimators.
- `run_experiment.py`: launch the root scientific scripts sequentially with one numerical thread per worker and a separate log for each script.
- `verify_optimizations.py`: compare complete native numerical outputs before and after the existing performance optimizations.
- `verify_ccm_neighbor_reuse.py`: compare the migrated full CCM pipeline with the local timestamped source freeze, including failed scans and two spawned workers.
- `benchmark_workers.py`: compare worker counts on the same prepared experiment-02 inputs.
- `explore_cuda.py`: run an optional CUDA matrix-operation probe in a separate environment; this does not move the causal estimators to the GPU.
- `ccm_gpu_assessment.md`: CCM/GPU feasibility review, complete-output native query-cache pilot, and computation planning for a larger 02_3 study.
- [`figure_review/`](figure_review/README.md): GRL/AGU style sources, figure plans and English captions for 45 current exports and 10 additional notebook figure definitions. Captions describe the restored research figures and the new 02_8/03/03_1/03_2 branches; the 10 missing notebook exports remain unexecuted.
- `requirements.txt`: pinned scientific dependencies, including the corrected LKIF commit.
- `requirements-notebook.txt`: optional notebook tools plus the scientific dependencies.
- `requirements-gpu.txt`: optional CUDA probe dependencies plus the scientific dependencies.
- `requirements-gpdc.txt`: isolated PyTorch CUDA/GPyTorch and notebook dependencies for the small root `02_5_compare_gpdc_backends.ipynb` experiment.

The notebook and GPU requirements include their sibling `requirements.txt` using `-r requirements.txt`; keep the three requirements files together.

Run these Windows commands from the project root:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r agent\requirements.txt
.\.venv\Scripts\python.exe agent\verify_optimizations.py
```

The current launcher selects only `04_1_check_tbk_excitation.py`, the accurate-age TBK excitation diagnostic described in its [plan](04_1_tbk_excitation_plan.md). The original 04 chronology pilot is complete. The earlier `03_2_check_driver_error_models.py` 732-case Z-only C/L/E/J comparison is complete, with 192 complete imports and 540 new conditional inputs using 32 workers. See the [03_2 plan](03_2_driver_error_models_plan.md) and [results](03_2_driver_error_models_results.md). The earlier 03_1 cross is complete; its [plan](03_1_process_chronology_cross_plan.md) and [results](03_1_process_chronology_cross_results.md) remain available. The original 03 run is also complete; its [plan](03_common_driver_plan.md) and [results](03_common_driver_results.md) remain unchanged. The [02_8 memory branch](02_8_memory_sensitivity_plan.md) and the main 02/02_3 chronology pilot are complete. The separate sampling branch is described in its [experiment plan](02_4_sampling_sensitivity_plan.md). For a manual launch from the project root, after confirming that no other instance is active:

```powershell
.\.venv\Scripts\python.exe agent\run_experiment.py
```

For a new chronology-study setup on different hardware, set `EXPERIMENT="02_check_causal_pipeline.py"`, `FOLLOWUP_EXPERIMENTS=[]` in the launcher and `PREPARE_ONLY=True` in 02. Prepare the baseline and chronology inputs, then benchmark:

```powershell
.\.venv\Scripts\python.exe agent\run_experiment.py
.\.venv\Scripts\python.exe agent\benchmark_workers.py
```

For the chronology production sequence, set `PREPARE_ONLY=False`, select the measured worker count in both 02 and 02_3, and set `FOLLOWUP_EXPERIMENTS=["02_3_check_uncertainty_escalation.py"]` in the launcher. Then run:

```powershell
.\.venv\Scripts\python.exe agent\run_experiment.py
```

Keep `EXPERIMENT="02_check_causal_pipeline.py"` for this sequence. The launcher stops on a nonzero exit. Read the project-root [AGENTS.md](../AGENTS.md) for scientific constraints, optional environments and resume rules.

For the GPDC notebook, select `.venv-gpdc/Scripts/python.exe` in the IDE. This environment was installed and checked separately from production. To recreate it with Python 3.11:

```powershell
py -3.11 -m venv .venv-gpdc
.\.venv-gpdc\Scripts\python.exe -m pip install -r agent\requirements-gpdc.txt
```

Run the notebook from the project root. It requires the existing 02 original-uniform CSV files and the matching native GPDC null cache. It evaluates six fixed conditional tests, preserves full regression arrays, and writes tables under `data/processed/02_5_compare_gpdc_backends/` and plots under the corresponding `figures/` directory. A verified executed copy is stored beside those tables; the root source notebook keeps generated outputs empty. The aligned GPyTorch regression is local to the notebook and is not used by production PCMCI+.

The root `02_6_check_ccm_library_length.ipynb` is a separate small CPU diagnostic. It reuses the three original-uniform CSV files, compares library lengths 350/700/1119 using 100 paired draws and 31 target shifts, and requires the production scientific packages plus notebook tools (no Torch dependency). Arrays, tables and an executed notebook are under `data/processed/02_6_check_ccm_library_length/`; the two plots are under the matching `figures/` directory. The 2026-11-07 run completed in about 188 seconds with all 55,800 rho values finite. Production CCM settings remain unchanged.

The root `02_7_check_ccm_neighbor_reuse.ipynb` defines `scan_ccm_reusing_neighbors` locally and compares it with the native lag scan at L=300. It keeps at most two exact KDTree query results per draw and restores the temporary native-module patch on exit. It requires the same CPU scientific and notebook packages, the three 02 original-uniform inputs, and one saved nonlinear level-4 C-J member from 02_3 for a small check. Run sequentially; the prototype is not thread-safe. Results and an executed copy are under `data/processed/02_7_check_ccm_neighbor_reuse/`, with figures under the matching `figures/` directory. Its scientific functions and settings remain unchanged after migration; the 2026-11-08 visual revision updates only style imports and the final plotting cell.

## Source freeze and optimization coverage

The pre-migration source archive is `data/archives/code_freeze_20261107T200542517660Z/source.zip`, with a verified per-file manifest and archive SHA-256. It contains the current code, including uncommitted notebooks; generated results and environments remain outside the archive. No old experiment outputs have been removed or rewritten.

| Entry point | Coverage and current status |
|---|---|
| `toolbox/causal.py`: `run_ccm` / `filter_ccm_by_lag` | Now uses bounded exact native neighbor-query reuse. Original Tp=0 and swap-null calculations retain their native implementation. |
| `02_check_causal_pipeline.py` | Both accurate-age baselines and all level-2 members call the shared CCM entry point. |
| `02_1_check_coupling_strength.py` | All gains and process seeds use the shared entry point. Preserve old outputs before launching its existing `RECOMPUTE=True` path. |
| `02_3_check_uncertainty_escalation.py` | New levels 1/3/4 benefit; level 2 must come from a newly completed, configuration-matching 02. |
| `02_4_check_sampling_sensitivity.py` | All 492 input pairs share the same estimator entry points. |
| `02_2_inspect_ccm_lags.ipynb` | Calls the toolbox and benefits indirectly; its scientific pilot settings are unchanged. |
| `02_5`, `02_6`, `02_7` | Remain standalone comparison notebooks; 02_7 retains its independent scientific reference; the later style-only cell edits are explicitly allowed by the verifier. |
| `agent/benchmark_workers.py`, `agent/verify_optimizations.py` | Use the shared optimized implementation and explicit library size for new inputs. No new 24-worker benchmark is needed. |
| `toolbox/causal.py`: `run_pcmciplus` | Future aligned GPyTorch regression entry point; production still uses CPU GPDC. |

The four production scripts now use `CCM_LIBRARY_SIZE=300` for observed scores, nulls and lag scans. Version guards were advanced consistently; old L=1119 metadata is not eligible for resume/import. Existing run flags, output paths, 10 chronology members and 32 CPU workers remain unchanged. Before the next production run, preserve the previous data/figure directories and prepare a fresh 02, then 02_3; do not edit old metadata to bypass the guards.

Migration acceptance passed in `data/benchmarks/ccm_migration_checks/20261107T201205780934Z/report.json`: five complete frozen-versus-migrated CCM comparisons, both L=300 and L=1119, two spawned-worker results, failure retention, repeated states and interruption cleanup all matched or passed. Full numerical results matched exactly; complete CCM timings improved by about 1.85--1.96x in these single-process comparisons. The 02_7 source hash was unchanged at migration acceptance; subsequent plotting-only edits retain the scientific-cell comparison against that freeze. Reproduce this local acceptance check with `.\.venv\Scripts\python.exe agent\verify_ccm_neighbor_reuse.py` while retaining its source archive and fixed input files. The portable `verify_optimizations.py` remains available for a newly prepared checkout.

## GPDC and the next experiment size

The promising GPU candidate is the **aligned float64, fixed-kernel GPyTorch regression** in 02_5. Stock GPDCtorch trains different kernel/noise/mean parameters and uses float32, so it is not an equivalent replacement. A production adapter belongs behind `run_pcmciplus`, keeping the native standardization sequence, distance correlation and cached null calculation. The current notebook checks only six fixed nonempty condition sets; it does not implement the native no-conditioning bypass or the complete `standardize`/`return_means`/`return_likelihood` interface. PCMCI deep-copies its conditional test, which also needs explicit handling.

Before enabling that backend, compare complete PCMCI+ searches on fixed accurate-age and uncertain inputs: visited conditions, residuals, scores/p values, full p/val matrices, native graph, Holm arrays and unresolved/conflicting states. GPU backend/device/dtype and package/runtime details must become saved metadata. The isolated `.venv-gpdc` is not a complete production environment; its requirements omit other causal packages. Test a controlled GPU worker arrangement and memory use rather than starting 32 independent CUDA clients on the single GPU.

Keep **N=10** for the first fresh configuration run after implementation acceptance. This means 480 level-2 input pairs in 02 and 1920 pairs across all four levels in 02_3, including reused level 2, rather than ten total fits. It is large enough to check storage, resume, level-2 reuse and 32-worker throughput/memory. GPDC should remain on the validated CPU path until its separate full-search/GPU checks pass. After the pilot, preselect a larger member count such as 50 or 100; increasing members improves within-record Monte Carlo precision but does not replace independently generated process histories. No production pilot or larger ensemble is started by the migration acceptance script.

The user subsequently authorized the N=10 pilot, started on 2026-11-07 at 20:24:52 UTC. Previous 02/02_3 data and figures are retained with verified file hashes under `data/archives/pre_ccm300_pilot_20261107T202152109Z/`. The original invocation in `data/logs/pilot_ccm300_20261107T202452271853Z/` stopped after a Windows status-file replacement error, with 02 complete and 02_3 at 1440/1920 pairs. A tested, bounded retry for atomic replacements leaves scientific calculations/settings unchanged. Recovery started at 23:13:49 UTC under `data/logs/pilot_ccm300_resume_20261107T231349346354Z/`: check its `pilot_status.json`, `launcher_console.log` and `resources.csv`. Recovery completed all 480 remaining level-4 pairs. The final totals are 480/480 for 02 and 1920/1920 for 02_3; the complete artifact audit and preservation check passed. The older failed invocation is retained as history, not current progress.

For a safe figure-only refresh, run `.\.venv\Scripts\python.exe agent\redraw_figures.py`. The 30 figures retain their original layouts, including all dense panels, and are exported as PNG/PDF pairs. Existing captions are under `figure_review/`; old figures/source remain in `data/archives/pre_figure_style_20261108T121553339241Z/`. The broad visual restyling was rolled back at the user's request: original annotations and appearance are retained, together with panel labels, corrected scientific labels and captions. The 01 and 02_2 notebook plotting cells retain their original content with added panel labels; they are not executed by this launcher because their saved artifacts are absent.

The memory branch preflight is `.\.venv\Scripts\python.exe agent\verify_memory_sensitivity.py`. It checks fixed-clock and innovation pairing, exact reproduction of 36 retained reference inputs, failure/resume semantics and two full-setting cases against direct native calls and Windows-spawned workers. The branch retains 432 complete per-case records with full native evidence; its figures use the original research styling and panel labels. Run status is under `data/processed/02_8_check_memory_sensitivity/run_status.json`.

The memory branch completed on 2026-11-09: 432/432 pairs, 1728 native archives and 3456 successful directional rows in 1458.999 seconds (24.317 minutes), using 32 workers. Minimum sampled available memory was 8.90 GiB. [Results and interpretation](02_8_memory_sensitivity_results.md) compare all four methods and distinguish accurate-age limitations from added chronology effects. Four PNG/PDF figures and their [English captions](figure_review/05_memory_sensitivity.md) are available. Run logs, source snapshot, artifact audit, preservation check and figure review are under `data/logs/memory_sensitivity_20261109T114908564938Z/`. The existing main-study scientific files and outputs were preserved.

The later 2026-11-09 memory-figure revisions change only presentation: graph recovery (01), CCM decision states (02), lag profiles and peak-sign diagnostics (03), paired score changes (04), with descriptive displays replacing seed-specific marker shapes. The latest version adds direct baseline-to-perturbed comparisons, explicit lag-rejected and unresolved states, saved per-input peak-sign fractions, and an expanded independent-system score scale. All three histories, numerical curves and outputs are retained. The source calculation guard and all 3616 existing branch data files are unchanged; see `data/benchmarks/memory_figure_checks/20261109T180042688693Z/report.json`. The preceding figures and source are archived in `data/archives/memory_figures_before_readability_20261109T180042688693Z/`.

The new [03 script](../03_check_common_driver.py) compares **Common driver only** with **Common driver + direct X to Y** at fixed Z/X/Y memories of 500/200/200 years. CPU GPDC, LKIF, equation-F GC and the separate pairwise CCM comparator retain full native evidence. Run `.\.venv\Scripts\python.exe agent\verify_common_driver.py` for the paired-input, native and Windows-spawn preflight. The six figures and [captions](figure_review/06_common_driver.md) focus on whether including Z removes misleading X/Y evidence and how wrong dates change that conclusion. The 10 November additions show conditioning gain and method comparisons using the existing exact XY recovery metric, with the original views retained. The run completed 192/192 inputs in 2046.07 seconds; production logs and audits are in `data/logs/common_driver_20261109T211815766544Z/`, and the figure revision checks are in `data/benchmarks/common_driver_figure_checks/20261110T105738358740Z/`.

The [03_1 follow-up](../03_1_check_process_chronology_cross.py) crosses the three saved realizations with all three saved Z age schemes, keeping X/Y accurate and assessing only their relation. It completed 192/192 cases in 799.423 seconds with 32 workers, including 360 new conditional calls and 984 exact native-archive copies. The [results](03_1_process_chronology_cross_results.md) and [two reviewed figure captions](figure_review/07_process_chronology_cross.md) explain the row/column comparisons and baseline exceptions. Complete native-output, summary, figure and preservation audits passed under `data/logs/process_chronology_cross_20261110T124145323892Z/`; all 26,464 pre-existing scientific files remain unchanged. The maintained preflight is `agent/verify_process_chronology_cross.py`.

The [03_2 comparison](../03_2_check_driver_error_models.py) extends that crossing to Z-only C/L/E/J errors and completed 732/732 cases. Following an accidental closure, it retained all 297 complete cases and finished the remaining 435 with 32 workers; 2535.030 seconds (42.25 minutes) is the recovery-stage runtime only. All native, summary, figure and preservation audits passed under `data/logs/driver_error_models_resume_20261110T135403003994Z/`. The 28595 pre-existing scientific files and 3267 files from completed pre-recovery cases remain unchanged. The [results report](03_2_driver_error_models_results.md) explains the unequal error scales and accurate-age exceptions; [reviewed captions](figure_review/08_driver_error_models.md) accompany three PNG/PDF figures. The maintained checker is `agent/verify_driver_error_models.py`.

The TBK model migration checker is `agent/verify_tbk_migration.py`; its exact same-input acceptance report is `data/benchmarks/tbk_migration_checks/20261110T215254006646Z/report.json`. The original notebook has been moved to the verified [timestamped archive](../data/archives/tbk_original_20261110T214805311811Z/One_D_model_TBK.ipynb). Use [04_0_inspect_tbk_records.ipynb](../04_0_inspect_tbk_records.ipynb) to inspect the model and chronology changes; [04_check_tbk_chronology.py](../04_check_tbk_chronology.py) is the completed chronology pilot. The launcher selects the subsequent 04_1 diagnostic only.

The TBK pilot completed 123/123 inputs in 1460.470 seconds (24.34 minutes) with 32 workers and at least 6.892 GiB sampled available physical memory. All 492 native archives, 984 directional rows and seven summaries passed audit; all 36,677 pre-existing scientific files remain unchanged. [Results and interpretation](04_tbk_chronology_results.md) distinguish accurate-age limitations from chronology-induced changes; [English captions](figure_review/09_tbk_chronology.md) explain the four production figures and three inspection figures. Logs, source capture and final audits are in `data/logs/tbk_chronology_20261111T171153343045Z/`. Do not refit to remove unresolved or unfavorable evidence.

The subsequent [04_1 excitation diagnostic](../04_1_check_tbk_excitation.py) completed twelve accurate-age inputs and three baseline ParCorr controls in 71.771 seconds, with 39 new native calls and twelve exact source copies. [Results](04_1_tbk_excitation_results.md) separate the forced cycle, stochastic response and fixed physical memory; [five reviewed figure captions](figure_review/10_tbk_excitation.md) cover directional, lag and residual evidence. All 51 archives, 102 directional rows and thirteen tables passed audit, and all 1,026 inventoried prior TBK files remain unchanged. Use `verify_tbk_excitation.py` for its isolated preflight. Run capture and final audit files are in `data/logs/tbk_excitation_20261113T121730984237Z/`. Further ParCorr, full-state, sampling and chronology comparisons are proposals, not additional completed experiments.
