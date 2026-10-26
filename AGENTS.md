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
| `toolbox/causal.py` | Thin calls to native causal packages, uniform resampling and CCM lag screening; no plotting. |
| `toolbox/evaluation.py` | Ensemble support, complete member states, native-score summaries and independent-realization evaluation. |
| `01_inspect_simulated_records.ipynb` | Inspect generation and simulated observations. |
| `02_check_causal_pipeline.py` | Two accurate-age baselines and 16 ordered type pairs at uncertainty level 2. |
| `02_1_check_coupling_strength.py` | Original-uniform-record native-score ranking under K=1/2/4, with three process seeds. Former 03. |
| `02_2_inspect_ccm_lags.ipynb` | CCM embedding/library/lag-window parameter investigation. Former 04. |
| `02_3_check_uncertainty_escalation.py` | All four levels and 16 ordered pairs; imports completed level 2 from 02. Former 05. |
| `run_experiment.py` | Portable launcher: one numerical thread per worker, console plus saved log. |
| `benchmark_workers.py` | Windows 16/24/32-worker throughput comparison on 96 fixed prepared02 inputs. |
| `verify_optimizations.py` | Portable same-input comparison of uncached/fresh native calls and optimized calls. |
| `explore_cuda.py` | Optional CUDA hardware/matrix-operation probe; **not** a GPU causal estimator. |

Branch numbering expresses a relationship to 02, not a requirement to run every notebook before the main experiment. Production calculations do not import code or data from `tmp/`.

## Scientific question and scope

The main question is how equal or different chronology-uncertainty types assigned to two records affect retained causal arrows and native method scores. The ensemble provides within-record support. One highly consistent ensemble does not establish causal correctness across independently generated climate histories.

Current controls deliberately keep sampling, proxy values, dynamics and estimator settings fixed. There is no main-study systematic age shift, feedback network, five-variable network, tipping experiment, measurement-noise sweep, or sampling-sensitivity sweep. Three-variable common forcing/confounding and a real-data case remain later work.

### Time and data conventions

- Age is in **kyr BP**, with 1950 CE as zero; positive values are older.
- Data rows advance forward in physical time: `time` increases and `age` decreases.
- `calendar_year = 1950 - 1000 * age`. Do not reverse the causal direction when plotting BP ages.
- Plots put older ages on the left and use `Age (kyr BP)`.
- Native process step and analysis step are 20 years (`0.020 kyr`). The retained process spans 30 kyr, 1501 points.
- Analysis uses the fixed 29–1 kyr window, 1401 points. Every accepted chronology member must cover the same grid. No member-specific cropping or extrapolation.
- The accurate-age references are: (1) original uniform data, cropped only; (2) irregular observations with their accurate ages, interpolated onto the same grid. There is no separate downsampled baseline.
- Main irregular sampling uses 27.5-year pre-rounding target gaps at both ends with 10% jitter; native retained gaps are 20/40 years. This mild sampling was chosen to separate chronology effects. Age-dependent thinning is deferred sensitivity work.

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
- **CCM:** PyEDM 2.5.7; E=3,tau=-1, temporal exclusion 10 steps, library fraction 0.8. Original Tp=0 skill averages 10 libraries; the directional swap-null test uses 100 surrogates. X→Y means reconstructing X from Y's embedding.
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

### GPDC null cache

`prepare_gpdc_null(filename, sample_sizes=[T], sig_samples=..., seed=...)` calls the native generator once in the parent. Pass its returned path to `run_pcmciplus(..., null_dist_filename=...)`.

Current effective T is 1401−2×5=1391. Use one T and seed per file to preserve the native random sequence. Cache metadata checks T, seed, sample count, relevant packages and platform. Do not reuse an unrelated null distribution. PCMCI deep-copies its conditional test, so simply retaining an external GPDC object is not sufficient.

Every member still receives its own conditional GP regressions. Do not cache residuals or conditional-test results across different series.

### CCM object reuse

Each lag uses one native PyEDM object, then calls its `Project(sample=1, seed=original_library_seed, parallel=False)` for the requested libraries. This avoids repeated embedding preparation but retains all neighbor searches and 100×31×2 rho values.

Do not replace this with 31 calls at `sample=100`: that does not preserve the required individual paired library curves and random streams. Do not change embedding, nulls, library count or lag range as a performance shortcut.

### Multiprocessing and scale

- 02 and 02_3 default to 24 outer workers as a **candidate** for the 14900KF, not a proven optimum.
- Use one BLAS/OpenMP thread per worker. The formal launcher sets these variables before importing numerical packages.
- Keep an outer process pool; do not substitute threads. The current PyEDM sequential path uses mutable module-level worker data.
- Keep PyEDM internal parallelism off. Repeated nested pools can waste resources and memory.
- Keep worker functions at module scope and preserve `if __name__ == "__main__"` for Windows spawn.
- Benchmark16/24/32 on the target PC with 96 identical inputs. Observe memory/commit and paging in Task Manager. Select the fastest configuration with acceptable memory pressure; do not pick 32 from the logical CPU count alone.
- Main case counts are derived from configuration:02 is 16×3×N;02_3 is 4×16×3×N, including imported level 2. N = 10 gives 480 and 1920 pairs respectively.

## Windows environment and launch procedure

Use 64-bit Python 3.11 and Git for Windows. Run commands from the project root. Do not copy `.venv` from Mac, and do not install unconstrained latest scientific packages.

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

`requirements.txt` fixes the scientific versions. NumPy 1.26.4 is intentional for LKIF compatibility. LKIF is installed from Git commit `95c088ef8cc5e5dab2c6279158e98338ee51124a`; do not replace it with the known problematic PyPI wheel. A fresh Windows installation still needs to be tested on the actual machine.

For notebooks, optionally install and register the project kernel:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements-notebook.txt
.\.venv\Scripts\python.exe -m ipykernel install --user --name chronology-causal --display-name "Python 3 (project .venv)"
```

Recommended first Windows workflow:

1. Run `.\.venv\Scripts\python.exe verify_optimizations.py`. This is a few-case numerical check, not a full experiment.
2. In 02 set `PREPARE_ONLY=True`; leave `RUN_BASELINES=True` and `RESUME=True`. Run the launcher below. It computes the small baseline set and prepares all chronology inputs, then stops before the ensemble fits.
3. Run `.\.venv\Scripts\python.exe benchmark_workers.py`. Choose the worker count from measured throughput and memory behavior. The96-case benchmark is real causal computation and may take several minutes per worker setting.
4. Set `PREPARE_ONLY=False` and the selected `N_WORKERS` in 02; launch it again. It reuses the matching prepared baseline and completes level 2.
5. Keep the same chronology counts in 02_3. Change the launcher's `EXPERIMENT` to `02_3_check_uncertainty_escalation.py`, then launch to calculate the other levels and import the completed local 02 level 2.

```powershell
.\.venv\Scripts\python.exe run_experiment.py
```

The launcher uses its own interpreter, creates a timestamped log under `data/logs/`, and streams output to the terminal. Keep that terminal open and run only one instance per experiment output folder. Ctrl+C requests termination of the child process group; completed member files remain. Check `run_status.json`, not an old log selection, for the latest saved progress. On Windows the launcher temporarily requests that automatic system sleep be prevented, then releases that request on exit. Display sleep remains allowed; no persistent power setting is changed. If this API fails, the launcher warns and continues.

The strength script and lag notebook are optional branches; they need not be rerun before every 02 experiment.

## Resume and fresh experiments

- `RESUME=True` resumes complete per-member outputs only when the scientific configuration, packages, implementation revision, numerical environment and regenerated inputs match.
- A completed member can include an explicitly recorded failed fit. Keep it as evidence; do not automatically rerun only unfavorable or failed scientific results.
- A member without its complete committed result table is incomplete and may be recalculated. Partial native files are not valid complete votes.
- Worker count and run-control flags can change without redefining the scientific experiment. Changing chronology count, dynamics, chronology settings, methods, package versions or platform requires a fresh run.
- To increase ensemble size, preserve/rename the old output folders, set matching counts in 02 and 02_3, and run 02 first. Do not edit saved metadata to force an old result through a guard.
- For 02_3 figures only, set `RUN_CAUSAL=False` with matching saved settings. This does not recompute estimates. Pending levels must remain explicitly partial.
- This is a simple same-run resume mechanism, not a general cache search across historical experiments.

## Files to synchronize

Git should contain scripts, source notebooks without embedded generated output, toolbox, requirements, this guide and the small directory README files. This repository does not require a transferred Mac dataset to start the Windows pipeline.

Keep `.venv/`,`.venv-gpu/`, `tmp/`, generated data, figures, caches, logs and historical archives out of Git. Ignoring a file is not deleting it. If a cloud folder is used instead of Git, configure its exclusions separately; `.gitignore` does not control a cloud-sync application.

The original repository had no tracked files or commits at migration time. Check `git status` before creating the first curated commit. Do not upload credentials, local environment directories or the legacy experiment archive by accident.

All temporary agent work belongs in `tmp/` locally. The Windows agent must be able to understand/run the project from this guide and checked-in source; no personal skill installation or `/Users/...` path is a runtime requirement. Readable comments and notebook prose should be English; discussion with the user can use their chosen language.

Use plain functions, explicit parameters and readable intermediate calculations. Do not introduce a class/configuration framework merely to remove a few repeated lines. Give every configurable scientific parameter a comment explaining its meaning and units. Notebook sections should have a title and one or two plain English sentences. Keep plotting in root analysis scripts/notebooks, and agent scratch calculations under `tmp/`. Future real-climate inputs will need an explicitly chosen storage/ignore policy; the current broad data exclusion is for this synthetic migration.

## Optional CUDA exploration

Production GPDC/sklearn GPR, PyEDM KDTree, LKIF and statsmodels currently run on CPU. Installing CuPy does not switch these packages to CUDA.

If a compatible NVIDIA driver and CUDA Toolkit 12.x are installed, use a separate environment. The pip wheel requires the CUDA dependencies described in the [CuPy 13.6 installation guide](https://docs.cupy.dev/en/v13.6.0/install.html); installing the Python wheel alone is not a complete CUDA setup.

```powershell
py -3.11 -m venv .venv-gpu
.\.venv-gpu\Scripts\python.exe -m pip install -r requirements-gpu.txt
.\.venv-gpu\Scripts\python.exe explore_cuda.py
```

The optional requirements pin CuPy 13.6 to remain compatible with the existing NumPy 1.x line. The probe checks CUDA availability, compares a float64 RBF-covariance/Cholesky solve and records timings including transfers. **Its speed ratio is not an end-to-end PCMCI/CCM speedup.** It does not implement a new causal estimator.

Keep this experiment independent of CPU production. A real GPU backend would need a separate validated implementation from an existing package, full numerical/decision comparisons and a measured batch-throughput benefit. Do not replace methods or lower precision just to use the GPU.

## Validation already performed and limits

The Mac implementation check used two fixed 1401-row inputs (linear and independent). Cached/uncached PCMCI arrays and graphs were exactly equal; fresh/reused-object CCM matched the full 100×31×2 rho array, nulls, seeds and arrow results exactly. Cache mismatches were rejected.

Observed single-process times were approximately 18.7→4.9 seconds and 17.7→3.8 seconds for PCMCI+, with a one-time 16.8-second null preparation. CCM changed 49.2→47.0 seconds and 51.2→49.9 seconds. These checks suggest about 1.3× combined improvement, **not** a 4× improvement of the whole workflow.

Fresh setup, resume, failure retention and metadata rejection were checked using small mocked workflows. A separate real-estimator integration check used two spawned processes and reduced test-only settings: both scripts saved eight successful directional rows. Launcher logging, exit codes and temporary sleep-control lifecycle were checked locally or with mocks. Source and notebook-cell syntax checks passed; notebook execution was not repeated. Windows throughput, 32-worker viability, fresh package/notebook installation, Windows process control and CUDA execution still require acceptance on the target PC.

Retain the original scientific outputs if a future optimization fails the numerical comparison. Investigate discrepancies rather than changing significance thresholds or dropping difficult cases.
