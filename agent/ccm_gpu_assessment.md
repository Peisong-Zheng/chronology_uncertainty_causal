# CCM acceleration assessment for experiment 02_3

Review date: 2026-11-06. This is a feasibility assessment, not a new production
backend. The current scientific scripts, estimators, package environment and
completed experiment outputs remain unchanged.

## Recommendation

CCM is a worthwhile acceleration target, but GPU use is not an automatic
consequence of installing CuPy. First remove repeated native neighbor searches
across the existing paired lag curves. Then evaluate a compatible GPU backend
against that improved CPU reference, using complete numerical results and total
batch throughput. Do not compare a GPU prototype only against one CPU worker.

The GPU route remains conditional on matching the project's CCM definition.
Current public GPU implementations are useful references, but are not verified
replacements for the installed PyEDM 2.5.7 pipeline.

## What the existing files and machine establish

- `agent/explore_cuda.py` tests a float64 RBF covariance construction and
  Cholesky solve at T=1391. It does not execute CCM or a complete GPDC test.
  Its speed ratio cannot establish an end-to-end causal-method speedup.
- `agent/requirements-gpu.txt` pins CuPy 13.6.0 in an optional separate
  environment. Installing this file does not redirect SciPy or PyEDM to CUDA.
- `agent/verify_optimizations.py` already provides a useful full-output
  comparison pattern. The earlier optimization reused native CCM objects,
  without reusing neighbor queries across lags.
- `nvidia-smi` detects an RTX 4080 with 16,376 MiB memory and driver 610.88.
  During inspection about 1,146 MiB was used by the desktop and applications.
- The project has `.venv`, but no `.venv-gpu`; CuPy and Torch are absent from
  `.venv`. No `nvcc`, `CUDA_PATH`, or default CUDA Toolkit directory was found.
  These checks do not prove that no CUDA library exists elsewhere, but do not
  establish a usable project GPU environment. The CUDA version displayed by
  the driver is not evidence of an installed Toolkit.
- Installed numerical versions are PyEDM 2.5.7, NumPy 1.26.4, SciPy 1.15.3,
  Tigramite 5.2.10.1 and scikit-learn 1.6.1.

CuPy's [versioned installation documentation](https://docs.cupy.dev/en/v13.6.0/install.html)
describes the CUDA requirements and isolated installation alternatives. No GPU
packages or system Toolkit were installed during this review, and no GPU
performance result is claimed.

## Measured production workload

The saved 02_3 results contain 1,440 newly computed input pairs at levels 1, 3
and 4. Level 2 is excluded here because its timings came from the earlier 02
run. Deduplicate the two directional rows before summing each method's time.

| Method | Mean saved seconds per pair | Share of summed method time |
|---|---:|---:|
| CCM | 68.972 | 54.68% |
| PCMCI+ | 56.685 | 44.94% |
| LKIF | 0.423 | 0.34% |
| GC | 0.048 | 0.04% |

These are worker wall-clock times under 32-worker load, including method-output
serialization. They are not isolated single-core timings or run elapsed time.
The source is `data/processed/02_3_check_uncertainty_escalation/member_results.csv`;
the calculation is recorded in
`tmp/gpu_ccm_review_20261106/saved_timing_review.json`.

If resource contention and scheduling stayed unchanged, the simple Amdahl
estimate for accelerating only CCM would be
`1 / (0.45316 + 0.54684 / CCM_speedup)`. Even eliminating all CCM time gives
about 2.21 times overall speedup in that model. This is a planning illustration,
not a physical bound on a different CPU/GPU schedule or a GPU prediction.
PCMCI+ would become the next major bottleneck.

## Why the paired lag scan repeats work

`toolbox/causal.py` uses the same common anchors, library size and per-library
seeds at every lag. With the current finite 1401-row inputs, E=3 and shifts
-15 through +15, the 1371 common anchors have the same embedding coordinates
at every target lag. Only the target values being reconstructed change.

The installed `pyEDM/CCM.py` builds a SciPy KDTree from each sampled library,
queries 25 candidate neighbors, excludes self and temporal neighbors, then
selects E+1=4 neighbors and calculates simplex weights and Pearson rho.
The current wrapper repeats this search for every lag:

- lag screen: 100 libraries x 31 lags x 2 directions = 6,200 searches;
- Tp=0 and swap nulls: 2 wrapper directions x 101 observed/null calls x
  10 libraries x 2 native directions = 4,040 searches;
- total: 10,240 native directional library searches per input pair.

The lag screen needs only 100 x 2 distinct neighbor geometries when these
invariants hold. All 6,200 rho values must still be computed and retained.
This does not justify replacing the paired seeds with one `sample=100` call.

Caching the exact native KDTree query output is a smaller change than porting
the estimator. Native random draws, tie behavior, exclusion filtering,
weighting, predictions and correlations can all remain in PyEDM. Any cache
must be local to matching input arrays and query settings; different records,
chronology members or shifted embedding surrogates cannot share geometry
merely because they have the same dimensions.

The swap-null stage also computes both native maps although each wrapper
call selects one. Removing the unused map is a possible later optimization,
but would need an explicit supported native interface and a review of saved
raw-output obligations. It is not included in the query-cache experiment.

## Isolated native-cache feasibility experiment

Two real, unchanged level-4 input files were tested at the complete production
CCM settings: linear C-J member 0 and independent E-L member 0. The scratch
experiment temporarily wraps only the installed PyEDM module's KDTree calls
and restores them after each run. Cache keys include the complete library and
query arrays (shape, dtype and SHA256 of their bytes) and query settings.
All misses execute the original SciPy query, and all later estimator operations
remain native PyEDM. This is an exploratory wrapper, not production code.

Both complete result dictionaries compare exactly in numerical values:
five DataFrames, eight arrays and 130 scalar leaves per input. This includes
the 100 x 31 x 2 lag rho array, 100 x 2 null array, original native tables,
cuts, seeds, scores, p values, flags and diagnostics. No tolerance was used.
All lag values are finite and both directional statuses are `ok`.

The first baseline used cProfile to identify the bottleneck. Native KDTree
queries took about 25.28 of its 33.37 seconds (about 76%); tree construction
took about 1.30 seconds. Its profile-instrumented time should not be used as
a fair speed ratio against an unprofiled call.

The second input used no profiler in either execution:

| Independent level-4 E-L input | Current native path | Native query reuse | Ratio |
|---|---:|---:|---:|
| Original Tp=0 and swap nulls | 13.263 s | 13.253 s | About 1.00 |
| Complete paired lag scan | 20.782 s | 4.359 s | 4.77 |
| Complete CCM | 34.045 s | 17.612 s | 1.93 |

For both inputs the lag stage executes 200 real queries and returns 6,000
cached query results, while the original/null stage still executes all 4,040
queries. This is a single-pass, single-process CPU feasibility result with
lightweight instrumentation. It is neither a GPU speedup nor measured
32-worker production throughput, and fixed execution order can affect timing.

The current scratch cache holds 109,680,000 bytes of query arrays and
5,371,200 bytes of tree data, about 109.72 MiB per worker before tree nodes,
Python objects and allocator overhead. At 32 workers that is at least about
3.43 GiB additional payload. A production design should bound cache lifetime
and measure memory pressure. Changing loop order to retain only the current
library's geometry is a possible further design, requiring its own validation.

After this optimization, the original/null stage is about 75% of complete CCM
time in the unprofiled pilot. A later GPU investigation should therefore also
address batched swap-null neighbor searches; speeding only the remaining lag
predictions would miss most of the residual cost. Swapping the embedding
changes both boundary embeddings and temporal positions: do not assume its
geometry is an unrestricted rotation of the original record.

Evidence and reproducible scratch code are in
`tmp/gpu_ccm_review_20261106/check_native_geometry_cache.py`,
`native_geometry_cache_report.json`, the profile text and complete local
baseline/candidate result archives. The first comparison validator incorrectly
used NaN comparison on a Unicode directions array. Its failed audit was
retained; the validator was corrected and the already-saved results compared
again, without recomputing any estimate.

## GPU options and compatibility

The official [kEDM repository](https://github.com/keichi/kEDM) implements CCM
with CPU/GPU backends. Its current [CCM API](https://kedm.readthedocs.io/en/latest/_generate/kedm.ccm.html)
documents float32 inputs and an approximation control. The published signature
does not expose equivalents of PyEDM's `exclusionRadius` or `validLib` mask.
`accuracy=1` avoids its approximate-search option, but does not establish
float64 arithmetic, identical random libraries, temporal exclusion, common
anchors or native tie behavior. Its
[installation guide](https://kedm.readthedocs.io/en/latest/setup.html) does not
document a prebuilt Windows wheel route. It is therefore a candidate for a
separate compatibility investigation, not a drop-in production replacement.

For this workload a plausible GPU design would batch exact neighbor geometry
and prediction work across independent input pairs, while retaining CPU-generated
library indices and swap cuts. For one direction a dense 1371 x 1371 float64
distance matrix is about 14.34 MiB; two directions are about 28.68 MiB before
selection buffers and other workspaces. This is a size estimate, not an
implementation or a timing result. Avoid materializing copies for every lag
and library at once.

A candidate backend must preserve:

- float64 inputs and calculations; no silent float32 or approximate k-NN;
- exact sampled library indices, sampling without replacement and existing
  independent random streams for cuts and libraries;
- self exclusion, the 10-step temporal exclusion and zero-distance handling;
- common anchors, target-shift sign and all 31 shifts in both directions;
- 100 paired lag libraries, 100 swap surrogates and 10 Tp=0 libraries;
- original rho/p, complete lag/null arrays, screened decisions and unresolved
  states, including boundary, zero, tied and unstable lag peaks.

GPU reductions, distance arithmetic and nearest-neighbor tie order can differ
even at float64. Record exact equality first, then absolute/relative errors and
the specific changed components. A numerical tolerance is a scientific
acceptance decision, not a way to hide a failed comparison. Matching only
significance arrows is insufficient.

The RTX 4080's graphics/AI throughput is not its float64 throughput. NVIDIA's
[Ada architecture description](https://images.nvidia.com/aem-dam/Solutions/geforce/ada/nvidia-ada-gpu-architecture.pdf)
shows the much smaller FP64 resources. Also, each embedding here is only
three-dimensional and about 1400 points: exact CPU KDTree search is already a
reasonable algorithm. GPU batching may help, but transfer, selection and
launch overhead must be measured rather than inferred from advertised FLOPS.

Use a dedicated GPU process with bounded batches if a compatible backend is
developed. Do not let all 32 outer workers create separate GPU contexts and
compete for one device. Compare complete pipeline throughput with the existing
32-worker CPU configuration, counting transfers, synchronization, serialization,
memory pressure and warm-up separately. No repeat of the 24-worker benchmark
is needed.

## What increasing the current 10 members means

The current 10 are chronology members for each setting, not 10 independent
process histories. Each process form still has one history (`r0`). For N
chronology members and R independent histories per process form, the four-level
study has `4 x 16 x 3 x N x R = 192 N R` input pairs.

| N | R per process form | Input pairs | Meaning |
|---:|---:|---:|---|
| 10 | 1 | 1,920 | Current pilot |
| 50 | 1 | 9,600 | More precise support for the same histories |
| 100 | 1 | 19,200 | Further within-record Monte Carlo precision |
| 50 | 10 | 96,000 | Both chronology and independent-history variation |

These are planning examples, not a declaration that a particular N or R is
sufficient. At conditional q=0.9, an independent Bernoulli approximation gives
Monte Carlo SE about 0.095 for N=10 and 0.030 for N=100. This describes
within-record simulation precision, not causal accuracy across histories.
The shared clocks make comparisons across settings paired.

The previous 02-to-02_3 sequence took approximately 130 minutes for 1920
pairs. A linear extrapolation to N=100, R=1 is about 22 hours on the existing
CPU path. This rough estimate excludes changes in workload, preparation,
resource competition and I/O; it is not a promised completion time.

Before a larger study, examine nested N=10/30/50/100 prefixes of a prespecified
full ensemble for q, native-score quantiles and graph-state stability. Treat
this as a pilot diagnostic, not a rule to stop as soon as a favorable q is
reached. Then allocate additional budget to independent histories and quantify
uncertainty at that independent-history level. More histories require explicit
seed and source handling; changing only `N_CAUSAL_MEMBERS` cannot provide them.

The current resume guards intentionally reject a changed member count and
require a matching completed 02 run before importing level 2. Preserve the old
outputs and create a fresh, consistently configured 02/02_3 run; do not edit
saved metadata to force compatibility. No larger production run was launched
as part of this review.

## Implementation order for a later production change

1. Validate native geometry reuse on representative accurate-age and uncertain
   inputs across all process forms and levels, including tie/zero-distance and
   failure cases. Compare complete arrays and nested outputs, not only arrows.
2. Implement bounded cache lifetime or process one library's lags together to
   control per-worker memory. Verify Windows spawn, restoration after errors,
   and same-input determinism. Keep the native reference available.
3. Benchmark the complete four-method workload at 32 workers against the
   original path on identical prepared inputs; inspect memory and paging.
4. Only then assess a compatible GPU backend on the same batch, with an
   explicit numerical acceptance contract. A newly written estimator requires
   a separate scientific decision under the project guide.
5. Freeze the accepted implementation and plan N and R separately before
   starting a fresh larger 02/02_3 study.

## Follow-up: GPDC notebook and paper check, 2026-11-07

The root `02_5_compare_gpdc_backends.ipynb` now compares six fixed conditional
tests from the three saved original uniform records. It uses an isolated
`.venv-gpdc` with CUDA-enabled Torch 2.5.1 and GPyTorch 1.13. All 24
case/backend combinations completed. Aligned CPU and GPU fits differ from the
current fixed-kernel GPDC by at most about 2.4e-14 in full regression vectors;
all six GPDC scores and p values agree exactly. The stock trained float32
GPDCtorch does not reproduce those regression vectors. Single-pass times were
about 0.28 seconds per fixed test for current GPDC and 0.10 seconds for aligned
CUDA. This is neither a full PCMCI+ validation nor production-throughput
evidence, and no production backend has been replaced.

Reading Van Nes et al. (2015), DOI 10.1038/nclimate2568, Methods and Fig. 3
on printed page 447 confirms **500 bootstrapped library sets of length 100**:
500 is the lag resampling count; 100 is the library length. Their separate
swap-null analysis uses 100 surrogate records. Their lag analysis uses shifts
up to 10 kyr, Gaussian smoothing (bandwidth 5), 5th/95th percentiles and a
bootstrapped paired-samples test. The four-page main paper does not explicitly
state replacement or whether every lag uses identical sampled indices; the
supplement was not obtained. Our 100 paired, without-replacement libraries
and 90% negative-peak rule are an adaptation, not a reproduction of that
bootstrap significance test. The notebook does not change the CCM procedure.

## Follow-up: bounded query reuse notebook, 2026-11-07

The user selected fixed L=300 for future work. The new root notebook
`02_7_check_ccm_neighbor_reuse.ipynb` defines `scan_ccm_reusing_neighbors`
locally; no production or toolbox function has been replaced. It changes the
loop order to library then lag, retaining at most two exact native KDTree query
results at once. All filtering, weights, predictions and correlations stay native.

Three saved original uniform records were compared at 100 libraries and 31 lags,
with two timing passes in opposite order. All 37,200 rho comparisons, full array
dtypes/metadata and 1,200 peak rows matched exactly, with maximum error zero.
Small checks also matched full prediction/target vectors and raw neighbor arrays,
including duplicate states and zero distances; a saved nonlinear level-4 input
and restoration after an injected query failure passed.

Each scan used 200 real searches and 6,000 cache hits. Median native/reuse times
were 21.21/5.58 s (linear), 21.82/5.46 s (nonlinear), and 23.00/5.58 s
(independent), approximately 3.8--4.1x for the lag stage. Peak cached query arrays
plus key bytes occupied 1.122 MiB, excluding the retained native models and other
overhead. This is a sequential CPU lag-stage experiment, not full CCM or
32-worker throughput evidence. Full arrays, checks, figures and an executed
notebook are saved under the corresponding `data/processed/` and `figures/`
directories. The temporary module patch remains unsuitable for concurrent
threads; production migration and its acceptance checks remain future work.

## Toolbox migration, 2026-11-07

After freezing 28 source/support files in
`data/archives/code_freeze_20261107T200542517660Z/source.zip`, the user requested
the CCM helper's migration. It now lives in `toolbox/causal.py` beside the lag
filter and is reached through the existing `run_ccm` call path. Unlike the
notebook prototype, it retains per-library/per-lag failures, preserves the
original/null evidence when a scan fails and restores the native factory on
interruption. The 02_7 notebook remains byte-identical to the freeze.

The four production scripts now explicitly set L=300, with updated implementation
and configuration guards; existing result directories have not been changed and
cannot be resumed as the new experiment. Five complete same-input comparisons
against the frozen code passed exactly, including original/null arrays, lag
arrays, tables, decisions and diagnostics. Two spawned workers and isolated
failure/repeated-state/interruption controls passed. The acceptance report is
`data/benchmarks/ccm_migration_checks/20261107T201205780934Z/report.json`.
Single-process complete CCM improved by approximately 1.85--1.96x; 32-worker
throughput has not been remeasured. GPDC remains on CPU. The recommended next
step is a fresh N=10 integration run before increasing chronology counts;
GPU regression still needs full PCMCI+ search and scheduling validation.
