# Experiment 02_4: sampling sensitivity of the bivariate causal pipeline

## Question and comparison contract

Measure how observation density, gap irregularity, asynchronous X/Y sampling,
and the age distribution of sparse sampling change the existing causal results.
This is a sensitivity analysis of the fixed sampling-plus-linear-interpolation
pipeline. An interpolated grid value is not an additional independent observation.

Use the same linear, nonlinear and independent process definitions as experiment
02: seed 20261009, 30 kyr retained history, 20 yr native step, 20 yr true delay,
500/200 yr X/Y memories, K=4 for coupled processes, and zero measurement noise.
Generate one realization per process form, with innovations paired across forms.
All ages are exact. There are no chronology ensembles in this branch.

All designs sample original native values within 29–1 kyr BP, include both
window endpoints as observations, and interpolate to the same 1401-point,
20-year grid without extrapolation. Endpoint inclusion is an explicit new
sampling-control condition, not a change to experiment 02's observation model.
The branch computes its own original-uniform reference; it does not modify or
force-reuse the completed 02/02_3 outputs.

## Fixed experiment matrix

| Stage | Designs per process | Repetitions | Input pairs across three processes |
|---|---|---|---:|
| Density | Synchronous regular 20, 40, 80, 160 yr intervals | One deterministic design per spacing | 12 |
| Irregularity and synchrony | Mean 80 yr; CV 0.25 and 0.75; AA, BB, AB, BA | 10 independent A/B mask groups | 240 |
| Age-dependent density | The same two gap multisets; old-sparse and young-sparse; AA and BB | The same 10 paired mask groups | 240 |
| Total | | | 492 |

For CV=0.25 use exactly 175 gaps of 60 yr and 175 of 100 yr. For CV=0.75
use exactly 175 gaps of 20 yr and 175 of 140 yr. Each irregular design has
350 intervals, 351 observations, a 28 kyr span and an 80 yr mean gap.
CV is the population SD of the intervals divided by their mean.
Regular designs have 1401, 701, 351 and 176 observations respectively.

Generate masks without looking at process values or causal outcomes. Share
each mask across the three process forms. A and B use independent seeds;
the same A/B seeds are reused across gradient modes and gap multisets.
Thus AA and BB are synchronous controls, while AB and BA exchange which
record receives each mask. Do not repeat a deterministic regular design
ten times and count the copies as different sampling realizations.

## Age-gradient construction

Place a fixed number of long intervals among the 350 interval positions.
For normalized position u from 0 (old end) to 1 (young end), assign log-weight
beta*(1-2*u), with beta=2.0 for old-sparse, beta=-2.0 for young-sparse and
beta=0 for the no-gradient control. Add the same per-position random Gumbel
draws within a paired A/B seed and select the required number of largest
scores. Put short intervals in the remaining positions. Counts and the full
gap histogram are therefore identical in all three modes.

The gradient is specified in interval order, not as an exact linear density
function of BP age. Save and plot actual density in equal-width age bins.
Random variation is retained; do not reject a mask because its causal result
or a small-sample gradient diagnostic is unfavorable. No observed ages or
proxy values are sorted or reassigned: only the proposed gap arrangement is
constructed before taking observations from the original process.

Mixed X/Y age gradients and additional independent process histories are
later extensions. The present age-gradient comparison keeps X/Y synchronous.

## Estimators and execution

Keep PCMCI+ GPDC (max lag 5, alpha/pc_alpha 0.05, 10000 null draws,
optimizer=None), corrected native LKIF (max lag 5), Granger (fixed order 1),
and PyEDM 2.5.7 CCM unchanged. CCM uses E=3, tau=-1, exclusion=10,
library fraction=0.8, 10 original libraries, 100 swap surrogates, and the
complete 100-library by 31-shift by two-direction lag screen.

Use method seed 20261309 throughout the paired sampling comparisons. Retain
the existing native GPDC cache for T=1391 with matching metadata and the
existing native CCM object reuse. Every input still receives its own fits.
Use 32 outer spawned processes, one BLAS/OpenMP thread per process, and no
internal PyEDM process pool. Prepare inputs and the null in the parent.
Use separate per-case output files; commit the eight-direction edge table
only after its native results and method diagnostics have been written.

Resume complete cases, including explicitly recorded failed fits. Reject
scientific-parameter, package, implementation-revision or numerical-platform
mismatches; compare regenerated inputs before resume. Worker count and
run-control flags may change. Do not create a historical cross-run cache.

## Saved evidence and summaries

Save the original processes, the canonical sampled indices/ages for every
mask, case-to-mask assignments, and every interpolated analysis input. These
are sufficient to reconstruct the exact observed values without duplicating
the same record for each AA/AB combination. Save per-record and pair sampling
diagnostics, native outputs, all directional rows, and explicit failures.

Report actual observation count, mean/CV/max gap, age-binned density, shared
X/Y observation fraction and nearest observation-time differences. Assess
interpolation using RMSE/original SD, correlation, SD ratio, increment
variance ratio, and lag-one autocorrelation.

Keep both directional native scores and all finite successful scores, whether
significant or not. Preserve CCM raw Tp=0 rho/p/votes, screened votes and all
lag arrays/profiles/seeds/flags. Preserve native PCMCI matrices/graphs, LKIF
arrays and GC regression arrays/statistics. A failed lag screen can coexist
with a usable native Tp=0 score.

Within each process/scenario/method report usable paired decisions, forward
and reverse retention, unresolved evidence, failed/pending counts, and exact
graph recovery divided by the requested case count. Correct graphs are X→Y
only for coupled forms and no arrows for the independent form, with no
unresolved evidence. These are conditional sampling frequencies for one
process history, not independently calibrated error rates.

For each mask group and each native directional score S, save

- asynchronous contrast: (S_AB + S_BA - S_AA - S_BB)/2;
- old-sparse contrast: mean(S_old_AA, S_old_BB) - mean(S_none_AA, S_none_BB);
- young-sparse contrast: mean(S_young_AA, S_young_BB) - mean(S_none_AA, S_none_BB).

Require every component of a contrast to be finite and successful, and
report the number of usable contrasts out of the ten requested groups.
Keep the individual combinations; their effects can cancel in an average.
Do not average p-values or interpret the ten mask groups as ten climate
histories. Use paired differences, medians and IQRs rather than ratios to a
near-zero baseline or confidence intervals across dependent combinations.

## Figures and acceptance

Produce sampling-layout and age-density diagnostics, a reconstruction example,
density-versus-native-score panels, full-graph recovery fractions, paired
native-score contrasts, CCM raw-versus-screened retention, and the individual
AA/BB/AB/BA native-score distributions against the regular 80-year reference. Use older
ages on the left and label Age (kyr BP). Mark missing/failed/partial cases.
Keep plotting in the root analysis script, with PNG and PDF outputs under
figures/02_4_check_sampling_sensitivity/ and tables/native evidence under
data/processed/02_4_check_sampling_sensitivity/.

Before the full run, check gap multiset/count/span invariants, reproducibility,
paired mask reuse, zero-error full-grid interpolation, complete contrast
arithmetic and failure retention. Run a small real-estimator integration
check on identical inputs and compare full native numerical outputs with
direct existing-wrapper calls. Then run all 492 cases with frozen settings,
inspect method statuses and figures, and verify same-run resume without
recomputing committed cases. Temporary checks belong in tmp/.

Reference: Anderson et al., On the Sensitivity of Granger Causality to
Errors-In-Variables, Linear Transformations and Subsampling,
https://doi.org/10.1111/jtsa.12430. This motivates explicit sampling controls;
the outcomes of the current four-method experiment remain empirical.

## Execution record: 2026-11-04

All 492 input pairs completed on Windows using 32 workers in 3739.757 seconds
(62.329 minutes; 473.614 pairs/hour, including preparation and figures).
There are 1968 method fits, 3936 directional rows and 1968 native archives.
All fits succeeded; every CCM direction has 100 valid lag libraries, and all
1440 paired score contrasts are usable. Full artifact validation checked
native array shapes, GC fit fields, summary denominators and contrast arithmetic.
The production duration is retained in the first artifact-audit report under
`data/benchmarks/sampling_checks/20261104/`; `run_status.json` reports the most
recent invocation and its elapsed time, including a later resume-only run.

System CPU use averaged 93.98% with a median of 100%. Minimum available
physical memory was 9.67 GiB; maximum sampled commit was 29.34/63.82 GiB.
No sustained memory-pressure pattern was observed. This workload is different
from the earlier worker benchmark, so its throughput is not a speedup comparison.

The full native-output comparison on two cases passed exactly before production.
Resume validation identified last-digit drift from the default CSV result parser;
committed edge tables now use `float_precision="round_trip"` when rebuilding
summaries. Native estimator inputs and per-case numerical outputs are unchanged.
The final resume audit (`artifact_audit_20261104T215527447629Z.json`) passed:
all 2460 committed member/native files retain identical byte hashes and
modification times, and rebuilt aggregates exactly match committed edge values.

### Conditional findings from the fixed histories

- Sampling density and interpolation already affect direction recovery. All
  methods recover the true graph on the original 20-year controls; coarser
  controls introduce reverse GC arrows and unresolved PCMCI+ evidence in
  coupled systems. Native scores can increase without improving graph recovery.
- Asynchrony is not uniformly detrimental. In the independent history at
  CV 0.75, PCMCI+ correct complete graphs increase from 4/20 synchronous cases
  to 18/20 asynchronous cases. These counts pool AA/BB or AB/BA within the same
  ten paired mask groups and are not independent-history error rates.
- Relocating the same intervals changes results. LKIF correct complete graphs
  for synchronous masks are shown below; every cell has 20 requested cases.

| Process and gap CV | No gradient | Old-sparse | Young-sparse |
|---|---:|---:|---:|
| Linear, 0.25 | 10/20 | 11/20 | 7/20 |
| Linear, 0.75 | 10/20 | 8/20 | 8/20 |
| Nonlinear, 0.25 | 2/20 | 3/20 | 1/20 |
| Nonlinear, 0.75 | 12/20 | 11/20 | 9/20 |
| Independent, 0.25 | 14/20 | 6/20 | 16/20 |
| Independent, 0.75 | 11/20 | 10/20 | 1/20 |

- Screened CCM recovers the complete true graph in all 492 cases here. Raw
  Tp=0 CCM retains both directions in all 328 coupled cases; the existing lag
  screen removes the reverse arrow. Raw scores and screened decisions must
  therefore remain separate. The 164 independent cases have no raw CCM arrows.
- PCMCI+ has 361 unresolved case states. These remain explicit evidence and
  do not count as complete-graph recovery, even when X to Y is retained.

These observations are conditional on one history per process form. They do
not establish general method rankings or calibrated false-positive rates.
Seven PNG/PDF figures and all numerical tables are saved in the experiment's
figure and data directories described above.
