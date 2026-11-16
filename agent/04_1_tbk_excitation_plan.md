# TBK excitation and forced-response diagnostic

Authorized on 13 November 2026. Entry point: `04_1_check_tbk_excitation.py`.

## Question and scope

The completed 04 chronology pilot has unresolved contemporaneous PCMCI+ links, dense retained lagged graphs, structured GC residuals, and high baseline CCM reconstruction skill without equally clear timing evidence. This diagnostic asks whether strong prescribed forcing and weak independent stochastic excitation contribute to these outcomes while the physical memory and interaction coefficients remain fixed.

Use accurate ages only. Do not add chronology perturbations, change diffusion or restoring parameters, introduce measurement error, tune a method after seeing its result, or expand this run into full-network recovery. The subsequent information-set, sampling-resolution and renewed chronology comparisons remain later decisions. The high-noise condition is an identification diagnostic, not a physically calibrated parameter choice.

## Paired inputs

Retain the existing three seeds, 20261009, 20261301 and 20261302. Integrate the unchanged four-box core for 39,000 annual updates, discard the first 9 kyr, sample every 20 years, and analyze the same 29–1 kyr BP grid of 1401 points. X is southern T1 and Y is northern T4. Each seed retains estimator seed `process_seed + 400` for all conditions.

| Condition | Annual increment SD | Analyzed trajectory | Purpose |
|---|---:|---|---|
| `baseline` | 0.005 | Original forced temperature | Exact reference from completed 04 |
| `noise_x10` | 0.05 | Independently integrated forced temperature | Tenfold stochastic excitation |
| `noise_x100` | 0.5 | Independently integrated forced temperature | Stronger diagnostic excitation |
| `random_response` | 0.005 | Baseline temperature minus same-parameter noise-free temperature | Isolate the random perturbation response using model knowledge |

Use the same local RandomState seed and array layout at every amplitude. Reintegrate the model rather than scaling plotted lines. The initial conditions, prescribed q, diffusion, global/area restoration, update order and duration are identical. For this linear temperature model, differences from the noise-free trajectory should scale with the noise amplitude to floating-point accuracy; verify this invariant before fitting.

The random-response condition is a model-assisted diagnostic, not a proposed preprocessing step for real proxy data. It retains the same temperature transition operator. The noise-free reference is not a thirteenth causal input. Three seeds provide paired descriptive comparisons, not population-level false-positive calibration. No chronology members are present in this experiment.

## Estimators and reuse

Run all four existing pairwise methods on all twelve inputs, preserving the complete 04 settings:

- PCMCI+/GPDC: lags 0–5, alpha and pc_alpha 0.05, 10,000 native reference null draws, `optimizer=None`, existing directional Holm families and native-graph screening.
- LKIF: corrected native package, maximum lag five and unchanged analytical significance/normalization.
- GC: fixed VAR order one, existing SSR-F test and native score; retain fitted coefficients, residuals, stability and whiteness diagnostics.
- CCM: E=3, tau=-1, temporal exclusion ten, L=300, ten original libraries, 100 swap surrogates, 61 target shifts from -30 to +30 steps, 100 paired lag libraries and the existing 90% sign rule. Preserve raw and screened evidence separately.

Copy the three baseline input files and their twelve native archives only after checking regenerated inputs, the actual CSV-loaded estimator arrays, seeds, packages, numerical platform, source configuration and archive completeness. Do not refit accurate baselines or selectively discard unfavorable evidence. Reuse the existing matching GPDC reference-null caches with their normal metadata checks.

In addition, run native PCMCI+ with ParCorr on the three original baseline inputs only. Keep the same search limits, thresholds, Holm handling and unresolved-state semantics. This is a separately labeled conditional-test diagnostic, not a replacement for GPDC and not a fifth algorithm ranked against the other four. Compare the independent native call with the wrapper's full p/value/graph arrays and derived decisions.

The primary experiment contains 48 method/input records and 96 directional rows. The ParCorr comparison adds three records and six directional rows. Twelve native records are source copies; 39 are new fits. Keep all 51 full native archives and all 102 directional rows.

## Interpretation and outputs

Do not use “both directions without unresolved evidence” as a cross-method accuracy score. PCMCI+ and CCM have explicit unresolved categories that current GC/LKIF wrappers do not share. Keep directional flags and unresolved evidence visible simultaneously.

The comparisons should expose:

1. The actual temperature traces and random-response amplitude at fixed forcing, plus exact dynamics/decay and signal-versus-random-response diagnostics.
2. South-to-North and North-to-South evidence for every history and condition, with separate contemporaneous/lag-screen uncertainty.
3. PCMCI+ lag-specific p values, scores and graph symbols from lag zero through five, showing whether all-lag retention changes; compare baseline GPDC with ParCorr without combining score scales.
4. Original CCM rho/significance, lag curves and per-library peak categories, including boundary/zero/tied/unstable peaks.
5. GC residual serial and contemporaneous dependence, together with its directional results. Residual diagnostics qualify interpretation; they do not delete fitted results.

Keep original native score definitions. In particular, maximum GPDC may occur at an unoriented zero-lag link, and raw CCM skill is Tp=0 rather than the maximum scanned value. Retain per-history results and all failures/unresolved states. Any medians or ranges across three histories are descriptive, not confidence intervals.

Save inputs, prepared model trajectories, native archives, complete case commits and scientific summary tables under `data/processed/04_1_check_tbk_excitation/`. Export labeled PNG/PDF figures under the matching `figures/` directory. English captions belong in `agent/figure_review/10_tbk_excitation.md`; the results report belongs in `agent/04_1_tbk_excitation_results.md`. A Chinese reading note may be placed only in the Git-ignored `Note/` directory.

## Execution and acceptance

Use the Windows project virtual environment and the existing 32-worker limit with one numerical thread per worker. Schedule only the available independent jobs; do not repeat worker benchmarks. Native method jobs may be scheduled independently to use the available workers. Keep native internal CCM pools disabled.

Before production, validate all prepared inputs, source-reuse equality, paired noise/forcing invariants, direct native ParCorr equality and Windows-spawn behavior. Preserve the old 04/04_0 sources, inputs, native results and figures by recording their hashes, sizes and modification times. Source-result guards must reject changed input/settings or incomplete native evidence on resume.

After production, audit all 51 native archives, requested directional rows, CCM library arrays and GC scientific fits. Reconstruct saved summaries from native evidence, review the actual figures, and verify the old scientific artifacts remain unchanged. A completed computation is not a resolved causal conclusion. Report any remaining baseline ambiguities and estimator limitations, including outcomes that contradict the motivating hypothesis.
