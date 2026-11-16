# Memory sensitivity of chronology effects

Approved implementation scope: 9 November 2026. The user requested all four causal methods after reviewing the compact CCM-focused design. Main experiments and their saved results remain unchanged.

## Question and discriminating comparison

Does the excellent CCM result in 02_3 depend on long process memory? More generally, does shorter memory increase chronology-induced loss for PCMCI+, LKIF, CCM and Granger causality, or does it already reduce performance with accurate ages?

The existing CCM result is close to a ceiling: 1918/1920 complete graphs are correct. All 1280 coupled inputs have significant raw Tp=0 tests in both directions; the lag screen rejects all reverse arrows and retains 1278 forward arrows. Accurate-age mean-curve optima are -180 yr for forward linear, -160 yr for forward nonlinear, and +120 yr for both reverse directions. This motivates examining the original test and the lag screen separately. It does not establish memory as the cause of these outcomes.

For every memory condition, compare perturbed records with their own irregular-sampling accurate-age reference. Preserve the original-uniform reference to identify sampling/interpolation effects separately. Keep unfavorable or unresolved accurate-age results rather than selecting only successful baselines.

## Fixed experimental matrix

| Memory label | X intrinsic memory (yr) | Y intrinsic memory (yr) |
|---|---:|---:|
| reference | 500 | 200 |
| short_x | 100 | 200 |
| short_y | 500 | 40 |
| short_xy | 100 | 40 |

Use process forms linear, nonlinear and independent, each with process seeds 20261009, 20261301 and 20261302. Within a seed, reuse standardized innovations across memory conditions and process forms. Fixed 125-step warm-up ensures that the retained innovations do not shift when memory changes. Each form/seed/memory has:

- One original-uniform accurate-age input.
- One irregular-sampling accurate-age input.
- Ten chronology members using the existing level-4 C-C clocks, with 80-year per-record proposal SD.

Total: 4 memory conditions × 3 forms × 3 seeds × (2 references + 10 members) = **432 input pairs**. All four methods analyze every pair: **1728 method calls and 3456 directional rows**. No intermediate uncertainty level or additional uncertainty type is introduced in this first run.

Three process seeds supply three process realizations, with paired conditions within each seed. Ten chronology members are within-record sensitivity draws, not ten independent process histories. The same clock set is deliberately shared across seeds, so inference is conditional on that set of chronology perturbations.

## Process and chronology controls

Keep the 20-year native and analysis grid, 30-kyr retained process, 29-1 kyr common analysis window, one-step true delay, theoretical background SDs, normalized driver transformations and zero measurement noise. Keep long-run K=4 for coupled forms and K=0 for independent controls. For each branch memory, derive

$$
a_i=\exp(-\Delta t/T_i),\qquad c=(1-a_Y)K,
$$

$$
\sigma_u=S_X\sqrt{1-a_X^2},\qquad
\sigma_v=S_Y\sqrt{1-a_Y^2}.
$$

These changes belong to the new branch only. Fixing K and unforced background SD does not fix the total variance of coupled Y. Save total variability, empirical lag-one correlation and the forced/background SD ratio. With paired innovations, the forced component can be diagnosed as coupled Y minus independent Y at the same seed and memory, without extra causal fits.

Read the existing level-4 C clock NPZ files from 02_3, together with their saved metadata. Reuse the exact sample identities, accurate ages, uncertain observed centers, ordered chronology members and proposal SDs. Attach each new process's proxy values to those fixed sample positions. Do not regenerate clocks around true ages. X/Y clocks remain independent of one another. The chronology error correlation length stays fixed at its current normalized-sample-order setting, approximately 500 years; it does not shrink with process memory.

All inputs cover the same 1401-point grid without extrapolation. Save realized clock diagnostics as well as the nominal 80-year SD; the existing C-C relative error RMS is approximately 137 years. Acceptance, coverage and effective spread remain explicit.

## Frozen methods

| Method | Settings and retained evidence |
|---|---|
| PCMCI+ | Native CPU Tigramite GPDC, max lag 5, alpha and pc_alpha 0.05, optimizer=None, 10000 null samples; native graph, full p/val arrays, directional Holm decisions and unresolved contemporaneous evidence |
| LKIF | Existing corrected author package, max lag 5, analytical significance and native normalization |
| CCM | Existing exact neighbor-reuse toolbox implementation; E=3, tau=-1, exclusion 10 steps, L=300, 10 original libraries, 100 swap surrogates, 31 shifts from -15 to +15, 100 paired libraries, current sign rule |
| GC | Existing statsmodels wrapper, fixed VAR order 1, same score and test; scientific OLS fit fields retained |

Use estimator seed = process seed + 300, shared across paired memory/form/chronology conditions. Prepare one native GPDC null cache per estimator seed at effective T=1391 and 10000 samples. Every input retains its own conditional regressions. Do not switch conditional tests, introduce GPDCtorch, tune estimator parameters by memory, reduce Monte Carlo counts or repeat worker benchmarks.

Keep 32 outer workers, one numerical thread per worker and no PyEDM internal multiprocessing. Reuse the established Windows launcher and temporary sleep prevention. Record resource pressure during execution.

Because embedding and exclusion times remain fixed, their ratios to process memory change. Interpret the result as sensitivity of the current analysis configurations, rather than a comparison of separately optimized methods at each timescale.

## Results and paired contrasts

For each method/direction, record all finite successful native scores, original significance where applicable, final arrows, unresolved evidence and failures. Do not compare absolute native-score magnitudes between methods.

The primary score contrast is

$$
\Delta s_{r,m,j}=s_{r,m,j,\mathrm{chronology}}
-s_{r,m,\mathrm{irregular\ accurate}},
$$

where r is process seed, m is memory condition, and j is chronology member. A contrast is usable only when both scores are finite and successful. Save individual paired differences, counts, seed-specific medians and IQRs. A matched contrast of these changes between each short-memory condition and the reference isolates differences in chronology sensitivity; all four component scores must be usable.

For every method, report forward/reverse support, complete graph recovery over the requested member count, unresolved states and failures. Support uses members with usable decisions in both directions. Correct coupled graphs contain only X to Y; independent graphs contain no cross-arrow; both require no unresolved evidence. Never treat a pending or failed fit as an absent arrow.

For CCM additionally retain raw versus screened votes, the complete 100×31×2 rho arrays, original nulls, per-library peaks, mean lag profiles, and past/future/zero/boundary/tied/failure fractions. Examine whether shorter memory reduces the accurate-age peak separation from zero and whether chronology perturbations destabilize the signs.

Preserve separate results for the three seeds in saved tables. Following the user's 9 November presentation revision, the main graphs show pooled descriptive fractions and paired-score distributions without seed-specific marker shapes. Per-history quantiles across ten chronology members describe chronology sensitivity. The revised pooled-score IQR includes both within-history chronology variation and between-history variation; neither is a confidence interval across independent histories. Do not calibrate population false-positive rates from three seeds.

## Script, outputs and figures

Implement `02_8_check_memory_sensitivity.py` at the project root. Save outputs under `data/processed/02_8_check_memory_sensitivity/` and figures under `figures/02_8_check_memory_sensitivity/`.

Retain case metadata and prepared inputs; process/clock diagnostics; four native JSON-gzip payloads per case; a committed eight-row edge table per case; aggregate member results; case states and scenario summaries; paired score changes and cross-memory contrasts; CCM lag summaries and profiles. The edge table is written last so partial native files cannot become complete votes. Resume checks scientific settings, source input fingerprints, packages and numerical platform. Atomic file replacement uses bounded retries for transient Windows sharing errors. Result CSVs are read with round-trip precision.

Figures should expose all four methods, both directions and reference-versus-chronology differences. Include a separate CCM lag/screening diagnostic and a process/clock diagnostic where useful. Keep informative numerical annotations, units, legends and panel labels. Preserve the established research appearance rather than imposing the rejected restyle. Detailed English captions belong in `agent/figure_review/05_memory_sensitivity.md`.

## Acceptance and interpretation

Before production, verify all 432 identities, fixed windows, matched clocks/masks/innovations and coefficient calculations. Compare the reference-seed inputs with their saved 02/02_3 counterparts, reporting any CSV serialization precision difference rather than silently calling it exact. Check two full-setting cases against direct native calls and two spawned workers using complete scientific payloads, excluding runtime only. Test failure retention, final commits, same-input resume and mismatch rejection. No main-study scientific files should change.

After production, audit all expected directional rows and native archives, complete lag-array shapes, recorded failures, summary denominators, paired contrasts and rendered figures. A failed scientific fit remains a recorded outcome. A serialization or bookkeeping failure requires repair and resume rather than a silent scientific failure label.

Interpret a chronology-specific memory effect only alongside each condition's accurate-age performance. If a short-memory reference already fails, it identifies a baseline/configuration limitation. Stable native scores with unstable lag signs implicate the CCM screen. Differences confined to short_x or short_y help locate the contribution of driver or response memory. A flat result limits the hypothesis only for these four memories and this C-C clock set.

The optional 20-year proposal-SD condition remains a separately considered follow-up; it is not part of this run.

## Execution outcome

The authorized experiment completed all 432 pairs on 9 November 2026 in 1458.999 seconds using 32 workers. All 3456 directional rows succeeded. See [the results and interpretation](02_8_memory_sensitivity_results.md), including accurate-age exceptions, all four methods, independent controls and the limits of the memory comparison. Detailed run and acceptance records are under `data/logs/memory_sensitivity_20261109T114908564938Z/`.
