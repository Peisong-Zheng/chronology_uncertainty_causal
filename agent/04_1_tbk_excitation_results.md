# TBK excitation diagnostic: completed results

Completed on 13 November 2026 with `04_1_check_tbk_excitation.py`. The [plan](04_1_tbk_excitation_plan.md) specifies the comparison; the [figure captions](figure_review/10_tbk_excitation.md) describe all five exports.

The result supports a substantial role for prescribed forcing and the amount of independent stochastic excitation in the endpoint evidence. It does **not** show that increasing process noise recovers the physical feedback, or establish intrinsic memory as the cause of the original ambiguities.

## What was held fixed

All inputs have accurate ages, the original annual four-box temperature equations, the same prescribed 1500-year transport cycle, 39 kyr integration, 9 kyr startup, and the same 1401-point, 20-year analysis grid. The three histories share forcing and differ in stochastic innovations. Each history uses paired standardized innovations across noise amplitudes and the same estimator seed across conditions.

| Figure code | Saved scenario | Input |
|---|---|---|
| A | `baseline` | Original forced trajectory, annual process-noise SD 0.005 |
| B | `noise_x10` | Reintegrated forced trajectory, SD 0.05 |
| C | `noise_x100` | Reintegrated forced trajectory, SD 0.5 |
| D | `random_response` | A minus its exactly known noise-free trajectory |

The three original A inputs and twelve native archives were copied exactly from completed 04. The nine other inputs received the same four primary methods. Three additional native PCMCI+/ParCorr fits use A only; they are a conditional-test diagnostic, not a production replacement.

For the unchanged temperature model, write

\[
\mathbf T_{t+1}=A\mathbf T_t+\mathbf bq_t+s\boldsymbol\epsilon_t,
\qquad
\mathbf T^{(0)}_{t+1}=A\mathbf T^{(0)}_t+\mathbf bq_t.
\]

Subtracting the same-initial-state noise-free trajectory gives

\[
\mathbf R_t=\mathbf T_t-\mathbf T_t^{(0)},\qquad
\mathbf R_{t+1}=A\mathbf R_t+s\boldsymbol\epsilon_t.
\]

Thus neither increasing s nor forming D changes A or its memory eigenmodes. Independent verification recovered e-folding times of approximately 154, 99, 9.2 and 2.8 years. The prescribed 1500-year recurrence is a different timescale. D is available because the generating model is known; it is not an ordinary detrending prescription for observed proxies.

## Input changes explain why this comparison matters

Across the three histories, SD(random response)/SD(noise-free response) is only 0.00937–0.00980 in the South and 0.00332–0.00348 in the North at baseline. B multiplies those ratios by ten. In C the ratios reach 0.937–0.980 and 0.332–0.348, respectively: even the hundredfold-noise North record still contains a strong forced response.

The baseline 1500-year autocorrelation is 0.99990–0.99991 in the South and approximately 0.99999 in the North. In D it falls to −0.064–0.007 and −0.066–0.044. The 20-year autocorrelation remains 0.855–0.875 and 0.750–0.766. Removing the known forced cycle therefore removes most of the recurrence at its period while retaining substantial short-term physical persistence. These are empirical autocorrelations, not estimates of the model's intrinsic e-folding times.

## Directional evidence

Counts below are retained directions out of three histories. Each cell lists **South to North / North to South**. CCM is the lag-screened result. These counts describe the fitted evidence; they are not calibrated accuracy scores or three independent forcing histories.

| Condition | PCMCI+/GPDC | LKIF | GC | CCM |
|---|---|---|---|---|
| A: original | 3/3 ; 3/3 | 0/3 ; 3/3 | 3/3 ; 3/3 | 0/3 ; 2/3 |
| B: noise ×10 | 0/3 ; 3/3 | 0/3 ; 3/3 | 3/3 ; 3/3 | 0/3 ; 3/3 |
| C: noise ×100 | 0/3 ; 0/3 | 0/3 ; 3/3 | 2/3 ; 3/3 | 0/3 ; 3/3 |
| D: random response | 0/3 ; 1/3 | 1/3 ; 1/3 | 0/3 ; 2/3 | 0/3 ; 0/3 |

Unresolved evidence is separate. GPDC has a contemporaneous `o-o` link for all three A histories and `x-x` conflicts for all three B histories; C and D have none. CCM has raw-significant directions with unresolved timing in all three A and D histories. A retained arrow can coexist with unresolved evidence. GC and LKIF wrappers do not assess the same unresolved category.

### PCMCI+ and the conditional-test control

In A, GPDC retains **every tested positive cross-lag, 20–100 years, in both directions**, in all three histories, while leaving lag zero unoriented. It is incorrect to summarize that baseline as no directional evidence.

In B, every South-to-North lag disappears. North-to-South retains 20 years in all three histories and additionally 60 years in one; lag-zero conflicts remain. C retains no cross-links. D retains only North-to-South at 20 years in history R2. Increasing noise therefore changes the dense baseline graph without simply restoring both physical directions.

Using ParCorr on exactly the A inputs retains North-to-South in all three histories and South-to-North in none. At lag zero, two histories have an oriented North-to-South link and one has a conflict. This demonstrates sensitivity to the conditional-independence test at fixed data and physical memory. It does not establish that the ParCorr graph is the correct physical graph. Native GPDC statistics and signed partial correlations have separate scales in Figure 02.

### CCM

Original Tp=0 reconstruction skill, across both directions and all histories, changes as follows:

| Condition | Native rho range |
|---|---|
| A | 0.99883–0.99984 |
| B | 0.98043–0.99071 |
| C | 0.43138–0.62882 |
| D | 0.00353–0.09777 |

B and C retain North-to-South after the lag screen in all histories and reject South-to-North as temporally incompatible. D has raw significant evidence in one South-to-North and all three North-to-South tests, but no direction passes the lag screen. Low reconstruction skill can still be significant relative to this particular surrogate null; significance, skill and timing are different quantities.

The collapse in skill after exact forced-response subtraction is evidence that the almost perfect baseline reconstruction relies strongly on the prescribed response. It is not proof that CCM always measures common forcing or that the remaining physical coupling is absent. The random response is less predictable under the unchanged embedding/library settings.

### GC and LKIF

GC retains both directions in all A and B histories. However, its residual contemporaneous correlation is −0.946 to −0.943 in A, and all three residual-whiteness tests reject at extremely small p values. In B the correlation magnitude falls to 0.317–0.367, but all three whiteness tests still reject. In C and D, residual correlation falls to −0.082 to −0.050 and −0.044 to −0.021, respectively; whiteness still rejects in two of three histories for each condition. Directional significance alone therefore does not demonstrate adequate VAR(1) fit.

LKIF consistently retains only North-to-South for A/B/C. For D, R1 retains South-to-North, R2 retains North-to-South, and R3 retains neither. The fitted direction is less stable once the forced response is removed. These are descriptive endpoint outcomes, not a calibrated false-positive assessment.

## What to do next

The present evidence narrows the question: how much of the endpoint interpretation comes from the forced cycle, the observed variables, the 20-year sampling and the estimator's conditional test? A generic statement that PCMCI+ fails because the model has memory is not supported by these controls.

The smallest useful next check would apply the same native ParCorr diagnostic to the three D inputs. If further separation is needed, compare an explicitly specified full-state linear reference with the endpoint-only analysis using the saved four-box trajectories and known forcing. The endpoint pair need not form a closed VAR(1), even though the full annual temperature model is linear. A sampling-resolution comparison should remain a separate factor. None of these follow-ups, new chronology ensembles or memory sweeps was launched in this run.

Only three paired innovation histories and one prescribed forcing schedule were used. The high-noise condition is diagnostic rather than physically calibrated. Do not select a future noise amplitude merely because it produces preferred arrows. Before extending chronology uncertainty, specify the mechanism contrast and the accurate-age observable benchmark against which a change in interpretation will be judged.

## Files, runtime and verification

- Data: `data/processed/04_1_check_tbk_excitation/`; primary tables include `member_results.csv`, `parcorr_results.csv`, `scenario_summary.csv`, `pcmci_lag_diagnostics.csv`, `ccm_lag_summary.csv`, and the input/GC diagnostics.
- Figures: `figures/04_1_check_tbk_excitation/00` through `04`, as PNG/PDF pairs. Start with Figure 00, then Figure 02 for PCMCI+, Figure 03 for CCM, Figure 04 for GC, and Figure 01 for the combined overview.
- Preflight: `data/benchmarks/tbk_excitation_checks/20261113T121452115736Z/report.json` and `final_source_acceptance.json`.
- Production: `data/logs/tbk_excitation_20261113T121730984237Z/`, containing the source capture, resource log, numerical audit and preservation report.
- Completed 12/12 primary inputs plus 3/3 ParCorr controls: 51 native archives, 102 successful directional rows, 39 new method calls. Main unresolved evidence is retained in 20 directional rows; it is not a failed fit.
- Experiment elapsed time: 71.771 seconds with a 32-worker limit and one numerical thread per worker. Minimum sampled available physical memory: 10.004 GiB. This is a small-run timing with existing GPDC caches and reused baselines, not a new worker-throughput benchmark.
- Audit passed complete CCM 100×61×2 lag arrays and 100×2 null arrays, native PCMCI graphs/Holm decisions, GC regression fits/residuals, and exact reconstruction of six scientific plus seven figure tables. Twelve copied baseline archives remain byte-identical. All 1,026 inventoried pre-existing TBK source/data/figure files retain hashes, sizes and modification times.
