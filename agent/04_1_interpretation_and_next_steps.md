# Interpreting 04_1 and choosing the next discriminating comparison

Discussion record, 13 November 2026. This note reinterprets completed 04_1 outputs and the unchanged model equations. No new causal fits, chronology members or model variants were run. The proposed cross-realization experiment below is not yet authorized or implemented.

## What the excitation comparison actually changes

The annual temperature model is linear in its state with prescribed external transport:

\[
T_{t+1}=A T_t+bq_t+s\epsilon_{t+1}.
\]

For the same initial condition and standardized innovations, its trajectory can be decomposed exactly as

\[
T^{(s)}_t=F_t+sR^{(1)}_t,
\qquad R^{(1)}_{t+1}=A R^{(1)}_t+\epsilon_{t+1}.
\]

Conditions A/B/C in 04_1 change the mixture of a fixed forced response and a random response, rather than changing the random system's coupling matrix or memory. Condition D removes F exactly and preserves A. In standardized variables, sufficiently large s approaches the random-response limit represented by D; this does not imply monotonic estimator behavior at finite mixtures. Repeating D at different uniform noise amplitudes would add little to scale-invariant diagnostics.

The model eigenmode e-folding times are 2.794, 9.231, 99.449 and 153.899 years. The 1500-year event recurrence is prescribed forcing, not an intrinsic decay time. D loses the nearly perfect 1500-year recurrence but retains 20-year autocorrelation of 0.855-0.875 in the South and 0.750-0.766 in the North. Weak results in D cannot by themselves be explained by its small absolute amplitude.

## A newly explicit physical benchmark

The annual transition matrix saved in the independent preflight is exactly symmetric, including diffusion and the global/area restoring terms:

\[
A=A^\mathsf{T},\qquad
(A^k)_{1,4}=(A^k)_{4,1}.
\]

Consequently, with forcing and subsequent innovations held fixed, equal endpoint temperature impulses in model units produce identical opposite-endpoint response curves at every horizon. Direct read-only evaluation gives both 20-year cross-response coefficients as -0.03207167671382581. This is a property of this specified model, not a general statement about real bipolar climate dynamics.

The forcing loading is asymmetric, b=(0,-1,0.5,0.5): northern T4 receives a direct forcing term, whereas southern T1 does not. Event-response timing, predictive information and equal-unit physical impulse response are therefore distinct objects. Repeated North-to-South evidence across methods is not evidence that the underlying equal-unit coupling is stronger in that direction. Standardization, different state/innovation variances, omitted middle boxes and external forcing all matter for detectability.

Temperature does not feed back to prescribed q. This model can test interpretation of temperature responses to imposed transport changes; it cannot establish that temperature feedback generates the imposed DO-event timing or AMOC transitions.

## What the method results support

- **PCMCI+/GPDC:** A retains all five positive cross-lags in both directions while lag zero is unoriented. B retains only North-to-South positive lags and has lag-zero conflicts; C retains no cross-links; D retains North-to-South in one history. This changes substantially at fixed physical memory, but does not show that adding noise restores feedback. Baseline ParCorr differs from GPDC, implicating conditional-test/search sensitivity without validating either physical graph.
- **Scope of the PCMCI claim:** the production configuration is PCMCI+ with fixed GPDC regression hyperparameters (`optimizer=None`). It is not all possible PCMCI+ or GPDC configurations. GP-regression underfitting is a possible explanation, not an observed finding. The [PCMCI+ paper](https://proceedings.mlr.press/v124/runge20a.html) explicitly addresses strong autocorrelation under causal sufficiency; memory alone is not an appropriate blanket explanation for failure.
- **CCM:** A has nearly perfect raw reconstruction, but the saved 100-library lag profiles improve over their own zero-lag values by only about 0.00005-0.00014 for South-to-North and 0.00036-0.00085 for North-to-South. High reconstruction skill therefore coexists with weak peak separation. B/C have lower skill but more stable lag-screen direction. D has low skill and no screened arrow despite some raw-significant results. Reconstruction, surrogate significance and temporal compatibility remain separate quantities.
- **CCM null:** cut-and-swap shifts one record relative to the other; it does not implement a physical null of no interaction conditional on shared q. Periodicity and circular splice structure can affect that test. The proposed control below would test a different and more relevant mechanism distinction. It has not yet established that the original high rho is wholly attributable to shared forcing. Related methodological evidence is discussed by [Cobey and Baskerville (2016)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169050).
- **GC/LKIF:** GC's baseline bidirectional retention coexists with very strongly correlated and non-white residuals. Removing F markedly reduces contemporaneous residual correlation, yet two of three D fits still reject whiteness. LKIF's stable North-to-South result for A/B/C becomes history-dependent in D. Neither stable arrows nor agreement across methods is a physical-effect strength measure.

The endpoint pair is also not generally a closed VAR(1). For all four D states sampled every 20 years,

\[
R_{n+1}=A^{20}R_n+\eta_n,
\qquad Q_{20}=s^2\sum_{j=0}^{19}A^j(A^j)^\mathsf{T}.
\]

Marginalizing T2/T3 changes the effective predictor history. Q20 is non-diagonal even though annual innovations are independent. Full-state 20-year VAR coefficients and residual covariance can be checked against A20 and Q20, but that is not automatically a causal-sufficiency benchmark for instantaneous PCMCI+ orientation. [Hyttinen et al. (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5944866/) distinguish system-timescale and measurement-timescale causal graphs.

## Recommended next core experiment: same versus different model realizations

Use the three existing histories without changing any single-record values, times, physical parameters or method settings. Pair southern history r with northern history s. Diagonal pairs (r=s) are the completed physical-system inputs; off-diagonal pairs (r!=s) join records from independent four-box copies that share the same prescribed forcing.

| Response retained | Same model realization | Different model realizations |
|---|---|---|
| A: forced + random | Shared forcing and actual within-system interaction | Shared forcing, no physical path between the selected records |
| D: random response | Interaction carried by stochastic responses | Independent stochastic-response control |

The selected off-diagonal endpoints cannot physically affect one another, although each endpoint still responds to interactions inside its own four-box copy. This is a control built from independent model copies, not a modification that deletes coupling terms from the original model. In A the records remain associated through common forcing; GC predictive evidence on such a pair is not automatically a computational error.

Use all six ordered off-diagonal pairs for each of A and D. The design has eighteen input pairs in total: six completed diagonal references and twelve new inputs, requiring 48 new calls to the original four methods. Retain the original estimator seed of each southern-history row so that diagonal reuse and row-wise contrasts preserve estimator randomization. No new process histories or chronology arrays are needed.

Display two 3-by-3 matrices per relevant method/outcome, distinguishing diagonal and off-diagonal cells. Keep both directional decisions, unresolved evidence and native scores. Compare each same-realization result with its row's two cross-realization controls. The cells are dependent recombinations of only three histories, not nine independent climate histories; do not calculate binomial uncertainty or claim general false-positive calibration from them.

## Outcomes that would guide the following step

| Possible result | Permitted interpretation and next action |
|---|---|
| A cross-realization evidence closely matches A diagonal evidence | The present statistic does not clearly distinguish the selected physical interaction from shared forced responses; chronology robustness of that statistic is not enough to validate the mechanism. |
| D diagonal differs clearly from D cross-realization | Random-response propagation carries observable interaction information; compare whether forcing or chronology obscures that distinction. |
| D diagonal and off-diagonal evidence are both weak and similar | Interaction is difficult to identify with the present endpoint information, resolution and method settings; prioritize a bounded ParCorr or full-state reference diagnostic. |
| D off-diagonal shows persistent evidence | Investigate finite-sample dependence, null construction and method assumptions; retain the result and avoid claiming calibrated error rates from the small dependent design. |
| A diagonal already differs from A off-diagonal | The method uses some information beyond the common forced waveform; a subsequent chronology comparison can test preservation of that distinction. |

A three-call ParCorr diagnostic on the existing D inputs remains inexpensive and useful, but it does not by itself distinguish common forcing from actual selected-record interaction. Further kernel tuning, a full-state/sampling matrix, and larger chronology ensembles should follow a specific unresolved explanation rather than be launched together.

The research target becomes: **does chronology uncertainty change the ability to distinguish actual interaction between observed records from common forced responses, and can it make an unsupported mechanism interpretation appear more consistent?** The original 04 remains a valid chronology-sensitivity comparison; this proposed control would provide a more explicit mechanism reference for interpreting it.
