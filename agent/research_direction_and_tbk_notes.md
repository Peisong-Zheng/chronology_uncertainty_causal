# Research direction and TBK model discussion

Discussion record initiated 10 November 2026; implementation status updated 11 November 2026. This file preserves the research-architecture discussion and distinguishes earlier proposals from accepted decisions. It does not replace the completed 02/03 experiments.

Current specification: [04 TBK chronology plan](04_tbk_chronology_plan.md). The user accepted graded uncertainty in both endpoint records, archival and near-verbatim model extraction, separate chronology preparation, an inspection notebook and a root experiment. The subsequent duration decision retains 30 kyr after 9 kyr startup, with 20-year observations and the common 29-1 kyr analysis grid. Model migration and complete native-output preflight have passed. The 123-input pilot started on 11 November 2026; its current completion state and scientific findings are recorded in [04 results](04_tbk_chronology_results.md). Earlier planning-only statements below refer to their historical discussion stage.

The central question accepted by the user is: **Can chronology uncertainty change our interpretation of climate mechanisms?** In the user's words: “它可以检验年代误差会不会改变我们对气候机制的解释”. For a system with reciprocal influences, one concrete concern is that uncertain dates could lead to an apparent one-way relationship and an unsupported claim of hemispheric dominance.

## Current decisions

| Item | Current status |
|---|---|
| Scientific framing | Focus on what chronology uncertainty does to causal interpretation, including directional ambiguity, mistaken direct-effect interpretations and stable but incorrect evidence. |
| Three-variable mistiming follow-up | Deferred by the user. The proposed simplified analytical explanation and further mistiming controls remain future sensitivity work. Existing 03/03_1/03_2 results are retained. |
| TBK physical structure | Keep the original four physical boxes for now. The proposed two-box replacement is not adopted. Do not add a temperature-to-circulation feedback equation. |
| Records of interest | Use the southernmost physical box, T1, and northernmost physical box, T4, as the initial candidate records. Keep the other states and forcing available as model diagnostics. |
| Immediate task | Complete the authorized 123-input both-records chronology-amplitude pilot, audit full outputs and inspect the resulting figures. |
| Current implementation | Original notebook archived; model and observation modules extracted; inspection notebook executed; exact migration, native-output and Windows-spawn preflight passed; production complete; full native-output and preservation audits passed. |

## How the completed experiments fit the research question

The completed experiments support two related mechanisms. In the bivariate branch, chronology perturbations can change whether a relationship is detected and whether its direction is resolved. In the common-driver branch, a measured driver with inaccurate dates may fail to remove shared information from the two target records. The latter is more specific than simply increasing the number of variables.

The [02_8 memory results](02_8_memory_sensitivity_results.md) show why memory belongs in the mechanism discussion. All 240 coupled chronology inputs remain significant in both directions in the original Tp=0 CCM stage. Much of the decline in complete relation recovery under short Y memory arises from unresolved lag evidence. High cross-map skill, retained forward evidence and successful exclusion of the reverse direction are different outcomes. The reconstruction optimum is not the imposed 20-year equation delay.

The [03_1 crossing results](03_1_process_chronology_cross_results.md) provide a concrete example of high chronology support for an incorrect direct-edge interpretation: GC reports Y-to-X in 9/10 members for Common driver only, realization R1 and Z scheme A2. This does not estimate a population false-positive rate. Members around one uncertain chronology center are not independent climate histories, and increasing their count does not correct a biased center or omitted physical variables.

The [03_2 type comparison](03_2_driver_error_models_results.md) is complete at 732/732 cases. It compares complete prescribed error settings: C/L/E/J have different realized overall errors. It does not isolate envelope shape at equal RMS. An equal-magnitude control is only needed if a future claim specifically attributes differences to shape independent of magnitude.

The proposed manuscript structure is: observation and chronology framework; loss of directional information with a physical-timescale explanation; degradation of common-driver adjustment; an optional TBK mechanism case; practical interpretation limits. Full parameter grids, sampling variants, engineering benchmarks and exhaustive type combinations can remain in supporting material or the repository. This is an organizational proposal, not a decision to discard results.

Future confirmation should prioritize independent simulation histories and independent uncertain chronology centers in a few selected comparisons. More members improve within-record Monte Carlo stability but do not supply independent histories. Method configurations must be specified before evaluating chronology effects; difficult accurate-age cases must not be selectively removed or tuned until successful.

## What the current TBK notebook represents

Source: [archived One_D_model_TBK.ipynb](../data/archives/tbk_original_20261110T214805311811Z/One_D_model_TBK.ipynb), particularly zero-based code cells 3, 8, 12 and 40. The notebook describes the model as a rough draft adapted from a diffusion/advection model. It contains repeated implementations and exploratory analyses rather than one designated production configuration.

There are four physical temperature boxes ordered from south to north, plus two ghost boxes imposing no-flux boundaries. The model step is one nominal year. An externally prescribed target u switches between -1 and +1, with approximately 1000 years off and 500 years on per complete 1500-year cycle. Actual transport forcing q approaches this target with a 10-year relaxation:

$$
q_{n+1}=q_n+\frac{\Delta t}{10}(u_n-q_n),\qquad \Delta t=1\ \mathrm{yr}.
$$

The current code updates the target at the end of each step; its indexing must be documented explicitly before assigning an exact physical lag. The temperature equations, with unit spatial grid spacing and zero restoring targets, are

$$
T_{i,n+1}=T_{i,n}+\Delta t\left[
D_i(T_{i+1,n}-T_{i,n})+D_{i-1}(T_{i-1,n}-T_{i,n})
+b_iq_n-\frac{\overline T_n}{100}
-\mathbb{1}_{i\in\{2,3\}}\frac{\overline T_{23,n}}{100000}
\right]+\eta_{i,n+1},
$$

$$
b=(0,-1,1/2,1/2)^\mathsf{T},\quad
\overline T=(T_1+T_2+T_3+T_4)/4,\quad
\overline T_{23}=(T_2+T_3)/2,
$$

with T0=T1 and T5=T4 at the boundaries. The implemented gain coefficients are zero. Temperature and noise amplitudes should retain the notebook's model units unless a physical calibration is supplied.

Transport removes an imposed amount from box 2 and adds half to each of boxes 3 and 4. This forcing sums to zero across the physical boxes. It does not depend on donor temperature, salinity or density, and temperature does not feed back to q. The model therefore simulates the response to prescribed transport changes; it does not generate AMOC or DO transitions through its own temperature/salinity feedback.

Temperature feedback already exists. Diffusion is reciprocal, and the global restoring term is

$$
-\overline T/100=-(T_1+T_2+T_3+T_4)/400.
$$

Consequently, the mathematical one-step temperature graph includes all four boxes, including reciprocal endpoint dependencies. Some of those links arise from the global restoring parameterization, so they must not all be described as literal local ocean transport. The two-endpoint records retain unobserved middle boxes and shared prescribed forcing. Their pairwise statistical arrows need not equal the full-state annual equation graph, particularly after coarser observation or interpolation.

Use the model to ask what mechanism interpretations the observations support. Do not assume that northern temperature must cause southern temperature alone, that the endpoint lead-lag is a physical causal delay, or that two significant arrows measure the strength of a physical feedback loop. The asymmetry q influences temperature while temperature does not determine q remains a useful model diagnostic, without making q a third analyzed record in the initial endpoint experiment.

## Verified origin and propagation of the existing noise

The notebook generates a Gaussian array with `np.random.randn`, scales it using `wnoise_sigma`, and interpolates it to the integration grid. In the current configuration both grids are annual, so this interpolation does not average neighboring annual noise draws. Across the four physical boxes and annual innovations, the input draws are independent:

$$
\eta_{i,n}\overset{\mathrm{iid}}{\sim}\mathcal N(0,0.005^2)
$$

for simulations 1 through 9. Simulation 10 sets the amplitude to zero. Multiple later analysis cells use `idx=9`, which is the tenth, noise-free trajectory, not an average of ten noisy realizations.

The addition `+ wnoiseint[:, :, t + 1]` is inside the temperature update. Noise is therefore **process noise**, entering the evolving physical state. Generating the whole random array before the loop does not make it post-processing noise. Every physical box receives a new perturbation each year; the experiment is not a single source perturbing one box only. The ghost-box random draws are overwritten by boundary copying and do not form two additional physical forcing sources.

For the present linear model, write

$$
\mathbf T_{n+1}=A\mathbf T_n+bq_n+\boldsymbol\eta_{n+1}.
$$

A perturbation inserted into box j at step n+1 changes later states, with the same subsequent forcing and innovations, by

$$
\delta\mathbf T_{n+1+k}=A^k\mathbf e_j\,\delta.
$$

It persists and spreads through diffusion and restoring terms. The global restoring term allows nonlocal effects at the next update; this is not propagation only between neighboring boxes. The resulting temperature variability becomes correlated even though the new noise inputs were independent. Noise in temperature does not propagate back into q, because that feedback is absent.

The code adds noise as a temperature increment, with no explicit sqrt(dt) scaling. Keep the annual integration step when preserving this implementation; changing the numerical step while holding 0.005 fixed would change the stochastic model. The noise is an imposed representation of unresolved process disturbances, not evidence that the deterministic model generates stochastic variability by itself.

Relevant source locations in the current notebook text: noise selection at lines 882-888; random draws at 948; noise-grid interpolation at 953; global restoring at 981; temperature update and noise addition at 984-990; boundary replacement at 1000-1001; example noise-free selection at 1852. Equivalent code occurs in the earlier and final core cells.

## Separate the three error layers

| Layer | Mathematical role | Does it change the physical trajectory? | Initial treatment proposed here |
|---|---|---|---|
| Process noise | eta in the model update | Yes; its effects persist and propagate | Keep the existing construction. Freeze each realization across all chronology comparisons. |
| Measurement or proxy-value error | y[k]=T(t[k])+epsilon[k] after integration | No, unless explicitly fed back into a different model | Set added value error to zero initially. Do not introduce a new noise sweep. |
| Chronology error | Assign the sampled value to an uncertain age | No; changes observed timing and subsequent resampling | Perturb only the two extracted records' age labels, using ordered, temporally correlated chronology curves. |

“Accurate-age data” may contain process noise. They are the true simulated realization at correct dates. The deterministic trajectory is a separate reference useful for illustrating forced response, not a definition of all accurate-age data. Post hoc value noise cannot reproduce the dynamical propagation of process disturbances.

## Earlier chronology proposal and current revision

The following records the earlier discussion. Its one-sided conditions in item 5 are superseded by the user's later decision: both endpoints are uncertain at every nonzero level, and amplitude is varied. The current proposed settings and file architecture are in [04_tbk_chronology_plan.md](04_tbk_chronology_plan.md). The four-box equations remain unchanged.

1. Generate and save a four-box realization, its forcing and innovations. Extract X=T1 and Y=T4 in forward physical-time order. Reuse that exact realization for every chronology condition. Use the noise-free run to explain the response shape, and noisy runs to assess the inference on stochastic trajectories; do not silently replace one with the other.
2. Fix the observation locations before chronology perturbation. First use a simple common observation schedule so as not to introduce asynchronous sampling as another experimental factor. Keep both an original accurate-grid reference and the same sampled records with accurate ages when interpolation is involved.
3. For each record, generate a smooth uncertain observed chronology around the known synthetic ages, using the existing observation-model logic. Keep sampled values fixed. Draw members around that observed chronology, rather than recentring every member on truth.
4. Reject whole curves that violate stratigraphic order. Every accepted member must cover the same fixed analysis window, with no sorting, clipping, moving values or extrapolation. Preserve the BP-age convention if the model's elapsed years are mapped to synthetic ages; choose an explicit reference origin rather than implying a historical date for the simulation.
5. Apply four paired conditions: accurate ages, southern record uncertain only, northern record uncertain only, both uncertain. Reuse the same southern curve in southern-only and both, and the same northern curve in northern-only and both. Initially use independent record chronologies, consistent with the existing main-study convention. Shared dating constraints would require a separate joint chronology model later.
6. Start with one envelope shape and a small number of nonzero amplitudes selected after reviewing actual response scales and sampling. C is a candidate for this first comparison, not a finalized choice. Ten chronology members remain the default pilot size. No amplitude, observation interval, history count or analysis window has been selected in this note.

The existing ordered chronology framework can be reused, but its old 30-kyr-specific settings cannot simply be copied into this 10-kyr notebook. Its correlation parameter is a fraction of normalized sample order, not a physical duration. Record the implied physical smoothness, proposed SD, realized RMS, relative endpoint timing errors, rejection/coverage counts and effective spread. Do not reuse the 29-1 kyr grid or existing age arrays without a new, explicit mapping. Sampling and age-model choices must remain fixed across the paired chronology conditions.

For samples at the same true age a[k], a schematic observation is

$$
x_k=T_1(t_k),\quad y_k=T_4(t_k),\qquad
\widehat a_{S,k}=a_k+e_{S,k},\quad
\widehat a_{N,k}=a_k+e_{N,k}.
$$

Changing e changes which temperatures are aligned during analysis, not the simulated temperatures. The difference e_N-e_S diagnoses relative chronology mismatch; it is distinct from the chronological uncertainty of either record alone. A 20-year plotting or analysis grid, if later chosen, would not create 20-year independent information from more sparsely observed data.

Do not inject chronology error into the model clock, the physical forcing schedule or a coupling term. Doing so would change the physical experiment. Do not draw a different process-noise realization for each age member, because this would mix chronology sensitivity with changing climate histories. Do not interpret a smooth age perturbation as independent white measurement error at each sample.

## Earlier evaluation discussion and choices subsequently resolved

The primary comparison is paired: inference for one stochastic physical trajectory with accurate ages versus inference for that same trajectory with uncertain ages. Retain both directional evidence, missing effects, unresolved states and native scores. Distinguish bidirectional retained evidence from an apparent one-way result, and never interpret a nonsignificant reverse test as proof that physical feedback is absent.

An accurate-age one-way result is already an observation/estimator limitation, not damage caused by chronology. Report this baseline explicitly rather than selecting histories that recover both directions. A future mechanism audit can perturb one physical state while keeping subsequent forcing and noise fixed, to document reciprocal responses without changing the model equations. Such an audit has not been run here and does not imply that estimator scores measure those intervention responses.

Before specifying a production script, resolve which existing diffusion variant is authoritative: cells 3/40 use the southern coefficient 0.0005, whereas cells 8/12 use 0.005. Also inspect the zero-state start and required spin-up, document the global restoring interpretation and the selected record/analysis resolution, and add explicit seeds and saved innovations when constructing a reproducible implementation. These were pending choices at that discussion stage; the accepted 04 plan subsequently selected cell 12 (0.005), 9 kyr startup and 20-year observations, with explicit process and chronology seeds. Do not infer a single authoritative physical version from notebook execution order or the last visible output.

## Literature context retained from the previous discussion

- [Runge et al. (2023), Causal inference for time series](https://www.nature.com/articles/s43017-023-00431-y): distinguishes causal discovery and effect estimation and discusses temporal resolution, hidden confounding and nonstationarity. This motivates an explicit observation-and-assumption framework.
- [Correia et al. (2026), Best practices for moving from correlation to causation in ecological research](https://www.nature.com/articles/s41467-026-69878-z): connects causal tasks, frameworks, methods and mechanistic knowledge. Agreement or predictive fit alone does not validate a physical interpretation.
- [Suzuki et al. (2026), Advances in causal discovery methods for ecological time series](https://pmc.ncbi.nlm.nih.gov/articles/PMC13542360/): discusses GC/CCM and statistical approaches, and the need to define simulation truth and evaluation targets carefully. This supports separating direct edges, transmitted effects and reconstruction information.
- [Lohmann et al. (2026), Quantitative comparison of causal inference methods for climate tipping points](https://npg.copernicus.org/articles/33/313/2026/): already combines multivariate method comparisons with an AMOC/sea-ice application. Our proposed distinction is chronology-specific degradation of interpretation, not simply adding a physical case to an algorithm comparison.
- [Saberski et al. (2024), The impact of data resolution on dynamic causal inference in multiscale ecological networks](https://www.nature.com/articles/s42003-024-07054-z): shows scale dependence of reconstructed dynamical relationships, motivating care when comparing model-equation and observed-record graphs.

These sources provide context; they do not establish an exhaustive novelty claim for the proposed research. A possible manuscript framing is: “We test whether uncertainty in the timing of climate records changes the inferred direction of influence and the resulting interpretation of known climate mechanisms.” Claims about observed effects will require the corresponding experiment.

## Preservation record

This discussion changed documentation only. At inspection, the source notebook SHA256 was `7F8E2A5CF055EA56FC5D7F90A8F76D88C6BCC20CDEA04A3A182331E253FD51EE`, length 5024031 bytes. Existing completed results and the notebook remain unchanged. Update this note when the user selects a model version or approves an experiment; keep proposed choices distinct from accepted decisions.

## Accepted implementation update: 10 November 2026

The user accepted the TBK proposal and asked whether a longer record could align with previous experiments. The implementation retains 30 kyr after 9 kyr startup, uses annual model integration and 20-year synchronous observations, and analyzes 29-1 kyr (1401 points). Alignment is useful for sample-count and CCM-library comparability, but does not make the physical system identical to the synthetic AR systems. Extend only the original forcing event table to cover 39 kyr. The original 10 kyr default remains available for exact migration tests. Both endpoints receive independent chronology errors at 20/50/100/200 yr proposal SD, with 1000 yr error correlation length (1/30 of support), three histories and ten members per level. Full specification and acceptance checks are in `04_tbk_chronology_plan.md`; earlier shorter-record settings are superseded.
