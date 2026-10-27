"""Inspect accurate-age baselines and chronology-error ensembles.

Run with python 02_check_causal_pipeline.py in the project environment.
This is a functional inspection of all sixteen ordered chronology-type pairs.
Figures and their numerical inputs are saved under this script's name.
Both baselines and every chronology member use the same 20-year analysis
grid. Ensemble support describes this record, not a calibrated accuracy.
"""

from pathlib import Path
import json
import time
import platform
from importlib.metadata import version
from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from toolbox.process import generate_process
from toolbox.observation import (sample_irregular_indices, simulate_observations,
                                 make_age_uncertainty, draw_age_ensemble)
from toolbox.causal import resample_pair, run_granger, run_pcmciplus, run_lkif, run_ccm, prepare_gpdc_null
from toolbox.evaluation import summarize_ensemble, evaluate_replicates


# Simulation settings. All physical times use kyr; BP zero is 1950 CE.
SEED = 20261009                 # Base seed, separate from the parameter-pilot seeds.
DURATION = 30.0                # Retained record length: 30 kyr.
DT = 0.020                     # Native generator step: 20 years; 1501 retained points.
TRUE_LAG = 0.020               # X at the previous step drives Y now: one native step, 20 years.
MEMORY_X = 0.5                 # Intrinsic X decay time: 500 years.
MEMORY_Y = 0.2                 # Intrinsic Y decay time: 200 years.
BACKGROUND_X = 1 / np.sqrt(1 - 0.95**2)  # Stationary driver SD, in arbitrary units.
BACKGROUND_Y = 1 / np.sqrt(1 - 0.8**2)   # Y SD without coupling, in arbitrary units.
COUPLING_GAIN = 4.0             # Long-run Y response to a unit held g(X); checked for the one-step benchmark.
# K is a model gain, not an estimator's score; reverse detections remain in all outputs.
BURN_TIMESCALES = 5             # Warm-up duration in units of the longest memory.
YOUNG_SPACING = 0.0275          # Pre-rounding target gap, kyr; this is not the actual mean spacing.
OLD_SPACING = 0.0275            # Same target at both ends: no age-dependent thinning in the main control.
SPACING_JITTER = 0.10           # Rounded gaps are 20/40 years; about 95.7% of native points are retained.
# The earlier 40-to-60-year targets with jitter=0.25 remain a sampling sensitivity scenario.
MEASUREMENT_NOISE_RATIO = 0.0   # Measurement SD / full-process SD; zero in the main design.

# Analysis settings. Both baselines use the same window and native step.
ANALYSIS_OLD = 29.0            # Old end of the shared analysis window, kyr BP.
ANALYSIS_YOUNG = 1.0           # Young end; original data outside the window bracket interpolation.
GRID_STEP = DT                 # Match the generator step; no separate downsampling baseline.
MAX_LAG = 5                   # PCMCI+/LKIF history up to 100 years; wider than the true 20-year lag.
GC_ORDER = 1                  # Fixed VAR order in both baselines; matches the known one-step process.
ALPHA = 0.05                  # Per-direction significance level; lag families use Holm.
N_REPLICATES = 1              # Independent process realizations; 1 is inspection, not an error-rate study.
GPDC_NULL_SAMPLES = 10000      # Reference-null draws; fine p resolution before the 6-lag Holm correction.
GPDC_OPTIMIZE_KERNEL = False   # False fixes the package's RBF + white-noise kernel at its initial values.

# CCM uses the frozen notebook configuration; the original test remains at Tp=0.
CCM_E = 3                     # Embedding coordinates; chosen in the separate accurate-age pilot.
CCM_TAU = -1                  # Backward embedding interval: 1 step = 20 years; history spans 40 years.
CCM_EXCLUSION_RADIUS = 10      # Exclude temporally adjacent neighbors within 200 years.
CCM_LIBRARY_FRACTION = 0.80    # Fraction of usable embedded rows in the single fixed library.
CCM_SAMPLES = 10               # Random libraries averaged at that one library size.
CCM_SURROGATES = 100           # Random cut-and-swap records for each direction; p resolution 1/101.
CCM_LAG_STEPS = list(range(-15, 16))  # Target shifts from -300 to +300 years, including zero.
CCM_LAG_LIBRARIES = 100        # Matched random libraries for the observed lag curves, not new null records.
CCM_LAG_SIGN_THRESHOLD = 0.90  # Required stable peak-sign fraction; an engineering rule, not confidence.

# Chronology settings. Scales are per-record proposal SDs, not hard error bounds.
RUN_BASELINES = True           # Build exact-age baselines on a fresh Windows checkout.
RESUME = True                 # Keep complete members from this same configuration; never discard failed fits.
PREPARE_ONLY = False          # True stops after baseline fits and chronology inputs, for worker benchmarking.
RUN_AGE_ENSEMBLES = True       # Run the restored age-error comparison after the exact-age checks.
N_AGE_MEMBERS = 10             # Generate ten ordered chronology realizations for each record type.
N_CAUSAL_MEMBERS = 10          # Analyze all ten; complete-set support changes in steps of 0.1.
N_WORKERS = 32                # Fastest tested on this Windows PC; one numerical thread per worker.
MAX_AGE_DRAWS = 10000          # Whole-curve proposals; missing members remain explicit failures.
AGE_CORRELATION_LENGTH = 0.05 * 10 / DURATION  # Normalized sample-order kernel length, about 500 years.
AGE_EXP_SHAPE = 2.0           # Fixed E-envelope curvature; each type's peak SD is set in AGE_SCALES.
UNCERTAINTY_LEVEL = 2         # Select one of the four planned uncertainty levels for every type pair.
AGE_SCALE_LEVELS = [0.010, 0.020, 0.040, 0.080]  # C/L/E peak SDs in kyr, for levels 1–4.
JUMP_COUNT_LEVELS = [1, 3, 5, 7]  # J increases the number of local error peaks, not their amplitude.
JUMP_PEAK_SD = 0.010          # Fixed J peak SD: 10 years at every level.
JUMP_INTERVAL_WIDTH = 2.8     # Width of each compact Gaussian-like J interval, in kyr.
JUMP_CENTER_FRACTIONS = [0.5, 0.2, 0.8, 0.35, 0.65, 0.05, 0.95]  # Nested centers within the analysis window.
SUPPORT_THRESHOLD = 0.90       # Reporting threshold, not the probability that a direction is correct.
MIN_VALID_FRACTION = 0.90      # Require at least 90% of requested member pairs for a summary decision.
ACTIVE_AGE_SD_TOLERANCE = 1e-12  # kyr; exclude numerical-zero J tails only from age diagnostics.
if UNCERTAINTY_LEVEL not in [1, 2, 3, 4]:
    raise ValueError("UNCERTAINTY_LEVEL must be 1, 2, 3, or 4.")
AGE_KINDS = ["C", "L", "E", "J"]
AGE_SCALES = {kind: AGE_SCALE_LEVELS[UNCERTAINTY_LEVEL - 1] for kind in ["C", "L", "E"]}
AGE_SCALES["J"] = JUMP_PEAK_SD
JUMP_COUNT = JUMP_COUNT_LEVELS[UNCERTAINTY_LEVEL - 1]
JUMP_COVERAGE = JUMP_COUNT * JUMP_INTERVAL_WIDTH / (ANALYSIS_OLD - ANALYSIS_YOUNG)
# Each tuple names the X type and Y type; C-J and J-C are distinct experiments.
AGE_SCENARIOS = [(f"{kind_x}-{kind_y}", kind_x, kind_y)
                 for kind_x in AGE_KINDS for kind_y in AGE_KINDS]

ROOT = Path(__file__).resolve().parent
STEM = Path(__file__).stem
DATA_DIR = ROOT / "data" / "processed" / STEM
FIGURE_DIR = ROOT / "figures" / STEM
METHODS = ["PCMCI+", "LKIF", "CCM", "GC"]
CACHE_DIR = ROOT / "data" / "cache" / "gpdc"


def execution_environment():
    """Keep resumed calculations on the same numerical platform."""
    blas = np.__config__.CONFIG.get("Build Dependencies", {}).get("blas", {})
    return {"system": platform.system(), "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy_blas": {name: blas.get(name) for name in
                           ["name", "version", "openblas configuration"]}}


def check_resume_parameters(path, parameters):
    """Resume only this exact experiment; changing worker count is allowed."""
    if not path.exists():
        return
    previous = json.loads(path.read_text())
    ignored = {"RESUME", "N_WORKERS", "RUN_BASELINES", "RUN_AGE_ENSEMBLES", "RUN_CAUSAL", "PREPARE_ONLY"}
    names = [name for name in parameters if name.isupper() and name not in ignored]
    names = names + ["execution_environment", "analysis_revision", "packages"]
    changed = [name for name in names if previous.get(name) != parameters.get(name)]
    if changed:
        raise ValueError(f"Saved run differs in {changed}. Preserve it in another folder before starting a new experiment.")


def read_completed_member(job):
    """Reuse complete two-direction results, including recorded method failures."""
    path = DATA_DIR / f"{job['label']}_edges.csv"
    if not RESUME or not path.exists():
        return None
    edges = pd.read_csv(path, dtype={"edge_detected": "boolean", "unoriented": "boolean"})
    keys = ["method", "source", "target"]
    expected = {(method, source, target) for method in METHODS
                for source, target in [("x", "y"), ("y", "x")]}
    if not set(keys + ["status", "input_status", "edge_detected", "score", "p_value", "unoriented"]).issubset(edges):
        return None
    if len(edges) != len(expected) or edges.duplicated(keys).any():
        return None
    if set(edges[keys].itertuples(index=False, name=None)) != expected or edges["status"].isna().any():
        return None
    for name in ["replicate", "process", "scenario", "member"]:
        if name not in edges or not edges[name].eq(job[name]).all():
            raise ValueError(f"Saved member identity differs: {path.name}")
    if not edges["input_status"].eq(job["status"]).all():
        raise ValueError(f"Saved input preparation status differs: {path.name}")
    return edges


def save_input(series, path):
    """Check regenerated inputs before reusing any fitted member."""
    if RESUME and path.exists():
        previous = pd.read_csv(path)
        if (list(previous.columns) != list(series.columns) or previous.shape != series.shape
                or not np.allclose(previous.to_numpy(), series.to_numpy(), rtol=0, atol=1e-12)):
            raise RuntimeError(f"Input changed for {path.name}; preserve the old run before recomputing.")
        return
    series.to_csv(path, index=False)


def gpdc_null_path(seed):
    n_points = int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / GRID_STEP)) + 1
    sample_size = n_points - 2 * MAX_LAG
    return CACHE_DIR / f"null_T{sample_size}_n{GPDC_NULL_SAMPLES}_seed{seed}.npz"


def prepare_gpdc_caches():
    """Prepare native reference nulls once in the parent, before worker creation."""
    sample_size = int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / GRID_STEP)) + 1 - 2 * MAX_LAG
    with threadpool_limits(limits=1):
        for replicate in range(N_REPLICATES):
            seed = SEED + 1000 * replicate + 300
            prepare_gpdc_null(gpdc_null_path(seed), sample_sizes=[sample_size],
                              sig_samples=GPDC_NULL_SAMPLES, seed=seed)


def analyze_pair(series, *, seed, label, diagnostics_name="method_diagnostics.jsonl"):
    """Call each existing estimator and retain failures explicitly."""
    edge_tables = []
    details = {}
    effective_length = len(series) - (CCM_E - 1) * abs(CCM_TAU)
    library_size = int(CCM_LIBRARY_FRACTION * effective_length)
    for method in METHODS:
        started = time.perf_counter()
        try:
            if method == "PCMCI+":
                # This is a native GPDC option; the default optimizer remains available.
                gp_params = None if GPDC_OPTIMIZE_KERNEL else {"optimizer": None}
                result = run_pcmciplus(series, tau_max=MAX_LAG, pc_alpha=ALPHA,
                                      alpha=ALPHA, sig_samples=GPDC_NULL_SAMPLES, seed=seed,
                                      gp_params=gp_params,
                                      null_dist_filename=str(gpdc_null_path(seed)))
            elif method == "LKIF":
                result = run_lkif(series, max_lag=MAX_LAG, alpha=ALPHA)
            elif method == "CCM":
                result = run_ccm(series, E=CCM_E, tau=CCM_TAU,
                                 exclusion_radius=CCM_EXCLUSION_RADIUS, lib_size=library_size,
                                 sample=CCM_SAMPLES, n_surrogates=CCM_SURROGATES,
                                 seed=seed, alpha=ALPHA, lag_steps=CCM_LAG_STEPS,
                                 n_libraries=CCM_LAG_LIBRARIES,
                                 sign_threshold=CCM_LAG_SIGN_THRESHOLD)
            else:
                result = run_granger(series, max_lag=GC_ORDER, alpha=ALPHA)
            edges = result["edges"].copy()
            edges["method"] = method
            diagnostics = result["diagnostics"]
            with (DATA_DIR / diagnostics_name).open("a", encoding="utf-8") as handle:
                handle.write(json.dumps({"label": label, "method": method, **diagnostics}) + "\n")
            if "unoriented" not in edges:
                edges["unoriented"] = False
            if (diagnostics.get("undirected_contemporaneous", False)
                    or diagnostics.get("conflicting_contemporaneous", False)):
                edges["unoriented"] = True
            details[method] = result
            # Save native numerical outputs, including lag arrays, without evaluating lag accuracy.
            arrays = {}
            for key, value in result["raw"].items():
                if isinstance(value, np.ndarray):
                    arrays[key] = value.astype(str) if value.dtype == object else value
                elif isinstance(value, pd.DataFrame):
                    value.to_csv(DATA_DIR / f"{label}_{method.replace('+', 'plus')}_{key}.csv", index=False)
            if arrays:
                np.savez_compressed(DATA_DIR / f"{label}_{method.replace('+', 'plus')}_raw.npz", **arrays)
        except Exception as error:
            edges = pd.DataFrame({
                "method": [method, method], "source": ["x", "y"], "target": ["y", "x"],
                "score": [np.nan, np.nan], "score_name": ["unavailable", "unavailable"],
                "edge_detected": [None, None], "p_value": [np.nan, np.nan],
                "p_kind": ["unavailable", "unavailable"], "status": ["error", "error"],
                "unoriented": [False, False], "error": [str(error), str(error)],
            })
        elapsed = time.perf_counter() - started
        edges["seconds"] = elapsed
        edge_tables.append(edges)
        lag_message = f", lag={edges['lag_status'].tolist()}" if "lag_status" in edges else ""
        print(f"{label}: {method}, {elapsed:.1f} s, {edges['status'].tolist()}{lag_message}", flush=True)
    return pd.concat(edge_tables, ignore_index=True), details


def run_baselines():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    parameters = {name: value for name, value in globals().items()
                  if name.isupper() and isinstance(value, (str, int, float, list, dict))}
    parameters = json.loads(json.dumps(parameters))
    parameters["packages"] = {name: version(name) for name in
                              ["numpy", "scipy", "pandas", "scikit-learn", "tigramite", "pyEDM", "LK-Info-Flow", "statsmodels", "dcor"]}
    parameters["execution_environment"] = execution_environment()
    parameters["analysis_revision"] = "native_null_cache_and_ccm_object_reuse_v1"
    parameter_path = DATA_DIR / "baseline_parameters.json"
    # A completed baseline may be reused when only chronology settings change.
    if RESUME and parameter_path.exists():
        previous = json.loads(parameter_path.read_text())
        check_saved_baselines(require_complete=False)
        if previous.get("complete", False):
            check_saved_baselines()
            print("Reusing both complete accurate-age baselines.", flush=True)
            return
    elif parameter_path.exists():
        raise ValueError("Baseline outputs already exist. Set RESUME=True or preserve them in another folder.")
    parameters["complete"] = False
    parameter_path.write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    (DATA_DIR / "baseline_method_diagnostics.jsonl").write_text("", encoding="utf-8")
    if N_REPLICATES < 1 or DT <= 0:
        raise ValueError("Require a positive replicate count and generator step.")
    if GRID_STEP != DT:
        raise ValueError("Both baselines must use the native generator step.")
    if MEASUREMENT_NOISE_RATIO != 0:
        raise ValueError("The paired sampling comparison requires zero measurement noise.")
    if not np.isclose(DURATION / DT, round(DURATION / DT), rtol=0, atol=1e-8):
        raise ValueError("DURATION / DT must be an integer.")
    grid_intervals = (ANALYSIS_OLD - ANALYSIS_YOUNG) / DT
    if grid_intervals <= 0 or not np.isclose(grid_intervals, round(grid_intervals), rtol=0, atol=1e-8):
        raise ValueError("The analysis span must be a positive integer number of native steps.")

    n_steps = int(round(DURATION / DT)) + 1
    age = DURATION - np.arange(n_steps) * DT
    age_grid = np.linspace(ANALYSIS_OLD, ANALYSIS_YOUNG, int(round(grid_intervals)) + 1)
    # Crop the original process to the common window without skipping any points.
    analysis_indices = np.rint((DURATION - age_grid) / DT).astype(int)
    if np.any(analysis_indices < 0) or np.any(analysis_indices >= n_steps):
        raise ValueError("The analysis window must lie within the generated process.")
    if not np.allclose(age[analysis_indices], age_grid, rtol=0, atol=1e-10):
        raise ValueError("The analysis ages must coincide with original generator ages.")
    if not np.all(np.diff(analysis_indices) == 1):
        raise ValueError("The original baseline must retain every native point in the window.")

    # Recalculate discrete coefficients while retaining memory, background SD and long-run gain.
    # This defines a new 20-year discrete model, not a subsample of the previous 2-year model.
    ar_x = np.exp(-DT / MEMORY_X)
    ar_y = np.exp(-DT / MEMORY_Y)
    noise = (BACKGROUND_X * np.sqrt(1 - ar_x**2), BACKGROUND_Y * np.sqrt(1 - ar_y**2))
    burn_in = int(np.ceil(BURN_TIMESCALES * max(MEMORY_X, MEMORY_Y) / DT))

    uniform_tables = []
    irregular_tables = []
    sampling_diagnostics = []
    sampling_examples = {}
    for replicate in range(N_REPLICATES):
        # Reuse the same innovations and sample locations across the three process types.
        replicate_seed = SEED + 1000 * replicate
        indices = {}
        for record, offset in [("x", 1), ("y", 2)]:
            indices[record] = sample_irregular_indices(
                age, young_spacing=YOUNG_SPACING, old_spacing=OLD_SPACING,
                jitter=SPACING_JITTER, rng=replicate_seed + offset,
            )
        processes = {}
        for process_name in ["linear", "nonlinear", "independent"]:
            gain = 0.0 if process_name == "independent" else COUPLING_GAIN
            form = "nonlinear" if process_name == "nonlinear" else "linear"
            processes[process_name] = generate_process(
                n_steps=n_steps, dt=DT, age_start=DURATION, lag=TRUE_LAG,
                ar_x=ar_x, ar_y=ar_y, coupling=(1 - ar_y) * gain,
                coupling_form=form, noise_std=noise, burn_in=burn_in, rng=replicate_seed,
            )

        if replicate == 0:
            fig, _ = plot_process_truth(processes)
            fig.savefig(FIGURE_DIR / "00_model_truth.png", dpi=180)
            fig.savefig(FIGURE_DIR / "00_model_truth.pdf")
            plt.close(fig)

        for process_name, process in processes.items():
            process["series"].to_csv(DATA_DIR / f"r{replicate}_{process_name}_generated_series.csv", index=False)
            uniform = process["series"].iloc[analysis_indices].copy().reset_index(drop=True)
            uniform["time"] = age_grid[0] - age_grid
            uniform["age"] = age_grid
            label = f"r{replicate}_{process_name}_original_uniform"
            uniform.to_csv(DATA_DIR / f"{label}_series.csv", index=False)
            uniform_result, _ = analyze_pair(
                uniform, seed=replicate_seed + 300, label=label,
                diagnostics_name="baseline_method_diagnostics.jsonl",
            )
            uniform_result["replicate"] = replicate
            uniform_result["process"] = process_name
            uniform_tables.append(uniform_result)
            pd.concat(uniform_tables, ignore_index=True).to_csv(DATA_DIR / "baseline_original_uniform_results.csv", index=False)

            # Sample the same original process directly, with exact ages and zero measurement error.
            # Save the irregular observations before interpolation so the information loss is visible.
            observed = {}
            for record, offset in [("x", 31), ("y", 32)]:
                observed[record] = simulate_observations(
                    age, process["series"][record], indices[record], np.zeros(len(indices[record])),
                    measurement_noise_ratio=MEASUREMENT_NOISE_RATIO,
                    rng_age=replicate_seed + offset, rng_measurement=replicate_seed + offset,
                )["observed"]
                observed[record].to_csv(DATA_DIR / f"r{replicate}_{process_name}_{record}_irregular_observed.csv", index=False)
            irregular = resample_pair(
                observed["x"]["age_obs"], observed["x"]["value_obs"],
                observed["y"]["age_obs"], observed["y"]["value_obs"], age_grid=age_grid,
            )["series"]
            label = f"r{replicate}_{process_name}_irregular_accurate"
            irregular.to_csv(DATA_DIR / f"{label}_series.csv", index=False)
            irregular_result, _ = analyze_pair(
                irregular, seed=replicate_seed + 300, label=label,
                diagnostics_name="baseline_method_diagnostics.jsonl",
            )
            irregular_result["replicate"] = replicate
            irregular_result["process"] = process_name
            irregular_tables.append(irregular_result)
            pd.concat(irregular_tables, ignore_index=True).to_csv(DATA_DIR / "baseline_irregular_accurate_results.csv", index=False)

            for record in ["x", "y"]:
                direct = uniform[record].to_numpy()
                interpolated = irregular[record].to_numpy()
                sample_age = observed[record]["age_obs"].to_numpy()
                gaps = -np.diff(sample_age)
                in_window = (sample_age >= ANALYSIS_YOUNG) & (sample_age <= ANALYSIS_OLD)
                sampling_diagnostics.append({
                    "replicate": replicate, "process": process_name, "record": record,
                    "n_generated": n_steps, "n_analysis_grid": len(uniform),
                    "n_irregular_full_record": len(sample_age), "n_irregular_in_window": int(in_window.sum()),
                    "fraction_native_points_retained": len(sample_age) / n_steps,
                    "fraction_grid_points_interpolated": 1 - int(in_window.sum()) / len(uniform),
                    "fraction_gaps_above_native_step": float(np.mean(gaps > DT * (1 + 1e-8))),
                    "minimum_gap_kyr": gaps.min(), "mean_gap_kyr": gaps.mean(), "maximum_gap_kyr": gaps.max(),
                    "rmse_over_original_sd": np.sqrt(np.mean((interpolated - direct)**2)) / np.std(direct),
                    "correlation": np.corrcoef(direct, interpolated)[0, 1],
                    "sd_ratio": np.std(interpolated) / np.std(direct),
                    "increment_variance_ratio": np.var(np.diff(interpolated)) / np.var(np.diff(direct)),
                    "lag1_original": np.corrcoef(direct[:-1], direct[1:])[0, 1],
                    "lag1_interpolated": np.corrcoef(interpolated[:-1], interpolated[1:])[0, 1],
                })
            if replicate == 0:
                sampling_examples[process_name] = {"uniform": uniform, "irregular": irregular, "observed": observed}

    uniform_results = pd.concat(uniform_tables, ignore_index=True)
    irregular_results = pd.concat(irregular_tables, ignore_index=True)
    pd.DataFrame(sampling_diagnostics).to_csv(DATA_DIR / "sampling_diagnostics.csv", index=False)
    keys = ["replicate", "process", "method", "source", "target"]
    columns = keys + ["score", "p_value", "edge_detected", "status"]
    comparison = uniform_results[columns].merge(irregular_results[columns], on=keys,
                                               suffixes=("_uniform", "_irregular_accurate"), validate="one_to_one")
    comparison["score_change"] = comparison["score_irregular_accurate"] - comparison["score_uniform"]
    comparison.to_csv(DATA_DIR / "sampling_method_comparison.csv", index=False)
    ccm_baselines = pd.concat([
        uniform_results[uniform_results["method"] == "CCM"].assign(baseline="original_uniform"),
        irregular_results[irregular_results["method"] == "CCM"].assign(baseline="irregular_accurate"),
    ], ignore_index=True)
    ccm_columns = ["replicate", "process", "baseline", "source", "target", "score", "p_value",
                   "raw_edge_detected", "edge_detected", "unoriented", "status", "lag_status",
                   "lag_reason", "optimal_lag_steps", "lag_past_fraction", "lag_future_fraction"]
    ccm_baselines.reindex(columns=ccm_columns).to_csv(DATA_DIR / "ccm_baseline_filter_comparison.csv", index=False)

    parameters["complete"] = True
    parameters["LKIF_author_commit"] = "95c088ef8cc5e5dab2c6279158e98338ee51124a"
    parameters["process_truth"] = {name: process["truth"] for name, process in processes.items()}
    parameters["n_analysis_points"] = len(age_grid)
    parameters["interpretation"] = "One realization per process; no ensemble, detection-rate estimate or common causal-strength scale."
    parameters["original_uniform_baseline"] = "Native 20-year values cropped to the common window; no decimation or interpolation."
    parameters["irregular_accurate_baseline"] = "Irregular samples directly from the same process; exact ages; linear interpolation to the same 20-year grid."
    parameters["sampling_control"] = "Gentle main control: mostly 20-year gaps with occasional 40-year gaps; report actual retention, not the pre-rounding target."
    parameters["GC_order_rule"] = "Fixed GC_ORDER for both baselines; no data-driven selection."
    parameters["CCM_lag_rule"] = "Tp=0 rho and swap p retained; both directions independently screened for stable interior negative peaks. Lag seed equals the CCM seed."
    parameters["sampling_comparison"] = "Combined effect of fewer observations, irregular gaps and interpolation; no age or measurement error."
    (DATA_DIR / "baseline_parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")

    plot_checks(uniform_results, irregular_results, sampling_examples)
    print(comparison[["process", "method", "source", "target", "p_value_uniform", "p_value_irregular_accurate",
                      "edge_detected_uniform", "edge_detected_irregular_accurate"]].to_string(index=False))


def plot_process_truth(processes):
    """Plot true cross-variable edges and the mean response to a held driver.

    Read parameters from generate_process's truth, without fitting the data.
    The curves describe equilibrium means after holding X constant, not
    correlations or the realized noisy Y series. K is a model gain, not a
    target value shared by the four causal estimators.
    """
    names = ["linear", "nonlinear", "independent"]
    z = np.linspace(-3, 3, 301)  # Held X in units of its theoretical stationary SD.
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), height_ratios=[1, 1.3])
    response_limit = 0.0
    for column, name in enumerate(names):
        truth = processes[name]["truth"]
        gain = truth["coupling_gain"]
        coupling = truth["coupling"]
        ax = axes[0, column]
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.set_title(name.capitalize(), fontsize=13, fontweight="bold")
        for x_position, label, color in [(0.17, "X", "tab:blue"), (0.83, "Y", "tab:orange")]:
            ax.text(x_position, 0.60, label, ha="center", va="center", fontsize=16,
                    bbox={"boxstyle": "circle,pad=0.5", "fc": "white", "ec": color, "lw": 2})
        if ["x", "y"] in truth["edges"]:
            ax.annotate("", xy=(0.73, 0.60), xytext=(0.27, 0.60),
                        arrowprops={"arrowstyle": "-|>", "lw": 2, "color": "0.2"})
            ax.text(0.5, 0.82, f"K = {gain:g}", ha="center", fontsize=12)
            ax.text(0.5, 0.35, f"Delay = {1000 * truth['lag']:g} yr", ha="center")
            ax.text(0.5, 0.17, "No Y → X link", ha="center", color="0.35")
        else:
            ax.text(0.5, 0.82, "K = 0", ha="center", fontsize=12)
            ax.text(0.5, 0.35, "No link in either direction", ha="center")
        ax.text(0.5, 0.01, f"c = {coupling:.5f} per {1000 * truth['dt']:g}-yr step", ha="center", fontsize=9)

        driver = z
        formula = r"$g(X) = z$"
        if truth["coupling_form"] == "nonlinear":
            scale = truth["coupling_transform_std"]
            driver = np.tanh(2 * z) / scale
            formula = rf"$g(X) = \tanh(2z) / {scale:.3f}$"
        response = gain * driver
        response_limit = max(response_limit, float(np.max(np.abs(response))))
        ax = axes[1, column]
        ax.axhline(0, color="0.8", lw=0.8)
        ax.axvline(0, color="0.8", lw=0.8)
        ax.plot(z, response, color="tab:blue" if gain != 0 else "0.4", lw=2)
        ax.set_title(formula if gain != 0 else "No forced response", fontsize=11)
        ax.set_xlabel(r"Held driver $z = X / \mathrm{SD}(X)$")
        ax.set_xlim(z[0], z[-1])
        if column == 0:
            ax.set_ylabel("Equilibrium mean of Y (a.u.)")
    for ax in axes[1]:
        ax.set_ylim(-1.1 * max(response_limit, 1), 1.1 * max(response_limit, 1))

    fig.suptitle("Known model relationships and held-driver responses", fontsize=15)
    fig.text(0.5, 0.075,
             r"$Y_t = a_Y Y_{t-1} + c\,g(X_{t-L}) + v_t$; "
             r"$K = c/(1-a_Y)$; held-driver mean $\overline{Y} = K\,g(X)$.",
             ha="center", fontsize=11)
    fig.text(0.5, 0.025,
             "Independent, zero-mean process noise is present. Model K and causal-method scores have different meanings.",
             ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.12, 1, 0.95))
    return fig, axes


def plot_baseline_results(baselines, reference, title):
    """Use the same score limits for the two paired baseline figures."""
    fig, axes = plt.subplots(1, 4, figsize=(16, 5.2))
    processes_order = ["linear", "nonlinear", "independent"]
    for ax, method in zip(axes, METHODS):
        p_labels = {}
        for direction, marker, color in [("x", "o", "tab:blue"), ("y", "s", "tab:orange")]:
            subset = baselines[(baselines["method"] == method) & (baselines["source"] == direction) & (baselines["replicate"] == 0)]
            for _, result in subset.iterrows():
                x = processes_order.index(result["process"]) + (-0.1 if direction == "x" else 0.1)
                if result["status"] == "ok" and np.isfinite(result["score"]):
                    detected = result["edge_detected"] == True
                    if method == "CCM" and result.get("lag_status") == "failed":
                        ax.scatter(x, result["score"], marker="X", color="0.4", s=55)
                    elif result.get("unoriented", False) == True:
                        ax.scatter(x, result["score"], marker="D", edgecolor=color,
                                   facecolor="#7FC97F", s=45)
                    else:
                        ax.scatter(x, result["score"], marker=marker, edgecolor=color,
                                   facecolor=color if detected else "none", s=45)
                    p_labels[(result["process"], direction)] = f"{result['p_value']:.3g}"
                else:
                    ax.text(x, 0.15, "Failed", transform=ax.get_xaxis_transform(), ha="center", fontsize=7, rotation=90)
            ax.scatter([], [], marker=marker, color=color, label="X→Y" if direction == "x" else "Y→X")
        labels = []
        for process_name in processes_order:
            forward_p = p_labels.get((process_name, "x"), "NA")
            reverse_p = p_labels.get((process_name, "y"), "NA")
            p_name = "p0" if method == "CCM" else "p"
            labels.append(f"{process_name}\nX→Y {p_name}={forward_p}\nY→X {p_name}={reverse_p}")
        ax.set_xticks(range(3), labels, fontsize=7)
        ax.set_xlim(-0.65, 2.65)
        values = reference.loc[reference["method"] == method, "score"].to_numpy()
        values = values[np.isfinite(values)]
        if len(values):
            span = max(float(np.ptp(values)), 0.01)
            ax.set_ylim(float(np.min(values)) - 0.16 * span, float(np.max(values)) + 0.16 * span)
        ax.set_title(method)
        ax.set_ylabel("Native score")
    axes[0].legend(fontsize=8)
    fig.legend(handles=[
        Line2D([], [], marker="o", color="black", ls="none", label="Retained direction"),
        Line2D([], [], marker="o", color="black", markerfacecolor="none", ls="none", label="Direction not retained"),
        Line2D([], [], marker="D", color="black", markerfacecolor="#7FC97F", ls="none", label="Unresolved evidence"),
        Line2D([], [], marker="X", color="0.4", ls="none", label="CCM lag failed"),
    ], loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.035), frameon=False, fontsize=8)
    fig.text(0.5, 0.015, "Hollow markers include CCM lag-screen rejection. CCM score and p0 remain the original Tp=0 results.",
             ha="center", fontsize=8)
    fig.suptitle(title + f" ({1000 * GRID_STEP:g}-yr grid)"
                 + "\nReplicate 0; scores differ by method")
    return fig


def plot_ccm_baseline_lags(uniform_results, irregular_results):
    """Inspect saved lag curves without repeating CCM or changing its Tp=0 test."""
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True, sharey="col")
    baseline_tables = [("original_uniform", "Original uniform", uniform_results),
                       ("irregular_accurate", "Irregular, accurate ages", irregular_results)]
    for row, (baseline_name, title, results) in enumerate(baseline_tables):
        for column, process in enumerate(["linear", "nonlinear", "independent"]):
            ax = axes[row, column]
            path = DATA_DIR / f"r0_{process}_{baseline_name}_CCM_lag_profile.csv"
            selected = results[(results["replicate"] == 0) & (results["process"] == process)
                               & (results["method"] == "CCM")]
            has_scan = "lag_status" in selected and selected["lag_status"].notna().any()
            if path.exists() and has_scan:
                profile = pd.read_csv(path)
                notes = []
                for source, color, label in [("x", "tab:blue", "X → Y"), ("y", "tab:orange", "Y → X")]:
                    curve = profile[profile["source"] == source]
                    lag_years = 1000 * curve["lag_kyr"].to_numpy()
                    mean = curve["rho_mean"].to_numpy()
                    spread = curve["rho_std"].to_numpy()
                    ax.plot(lag_years, mean, color=color, label=label)
                    ax.fill_between(lag_years, mean - spread, mean + spread, color=color, alpha=0.15)
                    edge_rows = selected[selected["source"] == source]
                    if len(edge_rows) == 1:
                        edge = edge_rows.iloc[0]
                        peak = edge.get("optimal_lag_kyr", np.nan)
                        if np.isfinite(peak):
                            ax.axvline(1000 * peak, color=color, ls=":", lw=0.8)
                        notes.append(f"{label}: p0={edge['p_value']:.3g}, {edge.get('lag_status', 'unavailable')}"
                                     + f"\npast={edge.get('lag_past_fraction', np.nan):.0%}, future={edge.get('lag_future_fraction', np.nan):.0%}")
                ax.text(0.03, 0.04, "\n".join(notes), transform=ax.transAxes, fontsize=8,
                        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"})
            else:
                ax.text(0.5, 0.5, "Lag scan unavailable", transform=ax.transAxes, ha="center")
            history_years = (CCM_E - 1) * abs(CCM_TAU) * GRID_STEP * 1000
            ax.axvspan(-history_years, 0, color="0.7", alpha=0.15)
            ax.axvline(0, color="black", ls="--", lw=0.8)
            ax.set_xlim(1000 * GRID_STEP * CCM_LAG_STEPS[0], 1000 * GRID_STEP * CCM_LAG_STEPS[-1])
            ax.set_title(f"{title}: {process}", fontsize=10)
            ax.set_xlabel("Cross-map lag (yr)")
            if column == 0:
                ax.set_ylabel("Cross-map skill (rho)")
    fig.legend(handles=[Line2D([], [], color="tab:blue", label="X → Y"),
                        Line2D([], [], color="tab:orange", label="Y → X"),
                        Patch(color="0.7", alpha=0.3, label="Historical embedding span")],
               loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.035), frameon=False)
    fig.suptitle(f"CCM temporal screen: {CCM_LAG_LIBRARIES} matched libraries; both directions checked independently", fontsize=13)
    fig.text(0.5, 0.015, "Bands: library SD, not confidence. Scan uses common anchors; original p0 does not test the selected peak.",
             ha="center", fontsize=9)
    return fig


def plot_checks(uniform_results, irregular_results, sampling_examples):
    """Plot the inspection outputs; keep every plotting step outside toolbox."""
    reference = pd.concat([uniform_results, irregular_results], ignore_index=True)
    original_figure = plot_baseline_results(uniform_results, reference, "Baseline 1: original uniform process")
    accurate_figure = plot_baseline_results(irregular_results, reference, "Baseline 2: irregular samples with accurate ages, then interpolation")
    figures = [("01_baseline_original_uniform", original_figure), ("02_baseline_irregular_accurate", accurate_figure)]
    processes_order = ["linear", "nonlinear", "independent"]

    fig, axes = plt.subplots(2, 3, figsize=(13, 6), sharex=True)
    middle = (ANALYSIS_OLD + ANALYSIS_YOUNG) / 2
    for column, process_name in enumerate(processes_order):
        original = sampling_examples[process_name]["uniform"]
        accurate = sampling_examples[process_name]["irregular"]
        for row, record in enumerate(["x", "y"]):
            ax = axes[row, column]
            direct = original[record].to_numpy()
            interpolated = accurate[record].to_numpy()
            rmse_ratio = np.sqrt(np.mean((interpolated - direct)**2)) / np.std(direct)
            increment_ratio = np.var(np.diff(interpolated)) / np.var(np.diff(direct))
            ax.plot(original["age"], direct, color="black", lw=1.2, label="Original uniform")
            ax.plot(accurate["age"], interpolated, color="tab:orange", lw=1, ls="--", label="Irregular + interpolation")
            observation = sampling_examples[process_name]["observed"][record]
            ax.scatter(observation["age_obs"], observation["value_obs"], s=10,
                       color="tab:orange", label="Irregular observations", zorder=3)
            ax.set_xlim(middle + 0.6, middle - 0.6)
            ax.set_ylabel(record.upper())
            ax.set_title(f"{process_name}: RMSE/SD = {rmse_ratio:.1%}\nIncrement variance ratio = {increment_ratio:.3f}", fontsize=10)
            ax.set_xlabel("Age (kyr BP)")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Sampling and interpolation: detail view; diagnostics use the full analysis window")
    figures.append(("03_irregular_sampling_effect", fig))
    figures.append(("09_ccm_baseline_lags", plot_ccm_baseline_lags(uniform_results, irregular_results)))

    for name, fig in figures:
        if name.startswith(("01_", "02_", "09_")):
            fig.tight_layout(rect=(0, 0.12, 1, 0.94))
        else:
            fig.tight_layout()
        fig.savefig(FIGURE_DIR / f"{name}.png", dpi=180)
        fig.savefig(FIGURE_DIR / f"{name}.pdf")
        plt.close(fig)
    print(f"Saved {len(figures) + 1} figures (including model truth) to {FIGURE_DIR}", flush=True)


def check_saved_baselines(*, require_complete=True):
    """Require the same data and estimator settings before reusing exact-age fits."""
    saved = json.loads((DATA_DIR / "baseline_parameters.json").read_text())
    names = [
        "SEED", "DURATION", "DT", "TRUE_LAG", "MEMORY_X", "MEMORY_Y", "BACKGROUND_X",
        "BACKGROUND_Y", "COUPLING_GAIN", "BURN_TIMESCALES", "YOUNG_SPACING", "OLD_SPACING",
        "SPACING_JITTER", "MEASUREMENT_NOISE_RATIO", "ANALYSIS_OLD", "ANALYSIS_YOUNG",
        "GRID_STEP", "MAX_LAG", "GC_ORDER", "ALPHA", "N_REPLICATES", "GPDC_NULL_SAMPLES",
        "GPDC_OPTIMIZE_KERNEL", "CCM_E", "CCM_TAU", "CCM_EXCLUSION_RADIUS",
        "CCM_LIBRARY_FRACTION", "CCM_SAMPLES", "CCM_SURROGATES",
        "CCM_LAG_STEPS", "CCM_LAG_LIBRARIES", "CCM_LAG_SIGN_THRESHOLD",
    ]
    for name in names:
        if saved.get(name) != globals()[name]:
            raise ValueError(f"Saved baseline setting {name} differs; run the baselines again.")
    for package, saved_version in saved["packages"].items():
        if version(package) != saved_version:
            raise ValueError(f"{package} changed since the baseline run; rerun both stages.")
    if (saved.get("analysis_revision") != "native_null_cache_and_ccm_object_reuse_v1"
            or saved.get("execution_environment") != execution_environment()):
        raise ValueError("Baseline numerical implementation/platform changed; start a fresh local run.")
    if require_complete:
        if not saved.get("complete", False):
            raise ValueError("Both accurate-age baselines must finish before chronology analysis.")
        expected = {(replicate, process, method, source, target)
                    for replicate in range(N_REPLICATES)
                    for process in ["linear", "nonlinear", "independent"] for method in METHODS
                    for source, target in [("x", "y"), ("y", "x")]}
        for name in ["original_uniform", "irregular_accurate"]:
            table = pd.read_csv(DATA_DIR / f"baseline_{name}_results.csv")
            keys = ["replicate", "process", "method", "source", "target"]
            if table.duplicated(keys).any() or set(table[keys].itertuples(index=False, name=None)) != expected:
                raise ValueError("Saved baseline coverage is incomplete.")
        for replicate in range(N_REPLICATES):
            for process in ["linear", "nonlinear", "independent"]:
                required = [f"r{replicate}_{process}_generated_series.csv",
                            f"r{replicate}_{process}_irregular_accurate_series.csv"]
                required.extend(f"r{replicate}_{process}_{record}_irregular_observed.csv" for record in ["x", "y"])
                if any(not (DATA_DIR / name).exists() for name in required):
                    raise ValueError("Saved baseline source observations are missing; restore them or start a fresh run.")


def prepare_age_comparison():
    """Prepare shared chronologies and paired member inputs before any causal fit."""
    check_saved_baselines()
    if not 1 <= N_CAUSAL_MEMBERS <= N_AGE_MEMBERS or N_WORKERS < 1:
        raise ValueError("Require 1 <= N_CAUSAL_MEMBERS <= N_AGE_MEMBERS and positive N_WORKERS.")
    grid_size = int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / GRID_STEP)) + 1
    age_grid = np.linspace(ANALYSIS_OLD, ANALYSIS_YOUNG, grid_size)
    jump_centers = ANALYSIS_YOUNG + (ANALYSIS_OLD - ANALYSIS_YOUNG) * np.array(JUMP_CENTER_FRACTIONS[:JUMP_COUNT])
    jump_intervals = [(float(center - JUMP_INTERVAL_WIDTH / 2), float(center + JUMP_INTERVAL_WIDTH / 2))
                      for center in jump_centers]
    scenario_rows = []
    for scenario, kind_x, kind_y in AGE_SCENARIOS:
        scenario_rows.append({
            "scenario": scenario, "uncertainty_level": UNCERTAINTY_LEVEL,
            "x_kind": kind_x, "y_kind": kind_y,
            "x_peak_sd_years": 1000 * AGE_SCALES[kind_x], "y_peak_sd_years": 1000 * AGE_SCALES[kind_y],
            "x_jump_count": JUMP_COUNT if kind_x == "J" else 0,
            "y_jump_count": JUMP_COUNT if kind_y == "J" else 0,
            "x_jump_width_kyr": JUMP_INTERVAL_WIDTH if kind_x == "J" else 0,
            "y_jump_width_kyr": JUMP_INTERVAL_WIDTH if kind_y == "J" else 0,
            "x_jump_coverage": JUMP_COVERAGE if kind_x == "J" else 0,
            "y_jump_coverage": JUMP_COVERAGE if kind_y == "J" else 0,
            "true_lag_years": 1000 * TRUE_LAG,
            "n_age_members": N_AGE_MEMBERS, "n_causal_members": N_CAUSAL_MEMBERS,
        })
    pd.DataFrame(scenario_rows).to_csv(DATA_DIR / "scenarios.csv", index=False)
    jobs, age_rows, relative_rows, resampling_rows = [], [], [], []
    chronology_example, examples = {}, {}
    baseline = pd.read_csv(DATA_DIR / "baseline_irregular_accurate_results.csv")

    for replicate in range(N_REPLICATES):
        seed = SEED + 1000 * replicate
        full = pd.read_csv(DATA_DIR / f"r{replicate}_linear_generated_series.csv")
        age = full["age"].to_numpy()
        indices = {}
        for record, offset in [("x", 1), ("y", 2)]:
            indices[record] = sample_irregular_indices(
                age, young_spacing=YOUNG_SPACING, old_spacing=OLD_SPACING,
                jitter=SPACING_JITTER, rng=seed + offset,
            )
            exact = pd.read_csv(DATA_DIR / f"r{replicate}_linear_{record}_irregular_observed.csv")
            if not np.allclose(age[indices[record]], exact["age_obs"], rtol=0, atol=1e-12):
                raise ValueError("Age scenarios must reuse the exact baseline sample positions.")

        # X and Y have independent clocks. Within each clock, reuse random draws
        # across envelope shapes/scales and reuse the resulting ages across models.
        specifications = {}
        for _, kind_x, kind_y in AGE_SCENARIOS:
            for record, kind in [("x", kind_x), ("y", kind_y)]:
                key = f"{record}_{kind}"
                specifications[key] = (record, kind, AGE_SCALES[kind])
        chronology = {}
        for key, (record, kind, scale) in specifications.items():
            offset = 10 if record == "x" else 20
            age_std = make_age_uncertainty(
                age[indices[record]], kind=kind, scale=scale, age_range=(0.0, DURATION),
                exp_shape=AGE_EXP_SHAPE, jump_intervals=jump_intervals if kind == "J" else None,
            )
            observation = simulate_observations(
                age, full[record], indices[record], age_std,
                age_correlation_length=AGE_CORRELATION_LENGTH,
                measurement_noise_ratio=0, rng_age=seed + offset,
                rng_measurement=seed + offset + 1,
            )
            ensemble = draw_age_ensemble(
                observation["observed"], n_members=N_AGE_MEMBERS,
                age_correlation_length=AGE_CORRELATION_LENGTH,
                max_draws=MAX_AGE_DRAWS, rng=seed + offset + 100,
            )
            chronology[key] = {"observation": observation, "ensemble": ensemble}
            truth_age = observation["truth"]["age_true"].to_numpy()
            active = age_std > ACTIVE_AGE_SD_TOLERANCE
            accepted = ensemble["ages"]
            diagnostic = {
                "replicate": replicate, "record": key, "kind": kind, "scale_years": 1000 * scale,
                "n_samples": len(age_std), "n_active_age_samples": int(active.sum()),
                "input_rms_sd_years": 1000 * float(np.sqrt(np.mean(age_std**2))),
                **ensemble["diagnostics"],
                "sample_fraction_inside_90_interval": np.nan,
                "active_sample_fraction_inside_90_interval": np.nan,
                "median_sd_ratio": np.nan, "fraction_covering_analysis_window": np.nan,
            }
            if len(accepted):
                low, high = np.percentile(accepted, [5, 95], axis=0)
                inside = (truth_age >= low) & (truth_age <= high)
                diagnostic["sample_fraction_inside_90_interval"] = float(np.mean(inside))
                if np.any(active):
                    diagnostic["active_sample_fraction_inside_90_interval"] = float(np.mean(inside[active]))
                    diagnostic["median_sd_ratio"] = float(np.median(np.std(accepted[:, active], axis=0) / age_std[active]))
                diagnostic["fraction_covering_analysis_window"] = float(np.mean(
                    (accepted[:, 0] >= ANALYSIS_OLD) & (accepted[:, -1] <= ANALYSIS_YOUNG)))
            age_rows.append(diagnostic)
            np.savez_compressed(
                DATA_DIR / f"r{replicate}_{key}_chronology.npz", ages=accepted,
                sample_id=ensemble["sample_id"], age_true=truth_age,
                age_obs=observation["observed"]["age_obs"], age_std=age_std,
            )
            print(f"Chronology r{replicate} {key}: {len(accepted)}/{N_AGE_MEMBERS} accepted", flush=True)

        if replicate == 0:
            chronology_example = chronology
            examples = {
                "accurate": pd.read_csv(DATA_DIR / "r0_linear_irregular_accurate_series.csv"),
                "resampled": [], "observed": {},
                "baseline": baseline[baseline["replicate"] == 0].copy(),
            }

        for scenario, kind_x, kind_y in AGE_SCENARIOS:
            clocks = {record: chronology[f"{record}_{kind}"]
                      for record, kind in [("x", kind_x), ("y", kind_y)]}
            available = min(len(clocks["x"]["ensemble"]["ages"]), len(clocks["y"]["ensemble"]["ages"]))
            # Truth is used only here for diagnostics, never supplied to the causal methods.
            relative_errors = []
            for member in range(available):
                errors = {}
                for record in ["x", "y"]:
                    true_age = clocks[record]["observation"]["truth"]["age_true"].to_numpy()
                    error = clocks[record]["ensemble"]["ages"][member] - true_age
                    errors[record] = np.interp(-age_grid, -true_age, error)
                relative_errors.append(errors["y"] - errors["x"])
            relative = np.asarray(relative_errors)
            relative_rows.append({
                "replicate": replicate, "scenario": scenario, "n_age_pairs": available,
                "relative_error_rms_years": 1000 * float(np.sqrt(np.mean(relative**2))) if available else np.nan,
                "fraction_relative_error_above_true_lag": float(np.mean(np.abs(relative) > TRUE_LAG)) if available else np.nan,
            })

            for process_name in ["linear", "nonlinear", "independent"]:
                measured = {}
                for record in ["x", "y"]:
                    exact = pd.read_csv(DATA_DIR / f"r{replicate}_{process_name}_{record}_irregular_observed.csv")
                    expected_ids = clocks[record]["observation"]["observed"]["sample_id"]
                    if not np.array_equal(exact["sample_id"], expected_ids):
                        raise ValueError("Every process must retain the same sample identities.")
                    if not np.allclose(exact["age_obs"], age[indices[record]], rtol=0, atol=1e-12):
                        raise ValueError("Every process must use the baseline's exact sample ages.")
                    measured[record] = exact["value_obs"].to_numpy()
                    observed = clocks[record]["observation"]["observed"].copy()
                    observed["value_obs"] = measured[record]
                    observed["measurement_std"] = exact["measurement_std"].to_numpy()
                    observed.to_csv(DATA_DIR / f"r{replicate}_{process_name}_{scenario}_{record}_observed.csv", index=False)
                    if replicate == 0 and process_name == "linear" and scenario == "L-J":
                        examples["observed"][record] = observed

                for member in range(N_CAUSAL_MEMBERS):
                    label = f"r{replicate}_{process_name}_{scenario}_m{member}"
                    job = {"replicate": replicate, "process": process_name, "scenario": scenario,
                           "member": member, "label": label, "seed": seed + 300,
                           "status": "ready", "error": ""}
                    if member >= available:
                        job["status"] = "chronology_error"
                        job["error"] = "Requested chronology member was not accepted."
                    else:
                        try:
                            sampled = resample_pair(
                                clocks["x"]["ensemble"]["ages"][member], measured["x"],
                                clocks["y"]["ensemble"]["ages"][member], measured["y"], age_grid=age_grid,
                            )
                            save_input(sampled["series"], DATA_DIR / f"{label}_series.csv")
                            resampling_rows.append({**job, **sampled["diagnostics"]})
                            if replicate == 0 and process_name == "linear" and scenario == "L-J":
                                examples["resampled"].append(sampled["series"])
                        except ValueError as error:
                            job["status"] = "resampling_error"
                            job["error"] = str(error)
                    jobs.append(job)
    pd.DataFrame(age_rows).to_csv(DATA_DIR / "age_diagnostics.csv", index=False)
    pd.DataFrame(relative_rows).to_csv(DATA_DIR / "relative_age_diagnostics.csv", index=False)
    pd.DataFrame(resampling_rows).to_csv(DATA_DIR / "resampling_diagnostics.csv", index=False)
    pd.DataFrame(jobs).to_csv(DATA_DIR / "age_member_inputs.csv", index=False)
    return jobs, pd.DataFrame(age_rows), chronology_example, examples, baseline


def analyze_age_member(job):
    """Run the fixed methods on one prepared member; all four see identical data."""
    label = job["label"]
    if job["status"] != "ready":
        edges = pd.DataFrame([
            {"method": method, "source": source, "target": target, "score": np.nan,
             "score_name": "unavailable", "edge_detected": None, "p_value": np.nan,
             "p_kind": "unavailable", "status": job["status"], "unoriented": False,
             "error": job["error"]}
            for method in METHODS for source, target in [("x", "y"), ("y", "x")]
        ])
    else:
        series = pd.read_csv(DATA_DIR / f"{label}_series.csv")
        diagnostic_name = f"{label}_method_diagnostics.jsonl"
        (DATA_DIR / diagnostic_name).write_text("", encoding="utf-8")
        with threadpool_limits(limits=1):
            edges, _ = analyze_pair(series, seed=job["seed"], label=label,
                                    diagnostics_name=diagnostic_name)
    for column in ["replicate", "process", "scenario", "member"]:
        edges[column] = job[column]
    edges["input_status"] = job["status"]
    edges["edge_detected"] = edges["edge_detected"].astype("boolean")
    edges["unoriented"] = edges["unoriented"].astype("boolean")
    path = DATA_DIR / f"{label}_edges.csv"
    temporary = path.with_suffix(".partial.csv")
    edges.to_csv(temporary, index=False)
    temporary.replace(path)
    return edges


def summarize_ccm_filter(members):
    """Compare original and retained votes on the same usable member pairs."""
    ccm = members[members["method"] == "CCM"]
    rows = []
    for (replicate, process, scenario), group in ccm.groupby(["replicate", "process", "scenario"], sort=False):
        raw_all, raw_paired, filtered, unresolved = [], [], [], []
        incompatible = []
        for member, pair in group.groupby("member", sort=False):
            forward = pair[(pair["source"] == "x") & (pair["target"] == "y")]
            reverse = pair[(pair["source"] == "y") & (pair["target"] == "x")]
            if len(forward) != 1 or len(reverse) != 1:
                continue
            edges = [forward.iloc[0], reverse.iloc[0]]
            raw_valid = all(edge["status"] == "ok" and pd.notna(edge.get("raw_edge_detected")) for edge in edges)
            if not raw_valid:
                continue
            original_votes = [bool(edge["raw_edge_detected"]) for edge in edges]
            raw_all.append(original_votes)
            if any(pd.isna(edge["edge_detected"]) for edge in edges):
                continue
            final_votes = [bool(edge["edge_detected"]) for edge in edges]
            if any(final and not original for final, original in zip(final_votes, original_votes)):
                raise ValueError("The CCM lag screen cannot promote an originally absent arrow.")
            raw_paired.append(original_votes)
            filtered.append(final_votes)
            unresolved.append([edge.get("unoriented", False) == True for edge in edges])
            incompatible.append([original and edge.get("lag_status") == "incompatible"
                                 for original, edge in zip(original_votes, edges)])
        n_raw = len(raw_all)
        n_paired = len(raw_paired)
        for direction, (source, target) in enumerate([("x", "y"), ("y", "x")]):
            row = {"replicate": replicate, "process": process, "scenario": scenario,
                   "source": source, "target": target, "n_requested": N_CAUSAL_MEMBERS,
                   "n_raw_valid": n_raw, "n_paired_valid": n_paired,
                   "q_raw_all": np.nan, "q_raw": np.nan, "q_filtered": np.nan,
                   "unoriented_fraction": np.nan, "incompatible_fraction": np.nan}
            if n_raw:
                row["q_raw_all"] = float(np.mean(np.asarray(raw_all)[:, direction]))
            if n_paired:
                row["q_raw"] = float(np.mean(np.asarray(raw_paired)[:, direction]))
                row["q_filtered"] = float(np.mean(np.asarray(filtered)[:, direction]))
                row["unoriented_fraction"] = float(np.mean(np.asarray(unresolved)[:, direction]))
                row["incompatible_fraction"] = float(np.mean(np.asarray(incompatible)[:, direction]))
                if row["q_filtered"] > row["q_raw"]:
                    raise ValueError("Paired CCM support cannot increase after filtering.")
            rows.append(row)
    return pd.DataFrame(rows)


def plot_ccm_filter_comparison(comparison):
    """Display paired-denominator CCM support and unresolved original responses."""
    names = [item[0] for item in AGE_SCENARIOS]
    positions = np.arange(len(names))
    fig, axes = plt.subplots(3, 2, figsize=(17, 10), sharex=True, sharey=True)
    for row, process in enumerate(["linear", "nonlinear", "independent"]):
        for column, source in enumerate(["x", "y"]):
            ax = axes[row, column]
            subset = comparison[(comparison["replicate"] == 0) & (comparison["process"] == process)
                                & (comparison["source"] == source)]
            for position, scenario in enumerate(names):
                selected = subset[subset["scenario"] == scenario]
                if selected.empty:
                    ax.text(position, 0.5, "Missing", rotation=90, ha="center", fontsize=8)
                    continue
                result = selected.iloc[0]
                complete = result["n_paired_valid"] >= MIN_VALID_FRACTION * result["n_requested"]
                for field, offset, color in [("q_raw", -0.25, "0.65"),
                                             ("q_filtered", 0, "tab:blue" if source == "x" else "tab:orange"),
                                             ("unoriented_fraction", 0.25, "#7FC97F")]:
                    value = result[field]
                    if np.isfinite(value):
                        ax.bar(position + offset, value, width=0.24, color=color,
                               alpha=1 if complete else 0.35, hatch=None if complete else "//")
                ax.text(position, 1.06, f"{int(result['n_paired_valid'])}/{int(result['n_requested'])}",
                        ha="center", fontsize=8, color="0.3" if complete else "tab:red")
            ax.axhline(SUPPORT_THRESHOLD, color="black", ls="--", lw=0.8)
            direction = "X → Y" if source == "x" else "Y → X"
            ax.set_title(f"{process.capitalize()}: {direction}", fontsize=11)
            ax.set_ylim(0, 1.15)
            ax.set_xlim(-0.6, len(names) - 0.4)
            ax.set_xticks(positions, names, fontsize=8)
            ax.tick_params(axis="x", labelbottom=True)
            ax.set_xlabel("Age-error pair: X type – Y type")
            if column == 0:
                ax.set_ylabel("Fraction of paired usable members")
    fig.legend(handles=[Patch(color="0.65", label="Original Tp=0 significant"),
                        Patch(color="tab:blue", label="Retained X → Y"),
                        Patch(color="tab:orange", label="Retained Y → X"),
                        Patch(color="#7FC97F", label="Significant response, lag unresolved")],
               loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.025), frameon=False, fontsize=9)
    fig.suptitle("CCM before and after the lag screen: identical member denominators\n"
                 f"Level {UNCERTAINTY_LEVEL}: C/L/E peak SD {1000 * AGE_SCALES['C']:g} yr; "
                 f"J {JUMP_COUNT} peaks at {1000 * JUMP_PEAK_SD:g} yr, {JUMP_COVERAGE:.0%} coverage; "
                 f"coupled-model true lag {1000 * TRUE_LAG:g} yr", fontsize=12)
    fig.text(0.5, 0.008, "Counts are paired usable / requested. Lag failures are excluded from both q values; original-all support is also saved in the table.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.10, 1, 0.93))
    for extension in ["png", "pdf"]:
        fig.savefig(FIGURE_DIR / f"10_ccm_filter_support.{extension}", dpi=180)
    plt.close(fig)


def summarize_age_results(members, baseline):
    """Separate within-record support from truth-based inspection of that record."""
    summaries = []
    for (replicate, process_name, scenario), group in members.groupby(["replicate", "process", "scenario"], sort=False):
        reference = baseline[(baseline["replicate"] == replicate) & (baseline["process"] == process_name)]
        summary = summarize_ensemble(
            group, n_members=N_CAUSAL_MEMBERS, baseline=reference,
            support_threshold=SUPPORT_THRESHOLD, min_valid_fraction=MIN_VALID_FRACTION,
        )
        summary["replicate"] = replicate
        summary["process"] = process_name
        summary["scenario"] = scenario
        summaries.append(summary)
    summaries = pd.concat(summaries, ignore_index=True)
    summaries.to_csv(DATA_DIR / "ensemble_summary.csv", index=False)
    ccm_comparison = summarize_ccm_filter(members)
    if not ccm_comparison.empty:
        ccm_comparison.to_csv(DATA_DIR / "ccm_filter_comparison.csv", index=False)

    # With R=1 these are truth checks for this record, not estimates of error rates.
    truth_checks = []
    for row in summaries.to_dict("records"):
        valid = row["decision"] != "insufficient_valid_members"
        forward = row["q_x_to_y"] >= SUPPORT_THRESHOLD
        reverse = row["q_y_to_x"] >= SUPPORT_THRESHOLD
        unresolved = row["high_unoriented"]
        independent = row["process"] == "independent"
        if independent:
            exact = not forward and not reverse and not unresolved
            false_direction = forward or reverse
            correct_fraction = row["fraction_neither"]
        else:
            exact = forward and not reverse and not unresolved
            false_direction = reverse
            correct_fraction = row["fraction_x_only"]
        truth_checks.append({
            "replicate": row["replicate"], "process": row["process"], "scenario": row["scenario"],
            "method": row["method"], "valid": valid,
            "exact_thresholded_graph": bool(exact) if valid else None,
            "high_support_false_direction": bool(false_direction) if valid else None,
            "fraction_correct_member_graph": correct_fraction,
        })
    pd.DataFrame(truth_checks).to_csv(DATA_DIR / "ensemble_truth_checks.csv", index=False)
    if N_REPLICATES > 1:
        performance = []
        for process_name, group in summaries.groupby("process", sort=False):
            truth = [] if process_name == "independent" else [("x", "y")]
            table = evaluate_replicates(group, truth_edges=truth)
            table["process"] = process_name
            performance.append(table)
        pd.concat(performance, ignore_index=True).to_csv(DATA_DIR / "replicate_metrics.csv", index=False)
    return summaries


def save_member_checkpoint(tables):
    """Replace the aggregate only after a complete checkpoint has been written."""
    path = DATA_DIR / "member_results.csv"
    temporary = path.with_suffix(".partial.csv")
    pd.concat(tables, ignore_index=True).to_csv(temporary, index=False)
    temporary.replace(path)


def write_run_status(state, completed, total, error=""):
    """Keep a small readable progress file for the current fresh run."""
    status = {
        "state": state, "completed_member_pairs": completed, "requested_member_pairs": total,
        "uncertainty_level": UNCERTAINTY_LEVEL, "n_age_members": N_AGE_MEMBERS,
        "n_causal_members": N_CAUSAL_MEMBERS, "n_workers": N_WORKERS,
        "updated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "error": error,
    }
    path = DATA_DIR / "run_status.json"
    temporary = path.with_suffix(".partial.json")
    temporary.write_text(json.dumps(status, indent=2), encoding="utf-8")
    temporary.replace(path)


def run_age_comparison():
    """Resume complete members of this run, preserving failures and reverse arrows."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    total = N_REPLICATES * 3 * len(AGE_SCENARIOS) * N_CAUSAL_MEMBERS
    parameters = {name: value for name, value in globals().items()
                  if name.isupper() and isinstance(value, (str, int, float, list, dict))}
    parameters = json.loads(json.dumps(parameters))
    parameters["interpretation"] = "Within-record working-age support; R=1 does not calibrate causal correctness or false-positive rates."
    parameters["reference"] = "The same proxy observations resampled with accurate ages."
    parameters["scale_definition"] = "Per-record Gaussian proposal SD: C constant, L/E old-end maximum; J peak SD fixed while interval count changes. Not matched RMS."
    parameters["pairing"] = "Independent X/Y clocks; shared standard draws across envelope shapes and process forms."
    parameters["packages"] = json.loads((DATA_DIR / "baseline_parameters.json").read_text())["packages"]
    parameters["packages"]["threadpoolctl"] = version("threadpoolctl")
    parameters["numerical_threads_per_worker"] = 1
    parameters["execution_environment"] = execution_environment()
    parameters["analysis_revision"] = "native_null_cache_and_ccm_object_reuse_v1"
    parameters["member_run"] = "Resume only complete members from this configuration; scientific failures are retained."
    parameter_path = DATA_DIR / "age_comparison_parameters.json"
    if parameter_path.exists():
        if not RESUME:
            raise ValueError("Chronology outputs exist. Set RESUME=True or preserve them in another folder.")
        check_resume_parameters(parameter_path, parameters)
    elif any(DATA_DIR.glob("r*_m*_edges.csv")):
        raise ValueError("Member files lack their run parameters; do not reuse them without provenance.")
    parameter_path.write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    tables = []
    write_run_status("preparing", 0, total)
    try:
        with threadpool_limits(limits=1):
            jobs, age_diagnostics, chronology, examples, baseline = prepare_age_comparison()
        if PREPARE_ONLY:
            write_run_status("prepared", 0, total)
            print(f"Prepared {len(jobs)} chronology inputs; no age-member methods were run.", flush=True)
            return
        pending = []
        for job in jobs:
            completed = read_completed_member(job)
            if completed is None:
                pending.append(job)
            else:
                tables.append(completed)
        print(f"Age members: {len(tables)} complete; {len(pending)} pending with {N_WORKERS} workers.", flush=True)
        write_run_status("running", len(tables), total)
        if N_WORKERS == 1:
            for job in pending:
                tables.append(analyze_age_member(job))
                if len(tables) % 10 == 0:
                    save_member_checkpoint(tables)
                    write_run_status("running", len(tables), total)
                    print(f"Age members completed: {len(tables)}/{total}", flush=True)
        elif pending:
            with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
                futures = [pool.submit(analyze_age_member, job) for job in pending]
                for future in as_completed(futures):
                    tables.append(future.result())
                    if len(tables) % 10 == 0:
                        save_member_checkpoint(tables)
                        write_run_status("running", len(tables), total)
                        print(f"Age members completed: {len(tables)}/{total}", flush=True)
        members = pd.concat(tables, ignore_index=True)
        members = members.sort_values(["replicate", "process", "scenario", "member", "method", "source"])
        save_member_checkpoint([members])
        write_run_status("summarizing", len(tables), total)
        summaries = summarize_age_results(members, baseline)
        plot_age_checks(
            summaries, age_diagnostics, chronology, examples,
            scenario_names=[item[0] for item in AGE_SCENARIOS], methods=METHODS,
            processes_order=["linear", "nonlinear", "independent"], figure_dir=FIGURE_DIR,
            n_causal_members=N_CAUSAL_MEMBERS, support_threshold=SUPPORT_THRESHOLD,
        )
        plot_ccm_filter_comparison(pd.read_csv(DATA_DIR / "ccm_filter_comparison.csv"))
        write_run_status("complete", len(tables), total)
        print(summaries[["process", "scenario", "method", "n_valid", "q_x_to_y", "q_y_to_x", "decision"]].to_string(index=False))
    except BaseException as error:
        if tables:
            save_member_checkpoint(tables)
        write_run_status("interrupted_or_failed", len(tables), total, f"{type(error).__name__}: {error}")
        raise



def plot_age_checks(
    summaries, age_diagnostics, chronology, examples, *, scenario_names,
    methods, processes_order, figure_dir, n_causal_members, support_threshold,
):
    """Save figures 04–08 for replicate zero, with all failures visible.

    Chronology keys are x_C, y_C, and the corresponding L/E/J
    keys. Each value contains observation and ensemble dictionaries.
    Summary member-state fractions are conditional on valid paired fits;
    figure 07 converts them to the requested-member denominator.
    """
    figure_dir = Path(figure_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    summaries = summaries.loc[summaries["replicate"] == 0].copy()
    age_diagnostics = age_diagnostics.copy()
    if "replicate" in age_diagnostics:
        age_diagnostics = age_diagnostics.loc[age_diagnostics["replicate"] == 0]
    if summaries.duplicated(["process", "scenario", "method"]).any():
        raise ValueError("Age plots require one summary per process, scenario and method.")
    baseline = examples["baseline"]
    if "replicate" in baseline:
        baseline = baseline.loc[baseline["replicate"] == 0]
    labels = [name.replace("-", "–") for name in scenario_names]
    positions = np.arange(len(scenario_names))
    figures = []
    level_label = (f"Level {UNCERTAINTY_LEVEL}: C/L/E peak SD {1000 * AGE_SCALES['C']:g} yr; "
                   f"J peak SD {1000 * AGE_SCALES['J']:g} yr, {JUMP_COUNT} peaks, {JUMP_COVERAGE:.0%} coverage")
    lag_label = f"coupled-model true lag = {1000 * TRUE_LAG:g} yr"

    # Figure 04: observed center errors and accepted chronology intervals.
    fig, axes = plt.subplots(4, 2, figsize=(12, 10), sharex=True, sharey=True)
    for row, kind in enumerate(["C", "L", "E", "J"]):
        for column, record in enumerate(["x", "y"]):
            ax = axes[row, column]
            item = chronology[f"{record}_{kind}"]
            observation = item["observation"]
            ensemble = item["ensemble"]
            true_age = observation["truth"]["age_true"].to_numpy()
            center_age = observation["observed"]["age_obs"].to_numpy()
            accepted_ages = np.asarray(ensemble["ages"])
            if len(accepted_ages):
                errors = 1000 * (accepted_ages - true_age)
                low, high = np.percentile(errors, [5, 95], axis=0)
                ax.fill_between(true_age, low, high, color="tab:blue", alpha=0.2,
                                label="5–95% age-member interval")
                for index, error in enumerate(errors[:5]):
                    ax.plot(true_age, error, color="0.6", lw=0.5, alpha=0.7,
                            label="First five age members" if index == 0 else None)
            else:
                ax.text(0.5, 0.7, "No accepted chronologies", transform=ax.transAxes,
                        ha="center", color="tab:red")
            ax.plot(true_age, 1000 * (center_age - true_age), color="black", lw=0.8,
                    label="Observed center error")
            ax.axhline(0, color="tab:red", ls="--", lw=0.8, label="True age")
            rate = ensemble["diagnostics"].get("acceptance_rate", np.nan)
            rate_text = f"{rate:.0%}" if np.isfinite(rate) else "unavailable"
            ax.set_title(f"{kind}, {record.upper()} — acceptance {rate_text}; {len(accepted_ages)} members", fontsize=10)
            ax.set_xlim(float(np.max(true_age)), float(np.min(true_age)))
            if column == 0:
                ax.set_ylabel("Age error (yr)")
            if row == 3:
                ax.set_xlabel("Age (kyr BP)")
    handles, legend_labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="lower center", ncol=2,
               bbox_to_anchor=(0.5, 0.025), frameon=False, fontsize=9)
    fig.suptitle("Working chronology ensembles\n" + level_label, fontsize=12)
    fig.text(0.5, 0.01, "Intervals describe accepted age members around one observed chronology; they are not calibrated causal confidence intervals.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.09, 1, 0.93))
    figures.append(("04_age_ensemble", fig))

    # Figure 05: every member uses the same full analysis window and grid.
    accurate = examples["accurate"]
    limits = [(float(accurate["age"].max()), float(accurate["age"].min())), (15.6, 14.4)]
    fig, axes = plt.subplots(2, 2, figsize=(13, 6), sharex="col", sharey="row")
    for row, record in enumerate(["x", "y"]):
        for column, age_limits in enumerate(limits):
            ax = axes[row, column]
            for member, series in enumerate(examples["resampled"]):
                ax.plot(series["age"], series[record], color="tab:blue", alpha=0.3, lw=0.8,
                        label="Chronology members" if member == 0 else None)
            ax.plot(accurate["age"], accurate[record], color="black", lw=0.9,
                    label="Accurate-age reference", zorder=3)
            observed = examples["observed"][record]
            ax.scatter(observed["age_obs"], observed["value_obs"], color="tab:orange",
                       s=3 if column == 0 else 12, alpha=0.65, label="Observed ages and values", zorder=4)
            ax.set_xlim(*age_limits)
            ax.set_title("Full common window" if column == 0 else "Fixed detail window", fontsize=11)
            ax.set_xlabel("Age (kyr BP)")
            ax.set_ylabel(f"{record.upper()} (a.u.)")
    axes[0, 1].legend(fontsize=8, loc="best")
    fig.suptitle(f"L–J: uniform resampling of {n_causal_members} chronology members\n"
                 + level_label + "; " + lag_label, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    figures.append(("05_resampling", fig))

    # Figure 06: q is a frequency among valid paired fits, not correctness.
    fig, axes = plt.subplots(len(processes_order), len(methods), figsize=(18, 10),
                             sharex=True, sharey=True, squeeze=False)
    for row, process in enumerate(processes_order):
        for column, method in enumerate(methods):
            ax = axes[row, column]
            subset = summaries[(summaries["process"] == process) & (summaries["method"] == method)]
            for position, scenario in enumerate(scenario_names):
                selected = subset.loc[subset["scenario"] == scenario]
                if selected.empty:
                    ax.text(position, 0.5, "Missing", ha="center", rotation=90, fontsize=8)
                    continue
                result = selected.iloc[0]
                incomplete = result["decision"] == "insufficient_valid_members"
                for direction, offset, color in [("x_to_y", -0.17, "tab:blue"),
                                                  ("y_to_x", 0.17, "tab:orange")]:
                    q = result[f"q_{direction}"]
                    if np.isfinite(q):
                        ax.bar(position + offset, q, width=0.32, color=color,
                               alpha=0.3 if incomplete else 0.85,
                               hatch="//" if incomplete else None)
                if result["n_valid"] < result["n_requested"]:
                    ax.text(position, 1.07, f"{int(result['n_valid'])}/{int(result['n_requested'])}",
                            ha="center", fontsize=7, color="tab:red" if incomplete else "0.25")
                if result["high_unoriented"]:
                    ax.text(position, 1.01, "U", ha="center", fontsize=8, color="0.25")
            flags = []
            for source in ["x", "y"]:
                selected = baseline[(baseline["process"] == process) & (baseline["method"] == method)
                                    & (baseline["source"] == source)]
                if len(selected) != 1 or selected.iloc[0]["status"] != "ok" or pd.isna(selected.iloc[0]["edge_detected"]):
                    flags.append("unavailable")
                else:
                    flags.append("yes" if bool(selected.iloc[0]["edge_detected"]) else "no")
            ax.axhline(support_threshold, color="black", ls="--", lw=0.8)
            ax.set_xticks(positions, labels, fontsize=7)
            ax.tick_params(axis="x", labelbottom=True)
            for boundary in [3.5, 7.5, 11.5]:
                ax.axvline(boundary, color="0.8", lw=0.6)
            ax.set_xlabel("Chronology pair: X model–Y model", fontsize=8)
            ax.set_xlim(-0.6, len(scenario_names) - 0.4)
            ax.set_ylim(0, 1.16)
            ax.set_title(f"{process.capitalize()}: {method}\nAccurate ages: X→Y {flags[0]}; Y→X {flags[1]}", fontsize=9)
            if column == 0:
                ax.set_ylabel("Direction support q")
    fig.legend(handles=[Patch(color="tab:blue", label="X → Y"), Patch(color="tab:orange", label="Y → X")],
               loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.035), frameon=False)
    fig.suptitle(f"Direction support: {lag_label}; {n_causal_members} requested chronology members\n"
                 + level_label, fontsize=12)
    fig.text(0.5, 0.015, "Dashed line: reporting threshold, not calibrated confidence. U: high unresolved support. Counts/hatching mark incomplete members.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.08, 1, 0.91))
    figures.append(("06_direction_support", fig))

    # Figure 07: include unavailable members in the requested denominator.
    state_columns = ["fraction_x_only", "fraction_y_only", "fraction_both",
                     "fraction_neither", "fraction_unoriented"]
    state_labels = ["X → Y only", "Y → X only", "Both directions", "No retained arrows",
                    "Unoriented contemporaneous link (PCMCI+)\nor inconclusive lag screen (CCM)", "Failed / missing"]
    state_colors = ["#4682B4", "#FDC086", "#B084CC", "#D9D9D9", "#7FC97F", "#666666"]
    fig, axes = plt.subplots(len(processes_order), len(methods), figsize=(18, 10),
                             sharex=True, sharey=True, squeeze=False)
    for row, process in enumerate(processes_order):
        for column, method in enumerate(methods):
            ax = axes[row, column]
            subset = summaries[(summaries["process"] == process) & (summaries["method"] == method)]
            for position, scenario in enumerate(scenario_names):
                selected = subset.loc[subset["scenario"] == scenario]
                fractions = np.zeros(6)
                if selected.empty:
                    fractions[-1] = 1
                else:
                    result = selected.iloc[0]
                    valid_fraction = float(result["valid_fraction"])
                    if result["n_valid"] > 0:
                        fractions[:5] = result[state_columns].to_numpy(dtype=float) * valid_fraction
                    fractions[-1] = 1 - valid_fraction
                if not np.all(np.isfinite(fractions)) or not np.isclose(np.sum(fractions), 1, atol=1e-10):
                    raise ValueError("Requested-member state fractions must be finite and sum to one.")
                bottom = 0.0
                for fraction, color in zip(fractions, state_colors):
                    ax.bar(position, fraction, bottom=bottom, width=0.78, color=color, edgecolor="white", linewidth=0.3)
                    bottom = bottom + fraction
            ax.set_xticks(positions, labels, fontsize=7)
            ax.tick_params(axis="x", labelbottom=True)
            for boundary in [3.5, 7.5, 11.5]:
                ax.axvline(boundary, color="0.8", lw=0.6)
            ax.set_xlabel("Chronology pair: X model–Y model", fontsize=8)
            ax.set_xlim(-0.6, len(scenario_names) - 0.4)
            ax.set_ylim(0, 1.02)
            ax.set_title(f"{process.capitalize()}: {method}", fontsize=11)
            if column == 0:
                ax.set_ylabel("Fraction of requested members")
    fig.legend(handles=[Patch(color=color, label=label) for color, label in zip(state_colors, state_labels)],
               loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.02), frameon=False, fontsize=10)
    fig.suptitle(f"Complete member outcomes: {lag_label}; {n_causal_members} requested chronology members\n"
                 + level_label, fontsize=12)
    fig.text(0.5, 0.008, "Green: PCMCI+ lag=0 links with unresolved/conflicting direction, or significant CCM responses with inconclusive lag screening. Arrows may coexist.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.10, 1, 0.91))
    figures.append(("07_member_states", fig))

    # Figure 08: coverage along correlated sample positions is not calibration.
    diagnostic_labels = age_diagnostics["record"].astype(str).str.replace("_", "\n", regex=False).tolist()
    diagnostic_positions = np.arange(len(age_diagnostics))
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    axes[0].bar(diagnostic_positions, age_diagnostics["acceptance_rate"], color="tab:blue", width=0.65)
    axes[0].set_ylim(0, 1.08)
    axes[0].set_ylabel("Accepted / proposed")
    axes[0].set_title("Whole-curve order constraint: proposal acceptance", fontsize=11)
    axes[1].bar(diagnostic_positions - 0.17, age_diagnostics["sample_fraction_inside_90_interval"],
                width=0.32, color="0.65", label="All sample ages")
    axes[1].bar(diagnostic_positions + 0.17, age_diagnostics["active_sample_fraction_inside_90_interval"],
                width=0.32, color="tab:blue", label="Only samples with nonzero input SD")
    axes[1].axhline(0.9, color="black", ls="--", lw=0.8, label="Nominal interval content (90%)")
    axes[1].set_ylim(0, 1.08)
    axes[1].set_ylabel("Fraction of true sample ages\ninside the 5–95% member interval")
    axes[1].set_title("Along-record inclusion, not coverage across independent records", fontsize=11)
    axes[1].legend(fontsize=9, loc="lower left")
    axes[2].scatter(diagnostic_positions, age_diagnostics["median_sd_ratio"], color="tab:blue", s=30)
    axes[2].axhline(1, color="black", ls="--", lw=0.8)
    axes[2].set_ylabel("Median accepted SD / input SD\n(nonzero-SD samples)")
    axes[2].set_title("Accepted chronology spread relative to the unconstrained input", fontsize=11)
    axes[2].set_xticks(diagnostic_positions, diagnostic_labels, fontsize=9)
    axes[2].set_xlim(-0.6, len(age_diagnostics) - 0.4)
    for ax in axes:
        ax.grid(axis="y", alpha=0.15)
    fig.suptitle(f"Chronology diagnostics: {N_AGE_MEMBERS} members per observed record\n" + level_label, fontsize=12)
    fig.text(0.5, 0.015, f"{N_AGE_MEMBERS}-member intervals and SDs describe this record. Correlated ages are not independent coverage trials; zero-SD sections can inflate inclusion.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 0.91))
    figures.append(("08_age_diagnostics", fig))

    saved_paths = []
    for name, fig in figures:
        for extension in ["png", "pdf"]:
            path = figure_dir / f"{name}.{extension}"
            fig.savefig(path, dpi=180)
            saved_paths.append(path)
        plt.close(fig)
    return saved_paths


def main():
    prepare_gpdc_caches()
    if RUN_BASELINES:
        run_baselines()
    if RUN_AGE_ENSEMBLES:
        run_age_comparison()


if __name__ == "__main__":
    main()
