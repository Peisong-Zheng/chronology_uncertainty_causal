"""Compare all sixteen ordered chronology-type pairs across four uncertainty levels.

C/L/E change nominal peak SD; J changes affected interval coverage at fixed
peak SD and width. All observations and exact-age references come from 02.
This is one process realization with 10 analyzed age members, not an estimate
of error rates across independent records. Run with python in the project environment.
"""

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from importlib.metadata import version
import json
import time
import platform

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from toolbox.observation import make_age_uncertainty, simulate_observations, draw_age_ensemble
from toolbox.causal import resample_pair, run_pcmciplus, run_lkif, run_ccm, run_granger, prepare_gpdc_null
from toolbox.evaluation import summarize_ensemble


RUN_CAUSAL = True             # False redraws existing summaries without preparing or fitting data.
RESUME = True                 # Keep completed members of this configuration, including failed fits.
SEED = 20261009               # Same independent process realization and random streams as 02.
DURATION = 30.0               # Full source record, kyr; uncertainty envelopes use this age range.
DT = 0.020                    # Native and analysis step, kyr: 20 years.
TRUE_LAG = 0.020              # Known model delay, kyr; used for diagnostics only.
ANALYSIS_OLD = 29.0           # Fixed old boundary of the analysis window, kyr BP.
ANALYSIS_YOUNG = 1.0          # Fixed young boundary; no member-specific shortening or extrapolation.
N_AGE_MEMBERS = 10            # Generate ten ordered working chronologies per record/type/level.
N_CAUSAL_MEMBERS = 10         # Analyze all ten; one process realization remains a pilot.
N_WORKERS = 32               # Fastest tested on this Windows PC; one numerical thread per worker.
MAX_AGE_DRAWS = 10000         # Maximum whole-curve proposals; missing members remain failures.
AGE_CORRELATION_LENGTH = 1 / 60  # Correlation length in normalized sample order, as in 02.
AGE_EXP_SHAPE = 2.0           # Fixed exponential curvature; only its SD multiplier changes.
AMPLITUDES_YEARS = [10, 20, 40, 80]  # Per-record proposal SD: C constant, L/E old-end maximum.
J_PEAK_COUNTS = [1, 3, 5, 7]  # Levels 1–4: nested intervals cover 10%, 30%, 50%, 70% of the window.
J_PEAK_SD_YEARS = 10.0       # Fixed J maximum proposal SD, years; not an absolute error bound.
J_INTERVAL_WIDTH = 2.8       # Fixed width of each smooth J interval, kyr; intervals do not overlap.
J_CENTER_FRACTIONS = [0.5, 0.2, 0.8, 0.35, 0.65, 0.05, 0.95]  # Add centers in this order.
ACTIVE_AGE_SD_TOLERANCE = 1e-12  # kyr; omit numerical-zero tails from active-age diagnostics.
SUPPORT_THRESHOLD = 0.90     # Reporting rule for q; not a calibrated probability of correctness.
MIN_VALID_FRACTION = 0.90    # At least 9/10 usable paired member decisions for a summary decision.

# Frozen estimator settings, matching 02; no tuning at higher uncertainty.
MAX_LAG = 5                 # PCMCI+/LKIF history, 20-year steps: up to 100 years.
GC_ORDER = 1                # Fixed Granger order, one 20-year step.
ALPHA = 0.05                # Member-level significance threshold; wrappers retain their Holm rules.
GPDC_NULL_SAMPLES = 10000    # GPDC native reference-null draws.
GPDC_OPTIMIZE_KERNEL = False  # Keep the fixed native kernel used by the exact-age reference.
CCM_E = 3                   # Embedding dimension; three coordinates.
CCM_TAU = -1                # Backward embedding interval: 20 years, 40-year total history.
CCM_EXCLUSION_RADIUS = 10   # Exclude temporal neighbors within 200 years.
CCM_LIBRARY_FRACTION = 0.80  # Fraction of usable embedded points in each random library.
CCM_SAMPLES = 10            # Libraries averaged for the original Tp=0 CCM test.
CCM_SURROGATES = 100         # Random cut-and-swap null records per direction.
CCM_LAG_STEPS = list(range(-15, 16))  # Fixed -300 to +300-year scan; boundary peaks remain unresolved.
CCM_LAG_LIBRARIES = 100      # Matched libraries for the observed lag curves, separate from surrogates.
CCM_LAG_SIGN_THRESHOLD = 0.90  # Peak-sign stability required to retain a significant direction.

ROOT = Path(__file__).resolve().parent
STEM = Path(__file__).stem
SOURCE_DIR = ROOT / "data" / "processed" / "02_check_causal_pipeline"
DATA_DIR = ROOT / "data" / "processed" / STEM
FIGURE_DIR = ROOT / "figures" / STEM
CACHE_DIR = ROOT / "data" / "cache" / "gpdc"
METHODS = ["PCMCI+", "LKIF", "CCM", "GC"]
PROCESSES = ["linear", "nonlinear", "independent"]
LEVELS = [1, 2, 3, 4]          # X and Y always share one level index; no cross-level pairing.
AGE_KINDS = ["C", "L", "E", "J"]
PAIR_NAMES = [f"{x}-{y}" for x in AGE_KINDS for y in AGE_KINDS]


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
    ignored = {"RESUME", "N_WORKERS", "RUN_BASELINES", "RUN_AGE_ENSEMBLES", "RUN_CAUSAL"}
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


def gpdc_null_path():
    sample_size = int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / DT)) + 1 - 2 * MAX_LAG
    return CACHE_DIR / f"null_T{sample_size}_n{GPDC_NULL_SAMPLES}_seed{SEED + 300}.npz"


def prepare_gpdc_cache():
    """Prepare the shared native null once before creating workers."""
    sample_size = int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / DT)) + 1 - 2 * MAX_LAG
    with threadpool_limits(limits=1):
        prepare_gpdc_null(gpdc_null_path(), sample_sizes=[sample_size],
                          sig_samples=GPDC_NULL_SAMPLES, seed=SEED + 300)


def make_scenarios():
    """Combine all sixteen ordered pairs within each of four shared levels."""
    if LEVELS != [1, 2, 3, 4] or len(AMPLITUDES_YEARS) != 4 or len(J_PEAK_COUNTS) != 4:
        raise ValueError("Provide four matching amplitude and J-count levels.")
    rows = []
    for level in LEVELS:
        for pair in PAIR_NAMES:
            kind_x, kind_y = pair.split("-")
            row = {"level": level, "pair": pair, "scenario": f"level{level}_{pair}",
                   "x_kind": kind_x, "y_kind": kind_y, "true_lag_years": 1000 * TRUE_LAG,
                   "n_age_members": N_AGE_MEMBERS, "n_causal_members": N_CAUSAL_MEMBERS}
            for record, kind in [("x", kind_x), ("y", kind_y)]:
                count = J_PEAK_COUNTS[level - 1] if kind == "J" else 0
                row[f"{record}_peak_sd_years"] = J_PEAK_SD_YEARS if kind == "J" else AMPLITUDES_YEARS[level - 1]
                row[f"{record}_jump_count"] = count
                row[f"{record}_jump_width_kyr"] = J_INTERVAL_WIDTH if kind == "J" else 0.0
                row[f"{record}_jump_coverage"] = count * J_INTERVAL_WIDTH / (ANALYSIS_OLD - ANALYSIS_YOUNG)
            rows.append(row)
    return pd.DataFrame(rows)


def jump_intervals(count):
    """Add fixed-width peaks without moving or overlapping the earlier peaks."""
    if count not in J_PEAK_COUNTS or count > len(J_CENTER_FRACTIONS):
        raise ValueError("Unsupported J peak count.")
    span = ANALYSIS_OLD - ANALYSIS_YOUNG
    centers = ANALYSIS_YOUNG + span * np.asarray(J_CENTER_FRACTIONS[:count])
    intervals = [(float(center - J_INTERVAL_WIDTH / 2), float(center + J_INTERVAL_WIDTH / 2))
                 for center in centers]
    ordered = np.asarray(sorted(intervals))
    if ordered[0, 0] < ANALYSIS_YOUNG or ordered[-1, 1] > ANALYSIS_OLD:
        raise ValueError("J intervals must remain inside the analysis window.")
    if np.any(ordered[1:, 0] < ordered[:-1, 1]):
        raise ValueError("J intervals must not overlap.")
    return intervals


def check_reference():
    """Check the single completed level-2 source before importing its results."""
    saved = json.loads((SOURCE_DIR / "baseline_parameters.json").read_text())
    source_age = json.loads((SOURCE_DIR / "age_comparison_parameters.json").read_text())
    status = json.loads((SOURCE_DIR / "run_status.json").read_text())
    source_pairs = len(PAIR_NAMES) * len(PROCESSES) * N_CAUSAL_MEMBERS
    if status["state"] != "complete" or status["completed_member_pairs"] != source_pairs:
        raise ValueError(f"Run experiment02 to completion with the same member count first ({source_pairs} level-2 pairs).")
    if (source_age.get("analysis_revision") != "native_null_cache_and_ccm_object_reuse_v1"
            or source_age.get("execution_environment") != execution_environment()):
        raise ValueError("Source02 was generated on a different platform/implementation; regenerate it locally.")
    names = ["SEED", "DURATION", "DT", "TRUE_LAG", "ANALYSIS_OLD", "ANALYSIS_YOUNG",
             "MAX_LAG", "GC_ORDER", "ALPHA", "GPDC_NULL_SAMPLES", "GPDC_OPTIMIZE_KERNEL",
             "CCM_E", "CCM_TAU", "CCM_EXCLUSION_RADIUS", "CCM_LIBRARY_FRACTION",
             "CCM_SAMPLES", "CCM_SURROGATES", "CCM_LAG_STEPS", "CCM_LAG_LIBRARIES",
             "CCM_LAG_SIGN_THRESHOLD"]
    for name in names:
        if saved[name] != globals()[name] or source_age[name] != globals()[name]:
            raise ValueError(f"Reference setting {name} differs from this experiment.")
    fixed = {"N_REPLICATES": 1, "MEASUREMENT_NOISE_RATIO": 0, "GRID_STEP": DT,
             "COUPLING_GAIN": 4.0, "MEMORY_X": 0.5, "MEMORY_Y": 0.2,
             "YOUNG_SPACING": 0.0275, "OLD_SPACING": 0.0275, "SPACING_JITTER": 0.10}
    for name, value in fixed.items():
        if saved[name] != value or source_age[name] != value:
            raise ValueError(f"Source process/sampling setting {name} differs.")
    chronology = {"N_AGE_MEMBERS": N_AGE_MEMBERS, "N_CAUSAL_MEMBERS": N_CAUSAL_MEMBERS,
                  "UNCERTAINTY_LEVEL": 2, "MAX_AGE_DRAWS": MAX_AGE_DRAWS,
                  "AGE_CORRELATION_LENGTH": AGE_CORRELATION_LENGTH, "AGE_EXP_SHAPE": AGE_EXP_SHAPE,
                  "AGE_SCALE_LEVELS": [value / 1000 for value in AMPLITUDES_YEARS],
                  "JUMP_COUNT_LEVELS": J_PEAK_COUNTS, "JUMP_PEAK_SD": J_PEAK_SD_YEARS / 1000,
                  "JUMP_INTERVAL_WIDTH": J_INTERVAL_WIDTH, "JUMP_CENTER_FRACTIONS": J_CENTER_FRACTIONS,
                  "SUPPORT_THRESHOLD": SUPPORT_THRESHOLD, "MIN_VALID_FRACTION": MIN_VALID_FRACTION}
    for name, value in chronology.items():
        if source_age[name] != value:
            raise ValueError(f"Source chronology setting {name} differs.")
    expected_scenarios = make_scenarios().query("level == 2").set_index("pair")
    actual_scenarios = pd.read_csv(SOURCE_DIR / "scenarios.csv").set_index("scenario")
    if set(actual_scenarios.index) != set(PAIR_NAMES):
        raise ValueError("Experiment 02 does not contain all sixteen ordered pairs.")
    for column in ["x_kind", "y_kind", "x_peak_sd_years", "y_peak_sd_years",
                   "x_jump_count", "y_jump_count", "x_jump_width_kyr", "y_jump_width_kyr",
                   "x_jump_coverage", "y_jump_coverage", "true_lag_years",
                   "n_age_members", "n_causal_members"]:
        expected = expected_scenarios.loc[PAIR_NAMES, column].to_numpy()
        actual = actual_scenarios.loc[PAIR_NAMES, column].to_numpy()
        matches = np.array_equal(actual, expected) if column.endswith("kind") else np.allclose(actual, expected, rtol=0, atol=1e-12)
        if not matches:
            raise ValueError(f"Experiment 02 scenario column {column} differs.")
    for package, saved_version in source_age["packages"].items():
        if version(package) != saved_version:
            raise ValueError(f"Package {package} differs from the completed 02 run.")
    saved["age_settings"] = source_age
    return saved


def read_level2_results():
    """Read one completed result table; keep detailed files in experiment 02."""
    members = pd.read_csv(SOURCE_DIR / "member_results.csv",
                          dtype={"edge_detected": "boolean", "unoriented": "boolean"})
    keys = ["replicate", "process", "scenario", "member", "method", "source", "target"]
    expected = {(0, process, pair, member, method, source, target)
                for process in PROCESSES for pair in PAIR_NAMES
                for member in range(N_CAUSAL_MEMBERS) for method in METHODS
                for source, target in [("x", "y"), ("y", "x")]}
    actual = set(members[keys].itertuples(index=False, name=None))
    if members.duplicated(keys).any() or actual != expected:
        raise ValueError("Experiment 02 member coverage is incomplete or duplicated.")
    required = {"raw_edge_detected", "lag_status", "lag_n_requested_libraries", "unoriented"}
    if not required.issubset(members.columns):
        raise ValueError("Experiment 02 results lack the CCM lag-screen fields.")
    ccm = members[(members["method"] == "CCM") & (members["status"] == "ok")]
    if not ccm["lag_n_requested_libraries"].eq(CCM_LAG_LIBRARIES).all():
        raise ValueError("Experiment 02 did not request the frozen CCM lag screen.")
    members["level"] = 2
    members["pair"] = members["scenario"]
    members["scenario"] = "level2_" + members["pair"]
    members["source_experiment"] = SOURCE_DIR.name
    members["result_origin"] = str((SOURCE_DIR / "member_results.csv").relative_to(ROOT))
    return members


def prepare_inputs(scenarios):
    """Generate only levels 1/3/4, sharing clocks across pairs and process forms."""
    if not 1 <= N_CAUSAL_MEMBERS <= N_AGE_MEMBERS:
        raise ValueError("Causal members must be a fixed prefix of the chronology ensemble.")
    full = pd.read_csv(SOURCE_DIR / "r0_linear_generated_series.csv")
    age = full["age"].to_numpy()
    age_grid = np.linspace(ANALYSIS_OLD, ANALYSIS_YOUNG,
                           int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / DT)) + 1)
    observed, indices = {}, {}
    for process in PROCESSES:
        for record in ["x", "y"]:
            exact = pd.read_csv(SOURCE_DIR / f"r0_{process}_{record}_irregular_observed.csv")
            source_indices = np.rint((DURATION - exact["age_obs"].to_numpy()) / DT).astype(int)
            if not np.allclose(age[source_indices], exact["age_obs"], rtol=0, atol=1e-12):
                raise ValueError("Accurate-age observations must match the native grid.")
            if record in indices:
                if not np.array_equal(indices[record], source_indices):
                    raise ValueError("Process forms must retain the same sample ages.")
                if not np.array_equal(observed["linear", record]["sample_id"], exact["sample_id"]):
                    raise ValueError("Process forms must retain the same sample identities.")
            indices[record] = source_indices
            observed[process, record] = exact

    # Import the level-2 diagnostics without copying its canonical chronology files.
    age_source = pd.read_csv(SOURCE_DIR / "age_diagnostics.csv")
    age_source["level"] = 2
    age_source["record"] = age_source["record"].str[0]
    age_source["source_experiment"] = SOURCE_DIR.name
    age_rows = age_source.to_dict("records")
    relative_source = pd.read_csv(SOURCE_DIR / "relative_age_diagnostics.csv")
    relative_source["level"] = 2
    relative_source["pair"] = relative_source["scenario"]
    relative_source["scenario"] = "level2_" + relative_source["pair"]
    relative_rows = relative_source.to_dict("records")
    resampling_source = pd.read_csv(SOURCE_DIR / "resampling_diagnostics.csv")
    resampling_source["level"] = 2
    resampling_source["pair"] = resampling_source["scenario"]
    resampling_source["scenario"] = "level2_" + resampling_source["pair"]
    resampling_source["source_experiment"] = SOURCE_DIR.name
    resampling_rows = resampling_source.to_dict("records")
    source_inputs = pd.read_csv(SOURCE_DIR / "age_member_inputs.csv")
    source_inputs["level"] = 2
    source_inputs["pair"] = source_inputs["scenario"]
    source_inputs["scenario"] = "level2_" + source_inputs["pair"]
    source_inputs["source_experiment"] = SOURCE_DIR.name
    source_inputs["input_path"] = [str((SOURCE_DIR / f"{label}_series.csv").relative_to(ROOT))
                                   for label in source_inputs["label"]]
    all_inputs = source_inputs.to_dict("records")
    jobs = []

    for level in [value for value in LEVELS if value != 2]:
        clocks = {}
        for record, offset in [("x", 10), ("y", 20)]:
            for kind in AGE_KINDS:
                scale = J_PEAK_SD_YEARS if kind == "J" else AMPLITUDES_YEARS[level - 1]
                intervals = jump_intervals(J_PEAK_COUNTS[level - 1]) if kind == "J" else None
                sd = make_age_uncertainty(age[indices[record]], kind=kind, scale=scale / 1000,
                    age_range=(0, DURATION), exp_shape=AGE_EXP_SHAPE, jump_intervals=intervals)
                simulated = simulate_observations(age, full[record].to_numpy(), indices[record], sd,
                    age_correlation_length=AGE_CORRELATION_LENGTH, measurement_noise_ratio=0,
                    rng_age=SEED + offset, rng_measurement=SEED + offset + 1)
                if not np.array_equal(simulated["observed"]["sample_id"], observed["linear", record]["sample_id"]):
                    raise ValueError("The simulated clock changed sample identities.")
                ensemble = draw_age_ensemble(simulated["observed"], n_members=N_AGE_MEMBERS,
                    age_correlation_length=AGE_CORRELATION_LENGTH, max_draws=MAX_AGE_DRAWS,
                    rng=SEED + offset + 100)
                ages = ensemble["ages"]
                true_age = simulated["truth"]["age_true"].to_numpy()
                clocks[record, kind] = {"ages": ages, "true_age": true_age,
                                        "observed": simulated["observed"]}
                np.savez_compressed(DATA_DIR / f"r0_level{level}_{record}_{kind}_chronology.npz",
                    ages=ages, age_true=true_age, age_obs=simulated["observed"]["age_obs"],
                    age_std=sd, sample_id=ensemble["sample_id"])
                active = sd > ACTIVE_AGE_SD_TOLERANCE
                row = {"replicate": 0, "level": level, "record": record, "kind": kind,
                       "scale_years": scale, "n_samples": len(sd), "n_active_age_samples": int(active.sum()),
                       "source_experiment": STEM, **ensemble["diagnostics"],
                       "input_rms_sd_years": 1000 * float(np.sqrt(np.mean(sd**2))),
                       "median_sd_ratio": np.nan, "sample_fraction_inside_90_interval": np.nan,
                       "active_sample_fraction_inside_90_interval": np.nan,
                       "fraction_covering_analysis_window": np.nan}
                if len(ages):
                    lo, hi = np.percentile(ages, [5, 95], axis=0)
                    inside = (true_age >= lo) & (true_age <= hi)
                    row["sample_fraction_inside_90_interval"] = float(np.mean(inside))
                    if active.any():
                        row["active_sample_fraction_inside_90_interval"] = float(np.mean(inside[active]))
                        row["median_sd_ratio"] = float(np.median(np.std(ages[:, active], axis=0) / sd[active]))
                    row["fraction_covering_analysis_window"] = float(np.mean(
                        (ages[:, 0] >= ANALYSIS_OLD) & (ages[:, -1] <= ANALYSIS_YOUNG)))
                age_rows.append(row)
                print(f"Chronology level{level} {record}_{kind}: {len(ages)}/{N_AGE_MEMBERS} accepted", flush=True)

        for scenario in scenarios[scenarios["level"] == level].to_dict("records"):
            pair, name = scenario["pair"], scenario["scenario"]
            clock_x = clocks["x", scenario["x_kind"]]
            clock_y = clocks["y", scenario["y_kind"]]
            available = min(len(clock_x["ages"]), len(clock_y["ages"]))
            relative = []
            for member in range(available):
                error_x = np.interp(-age_grid, -clock_x["true_age"],
                                   clock_x["ages"][member] - clock_x["true_age"])
                error_y = np.interp(-age_grid, -clock_y["true_age"],
                                   clock_y["ages"][member] - clock_y["true_age"])
                relative.append(error_y - error_x)
            relative = np.asarray(relative)
            relative_rows.append({"replicate": 0, "level": level, "pair": pair, "scenario": name,
                "n_age_pairs": available,
                "relative_error_rms_years": 1000 * float(np.sqrt(np.mean(relative**2))) if available else np.nan,
                "fraction_relative_error_above_true_lag": float(np.mean(np.abs(relative) > TRUE_LAG)) if available else np.nan})
            for process in PROCESSES:
                for member in range(N_CAUSAL_MEMBERS):
                    label = f"r0_{process}_{name}_m{member}"
                    path = DATA_DIR / f"{label}_series.csv"
                    job = {"replicate": 0, "level": level, "pair": pair, "process": process,
                           "scenario": name, "member": member, "label": label,
                           "status": "ready", "error": "", "source_experiment": STEM,
                           "input_path": str(path.relative_to(ROOT))}
                    if member >= available:
                        job["status"], job["error"] = "chronology_error", "Requested chronology member was not accepted."
                    else:
                        try:
                            sampled = resample_pair(clock_x["ages"][member], observed[process, "x"]["value_obs"].to_numpy(),
                                clock_y["ages"][member], observed[process, "y"]["value_obs"].to_numpy(), age_grid=age_grid)
                            save_input(sampled["series"], path)
                            resampling_rows.append({**job, **sampled["diagnostics"]})
                        except ValueError as error:
                            job["status"], job["error"] = "resampling_error", str(error)
                    jobs.append(job)
                    all_inputs.append(job)
    pd.DataFrame(age_rows).to_csv(DATA_DIR / "age_diagnostics.csv", index=False)
    pd.DataFrame(relative_rows).to_csv(DATA_DIR / "relative_age_diagnostics.csv", index=False)
    pd.DataFrame(resampling_rows).to_csv(DATA_DIR / "resampling_diagnostics.csv", index=False)
    pd.DataFrame(all_inputs).to_csv(DATA_DIR / "age_member_inputs.csv", index=False)
    return jobs


def analyze_member(job):
    """Call native package wrappers and keep failed or reversed results visible."""
    label = job["label"]
    tables, diagnostics = [], []
    series = None
    input_error = job["error"]
    input_status = job["status"]
    if job["status"] == "ready":
        try:
            series = pd.read_csv(DATA_DIR / f"{label}_series.csv")
        except (OSError, ValueError) as error:
            input_error, input_status = str(error), "input_error"
    for method in METHODS:
        start = time.perf_counter()
        try:
            if series is None:
                raise ValueError(input_error)
            with threadpool_limits(limits=1):
                if method == "PCMCI+":
                    result = run_pcmciplus(series, tau_max=MAX_LAG, pc_alpha=ALPHA, alpha=ALPHA,
                        sig_samples=GPDC_NULL_SAMPLES, seed=SEED + 300,
                        gp_params=None if GPDC_OPTIMIZE_KERNEL else {"optimizer": None},
                        null_dist_filename=str(gpdc_null_path()))
                elif method == "LKIF":
                    result = run_lkif(series, max_lag=MAX_LAG, alpha=ALPHA)
                elif method == "CCM":
                    library_size = int(CCM_LIBRARY_FRACTION * (len(series) - (CCM_E - 1) * abs(CCM_TAU)))
                    result = run_ccm(series, E=CCM_E, tau=CCM_TAU,
                        exclusion_radius=CCM_EXCLUSION_RADIUS, lib_size=library_size,
                        sample=CCM_SAMPLES, n_surrogates=CCM_SURROGATES, seed=SEED + 300, alpha=ALPHA,
                        lag_steps=CCM_LAG_STEPS, n_libraries=CCM_LAG_LIBRARIES,
                        sign_threshold=CCM_LAG_SIGN_THRESHOLD)
                else:
                    result = run_granger(series, max_lag=GC_ORDER, alpha=ALPHA)
            edges = result["edges"].copy()
            diagnostic = result["diagnostics"]
            if "unoriented" not in edges:
                edges["unoriented"] = False
            if diagnostic.get("undirected_contemporaneous", False) or diagnostic.get("conflicting_contemporaneous", False):
                edges["unoriented"] = True
            arrays = {}
            for key, value in result["raw"].items():
                if isinstance(value, np.ndarray):
                    arrays[key] = value.astype(str) if value.dtype == object else value
                elif isinstance(value, pd.DataFrame):
                    value.to_csv(DATA_DIR / f"{label}_{method.replace('+', 'plus')}_{key}.csv", index=False)
            if arrays:
                np.savez_compressed(DATA_DIR / f"{label}_{method.replace('+', 'plus')}_raw.npz", **arrays)
        except Exception as error:
            edges = pd.DataFrame([
                {"source": source, "target": target, "score": np.nan, "p_value": np.nan,
                 "score_name": "unavailable", "p_kind": "unavailable", "edge_detected": None,
                 "unoriented": False, "status": "error" if series is not None else input_status, "error": str(error)}
                for source, target in [("x", "y"), ("y", "x")]])
            diagnostic = {"error": str(error)}
        edges["method"] = method
        edges["seconds"] = time.perf_counter() - start
        tables.append(edges)
        diagnostics.append({"method": method, **diagnostic})
        lag_message = f", lag={edges['lag_status'].tolist()}" if "lag_status" in edges else ""
        print(f"{label}: {method}, {edges['seconds'].iloc[0]:.1f} s, {edges['status'].tolist()}{lag_message}", flush=True)
    edges = pd.concat(tables, ignore_index=True)
    for column in ["replicate", "level", "pair", "process", "scenario", "member"]:
        edges[column] = job[column]
    edges["source_experiment"] = STEM
    edges["input_status"] = job["status"]
    edges["edge_detected"] = edges["edge_detected"].astype("boolean")
    edges["unoriented"] = edges["unoriented"].astype("boolean")
    path = DATA_DIR / f"{label}_edges.csv"
    edges["result_origin"] = str(path.relative_to(ROOT))
    temporary = path.with_suffix(".partial.csv")
    edges.to_csv(temporary, index=False)
    temporary.replace(path)
    (DATA_DIR / f"{label}_method_diagnostics.json").write_text(json.dumps(diagnostics, indent=2), encoding="utf-8")
    return edges


def summarize_results(members, baseline, scenarios):
    """Summarize complete requested groups; pending groups are never called failures."""
    tables = []
    keys = ["level", "pair", "process", "member", "method", "source", "target"]
    if members.duplicated(keys).any():
        raise ValueError("Duplicate member/method/direction results cannot be summarized.")
    expected = {(member, method, source, target) for member in range(N_CAUSAL_MEMBERS)
                for method in METHODS for source, target in [("x", "y"), ("y", "x")]}
    for (level, pair, process), group in members.groupby(["level", "pair", "process"], sort=False):
        actual = set(group[["member", "method", "source", "target"]].itertuples(index=False, name=None))
        if actual != expected:
            continue
        reference = baseline[baseline["process"] == process]
        summary = summarize_ensemble(group, n_members=N_CAUSAL_MEMBERS, baseline=reference,
            support_threshold=SUPPORT_THRESHOLD, min_valid_fraction=MIN_VALID_FRACTION)
        summary["replicate"] = 0
        summary["level"], summary["pair"], summary["process"] = level, pair, process
        summary["scenario"] = f"level{level}_{pair}"
        fraction = summary["fraction_neither"] if process == "independent" else summary["fraction_x_only"]
        summary["n_correct_graph"] = (fraction.fillna(0) * summary["n_valid"]).round().astype(int)
        summary["fraction_correct_requested"] = summary["n_correct_graph"] / N_CAUSAL_MEMBERS
        summary["n_unresolved"] = (summary["fraction_unoriented"].fillna(0) * summary["n_valid"]).round().astype(int)
        summary["n_failed"] = N_CAUSAL_MEMBERS - summary["n_valid"]
        tables.append(summary)
    if not tables:
        raise ValueError("No complete scenario/process groups are ready to summarize.")
    summary = pd.concat(tables, ignore_index=True)
    summary = summary.merge(scenarios.drop(columns=["level", "pair"]), on="scenario", validate="many_to_one")
    summary.to_csv(DATA_DIR / "ensemble_summary.csv", index=False)
    ccm = members[members["method"] == "CCM"]
    if "lag_reason" in ccm:
        counts = ccm.groupby(["level", "pair", "process", "source", "lag_status", "lag_reason"], dropna=False).size()
        counts.rename("n_members").reset_index().to_csv(DATA_DIR / "ccm_lag_reason_counts.csv", index=False)
    return summary


def save_progress(tables, state, current_level=None, error=""):
    """Atomically save the combined results and a small progress summary."""
    members = pd.concat(tables, ignore_index=True).sort_values(
        ["level", "pair", "process", "member", "method", "source"])
    path = DATA_DIR / "member_results.csv"
    temporary = path.with_suffix(".partial.csv")
    members.to_csv(temporary, index=False)
    temporary.replace(path)
    pair_counts = members[["level", "pair", "process", "member"]].drop_duplicates().groupby("level").size()
    per_level = len(PAIR_NAMES) * len(PROCESSES) * N_CAUSAL_MEMBERS
    completed = int(pair_counts.sum())
    reused = int(pair_counts.get(2, 0))
    status = {"state": state, "completed_pairs": completed, "requested_pairs": len(LEVELS) * per_level,
              "reused_pairs": reused, "new_completed_pairs": completed - reused,
              "current_level": current_level, "completed_levels": [int(level) for level, count in pair_counts.items() if count == per_level],
              "causal_members_per_setting": N_CAUSAL_MEMBERS, "generated_members_per_clock": N_AGE_MEMBERS,
              "n_workers": N_WORKERS, "unusable_direction_rows": int(members["edge_detected"].isna().sum()),
              "updated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "error": error}
    path = DATA_DIR / "run_status.json"
    temporary = path.with_suffix(".partial.json")
    temporary.write_text(json.dumps(status, indent=2), encoding="utf-8")
    temporary.replace(path)
    return members


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    scenarios = make_scenarios()
    reference = check_reference()
    parameters = {name: value for name, value in globals().items()
                  if name.isupper() and isinstance(value, (str, int, float, list))}
    parameters = json.loads(json.dumps(parameters))
    parameters["reference_settings"] = reference
    parameters["packages"] = reference["age_settings"]["packages"]
    parameters["execution_environment"] = execution_environment()
    parameters["analysis_revision"] = "native_null_cache_and_ccm_object_reuse_v1"
    parameters["interpretation"] = "One process realization per form; observation-centered chronologies. Within-record support is not calibrated accuracy."
    parameters["reuse"] = "Read completed level2 from local experiment02; resume complete same-configuration members in other levels."
    parameter_path = DATA_DIR / "parameters.json"
    if parameter_path.exists():
        check_resume_parameters(parameter_path, parameters)
        previous = json.loads(parameter_path.read_text())
        if previous.get("reference_settings", {}).get("process_truth") != reference.get("process_truth"):
            raise ValueError("The source process definition changed; preserve the old run before recomputing.")
        if RUN_CAUSAL and not RESUME:
            raise ValueError("Outputs already exist. Set RESUME=True or preserve them in another folder.")
    elif any(DATA_DIR.glob("r*_m*_edges.csv")):
        raise ValueError("Member files lack their run parameters; do not reuse them without provenance.")
    if not RUN_CAUSAL:
        if not parameter_path.exists():
            raise ValueError("No saved experiment is available to redraw.")
        plot_uncertainty_envelopes(scenarios)
        plot_results(pd.read_csv(DATA_DIR / "ensemble_summary.csv"), scenarios)
        return
    prepare_gpdc_cache()
    parameter_path.write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    scenarios.to_csv(DATA_DIR / "scenarios.csv", index=False)
    baseline = pd.read_csv(SOURCE_DIR / "baseline_irregular_accurate_results.csv")
    baseline = baseline[baseline["replicate"] == 0].copy()
    source_members = read_level2_results()
    source_pairs = source_members[["pair", "process", "member"]].drop_duplicates().shape[0]
    per_level = len(PAIR_NAMES) * len(PROCESSES) * N_CAUSAL_MEMBERS
    total = len(LEVELS) * per_level
    tables = [source_members]
    try:
        save_progress(tables, "preparing")
        plot_uncertainty_envelopes(scenarios)
        with threadpool_limits(limits=1):
            jobs = prepare_inputs(scenarios)
        pending = []
        for job in jobs:
            completed = read_completed_member(job)
            if completed is None:
                pending.append(job)
            else:
                tables.append(completed)
        resumed = len(tables) - 1
        members = save_progress(tables, "running")
        plot_results(summarize_results(members, baseline, scenarios), scenarios)
        print(f"Reused {source_pairs} level-2 pairs, resumed {resumed} local pairs; "
              f"computing {len(pending)} with {N_WORKERS} workers.", flush=True)
        completed_new = resumed
        with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
            for level in [value for value in LEVELS if value != 2]:
                save_progress(tables, "running", level)
                level_jobs = [job for job in pending if job["level"] == level]
                futures = [pool.submit(analyze_member, job) for job in level_jobs]
                for future in as_completed(futures):
                    tables.append(future.result())
                    completed_new = completed_new + 1
                    if completed_new % 10 == 0:
                        save_progress(tables, "running", level)
                        print(f"Completed {source_pairs + completed_new}/{total} input pairs; {completed_new} outside level2.", flush=True)
                members = save_progress(tables, "summarizing", level)
                plot_results(summarize_results(members, baseline, scenarios), scenarios)
                print(f"Level {level} complete; figures updated.", flush=True)
        save_progress(tables, "complete")
    except BaseException as error:
        save_progress(tables, "interrupted_or_failed", error=f"{type(error).__name__}: {error}")
        raise



def plot_uncertainty_envelopes(scenarios):
    """Show the requested SD envelopes, before interpreting any causal result."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    ages = np.linspace(0, DURATION, 1501)
    colors = {1: "#332288", 2: "#4477AA", 3: "#EE7733", 4: "#66AA55"}
    widths = {1: 1.1, 2: 2.2, 3: 3.6, 4: 5.0}
    titles = {"C": "C: constant SD", "L": "L: linear increase",
              "E": "E: fixed exponential shape", "J": "J: more fixed-height intervals"}
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    for ax, kind in zip(axes.flat, ["C", "L", "E", "J"]):
        selected = scenarios[scenarios["pair"] == f"{kind}-{kind}"].sort_values("level")
        if kind == "J":
            # Draw the widest strokes first so all levels remain visible on shared peaks.
            selected = selected.sort_values("level", ascending=False)
        for row in selected.to_dict("records"):
            level = int(row["level"])
            intervals = jump_intervals(int(row["x_jump_count"])) if kind == "J" else None
            age_sd = make_age_uncertainty(
                ages, kind=kind, scale=row["x_peak_sd_years"] / 1000,
                age_range=(0, DURATION), exp_shape=AGE_EXP_SHAPE,
                jump_intervals=intervals,
            )
            if kind == "J":
                peak_word = "peak" if row["x_jump_count"] == 1 else "peaks"
                label = f"Level {level}: {row['x_jump_count']:g} {peak_word} ({row['x_jump_coverage']:.0%})"
                width = widths[level]
            else:
                label = f"Level {level}: peak SD {row['x_peak_sd_years']:g} yr"
                width = 1.8
            ax.plot(ages, 1000 * age_sd, color=colors[level], lw=width, label=label)
        ax.set_xlim(DURATION, 0)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("Age (kyr BP)")
        ax.set_ylabel("Input age SD (years)")
        ax.set_title(titles[kind], fontsize=12)
        handles, labels = ax.get_legend_handles_labels()
        if kind == "J":
            handles, labels = handles[::-1], labels[::-1]
        ax.legend(handles, labels, fontsize=9, loc="best")
        ax.grid(alpha=0.15)
    fig.suptitle("Four chronology-uncertainty levels\n"
                 "C/L/E increase SD; J increases affected duration at fixed 10-year peak SD", fontsize=14)
    fig.text(0.5, 0.018,
             "J intervals are 2.8 kyr wide; coverage is their total width / 28 kyr. "
             "SD is not an error bound, and levels do not imply equal total uncertainty across types.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.055, 1, 0.92))
    for extension in ["png", "pdf"]:
        fig.savefig(FIGURE_DIR / f"00_uncertainty_envelopes.{extension}", dpi=180)
    plt.close(fig)


def plot_results(summaries, scenarios, *, note=""):
    """Compare complete scenario groups; leave uncalculated cells visibly pending."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    keys = ["level", "process", "pair", "method"]
    if summaries.duplicated(keys).any():
        raise ValueError("Plot summaries must contain one row per level, process, pair and method.")
    if not summaries["level"].isin(LEVELS).all() or not summaries["pair"].isin(PAIR_NAMES).all():
        raise ValueError("Plot summaries contain an unexpected uncertainty level or type pair.")

    expected = {(process, pair, method) for process in PROCESSES
                for pair in PAIR_NAMES for method in METHODS}
    complete_levels = []
    for level in LEVELS:
        subset = summaries[summaries["level"] == level]
        observed = set(zip(subset["process"], subset["pair"], subset["method"]))
        if observed == expected:
            complete_levels.append(level)
    if complete_levels == LEVELS:
        completion_label = "All four levels complete"
    else:
        available = ", ".join(str(level) for level in complete_levels) or "none"
        completion_label = f"Partial results — completed levels: {available}; uncalculated cells remain pending"

    level_labels = []
    for level in LEVELS:
        constant = scenarios[(scenarios["level"] == level) & (scenarios["pair"] == "C-C")].iloc[0]
        jitter = scenarios[(scenarios["level"] == level) & (scenarios["pair"] == "J-J")].iloc[0]
        peak_word = "peak" if jitter["x_jump_count"] == 1 else "peaks"
        level_labels.append(f"{level}: {constant['x_peak_sd_years']:g} yr / "
                            f"{jitter['x_jump_count']:g} {peak_word} ({jitter['x_jump_coverage']:.0%})")
    level_text = "Level: C/L/E peak SD / J coverage — " + "; ".join(level_labels)

    figures = [
        ("01_x_to_y_support", "q_x_to_y", "Retained X → Y support",
         "Fraction of usable members"),
        ("02_y_to_x_support", "q_y_to_x", "Retained Y → X support",
         "Fraction of usable members"),
        ("03_correct_member_graph", "fraction_correct_requested", "Complete true-graph recovery",
         "Fraction of requested members"),
    ]
    for filename, value_column, title, colorbar_label in figures:
        is_graph = value_column == "fraction_correct_requested"
        fig, axes = plt.subplots(len(PROCESSES), len(METHODS), figsize=(17, 16), squeeze=False)
        color_map = plt.get_cmap("Blues").copy()
        color_map.set_bad("white")
        for process_index, process in enumerate(PROCESSES):
            for method_index, method in enumerate(METHODS):
                ax = axes[process_index, method_index]
                selected = summaries[(summaries["process"] == process) & (summaries["method"] == method)]
                lookup = {(int(row["level"]), row["pair"]): row for row in selected.to_dict("records")}
                values = np.full((len(PAIR_NAMES), len(LEVELS)), np.nan)
                for pair_index, pair in enumerate(PAIR_NAMES):
                    for level_index, level in enumerate(LEVELS):
                        row = lookup.get((level, pair))
                        if row is not None:
                            values[pair_index, level_index] = row[value_column]
                picture = ax.imshow(np.ma.masked_invalid(values), cmap=color_map, vmin=0, vmax=1,
                                    interpolation="nearest", aspect="auto")
                for pair_index, pair in enumerate(PAIR_NAMES):
                    for level_index, level in enumerate(LEVELS):
                        row = lookup.get((level, pair))
                        if row is None:
                            ax.text(level_index, pair_index, "—", ha="center", va="center",
                                    fontsize=8, color="0.65")
                            continue
                        value = row[value_column]
                        n_valid = int(row["n_valid"])
                        n_requested = int(row["n_requested"])
                        color = "white" if np.isfinite(value) and value >= 0.60 else "black"
                        if is_graph:
                            label = f"{int(row['n_correct_graph'])}/{n_requested}"
                        elif np.isfinite(value):
                            label = f"{value:.2f}"
                            if value >= SUPPORT_THRESHOLD and n_valid / n_requested >= MIN_VALID_FRACTION:
                                label = label + "*"
                        else:
                            label = "NA"
                        ax.text(level_index, pair_index, label, ha="center", va="center",
                                fontsize=7.5, color=color)
                        if int(row["n_unresolved"]) > 0:
                            ax.scatter(level_index + 0.35, pair_index - 0.30, marker="^", s=13,
                                       color="#33AA66", edgecolor="black", linewidths=0.3, zorder=3)
                        if int(row["n_failed"]) > 0:
                            ax.scatter(level_index + 0.35, pair_index + 0.30, marker="s", s=11,
                                       color="#EE9933", edgecolor="black", linewidths=0.3, zorder=3)
                            if not is_graph:
                                ax.text(level_index - 0.44, pair_index + 0.32,
                                        f"{n_valid}/{n_requested}", ha="left", va="center",
                                        fontsize=5.5, color=color)
                ax.set_xticks(np.arange(len(LEVELS)), [str(level) for level in LEVELS], fontsize=9)
                ax.set_yticks(np.arange(len(PAIR_NAMES)), [pair.replace("-", "–") for pair in PAIR_NAMES], fontsize=8)
                ax.set_xticks(np.arange(-0.5, len(LEVELS), 1), minor=True)
                ax.set_yticks(np.arange(-0.5, len(PAIR_NAMES), 1), minor=True)
                ax.grid(which="minor", color="0.88", linewidth=0.35)
                ax.tick_params(which="minor", bottom=False, left=False)
                for boundary in [3.5, 7.5, 11.5]:
                    ax.axhline(boundary, color="0.5", linewidth=0.8)
                ax.set_xlabel("Uncertainty level", fontsize=9)
                if method_index == 0:
                    ax.set_ylabel("Chronology pair: X model – Y model", fontsize=9)
                truth = "True graph: no arrows" if process == "independent" else "True graph: X → Y only"
                ax.set_title(f"{process.capitalize()} — {method}\n{truth}", fontsize=10)

        colorbar_axis = fig.add_axes([0.925, 0.20, 0.015, 0.58])
        colorbar = fig.colorbar(picture, cax=colorbar_axis)
        colorbar.set_label(colorbar_label, fontsize=10)
        figure_title = f"{title}\n{completion_label}"
        if note:
            figure_title = f"{note}\n{figure_title}"
        fig.suptitle(figure_title, fontsize=14, y=0.99 if note else 0.98)
        fig.text(0.5, 0.935, level_text, ha="center", fontsize=9)
        fig.text(0.5, 0.918,
                 f"J peak SD stays 10 yr; coupled-model true lag = {1000 * TRUE_LAG:g} yr; "
                 f"{N_CAUSAL_MEMBERS} requested chronology members per scenario; CCM lag screening enabled.",
                 ha="center", fontsize=9)
        handles = [
            Line2D([], [], marker="^", color="none", markerfacecolor="#33AA66",
                   markeredgecolor="black", markersize=5,
                   label="At least one member has unresolved evidence"),
            Line2D([], [], marker="s", color="none", markerfacecolor="#EE9933",
                   markeredgecolor="black", markersize=5,
                   label="Some requested members are unusable or missing"),
        ]
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.057),
                   frameon=False, ncol=2, fontsize=9)
        if is_graph:
            explanation = ("Cells show correct / requested members. Unresolved and unusable members do not count as recovered graphs; "
                           "their separate counts are saved in the summary table.")
        else:
            explanation = ("Cells show support among usable paired fits. * marks support ≥90% with ≥90% usable members. "
                           "Small fractions show usable / requested when incomplete; NA means no usable decision.")
        fig.text(0.5, 0.043, explanation, ha="center", fontsize=8.5)
        fig.text(0.5, 0.027,
                 "Unresolved means an unoriented contemporaneous link (PCMCI+) or an inconclusive lag screen (CCM). "
                 "A dash means not calculated, not absence of an arrow.",
                 ha="center", fontsize=8.5)
        fig.text(0.5, 0.012,
                 "These are chronology-ensemble summaries for one record per process, not calibrated accuracy estimates across independent climate histories.",
                 ha="center", fontsize=8.5)
        fig.subplots_adjust(left=0.065, right=0.90, bottom=0.105, top=0.885, hspace=0.37, wspace=0.30)
        for extension in ["png", "pdf"]:
            fig.savefig(FIGURE_DIR / f"{filename}.{extension}", dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    main()
