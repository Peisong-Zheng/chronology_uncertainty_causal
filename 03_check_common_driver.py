"""Fixed-memory linear common driving: can chronology error defeat adjustment?

The common-driver-only system has Z -> X and Z -> Y; the second system adds
direct X -> Y. Historical h0/h1 storage keys are retained. Compare pairwise X/Y calls with
conditional X/Y/Z calls while perturbing only Z, only X/Y, or every clock.
See agent/03_common_driver_plan.md. Existing 02 outputs are not inputs or targets.
"""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
for thread_name in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"]:
    os.environ[thread_name] = "1"
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "data/cache/matplotlib"))

from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import ast
import gzip
import hashlib
from importlib.metadata import version
import json
import multiprocessing
import platform
import shutil
import time

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch
import numpy as np
import pandas as pd
from scipy.linalg import solve_discrete_lyapunov
from threadpoolctl import threadpool_limits

from toolbox.observation import sample_irregular_indices, simulate_observations, draw_age_ensemble
from toolbox.causal import prepare_gpdc_null, run_pcmciplus, run_lkif, run_granger, run_ccm
from toolbox.multivariate import run_pcmciplus_multivariate, run_lkif_multivariate, run_granger_multivariate
from toolbox.plotting import add_panel_labels, save_figure

RUN_CAUSAL = True             # False summarizes saved cases without fitting.
PREPARE_ONLY = False          # True generates inputs and stops before fitting.
RESUME = True                 # Retain complete cases including scientific failures.
N_WORKERS = 32                # Windows outer processes; one numerical thread each.
PROCESS_SEEDS = [20261009, 20261301, 20261302]  # Three independent process/clock blocks.
STRUCTURES = ["h0_common_driver", "h1_direct"]
MEMORY_Z_YEARS = 500.0        # Autonomous driver's intrinsic decay time, years.
MEMORY_X_YEARS = 200.0        # X intrinsic response memory, years.
MEMORY_Y_YEARS = 200.0        # Y intrinsic response memory, years.
COMMON_GAIN = 2.0            # Long-run response to a unit Z, fixed background units.
DIRECT_GAIN = 1.0            # Long-run response to a unit X; enabled in the direct-link system.
DURATION = 30.0              # Retained process duration, kyr.
DT = 0.020                   # Native/analysis step and direct edge delay, kyr.
BURN_IN_STEPS = 250          # Fixed 5-kyr warm-up, ten longest memory times.
ANALYSIS_OLD = 29.0          # Old analysis boundary, kyr BP.
ANALYSIS_YOUNG = 1.0         # Young analysis boundary, kyr BP.
SAMPLING_SEED = 20261110     # One synchronous mask, fixed across variables/blocks.
SAMPLING_GAP = 0.0275        # Pre-rounding target gap, kyr; native gaps 20/40 years.
SAMPLING_JITTER = 0.10       # Fractional gap jitter before native-grid rounding.
N_AGE_MEMBERS = 10           # Conditional chronology members per record/block.
AGE_PROPOSAL_SD_YEARS = 40.0 # C-envelope proposal SD, not realized chronology error.
AGE_CORRELATION_LENGTH = 1 / 60  # Fixed normalized-sample-order covariance length.
AGE_MAX_DRAWS = 20000        # Whole-curve proposal budget; no age sorting/repair.
METHOD_SEED_OFFSET = 300     # Paired estimator seeds within each process/clock block.
CONDITIONS = ["original_uniform", "irregular_accurate", "z_only", "xy_only", "all_uncertain"]
METHODS = ["PCMCI+", "LKIF", "GC", "CCM"]
ANALYSES = {"pairwise": METHODS, "conditional": ["PCMCI+", "LKIF", "GC"]}
MAX_LAG = 5                 # PCMCI+/LKIF maximum lag in 20-year steps.
GC_ORDER = 1                # Fixed GC order; no order selection by chronology.
ALPHA = 0.05                # Native significance threshold and PCMCI pc_alpha.
GPDC_NULL_SAMPLES = 10000    # Native GPDC null samples at T=1391.
GPDC_GP_PARAMS = {"optimizer": None}  # Existing native CPU GPDC kernel settings.
CCM_E = 3                   # Existing embedding dimension.
CCM_TAU = -1                # Backward embedding interval, grid steps.
CCM_EXCLUSION_RADIUS = 10   # Temporal neighbor exclusion, grid steps.
CCM_LIBRARY_SIZE = 300      # States per original, surrogate and lag library.
CCM_SAMPLES = 10            # Original Tp=0 libraries.
CCM_SURROGATES = 100        # Native swap-null repetitions per direction.
CCM_LAG_STEPS = list(range(-15, 16))  # +/-300-year target shifts.
CCM_LAG_LIBRARIES = 100     # Paired library draws across all target shifts.
CCM_LAG_SIGN_THRESHOLD = 0.90  # Existing diagnostic sign rule.
SUPPORT_THRESHOLD = 0.90    # Per-history reporting threshold, not causal probability.
MIN_VALID_FRACTION = 0.90   # Minimum jointly usable requested X/Y decisions.
IMPLEMENTATION_REVISION = "fixed_linear_common_driver_cc40_v1"
STEM = Path(__file__).stem
DATA_DIR = ROOT / "data/processed" / STEM
FIGURE_DIR = ROOT / "figures" / STEM
CACHE_DIR = ROOT / "data/cache/gpdc"


def calculation_digest():
    """Exclude presentation and permissible run controls from the resume guard."""
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    controls = {"RUN_CAUSAL", "PREPARE_ONLY", "RESUME", "N_WORKERS"}
    retained = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id in controls for t in node.targets):
            continue
        if isinstance(node, ast.FunctionDef) and (node.name.startswith("_plot_") or node.name in {"plot_results", "_save_common_figure"}):
            continue
        if isinstance(node, ast.ImportFrom) and ((node.module or "").startswith("matplotlib") or node.module == "toolbox.plotting"):
            continue
        if isinstance(node, ast.Import) and all(a.name.startswith("matplotlib") for a in node.names):
            continue
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
            continue
        retained.append(node)
    tree.body = retained
    return hashlib.sha256(ast.dump(tree, include_attributes=False).encode()).hexdigest()


def execution_parameters():
    parameters = {name: value for name, value in globals().items()
                  if name.isupper() and isinstance(value, (str, int, float, list, dict))}
    parameters["packages"] = {name: version(name) for name in ["numpy", "scipy", "pandas", "scikit-learn",
        "tigramite", "pyEDM", "LK-Info-Flow", "statsmodels", "dcor", "threadpoolctl", "numba", "llvmlite"]}
    parameters["environment"] = {"system": platform.system(), "machine": platform.machine(),
        "python": platform.python_version(), "numpy_blas": np.__config__.CONFIG.get("Build Dependencies", {}).get("blas", {})}
    parameters["implementation_sha256"] = {name: file_digest(ROOT / name) for name in
        ["toolbox/observation.py", "toolbox/causal.py", "toolbox/multivariate.py"]}
    parameters["implementation_sha256"][Path(__file__).name] = calculation_digest()
    parameters["interpretation"] = "Three independent process/clock blocks; ten within-record chronology members; fixed memory."
    return json.loads(json.dumps(parameters))


def check_parameters(parameters, output_dir=DATA_DIR):
    output_dir = Path(output_dir)
    path = output_dir / "parameters.json"
    if path.exists():
        if RUN_CAUSAL and not RESUME:
            raise ValueError("Preserve existing output before starting a new experiment.")
        old = json.loads(path.read_text(encoding="utf-8"))
        controls = {"RUN_CAUSAL", "PREPARE_ONLY", "RESUME", "N_WORKERS"}
        changed = [key for key in set(old) | set(parameters) if key not in controls and old.get(key) != parameters.get(key)]
        if changed:
            raise ValueError(f"Saved scientific settings or implementation differ: {sorted(changed)}")
    elif any((output_dir / "members").glob("*_edges.csv")):
        raise ValueError("Committed results lack their scientific parameters.")


def generate_common_driver(process_seed, structure):
    """Explicit paired linear recurrence; column order in A/Q is Z, X, Y."""
    if structure not in STRUCTURES:
        raise ValueError("Unknown structure.")
    az, ax, ay = np.exp(-DT * 1000 / np.array([MEMORY_Z_YEARS, MEMORY_X_YEARS, MEMORY_Y_YEARS]))
    direct = DIRECT_GAIN if structure == "h1_direct" else 0.0
    transition = np.array([[az, 0, 0], [(1 - ax) * COMMON_GAIN, ax, 0],
                           [(1 - ay) * COMMON_GAIN, (1 - ay) * direct, ay]])
    innovation_sd = np.sqrt(1 - np.array([az, ax, ay])**2)
    covariance = solve_discrete_lyapunov(transition, np.diag(innovation_sd**2))
    n_steps = int(round(DURATION / DT)) + 1
    start = BURN_IN_STEPS + 1
    random = np.random.default_rng(process_seed)
    innovations = random.normal(size=(start + n_steps, 3))
    values = np.zeros_like(innovations)
    for step in range(1, len(values)):
        z, x, y = values[step - 1]
        values[step, 0] = az * z + innovation_sd[0] * innovations[step, 0]
        values[step, 1] = ax * x + (1 - ax) * COMMON_GAIN * z + innovation_sd[1] * innovations[step, 1]
        values[step, 2] = ay * y + (1 - ay) * (COMMON_GAIN * z + direct * x) + innovation_sd[2] * innovations[step, 2]
    values = values[start:]
    series = pd.DataFrame({"time": np.arange(n_steps) * DT, "age": DURATION - np.arange(n_steps) * DT,
                           "x": values[:, 1], "y": values[:, 2], "z": values[:, 0]})
    others = [0, 2]
    conditional_x_variance = covariance[1, 1] - covariance[1, others] @ np.linalg.solve(
        covariance[np.ix_(others, others)], covariance[others, 1])
    truth = {"structure": structure, "process_seed": process_seed, "order": ["z", "x", "y"],
        "transition": transition.tolist(), "innovation_sd": innovation_sd.tolist(),
        "stationary_covariance": covariance.tolist(), "eigenvalues": np.linalg.eigvals(transition).tolist(),
        "edges": [["z", "x"], ["z", "y"]] + ([["x", "y"]] if direct else []),
        "conditional_x_variance": float(conditional_x_variance),
        "conditional_direct_signal_ratio": float(((1 - ay) * direct)**2 * conditional_x_variance / innovation_sd[2]**2),
        "memory_years": [MEMORY_Z_YEARS, MEMORY_X_YEARS, MEMORY_Y_YEARS],
        "direct_delay_years": DT * 1000, "value_units": "fixed unforced background SD = 1"}
    return series, truth


def save_clock(clock, path):
    """Compare regenerated arrays before reusing an existing chronology file."""
    path = Path(path)
    if path.exists():
        with np.load(path, allow_pickle=False) as saved:
            if set(saved.files) != set(clock) or any(not np.array_equal(saved[k], v) for k, v in clock.items()):
                raise ValueError(f"Regenerated clock differs: {path.name}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".partial.npz")
    np.savez_compressed(temporary, **clock)
    replace_output(temporary, path)


def resample_three(full, indices, ages, grid):
    """Linear interpolation in sample order, with common support and no extrapolation."""
    result = pd.DataFrame({"time": grid[0] - grid, "age": grid})
    for record in ["x", "y", "z"]:
        age = np.asarray(ages[record])
        if not np.isfinite(age).all() or np.any(np.diff(age) >= 0):
            raise ValueError(f"{record}: unordered or nonfinite chronology.")
        if age[0] < grid[0] or age[-1] > grid[-1]:
            raise ValueError(f"{record}: chronology does not cover common window.")
        result[record] = np.interp(-grid, -age, full[record].to_numpy()[indices])
    return result


def case_name(seed, structure, condition, member=-1):
    return f"s{seed}_{structure}_{condition}_m{member}"


def prepare_inputs(output_dir=DATA_DIR):
    """Generate 192 inputs, with independent clock blocks and paired assignments."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    full_age = DURATION - np.arange(int(round(DURATION / DT)) + 1) * DT
    grid = np.linspace(ANALYSIS_OLD, ANALYSIS_YOUNG, int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / DT)) + 1)
    native_indices = np.rint((DURATION - grid) / DT).astype(int)
    indices = sample_irregular_indices(full_age, young_spacing=SAMPLING_GAP,
        old_spacing=SAMPLING_GAP, jitter=SAMPLING_JITTER, rng=SAMPLING_SEED)
    age_true = full_age[indices]
    save_input(pd.DataFrame({"sample_id": np.arange(len(indices)), "native_index": indices, "age_true": age_true}),
               output_dir / "sampling.csv")
    jobs, process_rows, clock_rows, acceptance_rows = [], [], [], []
    for seed in PROCESS_SEEDS:
        clocks = {}
        for number, record in enumerate(["x", "y", "z"]):
            age_seed = int(np.random.SeedSequence([seed, 3103, number]).generate_state(1)[0])
            ensemble_seed = int(np.random.SeedSequence([seed, 3104, number]).generate_state(1)[0])
            observation = simulate_observations(full_age, np.zeros(len(full_age)), indices,
                np.full(len(indices), AGE_PROPOSAL_SD_YEARS / 1000),
                age_correlation_length=AGE_CORRELATION_LENGTH, measurement_noise_ratio=0,
                rng_age=age_seed, rng_measurement=0)
            observed = observation["observed"]
            ensemble = draw_age_ensemble(observed, n_members=N_AGE_MEMBERS,
                age_correlation_length=AGE_CORRELATION_LENGTH, max_draws=AGE_MAX_DRAWS, rng=ensemble_seed)
            if not ensemble["diagnostics"]["complete"]:
                raise RuntimeError(f"Incomplete chronology generation: {seed}/{record}")
            clock = {"age_true": age_true, "age_obs": observed.age_obs.to_numpy(),
                "age_std": observed.age_std.to_numpy(), "ages": ensemble["ages"],
                "sample_id": ensemble["sample_id"], "native_indices": indices}
            save_clock(clock, output_dir / "clocks" / f"s{seed}_{record}_C.npz")
            clocks[record] = clock
            acceptance_rows.append({"process_seed": seed, "record": record, "age_seed": age_seed,
                "ensemble_seed": ensemble_seed, **ensemble["diagnostics"],
                "observed_rms_years": float(np.sqrt(np.mean((clock["age_obs"] - age_true)**2))) * 1000,
                "effective_member_sd_years": float(np.mean(np.std(clock["ages"], axis=0))) * 1000})
        for condition in CONDITIONS[1:]:
            members = [-1] if condition == "irregular_accurate" else range(N_AGE_MEMBERS)
            for member in members:
                errors = {}
                for record in ["x", "y", "z"]:
                    uncertain = condition == "all_uncertain" or (condition == "z_only" and record == "z") or (condition == "xy_only" and record != "z")
                    ages = clocks[record]["ages"][member] if uncertain else age_true
                    error = np.interp(-grid, -age_true, ages - age_true) * 1000
                    errors[record] = error
                    clock_rows.append({"process_seed": seed, "condition": condition, "member": member,
                        "record": record, "rms_error_years": float(np.sqrt(np.mean(error**2))),
                        "mean_error_years": float(np.mean(error)), "minimum_gap_years": float(np.min(-np.diff(ages))) * 1000,
                        "covers_analysis_window": bool(ages[0] >= grid[0] and ages[-1] <= grid[-1])})
                for row in clock_rows[-3:]:
                    for first, second in [("x", "y"), ("x", "z"), ("y", "z")]:
                        row[f"relative_{first}_{second}_rms_years"] = float(np.sqrt(np.mean((errors[first] - errors[second])**2)))
        for structure in STRUCTURES:
            full, truth = generate_common_driver(seed, structure)
            save_input(full, output_dir / "generated" / f"s{seed}_{structure}.csv")
            write_json(truth, output_dir / "generated" / f"s{seed}_{structure}_truth.json")
            for record in ["x", "y", "z"]:
                values = full[record].to_numpy()[native_indices]
                process_rows.append({"process_seed": seed, "structure": structure, "record": record,
                    "sd": float(np.std(values)), "lag1_correlation": float(np.corrcoef(values[:-1], values[1:])[0, 1]),
                    "conditional_direct_signal_ratio": truth["conditional_direct_signal_ratio"]})
            for condition in CONDITIONS:
                members = [-1] if condition in CONDITIONS[:2] else range(N_AGE_MEMBERS)
                for member in members:
                    identity = case_name(seed, structure, condition, member)
                    source = identity
                    if condition == "z_only":
                        source = case_name(seed, structure, "irregular_accurate")
                    elif condition == "all_uncertain":
                        source = case_name(seed, structure, "xy_only", member)
                    job = {"case_id": identity, "process_seed": seed, "structure": structure,
                        "condition": condition, "member": member, "estimator_seed": seed + METHOD_SEED_OFFSET,
                        "pairwise_source_case": source, "input_status": "ok"}
                    try:
                        if condition == "original_uniform":
                            series = full.iloc[native_indices].copy().reset_index(drop=True)
                            series["time"], series["age"] = grid[0] - grid, grid
                        else:
                            ages = {}
                            for record in ["x", "y", "z"]:
                                uncertain = condition == "all_uncertain" or (condition == "z_only" and record == "z") or (condition == "xy_only" and record != "z")
                                ages[record] = clocks[record]["ages"][member] if uncertain else age_true
                            series = resample_three(full, indices, ages, grid)
                    except ValueError:
                        job["input_status"] = "invalid_coverage_or_order"
                        series = pd.DataFrame({"time": grid[0] - grid, "age": grid, "x": np.nan, "y": np.nan, "z": np.nan})
                    save_input(series, output_dir / "inputs" / f"{identity}.csv")
                    jobs.append(job)
    write_table(pd.DataFrame(jobs), output_dir / "cases.csv")
    process_diagnostics, clock_diagnostics = pd.DataFrame(process_rows), pd.DataFrame(clock_rows)
    write_table(process_diagnostics, output_dir / "process_diagnostics.csv")
    write_table(clock_diagnostics, output_dir / "clock_diagnostics.csv")
    write_table(pd.DataFrame(acceptance_rows), output_dir / "chronology_acceptance.csv")
    return jobs, process_diagnostics, clock_diagnostics


def gpdc_null_path(estimator_seed):
    effective_n = int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / DT)) + 1 - 2 * MAX_LAG
    return CACHE_DIR / f"null_T{effective_n}_n{GPDC_NULL_SAMPLES}_seed{estimator_seed}.npz"


def call_method(series, analysis, method, estimator_seed, null_file):
    if analysis == "conditional":
        if method == "PCMCI+":
            return run_pcmciplus_multivariate(series, tau_max=MAX_LAG, pc_alpha=ALPHA, alpha=ALPHA,
                sig_samples=GPDC_NULL_SAMPLES, seed=estimator_seed, gp_params=GPDC_GP_PARAMS, null_dist_filename=str(null_file))
        if method == "LKIF":
            return run_lkif_multivariate(series, max_lag=MAX_LAG, alpha=ALPHA)
        if method == "GC":
            return run_granger_multivariate(series, order=GC_ORDER, alpha=ALPHA)
    if analysis == "pairwise":
        pair = series[["time", "age", "x", "y"]]
        if method == "PCMCI+":
            return run_pcmciplus(pair, tau_max=MAX_LAG, pc_alpha=ALPHA, alpha=ALPHA,
                sig_samples=GPDC_NULL_SAMPLES, seed=estimator_seed, gp_params=GPDC_GP_PARAMS, null_dist_filename=str(null_file))
        if method == "LKIF":
            return run_lkif(pair, max_lag=MAX_LAG, alpha=ALPHA)
        if method == "GC":
            return run_granger(pair, max_lag=GC_ORDER, alpha=ALPHA)
        if method == "CCM":
            return run_ccm(pair, E=CCM_E, tau=CCM_TAU, exclusion_radius=CCM_EXCLUSION_RADIUS,
                lib_size=CCM_LIBRARY_SIZE, sample=CCM_SAMPLES, n_surrogates=CCM_SURROGATES,
                seed=estimator_seed, alpha=ALPHA, lag_steps=CCM_LAG_STEPS,
                n_libraries=CCM_LAG_LIBRARIES, sign_threshold=CCM_LAG_SIGN_THRESHOLD)
    raise ValueError(f"Unknown analysis/method: {analysis}/{method}")


def edge_pairs(analysis):
    return [(a, b) for a in ["x", "y", "z"] for b in ["x", "y", "z"] if a != b] if analysis == "conditional" else [("x", "y"), ("y", "x")]


def empty_edges(analysis, method, status, error=""):
    return pd.DataFrame([{"analysis": analysis, "method": method, "source": source, "target": target,
        "score": np.nan, "p_value": np.nan, "edge_detected": None, "unoriented": False,
        "status": status, "error": error} for source, target in edge_pairs(analysis)])


def native_path(output_dir, case_id, analysis, method):
    return Path(output_dir) / "native" / f"{case_id}_{analysis}_{method.replace('+', 'plus')}.json.gz"


def run_case(job, output_dir, null_file):
    """One case commits 26 rows after seven complete native scientific archives."""
    output_dir = Path(output_dir)
    series = pd.read_csv(output_dir / "inputs" / f"{job['case_id']}.csv")
    source_case = job["pairwise_source_case"]
    reuse = source_case != job["case_id"] and job["input_status"] == "ok"
    source_rows = None
    if reuse:
        old_input = pd.read_csv(output_dir / "inputs" / f"{source_case}.csv")
        if not np.array_equal(series[["time", "age", "x", "y"]].to_numpy(), old_input[["time", "age", "x", "y"]].to_numpy()):
            raise ValueError("Pairwise reuse requested for different numerical inputs.")
        source_rows = pd.read_csv(output_dir / "members" / f"{source_case}_edges.csv", float_precision="round_trip",
                                 dtype={"edge_detected": "boolean", "unoriented": "boolean", "raw_edge_detected": "boolean"})
        if len(source_rows) != 26 or not source_rows.estimator_seed.eq(job["estimator_seed"]).all():
            raise ValueError("Pairwise source is incomplete or uses different estimator settings.")
    tables = []
    with threadpool_limits(limits=1):
        for analysis, methods in ANALYSES.items():
            for method in methods:
                started = time.perf_counter()
                destination = native_path(output_dir, job["case_id"], analysis, method)
                destination.parent.mkdir(parents=True, exist_ok=True)
                if reuse and analysis == "pairwise":
                    edges = source_rows[source_rows.analysis.eq(analysis) & source_rows.method.eq(method)].copy()
                    temporary = destination.with_suffix(".partial.gz")
                    shutil.copyfile(native_path(output_dir, source_case, analysis, method), temporary)
                    replace_output(temporary, destination)
                    if method == "CCM":
                        for suffix in ["profile", "peaks"]:
                            source = output_dir / "lag" / f"{source_case}_{suffix}.csv"
                            if source.exists():
                                target = output_dir / "lag" / f"{job['case_id']}_{suffix}.csv"
                                target.parent.mkdir(parents=True, exist_ok=True)
                                shutil.copyfile(source, target)
                else:
                    try:
                        if job["input_status"] != "ok":
                            raise ValueError(job["input_status"])
                        result = call_method(series, analysis, method, job["estimator_seed"], null_file)
                    except Exception as error:
                        edges = empty_edges(analysis, method, "error", f"{type(error).__name__}: {error}")
                        result = {"edges": edges, "raw": {}, "diagnostics": {"error": str(error)}}
                    else:
                        edges = result["edges"].copy()
                        if "unoriented" not in edges:
                            edges["unoriented"] = False
                        if analysis == "pairwise" and (result["diagnostics"].get("undirected_contemporaneous", False)
                                or result["diagnostics"].get("conflicting_contemporaneous", False)):
                            edges["unoriented"] = True
                    temporary = destination.with_suffix(".partial.gz")
                    with gzip.open(temporary, "wt", encoding="utf-8") as handle:
                        json.dump(native_snapshot(result), handle, separators=(",", ":"))
                    replace_output(temporary, destination)
                    if method == "CCM":
                        for key, suffix in [("lag_profile", "profile"), ("lag_peaks", "peaks")]:
                            if isinstance(result["raw"].get(key), pd.DataFrame):
                                write_table(result["raw"][key], output_dir / "lag" / f"{job['case_id']}_{suffix}.csv")
                edges["analysis"], edges["method"] = analysis, method
                edges["pairwise_reused"] = reuse and analysis == "pairwise"
                edges["seconds"] = 0.0 if reuse and analysis == "pairwise" else time.perf_counter() - started
                for key, value in job.items():
                    edges[key] = value
                tables.append(edges)
    result = pd.concat(tables, ignore_index=True)
    write_table(result, output_dir / "members" / f"{job['case_id']}_edges.csv")
    return result


def read_completed_case(job, output_dir=DATA_DIR):
    output_dir = Path(output_dir)
    path = output_dir / "members" / f"{job['case_id']}_edges.csv"
    if not path.exists():
        return None
    result = pd.read_csv(path, float_precision="round_trip", dtype={
        "edge_detected": "boolean", "unoriented": "boolean", "raw_edge_detected": "boolean"})
    keys = ["analysis", "method", "source", "target"]
    expected = {(analysis, method, source, target) for analysis, methods in ANALYSES.items()
                for method in methods for source, target in edge_pairs(analysis)}
    required = set(keys + ["status", "edge_detected", "unoriented", "score", "p_value", *job])
    if not required.issubset(result) or len(result) != len(expected) or result.duplicated(keys).any():
        return None
    if set(result[keys].itertuples(index=False, name=None)) != expected or result.status.isna().any() or result.status.eq("pending").any():
        return None
    for key, value in job.items():
        if not result[key].eq(value).all():
            raise ValueError(f"Saved case identity differs: {key}/{path.name}")
    for analysis, methods in ANALYSES.items():
        for method in methods:
            if not native_path(output_dir, job["case_id"], analysis, method).is_file():
                return None
    return result


def collect_results(jobs, completed):
    tables = []
    for job in jobs:
        result = completed.get(job["case_id"])
        if result is None:
            result = pd.concat([empty_edges(analysis, method, "pending")
                for analysis, methods in ANALYSES.items() for method in methods], ignore_index=True)
            for key, value in job.items():
                result[key] = value
        tables.append(result)
    combined = pd.concat(tables, ignore_index=True)
    for name in ["edge_detected", "unoriented", "raw_edge_detected"]:
        if name in combined:
            combined[name] = combined[name].astype("boolean")
    return combined


def summarize_results(results):
    """Save pair states, per-block q and matched within-analysis score changes."""
    states, networks = [], []
    keys = ["case_id", "structure", "process_seed", "condition", "member", "analysis", "method"]
    for identity, group in results.groupby(keys, sort=False):
        meta = dict(zip(keys, identity))
        pair = group[group.source.isin(["x", "y"]) & group.target.isin(["x", "y"])]
        usable = len(pair) == 2 and pair.status.eq("ok").all() and pair.edge_detected.notna().all()
        unresolved = bool(pair.unoriented.fillna(False).any())
        forward = bool(pair.loc[pair.source.eq("x"), "edge_detected"].iloc[0]) if usable else None
        reverse = bool(pair.loc[pair.source.eq("y"), "edge_detected"].iloc[0]) if usable else None
        state = "unavailable" if not usable else "unresolved" if unresolved else (
            "both" if forward and reverse else "x_only" if forward else "y_only" if reverse else "neither")
        correct = usable and not unresolved and (forward == (meta["structure"] == "h1_direct")) and not reverse
        states.append({**meta, "usable": usable, "unresolved": unresolved, "state": state,
            "forward": forward, "reverse": reverse, "correct_xy": bool(correct),
            "failed": bool(not usable and not pair.status.eq("pending").any()),
            "pending": bool(pair.status.eq("pending").any())})
        if meta["analysis"] == "conditional":
            truth = {("z", "x"), ("z", "y")} | ({("x", "y")} if meta["structure"] == "h1_direct" else set())
            complete = len(group) == 6 and group.status.eq("ok").all() and group.edge_detected.notna().all()
            detected = set(zip(group.loc[group.edge_detected.fillna(False), "source"], group.loc[group.edge_detected.fillna(False), "target"]))
            networks.append({**meta, "usable": complete, "unresolved": bool(group.unoriented.fillna(False).any()),
                "correct_graph": bool(complete and not group.unoriented.fillna(False).any() and detected == truth)})
    states = pd.DataFrame(states)
    support = []
    group_keys = ["structure", "process_seed", "condition", "analysis", "method"]
    for identity, rows in states.groupby(group_keys, sort=False):
        meta = dict(zip(group_keys, identity))
        valid = rows[rows.usable]
        q_forward = float(valid.forward.mean()) if len(valid) else np.nan
        q_reverse = float(valid.reverse.mean()) if len(valid) else np.nan
        eligible = len(valid) >= MIN_VALID_FRACTION * len(rows)
        support.append({**meta, "requested": len(rows), "usable": len(valid),
            "q_x_to_y": q_forward, "q_y_to_x": q_reverse, "support_eligible": eligible,
            "supported_x_to_y": bool(eligible and q_forward >= SUPPORT_THRESHOLD),
            "supported_y_to_x": bool(eligible and q_reverse >= SUPPORT_THRESHOLD),
            "correct_xy": int(rows.correct_xy.sum()), "unresolved": int(rows.unresolved.sum()),
            "failed": int(rows.failed.sum()), "pending": int(rows.pending.sum())})
    reference = results[results.condition.eq("irregular_accurate")]
    key = ["structure", "process_seed", "analysis", "method", "source", "target"]
    changes = results[results.condition.isin(CONDITIONS[2:])].merge(reference[key + ["score", "status"]],
        on=key, how="left", suffixes=("", "_reference"), validate="many_to_one")
    changes["usable_score_pair"] = changes.status.eq("ok") & changes.status_reference.eq("ok") & np.isfinite(changes.score) & np.isfinite(changes.score_reference)
    changes["score_change"] = np.where(changes.usable_score_pair, changes.score - changes.score_reference, np.nan)
    return {"xy_case_states": states, "direction_support": pd.DataFrame(support),
            "conditional_graph_states": pd.DataFrame(networks), "paired_score_changes": changes}


def write_status(state, completed, requested, started, output_dir=DATA_DIR, error=""):
    write_json({"state": state, "completed_cases": completed, "requested_cases": requested,
        "n_workers": N_WORKERS, "elapsed_seconds": time.perf_counter() - started,
        "updated_utc": datetime.now(timezone.utc).isoformat(), "error": error}, Path(output_dir) / "run_status.json")


def save_outputs(results, jobs, process_diagnostics, clock_diagnostics, output_dir=DATA_DIR, figure_dir=FIGURE_DIR):
    tables = summarize_results(results)
    write_table(results, Path(output_dir) / "member_results.csv")
    for name, table in tables.items():
        write_table(table, Path(output_dir) / f"{name}.csv")
    if "plot_results" in globals():
        plot_results(results, jobs, process_diagnostics, clock_diagnostics, Path(output_dir), Path(figure_dir))
    return tables


def main():
    parameters = execution_parameters()
    check_parameters(parameters)
    write_json(parameters, DATA_DIR / "parameters.json")
    started = time.perf_counter()
    completed, jobs = {}, []
    try:
        jobs, process_diagnostics, clock_diagnostics = prepare_inputs(DATA_DIR)
        for job in jobs:
            result = read_completed_case(job, DATA_DIR)
            if result is not None:
                completed[job["case_id"]] = result
        if PREPARE_ONLY:
            write_status("prepared", len(completed), len(jobs), started)
            print(f"Prepared {len(jobs)} common-driver inputs.", flush=True)
            return
        pending = [job for job in jobs if job["case_id"] not in completed]
        if RUN_CAUSAL and pending:
            with threadpool_limits(limits=1):
                for seed in sorted({job["estimator_seed"] for job in pending}):
                    prepare_gpdc_null(gpdc_null_path(seed), sample_sizes=[1391], sig_samples=GPDC_NULL_SAMPLES, seed=seed)
            write_status("running", len(completed), len(jobs), started)
            context = multiprocessing.get_context("spawn")
            # Source cases precede exact-input reuse; both stages use the full outer pool.
            for dependent in [False, True]:
                stage = [job for job in pending if (job["pairwise_source_case"] != job["case_id"]) == dependent]
                if not stage:
                    continue
                print(f"Starting {len(stage)} cases; exact-input pair reuse={dependent}; workers={N_WORKERS}.", flush=True)
                with ProcessPoolExecutor(max_workers=min(N_WORKERS, len(stage)), mp_context=context) as pool:
                    futures = [pool.submit(run_case, job, str(DATA_DIR), str(gpdc_null_path(job["estimator_seed"]))) for job in stage]
                    for future in as_completed(futures):
                        result = future.result()
                        completed[result.case_id.iloc[0]] = result
                        if len(completed) % 6 == 0 or len(completed) == len(jobs):
                            write_status("running", len(completed), len(jobs), started)
                            print(f"Completed {len(completed)}/{len(jobs)} common-driver cases.", flush=True)
        results = collect_results(jobs, completed)
        save_outputs(results, jobs, process_diagnostics, clock_diagnostics)
        state = "complete" if len(completed) == len(jobs) else "partial"
        write_status(state, len(completed), len(jobs), started)
        print(f"Common-driver experiment {state}: {len(completed)}/{len(jobs)} cases.", flush=True)
    except BaseException as error:
        if jobs:
            write_table(collect_results(jobs, completed), DATA_DIR / "member_results.csv")
        write_status("interrupted_or_failed", len(completed), len(jobs), started, error=f"{type(error).__name__}: {error}")
        raise


def replace_output(temporary, target, timeout=5.0):
    """Retry only atomic file replacement after transient Windows sharing errors."""
    deadline = time.monotonic() + timeout
    while True:
        try:
            Path(temporary).replace(target)
            return
        except PermissionError as error:
            remaining = deadline - time.monotonic()
            if getattr(error, "winerror", None) not in {5, 32, 33} or remaining <= 0:
                raise
            time.sleep(min(0.05, remaining))


def write_table(table, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".partial.csv")
    table.to_csv(temporary, index=False)
    replace_output(temporary, path)


def write_json(value, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".partial.json")
    temporary.write_text(json.dumps(value, indent=2), encoding="utf-8")
    replace_output(temporary, path)


def file_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_input(series, path):
    """Check regenerated numerical values; never overwrite a matching committed input."""
    path = Path(path)
    if path.exists():
        previous = pd.read_csv(path, float_precision="round_trip")
        if (list(previous.columns) != list(series.columns) or previous.shape != series.shape
                or not np.array_equal(previous.to_numpy(), series.to_numpy(), equal_nan=True)):
            raise ValueError(f"Regenerated input changed: {path.name}")
        return
    write_table(series, path)


def native_snapshot(value):
    """Preserve complete native arrays, tables, nested results and GC fit fields."""
    if value is pd.NA:
        return None
    if isinstance(value, pd.DataFrame):
        return {"columns": list(value.columns), "index": native_snapshot(value.index.to_numpy()),
                "data": {str(column): native_snapshot(value[column].to_numpy()) for column in value.columns}}
    if isinstance(value, np.ndarray):
        return {"dtype": str(value.dtype), "shape": list(value.shape), "values": native_snapshot(value.tolist())}
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        converted = {str(key): native_snapshot(item) for key, item in value.items()}
        if len(converted) != len(value):
            raise ValueError("Native result keys collide after text serialization.")
        return converted
    if isinstance(value, (list, tuple, set)):
        items = sorted(value, key=repr) if isinstance(value, set) else value
        return [native_snapshot(item) for item in items]
    if hasattr(value, "model") and hasattr(value, "params"):
        fields = ["params", "resid", "fittedvalues", "normalized_cov_params", "bse", "tvalues",
                  "pvalues", "ssr", "scale", "nobs", "df_model", "df_resid", "rsquared",
                  "rsquared_adj", "centered_tss", "uncentered_tss", "ess", "mse_model", "mse_resid",
                  "mse_total", "fvalue", "f_pvalue", "llf", "aic", "bic", "condition_number"]
        fit = {name: getattr(value, name) for name in fields}
        fit.update({"cov_params": value.cov_params(), "conf_int": value.conf_int(),
                    "endog": value.model.endog, "exog": value.model.exog, "design_rank": value.model.rank})
        return native_snapshot(fit)
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError(f"Unsupported native result type: {type(value).__name__}")


def _plot_03_bool(value):
    """Read saved decisions without converting missing decisions into False."""
    if value is None or pd.isna(value):
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if value in (0, 1):
        return bool(value)
    if isinstance(value, str) and value.lower() in ("true", "false"):
        return value.lower() == "true"
    raise ValueError(f"Unknown saved decision: {value!r}")


def _plot_03_conditions():
    return [("irregular_accurate", "Accurate"), ("z_only", "Z only"),
            ("xy_only", "X and Y"), ("all_uncertain", "All three")]


def _plot_03_structures():
    return [("h0_common_driver", "Common driver only"),
            ("h1_direct", "Common driver + direct X → Y")]


def _plot_03_rows(edges, case_id, method, analysis):
    selected = edges[(edges.case_id == case_id) & (edges.method == method)
                     & (edges.analysis == analysis)]
    selected = selected[((selected.source == "x") & (selected.target == "y"))
                        | ((selected.source == "y") & (selected.target == "x"))]
    if selected.duplicated(["source", "target"]).any():
        raise ValueError(f"Duplicate XY rows: {case_id}/{method}/{analysis}")
    return selected


def _plot_03_case_state(rows):
    """Classify whole XY decisions, retaining unresolved and unavailable states."""
    if len(rows) != 2 or set(zip(rows.source, rows.target)) != {("x", "y"), ("y", "x")}:
        return "unavailable"
    if not rows.status.eq("ok").all():
        return "unavailable"
    decisions = [_plot_03_bool(value) for value in rows.edge_detected]
    if any(value is None for value in decisions):
        return "unavailable"
    if "unoriented" in rows and any(_plot_03_bool(value) is True for value in rows.unoriented):
        return "unresolved"
    forward = _plot_03_bool(rows.loc[rows.source.eq("x"), "edge_detected"].iloc[0])
    reverse = _plot_03_bool(rows.loc[rows.source.eq("y"), "edge_detected"].iloc[0])
    return "both" if forward and reverse else "x_only" if forward else "y_only" if reverse else "neither"


def _plot_03_counts(edges, cases, method, analysis, raw=False):
    """Use all requested cases as the denominator, separately for each direction."""
    counts = {"requested": len(cases), "forward": 0, "reverse": 0,
              "unresolved": 0, "unavailable": 0}
    field = "raw_edge_detected" if raw else "edge_detected"
    for case_id in cases.case_id:
        rows = _plot_03_rows(edges, case_id, method, analysis)
        state = _plot_03_case_state(rows)
        if state in ("unresolved", "unavailable"):
            counts[state] = counts[state] + 1
        for source, key in [("x", "forward"), ("y", "reverse")]:
            row = rows[rows.source.eq(source)]
            if len(row) == 1 and row.status.iloc[0] == "ok" and field in row:
                if _plot_03_bool(row[field].iloc[0]) is True:
                    counts[key] = counts[key] + 1
    return counts


def _plot_03_percent_axis(ax):
    ax.set_ylim(-0.04, 1.13)
    ax.set_yticks(np.linspace(0, 1, 6), ["0%", "20%", "40%", "60%", "80%", "100%"])
    ax.grid(axis="y", color="0.90", lw=0.7)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=10)


def _plot_03_save(fig, axes, output_dir, stem, rect):
    add_panel_labels(axes, x=-0.095, y=1.035)
    fig.tight_layout(rect=rect, h_pad=3.0, w_pad=2.5)
    save_figure(fig, output_dir, stem)
    plt.close(fig)


def _plot_03_design(jobs, data_dir, output_dir):
    """Display the prescribed graphs and a fixed input diagnostic window."""
    cases = pd.DataFrame(jobs).copy()
    seed = int(cases.process_seed.min())
    data_dir = Path(data_dir)
    fig, axes = plt.subplots(3, 2, figsize=(18, 13),
                             gridspec_kw={"height_ratios": [1.05, 1.35, 1.25]})
    variables = [("z", "Z: common driver", "#555555"),
                 ("x", "X", "#0072B2"), ("y", "Y", "#D55E00")]
    selected_series = []
    for column, (structure, title) in enumerate(_plot_03_structures()):
        ax = axes[0, column]
        positions = {"Z": (0.5, 0.78), "X": (0.22, 0.24), "Y": (0.78, 0.24)}
        for source, target in [("Z", "X"), ("Z", "Y")]:
            ax.annotate("", positions[target], positions[source],
                        arrowprops={"arrowstyle": "-|>", "lw": 2, "color": "0.3",
                                    "shrinkA": 23, "shrinkB": 23})
        if structure == "h1_direct":
            ax.annotate("", (0.71, 0.24), (0.29, 0.24),
                        arrowprops={"arrowstyle": "-|>", "lw": 2, "color": "#0072B2",
                                    "shrinkA": 0, "shrinkB": 0})
            ax.text(0.5, 0.29, f"K = {DIRECT_GAIN:g}; {DT * 1000:g} yr", ha="center", fontsize=11, color="#0072B2")
        else:
            ax.text(0.5, 0.25, "No X–Y direct link", ha="center", fontsize=11, color="0.35")
        for name, position in positions.items():
            ax.add_patch(Circle(position, radius=0.065, facecolor="white", edgecolor="0.25", lw=1.5))
            ax.text(*position, name, ha="center", va="center", fontsize=16, fontweight="bold")
            memory = {"Z": MEMORY_Z_YEARS, "X": MEMORY_X_YEARS, "Y": MEMORY_Y_YEARS}[name]
            memory_y = position[1]+0.115 if name == "Z" else position[1]-0.125
            ax.text(position[0], memory_y, f"Memory = {memory:g} yr", ha="center", fontsize=10)
        ax.text(0.28, 0.56, f"K = {COMMON_GAIN:g}\n{DT * 1000:g} yr", ha="center", fontsize=11)
        ax.text(0.72, 0.56, f"K = {COMMON_GAIN:g}\n{DT * 1000:g} yr", ha="center", fontsize=11)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_axis_off()
        ax.set_title(title, fontsize=14)
        case = cases[(cases.structure == structure) & (cases.process_seed == seed)
                     & (cases.condition == "irregular_accurate")].iloc[0]
        series = pd.read_csv(data_dir / "inputs" / f"{case.case_id}.csv", float_precision="round_trip")
        selected_series.append(series)
        window = series.age.between(13, 17)
        for variable, label, color in variables:
            axes[1, column].plot(series.loc[window, "age"], series.loc[window, variable],
                                 color=color, lw=1.3, label=label)
        axes[1, column].set_xlim(17, 13)
        axes[1, column].set_xlabel("Age (kyr BP)")
        axes[1, column].set_ylabel("Value (unforced background SD units)")
        axes[1, column].set_title(f"Accurate-age inputs | process seed {seed}", fontsize=12)
        axes[1, column].grid(color="0.92", lw=0.6)
        axes[1, column].legend(frameon=False, ncol=3, fontsize=10, loc="upper right")
    common_values = np.concatenate([frame.loc[frame.age.between(13, 17), ["x", "y", "z"]].to_numpy().ravel()
                                    for frame in selected_series])
    margin = 0.09 * np.ptp(common_values)
    for ax in axes[1]:
        ax.set_ylim(common_values.min()-margin, common_values.max()+margin)

    accurate = selected_series[0]
    perturbed_case = cases[(cases.structure == "h0_common_driver") & (cases.process_seed == seed)
                           & (cases.condition == "z_only")].sort_values("member").iloc[0]
    perturbed = pd.read_csv(data_dir / "inputs" / f"{perturbed_case.case_id}.csv", float_precision="round_trip")
    np.testing.assert_array_equal(accurate.age.to_numpy(), perturbed.age.to_numpy())
    window = accurate.age.between(13, 17)
    difference = perturbed[["x", "y", "z"]] - accurate[["x", "y", "z"]]
    ax = axes[2, 0]
    for variable, label, color in variables:
        style = "-" if variable == "z" else "--" if variable == "x" else ":"
        ax.plot(accurate.loc[window, "age"], difference.loc[window, variable],
                color=color, lw=1.5, ls=style, label=label)
    ax.set_xlim(17, 13)
    ax.set_xlabel("Age (kyr BP)")
    ax.set_ylabel("Perturbed input − accurate input")
    ax.set_title("Only Z is dated incorrectly: X and Y stay identical", fontsize=12)
    ax.grid(color="0.92", lw=0.6)
    ax.legend(frameon=False, ncol=3, fontsize=10, loc="upper right")
    maximum = difference[["x", "y"]].abs().max()
    ax.text(0.02, 0.04, f"Full-window maximum |change|: X={maximum.x:.2g}, Y={maximum.y:.2g}",
            transform=ax.transAxes, fontsize=10, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8})
    ax = axes[2, 1]
    table = ax.table(cellText=[["Accurate", "Accurate", "Accurate"],
                              ["Perturbed", "Accurate", "Accurate"],
                              ["Accurate", "Perturbed", "Perturbed"],
                              ["Perturbed", "Perturbed", "Perturbed"]],
                     rowLabels=["Accurate ages", "Z only", "X and Y", "All three"],
                     colLabels=["Z", "X", "Y"], cellLoc="center", rowLoc="center",
                     bbox=[0.18, 0.28, 0.8, 0.63])
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    for (row, column), cell in table.get_celld().items():
        cell.set_edgecolor("0.8")
        if row == 0:
            cell.set_facecolor("0.9")
        elif cell.get_text().get_text() == "Perturbed":
            cell.set_facecolor("#FFF0CF")
    ax.text(0.5, 0.12, f"Independent C-type clocks; proposal SD = {AGE_PROPOSAL_SD_YEARS:g} yr\n"
            "Same clocks and process innovations paired across the two systems",
            transform=ax.transAxes, ha="center", va="center", fontsize=11)
    ax.set_axis_off()
    ax.set_title("Chronology contrasts", fontsize=12)
    fig.suptitle("0. Can uncertain dating undermine control of a common driver?", fontsize=18, y=0.99)
    fig.text(0.5, 0.015, "The graphs show prescribed cross-variable effects; self-memory is written at each node. "
             "Curves use a fixed 17–13 kyr BP illustration window, not a selected event.", ha="center", fontsize=10)
    _plot_03_save(fig, axes, output_dir, "00_design_and_inputs", (0.02, 0.04, 0.99, 0.955))


def _plot_03_conditioning(edges, jobs, output_dir):
    """Compare directional detection before and after adding observed Z."""
    cases = pd.DataFrame(jobs)
    fig, axes = plt.subplots(2, 3, figsize=(21, 12), sharey=True)
    directions = [("forward", "X → Y", "#0072B2", -0.19),
                  ("reverse", "Y → X", "#D55E00", 0.19)]
    for row, (structure, title) in enumerate(_plot_03_structures()):
        for column, method in enumerate(["PCMCI+", "LKIF", "GC"]):
            ax = axes[row, column]
            tick_labels = []
            for position, (condition, name) in enumerate(_plot_03_conditions()):
                selected = cases[(cases.structure == structure) & (cases.condition == condition)]
                pair = _plot_03_counts(edges, selected, method, "pairwise")
                adjusted = _plot_03_counts(edges, selected, method, "conditional")
                total = len(selected)
                for key, label, color, offset in directions:
                    xx = np.array([position+offset-0.065, position+offset+0.065])
                    values = np.array([pair[key], adjusted[key]], dtype=float) / total
                    ax.plot(xx, values, color=color, alpha=0.60, lw=1.5)
                    ax.scatter(xx[0], values[0], s=65, facecolor="white", edgecolor=color, lw=1.5, zorder=3)
                    ax.scatter(xx[1], values[1], s=65, facecolor=color, edgecolor="white", lw=0.6, zorder=4)
                tick_labels.append(f"{name}\nn={total}\nX→Y {pair['forward']}→{adjusted['forward']}\n"
                                   f"Y→X {pair['reverse']}→{adjusted['reverse']}\n"
                                   f"U {pair['unresolved']}/{adjusted['unresolved']}; "
                                   f"F {pair['unavailable']}/{adjusted['unavailable']}")
            _plot_03_percent_axis(ax)
            ax.set_xticks(range(4), tick_labels, fontsize=9.2)
            ax.set_xlim(-0.55, 3.55)
            ax.set_title(f"{title}\n{method}", fontsize=13)
            if column == 0:
                ax.set_ylabel("Detected direction / requested inputs", fontsize=11)
    handles = [Line2D([], [], color=color, marker="o", ls="none", label=label) for _, label, color, _ in directions]
    handles = handles + [Line2D([], [], color="0.3", marker="o", markerfacecolor="white", ls="none", label="P: pairwise X,Y"),
                         Line2D([], [], color="0.3", marker="o", ls="none", label="C: conditional X,Y,Z")]
    fig.suptitle("1. Does adding Z remove spurious links and preserve the direct effect?", fontsize=18, y=0.99)
    fig.legend(handles=handles, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 0.952), frameon=False, fontsize=11)
    fig.text(0.5, 0.027, "Tick counts read pairwise → conditional. Correct XY state: neither arrow without a direct link; X → Y only with it. "
             "U/F: unresolved / unavailable case counts, each shown P/C.", ha="center", fontsize=10)
    fig.text(0.5, 0.008, "Accurate = irregular-sampling accurate-age reference. Denominators include all requested inputs; "
             "3 histories × 10 members give descriptive fractions. Retained arrows may coexist with unresolved evidence.", ha="center", fontsize=10)
    _plot_03_save(fig, axes, output_dir, "01_effect_of_conditioning", (0.02, 0.065, 0.995, 0.915))


def _plot_03_graph_states(edges, jobs, output_dir):
    """Expose how whole XY decisions change, including unresolved evidence."""
    cases = pd.DataFrame(jobs)
    states = [("neither", "Neither direction", "#D9D9D9", None),
              ("x_only", "X → Y only", "#0072B2", None),
              ("y_only", "Y → X only", "#D55E00", None),
              ("both", "Both directions", "#8064A2", None),
              ("unresolved", "Unresolved evidence", "#E6C84F", "///"),
              ("unavailable", "Unavailable", "#8B2635", "xx")]
    fig, axes = plt.subplots(2, 3, figsize=(21, 11.5), sharey=True)
    for row, (structure, title) in enumerate(_plot_03_structures()):
        for column, method in enumerate(["PCMCI+", "LKIF", "GC"]):
            ax = axes[row, column]
            labels = []
            for position, (condition, name) in enumerate(_plot_03_conditions()):
                selected = cases[(cases.structure == structure) & (cases.condition == condition)]
                total = len(selected)
                labels.append(f"{name}\nn={total}")
                for analysis, offset, short in [("pairwise", -0.20, "P"), ("conditional", 0.20, "C")]:
                    counts = {key: 0 for key, _, _, _ in states}
                    for case_id in selected.case_id:
                        state = _plot_03_case_state(_plot_03_rows(edges, case_id, method, analysis))
                        counts[state] = counts[state] + 1
                    bottom = 0.0
                    for key, label, color, hatch in states:
                        count = counts[key]
                        height = count / total
                        ax.bar(position+offset, height, bottom=bottom, width=0.34,
                               color=color, edgecolor="white", lw=0.5, hatch=hatch)
                        if count:
                            y = bottom + height/2
                            text_color = "white" if key in ("x_only", "y_only", "both", "unavailable") else "black"
                            if height >= 0.065:
                                ax.text(position+offset, y, str(count), ha="center", va="center",
                                        fontsize=10, fontweight="bold", color=text_color)
                            else:
                                ax.annotate(str(count), (position+offset+0.17, y), xytext=(3, 0),
                                            textcoords="offset points", fontsize=8, va="center", ha="left",
                                            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.2})
                        bottom = bottom + height
                    ax.text(position+offset, -0.025, short, ha="center", va="top", fontsize=10, fontweight="bold")
            _plot_03_percent_axis(ax)
            ax.set_xticks(range(4), labels, fontsize=10)
            ax.tick_params(axis="x", pad=20, length=0)
            ax.set_xlim(-0.65, 3.65)
            truth = "Correct XY state: neither direction" if row == 0 else "Correct XY state: X → Y only"
            ax.set_title(f"{title} | {method}\n{truth}", fontsize=12)
            if column == 0:
                ax.set_ylabel("Share of requested inputs", fontsize=11)
    handles = [Patch(facecolor=color, edgecolor="0.4", hatch=hatch, label=label) for _, label, color, hatch in states]
    fig.suptitle("2. What replaces the correct X–Y decision?", fontsize=18, y=0.99)
    fig.legend(handles=handles, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 0.945), frameon=False, fontsize=11)
    fig.text(0.5, 0.051, "P = pairwise X,Y; C = conditional X,Y,Z. Segment labels are counts. "
             "All requested inputs remain in the denominator.", ha="center", fontsize=11)
    fig.text(0.5, 0.028, "Unresolved evidence takes precedence over arrow combinations in this state display; "
             "a retained arrow can coexist with it (Figure 1). Unavailable includes failed or pending decisions.", ha="center", fontsize=10)
    fig.text(0.5, 0.008, "This is the X–Y subgraph only. Conditional estimates for Z → X, Z → Y, X → Z and Y → Z remain in the saved results.",
             ha="center", fontsize=10)
    _plot_03_save(fig, axes, output_dir, "02_xy_decision_states", (0.02, 0.085, 0.995, 0.89))


def _plot_03_ccm(edges, jobs, output_dir):
    """Keep pairwise CCM evidence distinct from conditional direct-edge tests."""
    cases = pd.DataFrame(jobs)
    fig, axes = plt.subplots(2, 2, figsize=(18, 11.5))
    directions = [("forward", "x", "y", "X → Y", "#0072B2", -0.17),
                  ("reverse", "y", "x", "Y → X", "#D55E00", 0.17)]
    ccm = edges[(edges.method == "CCM") & (edges.analysis == "pairwise")]
    finite_scores = pd.to_numeric(ccm.loc[ccm.status.eq("ok"), "score"], errors="coerce")
    finite_scores = finite_scores[np.isfinite(finite_scores)]
    lower = min(-0.05, float(finite_scores.min())-0.05) if len(finite_scores) else -0.05
    upper = max(1.0, float(finite_scores.max())+0.05) if len(finite_scores) else 1.0
    for column, (structure, title) in enumerate(_plot_03_structures()):
        ax = axes[0, column]
        labels = []
        score_labels = []
        for position, (condition, name) in enumerate(_plot_03_conditions()):
            selected = cases[(cases.structure == structure) & (cases.condition == condition)]
            raw = _plot_03_counts(edges, selected, "CCM", "pairwise", raw=True)
            screened = _plot_03_counts(edges, selected, "CCM", "pairwise")
            total = len(selected)
            labels.append(f"{name}\nn={total}\nX→Y {raw['forward']}→{screened['forward']}\n"
                          f"Y→X {raw['reverse']}→{screened['reverse']}\n"
                          f"U {screened['unresolved']}; F/P {screened['unavailable']}")
            counts = []
            for key, source, target, label, color, offset in directions:
                xx = np.array([position+offset-0.06, position+offset+0.06])
                values = np.array([raw[key], screened[key]], dtype=float) / total
                ax.plot(xx, values, color=color, lw=1.5, alpha=0.6)
                ax.scatter(xx[0], values[0], s=65, facecolor="white", edgecolor=color, lw=1.5, zorder=3)
                ax.scatter(xx[1], values[1], s=65, color=color, edgecolor="white", lw=0.6, zorder=4)
                subset = ccm[(ccm.structure == structure) & (ccm.condition == condition)
                             & (ccm.source == source) & (ccm.target == target)].sort_values(["process_seed", "member"])
                scores = pd.to_numeric(subset.loc[subset.status.eq("ok"), "score"], errors="coerce").to_numpy(dtype=float)
                scores = scores[np.isfinite(scores)]
                counts.append(f"{len(scores)}/{total}")
                if len(scores):
                    jitter = np.random.default_rng(0).uniform(-0.06, 0.06, len(scores))
                    axes[1, column].scatter(position+offset+jitter, scores, s=20, color=color, alpha=0.35, linewidths=0)
                    # Empirical-CDF quantiles are invariant to the ten exact
                    # repeats of each XY baseline in the Z-only condition.
                    q25, q75 = np.percentile(scores, [25, 75], method="inverted_cdf")
                    median = np.median(scores)
                    axes[1, column].errorbar(position+offset, median, yerr=[[median-q25], [q75-median]],
                                             fmt="o", color=color, markersize=7, capsize=4, lw=2,
                                             markeredgecolor="white", markeredgewidth=0.7)
            score_labels.append(f"{name}\nX→Y {counts[0]}\nY→X {counts[1]}")
        _plot_03_percent_axis(ax)
        ax.set_xlim(-0.55, 3.55)
        ax.set_xticks(range(4), labels, fontsize=10)
        ax.set_ylabel("Detected direction / requested inputs")
        ax.set_title(f"{title}\nRaw significance → lag-screened retention", fontsize=13)
        score_ax = axes[1, column]
        score_ax.set_xlim(-0.55, 3.55)
        score_ax.set_ylim(lower, upper)
        score_ax.set_xticks(range(4), score_labels, fontsize=10)
        score_ax.set_ylabel("Reconstruction correlation ρ (Tp = 0)")
        score_ax.set_title("Reconstruction skill before lag screening\n"
                           "Small dots: inputs; large dot: median; whiskers: middle 50%", fontsize=11.5)
        score_ax.grid(axis="y", color="0.9", lw=0.7)
        score_ax.set_axisbelow(True)
        score_ax.spines["top"].set_visible(False)
        score_ax.spines["right"].set_visible(False)
    handles = [Line2D([], [], color=color, marker="o", ls="none", label=label) for _, _, _, label, color, _ in directions]
    handles = handles + [Line2D([], [], color="0.3", marker="o", markerfacecolor="white", ls="none", label="Open: raw significant"),
                         Line2D([], [], color="0.3", marker="o", ls="none", label="Filled: lag-screened retained")]
    fig.suptitle("3. CCM is a pairwise diagnostic: Z is not conditioned on", fontsize=18, y=0.99)
    fig.legend(handles=handles, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 0.947), frameon=False, fontsize=11)
    fig.text(0.5, 0.069, "Upper tick counts read raw → screened; U = input with unresolved evidence; F/P = failed or pending decision. "
             "Raw significance is retained when a later lag screen fails.", ha="center", fontsize=10)
    fig.text(0.5, 0.049, "Lower panels: blue X → Y reconstructs X from Y's embedding; orange Y → X reconstructs Y from X's embedding.",
             ha="center", fontsize=10)
    fig.text(0.5, 0.029, "Lower counts are finite successful scores / requested inputs. Whiskers are the 25th–75th percentiles, not confidence intervals. "
             "A higher ρ means better reconstruction, not a larger causal effect.", ha="center", fontsize=10)
    fig.text(0.5, 0.009, "Changing only Z's chronology leaves X,Y inputs unchanged. Pairwise results should therefore agree with their accurate-age counterparts.",
             ha="center", fontsize=10)
    _plot_03_save(fig, axes, output_dir, "03_pairwise_ccm_evidence", (0.02, 0.105, 0.995, 0.91))


def _plot_03_conditioning_gain(edges, jobs, output_dir):
    """Plot the change in exact XY recovery after including observed Z."""
    cases = pd.DataFrame(jobs)
    systems = [("h0_common_driver", "Common driver only", "neither", "Correct X–Y state: neither direction"),
               ("h1_direct", "Common driver + direct X → Y", "x_only", "Correct X–Y state: X → Y only")]
    conditions = [("irregular_accurate", "Accurate\nages"), ("z_only", "Only Z\nuncertain"),
                  ("xy_only", "Only X,Y\nuncertain"), ("all_uncertain", "All three\nuncertain")]
    method_colors = [("PCMCI+", "#0072B2"), ("LKIF", "#009E73"), ("GC", "#D55E00")]
    figure_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(18, 8.5), sharey=True, squeeze=False)
    for column, (structure, system_name, correct_state, truth_label) in enumerate(systems):
        ax = axes[0, column]
        tick_labels = []
        for position, (condition, condition_label) in enumerate(conditions):
            selected = cases[(cases.structure == structure) & (cases.condition == condition)]
            requested = len(selected)
            if requested == 0:
                raise ValueError(f"No requested cases for {structure}/{condition}.")
            tick_labels.append(f"{condition_label}\nn={requested}")
            for method_index, (method, color) in enumerate(method_colors):
                counts = {}
                for analysis in ["pairwise", "conditional"]:
                    correct = 0
                    for case_id in selected.case_id:
                        state = _plot_03_case_state(_plot_03_rows(edges, case_id, method, analysis))
                        correct = correct + int(state == correct_state)
                    counts[analysis] = correct
                pairwise_percent = 100 * counts["pairwise"] / requested
                conditional_percent = 100 * counts["conditional"] / requested
                gain = conditional_percent - pairwise_percent
                figure_rows.append({"structure": structure, "condition": condition, "method": method,
                    "requested": requested, "pairwise_correct": counts["pairwise"],
                    "conditional_correct": counts["conditional"],
                    "pairwise_recovery_percent": pairwise_percent,
                    "conditional_recovery_percent": conditional_percent, "conditioning_gain_pp": gain})
                x = position + (method_index-1) * 0.24
                ax.bar(x, gain, width=0.22, color=color, edgecolor="white", linewidth=0.6, zorder=3)
                value_label = f"{gain:+.0f}" if np.isclose(gain, round(gain), rtol=0, atol=1e-10) else f"{gain:+.1f}"
                if abs(gain) >= 86:
                    text_y = gain - np.sign(gain) * 5
                    alignment = "top" if gain > 0 else "bottom"
                    text_color = "white"
                else:
                    text_y = gain + (4 if gain >= 0 else -4)
                    alignment = "bottom" if gain >= 0 else "top"
                    text_color = color
                ax.text(x, text_y, value_label, ha="center", va=alignment,
                        fontsize=10, fontweight="bold", color=text_color)
        ax.set_xticks(range(4), tick_labels, fontsize=11)
        ax.set_xlim(-0.55, 3.55)
        ax.set_ylim(-100, 100)
        ax.set_yticks(np.arange(-100, 101, 25))
        ax.tick_params(axis="y", labelsize=11)
        ax.axhline(0, color="0.2", linewidth=1.1, zorder=2)
        ax.grid(axis="y", color="0.90", linewidth=0.7)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_title(f"{system_name}\n{truth_label}", fontsize=13, pad=14)
        if column == 0:
            ax.set_ylabel("Conditional − pairwise exact X–Y recovery\n(percentage points)", fontsize=12)
    handles = [Patch(facecolor=color, label=method) for method, color in method_colors]
    fig.suptitle("4. How much does conditioning on Z improve the complete X–Y decision?", fontsize=17, y=0.99)
    fig.legend(handles=handles, ncol=3, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 0.943), fontsize=12)
    fig.text(0.5, 0.070, "Positive: more correct X–Y decisions after adding Z. Negative: fewer correct decisions. "
             "Labels are percentage-point changes, not changes in single-arrow detection.", ha="center", fontsize=10.5)
    fig.text(0.5, 0.043, "Correct states require both directions to match the prescribed graph with no unresolved evidence. "
             "Unavailable inputs remain in the requested denominator.", ha="center", fontsize=10.5)
    fig.text(0.5, 0.016, "Same requested inputs are compared within each condition. CCM has no conditional counterpart here. "
             "Fractions pool three histories; chronology members are not independent histories.", ha="center", fontsize=10)
    _plot_03_save(fig, axes, output_dir, "04_conditioning_gain", (0.02, 0.11, 0.995, 0.895))
    return pd.DataFrame(figure_rows)


def _plot_03_xy_recovery_comparison(edges, jobs, output_dir):
    """Compare method-wise exact XY recovery under the same requested conditions."""
    cases = pd.DataFrame(jobs)
    systems = [("h0_common_driver", "Common driver only", "neither", "Correct X–Y state: neither direction"),
               ("h1_direct", "Common driver + direct X → Y", "x_only", "Correct X–Y state: X → Y only")]
    conditions = [("irregular_accurate", "Accurate\nages"), ("z_only", "Only Z\nuncertain"),
                  ("xy_only", "Only X,Y\nuncertain"), ("all_uncertain", "All three\nuncertain")]
    method_styles = [("PCMCI+", "#0072B2", "o", "-"), ("LKIF", "#009E73", "s", "--"),
                     ("GC", "#D55E00", "^", "-."), ("CCM", "#8064A2", "D", ":")]
    figure_rows = []
    fig, axes = plt.subplots(2, 2, figsize=(20, 14), sharey=True, squeeze=False)
    for row, (structure, system_name, correct_state, truth_label) in enumerate(systems):
        for column, (analysis, analysis_label) in enumerate([("pairwise", "Pairwise X,Y"),
                                                          ("conditional", "Conditional X,Y,Z")]):
            ax = axes[row, column]
            styles = method_styles if analysis == "pairwise" else method_styles[:3]
            requested_counts = []
            for condition, _ in conditions:
                requested = len(cases[(cases.structure == structure) & (cases.condition == condition)])
                if requested == 0:
                    raise ValueError(f"No requested cases for {structure}/{condition}.")
                requested_counts.append(requested)
            count_table = []
            for method_index, (method, color, marker, linestyle) in enumerate(styles):
                rates = []
                labels = []
                for position, (condition, _) in enumerate(conditions):
                    selected = cases[(cases.structure == structure) & (cases.condition == condition)]
                    requested = requested_counts[position]
                    correct = 0
                    for case_id in selected.case_id:
                        state = _plot_03_case_state(_plot_03_rows(edges, case_id, method, analysis))
                        correct = correct + int(state == correct_state)
                    rate = 100 * correct / requested
                    rates.append(rate)
                    labels.append(f"{correct}/{requested}")
                    figure_rows.append({"structure": structure, "condition": condition, "analysis": analysis,
                        "method": method, "correct_xy": correct, "requested": requested,
                        "recovery_percent": rate})
                count_table.append(labels)
                # Tiny horizontal offsets expose coincident methods; all vertical
                # coordinates remain the exact requested-denominator percentages.
                offset = (method_index-(len(styles)-1)/2) * 0.045
                ax.plot(np.arange(4)+offset, rates, color=color, linestyle=linestyle, linewidth=1.8,
                        marker=marker, markersize=7.5, markeredgecolor="white", markeredgewidth=0.6,
                        label=method, zorder=3)
            ax.set_xticks(range(4), [f"{label}\nn={requested}" for (_, label), requested in zip(conditions, requested_counts)], fontsize=11)
            ax.set_xlim(-0.5, 3.5)
            ax.set_ylim(-3, 103)
            ax.set_yticks(np.arange(0, 101, 20), ["0%", "20%", "40%", "60%", "80%", "100%"])
            ax.tick_params(axis="y", labelsize=11)
            ax.grid(axis="y", color="0.90", linewidth=0.7)
            ax.set_axisbelow(True)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.set_title(f"{system_name}\n{analysis_label} | {truth_label}", fontsize=12, pad=14)
            if column == 0:
                ax.set_ylabel("Exact X–Y recovery (%)", fontsize=12)
            # Rows are methods; columns align with the four x-axis conditions.
            table_height = 0.075 * len(styles)
            table_top = -0.37
            table = ax.table(cellText=count_table, rowLabels=[style[0] for style in styles],
                             cellLoc="center", rowLoc="center",
                             bbox=[0.0, table_top-table_height, 1.0, table_height])
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            for (table_row, table_column), cell in table.get_celld().items():
                cell.set_edgecolor("0.84")
                cell.set_linewidth(0.6)
                if table_column == -1:
                    cell.get_text().set_color(styles[table_row][1])
                    cell.get_text().set_fontweight("bold")
                else:
                    cell.set_facecolor("0.97" if table_row % 2 == 0 else "white")
            ax.text(0.5, -0.31, "Correct / requested", transform=ax.transAxes,
                    ha="center", va="center", fontsize=10)
    handles = [Line2D([], [], color=color, marker=marker, linestyle=linestyle, linewidth=1.8,
                      markersize=7, label=method) for method, color, marker, linestyle in method_styles]
    fig.suptitle("5. How often is the complete X–Y decision correct?", fontsize=18, y=0.99)
    fig.legend(handles=handles, ncol=4, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 0.952), fontsize=12)
    fig.text(0.5, 0.058, "All panels use the same recovery scale and requested denominators. "
             "Unresolved or unavailable decisions are not counted as correct. Tables retain exact counts.", ha="center", fontsize=10.5)
    fig.text(0.5, 0.035, "Conditions are categorical: connecting lines guide the eye, not an increasing uncertainty scale. "
             "Method positions are slightly offset horizontally to expose overlapping values.", ha="center", fontsize=10.5)
    fig.text(0.5, 0.012, "CCM is evaluated as a pairwise workflow outcome; it does not condition on Z. "
             "Three histories and their chronology members give descriptive fractions, without confidence intervals.", ha="center", fontsize=10)
    # Explicit spacing reserves a separate band for the count tables. An
    # automatic tight layout would shrink the axes around these outside tables.
    for row in range(2):
        for column in range(2):
            axes[row, column].set_position([0.075 + 0.48 * column,
                                           0.65 - 0.395 * row, 0.425, 0.21])
    add_panel_labels(axes, x=-0.095, y=1.035)
    save_figure(fig, output_dir, "05_xy_recovery_comparison")
    plt.close(fig)
    return pd.DataFrame(figure_rows)


def plot_results(edges, jobs, process_diagnostics, clock_diagnostics, data_dir, output_dir):
    """Render six figures from saved inputs and decisions without causal refitting."""
    _plot_03_design(jobs, data_dir, output_dir)
    _plot_03_conditioning(edges, jobs, output_dir)
    _plot_03_graph_states(edges, jobs, output_dir)
    _plot_03_ccm(edges, jobs, output_dir)
    _plot_03_conditioning_gain(edges, jobs, output_dir)
    _plot_03_xy_recovery_comparison(edges, jobs, output_dir)


if __name__ == "__main__":
    main()
