"""Paired memory sensitivity of chronology effects for all four causal methods.

Four memory combinations, three process seeds and three process forms use two
accurate-age references and the same ten saved level-4 C-C chronologies.
See agent/02_8_memory_sensitivity_plan.md. Main-experiment files are read only.
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
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from toolbox.plotting import add_panel_labels
from toolbox.process import generate_process
from toolbox.causal import resample_pair, prepare_gpdc_null, run_pcmciplus, run_lkif, run_ccm, run_granger


RUN_CAUSAL = True             # False summarizes committed results without causal fits.
PREPARE_ONLY = False          # True prepares and validates inputs, then stops before fits.
RESUME = True                 # Keep all committed cases, including explicit scientific failures.
N_WORKERS = 32                # Measured Windows setting; one numerical thread per worker.
PROCESS_SEEDS = [20261009, 20261301, 20261302]  # Independent histories, paired across memory and form.
MEMORY_SETTINGS = [
    {"memory": "reference", "memory_x_years": 500.0, "memory_y_years": 200.0},
    {"memory": "short_x", "memory_x_years": 100.0, "memory_y_years": 200.0},
    {"memory": "short_y", "memory_x_years": 500.0, "memory_y_years": 40.0},
    {"memory": "short_xy", "memory_x_years": 100.0, "memory_y_years": 40.0},
]
DURATION = 30.0              # Retained process duration, kyr, including both endpoints.
DT = 0.020                   # Native and analysis spacing, kyr: 20 years.
TRUE_LAG = 0.020             # Direct X->Y equation delay, kyr: one step.
BACKGROUND_X = 1 / np.sqrt(1 - 0.95**2)  # Theoretical stationary driver SD.
BACKGROUND_Y = 1 / np.sqrt(1 - 0.8**2)   # Unforced response SD, not total coupled Y SD.
COUPLING_GAIN = 4.0          # Fixed long-run response to a held unit normalized driver.
BURN_IN_STEPS = 125          # Same 2.5-kyr warm-up across all conditions preserves innovation pairing.
ANALYSIS_OLD = 29.0          # Fixed analysis boundary, kyr BP.
ANALYSIS_YOUNG = 1.0         # Fixed analysis boundary, kyr BP; no extrapolation.
N_AGE_MEMBERS = 10           # Reused chronology members, not independent process histories.
SOURCE_LEVEL = 4             # Saved C-C clocks; proposal SD is 80 years per record.
AGE_PROPOSAL_SD_YEARS = 80.0 # Nominal proposal scale, not realized relative X/Y error.
AGE_CORRELATION_LENGTH = 1 / 60  # Original normalized sample-order error correlation, held fixed.
SOURCE_PROCESS_SEED = 20261009  # Identifies the fixed sampling/chronology source.
METHOD_SEED_OFFSET = 300     # A paired estimator random stream for each process seed.
SUPPORT_THRESHOLD = 0.90    # Reporting threshold, not a calibrated causal probability.
MIN_VALID_FRACTION = 0.90   # Minimum paired usable fraction for thresholded support.

MAX_LAG = 5                 # PCMCI+/LKIF maximum lag, 20-year steps.
GC_ORDER = 1                # Fixed statsmodels Granger order.
ALPHA = 0.05                # Existing native/wrapper significance thresholds.
GPDC_NULL_SAMPLES = 10000    # Native reference draws for effective sample size 1391.
GPDC_GP_PARAMS = {"optimizer": None}  # CPU GPDC with fixed kernel hyperparameters.
CCM_E = 3                   # Embedding dimension.
CCM_TAU = -1                # Backward embedding spacing, native grid steps.
CCM_EXCLUSION_RADIUS = 10   # Neighbor exclusion: 200 years, held fixed across memory.
CCM_LIBRARY_SIZE = 300      # Embedded states per original/null/lag reference library.
CCM_SAMPLES = 10            # Libraries averaged for the original Tp=0 score.
CCM_SURROGATES = 100        # Cut-and-swap null records per direction.
CCM_LAG_STEPS = list(range(-15, 16))  # Full +/-300-year target-shift scan.
CCM_LAG_LIBRARIES = 100     # Paired random libraries shared across target shifts.
CCM_LAG_SIGN_THRESHOLD = 0.90  # Existing lag-sign stability rule, not a probability.

STEM = Path(__file__).stem
DATA_DIR = ROOT / "data/processed" / STEM
FIGURE_DIR = ROOT / "figures" / STEM
CACHE_DIR = ROOT / "data/cache/gpdc"
BASELINE_DIR = ROOT / "data/processed/02_check_causal_pipeline"
CLOCK_DIR = ROOT / "data/processed/02_3_check_uncertainty_escalation"
PROCESSES = ["linear", "nonlinear", "independent"]
METHODS = ["PCMCI+", "LKIF", "CCM", "GC"]
IMPLEMENTATION_REVISION = "paired_memory_cc80_all_methods_v1"
IDENTITY_COLUMNS = ["case_id", "process_seed", "process", "memory", "memory_x_years",
                    "memory_y_years", "condition", "member", "estimator_seed", "input_status"]


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


def calculation_digest(path=None):
    """Ignore presentation and permitted run-control edits in the source guard."""
    path = Path(__file__) if path is None else Path(path)
    tree = ast.parse(path.read_text(encoding="utf-8-sig"))
    controls = {"RUN_CAUSAL", "PREPARE_ONLY", "RESUME", "N_WORKERS"}
    plotting_functions = {"plot_results", "_memory_plot_labels", "_save_memory_figure"}
    nodes = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id in controls
                                              for target in node.targets):
            continue
        if isinstance(node, ast.FunctionDef) and (node.name in plotting_functions or node.name.startswith("_plot_")):
            continue
        if isinstance(node, ast.ImportFrom) and (node.module == "toolbox.plotting" or (node.module or "").startswith("matplotlib")):
            continue
        if isinstance(node, ast.Import) and all(alias.name.startswith("matplotlib") for alias in node.names):
            continue
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
            continue
        nodes.append(node)
    tree.body = nodes
    return hashlib.sha256(ast.dump(tree, include_attributes=False).encode("utf-8")).hexdigest()


def source_files():
    return [BASELINE_DIR / "baseline_parameters.json", CLOCK_DIR / "parameters.json",
            CLOCK_DIR / "age_diagnostics.csv",
            *[CLOCK_DIR / f"r0_level{SOURCE_LEVEL}_{record}_C_chronology.npz" for record in ["x", "y"]]]


def execution_parameters():
    """Keep source-clock identity, numerical implementation and science in resume checks."""
    parameters = {name: value for name, value in globals().items()
                  if name.isupper() and isinstance(value, (str, int, float, list, dict))}
    parameters["packages"] = {name: version(name) for name in [
        "numpy", "scipy", "pandas", "scikit-learn", "tigramite", "pyEDM",
        "LK-Info-Flow", "statsmodels", "dcor", "threadpoolctl", "numba", "llvmlite"]}
    parameters["environment"] = {
        "system": platform.system(), "machine": platform.machine(), "python": platform.python_version(),
        "numpy_blas": np.__config__.CONFIG.get("Build Dependencies", {}).get("blas", {})}
    implementation = [ROOT / "toolbox/process.py", ROOT / "toolbox/causal.py"]
    parameters["source_sha256"] = {p.relative_to(ROOT).as_posix(): file_digest(p) for p in source_files()}
    parameters["implementation_sha256"] = {p.relative_to(ROOT).as_posix(): file_digest(p) for p in implementation}
    parameters["implementation_sha256"][Path(__file__).relative_to(ROOT).as_posix()] = calculation_digest()
    parameters["interpretation"] = "Three process histories; ten common clocks reused across histories, memories and forms."
    parameters["input_loading"] = "Default CSV parser for estimator inputs; round_trip only for committed result tables."
    return json.loads(json.dumps(parameters))


def check_parameters(parameters, output_dir=DATA_DIR):
    output_dir = Path(output_dir)
    path = output_dir / "parameters.json"
    if path.exists():
        if not RESUME and RUN_CAUSAL:
            raise ValueError("Outputs already exist; preserve them before a fresh experiment.")
        previous = json.loads(path.read_text(encoding="utf-8"))
        controls = {"RUN_CAUSAL", "PREPARE_ONLY", "RESUME", "N_WORKERS"}
        changed = [key for key in set(parameters) | set(previous)
                   if key not in controls and parameters.get(key) != previous.get(key)]
        if changed:
            raise ValueError(f"Saved experiment differs in {sorted(changed)}; preserve it before a fresh run.")
    elif any((output_dir / "members").glob("*_edges.csv")):
        raise ValueError("Committed results lack their scientific parameters.")


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


def load_clock_inputs():
    """Read saved ordered clocks without changing or regenerating their error fields."""
    baseline = json.loads((BASELINE_DIR / "baseline_parameters.json").read_text(encoding="utf-8"))
    source = json.loads((CLOCK_DIR / "parameters.json").read_text(encoding="utf-8"))
    for name in ["DURATION", "DT", "TRUE_LAG", "ANALYSIS_OLD", "ANALYSIS_YOUNG"]:
        if baseline[name] != globals()[name] or source[name] != globals()[name]:
            raise ValueError(f"Source chronology has a different {name}.")
    if baseline["SEED"] != SOURCE_PROCESS_SEED or source["SEED"] != SOURCE_PROCESS_SEED:
        raise ValueError("Source chronology seed differs from the specified fixed clocks.")
    if source["AGE_CORRELATION_LENGTH"] != AGE_CORRELATION_LENGTH:
        raise ValueError("Source chronology error correlation differs.")
    if source["AMPLITUDES_YEARS"][SOURCE_LEVEL - 1] != AGE_PROPOSAL_SD_YEARS:
        raise ValueError("Source proposal amplitude differs.")
    if source["execution_environment"]["system"] != platform.system():
        raise ValueError("Use authoritative local-platform source clocks.")
    clocks = {}
    n_steps = int(round(DURATION / DT)) + 1
    for record in ["x", "y"]:
        path = CLOCK_DIR / f"r0_level{SOURCE_LEVEL}_{record}_C_chronology.npz"
        with np.load(path, allow_pickle=False) as saved:
            clock = {name: saved[name].copy() for name in saved.files}
        if clock["ages"].shape != (N_AGE_MEMBERS, len(clock["age_true"])):
            raise ValueError("Saved clock count or sample support differs.")
        indices = np.rint((DURATION - clock["age_true"]) / DT).astype(int)
        if (np.any(indices < 0) or np.any(indices >= n_steps) or np.any(np.diff(indices) <= 0)
                or not np.allclose(DURATION - indices * DT, clock["age_true"], rtol=0, atol=1e-12)):
            raise ValueError("Clock truth does not identify ordered native sample locations.")
        for ages in [clock["age_true"], clock["age_obs"], *clock["ages"]]:
            if (not np.isfinite(ages).all() or np.any(np.diff(ages) >= 0)
                    or ages[0] < ANALYSIS_OLD or ages[-1] > ANALYSIS_YOUNG):
                raise ValueError("Every saved chronology must be ordered and cover the fixed grid.")
        if not np.allclose(clock["age_std"], AGE_PROPOSAL_SD_YEARS / 1000, rtol=0, atol=1e-15):
            raise ValueError("Saved clock is not the specified constant-SD chronology.")
        clock["native_indices"] = indices
        clocks[record] = clock
    return clocks


def prepare_inputs(output_dir=DATA_DIR):
    """Create 432 paired inputs; the same sampling and clocks apply to every history."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    clocks = load_clock_inputs()
    age_grid = np.linspace(ANALYSIS_OLD, ANALYSIS_YOUNG,
                           int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / DT)) + 1)
    full_age = np.linspace(DURATION, 0, int(round(DURATION / DT)) + 1)
    analysis_indices = np.rint((DURATION - age_grid) / DT).astype(int)
    if BURN_IN_STEPS * DT < 5 * max(max(s["memory_x_years"], s["memory_y_years"]) for s in MEMORY_SETTINGS) / 1000:
        raise ValueError("The fixed warm-up must cover at least five longest memory times.")
    source_diagnostics = pd.read_csv(CLOCK_DIR / "age_diagnostics.csv", float_precision="round_trip")
    selected = source_diagnostics[source_diagnostics.level.eq(SOURCE_LEVEL) & source_diagnostics.kind.eq("C")]
    if len(selected) != 2 or not selected.accepted_members.eq(N_AGE_MEMBERS).all():
        raise ValueError("Saved source diagnostics do not describe two complete C clocks.")
    write_table(selected, output_dir / "source_clock_diagnostics.csv")
    for record in ["x", "y"]:
        path = output_dir / "clocks" / f"{record}_C_level{SOURCE_LEVEL}.npz"
        source_path = CLOCK_DIR / f"r0_level{SOURCE_LEVEL}_{record}_C_chronology.npz"
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and file_digest(path) != file_digest(source_path):
            raise ValueError("Previously copied source clock changed.")
        if not path.exists():
            shutil.copyfile(source_path, path)

    clock_rows = []
    for member in range(N_AGE_MEMBERS):
        errors = {}
        for record in ["x", "y"]:
            clock = clocks[record]
            errors[record] = np.interp(-age_grid, -clock["age_true"],
                                       clock["ages"][member] - clock["age_true"]) * 1000
            clock_rows.append({"member": member, "record": record,
                               "proposal_sd_years": AGE_PROPOSAL_SD_YEARS,
                               "rms_error_years": float(np.sqrt(np.mean(errors[record]**2))),
                               "mean_error_years": float(np.mean(errors[record])),
                               "error_sd_years": float(np.std(errors[record])),
                               "minimum_gap_years": float(np.min(-np.diff(clock["ages"][member]))) * 1000,
                               "covers_analysis_window": True})
        relative = errors["y"] - errors["x"]
        for row in clock_rows[-2:]:
            row["relative_error_rms_years"] = float(np.sqrt(np.mean(relative**2)))
            row["fraction_relative_error_above_true_lag"] = float(np.mean(np.abs(relative) > TRUE_LAG * 1000))
    clock_diagnostics = pd.DataFrame(clock_rows)
    write_table(clock_diagnostics, output_dir / "clock_diagnostics.csv")

    jobs, process_rows, reference_rows = [], [], []
    for setting in MEMORY_SETTINGS:
        ar_x = np.exp(-DT / (setting["memory_x_years"] / 1000))
        ar_y = np.exp(-DT / (setting["memory_y_years"] / 1000))
        noise = (BACKGROUND_X * np.sqrt(1 - ar_x**2), BACKGROUND_Y * np.sqrt(1 - ar_y**2))
        for process_seed in PROCESS_SEEDS:
            # Equal lengths, initial draws and warm-up retain the same innovation positions.
            backgrounds = generate_process(n_steps=len(full_age), dt=DT, age_start=DURATION,
                lag=TRUE_LAG, ar_x=ar_x, ar_y=ar_y, coupling=0, coupling_form="linear",
                noise_std=noise, burn_in=BURN_IN_STEPS, rng=process_seed)["series"]
            for process in PROCESSES:
                gain = 0.0 if process == "independent" else COUPLING_GAIN
                generated = generate_process(n_steps=len(full_age), dt=DT, age_start=DURATION,
                    lag=TRUE_LAG, ar_x=ar_x, ar_y=ar_y, coupling=(1 - ar_y) * gain,
                    coupling_form="nonlinear" if process == "nonlinear" else "linear",
                    noise_std=noise, burn_in=BURN_IN_STEPS, rng=process_seed)
                full = generated["series"]
                prefix = f"s{process_seed}_{setting['memory']}_{process}"
                save_input(full, output_dir / "generated" / f"{prefix}.csv")
                write_json(generated["truth"], output_dir / "generated" / f"{prefix}_truth.json")
                original = full.iloc[analysis_indices].copy().reset_index(drop=True)
                original["time"], original["age"] = age_grid[0] - age_grid, age_grid
                observed, loaded_values = {}, {}
                for record in ["x", "y"]:
                    indices = clocks[record]["native_indices"]
                    observed[record] = pd.DataFrame({"sample_id": np.arange(len(indices), dtype=int),
                        "age_obs": full_age[indices], "value_obs": full[record].to_numpy()[indices],
                        "age_std": np.zeros(len(indices)), "measurement_std": np.zeros(len(indices))})
                    path = output_dir / "observations" / f"{prefix}_{record}.csv"
                    save_input(observed[record], path)
                    # Follow 02_3: chronology inputs attach values loaded from observed CSVs.
                    loaded_values[record] = pd.read_csv(path)["value_obs"].to_numpy()
                accurate = resample_pair(observed["x"].age_obs, observed["x"].value_obs,
                    observed["y"].age_obs, observed["y"].value_obs, age_grid=age_grid)["series"]
                series_by_condition = [("original_uniform", -1, original), ("irregular_accurate", -1, accurate)]
                for member in range(N_AGE_MEMBERS):
                    series = resample_pair(clocks["x"]["ages"][member], loaded_values["x"],
                        clocks["y"]["ages"][member], loaded_values["y"], age_grid=age_grid)["series"]
                    series_by_condition.append(("chronology", member, series))
                for condition, member, series in series_by_condition:
                    case_id = f"{prefix}_{condition}" + (f"_m{member:02d}" if member >= 0 else "")
                    job = {"case_id": case_id, "process_seed": process_seed, "process": process,
                           **setting, "condition": condition, "member": member,
                           "estimator_seed": process_seed + METHOD_SEED_OFFSET, "input_status": "ok"}
                    save_input(series, output_dir / "inputs" / f"{case_id}.csv")
                    jobs.append(job)
                    if setting["memory"] == "reference" and process_seed == SOURCE_PROCESS_SEED:
                        reference_path = (BASELINE_DIR / f"r0_{process}_{condition}_series.csv" if member < 0
                                          else CLOCK_DIR / f"r0_{process}_level4_C-C_m{member}_series.csv")
                        reference = pd.read_csv(reference_path)
                        prepared = pd.read_csv(output_dir / "inputs" / f"{case_id}.csv")
                        difference = np.max(np.abs(prepared.to_numpy() - reference.to_numpy()))
                        if difference > 1e-12:
                            raise ValueError(f"Reference dynamics do not reproduce the main input: {case_id}")
                        reference_rows.append({"case_id": case_id, "reference_path": reference_path.relative_to(ROOT).as_posix(),
                            "arrays_exactly_equal": bool(np.array_equal(prepared.to_numpy(), reference.to_numpy())),
                            "max_abs_difference": float(difference)})
                background_y = backgrounds.y.to_numpy()[analysis_indices]
                forced_y = original.y.to_numpy() - background_y
                process_rows.append({"process_seed": process_seed, "process": process, **setting,
                    **generated["truth"], "x_sd": float(np.std(original.x)), "y_sd": float(np.std(original.y)),
                    "x_lag1": float(np.corrcoef(original.x[:-1], original.x[1:])[0, 1]),
                    "y_lag1": float(np.corrcoef(original.y[:-1], original.y[1:])[0, 1]),
                    "background_y_sd": float(np.std(background_y)), "forced_y_sd": float(np.std(forced_y)),
                    "forced_to_background_sd": float(np.std(forced_y) / np.std(background_y)),
                    "x_innovation_sd": float(noise[0]), "y_innovation_sd": float(noise[1])})
    expected = len(MEMORY_SETTINGS) * len(PROCESS_SEEDS) * len(PROCESSES) * (2 + N_AGE_MEMBERS)
    if len(jobs) != expected or len({job["case_id"] for job in jobs}) != expected:
        raise ValueError("Memory input identities are incomplete or duplicated.")
    manifest = pd.DataFrame(jobs)
    manifest["input_sha256"] = [file_digest(output_dir / "inputs" / f"{job['case_id']}.csv") for job in jobs]
    prior = output_dir / "cases.csv"
    if prior.exists():
        pd.testing.assert_frame_equal(pd.read_csv(prior, float_precision="round_trip"), manifest,
                                      check_dtype=False, check_exact=True)
    write_table(manifest, prior)
    process_diagnostics = pd.DataFrame(process_rows)
    write_table(process_diagnostics, output_dir / "process_diagnostics.csv")
    write_table(pd.DataFrame(reference_rows), output_dir / "reference_input_comparison.csv")
    return jobs, process_diagnostics, clock_diagnostics


def gpdc_null_path(estimator_seed):
    effective_n = int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / DT)) + 1 - 2 * MAX_LAG
    return CACHE_DIR / f"null_T{effective_n}_n{GPDC_NULL_SAMPLES}_seed{estimator_seed}.npz"


def call_method(series, method, estimator_seed, null_file):
    if method == "PCMCI+":
        return run_pcmciplus(series, tau_max=MAX_LAG, pc_alpha=ALPHA, alpha=ALPHA,
            sig_samples=GPDC_NULL_SAMPLES, seed=estimator_seed, gp_params=GPDC_GP_PARAMS,
            null_dist_filename=str(null_file))
    if method == "LKIF":
        return run_lkif(series, max_lag=MAX_LAG, alpha=ALPHA)
    if method == "CCM":
        return run_ccm(series, E=CCM_E, tau=CCM_TAU, exclusion_radius=CCM_EXCLUSION_RADIUS,
            lib_size=CCM_LIBRARY_SIZE, sample=CCM_SAMPLES, n_surrogates=CCM_SURROGATES,
            seed=estimator_seed, alpha=ALPHA, lag_steps=CCM_LAG_STEPS,
            n_libraries=CCM_LAG_LIBRARIES, sign_threshold=CCM_LAG_SIGN_THRESHOLD)
    if method == "GC":
        return run_granger(series, max_lag=GC_ORDER, alpha=ALPHA)
    raise ValueError(f"Unknown method: {method}")


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


def empty_edges(method, status, error=""):
    return pd.DataFrame([{"method": method, "source": source, "target": target,
        "score": np.nan, "p_value": np.nan, "edge_detected": None,
        "unoriented": False, "status": status, "error": error}
        for source, target in [("x", "y"), ("y", "x")]])


def run_case(job, output_dir, null_file):
    """One worker owns each case; eight directional rows are its final commit."""
    output_dir = Path(output_dir)
    series = pd.read_csv(output_dir / "inputs" / f"{job['case_id']}.csv")
    tables = []
    with threadpool_limits(limits=1):
        for method in METHODS:
            started = time.perf_counter()
            try:
                result = call_method(series, method, job["estimator_seed"], null_file)
            except Exception as error:
                edges = empty_edges(method, "error", f"{type(error).__name__}: {error}")
                result = {"edges": edges, "raw": {}, "diagnostics": {"error": str(error)}}
            else:
                edges = result["edges"].copy()
                if "unoriented" not in edges:
                    edges["unoriented"] = False
                if (result["diagnostics"].get("undirected_contemporaneous", False)
                        or result["diagnostics"].get("conflicting_contemporaneous", False)):
                    edges["unoriented"] = True
            # Storage failures stop execution; they are not estimator failures.
            native_path = output_dir / "native" / f"{job['case_id']}_{method.replace('+', 'plus')}.json.gz"
            native_path.parent.mkdir(parents=True, exist_ok=True)
            temporary = native_path.with_suffix(".partial.gz")
            with gzip.open(temporary, "wt", encoding="utf-8") as handle:
                json.dump(native_snapshot(result), handle, separators=(",", ":"))
            replace_output(temporary, native_path)
            if method == "CCM":
                for key, suffix in [("lag_profile", "profile"), ("lag_peaks", "peaks")]:
                    if isinstance(result["raw"].get(key), pd.DataFrame):
                        write_table(result["raw"][key], output_dir / "lag" / f"{job['case_id']}_{suffix}.csv")
            edges["method"] = method
            edges["seconds"] = time.perf_counter() - started
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
    keys = ["method", "source", "target"]
    expected = {(method, source, target) for method in METHODS for source, target in [("x", "y"), ("y", "x")]}
    required = set(keys + ["status", "edge_detected", "unoriented", "score", "p_value", *job])
    if not required.issubset(result) or len(result) != len(expected) or result.duplicated(keys).any():
        return None
    if set(result[keys].itertuples(index=False, name=None)) != expected or result.status.isna().any():
        return None
    if result.status.eq("pending").any():
        return None
    for key, value in job.items():
        if not result[key].eq(value).all():
            raise ValueError(f"Saved case identity differs in {key}: {path.name}")
    for method in METHODS:
        if not (output_dir / "native" / f"{job['case_id']}_{method.replace('+', 'plus')}.json.gz").is_file():
            return None
    return result


def collect_results(jobs, completed):
    tables = []
    for job in jobs:
        result = completed.get(job["case_id"])
        if result is None:
            result = pd.concat([empty_edges(method, "pending") for method in METHODS], ignore_index=True)
            for key, value in job.items():
                result[key] = value
        tables.append(result)
    result = pd.concat(tables, ignore_index=True)
    for column in ["edge_detected", "unoriented", "raw_edge_detected"]:
        if column in result:
            result[column] = result[column].astype("boolean")
    return result


def summarize_results(results):
    """Keep per-seed conditional support separate from paired finite-score changes."""
    metadata = IDENTITY_COLUMNS + ["method"]
    required = metadata + ["source", "target", "score", "status", "edge_detected", "unoriented"]
    if results.empty or not set(required).issubset(results.columns):
        raise ValueError("Memory results require all requested cases and explicit statuses.")
    if results.duplicated(["case_id", "method", "source", "target"]).any() or results.status.isna().any():
        raise ValueError("Directional rows must be unique and have explicit status.")
    for column in ["edge_detected", "unoriented", "raw_edge_detected"]:
        if column in results and any(not isinstance(value, (bool, np.bool_)) for value in results[column].dropna()):
            raise ValueError(f"{column} must be boolean or missing.")
    state_rows = []
    for (_, method), pair in results.groupby(["case_id", "method"], sort=False):
        if len(pair) != 2 or set(zip(pair.source, pair.target)) != {("x", "y"), ("y", "x")}:
            raise ValueError("Every case and method needs both directions.")
        if any(pair[column].nunique(dropna=False) != 1 for column in metadata):
            raise ValueError("The two directions disagree on case identity.")
        row = pair.iloc[0][metadata].to_dict()
        row["pending"] = bool(pair.status.eq("pending").any())
        row["unresolved"] = bool(pair.unoriented.fillna(False).any())
        for source, target in [("x", "y"), ("y", "x")]:
            edge = pair[pair.source.eq(source)].iloc[0]
            direction = f"{source}_to_{target}"
            row[direction] = edge.edge_detected if edge.status == "ok" else pd.NA
            row[f"raw_{direction}"] = (edge.get("raw_edge_detected", pd.NA)
                                        if method == "CCM" and edge.status == "ok" else pd.NA)
        row["valid"] = all(pd.notna(row[name]) for name in ["x_to_y", "y_to_x"])
        row["raw_valid"] = all(pd.notna(row[f"raw_{name}"]) for name in ["x_to_y", "y_to_x"])
        row["failed"] = not row["pending"] and not row["valid"]
        row["correct"] = False
        if row["valid"] and not row["unresolved"]:
            forward, reverse = bool(row["x_to_y"]), bool(row["y_to_x"])
            row["correct"] = (not forward and not reverse if row["process"] == "independent" else forward and not reverse)
        row["ccm_lag_failed"] = bool(method == "CCM" and "lag_status" in pair and pair.lag_status.eq("failed").any())
        state_rows.append(row)
    states = pd.DataFrame(state_rows)
    for column in ["x_to_y", "y_to_x", "raw_x_to_y", "raw_y_to_x"]:
        states[column] = states[column].astype("boolean")
    group_keys = ["process_seed", "process", "memory", "condition", "method"]
    state_groups = dict(tuple(states.groupby(group_keys, sort=False)))
    summary_rows = []
    for key, group in results.groupby(group_keys, sort=False):
        paired = state_groups[key]
        expected_n = N_AGE_MEMBERS if key[3] == "chronology" else 1
        if len(paired) != expected_n:
            raise ValueError("Summary requires placeholders for all requested members.")
        valid, raw_valid = paired[paired.valid], paired[paired.raw_valid]
        row = group.iloc[0][[name for name in metadata if name not in ["case_id", "member", "input_status"]]].to_dict()
        row.update({"n_requested": expected_n, "n_valid": len(valid), "n_raw_valid": len(raw_valid),
                    "n_pending": int(paired.pending.sum()), "n_failed": int(paired.failed.sum()),
                    "n_unresolved": int(paired.unresolved.sum()), "n_correct_graph": int(paired.correct.sum()),
                    "n_ccm_lag_failed": int(paired.ccm_lag_failed.sum())})
        row["valid_fraction"] = len(valid) / expected_n
        row["fraction_correct_requested"] = row["n_correct_graph"] / expected_n
        row["fraction_unoriented"] = float(valid.unresolved.mean()) if len(valid) else np.nan
        row["support_threshold"] = SUPPORT_THRESHOLD
        for direction in ["x_to_y", "y_to_x"]:
            row[f"n_retained_{direction}"] = int(valid[direction].sum())
            row[f"q_{direction}"] = float(valid[direction].mean()) if len(valid) else np.nan
            row[f"n_raw_retained_{direction}"] = int(raw_valid[f"raw_{direction}"].sum())
            row[f"raw_q_{direction}"] = float(raw_valid[f"raw_{direction}"].mean()) if len(raw_valid) else np.nan
            scores = pd.to_numeric(group.loc[group.source.eq(direction[0]) & group.status.eq("ok"), "score"], errors="coerce")
            scores = scores[np.isfinite(scores)]
            quartiles = np.percentile(scores, [25, 50, 75]) if len(scores) else [np.nan] * 3
            row[f"n_scores_{direction}"] = len(scores)
            row[f"q25_{direction}"], row[f"median_{direction}"], row[f"q75_{direction}"] = quartiles
            row[f"iqr_{direction}"] = quartiles[2] - quartiles[0]
        row["decision"] = "insufficient_valid_members"
        row["high_unoriented"] = bool(len(valid) and row["fraction_unoriented"] >= SUPPORT_THRESHOLD)
        if len(valid) and row["valid_fraction"] >= MIN_VALID_FRACTION:
            forward = row["q_x_to_y"] >= SUPPORT_THRESHOLD
            reverse = row["q_y_to_x"] >= SUPPORT_THRESHOLD
            row["decision"] = ("both" if forward and reverse else "x_to_y" if forward else
                               "y_to_x" if reverse else "unoriented" if row["high_unoriented"] else "no_high_support")
        summary_rows.append(row)
    summaries = pd.DataFrame(summary_rows)
    keys = ["process_seed", "process", "memory", "method", "source", "target"]
    reference = results.loc[results.condition.eq("irregular_accurate"), keys + ["score", "status"]]
    reference = reference.rename(columns={"score": "baseline_score", "status": "baseline_status"})
    changes = results.loc[results.condition.eq("chronology")].merge(reference, on=keys, how="left", validate="many_to_one")
    usable = (changes.status.eq("ok") & changes.baseline_status.eq("ok")
              & np.isfinite(changes.score) & np.isfinite(changes.baseline_score))
    changes["score_change_usable"] = usable
    changes["delta_score"] = np.where(usable, changes.score - changes.baseline_score, np.nan)
    change_rows = []
    for key, group in changes.groupby(keys, sort=False):
        values = group.loc[group.score_change_usable, "delta_score"].to_numpy()
        quartiles = np.percentile(values, [25, 50, 75]) if len(values) else [np.nan] * 3
        change_rows.append({**dict(zip(keys, key)), "n_requested": len(group), "n_usable": len(values),
            "q25": quartiles[0], "median": quartiles[1], "q75": quartiles[2]})
    contrast_keys = ["process_seed", "process", "method", "source", "target", "member"]
    baseline_change = changes.loc[changes.memory.eq("reference"), contrast_keys + ["delta_score", "score_change_usable"]]
    baseline_change = baseline_change.rename(columns={"delta_score": "reference_delta_score", "score_change_usable": "reference_change_usable"})
    contrasts = changes.loc[~changes.memory.eq("reference")].merge(baseline_change, on=contrast_keys,
                                                                  how="left", validate="many_to_one")
    valid_contrast = contrasts.score_change_usable & contrasts.reference_change_usable.fillna(False)
    contrasts["contrast_usable"] = valid_contrast
    contrasts["difference_in_score_change"] = np.where(valid_contrast,
        contrasts.delta_score - contrasts.reference_delta_score, np.nan)
    return {"scenario_summary": summaries, "case_states": states,
            "paired_score_changes": changes, "paired_change_summary": pd.DataFrame(change_rows),
            "memory_contrasts": contrasts, "ccm_lag_summary": results.loc[results.method.eq("CCM")].copy()}


def write_status(state, completed, requested, started, output_dir=DATA_DIR, error=""):
    write_json({"state": state, "completed_pairs": completed, "requested_pairs": requested,
        "n_workers": N_WORKERS, "n_process_histories": len(PROCESS_SEEDS), "n_age_members": N_AGE_MEMBERS,
        "n_methods": len(METHODS), "elapsed_seconds": time.perf_counter() - started,
        "updated_utc": datetime.now(timezone.utc).isoformat(), "error": error}, Path(output_dir) / "run_status.json")


def save_outputs(results, process_diagnostics, output_dir=DATA_DIR):
    tables = summarize_results(results)
    write_table(results, Path(output_dir) / "member_results.csv")
    for name, table in tables.items():
        write_table(table, Path(output_dir) / f"{name}.csv")
    if "plot_results" in globals():
        plot_results(tables, process_diagnostics, output_dir=FIGURE_DIR, data_dir=Path(output_dir))
    return tables


def main():
    if N_WORKERS < 1 or N_AGE_MEMBERS < 1 or len(set(PROCESS_SEEDS)) != len(PROCESS_SEEDS):
        raise ValueError("Provide positive counts and distinct process seeds.")
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
            print(f"Prepared {len(jobs)} memory inputs; {len(completed)} complete.", flush=True)
            return
        pending = [job for job in jobs if job["case_id"] not in completed]
        if RUN_CAUSAL and pending:
            effective_n = int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / DT)) + 1 - 2 * MAX_LAG
            with threadpool_limits(limits=1):
                for estimator_seed in sorted({job["estimator_seed"] for job in pending}):
                    prepare_gpdc_null(gpdc_null_path(estimator_seed), sample_sizes=[effective_n],
                                      sig_samples=GPDC_NULL_SAMPLES, seed=estimator_seed)
            print(f"Memory experiment: {len(completed)} complete, {len(pending)} pending; {N_WORKERS} workers.", flush=True)
            write_status("running", len(completed), len(jobs), started)
            context = multiprocessing.get_context("spawn")
            with ProcessPoolExecutor(max_workers=min(N_WORKERS, len(pending)), mp_context=context) as pool:
                futures = [pool.submit(run_case, job, str(DATA_DIR), str(gpdc_null_path(job["estimator_seed"]))) for job in pending]
                for future in as_completed(futures):
                    result = future.result()
                    completed[result.case_id.iloc[0]] = result
                    if len(completed) % 10 == 0 or len(completed) == len(jobs):
                        write_status("running", len(completed), len(jobs), started)
                        print(f"Completed {len(completed)}/{len(jobs)} memory pairs.", flush=True)
        results = collect_results(jobs, completed)
        save_outputs(results, process_diagnostics)
        state = "complete" if len(completed) == len(jobs) else "partial"
        write_status(state, len(completed), len(jobs), started)
        print(f"Memory experiment {state}: {len(completed)}/{len(jobs)} committed pairs.", flush=True)
    except BaseException as error:
        if jobs:
            write_table(collect_results(jobs, completed), DATA_DIR / "member_results.csv")
        write_status("interrupted_or_failed", len(completed), len(jobs), started,
                     error=f"{type(error).__name__}: {error}")
        raise


def _memory_plot_labels():
    """Keep the four factorial conditions in their scientific reading order."""
    memories = [item["memory"] for item in MEMORY_SETTINGS]
    labels = [f"({item['memory_x_years']:g}, {item['memory_y_years']:g})"
              for item in MEMORY_SETTINGS]
    return memories, labels


def _save_memory_figure(fig, output_dir, name):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{name}.png", dpi=180)
    fig.savefig(output_dir / f"{name}.pdf")
    plt.close(fig)


def _plot_graph_recovery(summaries, output_dir):
    """Compare pooled graph recovery against each accurate-age reference."""
    memories, memory_labels = _memory_plot_labels()
    names = ["Reference", "Short X", "Short Y", "Short X + Y"]
    conditions = [("original_uniform", "Original uniform, accurate age", "0.60", -0.23),
                  ("irregular_accurate", "Irregular sampling, accurate age", "0.22", 0),
                  ("chronology", "Irregular sampling + age uncertainty", "#009E73", 0.23)]
    fig, axes = plt.subplots(3, 4, figsize=(22, 14), sharey=True, squeeze=False)
    for row, process in enumerate(PROCESSES):
        for column, method in enumerate(METHODS):
            ax = axes[row, column]
            selected = summaries[(summaries["process"] == process) & (summaries["method"] == method)]
            labels = []
            for position, memory in enumerate(memories):
                memory_rows = selected[selected["memory"] == memory]
                counts, values = {}, {}
                for condition, label, color, offset in conditions:
                    subset = memory_rows[memory_rows["condition"] == condition]
                    requested = int(subset["n_requested"].sum())
                    correct = int(subset["n_correct_graph"].sum())
                    value = correct / requested if requested else np.nan
                    counts[condition] = f"{correct}/{requested}"
                    values[condition] = value
                    if np.isfinite(value):
                        ax.scatter(position + offset, value, s=58, marker="o",
                                   facecolor="white" if condition == "original_uniform" else color,
                                   edgecolor=color, linewidth=1.2, zorder=4)
                        if condition == "chronology":
                            ax.annotate(counts[condition], (position + offset, value),
                                        xytext=(0, 8), textcoords="offset points", ha="center",
                                        fontsize=9, fontweight="bold", color="#006B50")
                if np.isfinite(values["irregular_accurate"]) and np.isfinite(values["chronology"]):
                    ax.plot([position, position + 0.23],
                            [values["irregular_accurate"], values["chronology"]],
                            color="0.72", lw=1.5, zorder=2)
                age = memory_rows[memory_rows["condition"] == "chronology"]
                unresolved, failed, pending = [int(age[key].sum())
                                              for key in ["n_unresolved", "n_failed", "n_pending"]]
                labels.append(f"{names[position]}\n{memory_labels[position]}\n"
                              f"O {counts['original_uniform']} | I {counts['irregular_accurate']}\n"
                              f"U{unresolved} F{failed} P{pending}")
            ax.set_xticks(range(4), labels, fontsize=8.5)
            ax.set_xlim(-0.5, 3.55)
            ax.set_ylim(-0.05, 1.16)
            ax.set_yticks(np.linspace(0, 1, 6), ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=10)
            ax.grid(axis="y", color="0.90", lw=0.7)
            ax.set_axisbelow(True)
            ax.set_xlabel(r"Memory $(T_X, T_Y)$, years", fontsize=10)
            truth = "X -> Y only" if process != "independent" else "Neither direction"
            ax.set_title(f"{process.capitalize()} | {method}\nTrue graph: {truth}", fontsize=12)
            if column == 0:
                ax.set_ylabel("Complete graph recovered", fontsize=11)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
    handles = [Line2D([], [], marker="o", color=color, ls="none", markersize=7,
                      markerfacecolor="white" if condition == "original_uniform" else color, label=label)
               for condition, label, color, _ in conditions]
    fig.suptitle("1. Does the method recover the complete causal graph?\n"
                 "Compare the accurate-age baseline with the same data under chronology uncertainty", fontsize=17, y=0.99)
    fig.legend(handles=handles, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 0.935),
               frameon=False, fontsize=11)
    fig.text(0.5, 0.044,
             "Green labels: correct / requested chronology inputs. Gray links show accurate -> perturbed change. "
             "O / I: original / irregular accurate-age counts.", ha="center", fontsize=10)
    fig.text(0.5, 0.025,
             "U / F / P: unresolved / failed / pending chronology cases. A graph is correct only when both directions are correct and no evidence is unresolved.",
             ha="center", fontsize=10)
    fig.text(0.5, 0.008,
             "3 histories x 10 shared clock pairs; descriptive pooled fractions. K=4 for coupled systems; true delay=20 years; C-C proposal SD=80 years.",
             ha="center", fontsize=10)
    add_panel_labels(axes, x=-0.10)
    fig.tight_layout(rect=(0.015, 0.075, 1, 0.90), h_pad=2.5, w_pad=1.8)
    _save_memory_figure(fig, output_dir, "01_graph_recovery")


def _plot_paired_scores(changes, output_dir):
    """Display pooled paired changes, with a clearly expanded independent row."""
    memories, memory_labels = _memory_plot_labels()
    names = ["Reference", "Short X", "Short Y", "Short X + Y"]
    directions = [("x", "y", "X -> Y", "#0072B2", -0.18),
                  ("y", "x", "Y -> X", "#D55E00", 0.18)]
    score_labels = {"PCMCI+": "Change in max. GPDC statistic",
                    "LKIF": "Change in max. |normalized flow|",
                    "CCM": "Change in original Tp=0 rho",
                    "GC": "Change in log variance ratio"}
    fig, axes = plt.subplots(3, 4, figsize=(22, 14), squeeze=False)
    for row, process in enumerate(PROCESSES):
        for column, method in enumerate(METHODS):
            ax = axes[row, column]
            selected = changes[(changes["process"] == process) & (changes["method"] == method)]
            count_labels = []
            for position, memory in enumerate(memories):
                counts = []
                for source, target, label, color, offset in directions:
                    subset = selected[(selected["memory"] == memory) & (selected["source"] == source)
                                      & (selected["target"] == target)]
                    usable = subset["score_change_usable"].fillna(False).astype(bool)
                    paired = subset.loc[usable].sort_values(["process_seed", "member"])
                    values = paired["delta_score"].to_numpy(dtype=float)
                    values = values[np.isfinite(values)]
                    counts.append(f"{len(values)}/{len(subset)}")
                    x = position + offset
                    if not len(values):
                        continue
                    jitter = np.random.default_rng(0).uniform(-0.065, 0.065, len(values))
                    ax.scatter(x + jitter, values, marker="o", s=17, color=color, alpha=0.30,
                               linewidths=0, zorder=2)
                    q25, median, q75 = np.percentile(values, [25, 50, 75])
                    ax.errorbar(x, median, yerr=[[median-q25], [q75-median]], marker="o", color=color,
                                markersize=7, markeredgecolor="white", markeredgewidth=0.8,
                                linewidth=2, capsize=4, linestyle="none", zorder=4)
                count_labels.append(f"{names[position]}\n{memory_labels[position]}\n"
                                    f"X->Y {counts[0]}\nY->X {counts[1]}")
            # Keep coupled rows comparable within a method. The independent row
            # gets its own explicitly labeled range so small changes stay visible.
            limit_rows = changes[(changes["method"] == method)
                                 & (changes["process"].eq("independent") if process == "independent"
                                    else ~changes["process"].eq("independent"))]
            usable_limits = limit_rows["score_change_usable"].fillna(False).astype(bool)
            finite = limit_rows.loc[usable_limits, "delta_score"].to_numpy(dtype=float)
            finite = finite[np.isfinite(finite)]
            if len(finite):
                low, high = min(0.0, float(finite.min())), max(0.0, float(finite.max()))
                span = max(high-low, 0.001)
                ax.set_ylim(low-0.10*span, high+0.12*span)
            ax.axhline(0, color="0.25", lw=1.1, zorder=1)
            ax.grid(axis="y", color="0.92", lw=0.7)
            ax.set_axisbelow(True)
            ax.set_xticks(range(4), count_labels, fontsize=8.5)
            ax.set_xlim(-0.5, 3.5)
            ax.tick_params(axis="y", labelsize=10)
            ax.set_xlabel(r"Memory $(T_X, T_Y)$, years", fontsize=10)
            title = f"{process.capitalize()} | {method}"
            if process == "independent":
                title = title + "\nExpanded y-scale"
            ax.set_title(title, fontsize=12)
            ax.set_ylabel(score_labels[method], fontsize=10)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
    handles = [Line2D([], [], color=color, marker="o", ls="none", markersize=7, label=label)
               for _, _, label, color, _ in directions]
    handles.append(Line2D([], [], color="0.3", marker="o", ls="-", label="Large point / bar: median / IQR"))
    fig.suptitle("4. How much does chronology change the native score?\n"
                 "Change = chronology member score - its own irregular accurate-age baseline score", fontsize=17, y=0.99)
    fig.legend(handles=handles, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 0.935),
               frameon=False, fontsize=11)
    fig.text(0.5, 0.044,
             "Small points: individual paired changes. Large points and bars: pooled median and IQR. "
             "Tick counts: usable / requested paired changes.", ha="center", fontsize=10)
    fig.text(0.5, 0.025,
             "Linear and nonlinear rows share limits within each method; the independent row uses an expanded scale. "
             "Method units differ; score loss is not a graph-error rate.", ha="center", fontsize=10)
    fig.text(0.5, 0.008,
             "3 histories x 10 shared clocks; descriptive spread, not a confidence interval. "
             "CCM uses original Tp=0 rho, including successful scores with an unresolved later lag screen.", ha="center", fontsize=10)
    add_panel_labels(axes, x=-0.10)
    fig.tight_layout(rect=(0.015, 0.075, 1, 0.90), h_pad=2.5, w_pad=2.2)
    _save_memory_figure(fig, output_dir, "04_paired_native_score_changes")


def _plot_ccm_lag_profiles(lag_summary, data_dir, output_dir):
    """Keep the existing curve aggregation and show saved peak-sign diagnostics."""
    memories, memory_labels = _memory_plot_labels()
    seeds = sorted(lag_summary["process_seed"].unique())
    coupled = [name for name in PROCESSES if name != "independent"]
    directions = [("x", "y", "X → Y", "#0072B2"), ("y", "x", "Y → X", "#D55E00")]
    metadata = ["case_id", "process_seed", "process", "memory", "condition", "member"]
    rows = []
    for item in lag_summary[metadata].drop_duplicates().to_dict("records"):
        if item["process"] not in coupled or item["condition"] == "original_uniform":
            continue
        path = Path(data_dir) / "lag" / f"{item['case_id']}_profile.csv"
        if not path.exists():
            continue
        profile = pd.read_csv(path, float_precision="round_trip")
        if not {"source", "target", "lag_steps", "rho_mean"}.issubset(profile.columns):
            raise ValueError(f"CCM profile columns are incomplete: {path.name}")
        for key, value in item.items():
            profile[key] = value
        rows.append(profile)
    if rows:
        profiles = pd.concat(rows, ignore_index=True)
    else:
        profiles = pd.DataFrame(columns=metadata + ["source", "target", "lag_steps", "rho_mean"])
    finite_rho = pd.to_numeric(profiles["rho_mean"], errors="coerce")
    finite_rho = finite_rho[np.isfinite(finite_rho)]
    rho_lower = min(-0.15, float(finite_rho.min()) - 0.05) if len(finite_rho) else -0.15

    fig = plt.figure(figsize=(18.5, 21.5))
    grid = fig.add_gridspec(len(memories), len(coupled), left=0.085, right=0.975,
                            bottom=0.115, top=0.902, hspace=0.24, wspace=0.27)
    axes = np.empty((len(memories), len(coupled)), dtype=object)
    for row, memory in enumerate(memories):
        for column, process in enumerate(coupled):
            # The strip has its own fraction axis, independent of target shift.
            panel = grid[row, column].subgridspec(3, 1, height_ratios=[4.2, 1.6, 0.62], hspace=0.75)
            ax = fig.add_subplot(panel[0])
            strip = fig.add_subplot(panel[1])
            note_ax = fig.add_subplot(panel[2])
            note_ax.set_axis_off()
            axes[row, column] = ax
            notes = []
            for direction_index, (source, target, label, color) in enumerate(directions):
                subset = profiles[(profiles["memory"] == memory) & (profiles["process"] == process)
                                  & (profiles["source"] == source) & (profiles["target"] == target)]
                peak_rows = lag_summary[(lag_summary["memory"] == memory)
                                        & (lag_summary["process"] == process)
                                        & (lag_summary["source"] == source) & (lag_summary["target"] == target)]
                counts = []
                peak_labels = []
                diagnostic_counts = []
                for condition, linestyle in [("irregular_accurate", "--"), ("chronology", "-")]:
                    condition_rows = subset[subset["condition"] == condition]
                    requested = len(peak_rows[peak_rows["condition"] == condition])
                    available = condition_rows.loc[np.isfinite(condition_rows["rho_mean"].astype(float)), "case_id"].nunique()
                    counts.append(f"{available}/{requested}")
                    seed_curves = []
                    for seed in seeds:
                        one_seed = condition_rows[condition_rows["process_seed"] == seed]
                        if not one_seed.empty:
                            # Each saved rho_mean averages that case's paired libraries.
                            # The next median is across chronology members, not libraries.
                            curve = one_seed.groupby("lag_steps")["rho_mean"].median().rename(seed)
                            seed_curves.append(curve)
                    if seed_curves:
                        wide = pd.concat(seed_curves, axis=1).sort_index()
                        values = wide.to_numpy(dtype=float)
                        finite = np.isfinite(values).any(axis=1)
                        if finite.any():
                            lags = wide.index.to_numpy(dtype=float)[finite] * DT * 1000
                            middle = np.nanmedian(values[finite], axis=1)
                            ax.plot(lags, middle, color=color, ls=linestyle, lw=2.0)
                            if condition == "chronology":
                                lower = np.nanmin(values[finite], axis=1)
                                upper = np.nanmax(values[finite], axis=1)
                                ax.fill_between(lags, lower, upper, color=color, alpha=0.13, linewidth=0)
                    peaks = peak_rows[peak_rows["condition"] == condition].groupby("process_seed")["optimal_lag_steps"].median()
                    finite_peaks = peaks[np.isfinite(peaks)]
                    peak_labels.append(f"{float(finite_peaks.median()) * DT * 1000:+.0f}" if len(finite_peaks) else "NA")

                    # Display saved per-input fractions; no peak classification or
                    # significance test is recomputed here. Jitter changes y only.
                    diagnostic = peak_rows[peak_rows["condition"] == condition].sort_values(["process_seed", "member"])
                    fraction_column = "lag_past_fraction" if source == "x" else "lag_future_fraction"
                    fractions = diagnostic[fraction_column].to_numpy(dtype=float)
                    fractions = fractions[np.isfinite(fractions)]
                    diagnostic_counts.append(f"{len(fractions)}/{requested}")
                    filled = condition == "chronology"
                    level = (1.3 if source == "x" else 0.0) + (-0.16 if filled else 0.20)
                    jitter = np.random.default_rng(0).uniform(-0.12, 0.12, len(fractions))
                    strip.scatter(fractions, level + jitter, marker="o", s=20 if filled else 38,
                                  facecolors=color if filled else "white", edgecolors=color,
                                  linewidths=0.65 if filled else 1.25, alpha=0.65 if filled else 1,
                                  zorder=3 if filled else 4)
                notes.append((color, f"{label}: profiles A {counts[0]}, C {counts[1]}; "
                              f"peak medians A {peak_labels[0]}, C {peak_labels[1]} yr"))
                if diagnostic_counts != counts:
                    notes[-1] = (color, notes[-1][1] + "; diagnostic n " + ", ".join(diagnostic_counts))

            ax.axvline(0, color="0.3", ls=":", lw=1.0)
            ax.axvline(-TRUE_LAG * 1000, color="0.65", ls="--", lw=1.0)
            ax.set_title(f"{process.capitalize()}  ·  (T_X, T_Y) = {memory_labels[row]} years", fontsize=12, pad=8)
            ax.set_ylabel("Cross-map ρ\n(library mean)", fontsize=11)
            ax.set_xlabel("Target shift (years; negative = past)", fontsize=10, labelpad=3)
            ax.set_xlim(min(CCM_LAG_STEPS) * DT * 1000, max(CCM_LAG_STEPS) * DT * 1000)
            ax.set_ylim(rho_lower, 1.02)
            ax.set_xticks(np.arange(-300, 301, 100))
            ax.tick_params(axis="both", labelsize=10)
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(axis="y", color="0.92", lw=0.6)
            ax.set_axisbelow(True)

            strip.axvline(CCM_LAG_SIGN_THRESHOLD, color="0.35", ls=":", lw=1.0, zorder=1)
            strip.set_xlim(-0.035, 1.035)
            strip.set_ylim(-0.45, 1.75)
            strip.set_xticks([0, 0.5, 0.9, 1.0], ["0", "0.5", "0.9", "1"])
            strip.set_yticks([1.3, 0], ["X → Y: past", "Y → X: future"])
            for tick, (_, _, _, color) in zip(strip.get_yticklabels(), directions):
                tick.set_color(color)
            strip.set_xlabel("Peak-sign diagnostic: fraction of 100 libraries per input", fontsize=9, labelpad=2)
            strip.tick_params(axis="x", labelsize=9, length=3, pad=2)
            strip.tick_params(axis="y", labelsize=10, length=0, pad=6)
            strip.spines[["top", "right", "left"]].set_visible(False)
            strip.spines["bottom"].set_color("0.6")
            for index, (color, note) in enumerate(notes):
                note_ax.text(0, 1 - index * 0.62, note, color=color, fontsize=9,
                             ha="left", va="top", transform=note_ax.transAxes)

    handles = [Line2D([], [], color=color, lw=2, label=label) for _, _, label, color in directions]
    handles.extend([Line2D([], [], color="0.3", ls="--", marker="o", markerfacecolor="white",
                           label="A: irregular accurate age (3 inputs)"),
                    Line2D([], [], color="0.3", ls="-", marker="o",
                           label="C: C-C chronology (30 inputs)")])
    fig.suptitle("3. CCM lag profiles and per-input peak-sign diagnostics\n"
                 "Four memory conditions · L = 300 · 100 paired libraries · 31 target shifts", fontsize=16, y=0.991)
    fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.955),
               frameon=False, fontsize=11, handlelength=3, columnspacing=3)
    fig.text(0.5, 0.083,
             "Curves: library mean within an input → median members within a history → median of 3 histories. "
             "Shading: range of the 3 history curves, not a confidence interval.", ha="center", fontsize=10)
    fig.text(0.5, 0.064,
             "Open circles: 3 accurate inputs. Filled circles: 30 chronology inputs (3 histories × 10 shared clocks). "
             "Strip positions are saved fractions; vertical jitter is only for display.", ha="center", fontsize=10)
    fig.text(0.5, 0.045,
             "The 0.9 strip guide applies to each input's 100 libraries, not to a pooled retained fraction. "
             "Sign fractions alone do not determine the final arrow.", ha="center", fontsize=10)
    fig.text(0.5, 0.026,
             "Retention also requires original significance, complete valid curves, and a negative interior mean-curve peak. "
             "Gray −20-year guide: imposed X → Y delay; dotted line: zero shift.", ha="center", fontsize=10)
    fig.text(0.5, 0.009,
             "Notes give available/requested profiles and median peak shifts (A = accurate, C = chronology). "
             "A reconstruction peak need not equal the imposed delay; native strength remains Tp = 0 ρ.", ha="center", fontsize=10)
    add_panel_labels(axes, x=-0.10, y=1.03)
    _save_memory_figure(fig, output_dir, "03_ccm_lag_profiles")



def _plot_ccm_screening_state(row, screened):
    """Classify one requested CCM direction without discarding raw evidence."""
    raw = row.get("raw_edge_detected", pd.NA)
    edge = row.get("edge_detected", pd.NA)
    unresolved = row.get("unoriented", pd.NA)
    raw_known = pd.notna(raw) and isinstance(raw, (bool, np.bool_))
    edge_known = pd.notna(edge) and isinstance(edge, (bool, np.bool_))
    if row.get("status") != "ok" or not raw_known:
        return "unavailable"
    if not screened:
        return "retained" if bool(raw) else "nonsignificant"
    if not edge_known or row.get("lag_status") == "failed":
        return "unavailable"
    if not bool(raw):
        if bool(edge):
            raise ValueError("A lag screen cannot upgrade a nonsignificant CCM direction.")
        return "nonsignificant"
    if bool(edge) and row.get("lag_status") == "compatible":
        return "retained"
    if not bool(edge) and row.get("lag_status") == "unresolved" and pd.notna(unresolved) and bool(unresolved):
        return "unresolved"
    if not bool(edge) and row.get("lag_status") == "incompatible":
        return "rejected"
    raise ValueError(f"Unexpected CCM screening state: {row.get('case_id')}, "
                     f"{row.get('source')} -> {row.get('target')}.")


def _plot_ccm_case_counts(rows):
    """Preserve input-level U/F/P counts rather than adding directional flags."""
    unresolved, failed, pending = 0, 0, 0
    for case_id, pair in rows.groupby("case_id", sort=False):
        if len(pair) != 2 or set(zip(pair["source"], pair["target"])) != {("x", "y"), ("y", "x")}:
            raise ValueError(f"Both CCM directions are required for {case_id}.")
        case_pending = bool(pair["status"].eq("pending").any())
        case_unresolved = bool(pair["unoriented"].fillna(False).any())
        case_valid = bool(pair["status"].eq("ok").all() and pair["edge_detected"].notna().all())
        unresolved = unresolved + int(case_unresolved)
        pending = pending + int(case_pending)
        failed = failed + int(not case_pending and not case_valid)
    return unresolved, failed, pending


def _plot_ccm_screening(ccm_rows, output_dir):
    """Show all requested directional states before and after CCM lag screening."""
    from matplotlib.patches import Patch

    memories, memory_labels = _memory_plot_labels()
    required = {"case_id", "process", "memory", "condition", "source", "target", "status",
                "raw_edge_detected", "edge_detected", "unoriented", "lag_status"}
    if not required.issubset(ccm_rows.columns):
        raise ValueError("The CCM state figure requires the complete directional result table.")
    if ccm_rows.duplicated(["case_id", "source", "target"]).any():
        raise ValueError("Duplicate CCM directional rows cannot be pooled in a state bar.")

    colors = {"x": "#0072B2", "y": "#D55E00", "nonsignificant": "#D9D9D9",
              "rejected": "#8C8C8C", "unresolved": "#E6C84F", "unavailable": "#8B2635"}
    states = ["retained", "nonsignificant", "rejected", "unresolved", "unavailable"]
    hatches = {"retained": "", "nonsignificant": "", "rejected": "",
               "unresolved": "///", "unavailable": "xx"}
    groups = [("x", "y", False, "Raw\nX → Y"),
              ("x", "y", True, "Screened\nX → Y"),
              ("y", "x", False, "Raw\nY → X"),
              ("y", "x", True, "Screened\nY → X")]
    conditions = [("irregular_accurate", "A", -0.18), ("chronology", "C", 0.18)]
    bar_width = 0.27
    fig, axes = plt.subplots(len(PROCESSES), len(memories), figsize=(22, 13),
                             sharey=True, squeeze=False)
    for row_index, process in enumerate(PROCESSES):
        for column, memory in enumerate(memories):
            ax = axes[row_index, column]
            selected = ccm_rows[(ccm_rows["process"] == process) & (ccm_rows["memory"] == memory)]
            case_notes = []
            for condition, short_label, offset in conditions:
                condition_rows = selected[selected["condition"] == condition]
                requested = len(PROCESS_SEEDS) * (N_AGE_MEMBERS if condition == "chronology" else 1)
                if condition_rows["case_id"].nunique() != requested:
                    raise ValueError(f"Missing requested CCM inputs: {process}, {memory}, {condition}.")
                n_unresolved, n_failed, n_pending = _plot_ccm_case_counts(condition_rows)
                case_notes.append(f"{short_label}: n={requested}; U{n_unresolved} F{n_failed} P{n_pending}")
                for position, (source, target, screened, label) in enumerate(groups):
                    directions = condition_rows[(condition_rows["source"] == source)
                                                & (condition_rows["target"] == target)]
                    if len(directions) != requested:
                        raise ValueError("Every requested input must supply a directional row, including failures.")
                    counts = {state: 0 for state in states}
                    for record in directions.to_dict("records"):
                        state = _plot_ccm_screening_state(record, screened)
                        counts[state] = counts[state] + 1
                    if sum(counts.values()) != requested:
                        raise ValueError("CCM state counts must sum to the requested denominator.")
                    x = position + offset
                    bottom = 0.0
                    for state in states:
                        count = counts[state]
                        if not count:
                            continue
                        fraction = count / requested
                        color = colors[source] if state == "retained" else colors[state]
                        ax.bar(x, fraction, width=bar_width, bottom=bottom, color=color,
                               edgecolor="white" if not hatches[state] else "0.25", linewidth=0.6,
                               hatch=hatches[state], zorder=3)
                        middle = bottom + fraction / 2
                        if fraction >= 0.10:
                            text_color = "white" if state in {"retained", "unavailable"} else "black"
                            ax.text(x, middle, str(count), color=text_color, fontsize=10,
                                    fontweight="bold", ha="center", va="center", zorder=4)
                        else:
                            # Keep small nonzero segments legible without enlarging their fractions.
                            ax.annotate(str(count), xy=(x + bar_width / 2, middle), xytext=(6, 0),
                                        textcoords="offset points", ha="left", va="center", fontsize=9,
                                        arrowprops={"arrowstyle": "-", "color": "0.25", "lw": 0.7},
                                        zorder=4)
                        bottom = bottom + fraction
                    ax.text(x, -0.025, short_label, transform=ax.get_xaxis_transform(),
                            ha="center", va="top", fontsize=10, fontweight="bold")
            ax.set_xticks(range(len(groups)), [label for _, _, _, label in groups], fontsize=10)
            ax.tick_params(axis="x", length=0, pad=20)
            ax.set_xlim(-0.53, len(groups) - 0.47)
            ax.set_ylim(0, 1.05)
            ax.set_yticks([0, 0.25, 0.5, 0.75, 1], ["0%", "25%", "50%", "75%", "100%"], fontsize=10)
            ax.set_axisbelow(True)
            ax.grid(axis="y", color="0.9", linewidth=0.7)
            ax.set_title(f"{process.capitalize()} | {memory_labels[column]} yr", fontsize=12, pad=10)
            ax.set_xlabel("   ".join(case_notes), fontsize=9, labelpad=8)
            if column == 0:
                ax.set_ylabel("Share of requested inputs", fontsize=11)

    handles = [Patch(facecolor=colors["x"], label="Significant / retained: X → Y"),
               Patch(facecolor=colors["y"], label="Significant / retained: Y → X"),
               Patch(facecolor=colors["nonsignificant"], label="Not significant in raw test"),
               Patch(facecolor=colors["rejected"], label="Rejected by lag screen"),
               Patch(facecolor=colors["unresolved"], edgecolor="0.25", hatch=hatches["unresolved"],
                     label="Unresolved after lag screen"),
               Patch(facecolor=colors["unavailable"], edgecolor="0.25", hatch=hatches["unavailable"],
                     label="Unavailable")]
    fig.suptitle("2. CCM: raw significance and lag-screened states", fontsize=18, y=0.995)
    fig.legend(handles=handles, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 0.967),
               frameon=False, fontsize=11, columnspacing=2.2, handlelength=2.1)
    fig.text(0.5, 0.895,
             "A = irregular accurate age (3 inputs); C = chronology perturbed (3 histories × 10 shared clocks). "
             "Segment labels are counts; titles give (T_X, T_Y).",
             ha="center", fontsize=11)
    fig.text(0.5, 0.035,
             "Bars include every requested input. U/F/P count unresolved/failed/pending cases once per input, across both directions.",
             ha="center", fontsize=11)
    fig.text(0.5, 0.015,
             "A nonsignificant test does not establish absence. Raw evidence is preserved if a lag scan fails; the screened result is then unavailable.",
             ha="center", fontsize=11)
    add_panel_labels(axes, x=-0.08)
    fig.tight_layout(rect=(0, 0.07, 1, 0.875), h_pad=2.5, w_pad=2.0)
    _save_memory_figure(fig, output_dir, "02_ccm_raw_screened_retention")



def plot_results(tables, process_diagnostics, output_dir=None, data_dir=None):
    """Plot all-method memory contrasts and the full CCM-screen explanation."""
    output_dir = FIGURE_DIR if output_dir is None else Path(output_dir)
    data_dir = DATA_DIR if data_dir is None else Path(data_dir)
    # Process diagnostics remain a separate numerical table, including actual
    # autocorrelation and forced/background amplitudes; no data are re-fitted.
    _plot_graph_recovery(tables["scenario_summary"], output_dir)
    _plot_ccm_screening(tables["ccm_lag_summary"], output_dir)
    _plot_ccm_lag_profiles(tables["ccm_lag_summary"], data_dir, output_dir)
    _plot_paired_scores(tables["paired_score_changes"], output_dir)


if __name__ == "__main__":
    main()
