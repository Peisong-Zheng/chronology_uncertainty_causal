"""Chronology sensitivity of causal evidence between the TBK southern/northern endpoints.

Three independent noisy four-box histories, one accurate-age reference and ten
members at each of four C-C chronology levels: 123 inputs, four native methods.
Bidirectional evidence is a pipeline outcome, not recovery of the hidden network.
See agent/04_tbk_chronology_plan.md for the scientific scope and interpretation.
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
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from toolbox.plotting import save_figure
from threadpoolctl import threadpool_limits

from toolbox.tbk_observation import prepare_tbk_records
from toolbox.causal import prepare_gpdc_null, run_pcmciplus, run_lkif, run_ccm, run_granger


RUN_CAUSAL = True             # False rebuilds summaries from committed results without fitting.
PREPARE_ONLY = False          # True validates and saves inputs, then stops before fits.
RESUME = True                 # Preserve completed cases, including unfavorable/failed evidence.
N_WORKERS = 32                # Existing measured Windows setting; one numerical thread each.
PROCESS_SEEDS = [20261009, 20261301, 20261302]
TOTAL_YEARS = 39000           # Original annual update; nine kyr startup plus 30 kyr retained.
SPINUP_YEARS = 9000
MODEL_STEP_YEARS = 1
D_SOUTH = 0.005               # Original notebook cell-12 diffusion variant.
PROCESS_NOISE_SD = 0.005      # Original intrinsic process noise; no measurement noise added.
DT = 0.020                   # Twenty-year observation and common analysis grid, kyr.
ANALYSIS_OLD = 29.0           # Fixed ages, kyr BP; no member-specific cropping/extrapolation.
ANALYSIS_YOUNG = 1.0
AGE_SIGMA_YEARS = [20.0, 50.0, 100.0, 200.0]
AGE_CORRELATION_YEARS = 1000.0  # One kyr on the uniformly sampled 30-kyr record (1/30).
MAX_AGE_DRAWS = 10000         # Whole-curve rejection budget per observed center or ensemble.
N_AGE_MEMBERS = 10            # Around each uncertain observed center, not around known truth.
METHOD_SEED_OFFSET = 400      # One estimator stream per history, paired across age levels.
SUPPORT_THRESHOLD = 0.90      # Reporting rule; not a calibrated causal probability.
MIN_VALID_FRACTION = 0.90

MAX_LAG = 5                  # PCMCI+/LKIF history: 100 years.
GC_ORDER = 1                 # Fixed 20-year history, unchanged across chronology members.
ALPHA = 0.05
GPDC_NULL_SAMPLES = 10000     # Effective PCMCI null sample size: 1401 - 2*5 = 1391.
GPDC_GP_PARAMS = {"optimizer": None}
CCM_E = 3
CCM_TAU = -1                 # Backward 20-year embedding spacing; 40-year span.
CCM_EXCLUSION_RADIUS = 10    # Two hundred years.
CCM_LIBRARY_SIZE = 300       # Absolute size, used for Tp=0, nulls and lag scans.
CCM_SAMPLES = 10
CCM_SURROGATES = 100
CCM_LAG_STEPS = list(range(-30, 31))  # Sixty-one shifts, -600 to +600 years.
CCM_LAG_LIBRARIES = 100      # Same random libraries across all shifts on common anchors.
CCM_LAG_SIGN_THRESHOLD = 0.90

STEM = Path(__file__).stem
DATA_DIR = ROOT / "data/processed" / STEM
FIGURE_DIR = ROOT / "figures" / STEM
CACHE_DIR = ROOT / "data/cache/gpdc"
METHODS = ["PCMCI+", "LKIF", "CCM", "GC"]
STATES = ["both", "x_to_y", "y_to_x", "neither", "unresolved", "failed"]
IMPLEMENTATION_REVISION = "tbk_cc_30k_20yr_paired_native_v1"
IDENTITY_COLUMNS = ["case_id", "process_seed", "condition", "level", "sigma_years",
                    "member", "estimator_seed", "input_status", "input_sha256"]


def replace_output(temporary, target, timeout=5.0):
    """Retry atomic replacement only for transient Windows sharing failures."""
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
    """Presentation and run controls may change without invalidating science."""
    path = Path(__file__) if path is None else Path(path)
    tree = ast.parse(path.read_text(encoding="utf-8-sig"))
    controls = {"RUN_CAUSAL", "PREPARE_ONLY", "RESUME", "N_WORKERS"}
    nodes = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id in controls
                                              for target in node.targets):
            continue
        if isinstance(node, ast.FunctionDef) and (node.name == "plot_results" or node.name.startswith("_plot_")):
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


def execution_parameters():
    parameters = {name: value for name, value in globals().items()
                  if name.isupper() and isinstance(value, (str, int, float, list, dict))}
    parameters["packages"] = {name: version(name) for name in [
        "numpy", "scipy", "pandas", "scikit-learn", "tigramite", "pyEDM",
        "LK-Info-Flow", "statsmodels", "dcor", "threadpoolctl", "numba", "llvmlite"]}
    parameters["environment"] = {
        "system": platform.system(), "machine": platform.machine(), "python": platform.python_version(),
        "numpy_blas": np.__config__.CONFIG.get("Build Dependencies", {}).get("blas", {})}
    implementation = [ROOT / "toolbox" / name for name in
                      ["tbk_model.py", "tbk_observation.py", "observation.py", "causal.py"]]
    parameters["implementation_sha256"] = {p.relative_to(ROOT).as_posix(): file_digest(p) for p in implementation}
    parameters["implementation_sha256"][Path(__file__).relative_to(ROOT).as_posix()] = calculation_digest()
    parameters["interpretation"] = "Three independent noisy histories and chronology centers; ten members per nonzero level within each history."
    parameters["record_mapping"] = {"x": "South / T1", "y": "North / T4"}
    parameters["input_loading"] = "Default CSV parser for estimator inputs; round_trip only for result and equality checks."
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
    path = Path(path)
    if path.exists():
        previous = pd.read_csv(path, float_precision="round_trip")
        if (list(previous.columns) != list(series.columns) or previous.shape != series.shape
                or not np.array_equal(previous.to_numpy(), series.to_numpy(), equal_nan=True)):
            raise ValueError(f"Regenerated input changed: {path.name}")
        return
    write_table(series, path)


def native_snapshot(value):
    """Serialize complete native arrays/tables and scientifically relevant OLS fits."""
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
            raise ValueError("Native keys collide after conversion to text.")
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
    raise TypeError(f"Unsupported native type: {type(value).__name__}")


def save_snapshot(value, path, check_existing=False):
    """Keep the same saved prepared arrays when regenerating a resumed input."""
    path = Path(path)
    text = json.dumps(native_snapshot(value), sort_keys=True, separators=(",", ":"))
    if check_existing and path.exists():
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            previous = json.dumps(json.load(handle), sort_keys=True, separators=(",", ":"))
        if previous != text:
            raise ValueError(f"Regenerated scientific arrays changed: {path.name}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".partial.gz")
    with gzip.open(temporary, "wt", encoding="utf-8") as handle:
        handle.write(text)
    replace_output(temporary, path)


def prepare_inputs(output_dir=DATA_DIR):
    """Generate each history once; dates change while sampled proxy values stay fixed."""
    output_dir = Path(output_dir)
    jobs, process_rows, chronology_tables, relative_tables, acceptance_tables = [], [], [], [], []
    expected_grid = np.linspace(ANALYSIS_OLD, ANALYSIS_YOUNG,
                                int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / DT)) + 1)
    for process_seed in PROCESS_SEEDS:
        sample_step_years = int(round(DT * 1000))
        if MODEL_STEP_YEARS != 1 or sample_step_years != DT * 1000:
            raise ValueError("The original annual model and integer-year observation spacing are required.")
        prepared = prepare_tbk_records(process_seed, sigma_years=AGE_SIGMA_YEARS,
            n_members=N_AGE_MEMBERS, runlength=TOTAL_YEARS, spinup_years=SPINUP_YEARS,
            south_diffusion=D_SOUTH, noise_sd=PROCESS_NOISE_SD, sample_step_years=sample_step_years,
            analysis_old=ANALYSIS_OLD, analysis_young=ANALYSIS_YOUNG,
            correlation_years=AGE_CORRELATION_YEARS, max_draws=MAX_AGE_DRAWS)
        if not np.allclose(prepared["age_grid"], expected_grid, rtol=0, atol=1e-12):
            raise ValueError("TBK preparation does not use the requested fixed analysis grid.")
        expected_n = 1 + len(AGE_SIGMA_YEARS) * N_AGE_MEMBERS
        if len(prepared["cases"]) != expected_n:
            raise ValueError("Every history requires its baseline and all requested chronology members.")
        scientific_arrays = {key: value for key, value in prepared.items() if key != "cases"}
        save_snapshot(scientific_arrays, output_dir / "prepared" / f"history_{process_seed}.json.gz", check_existing=True)
        save_input(prepared["sampled"], output_dir / "prepared" / f"sampled_{process_seed}.csv")
        process_rows.append({"process_seed": process_seed, **prepared["diagnostics"]["process"]})
        for name, tables in [("chronology", chronology_tables), ("relative", relative_tables), ("acceptance", acceptance_tables)]:
            table = prepared["diagnostics"][name].copy()
            table["process_seed"] = process_seed
            tables.append(table)
        for case in prepared["cases"]:
            series = case["series"]
            if (list(series.columns) != ["time", "age", "x", "y"] or len(series) != len(expected_grid)
                    or not np.isfinite(series.to_numpy()).all()
                    or not np.allclose(series.age, expected_grid, rtol=0, atol=1e-12)
                    or not np.all(np.diff(series.time) > 0)):
                raise ValueError(f"Invalid fixed-grid input: {case['case_id']}")
            job = {key: value for key, value in case.items() if key != "series"}
            job.setdefault("condition", "accurate" if job["level"] == 0 else "chronology")
            job.setdefault("input_status", "ok")
            if job["estimator_seed"] != process_seed + METHOD_SEED_OFFSET:
                raise ValueError("Estimator seed must be shared across all chronology levels of one history.")
            path = output_dir / "inputs" / f"{job['case_id']}.csv"
            save_input(series, path)
            job["input_sha256"] = file_digest(path)
            jobs.append(job)
    if len({job["case_id"] for job in jobs}) != len(jobs):
        raise ValueError("Prepared case identifiers are not unique.")
    process_diagnostics = pd.DataFrame(process_rows)
    chronology_diagnostics = pd.concat(chronology_tables, ignore_index=True)
    write_table(pd.DataFrame(jobs), output_dir / "input_cases.csv")
    write_table(process_diagnostics, output_dir / "process_diagnostics.csv")
    write_table(chronology_diagnostics, output_dir / "chronology_diagnostics.csv")
    write_table(pd.concat(relative_tables, ignore_index=True), output_dir / "relative_error_diagnostics.csv")
    write_table(pd.concat(acceptance_tables, ignore_index=True), output_dir / "chronology_acceptance.csv")
    return jobs, process_diagnostics, chronology_diagnostics


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


def empty_edges(method, status, error=""):
    return pd.DataFrame([{"method": method, "source": source, "target": target,
        "score": np.nan, "p_value": np.nan, "edge_detected": None,
        "unoriented": False, "status": status, "error": error}
        for source, target in [("x", "y"), ("y", "x")]])


def run_case(job, output_dir, null_file):
    """One worker writes four full archives; eight directional rows are its commit."""
    output_dir = Path(output_dir)
    input_path = output_dir / "inputs" / f"{job['case_id']}.csv"
    if file_digest(input_path) != job["input_sha256"]:
        raise ValueError(f"Input changed before calculation: {job['case_id']}")
    series = pd.read_csv(input_path)
    tables = []
    with threadpool_limits(limits=1):
        for method in METHODS:
            started = time.perf_counter()
            try:
                result = call_method(series, method, job["estimator_seed"], null_file)
            except Exception as error:
                edges = empty_edges(method, "error", f"{type(error).__name__}: {error}")
                result = {"edges": edges.copy(), "raw": {}, "diagnostics": {"error": str(error)}}
            else:
                edges = result["edges"].copy()
                if "unoriented" not in edges:
                    edges["unoriented"] = False
                if (result["diagnostics"].get("undirected_contemporaneous", False)
                        or result["diagnostics"].get("conflicting_contemporaneous", False)):
                    edges["unoriented"] = True
            # Storage errors stop the run. They must not masquerade as scientific failures.
            edges["method"] = method
            # The archive retains the package-wrapper output verbatim; only the
            # separate case table uses display names and joint unresolved flags.
            native_path = output_dir / "native" / f"{job['case_id']}_{method.replace('+', 'plus')}.json.gz"
            archive = {**result, "case_identity": job, "method": method}
            save_snapshot(archive, native_path)
            if method == "CCM":
                for key, suffix in [("lag_profile", "profile"), ("lag_peaks", "peaks")]:
                    if isinstance(result["raw"].get(key), pd.DataFrame):
                        write_table(result["raw"][key], output_dir / "lag" / f"{job['case_id']}_{suffix}.csv")
            edges["seconds"] = time.perf_counter() - started
            edges["native_sha256"] = file_digest(native_path)
            for key, value in job.items():
                edges[key] = value
            tables.append(edges)
    result = pd.concat(tables, ignore_index=True)
    write_table(result, output_dir / "members" / f"{job['case_id']}_edges.csv")
    return result


def check_native_archive(archive, method, job):
    """Check complete successful CCM arrays/GC fits; keep explicitly failed outputs."""
    if archive.get("case_identity") != job or archive.get("method") != method:
        raise ValueError("Native archive identity differs from its committed case.")
    if not {"edges", "raw", "diagnostics"}.issubset(archive):
        raise ValueError("Incomplete native archive.")
    edge_data = archive["edges"]["data"]
    sources = edge_data["source"]["values"]
    targets = edge_data["target"]["values"]
    if set(zip(sources, targets)) != {("x", "y"), ("y", "x")} or len(sources) != 2:
        raise ValueError("Native archive is missing directional output.")
    if archive["diagnostics"].get("error") is not None:
        if any(status == "ok" for status in edge_data["status"]["values"]):
            raise ValueError("A failed native archive contains successful directional rows.")
        return
    raw = archive["raw"]
    if method == "PCMCI+" and not {"p_matrix", "val_matrix", "graph", "p_holm", "graph_holm"}.issubset(raw):
        raise ValueError("Native PCMCI+ archive lacks complete search and screening arrays.")
    if method == "LKIF" and not {"IF", "nIF", "p", "p_holm"}.issubset(raw):
        raise ValueError("Native LKIF archive lacks information flow and significance arrays.")
    if method == "CCM":
        expected = {"rho_null": [CCM_SURROGATES, 2], "cuts": [CCM_SURROGATES, 2],
            "library_seeds": [CCM_SURROGATES + 1, 2],
            "lag_rho": [CCM_LAG_LIBRARIES, len(CCM_LAG_STEPS), 2],
            "lag_steps": [len(CCM_LAG_STEPS)], "lag_library_seeds": [CCM_LAG_LIBRARIES]}
        for name, shape in expected.items():
            if name not in raw or raw[name].get("shape") != shape:
                raise ValueError(f"Missing or incomplete native CCM array: {name}")
            if list(np.asarray(raw[name]["values"]).shape) != shape:
                raise ValueError(f"CCM array values disagree with recorded shape: {name}")
        if raw["lag_steps"]["values"] != CCM_LAG_STEPS:
            raise ValueError("Native CCM lag grid changed.")
        if not {"lag_profile", "lag_peaks", "lag_common_anchor_indices", "observed"}.issubset(raw):
            raise ValueError("Native CCM archive lacks complete lag/observed evidence.")
        expected_n = int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / DT)) + 1
        anchors = list(range(max((CCM_E - 1) * abs(CCM_TAU), -CCM_LAG_STEPS[0]),
                             expected_n - CCM_LAG_STEPS[-1]))
        if raw["lag_common_anchor_indices"]["values"] != anchors:
            raise ValueError("Native CCM common anchors differ from the prescribed scan.")
        if len(raw["lag_profile"]["index"]["values"]) != 2 * len(CCM_LAG_STEPS):
            raise ValueError("Native CCM lag profile is incomplete.")
        if len(raw["lag_peaks"]["index"]["values"]) != 2 * CCM_LAG_LIBRARIES:
            raise ValueError("Native CCM library peak table is incomplete.")
    if method == "GC":
        required_fit = {"params", "resid", "fittedvalues", "normalized_cov_params", "bse", "tvalues",
                        "pvalues", "ssr", "scale", "nobs", "df_model", "df_resid", "cov_params",
                        "conf_int", "endog", "exog", "design_rank"}
        for direction in ["x_to_y", "y_to_x"]:
            for fit in ["restricted", "full"]:
                if not required_fit.issubset(raw.get(direction, {}).get(fit, {})):
                    raise ValueError(f"Native GC archive lacks complete {direction} {fit} fit.")
        if not {"var_coefficients", "var_intercept", "var_residuals"}.issubset(raw):
            raise ValueError("Native GC archive lacks VAR diagnostics.")


def read_completed_case(job, output_dir=DATA_DIR):
    output_dir = Path(output_dir)
    path = output_dir / "members" / f"{job['case_id']}_edges.csv"
    if not path.exists():
        return None
    result = pd.read_csv(path, float_precision="round_trip", dtype={
        "edge_detected": "boolean", "unoriented": "boolean", "raw_edge_detected": "boolean"})
    keys = ["method", "source", "target"]
    expected = {(method, source, target) for method in METHODS for source, target in [("x", "y"), ("y", "x")]}
    required = set(keys + ["status", "edge_detected", "unoriented", "score", "p_value", "native_sha256", *job])
    if not required.issubset(result) or len(result) != 8 or result.duplicated(keys).any():
        raise ValueError(f"Malformed committed case: {path.name}")
    if set(result[keys].itertuples(index=False, name=None)) != expected or result.status.isna().any() or result.status.eq("pending").any():
        raise ValueError(f"Incomplete committed directions: {path.name}")
    for key, value in job.items():
        if not result[key].eq(value).all():
            raise ValueError(f"Saved case identity differs in {key}: {path.name}")
    if file_digest(output_dir / "inputs" / f"{job['case_id']}.csv") != job["input_sha256"]:
        raise ValueError(f"Committed input bytes changed: {path.name}")
    for method in METHODS:
        native_path = output_dir / "native" / f"{job['case_id']}_{method.replace('+', 'plus')}.json.gz"
        if not native_path.is_file():
            raise ValueError(f"Committed case lacks native {method} archive: {path.name}")
        if not result.loc[result.method.eq(method), "native_sha256"].eq(file_digest(native_path)).all():
            raise ValueError(f"Committed native {method} bytes changed: {path.name}")
        with gzip.open(native_path, "rt", encoding="utf-8") as handle:
            archive = json.load(handle)
        check_native_archive(archive, method, job)
        archived_edges = pd.DataFrame({key: column["values"] for key, column in archive["edges"]["data"].items()})
        # Reconstruct the same reporting transformations without modifying the
        # archived native evidence or its complete numerical arrays.
        archived_edges["method"] = method
        if "unoriented" not in archived_edges:
            archived_edges["unoriented"] = False
        if (archive["diagnostics"].get("undirected_contemporaneous", False)
                or archive["diagnostics"].get("conflicting_contemporaneous", False)):
            archived_edges["unoriented"] = True
        committed = result.loc[result.method.eq(method)].set_index(["source", "target"])
        for _, edge in archived_edges.iterrows():
            saved = committed.loc[(edge.source, edge.target)]
            for column, value in edge.items():
                if column in ["source", "target"]:
                    continue
                if column not in committed:
                    raise ValueError(f"Committed table lacks native field {column}: {path.name}")
                actual = saved[column]
                both_missing = pd.isna(value) and pd.isna(actual)
                empty_error = column == "error" and value == "" and pd.isna(actual)
                if not both_missing and not empty_error and (pd.isna(value) or pd.isna(actual) or actual != value):
                    raise ValueError(f"Committed native field differs: {method}/{column}, {path.name}")
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


def evidence_state(forward, reverse, valid, unresolved=False):
    if not valid:
        return "failed"
    if unresolved:
        return "unresolved"
    return "both" if forward and reverse else "x_to_y" if forward else "y_to_x" if reverse else "neither"


def summarize_results(results):
    """Retain complete states, history-specific support and paired native-score changes."""
    metadata = IDENTITY_COLUMNS + ["method"]
    required = metadata + ["source", "target", "score", "status", "edge_detected", "unoriented"]
    if results.empty or not set(required).issubset(results.columns):
        raise ValueError("TBK summaries require all requested cases and explicit statuses.")
    if results.duplicated(["case_id", "method", "source", "target"]).any() or results.status.isna().any():
        raise ValueError("Directional rows must be unique with explicit status.")
    for column in ["edge_detected", "unoriented", "raw_edge_detected"]:
        if column in results and any(not isinstance(value, (bool, np.bool_)) for value in results[column].dropna()):
            raise ValueError(f"{column} must be boolean or missing.")
    state_rows = []
    for (_, method), pair in results.groupby(["case_id", "method"], sort=False):
        if len(pair) != 2 or set(zip(pair.source, pair.target)) != {("x", "y"), ("y", "x")}:
            raise ValueError("Every method/input needs both directions.")
        if any(pair[column].nunique(dropna=False) != 1 for column in metadata):
            raise ValueError("Directions disagree on case identity.")
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
        row["state"] = evidence_state(row["x_to_y"], row["y_to_x"], row["valid"], row["unresolved"])
        row["raw_state"] = (evidence_state(row["raw_x_to_y"], row["raw_y_to_x"], row["raw_valid"])
                            if method == "CCM" else None)
        row["ccm_lag_failed"] = bool(method == "CCM" and "lag_status" in pair and pair.lag_status.eq("failed").any())
        state_rows.append(row)
    states = pd.DataFrame(state_rows)
    for column in ["x_to_y", "y_to_x", "raw_x_to_y", "raw_y_to_x"]:
        states[column] = states[column].astype("boolean")
    reference_states = states.loc[states.level.eq(0), ["process_seed", "method", "state"]]
    reference_states = reference_states.rename(columns={"state": "baseline_state"})
    states = states.merge(reference_states, on=["process_seed", "method"], how="left", validate="many_to_one")
    if states.baseline_state.isna().any():
        raise ValueError("Every history/method requires its accurate-age reference.")

    group_keys = ["process_seed", "level", "sigma_years", "method"]
    summary_rows = []
    for key, group in states.groupby(group_keys, sort=False):
        expected_n = 1 if key[1] == 0 else N_AGE_MEMBERS
        if len(group) != expected_n:
            raise ValueError("Summary requires all requested members, including unavailable placeholders.")
        valid, raw_valid = group[group.valid], group[group.raw_valid]
        row = dict(zip(group_keys, key))
        row.update({"n_requested": expected_n, "n_valid": len(valid), "n_raw_valid": len(raw_valid),
                    "n_pending": int(group.pending.sum()), "n_failed": int(group.failed.sum()),
                    "n_unresolved": int(group.unresolved.sum()), "valid_fraction": len(valid) / expected_n,
                    "n_ccm_lag_failed": int(group.ccm_lag_failed.sum()), "baseline_state": group.baseline_state.iloc[0]})
        for state in STATES:
            row[f"n_state_{state}"] = int(group.state.eq(state).sum())
        row["fraction_both_requested"] = row["n_state_both"] / expected_n
        row["fraction_unoriented"] = float(valid.unresolved.mean()) if len(valid) else np.nan
        row["high_unoriented"] = bool(len(valid) and row["fraction_unoriented"] >= SUPPORT_THRESHOLD)
        for direction in ["x_to_y", "y_to_x"]:
            row[f"n_retained_{direction}"] = int(valid[direction].sum())
            row[f"q_{direction}"] = float(valid[direction].mean()) if len(valid) else np.nan
            row[f"n_raw_retained_{direction}"] = int(raw_valid[f"raw_{direction}"].sum())
            row[f"raw_q_{direction}"] = float(raw_valid[f"raw_{direction}"].mean()) if len(raw_valid) else np.nan
            # Marginal counts retain a usable direction even when the other direction is unavailable.
            row[f"marginal_retained_{direction}"] = int(group[direction].fillna(False).sum())
            row[f"marginal_fraction_{direction}"] = row[f"marginal_retained_{direction}"] / expected_n
        row["support_threshold"] = SUPPORT_THRESHOLD
        row["decision"] = "insufficient_valid_members"
        if len(valid) and row["valid_fraction"] >= MIN_VALID_FRACTION:
            forward = row["q_x_to_y"] >= SUPPORT_THRESHOLD
            reverse = row["q_y_to_x"] >= SUPPORT_THRESHOLD
            row["decision"] = ("both" if forward and reverse else "x_to_y" if forward else
                               "y_to_x" if reverse else "unoriented" if row["high_unoriented"] else "no_high_support")
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)

    fraction_rows = []
    for key, group in states.groupby(["level", "sigma_years", "method"], sort=False):
        for state in STATES:
            count = int(group.state.eq(state).sum())
            fraction_rows.append({**dict(zip(["level", "sigma_years", "method"], key)),
                "state": state, "n_requested": len(group), "n_state": count,
                "fraction_requested": count / len(group), "n_histories": group.process_seed.nunique()})

    baseline_keys = ["process_seed", "method", "source", "target"]
    reference = results.loc[results.level.eq(0), baseline_keys + ["score", "status"]]
    reference = reference.rename(columns={"score": "baseline_score", "status": "baseline_status"})
    changes = results.loc[results.level.gt(0)].merge(reference, on=baseline_keys, how="left", validate="many_to_one")
    usable = (changes.status.eq("ok") & changes.baseline_status.eq("ok")
              & np.isfinite(changes.score) & np.isfinite(changes.baseline_score))
    changes["score_change_usable"] = usable
    changes["delta_score"] = np.where(usable, changes.score - changes.baseline_score, np.nan)
    change_keys = ["process_seed", "level", "sigma_years", "method", "source", "target"]
    change_rows = []
    for key, group in changes.groupby(change_keys, sort=False):
        values = group.loc[group.score_change_usable, "delta_score"].to_numpy()
        quartiles = np.percentile(values, [25, 50, 75]) if len(values) else [np.nan] * 3
        change_rows.append({**dict(zip(change_keys, key)), "n_requested": len(group), "n_usable": len(values),
                           "q25": quartiles[0], "median": quartiles[1], "q75": quartiles[2]})
    transition_rows = []
    for key, group in states.loc[states.level.gt(0)].groupby(group_keys, sort=False):
        baseline = group.baseline_state.iloc[0]
        for state in STATES:
            count = int(group.state.eq(state).sum())
            transition_rows.append({**dict(zip(group_keys, key)), "baseline_state": baseline,
                "state": state, "n_requested": len(group), "n_transition": count,
                "fraction_requested": count / len(group),
                "baseline_both_denominator": len(group) if baseline == "both" else 0,
                "n_both_to_one": count if baseline == "both" and state in ["x_to_y", "y_to_x"] else 0})
    return {"scenario_summary": summary, "case_states": states, "state_fractions": pd.DataFrame(fraction_rows),
            "paired_score_changes": changes, "paired_change_summary": pd.DataFrame(change_rows),
            "baseline_transitions": pd.DataFrame(transition_rows),
            "ccm_lag_summary": results.loc[results.method.eq("CCM")].copy()}


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
        jobs, process_diagnostics, chronology_diagnostics = prepare_inputs(DATA_DIR)
        for job in jobs:
            result = read_completed_case(job, DATA_DIR)
            if result is not None:
                completed[job["case_id"]] = result
        if PREPARE_ONLY:
            write_status("prepared", len(completed), len(jobs), started)
            print(f"Prepared {len(jobs)} TBK inputs; {len(completed)} complete.", flush=True)
            return
        pending = [job for job in jobs if job["case_id"] not in completed]
        if RUN_CAUSAL and pending:
            effective_n = int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / DT)) + 1 - 2 * MAX_LAG
            with threadpool_limits(limits=1):
                for estimator_seed in sorted({job["estimator_seed"] for job in pending}):
                    prepare_gpdc_null(gpdc_null_path(estimator_seed), sample_sizes=[effective_n],
                                      sig_samples=GPDC_NULL_SAMPLES, seed=estimator_seed)
            print(f"TBK experiment: {len(completed)} complete, {len(pending)} pending; {N_WORKERS} workers.", flush=True)
            write_status("running", len(completed), len(jobs), started)
            context = multiprocessing.get_context("spawn")
            with ProcessPoolExecutor(max_workers=min(N_WORKERS, len(pending)), mp_context=context) as pool:
                futures = [pool.submit(run_case, job, str(DATA_DIR), str(gpdc_null_path(job["estimator_seed"]))) for job in pending]
                for future in as_completed(futures):
                    result = future.result()
                    completed[result.case_id.iloc[0]] = result
                    if len(completed) % 10 == 0 or len(completed) == len(jobs):
                        write_status("running", len(completed), len(jobs), started)
                        print(f"Completed {len(completed)}/{len(jobs)} TBK inputs.", flush=True)
        results = collect_results(jobs, completed)
        save_outputs(results, process_diagnostics)
        state = "complete" if len(completed) == len(jobs) else "partial"
        write_status(state, len(completed), len(jobs), started)
        print(f"TBK experiment {state}: {len(completed)}/{len(jobs)} committed inputs.", flush=True)
    except BaseException as error:
        if jobs:
            write_table(collect_results(jobs, completed), DATA_DIR / "member_results.csv")
        write_status("interrupted_or_failed", len(completed), len(jobs), started,
                     error=f"{type(error).__name__}: {error}")
        raise


def _plot_tbk_style(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=9)


def _plot_read_tbk_archive(path):
    """Read selected saved arrays without regenerating a trajectory or chronology."""
    import gzip
    import json

    with gzip.open(path, 'rt', encoding='utf-8') as handle:
        saved = json.load(handle)
    def array(item):
        return np.asarray(item['values'], dtype=item['dtype']).reshape(item['shape'])
    result = {'model_time': array(saved['model']['time']),
              'temperature': array(saved['model']['temperature']),
              'age_true': array(saved['sampled']['data']['age_true'])}
    for key, package in saved['chronologies'].items():
        result[f'{key}_ages'] = array(package['ages'])
        result[f'{key}_obs'] = array(package['age_obs'])
    return result


def _plot_tbk_design(data_dir, output_dir):
    from matplotlib.patches import FancyArrowPatch, Rectangle

    archives = sorted((Path(data_dir) / 'prepared').glob('history_*.json.gz'))
    if not archives:
        raise FileNotFoundError('TBK figure 00 requires the saved prepared/history_*.json.gz inputs.')
    example = _plot_read_tbk_archive(archives[0])
    age = example['age_true']
    model_age = (example['model_time'][-1] - example['model_time']) / 1000
    temperature = example['temperature']
    if temperature.ndim == 3:
        temperature = temperature[0]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10.5))
    ax = axes[0, 0]
    ax.set_axis_off()
    box_colors = ['#0072B2', '#56B4E9', '#E69F00', '#D55E00']
    for box, x in enumerate([0.05, 0.29, 0.53, 0.77], start=1):
        ax.add_patch(Rectangle((x, 0.42), 0.17, 0.19, facecolor=box_colors[box - 1], alpha=0.20,
                               edgecolor=box_colors[box - 1], lw=1.7, transform=ax.transAxes))
        ax.text(x + 0.085, 0.515, f'T{box}', ha='center', va='center', transform=ax.transAxes, fontsize=13)
        if box < 4:
            ax.add_patch(FancyArrowPatch((x + 0.17, 0.515), (x + 0.24, 0.515),
                                         arrowstyle='<->', mutation_scale=11, transform=ax.transAxes))
    ax.text(0.135, 0.34, 'South / X', ha='center', color=box_colors[0], transform=ax.transAxes)
    ax.text(0.855, 0.34, 'North / Y', ha='center', color=box_colors[3], transform=ax.transAxes)
    ax.text(0.50, 0.79, 'Shared prescribed transport q(t)', ha='center', transform=ax.transAxes, fontsize=11)
    for x in [0.375, 0.615, 0.855]:
        ax.annotate('', xy=(x, 0.62), xytext=(0.5, 0.74), xycoords='axes fraction',
                    arrowprops={'arrowstyle': '->', 'color': '0.45', 'lw': 1})
    ax.text(0.5, 0.18, 'Neighbour diffusion + global mean restoring\nProcess noise enters every box update',
            ha='center', va='center', transform=ax.transAxes, fontsize=10)
    ax.text(0.5, 0.03, 'Endpoint statistics also contain hidden-box and shared-forcing effects',
            ha='center', transform=ax.transAxes, fontsize=9)
    ax.set_title('Four-box model; endpoint records used for causal comparison', loc='left', pad=15)
    zoom = (model_age <= 15) & (model_age >= 12)
    for box in range(1, 5):
        axes[0, 1].plot(model_age[zoom], temperature[box, zoom], color=box_colors[box - 1],
                        lw=1.5 if box in (1, 4) else 0.9,
                        label=f'T{box}' + (' South' if box == 1 else ' North' if box == 4 else ''))
    axes[0, 1].set(xlim=(15, 12), xlabel='Age (kyr BP)', ylabel='Temperature (model units)',
                    title='Accurate ages: fixed two-cycle illustration')
    axes[0, 1].legend(ncol=2, fontsize=9)
    colors = ['#4477AA', '#228833', '#CCBB44', '#AA3377']
    for variable, color, label in [('x', box_colors[0], 'South'), ('y', box_colors[3], 'North')]:
        offsets = (example[f'level4_{variable}_ages'][0] - age) * 1000
        centers = (example[f'level4_{variable}_obs'] - age) * 1000
        axes[1, 0].plot(age, offsets, color=color, label=label + ', member 1')
        axes[1, 0].plot(age, centers, color=color, ls=':', lw=1.1, label=label + ', centre')
    axes[1, 0].axhline(0, color='0.5', lw=0.8)
    axes[1, 0].set(xlim=(30, 0), xlabel='Age (kyr BP)', ylabel='Age error (yr)',
                    title='Both records uncertain: example at proposal SD 200 yr')
    axes[1, 0].legend(ncol=2, fontsize=8)
    relative_rows = []
    for archive in archives:
        saved = _plot_read_tbk_archive(archive)
        true_age = saved['age_true']
        window = (true_age <= 29) & (true_age >= 1)
        for level, sigma in enumerate([20, 50, 100, 200], start=1):
            error = (saved[f'level{level}_x_ages'] - saved[f'level{level}_y_ages'])[:, window] * 1000
            for rms in np.sqrt(np.mean(error ** 2, axis=1)):
                relative_rows.append((level, rms))
    for level, sigma in enumerate([20, 50, 100, 200], start=1):
        values = np.array([value for group, value in relative_rows if group == level])
        axes[1, 1].scatter(level + np.linspace(-0.16, 0.16, len(values)), values,
                           color=colors[level - 1], s=25, alpha=0.65)
        axes[1, 1].plot([level - 0.22, level + 0.22], [2 * sigma, 2 * sigma], color='0.3', ls='--')
        axes[1, 1].text(level, -0.11, f'{sigma}\n(n={len(values)})', ha='center',
                        transform=axes[1, 1].get_xaxis_transform(), fontsize=9)
    axes[1, 1].set(xticks=[], xlabel='Per-layer proposal SD (yr)', ylabel='Inter-record RMS age error (yr)',
                    title='Realized mismatch across 3 histories × 10 members')
    axes[1, 1].xaxis.labelpad = 36
    axes[1, 1].set_ylim(bottom=0)
    axes[1, 1].text(0.02, 0.96, 'Dashed: idealized pointwise SD = 2σ\nDots: actual 29–1 kyr temporal RMS',
                    transform=axes[1, 1].transAxes, va='top', fontsize=9)
    for ax in axes.flat:
        _plot_tbk_style(ax)
    fig.suptitle('TBK: fixed dynamics and values, increasingly uncertain endpoint dates', fontsize=15, y=0.98)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.91, bottom=0.11, hspace=0.39, wspace=0.28)
    save_figure(fig, output_dir, '00_model_and_records', axes=axes)
    plt.close(fig)


def _plot_tbk_states(tables, output_dir):
    from matplotlib.patches import Patch
    from toolbox.plotting import add_panel_labels

    states = tables['case_states']
    methods = ['PCMCI+', 'LKIF', 'GC', 'CCM']
    state_names = ['both', 'x_to_y', 'y_to_x', 'neither', 'unresolved', 'failed']
    state_labels = ['Both directions', 'South → North only', 'North → South only', 'Neither', 'Unresolved', 'Failed / unavailable']
    colors = ['#4477AA', '#EE7733', '#AA3377', '#CCCCCC', '#DDCC77', '#333333']
    fig = plt.figure(figsize=(14, 13))
    grid = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.95], left=0.08, right=0.98,
                            bottom=0.10, top=0.87, hspace=0.43, wspace=0.25)
    axes = [fig.add_subplot(grid[i // 2, i % 2]) for i in range(4)]
    comparison = fig.add_subplot(grid[2, :])
    method_colors = ['#0072B2', '#D55E00', '#009E73', '#CC79A7']
    for ax, method, color in zip(axes, methods, method_colors):
        subset = states[states.method.eq(method)]
        denominators = []
        both_fraction = []
        bottom = np.zeros(5)
        for state, label, state_color in zip(state_names, state_labels, colors):
            heights, counts = [], []
            for level in range(5):
                group = subset[subset.level.eq(level)]
                counts.append(int(group.state.eq(state).sum()))
                heights.append(100 * counts[-1] / len(group) if len(group) else 0)
            ax.bar(range(5), heights, bottom=bottom, color=state_color, edgecolor='white', linewidth=0.7)
            for level, (height, count) in enumerate(zip(heights, counts)):
                if height >= 11:
                    ax.text(level, bottom[level] + height / 2, str(count), ha='center', va='center',
                            fontsize=10, color='white' if state in ('both', 'y_to_x', 'failed') else 'black')
            bottom = bottom + np.array(heights)
        for level in range(5):
            group = subset[subset.level.eq(level)]
            denominators.append(len(group))
            both_fraction.append(100 * group.state.eq('both').mean() if len(group) else np.nan)
        ax.set(ylim=(0, 100), ylabel='Requested inputs (%)', title=method,
               xticks=range(5), xticklabels=[f'{label}\nn={n}' for label, n in zip(['Accurate', '20', '50', '100', '200'], denominators)])
        if method == 'PCMCI+' and subset.unresolved.any():
            requested = len(subset)
            unresolved = int(subset.unresolved.sum())
            south_to_north = int(subset.x_to_y.astype('boolean').fillna(False).sum())
            north_to_south = int(subset.y_to_x.astype('boolean').fillna(False).sum())
            note = (f'Unresolved contemporaneous links: {unresolved}/{requested}\n'
                    f'Retained lagged South → North: {south_to_north}/{requested}\n'
                    f'Retained lagged North → South: {north_to_south}/{requested}')
            ax.text(0.04, 0.96, note, transform=ax.transAxes, va='top', fontsize=8,
                    bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.88, 'pad': 3})
        comparison.plot(range(5), both_fraction, color=color, marker='o', label=method, lw=2)
        _plot_tbk_style(ax)
    comparison.set(ylim=(-3, 103), xticks=range(5), xticklabels=['Accurate (n=3)', '20 (n=30)', '50 (n=30)', '100 (n=30)', '200 (n=30)'],
                    xlabel='Per-layer proposal SD (yr)', ylabel='Bidirectional evidence retained (%)',
                    title='Both directions retained without unresolved evidence')
    comparison.grid(axis='y', alpha=0.2)
    comparison.legend(ncol=4, loc='best')
    _plot_tbk_style(comparison)
    fig.legend(handles=[Patch(facecolor=c, label=l) for c, l in zip(colors, state_labels)],
                loc='upper center', bbox_to_anchor=(0.53, 0.945), ncol=3, frameon=False, fontsize=10)
    fig.suptitle('How dating uncertainty changes directional evidence', fontsize=16, y=0.98)
    fig.text(0.08, 0.035, 'Numbers in bars are counts. Accurate: 3 histories; each nonzero level: 3 histories × 10 members.\n'
             'Unresolved has priority in the six-state display and may coexist with a retained arrow. CCM uses its lag screen.\n'
             'Bidirectional evidence retention is a descriptive statistic, not network accuracy or a feedback-detection probability.', fontsize=10)
    add_panel_labels(axes + [comparison], x=-0.07)
    save_figure(fig, output_dir, '01_evidence_vs_uncertainty')
    plt.close(fig)


def _plot_tbk_score_changes(tables, output_dir):
    changes = tables['paired_score_changes']
    methods = ['PCMCI+', 'LKIF', 'GC', 'CCM']
    score_labels = {'PCMCI+': 'Δ max GPDC',
                    'LKIF': 'Δ max |normalized IF|',
                    'GC': 'Δ log residual variance ratio',
                    'CCM': 'Δ Tp=0 reconstruction skill'}
    colors = ['#4477AA', '#228833', '#CCBB44', '#AA3377']
    fig, axes = plt.subplots(2, 4, figsize=(17.5, 9.5), sharex=False)
    for row, (source, target, direction) in enumerate([('x', 'y', 'South → North'), ('y', 'x', 'North → South')]):
        for column, method in enumerate(methods):
            ax = axes[row, column]
            subset = changes[changes.method.eq(method) & changes.source.eq(source) & changes.target.eq(target)]
            ticklabels = []
            for level, (sigma, color) in enumerate(zip([20, 50, 100, 200], colors), start=1):
                group = subset[subset.level.eq(level)]
                usable = group[group.score_change_usable.astype(bool)]
                values = usable.delta_score.to_numpy(dtype=float)
                values = values[np.isfinite(values)]
                if len(values):
                    ax.scatter(level + np.linspace(-0.16, 0.16, len(values)), values, s=20, color=color, alpha=0.5)
                    q25, median, q75 = np.percentile(values, [25, 50, 75])
                    ax.plot([level, level], [q25, q75], color='black', lw=2)
                    ax.plot([level - 0.17, level + 0.17], [median, median], color='black', lw=2)
                ticklabels.append(f'{sigma}\n{len(values)}/{len(group)}')
            ax.axhline(0, color='0.4', ls='--', lw=1)
            ax.set(xticks=range(1, 5), xticklabels=ticklabels, ylabel=score_labels[method],
                   title=f'{method} | {direction}')
            ax.tick_params(labelbottom=True)
            if row == 1:
                ax.set_xlabel('Proposal SD (yr)\nusable / requested')
            _plot_tbk_style(ax)
    fig.suptitle('Each uncertain-date score minus the same history’s accurate-date score', fontsize=15, y=0.98)
    fig.text(0.07, 0.025, 'Dots: all finite paired member scores, including nonsignificant tests. Black: median and interquartile range.\n'
             'Three histories × ten members per level; spreads are descriptive. Score definitions differ across methods. CCM uses raw Tp=0 reconstruction skill ρ.', fontsize=10)
    fig.subplots_adjust(left=0.065, right=0.985, top=0.90, bottom=0.14, hspace=0.35, wspace=0.43)
    save_figure(fig, output_dir, '02_paired_score_changes', axes=axes)
    plt.close(fig)


def _plot_tbk_ccm(tables, data_dir, output_dir):
    from matplotlib.patches import Patch

    lag = tables['ccm_lag_summary'].copy()
    # A fully unavailable CCM table has no native lag fields; keep every
    # requested case in the figure instead of failing during presentation.
    for column in ['raw_edge_detected', 'edge_detected']:
        if column not in lag:
            lag[column] = pd.Series(pd.NA, index=lag.index, dtype='boolean')
    if 'lag_status' not in lag:
        lag['lag_status'] = 'unavailable'
    levels = ['Accurate', '20', '50', '100', '200']
    level_colors = ['#333333', '#4477AA', '#228833', '#CCBB44', '#AA3377']
    profiles = []
    metadata = lag[['case_id', 'level']].drop_duplicates()
    for item in metadata.to_dict('records'):
        path = Path(data_dir) / 'lag' / f"{item['case_id']}_profile.csv"
        if path.exists():
            profile = pd.read_csv(path, float_precision='round_trip')
            profile['case_id'] = item['case_id']
            profile['level'] = item['level']
            profiles.append(profile)
    profiles = pd.concat(profiles, ignore_index=True) if profiles else pd.DataFrame(
        columns=['source', 'target', 'level', 'lag_kyr', 'rho_mean'])
    fig, axes = plt.subplots(3, 2, figsize=(14, 14.5))
    diagnostic_kinds = ['past', 'future', 'boundary', 'zero', 'tied', 'failed']
    diagnostic_colors = ['#4477AA', '#EE7733', '#AA3377', '#CCCCCC', '#DDCC77', '#333333']
    for column, (source, target, direction) in enumerate([('x', 'y', 'South → North'), ('y', 'x', 'North → South')]):
        subset = lag[lag.source.eq(source) & lag.target.eq(target)]
        diagnostic = np.zeros((6, 5))
        denominator_labels = []
        for key, label, color, linestyle in [
                ('raw_edge_detected', 'Raw Tp=0 significance', '#333333', '--'),
                ('edge_detected', 'Retained after lag screen', '#0072B2', '-'),
                ('screen_unresolved', 'Raw significant + lag unresolved', '#AA3377', ':'),
                ('screen_rejected', 'Raw significant + lag incompatible', '#D55E00', '-.')]:
            fractions = []
            for level in range(5):
                group = subset[subset.level.eq(level)]
                raw = group.raw_edge_detected.astype('boolean').fillna(False)
                if key == 'screen_unresolved':
                    values = raw & group.lag_status.eq('unresolved')
                elif key == 'screen_rejected':
                    values = raw & group.lag_status.eq('incompatible')
                else:
                    values = group[key].astype('boolean').fillna(False)
                fractions.append(100 * values.sum() / len(group) if len(group) else np.nan)
            axes[0, column].plot(range(5), fractions, color=color, ls=linestyle, marker='o', label=label)
        axes[0, column].set(ylim=(-3, 103), ylabel='Requested directional inputs (%)',
                            title=direction, xticks=range(5), xticklabels=levels,
                            xlabel='Per-layer proposal SD (yr)')
        axes[0, column].legend(fontsize=8, loc='best')
        for level, color in enumerate(level_colors):
            curves = profiles[profiles.source.eq(source) & profiles.target.eq(target) & profiles.level.eq(level)]
            if len(curves):
                summary = curves.groupby('lag_kyr').rho_mean.agg(['median', 'count'])
                axes[1, column].plot(summary.index.to_numpy() * 1000, summary['median'], color=color,
                                     label=f'{levels[level]} (n={curves.case_id.nunique()})', lw=1.8)
            group = subset[subset.level.eq(level)]
            denominator_labels.append(f'{levels[level]}\nn={len(group)}')
            for index, kind in enumerate(diagnostic_kinds):
                column_name = f'lag_{kind}_fraction'
                values = pd.to_numeric(group[column_name], errors='coerce') if column_name in group else pd.Series(np.nan, index=group.index)
                diagnostic[index, level] = values.fillna(1.0 if kind == 'failed' else 0.0).sum() / len(group) if len(group) else 0
        axes[1, column].axvline(0, color='0.5', ls=':', lw=1)
        axes[1, column].set(xlim=(-600, 600), xlabel='Target shift (yr)', ylabel='Mean-library reconstruction skill ρ',
                            title='Median of saved case lag profiles')
        if axes[1, column].get_legend_handles_labels()[0]:
            axes[1, column].legend(ncol=2, fontsize=8, loc='center left',
                                      bbox_to_anchor=(0.01, 0.30))
        else:
            axes[1, column].text(0.5, 0.5, 'No available lag profiles',
                                 ha='center', transform=axes[1, column].transAxes)
        bottom = np.zeros(5)
        for index, (kind, color) in enumerate(zip(diagnostic_kinds, diagnostic_colors)):
            axes[2, column].bar(range(5), 100 * diagnostic[index], bottom=bottom, color=color, edgecolor='white', lw=0.7)
            bottom = bottom + 100 * diagnostic[index]
        axes[2, column].set(ylim=(0, 100), xticks=range(5), xticklabels=denominator_labels,
                            ylabel='Mean library peak-category fraction (%)',
                            title='All 100 library peaks per input, regardless of raw significance',
                            xlabel='Per-layer proposal SD (yr)')
    for ax in axes.flat:
        _plot_tbk_style(ax)
    fig.legend(handles=[Patch(facecolor=c, label=k.capitalize()) for k, c in zip(diagnostic_kinds, diagnostic_colors)],
               loc='lower center', bbox_to_anchor=(0.53, 0.064), ncol=6, frameon=False)
    fig.suptitle('CCM: original significance, lag-screen decisions and timing evidence', fontsize=15, y=0.985)
    fig.text(0.08, 0.017, 'Accurate: 3 histories; each nonzero level: 30 member inputs from the same 3 histories.\n'
             'Raw significance and screened decisions are distinct; a lag-rejected or unresolved arrow does not establish absence of physical feedback.\n'
             'Profiles average 100 paired L=300 libraries within each case, then take the median across cases. No causal probability or confidence interval is implied.', fontsize=9)
    fig.subplots_adjust(left=0.085, right=0.98, top=0.93, bottom=0.14, hspace=0.40, wspace=0.25)
    save_figure(fig, output_dir, '03_ccm_evidence_and_timing', axes=axes)
    plt.close(fig)


def plot_results(tables, process_diagnostics, output_dir=None, data_dir=None):
    """Render the saved inputs/results only; no model or causal calls occur here."""
    output_dir = FIGURE_DIR if output_dir is None else Path(output_dir)
    data_dir = DATA_DIR if data_dir is None else Path(data_dir)
    _plot_tbk_design(data_dir, output_dir)
    _plot_tbk_states(tables, output_dir)
    _plot_tbk_score_changes(tables, output_dir)
    _plot_tbk_ccm(tables, data_dir, output_dir)



if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
