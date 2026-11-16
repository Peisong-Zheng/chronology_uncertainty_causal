"""Accurate-age TBK excitation diagnostic with fixed physical memory.

Three paired histories compare the original forcing response at three process
noise amplitudes with its known-forcing-subtracted random response. The twelve
inputs retain complete native evidence from four methods. Three original inputs
also receive a separate PCMCI+/ParCorr diagnostic, without replacing GPDC.
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
import gzip
import importlib
from io import StringIO
import json
import multiprocessing
import shutil
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits

from toolbox import causal
from toolbox.plotting import save_figure
from toolbox.tbk_model import simulate_tbk

core = importlib.import_module("04_check_tbk_chronology")

RUN_CAUSAL = True             # Run the authorized diagnostic after validated preflight.
PREPARE_ONLY = False          # False completes all native fits and saved-result figures.
RESUME = True                 # Keep complete cases, including unfavorable evidence.
N_WORKERS = 32                # One numerical thread; one native method per pool task.
PROCESS_SEEDS = [20261009, 20261301, 20261302]
TOTAL_YEARS = 39000           # Annual integration, unchanged from completed 04.
SPINUP_YEARS = 9000           # Discarded startup; retained support is 30 kyr.
MODEL_STEP_YEARS = 1
D_SOUTH = 0.005               # Original cell-12 southern diffusion coefficient.
DT = 0.020                   # Sampling and analysis spacing, kyr.
ANALYSIS_OLD = 29.0
ANALYSIS_YOUNG = 1.0
METHOD_SEED_OFFSET = 400
MAX_LAG = 5                  # PCMCI+/LKIF: 100 years at the fixed sampling step.
GC_ORDER = 1                 # Fixed twenty-year history; report residual checks.
ALPHA = 0.05
GPDC_NULL_SAMPLES = 10000
GPDC_GP_PARAMS = {"optimizer": None}
CCM_E = 3
CCM_TAU = -1
CCM_EXCLUSION_RADIUS = 10
CCM_LIBRARY_SIZE = 300
CCM_SAMPLES = 10
CCM_SURROGATES = 100
CCM_LAG_STEPS = list(range(-30, 31))
CCM_LAG_LIBRARIES = 100
CCM_LAG_SIGN_THRESHOLD = 0.90
SCENARIOS = [
    {"scenario": "baseline", "noise_sd": 0.005, "label": "Original (SD 0.005)"},
    {"scenario": "noise_x10", "noise_sd": 0.05, "label": "Tenfold noise (SD 0.05)"},
    {"scenario": "noise_x100", "noise_sd": 0.5, "label": "Hundredfold noise (SD 0.5)"},
    {"scenario": "random_response", "noise_sd": 0.005, "label": "Known-forcing-subtracted response"},
]
METHODS = ["PCMCI+", "LKIF", "CCM", "GC"]
PARCORR_METHOD = "PCMCI+_ParCorr"
STEM = Path(__file__).stem
DATA_DIR = ROOT / "data/processed" / STEM
FIGURE_DIR = ROOT / "figures" / STEM
SOURCE_DIR = ROOT / "data/processed/04_check_tbk_chronology"
CACHE_DIR = ROOT / "data/cache/gpdc"
IMPLEMENTATION_REVISION = "tbk_excitation_accurate_paired_v1"


def source_job(case_id):
    """Read the exact original case identity retained in copied native archives."""
    cases = pd.read_csv(SOURCE_DIR / "input_cases.csv", float_precision="round_trip")
    selected = cases[cases.case_id.eq(case_id)]
    if len(selected) != 1 or selected.iloc[0].level != 0:
        raise ValueError(f"Missing accurate-age source case: {case_id}")
    return selected.iloc[0].to_dict()


def check_source_parameters():
    """Require the completed source to match the current native implementation."""
    previous = json.loads((SOURCE_DIR / "parameters.json").read_text(encoding="utf-8"))
    current = core.execution_parameters()
    ignored = {"RUN_CAUSAL", "PREPARE_ONLY", "RESUME", "N_WORKERS"}
    changed = [key for key in set(previous) | set(current)
               if key not in ignored and previous.get(key) != current.get(key)]
    if changed:
        raise ValueError(f"Source 04 settings/implementation differ: {sorted(changed)}")
    names = ["PROCESS_SEEDS", "TOTAL_YEARS", "SPINUP_YEARS", "MODEL_STEP_YEARS",
             "D_SOUTH", "DT", "ANALYSIS_OLD", "ANALYSIS_YOUNG", "METHOD_SEED_OFFSET",
             "MAX_LAG", "GC_ORDER", "ALPHA", "GPDC_NULL_SAMPLES", "GPDC_GP_PARAMS",
             "CCM_E", "CCM_TAU", "CCM_EXCLUSION_RADIUS", "CCM_LIBRARY_SIZE", "CCM_SAMPLES",
             "CCM_SURROGATES", "CCM_LAG_STEPS", "CCM_LAG_LIBRARIES", "CCM_LAG_SIGN_THRESHOLD"]
    if any(previous[name] != globals()[name] for name in names):
        raise ValueError("Excitation diagnostic must retain source 04 method/grid settings.")
    if previous["PROCESS_NOISE_SD"] != SCENARIOS[0]["noise_sd"]:
        raise ValueError("The reusable original-noise amplitude changed.")
    return previous


def source_hashes():
    paths = [SOURCE_DIR / "parameters.json", SOURCE_DIR / "input_cases.csv"]
    for seed in PROCESS_SEEDS:
        case_id = f"tbk_s{seed}_accurate"
        paths.extend([SOURCE_DIR / "inputs" / f"{case_id}.csv",
                      SOURCE_DIR / "members" / f"{case_id}_edges.csv"])
        paths.extend(SOURCE_DIR / "native" / f"{case_id}_{method.replace('+', 'plus')}.json.gz"
                     for method in METHODS)
        paths.extend(SOURCE_DIR / "lag" / f"{case_id}_{suffix}.csv" for suffix in ["profile", "peaks"])
    return {p.relative_to(ROOT).as_posix(): core.file_digest(p) for p in paths}


def execution_parameters():
    source = check_source_parameters()
    parameters = {name: value for name, value in globals().items()
                  if name.isupper() and isinstance(value, (str, int, float, list, dict))}
    parameters["packages"] = source["packages"]
    parameters["environment"] = source["environment"]
    parameters["source_sha256"] = source_hashes()
    parameters["implementation_sha256"] = {
        "04_check_tbk_chronology.py": core.calculation_digest(),
        "04_1_check_tbk_excitation.py": core.calculation_digest(Path(__file__)),
        "toolbox/causal.py": core.file_digest(ROOT / "toolbox/causal.py"),
        "toolbox/tbk_model.py": core.file_digest(ROOT / "toolbox/tbk_model.py"),
    }
    parameters["interpretation"] = (
        "Twelve accurate-age inputs from three paired noise histories, not chronology members. "
        "The SD 0.5 and exact forcing subtraction are diagnostic controls, not calibrated observations.")
    parameters["input_loading"] = "Default CSV parser for estimator inputs; round_trip for saved-result audit."
    return json.loads(json.dumps(parameters))


def check_parameters(parameters, output_dir=DATA_DIR):
    path = Path(output_dir) / "parameters.json"
    if path.exists():
        if not RESUME and RUN_CAUSAL:
            raise ValueError("Preserve existing outputs before a new experiment.")
        previous = json.loads(path.read_text(encoding="utf-8"))
        ignored = {"RUN_CAUSAL", "PREPARE_ONLY", "RESUME", "N_WORKERS"}
        changed = [key for key in set(previous) | set(parameters)
                   if key not in ignored and previous.get(key) != parameters.get(key)]
        if changed:
            raise ValueError(f"Saved excitation experiment differs in {sorted(changed)}")
    elif any((Path(output_dir) / "members").glob("*_edges.csv")):
        raise ValueError("Committed excitation results lack scientific parameters.")


def save_arrays(arrays, path):
    path = Path(path)
    if path.exists():
        with np.load(path, allow_pickle=False) as previous:
            if set(previous.files) != set(arrays):
                raise ValueError(f"Prepared array fields changed: {path.name}")
            for key, value in arrays.items():
                np.testing.assert_array_equal(previous[key], value, err_msg=f"{path.name}/{key}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".partial.npz")
    np.savez_compressed(temporary, **arrays)
    core.replace_output(temporary, path)


def acf_at(values, steps):
    values = np.asarray(values, dtype=float)
    if steps == 0:
        return 1.0
    if len(values) <= steps or np.std(values) == 0:
        return np.nan
    return float(np.corrcoef(values[:-steps], values[steps:])[0, 1])


def prepare_inputs(output_dir=DATA_DIR):
    """Use fixed standard-normal draws across amplitudes; never draw age errors."""
    check_source_parameters()
    output_dir = Path(output_dir)
    jobs, diagnostic_rows = [], []
    sampled_time = np.arange(SPINUP_YEARS, TOTAL_YEARS + 1, int(DT * 1000))
    indices = sampled_time - 1
    age_true = (TOTAL_YEARS - sampled_time) / 1000.
    selected = (age_true <= ANALYSIS_OLD) & (age_true >= ANALYSIS_YOUNG)
    grid = age_true[selected]
    for seed in PROCESS_SEEDS:
        clean = simulate_tbk(runlength=TOTAL_YEARS, south_diffusion=D_SOUTH, noise_sd=0, seed=seed)
        models = {item["scenario"]: simulate_tbk(runlength=TOTAL_YEARS, south_diffusion=D_SOUTH,
                   noise_sd=item["noise_sd"], seed=seed) for item in SCENARIOS[:3]}
        temperature = {name: model["temperature"][0, 1:5] for name, model in models.items()}
        clean_temperature = clean["temperature"][0, 1:5]
        temperature["random_response"] = temperature["baseline"] - clean_temperature
        for name, multiplier in [("noise_x10", 10), ("noise_x100", 100)]:
            np.testing.assert_allclose(models[name]["noise"], models["baseline"]["noise"] * multiplier,
                                       rtol=1e-14, atol=1e-15)
            np.testing.assert_array_equal(models[name]["transport_actual"], clean["transport_actual"])
        arrays = {"time": clean["time"], "age": (TOTAL_YEARS - clean["time"]) / 1000.,
                  "transport": clean["transport_actual"], "clean_temperature": clean_temperature,
                  "sampled_time": sampled_time, "age_true": age_true}
        arrays.update({f"{name}_temperature": values for name, values in temperature.items()})
        save_arrays(arrays, output_dir / "prepared" / f"history_{seed}.npz")
        for item in SCENARIOS:
            scenario = item["scenario"]
            case_id = f"tbk_s{seed}_accurate" if scenario == "baseline" else f"tbk_s{seed}_{scenario}"
            sample = temperature[scenario][:, indices]
            series = causal.resample_pair(age_true, sample[0], age_true, sample[3], age_grid=grid)["series"]
            if len(series) != 1401 or not np.isfinite(series.to_numpy()).all():
                raise ValueError("The diagnostic requires 1401 finite fixed-grid observations.")
            input_path = output_dir / "inputs" / f"{case_id}.csv"
            if scenario == "baseline":
                original_path = SOURCE_DIR / "inputs" / f"{case_id}.csv"
                original_exact = pd.read_csv(original_path, float_precision="round_trip")
                pd.testing.assert_frame_equal(series, original_exact, check_exact=True)
                expected_parsed = pd.read_csv(StringIO(series.to_csv(index=False)))
                pd.testing.assert_frame_equal(expected_parsed, pd.read_csv(original_path), check_exact=True)
                copy_verified(original_path, input_path)
                job = source_job(case_id)
            else:
                core.save_input(series, input_path)
                job = dict(case_id=case_id, process_seed=seed, condition=scenario,
                    level=0, sigma_years=0., member=-1, estimator_seed=seed + METHOD_SEED_OFFSET,
                    input_status="ok", input_sha256=core.file_digest(input_path))
            job.update(scenario=scenario, noise_sd=item["noise_sd"])
            jobs.append(job)
            for record, box in [("x", 0), ("y", 3)]:
                values = sample[box, selected]
                clean_values = clean_temperature[box, indices[selected]]
                random_values = (temperature[scenario][box, indices[selected]] - clean_values
                                 if scenario != "random_response" else values)
                diagnostic_rows.append(dict(case_id=case_id, scenario=scenario, process_seed=seed,
                    noise_sd=item["noise_sd"], record=record, n_samples=len(values),
                    sd=float(np.std(values)), clean_sd=float(np.std(clean_values)),
                    random_response_sd=float(np.std(random_values)),
                    random_to_total_variance=float(np.var(random_values) / np.var(values)),
                    acf_20yr=acf_at(values, 1), acf_100yr=acf_at(values, 5), acf_1500yr=acf_at(values, 75),
                    random_acf_20yr=acf_at(random_values, 1), random_acf_1500yr=acf_at(random_values, 75)))
    if len(jobs) != 12 or len({job["case_id"] for job in jobs}) != 12:
        raise ValueError("Exactly twelve distinct diagnostic inputs are required.")
    diagnostics = pd.DataFrame(diagnostic_rows)
    core.write_table(pd.DataFrame(jobs), output_dir / "input_cases.csv")
    core.write_table(diagnostics, output_dir / "input_diagnostics.csv")
    return jobs, diagnostics


def copy_verified(source, target):
    source, target = Path(source), Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if core.file_digest(target) != core.file_digest(source):
            raise ValueError(f"Existing reused output differs: {target.name}")
    else:
        shutil.copy2(source, target)


def restore_native(value):
    if isinstance(value, dict) and set(value) == {"dtype", "shape", "values"}:
        result = np.asarray(value["values"], dtype=value["dtype"])
        if list(result.shape) != value["shape"]:
            raise ValueError("Native array shape is incomplete.")
        return result
    if isinstance(value, dict) and set(value) == {"columns", "index", "data"}:
        return pd.DataFrame({key: restore_native(value["data"][str(key)]) for key in value["columns"]},
                            index=restore_native(value["index"]))
    if isinstance(value, dict):
        return {key: restore_native(item) for key, item in value.items()}
    if isinstance(value, list):
        return [restore_native(item) for item in value]
    return value


def native_path(case_id, method, output_dir=DATA_DIR):
    return Path(output_dir) / "native" / f"{case_id}_{method.replace('+', 'plus')}.json.gz"


def load_native(case_id, method, data_dir=DATA_DIR):
    with gzip.open(native_path(case_id, method, data_dir), "rt", encoding="utf-8") as handle:
        return restore_native(json.load(handle))


def validate_native(job, method, output_dir=DATA_DIR):
    path = native_path(job["case_id"], method, output_dir)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        archive = json.load(handle)
    identity = source_job(job["case_id"]) if job["scenario"] == "baseline" and method != PARCORR_METHOD else job
    if method == PARCORR_METHOD:
        core.check_native_archive({**archive, "method": "PCMCI+"}, "PCMCI+", identity)
        if archive["method"] != PARCORR_METHOD or archive["diagnostics"].get("conditional_test") != "ParCorr":
            raise ValueError("ParCorr archive is missing its explicit test identity.")
    else:
        core.check_native_archive(archive, method, identity)
    return restore_native(archive)


def report_edges(result, method, job, archive_path):
    edges = result["edges"].copy()
    edges["method"] = method
    if "unoriented" not in edges:
        edges["unoriented"] = False
    if result["diagnostics"].get("undirected_contemporaneous", False) or result["diagnostics"].get("conflicting_contemporaneous", False):
        edges["unoriented"] = True
    for key, value in job.items():
        edges[key] = value
    edges["native_sha256"] = core.file_digest(archive_path)
    return edges


def reuse_baselines(jobs, output_dir=DATA_DIR):
    """Import complete original native evidence without re-fitting an estimator."""
    check_source_parameters()
    output_dir = Path(output_dir)
    completed = {}
    for job in jobs:
        if job["scenario"] != "baseline":
            continue
        original_job = source_job(job["case_id"])
        original = core.read_completed_case(original_job, SOURCE_DIR)
        if original is None:
            raise ValueError("Source baseline is not fully committed.")
        pd.testing.assert_frame_equal(pd.read_csv(output_dir / "inputs" / f"{job['case_id']}.csv"),
            pd.read_csv(SOURCE_DIR / "inputs" / f"{job['case_id']}.csv"), check_exact=True)
        for method in METHODS:
            copy_verified(native_path(job["case_id"], method, SOURCE_DIR), native_path(job["case_id"], method, output_dir))
        for suffix in ["profile", "peaks"]:
            copy_verified(SOURCE_DIR / "lag" / f"{job['case_id']}_{suffix}.csv",
                          output_dir / "lag" / f"{job['case_id']}_{suffix}.csv")
        existing = read_completed_case(job, output_dir)
        if existing is None:
            result = original.copy()
            for key, value in job.items():
                result[key] = value
            core.write_table(result, output_dir / "members" / f"{job['case_id']}_edges.csv")
            existing = read_completed_case(job, output_dir)
        completed[job["case_id"]] = existing
    return completed


def run_parcorr(series, estimator_seed):
    """Native linear CI diagnostic; same search, Holm families and graph rules."""
    from tigramite import data_processing
    from tigramite.independence_tests.parcorr import ParCorr
    from tigramite.pcmci import PCMCI

    values, dt = causal._series_values(series)
    test = ParCorr(significance="analytic", seed=int(estimator_seed))
    dataframe = data_processing.DataFrame(values, var_names=["x", "y"])
    pcmci = PCMCI(dataframe=dataframe, cond_ind_test=test, verbosity=0)
    with warnings.catch_warnings(record=True) as fit_warnings:
        warnings.simplefilter("always")
        result = pcmci.run_pcmciplus(tau_min=0, tau_max=MAX_LAG, pc_alpha=ALPHA, fdr_method="none")
    p_matrix, values = np.asarray(result["p_matrix"]), np.asarray(result["val_matrix"])
    graph = np.asarray(result["graph"])
    p_holm = np.full(p_matrix.shape, np.nan)
    graph_holm = graph.copy()
    rows = []
    for source, target in [(0, 1), (1, 0)]:
        p, score = p_matrix[source, target], values[source, target]
        valid = np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all() and np.isfinite(score).all()
        if valid:
            adjusted = multipletests(p, alpha=ALPHA, method="holm")[1]
            p_holm[source, target] = adjusted
            detected = bool(np.any((graph[source, target] == "-->") & (adjusted < ALPHA)))
            graph_holm[source, target, adjusted >= ALPHA] = ""
            magnitude, p_summary = float(np.max(np.abs(score))), float(np.min(adjusted))
        else:
            detected, magnitude, p_summary = False, np.nan, np.nan
            graph_holm[source, target] = ""
        rows.append(causal._edge_row("pcmci+_parcorr", ["x", "y"][source], ["x", "y"][target],
            magnitude, "max_abs_partial_correlation", detected, p_summary,
            "native_analytic_holm_within_direction"))
    zero_links = [str(graph_holm[0, 1, 0]), str(graph_holm[1, 0, 0])]
    return {"edges": pd.DataFrame(rows), "raw": {**result, "p_holm": p_holm, "graph_holm": graph_holm},
        "diagnostics": {"conditional_test": "ParCorr", "tau_max": MAX_LAG, "pc_alpha": ALPHA,
            "seed": int(estimator_seed), "package_version": core.execution_parameters()["packages"]["tigramite"],
            "significance": "native_analytic",
            "alpha": ALPHA, "lag_steps": list(range(MAX_LAG + 1)),
            "lag_kyr": (np.arange(MAX_LAG + 1) * dt).tolist(),
            "undirected_contemporaneous": "o-o" in zero_links,
            "conflicting_contemporaneous": "x-x" in zero_links,
            "native_contemporaneous_link": str(graph[0, 1, 0]),
            "n_fit_warnings": len(fit_warnings), "fit_warnings": sorted(set(str(w.message) for w in fit_warnings))}}


def run_method_job(job, method, output_dir, null_file=None):
    """One spawned process fits and saves one complete native method result."""
    output_dir = Path(output_dir)
    path = output_dir / "inputs" / f"{job['case_id']}.csv"
    if core.file_digest(path) != job["input_sha256"]:
        raise ValueError("Input bytes changed before calculation.")
    series = pd.read_csv(path)
    started = time.perf_counter()
    with threadpool_limits(limits=1):
        try:
            result = run_parcorr(series, job["estimator_seed"]) if method == PARCORR_METHOD else core.call_method(
                series, method, job["estimator_seed"], null_file)
        except Exception as error:
            result = {"edges": core.empty_edges(method, "error", f"{type(error).__name__}: {error}"),
                "raw": {}, "diagnostics": {"error": str(error)}}
            if method == PARCORR_METHOD:
                result["diagnostics"]["conditional_test"] = "ParCorr"
    archive_path = native_path(job["case_id"], method, output_dir)
    core.save_snapshot({**result, "case_identity": job, "method": method}, archive_path)
    if method == "CCM":
        for key, suffix in [("lag_profile", "profile"), ("lag_peaks", "peaks")]:
            if isinstance(result["raw"].get(key), pd.DataFrame):
                core.write_table(result["raw"][key], output_dir / "lag" / f"{job['case_id']}_{suffix}.csv")
    edges = report_edges(result, method, job, archive_path)
    edges["seconds"] = time.perf_counter() - started
    return edges


def commit_case(job, method_tables, output_dir=DATA_DIR):
    result = pd.concat([method_tables[method] for method in METHODS], ignore_index=True)
    core.write_table(result, Path(output_dir) / "members" / f"{job['case_id']}_edges.csv")
    return read_completed_case(job, output_dir)


def run_case(job, output_dir, null_file):
    tables = {method: run_method_job(job, method, output_dir, null_file) for method in METHODS}
    return commit_case(job, tables, output_dir)


def run_parcorr_case(job, output_dir):
    if job["scenario"] != "baseline":
        raise ValueError("ParCorr is only authorized for the three original inputs.")
    result = run_method_job(job, PARCORR_METHOD, output_dir)
    core.write_table(result, Path(output_dir) / "parcorr" / f"{job['case_id']}_edges.csv")
    return read_completed_parcorr(job, output_dir)


def read_committed(job, methods, subfolder, output_dir):
    output_dir = Path(output_dir)
    path = output_dir / subfolder / f"{job['case_id']}_edges.csv"
    if not path.exists():
        return None
    table = pd.read_csv(path, float_precision="round_trip", dtype={
        "edge_detected": "boolean", "unoriented": "boolean", "raw_edge_detected": "boolean"})
    required = {"method", "source", "target", "status", "edge_detected", "unoriented",
                "score", "p_value", "native_sha256", *job}
    if not required.issubset(table) or table.status.isna().any() or table.status.eq("pending").any():
        raise ValueError("Committed case is missing required fields or retains pending directions.")
    expected = {(method, source, target) for method in methods for source, target in [("x", "y"), ("y", "x")]}
    if len(table) != len(expected) or table.duplicated(["method", "source", "target"]).any():
        raise ValueError("Committed case has incomplete or duplicate directional rows.")
    if set(table[["method", "source", "target"]].itertuples(index=False, name=None)) != expected:
        raise ValueError("Committed case methods/directions changed.")
    for key, value in job.items():
        if key not in table or not table[key].eq(value).all():
            raise ValueError(f"Committed identity changed: {key}")
    if core.file_digest(output_dir / "inputs" / f"{job['case_id']}.csv") != job["input_sha256"]:
        raise ValueError("Committed input bytes changed.")
    for method in methods:
        archive = validate_native(job, method, output_dir)
        expected_rows = report_edges(archive, method, job, native_path(job["case_id"], method, output_dir))
        actual = table[table.method.eq(method)].set_index(["source", "target"])
        for _, row in expected_rows.iterrows():
            saved = actual.loc[(row.source, row.target)]
            for column, value in row.items():
                if column in {"source", "target"}:
                    continue
                actual_value = saved[column]
                missing = pd.isna(value) and pd.isna(actual_value)
                empty_error = column == "error" and value == "" and pd.isna(actual_value)
                if not missing and not empty_error and (pd.isna(value) or pd.isna(actual_value) or value != actual_value):
                    raise ValueError(f"Native versus commit mismatch: {method}/{column}")
    return table


def read_completed_case(job, output_dir=DATA_DIR):
    return read_committed(job, METHODS, "members", output_dir)


def read_completed_parcorr(job, output_dir=DATA_DIR):
    return read_committed(job, [PARCORR_METHOD], "parcorr", output_dir)


def collect_results(jobs, completed, methods=None):
    methods = METHODS if methods is None else methods
    frames = []
    for job in jobs:
        frame = completed.get(job["case_id"])
        if frame is None:
            frame = pd.concat([core.empty_edges(method, "pending") for method in methods], ignore_index=True)
            for key, value in job.items():
                frame[key] = value
        frames.append(frame)
    result = pd.concat(frames, ignore_index=True)
    for column in ["edge_detected", "unoriented", "raw_edge_detected"]:
        if column in result:
            result[column] = result[column].astype("boolean")
    return result


def summarize_results(results, parcorr_results, data_dir=DATA_DIR):
    """Keep directional retention separate from contemporaneous uncertainty."""
    combined = pd.concat([results, parcorr_results], ignore_index=True)
    state_rows, pcmci_rows, gc_rows = [], [], []
    for (case_id, method), pair in combined.groupby(["case_id", "method"], sort=False):
        row = pair.iloc[0]
        meta = {name: row[name] for name in ["case_id", "scenario", "process_seed", "noise_sd", "method"]}
        valid = len(pair) == 2 and pair.status.eq("ok").all() and pair.edge_detected.notna().all()
        forward = bool(pair.loc[pair.source.eq("x"), "edge_detected"].fillna(False).iloc[0])
        reverse = bool(pair.loc[pair.source.eq("y"), "edge_detected"].fillna(False).iloc[0])
        unresolved = bool(pair.unoriented.fillna(False).any())
        state_rows.append({**meta, "valid": valid, "pending": bool(pair.status.eq("pending").any()),
            "x_to_y": forward if valid else None, "y_to_x": reverse if valid else None,
            "unresolved": unresolved, "both_directions": bool(valid and forward and reverse),
            "state": core.evidence_state(forward, reverse, valid, unresolved)})
        path = native_path(case_id, method, data_dir)
        if not path.exists() or not valid:
            continue
        archive = load_native(case_id, method, data_dir)
        raw, diagnostics = archive["raw"], archive["diagnostics"]
        if method in ["PCMCI+", PARCORR_METHOD]:
            for source, target in [(0, 1), (1, 0)]:
                for lag in range(MAX_LAG + 1):
                    symbol = str(raw["graph_holm"][source, target, lag])
                    pcmci_rows.append({**meta, "source": ["x", "y"][source], "target": ["x", "y"][target],
                        "lag_steps": lag, "lag_years": lag * DT * 1000,
                        "score": raw["val_matrix"][source, target, lag],
                        "p_value": raw["p_matrix"][source, target, lag],
                        "p_holm": raw["p_holm"][source, target, lag],
                        "graph": str(raw["graph"][source, target, lag]), "graph_holm": symbol,
                        "retained_arrow": symbol == "-->", "unresolved": symbol in ["o-o", "x-x"]})
        if method == "GC":
            residuals = raw["var_residuals"]
            for column, record in enumerate(["x", "y"]):
                for lag in range(11):
                    gc_rows.append({**meta, "record": record, "lag_steps": lag,
                        "lag_years": lag * DT * 1000, "acf": acf_at(residuals[:, column], lag),
                        "contemporaneous_correlation": float(np.corrcoef(residuals.T)[0, 1]),
                        "whiteness_p": diagnostics.get("residual_whiteness_p"),
                        "whiteness_status": diagnostics.get("residual_whiteness_status"),
                        "var_is_stable": diagnostics.get("var_is_stable")})
    states = pd.DataFrame(state_rows)
    summary_rows = []
    for (scenario, method), group in states.groupby(["scenario", "method"], sort=False):
        requested = len(PROCESS_SEEDS)
        summary_rows.append(dict(scenario=scenario, method=method, n_requested=requested,
            n_present=len(group), n_valid=int(group.valid.sum()), n_pending=int(group.pending.sum()),
            n_x_to_y=int(group.x_to_y.fillna(False).sum()), n_y_to_x=int(group.y_to_x.fillna(False).sum()),
            n_both=int(group.both_directions.sum()), n_unresolved=int(group.unresolved.sum()),
            n_both_resolved=int(group.state.eq("both").sum())))
    reference = results[results.scenario.eq("baseline")][["process_seed", "method", "source", "target", "score", "status"]]
    reference = reference.rename(columns={"score": "baseline_score", "status": "baseline_status"})
    changes = results[~results.scenario.eq("baseline")].merge(reference,
        on=["process_seed", "method", "source", "target"], how="left", validate="many_to_one")
    usable = changes.status.eq("ok") & changes.baseline_status.eq("ok") & np.isfinite(changes.score) & np.isfinite(changes.baseline_score)
    changes["score_change_usable"] = usable
    changes["delta_score"] = np.where(usable, changes.score - changes.baseline_score, np.nan)
    return {"case_states": states, "scenario_summary": pd.DataFrame(summary_rows),
        "pcmci_lag_diagnostics": pd.DataFrame(pcmci_rows), "gc_residual_diagnostics": pd.DataFrame(gc_rows),
        "ccm_lag_summary": results[results.method.eq("CCM")].copy(), "paired_score_changes": changes}


def save_outputs(results, parcorr_results, output_dir=DATA_DIR, make_figures=True):
    output_dir = Path(output_dir)
    core.write_table(results, output_dir / "member_results.csv")
    core.write_table(parcorr_results, output_dir / "parcorr_results.csv")
    tables = summarize_results(results, parcorr_results, output_dir)
    for name, table in tables.items():
        core.write_table(table, output_dir / f"{name}.csv")
    if make_figures:
        plot_results(data_dir=output_dir, output_dir=FIGURE_DIR)
    return tables


def write_status(state, completed, parcorr_complete, started, output_dir=DATA_DIR, error=""):
    core.write_json({"state": state, "completed_pairs": completed, "requested_pairs": 12,
        "completed_parcorr": parcorr_complete, "requested_parcorr": 3,
        "configured_workers": N_WORKERS, "new_native_calls_requested": 39,
        "elapsed_seconds": time.perf_counter() - started,
        "updated_utc": datetime.now(timezone.utc).isoformat(), "error": error}, Path(output_dir) / "run_status.json")


def _plot_build_diagnostic_tables(data_dir):
    """Extract directional, lag and residual evidence without new model fits."""
    data_dir = Path(data_dir)
    jobs = pd.read_csv(data_dir / 'input_cases.csv', float_precision='round_trip')
    direction_rows, pcmci_rows, profile_tables, peak_rows = [], [], [], []
    gc_rows, acf_rows, variability_rows = [], [], []
    methods = ['PCMCI+', 'LKIF', 'GC', 'CCM']
    for job in jobs.to_dict('records'):
        meta = {key: job[key] for key in ['case_id', 'scenario', 'process_seed', 'noise_sd']}
        series = pd.read_csv(data_dir / 'inputs' / f"{job['case_id']}.csv")
        with np.load(data_dir / 'prepared' / f"history_{job['process_seed']}.npz") as history:
            sample = ((history['age'] <= 29.) & (history['age'] >= 1.)
                      & (history['time'] % 20 == 0))
            clean = history['clean_temperature'][:, sample]
            temperature = history[f"{job['scenario']}_temperature"][:, sample]
            fluctuation = temperature if job['scenario'] == 'random_response' else temperature - clean
            for record, box in [('x', 0), ('y', 3)]:
                values = series[record].to_numpy(dtype=float)
                forced_sd = float(np.std(clean[box]))
                fluctuation_sd = float(np.std(fluctuation[box]))
                variability_rows.append(dict(meta, record=record, input_sd=float(np.std(values)),
                    forced_sd=forced_sd, fluctuation_sd=fluctuation_sd,
                    fluctuation_to_forced_sd=fluctuation_sd / forced_sd,
                    input_acf_20yr=float(np.corrcoef(values[:-1], values[1:])[0, 1]),
                    input_acf_1500yr=float(np.corrcoef(values[:-75], values[75:])[0, 1]),
                    forced_response_removed=job['scenario'] == 'random_response'))
        case_methods = methods + (['PCMCI+_ParCorr'] if job['scenario'] == 'baseline' else [])
        committed_methods = set()
        if (data_dir / 'members' / f"{job['case_id']}_edges.csv").exists():
            read_completed_case(job, data_dir)
            committed_methods.update(methods)
        if (data_dir / 'parcorr' / f"{job['case_id']}_edges.csv").exists():
            read_completed_parcorr(job, data_dir)
            committed_methods.add('PCMCI+_ParCorr')
        for method in case_methods:
            native = (load_native(job['case_id'], method, data_dir) if method in committed_methods
                      else {'edges': pd.DataFrame(), 'raw': {}, 'diagnostics': {}})
            edges = native.get('edges', pd.DataFrame())
            raw, diagnostics = native.get('raw', {}), native.get('diagnostics', {})
            for source, target in [('x', 'y'), ('y', 'x')]:
                selected = edges[edges.source.eq(source) & edges.target.eq(target)] if len(edges) else pd.DataFrame()
                row = selected.iloc[0].to_dict() if len(selected) else {}
                available = row.get('status') == 'ok' and pd.notna(row.get('edge_detected'))
                contemporaneous = method in ['PCMCI+', 'PCMCI+_ParCorr']
                unresolved = (bool(diagnostics.get('undirected_contemporaneous', False))
                              or bool(diagnostics.get('conflicting_contemporaneous', False))) if contemporaneous else (
                              bool(row.get('unoriented', False)) if available and method == 'CCM' else False)
                direction_rows.append(dict(meta, method=method, source=source, target=target,
                    status=row.get('status', 'unavailable'), available=available,
                    retained=bool(row['edge_detected']) if available else False,
                    unresolved=unresolved, unresolved_assessed=method in ['PCMCI+', 'PCMCI+_ParCorr', 'CCM'],
                    native_score=row.get('score', np.nan), p_value=row.get('p_value', np.nan),
                    raw_significant=bool(row.get('raw_edge_detected', False)) if available and method == 'CCM' else False,
                    lag_status=row.get('lag_status', 'not_applicable')))
                if contemporaneous:
                    source_index, target_index = (0, 1) if source == 'x' else (1, 0)
                    for lag in range(6):
                        present = all(key in raw for key in ['graph_holm', 'val_matrix', 'p_matrix', 'p_holm'])
                        pcmci_rows.append(dict(meta, method=method, source=source, target=target,
                            lag_steps=lag, lag_years=20 * lag, available=present,
                            score=float(raw['val_matrix'][source_index, target_index, lag]) if present else np.nan,
                            p_raw=float(raw['p_matrix'][source_index, target_index, lag]) if present else np.nan,
                            p_holm=float(raw['p_holm'][source_index, target_index, lag]) if present else np.nan,
                            graph=str(raw['graph_holm'][source_index, target_index, lag]) if present else 'failed'))
                if method == 'CCM':
                    peaks = raw.get('lag_peaks', pd.DataFrame())
                    peaks = peaks[peaks.source.eq(source) & peaks.target.eq(target)] if len(peaks) else peaks
                    requested = 100
                    for kind in ['past', 'future', 'boundary', 'zero', 'tied', 'failed']:
                        count = int(peaks.peak_kind.eq(kind).sum()) if len(peaks) else 0
                        if kind == 'failed':
                            count = count + max(0, requested - len(peaks))
                        peak_rows.append(dict(meta, source=source, target=target, peak_kind=kind,
                                              count=count, requested_libraries=requested, fraction=count / requested))
            if method == 'CCM' and isinstance(raw.get('lag_profile'), pd.DataFrame):
                profile = raw['lag_profile'].copy()
                for key, value in meta.items():
                    profile[key] = value
                profile['lag_years'] = profile.lag_steps * 20
                profile_tables.append(profile)
            if method == 'GC':
                residuals = np.asarray(raw.get('var_residuals', []), dtype=float)
                usable = residuals.ndim == 2 and residuals.shape[1] == 2 and np.isfinite(residuals).all()
                gc_rows.append(dict(meta, available=usable,
                    contemporaneous_residual_correlation=float(np.corrcoef(residuals.T)[0, 1]) if usable else np.nan,
                    whiteness_p=diagnostics.get('residual_whiteness_p', np.nan),
                    whiteness_status=diagnostics.get('residual_whiteness_status', 'unavailable'),
                    var_is_stable=diagnostics.get('var_is_stable', False)))
                for column, record in enumerate(['x', 'y']):
                    for lag in range(26):
                        value = np.nan
                        if usable:
                            value = 1. if lag == 0 else float(np.corrcoef(residuals[:-lag, column], residuals[lag:, column])[0, 1])
                        acf_rows.append(dict(meta, record=record, lag_steps=lag,
                                             lag_years=20 * lag, acf=value, available=usable))
    profiles = pd.concat(profile_tables, ignore_index=True) if profile_tables else pd.DataFrame(
        columns=['case_id', 'scenario', 'process_seed', 'source', 'target', 'lag_steps', 'lag_years', 'rho_mean'])
    return dict(direction_evidence=pd.DataFrame(direction_rows), pcmci_lag_evidence=pd.DataFrame(pcmci_rows),
                ccm_lag_profiles=profiles, ccm_peak_fractions=pd.DataFrame(peak_rows),
                gc_residual_checks=pd.DataFrame(gc_rows), gc_residual_acf=pd.DataFrame(acf_rows),
                input_variability=pd.DataFrame(variability_rows))


def _plot_excitation_specs():
    scenarios = ['baseline', 'noise_x10', 'noise_x100', 'random_response']
    names = ['A: noise SD 0.005', 'B: noise SD 0.05', 'C: noise SD 0.5', 'D: random response']
    short = ['A\n0.005', 'B\n0.05', 'C\n0.5', 'D\nresponse']
    colors = ['#333333', '#4477AA', '#AA3377', '#228833']
    return scenarios, names, short, colors


def _plot_excitation_order(table):
    scenarios, _, _, _ = _plot_excitation_specs()
    result = table[['case_id', 'scenario', 'process_seed']].drop_duplicates().copy()
    result['scenario_order'] = result.scenario.map({key: i for i, key in enumerate(scenarios)})
    return result.sort_values(['scenario_order', 'process_seed']).reset_index(drop=True)


def _plot_excitation_style(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=10)


def _plot_excitation_inputs(tables, data_dir, output_dir):
    scenarios, names, short, colors = _plot_excitation_specs()
    variability = tables['input_variability']
    seed = int(variability.process_seed.min())
    with np.load(Path(data_dir) / 'prepared' / f'history_{seed}.npz') as saved:
        history = {key: saved[key] for key in saved.files}
    age = history['age']
    zoom = (age <= 15.) & (age >= 12.)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10.5))
    for column, (box, record) in enumerate([(0, 'South / T1'), (3, 'North / T4')]):
        for scenario, name, color in zip(scenarios[:3], names[:3], colors[:3]):
            axes[0, column].plot(age[zoom], history[f'{scenario}_temperature'][box, zoom],
                                color=color, lw=1.1, alpha=0.9, label=name)
        axes[0, column].plot(age[zoom], history['clean_temperature'][box, zoom],
                            color='0.55', ls='--', lw=1.5, label='Known noise-free response')
        axes[0, column].set(xlim=(15, 12), xlabel='Age (kyr BP)', ylabel='Temperature (model units)', title=record)
    axes[0, 0].legend(fontsize=8, loc='lower left')
    for box, label, color in [(0, 'South / T1', '#0072B2'), (3, 'North / T4', '#D55E00')]:
        axes[1, 0].plot(age[zoom], history['random_response_temperature'][box, zoom], color=color, lw=1.1, label=label)
    axes[1, 0].axhline(0, color='0.5', lw=0.7)
    axes[1, 0].set(xlim=(15, 12), xlabel='Age (kyr BP)', ylabel='Temperature difference (model units)',
                   title='D: baseline minus the known noise-free response')
    axes[1, 0].legend(fontsize=9)
    for record, offset, color, name in [('x', -0.12, '#0072B2', 'South'), ('y', 0.12, '#D55E00', 'North')]:
        for position, scenario in enumerate(scenarios[:3]):
            values = variability[variability.record.eq(record) & variability.scenario.eq(scenario)].sort_values('process_seed')
            x = position + offset + np.linspace(-0.035, 0.035, len(values))
            axes[1, 1].scatter(x, values.fluctuation_to_forced_sd, color=color, s=42,
                               label=name if position == 0 else None)
    axes[1, 1].axhline(1, color='0.5', ls='--', lw=1)
    axes[1, 1].set(yscale='log', xticks=range(3), xticklabels=short[:3],
                   ylabel='SD(random response) / SD(forced response)',
                       title='Random-response amplitude relative to the forced cycle')
    axes[1, 1].legend(fontsize=9)
    for ax in axes.flat:
        _plot_excitation_style(ax)
    fig.suptitle('TBK: fixed dynamics and accurate ages, different stochastic excitation', fontsize=15, y=0.98)
    fig.text(0.08, 0.025, 'Curves: preselected history R1. Points: all three histories, over the 29–1 kyr analysis window.\n'
             'D removes an exactly known model response; it is a diagnostic, not an observational preprocessing prescription.', fontsize=10)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.92, bottom=0.12, hspace=0.40, wspace=0.28)
    save_figure(fig, output_dir, '00_excitation_and_inputs', axes=axes)
    plt.close(fig)


def _plot_excitation_directions(tables, output_dir):
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    evidence = tables['direction_evidence']
    order = _plot_excitation_order(evidence)
    scenarios, names, _, _ = _plot_excitation_specs()
    seeds = sorted(order.process_seed.unique())
    seed_labels = {seed: f'R{i + 1}' for i, seed in enumerate(seeds)}
    labels = [f"{chr(65 + scenarios.index(row.scenario))}  {seed_labels[row.process_seed]}" for row in order.itertuples()]
    colors = ['#F0F0F0', '#66A182', '#E5CF77', '#555555', '#CCCCCC']
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, 5.5), cmap.N)
    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    for ax, method in zip(axes.flat, ['PCMCI+', 'LKIF', 'GC', 'CCM']):
        matrix = np.zeros((len(order), 3))
        text = np.full((len(order), 3), 'No', dtype=object)
        for i, case in order.iterrows():
            subset = evidence[evidence.case_id.eq(case.case_id) & evidence.method.eq(method)]
            for j, source in enumerate(['x', 'y']):
                row = subset[subset.source.eq(source)]
                available = len(row) == 1 and bool(row.iloc[0].available)
                matrix[i, j] = 1 if available and row.iloc[0].retained else 0 if available else 3
                text[i, j] = 'Yes' if matrix[i, j] == 1 else 'No' if available else 'NA'
            if method in ['GC', 'LKIF']:
                matrix[i, 2], text[i, 2] = 4, 'N/A'
            elif len(subset) != 2 or not subset.available.all():
                matrix[i, 2], text[i, 2] = 3, 'NA'
            elif subset.unresolved.any():
                matrix[i, 2], text[i, 2] = 2, 'Yes'
        ax.imshow(matrix, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')
        for i in range(len(order)):
            for j in range(3):
                ax.text(j, i, text[i, j], ha='center', va='center', fontsize=10,
                        color='white' if matrix[i, j] == 3 else 'black')
        for boundary in [2.5, 5.5, 8.5]:
            ax.axhline(boundary, color='white', lw=3)
        ax.set(xticks=range(3), xticklabels=['South → North\nretained', 'North → South\nretained', 'Unresolved\nevidence'],
               yticks=range(len(order)), yticklabels=labels, title=method + (' (lag screened)' if method == 'CCM' else ''))
        ax.tick_params(length=0, labelsize=10)
    fig.suptitle('Directional retention and unresolved evidence are separate outcomes', fontsize=15, y=0.98)
    fig.legend(handles=[Patch(facecolor=colors[1], label='Direction retained'), Patch(facecolor=colors[2], label='Unresolved present'),
                        Patch(facecolor=colors[4], label='Not assessed by this wrapper'), Patch(facecolor=colors[3], label='Failed / unavailable')],
               ncol=2, loc='upper center', bbox_to_anchor=(0.52, 0.947), frameon=False, fontsize=10)
    fig.text(0.08, 0.02, 'A: SD 0.005; B: SD 0.05; C: SD 0.5; D: baseline minus known noise-free response. R1–R3 are paired histories.\n'
             'PCMCI: unresolved contemporaneous link. CCM: raw-significant direction with unresolved lag evidence.\n'
             'Unresolved evidence does not hide retained arrows. These matrices are not a ranking of physical-network accuracy.', fontsize=9.5)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.85, bottom=0.15, hspace=0.32, wspace=0.28)
    save_figure(fig, output_dir, '01_directional_evidence', axes=axes)
    plt.close(fig)


def _plot_excitation_pcmci(tables, output_dir):
    from matplotlib.colors import Normalize
    lag_table = tables['pcmci_lag_evidence']
    order = _plot_excitation_order(lag_table)
    scenarios, _, _, _ = _plot_excitation_specs()
    seeds = sorted(order.process_seed.unique())
    seed_labels = {seed: f'R{i + 1}' for i, seed in enumerate(seeds)}
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), gridspec_kw={'height_ratios': [3, 1.1]})
    for row_index, method in enumerate(['PCMCI+', 'PCMCI+_ParCorr']):
        selected_order = order if row_index == 0 else order[order.scenario.eq('baseline')].reset_index(drop=True)
        for column, (source, direction) in enumerate([('x', 'South → North'), ('y', 'North → South')]):
            ax = axes[row_index, column]
            matrix = np.full((len(selected_order), 6), np.nan)
            symbols = np.full(matrix.shape, '', dtype=object)
            for i, case in selected_order.iterrows():
                selected = lag_table[lag_table.case_id.eq(case.case_id) & lag_table.method.eq(method) & lag_table.source.eq(source)]
                for item in selected.itertuples():
                    matrix[i, item.lag_steps] = item.score
                    symbols[i, item.lag_steps] = {'-->': '→', '<--': '←', 'o-o': 'o–o', 'x-x': 'x–x', 'failed': 'NA'}.get(item.graph, '·')
            image = ax.imshow(matrix, aspect='auto', cmap='viridis' if row_index == 0 else 'RdBu_r',
                              norm=Normalize(0, 1) if row_index == 0 else Normalize(-1, 1), interpolation='nearest')
            for i in range(len(selected_order)):
                for j in range(6):
                    value = matrix[i, j]
                    color = 'white' if np.isfinite(value) and ((row_index == 0 and value < .55) or (row_index == 1 and abs(value) > .5)) else 'black'
                    ax.text(j, i, symbols[i, j], ha='center', va='center', fontsize=10, color=color)
            labels = [f"{chr(65 + scenarios.index(case.scenario))}  {seed_labels[case.process_seed]}" for case in selected_order.itertuples()]
            ax.set(xticks=range(6), xticklabels=[0, 20, 40, 60, 80, 100], yticks=range(len(labels)), yticklabels=labels,
                   xlabel='Lag (yr)', title=('GPDC | ' if row_index == 0 else 'ParCorr baseline | ') + direction)
            if row_index == 0:
                for boundary in [2.5, 5.5, 8.5]:
                    ax.axhline(boundary, color='white', lw=3)
            ax.tick_params(length=0, labelsize=10)
            bar = fig.colorbar(image, ax=ax, fraction=0.04, pad=0.025)
            bar.set_label('GPDC adjacency statistic' if row_index == 0 else 'Partial correlation', fontsize=9)
    fig.suptitle('PCMCI+: lag-specific evidence, including the baseline ParCorr diagnostic', fontsize=15, y=0.98)
    fig.text(0.08, 0.023, 'Cell colour: native statistic, including nonsignificant entries. Symbols: Holm-filtered native graph.\n'
             '→ indicated direction; ← opposite direction at lag 0; o–o unoriented; x–x conflict; · no retained link; NA unavailable.\n'
             'A/B/C/D are the four input conditions. GPDC and ParCorr use different score scales; ParCorr is tested on A only.', fontsize=10)
    fig.subplots_adjust(left=0.08, right=0.96, top=0.91, bottom=0.16, hspace=0.42, wspace=0.28)
    save_figure(fig, output_dir, '02_pcmci_lag_evidence', axes=axes)
    plt.close(fig)


def _plot_excitation_ccm(tables, output_dir):
    from matplotlib.patches import Patch
    scenarios, names, _, colors = _plot_excitation_specs()
    profiles, peaks = tables['ccm_lag_profiles'], tables['ccm_peak_fractions']
    evidence = tables['direction_evidence']
    order = _plot_excitation_order(peaks)
    seeds = sorted(order.process_seed.unique())
    seed_labels = {seed: f'R{i + 1}' for i, seed in enumerate(seeds)}
    kinds = ['past', 'future', 'boundary', 'zero', 'tied', 'failed']
    peak_colors = ['#4477AA', '#EE7733', '#AA3377', '#CCCCCC', '#DDCC77', '#333333']
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    for column, (source, direction) in enumerate([('x', 'South → North'), ('y', 'North → South')]):
        for scenario, name, color in zip(scenarios, names, colors):
            selected = profiles[profiles.scenario.eq(scenario) & profiles.source.eq(source)]
            for case_id, case in selected.groupby('case_id'):
                axes[0, column].plot(case.lag_years, case.rho_mean, color=color, alpha=.3, lw=.8)
            votes = evidence[evidence.method.eq('CCM') & evidence.source.eq(source) & evidence.scenario.eq(scenario)]
            label = f"{name}: raw {int(votes.raw_significant.sum())}/3, kept {int(votes.retained.sum())}/3"
            if len(selected):
                summary = selected.groupby('lag_years').rho_mean.median()
                axes[0, column].plot(summary.index, summary.values, color=color, lw=2, label=label)
        axes[0, column].axvline(0, color='.5', ls=':', lw=1)
        axes[0, column].set(xlim=(-600, 600), xlabel='Target shift (yr)', ylabel='Mean-library reconstruction skill ρ', title=direction)
        axes[0, column].legend(fontsize=8, loc='upper center', bbox_to_anchor=(.5, -.18),
                               ncol=2, columnspacing=1, handlelength=1.7)
        bottom = np.zeros(len(order))
        for kind, color in zip(kinds, peak_colors):
            values = []
            for case in order.itertuples():
                subset = peaks[peaks.case_id.eq(case.case_id) & peaks.source.eq(source) & peaks.peak_kind.eq(kind)]
                values.append(float(subset.fraction.iloc[0]) if len(subset) else (1 if kind == 'failed' else 0))
            heights = 100 * np.asarray(values)
            axes[1, column].bar(range(len(order)), heights, bottom=bottom, color=color, edgecolor='white', lw=.5)
            bottom = bottom + heights
        labels = [f"{chr(65 + scenarios.index(case.scenario))}\n{seed_labels[case.process_seed]}" for case in order.itertuples()]
        axes[1, column].set(ylim=(0, 100), xticks=range(len(order)), xticklabels=labels,
                            ylabel='Library peak fraction (%)', title='All 100 library peaks for each individual input')
        for boundary in [2.5, 5.5, 8.5]:
            axes[1, column].axvline(boundary, color='white', lw=4)
    for ax in axes.flat:
        _plot_excitation_style(ax)
    fig.suptitle('CCM: reconstruction skill and timing consistency', fontsize=15, y=.98)
    fig.legend(handles=[Patch(facecolor=c, label=k.capitalize()) for k, c in zip(kinds, peak_colors)],
               loc='lower center', bbox_to_anchor=(.52, .085), ncol=6, frameon=False)
    fig.text(.08, .025, 'Thin curves: three histories; thick curves: their descriptive median. Each case curve averages 100 paired L=300 libraries.\n'
             'The 90% sign rule applies to each input separately. Peak bars include raw-nonsignificant cases; no causal probability is implied.', fontsize=10)
    fig.subplots_adjust(left=.08, right=.98, top=.91, bottom=.17, hspace=.65, wspace=.25)
    save_figure(fig, output_dir, '03_ccm_reconstruction_and_timing', axes=axes)
    plt.close(fig)


def _plot_excitation_gc(tables, output_dir):
    scenarios, names, short, colors = _plot_excitation_specs()
    acf = tables['gc_residual_acf']
    diagnostics = tables['gc_residual_checks']
    fig, axes = plt.subplots(2, 2, figsize=(14, 10.5))
    for column, (record, title) in enumerate([('x', 'South residuals'), ('y', 'North residuals')]):
        for scenario, name, color in zip(scenarios, names, colors):
            selected = acf[acf.scenario.eq(scenario) & acf.record.eq(record) & acf.lag_steps.gt(0)]
            for case_id, case in selected.groupby('case_id'):
                axes[0, column].plot(case.lag_years, case.acf, color=color, alpha=.3, lw=.8)
            summary = selected.groupby('lag_years').acf.median()
            axes[0, column].plot(summary.index, summary.values, color=color, lw=2, label=name)
        axes[0, column].axhline(0, color='.5', lw=.8)
        axes[0, column].set(xlabel='Lag (yr)', ylabel='Residual autocorrelation', title=title)
    axes[0, 0].legend(fontsize=8)
    for position, (scenario, color) in enumerate(zip(scenarios, colors)):
        group = diagnostics[diagnostics.scenario.eq(scenario)].sort_values('process_seed')
        x = position + np.linspace(-.08, .08, len(group))
        axes[1, 0].scatter(x, group.contemporaneous_residual_correlation, color=color, s=50)
        p = pd.to_numeric(group.whiteness_p, errors='coerce').to_numpy()
        values = -np.log10(np.maximum(p, 1e-300))
        axes[1, 1].scatter(x, values, color=color, s=50)
    axes[1, 0].axhline(0, color='.5', lw=.8)
    axes[1, 0].set(ylim=(-1.03, 1.03), ylabel='South–North residual correlation', title='Unexplained contemporaneous association')
    axes[1, 1].axhline(-np.log10(.05), color='.5', ls='--', lw=1, label='p = 0.05')
    axes[1, 1].set(ylabel='−log10(whiteness p)', title='Native VAR(1) residual whiteness test')
    axes[1, 1].legend(fontsize=9)
    for ax in axes[1]:
        ax.set_xticks(range(4), short)
        ax.set_xlim(-.4, 3.4)
    for ax in axes.flat:
        _plot_excitation_style(ax)
    fig.suptitle('GC: retained arrows must be read alongside residual diagnostics', fontsize=15, y=.98)
    fig.text(.08, .025, 'All diagnostics use saved native VAR(1) residuals; no new regression is fitted. Thin curves and dots: each of three histories.\n'
             'Thick curves: descriptive medians. Higher −log10(p) means stronger evidence against white residuals; p=0 is displayed at 300.', fontsize=10)
    fig.subplots_adjust(left=.09, right=.98, top=.91, bottom=.13, hspace=.40, wspace=.28)
    save_figure(fig, output_dir, '04_gc_residual_diagnostics', axes=axes)
    plt.close(fig)


def plot_results(data_dir=DATA_DIR, output_dir=FIGURE_DIR):
    """Build diagnostic tables and figures solely from saved model/native outputs."""
    data_dir, output_dir = Path(data_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tables = _plot_build_diagnostic_tables(data_dir)
    for name, table in tables.items():
        core.write_table(table, data_dir / f'{name}.csv')
    _plot_excitation_inputs(tables, data_dir, output_dir)
    _plot_excitation_directions(tables, output_dir)
    _plot_excitation_pcmci(tables, output_dir)
    _plot_excitation_ccm(tables, output_dir)
    _plot_excitation_gc(tables, output_dir)
    return tables


def main():
    parameters = execution_parameters()
    check_parameters(parameters)
    core.write_json(parameters, DATA_DIR / "parameters.json")
    started = time.perf_counter()
    completed, parcorr_complete, jobs = {}, {}, []
    try:
        jobs, diagnostics = prepare_inputs(DATA_DIR)
        completed = reuse_baselines(jobs, DATA_DIR)
        baseline_jobs = [job for job in jobs if job["scenario"] == "baseline"]
        for job in jobs:
            result = read_completed_case(job, DATA_DIR)
            if result is not None:
                completed[job["case_id"]] = result
        for job in baseline_jobs:
            result = read_completed_parcorr(job, DATA_DIR)
            if result is not None:
                parcorr_complete[job["case_id"]] = result
        if PREPARE_ONLY:
            write_status("prepared", len(completed), len(parcorr_complete), started)
            print(f"Prepared 12 inputs; {len(completed)} complete primary cases, {len(parcorr_complete)} ParCorr.", flush=True)
            return
        pending = [job for job in jobs if job["case_id"] not in completed]
        pending_parcorr = [job for job in baseline_jobs if job["case_id"] not in parcorr_complete]
        tasks = [(job, method) for job in pending for method in METHODS]
        tasks.extend((job, PARCORR_METHOD) for job in pending_parcorr)
        if RUN_CAUSAL and tasks:
            with threadpool_limits(limits=1):
                for seed in sorted({job["estimator_seed"] for job in pending}):
                    causal.prepare_gpdc_null(core.gpdc_null_path(seed), sample_sizes=[1401 - 2 * MAX_LAG],
                        sig_samples=GPDC_NULL_SAMPLES, seed=seed)
            workers = min(N_WORKERS, len(tasks))
            print(f"Excitation diagnostic: {len(tasks)} native tasks; {workers} workers, one numerical thread each.", flush=True)
            write_status("running", len(completed), len(parcorr_complete), started)
            case_tables = {job["case_id"]: {} for job in pending}
            context = multiprocessing.get_context("spawn")
            with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
                futures = {}
                for job, method in tasks:
                    future = pool.submit(run_parcorr_case, job, str(DATA_DIR)) if method == PARCORR_METHOD else pool.submit(
                        run_method_job, job, method, str(DATA_DIR), str(core.gpdc_null_path(job["estimator_seed"])))
                    futures[future] = (job, method)
                for future in as_completed(futures):
                    job, method = futures[future]
                    result = future.result()
                    if method == PARCORR_METHOD:
                        parcorr_complete[job["case_id"]] = result
                    else:
                        case_tables[job["case_id"]][method] = result
                        if len(case_tables[job["case_id"]]) == len(METHODS):
                            completed[job["case_id"]] = commit_case(job, case_tables[job["case_id"]], DATA_DIR)
                    print(f"Saved {job['case_id']}/{method}; primary {len(completed)}/12, ParCorr {len(parcorr_complete)}/3.", flush=True)
                    write_status("running", len(completed), len(parcorr_complete), started)
        results = collect_results(jobs, completed)
        parcorr = collect_results(baseline_jobs, parcorr_complete, [PARCORR_METHOD])
        save_outputs(results, parcorr)
        state = "complete" if len(completed) == 12 and len(parcorr_complete) == 3 else "partial"
        write_status(state, len(completed), len(parcorr_complete), started)
        print(f"Excitation diagnostic {state}: {len(completed)}/12 primary cases and {len(parcorr_complete)}/3 ParCorr.", flush=True)
    except BaseException as error:
        write_status("interrupted_or_failed", len(completed), len(parcorr_complete), started,
                     error=f"{type(error).__name__}: {error}")
        raise


if __name__ == "__main__":
    main()
