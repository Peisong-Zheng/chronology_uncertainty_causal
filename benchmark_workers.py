"""Compare Windows process-pool sizes using identical prepared experiment 02 inputs.

First run 02 with PREPARE_ONLY=True to prepare its chronology-member inputs.
After benchmarking, set the chosen N_WORKERS and PREPARE_ONLY=False in 02.
This script does not
write into any experiment output folder or select scientific parameters.
Close competing numerical jobs before timing. Inspect Windows memory pressure
alongside the saved throughput; the fastest run is not useful if it causes paging.
"""

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from importlib.metadata import version
import hashlib
import json
import multiprocessing
import os
import platform
import time

ROOT = Path(__file__).resolve().parent
# Set these before scientific imports, including imports made by spawned workers.
for thread_name in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"]:
    os.environ[thread_name] = "1"
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "data" / "cache" / "matplotlib"))

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from toolbox.causal import prepare_gpdc_null, run_pcmciplus, run_lkif, run_ccm, run_granger


WORKER_COUNTS = [16, 24, 32]  # Candidate counts for this CPU; 32 was fastest in the local Windows benchmark.
MEMBERS_PER_SETTING = 2      # Two of each 16 pairs × 3 processes = 96 input pairs, enough to occupy 32 workers.
REPEATS = 1                  # Repeat a close comparison on the target PC before choosing its production setting.
SCORE_RTOL = 1e-10           # Tolerance for numerical results between worker counts on the same machine.
SCORE_ATOL = 1e-12           # Absolute tolerance near zero; directions and lag classifications must match exactly.

SOURCE_DIR = ROOT / "data" / "processed" / "02_check_causal_pipeline"
OUTPUT_ROOT = ROOT / "data" / "benchmarks" / "workers"
METHODS = ["PCMCI+", "LKIF", "CCM", "GC"]


def fingerprint_value(value, component, rows):
    """Record exact native values and their structure without copying full arrays."""
    if isinstance(value, pd.DataFrame):
        fingerprint_value(list(value.columns), component + "/columns", rows)
        fingerprint_value(value.index.to_numpy(), component + "/index", rows)
        for position, column in enumerate(value.columns):
            fingerprint_value(value.iloc[:, position].to_numpy(),
                              component + f"/column[{position}]={column!r}", rows)
        return
    if isinstance(value, dict):
        keys = sorted(value, key=repr)
        fingerprint_value([repr(key) for key in keys], component + "/keys", rows)
        for key in keys:
            fingerprint_value(value[key], component + f"/[{key!r}]", rows)
        return
    if isinstance(value, (list, tuple, set)):
        items = sorted(value, key=repr) if isinstance(value, set) else value
        # Retain sequence type and length, including empty native collections.
        fingerprint_value(f"{type(value).__name__}:{len(items)}", component + "/structure", rows)
        for position, item in enumerate(items):
            fingerprint_value(item, component + f"/[{position}]", rows)
        return
    if value is None:
        array = np.asarray("None")
        dtype = "NoneType"
    elif isinstance(value, (np.ndarray, np.generic, str, bool, int, float)):
        array = np.asarray(value)
        dtype = str(array.dtype)
    else:
        raise TypeError(f"Unsupported native value at {component}: {type(value).__name__}")
    if array.dtype.hasobject:
        # Object columns contain labels and optional decisions, not object addresses.
        def json_scalar(item):
            if item is pd.NA:
                return {"missing": "pd.NA"}
            if isinstance(item, np.generic):
                return item.item()
            raise TypeError(f"Unsupported object at {component}: {type(item).__name__}")
        payload = json.dumps(array.tolist(), default=json_scalar, sort_keys=True,
                             ensure_ascii=True, separators=(",", ":")).encode("utf-8")
    else:
        payload = np.ascontiguousarray(array).tobytes()
    rows.append({"component": component, "dtype": dtype,
                 "shape": json.dumps(list(array.shape)),
                 "sha256": hashlib.sha256(payload).hexdigest()})


def native_fingerprints(result, method, case_label):
    """Cover all raw results and edges, including statsmodels' fitted models."""
    raw = dict(result["raw"])
    if method == "GC":
        for direction in ["x_to_y", "y_to_x"]:
            native = dict(raw[direction])
            for name in ["restricted", "full"]:
                fit = native[name]
                # Keep the complete regression vectors and native uncertainty outputs.
                # Object identity and lazy-cache layout vary across processes.
                fitted_values = {attribute: getattr(fit, attribute) for attribute in [
                    "params", "resid", "fittedvalues", "normalized_cov_params",
                    "bse", "tvalues", "pvalues", "ssr", "scale", "nobs",
                    "df_model", "df_resid", "rsquared", "rsquared_adj",
                    "centered_tss", "uncentered_tss", "ess", "mse_model",
                    "mse_resid", "mse_total", "fvalue", "f_pvalue", "llf",
                    "aic", "bic", "condition_number",
                ]}
                fitted_values["cov_params"] = fit.cov_params()
                fitted_values["conf_int"] = fit.conf_int()
                fitted_values["endog"] = fit.model.endog
                fitted_values["exog"] = fit.model.exog
                fitted_values["design_rank"] = fit.model.rank
                native[name] = fitted_values
            raw[direction] = native
    rows = []
    fingerprint_value(raw, "raw", rows)
    fingerprint_value(result["edges"], "edges", rows)
    fingerprint_value(result["diagnostics"], "diagnostics", rows)
    for row in rows:
        row["case"] = case_label
        row["method"] = method
    return rows


def select_cases():
    """Select a fixed prefix per setting without looking at causal outcomes."""
    settings = json.loads((SOURCE_DIR / "age_comparison_parameters.json").read_text())
    for package, saved_version in settings["packages"].items():
        if version(package) != saved_version:
            raise ValueError(f"Prepared experiment package {package} differs from this environment.")
    inputs = pd.read_csv(SOURCE_DIR / "age_member_inputs.csv")
    keys = ["replicate", "process", "scenario"]
    selected = inputs.sort_values(keys + ["member"]).groupby(keys, sort=False).head(MEMBERS_PER_SETTING)
    counts = selected.groupby(keys).size()
    if counts.empty or not counts.eq(MEMBERS_PER_SETTING).all():
        raise ValueError("Every setting needs the requested number of prepared members.")
    if not selected["status"].eq("ready").all():
        raise ValueError("Benchmark inputs contain chronology/resampling failures; do not silently select replacements.")
    if len(selected) < max(WORKER_COUNTS):
        raise ValueError("There must be at least as many benchmark cases as workers.")
    cases = []
    for row in selected.to_dict("records"):
        path = SOURCE_DIR / f"{row['label']}_series.csv"
        if not path.is_file():
            raise FileNotFoundError(f"Missing prepared input: {path}. Prepare 02 on this computer first.")
        seed = int(row["seed"])
        n_rows = len(pd.read_csv(path, usecols=["time"]))
        effective_n = n_rows - 2 * settings["MAX_LAG"]
        cache = ROOT / "data" / "cache" / "gpdc" / f"null_T{effective_n}_n{settings['GPDC_NULL_SAMPLES']}_seed{seed}.npz"
        cases.append({"label": row["label"], "input_path": str(path), "seed": seed,
                      "effective_n": effective_n, "null_path": str(cache),
                      "process": row["process"], "scenario": row["scenario"], "member": int(row["member"])})
    # Prepare each native null once in the parent, before worker processes start.
    prepared = set()
    for case in cases:
        if case["null_path"] not in prepared:
            prepare_gpdc_null(case["null_path"], sample_sizes=[case["effective_n"]],
                              sig_samples=settings["GPDC_NULL_SAMPLES"], seed=case["seed"])
            prepared.add(case["null_path"])
    return cases, settings


def fit_case(case, settings):
    """Run the same four native methods without writing production files."""
    series = pd.read_csv(case["input_path"])
    rows = []
    fingerprints = []
    with threadpool_limits(limits=1):
        for method in METHODS:
            started = time.perf_counter()
            try:
                if method == "PCMCI+":
                    result = run_pcmciplus(
                        series, tau_max=settings["MAX_LAG"], pc_alpha=settings["ALPHA"],
                        alpha=settings["ALPHA"], sig_samples=settings["GPDC_NULL_SAMPLES"],
                        seed=case["seed"],
                        gp_params=None if settings["GPDC_OPTIMIZE_KERNEL"] else {"optimizer": None},
                        null_dist_filename=case["null_path"],
                    )
                elif method == "LKIF":
                    result = run_lkif(series, max_lag=settings["MAX_LAG"], alpha=settings["ALPHA"])
                elif method == "CCM":
                    embedded_rows = len(series) - (settings["CCM_E"] - 1) * abs(settings["CCM_TAU"])
                    result = run_ccm(
                        series, E=settings["CCM_E"], tau=settings["CCM_TAU"],
                        exclusion_radius=settings["CCM_EXCLUSION_RADIUS"],
                        lib_size=int(settings["CCM_LIBRARY_FRACTION"] * embedded_rows),
                        sample=settings["CCM_SAMPLES"], n_surrogates=settings["CCM_SURROGATES"],
                        seed=case["seed"], alpha=settings["ALPHA"], lag_steps=settings["CCM_LAG_STEPS"],
                        n_libraries=settings["CCM_LAG_LIBRARIES"], sign_threshold=settings["CCM_LAG_SIGN_THRESHOLD"],
                    )
                else:
                    result = run_granger(series, max_lag=settings["GC_ORDER"], alpha=settings["ALPHA"])
                edges = result["edges"].copy()
                diagnostics = result["diagnostics"]
                if "unoriented" not in edges:
                    edges["unoriented"] = False
                if diagnostics.get("undirected_contemporaneous", False) or diagnostics.get("conflicting_contemporaneous", False):
                    edges["unoriented"] = True
            except Exception as error:
                edges = pd.DataFrame([
                    {"source": source, "target": target, "score": np.nan, "p_value": np.nan,
                     "edge_detected": None, "unoriented": False, "status": "error", "error": str(error)}
                    for source, target in [("x", "y"), ("y", "x")]
                ])
            else:
                try:
                    fingerprints.extend(native_fingerprints(result, method, case["label"]))
                except Exception as error:
                    raise RuntimeError(
                        f"Native-output validation failed for {case['label']}, {method}: {error}. "
                        "The estimator succeeded; this is a benchmark validation error, not a failed scientific fit."
                    ) from error
            edges["case"] = case["label"]
            edges["method"] = method
            edges["seconds"] = time.perf_counter() - started
            rows.append(edges)
    return pd.concat(rows, ignore_index=True), fingerprints


def check_equal_results(reference, result):
    """Require the same scientific decisions and numerical summaries at every count."""
    keys = ["case", "method", "source", "target"]
    if reference.duplicated(keys).any() or result.duplicated(keys).any():
        raise ValueError("Duplicate benchmark results.")
    first = reference.sort_values(keys).reset_index(drop=True)
    second = result.sort_values(keys).reset_index(drop=True)
    decision_columns = keys + ["status", "edge_detected", "unoriented", "lag_status", "lag_reason"]
    pd.testing.assert_frame_equal(first[decision_columns], second[decision_columns], check_dtype=False)
    for column in ["score", "p_value", "optimal_lag_steps", "lag_past_fraction", "lag_future_fraction"]:
        np.testing.assert_allclose(first[column], second[column], rtol=SCORE_RTOL, atol=SCORE_ATOL, equal_nan=True)


def check_equal_native(reference, result, mismatch_path):
    """Reject any changed native component and save the exact mismatch locations."""
    keys = ["case", "method", "component"]
    if reference.duplicated(keys).any() or result.duplicated(keys).any():
        raise ValueError("Duplicate native benchmark fingerprint components.")
    compared = reference.merge(result, on=keys, how="outer", suffixes=("_reference", "_candidate"),
                               indicator=True, validate="one_to_one")
    different = compared["_merge"].ne("both")
    for column in ["dtype", "shape", "sha256"]:
        different = different | compared[column + "_reference"].ne(compared[column + "_candidate"])
    if different.any():
        mismatches = compared.loc[different].sort_values(keys)
        mismatches.to_csv(mismatch_path, index=False)
        first = mismatches.iloc[0]
        raise AssertionError(
            f"Native outputs differ for {first['case']}, {first['method']}, {first['component']}; "
            f"{len(mismatches)} changed components saved to {mismatch_path}. "
            "Do not select this worker count from its timing."
        )


def main():
    if MEMBERS_PER_SETTING < 1 or REPEATS < 1 or any(count < 1 for count in WORKER_COUNTS):
        raise ValueError("Worker counts, selected members and repeats must be positive.")
    cases, settings = select_cases()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    output = OUTPUT_ROOT / stamp
    output.mkdir(parents=True, exist_ok=False)
    pd.DataFrame(cases).to_csv(output / "selected_cases.csv", index=False)
    parameters = {"worker_counts": WORKER_COUNTS, "cases": len(cases), "repeats": REPEATS,
                  "platform": platform.platform(), "python": platform.python_version(),
                  "logical_cpus": os.cpu_count(), "source_settings": settings,
                  "timing": "Includes process startup, imports, input reads, method calls, native fingerprinting and result transfer; excludes null precomputation and result-file writing.",
                  "validation": "Exact SHA-256 equality of every raw array/value, edge column and diagnostic, with shape and dtype checks. GC fitted models are compared through native regression arrays and statistics. No tolerance is applied to these fingerprints.",
                  "memory": "Inspect Task Manager memory/commit and paging on Windows; this benchmark does not automatically approve a worker count."}
    (output / "parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    reference = None
    reference_native = None
    timings = []
    context = multiprocessing.get_context("spawn")
    print(f"Benchmarking {len(cases)} fixed input pairs. Outputs: {output}", flush=True)
    for repeat in range(REPEATS):
        for workers in WORKER_COUNTS:
            started_utc = datetime.now(timezone.utc).isoformat()
            print(f"Workers {workers}, repeat {repeat + 1}: started {started_utc}", flush=True)
            started = time.perf_counter()
            tables = []
            fingerprint_rows = []
            with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
                futures = [pool.submit(fit_case, case, settings) for case in cases]
                for future in as_completed(futures):
                    table, native = future.result()
                    tables.append(table)
                    fingerprint_rows.extend(native)
                    if len(tables) % 16 == 0 or len(tables) == len(cases):
                        print(f"Workers {workers}, repeat {repeat + 1}: {len(tables)}/{len(cases)} pairs", flush=True)
            seconds = time.perf_counter() - started
            finished_utc = datetime.now(timezone.utc).isoformat()
            result = pd.concat(tables, ignore_index=True)
            result.to_csv(output / f"workers{workers}_repeat{repeat + 1}_results.csv", index=False)
            native = pd.DataFrame(fingerprint_rows, columns=["case", "method", "component", "dtype", "shape", "sha256"])
            native = native.sort_values(["case", "method", "component"])
            native.to_csv(output / f"workers{workers}_repeat{repeat + 1}_native_fingerprints.csv", index=False)
            if len(result) != len(cases) * 8 or not result["status"].eq("ok").all():
                raise RuntimeError("Some benchmark fits failed; inspect the saved results before comparing speeds.")
            ccm = result[result["method"] == "CCM"]
            if not ccm["lag_n_valid_libraries"].eq(settings["CCM_LAG_LIBRARIES"]).all():
                raise RuntimeError("Some CCM lag libraries failed; their shorter runtime is not a valid speed result.")
            if reference is None:
                reference = result
                reference_native = native
            else:
                check_equal_native(reference_native, native,
                                   output / f"workers{workers}_repeat{repeat + 1}_native_mismatches.csv")
                check_equal_results(reference, result)
            timings.append({"workers": workers, "repeat": repeat + 1, "n_pairs": len(cases),
                            "started_utc": started_utc, "finished_utc": finished_utc,
                            "wall_seconds": seconds, "pairs_per_hour": 3600 * len(cases) / seconds,
                            "results_match": True, "native_outputs_match_exactly": True})
            pd.DataFrame(timings).to_csv(output / "throughput.csv", index=False)
            print(f"Workers {workers}: {seconds / 60:.2f} min, {3600 * len(cases) / seconds:.1f} pairs/hour; complete native outputs match exactly.", flush=True)
    print("Choose a count only after reviewing throughput and Windows memory pressure.", flush=True)


if __name__ == "__main__":
    main()
