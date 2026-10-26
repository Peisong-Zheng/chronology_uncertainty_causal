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


WORKER_COUNTS = [16, 24, 32]  # Candidate counts for the 14900KF; this is a benchmark, not a claim that 32 is optimal.
MEMBERS_PER_SETTING = 2      # Two of each 16 pairs × 3 processes = 96 input pairs, enough to occupy 32 workers.
REPEATS = 1                  # Repeat a close comparison on the target PC before choosing its production setting.
SCORE_RTOL = 1e-10           # Tolerance for numerical results between worker counts on the same machine.
SCORE_ATOL = 1e-12           # Absolute tolerance near zero; directions and lag classifications must match exactly.

SOURCE_DIR = ROOT / "data" / "processed" / "02_check_causal_pipeline"
OUTPUT_ROOT = ROOT / "data" / "benchmarks" / "workers"
METHODS = ["PCMCI+", "LKIF", "CCM", "GC"]


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
            edges["case"] = case["label"]
            edges["method"] = method
            edges["seconds"] = time.perf_counter() - started
            rows.append(edges)
    return pd.concat(rows, ignore_index=True)


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
                  "timing": "Includes process startup, imports, input reads and method calls; excludes null precomputation and result-file writing.",
                  "memory": "Inspect Task Manager memory/commit and paging on Windows; this benchmark does not automatically approve a worker count."}
    (output / "parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    reference = None
    timings = []
    context = multiprocessing.get_context("spawn")
    print(f"Benchmarking {len(cases)} fixed input pairs. Outputs: {output}", flush=True)
    for repeat in range(REPEATS):
        for workers in WORKER_COUNTS:
            started = time.perf_counter()
            tables = []
            with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
                futures = [pool.submit(fit_case, case, settings) for case in cases]
                for future in as_completed(futures):
                    tables.append(future.result())
                    if len(tables) % 16 == 0 or len(tables) == len(cases):
                        print(f"Workers {workers}, repeat {repeat + 1}: {len(tables)}/{len(cases)} pairs", flush=True)
            seconds = time.perf_counter() - started
            result = pd.concat(tables, ignore_index=True)
            result.to_csv(output / f"workers{workers}_repeat{repeat + 1}_results.csv", index=False)
            if len(result) != len(cases) * 8 or not result["status"].eq("ok").all():
                raise RuntimeError("Some benchmark fits failed; inspect the saved results before comparing speeds.")
            ccm = result[result["method"] == "CCM"]
            if not ccm["lag_n_valid_libraries"].eq(settings["CCM_LAG_LIBRARIES"]).all():
                raise RuntimeError("Some CCM lag libraries failed; their shorter runtime is not a valid speed result.")
            if reference is None:
                reference = result
            else:
                check_equal_results(reference, result)
            timings.append({"workers": workers, "repeat": repeat + 1, "n_pairs": len(cases),
                            "wall_seconds": seconds, "pairs_per_hour": 3600 * len(cases) / seconds,
                            "results_match": True})
            pd.DataFrame(timings).to_csv(output / "throughput.csv", index=False)
            print(f"Workers {workers}: {seconds / 60:.2f} min, {3600 * len(cases) / seconds:.1f} pairs/hour; results match.", flush=True)
    print("Choose a count only after reviewing throughput and Windows memory pressure.", flush=True)


if __name__ == "__main__":
    main()
