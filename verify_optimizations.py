"""Check native GPDC caching and CCM object reuse on this computer.

Run ``python verify_optimizations.py`` before a large Windows experiment.
This uses the complete 100-library, 31-lag CCM screen on two 1401-point
records. Allow several minutes; it does not run any chronology ensembles.
Outputs go to data/benchmarks/optimization_checks/<UTC timestamp>/.

The reference and optimized calculations must agree exactly on the same
machine. This is not a test of bitwise agreement between operating systems.
"""

import os

# Set numerical thread limits before importing NumPy or any causal package.
for variable in [
    "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
]:
    os.environ[variable] = "1"

from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
import json
import platform
import sys
import time

import numpy as np
import pandas as pd
import pyEDM
from threadpoolctl import threadpool_limits

from toolbox.causal import (
    filter_ccm_by_lag, prepare_gpdc_null, run_ccm, run_pcmciplus,
    summarize_ccm_lag,
)
from toolbox.process import generate_process


# Fixed acceptance inputs, matching the main experiment's original baseline.
SEED = 20261009                 # Process seed; both cases share the same innovations.
DURATION = 30.0                 # Retained simulated record length, kyr.
DT = 0.020                      # Native step, kyr: 1501 points over 30 kyr.
TRUE_LAG = 0.020                # X drives Y one step later in the coupled case.
MEMORY_X = 0.5                  # Driver's intrinsic decay time, kyr.
MEMORY_Y = 0.2                  # Response's intrinsic decay time, kyr.
BACKGROUND_X = 1 / np.sqrt(1 - 0.95**2)  # Stationary X standard deviation.
BACKGROUND_Y = 1 / np.sqrt(1 - 0.8**2)   # Y standard deviation without coupling.
COUPLING_GAIN = 4.0             # Long-run response to a held unit driver.
BURN_TIMESCALES = 5             # Warm-up in units of the longest decay time.
ANALYSIS_OLD = 29.0             # Old analysis boundary, kyr BP.
ANALYSIS_YOUNG = 1.0            # Young analysis boundary: 1401 rows retained.

# Keep scientific settings unchanged while comparing execution strategies.
METHOD_SEED = SEED + 300        # Native GPDC null and CCM random-library seed.
MAX_LAG = 5                    # PCMCI+ maximum lag in 20-year steps.
ALPHA = 0.05                   # Per-direction significance threshold.
GPDC_NULL_SAMPLES = 10000       # Native GPDC reference-null draws.
CCM_E = 3                      # Embedding dimension.
CCM_TAU = -1                   # Historical embedding interval, in steps.
CCM_EXCLUSION_RADIUS = 10       # Neighbor exclusion radius, in steps.
CCM_LIBRARY_FRACTION = 0.80     # Fraction of embedded rows in each random library.
CCM_SAMPLES = 10               # Libraries averaged by the original Tp=0 test.
CCM_SURROGATES = 100            # Cut-and-swap surrogates for each direction.
CCM_LAG_STEPS = list(range(-15, 16))  # All 31 target shifts, from -300 to +300 yr.
CCM_LAG_LIBRARIES = 100         # Matched random libraries at every target shift.
CCM_LAG_SIGN_THRESHOLD = 0.90   # Existing temporal-consistency decision rule.

ROOT = Path(__file__).resolve().parent


def generate_inputs():
    """Generate coupled and independent records without external data files."""
    ar_x = np.exp(-DT / MEMORY_X)
    ar_y = np.exp(-DT / MEMORY_Y)
    noise = (
        BACKGROUND_X * np.sqrt(1 - ar_x**2),
        BACKGROUND_Y * np.sqrt(1 - ar_y**2),
    )
    burn_in = int(np.ceil(BURN_TIMESCALES * max(MEMORY_X, MEMORY_Y) / DT))
    first = int(round((DURATION - ANALYSIS_OLD) / DT))
    last = int(round((DURATION - ANALYSIS_YOUNG) / DT))
    age_grid = np.linspace(ANALYSIS_OLD, ANALYSIS_YOUNG, last - first + 1)
    cases = {}
    for name in ["linear", "independent"]:
        gain = COUPLING_GAIN if name == "linear" else 0.0
        process = generate_process(
            n_steps=int(round(DURATION / DT)) + 1, dt=DT,
            age_start=DURATION, lag=TRUE_LAG, ar_x=ar_x, ar_y=ar_y,
            coupling=(1 - ar_y) * gain, coupling_form="linear",
            noise_std=noise, burn_in=burn_in, rng=SEED,
        )
        series = process["series"].iloc[first:last + 1].copy().reset_index(drop=True)
        series["time"] = age_grid[0] - age_grid
        series["age"] = age_grid
        if len(series) != 1401:
            raise AssertionError("The acceptance inputs must contain 1401 rows.")
        cases[name] = series
    return cases


def reference_ccm_scan(series, original):
    """Use a fresh native CCM call for every library and lag, without reuse.

    This deliberately slow reference only verifies the optimized scan. Both
    paths use the production summary function, so decision rules are shared.
    The original Tp=0 scores, swap p values, and significance flags are reused.
    """
    settings = original["diagnostics"]
    lags = np.asarray(CCM_LAG_STEPS, dtype=int)
    history = (settings["E"] - 1) * abs(settings["tau"])
    first = max(history, -int(lags[0]))
    last = len(series) - 1 - int(lags[-1])
    anchors = np.arange(first, last + 1, dtype=int)
    valid_lib = np.zeros(len(series), dtype=bool)
    valid_lib[anchors] = True
    seeds = np.random.default_rng(METHOD_SEED).integers(
        1, 2**31 - 1, size=CCM_LAG_LIBRARIES,
    )
    rho = np.full((CCM_LAG_LIBRARIES, len(lags), 2), np.nan)
    data = series[["time", "x", "y"]]
    for library, seed in enumerate(seeds):
        for position, lag in enumerate(lags):
            native = pyEDM.CCM(
                dataFrame=data, columns="y", target="x",
                E=settings["E"], tau=settings["tau"], Tp=int(lag),
                exclusionRadius=settings["exclusion_radius"],
                libSizes=[settings["lib_size"]], sample=1, seed=int(seed),
                validLib=valid_lib, parallel=False, includeData=False, showPlot=False,
            )
            if int(native["LibSize"].iloc[0]) != settings["lib_size"]:
                raise AssertionError("Reference CCM changed the requested library size.")
            rho[library, position, :] = [native["y:x"].iloc[0], native["x:y"].iloc[0]]
    if not np.isfinite(rho).all():
        raise AssertionError("Reference CCM produced incomplete lag curves.")

    raw = dict(original["raw"])
    raw.update({
        "lag_rho": rho, "lag_steps": lags, "lag_library_seeds": seeds,
        "lag_directions": np.array(["x_to_y", "y_to_x"]),
        "lag_common_anchor_indices": anchors,
    })
    diagnostics = dict(settings)
    diagnostics["lag_filter"] = {"dt_kyr": settings["dt_kyr"]}
    reference = {"edges": original["edges"], "diagnostics": diagnostics, "raw": raw}
    return summarize_ccm_lag(reference, sign_threshold=CCM_LAG_SIGN_THRESHOLD)


def check_pcmci(series, cache_path, settings):
    """Compare complete native matrices, corrected matrices, and edge tables."""
    start = time.perf_counter()
    expected = run_pcmciplus(series, **settings)
    uncached_seconds = time.perf_counter() - start
    start = time.perf_counter()
    actual = run_pcmciplus(series, null_dist_filename=cache_path, **settings)
    cached_seconds = time.perf_counter() - start
    arrays = ["val_matrix", "p_matrix", "graph", "p_holm", "graph_holm"]
    for name in arrays:
        np.testing.assert_array_equal(actual["raw"][name], expected["raw"][name], err_msg=name)
    pd.testing.assert_frame_equal(actual["edges"], expected["edges"], check_exact=True)
    return {
        "method": "PCMCI+", "compared_arrays": arrays,
        "arrays_exactly_equal": True, "edge_tables_exactly_equal": True,
        "uncached_seconds": uncached_seconds, "cached_seconds": cached_seconds,
        "warm_speedup": uncached_seconds / cached_seconds,
    }


def check_ccm(series):
    """Compare all 100 x 31 x 2 lag scores and resulting directional decisions."""
    embedded_rows = len(series) - (CCM_E - 1) * abs(CCM_TAU)
    start = time.perf_counter()
    original = run_ccm(
        series, E=CCM_E, tau=CCM_TAU, exclusion_radius=CCM_EXCLUSION_RADIUS,
        lib_size=int(CCM_LIBRARY_FRACTION * embedded_rows), sample=CCM_SAMPLES,
        n_surrogates=CCM_SURROGATES, seed=METHOD_SEED, alpha=ALPHA, lag_steps=None,
    )
    original_seconds = time.perf_counter() - start
    start = time.perf_counter()
    expected = reference_ccm_scan(series, original)
    reference_seconds = time.perf_counter() - start
    start = time.perf_counter()
    actual = filter_ccm_by_lag(
        series, original, lag_steps=CCM_LAG_STEPS, n_libraries=CCM_LAG_LIBRARIES,
        sign_threshold=CCM_LAG_SIGN_THRESHOLD, seed=METHOD_SEED,
    )
    optimized_seconds = time.perf_counter() - start
    if actual["diagnostics"]["lag_filter"]["scan_errors"]:
        raise AssertionError("Optimized CCM reported failed lag scans.")
    arrays = [
        "rho_null", "cuts", "library_seeds", "lag_rho", "lag_steps",
        "lag_library_seeds", "lag_common_anchor_indices",
    ]
    for name in arrays:
        np.testing.assert_array_equal(actual["raw"][name], expected["raw"][name], err_msg=name)
    pd.testing.assert_frame_equal(actual["edges"], expected["edges"], check_exact=True)
    for name in ["lag_profile", "lag_peaks"]:
        pd.testing.assert_frame_equal(actual["raw"][name], expected["raw"][name], check_exact=True)
    return {
        "method": "CCM", "compared_arrays": arrays,
        "lag_array_shape": list(actual["raw"]["lag_rho"].shape),
        "arrays_exactly_equal": True, "edge_tables_exactly_equal": True,
        "lag_profiles_and_peaks_exactly_equal": True,
        "tp0_shared_by_construction": True,
        "tp0_seconds": original_seconds,
        "reference_scan_seconds": reference_seconds,
        "optimized_scan_seconds": optimized_seconds,
        "scan_speedup": reference_seconds / optimized_seconds,
        "total_speedup": (original_seconds + reference_seconds) / (original_seconds + optimized_seconds),
    }


def check_cache_identity(series, cache_path, settings):
    """Ensure both preparation and consumption reject mismatched cache keys."""
    sample_size = len(series) - 2 * settings["tau_max"]
    cache_settings = {
        "sample_sizes": [sample_size], "sig_samples": GPDC_NULL_SAMPLES, "seed": METHOD_SEED,
    }
    start = time.perf_counter()
    prepare_gpdc_null(cache_path, **cache_settings)
    reuse_seconds = time.perf_counter() - start
    checks = {}
    for key, value in [
        ("seed", METHOD_SEED + 1), ("sig_samples", GPDC_NULL_SAMPLES - 1),
        ("sample_sizes", [sample_size - 1]),
    ]:
        changed = dict(cache_settings)
        changed[key] = value
        try:
            prepare_gpdc_null(cache_path, **changed)
        except ValueError:
            checks[f"prepare_rejects_{key}"] = True
        else:
            raise AssertionError(f"Preparation accepted a mismatched {key}.")
    for key, value in [
        ("seed", METHOD_SEED + 1), ("sig_samples", GPDC_NULL_SAMPLES - 1),
        ("tau_max", MAX_LAG + 1),
    ]:
        changed = dict(settings)
        changed[key] = value
        try:
            run_pcmciplus(series, null_dist_filename=cache_path, **changed)
        except ValueError:
            checks[f"consumer_rejects_{key}"] = True
        else:
            raise AssertionError(f"PCMCI+ accepted a mismatched {key}.")
    return {"matching_cache_reuse_seconds": reuse_seconds, "negative_checks": checks}


def main():
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output = ROOT / "data" / "benchmarks" / "optimization_checks" / timestamp
    output.mkdir(parents=True, exist_ok=False)
    report_path = output / "report.json"
    report = {
        "status": "running", "created_utc": timestamp,
        "platform": platform.platform(), "machine": platform.machine(),
        "python": sys.version, "numeric_threads": 1,
        "packages": {name: version(name) for name in [
            "numpy", "pandas", "scipy", "scikit-learn", "tigramite", "dcor", "pyEDM",
        ]},
        "parameters": {name: value for name, value in globals().items()
                       if name.isupper() and name != "ROOT"},
        "comparison": "Exact array and DataFrame equality; no numerical tolerance.",
        "timing_note": "Single sequential pass; cache preparation is excluded from warm timings.",
        "cases": [],
    }
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Acceptance report: {report_path}", flush=True)
    settings = {
        "tau_max": MAX_LAG, "pc_alpha": ALPHA, "alpha": ALPHA,
        "sig_samples": GPDC_NULL_SAMPLES, "seed": METHOD_SEED,
        "gp_params": {"optimizer": None},
    }
    start_all = time.perf_counter()
    try:
        with threadpool_limits(limits=1):
            cases = generate_inputs()
            for name, series in cases.items():
                series.to_csv(output / f"{name}_input.csv", index=False)
            sample_size = len(cases["linear"]) - 2 * MAX_LAG
            print("Preparing the native GPDC reference null...", flush=True)
            start = time.perf_counter()
            cache_path = prepare_gpdc_null(
                output / "gpdc_null.npz", sample_sizes=[sample_size],
                sig_samples=GPDC_NULL_SAMPLES, seed=METHOD_SEED,
            )
            report["gpdc_preparation_seconds"] = time.perf_counter() - start
            report["cache_identity"] = check_cache_identity(cases["linear"], cache_path, settings)
            for name, series in cases.items():
                print(f"{name}: comparing complete PCMCI+ outputs...", flush=True)
                pcmci_check = check_pcmci(series, cache_path, settings)
                report["cases"].append({"process": name, **pcmci_check})
                report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
                print(f"{name}: comparing the full 100-library x 31-lag CCM scan...", flush=True)
                ccm_check = check_ccm(series)
                report["cases"].append({"process": name, **ccm_check})
                report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
                print(f"{name}: all comparisons passed exactly.", flush=True)
        report["status"] = "passed"
    except Exception as error:
        report["status"] = "failed"
        report["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        report["elapsed_seconds"] = time.perf_counter() - start_all
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"All optimization checks passed. Report: {report_path}", flush=True)


if __name__ == "__main__":
    main()
