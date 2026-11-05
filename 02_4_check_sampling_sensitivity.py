"""Paired sampling controls with exact ages and the fixed native causal pipeline.

See agent/02_4_sampling_sensitivity_plan.md for the complete comparison contract.
Run through agent/run_experiment.py. Results describe sampling sensitivity of
one history per process form, not error rates across independent histories.
"""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# Apply these before imports in both the parent and Windows spawned workers.
for thread_name in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"]:
    os.environ[thread_name] = "1"
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "data/cache/matplotlib"))

from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
from importlib.metadata import version
import json
import multiprocessing
import platform
import time

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from toolbox.process import generate_process
from toolbox.observation import sampling_indices_from_gaps
from toolbox.causal import resample_pair, prepare_gpdc_null, run_pcmciplus, run_lkif, run_ccm, run_granger


RUN_CAUSAL = True              # False only summarizes and plots existing committed results.
PREPARE_ONLY = False           # True writes inputs and sampling diagnostics, without causal fits.
RESUME = True                  # Preserve all complete cases, including explicitly failed fits.
N_WORKERS = 32                 # Measured fastest on this Windows PC; one numerical thread each.
N_MASK_GROUPS = 10             # Independent A/B sampling pairs, not independent climate histories.
SEED = 20261009                # One process realization per form, paired innovations as in 02.
MASK_SEED = 20261104           # Separate random stream for sampling arrangements.
METHOD_SEED = SEED + 300       # Shared native-null and CCM library seed across paired designs.
DURATION = 30.0               # Full retained process duration, kyr.
DT = 0.020                    # Native and analysis grid step, kyr: 20 years.
TRUE_LAG = 0.020              # Model X->Y delay, kyr; independent process has no coupling.
MEMORY_X = 0.5                # Driver memory, kyr.
MEMORY_Y = 0.2                # Response memory, kyr.
BACKGROUND_X = 1 / np.sqrt(1 - 0.95**2)  # Unforced stationary driver SD.
BACKGROUND_Y = 1 / np.sqrt(1 - 0.8**2)   # Unforced stationary response SD.
COUPLING_GAIN = 4.0           # Long-run response gain K, not an estimator score.
BURN_TIMESCALES = 5           # Warm-up in units of the longest memory time.
ANALYSIS_OLD = 29.0           # Fixed old endpoint, kyr BP; observed in every sampling design.
ANALYSIS_YOUNG = 1.0          # Fixed young endpoint, kyr BP; observed in every design.
REGULAR_GAPS_YEARS = [20, 40, 80, 160]  # Synchronous density controls, years.
IRREGULAR_GAPS_YEARS = {"cv025": [60, 100], "cv075": [20, 140]}  # Equal counts; mean 80 yr.
GRADIENT_STRENGTH = 2.0       # Dimensionless log-weight slope; positive favors long gaps at old end.
AGE_BIN_WIDTH = 4.0           # Width of diagnostic age bins, kyr; seven bins cover the window.

MAX_LAG = 5                  # PCMCI+/LKIF maximum history, 20-year grid steps.
GC_ORDER = 1                 # Fixed Granger order, 20-year grid steps.
ALPHA = 0.05                 # Native/significance threshold; wrapper correction rules unchanged.
GPDC_NULL_SAMPLES = 10000     # Native GPDC null draws, prepared once for T=1391.
GPDC_GP_PARAMS = {"optimizer": None}  # Disable kernel hyperparameter optimization, as in 02.
CCM_E = 3                    # Embedding dimension.
CCM_TAU = -1                 # Backward embedding interval in grid steps.
CCM_EXCLUSION_RADIUS = 10    # Neighbor exclusion, grid steps: 200 years.
CCM_LIBRARY_FRACTION = 0.8    # Fraction of usable embedded rows in each library.
CCM_SAMPLES = 10             # Libraries for the original Tp=0 score.
CCM_SURROGATES = 100         # Cut-and-swap null records per direction.
CCM_LAG_STEPS = list(range(-15, 16))  # Full -300..300yr scan, grid steps.
CCM_LAG_LIBRARIES = 100      # Paired random libraries at every scanned shift.
CCM_LAG_SIGN_THRESHOLD = 0.9 # Existing temporal compatibility rule, not a calibrated probability.

STEM = Path(__file__).stem
DATA_DIR = ROOT / "data/processed" / STEM
FIGURE_DIR = ROOT / "figures" / STEM
CACHE_DIR = ROOT / "data/cache/gpdc"
PROCESSES = ["linear", "nonlinear", "independent"]
METHODS = ["PCMCI+", "LKIF", "CCM", "GC"]
IMPLEMENTATION_REVISION = "fixed_gap_paired_sampling_v1"


def execution_parameters():
    """Record scientific settings, packages and the local numerical platform."""
    parameters = {name: value for name, value in globals().items()
                  if name.isupper() and isinstance(value, (str, int, float, list, dict))}
    parameters["packages"] = {name: version(name) for name in [
        "numpy", "scipy", "pandas", "scikit-learn", "tigramite", "pyEDM",
        "LK-Info-Flow", "statsmodels", "dcor", "threadpoolctl", "numba", "llvmlite"]}
    blas = np.__config__.CONFIG.get("Build Dependencies", {}).get("blas", {})
    parameters["environment"] = {"system": platform.system(), "machine": platform.machine(),
                                  "python": platform.python_version(), "numpy_blas": blas}
    parameters["interpretation"] = "One history per form; mask groups measure conditional sampling sensitivity."
    parameters["boundary_condition"] = "Both exact-age analysis-window endpoints are observed."
    return json.loads(json.dumps(parameters))


def check_parameters(parameters):
    path = DATA_DIR / "parameters.json"
    if path.exists():
        if not RESUME and RUN_CAUSAL:
            raise ValueError("Outputs already exist. Resume this run or preserve it before changing settings.")
        previous = json.loads(path.read_text(encoding="utf-8"))
        controls = {"RUN_CAUSAL", "PREPARE_ONLY", "RESUME", "N_WORKERS"}
        changed = [key for key in set(parameters) | set(previous)
                   if key not in controls and parameters.get(key) != previous.get(key)]
        if changed:
            raise ValueError(f"Saved experiment differs in {sorted(changed)}; preserve it before a fresh run.")
    elif any((DATA_DIR / "members").glob("*_edges.csv")):
        raise ValueError("Committed results lack their scientific parameters.")


def write_table(table, path):
    """Commit a complete table atomically so a partial write is never a vote."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".partial.csv")
    table.to_csv(temporary, index=False)
    temporary.replace(path)


def save_input(series, path):
    if path.exists():
        previous = pd.read_csv(path)
        if (list(previous.columns) != list(series.columns) or previous.shape != series.shape
                or not np.allclose(previous.to_numpy(), series.to_numpy(), rtol=0, atol=1e-12)):
            raise ValueError(f"Regenerated input changed: {path.name}")
        return
    write_table(series, path)


def generate_originals():
    """Generate the same physical processes as 02, then crop without resampling."""
    n_steps = int(round(DURATION / DT)) + 1
    age_grid = np.linspace(ANALYSIS_OLD, ANALYSIS_YOUNG,
                           int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / DT)) + 1)
    indices = np.rint((DURATION - age_grid) / DT).astype(int)
    ar_x = np.exp(-DT / MEMORY_X)
    ar_y = np.exp(-DT / MEMORY_Y)
    noise = (BACKGROUND_X * np.sqrt(1 - ar_x**2), BACKGROUND_Y * np.sqrt(1 - ar_y**2))
    burn_in = int(np.ceil(BURN_TIMESCALES * max(MEMORY_X, MEMORY_Y) / DT))
    originals = {}
    for process in PROCESSES:
        gain = 0.0 if process == "independent" else COUPLING_GAIN
        generated = generate_process(
            n_steps=n_steps, dt=DT, age_start=DURATION, lag=TRUE_LAG,
            ar_x=ar_x, ar_y=ar_y, coupling=(1 - ar_y) * gain,
            coupling_form="nonlinear" if process == "nonlinear" else "linear",
            noise_std=noise, burn_in=burn_in, rng=SEED)
        save_input(generated["series"], DATA_DIR / f"{process}_generated.csv")
        original = generated["series"].iloc[indices].copy().reset_index(drop=True)
        original["time"] = age_grid[0] - age_grid
        original["age"] = age_grid
        save_input(original, DATA_DIR / f"{process}_original.csv")
        originals[process] = original
    return originals


def prepare_sampling_designs(age):
    """Build one canonical copy of each sampling mask, shared across processes."""
    span_steps = len(age) - 1
    clocks, metadata, point_tables = {}, [], []

    def add_clock(clock_id, indices, gap_label, group, mask, gradient, seed):
        gaps = -np.diff(age[indices]) * 1000
        clocks[clock_id] = indices
        metadata.append({"clock_id": clock_id, "gap_label": gap_label, "mask_group": group,
                         "mask": mask, "gradient": gradient, "seed": seed,
                         "n_observations": len(indices), "mean_gap_years": gaps.mean(),
                         "gap_cv": gaps.std() / gaps.mean(), "minimum_gap_years": gaps.min(),
                         "maximum_gap_years": gaps.max()})
        point_tables.append(pd.DataFrame({"clock_id": clock_id, "native_index": indices,
                                          "age": age[indices]}))

    for years in REGULAR_GAPS_YEARS:
        step = years / (1000 * DT)
        if not float(step).is_integer() or step < 1 or span_steps % int(step):
            raise ValueError("Regular spacings must divide the analysis span in native steps.")
        gaps = np.full(span_steps // int(step), int(step), dtype=int)
        indices = sampling_indices_from_gaps(age, gap_steps=gaps, order="given")
        add_clock(f"regular_{years}", indices, f"regular_{years}", -1, "R", "none", -1)
    for gap_label, gap_years in IRREGULAR_GAPS_YEARS.items():
        steps = np.asarray(gap_years) / (1000 * DT)
        if (len(steps) != 2 or not np.all(np.isfinite(steps)) or np.any(steps < 1)
                or len(np.unique(steps)) != 2 or np.any(steps != np.rint(steps))
                or not np.isclose(np.mean(gap_years), 80.0)
                or span_steps % int(steps.sum())):
            raise ValueError("Two distinct positive native-step gaps must average 80 yr and cover the window.")
        copies = span_steps // int(steps.sum())
        gap_steps = np.repeat(steps.astype(int), copies)
        for group in range(N_MASK_GROUPS):
            for offset, mask in enumerate(["A", "B"]):
                seed = MASK_SEED + 100 * group + offset
                for gradient, order in [("none", "random"), ("old_sparse", "old_sparse"),
                                        ("young_sparse", "young_sparse")]:
                    indices = sampling_indices_from_gaps(
                        age, gap_steps=gap_steps, order=order,
                        gradient_strength=GRADIENT_STRENGTH, rng=seed)
                    clock_id = f"{gap_label}_g{group:02d}_{mask}_{gradient}"
                    add_clock(clock_id, indices, gap_label, group, mask, gradient, seed)
    return clocks, pd.DataFrame(metadata), pd.concat(point_tables, ignore_index=True)


def build_jobs(clock_metadata):
    jobs = []
    for process in PROCESSES:
        for years in REGULAR_GAPS_YEARS:
            scenario = f"density_{years}"
            jobs.append({"case_id": f"r0_{process}_{scenario}", "replicate": 0, "process": process,
                         "family": "density", "scenario": scenario, "gap_label": f"regular_{years}",
                         "mean_gap_years": float(years), "gap_cv": 0.0, "gradient": "none",
                         "pairing": "RR", "mask_group": -1, "clock_x": f"regular_{years}",
                         "clock_y": f"regular_{years}"})
        for gap_label in IRREGULAR_GAPS_YEARS:
            cv = float(clock_metadata.loc[clock_metadata.gap_label.eq(gap_label), "gap_cv"].iloc[0])
            for gradient in ["none", "old_sparse", "young_sparse"]:
                pairings = ["AA", "BB", "AB", "BA"] if gradient == "none" else ["AA", "BB"]
                for pairing in pairings:
                    family = "irregularity" if gradient == "none" else "gradient"
                    scenario = (f"irregular_{gap_label}_{pairing}" if gradient == "none"
                                else f"gradient_{gap_label}_{gradient}_{pairing}")
                    for group in range(N_MASK_GROUPS):
                        jobs.append({"case_id": f"r0_{process}_{scenario}_g{group:02d}", "replicate": 0,
                                     "process": process, "family": family, "scenario": scenario,
                                     "gap_label": gap_label, "mean_gap_years": 80.0, "gap_cv": cv,
                                     "gradient": gradient, "pairing": pairing, "mask_group": group,
                                     "clock_x": f"{gap_label}_g{group:02d}_{pairing[0]}_{gradient}",
                                     "clock_y": f"{gap_label}_g{group:02d}_{pairing[1]}_{gradient}"})
    return jobs


def prepare_inputs():
    originals = generate_originals()
    age = originals[PROCESSES[0]]["age"].to_numpy()
    clocks, clock_metadata, points = prepare_sampling_designs(age)
    jobs = build_jobs(clock_metadata)
    existing = DATA_DIR / "cases.csv"
    manifest = pd.DataFrame(jobs)
    if existing.exists():
        pd.testing.assert_frame_equal(pd.read_csv(existing), manifest, check_dtype=False,
                                      rtol=0, atol=1e-12)
    write_table(manifest, existing)
    write_table(clock_metadata, DATA_DIR / "sampling_clocks.csv")
    write_table(points, DATA_DIR / "sampling_points.csv")
    diagnostics, bin_rows = [], []
    age_bins = np.arange(ANALYSIS_YOUNG, ANALYSIS_OLD + AGE_BIN_WIDTH / 2, AGE_BIN_WIDTH)
    if not np.isclose(age_bins[-1], ANALYSIS_OLD):
        raise ValueError("Diagnostic age bins must exactly cover the analysis window.")
    for clock_id, indices in clocks.items():
        counts, _ = np.histogram(age[indices], bins=age_bins)
        for left, right, count in zip(age_bins[:-1], age_bins[1:], counts):
            bin_rows.append({"clock_id": clock_id, "young_age": left, "old_age": right,
                             "age_midpoint": (left + right) / 2, "n_observations": int(count),
                             "observations_per_kyr": count / (right - left)})
    for job in jobs:
        original = originals[job["process"]]
        ix, iy = clocks[job["clock_x"]], clocks[job["clock_y"]]
        series = resample_pair(age[ix], original.x.to_numpy()[ix],
                               age[iy], original.y.to_numpy()[iy], age_grid=age)["series"]
        save_input(series, DATA_DIR / "inputs" / f"{job['case_id']}.csv")
        shared = len(np.intersect1d(ix, iy)) / len(np.union1d(ix, iy))
        for record, indices, other_indices in [("x", ix, iy), ("y", iy, ix)]:
            direct, interpolated = original[record].to_numpy(), series[record].to_numpy()
            observed_age = age[indices]
            gaps = -np.diff(observed_age) * 1000
            nearest_years = np.min(np.abs(observed_age[:, None] - age[other_indices][None, :]), axis=1) * 1000
            diagnostics.append({**job, "record": record, "n_observations": len(indices),
                                "fraction_grid_interpolated": 1 - len(indices) / len(age),
                                "actual_mean_gap_years": gaps.mean(), "actual_gap_cv": gaps.std() / gaps.mean(),
                                "maximum_gap_years": gaps.max(), "shared_observation_jaccard": shared,
                                "mean_nearest_other_years": nearest_years.mean(),
                                "max_nearest_other_years": nearest_years.max(),
                                "rmse_over_original_sd": np.sqrt(np.mean((interpolated - direct)**2)) / np.std(direct),
                                "correlation": np.corrcoef(direct, interpolated)[0, 1],
                                "sd_ratio": np.std(interpolated) / np.std(direct),
                                "increment_variance_ratio": np.var(np.diff(interpolated)) / np.var(np.diff(direct)),
                                "lag1_original": np.corrcoef(direct[:-1], direct[1:])[0, 1],
                                "lag1_interpolated": np.corrcoef(interpolated[:-1], interpolated[1:])[0, 1]})
    diagnostics = pd.DataFrame(diagnostics)
    write_table(diagnostics, DATA_DIR / "sampling_diagnostics.csv")
    write_table(pd.DataFrame(bin_rows), DATA_DIR / "sampling_age_bins.csv")
    return jobs, originals, clock_metadata, points, diagnostics


def gpdc_null_path():
    effective_n = int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / DT)) + 1 - 2 * MAX_LAG
    return CACHE_DIR / f"null_T{effective_n}_n{GPDC_NULL_SAMPLES}_seed{METHOD_SEED}.npz"


def call_method(series, method, settings, null_path):
    """Thin native calls with explicit frozen settings; no conditional-fit cache."""
    if method == "PCMCI+":
        return run_pcmciplus(series, tau_max=settings["MAX_LAG"], pc_alpha=settings["ALPHA"],
                             alpha=settings["ALPHA"], sig_samples=settings["GPDC_NULL_SAMPLES"],
                             seed=settings["METHOD_SEED"], gp_params=settings["GPDC_GP_PARAMS"],
                             null_dist_filename=str(null_path))
    if method == "LKIF":
        return run_lkif(series, max_lag=settings["MAX_LAG"], alpha=settings["ALPHA"])
    if method == "CCM":
        embedded_rows = len(series) - (settings["CCM_E"] - 1) * abs(settings["CCM_TAU"])
        return run_ccm(series, E=settings["CCM_E"], tau=settings["CCM_TAU"],
                       exclusion_radius=settings["CCM_EXCLUSION_RADIUS"],
                       lib_size=int(settings["CCM_LIBRARY_FRACTION"] * embedded_rows),
                       sample=settings["CCM_SAMPLES"], n_surrogates=settings["CCM_SURROGATES"],
                       seed=settings["METHOD_SEED"], alpha=settings["ALPHA"],
                       lag_steps=settings["CCM_LAG_STEPS"], n_libraries=settings["CCM_LAG_LIBRARIES"],
                       sign_threshold=settings["CCM_LAG_SIGN_THRESHOLD"])
    if method == "GC":
        return run_granger(series, max_lag=settings["GC_ORDER"], alpha=settings["ALPHA"])
    raise ValueError(f"Unknown method: {method}")


def native_snapshot(value):
    """Make native arrays, nested results and OLS fits portable without pickling."""
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
        # Native GC OLS objects contain lazy caches; save their scientific values.
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


def fit_case(job, settings, output_dir, null_path):
    """One worker owns all files for a case; the edge table is its final commit."""
    output_dir = Path(output_dir)
    series = pd.read_csv(output_dir / "inputs" / f"{job['case_id']}.csv")
    tables = []
    with threadpool_limits(limits=1):
        for method in METHODS:
            started = time.perf_counter()
            try:
                result = call_method(series, method, settings, null_path)
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
            # Serialization/I/O failures stop the run, rather than becoming scientific failures.
            native_path = output_dir / "native" / f"{job['case_id']}_{method.replace('+', 'plus')}.json.gz"
            native_path.parent.mkdir(parents=True, exist_ok=True)
            temporary = native_path.with_suffix(".partial.gz")
            with gzip.open(temporary, "wt", encoding="utf-8") as handle:
                json.dump(native_snapshot(result), handle, separators=(",", ":"))
            temporary.replace(native_path)
            edges["method"] = method
            edges["seconds"] = time.perf_counter() - started
            for key, value in job.items():
                edges[key] = value
            tables.append(edges)
    result = pd.concat(tables, ignore_index=True)
    write_table(result, output_dir / "members" / f"{job['case_id']}_edges.csv")
    return result


def read_completed_case(job):
    path = DATA_DIR / "members" / f"{job['case_id']}_edges.csv"
    if not path.exists():
        return None
    # Preserve the committed float64 values exactly when rebuilding aggregate tables.
    result = pd.read_csv(path, float_precision="round_trip", dtype={
        "edge_detected": "boolean", "unoriented": "boolean", "raw_edge_detected": "boolean"})
    keys = ["method", "source", "target"]
    expected = {(method, source, target) for method in METHODS for source, target in [("x", "y"), ("y", "x")]}
    required = set(keys + ["status", "edge_detected", "unoriented", "score", "p_value", *job])
    if not required.issubset(result) or len(result) != 8 or result.duplicated(keys).any():
        return None
    if set(result[keys].itertuples(index=False, name=None)) != expected or result.status.isna().any():
        return None
    if result.status.eq("pending").any():
        return None
    for key, value in job.items():
        matches = (np.allclose(result[key], value, rtol=0, atol=1e-12)
                   if isinstance(value, float) else result[key].eq(value).all())
        if not matches:
            raise ValueError(f"Saved case identity differs in {key}: {path.name}")
    for method in METHODS:
        if not (DATA_DIR / "native" / f"{job['case_id']}_{method.replace('+', 'plus')}.json.gz").is_file():
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


def write_status(state, completed, requested, started, error=""):
    status = {"state": state, "completed_pairs": completed, "requested_pairs": requested,
              "n_workers": N_WORKERS, "n_mask_groups": N_MASK_GROUPS, "n_process_histories_per_form": 1,
              "elapsed_seconds": time.perf_counter() - started,
              "updated_utc": datetime.now(timezone.utc).isoformat(), "error": error}
    path = DATA_DIR / "run_status.json"
    temporary = path.with_suffix(".partial.json")
    temporary.write_text(json.dumps(status, indent=2), encoding="utf-8")
    temporary.replace(path)


def plot_sampling_designs(originals, clock_metadata, points):
    """Display the actual BP-age density and one reconstruction, not target jitter."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    bins = pd.read_csv(DATA_DIR / "sampling_age_bins.csv")
    fig, axes = plt.subplots(2, 3, figsize=(13, 6), sharex=True, sharey=True)
    for row, gap_label in enumerate(IRREGULAR_GAPS_YEARS):
        for column, gradient in enumerate(["none", "old_sparse", "young_sparse"]):
            selected = clock_metadata[(clock_metadata.gap_label == gap_label) & (clock_metadata.gradient == gradient)]
            densities = bins[bins.clock_id.isin(selected.clock_id)]
            ax = axes[row, column]
            for _, group in densities.groupby("clock_id", sort=False):
                ax.plot(group.age_midpoint, group.observations_per_kyr, color="tab:blue", alpha=0.12, lw=0.7)
            mean = densities.groupby("age_midpoint").observations_per_kyr.mean()
            ax.plot(mean.index, mean.values, color="tab:blue", lw=2)
            ax.axhline(1000 / 80, color="0.4", ls="--", lw=0.8)
            ax.set_xlim(ANALYSIS_OLD, ANALYSIS_YOUNG)
            ax.set_title(f"{gap_label}: {gradient.replace('_', ' ')}")
            ax.set_xlabel("Age (kyr BP)")
            if column == 0:
                ax.set_ylabel("Observations / kyr")
    fig.suptitle("Fixed 351 observations and gap multiset; only their age distribution changes")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    for extension in ["png", "pdf"]:
        fig.savefig(FIGURE_DIR / f"01_sampling_age_density.{extension}", dpi=180)
    plt.close(fig)
    original = originals["linear"]
    fig, axes = plt.subplots(2, 3, figsize=(14, 6), sharex=True)
    for column, (gradient, pairing) in enumerate([("none", "AA"), ("none", "AB"), ("old_sparse", "AA")]):
        for row, record in enumerate(["x", "y"]):
            clock_id = f"cv075_g00_{pairing[row]}_{gradient}"
            indices = points.loc[points.clock_id.eq(clock_id), "native_index"].to_numpy(dtype=int)
            age = original.age.to_numpy()
            values = original[record].to_numpy()
            interpolated = np.interp(-age, -age[indices], values[indices])
            ax = axes[row, column]
            ax.plot(age, values, color="0.3", lw=0.9, label="Original uniform")
            ax.plot(age, interpolated, color="tab:orange", lw=1.0, label="Linear interpolation")
            ax.scatter(age[indices], values[indices], color="tab:orange", s=12, label="Actual observations")
            ax.set_xlim(24, 22)
            ax.set_title(f"{record.upper()}: {gradient.replace('_', ' ')}, {pairing}")
            ax.set_xlabel("Age (kyr BP)")
            ax.set_ylabel("Proxy value")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Linear process, CV=0.75, mask group 0; exact ages and a common 20-year grid")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    for extension in ["png", "pdf"]:
        fig.savefig(FIGURE_DIR / f"02_sampling_reconstruction.{extension}", dpi=180)
    plt.close(fig)


def summarize_sampling_results(results):
    """Describe sampling masks conditional on one history, retaining all requests.

    Results must include pending rows for every requested case and method.
    Decision validity requires both directions; finite successful native scores
    remain usable separately, including CCM scores whose lag screen failed.
    """
    metadata = ["case_id", "replicate", "process", "family", "scenario", "gap_label",
                "mean_gap_years", "gap_cv", "gradient", "pairing", "mask_group",
                "clock_x", "clock_y", "method"]
    keys = ["case_id", "method", "source", "target"]
    required = metadata + ["source", "target", "score", "status", "edge_detected", "unoriented"]
    if not set(required).issubset(results.columns) or results.empty:
        raise ValueError("Sampling results require metadata and all requested directional rows.")
    if results.duplicated(keys).any() or results["status"].isna().any():
        raise ValueError("Directional rows must be unique and have explicit status.")
    if not results.groupby("case_id")["method"].nunique().eq(results["method"].nunique()).all():
        raise ValueError("Every requested case needs all methods, including pending placeholders.")
    for column in ["edge_detected", "unoriented", "raw_edge_detected"]:
        if column in results and any(not isinstance(value, (bool, np.bool_))
                                     for value in results[column].dropna()):
            raise ValueError(f"{column} must contain booleans or missing values.")
    state_rows = []
    for (_, method), pair in results.groupby(["case_id", "method"], sort=False):
        if set(zip(pair["source"], pair["target"])) != {("x", "y"), ("y", "x")}:
            raise ValueError("Every requested case and method needs both directional rows.")
        if any(pair[column].nunique(dropna=False) != 1 for column in metadata):
            raise ValueError("The two directions disagree on case metadata.")
        row = pair.iloc[0][metadata].to_dict()
        row["pending"] = bool(pair["status"].eq("pending").any())
        row["unresolved"] = bool(pair["unoriented"].fillna(False).any())
        for source, target in [("x", "y"), ("y", "x")]:
            edge = pair[pair["source"] == source].iloc[0]
            direction = f"{source}_to_{target}"
            row[direction] = edge["edge_detected"] if edge["status"] == "ok" else pd.NA
            row[f"raw_{direction}"] = (edge.get("raw_edge_detected", pd.NA)
                                        if method == "CCM" and edge["status"] == "ok" else pd.NA)
        row["valid"] = all(pd.notna(row[name]) for name in ["x_to_y", "y_to_x"])
        row["raw_valid"] = all(pd.notna(row[f"raw_{name}"]) for name in ["x_to_y", "y_to_x"])
        row["failed"] = not row["pending"] and not row["valid"]
        row["correct"] = False
        if row["valid"] and not row["unresolved"]:
            forward, reverse = bool(row["x_to_y"]), bool(row["y_to_x"])
            row["correct"] = (not forward and not reverse if row["process"] == "independent"
                              else forward and not reverse)
        row["ccm_lag_failed"] = bool(method == "CCM" and "lag_status" in pair
                                     and pair["lag_status"].eq("failed").any())
        state_rows.append(row)
    states = pd.DataFrame(state_rows)
    for column in ["x_to_y", "y_to_x", "raw_x_to_y", "raw_y_to_x"]:
        states[column] = states[column].astype("boolean")

    group_keys = ["replicate", "process", "scenario", "method"]
    state_groups = dict(tuple(states.groupby(group_keys, sort=False)))
    summary_rows = []
    for key, group in results.groupby(group_keys, sort=False):
        paired = state_groups[key]
        valid = paired[paired["valid"]]
        raw_valid = paired[paired["raw_valid"]]
        row = group.iloc[0][[name for name in metadata if name not in
                            ["case_id", "mask_group", "clock_x", "clock_y"]]].to_dict()
        row.update({"n_requested": len(paired), "n_valid": len(valid),
                    "n_pending": int(paired["pending"].sum()),
                    "n_failed": int(paired["failed"].sum()),
                    "n_unresolved": int(paired["unresolved"].sum()),
                    "n_correct_graph": int(paired["correct"].sum()),
                    "n_raw_valid": len(raw_valid),
                    "n_ccm_lag_failed": int(paired["ccm_lag_failed"].sum())})
        row["valid_fraction"] = len(valid) / len(paired)
        row["fraction_correct_requested"] = row["n_correct_graph"] / len(paired)
        for direction in ["x_to_y", "y_to_x"]:
            row[f"n_retained_{direction}"] = int(valid[direction].sum())
            row[f"q_{direction}"] = float(valid[direction].mean()) if len(valid) else np.nan
            row[f"n_raw_retained_{direction}"] = int(raw_valid[f"raw_{direction}"].sum())
            row[f"raw_q_{direction}"] = (float(raw_valid[f"raw_{direction}"].mean())
                                          if len(raw_valid) else np.nan)
            source = direction[0]
            scores = pd.to_numeric(group.loc[(group["source"] == source)
                                             & group["status"].eq("ok"), "score"], errors="coerce")
            scores = scores[np.isfinite(scores)]
            quartiles = np.percentile(scores, [25, 50, 75]) if len(scores) else [np.nan] * 3
            row[f"n_scores_{direction}"] = len(scores)
            row[f"q25_{direction}"], row[f"median_{direction}"], row[f"q75_{direction}"] = quartiles
            row[f"iqr_{direction}"] = quartiles[2] - quartiles[0]
        summary_rows.append(row)
    summaries = pd.DataFrame(summary_rows)

    # A/B permutations are matched within a mask group, never independent histories.
    contrast_specs = [
        ("asynchrony", [("irregularity", "none", "AB", 0.5), ("irregularity", "none", "BA", 0.5),
                        ("irregularity", "none", "AA", -0.5), ("irregularity", "none", "BB", -0.5)]),
        ("old_sparse", [("gradient", "old_sparse", "AA", 0.5), ("gradient", "old_sparse", "BB", 0.5),
                        ("irregularity", "none", "AA", -0.5), ("irregularity", "none", "BB", -0.5)]),
        ("young_sparse", [("gradient", "young_sparse", "AA", 0.5), ("gradient", "young_sparse", "BB", 0.5),
                          ("irregularity", "none", "AA", -0.5), ("irregularity", "none", "BB", -0.5)]),
    ]
    contrast_keys = ["replicate", "process", "gap_label", "mask_group", "method", "source", "target"]
    contrast_rows = []
    masks = results[results["family"].isin(["irregularity", "gradient"])]
    for key, group in masks.groupby(contrast_keys, sort=False):
        for name, terms in contrast_specs:
            row = dict(zip(contrast_keys, key))
            row.update({"contrast": name, "n_terms_requested": 4, "n_terms_present": 0,
                        "n_terms_valid": 0, "n_terms_pending": 0, "delta_score": np.nan})
            weighted_scores = []
            for family, gradient, pairing, weight in terms:
                term = group[(group["family"] == family) & (group["gradient"] == gradient)
                             & (group["pairing"] == pairing)]
                if len(term) > 1:
                    raise ValueError("A paired contrast contains duplicate sampling cases.")
                if term.empty:
                    continue
                row["n_terms_present"] = row["n_terms_present"] + 1
                edge = term.iloc[0]
                row["n_terms_pending"] = row["n_terms_pending"] + int(edge["status"] == "pending")
                score = pd.to_numeric(edge["score"], errors="coerce")
                if edge["status"] == "ok" and pd.notna(score) and np.isfinite(score):
                    row["n_terms_valid"] = row["n_terms_valid"] + 1
                    weighted_scores.append(weight * float(score))
            row["n_terms_missing"] = 4 - row["n_terms_present"]
            row["n_terms_failed"] = row["n_terms_present"] - row["n_terms_valid"] - row["n_terms_pending"]
            row["status"] = ("ok" if row["n_terms_valid"] == 4 else "pending" if row["n_terms_pending"]
                             else "missing" if row["n_terms_missing"] else "failed")
            if row["status"] == "ok":
                row["delta_score"] = sum(weighted_scores)
            contrast_rows.append(row)
    return summaries, states, pd.DataFrame(contrast_rows)


def _save_sampling_figure(fig, output_dir, name):
    output_dir.mkdir(parents=True, exist_ok=True)
    for extension in ["png", "pdf"]:
        fig.savefig(output_dir / f"{name}.{extension}", dpi=180)
    plt.close(fig)


def plot_sampling_results(results, summaries, contrasts, output_dir):
    """Plot native scores and paired mask contrasts; spreads are not confidence intervals."""
    if results["replicate"].nunique() != 1:
        raise ValueError("Plot one independent process realization at a time.")
    processes = ["linear", "nonlinear", "independent"]
    methods = ["PCMCI+", "LKIF", "CCM", "GC"]
    colors = {"x": "tab:blue", "y": "tab:orange"}
    direction_labels = {"x": "X → Y", "y": "Y → X"}
    score_labels = ["Max. GPDC", "Max. |normalized flow|", "Original Tp=0 rho", "Log residual variance ratio"]
    direction_legend = [Line2D([], [], color=colors[source], marker="o", label=direction_labels[source])
                        for source in ["x", "y"]]

    fig, axes = plt.subplots(3, 4, figsize=(14, 8.5), sharex=True, sharey="col")
    density = results[results["family"] == "density"]
    for row, process in enumerate(processes):
        for column, method in enumerate(methods):
            ax = axes[row, column]
            group = density[(density["process"] == process) & (density["method"] == method)]
            for source in ["x", "y"]:
                edges = group[group["source"] == source].sort_values("mean_gap_years")
                scores = pd.to_numeric(edges["score"], errors="coerce").to_numpy(dtype=float)
                usable = edges["status"].eq("ok").to_numpy() & np.isfinite(scores)
                ax.plot(edges["mean_gap_years"], np.where(usable, scores, np.nan),
                        color=colors[source], marker="o", markersize=4)
                for _, edge in edges.loc[~usable].iterrows():
                    label = "pending" if edge["status"] == "pending" else "failed"
                    ax.text(edge["mean_gap_years"], 0.04 if source == "x" else 0.18, label,
                            transform=ax.get_xaxis_transform(), color=colors[source], fontsize=6, rotation=60)
            ax.set_xscale("log", base=2)
            ax.set_xticks([20, 40, 80, 160], labels=["20", "40", "80", "160"])
            ax.grid(alpha=0.2)
            if row == 0:
                ax.set_title(f"{method}: {score_labels[column]}", fontsize=9)
            if column == 0:
                ax.set_ylabel(process.capitalize())
            if row == 2:
                ax.set_xlabel("Observation gap (years)")
    fig.legend(handles=direction_legend, loc="lower center", ncol=2, frameon=False)
    fig.suptitle("Regular sampling density: native scores on the common 20-year analysis grid")
    fig.tight_layout(rect=(0, 0.045, 1, 0.96))
    _save_sampling_figure(fig, output_dir, "03_density_native_scores")

    scenarios = [f"density_{gap}" for gap in [20, 40, 80, 160]]
    scenarios = scenarios + [f"irregular_{gap}_{pair}" for gap in ["cv025", "cv075"]
                             for pair in ["AA", "BB", "AB", "BA"]]
    scenarios = scenarios + [f"gradient_{gap}_{gradient}_{pair}" for gap in ["cv025", "cv075"]
                             for gradient in ["old_sparse", "young_sparse"] for pair in ["AA", "BB"]]
    columns = [(process, method) for process in processes for method in methods]
    matrix = np.full((len(scenarios), len(columns)), np.nan)
    cell_labels = {}
    for row, scenario in enumerate(scenarios):
        for column, (process, method) in enumerate(columns):
            match = summaries[(summaries["scenario"] == scenario) & (summaries["process"] == process)
                              & (summaries["method"] == method)]
            if match.empty:
                cell_labels[row, column] = "not reported"
                continue
            item = match.iloc[0]
            label = f"{int(item['n_correct_graph'])}/{int(item['n_requested'])}"
            for count, prefix in [("n_pending", "P"), ("n_failed", "F"), ("n_unresolved", "U")]:
                if item[count]:
                    label = label + f"\n{prefix}:{int(item[count])}"
            cell_labels[row, column] = label
            if item["n_pending"] == 0:
                matrix[row, column] = item["fraction_correct_requested"]
    fig, ax = plt.subplots(figsize=(15, 13))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("0.9")
    plotted = ax.imshow(np.ma.masked_invalid(matrix), vmin=0, vmax=1, cmap=cmap, aspect="auto")
    for (row, column), label in cell_labels.items():
        color = "white" if np.isfinite(matrix[row, column]) and matrix[row, column] < 0.55 else "black"
        ax.text(column, row, label, ha="center", va="center", fontsize=6.5, color=color)
    ax.set_xticks(range(len(columns)), labels=[f"{process}\n{method}" for process, method in columns], fontsize=8)
    ax.set_yticks(range(len(scenarios)), labels=[f"{index + 1:02d}  {name.replace('_', ' ')}"
                                               for index, name in enumerate(scenarios)], fontsize=8)
    ax.set_title("Complete true graph recovery / requested sampling cases")
    fig.colorbar(plotted, ax=ax, label="Fraction of requested cases", shrink=0.7)
    fig.text(0.5, 0.01, "P: pending; F: unusable decisions; U: unresolved evidence. Pending cells stay gray. "
             "Unresolved and unusable cases do not count as correct; these are not independent climate histories.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    _save_sampling_figure(fig, output_dir, "04_correct_graph_requested")

    categories = [(name, gap) for name in ["asynchrony", "old_sparse", "young_sparse"]
                  for gap in ["cv025", "cv075"]]
    contrast_labels = {"asynchrony": "Async.", "old_sparse": "Old\nsparse", "young_sparse": "Young\nsparse"}
    labels = [f"{contrast_labels[name]}\nCV {0.25 if gap == 'cv025' else 0.75:g}" for name, gap in categories]
    fig, axes = plt.subplots(3, 4, figsize=(17, 9), sharex=True, sharey="col")
    for row, process in enumerate(processes):
        for column, method in enumerate(methods):
            ax = axes[row, column]
            ax.axhline(0, color="0.6", lw=0.8)
            for position, (name, gap) in enumerate(categories):
                for source, offset in [("x", -0.13), ("y", 0.13)]:
                    group = contrasts[(contrasts["process"] == process) & (contrasts["method"] == method)
                                      & (contrasts["contrast"] == name) & (contrasts["gap_label"] == gap)
                                      & (contrasts["source"] == source)].sort_values("mask_group")
                    values = pd.to_numeric(group.loc[group["status"] == "ok", "delta_score"], errors="coerce")
                    values = values[np.isfinite(values)].to_numpy()
                    center = position + offset
                    if len(values):
                        ax.scatter(center + np.linspace(-0.045, 0.045, len(values)), values,
                                   s=13, color=colors[source], alpha=0.45)
                        q25, median, q75 = np.percentile(values, [25, 50, 75])
                        ax.vlines(center, q25, q75, color=colors[source], lw=3)
                        ax.plot(center, median, marker="D", color=colors[source], markersize=4)
                    if len(values) < len(group):
                        ax.text(center, 0.97, f"{len(values)}/{len(group)}", color=colors[source], fontsize=6,
                                rotation=90, ha="center", va="top", transform=ax.get_xaxis_transform())
            ax.set_xticks(range(len(categories)), labels=labels, fontsize=7)
            ax.grid(axis="y", alpha=0.2)
            if row == 0:
                ax.set_title(method)
            if column == 0:
                ax.set_ylabel(f"{process.capitalize()}\nNative-score difference")
    # Give constant differences a readable display range without changing their values.
    for column, method in enumerate(methods):
        values = contrasts.loc[(contrasts["method"] == method) & contrasts["status"].eq("ok"), "delta_score"]
        values = pd.to_numeric(values, errors="coerce")
        values = values[np.isfinite(values)].to_numpy()
        if len(values) and np.ptp(values) == 0:
            center = float(values[0])
            margin = 0.05 if center == 0 else abs(center) * 0.05
            axes[0, column].set_ylim(center - margin, center + margin)
    fig.legend(handles=direction_legend, loc="lower center", ncol=2, frameon=False)
    fig.suptitle("Paired sampling effects: points = mask groups; diamonds = median; bars = IQR")
    fig.text(0.5, 0.045, "Asynchrony: (AB + BA − AA − BB)/2. Gradients: synchronous gradient mean − synchronous no-gradient mean. "
             "IQR is not a confidence interval; labels show usable/requested groups when incomplete.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.075, 1, 0.96))
    _save_sampling_figure(fig, output_dir, "05_paired_native_score_contrasts")

    fig, axes = plt.subplots(3, 2, figsize=(16, 10), sharex=True, sharey=True)
    ccm = summaries[summaries["method"] == "CCM"]
    positions = np.arange(len(scenarios))
    for row, process in enumerate(processes):
        for column, direction in enumerate(["x_to_y", "y_to_x"]):
            ax = axes[row, column]
            group = ccm[ccm["process"] == process].set_index("scenario").reindex(scenarios)
            ax.bar(positions - 0.18, group[f"n_raw_retained_{direction}"] / group["n_requested"],
                   width=0.36, facecolor="white", edgecolor="tab:blue", label="Raw Tp=0")
            ax.bar(positions + 0.18, group[f"n_retained_{direction}"] / group["n_requested"],
                   width=0.36, color="tab:blue", label="After lag screen")
            for position in np.flatnonzero(group["n_pending"].fillna(0).to_numpy() > 0):
                ax.text(position, 1.02, "P", ha="center", fontsize=7, color="0.4")
            ax.set_ylim(0, 1.12)
            ax.set_title(f"{process.capitalize()}: {direction.replace('_to_', ' → ').upper()}")
            ax.set_xticks(positions, labels=[str(index + 1) for index in positions], fontsize=7)
            ax.grid(axis="y", alpha=0.2)
            if column == 0:
                ax.set_ylabel("Retained / requested cases")
            if row == 2:
                ax.set_xlabel("Scenario row in the recovery heatmap")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.965),
               ncol=2, frameon=False, fontsize=9)
    fig.suptitle("CCM raw and screened retention; only cases with usable paired directions count")
    fig.text(0.5, 0.01, "The denominator includes every requested case. P marks pending cases. "
             "Original rho0 and p values remain unchanged by the lag screen.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.035, 1, 0.93))
    _save_sampling_figure(fig, output_dir, "06_ccm_raw_screened_retention")



def plot_sampling_individual_scores(results, output_dir):
    """Show all sampling pairings separately; IQR describes masks, not uncertainty across histories."""
    if results["replicate"].nunique() != 1:
        raise ValueError("Plot one independent process realization at a time.")
    processes = ["linear", "nonlinear", "independent"]
    methods = ["PCMCI+", "LKIF", "CCM", "GC"]
    score_labels = ["Max. GPDC", "Max. |normalized flow|", "Original Tp=0 rho", "Log residual variance ratio"]
    categories = [(gap, pair) for gap in ["cv025", "cv075"] for pair in ["AA", "BB", "AB", "BA"]]
    labels = [f"{pair}\nCV {0.25 if gap == 'cv025' else 0.75:g}" for gap, pair in categories]
    colors = {"x": "tab:blue", "y": "tab:orange"}
    irregular = results[results["family"] == "irregularity"]
    baseline = results[(results["family"] == "density") & (results["scenario"] == "density_80")]
    fig, axes = plt.subplots(3, 4, figsize=(17, 9), sharex=True, sharey="col")
    for row, process in enumerate(processes):
        for column, method in enumerate(methods):
            ax = axes[row, column]
            for source, offset in [("x", -0.14), ("y", 0.14)]:
                reference = baseline[(baseline["process"] == process) & (baseline["method"] == method)
                                     & (baseline["source"] == source)]
                if len(reference) == 1 and reference["status"].iloc[0] == "ok":
                    value = pd.to_numeric(reference["score"].iloc[0], errors="coerce")
                    if pd.notna(value) and np.isfinite(value):
                        ax.axhline(value, color=colors[source], linestyle="--", linewidth=1, alpha=0.8)
                for position, (gap, pairing) in enumerate(categories):
                    group = irregular[(irregular["process"] == process) & (irregular["method"] == method)
                                      & (irregular["source"] == source) & (irregular["gap_label"] == gap)
                                      & (irregular["pairing"] == pairing)].sort_values("mask_group")
                    values = pd.to_numeric(group.loc[group["status"].eq("ok"), "score"], errors="coerce")
                    values = values[np.isfinite(values)].to_numpy()
                    center = position + offset
                    if len(values):
                        ax.scatter(center + np.linspace(-0.045, 0.045, len(values)), values,
                                   color=colors[source], s=13, alpha=0.45)
                        q25, median, q75 = np.percentile(values, [25, 50, 75])
                        ax.vlines(center, q25, q75, color=colors[source], linewidth=3)
                        ax.plot(center, median, color=colors[source], marker="D", markersize=4)
                    if len(values) < len(group):
                        ax.text(center, 0.97, f"{len(values)}/{len(group)}", color=colors[source], fontsize=6,
                                rotation=90, ha="center", va="top", transform=ax.get_xaxis_transform())
            ax.axvline(3.5, color="0.8", linewidth=0.8)
            ax.set_xticks(range(len(categories)), labels=labels, fontsize=7)
            ax.grid(axis="y", alpha=0.2)
            if row == 0:
                ax.set_title(f"{method}: {score_labels[column]}", fontsize=9)
            if column == 0:
                ax.set_ylabel(process.capitalize())
    for column, method in enumerate(methods):
        plotted = results[(results["method"] == method) & results["status"].eq("ok")
                          & (results["family"].eq("irregularity") | results["scenario"].eq("density_80"))]
        values = pd.to_numeric(plotted["score"], errors="coerce")
        values = values[np.isfinite(values)].to_numpy()
        if len(values) and np.ptp(values) == 0:
            center = float(values[0])
            margin = 0.05 if center == 0 else abs(center) * 0.05
            axes[0, column].set_ylim(center - margin, center + margin)
    handles = [Line2D([], [], color=colors[source], marker="o", label=label)
               for source, label in [("x", "X → Y"), ("y", "Y → X")]]
    handles.append(Line2D([], [], color="0.4", linestyle="--", label="Regular 80-year baseline"))
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False)
    fig.suptitle("Individual irregular sampling pairings: points = mask groups; diamonds = median; bars = IQR")
    fig.text(0.5, 0.045, "AA/BB use synchronous masks; AB/BA use different X/Y masks. "
             "Labels show usable/requested scores when pending or failed. Successful CCM rho0 survives lag-screen failure. "
             "IQR is not a confidence interval.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.075, 1, 0.96))
    _save_sampling_figure(fig, output_dir, "07_irregular_sampling_native_scores")


def save_summaries_and_figures(results):
    summaries, states, contrasts = summarize_sampling_results(results)
    write_table(summaries, DATA_DIR / "scenario_summary.csv")
    write_table(states, DATA_DIR / "case_states.csv")
    write_table(contrasts, DATA_DIR / "paired_score_contrasts.csv")

    # Differences are defined only when both the case and its original reference succeeded.
    keys = ["replicate", "process", "method", "source", "target"]
    reference = results.loc[results["scenario"].eq("density_20"), keys + ["score", "status"]]
    reference = reference.rename(columns={"score": "original_score", "status": "original_status"})
    changes = results.merge(reference, on=keys, how="left", validate="many_to_one")
    usable = (changes["status"].eq("ok") & changes["original_status"].eq("ok")
              & np.isfinite(changes["score"]) & np.isfinite(changes["original_score"]))
    changes["score_change_usable"] = usable
    changes["score_change_from_original"] = np.where(
        usable, changes["score"] - changes["original_score"], np.nan)
    write_table(changes, DATA_DIR / "score_changes_from_original.csv")

    contrast_summary = []
    keys = ["replicate", "process", "gap_label", "method", "source", "target", "contrast"]
    for key, group in contrasts.groupby(keys, sort=False):
        values = group.loc[group["status"].eq("ok"), "delta_score"].to_numpy(dtype=float)
        q25, median, q75 = np.percentile(values, [25, 50, 75]) if len(values) else [np.nan] * 3
        contrast_summary.append({**dict(zip(keys, key)), "n_requested": len(group),
                                 "n_usable": len(values), "n_pending": int(group.status.eq("pending").sum()),
                                 "q25": q25, "median": median, "q75": q75})
    write_table(pd.DataFrame(contrast_summary), DATA_DIR / "paired_contrast_summary.csv")
    plot_sampling_results(results, summaries, contrasts, FIGURE_DIR)
    plot_sampling_individual_scores(results, FIGURE_DIR)


def main():
    if N_WORKERS < 1 or N_MASK_GROUPS < 1:
        raise ValueError("Worker and sampling-group counts must be positive.")
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    parameters = execution_parameters()
    check_parameters(parameters)
    (DATA_DIR / "parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    started = time.perf_counter()
    completed, jobs = {}, []
    try:
        jobs, originals, clock_metadata, points, diagnostics = prepare_inputs()
        for job in jobs:
            result = read_completed_case(job)
            if result is not None:
                completed[job["case_id"]] = result
        plot_sampling_designs(originals, clock_metadata, points)
        if PREPARE_ONLY:
            write_status("prepared", len(completed), len(jobs), started)
            print(f"Prepared {len(jobs)} fixed sampling inputs; {len(completed)} complete.", flush=True)
            return
        pending = [job for job in jobs if job["case_id"] not in completed]
        if RUN_CAUSAL and pending:
            effective_n = len(originals[PROCESSES[0]]) - 2 * MAX_LAG
            prepare_gpdc_null(gpdc_null_path(), sample_sizes=[effective_n],
                              sig_samples=GPDC_NULL_SAMPLES, seed=METHOD_SEED)
            print(f"Sampling experiment: {len(completed)} complete, {len(pending)} pending; {N_WORKERS} workers.", flush=True)
            write_status("running", len(completed), len(jobs), started)
            # Complete stages in order; already-computed synchronous controls are shared.
            for family in ["density", "irregularity", "gradient"]:
                stage_jobs = [job for job in pending if job["family"] == family]
                if not stage_jobs:
                    continue
                context = multiprocessing.get_context("spawn")
                with ProcessPoolExecutor(max_workers=min(N_WORKERS, len(stage_jobs)), mp_context=context) as pool:
                    futures = [pool.submit(fit_case, job, parameters, str(DATA_DIR), str(gpdc_null_path()))
                               for job in stage_jobs]
                    for future in as_completed(futures):
                        result = future.result()
                        completed[result.case_id.iloc[0]] = result
                        if len(completed) % 10 == 0 or len(completed) == len(jobs):
                            write_status("running", len(completed), len(jobs), started)
                            print(f"Completed {len(completed)}/{len(jobs)} sampling pairs.", flush=True)
                write_status("running", len(completed), len(jobs), started)
                if len(completed) < len(jobs):
                    results = collect_results(jobs, completed)
                    write_table(results, DATA_DIR / "member_results.csv")
                    save_summaries_and_figures(results)
                print(f"Stage {family} calculations complete.", flush=True)
        results = collect_results(jobs, completed)
        write_table(results, DATA_DIR / "member_results.csv")
        save_summaries_and_figures(results)
        state = "complete" if len(completed) == len(jobs) else "partial"
        write_status(state, len(completed), len(jobs), started)
        print(f"Sampling experiment {state}: {len(completed)}/{len(jobs)} committed pairs.", flush=True)
    except BaseException as error:
        if jobs:
            write_table(collect_results(jobs, completed), DATA_DIR / "member_results.csv")
        write_status("interrupted_or_failed", len(completed), len(jobs), started, f"{type(error).__name__}: {error}")
        raise


if __name__ == "__main__":
    main()
