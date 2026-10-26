"""Check whether native method scores rank paired coupling strengths.

Use only original uniform records. K changes while innovations, memory,
record length and estimator settings remain fixed. Scores are not estimates
of K. Run with python 02_1_check_coupling_strength.py in the project environment.
"""

from pathlib import Path
from importlib.metadata import version
import json

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from toolbox.process import generate_process
from toolbox.causal import prepare_gpdc_null, run_pcmciplus, run_lkif, run_ccm, run_granger


RECOMPUTE = True               # False redraws the saved result table without rerunning estimators.
SEEDS = [20261009, 20261301, 20261302]  # Original example and two fresh process realizations.
GAINS = [1.0, 2.0, 4.0]        # Long-run Y response to a held unit g(X); main baseline uses K=4.
DURATION = 30.0                # Total retained duration, kyr.
DT = 0.020                     # Original uniform step, kyr; no sampling or interpolation.
TRUE_LAG = 0.020               # One-step X -> Y coupling; no reverse term.
MEMORY_X = 0.5                 # Driver decay time, kyr.
MEMORY_Y = 0.2                 # Response decay time, kyr.
BACKGROUND_X = 1 / np.sqrt(1 - 0.95**2)  # Driver stationary SD in arbitrary units.
BACKGROUND_Y = 1 / np.sqrt(1 - 0.8**2)   # Unforced Y stationary SD, fixed across K.
BURN_TIMESCALES = 5            # Warm-up in units of the longest memory time.
ANALYSIS_OLD = 29.0            # Old edge of the shared analysis window, kyr BP.
ANALYSIS_YOUNG = 1.0           # Young edge; all original samples inside are retained.
MAX_LAG = 5                    # Method history bound, in grid steps.
GC_ORDER = 1                   # Fixed VAR order, matching the current experiment 02 baseline.
ALPHA = 0.05                   # Direction threshold, identical to the two-baseline inspection.
GPDC_NULL_SAMPLES = 10000      # Native GPDC reference-null draws.
CCM_E = 3                     # Number of embedding coordinates.
CCM_TAU = -1                  # Backward embedding interval: 20 years; total history spans 40 years.
CCM_EXCLUSION_RADIUS = 10      # Exclude neighbors within this many time steps.
CCM_LIBRARY_FRACTION = 0.80    # Fraction of usable embedded rows in each library.
CCM_SAMPLES = 10               # Number of random libraries averaged.
CCM_SURROGATES = 100           # Cut-and-swap null records per tested direction.
CCM_LAG_STEPS = list(range(-15, 16))  # Target shifts from -300 to +300 years; zero remains in the scan.
CCM_LAG_LIBRARIES = 100        # Matched random libraries for observed lag curves; distinct from null records.
CCM_LAG_SIGN_THRESHOLD = 0.90  # Stable peak-sign fraction; not a calibrated causal probability.

ROOT = Path(__file__).resolve().parent
STEM = Path(__file__).stem
DATA_DIR = ROOT / "data" / "processed" / STEM
FIGURE_DIR = ROOT / "figures" / STEM
METHODS = ["PCMCI+", "LKIF", "CCM", "GC"]


def generate_pair(seed, form, gain):
    """Change coupling only; reuse the seed and all innovations within a triplet."""
    ar_x = np.exp(-DT / MEMORY_X)
    ar_y = np.exp(-DT / MEMORY_Y)
    noise = (BACKGROUND_X * np.sqrt(1 - ar_x**2), BACKGROUND_Y * np.sqrt(1 - ar_y**2))
    process = generate_process(
        n_steps=int(round(DURATION / DT)) + 1, dt=DT, age_start=DURATION,
        lag=TRUE_LAG, ar_x=ar_x, ar_y=ar_y, coupling=(1 - ar_y) * gain,
        coupling_form=form, noise_std=noise,
        burn_in=int(np.ceil(BURN_TIMESCALES * max(MEMORY_X, MEMORY_Y) / DT)), rng=seed,
    )
    first = int(round((DURATION - ANALYSIS_OLD) / DT))
    last = int(round((DURATION - ANALYSIS_YOUNG) / DT))
    series = process["series"].iloc[first:last + 1].copy().reset_index(drop=True)
    series["age"] = np.linspace(ANALYSIS_OLD, ANALYSIS_YOUNG, len(series))
    series["time"] = ANALYSIS_OLD - series["age"]
    return series


def run_experiment():
    """Retain scores and p values for both directions, including non-detections."""
    tables = []
    diagnostics = []
    for seed in SEEDS:
        n_grid = int(round((ANALYSIS_OLD - ANALYSIS_YOUNG) / DT)) + 1
        effective_size = n_grid - 2 * MAX_LAG
        cache_file = ROOT / "data" / "cache" / f"gpdc_T{effective_size}_n{GPDC_NULL_SAMPLES}_seed{seed + 300}.npz"
        null_file = prepare_gpdc_null(cache_file, sample_sizes=[effective_size],
                                     sig_samples=GPDC_NULL_SAMPLES, seed=seed + 300)
        for form in ["linear", "nonlinear"]:
            for gain in GAINS:
                series = generate_pair(seed, form, gain)
                library_size = int(CCM_LIBRARY_FRACTION * (len(series) - (CCM_E - 1) * abs(CCM_TAU)))
                for method in METHODS:
                    print(f"seed={seed}, {form}, K={gain:g}: {method}", flush=True)
                    label = f"s{seed}_{form}_K{gain:g}_{method.replace('+', 'plus')}"
                    try:
                        if method == "PCMCI+":
                            result = run_pcmciplus(
                                series, tau_max=MAX_LAG, pc_alpha=ALPHA, alpha=ALPHA,
                                sig_samples=GPDC_NULL_SAMPLES, seed=seed + 300,
                                gp_params={"optimizer": None}, null_dist_filename=null_file,
                            )
                        elif method == "LKIF":
                            result = run_lkif(series, max_lag=MAX_LAG, alpha=ALPHA)
                        elif method == "CCM":
                            result = run_ccm(
                                series, E=CCM_E, tau=CCM_TAU, exclusion_radius=CCM_EXCLUSION_RADIUS,
                                lib_size=library_size, sample=CCM_SAMPLES, n_surrogates=CCM_SURROGATES,
                                seed=seed + 300, alpha=ALPHA, lag_steps=CCM_LAG_STEPS,
                                n_libraries=CCM_LAG_LIBRARIES,
                                sign_threshold=CCM_LAG_SIGN_THRESHOLD,
                            )
                        else:
                            result = run_granger(series, max_lag=GC_ORDER, alpha=ALPHA)
                        edges = result["edges"].copy()
                        diagnostic = result["diagnostics"]
                        if "unoriented" not in edges:
                            edges["unoriented"] = False
                        if (diagnostic.get("undirected_contemporaneous", False)
                                or diagnostic.get("conflicting_contemporaneous", False)):
                            edges["unoriented"] = True
                        arrays = {}
                        for key, value in result["raw"].items():
                            if isinstance(value, np.ndarray):
                                arrays[key] = value.astype(str) if value.dtype == object else value
                            elif isinstance(value, pd.DataFrame):
                                value.to_csv(DATA_DIR / f"{label}_{key}.csv", index=False)
                        if arrays:
                            np.savez_compressed(DATA_DIR / f"{label}_raw.npz", **arrays)
                    except Exception as error:
                        edges = pd.DataFrame([
                            {"source": source, "target": target, "score": np.nan,
                             "score_name": "unavailable", "p_value": np.nan, "p_kind": "unavailable",
                             "edge_detected": None, "unoriented": False, "status": "error", "error": str(error)}
                            for source, target in [("x", "y"), ("y", "x")]
                        ])
                        diagnostic = {"status": "error", "error": str(error)}
                    edges["method"] = method
                    edges["seed"] = seed
                    edges["process"] = form
                    edges["gain"] = gain
                    tables.append(edges)
                    diagnostics.append({"process_seed": seed, "process": form, "gain": gain,
                                        "method": method, **diagnostic})
                    lag_message = f", lag={edges['lag_status'].tolist()}" if "lag_status" in edges else ""
                    print(f"  status={edges['status'].tolist()}{lag_message}", flush=True)
    results = pd.concat(tables, ignore_index=True)
    results.to_csv(DATA_DIR / "strength_results.csv", index=False)
    ccm_columns = ["seed", "process", "gain", "source", "target", "score", "p_value",
                   "raw_edge_detected", "edge_detected", "unoriented", "status", "lag_status",
                   "lag_reason", "optimal_lag_steps", "lag_past_fraction", "lag_future_fraction"]
    results[results["method"] == "CCM"].reindex(columns=ccm_columns).to_csv(DATA_DIR / "ccm_filter_comparison.csv", index=False)
    (DATA_DIR / "method_diagnostics.json").write_text(json.dumps(diagnostics, indent=2), encoding="utf-8")
    return results


def summarize_rankings(results):
    """Count complete increasing triplets without filtering on significance."""
    rows = []
    for (method, form, source), group in results.groupby(["method", "process", "source"]):
        valid_triplets = 0
        ordered_triplets = 0
        for _, triplet in group.groupby("seed"):
            triplet = triplet.sort_values("gain")
            complete = triplet["gain"].tolist() == sorted(GAINS)
            valid = complete and triplet["status"].eq("ok").all() and np.isfinite(triplet["score"]).all()
            if valid:
                valid_triplets = valid_triplets + 1
                ordered_triplets = ordered_triplets + int(np.all(np.diff(triplet["score"]) > 0))
        rows.append({"method": method, "process": form, "source": source,
                     "ordered_triplets": ordered_triplets, "valid_triplets": valid_triplets,
                     "expected_triplets": len(SEEDS)})
    return pd.DataFrame(rows)


def plot_strength_scores(results):
    """Show each paired record, with the same vertical scale within a method."""
    fig, axes = plt.subplots(2, 4, figsize=(14, 8), sharex=True, sharey="col")
    score_labels = ["Max. GPDC statistic", "Max. |normalized information flow|",
                    "Original Tp=0 skill (rho)", "Log residual variance ratio"]
    for column, method in enumerate(METHODS):
        for row, form in enumerate(["linear", "nonlinear"]):
            ax = axes[row, column]
            for source, color, style, marker in [("x", "tab:blue", "-", "o"),
                                                  ("y", "tab:orange", "--", "s")]:
                subset = results[(results["method"] == method) & (results["process"] == form)
                                 & (results["source"] == source)]
                for _, triplet in subset.groupby("seed"):
                    triplet = triplet.sort_values("gain")
                    ax.plot(triplet["gain"], triplet["score"], color=color, linestyle=style, alpha=0.45)
                    for _, edge in triplet.iterrows():
                        if edge["status"] != "ok" or not np.isfinite(edge["score"]):
                            ax.text(edge["gain"], 0.04, "Fit failed", transform=ax.get_xaxis_transform(),
                                    ha="center", rotation=90, fontsize=7, color="0.4")
                        elif method == "CCM" and edge.get("lag_status") == "failed":
                            ax.scatter(edge["gain"], edge["score"], marker="X", color="0.4", s=45, zorder=4)
                        elif edge.get("unoriented", False) == True:
                            ax.scatter(edge["gain"], edge["score"], marker="D", edgecolor=color,
                                       facecolor="#7FC97F", s=40, zorder=4)
                        elif method == "CCM" and edge.get("raw_edge_detected") == True and edge["edge_detected"] == False:
                            ax.scatter(edge["gain"], edge["score"], marker="x", color=color, s=40, zorder=4)
                        else:
                            ax.scatter(edge["gain"], edge["score"], marker=marker, s=30, edgecolor=color,
                                       facecolor=color if edge["edge_detected"] == True else "white", zorder=3)
            ax.set_title(f"{method} — {form.capitalize()}", fontsize=11)
            ax.set_xticks(GAINS)
            ax.grid(alpha=0.15)
            if row == 1:
                ax.set_xlabel("True model gain K")
            if row == 0:
                ax.set_ylabel(score_labels[column], fontsize=9)
    handles = [Line2D([], [], color="tab:blue", marker="o", label="True direction: X → Y"),
               Line2D([], [], color="tab:orange", marker="s", linestyle="--", label="Reverse: Y → X"),
               Line2D([], [], color="black", marker="o", ls="none", label="Retained direction"),
               Line2D([], [], color="black", marker="o", markerfacecolor="white", ls="none", label="Original test not significant"),
               Line2D([], [], color="black", marker="x", ls="none", label="CCM lag incompatible"),
               Line2D([], [], color="black", marker="D", markerfacecolor="#7FC97F", ls="none", label="Unresolved evidence"),
               Line2D([], [], color="0.4", marker="X", ls="none", label="CCM lag failed; score retained")]
    fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.035), frameon=False, fontsize=8)
    fig.suptitle("Native scores versus known coupling: original uniform records", fontsize=14)
    fig.text(0.5, 0.015, "Each line follows one paired realization. CCM keeps rho0 for every successful original fit, including unresolved or failed lag checks. Scores do not estimate K.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.13, 1, 0.95))
    return fig


def check_saved_results(results):
    """Reject earlier unscreened scores or changed scientific settings on redraw."""
    saved = json.loads((DATA_DIR / "parameters.json").read_text())
    names = ["SEEDS", "GAINS", "DURATION", "DT", "TRUE_LAG", "MEMORY_X", "MEMORY_Y",
             "BACKGROUND_X", "BACKGROUND_Y", "BURN_TIMESCALES", "ANALYSIS_OLD", "ANALYSIS_YOUNG",
             "MAX_LAG", "GC_ORDER", "ALPHA", "GPDC_NULL_SAMPLES", "CCM_E", "CCM_TAU",
             "CCM_EXCLUSION_RADIUS", "CCM_LIBRARY_FRACTION", "CCM_SAMPLES", "CCM_SURROGATES",
             "CCM_LAG_STEPS", "CCM_LAG_LIBRARIES", "CCM_LAG_SIGN_THRESHOLD"]
    for name in names:
        if saved.get(name) != globals()[name]:
            raise ValueError(f"Saved strength setting {name} differs; set RECOMPUTE=True.")
    for package, saved_version in saved["packages"].items():
        if version(package) != saved_version:
            raise ValueError(f"{package} changed; set RECOMPUTE=True.")
    ccm = results[(results["method"] == "CCM") & (results["status"] == "ok")]
    if len(ccm):
        required = ["raw_edge_detected", "lag_status", "unoriented"]
        if any(column not in ccm for column in required) or ccm[required].isna().any().any():
            raise ValueError("Saved successful CCM fits lack lag-screen results; set RECOMPUTE=True.")


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    if RECOMPUTE:
        results = run_experiment()
        parameters = {name: value for name, value in globals().items()
                      if name.isupper() and isinstance(value, (str, int, float, list))}
        parameters["packages"] = {name: version(name) for name in
                                  ["numpy", "scipy", "pandas", "tigramite", "pyEDM", "LK-Info-Flow", "statsmodels"]}
        parameters["interpretation"] = "Paired strength ranking in three records; not calibration of effects or error rates. GPDC kernel optimizer=None."
        parameters["CCM_lag_rule"] = "Preserve Tp=0 rho and swap p; independently screen directions with the same settings as experiment 02. Lag seed equals CCM seed."
        parameters["GC_order_rule"] = "Fixed GC_ORDER=1, aligned with experiment 02; replaces the earlier BIC configuration."
        (DATA_DIR / "parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    else:
        results = pd.read_csv(DATA_DIR / "strength_results.csv")
        check_saved_results(results)
    summary = summarize_rankings(results)
    summary.to_csv(DATA_DIR / "ranking_summary.csv", index=False)
    fig = plot_strength_scores(results)
    fig.savefig(FIGURE_DIR / "01_strength_ranking.png", dpi=180)
    fig.savefig(FIGURE_DIR / "01_strength_ranking.pdf")
    plt.close(fig)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
