"""Regular-grid pairwise causal analyses using the authors' Python packages.

Coordinates are kyr BP (age decreases) and elapsed kyr (time increases).
No function draws chronologies or interprets simulation truth. The explicit
GPDC cache preparation helper writes the package's reusable null distribution.
"""

from importlib.metadata import version
from pathlib import Path
import json
import platform
import warnings

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def _series_values(series):
    """Check the shared regular, forward-time input without changing its values."""
    required = ["time", "age", "x", "y"]
    if not isinstance(series, pd.DataFrame) or not set(required).issubset(series.columns):
        raise ValueError("series must be a DataFrame containing time, age, x, and y.")
    data = series[required].to_numpy(dtype=float)
    if len(data) < 3 or not np.all(np.isfinite(data)):
        raise ValueError("series must have at least three finite rows.")
    spacing = np.diff(data[:, 0])
    if np.any(spacing <= 0) or not np.allclose(spacing, spacing[0], rtol=1e-8, atol=1e-12):
        raise ValueError("time must increase on a regular grid.")
    if not np.allclose(np.diff(data[:, 1]), -spacing, rtol=1e-8, atol=1e-12):
        raise ValueError("age must decrease by the same intervals as time increases.")
    values = data[:, 2:]
    if np.any(np.ptp(values, axis=0) == 0):
        raise ValueError("Causal analysis requires two nonconstant series.")
    return values, float(spacing[0])


def _edge_row(method, source, target, score, score_name, detected, p_value, p_kind):
    valid = np.isfinite(score) and np.isfinite(p_value) and 0 <= p_value <= 1
    return {
        "method": method,
        "source": source,
        "target": target,
        "score": float(score),
        "score_name": score_name,
        "edge_detected": bool(detected) if valid else None,
        "p_value": float(p_value),
        "p_kind": p_kind,
        "status": "ok" if valid else "invalid_output",
    }


def resample_pair(age_x, value_x, age_y, value_y, *, age_grid):
    """Interpolate two ordered irregular records onto one fixed age grid.

    All ages are in kyr BP, strictly decreasing in existing sample order.
    The regular age_grid must be covered by both records; extrapolation,
    sorting, smoothing, and value imputation are not performed.
    """
    age_grid = np.asarray(age_grid, dtype=float)
    if age_grid.ndim != 1 or len(age_grid) < 2 or not np.all(np.isfinite(age_grid)):
        raise ValueError("age_grid must contain at least two finite ages.")
    spacing = -np.diff(age_grid)
    if np.any(spacing <= 0) or not np.allclose(spacing, spacing[0], rtol=1e-8, atol=1e-12):
        raise ValueError("age_grid must be regular and strictly decreasing.")

    time_grid = -age_grid
    interpolated = {}
    ranges = []
    max_gaps = {}
    for name, age, value in [("x", age_x, value_x), ("y", age_y, value_y)]:
        age = np.asarray(age, dtype=float)
        value = np.asarray(value, dtype=float)
        if age.ndim != 1 or len(age) < 2 or value.shape != age.shape:
            raise ValueError(f"{name}: age and value must be matching one-dimensional arrays.")
        if not np.all(np.isfinite(age)) or not np.all(np.isfinite(value)):
            raise ValueError(f"{name}: ages and values must be finite.")
        if np.any(np.diff(age) >= 0):
            raise ValueError(f"{name}: ages must strictly decrease in sample order.")
        if age_grid[0] > age[0] or age_grid[-1] < age[-1]:
            raise ValueError(f"{name}: record does not cover the fixed grid; extrapolation is forbidden.")

        # Negating BP ages preserves sample order and gives forward-time coordinates.
        sample_time = -age
        interpolated[name] = np.interp(time_grid, sample_time, value)
        right = np.searchsorted(sample_time, time_grid, side="left")
        exact = sample_time[right] == time_grid
        left = np.maximum(right - 1, 0)
        bracket_gaps = sample_time[right] - sample_time[left]
        bracket_gaps[exact] = 0.0
        max_gaps[name] = float(np.max(bracket_gaps))
        ranges.append((float(age[-1]), float(age[0])))

    series = pd.DataFrame({
        "time": time_grid - time_grid[0],
        "age": age_grid,
        "x": interpolated["x"],
        "y": interpolated["y"],
    })
    return {
        "series": series,
        "diagnostics": {
            "status": "ok",
            "n_grid": len(age_grid),
            "grid_step_kyr": float(spacing[0]),
            "overlap_age_range": (max(r[0] for r in ranges), min(r[1] for r in ranges)),
            "analysis_age_range": (float(age_grid[-1]), float(age_grid[0])),
            "max_interpolation_gap_kyr": max_gaps,
        },
    }


def _gpdc_null_settings(sample_sizes, sig_samples, seed):
    """Identify one fixed-length native null distribution and its random stream."""
    sizes = list(sample_sizes)
    if len(sizes) != 1 or not isinstance(sizes[0], (int, np.integer)) or sizes[0] < 3:
        raise ValueError("Use one integer effective sample size >=3 per GPDC cache file.")
    if not isinstance(sig_samples, (int, np.integer)) or sig_samples < 1:
        raise ValueError("sig_samples must be a positive integer.")
    if not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer.")
    return {
        "sample_sizes": [int(sizes[0])], "sig_samples": int(sig_samples), "seed": int(seed),
        "packages": {name: version(name) for name in ["tigramite", "numpy", "scipy", "dcor"]},
        "platform": platform.system(), "machine": platform.machine(),
    }


def _check_gpdc_null(filename, settings):
    """Reject a stale or incomplete cache instead of changing its p-value reference."""
    metadata = json.loads(filename.with_suffix(".json").read_text(encoding="utf-8"))
    if metadata != settings:
        raise ValueError("GPDC null cache settings differ; prepare a different cache file.")
    with np.load(filename, allow_pickle=False) as saved:
        if not np.array_equal(saved["T"], settings["sample_sizes"]):
            raise ValueError("GPDC cache has different effective sample sizes.")
        distribution = saved["exact_dist"]
        if distribution.shape != (1, settings["sig_samples"]) or not np.isfinite(distribution).all():
            raise ValueError("GPDC cache has an incomplete null distribution.")


def prepare_gpdc_null(filename, *, sample_sizes, sig_samples, seed):
    """Prepare one native GPDC null distribution before starting worker processes.

    Use a separate file for each effective sample size and seed. This keeps
    the random stream identical to an uncached fixed-length PCMCI+ call;
    generating several lengths sequentially would advance that stream.
    Files are small and reusable, but are not portable across environments.
    """
    from tigramite.independence_tests.gpdc import GPDC

    filename = Path(filename)
    if filename.suffix != ".npz":
        raise ValueError("The GPDC null cache filename must end in .npz.")
    settings = _gpdc_null_settings(sample_sizes, sig_samples, seed)
    metadata_path = filename.with_suffix(".json")
    if filename.exists() or metadata_path.exists():
        _check_gpdc_null(filename, settings)
        return str(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)
    test = GPDC(significance="analytic", sig_samples=int(sig_samples), seed=int(seed))
    temporary = filename.with_suffix(".partial.npz")
    test.generate_and_save_nulldists(settings["sample_sizes"], str(temporary))
    temporary.replace(filename)
    metadata_temporary = metadata_path.with_suffix(".partial.json")
    metadata_temporary.write_text(json.dumps(settings, indent=2), encoding="utf-8")
    metadata_temporary.replace(metadata_path)
    return str(filename)


def run_pcmciplus(
    series, *, tau_max=30, pc_alpha=0.05, alpha=0.05, sig_samples=1000, seed=1, gp_params=None,
    null_dist_filename=None,
):
    """Run Tigramite PCMCI+ / GPDC with its native reference null distribution.

    Each directional Holm family contains all cross-lags 0..tau_max. A vote
    additionally requires an arrow in the native graph. Zero-lag unoriented
    and conflicting links are retained separately, never made directional.
    If either directional family retains such a zero-lag symbol, preserve it
    as unresolved evidence; this rule does not depend on the column order.
    gp_params is passed to Tigramite's native GP regression. With None, its
    RBF+white-noise kernel is optimized. {'optimizer': None} fixes the native
    unit length scale and noise level; the conditional regression still fits
    the data. This is an explicit model setting, not equivalent estimation.
    An optional prepare_gpdc_null() file reuses only the native reference
    distribution; every member's conditional regressions are fitted again.
    """
    from tigramite import data_processing
    from tigramite.independence_tests.gpdc import GPDC
    from tigramite.pcmci import PCMCI

    values, dt = _series_values(series)
    if not isinstance(tau_max, (int, np.integer)) or tau_max < 1 or 2 * tau_max >= len(values):
        raise ValueError("tau_max must be a positive integer shorter than half the record.")
    if not 0 < alpha < 1 or not 0 < pc_alpha < 1 or sig_samples < 1:
        raise ValueError("Require alpha and pc_alpha in (0, 1), and positive sig_samples.")

    # GPDC performs its own standardization, Gaussian regression, and null sampling.
    if null_dist_filename is not None:
        # All input rows are finite; Tigramite's default cut_off='2xtau_max'
        # leaves the same effective T for every conditional test in this run.
        settings = _gpdc_null_settings([len(values) - 2 * tau_max], sig_samples, seed)
        _check_gpdc_null(Path(null_dist_filename), settings)
    test = GPDC(significance="analytic", sig_samples=sig_samples, seed=seed, gp_params=gp_params,
                null_dist_filename=None if null_dist_filename is None else str(null_dist_filename))
    dataframe = data_processing.DataFrame(values, var_names=["x", "y"])
    pcmci = PCMCI(dataframe=dataframe, cond_ind_test=test, verbosity=0)
    with warnings.catch_warnings(record=True) as fit_warnings:
        warnings.simplefilter("always")
        result = pcmci.run_pcmciplus(
            tau_min=0, tau_max=tau_max, pc_alpha=pc_alpha, fdr_method="none",
        )
    p_matrix = np.asarray(result["p_matrix"], dtype=float)
    val_matrix = np.asarray(result["val_matrix"], dtype=float)
    graph = np.asarray(result["graph"])
    p_holm = np.full(p_matrix.shape, np.nan)
    graph_holm = graph.copy()
    rows = []
    for source, target in [(0, 1), (1, 0)]:
        p = p_matrix[source, target, :]
        scores = val_matrix[source, target, :]
        valid = np.all(np.isfinite(p)) and np.all((p >= 0) & (p <= 1))
        valid = valid and np.all(np.isfinite(scores))
        if valid:
            adjusted = multipletests(p, alpha=alpha, method="holm")[1]
            p_holm[source, target, :] = adjusted
            arrows = graph[source, target, :] == "-->"
            detected = np.any(arrows & (adjusted < alpha))
            graph_holm[source, target, adjusted >= alpha] = ""
            p_summary = float(np.min(adjusted))
            score = float(np.max(scores))
        else:
            detected, p_summary, score = False, np.nan, np.nan
            graph_holm[source, target, :] = ""
        rows.append(_edge_row(
            "pcmci+", ["x", "y"][source], ["x", "y"][target], score,
            "max_gpdc", detected, p_summary, "gpdc_reference_holm_within_direction",
        ))

    # Directional Holm families can retain different halves of a zero-lag link.
    zero_links = [str(graph_holm[0, 1, 0]), str(graph_holm[1, 0, 0])]
    diagnostics = {
        "package_version": version("tigramite"),
        "tau_max": tau_max,
        "lag_steps": np.arange(tau_max + 1).tolist(),
        "lag_kyr": (np.arange(tau_max + 1) * dt).tolist(),
        "pc_alpha": pc_alpha,
        "alpha": alpha,
        "sig_samples": sig_samples,
        "seed": seed,
        "gp_params": gp_params,
        "null_distribution": "generated" if null_dist_filename is None else "precomputed",
        "undirected_contemporaneous": "o-o" in zero_links,
        "conflicting_contemporaneous": "x-x" in zero_links,
        "native_contemporaneous_link": str(graph[0, 1, 0]),
        "n_fit_warnings": len(fit_warnings),
        "fit_warnings": sorted(set(str(w.message) for w in fit_warnings)),
    }
    return {
        "edges": pd.DataFrame(rows), "diagnostics": diagnostics,
        "raw": {**result, "p_holm": p_holm, "graph_holm": graph_holm},
    }


def run_lkif(series, *, max_lag=30, np_step=1, alpha=0.05):
    """Call LK_Info_Flow's multivariate estimator and native analytic p values.

    Information-flow arrays use [target, source, lag index]. The package's
    native lag indices are retained; they are not reported as propagation
    delays. IF has units kyr^-1; normalized IF is dimensionless.
    """
    from LK_Info_Flow.multi_causality_est import multi_causality_est

    values, dt = _series_values(series)
    if not isinstance(max_lag, (int, np.integer)) or max_lag < 1:
        raise ValueError("max_lag must be a positive integer.")
    if len(values) - max_lag <= 2 * max_lag + 1:
        raise ValueError("Too few rows for the requested LKIF lag model.")
    if np_step != 1 or not 0 < alpha < 1:
        raise ValueError("This workflow uses the author's recommended np_step=1 and alpha in (0, 1).")

    # Empty temporal order invokes the author's contiguous-integer-sample branch;
    # supplying floating kyr timestamps can fail its exact unit-spacing check.
    result = multi_causality_est(
        values.T, max_lag=max_lag, np=np_step, dt=dt,
        series_temporal_order=np.array([]), significance_test=1,
    )
    p_matrix = np.asarray(result["p"], dtype=float)
    normalized_flow = np.asarray(result["nIF"], dtype=float)
    if p_matrix.shape != (2, 2, max_lag) or normalized_flow.shape != p_matrix.shape:
        raise ValueError("Unexpected LK_Info_Flow output shape; direction mapping needs review.")
    p_holm = np.full(p_matrix.shape, np.nan)
    signed_max = {}
    rows = []
    for source, target in [(0, 1), (1, 0)]:
        p = p_matrix[target, source, :]
        flow = normalized_flow[target, source, :]
        valid = np.all(np.isfinite(p)) and np.all((p >= 0) & (p <= 1))
        valid = valid and np.all(np.isfinite(flow))
        direction = f"{['x', 'y'][source]}_to_{['x', 'y'][target]}"
        if valid:
            adjusted = multipletests(p, alpha=alpha, method="holm")[1]
            p_holm[target, source, :] = adjusted
            detected = np.any(adjusted < alpha)
            p_summary = float(np.min(adjusted))
            signed_max[direction] = float(flow[np.argmax(np.abs(flow))])
            score = abs(signed_max[direction])
        else:
            detected, p_summary, score = False, np.nan, np.nan
            signed_max[direction] = np.nan
        rows.append(_edge_row(
            "lkif", ["x", "y"][source], ["x", "y"][target], score,
            "max_abs_normalized_information_flow", detected, p_summary,
            "native_normal_approximation_holm_within_direction",
        ))
    return {
        "edges": pd.DataFrame(rows),
        "diagnostics": {
            "package_version": version("LK_Info_Flow"), "max_lag": max_lag,
            "native_lag_indices": np.arange(max_lag).tolist(), "dt_kyr": dt,
            "np_step": np_step, "alpha": alpha, "signed_max_normalized_flow": signed_max,
            "temporal_order": "contiguous_integer_samples", "n_fitted_samples": len(values) - max_lag,
        },
        "raw": {**result, "p_holm": p_holm},
    }


def run_ccm(
    series, *, E, tau, exclusion_radius, lib_size, sample, n_surrogates, seed, alpha=0.05,
    lag_steps=None, n_libraries=100, sign_threshold=0.90, lag_seed=None,
):
    """Compute fixed-Tp CCM and directional cut-and-swap empirical p values.

    For X->Y, the Y embedding reconstructs the original X target. Each null
    moves Y at a uniform random cut in 1..N-1; X stays fixed (and vice versa).
    The original test does not compare directions or search lags. When
    lag_steps is provided, apply the separate temporal-consistency screen
    after this test, preserving its Tp=0 score and p value.
    """
    import pyEDM

    values, dt = _series_values(series)
    integer_settings = [E, tau, exclusion_radius, lib_size, sample, n_surrogates, seed]
    if not all(isinstance(v, (int, np.integer)) for v in integer_settings):
        raise ValueError("CCM dimensions, offsets, library settings, counts, and seed must be integers.")
    if E < 1 or tau >= 0 or exclusion_radius < 0 or sample < 1 or n_surrogates < 1 or seed < 1:
        raise ValueError("Require E>=1, tau<0, nonnegative exclusion radius, positive counts and seed.")
    n_embedded = len(values) - (E - 1) * abs(tau)
    if not E + 2 <= lib_size <= n_embedded or not 0 < alpha < 1:
        raise ValueError("Library must fit the embedded record and contain at least E+2 points.")

    # Independent RNG streams keep the swap cuts separate from random libraries.
    cut_sequence, library_sequence = np.random.SeedSequence(seed).spawn(2)
    cut_rng = np.random.default_rng(cut_sequence)
    library_rng = np.random.default_rng(library_sequence)
    cuts = cut_rng.integers(1, len(values), size=(n_surrogates, 2))
    library_seeds = library_rng.integers(1, 2**31 - 1, size=(n_surrogates + 1, 2))
    null_rho = np.full((n_surrogates, 2), np.nan)
    observed_results = {}
    rows = []
    for direction, (source, target) in enumerate([("x", "y"), ("y", "x")]):
        observed = pyEDM.CCM(
            dataFrame=series[["time", "x", "y"]], columns=target, target=source,
            E=E, tau=tau, Tp=0, exclusionRadius=exclusion_radius,
            libSizes=[lib_size], sample=sample, seed=int(library_seeds[0, direction]),
            parallel=False, showPlot=False,
        )
        column = f"{target}:{source}"
        rho = float(observed[column].iloc[0])
        observed_results[f"{source}_to_{target}"] = observed
        for surrogate in range(n_surrogates):
            shifted = series[["time", "x", "y"]].copy()
            embedding_values = series[target].to_numpy()
            cut = cuts[surrogate, direction]
            shifted[target] = np.concatenate((embedding_values[cut:], embedding_values[:cut]))
            comparison = pyEDM.CCM(
                dataFrame=shifted, columns=target, target=source,
                E=E, tau=tau, Tp=0, exclusionRadius=exclusion_radius,
                libSizes=[lib_size], sample=sample,
                seed=int(library_seeds[surrogate + 1, direction]),
                parallel=False, showPlot=False,
            )
            null_rho[surrogate, direction] = float(comparison[column].iloc[0])
        complete = np.all(np.isfinite(null_rho[:, direction])) and np.isfinite(rho)
        if complete:
            p_value = (1 + np.count_nonzero(null_rho[:, direction] >= rho)) / (n_surrogates + 1)
        else:
            p_value = np.nan
        rows.append(_edge_row(
            "ccm", source, target, rho, "cross_map_rho", p_value < alpha,
            p_value, "embedding_swap_upper_tail",
        ))

    directions = ["x_to_y", "y_to_x"]
    result = {
        "edges": pd.DataFrame(rows),
        "diagnostics": {
            "package_version": version("pyEDM"), "E": E, "tau": tau, "Tp": 0,
            "exclusion_radius": exclusion_radius, "lib_size": lib_size,
            "sample": sample, "n_surrogates": n_surrogates, "seed": seed,
            "alpha": alpha, "dt_kyr": dt, "n_embedded": n_embedded,
            "swap_sequence": "embedding", "cut_range": (1, len(values) - 1),
        },
        "raw": {
            "observed": observed_results,
            "rho_null": null_rho,
            "directions": directions,
            "cuts": cuts,
            "library_seeds": library_seeds,
        },
    }
    if lag_steps is not None:
        result = filter_ccm_by_lag(
            series, result, lag_steps=lag_steps, n_libraries=n_libraries,
            sign_threshold=sign_threshold, seed=seed if lag_seed is None else lag_seed,
        )
    return result


def _ccm_lag_grid(lag_steps):
    """Require an ordered integer scan containing past, zero, and future."""
    lags = np.asarray(lag_steps)
    if lags.ndim != 1 or len(lags) < 3 or not np.issubdtype(lags.dtype, np.integer):
        raise ValueError("lag_steps must be a one-dimensional integer grid.")
    if np.any(np.diff(lags) <= 0) or lags[0] >= 0 or lags[-1] <= 0 or 0 not in lags:
        raise ValueError("lag_steps must increase strictly and include negative, zero, and positive lags.")
    return lags.astype(int)


def _ccm_lag_peak(rho, lag_steps):
    """Classify a peak using only floating-point tolerance for exact ties."""
    if not np.all(np.isfinite(rho)):
        return np.nan, np.nan, "failed"
    peak = float(np.max(rho))
    tolerance = 32 * np.finfo(float).eps * max(1.0, abs(peak))
    positions = np.flatnonzero(np.abs(rho - peak) <= tolerance)
    if len(positions) != 1:
        return np.nan, peak, "tied"
    position = int(positions[0])
    lag = int(lag_steps[position])
    if position in (0, len(lag_steps) - 1):
        return lag, peak, "boundary"
    if lag == 0:
        return lag, peak, "zero"
    return lag, peak, "past" if lag < 0 else "future"


def summarize_ccm_lag(
    ccm_result, *, n_libraries=None, lag_steps=None, sign_threshold=0.90,
):
    """Reclassify stored lag curves without fitting CCM again.

    Use a prefix of the random libraries and/or a subset of the saved lag
    grid. Subsetting retains the original common anchors and library length;
    it is not a refit with the extra data available in a narrower window.
    Returned raw lag arrays contain the selected subset. Keep the original
    result to compare other prefixes or windows. Fractions describe random
    library stability conditional on this record, not causal probabilities.
    """
    raw = dict(ccm_result["raw"])
    saved_lags = _ccm_lag_grid(raw["lag_steps"])
    rho = np.asarray(raw["lag_rho"], dtype=float)
    seeds = np.asarray(raw["lag_library_seeds"])
    directions = ["x_to_y", "y_to_x"]
    if rho.ndim != 3 or rho.shape[1:] != (len(saved_lags), 2) or len(rho) < 1:
        raise ValueError("lag_rho must have shape (n_libraries, n_lags, 2).")
    if list(raw["lag_directions"]) != directions or seeds.shape != (len(rho),):
        raise ValueError("Stored lag directions or random-library seeds are inconsistent.")
    if n_libraries is None:
        n_libraries = len(rho)
    if not isinstance(n_libraries, (int, np.integer)) or not 1 <= n_libraries <= len(rho):
        raise ValueError("n_libraries must select a nonempty prefix of the saved libraries.")
    if not 0.5 < sign_threshold <= 1:
        raise ValueError("sign_threshold must be in (0.5, 1].")
    lags = saved_lags if lag_steps is None else _ccm_lag_grid(lag_steps)
    if not np.all(np.isin(lags, saved_lags)):
        raise ValueError("The requested lag grid must be contained in the saved scan.")
    rho = rho[:n_libraries, np.searchsorted(saved_lags, lags), :].copy()
    edges = ccm_result["edges"].copy(deep=True).reset_index(drop=True)
    pairs = list(zip(edges["source"], edges["target"]))
    if len(pairs) != 2 or set(pairs) != {("x", "y"), ("y", "x")}:
        raise ValueError("CCM edges must contain exactly x->y and y->x.")
    if "raw_edge_detected" not in edges:
        edges["raw_edge_detected"] = edges["edge_detected"].copy()
    edges["edge_detected"] = edges["edge_detected"].astype(object)
    diagnostics = dict(ccm_result["diagnostics"])
    lag_diagnostics = dict(diagnostics["lag_filter"])
    dt = float(lag_diagnostics["dt_kyr"])
    peaks, profiles, n_valid = [], [], {}

    for direction, (source, target) in enumerate([("x", "y"), ("y", "x")]):
        curves = rho[:, :, direction]
        library_peaks = [_ccm_lag_peak(curve, lags) for curve in curves]
        peak_kinds = np.array([peak[2] for peak in library_peaks])
        valid_count = int(np.count_nonzero(peak_kinds != "failed"))
        n_valid[directions[direction]] = valid_count
        for library, (lag, peak, kind) in enumerate(library_peaks):
            peaks.append({
                "library": library, "library_seed": int(seeds[library]),
                "source": source, "target": target, "optimal_lag_steps": lag,
                "optimal_lag_kyr": lag * dt, "rho_peak": peak, "peak_kind": kind,
            })
        # A failed library is never omitted from a mean or fraction denominator.
        mean_curve = np.mean(curves, axis=0)
        spread = np.std(curves, axis=0, ddof=1) if n_libraries > 1 else np.full(len(lags), np.nan)
        mean_lag, mean_peak, mean_kind = _ccm_lag_peak(mean_curve, lags)
        fractions = {
            kind: float(np.count_nonzero(peak_kinds == kind) / n_libraries)
            for kind in ["past", "future", "boundary", "zero", "tied", "failed"]
        }
        if valid_count != n_libraries:
            lag_status, reason = "failed", "incomplete_library_curves"
        elif mean_kind == "past" and fractions["past"] >= sign_threshold:
            lag_status, reason = "compatible", "stable_past_peak"
        elif mean_kind == "future" and fractions["future"] >= sign_threshold:
            lag_status, reason = "incompatible", "stable_future_peak"
        else:
            lag_status = "unresolved"
            reason = f"mean_{mean_kind}_peak" if mean_kind not in ("past", "future") else "unstable_peak_sign"

        row = pairs.index((source, target))
        raw_vote = edges.at[row, "raw_edge_detected"]
        original_valid = edges.at[row, "status"] == "ok" and pd.notna(raw_vote)
        if not original_valid:
            lag_status, reason = "failed", "original_ccm_invalid"
        detected = bool(raw_vote) and lag_status == "compatible" if original_valid else None
        if lag_status == "failed":
            detected = None
        edges.at[row, "edge_detected"] = detected
        edges.at[row, "unoriented"] = bool(original_valid and raw_vote and lag_status == "unresolved")
        edges.at[row, "lag_status"] = lag_status
        edges.at[row, "lag_reason"] = reason
        edges.at[row, "optimal_lag_steps"] = mean_lag
        edges.at[row, "optimal_lag_kyr"] = mean_lag * dt
        edges.at[row, "rho_peak"] = mean_peak
        edges.at[row, "lag_n_valid_libraries"] = valid_count
        edges.at[row, "lag_n_requested_libraries"] = int(n_libraries)
        for kind, fraction in fractions.items():
            edges.at[row, f"lag_{kind}_fraction"] = fraction
        for position, lag in enumerate(lags):
            profiles.append({
                "source": source, "target": target, "lag_steps": int(lag),
                "lag_kyr": float(lag * dt), "rho_mean": float(mean_curve[position]),
                "rho_std": float(spread[position]),
            })

    lag_diagnostics.update({
        "n_libraries": int(n_libraries), "sign_threshold": float(sign_threshold),
        "lag_steps": lags.tolist(), "lag_kyr": (lags * dt).tolist(),
        "n_valid_libraries": n_valid, "fraction_denominator": int(n_libraries),
        "all_directions_complete": all(count == n_libraries for count in n_valid.values()),
    })
    edges["unoriented"] = edges["unoriented"].astype(bool)
    for column in ["lag_n_valid_libraries", "lag_n_requested_libraries"]:
        edges[column] = edges[column].astype(int)
    diagnostics["lag_filter"] = lag_diagnostics
    raw.update({
        "lag_rho": rho, "lag_steps": lags.copy(),
        "lag_library_seeds": seeds[:n_libraries].copy(),
        "lag_profile": pd.DataFrame(profiles), "lag_peaks": pd.DataFrame(peaks),
    })
    return {"edges": edges, "diagnostics": diagnostics, "raw": raw}


def filter_ccm_by_lag(
    series, ccm_result, *, lag_steps, n_libraries=100, sign_threshold=0.90, seed=1,
):
    """Screen both original CCM arrows for temporal consistency with pyEDM.

    ccm_result must belong to this same series; length and grid are checked,
    but the original result does not store the input values for comparison.
    For X->Y, Y's historical embedding reconstructs X(t+Tp); negative Tp
    refers to past X in forward time, even though BP ages decrease. A unique
    interior negative peak must occur in the mean curve and the requested
    fraction of random libraries. Positive peaks reject the original arrow;
    zero, tied, boundary, and unstable peaks leave its direction unresolved.
    This screen cannot create significance or replace the original Tp=0 p.

    The library subsampling follows native pyEDM (without replacement within
    each library). Each repeat uses one seed at every Tp and fixed prediction
    and embedding anchors. In pyEDM 2.5.7, both maps also use the same draws.
    It is inspired by Ye et al. (2015), doi:10.1038/srep14750, and Van Nes et
    al. (2015), doi:10.1038/nclimate2568, not their complete model-null test.
    """
    import pyEDM

    values, dt = _series_values(series)
    lags = _ccm_lag_grid(lag_steps)
    if not isinstance(n_libraries, (int, np.integer)) or n_libraries < 1:
        raise ValueError("n_libraries must be a positive integer.")
    if not isinstance(seed, (int, np.integer)) or seed < 1 or not 0.5 < sign_threshold <= 1:
        raise ValueError("Require a positive integer seed and sign_threshold in (0.5, 1].")
    diagnostics = dict(ccm_result["diagnostics"])
    E, tau = diagnostics["E"], diagnostics["tau"]
    library_size, exclusion_radius = diagnostics["lib_size"], diagnostics["exclusion_radius"]
    if E < 1 or tau >= 0 or not np.isclose(dt, diagnostics["dt_kyr"], rtol=1e-8, atol=1e-12):
        raise ValueError("The original CCM configuration and current series are inconsistent.")
    history = (E - 1) * abs(tau)
    if diagnostics["n_embedded"] != len(values) - history:
        raise ValueError("The current record length differs from the original CCM input.")
    anchor_start = max(history, -int(lags[0]))
    anchor_stop = len(values) - 1 - int(lags[-1])
    anchors = np.arange(anchor_start, anchor_stop + 1, dtype=int)
    valid_lib = np.zeros(len(values), dtype=bool)
    valid_lib[anchors] = True
    library_seeds = np.random.default_rng(seed).integers(1, 2**31 - 1, size=n_libraries)
    rho = np.full((n_libraries, len(lags), 2), np.nan)
    errors = []
    if not E + 2 <= library_size <= len(anchors):
        errors.append({"reason": "library_does_not_fit_common_anchors", "n_common_anchors": len(anchors)})
    else:
        data = series[["time", "x", "y"]]
        # A fixed lag has fixed embeddings and target alignment. Reuse the
        # native object while keeping every original library seed and curve.
        # Keep process-level parallelism outside this function: pyEDM's
        # sequential Project uses module-level worker arrays, not thread-local state.
        for position, lag in enumerate(lags):
            model = None
            for library, library_seed in enumerate(library_seeds):
                try:
                    if model is None:
                        model = pyEDM.CCM(
                            dataFrame=data, columns="y", target="x", E=E, tau=tau,
                            Tp=int(lag), exclusionRadius=exclusion_radius,
                            libSizes=[library_size], sample=1, seed=int(library_seed),
                            validLib=valid_lib, parallel=False, includeData=False, showPlot=False,
                            returnObject=True,
                        )
                    else:
                        model.Project(sample=1, seed=int(library_seed), parallel=False, includeData=False)
                    native = model.libMeans
                    if int(native["LibSize"].iloc[0]) != library_size:
                        raise ValueError("Native CCM changed the requested library size.")
                    rho[library, position, :] = [native["y:x"].iloc[0], native["x:y"].iloc[0]]
                except Exception as error:
                    # Preserve successful Tp=0 fits even if this extra scan fails.
                    model = None
                    errors.append({
                        "library": library, "lag_steps": int(lag),
                        "reason": f"{type(error).__name__}: {error}",
                    })
        errors.sort(key=lambda item: (item["library"], item["lag_steps"]))
    diagnostics["lag_filter"] = {
        "package_version": version("pyEDM"), "seed": int(seed), "dt_kyr": dt,
        "n_libraries": int(n_libraries), "lib_size": int(library_size),
        "E": int(E), "tau": int(tau), "exclusion_radius": int(exclusion_radius),
        "embedding_span_steps": int(history), "embedding_span_kyr": float(history * dt),
        "common_anchor_start": int(anchor_start), "common_anchor_stop": int(anchor_stop),
        "n_common_anchors": len(anchors), "anchor_scan_lag_steps": lags.tolist(),
        "native_sample": 1, "sampling": "without_replacement_within_library",
        "scan_errors": errors,
    }
    raw = dict(ccm_result["raw"])
    raw.update({
        "lag_rho": rho, "lag_steps": lags.copy(), "lag_library_seeds": library_seeds,
        "lag_directions": np.array(["x_to_y", "y_to_x"]),
        "lag_common_anchor_indices": anchors,
    })
    scanned = {"edges": ccm_result["edges"], "diagnostics": diagnostics, "raw": raw}
    return summarize_ccm_lag(scanned, sign_threshold=sign_threshold)


def run_granger(series, *, max_lag=30, alpha=0.05, order_selection=None):
    """Classical bivariate GC with fixed order or BIC selection in 1..max_lag.

    statsmodels returns both regressions on the same response rows. The
    score log(SSR_restricted/SSR_full) uses their common variance denominator,
    not the different degrees-of-freedom-corrected residual variances.

    None uses exactly max_lag, preserving the fixed-order comparison. 'bic'
    uses the package's BIC on common response rows to select one positive
    order for both directions, then refits at that order. Order zero is
    excluded because this analysis tests predictive effects of past values.
    The selected-order F p value is approximate: it does not account for
    order selection. Residual whiteness is reported, not used to hide edges.
    """
    from statsmodels.tsa.api import VAR
    from statsmodels.tsa.stattools import grangercausalitytests

    values, dt = _series_values(series)
    if not isinstance(max_lag, (int, np.integer)) or max_lag < 1 or not 0 < alpha < 1:
        raise ValueError("Require a positive integer max_lag and alpha in (0, 1).")
    if order_selection not in (None, "bic"):
        raise ValueError("order_selection must be None (fixed order) or 'bic'.")
    order = int(max_lag)
    bic_by_order = None
    if order_selection == "bic":
        selection = VAR(values).select_order(maxlags=max_lag, trend="c")
        bic_by_order = [float(value) for value in selection.ics["bic"][1:]]
        order = int(np.argmin(bic_by_order)) + 1

    # Diagnose the joint unrestricted model, using the same order as the tests.
    var_fit = VAR(values).fit(order, trend="c")
    whiteness_lags = min(max(10, order + 1), var_fit.nobs - 1)
    try:
        whiteness = var_fit.test_whiteness(nlags=whiteness_lags, adjusted=True)
        whiteness_p = float(whiteness.pvalue)
        whiteness_status = "ok"
    except np.linalg.LinAlgError:
        # Shared innovations can prevent this joint diagnostic, while the
        # two directional regressions remain estimable. Keep their results.
        whiteness_p = None
        whiteness_status = "singular_residual_covariance"
    rows, raw, fit_diagnostics = [], {}, {}
    for source, target in [(0, 1), (1, 0)]:
        # The second input column causes the first in this statsmodels interface.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="verbose is deprecated", category=FutureWarning)
            result = grangercausalitytests(values[:, [target, source]], maxlag=[order], verbose=False)
        tests, (restricted, full, contrast) = result[order]
        f_statistic, p_value, df_denom, df_num = tests["ssr_ftest"]
        if restricted.nobs != full.nobs:
            raise ValueError("GC restricted and full fits use different response samples.")
        if full.model.rank < len(full.params) or restricted.model.rank < len(restricted.params):
            raise ValueError("GC regression is rank deficient; its non-detection is not a valid null result.")
        score = np.log(restricted.ssr / full.ssr)
        source_name, target_name = ["x", "y"][source], ["x", "y"][target]
        direction = f"{source_name}_to_{target_name}"
        rows.append(_edge_row(
            "granger", source_name, target_name, score, "log_residual_variance_ratio",
            p_value < alpha, p_value,
            "native_joint_F_after_BIC_selection" if order_selection == "bic" else "native_joint_F",
        ))
        rows[-1]["selected_order"] = order
        raw[direction] = {"tests": tests, "restricted": restricted, "full": full, "contrast": contrast}
        fit_diagnostics[direction] = {
            "f_statistic": float(f_statistic), "df_num": float(df_num),
            "df_denom": float(df_denom), "nobs": int(full.nobs),
            "ssr_restricted": float(restricted.ssr), "ssr_full": float(full.ssr),
            "condition_number": float(full.condition_number),
            "design_rank": int(full.model.rank), "n_parameters": len(full.params),
        }
    return {
        "edges": pd.DataFrame(rows),
        "diagnostics": {
            "package_version": version("statsmodels"), "max_lag": max_lag,
            "order_selection": order_selection, "selected_order": order,
            "bic_by_positive_order": bic_by_order,
            "history_kyr": order * dt, "alpha": alpha, "fits": fit_diagnostics,
            "var_is_stable": bool(var_fit.is_stable()),
            "residual_whiteness_p": whiteness_p,
            "residual_whiteness_status": whiteness_status,
            "residual_whiteness_lags": whiteness_lags,
            "residual_whiteness_adjusted": True,
        },
        "raw": {**raw, "var_coefficients": var_fit.coefs,
                "var_intercept": var_fit.intercept, "var_residuals": var_fit.resid},
    }
