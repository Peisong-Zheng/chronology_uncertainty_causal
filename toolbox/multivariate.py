"""Three-variable conditional analyses using the existing native packages.

Inputs contain x, y, z on the same regular forward-time grid. These functions
fit all three variables together and return all six cross-variable directions.
They do not change the established pairwise wrappers in ``toolbox.causal``.
"""

from importlib.metadata import version
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from toolbox.causal import _check_gpdc_null, _edge_row, _gpdc_null_settings


VARIABLES = ("x", "y", "z")
DIRECTIONS = tuple((source, target) for source in range(3)
                   for target in range(3) if source != target)


def _series_values(series):
    """Validate finite regular inputs without sorting or rescaling variables."""
    required = ["time", "age", *VARIABLES]
    if not isinstance(series, pd.DataFrame) or not set(required).issubset(series.columns):
        raise ValueError("series must contain time, age, x, y, and z.")
    data = series[required].to_numpy(dtype=float)
    if len(data) < 3 or not np.isfinite(data).all():
        raise ValueError("series must contain at least three finite rows.")
    spacing = np.diff(data[:, 0])
    if np.any(spacing <= 0) or not np.allclose(spacing, spacing[0], rtol=1e-8, atol=1e-12):
        raise ValueError("time must increase on a regular grid.")
    if not np.allclose(np.diff(data[:, 1]), -spacing, rtol=1e-8, atol=1e-12):
        raise ValueError("BP ages must decrease at the forward-time grid interval.")
    values = data[:, 2:]
    if np.any(np.ptp(values, axis=0) == 0):
        raise ValueError("All three causal input series must be nonconstant.")
    return values, float(spacing[0])


def run_pcmciplus_multivariate(
    series, *, tau_max=5, pc_alpha=0.05, alpha=0.05, sig_samples=10000,
    seed=1, gp_params=None, null_dist_filename=None,
):
    """Native PCMCI+/GPDC with all three variables in its conditioning search.

    tau_max is the largest tested lag in grid steps; alpha and pc_alpha are
    decision and parent-search thresholds. sig_samples and seed specify the
    native GPDC null distribution. gp_params follows native GPDC: callers use
    {'optimizer': None} for the frozen unit-kernel production configuration.

    Holm adjustment is within each direction across lags 0..tau_max, as in
    the pairwise pipeline; it is not an additional graph-wide correction.
    Retention also requires a native arrow. An unresolved contemporaneous
    link marks its two directed rows only, leaving the third pair unaffected.
    """
    from tigramite import data_processing
    from tigramite.independence_tests.gpdc import GPDC
    from tigramite.pcmci import PCMCI

    values, dt = _series_values(series)
    if not isinstance(tau_max, (int, np.integer)) or tau_max < 1 or 2 * tau_max >= len(values):
        raise ValueError("tau_max must be a positive integer shorter than half the record.")
    if not 0 < alpha < 1 or not 0 < pc_alpha < 1 or sig_samples < 1:
        raise ValueError("Require alpha/pc_alpha in (0, 1) and positive sig_samples.")
    if null_dist_filename is not None:
        settings = _gpdc_null_settings([len(values) - 2 * tau_max], sig_samples, seed)
        _check_gpdc_null(Path(null_dist_filename), settings)
    test = GPDC(
        significance="analytic", sig_samples=sig_samples, seed=seed, gp_params=gp_params,
        null_dist_filename=None if null_dist_filename is None else str(null_dist_filename),
    )
    dataframe = data_processing.DataFrame(values, var_names=list(VARIABLES))
    pcmci = PCMCI(dataframe=dataframe, cond_ind_test=test, verbosity=0)
    with warnings.catch_warnings(record=True) as fit_warnings:
        warnings.simplefilter("always")
        result = pcmci.run_pcmciplus(
            tau_min=0, tau_max=tau_max, pc_alpha=pc_alpha, fdr_method="none",
        )
    p_matrix = np.asarray(result["p_matrix"], dtype=float)
    val_matrix = np.asarray(result["val_matrix"], dtype=float)
    graph = np.asarray(result["graph"])
    expected_shape = (3, 3, tau_max + 1)
    if p_matrix.shape != expected_shape or val_matrix.shape != expected_shape or graph.shape != expected_shape:
        raise ValueError("Unexpected native PCMCI+ shape for three variables.")
    p_holm = np.full(p_matrix.shape, np.nan)
    graph_holm = graph.copy()
    rows = []
    for source, target in DIRECTIONS:
        p = p_matrix[source, target, :]
        scores = val_matrix[source, target, :]
        valid = np.isfinite(p).all() and np.all((p >= 0) & (p <= 1))
        valid = valid and np.isfinite(scores).all()
        if valid:
            adjusted = multipletests(p, alpha=alpha, method="holm")[1]
            p_holm[source, target, :] = adjusted
            detected = np.any((graph[source, target, :] == "-->") & (adjusted < alpha))
            graph_holm[source, target, adjusted >= alpha] = ""
            score = float(np.max(scores))
            p_summary = float(np.min(adjusted))
        else:
            detected, score, p_summary = False, np.nan, np.nan
            graph_holm[source, target, :] = ""
        rows.append(_edge_row(
            "pcmci+", VARIABLES[source], VARIABLES[target], score, "max_gpdc",
            detected, p_summary, "gpdc_reference_holm_within_direction",
        ))

    pair_diagnostics = {}
    for row, (source, target) in zip(rows, DIRECTIONS):
        zero_links = [str(graph_holm[source, target, 0]), str(graph_holm[target, source, 0])]
        row["unoriented"] = "o-o" in zero_links or "x-x" in zero_links
        if source < target:
            pair_diagnostics[f"{VARIABLES[source]}_{VARIABLES[target]}"] = {
                "undirected_contemporaneous": "o-o" in zero_links,
                "conflicting_contemporaneous": "x-x" in zero_links,
                "native_contemporaneous_link": str(graph[source, target, 0]),
                "retained_directional_zero_links": zero_links,
            }
    return {
        "edges": pd.DataFrame(rows),
        "diagnostics": {
            "package_version": version("tigramite"), "variables": list(VARIABLES),
            "tau_max": tau_max, "lag_steps": np.arange(tau_max + 1).tolist(),
            "lag_kyr": (np.arange(tau_max + 1) * dt).tolist(), "pc_alpha": pc_alpha,
            "alpha": alpha, "sig_samples": sig_samples, "seed": seed, "gp_params": gp_params,
            "null_distribution": "generated" if null_dist_filename is None else "precomputed",
            "contemporaneous_by_pair": pair_diagnostics,
            "n_fit_warnings": len(fit_warnings),
            "fit_warnings": sorted(set(str(w.message) for w in fit_warnings)),
        },
        "raw": {**result, "p_holm": p_holm, "graph_holm": graph_holm},
    }


def run_lkif_multivariate(series, *, max_lag=5, np_step=1, alpha=0.05):
    """Native three-variable LKIF, retaining its analytic significance and units.

    max_lag is the native lag-model size; np_step=1 is the package's recommended
    difference step. alpha is applied after Holm adjustment within each ordered
    direction. The package stores information flow as [target, source, lag].
    Native lag indices are not interpreted as physical propagation delays.
    """
    from LK_Info_Flow.multi_causality_est import multi_causality_est

    values, dt = _series_values(series)
    if not isinstance(max_lag, (int, np.integer)) or max_lag < 1:
        raise ValueError("max_lag must be a positive integer.")
    if len(values) - max_lag <= 3 * max_lag + 1:
        raise ValueError("Too few samples for the three-variable LKIF lag model.")
    if np_step != 1 or not 0 < alpha < 1:
        raise ValueError("Require np_step=1 and alpha in (0, 1).")
    with warnings.catch_warnings(record=True) as fit_warnings:
        warnings.simplefilter("always")
        result = multi_causality_est(
            values.T, max_lag=max_lag, np=np_step, dt=dt,
            series_temporal_order=np.array([]), significance_test=1,
        )
    p_matrix = np.asarray(result["p"], dtype=float)
    normalized_flow = np.asarray(result["nIF"], dtype=float)
    if p_matrix.shape != (3, 3, max_lag) or normalized_flow.shape != p_matrix.shape:
        raise ValueError("Unexpected three-variable LK_Info_Flow output shape.")
    p_holm = np.full(p_matrix.shape, np.nan)
    signed_max, rows = {}, []
    for source, target in DIRECTIONS:
        p = p_matrix[target, source, :]
        flow = normalized_flow[target, source, :]
        valid = np.isfinite(p).all() and np.all((p >= 0) & (p <= 1))
        valid = valid and np.isfinite(flow).all()
        direction = f"{VARIABLES[source]}_to_{VARIABLES[target]}"
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
        row = _edge_row(
            "lkif", VARIABLES[source], VARIABLES[target], score,
            "max_abs_normalized_information_flow", detected, p_summary,
            "native_normal_approximation_holm_within_direction",
        )
        row["unoriented"] = False
        rows.append(row)
    return {
        "edges": pd.DataFrame(rows),
        "diagnostics": {
            "package_version": version("LK_Info_Flow"), "variables": list(VARIABLES),
            "max_lag": max_lag, "native_lag_indices": np.arange(max_lag).tolist(),
            "dt_kyr": dt, "np_step": np_step, "alpha": alpha,
            "signed_max_normalized_flow": signed_max,
            "temporal_order": "contiguous_integer_samples",
            "n_fitted_samples": len(values) - max_lag,
            "n_fit_warnings": len(fit_warnings),
            "fit_warnings": sorted(set(str(w.message) for w in fit_warnings)),
        },
        "raw": {**result, "p_holm": p_holm},
    }


def _ols_fields(fit):
    """Retain scientific native regression fields without a pickled fit object."""
    fields = ["params", "resid", "fittedvalues", "normalized_cov_params", "bse", "tvalues",
              "pvalues", "ssr", "scale", "nobs", "df_model", "df_resid", "rsquared",
              "rsquared_adj", "centered_tss", "uncentered_tss", "ess", "mse_model", "mse_resid",
              "mse_total", "fvalue", "f_pvalue", "llf", "aic", "bic", "condition_number"]
    result = {name: getattr(fit, name) for name in fields}
    result.update({"cov_params": fit.cov_params(), "conf_int": fit.conf_int(),
                   "endog": fit.model.endog, "exog": fit.model.exog,
                   "design_rank": fit.model.rank})
    return result


def run_granger_multivariate(series, *, order=1, alpha=0.05):
    """Native conditional GC in a three-variable VAR of fixed grid-step order.

    Each native OLS F-test removes one source's lag coefficients from one
    target equation while retaining the histories of the other two variables.
    Its single-equation denominator degrees of freedom match the existing
    pairwise SSR F-test convention. The native VAR system F-test, with its
    different denominator degrees of freedom, is saved as a diagnostic only.
    No order selection or contemporaneous-effect test is performed.

    The score remains log(SSR_restricted / SSR_full) on identical response rows.
    Native OLS fits provide this score, full fit records and the F-test p value;
    no custom regression estimator or reference distribution is introduced.
    """
    from statsmodels.regression.linear_model import OLS
    from statsmodels.tsa.api import VAR

    values, dt = _series_values(series)
    if not isinstance(order, (int, np.integer)) or order < 1 or not 0 < alpha < 1:
        raise ValueError("Require a positive integer VAR order and alpha in (0, 1).")
    if len(values) - order <= 1 + 3 * order:
        raise ValueError("Too few rows for the three-variable VAR.")
    with warnings.catch_warnings(record=True) as fit_warnings:
        warnings.simplefilter("always")
        var_fit = VAR(pd.DataFrame(values, columns=VARIABLES)).fit(order, trend="c")
    design = np.asarray(var_fit.endog_lagged)
    responses = values[order:]
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("Conditional VAR design is rank deficient.")
    full_fits = [OLS(responses[:, target], design).fit() for target in range(3)]
    # The native VAR and its separate OLS equations must use the same samples.
    for target, fit in enumerate(full_fits):
        if fit.nobs != var_fit.nobs or fit.model.rank < len(fit.params):
            raise ValueError("Conditional GC equation has invalid sample count or rank.")
        if not np.allclose(fit.params, np.asarray(var_fit.params)[:, target], rtol=1e-10, atol=1e-12):
            raise ValueError("Native VAR and OLS coefficient layouts do not agree.")

    whiteness_lags = min(max(10, order + 1), var_fit.nobs - 1)
    whiteness = None
    try:
        whiteness = vars(var_fit.test_whiteness(whiteness_lags, adjusted=True)).copy()
        whiteness_status = "ok"
    except np.linalg.LinAlgError:
        whiteness_status = "singular_residual_covariance"

    rows, edge_fits, fit_diagnostics = [], {}, {}
    for source, target in DIRECTIONS:
        var_test = var_fit.test_causality(VARIABLES[target], [VARIABLES[source]], kind="f", signif=alpha)
        # Native VAR order is intercept, then x/y/z at lag 1, lag 2, and so on.
        removed = [1 + lag * 3 + source for lag in range(order)]
        keep = [column for column in range(design.shape[1]) if column not in removed]
        restricted = OLS(responses[:, target], design[:, keep]).fit()
        full = full_fits[target]
        if restricted.nobs != full.nobs or restricted.model.rank < len(restricted.params):
            raise ValueError("Restricted conditional GC fit has invalid row count or rank.")
        score = float(np.log(restricted.ssr / full.ssr))
        contrast = np.zeros((order, design.shape[1]))
        for lag, column in enumerate(removed):
            contrast[lag, column] = 1.0
        test = full.f_test(contrast)
        p_value = float(test.pvalue)
        direction = f"{VARIABLES[source]}_to_{VARIABLES[target]}"
        row = _edge_row(
            "granger", VARIABLES[source], VARIABLES[target], score, "log_residual_variance_ratio",
            p_value < alpha, p_value, "native_OLS_conditional_equation_F",
        )
        row.update({"unoriented": False, "selected_order": order})
        rows.append(row)
        edge_fits[direction] = {
            "test": {"f_statistic": float(test.fvalue), "p_value": p_value,
                     "df_num": float(test.df_num), "df_denom": float(test.df_denom)},
            "native_VAR_test": vars(var_test).copy(), "full": _ols_fields(full),
            "restricted": _ols_fields(restricted), "removed_design_columns": removed,
            "retained_design_columns": keep, "contrast": contrast,
        }
        fit_diagnostics[direction] = {
            "f_statistic": float(test.fvalue), "df_num": float(test.df_num),
            "df_denom": float(test.df_denom), "nobs": int(full.nobs),
            "ssr_restricted": float(restricted.ssr), "ssr_full": float(full.ssr),
            "condition_number": float(full.condition_number), "design_rank": int(full.model.rank),
            "n_parameters": len(full.params), "conditioned_variables": [
                VARIABLES[index] for index in range(3) if index not in (source, target)],
        }
    var_fields = ["endog", "endog_lagged", "params", "resid", "fittedvalues", "stderr", "tvalues",
                  "pvalues", "sigma_u", "sigma_u_mle", "coefs", "intercept", "nobs", "df_model",
                  "df_resid", "k_ar", "k_trend", "k_exog", "names", "exog_names", "llf",
                  "aic", "bic", "hqic", "fpe"]
    var_snapshot = {name: getattr(var_fit, name) for name in var_fields}
    var_snapshot.update({"cov_params": var_fit.cov_params(),
                         "roots_real": np.real(var_fit.roots), "roots_imag": np.imag(var_fit.roots),
                         "is_stable": bool(var_fit.is_stable())})
    return {
        "edges": pd.DataFrame(rows),
        "diagnostics": {
            "package_version": version("statsmodels"), "variables": list(VARIABLES),
            "selected_order": order, "order_selection": None, "history_kyr": order * dt,
            "alpha": alpha, "fits": fit_diagnostics, "var_is_stable": bool(var_fit.is_stable()),
            "test": "native OLS conditional F; single-equation denominator degrees of freedom",
            "residual_whiteness_p": None if whiteness is None else float(whiteness["pvalue"]),
            "residual_whiteness_status": whiteness_status,
            "n_fit_warnings": len(fit_warnings),
            "fit_warnings": sorted(set(str(w.message) for w in fit_warnings)),
        },
        "raw": {"var": var_snapshot, "directions": edge_fits, "whiteness": whiteness},
    }
