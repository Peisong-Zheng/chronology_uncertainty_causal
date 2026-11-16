"""Verify the fixed-memory common-driver experiment before production.

Inputs, explicit failures, resume guards and complete native outputs are checked
in an isolated timestamped benchmark directory. This script never launches the
production experiment or changes earlier experiment outputs.
"""

import os

for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"]:
    os.environ[name] = "1"
os.environ.setdefault("MPLBACKEND", "Agg")

from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
from datetime import datetime, timezone
import gzip
import hashlib
import importlib
import json
import multiprocessing
from pathlib import Path
import shutil
import sys
import time
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy.linalg import solve_discrete_lyapunov
from threadpoolctl import threadpool_limits

core = importlib.import_module("03_check_common_driver")

ANALYSES = [("pairwise", method) for method in ["PCMCI+", "LKIF", "GC", "CCM"]]
ANALYSES = ANALYSES + [("conditional", method) for method in ["PCMCI+", "LKIF", "GC"]]


def file_state(path):
    return {"sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "mtime_ns": path.stat().st_mtime_ns}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True, allow_nan=True).encode("utf-8")


def native_path(output, job, analysis, method):
    return Path(output) / "native" / f"{job['case_id']}_{analysis}_{method.replace('+', 'plus')}.json.gz"


def expected_edges(analysis):
    variables = ["x", "y"] if analysis == "pairwise" else ["x", "y", "z"]
    return [(source, target) for source in variables for target in variables if source != target]


def expect_rejection(function, *args):
    try:
        function(*args)
    except (ValueError, RuntimeError, AssertionError):
        return
    raise AssertionError(f"Scientific mismatch accepted by {function.__name__}.")


def check_inputs(jobs, prepared):
    """Check the physical equations and exact pairing of all prepared inputs."""
    cases = pd.DataFrame(jobs)
    assert len(cases) == 192 and cases.case_id.nunique() == 192
    assert cases.groupby(["structure", "process_seed"]).size().eq(32).all()
    assert cases.groupby("condition").size().to_dict() == {
        "all_uncertain": 60, "irregular_accurate": 6, "original_uniform": 6,
        "xy_only": 60, "z_only": 60}
    assert set(cases.structure) == {"h0_common_driver", "h1_direct"}
    assert set(cases.process_seed) == {20261009, 20261301, 20261302}
    grid = np.linspace(29, 1, 1401)
    loaded = {}
    for row in cases.to_dict("records"):
        frame = pd.read_csv(prepared / "inputs" / f"{row['case_id']}.csv", float_precision="round_trip")
        assert list(frame) == ["time", "age", "x", "y", "z"]
        assert frame.shape == (1401, 5) and np.isfinite(frame.to_numpy()).all()
        np.testing.assert_allclose(frame.age, grid, rtol=0, atol=1e-12)
        np.testing.assert_allclose(frame.time, grid[0] - grid, rtol=0, atol=1e-12)
        key = (row["process_seed"], row["structure"], row["condition"], row["member"])
        loaded[key] = frame
    for seed in cases.process_seed.unique():
        for condition in cases.condition.unique():
            members = [-1] if condition in {"original_uniform", "irregular_accurate"} else range(10)
            for member in members:
                h0 = loaded[seed, "h0_common_driver", condition, member]
                h1 = loaded[seed, "h1_direct", condition, member]
                np.testing.assert_array_equal(h0[["x", "z"]], h1[["x", "z"]])
        for structure in cases.structure.unique():
            accurate = loaded[seed, structure, "irregular_accurate", -1]
            for member in range(10):
                only_z = loaded[seed, structure, "z_only", member]
                only_xy = loaded[seed, structure, "xy_only", member]
                both = loaded[seed, structure, "all_uncertain", member]
                np.testing.assert_array_equal(only_z[["x", "y"]], accurate[["x", "y"]])
                np.testing.assert_array_equal(only_xy.z, accurate.z)
                np.testing.assert_array_equal(both[["x", "y"]], only_xy[["x", "y"]])
                np.testing.assert_array_equal(both.z, only_z.z)
                assert not np.array_equal(only_z.z, accurate.z)
                assert not np.array_equal(only_xy[["x", "y"]], accurate[["x", "y"]])

    clock_arrays = []
    sample_indices = None
    for seed in cases.process_seed.unique():
        for record in ["x", "y", "z"]:
            path = prepared / "clocks" / f"s{seed}_{record}_C.npz"
            with np.load(path, allow_pickle=False) as clock:
                age_true, age_obs, ages = clock["age_true"], clock["age_obs"], clock["ages"]
                assert ages.shape == (10, len(age_true))
                assert np.all(np.diff(age_true) < 0) and np.all(np.diff(age_obs) < 0)
                assert np.all(np.diff(ages, axis=1) < 0)
                assert np.all(ages[:, 0] >= 29) and np.all(ages[:, -1] <= 1)
                assert np.isfinite(ages).all() and not np.array_equal(age_obs, age_true)
                np.testing.assert_array_equal(clock["age_std"], np.full(len(age_true), .04))
                np.testing.assert_allclose(age_true, 30 - .02 * clock["native_indices"], rtol=0, atol=1e-12)
                np.testing.assert_array_equal(clock["sample_id"], np.arange(len(age_true)))
                if sample_indices is None:
                    sample_indices = clock["native_indices"].copy()
                else:
                    np.testing.assert_array_equal(clock["native_indices"], sample_indices)
                for prior in clock_arrays:
                    assert not np.array_equal(ages, prior)
                clock_arrays.append(ages.copy())

    # Independent calculations use the stated background units, not sample SDs.
    a = np.exp(-20 / np.array([500., 200., 200.]))  # State order: Z, X, Y.
    innovation_sd = np.sqrt(1 - a**2)
    theoretical = []
    recovered_innovations = {}
    for seed in cases.process_seed.unique():
        for structure, direct in [("h0_common_driver", 0), ("h1_direct", 1)]:
            matrix = np.array([[a[0], 0, 0], [(1-a[1])*2, a[1], 0],
                               [(1-a[2])*2, (1-a[2])*direct, a[2]]])
            covariance = solve_discrete_lyapunov(matrix, np.diag(innovation_sd**2))
            assert np.max(np.abs(np.linalg.eigvals(matrix))) < 1
            assert np.min(np.linalg.eigvalsh(covariance)) > 0
            np.testing.assert_allclose(covariance, matrix @ covariance @ matrix.T
                                       + np.diag(innovation_sd**2), rtol=0, atol=1e-13)
            path = prepared / "generated" / f"s{seed}_{structure}.csv"
            full = pd.read_csv(path, float_precision="round_trip")
            assert len(full) == 1501
            values = full[["z", "x", "y"]].to_numpy()
            innovations = (values[1:] - values[:-1] @ matrix.T) / innovation_sd
            # A shared 250-step warm-up plus one history step fixes retained
            # positions in the independent three-column innovation stream.
            expected = np.random.default_rng(int(seed)).normal(size=(250 + 1 + 1501, 3))[252:]
            np.testing.assert_allclose(innovations, expected, rtol=0, atol=2e-14)
            if seed in recovered_innovations:
                np.testing.assert_allclose(innovations, recovered_innovations[seed], rtol=0, atol=2e-14)
            else:
                recovered_innovations[seed] = innovations
            for first, second in [(0, 1), (0, 2), (1, 2)]:
                assert not np.array_equal(innovations[:, first], innovations[:, second])
            theoretical.append({"process_seed": int(seed), "structure": structure,
                                "spectral_radius": float(np.max(np.abs(np.linalg.eigvals(matrix)))),
                                "covariance": covariance.tolist()})
    return {"prepared_cases": len(cases), "paired_structure_inputs": True,
            "only_z_preserves_xy_exactly": True, "only_xy_preserves_z_exactly": True,
            "combined_clocks_match_component_conditions": True,
            "nine_distinct_ordered_clock_sets": len(clock_arrays),
            "common_sampling_mask_and_constant_proposal_sd": True,
            "paired_independent_innovations": True, "theoretical_systems": theoretical}


def fake_method(series, analysis, method, estimator_seed, null_file):
    """Return defined mixed outcomes to test storage without fitting anything."""
    if method == "LKIF":
        raise RuntimeError("Intentional verification-only method failure.")
    rows = []
    for source, target in expected_edges(analysis):
        row = {"source": source, "target": target, "score": 0.3,
               "p_value": 0.01, "edge_detected": source == "x", "status": "ok",
               "unoriented": analysis == "conditional" and method == "PCMCI+"
                             and {source, target} == {"x", "y"}}
        if method == "CCM":
            row.update(raw_edge_detected=source == "x", lag_status="compatible",
                       unoriented=False)
            if source == "x":
                row.update(edge_detected=None, lag_status="failed")
        rows.append(row)
    return {"edges": pd.DataFrame(rows), "diagnostics": {}, "raw": {"array": np.arange(6.).reshape(3, 2)}}


def check_guards(jobs, prepared, output):
    job = jobs[0]
    source = prepared / "inputs" / f"{job['case_id']}.csv"
    target = output / "failure_retention"
    (target / "inputs").mkdir(parents=True)
    shutil.copyfile(source, target / "inputs" / source.name)
    with patch.object(core, "call_method", side_effect=fake_method):
        result = core.run_case(job, target, "unused_mock_null.npz")
    assert len(result) == 26
    assert result.loc[result.method.eq("LKIF"), "status"].eq("error").all()
    assert result.loc[result.method.eq("LKIF"), "edge_detected"].isna().all()
    selected = result.loc[result.analysis.eq("conditional") & result.method.eq("PCMCI+")]
    assert selected.loc[selected.source.eq("z") | selected.target.eq("z"), "unoriented"].eq(False).all()
    assert selected.loc[~(selected.source.eq("z") | selected.target.eq("z")), "unoriented"].eq(True).all()
    before = {p.relative_to(target).as_posix(): file_state(p) for p in target.rglob("*") if p.is_file()}
    resumed = core.read_completed_case(job, target)
    after = {p.relative_to(target).as_posix(): file_state(p) for p in target.rglob("*") if p.is_file()}
    assert resumed is not None and len(resumed) == 26 and before == after
    assert resumed.loc[resumed.method.eq("CCM") & resumed.source.eq("x"), "score"].iloc[0] == 0.3
    assert resumed.loc[resumed.method.eq("CCM") & resumed.source.eq("x"), "edge_detected"].isna().all()

    for label, shortened, missing in [("short_commit", True, None),
                                      ("missing_native", False, ("conditional", "GC"))]:
        folder = output / label
        (folder / "members").mkdir(parents=True)
        (folder / "native").mkdir()
        table = result.iloc[:-1] if shortened else result
        table.to_csv(folder / "members" / f"{job['case_id']}_edges.csv", index=False)
        for analysis, method in ANALYSES:
            if (analysis, method) != missing:
                path = native_path(target, job, analysis, method)
                shutil.copyfile(path, folder / "native" / path.name)
        assert core.read_completed_case(job, folder) is None

    failed = output / "serialization_failure"
    (failed / "inputs").mkdir(parents=True)
    shutil.copyfile(source, failed / "inputs" / source.name)

    def unserializable(series, analysis, method, seed, null):
        result = fake_method(series, analysis, "GC" if method == "LKIF" else method, seed, null)
        if method == "GC":
            result["raw"]["unsupported"] = object()
        return result

    with patch.object(core, "call_method", side_effect=unserializable):
        try:
            core.run_case(job, failed, "unused_mock_null.npz")
        except TypeError:
            pass
        else:
            raise AssertionError("Serialization error was committed as scientific evidence.")
    assert not (failed / "members" / f"{job['case_id']}_edges.csv").exists()

    parameter_dir = output / "parameter_guards"
    parameter_dir.mkdir()
    parameters = core.execution_parameters()
    (parameter_dir / "parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    core.check_parameters(parameters, parameter_dir)
    changed = deepcopy(parameters)
    changed["N_WORKERS"] = 16
    core.check_parameters(changed, parameter_dir)
    checked = []
    for key in ["ALPHA", "CCM_LIBRARY_SIZE", "BURN_IN_STEPS", "N_AGE_MEMBERS"]:
        if key in parameters:
            changed = deepcopy(parameters)
            changed[key] = parameters[key] + 1
            expect_rejection(core.check_parameters, changed, parameter_dir)
            checked.append(key)
    for key in ["packages", "environment", "source_sha256", "implementation_sha256"]:
        if key in parameters:
            changed = deepcopy(parameters)
            changed[key] = {"verification": "changed"}
            expect_rejection(core.check_parameters, changed, parameter_dir)
            checked.append(key)
    original_state = file_state(source)
    changed_input = pd.read_csv(source, float_precision="round_trip")
    changed_input.loc[1, "x"] = changed_input.loc[1, "x"] + .01
    expect_rejection(core.save_input, changed_input, source)
    assert file_state(source) == original_state
    return {"failed_fits_retained": True, "unresolved_is_pair_specific": True,
            "ccm_raw_score_retained_after_lag_failure": True, "resume_no_writes": True,
            "partial_commit_rejected": True, "missing_native_rejected": True,
            "serialization_failure_not_committed": True, "parameter_guards": checked,
            "changed_regenerated_input_rejected_without_write": True}


def check_conditional_mapping(series, method, result):
    """Independently check direction indexing and native equation-F decisions."""
    from statsmodels.stats.multitest import multipletests

    variables = ["x", "y", "z"]
    rows = result["edges"].set_index(["source", "target"])
    assert set(rows.index) == set(expected_edges("conditional"))
    raw = result["raw"]
    if method in {"PCMCI+", "LKIF"}:
        for source, target in expected_edges("conditional"):
            first, second = variables.index(source), variables.index(target)
            index = (first, second) if method == "PCMCI+" else (second, first)
            p = raw["p_matrix" if method == "PCMCI+" else "p"][index]
            scores = raw["val_matrix" if method == "PCMCI+" else "nIF"][index]
            adjusted = multipletests(p, alpha=.05, method="holm")[1]
            expected_score = float(np.max(scores)) if method == "PCMCI+" else float(np.max(np.abs(scores)))
            eligible = raw["graph"][index] == "-->" if method == "PCMCI+" else np.ones(len(p), dtype=bool)
            assert rows.loc[(source, target), "score"] == expected_score
            assert rows.loc[(source, target), "p_value"] == float(np.min(adjusted))
            assert bool(rows.loc[(source, target), "edge_detected"]) == bool(np.any(eligible & (adjusted < .05)))
            np.testing.assert_array_equal(raw["p_holm"][index], adjusted)
        return {"method": method, "all_six_native_directions_checked": True}

    from statsmodels.regression.linear_model import OLS
    from statsmodels.tsa.api import VAR

    values = series[variables].to_numpy(dtype=float)
    fit = VAR(pd.DataFrame(values, columns=variables)).fit(1, trend="c")
    design = np.asarray(fit.endog_lagged)
    for source, target in expected_edges("conditional"):
        source_index, target_index = variables.index(source), variables.index(target)
        response = values[1:, target_index]
        full = OLS(response, design).fit()
        keep = [column for column in range(4) if column != 1 + source_index]
        restricted = OLS(response, design[:, keep]).fit()
        contrast = np.zeros((1, 4))
        contrast[0, 1 + source_index] = 1
        native = full.f_test(contrast)
        record = raw["directions"][f"{source}_to_{target}"]
        np.testing.assert_allclose(rows.loc[(source, target), "p_value"], float(native.pvalue), rtol=1e-10, atol=1e-14)
        np.testing.assert_allclose(rows.loc[(source, target), "score"],
                                   np.log(restricted.ssr / full.ssr), rtol=1e-10, atol=1e-14)
        assert record["test"]["df_denom"] == len(response) - 4
        assert record["test"]["df_num"] == 1
        for name, native_fit in [("full", full), ("restricted", restricted)]:
            saved = record[name]
            for field in ["params", "resid", "fittedvalues"]:
                np.testing.assert_allclose(saved[field], getattr(native_fit, field), rtol=1e-10, atol=1e-12)
            np.testing.assert_array_equal(saved["endog"], response)
            assert saved["nobs"] == len(response)
        system_test = fit.test_causality(target, [source], kind="f", signif=.05)
        np.testing.assert_allclose(record["native_VAR_test"]["pvalue"], system_test.pvalue,
                                   rtol=1e-10, atol=1e-14)
    return {"method": method, "all_six_native_directions_checked": True,
            "primary_test": "native single-equation OLS F",
            "secondary_test": "native VAR system F retained separately",
            "complete_fits_and_common_response_rows_checked": True}


def check_pairwise_reuse(jobs, prepared, output):
    """Identical X/Y evidence is copied exactly, including recorded failures."""
    source_job = next(job for job in jobs if job["structure"] == "h0_common_driver"
                      and job["process_seed"] == 20261009 and job["condition"] == "irregular_accurate")
    target_job = next(job for job in jobs if job["structure"] == "h0_common_driver"
                      and job["process_seed"] == 20261009 and job["condition"] == "z_only"
                      and job["member"] == 0)
    folder = output / "pairwise_reuse"
    (folder / "inputs").mkdir(parents=True)
    pd.DataFrame(jobs).to_csv(folder / "cases.csv", index=False)
    for job in [source_job, target_job]:
        path = prepared / "inputs" / f"{job['case_id']}.csv"
        shutil.copyfile(path, folder / "inputs" / path.name)
    with patch.object(core, "call_method", side_effect=fake_method):
        source = core.run_case(source_job, folder, "unused_mock_null.npz")
    calls = []

    def conditional_only(series, analysis, method, seed, null):
        assert analysis == "conditional", "Identical pairwise input was needlessly recomputed."
        calls.append((analysis, method))
        return fake_method(series, analysis, method, seed, null)

    with patch.object(core, "call_method", side_effect=conditional_only):
        target = core.run_case(target_job, folder, "unused_mock_null.npz")
    assert len(calls) == 3 and len(target) == 26
    assert target.loc[target.analysis.eq("pairwise") & target.method.eq("LKIF"), "status"].eq("error").all()
    compared = []
    for method in ["PCMCI+", "LKIF", "GC", "CCM"]:
        first = native_path(folder, source_job, "pairwise", method)
        second = native_path(folder, target_job, "pairwise", method)
        assert first.read_bytes() == second.read_bytes(), method
        compared.append(method)
    keys = ["method", "source", "target"]
    scientific = ["score", "p_value", "edge_detected", "status", "unoriented"]
    first = source.loc[source.analysis.eq("pairwise")].set_index(keys)[scientific]
    second = target.loc[target.analysis.eq("pairwise")].set_index(keys)[scientific]
    first = first.astype({"edge_detected": "boolean", "unoriented": "boolean"})
    second = second.astype({"edge_detected": "boolean", "unoriented": "boolean"})
    pd.testing.assert_frame_equal(first, second, check_dtype=False, check_exact=True)
    input_path = folder / "inputs" / f"{target_job['case_id']}.csv"
    original = input_path.read_bytes()
    changed = pd.read_csv(input_path, float_precision="round_trip")
    changed.loc[1, "x"] = changed.loc[1, "x"] + .01
    changed.to_csv(input_path, index=False)
    try:
        with patch.object(core, "call_method", side_effect=conditional_only):
            expect_rejection(core.run_case, target_job, folder, "unused_mock_null.npz")
    finally:
        input_path.write_bytes(original)
    changed_job = dict(target_job)
    changed_job["estimator_seed"] = changed_job["estimator_seed"] + 1
    with patch.object(core, "call_method", side_effect=conditional_only):
        expect_rejection(core.run_case, changed_job, folder, "unused_mock_null.npz")
    return {"byte_identical_native_archives": compared, "only_three_conditional_calls": True,
            "source_scientific_failure_retained": True, "changed_xy_input_rejected": True,
            "changed_estimator_seed_rejected": True}


def check_summaries(jobs):
    """Separate missing decisions, unresolved evidence and valid raw scores."""
    subset = [job for job in jobs if job["process_seed"] == 20261009]
    pending = core.summarize_results(core.collect_results(subset, {}))
    assert pending["xy_case_states"].pending.all()
    assert not pending["xy_case_states"].failed.any()
    assert not pending["xy_case_states"].correct_xy.any()
    assert not pending["direction_support"].support_eligible.any()
    rows = []
    for job in subset:
        truth = {("z", "x"), ("z", "y")}
        if job["structure"] == "h1_direct":
            truth.add(("x", "y"))
        for analysis, method in ANALYSES:
            for source, target in expected_edges(analysis):
                perturbed = job["condition"] in {"z_only", "xy_only", "all_uncertain"}
                row = {**job, "analysis": analysis, "method": method, "source": source, "target": target,
                       "score": .4 if perturbed else .6, "p_value": .01 if (source, target) in truth else .7,
                       "edge_detected": (source, target) in truth, "status": "ok", "unoriented": False,
                       "raw_edge_detected": True if method == "CCM" else pd.NA,
                       "lag_status": "compatible" if method == "CCM" else ""}
                member = job["member"]
                if perturbed:
                    if member == 0:
                        row.update(status="error", score=np.nan, edge_detected=None)
                    elif member == 1 and {source, target} == {"x", "y"}:
                        row["unoriented"] = True
                    elif member == 2 and source == "y" and target == "x":
                        row["edge_detected"] = True
                    elif member == 3 and source == "x" and target == "y":
                        row["score"] = np.nan
                    elif member == 4 and method == "CCM" and source == "x":
                        row.update(edge_detected=None, lag_status="failed")
                    elif member == 5 and analysis == "conditional" and source == "z" and target == "x":
                        row.update(status="error", score=np.nan, edge_detected=None)
                rows.append(row)
    results = pd.DataFrame(rows)
    for column in ["edge_detected", "unoriented", "raw_edge_detected"]:
        results[column] = results[column].astype("boolean")
    tables = core.summarize_results(results)
    for row in tables["direction_support"].query("condition in ['z_only', 'xy_only', 'all_uncertain']").itertuples():
        is_ccm = row.method == "CCM"
        assert row.requested == 10 and row.usable == (8 if is_ccm else 9)
        assert row.failed == (2 if is_ccm else 1) and row.pending == 0
        assert row.unresolved == 1 and row.correct_xy == (6 if is_ccm else 7)
        assert row.q_y_to_x == 1 / row.usable
        assert row.q_x_to_y == (1.0 if row.structure == "h1_direct" else 0.0)
        assert row.support_eligible == (not is_ccm)
    xy = tables["xy_case_states"]
    separate_failure = xy.loc[xy.analysis.eq("conditional") & xy.member.eq(5)]
    assert separate_failure.usable.all() and separate_failure.correct_xy.all()
    assert not separate_failure.failed.any()
    network = tables["conditional_graph_states"]
    assert not network.loc[network.member.eq(5), "usable"].any()
    changes = tables["paired_score_changes"]
    raw_survives = changes.loc[changes.method.eq("CCM") & changes.member.eq(4) & changes.source.eq("x")]
    assert raw_survives.usable_score_pair.all()
    np.testing.assert_allclose(raw_survives.score_change, -.2, rtol=0, atol=1e-15)
    broken = results.copy()
    selected = broken.condition.eq("irregular_accurate") & broken.method.eq("GC") & broken.source.eq("x") & broken.target.eq("y")
    broken.loc[selected, "status"] = "error"
    broken.loc[selected, "score"] = np.nan
    affected = core.summarize_results(broken)["paired_score_changes"]
    affected = affected.loc[affected.method.eq("GC") & affected.source.eq("x") & affected.target.eq("y")]
    assert not affected.usable_score_pair.any() and affected.score_change.isna().all()
    return {"requested_denominators_and_pending_separated": True,
            "paired_valid_q_with_coverage_threshold": True, "unresolved_not_correct": True,
            "xy_failure_not_contaminated_by_other_edges": True,
            "ccm_lag_failure_counted_but_raw_score_retained": True,
            "paired_scores_require_successful_finite_reference": True}


def check_native(jobs, prepared, output, report, save_report):
    """Match full native evidence between direct calls and Windows workers."""
    selected = [next(job for job in jobs if job["process_seed"] == 20261009
                     and job["structure"] == structure and job["condition"] == condition
                     and job["member"] == member)
                for structure, condition, member in [("h0_common_driver", "irregular_accurate", -1),
                                                      ("h1_direct", "xy_only", 0)]]
    null = output / "cache" / f"null_T1391_n10000_seed{selected[0]['estimator_seed']}.npz"
    core.prepare_gpdc_null(null, sample_sizes=[1391], sig_samples=10000,
                           seed=selected[0]["estimator_seed"])
    parallel = output / "parallel"
    (parallel / "inputs").mkdir(parents=True)
    reference = {}
    for job in selected:
        source = prepared / "inputs" / f"{job['case_id']}.csv"
        shutil.copyfile(source, parallel / "inputs" / source.name)
        series = pd.read_csv(source)
        with threadpool_limits(limits=1):
            for analysis, method in ANALYSES:
                print(f"Direct check: {job['case_id']} {analysis} {method}", flush=True)
                result = core.call_method(series, analysis, method, job["estimator_seed"], str(null))
                assert result["edges"].status.eq("ok").all(), result["edges"].to_string()
                if analysis == "conditional":
                    report.setdefault("conditional_mapping", []).append(
                        {"case_id": job["case_id"], **check_conditional_mapping(series, method, result)})
                reference[job["case_id"], analysis, method] = core.native_snapshot(result)
        report.setdefault("direct_completed", []).append(job["case_id"])
        save_report()
    with ProcessPoolExecutor(max_workers=2, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(core.run_case, job, str(parallel), str(null)) for job in selected]
        for future in as_completed(futures):
            result = future.result()
            assert len(result) == 26 and result.status.eq("ok").all()
            report.setdefault("spawn_completed", []).append(result.case_id.iloc[0])
            save_report()
    comparisons = []
    for job in selected:
        for analysis, method in ANALYSES:
            path = native_path(parallel, job, analysis, method)
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                saved = json.load(handle)
            expected = canonical(reference[job["case_id"], analysis, method])
            assert expected == canonical(saved), (job["case_id"], analysis, method)
            if method == "CCM":
                assert saved["raw"]["lag_rho"]["shape"] == [100, 31, 2]
                assert saved["raw"]["rho_null"]["shape"] == [100, 2]
            comparisons.append({"case_id": job["case_id"], "analysis": analysis,
                                "method": method, "exact_complete_native_match": True})
    return {"complete_native_comparisons": comparisons, "ignored_fields": []}


def frame_from_snapshot(snapshot):
    return pd.DataFrame({name: snapshot["data"][name]["values"] for name in snapshot["columns"]})


def audit_outputs(output_dir, reference_inputs=None):
    """Audit every committed artifact and rebuild summaries without any fit."""
    output = Path(output_dir)
    cases = pd.read_csv(output / "cases.csv", float_precision="round_trip")
    jobs = cases.to_dict("records")
    assert len(jobs) == 192 and cases.case_id.nunique() == 192
    reference_count = 0
    for case in cases.to_dict("records"):
        path = output / "inputs" / f"{case['case_id']}.csv"
        if "input_sha256" in case:
            assert hashlib.sha256(path.read_bytes()).hexdigest() == case["input_sha256"]
        if reference_inputs is not None:
            reference = Path(reference_inputs) / "inputs" / path.name
            assert path.read_bytes() == reference.read_bytes(), path.name
            reference_count = reference_count + 1
    records, failures = [], []
    native_count, gc_fit_count, ccm_count, copied_count = 0, 0, 0, 0
    for job in jobs:
        rows = core.read_completed_case(job, output)
        assert rows is not None and len(rows) == 26, job["case_id"]
        records.append(rows)
        for analysis, method in ANALYSES:
            path = native_path(output, job, analysis, method)
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                payload = json.load(handle)
            native_count = native_count + 1
            assert set(payload) >= {"edges", "raw", "diagnostics"}
            native_edges = frame_from_snapshot(payload["edges"])
            selected = rows.loc[rows.analysis.eq(analysis) & rows.method.eq(method)].set_index(["source", "target"])
            assert set(selected.index) == set(expected_edges(analysis))
            for _, edge in native_edges.iterrows():
                committed = selected.loc[(edge.source, edge.target)]
                for field in ["status", "score", "p_value", "edge_detected", "unoriented",
                              "raw_edge_detected", "optimal_lag_steps", "lag_status"]:
                    if field not in native_edges:
                        continue
                    first, second = edge[field], committed[field]
                    assert (pd.isna(first) and pd.isna(second)) or first == second, (path.name, field)
            source_id = job.get("pairwise_source_case", job["case_id"])
            if analysis == "pairwise" and source_id != job["case_id"]:
                source_job = {"case_id": source_id}
                assert path.read_bytes() == native_path(output, source_job, analysis, method).read_bytes()
                copied_count = copied_count + 1
            if not selected.status.eq("ok").all():
                failures.append({"case_id": job["case_id"], "analysis": analysis, "method": method,
                                 "statuses": selected.status.tolist()})
                continue
            if method == "CCM":
                assert payload["raw"]["lag_rho"]["shape"] == [100, 31, 2]
                assert payload["raw"]["rho_null"]["shape"] == [100, 2]
                assert len(payload["raw"]["lag_library_seeds"]["values"]) == 100
                ccm_count = ccm_count + 1
            elif method == "PCMCI+":
                dimension = 2 if analysis == "pairwise" else 3
                for field in ["p_matrix", "val_matrix", "graph", "p_holm", "graph_holm"]:
                    assert payload["raw"][field]["shape"] == [dimension, dimension, 6]
            elif method == "LKIF":
                dimension = 2 if analysis == "pairwise" else 3
                for field in ["p", "nIF", "p_holm"]:
                    assert payload["raw"][field]["shape"] == [dimension, dimension, 5]
            elif method == "GC":
                fits = payload["raw"] if analysis == "pairwise" else payload["raw"]["directions"]
                for source, target in expected_edges(analysis):
                    direction = f"{source}_to_{target}"
                    for model in ["restricted", "full"]:
                        fields = fits[direction][model]
                        assert {"params", "resid", "fittedvalues", "cov_params", "endog", "exog",
                                "df_resid", "llf", "ssr"}.issubset(fields)
                        gc_fit_count = gc_fit_count + 1
    aggregate = pd.concat(records, ignore_index=True)
    assert len(aggregate) == 4992 and native_count == 1344 and copied_count == 480
    assert len(list((output / "members").glob("*_edges.csv"))) == 192
    assert len(list((output / "native").glob("*.json.gz"))) == 1344
    saved = pd.read_csv(output / "member_results.csv", float_precision="round_trip",
                        dtype={"edge_detected": "boolean", "unoriented": "boolean", "raw_edge_detected": "boolean"})
    pd.testing.assert_frame_equal(aggregate.reset_index(drop=True), saved.reset_index(drop=True),
                                  check_dtype=False, check_exact=True)
    tables = core.summarize_results(aggregate)
    checked_tables = []
    for name, frame in tables.items():
        path = output / f"{name}.csv"
        assert path.exists(), name
        bool_columns = []
        for column in frame:
            present = frame[column].dropna()
            if pd.api.types.is_bool_dtype(frame[column]) or (len(present) and present.map(
                    lambda value: isinstance(value, (bool, np.bool_))).all()):
                bool_columns.append(column)
        expected = frame.astype({column: "boolean" for column in bool_columns})
        saved_table = pd.read_csv(path, float_precision="round_trip",
                                  dtype={column: "boolean" for column in bool_columns})
        pd.testing.assert_frame_equal(expected.reset_index(drop=True), saved_table.reset_index(drop=True),
                                      check_dtype=False, check_exact=True)
        checked_tables.append(name)
    return {"state": "passed", "checked_utc": datetime.now(timezone.utc).isoformat(),
            "cases": len(jobs), "directional_rows": len(aggregate), "native_archives": native_count,
            "inputs_exactly_matching_preflight": reference_count,
            "pairwise_native_archives_copied_exactly": copied_count,
            "retained_scientific_failures": failures, "failed_method_cases": len(failures),
            "directional_status_counts": aggregate.status.value_counts(dropna=False).to_dict(),
            "ccm_arrays_checked": ccm_count, "gc_scientific_fit_records_checked": gc_fit_count,
            "rebuilt_summary_tables": checked_tables, "scientific_calculations_rerun": False}


def main():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output = ROOT / "data/benchmarks/common_driver_checks" / stamp
    output.mkdir(parents=True, exist_ok=False)
    prepared = output / "prepared"
    report_path = output / "report.json"
    started = time.perf_counter()
    report = {"state": "running", "started_utc": stamp,
              "output_directory": str(output), "production_run_started": False}

    def save_report():
        report["elapsed_seconds"] = time.perf_counter() - started
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"Common-driver verification: {report_path}", flush=True)
    save_report()
    try:
        jobs, process_diagnostics, clock_diagnostics = core.prepare_inputs(prepared)
        report["inputs"] = check_inputs(jobs, prepared)
        report["retention_and_resume"] = check_guards(jobs, prepared, output)
        report["pairwise_reuse"] = check_pairwise_reuse(jobs, prepared, output)
        report["summaries"] = check_summaries(jobs)
        save_report()
        report["native"] = check_native(jobs, prepared, output, report, save_report)
        report["validated_calculation_digest"] = core.calculation_digest()
        report["validated_implementation_sha256"] = core.execution_parameters()["implementation_sha256"]
        report["state"] = "passed"
    except BaseException as error:
        report["state"] = "failed"
        report["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        save_report()
    print(f"Common-driver checks passed: {report_path}", flush=True)


if __name__ == "__main__":
    main()
