"""Check the memory branch's paired inputs, result retention and native outputs.

This checks all 432 prepared inputs, then fits two complete-setting cases both
directly and with two Windows-spawned workers. All new outputs stay under a
timestamped data/benchmarks/memory_checks directory. No production run starts.
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
import re
import shutil
import sys
import time
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

# A normal module import lets Windows spawn resolve the same worker function.
core = importlib.import_module("02_8_check_memory_sensitivity")


def file_state(path):
    """Check both content and modification time when a completed case resumes."""
    return {"sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "mtime_ns": path.stat().st_mtime_ns}


def canonical(payload):
    return json.dumps(payload, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True, allow_nan=True).encode("utf-8")


def expect_rejection(function, *args):
    try:
        function(*args)
    except (ValueError, RuntimeError):
        return
    raise AssertionError(f"Changed inputs/settings were accepted by {function.__name__}.")


def check_inputs(jobs, prepared):
    """Check the physical design and compare the reference to retained inputs."""
    cases = pd.DataFrame(jobs)
    assert len(cases) == 432 and cases.case_id.nunique() == 432
    assert cases.groupby(["memory", "process_seed", "process"]).size().eq(12).all()
    assert cases.groupby("condition").size().to_dict() == {
        "chronology": 360, "irregular_accurate": 36, "original_uniform": 36}
    reference_comparisons = []
    grid = np.linspace(29.0, 1.0, 1401)
    for job in jobs:
        path = prepared / "inputs" / f"{job['case_id']}.csv"
        # The established estimator input path uses pandas' ordinary parser.
        series = pd.read_csv(path)
        assert list(series) == ["time", "age", "x", "y"]
        assert series.shape == (1401, 4)
        assert np.isfinite(series.to_numpy()).all()
        np.testing.assert_allclose(series.age, grid, rtol=0, atol=1e-13)
        np.testing.assert_allclose(series.time, grid[0] - grid, rtol=0, atol=1e-13)
        if job["memory"] != "reference" or job["process_seed"] != 20261009:
            continue
        process = job["process"]
        if job["condition"] == "chronology":
            original = ROOT / "data/processed/02_3_check_uncertainty_escalation" / (
                f"r0_{process}_level4_C-C_m{job['member']}_series.csv")
        else:
            original = ROOT / "data/processed/02_check_causal_pipeline" / (
                f"r0_{process}_{job['condition']}_series.csv")
        expected = pd.read_csv(original)
        pd.testing.assert_frame_equal(series, expected, check_exact=True)
        reference_comparisons.append(original.relative_to(ROOT).as_posix())
    assert len(reference_comparisons) == 36

    # Reconstruct innovations from the saved full processes, rather than trusting
    # a reported pairing flag. Standardized disturbances must match across memory.
    recovered = {}
    for row in cases[["memory", "process_seed", "process", "memory_x_years", "memory_y_years"]].drop_duplicates().to_dict("records"):
        name = f"s{row['process_seed']}_{row['memory']}_{row['process']}"
        series = pd.read_csv(prepared / "generated" / f"{name}.csv", float_precision="round_trip")
        assert len(series) == 1501
        ax = np.exp(-20 / row["memory_x_years"])
        ay = np.exp(-20 / row["memory_y_years"])
        sx, sy = 1 / np.sqrt(1 - 0.95**2), 1 / np.sqrt(1 - 0.8**2)
        sigma_u, sigma_v = sx * np.sqrt(1 - ax**2), sy * np.sqrt(1 - ay**2)
        g = series.x.to_numpy() / sx
        if row["process"] == "nonlinear":
            nodes, weights = np.polynomial.hermite.hermgauss(256)
            scale = np.sqrt(np.sum(weights * np.tanh(2 * np.sqrt(2) * nodes)**2) / np.sqrt(np.pi))
            g = np.tanh(2 * g) / scale
        c = 0 if row["process"] == "independent" else (1 - ay) * 4
        x, y = series.x.to_numpy(), series.y.to_numpy()
        innovations = np.column_stack(((x[1:] - ax * x[:-1]) / sigma_u,
                                       (y[1:] - ay * y[:-1] - c * g[:-1]) / sigma_v))
        seed = row["process_seed"]
        if seed in recovered:
            np.testing.assert_allclose(innovations, recovered[seed], rtol=0, atol=2e-14)
        else:
            recovered[seed] = innovations
        # The common 125-step warm-up plus one lag leaves these exact positions
        # in the seeded normal streams. A memory-dependent warm-up would fail.
        random = np.random.default_rng(seed)
        n_total = 125 + 1 + 1501
        expected_u = random.normal(size=n_total)[127:]
        expected_v = random.normal(size=n_total)[127:]
        np.testing.assert_allclose(innovations, np.column_stack((expected_u, expected_v)), rtol=0, atol=2e-14)

    # The clocks are an external fixed input, not regenerated for each memory.
    for record in ["x", "y"]:
        path = ROOT / "data/processed/02_3_check_uncertainty_escalation" / f"r0_level4_{record}_C_chronology.npz"
        with np.load(path, allow_pickle=False) as clock:
            assert clock["ages"].shape[0] == 10
            assert np.all(np.diff(clock["ages"], axis=1) < 0)
            assert np.all(clock["ages"][:, 0] >= 29) and np.all(clock["ages"][:, -1] <= 1)
            positions = np.rint((30 - clock["age_true"]) / 0.02).astype(int)
            np.testing.assert_allclose(30 - positions * 0.02, clock["age_true"], rtol=0, atol=1e-12)
            assert np.all(np.diff(positions) > 0)
    return {"prepared_pairs": len(jobs), "reference_inputs_exact": reference_comparisons,
            "paired_standardized_innovations_and_fixed_warmup": True,
            "fixed_grid_and_ordered_source_clocks": True}


def prepare_null(output, seed):
    path = output / "cache" / f"null_T1391_n10000_seed{seed}.npz"
    with threadpool_limits(limits=1):
        core.prepare_gpdc_null(path, sample_sizes=[1391], sig_samples=10000, seed=seed)
    return path


def fake_method(series, method, estimator_seed, null_file):
    """Return explicit mixed outcomes to test storage, not scientific estimates."""
    if method == "LKIF":
        raise RuntimeError("Intentional verification-only estimator failure.")
    edges = pd.DataFrame({"source": ["x", "y"], "target": ["y", "x"],
                          "score": [0.3, 0.2], "p_value": [0.01, 0.8],
                          "edge_detected": [True, False], "status": ["ok", "ok"]})
    diagnostics = {"conflicting_contemporaneous": method == "PCMCI+"}
    if method == "CCM":
        edges["raw_edge_detected"] = [True, False]
        edges["edge_detected"] = [None, False]
        edges["lag_status"] = ["failed", "unresolved"]
        edges["unoriented"] = False
    return {"edges": edges, "diagnostics": diagnostics,
            "raw": {"array": np.arange(6.).reshape(3, 2)}}


def check_guards(jobs, prepared, output):
    """Retain scientific failures but reject partial or scientifically changed runs."""
    job = jobs[0]
    target = output / "failure_retention"
    (target / "inputs").mkdir(parents=True)
    source = prepared / "inputs" / f"{job['case_id']}.csv"
    shutil.copyfile(source, target / "inputs" / source.name)
    with patch.object(core, "call_method", side_effect=fake_method):
        result = core.run_case(job, target, "unused_mock_null.npz")
    assert len(result) == 8
    assert result.loc[result.method.eq("LKIF"), "status"].eq("error").all()
    assert result.loc[result.method.eq("LKIF"), "edge_detected"].isna().all()
    assert result.loc[result.method.eq("PCMCI+"), "unoriented"].all()
    resumed = core.read_completed_case(job, target)
    assert resumed is not None
    assert resumed.loc[resumed.method.eq("LKIF"), "status"].eq("error").all()
    assert resumed.loc[resumed.method.eq("CCM") & resumed.source.eq("x"), "score"].iloc[0] == 0.3
    assert resumed.loc[resumed.method.eq("CCM") & resumed.source.eq("x"), "edge_detected"].isna().all()
    checks = {"scientific_failure_retained": True,
              "pcmci_unresolved_evidence_retained": True,
              "ccm_score_retained_after_failed_screen": True}

    for label, shortened, omit in [("short_commit", True, None),
                                    ("missing_native", False, "_GC.json.gz")]:
        folder = output / label
        (folder / "members").mkdir(parents=True)
        (folder / "native").mkdir()
        table = result.iloc[:-1] if shortened else result
        table.to_csv(folder / "members" / f"{job['case_id']}_edges.csv", index=False)
        for path in (target / "native").glob("*.json.gz"):
            if omit is None or not path.name.endswith(omit):
                shutil.copyfile(path, folder / "native" / path.name)
        assert core.read_completed_case(job, folder) is None
        checks[f"{label}_not_resumable"] = True

    failed = output / "serialization_failure"
    (failed / "inputs").mkdir(parents=True)
    shutil.copyfile(source, failed / "inputs" / source.name)

    def unserializable(series, method, seed, null):
        payload = fake_method(series, "GC" if method == "LKIF" else method, seed, null)
        if method == "GC":
            payload["raw"]["unsupported_object"] = object()
        return payload

    with patch.object(core, "call_method", side_effect=unserializable):
        try:
            core.run_case(job, failed, "unused_mock_null.npz")
        except TypeError:
            pass
        else:
            raise AssertionError("Serialization errors must stop before the case commit.")
    assert not (failed / "members" / f"{job['case_id']}_edges.csv").exists()
    assert core.read_completed_case(job, failed) is None
    checks["serialization_failure_not_a_committed_vote"] = True

    parameter_dir = output / "parameter_guards"
    parameter_dir.mkdir()
    parameters = core.execution_parameters()
    (parameter_dir / "parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    core.check_parameters(parameters, parameter_dir)
    changed = deepcopy(parameters)
    changed["N_WORKERS"] = 16
    core.check_parameters(changed, parameter_dir)
    for key in ["ALPHA", "CCM_LIBRARY_SIZE", "BURN_IN_STEPS"]:
        if key in parameters:
            changed = deepcopy(parameters)
            changed[key] = parameters[key] + 1
            expect_rejection(core.check_parameters, changed, parameter_dir)
            checks[f"changed_{key}_rejected"] = True
    for key in ["packages", "environment", "source_sha256", "implementation_sha256"]:
        if key in parameters:
            changed = deepcopy(parameters)
            changed[key] = {"verification": "different scientific environment or source"}
            expect_rejection(core.check_parameters, changed, parameter_dir)
            checks[f"changed_{key}_rejected"] = True

    # Modify an isolated prepared input, confirm a regeneration mismatch is
    # rejected, then restore it before any numerical acceptance calculation.
    contents, stat = source.read_bytes(), source.stat()
    changed = pd.read_csv(source, float_precision="round_trip")
    changed.loc[1, "x"] = changed.loc[1, "x"] + 0.01
    changed.to_csv(source, index=False)
    try:
        expect_rejection(core.prepare_inputs, prepared)
    finally:
        source.write_bytes(contents)
        os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    checks["changed_prepared_input_rejected"] = True
    # Run controls and figures may change without redefining the experiment;
    # scientific settings must remain part of the implementation fingerprint.
    source_text = (ROOT / "02_8_check_memory_sensitivity.py").read_text(encoding="utf-8")
    original_digest = core.calculation_digest()
    control_file = output / "controls_only.py"
    control_file.write_text(re.sub(r"N_WORKERS = 32", "N_WORKERS = 16", source_text, count=1), encoding="utf-8")
    assert core.calculation_digest(control_file) == original_digest
    scientific_file = output / "scientific_change.py"
    scientific_file.write_text(re.sub(r"CCM_LIBRARY_SIZE = 300", "CCM_LIBRARY_SIZE = 301", source_text, count=1), encoding="utf-8")
    assert core.calculation_digest(scientific_file) != original_digest
    checks["implementation_guard_separates_controls_from_science"] = True
    return checks


def check_summaries(jobs):
    """Check requested denominators and score pairing with explicit mixed outcomes."""
    subset = [job for job in jobs if job["process_seed"] == 20261009
              and job["process"] == "linear" and job["memory"] in ["reference", "short_x"]]
    pending = core.summarize_results(core.collect_results(subset, {}))["scenario_summary"]
    assert pending.n_valid.eq(0).all() and pending.n_correct_graph.eq(0).all()
    assert pending.fraction_correct_requested.eq(0).all()
    assert pending.loc[pending.condition.eq("chronology"), "n_requested"].eq(10).all()
    rows = []
    for job in subset:
        for method in core.METHODS:
            for source, target in [("x", "y"), ("y", "x")]:
                chronology = job["condition"] == "chronology"
                score = (0.4 if source == "x" else 0.1) if chronology else (0.6 if source == "x" else 0.2)
                if chronology and job["memory"] == "short_x" and source == "x":
                    score = 0.2
                row = {**job, "method": method, "source": source, "target": target,
                       "score": score, "p_value": 0.01 if source == "x" else 0.5,
                       "edge_detected": source == "x", "status": "ok", "unoriented": False,
                       "raw_edge_detected": True if method == "CCM" else pd.NA,
                       "lag_status": "compatible" if method == "CCM" else ""}
                if chronology:
                    if job["member"] == 0:
                        row.update(status="error", score=np.nan, edge_detected=None)
                    elif job["member"] == 1:
                        row["unoriented"] = True
                    elif job["member"] == 2 and source == "y":
                        row["edge_detected"] = True
                    elif job["member"] == 3 and source == "x":
                        row["score"] = np.nan
                    elif job["member"] == 4 and method == "CCM" and source == "x":
                        row.update(edge_detected=None, lag_status="failed")
                rows.append(row)
    results = pd.DataFrame(rows)
    for name in ["edge_detected", "unoriented", "raw_edge_detected"]:
        results[name] = results[name].astype("boolean")
    tables = core.summarize_results(results)
    summaries = tables["scenario_summary"]
    for _, summary in summaries.loc[summaries.condition.eq("chronology")].iterrows():
        ccm = summary.method == "CCM"
        assert summary.n_requested == 10 and summary.n_unresolved == 1
        assert summary.n_valid == (8 if ccm else 9)
        assert summary.n_failed == (2 if ccm else 1)
        assert summary.n_correct_graph == (6 if ccm else 7)
        assert summary.fraction_correct_requested == (0.6 if ccm else 0.7)
        assert summary.q_x_to_y == 1
        assert summary.q_y_to_x == 1 / (8 if ccm else 9)
        assert summary.n_scores_x_to_y == 8 and summary.n_scores_y_to_x == 9
        assert summary.decision == ("insufficient_valid_members" if ccm else "x_to_y")
        if ccm:
            assert summary.n_raw_valid == 9
            assert summary.raw_q_x_to_y == 1 and summary.raw_q_y_to_x == 1
    changes = tables["paired_score_changes"]
    assert changes.loc[changes.source.eq("x")].groupby(["memory", "method"]).score_change_usable.sum().eq(8).all()
    # A successful raw score remains usable even if its CCM lag decision failed.
    failed_screen = changes.loc[changes.method.eq("CCM") & changes.member.eq(4) & changes.source.eq("x")]
    assert failed_screen.score_change_usable.all()
    contrasts = tables["memory_contrasts"]
    forward = contrasts.loc[contrasts.source.eq("x") & contrasts.contrast_usable]
    assert len(forward) == 4 * 8
    np.testing.assert_allclose(forward.difference_in_score_change, -0.2, rtol=0, atol=1e-15)
    # An unavailable reference must invalidate the paired change; neither zero
    # replacement nor one-sided finite filtering is permissible.
    broken = results.copy()
    baseline = broken.condition.eq("irregular_accurate") & broken.memory.eq("short_x") & broken.method.eq("GC") & broken.source.eq("x")
    broken.loc[baseline, "status"] = "error"
    broken.loc[baseline, "score"] = np.nan
    rebuilt = core.summarize_results(broken)
    affected = rebuilt["paired_score_changes"]
    affected = affected.loc[affected.memory.eq("short_x") & affected.method.eq("GC") & affected.source.eq("x")]
    assert not affected.score_change_usable.any() and affected.delta_score.isna().all()
    return {"requested_denominators_and_pending_states": True,
            "unresolved_and_failed_excluded_from_complete_graph": True,
            "q_uses_two_usable_directions": True,
            "finite_successful_scores_and_baselines_required": True,
            "ccm_native_score_survives_failed_lag_screen": True,
            "memory_contrasts_require_both_paired_changes": True}


def check_native(jobs, prepared, output, report, save_report):
    """Compare the complete native evidence directly and with two workers."""
    selected = [next(job for job in jobs if job["process_seed"] == 20261009
                     and job["process"] == process and job["memory"] == memory
                     and job["condition"] == "chronology" and job["member"] == 0)
                for process, memory in [("linear", "short_x"), ("nonlinear", "short_xy")]]
    null_file = prepare_null(output, selected[0]["estimator_seed"])
    parallel = output / "parallel"
    (parallel / "inputs").mkdir(parents=True)
    direct = {}
    for job in selected:
        source = prepared / "inputs" / f"{job['case_id']}.csv"
        shutil.copyfile(source, parallel / "inputs" / source.name)
        series = pd.read_csv(source)
        with threadpool_limits(limits=1):
            for method in core.METHODS:
                print(f"Direct native check: {job['case_id']} {method}", flush=True)
                result = core.call_method(series, method, job["estimator_seed"], str(null_file))
                assert result["edges"].status.eq("ok").all(), result["edges"].to_string()
                direct[job["case_id"], method] = core.native_snapshot(result)
        report.setdefault("direct_completed", []).append(job["case_id"])
        save_report()
    print("Checking the same inputs with two Windows-spawned workers.", flush=True)
    with ProcessPoolExecutor(max_workers=2, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(core.run_case, job, str(parallel), str(null_file)) for job in selected]
        for future in as_completed(futures):
            result = future.result()
            assert len(result) == 8 and result.status.eq("ok").all(), result.to_string()
            report.setdefault("spawn_completed", []).append(result.case_id.iloc[0])
            save_report()

    comparisons = []
    for job in selected:
        case_id = job["case_id"]
        before = {p.relative_to(parallel).as_posix(): file_state(p)
                  for p in parallel.rglob(f"{case_id}*") if p.is_file()}
        resumed = core.read_completed_case(job, parallel)
        assert resumed is not None and len(resumed) == 8
        after = {p.relative_to(parallel).as_posix(): file_state(p)
                 for p in parallel.rglob(f"{case_id}*") if p.is_file()}
        assert before == after
        for method in core.METHODS:
            path = parallel / "native" / f"{case_id}_{method.replace('+', 'plus')}.json.gz"
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                worker = json.load(handle)
            reference = direct[case_id, method]
            expected, actual = canonical(reference), canonical(worker)
            assert expected == actual, f"Complete native output changed: {case_id}, {method}"
            comparison = {"case_id": case_id, "method": method,
                          "exact_complete_payload_match": True,
                          "sha256": hashlib.sha256(actual).hexdigest()}
            if method == "CCM":
                assert worker["raw"]["lag_rho"]["shape"] == [100, 31, 2]
                assert worker["raw"]["rho_null"]["shape"] == [100, 2]
                comparison.update(lag_rho_shape=[100, 31, 2], rho_null_shape=[100, 2])
            comparisons.append(comparison)
    return {"complete_native_comparisons": comparisons,
            "recorded_fields_ignored": [], "resume_preserves_committed_files": True}


def frame_from_snapshot(snapshot):
    """Read the portable native DataFrame representation used by the branch."""
    return pd.DataFrame({name: snapshot["data"][name]["values"] for name in snapshot["columns"]})


def audit_outputs(output_dir):
    """Read all production artifacts without rerunning a method or writing data."""
    output = Path(output_dir)
    cases = pd.read_csv(output / "cases.csv", float_precision="round_trip")
    jobs = cases[core.IDENTITY_COLUMNS].to_dict("records")
    assert len(jobs) == 432 and len({job["case_id"] for job in jobs}) == 432
    for case in cases.to_dict("records"):
        path = output / "inputs" / f"{case['case_id']}.csv"
        assert hashlib.sha256(path.read_bytes()).hexdigest() == case["input_sha256"]
    records = []
    failures = []
    native_count = 0
    ccm_arrays_checked = 0
    gc_models_checked = 0
    for job in jobs:
        rows = core.read_completed_case(job, output)
        assert rows is not None and len(rows) == 8, job["case_id"]
        records.append(rows)
        for method in core.METHODS:
            path = output / "native" / f"{job['case_id']}_{method.replace('+', 'plus')}.json.gz"
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                payload = json.load(handle)
            native_count = native_count + 1
            assert set(payload) >= {"edges", "raw", "diagnostics"}, path.name
            native_edges = frame_from_snapshot(payload["edges"])
            selected = rows.loc[rows.method.eq(method)].set_index("source")
            for _, native_edge in native_edges.iterrows():
                committed = selected.loc[native_edge["source"]]
                for field in ["target", "status", "score", "p_value", "edge_detected",
                              "raw_edge_detected", "optimal_lag_steps", "lag_status"]:
                    if field not in native_edges:
                        continue
                    a, b = native_edge[field], committed[field]
                    assert (pd.isna(a) and pd.isna(b)) or a == b, (path.name, field, a, b)
            if not selected.status.eq("ok").all():
                failures.append({"case_id": job["case_id"], "method": method,
                                 "statuses": selected.status.tolist()})
                continue
            if method == "CCM":
                assert payload["raw"]["lag_rho"]["shape"] == [100, 31, 2]
                assert payload["raw"]["rho_null"]["shape"] == [100, 2]
                assert len(payload["raw"]["lag_library_seeds"]["values"]) == 100
                ccm_arrays_checked = ccm_arrays_checked + 1
            if method == "GC":
                for direction in ["x_to_y", "y_to_x"]:
                    for model in ["restricted", "full"]:
                        fit = payload["raw"][direction][model]
                        assert {"params", "resid", "fittedvalues", "cov_params",
                                "endog", "exog", "df_resid", "llf"}.issubset(fit)
                        gc_models_checked = gc_models_checked + 1
    aggregate = pd.concat(records, ignore_index=True)
    assert len(aggregate) == 3456 and native_count == 1728
    assert len(list((output / "members").glob("*_edges.csv"))) == 432
    assert len(list((output / "native").glob("*.json.gz"))) == 1728
    saved = pd.read_csv(output / "member_results.csv", float_precision="round_trip",
                        dtype={"edge_detected": "boolean", "unoriented": "boolean", "raw_edge_detected": "boolean"})
    pd.testing.assert_frame_equal(aggregate.reset_index(drop=True), saved.reset_index(drop=True),
                                  check_dtype=False, check_exact=True)
    tables = core.summarize_results(aggregate)
    checked_tables = []
    for name, frame in tables.items():
        path = output / f"{name}.csv"
        assert path.exists(), path.name
        # CSV stores nullable Boolean flags as blanks. Read those columns with
        # their explicit nullable type so pd.NA is not mistaken for changed data.
        bool_columns = []
        for column in frame:
            present = frame[column].dropna()
            if pd.api.types.is_bool_dtype(frame[column]) or (
                    len(present) and present.map(lambda value: isinstance(value, (bool, np.bool_))).all()):
                bool_columns.append(column)
        expected = frame.astype({column: "boolean" for column in bool_columns})
        saved_table = pd.read_csv(path, float_precision="round_trip",
                                  dtype={column: "boolean" for column in bool_columns})
        pd.testing.assert_frame_equal(expected.reset_index(drop=True), saved_table.reset_index(drop=True),
                                      check_dtype=False, check_exact=True)
        checked_tables.append(name)
    summaries = tables["scenario_summary"]
    assert summaries.loc[summaries.condition.eq("chronology"), "n_requested"].eq(10).all()
    assert summaries.loc[~summaries.condition.eq("chronology"), "n_requested"].eq(1).all()
    np.testing.assert_array_equal(summaries.fraction_correct_requested,
                                  summaries.n_correct_graph / summaries.n_requested)
    changes, contrasts = tables["paired_score_changes"], tables["memory_contrasts"]
    usable_by_method = []
    for method in core.METHODS:
        paired = changes.loc[changes.method.eq(method)]
        difference = contrasts.loc[contrasts.method.eq(method)]
        usable_by_method.append({"method": method, "paired_changes_requested": len(paired),
            "paired_changes_usable": int(paired.score_change_usable.sum()),
            "memory_contrasts_requested": len(difference),
            "memory_contrasts_usable": int(difference.contrast_usable.sum())})
    return {"state": "passed", "checked_utc": datetime.now(timezone.utc).isoformat(),
            "cases": 432, "directional_rows": 3456, "native_archives": 1728,
            "retained_scientific_failures": failures, "failed_method_cases": len(failures),
            "directional_status_counts": aggregate.status.value_counts(dropna=False).to_dict(),
            "ccm_arrays_checked": ccm_arrays_checked,
            "ccm_lag_rho_shape": [100, 31, 2], "ccm_rho_null_shape": [100, 2],
            "gc_scientific_models_checked": gc_models_checked,
            "rebuilt_summary_tables": checked_tables,
            "summary_comparison": "Exact values, with nullable Boolean CSV fields read explicitly",
            "paired_score_coverage": usable_by_method,
            "scientific_calculations_rerun": False}


def main():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output = ROOT / "data/benchmarks/memory_checks" / stamp
    output.mkdir(parents=True, exist_ok=False)
    prepared = output / "prepared"
    report_path = output / "report.json"
    started = time.perf_counter()
    report = {"state": "running", "started_utc": stamp,
              "output_directory": str(output), "production_run_started": False,
              "parameters": core.execution_parameters()}

    def save_report():
        report["elapsed_seconds"] = time.perf_counter() - started
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"Memory-branch verification: {report_path}", flush=True)
    save_report()
    try:
        jobs, process_diagnostics, clock_diagnostics = core.prepare_inputs(prepared)
        report["inputs"] = check_inputs(jobs, prepared)
        save_report()
        report["retention_and_resume"] = check_guards(jobs, prepared, output)
        report["summaries"] = check_summaries(jobs)
        save_report()
        report["native"] = check_native(jobs, prepared, output, report, save_report)
        report["state"] = "passed"
    except BaseException as error:
        report["state"] = "failed"
        report["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        save_report()
    print(f"Memory-branch checks passed: {report_path}", flush=True)


if __name__ == "__main__":
    main()
