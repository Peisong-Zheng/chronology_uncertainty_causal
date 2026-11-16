"""Scientific preflight and read-only artifact audit for the TBK excitation test.

The default execution checks isolated inputs and the new ParCorr call. It never
launches the twelve-input production experiment or repeats baseline GPDC/CCM.
"""

import os
for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"]:
    os.environ[name] = "1"
os.environ.setdefault("MPLBACKEND", "Agg")

from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
from datetime import datetime, timezone
import gzip
import importlib
from io import StringIO
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
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits

from toolbox import causal
from toolbox.tbk_model import simulate_tbk

source = importlib.import_module("04_check_tbk_chronology")
BRANCH = "04_1_check_tbk_excitation"
SEEDS = [20261009, 20261301, 20261302]
SCENARIOS = {"baseline": .005, "noise_x10": .05,
             "noise_x100": .5, "random_response": .005}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=True)


def read_archive(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def restore_arrays(value):
    if isinstance(value, dict) and set(value) == {"dtype", "shape", "values"}:
        array = np.asarray(value["values"], dtype=value["dtype"])
        assert list(array.shape) == value["shape"]
        return array
    if isinstance(value, dict) and set(value) == {"columns", "index", "data"}:
        return pd.DataFrame({key: restore_arrays(value["data"][str(key)])
                             for key in value["columns"]},
                            index=restore_arrays(value["index"]))
    if isinstance(value, dict):
        return {key: restore_arrays(item) for key, item in value.items()}
    if isinstance(value, list):
        return [restore_arrays(item) for item in value]
    return value


def snapshot_previous(path):
    """Record the old TBK scientific evidence without snapshotting unrelated runs."""
    files = {ROOT / name for name in [
        "04_check_tbk_chronology.py", "04_0_inspect_tbk_records.ipynb",
        "toolbox/tbk_model.py", "toolbox/tbk_observation.py", "toolbox/causal.py",
        "toolbox/observation.py",
        "data/archives/tbk_original_20261110T214805311811Z/One_D_model_TBK.ipynb"]}
    for folder in ["data/processed/04_check_tbk_chronology",
                   "data/processed/04_0_inspect_tbk_records",
                   "figures/04_check_tbk_chronology", "figures/04_0_inspect_tbk_records"]:
        files.update(item for item in (ROOT / folder).rglob("*") if item.is_file())
    records = []
    for item in sorted(files):
        if not item.is_file():
            raise FileNotFoundError(item)
        stat = item.stat()
        records.append({"path": item.relative_to(ROOT).as_posix(),
                        "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
                        "sha256": source.file_digest(item)})
    source.write_json({"created_utc": datetime.now(timezone.utc).isoformat(),
                       "files": records}, path)
    return {"files": len(records), "snapshot": str(Path(path).relative_to(ROOT))}


def check_preservation(path):
    snapshot = json.loads(Path(path).read_text(encoding="utf-8"))
    for item in snapshot["files"]:
        file = ROOT / item["path"]
        stat = file.stat()
        assert stat.st_size == item["size"], file
        assert stat.st_mtime_ns == item["mtime_ns"], file
        assert source.file_digest(file) == item["sha256"], file
    return {"passed": True, "files": len(snapshot["files"]),
            "checks": ["sha256", "size", "mtime_ns"]}


def physical_operator():
    """Independently write the annual four-temperature transition and q loading."""
    diffusion = np.array([[-.005, .005, 0., 0.],
                          [.005, -.105, .1, 0.],
                          [0., .1, -.2, .1],
                          [0., 0., .1, -.1]])
    transition = np.eye(4) + diffusion - np.ones((4, 4)) / 400.
    transition[1:3, 1:3] = transition[1:3, 1:3] - 1. / 200000.
    loading = np.array([0., -1., .5, .5])
    return transition, loading


def check_inputs(jobs, output_dir):
    """Reintegrate every condition and verify exact ages and paired innovations."""
    output_dir = Path(output_dir)
    cases = pd.DataFrame(jobs)
    assert len(cases) == cases.case_id.nunique() == 12
    assert cases.groupby("scenario").size().to_dict() == dict.fromkeys(SCENARIOS, 3)
    assert cases.groupby("process_seed").size().to_dict() == dict.fromkeys(SEEDS, 4)
    assert cases.input_status.eq("ok").all()
    assert cases.level.eq(0).all() and cases.sigma_years.eq(0).all()
    assert cases.member.eq(-1).all()
    np.testing.assert_array_equal(cases.estimator_seed, cases.process_seed + 400)
    transition, loading = physical_operator()
    analysis_time = np.arange(10000, 38001, 20)
    grid = (39000 - analysis_time) / 1000.
    indices = analysis_time - 1
    report = []
    for seed in SEEDS:
        saved = np.load(output_dir / "prepared" / f"history_{seed}.npz")
        clean_model = simulate_tbk(runlength=39000, south_diffusion=.005,
                                   noise_sd=0., seed=seed)
        clean = clean_model["temperature"][0, 1:5]
        np.testing.assert_array_equal(saved["clean_temperature"], clean)
        np.testing.assert_array_equal(saved["time"], np.arange(1, 39001))
        np.testing.assert_array_equal(saved["age"], (39000 - np.arange(1, 39001)) / 1000.)
        np.testing.assert_array_equal(saved["sampled_time"], np.arange(9000, 39001, 20))
        np.testing.assert_array_equal(saved["age_true"], (39000 - np.arange(9000, 39001, 20)) / 1000.)
        np.testing.assert_array_equal(saved["transport"], clean_model["transport_actual"])
        baseline = None
        baseline_noise = None
        for scenario, noise_sd in SCENARIOS.items():
            job = next(job for job in jobs if job["process_seed"] == seed
                       and job["scenario"] == scenario)
            assert job["noise_sd"] == noise_sd
            if scenario != "random_response":
                model = simulate_tbk(runlength=39000, south_diffusion=.005,
                                     noise_sd=noise_sd, seed=seed)
                temperature = model["temperature"][0, 1:5]
                np.testing.assert_array_equal(saved[f"{scenario}_temperature"], temperature)
                np.testing.assert_array_equal(model["transport_actual"], clean_model["transport_actual"])
                np.testing.assert_array_equal(model["event_intervals"], clean_model["event_intervals"])
                assert model["parameters"]["diffusion_coefficients"] == [.005, .005, .1, .1, .1]
                independent_update = (transition @ temperature[:, :-1]
                    + loading[:, None] * model["transport_actual"][None, :-1]
                    + model["noise"][0, 1:5, 1:])
                update_error = float(np.max(np.abs(temperature[:, 1:] - independent_update)))
                np.testing.assert_allclose(temperature[:, 1:], independent_update,
                                           rtol=0, atol=2e-13)
                if scenario == "baseline":
                    baseline = temperature.copy()
                    baseline_noise = model["noise"].copy()
                factor = noise_sd / .005
                innovation_error = float(np.max(np.abs(model["noise"] - factor * baseline_noise)))
                np.testing.assert_allclose(model["noise"], factor * baseline_noise,
                                           rtol=5e-16, atol=5e-16)
                response_error = float(np.max(np.abs(temperature - clean - factor * (baseline - clean))))
                # Integrating and subtracting forced trajectories changes floating-point
                # rounding; the underlying linear operator and innovations are fixed.
                np.testing.assert_allclose(temperature - clean, factor * (baseline - clean),
                                           rtol=0, atol=2e-10)
            else:
                temperature = baseline - clean
                np.testing.assert_array_equal(saved["random_response_temperature"], temperature)
                update_error = float(np.max(np.abs(temperature[:, 1:]
                    - transition @ temperature[:, :-1] - baseline_noise[0, 1:5, 1:])))
                np.testing.assert_allclose(temperature[:, 1:],
                    transition @ temperature[:, :-1] + baseline_noise[0, 1:5, 1:],
                    rtol=0, atol=3e-13)
                innovation_error, response_error = 0., 0.
            expected = pd.DataFrame({"time": grid[0] - grid, "age": grid,
                                     "x": temperature[0, indices], "y": temperature[3, indices]})
            path = output_dir / "inputs" / f"{job['case_id']}.csv"
            actual = pd.read_csv(path)
            assert len(actual) == 1401 and np.isfinite(actual.to_numpy()).all()
            assert np.all(np.diff(actual.time) > 0) and np.all(np.diff(actual.age) < 0)
            # Match the estimator's established default CSV parser, including its
            # last-bit conversion, rather than silently changing input precision.
            expected_parsed = pd.read_csv(StringIO(expected.to_csv(index=False)))
            pd.testing.assert_frame_equal(actual, expected_parsed, check_exact=True)
            assert source.file_digest(path) == job["input_sha256"]
            if scenario == "baseline":
                original = source.DATA_DIR / "inputs" / path.name
                assert source.file_digest(path) == source.file_digest(original)
                pd.testing.assert_frame_equal(actual, pd.read_csv(original), check_exact=True)
            report.append({"process_seed": seed, "scenario": scenario,
                           "annual_recurrence_max_abs_error": update_error,
                           "paired_innovation_max_abs_error": innovation_error,
                           "scaled_response_max_abs_error": response_error})
        saved.close()
    return {"cases": 12, "points_per_case": 1401,
            "all_saved_inputs_match_independent_reintegration": True,
            "baseline_input_bytes_exact": True, "transition_matrix": transition.tolist(),
            "forcing_loading": loading.tolist(), "per_input": report}


def check_baseline_copies(jobs, output_dir):
    core = importlib.import_module(BRANCH)
    source_jobs = pd.read_csv(source.DATA_DIR / "input_cases.csv").to_dict("records")
    count = 0
    for job in jobs:
        if job["scenario"] != "baseline":
            continue
        original = next(item for item in source_jobs if item["case_id"] == job["case_id"])
        assert original["estimator_seed"] == job["estimator_seed"]
        old_edges = source.read_completed_case(original, source.DATA_DIR)
        new_edges = core.read_completed_case(job, output_dir)
        assert new_edges is not None
        for method in source.METHODS:
            filename = f"{job['case_id']}_{method.replace('+', 'plus')}.json.gz"
            old = source.DATA_DIR / "native" / filename
            new = Path(output_dir) / "native" / filename
            assert source.file_digest(old) == source.file_digest(new), filename
            source.check_native_archive(read_archive(new), method, original)
            count = count + 1
        pd.testing.assert_frame_equal(new_edges[old_edges.columns], old_edges,
                                      check_dtype=False, check_exact=True)
    return {"baseline_native_archives_byte_identical": count,
            "all_baseline_rows_preserved": True}


def spawned_parcorr(series, seed):
    core = importlib.import_module(BRANCH)
    return core.run_parcorr(series, seed)


def check_parcorr(job, output_dir):
    """Compare the wrapper with an independently configured native Tigramite call."""
    from tigramite import data_processing
    from tigramite.independence_tests.parcorr import ParCorr
    from tigramite.pcmci import PCMCI
    core = importlib.import_module(BRANCH)
    series = pd.read_csv(Path(output_dir) / "inputs" / f"{job['case_id']}.csv")
    seed = job["estimator_seed"]
    frame = data_processing.DataFrame(series[["x", "y"]].to_numpy(), var_names=["x", "y"])
    native = PCMCI(dataframe=frame,
                   cond_ind_test=ParCorr(significance="analytic", seed=seed), verbosity=0)
    reference = native.run_pcmciplus(tau_min=0, tau_max=5, pc_alpha=.05, fdr_method="none")
    wrapped = core.run_parcorr(series, seed)
    for key, value in reference.items():
        assert canonical(source.native_snapshot(value)) == canonical(source.native_snapshot(wrapped["raw"][key])), key
    expected_holm = np.full((2, 2, 6), np.nan)
    expected_graph = reference["graph"].copy()
    for source_index, target_index in [(0, 1), (1, 0)]:
        corrected = multipletests(reference["p_matrix"][source_index, target_index],
                                  alpha=.05, method="holm")[1]
        expected_holm[source_index, target_index] = corrected
        expected_graph[source_index, target_index, corrected >= .05] = ""
        edge = wrapped["edges"].loc[wrapped["edges"].source.eq(["x", "y"][source_index])].iloc[0]
        assert edge.edge_detected == bool(np.any((reference["graph"][source_index, target_index] == "-->") & (corrected < .05)))
        assert edge.p_value == float(corrected.min())
        assert edge.score == float(np.max(np.abs(reference["val_matrix"][source_index, target_index])))
    np.testing.assert_array_equal(wrapped["raw"]["p_holm"], expected_holm)
    np.testing.assert_array_equal(wrapped["raw"]["graph_holm"], expected_graph)
    zero_links = [expected_graph[0, 1, 0], expected_graph[1, 0, 0]]
    assert wrapped["diagnostics"]["undirected_contemporaneous"] == ("o-o" in zero_links)
    assert wrapped["diagnostics"]["conflicting_contemporaneous"] == ("x-x" in zero_links)
    with ProcessPoolExecutor(max_workers=1, mp_context=multiprocessing.get_context("spawn")) as pool:
        spawned = pool.submit(spawned_parcorr, series, seed).result()
    assert canonical(source.native_snapshot(spawned)) == canonical(source.native_snapshot(wrapped))
    core.run_parcorr_case(job, output_dir)
    assert core.read_completed_parcorr(job, output_dir) is not None
    restored = core.load_native(job["case_id"], "PCMCI+_ParCorr", output_dir)
    for key in ["raw", "edges", "diagnostics"]:
        assert canonical(source.native_snapshot(restored[key])) == canonical(source.native_snapshot(wrapped[key]))
    return {"all_native_fields_exact": list(reference), "holm_and_directional_votes_exact": True,
            "unresolved_zero_lag_preserved": True, "windows_spawn_exact": True,
            "saved_native_parcorr_exact": True}


def check_resume_and_failure(jobs, output_dir):
    """Reject changed science/input and retain a failed native call explicitly."""
    core = importlib.import_module(BRANCH)
    output_dir = Path(output_dir)
    parameters = core.execution_parameters()
    source.write_json(parameters, output_dir / "parameters.json")
    core.check_parameters(parameters, output_dir)
    changed = deepcopy(parameters)
    changed["CCM_LIBRARY_SIZE"] = 301
    try:
        core.check_parameters(changed, output_dir)
    except ValueError:
        pass
    else:
        raise AssertionError("Resume accepted changed CCM library length.")
    baseline = next(job for job in jobs if job["scenario"] == "baseline")
    try:
        core.read_completed_case(dict(baseline, input_sha256="changed"), output_dir)
    except ValueError:
        pass
    else:
        raise AssertionError("Resume accepted changed estimator input.")
    archive = read_archive(core.native_path(baseline["case_id"], "CCM", output_dir))
    archive["raw"]["lag_rho"]["shape"] = [100, 31, 2]
    try:
        source.check_native_archive(archive, "CCM", core.source_job(baseline["case_id"]))
    except ValueError:
        pass
    else:
        raise AssertionError("Resume accepted an incomplete lag array.")
    failed_job = next(job for job in jobs if job["scenario"] == "noise_x10")
    failed_dir = output_dir.parent / "failure_check"
    (failed_dir / "inputs").mkdir(parents=True, exist_ok=True)
    shutil.copy2(output_dir / "inputs" / f"{failed_job['case_id']}.csv",
                 failed_dir / "inputs" / f"{failed_job['case_id']}.csv")
    with patch.object(core.core, "call_method", side_effect=RuntimeError("Controlled verification failure")):
        rows = core.run_method_job(failed_job, "GC", failed_dir, None)
    assert len(rows) == 2 and rows.status.eq("error").all()
    assert rows.edge_detected.isna().all() and rows.score.isna().all()
    failed_archive = core.validate_native(failed_job, "GC", failed_dir)
    assert "Controlled verification failure" in failed_archive["diagnostics"]["error"]
    assert core.read_completed_case(failed_job, failed_dir) is None
    return {"changed_scientific_setting_rejected": True, "changed_input_rejected": True,
            "incomplete_ccm_array_rejected": True, "failure_native_record_retained": True,
            "partial_case_not_treated_as_complete": True}


def check_native_calculations(archive, method):
    """Check derived numerical evidence against complete stored native arrays."""
    raw, edges = archive["raw"], archive["edges"]
    if method in ["PCMCI+", "PCMCI+_ParCorr"]:
        assert raw["p_matrix"].shape == raw["val_matrix"].shape == raw["graph"].shape == (2, 2, 6)
        for source_index, target_index in [(0, 1), (1, 0)]:
            row = edges.loc[edges.source.eq(["x", "y"][source_index])].iloc[0]
            if row.status != "ok":
                continue
            adjusted = multipletests(raw["p_matrix"][source_index, target_index], alpha=.05, method="holm")[1]
            np.testing.assert_array_equal(adjusted, raw["p_holm"][source_index, target_index])
            assert row.p_value == float(np.min(adjusted))
            retained = (raw["graph"][source_index, target_index] == "-->") & (adjusted < .05)
            assert row.edge_detected == bool(np.any(retained))
            scores = raw["val_matrix"][source_index, target_index]
            expected_score = np.max(np.abs(scores)) if method == "PCMCI+_ParCorr" else np.max(scores)
            assert row.score == float(expected_score)
    if method == "GC":
        for direction in ["x_to_y", "y_to_x"]:
            for fit_kind in ["restricted", "full"]:
                fit = raw[direction][fit_kind]
                fitted = np.asarray(fit["exog"]) @ np.asarray(fit["params"])
                np.testing.assert_allclose(fitted, fit["fittedvalues"], rtol=1e-13, atol=1e-13)
                np.testing.assert_allclose(np.asarray(fit["endog"]) - fitted, fit["resid"], rtol=1e-13, atol=1e-13)
                np.testing.assert_allclose(np.dot(fit["resid"], fit["resid"]), fit["ssr"], rtol=1e-13, atol=1e-13)
                assert fit["nobs"] == 1400
        assert raw["var_residuals"].shape == (1400, 2)
    if method == "CCM":
        for name, shape in {"lag_rho": (100, 61, 2), "rho_null": (100, 2),
                            "cuts": (100, 2), "library_seeds": (101, 2),
                            "lag_library_seeds": (100,)}.items():
            assert raw[name].shape == shape
        np.testing.assert_array_equal(raw["lag_steps"], np.arange(-30, 31))
        np.testing.assert_array_equal(raw["lag_common_anchor_indices"], np.arange(30, 1371))
        for index, label in enumerate(["x", "y"]):
            row = edges.loc[edges.source.eq(label)].iloc[0]
            if row.status != "ok":
                continue
            p_value = (1 + np.count_nonzero(raw["rho_null"][:, index] >= row.score)) / 101
            assert row.p_value == p_value
            assert row.raw_edge_detected == (p_value < .05)


def compare_saved_table(expected, path):
    """Compare complete CSV payloads, respecting blank-string CSV conversion."""
    boolean_types = {column: "boolean" for column in expected if pd.api.types.is_bool_dtype(expected[column])}
    expected_csv = pd.read_csv(StringIO(expected.to_csv(index=False)),
                               float_precision="round_trip", dtype=boolean_types)
    actual = pd.read_csv(path, float_precision="round_trip", dtype=boolean_types)
    pd.testing.assert_frame_equal(expected_csv, actual, check_dtype=False, check_exact=True)


def check_artifacts(output_dir=None):
    """Audit every committed native result without fitting or overwriting anything."""
    core = importlib.import_module(BRANCH)
    output_dir = core.DATA_DIR if output_dir is None else Path(output_dir)
    core.check_parameters(core.execution_parameters(), output_dir)
    jobs = pd.read_csv(output_dir / "input_cases.csv").to_dict("records")
    completed = {job["case_id"]: core.read_completed_case(job, output_dir) for job in jobs}
    assert len(completed) == 12 and all(table is not None for table in completed.values())
    baseline_jobs = [job for job in jobs if job["scenario"] == "baseline"]
    parcorr_tables = [core.read_completed_parcorr(job, output_dir) for job in baseline_jobs]
    assert all(table is not None for table in parcorr_tables)
    results = core.collect_results(jobs, completed)
    parcorr = core.collect_results(baseline_jobs,
        {job["case_id"]: table for job, table in zip(baseline_jobs, parcorr_tables)}, [core.PARCORR_METHOD])
    assert len(results) == 96 and len(parcorr) == 6
    assert not results.status.eq("pending").any() and not parcorr.status.eq("pending").any()
    compare_saved_table(results, output_dir / "member_results.csv")
    compare_saved_table(parcorr, output_dir / "parcorr_results.csv")
    copied = check_baseline_copies(jobs, output_dir)
    successful_native = 0
    for job in jobs:
        for method in source.METHODS:
            archive = core.load_native(job["case_id"], method, output_dir)
            if archive["diagnostics"].get("error") is not None:
                continue
            successful_native = successful_native + 1
            check_native_calculations(archive, method)
            if method == "CCM":
                rebuilt = causal.summarize_ccm_lag(archive, sign_threshold=.9)
                for key in ["lag_profile", "lag_peaks"]:
                    pd.testing.assert_frame_equal(rebuilt["raw"][key], archive["raw"][key], check_exact=True)
                pd.testing.assert_frame_equal(rebuilt["edges"], archive["edges"], check_exact=True)
                for key, suffix in [("lag_profile", "profile"), ("lag_peaks", "peaks")]:
                    disk = pd.read_csv(output_dir / "lag" / f"{job['case_id']}_{suffix}.csv", float_precision="round_trip")
                    pd.testing.assert_frame_equal(disk, archive["raw"][key], check_dtype=False, check_exact=True)
    for job in baseline_jobs:
        archive = core.load_native(job["case_id"], core.PARCORR_METHOD, output_dir)
        if archive["diagnostics"].get("error") is None:
            check_native_calculations(archive, core.PARCORR_METHOD)
    tables = core.summarize_results(results, parcorr, output_dir)
    for name, expected in tables.items():
        compare_saved_table(expected, output_dir / f"{name}.csv")
    plot_tables = core._plot_build_diagnostic_tables(output_dir)
    for name, expected in plot_tables.items():
        compare_saved_table(expected, output_dir / f"{name}.csv")
    assert len(list((output_dir / "native").glob("*.json.gz"))) == 51
    return {"passed": True, "cases": 12, "main_directional_rows": len(results),
            "parcorr_directional_rows": len(parcorr), "native_archives": 51,
            "successful_main_native_archives": successful_native,
            "failed_directional_rows": int(results.status.ne("ok").sum()),
            "failed_parcorr_directional_rows": int(parcorr.status.ne("ok").sum()),
            "unresolved_directional_rows": int(results.unoriented.fillna(False).sum()),
            "summary_tables_exact": list(tables), "figure_tables_exact": list(plot_tables),
            "ccm_lag_summaries_rebuilt_exactly": True,
            "pcmci_holm_and_votes_recomputed": True,
            "gc_full_regression_fits_and_residuals_checked": True,
            **copied}


def main(output_folder=None):
    core = importlib.import_module(BRANCH)
    started = time.perf_counter()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    folder = ROOT / "tmp/tbk_excitation_checks" / stamp if output_folder is None else Path(output_folder)
    folder.mkdir(parents=True, exist_ok=True)
    manifest = folder / "prior_tbk_manifest.json"
    report = {"passed": False, "started_utc": datetime.now(timezone.utc).isoformat()}
    try:
        report["snapshot"] = snapshot_previous(manifest)
        with threadpool_limits(limits=1):
            output_dir = folder / "prepared_run"
            jobs, _ = core.prepare_inputs(output_dir)
            report["inputs"] = check_inputs(jobs, output_dir)
            print("All twelve inputs and paired annual dynamics verified.", flush=True)
            core.reuse_baselines(jobs, output_dir)
            report["baseline_copies"] = check_baseline_copies(jobs, output_dir)
            baseline = next(job for job in jobs if job["scenario"] == "baseline")
            report["parcorr"] = check_parcorr(baseline, output_dir)
            report["resume_and_failure"] = check_resume_and_failure(jobs, output_dir)
        report["preservation"] = check_preservation(manifest)
        report["passed"] = True
    finally:
        report["elapsed_seconds"] = time.perf_counter() - started
        source.write_json(report, folder / "report.json")
    print(f"PASS: {folder / 'report.json'}", flush=True)
    return report


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
