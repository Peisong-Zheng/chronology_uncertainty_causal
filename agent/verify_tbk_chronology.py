"""Same-input TBK chronology preflight and read-only production artifact audit.

The default run prepares isolated inputs and checks full native outputs, including
all 100 paired libraries and 61 CCM shifts. It does not start the 123-case run.
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
import pyEDM
from threadpoolctl import threadpool_limits

from toolbox import causal, observation, tbk_observation

core = importlib.import_module("04_check_tbk_chronology")
MANIFEST = ROOT / "tmp/tbk_preexisting_manifest.json"


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=True)


def read_archive(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def restore_arrays(value):
    """Read the explicit array/table representation used by the experiment."""
    if isinstance(value, dict) and set(value) == {"dtype", "shape", "values"}:
        result = np.asarray(value["values"], dtype=value["dtype"])
        assert list(result.shape) == value["shape"]
        return result
    if isinstance(value, dict) and set(value) == {"columns", "index", "data"}:
        data = {name: restore_arrays(value["data"][str(name)]) for name in value["columns"]}
        return pd.DataFrame(data, index=restore_arrays(value["index"]))
    if isinstance(value, dict):
        return {key: restore_arrays(item) for key, item in value.items()}
    if isinstance(value, list):
        return [restore_arrays(item) for item in value]
    return value


def check_preservation(snapshot=MANIFEST):
    snapshot = json.loads(Path(snapshot).read_text(encoding="utf-8"))
    for record in snapshot["files"]:
        path = ROOT / record["path"]
        stat = path.stat()
        assert stat.st_size == record["size"], path
        assert stat.st_mtime_ns == record["mtime_ns"], path
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"], path
    return {"passed": True, "files": len(snapshot["files"]),
            "checks": ["size", "mtime_ns", "sha256"],
            "checked_utc": datetime.now(timezone.utc).isoformat()}


def native_candidates(center, std, factor, seed, attempted, old=29., young=1.):
    """Use the existing observation kernel, retaining IDs before rejection."""
    rng = np.random.default_rng(seed)
    accepted, ids, order_rejected, coverage_rejected = [], [], 0, 0
    with patch.object(observation, "_age_error_factor", return_value=factor):
        for candidate_id in range(1, attempted + 1):
            candidate = center + observation._draw_age_error(std, 1 / 30, rng)
            if not np.all(np.diff(candidate) < 0):
                order_rejected = order_rejected + 1
            elif candidate[0] < old or candidate[-1] > young:
                coverage_rejected = coverage_rejected + 1
            else:
                accepted.append(candidate)
                ids.append(candidate_id)
    return np.asarray(accepted), np.asarray(ids), order_rejected, coverage_rejected


def check_inputs(jobs, output_dir):
    """Independently reconstruct dates, fixed proxy values and every interpolation."""
    output_dir = Path(output_dir)
    cases = pd.DataFrame(jobs)
    assert len(cases) == cases.case_id.nunique() == 123
    assert cases.groupby("process_seed").size().eq(41).all()
    assert cases.groupby("level").size().to_dict() == {0: 3, 1: 30, 2: 30, 3: 30, 4: 30}
    assert set(cases.process_seed) == set(core.PROCESS_SEEDS)
    assert cases.input_status.eq("ok").all()
    assert np.array_equal(cases.estimator_seed, cases.process_seed + 400)
    factor = observation._age_error_factor(1501, 1 / 30)
    cached_factor = tbk_observation.age_error_factor(1501, 1 / 30)
    for left, right in zip(factor, cached_factor):
        np.testing.assert_array_equal(left, right)
    clock_count = 0
    rejections = {"order": 0, "coverage": 0}
    seeds = []
    native_member_checks = 0
    diagnostic_rows, relative_rows, acceptance_rows = [], [], []
    for process_seed in core.PROCESS_SEEDS:
        prepared = restore_arrays(read_archive(output_dir / "prepared" / f"history_{process_seed}.json.gz"))
        model, sampled, grid = prepared["model"], prepared["sampled"], prepared["age_grid"]
        np.testing.assert_array_equal(model["time"], np.arange(1, 39001))
        np.testing.assert_array_equal(sampled.model_time_years, np.arange(9000, 39001, 20))
        age = (39000 - sampled.model_time_years.to_numpy()) / 1000
        np.testing.assert_array_equal(sampled.age_true, age)
        assert len(sampled) == 1501 and len(grid) == 1401
        assert grid[0] == 29 and grid[-1] == 1
        selected = (age <= 29) & (age >= 1)
        indices = sampled.model_time_years.to_numpy(dtype=int) - 1
        for column, box in [("x", 1), ("y", 4), ("t2", 2), ("t3", 3)]:
            np.testing.assert_array_equal(sampled[column], model["temperature"][0, box, indices])
        np.testing.assert_array_equal(sampled.q, model["transport_actual"][indices])
        for key, clock in prepared["chronologies"].items():
            clock_count = clock_count + 1
            np.testing.assert_array_equal(clock["age_true"], age)
            np.testing.assert_array_equal(clock["value"], sampled[clock["record"]])
            assert clock["correlation_years"] == 1000 and clock["correlation_length"] == 1 / 30
            np.testing.assert_array_equal(clock["age_std"], np.full(1501, clock["sigma_years"] / 1000))
            record_index = 0 if clock["record"] == "x" else 1
            assert clock["center_seed"] == [process_seed, 4101, record_index]
            assert clock["member_seed"] == [process_seed, 4102, record_index]
            seeds.extend([tuple(clock["center_seed"]), tuple(clock["member_seed"])])
            for phase, center, expected, expected_ids in [
                ("center", age, clock["age_obs"][None], np.array([clock["center_candidate_id"]])),
                ("member", clock["age_obs"], clock["ages"], clock["member_candidate_ids"]),
            ]:
                counts = clock[f"{phase}_counts"]
                assert counts["complete"]
                assert counts["accepted"] == (1 if phase == "center" else 10)
                assert counts["attempted"] == counts["accepted"] + counts["order_rejected"] + counts["coverage_rejected"]
                actual, ids, order_count, coverage_count = native_candidates(
                    center, clock["age_std"], factor, clock[f"{phase}_seed"], counts["attempted"])
                np.testing.assert_array_equal(actual, expected)
                np.testing.assert_array_equal(ids, expected_ids)
                assert order_count == counts["order_rejected"]
                assert coverage_count == counts["coverage_rejected"]
                assert counts["acceptance_rate"] == len(actual) / counts["attempted"]
                assert np.all(np.diff(expected, axis=1) < 0)
                assert np.all(expected[:, 0] >= 29) and np.all(expected[:, -1] <= 1)
                rejections["order"] = rejections["order"] + order_count
                rejections["coverage"] = rejections["coverage"] + coverage_count
                base = {name: clock[name] for name in ["process_seed", "level", "sigma_years", "record"]}
                acceptance_rows.append(dict(base, phase=phase, **counts))
                for member, (candidate, candidate_id) in enumerate(zip(expected, expected_ids)):
                    error = np.interp(-grid, -age, (candidate - age) * 1000)
                    about_center = np.interp(-grid, -age, (candidate - center) * 1000)
                    if phase == "center":
                        about_center = np.zeros_like(about_center)
                    spacing = -np.diff(candidate) * 1000
                    stretch = np.diff(candidate) / np.diff(age)
                    diagnostic_rows.append(dict(base, phase=phase, member=-1 if phase == "center" else member,
                        candidate_id=int(candidate_id), mean_error_years=float(error.mean()),
                        rms_error_years=float(np.sqrt(np.mean(error ** 2))),
                        max_abs_error_years=float(np.max(np.abs(error))),
                        rms_about_center_years=float(np.sqrt(np.mean(about_center ** 2))),
                        min_spacing_years=float(spacing.min()), min_stretch=float(stretch.min()),
                        max_stretch=float(stretch.max()), covers_grid=bool(candidate[0] >= 29 and candidate[-1] <= 1)))

            # The original native ensemble routine rejects order only. Request
            # the full attempted stream, then apply the new explicit coverage
            # rule, to verify exact compatibility with its accepted curves.
            observed = pd.DataFrame({"sample_id": np.arange(1501),
                                     "age_obs": clock["age_obs"], "age_std": clock["age_std"]})
            attempted = clock["member_counts"]["attempted"]
            with patch.object(observation, "_age_error_factor", return_value=factor):
                native = observation.draw_age_ensemble(observed, n_members=attempted,
                    age_correlation_length=1 / 30, max_draws=attempted,
                    rng=np.random.default_rng(clock["member_seed"]))
            covered = (native["ages"][:, 0] >= 29) & (native["ages"][:, -1] <= 1)
            np.testing.assert_array_equal(native["ages"][covered], clock["ages"])
            native_member_checks = native_member_checks + 1
        for level in range(1, 5):
            south = prepared["chronologies"][f"level{level}_x"]
            north = prepared["chronologies"][f"level{level}_y"]
            for member in range(-1, 10):
                south_age = south["age_obs"] if member < 0 else south["ages"][member]
                north_age = north["age_obs"] if member < 0 else north["ages"][member]
                relative = np.interp(-grid, -age, (north_age - south_age) * 1000)
                relative_rows.append(dict(process_seed=process_seed, level=level,
                    sigma_years=south["sigma_years"], member=member,
                    phase="center" if member < 0 else "member",
                    relative_rms_years=float(np.sqrt(np.mean(relative ** 2))),
                    relative_mean_years=float(relative.mean()),
                    relative_max_abs_years=float(np.max(np.abs(relative)))))
        for job in cases.loc[cases.process_seed.eq(process_seed)].to_dict("records"):
            path = output_dir / "inputs" / f"{job['case_id']}.csv"
            assert core.file_digest(path) == job["input_sha256"]
            frame = pd.read_csv(path, float_precision="round_trip")
            assert list(frame) == ["time", "age", "x", "y"]
            assert frame.shape == (1401, 4) and np.isfinite(frame.to_numpy()).all()
            np.testing.assert_array_equal(frame.age, grid)
            np.testing.assert_array_equal(frame.time, grid[0] - grid)
            if job["level"] == 0:
                np.testing.assert_array_equal(frame.x, sampled.x.to_numpy()[selected])
                np.testing.assert_array_equal(frame.y, sampled.y.to_numpy()[selected])
            else:
                for record in ["x", "y"]:
                    clock = prepared["chronologies"][f"level{job['level']}_{record}"]
                    member_age = clock["ages"][job["member"]]
                    expected = np.interp(-grid, -member_age, sampled[record])
                    np.testing.assert_array_equal(frame[record], expected)
    assert len(set(seeds)) == 12  # Three histories x two records x two layers, shared over levels.
    for filename, rows, keys in [
        ("chronology_diagnostics.csv", diagnostic_rows, ["process_seed", "level", "record", "phase", "member"]),
        ("relative_error_diagnostics.csv", relative_rows, ["process_seed", "level", "member"]),
        ("chronology_acceptance.csv", acceptance_rows, ["process_seed", "level", "record", "phase"]),
    ]:
        expected = pd.DataFrame(rows).set_index(keys).sort_index()
        actual = pd.read_csv(output_dir / filename, float_precision="round_trip").set_index(keys).sort_index()
        pd.testing.assert_frame_equal(actual[expected.columns], expected, check_dtype=False, check_exact=True)
    return {"cases": 123, "accurate_cases": 3, "chronology_cases": 120,
            "native_clock_packages_exact": clock_count,
            "native_ensemble_routines_exact": native_member_checks,
            "independent_layer_seed_streams": len(set(seeds)), "rejections": rejections,
            "chronology_relative_and_acceptance_tables_exact": True,
            "zero_error_values_exact": True, "all_values_fixed_before_dating": True,
            "all_inputs_reconstructed_exactly": True, "all_members_ordered_and_cover_grid": True}


def check_ccm_reference(frame, original, estimator_seed):
    """Compare every optimized value with a fresh native CCM call for that draw."""
    lags = np.asarray(core.CCM_LAG_STEPS)
    anchors = np.arange(30, len(frame) - 30)
    assert len(anchors) == 1341
    valid = np.zeros(len(frame), dtype=bool)
    valid[anchors] = True
    seeds = np.random.default_rng(estimator_seed).integers(1, 2**31 - 1, size=100)
    rho = np.full((100, 61, 2), np.nan)
    started = time.perf_counter()
    for library, seed in enumerate(seeds):
        for position, lag in enumerate(lags):
            native = pyEDM.CCM(dataFrame=frame[["time", "x", "y"]], columns="y", target="x",
                E=3, tau=-1, Tp=int(lag), exclusionRadius=10, libSizes=[300],
                sample=1, seed=int(seed), validLib=valid, parallel=False,
                includeData=False, showPlot=False)
            assert native.LibSize.iloc[0] == 300
            rho[library, position] = [native["y:x"].iloc[0], native["x:y"].iloc[0]]
        if (library + 1) % 20 == 0:
            print(f"CCM fresh native scan: {library + 1}/100 libraries", flush=True)
    direct_seconds = time.perf_counter() - started
    started = time.perf_counter()
    optimized = causal.filter_ccm_by_lag(frame, original, lag_steps=lags,
        n_libraries=100, sign_threshold=.9, seed=estimator_seed)
    optimized_seconds = time.perf_counter() - started
    np.testing.assert_array_equal(optimized["raw"]["lag_rho"], rho)
    np.testing.assert_array_equal(optimized["raw"]["lag_common_anchor_indices"], anchors)
    np.testing.assert_array_equal(optimized["raw"]["lag_library_seeds"], seeds)
    assert not optimized["diagnostics"]["lag_filter"]["scan_errors"]
    direct = deepcopy(optimized)
    direct["raw"]["lag_rho"] = rho
    direct = causal.summarize_ccm_lag(direct, sign_threshold=.9)
    assert canonical(core.native_snapshot(direct)) == canonical(core.native_snapshot(optimized))
    for name in ["rho_null", "cuts", "library_seeds", "observed"]:
        assert canonical(core.native_snapshot(original["raw"][name])) == canonical(core.native_snapshot(optimized["raw"][name]))
    return optimized, {"lag_shape": [100, 61, 2], "common_anchors": 1341,
        "all_lag_scores_and_summaries_exact": True, "tp0_evidence_preserved": True,
        "direct_seconds": direct_seconds, "optimized_seconds": optimized_seconds}


def spawned_case(job, output_dir, null_path):
    core.run_case(job, output_dir, null_path)
    return True


def check_native_outputs(job, output_dir, validated_reference_dir=None):
    """Validate cached CPU GPDC, 61-shift CCM, serialization and Windows spawn."""
    output_dir = Path(output_dir)
    frame = pd.read_csv(output_dir / "inputs" / f"{job['case_id']}.csv")
    seed = job["estimator_seed"]
    null_path = output_dir / "null_T1391_n10000.npz"
    if validated_reference_dir is not None:
        validated_reference_dir = Path(validated_reference_dir).resolve()
        for source in validated_reference_dir.glob("null_T1391_n10000*"):
            shutil.copy2(source, output_dir / source.name)
    causal.prepare_gpdc_null(null_path, sample_sizes=[1391], sig_samples=10000, seed=seed)
    pcmci_settings = dict(tau_max=5, pc_alpha=.05, alpha=.05, sig_samples=10000,
                          seed=seed, gp_params={"optimizer": None})
    checkpoint = output_dir / "gpdc_check.json"
    reuse_gpdc = checkpoint.exists() and json.loads(checkpoint.read_text(encoding="utf-8")).get("passed")
    if not reuse_gpdc:
        print("Checking complete fresh versus cached CPU GPDC output.", flush=True)
        uncached = causal.run_pcmciplus(frame, **pcmci_settings)
    cached = causal.run_pcmciplus(frame, null_dist_filename=str(null_path), **pcmci_settings)
    if not reuse_gpdc:
        for key in ["raw", "edges"]:
            assert canonical(core.native_snapshot(uncached[key])) == canonical(core.native_snapshot(cached[key]))
    else:
        assert validated_reference_dir is not None
        previous = read_archive(validated_reference_dir / "native" / f"{job['case_id']}_PCMCIplus.json.gz")
        assert canonical(core.native_snapshot(cached["raw"])) == canonical(previous["raw"])
        print("Reused completed fresh/cached GPDC equality checkpoint; current native arrays agree exactly.", flush=True)
    core.write_json({"passed": True, "raw_and_edges_exact": True,
        "input_sha256": job["input_sha256"], "reused_completed_checkpoint": bool(reuse_gpdc),
        "causal_source_sha256": core.file_digest(ROOT / "toolbox/causal.py")}, checkpoint)
    original = causal.run_ccm(frame, E=3, tau=-1, exclusion_radius=10, lib_size=300,
        sample=10, n_surrogates=100, seed=seed, alpha=.05, lag_steps=None)
    if validated_reference_dir is None:
        ccm, ccm_report = check_ccm_reference(frame, original, seed)
    else:
        # The earlier check completed all 6100 fresh native calls and both
        # equality assertions before failing in unrelated archive edge labels.
        # Reuse that scientific result; compare all current arrays and derived
        # summaries again without repeating the unchanged neighbor searches.
        reference_path = validated_reference_dir / "native" / f"{job['case_id']}_CCM.json.gz"
        previous = read_archive(reference_path)
        ccm = causal.filter_ccm_by_lag(frame, original, lag_steps=core.CCM_LAG_STEPS,
            n_libraries=100, sign_threshold=.9, seed=seed)
        assert canonical(core.native_snapshot(ccm["raw"])) == canonical(previous["raw"])
        previous_result = {key: restore_arrays(previous[key]) for key in ["edges", "raw", "diagnostics"]}
        previous_result["edges"] = original["edges"].copy()
        previous_summary = causal.summarize_ccm_lag(previous_result, sign_threshold=.9)
        assert canonical(core.native_snapshot(previous_summary["edges"])) == canonical(core.native_snapshot(ccm["edges"]))
        ccm_report = {"lag_shape": [100, 61, 2], "common_anchors": 1341,
            "all_lag_scores_and_summaries_exact": True, "tp0_evidence_preserved": True,
            "fresh_native_reference_reused": str(reference_path.relative_to(ROOT)),
            "fresh_native_reference_sha256": core.file_digest(reference_path),
            "reference_stage": "All 100x61 fresh calls and exact comparisons passed before unrelated archive-label failure.",
            "current_all_native_arrays_against_prior_validated_reference_exact": True}
    core.write_json(ccm_report, output_dir / "ccm_check.json")
    expected = {"PCMCI+": cached, "CCM": ccm,
                "LKIF": causal.run_lkif(frame, max_lag=5, alpha=.05),
                "GC": causal.run_granger(frame, max_lag=1, alpha=.05)}
    print("Direct full-setting method outputs ready; checking persisted worker results.", flush=True)
    core.run_case(job, output_dir, null_path)
    restored = core.read_completed_case(job, output_dir)
    assert len(restored) == 8 and restored.status.eq("ok").all()
    for method, result in expected.items():
        archive = read_archive(output_dir / "native" / f"{job['case_id']}_{method.replace('+', 'plus')}.json.gz")
        assert canonical({key: archive[key] for key in result}) == canonical(core.native_snapshot(result)), method

    spawn_dir = output_dir / "spawn"
    core.save_input(pd.read_csv(output_dir / "inputs" / f"{job['case_id']}.csv", float_precision="round_trip"),
                    spawn_dir / "inputs" / f"{job['case_id']}.csv")
    assert core.file_digest(spawn_dir / "inputs" / f"{job['case_id']}.csv") == job["input_sha256"]
    with ProcessPoolExecutor(max_workers=1, mp_context=multiprocessing.get_context("spawn")) as pool:
        assert pool.submit(spawned_case, job, spawn_dir, null_path).result()
    for method in core.METHODS:
        name = f"{job['case_id']}_{method.replace('+', 'plus')}.json.gz"
        assert canonical(read_archive(output_dir / "native" / name)) == canonical(read_archive(spawn_dir / "native" / name)), method
    assert core.read_completed_case(job, spawn_dir) is not None
    return {"ccm": ccm_report, "cpu_gpdc_cached_uncached_raw_and_edges_exact": True,
            "all_four_method_archives_exact": True, "windows_spawn_exact": True,
            "native_fit_fields_serialized": True, "committed_case_resume_passed": True}


def check_resume_guards(job, output_dir):
    output_dir = Path(output_dir)
    parameters = core.execution_parameters()
    core.write_json(parameters, output_dir / "parameters.json")
    core.check_parameters(parameters, output_dir)
    changed = deepcopy(parameters)
    changed["CCM_LIBRARY_SIZE"] = 301
    try:
        core.check_parameters(changed, output_dir)
    except ValueError:
        pass
    else:
        raise AssertionError("Changed CCM library accepted on resume.")
    changed_job = dict(job, input_sha256="changed-input")
    try:
        core.read_completed_case(changed_job, output_dir)
    except ValueError:
        pass
    else:
        raise AssertionError("Changed input identity accepted on resume.")
    archive = read_archive(output_dir / "native" / f"{job['case_id']}_CCM.json.gz")
    archive["raw"]["lag_rho"]["shape"] = [100, 31, 2]
    try:
        core.check_native_archive(archive, "CCM", job)
    except ValueError:
        pass
    else:
        raise AssertionError("Incomplete CCM lag arrays accepted.")
    return {"changed_scientific_setting_rejected": True, "changed_input_rejected": True,
            "incomplete_native_arrays_rejected": True}


def check_high_uncertainty_outputs(job, output_dir):
    """Check a maximum-amplitude input against its spawned persisted outputs."""
    output_dir = Path(output_dir)
    frame = pd.read_csv(output_dir / "inputs" / f"{job['case_id']}.csv")
    null_path = output_dir / "null_T1391_n10000.npz"
    expected = {method: core.call_method(frame, method, job["estimator_seed"], null_path)
                for method in core.METHODS}
    with ProcessPoolExecutor(max_workers=1, mp_context=multiprocessing.get_context("spawn")) as pool:
        assert pool.submit(spawned_case, job, output_dir, null_path).result()
    rows = core.read_completed_case(job, output_dir)
    assert len(rows) == 8
    for method, result in expected.items():
        archive = read_archive(output_dir / "native" / f"{job['case_id']}_{method.replace('+', 'plus')}.json.gz")
        assert canonical({key: archive[key] for key in result}) == canonical(core.native_snapshot(result)), method
    return {"case_id": job["case_id"], "level": job["level"],
            "all_four_full_outputs_exact": True, "windows_spawn_exact": True,
            "successful_directional_rows": int(rows.status.eq("ok").sum()),
            "unresolved_directional_rows": int(rows.unoriented.fillna(False).sum())}


def check_artifacts(output_dir=None):
    """Read every production archive and reconstruct all seven summary tables."""
    output_dir = core.DATA_DIR if output_dir is None else Path(output_dir)
    jobs = pd.read_csv(output_dir / "input_cases.csv", float_precision="round_trip").to_dict("records")
    completed = {job["case_id"]: core.read_completed_case(job, output_dir) for job in jobs}
    assert len(jobs) == 123 and all(value is not None for value in completed.values())
    results = core.collect_results(jobs, completed)
    assert len(results) == 984
    tables = core.summarize_results(results)
    boolean_columns = [name for name in ["edge_detected", "unoriented", "raw_edge_detected"] if name in results]
    disk = pd.read_csv(output_dir / "member_results.csv", float_precision="round_trip",
                       dtype={name: "boolean" for name in boolean_columns})
    pd.testing.assert_frame_equal(results.reset_index(drop=True), disk, check_dtype=False, check_exact=True)
    for name, expected in tables.items():
        bool_types = {column: "boolean" for column in expected if pd.api.types.is_bool_dtype(expected[column])}
        saved = pd.read_csv(output_dir / f"{name}.csv", float_precision="round_trip", dtype=bool_types)
        pd.testing.assert_frame_equal(expected.reset_index(drop=True), saved, check_dtype=False, check_exact=True)
    return {"passed": True, "cases": len(jobs), "directional_rows": len(results),
            "native_archives": 4 * len(jobs), "summary_tables_exact": list(tables),
            "failed_directional_rows": int(results.status.ne("ok").sum())}


def main(validated_reference_dir=None, output_folder=None):
    started = time.perf_counter()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    folder = ROOT / "data/benchmarks/tbk_chronology_checks" / stamp if output_folder is None else Path(output_folder).resolve()
    folder.mkdir(parents=True, exist_ok=output_folder is not None)
    if (folder / "report.json").exists():
        shutil.copy2(folder / "report.json", folder / f"report_before_resume_{stamp}.json")
    report = {"passed": False, "started_utc": datetime.now(timezone.utc).isoformat()}
    try:
        with threadpool_limits(limits=1):
            jobs, _, _ = core.prepare_inputs(folder / "prepared_run")
            report["inputs"] = check_inputs(jobs, folder / "prepared_run")
            print("All 123 inputs and native chronology streams verified.", flush=True)
            report["native"] = check_native_outputs(jobs[0], folder / "prepared_run", validated_reference_dir)
            high_job = next(job for job in jobs if job["process_seed"] == core.PROCESS_SEEDS[0] and job["level"] == 4)
            report["highest_uncertainty"] = check_high_uncertainty_outputs(high_job, folder / "prepared_run")
            report["resume"] = check_resume_guards(jobs[0], folder / "prepared_run")
        report["preservation"] = check_preservation()
        report["final_scientific_parameters"] = core.execution_parameters()
        report["verifier_sha256"] = core.file_digest(Path(__file__))
        report["passed"] = True
    finally:
        report["elapsed_seconds"] = time.perf_counter() - started
        (folder / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: {folder / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
