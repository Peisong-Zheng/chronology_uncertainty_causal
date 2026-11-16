"""Independently verify the 03_2 Z dating-error-model comparison.

Checks use isolated benchmark folders. Existing 03 outputs are read-only. The
native check compares one new E-type input directly and in two Windows
workers; it never refits CCM or launches the production experiment.
"""

import os
for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"]:
    os.environ[name] = "1"
os.environ.setdefault("MPLBACKEND", "Agg")

from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import ExitStack
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
from threadpoolctl import threadpool_limits
from agent.verify_common_driver import canonical, check_conditional_mapping, expected_edges, frame_from_snapshot

source = importlib.import_module("03_check_common_driver")
cross = importlib.import_module("03_1_check_process_chronology_cross")
core = importlib.import_module("03_2_check_driver_error_models")
ANALYSES = [("pairwise", method) for method in ["PCMCI+", "LKIF", "GC", "CCM"]]
ANALYSES = ANALYSES + [("conditional", method) for method in ["PCMCI+", "LKIF", "GC"]]


def file_state(path):
    path = Path(path)
    stat = path.stat()
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    return {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns, "sha256": digest}


def expect_rejection(function, *args, **kwargs):
    try:
        function(*args, **kwargs)
    except (ValueError, RuntimeError, AssertionError):
        return
    raise AssertionError(f"Scientific mismatch accepted by {function.__name__}.")


def check_preservation(snapshot=ROOT / "tmp/03_2_existing_files.json"):
    """Read every protected file and compare all original metadata and bytes."""
    snapshot = json.loads(Path(snapshot).read_text(encoding="utf-8"))
    for record in snapshot["files"]:
        path = ROOT / record["path"]
        assert path.is_file(), f"Missing pre-existing file: {path}"
        actual = file_state(path)
        for key in ["size", "mtime_ns", "sha256"]:
            assert actual[key] == record[key], f"Changed pre-existing {key}: {path}"
    return {"state": "passed", "files": len(snapshot["files"]), "bytes": snapshot["bytes"],
            "checked_utc": datetime.now(timezone.utc).isoformat(),
            "checks": ["existence", "size", "mtime_ns", "sha256"]}


def independent_clocks(prepared, source_dir):
    """Regenerate all envelopes/centres/members through the unchanged toolbox."""
    from toolbox.observation import make_age_uncertainty, simulate_observations, draw_age_ensemble
    sampling = pd.read_csv(Path(source_dir) / "sampling.csv", float_precision="round_trip")
    indices = sampling.native_index.to_numpy(dtype=int)
    age_true = sampling.age_true.to_numpy()
    full_age = 30.0 - np.arange(1501) * .020
    np.testing.assert_array_equal(full_age[indices], age_true)
    centers = 1 + 28 * np.array([.5, .2, .8, .35, .65])
    intervals = [(float(center - 1.4), float(center + 1.4)) for center in centers]
    clocks, diagnostics, acceptance_records = {}, [], []
    grid = np.linspace(29, 1, 1401)
    for seed in source.PROCESS_SEEDS:
        age_seed = int(np.random.SeedSequence([seed, 3103, 2]).generate_state(1)[0])
        ensemble_seed = int(np.random.SeedSequence([seed, 3104, 2]).generate_state(1)[0])
        for kind in ["C", "L", "E", "J"]:
            std = make_age_uncertainty(age_true, kind=kind,
                scale=.010 if kind == "J" else .040, age_range=(0, 30), exp_shape=2,
                jump_intervals=intervals if kind == "J" else None)
            observation = simulate_observations(full_age, np.zeros(1501), indices, std,
                age_correlation_length=1 / 60, measurement_noise_ratio=0,
                rng_age=age_seed, rng_measurement=0)
            ensemble = draw_age_ensemble(observation["observed"], n_members=10,
                age_correlation_length=1 / 60, max_draws=20000, rng=ensemble_seed)
            assert ensemble["diagnostics"]["complete"]
            expected = {"age_true": age_true, "age_obs": observation["observed"].age_obs.to_numpy(),
                "age_std": std, "ages": ensemble["ages"], "sample_id": np.arange(len(indices)),
                "native_indices": indices}
            path = Path(prepared) / "clocks" / f"s{seed}_z_{kind}.npz"
            with np.load(path, allow_pickle=False) as data:
                actual = {key: data[key].copy() for key in data.files}
            assert set(expected).issubset(actual)
            for key, values in expected.items():
                np.testing.assert_array_equal(actual[key], values, err_msg=f"{seed}/{kind}/{key}")
            assert np.all(np.diff(actual["ages"], axis=1) < 0)
            assert np.all(actual["ages"][:, 0] >= 29) and np.all(actual["ages"][:, -1] <= 1)
            assert not np.array_equal(actual["age_obs"], age_true)
            if kind == "C":
                assert path.read_bytes() == (Path(source_dir) / "clocks" / path.name).read_bytes()
            if kind == "J":
                fixed = std == 0
                np.testing.assert_array_equal(actual["age_obs"][fixed], age_true[fixed])
                np.testing.assert_array_equal(actual["ages"][:, fixed], np.tile(age_true[fixed], (10, 1)))
            clocks[seed, kind] = actual
            grid_std = make_age_uncertainty(grid, kind=kind,
                scale=.010 if kind == "J" else .040, age_range=(0, 30), exp_shape=2,
                jump_intervals=intervals if kind == "J" else None)
            error_curves = np.array([np.interp(-grid, -age_true, ages - age_true) * 1000
                                     for ages in actual["ages"]])
            proposal_rms = np.sqrt(np.mean((1000 * grid_std)**2))
            effective_sd = np.mean(np.std(error_curves, axis=0, ddof=0))
            acceptance_records.append({"chronology_seed": seed, "age_type": kind,
                "age_seed": age_seed, "ensemble_seed": ensemble_seed, **ensemble["diagnostics"],
                "proposal_rms_years": proposal_rms, "effective_member_sd_years": effective_sd})
            for member, ages in enumerate(actual["ages"]):
                error = np.interp(-grid, -age_true, ages - age_true) * 1000
                center_error = np.interp(-grid, -age_true, actual["age_obs"] - age_true) * 1000
                diagnostics.append({"chronology_seed": seed, "age_type": kind, "member": member,
                    "rms_error_years": np.sqrt(np.mean(error**2)), "mean_error_years": np.mean(error),
                    "observed_rms_years": np.sqrt(np.mean(center_error**2)),
                    "proposal_rms_years": proposal_rms, "effective_member_sd_years": effective_sd,
                    "minimum_gap_years": np.min(-np.diff(ages)) * 1000})
    saved_diagnostics = pd.read_csv(Path(prepared) / "clock_diagnostics.csv", float_precision="round_trip")
    assert len(saved_diagnostics) == 120
    keys = ["chronology_seed", "age_type", "member"]
    expected_diagnostics = pd.DataFrame(diagnostics).set_index(keys).sort_index()
    saved_diagnostics = saved_diagnostics.set_index(keys).sort_index()
    for name in expected_diagnostics:
        np.testing.assert_array_equal(saved_diagnostics[name], expected_diagnostics[name], err_msg=name)
    assert saved_diagnostics.covers_analysis_window.all()
    np.testing.assert_array_equal(saved_diagnostics.relative_x_z_rms_years, saved_diagnostics.rms_error_years)
    np.testing.assert_array_equal(saved_diagnostics.relative_y_z_rms_years, saved_diagnostics.rms_error_years)
    assert saved_diagnostics.relative_x_y_rms_years.eq(0).all()
    expected_acceptance = pd.DataFrame(acceptance_records).set_index(["chronology_seed", "age_type"]).sort_index()
    actual_acceptance = pd.read_csv(Path(prepared) / "chronology_acceptance.csv", float_precision="round_trip")
    assert len(actual_acceptance) == 12
    actual_acceptance = actual_acceptance.set_index(["chronology_seed", "age_type"]).sort_index()
    for name in expected_acceptance:
        np.testing.assert_array_equal(actual_acceptance[name], expected_acceptance[name], err_msg=name)
    assert actual_acceptance.loc[(slice(None), "C"), "source_c_verified"].all()
    return clocks


def check_inputs(jobs, prepared, source_dir=None):
    """Independently rebuild 720 Z inputs and verify 192 byte-identical imports."""
    prepared = Path(prepared)
    source_dir = core.SOURCE_DIR if source_dir is None else Path(source_dir)
    cases = pd.DataFrame(jobs)
    assert len(cases) == cases.case_id.nunique() == 732
    assert cases.groupby("condition").size().to_dict() == {
        "original_uniform": 6, "irregular_accurate": 6, "z_only": 720}
    crossed = cases[cases.condition.eq("z_only")]
    assert crossed.groupby("age_type").size().to_dict() == {kind: 180 for kind in ["C", "E", "J", "L"]}
    assert cases.loc[cases.condition.ne("z_only"), "age_type"].eq("accurate").all()
    assert crossed.groupby(["structure", "process_seed", "chronology_seed", "age_type"]).size().eq(10).all()
    assert len(crossed.groupby(["structure", "process_seed", "chronology_seed", "age_type"])) == 72
    expected_import = cases.condition.ne("z_only") | cases.age_type.eq("C")
    np.testing.assert_array_equal(cases.full_result_reused.to_numpy(dtype=bool), expected_import.to_numpy())
    assert int(expected_import.sum()) == 192
    for row in cases.itertuples():
        assert row.realization_id == source.PROCESS_SEEDS.index(row.process_seed) + 1
        assert row.age_scheme_id == (source.PROCESS_SEEDS.index(row.chronology_seed) + 1 if row.chronology_seed else 0)
        assert row.estimator_seed == row.process_seed + source.METHOD_SEED_OFFSET
    assert (prepared / "sampling.csv").read_bytes() == (source_dir / "sampling.csv").read_bytes()
    clocks = independent_clocks(prepared, source_dir)
    grid = np.linspace(29, 1, 1401)
    full, accurate, accurate_native = {}, {}, {}
    for seed in source.PROCESS_SEEDS:
        for structure in source.STRUCTURES:
            full[seed, structure] = pd.read_csv(source.DATA_DIR / "generated" / f"s{seed}_{structure}.csv", float_precision="round_trip")
            baseline_id = cross.case_name(seed, structure, "irregular_accurate")
            accurate[seed, structure] = pd.read_csv(source_dir / "inputs" / f"{baseline_id}.csv", float_precision="round_trip")
            accurate_native[seed, structure] = pd.read_csv(source_dir / "inputs" / f"{baseline_id}.csv")
    loaded, copied_inputs, reconstructed = {}, 0, 0
    for job in jobs:
        path = prepared / "inputs" / f"{job['case_id']}.csv"
        series = pd.read_csv(path, float_precision="round_trip")
        assert list(series) == ["time", "age", "x", "y", "z"]
        assert series.shape == (1401, 5) and np.isfinite(series.to_numpy()).all()
        np.testing.assert_array_equal(series.age, grid)
        np.testing.assert_array_equal(series.time, grid[0] - grid)
        seed, structure = job["process_seed"], job["structure"]
        if job["full_result_reused"]:
            assert path.read_bytes() == (source_dir / "inputs" / f"{job['source_case_id']}.csv").read_bytes()
            copied_inputs = copied_inputs + 1
        if job["condition"] == "z_only":
            np.testing.assert_array_equal(series[["time", "age", "x", "y"]], accurate[seed, structure][["time", "age", "x", "y"]])
            native = pd.read_csv(path)
            np.testing.assert_array_equal(native[["time", "age", "x", "y"]], accurate_native[seed, structure][["time", "age", "x", "y"]])
            clock = clocks[job["chronology_seed"], job["age_type"]]
            expected_z = np.interp(-grid, -clock["ages"][job["member"]], full[seed, structure].z.to_numpy()[clock["native_indices"]])
            np.testing.assert_array_equal(series.z, expected_z)
            loaded[seed, structure, job["chronology_seed"], job["age_type"], job["member"]] = series
            reconstructed = reconstructed + 1
    for seed in source.PROCESS_SEEDS:
        for chronology_seed in source.PROCESS_SEEDS:
            for kind in ["C", "L", "E", "J"]:
                for member in range(10):
                    first = loaded[seed, "h0_common_driver", chronology_seed, kind, member]
                    second = loaded[seed, "h1_direct", chronology_seed, kind, member]
                    np.testing.assert_array_equal(first[["time", "age", "x", "z"]], second[["time", "age", "x", "z"]])
    assert copied_inputs == 192 and reconstructed == 720
    return {"cases": 732, "baseline_cases": 12, "C_cases": 180, "new_type_cases": 540,
            "byte_identical_copied_inputs": copied_inputs, "independently_reconstructed_z_inputs": reconstructed,
            "independently_generated_clock_packages": 12, "C_clock_archives_byte_identical": 3,
            "xy_unchanged_under_native_parser": True, "clock_and_process_pairing_across_systems": True,
            "centres_and_members_reproduced_from_native_toolbox": True, "actual_error_diagnostics_exact": True}


def prepare_without_new_processes(prepared, source_dir=None):
    """New clocks are allowed, but dynamics and observation positions are fixed."""
    def forbidden(*args, **kwargs):
        raise AssertionError("03_2 must reuse saved processes and sampling positions.")
    with ExitStack() as stack:
        for module in [core, source, core.base, core.cross]:
            for name in ["generate_common_driver", "sample_irregular_indices"]:
                if hasattr(module, name):
                    stack.enter_context(patch.object(module, name, side_effect=forbidden))
        return core.prepare_inputs(prepared, core.SOURCE_DIR if source_dir is None else source_dir)


def fake_conditional(series, analysis, method, seed, null):
    """Deliberate failures/unresolved evidence; no numerical estimator executes."""
    assert analysis == "conditional", "Pairwise evidence must be copied."
    if method == "LKIF":
        raise RuntimeError("Deliberate isolated verification failure.")
    rows = []
    for first, second in expected_edges(analysis):
        xy = {first, second} == {"x", "y"}
        row = {"method": method, "source": first, "target": second, "status": "ok",
               "score": .4, "p_value": .7, "edge_detected": False,
               "unoriented": method == "PCMCI+" and xy}
        if method == "GC" and not xy:
            row.update(status="error", score=np.nan, edge_detected=None)
        rows.append(row)
    return {"edges": pd.DataFrame(rows), "raw": {"array": np.arange(6.).reshape(3, 2)}, "diagnostics": {}}


def copy_case_input(job, prepared, folder):
    destination = Path(folder) / "inputs" / f"{job['case_id']}.csv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(Path(prepared) / "inputs" / destination.name, destination)
    return destination


def check_guards(jobs, prepared, output):
    output = Path(output)
    job = next(job for job in jobs if job["structure"] == "h0_common_driver"
               and job["condition"] == "z_only" and not job["full_result_reused"])
    imported_job = next(job for job in jobs if job["condition"] == "z_only" and job["full_result_reused"])
    imported = output / "complete_import"
    copy_case_input(imported_job, prepared, imported)
    with patch.object(core, "call_method", side_effect=AssertionError("A full source import must not fit.")):
        inherited = core.run_case(imported_job, imported, None)
    assert len(inherited) == 26 and inherited.native_reused.all()
    for analysis, method in ANALYSES:
        assert core.native_path(imported, imported_job["case_id"], analysis, method).read_bytes() == source.native_path(
            core.SOURCE_DIR, imported_job["source_case_id"], analysis, method).read_bytes()

    folder = output / "failure_and_resume"
    input_path = copy_case_input(job, prepared, folder)
    with patch.object(core, "call_method", side_effect=fake_conditional) as calls:
        result = core.run_case(job, folder, "unused_mock_null")
    assert calls.call_count == 3 and len(result) == 26
    assert result.loc[result.analysis.eq("pairwise"), "native_reused"].all()
    assert not result.loc[result.analysis.eq("conditional"), "native_reused"].any()
    failed = result.analysis.eq("conditional") & result.method.eq("LKIF")
    assert result.loc[failed, "status"].eq("error").all() and result.loc[failed, "edge_detected"].isna().all()
    state_table = core.summarize_results(result)["xy_case_states"].set_index(["analysis", "method"])
    assert state_table.loc[("conditional", "PCMCI+"), "state"] == "unresolved"
    assert state_table.loc[("conditional", "LKIF"), "failed"]
    assert state_table.loc[("conditional", "GC"), "correct_xy"]
    before = {path.relative_to(folder).as_posix(): file_state(path) for path in folder.rglob("*") if path.is_file()}
    resumed = core.read_completed_case(job, folder)
    after = {path.relative_to(folder).as_posix(): file_state(path) for path in folder.rglob("*") if path.is_file()}
    assert resumed is not None and len(resumed) == 26 and before == after
    assert resumed.loc[failed, "status"].eq("error").all()

    for label, shortened in [("short_commit", True), ("missing_native", False)]:
        partial = output / label
        shutil.copytree(folder, partial)
        if shortened:
            result.iloc[:-1].to_csv(partial / "members" / f"{job['case_id']}_edges.csv", index=False)
        else:
            core.native_path(partial, job["case_id"], "conditional", "GC").unlink()
        assert core.read_completed_case(job, partial) is None

    original_input = input_path.read_bytes()
    changed = pd.read_csv(input_path, float_precision="round_trip")
    changed.loc[1, "x"] = changed.loc[1, "x"] + .01
    changed.to_csv(input_path, index=False)
    try:
        expect_rejection(core.read_completed_case, job, folder)
        expect_rejection(core.validate_case_source, job, folder)
        changed_job = dict(job, input_sha256=file_state(input_path)["sha256"])
        expect_rejection(core.validate_case_source, changed_job, folder)
    finally:
        input_path.write_bytes(original_input)
    expect_rejection(core.validate_case_source, dict(job, estimator_seed=job["estimator_seed"] + 1), folder)
    wrong_type = "E" if job["age_type"] != "E" else "J"
    expect_rejection(core.validate_case_source, dict(job, age_type=wrong_type), folder)
    expect_rejection(core.read_completed_case, dict(job, age_type=wrong_type), folder)
    before_input = file_state(input_path)
    expect_rejection(core.save_input, changed, input_path)
    assert file_state(input_path) == before_input

    serialization = output / "serialization_failure"
    copy_case_input(job, prepared, serialization)
    def unserializable(*args):
        value = fake_conditional(*args)
        value["raw"]["unsupported"] = object()
        return value
    with patch.object(core, "call_method", side_effect=unserializable):
        try:
            core.run_case(job, serialization, "unused_mock_null")
        except TypeError:
            pass
        else:
            raise AssertionError("Serialization failure produced a valid case commit.")
    assert not (serialization / "members" / f"{job['case_id']}_edges.csv").exists()

    # A small isolated source case is sufficient to exercise per-case imports.
    # Modify only this copy to represent a previously committed scientific failure.
    failed_source = output / "isolated_failed_source"
    source_id = job["source_case_id"]
    for path in core.source_case_paths(source_id, core.SOURCE_DIR):
        destination = failed_source / path.relative_to(core.SOURCE_DIR)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, destination)
    commit_path = failed_source / "members" / f"{source_id}_edges.csv"
    failed_commit = pd.read_csv(commit_path, float_precision="round_trip")
    selection = failed_commit.analysis.eq("pairwise") & failed_commit.method.eq("LKIF")
    failed_commit["edge_detected"] = failed_commit.edge_detected.astype("boolean")
    failed_commit.loc[selection, ["score", "p_value"]] = np.nan
    failed_commit.loc[selection, "edge_detected"] = pd.NA
    failed_commit.loc[selection, "status"] = "error"
    failed_commit.to_csv(commit_path, index=False)
    failed_native = source.native_path(failed_source, source_id, "pairwise", "LKIF")
    failed_edges = source.empty_edges("pairwise", "LKIF", "error", "Preserved source failure.")
    with gzip.open(failed_native, "wt", encoding="utf-8") as handle:
        json.dump(source.native_snapshot({"edges": failed_edges, "raw": {}, "diagnostics": {"error": "Preserved source failure."}}), handle)
    retained_job = dict(job, source_case_sha256=core.source_case_digest(source_id, failed_source))
    retained_folder = output / "retained_source_failure"
    copy_case_input(retained_job, prepared, retained_folder)
    with patch.object(core, "call_method", side_effect=fake_conditional) as calls:
        retained = core.run_case(retained_job, retained_folder, "unused_mock_null", failed_source)
    assert calls.call_count == 3
    selected = retained.analysis.eq("pairwise") & retained.method.eq("LKIF")
    assert retained.loc[selected, "status"].eq("error").all()
    assert retained.loc[selected, "edge_detected"].isna().all()
    assert core.native_path(retained_folder, retained_job["case_id"], "pairwise", "LKIF").read_bytes() == failed_native.read_bytes()
    assert core.read_completed_case(retained_job, retained_folder) is not None
    changed_source = failed_source / "inputs" / f"{source_id}.csv"
    with changed_source.open("ab") as handle:
        handle.write(b"\n")
    expect_rejection(core.validate_case_source, retained_job, retained_folder, failed_source)

    parameter_dir = output / "parameter_guards"
    parameter_dir.mkdir()
    parameters = core.execution_parameters()
    (parameter_dir / "parameters.json").write_text(json.dumps(parameters), encoding="utf-8")
    core.check_parameters(parameters, parameter_dir)
    permissible = deepcopy(parameters)
    permissible["N_WORKERS"] = 16
    core.check_parameters(permissible, parameter_dir)
    guarded = []
    for key in ["N_AGE_MEMBERS", "CHRONOLOGY_SEEDS", "AGE_TYPES", "TYPE_SCALES_YEARS",
                "AGE_EXP_SHAPE", "JUMP_COUNT", "JUMP_INTERVAL_WIDTH", "JUMP_CENTER_FRACTIONS",
                "AGE_CORRELATION_LENGTH", "age_range_kyr_bp", "source_sha256",
                "source_parameters", "calculation_sha256"]:
        modified = deepcopy(parameters)
        modified[key] = "deliberate mismatch"
        expect_rejection(core.check_parameters, modified, parameter_dir)
        guarded.append(key)
    for field in ["packages", "environment", "implementation_sha256", "CCM_LIBRARY_SIZE", "METHOD_SEED_OFFSET"]:
        modified = deepcopy(parameters)
        modified["source_parameters"]["source_parameters"][field] = "deliberate mismatch"
        expect_rejection(core.check_parameters, modified, parameter_dir)
        invalid_source = output / f"invalid_source_{field}"
        invalid_source.mkdir()
        metadata = json.loads((core.SOURCE_DIR / "parameters.json").read_text(encoding="utf-8"))
        metadata["source_parameters"][field] = "deliberate mismatch"
        (invalid_source / "parameters.json").write_text(json.dumps(metadata), encoding="utf-8")
        expect_rejection(core.validate_source, invalid_source)
        guarded.append(f"source_parameters.source_parameters.{field}")
    return {"complete_case_import_without_fitting": True, "only_three_new_conditional_calls": True,
            "new_and_source_failures_retained": True, "xy_state_not_contaminated_by_z_edge_failure": True,
            "unresolved_not_correct": True, "resume_no_writes": True, "partial_commits_rejected": True,
            "serialization_failure_not_committed": True, "input_bytes_and_actual_xy_guarded": True,
            "estimator_seed_guarded": True, "source_bytes_guarded": True,
            "age_type_identity_guarded": True,
            "changed_input_rejected_without_overwrite": True, "parameter_guards": guarded}


def native_worker(input_path, method, estimator_seed, null_file):
    """Module-level worker is imported under Windows spawn."""
    series = pd.read_csv(input_path)
    with threadpool_limits(limits=1):
        result = core.call_method(series, "conditional", method, estimator_seed, null_file)
    return os.getpid(), method, source.native_snapshot(result)


def check_native(jobs, prepared, output, report, save_report):
    job = next(job for job in jobs if job["structure"] == "h0_common_driver"
               and job["condition"] == "z_only" and job["realization_id"] == 1
               and job["age_scheme_id"] == 2 and job["member"] == 0 and job["age_type"] == "E")
    input_path = Path(prepared) / "inputs" / f"{job['case_id']}.csv"
    null = source.gpdc_null_path(job["estimator_seed"])
    assert null.is_file(), "Reuse the existing validated source03 GPDC null cache."
    series = pd.read_csv(input_path)
    reference = {}
    for method in ["PCMCI+", "LKIF", "GC"]:
        print(f"Direct conditional check: {job['case_id']} {method}", flush=True)
        with threadpool_limits(limits=1):
            # The reference is the frozen source03 wrapper, not the new caller.
            result = source.call_method(series, "conditional", method, job["estimator_seed"], str(null))
        assert result["edges"].status.eq("ok").all()
        reference[method] = source.native_snapshot(result)
        report.setdefault("conditional_mapping", []).append(check_conditional_mapping(series, method, result))
        save_report()
    workers = set()
    compared = []
    with ProcessPoolExecutor(max_workers=2, mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(native_worker, str(input_path), method, job["estimator_seed"], str(null))
                   for method in ["PCMCI+", "LKIF", "GC"]]
        for future in as_completed(futures):
            pid, method, payload = future.result()
            assert canonical(payload) == canonical(reference[method]), method
            workers.add(pid)
            compared.append(method)
            report["native_compared"] = sorted(compared)
            save_report()
    assert len(workers) == 2, f"Expected both Windows workers to execute: {workers}"
    return {"case_id": job["case_id"], "methods": sorted(compared), "windows_worker_pids": sorted(workers),
            "complete_native_payloads_exact": True, "ignored_fields": [], "native_calls": 6, "ccm_refits": 0}


def same_scalar(first, second):
    return bool((pd.isna(first) and pd.isna(second)) or first == second)


def check_serialized_arrays(value):
    """Check the stored values, shape and dtype of every nested native array."""
    count = 0
    if isinstance(value, dict):
        if {"dtype", "shape", "values"}.issubset(value):
            dtype = np.dtype(value["dtype"])
            if dtype.hasobject:
                # A native one-dimensional object array can contain tuples;
                # JSON preserves those as lists rather than extra array axes.
                def check_object_shape(items, dimensions):
                    if dimensions:
                        assert isinstance(items, list) and len(items) == dimensions[0]
                        for item in items:
                            check_object_shape(item, dimensions[1:])
                check_object_shape(value["values"], value["shape"])
            else:
                array = np.asarray(value["values"], dtype=dtype)
                assert list(array.shape) == value["shape"], (array.shape, value["shape"])
            return 1
        for nested in value.values():
            count = count + check_serialized_arrays(nested)
    elif isinstance(value, list):
        for nested in value:
            count = count + check_serialized_arrays(nested)
    return count


def check_summaries(jobs):
    """Test requested denominators and X-Y scope with one ten-member cell."""
    selected = [job for job in jobs if job["structure"] == "h0_common_driver" and job["realization_id"] == 1
                and ((job["age_scheme_id"] == 2 and job["age_type"] in ["E", "J"]) or job["condition"] == "irregular_accurate")]
    pending = core.summarize_results(core.collect_results(selected, {}))
    assert set(pending) == {"xy_case_states", "cell_summary", "paired_score_changes"}
    assert pending["xy_case_states"].pending.all() and not pending["xy_case_states"].failed.any()
    assert not pending["cell_summary"].correct_xy.any()
    assert pending["cell_summary"].query("condition == 'z_only'").requested.eq(10).all()
    rows = []
    for job in selected:
        for analysis, method in ANALYSES:
            for first, second in expected_edges(analysis):
                row = {**job, "analysis": analysis, "method": method, "source": first, "target": second,
                       "score": .4 if job["condition"] == "z_only" else .6,
                       "p_value": .7, "edge_detected": False, "unoriented": False, "status": "ok"}
                xy = {first, second} == {"x", "y"}
                if job["condition"] == "z_only":
                    if job["member"] == 0 and xy:
                        row.update(status="error", edge_detected=None, score=np.nan)
                    elif job["member"] == 1 and xy:
                        row["unoriented"] = True
                    elif job["member"] == 2 and first == "x" and second == "y":
                        row["edge_detected"] = True
                    elif job["member"] == 3 and not xy:
                        row.update(status="error", edge_detected=None, score=np.nan)
                    elif job["member"] == 4 and first == "x" and second == "y":
                        row["score"] = np.nan
                rows.append(row)
    results = pd.DataFrame(rows)
    tables = core.summarize_results(results)
    cells = tables["cell_summary"].query("condition == 'z_only'")
    assert len(cells) == 14
    assert set(cells.age_type) == {"E", "J"}
    for row in cells.itertuples():
        assert row.requested == 10 and row.usable == 9 and row.correct_xy == 7
        assert row.r_xy == .7 and row.unresolved == 1 and row.failed == 1 and row.pending == 0
        assert row.q_x_to_y == 1 / 9 and row.q_y_to_x == 0 and row.support_eligible
    changes = tables["paired_score_changes"]
    assert set(changes.source) <= {"x", "y"} and set(changes.target) <= {"x", "y"}
    assert len(changes) == 280
    assert changes.groupby("age_type").size().to_dict() == {"E": 140, "J": 140}
    assert changes.loc[changes.member.eq(4) & changes.source.eq("x"), "score_change"].isna().all()
    np.testing.assert_allclose(changes.loc[changes.usable_score_pair, "score_change"], -.2, atol=1e-15, rtol=0)
    return {"one_cell_requested_count": 10, "correct_count_including_unresolved_and_failed_denominator": 7,
            "baseline_requested_count": 1, "xy_scope_only": True, "pending_separate_from_failed": True,
            "unresolved_excluded_from_correct": True, "z_edge_failure_does_not_change_xy_state": True,
            "paired_scores_require_finite_success": True}


def audit_outputs(output_dir, reference_inputs=None, source_dir=None):
    """Read every archive and reconstruct saved summaries without any fit."""
    output = Path(output_dir)
    source_dir = core.SOURCE_DIR if source_dir is None else Path(source_dir)
    core.check_parameters(core.execution_parameters(source_dir), output)
    cases = pd.read_csv(output / "cases.csv", float_precision="round_trip")
    jobs = cases.to_dict("records")
    inputs = check_inputs(jobs, output, source_dir)
    source_copies, new_native, gc_fits, ccm_arrays, serialized_arrays = 0, 0, 0, 0, 0
    complete_imports, pairwise_imports, reference_count, lag_tables = 0, 0, 0, 0
    records, failures = [], []
    for job in jobs:
        if reference_inputs is not None:
            filename = f"{job['case_id']}.csv"
            assert (output / "inputs" / filename).read_bytes() == (Path(reference_inputs) / "inputs" / filename).read_bytes()
            reference_count = reference_count + 1
        rows = core.read_completed_case(job, output)
        assert rows is not None and len(rows) == 26, job["case_id"]
        records.append(rows)
        if job["full_result_reused"]:
            complete_imports = complete_imports + 1
        else:
            pairwise_imports = pairwise_imports + 1
        for suffix in ["profile", "peaks"]:
            original_lag = source_dir / "lag" / f"{job['source_case_id']}_{suffix}.csv"
            copied_lag = output / "lag" / f"{job['case_id']}_{suffix}.csv"
            if original_lag.exists():
                assert copied_lag.read_bytes() == original_lag.read_bytes()
                lag_tables = lag_tables + 1
        for analysis, method in ANALYSES:
            path = core.native_path(output, job["case_id"], analysis, method)
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                payload = json.load(handle)
            assert set(payload) >= {"edges", "raw", "diagnostics"}
            serialized_arrays = serialized_arrays + check_serialized_arrays(payload)
            native_edges = frame_from_snapshot(payload["edges"])
            native_label = {"PCMCI+": "pcmci+", "LKIF": "lkif", "GC": "granger", "CCM": "ccm"}[method]
            assert native_edges.loc[native_edges.status.eq("ok"), "method"].eq(native_label).all()
            assert native_edges.method.isin([native_label, method]).all()
            selected = rows.loc[rows.analysis.eq(analysis) & rows.method.eq(method)].set_index(["source", "target"])
            assert set(selected.index) == set(expected_edges(analysis))
            assert set(zip(native_edges.source, native_edges.target)) == set(expected_edges(analysis))
            for edge in native_edges.to_dict("records"):
                saved = selected.loc[edge["source"], edge["target"]]
                for field in native_edges.columns:
                    if field in {"source", "target", "method"} or field not in saved.index:
                        continue
                    assert same_scalar(edge[field], saved[field]), (path.name, field)
            if job["full_result_reused"] or analysis == "pairwise":
                source_case = job["source_case_id"] if job["full_result_reused"] else job["pairwise_source_case"]
                original = source.native_path(source_dir, source_case, analysis, method)
                assert path.read_bytes() == original.read_bytes(), path.name
                source_copies = source_copies + 1
            else:
                new_native = new_native + 1
            if not selected.status.eq("ok").all():
                failures.append({"case_id": job["case_id"], "analysis": analysis,
                                 "method": method, "statuses": selected.status.tolist()})
                continue
            raw = payload["raw"]
            if method == "CCM":
                assert raw["lag_rho"]["shape"] == [100, 31, 2]
                assert raw["rho_null"]["shape"] == [100, 2]
                assert len(raw["lag_library_seeds"]["values"]) == 100
                ccm_arrays = ccm_arrays + 1
            elif method == "PCMCI+":
                dimension = 2 if analysis == "pairwise" else 3
                for field in ["p_matrix", "val_matrix", "graph", "p_holm", "graph_holm"]:
                    assert raw[field]["shape"] == [dimension, dimension, 6]
            elif method == "LKIF":
                dimension = 2 if analysis == "pairwise" else 3
                for field in ["p", "nIF", "p_holm"]:
                    assert raw[field]["shape"] == [dimension, dimension, 5]
            elif method == "GC":
                fits = raw if analysis == "pairwise" else raw["directions"]
                for first, second in expected_edges(analysis):
                    for model in ["restricted", "full"]:
                        fit = fits[f"{first}_to_{second}"][model]
                        assert {"params", "resid", "fittedvalues", "cov_params", "endog", "exog",
                                "df_resid", "llf", "ssr"}.issubset(fit)
                        assert fit["resid"]["shape"] == [1400]
                        assert fit["endog"]["shape"] == [1400]
                        assert fit["exog"]["shape"][0] == 1400
                        gc_fits = gc_fits + 1
    assert complete_imports == 192 and pairwise_imports == 540
    assert source_copies == 3504 and new_native == 1620
    assert len(list((output / "members").glob("*_edges.csv"))) == 732
    assert len(list((output / "native").glob("*.json.gz"))) == 5124
    aggregate = pd.concat(records, ignore_index=True)
    assert len(aggregate) == 19032
    saved = pd.read_csv(output / "member_results.csv", float_precision="round_trip",
                        dtype={"edge_detected": "boolean", "unoriented": "boolean", "raw_edge_detected": "boolean"})
    pd.testing.assert_frame_equal(aggregate.reset_index(drop=True), saved.reset_index(drop=True),
                                  check_dtype=False, check_exact=True)
    tables = core.summarize_results(aggregate)
    assert set(tables) == {"xy_case_states", "cell_summary", "paired_score_changes"}
    summaries = []
    for name, frame in tables.items():
        bool_columns = [column for column in frame if pd.api.types.is_bool_dtype(frame[column])
                        or (len(frame[column].dropna()) and frame[column].dropna().map(
                            lambda value: isinstance(value, (bool, np.bool_))).all())]
        expected = frame.astype({column: "boolean" for column in bool_columns})
        saved_table = pd.read_csv(output / f"{name}.csv", float_precision="round_trip",
                                 dtype={column: "boolean" for column in bool_columns})
        pd.testing.assert_frame_equal(expected.reset_index(drop=True), saved_table.reset_index(drop=True),
                                      check_dtype=False, check_exact=True)
        summaries.append(name)
    assert not (output / "conditional_graph_states.csv").exists()
    state_rows = tables["xy_case_states"]
    assert len(state_rows) == 5124
    assert len(tables["cell_summary"]) == 588
    assert len(tables["paired_score_changes"]) == 10080
    assert tables["paired_score_changes"].groupby("age_type").size().eq(2520).all()
    pairwise_changes = tables["paired_score_changes"].query("analysis == 'pairwise'")
    assert len(pairwise_changes) == 5760
    assert pairwise_changes.loc[pairwise_changes.usable_score_pair, "score_change"].eq(0).all()
    assert set(state_rows.state) <= {"unavailable", "unresolved", "both", "x_only", "y_only", "neither"}
    for row in state_rows.itertuples():
        expected = row.state == ("x_only" if row.structure == "h1_direct" else "neither")
        assert expected == bool(row.correct_xy)
    return {"state": "passed", "checked_utc": datetime.now(timezone.utc).isoformat(),
            "inputs": inputs, "inputs_exactly_matching_preflight": reference_count,
            "committed_cases": len(jobs), "directional_rows": len(aggregate),
            "native_archives": source_copies + new_native, "exact_copied_native_archives": source_copies,
            "new_conditional_archives": new_native, "complete_imported_cases": complete_imports,
            "pairwise_imported_new_type_cases": pairwise_imports,
            "exact_copied_ccm_lag_tables": lag_tables,
            "ccm_complete_arrays_checked": ccm_arrays, "gc_complete_fit_records_checked": gc_fits,
            "serialized_native_arrays_checked": serialized_arrays,
            "retained_scientific_failures": failures, "directional_status_counts": aggregate.status.value_counts().to_dict(),
            "rebuilt_summary_tables": summaries, "xy_states_verified": len(state_rows),
            "pairwise_score_changes": len(pairwise_changes),
            "usable_pairwise_score_changes_exactly_zero": int(pairwise_changes.usable_score_pair.sum()),
            "xy_state_counts": state_rows.state.value_counts().to_dict(),
            "unresolved_xy_states_retained": int(state_rows.unresolved.sum()),
            "causal_fits_rerun": False, "chronologies_reconstructed_in_memory_for_verification": True}


def main():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output = ROOT / "data/benchmarks/driver_error_models_checks" / stamp
    output.mkdir(parents=True, exist_ok=False)
    prepared = output / "prepared"
    report_path = output / "report.json"
    started = time.perf_counter()
    report = {"state": "running", "started_utc": stamp, "output_directory": str(output),
              "production_run_started": False}
    def save_report():
        report["elapsed_seconds"] = time.perf_counter() - started
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Driver-error-model verification: {report_path}", flush=True)
    save_report()
    try:
        jobs, _, _ = prepare_without_new_processes(prepared)
        report["inputs"] = check_inputs(jobs, prepared)
        save_report()
        report["retention_and_guards"] = check_guards(jobs, prepared, output)
        report["summaries"] = check_summaries(jobs)
        save_report()
        report["native"] = check_native(jobs, prepared, output, report, save_report)
        report["validated_calculation_digest"] = core.calculation_digest()
        report["validated_source_calculation_digest"] = source.calculation_digest()
        report["validated_cross_calculation_digest"] = cross.calculation_digest()
        report["state"] = "passed"
    except BaseException as error:
        report["state"] = "failed"
        report["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        save_report()
    print(f"Driver-error-model checks passed: {report_path}", flush=True)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
