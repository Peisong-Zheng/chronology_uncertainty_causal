"""Validate the toolbox CCM migration against the timestamped source freeze.

Run from the project root with the production Python environment. This small
acceptance run compares full scientific outputs on five fixed local inputs,
then checks failures, zero distances and two Windows-spawned workers. It never
launches a production ensemble or overwrites prior experiment outputs.
"""
import os

for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"]:
    os.environ[name] = "1"

from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType
import ast
import gzip
import hashlib
import importlib
import json
import multiprocessing
import pickle
import platform
import sys
import time
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from toolbox import causal


# The pre-migration source freeze is the independent implementation reference.
FREEZE_DIR = ROOT / "data/archives/code_freeze_20261107T200542517660Z"
METHOD_SEED = 20261309       # Identical original/null/library random streams in both paths.
E = 3                      # Embedding dimension.
TAU = -1                   # Historical embedding spacing, 20-year steps.
EXCLUSION_RADIUS = 10      # Temporal neighbor exclusion, 200 years.
SAMPLES = 10               # Random libraries averaged for original/null Tp=0 skills.
N_SURROGATES = 100          # Cut-and-swap null records per direction.
N_LIBRARIES = 100           # Paired lag libraries.
LAG_STEPS = list(range(-15, 16))  # All 31 target shifts, +/-300 years.
ALPHA = 0.05               # Original CCM significance threshold.
SIGN_THRESHOLD = 0.90      # Existing lag-screen sign rule.
CASES = [
    ("linear_uniform_L300", "02_check_causal_pipeline/r0_linear_original_uniform_series.csv", 300),
    ("nonlinear_uniform_L300", "02_check_causal_pipeline/r0_nonlinear_original_uniform_series.csv", 300),
    ("independent_uniform_L300", "02_check_causal_pipeline/r0_independent_original_uniform_series.csv", 300),
    ("nonlinear_level4_L300", "02_3_check_uncertainty_escalation/r0_nonlinear_level4_C-J_m0_series.csv", 300),
    ("independent_level4_L1119", "02_3_check_uncertainty_escalation/r0_independent_level4_E-L_m0_series.csv", 1119),
]


def ccm_settings(library_size):
    return dict(E=E, tau=TAU, exclusion_radius=EXCLUSION_RADIUS, lib_size=library_size,
                sample=SAMPLES, n_surrogates=N_SURROGATES, seed=METHOD_SEED,
                alpha=ALPHA, lag_steps=LAG_STEPS, n_libraries=N_LIBRARIES,
                sign_threshold=SIGN_THRESHOLD)


def check_notebook_reference(frozen, current):
    """Allow the approved plot-only revision, keeping all scientific cells fixed."""
    assert len(frozen["cells"]) == len(current["cells"])
    for index, (before, after) in enumerate(zip(frozen["cells"], current["cells"])):
        assert before["cell_type"] == after["cell_type"], index
        if before["cell_type"] != "code":
            assert before["source"] == after["source"], index
            continue
        if index == 14:
            # This is the final figure-only cell, after the scientific results.
            # The independent CCM implementation lives in earlier checked cells.
            assert not any(isinstance(node, (ast.FunctionDef, ast.ClassDef))
                           for node in ast.walk(ast.parse("".join(after["source"]))))
            continue
        normalized = []
        for cell in [before, after]:
            tree = ast.parse("".join(cell["source"]))
            tree.body = [node for node in tree.body if not (
                isinstance(node, ast.ImportFrom) and node.module == "toolbox.plotting"
                or isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
                and isinstance(node.value.func, ast.Name) and node.value.func.id == "set_plot_style"
                and not node.value.args and not node.value.keywords)]
            normalized.append(ast.dump(tree, include_attributes=False))
        assert normalized[0] == normalized[1], f"02_7 scientific cell {index} changed"


def load_frozen_causal():
    manifest = json.loads((FREEZE_DIR / "manifest.json").read_text(encoding="utf-8"))
    archive = FREEZE_DIR / "source.zip"
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == manifest["archive_sha256"]
    with zipfile.ZipFile(archive) as zipped:
        source = zipped.read("toolbox/causal.py")
        # The 2026-11-08 style revision may alter imports and the last plot cell.
        check_notebook_reference(json.loads(zipped.read("02_7_check_ccm_neighbor_reuse.ipynb")),
            json.loads((ROOT / "02_7_check_ccm_neighbor_reuse.ipynb").read_text(encoding="utf-8")))
    frozen = ModuleType("frozen_causal_reference")
    frozen.__file__ = str(archive) + "/toolbox/causal.py"
    exec(compile(source, frozen.__file__, "exec"), frozen.__dict__)
    return frozen, manifest


def compare_complete(expected, actual, path="result", counts=None):
    """Require exact complete outputs, including dtype, NaNs and failure metadata."""
    if counts is None:
        counts = {"dataframes": 0, "arrays": 0, "scalar_leaves": 0}
    if isinstance(expected, pd.DataFrame):
        pd.testing.assert_frame_equal(expected, actual, check_exact=True)
        counts["dataframes"] = counts["dataframes"] + 1
    elif isinstance(expected, np.ndarray):
        assert expected.dtype == actual.dtype and expected.shape == actual.shape, path
        np.testing.assert_array_equal(expected, actual, err_msg=path)
        counts["arrays"] = counts["arrays"] + 1
    elif isinstance(expected, dict):
        assert list(expected) == list(actual), path
        for key in expected:
            compare_complete(expected[key], actual[key], f"{path}.{key}", counts)
    elif isinstance(expected, (list, tuple)):
        assert type(expected) is type(actual) and len(expected) == len(actual), path
        for index, (left, right) in enumerate(zip(expected, actual)):
            compare_complete(left, right, f"{path}[{index}]", counts)
    else:
        if isinstance(expected, (float, np.floating)) and np.isnan(expected):
            assert np.isnan(actual), path
        else:
            assert expected == actual, (path, expected, actual)
        counts["scalar_leaves"] = counts["scalar_leaves"] + 1
    return counts


def save_result(folder, name, result):
    with gzip.open(folder / f"{name}.pkl.gz", "wb") as handle:
        pickle.dump(result, handle, protocol=5)


def spawned_fit(case):
    """Module-scope worker; no closure or cache is passed between processes."""
    name, relative_path, library_size = case
    frame = pd.read_csv(ROOT / "data/processed" / relative_path)
    with threadpool_limits(limits=1):
        result = causal.run_ccm(frame, **ccm_settings(library_size))
    return name, result


def check_failure_controls(frozen, frame):
    """Use reduced test-only null counts for deterministic implementation controls."""
    native = importlib.import_module("pyEDM.CCM")
    original_tree = native.KDTree
    original_project = native.CCM.Project
    base_settings = ccm_settings(300)
    base_settings.update(lag_steps=None, sample=1, n_surrogates=2)
    original = frozen.run_ccm(frame, **base_settings)
    settings = dict(lag_steps=[-15, 0, 15], n_libraries=3,
                    sign_threshold=SIGN_THRESHOLD, seed=METHOD_SEED)
    seeds = np.random.default_rng(METHOD_SEED).integers(1, 2**31 - 1, 3)
    failed_seed = int(seeds[1])

    def fail_one_project(self, *args, **kwargs):
        seed = kwargs.get("seed", self.seed)
        if int(self.Tp) == 0 and int(seed) == failed_seed:
            raise RuntimeError("Injected one-library/one-lag failure")
        return original_project(self, *args, **kwargs)

    try:
        native.CCM.Project = fail_one_project
        reference = frozen.filter_ccm_by_lag(frame, original, **settings)
        candidate = causal.filter_ccm_by_lag(frame, original, **settings)
    finally:
        native.CCM.Project = original_project
    compare_complete(reference, candidate)
    errors = candidate["diagnostics"]["lag_filter"]["scan_errors"]
    assert len(errors) == 1 and errors[0]["library"] == 1 and errors[0]["lag_steps"] == 0
    assert np.isnan(candidate["raw"]["lag_rho"][1, 1]).all()
    assert np.isfinite(candidate["raw"]["lag_rho"][2, 1]).all()
    compare_complete(original["raw"]["rho_null"], candidate["raw"]["rho_null"])
    pd.testing.assert_series_equal(original["edges"]["score"], candidate["edges"]["score"], check_exact=True)
    pd.testing.assert_series_equal(original["edges"]["p_value"], candidate["edges"]["p_value"], check_exact=True)
    assert native.KDTree is original_tree

    # A native-valid Tp=0 library can be too large for the common lag anchors.
    oversized = frozen.run_ccm(frame, **dict(base_settings, lib_size=1390))
    compare_complete(frozen.filter_ccm_by_lag(frame, oversized, **settings),
                     causal.filter_ccm_by_lag(frame, oversized, **settings))

    repeated = frame.copy()
    repeated["x"] = np.resize(np.array([0., 1., 0., -1.]), len(frame))
    repeated["y"] = np.resize(np.array([1., 0., -1., 0.]), len(frame))
    repeated_original = frozen.run_ccm(repeated, **base_settings)
    compare_complete(frozen.filter_ccm_by_lag(repeated, repeated_original, **settings),
                     causal.filter_ccm_by_lag(repeated, repeated_original, **settings))

    def abort_tree(library):
        raise KeyboardInterrupt("Injected scan interruption")

    scan_settings = dict(E=E, tau=TAU, exclusion_radius=EXCLUSION_RADIUS,
                         library_size=300, lag_steps=settings["lag_steps"], library_seeds=seeds)
    try:
        native.KDTree = abort_tree
        try:
            causal.scan_ccm_reusing_neighbors(frame, **scan_settings)
        except KeyboardInterrupt as error:
            assert str(error) == "Injected scan interruption"
        else:
            raise AssertionError("Interruption was swallowed.")
        assert native.KDTree is abort_tree
    finally:
        native.KDTree = original_tree
    recovered = causal.scan_ccm_reusing_neighbors(frame, **scan_settings)
    assert not recovered["scan_errors"] and np.isfinite(recovered["rho"]).all()
    assert native.KDTree is original_tree
    assert recovered["cache_stats"]["native_searches"] == 2 * len(seeds)
    assert recovered["cache_stats"]["peak_cache_entries"] == 2
    return {"isolated_failure_and_recovery": True, "original_scores_and_nulls_retained": True,
            "oversized_library_failure_retained": True, "repeated_states_exact": True,
            "interrupt_restoration": True, "post_interrupt_native_recovery": True,
            "bounded_cache": recovered["cache_stats"]}


def main():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    folder = ROOT / "data/benchmarks/ccm_migration_checks" / stamp
    folder.mkdir(parents=True, exist_ok=False)
    frozen, manifest = load_frozen_causal()
    report = {"source_freeze": str(FREEZE_DIR.relative_to(ROOT)),
              "archive_sha256": manifest["archive_sha256"], "notebook_02_7_scientific_code_unchanged": True,
              "platform": platform.platform(), "cases": [], "passed": False}
    current_results = {}
    original_tree = importlib.import_module("pyEDM.CCM").KDTree
    started_all = time.perf_counter()
    try:
        with threadpool_limits(limits=1):
            for case_name, relative_path, library_size in CASES:
                path = ROOT / "data/processed" / relative_path
                frame = pd.read_csv(path)
                runs, elapsed = {}, {}
                for label, module in [("frozen", frozen), ("migrated", causal)]:
                    started = time.perf_counter()
                    runs[label] = module.run_ccm(frame, **ccm_settings(library_size))
                    elapsed[label] = time.perf_counter() - started
                    save_result(folder, f"{case_name}_{label}", runs[label])
                    print(f"{case_name} {label}: {elapsed[label]:.2f}s", flush=True)
                counts = compare_complete(runs["frozen"], runs["migrated"])
                assert np.isfinite(runs["migrated"]["raw"]["lag_rho"]).all()
                assert not runs["migrated"]["diagnostics"]["lag_filter"]["scan_errors"]
                assert importlib.import_module("pyEDM.CCM").KDTree is original_tree
                current_results[case_name] = runs["migrated"]
                report["cases"].append({"case": case_name, "library_size": library_size,
                                        "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                                        "complete_exact": True, "comparison_counts": counts,
                                        "seconds": elapsed})
                print(f"{case_name}: complete result EXACT", flush=True)
                (folder / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
            frame = pd.read_csv(ROOT / "data/processed" / CASES[0][1])
            report["controls"] = check_failure_controls(frozen, frame)
            print("Failure, zero-distance, cache and interruption controls passed.", flush=True)

        spawn_cases = [CASES[0], CASES[3]]
        with ProcessPoolExecutor(max_workers=2, mp_context=multiprocessing.get_context("spawn")) as pool:
            spawned = list(pool.map(spawned_fit, spawn_cases))
        report["spawned_cases"] = []
        for name, result in spawned:
            save_result(folder, f"{name}_spawned", result)
            counts = compare_complete(current_results[name], result)
            report["spawned_cases"].append({"case": name, "complete_exact": True,
                                            "comparison_counts": counts})
        report["passed"] = True
    finally:
        report["total_seconds"] = time.perf_counter() - started_all
        (folder / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"PASS: {folder / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
