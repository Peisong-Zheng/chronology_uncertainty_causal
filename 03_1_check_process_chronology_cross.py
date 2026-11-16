"""Cross saved process histories with saved Z chronologies, inspecting X-Y only.

The two fixed linear structures, three histories and three chronology schemes
come directly from completed experiment 03. Only off-diagonal conditional fits
are new. No process, sampling mask, observed chronology or ensemble is generated.
See agent/03_1_process_chronology_cross_plan.md.
"""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
for thread_name in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"]:
    os.environ[thread_name] = "1"
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "data/cache/matplotlib"))

from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import ast
import gzip
import hashlib
import importlib.util
import json
import multiprocessing
import shutil
import sys
import time

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
from toolbox.plotting import add_panel_labels, save_figure

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
source_spec = importlib.util.spec_from_file_location("_common_driver_source_03", ROOT / "03_check_common_driver.py")
base = importlib.util.module_from_spec(source_spec)
sys.modules[source_spec.name] = base
source_spec.loader.exec_module(base)

RUN_CAUSAL = True             # False reads committed outputs without new fits or imports.
PREPARE_ONLY = False          # True prepares/validates inputs and stops before case work.
RESUME = True                 # Retain committed results, including inherited fit failures.
N_WORKERS = 32                # Outer Windows processes; one numerical thread per worker.
PROCESS_SEEDS = [20261009, 20261301, 20261302]  # Existing process realization identifiers.
CHRONOLOGY_SEEDS = PROCESS_SEEDS.copy()       # Existing Z schemes, crossed with each history.
STRUCTURES = ["h0_common_driver", "h1_direct"]  # Without/with the existing direct X-to-Y link.
CONDITIONS = ["original_uniform", "irregular_accurate", "z_only"]
N_AGE_MEMBERS = 10            # All ten saved members around each observed Z chronology.
METHODS = ["PCMCI+", "LKIF", "GC", "CCM"]
ANALYSES = {"pairwise": METHODS, "conditional": ["PCMCI+", "LKIF", "GC"]}
SUPPORT_THRESHOLD = 0.90      # Descriptive within-cell directional support threshold.
MIN_VALID_FRACTION = 0.90     # Required fraction of jointly usable requested X/Y decisions.
IMPLEMENTATION_REVISION = "saved_process_z_clock_cross_xy_v1"
STEM = Path(__file__).stem
SOURCE_DIR = ROOT / "data/processed/03_check_common_driver"
DATA_DIR = ROOT / "data/processed" / STEM
FIGURE_DIR = ROOT / "figures" / STEM

# Use the original native calls and serialization unchanged; plotting stays here.
call_method = base.call_method
gpdc_null_path = base.gpdc_null_path
native_path = base.native_path
native_snapshot = base.native_snapshot
write_table = base.write_table
write_json = base.write_json
file_digest = base.file_digest
replace_output = base.replace_output
save_input = base.save_input


def calculation_digest():
    """Keep scientific code guarded while allowing presentation-only revisions."""
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    controls = {"RUN_CAUSAL", "PREPARE_ONLY", "RESUME", "N_WORKERS"}
    retained = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id in controls for t in node.targets):
            continue
        if isinstance(node, ast.FunctionDef) and (node.name.startswith("_plot_") or node.name == "plot_results"):
            continue
        if isinstance(node, ast.ImportFrom) and ((node.module or "").startswith("matplotlib") or node.module == "toolbox.plotting"):
            continue
        if isinstance(node, ast.Import) and all(a.name.startswith("matplotlib") for a in node.names):
            continue
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
            continue
        retained.append(node)
    tree.body = retained
    return hashlib.sha256(ast.dump(tree, include_attributes=False).encode()).hexdigest()


def source_case_paths(case_id, source_dir=SOURCE_DIR):
    """Files that define one complete imported result, including optional CCM tables."""
    source_dir = Path(source_dir)
    paths = [source_dir / "inputs" / f"{case_id}.csv", source_dir / "members" / f"{case_id}_edges.csv"]
    paths.extend(native_path(source_dir, case_id, analysis, method)
                 for analysis, methods in ANALYSES.items() for method in methods)
    paths.extend(path for suffix in ["profile", "peaks"]
                 if (path := source_dir / "lag" / f"{case_id}_{suffix}.csv").exists())
    return paths


def source_case_digest(case_id, source_dir=SOURCE_DIR):
    source_dir = Path(source_dir)
    hashes = {str(path.relative_to(source_dir)).replace("\\", "/"): file_digest(path)
              for path in source_case_paths(case_id, source_dir)}
    return hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()


def validate_source(source_dir=SOURCE_DIR):
    """Reject changed scientific settings/packages/source before importing any fit."""
    source_dir = Path(source_dir)
    current = base.execution_parameters()
    base.check_parameters(current, source_dir)
    saved = json.loads((source_dir / "parameters.json").read_text(encoding="utf-8"))
    if (saved["PROCESS_SEEDS"] != PROCESS_SEEDS or saved["STRUCTURES"] != STRUCTURES
            or saved["N_AGE_MEMBERS"] != N_AGE_MEMBERS or saved["ANALYSES"] != ANALYSES):
        raise ValueError("The cross design does not match completed experiment 03.")
    status = json.loads((source_dir / "run_status.json").read_text(encoding="utf-8"))
    if status["state"] != "complete" or status["completed_cases"] != status["requested_cases"]:
        raise ValueError("Source experiment 03 is not complete.")
    cases = pd.read_csv(source_dir / "cases.csv")
    if cases.case_id.duplicated().any():
        raise ValueError("Source case identifiers are duplicated.")
    required = {base.case_name(seed, structure, condition, member)
                for seed in PROCESS_SEEDS for structure in STRUCTURES for condition in CONDITIONS
                for member in ([-1] if condition != "z_only" else range(N_AGE_MEMBERS))}
    selected = cases[cases.case_id.isin(required)]
    if set(selected.case_id) != required or not selected.input_status.eq("ok").all():
        raise ValueError("Required source references/diagonals are missing or invalid.")
    for job in selected.to_dict("records"):
        if (job["estimator_seed"] != job["process_seed"] + saved["METHOD_SEED_OFFSET"]
                or base.read_completed_case(job, source_dir) is None):
            raise ValueError(f"Incomplete source case or altered estimator seed: {job['case_id']}")
    return saved, selected


def source_file_digests(source_dir=SOURCE_DIR):
    """Fingerprint only source artifacts used by this cross; never edit them."""
    source_dir = Path(source_dir)
    _, selected = validate_source(source_dir)
    paths = [source_dir / name for name in ["parameters.json", "cases.csv", "sampling.csv",
             "process_diagnostics.csv", "chronology_acceptance.csv"]]
    paths.extend(source_dir / "clocks" / f"s{seed}_z_C.npz" for seed in CHRONOLOGY_SEEDS)
    for seed in PROCESS_SEEDS:
        for structure in STRUCTURES:
            paths.extend(source_dir / "generated" / f"s{seed}_{structure}{suffix}"
                         for suffix in [".csv", "_truth.json"])
    for case_id in selected.case_id:
        paths.extend(source_case_paths(case_id, source_dir))
    return {str(path.relative_to(source_dir)).replace("\\", "/"): file_digest(path)
            for path in sorted(set(paths))}


def execution_parameters(source_dir=SOURCE_DIR):
    source_parameters, _ = validate_source(source_dir)
    controls = {"RUN_CAUSAL", "PREPARE_ONLY", "RESUME", "N_WORKERS"}
    return {"RUN_CAUSAL": RUN_CAUSAL, "PREPARE_ONLY": PREPARE_ONLY, "RESUME": RESUME,
        "N_WORKERS": N_WORKERS, "PROCESS_SEEDS": PROCESS_SEEDS, "CHRONOLOGY_SEEDS": CHRONOLOGY_SEEDS,
        "STRUCTURES": STRUCTURES, "CONDITIONS": CONDITIONS, "N_AGE_MEMBERS": N_AGE_MEMBERS,
        "ANALYSES": ANALYSES, "SUPPORT_THRESHOLD": SUPPORT_THRESHOLD,
        "MIN_VALID_FRACTION": MIN_VALID_FRACTION, "IMPLEMENTATION_REVISION": IMPLEMENTATION_REVISION,
        "source_parameters": {k: v for k, v in source_parameters.items() if k not in controls},
        "source_sha256": source_file_digests(source_dir), "calculation_sha256": calculation_digest(),
        "interpretation": "Three saved processes crossed with three saved Z schemes; ten conditional members per cell; X-Y outcomes only."}


def check_parameters(parameters, output_dir=DATA_DIR):
    output_dir = Path(output_dir)
    path = output_dir / "parameters.json"
    if path.exists():
        if RUN_CAUSAL and not RESUME:
            raise ValueError("Preserve existing outputs before a fresh cross experiment.")
        previous = json.loads(path.read_text(encoding="utf-8"))
        controls = {"RUN_CAUSAL", "PREPARE_ONLY", "RESUME", "N_WORKERS"}
        changed = [key for key in set(previous) | set(parameters)
                   if key not in controls and previous.get(key) != parameters.get(key)]
        if changed:
            raise ValueError(f"Saved settings, source artifacts or calculation changed: {sorted(changed)}")
    elif any((output_dir / "members").glob("*_edges.csv")):
        raise ValueError("Committed cross results have no scientific/source metadata.")


def copy_exact(source, destination):
    """Copy source bytes atomically; already-identical outputs retain their mtime."""
    source, destination = Path(source), Path(destination)
    digest = file_digest(source)
    if destination.exists() and file_digest(destination) == digest:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".partial")
    shutil.copyfile(source, temporary)
    if file_digest(temporary) != digest:
        raise ValueError(f"Source copy differs: {source.name}")
    replace_output(temporary, destination)


def case_name(process_seed, structure, condition, chronology_seed=0, member=-1):
    realization = PROCESS_SEEDS.index(process_seed) + 1
    scheme = CHRONOLOGY_SEEDS.index(chronology_seed) + 1 if chronology_seed else 0
    return f"r{realization}_a{scheme}_{structure}_{condition}_m{member}"


def prepare_inputs(output_dir=DATA_DIR, source_dir=SOURCE_DIR):
    """Prepare 12 references and 180 crossed inputs using saved arrays only."""
    output_dir, source_dir = Path(output_dir), Path(source_dir)
    source_parameters, source_cases = validate_source(source_dir)
    source_jobs = source_cases.set_index("case_id").to_dict("index")
    source_hashes = {case_id: source_case_digest(case_id, source_dir) for case_id in source_cases.case_id}
    sampling = pd.read_csv(source_dir / "sampling.csv", float_precision="round_trip")
    indices = sampling.native_index.to_numpy(dtype=int)
    age_true = sampling.age_true.to_numpy()
    if not np.array_equal(indices, sampling.native_index.to_numpy()) or np.any(np.diff(indices) <= 0):
        raise ValueError("Source sampling indices are not ordered integers.")
    copy_exact(source_dir / "sampling.csv", output_dir / "sampling.csv")
    clocks, clock_rows = {}, []
    grid = np.linspace(source_parameters["ANALYSIS_OLD"], source_parameters["ANALYSIS_YOUNG"],
        int(round((source_parameters["ANALYSIS_OLD"] - source_parameters["ANALYSIS_YOUNG"]) / source_parameters["DT"])) + 1)
    for scheme, seed in enumerate(CHRONOLOGY_SEEDS, start=1):
        path = source_dir / "clocks" / f"s{seed}_z_C.npz"
        with np.load(path, allow_pickle=False) as saved:
            clock = {name: saved[name].copy() for name in saved.files}
        if (not np.array_equal(clock["native_indices"], indices)
                or not np.array_equal(clock["age_true"], age_true)
                or clock["ages"].shape != (N_AGE_MEMBERS, len(indices))):
            raise ValueError("Source Z chronology and common sampling mask disagree.")
        copy_exact(path, output_dir / "clocks" / path.name)
        clocks[seed] = clock
        for member, ages in enumerate(clock["ages"]):
            if not np.isfinite(ages).all() or np.any(np.diff(ages) >= 0) or ages[0] < grid[0] or ages[-1] > grid[-1]:
                raise ValueError("A saved Z chronology violates order or analysis coverage.")
            error = np.interp(-grid, -age_true, ages - age_true) * 1000
            observed_error = np.interp(-grid, -age_true, clock["age_obs"] - age_true) * 1000
            rms = float(np.sqrt(np.mean(error**2)))
            clock_rows.append({"age_scheme_id": scheme, "chronology_seed": seed, "member": member,
                "record": "z", "condition": "z_only", "rms_error_years": rms,
                "mean_error_years": float(np.mean(error)), "observed_rms_years": float(np.sqrt(np.mean(observed_error**2))),
                "relative_x_y_rms_years": 0.0, "relative_x_z_rms_years": rms,
                "relative_y_z_rms_years": rms, "minimum_gap_years": float(np.min(-np.diff(ages))) * 1000,
                "covers_analysis_window": True})
    jobs = []
    for realization, seed in enumerate(PROCESS_SEEDS, start=1):
        for structure in STRUCTURES:
            full = pd.read_csv(source_dir / "generated" / f"s{seed}_{structure}.csv", float_precision="round_trip")
            if not np.array_equal(full.age.to_numpy()[indices], age_true):
                raise ValueError("Saved process and sampling ages disagree.")
            accurate_id = base.case_name(seed, structure, "irregular_accurate")
            accurate = pd.read_csv(source_dir / "inputs" / f"{accurate_id}.csv", float_precision="round_trip")
            accurate_native = pd.read_csv(source_dir / "inputs" / f"{accurate_id}.csv")
            if not np.array_equal(accurate.age.to_numpy(), grid):
                raise ValueError("Source reference uses a different analysis grid.")
            assignments = [(condition, 0, -1) for condition in CONDITIONS[:2]]
            assignments.extend(("z_only", chronology_seed, member)
                for chronology_seed in CHRONOLOGY_SEEDS for member in range(N_AGE_MEMBERS))
            for condition, chronology_seed, member in assignments:
                scheme = CHRONOLOGY_SEEDS.index(chronology_seed) + 1 if chronology_seed else 0
                identity = case_name(seed, structure, condition, chronology_seed, member)
                complete_reuse = condition != "z_only" or seed == chronology_seed
                source_id = base.case_name(seed, structure, condition, member) if complete_reuse else accurate_id
                if condition == "z_only":
                    series = accurate.copy()
                    series["z"] = np.interp(-grid, -clocks[chronology_seed]["ages"][member], full.z.to_numpy()[indices])
                    if complete_reuse:
                        original = pd.read_csv(source_dir / "inputs" / f"{source_id}.csv", float_precision="round_trip")
                        if not np.array_equal(series.to_numpy(), original.to_numpy()):
                            raise ValueError(f"Reconstructed diagonal differs from original 03: {identity}")
                else:
                    series = pd.read_csv(source_dir / "inputs" / f"{source_id}.csv", float_precision="round_trip")
                destination = output_dir / "inputs" / f"{identity}.csv"
                if destination.exists():
                    save_input(series, destination)
                if complete_reuse:
                    copy_exact(source_dir / "inputs" / f"{source_id}.csv", destination)
                else:
                    save_input(series, destination)
                actual_input = pd.read_csv(destination)
                reference = accurate_native if condition == "z_only" else pd.read_csv(source_dir / "inputs" / f"{source_id}.csv")
                if not np.array_equal(actual_input[["time", "age", "x", "y"]].to_numpy(), reference[["time", "age", "x", "y"]].to_numpy()):
                    raise ValueError("Default-parser X/Y estimator inputs changed during cross preparation.")
                job = {"case_id": identity, "process_seed": seed, "chronology_seed": chronology_seed,
                    "realization_id": realization, "age_scheme_id": scheme, "structure": structure,
                    "condition": condition, "member": member,
                    "estimator_seed": int(source_jobs[source_id]["estimator_seed"]),
                    "source_case_id": source_id, "pairwise_source_case": source_id if condition != "z_only" else accurate_id,
                    "input_origin": "source03_reference" if condition != "z_only" else "source03_diagonal" if complete_reuse else "saved_process_clock_cross",
                    "full_result_reused": complete_reuse, "input_status": "ok",
                    "source_case_sha256": source_hashes[source_id], "input_sha256": file_digest(destination)}
                jobs.append(job)
    write_table(pd.DataFrame(jobs), output_dir / "cases.csv")
    process_diagnostics = pd.read_csv(source_dir / "process_diagnostics.csv", float_precision="round_trip")
    process_diagnostics["realization_id"] = process_diagnostics.process_seed.map({seed: i + 1 for i, seed in enumerate(PROCESS_SEEDS)})
    clock_diagnostics = pd.DataFrame(clock_rows)
    write_table(process_diagnostics, output_dir / "process_diagnostics.csv")
    write_table(clock_diagnostics, output_dir / "clock_diagnostics.csv")
    return jobs, process_diagnostics, clock_diagnostics


def validate_case_source(job, output_dir, source_dir=SOURCE_DIR):
    """Check actual estimator inputs and complete source identity before reuse."""
    output_dir, source_dir = Path(output_dir), Path(source_dir)
    path = output_dir / "inputs" / f"{job['case_id']}.csv"
    if file_digest(path) != job["input_sha256"]:
        raise ValueError("Prepared input bytes changed.")
    source_id = job["source_case_id"]
    if source_case_digest(source_id, source_dir) != job["source_case_sha256"]:
        raise ValueError("Source case files changed after preparation.")
    rows = pd.read_csv(source_dir / "members" / f"{source_id}_edges.csv", float_precision="round_trip",
        dtype={"edge_detected": "boolean", "unoriented": "boolean", "raw_edge_detected": "boolean"})
    if (len(rows) != 26 or not rows.estimator_seed.eq(job["estimator_seed"]).all()
            or not rows.process_seed.eq(job["process_seed"]).all() or not rows.structure.eq(job["structure"]).all()
            or job["estimator_seed"] != job["process_seed"] + base.METHOD_SEED_OFFSET):
        raise ValueError("Source result identity or estimator seed changed.")
    series = pd.read_csv(path)
    original = pd.read_csv(source_dir / "inputs" / f"{source_id}.csv")
    columns = ["time", "age", "x", "y", "z"] if job["full_result_reused"] else ["time", "age", "x", "y"]
    if not np.array_equal(series[columns].to_numpy(), original[columns].to_numpy()):
        raise ValueError("Reused native result would have different numerical inputs.")
    if job["full_result_reused"] and file_digest(path) != file_digest(source_dir / "inputs" / f"{source_id}.csv"):
        raise ValueError("A fully imported case no longer has byte-identical source input.")
    return series, rows


def run_case(job, output_dir, null_file, source_dir=SOURCE_DIR):
    """Commit all 26 rows only after seven native archives, copied or newly fit."""
    output_dir, source_dir = Path(output_dir), Path(source_dir)
    series, source_rows = validate_case_source(job, output_dir, source_dir)
    tables = []
    with threadpool_limits(limits=1):
        for analysis, methods in ANALYSES.items():
            for method in methods:
                started = time.perf_counter()
                destination = native_path(output_dir, job["case_id"], analysis, method)
                reused = job["full_result_reused"] or analysis == "pairwise"
                if reused:
                    edges = source_rows[source_rows.analysis.eq(analysis) & source_rows.method.eq(method)].copy()
                    copy_exact(native_path(source_dir, job["source_case_id"], analysis, method), destination)
                    if method == "CCM":
                        for suffix in ["profile", "peaks"]:
                            source = source_dir / "lag" / f"{job['source_case_id']}_{suffix}.csv"
                            if source.exists():
                                copy_exact(source, output_dir / "lag" / f"{job['case_id']}_{suffix}.csv")
                    edges["source_seconds"] = edges["seconds"]
                else:
                    try:
                        result = call_method(series, analysis, method, job["estimator_seed"], null_file)
                    except Exception as error:
                        edges = base.empty_edges(analysis, method, "error", f"{type(error).__name__}: {error}")
                        result = {"edges": edges, "raw": {}, "diagnostics": {"error": str(error)}}
                    else:
                        edges = result["edges"].copy()
                        if "unoriented" not in edges:
                            edges["unoriented"] = False
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    temporary = destination.with_suffix(".partial.gz")
                    with gzip.open(temporary, "wt", encoding="utf-8") as handle:
                        json.dump(native_snapshot(result), handle, separators=(",", ":"))
                    replace_output(temporary, destination)
                    edges["source_seconds"] = np.nan
                edges["analysis"], edges["method"] = analysis, method
                edges["native_reused"] = bool(reused)
                edges["pairwise_reused"] = analysis == "pairwise"
                edges["result_origin"] = "source03_complete" if job["full_result_reused"] else "source03_pairwise" if reused else "new_conditional"
                edges["seconds"] = 0.0 if reused else time.perf_counter() - started
                for key, value in job.items():
                    edges[key] = value
                tables.append(edges)
    result = pd.concat(tables, ignore_index=True)
    write_table(result, output_dir / "members" / f"{job['case_id']}_edges.csv")
    return result


def read_completed_case(job, output_dir=DATA_DIR):
    path = Path(output_dir) / "inputs" / f"{job['case_id']}.csv"
    if not path.exists() or file_digest(path) != job["input_sha256"]:
        raise ValueError("Committed input bytes no longer match the prepared job.")
    return base.read_completed_case(job, output_dir)


def collect_results(jobs, completed):
    return base.collect_results(jobs, completed)


def summarize_results(results):
    """Inspect X/Y only; native Z directions remain archived without network scoring."""
    keys = ["case_id", "structure", "process_seed", "chronology_seed", "realization_id", "age_scheme_id",
            "condition", "member", "analysis", "method"]
    states = []
    xy = results[results.source.isin(["x", "y"]) & results.target.isin(["x", "y"])]
    for identity, pair in xy.groupby(keys, sort=False):
        meta = dict(zip(keys, identity))
        usable = len(pair) == 2 and pair.status.eq("ok").all() and pair.edge_detected.notna().all()
        unresolved = bool(pair.unoriented.fillna(False).any())
        forward = bool(pair.loc[pair.source.eq("x"), "edge_detected"].iloc[0]) if usable else None
        reverse = bool(pair.loc[pair.source.eq("y"), "edge_detected"].iloc[0]) if usable else None
        state = "unavailable" if not usable else "unresolved" if unresolved else (
            "both" if forward and reverse else "x_only" if forward else "y_only" if reverse else "neither")
        states.append({**meta, "usable": bool(usable), "forward": forward, "reverse": reverse,
            "unresolved": unresolved, "state": state,
            "correct_xy": bool(usable and not unresolved and forward == (meta["structure"] == "h1_direct") and not reverse),
            "failed": bool(not usable and not pair.status.eq("pending").any()),
            "pending": bool(pair.status.eq("pending").any())})
    states = pd.DataFrame(states)
    cell_keys = [key for key in keys if key not in {"case_id", "member"}]
    cells = []
    for identity, rows in states.groupby(cell_keys, sort=False):
        meta = dict(zip(cell_keys, identity))
        valid = rows[rows.usable]
        q_forward = float(valid.forward.mean()) if len(valid) else np.nan
        q_reverse = float(valid.reverse.mean()) if len(valid) else np.nan
        eligible = len(valid) >= MIN_VALID_FRACTION * len(rows)
        cells.append({**meta, "requested": len(rows), "usable": len(valid),
            "correct_xy": int(rows.correct_xy.sum()), "r_xy": float(rows.correct_xy.sum() / len(rows)),
            "forward": int(valid.forward.sum()), "reverse": int(valid.reverse.sum()),
            "unresolved": int(rows.unresolved.sum()), "failed": int(rows.failed.sum()), "pending": int(rows.pending.sum()),
            "q_x_to_y": q_forward, "q_y_to_x": q_reverse, "support_eligible": bool(eligible),
            "supported_x_to_y": bool(eligible and q_forward >= SUPPORT_THRESHOLD),
            "supported_y_to_x": bool(eligible and q_reverse >= SUPPORT_THRESHOLD)})
    reference = xy[xy.condition.eq("irregular_accurate")]
    reference_keys = ["structure", "process_seed", "analysis", "method", "source", "target"]
    changes = xy[xy.condition.eq("z_only")].merge(reference[reference_keys + ["score", "status"]],
        on=reference_keys, how="left", suffixes=("", "_reference"), validate="many_to_one")
    changes["usable_score_pair"] = changes.status.eq("ok") & changes.status_reference.eq("ok") & np.isfinite(changes.score) & np.isfinite(changes.score_reference)
    changes["score_change"] = np.where(changes.usable_score_pair, changes.score - changes.score_reference, np.nan)
    return {"xy_case_states": states, "cell_summary": pd.DataFrame(cells), "paired_score_changes": changes}


def write_status(state, completed, requested, started, output_dir=DATA_DIR, error=""):
    write_json({"state": state, "completed_cases": completed, "requested_cases": requested,
        "n_workers": N_WORKERS, "elapsed_seconds": time.perf_counter() - started,
        "updated_utc": datetime.now(timezone.utc).isoformat(), "error": error}, Path(output_dir) / "run_status.json")


def save_outputs(results, jobs, process_diagnostics, clock_diagnostics, output_dir=DATA_DIR, figure_dir=FIGURE_DIR):
    tables = summarize_results(results)
    write_table(results, Path(output_dir) / "member_results.csv")
    for name, table in tables.items():
        write_table(table, Path(output_dir) / f"{name}.csv")
    if "plot_results" in globals():
        plot_results(results, jobs, process_diagnostics, clock_diagnostics, Path(output_dir), Path(figure_dir))
    return tables


def main():
    parameters = execution_parameters()
    check_parameters(parameters)
    write_json(parameters, DATA_DIR / "parameters.json")
    started = time.perf_counter()
    jobs, completed = [], {}
    try:
        jobs, process_diagnostics, clock_diagnostics = prepare_inputs()
        for job in jobs:
            result = read_completed_case(job)
            if result is not None:
                completed[job["case_id"]] = result
        if PREPARE_ONLY:
            write_status("prepared", len(completed), len(jobs), started)
            print(f"Prepared {len(jobs)} saved-process/chronology cross inputs.", flush=True)
            return
        if RUN_CAUSAL:
            imports = [job for job in jobs if job["case_id"] not in completed and job["full_result_reused"]]
            for job in imports:
                result = run_case(job, DATA_DIR, None)
                completed[job["case_id"]] = result
            print(f"Complete source imports/resumed cases: {len(completed)}/{len(jobs)}.", flush=True)
            pending = [job for job in jobs if job["case_id"] not in completed]
            if pending:
                with threadpool_limits(limits=1):
                    effective_n = int(round((base.ANALYSIS_OLD - base.ANALYSIS_YOUNG) / base.DT)) + 1 - 2 * base.MAX_LAG
                    for seed in sorted({job["estimator_seed"] for job in pending}):
                        base.prepare_gpdc_null(gpdc_null_path(seed), sample_sizes=[effective_n],
                            sig_samples=base.GPDC_NULL_SAMPLES, seed=seed)
                write_status("running", len(completed), len(jobs), started)
                print(f"Starting {len(pending)} off-diagonal cases, three conditional methods each; workers={N_WORKERS}.", flush=True)
                context = multiprocessing.get_context("spawn")
                with ProcessPoolExecutor(max_workers=min(N_WORKERS, len(pending)), mp_context=context) as pool:
                    futures = [pool.submit(run_case, job, str(DATA_DIR), str(gpdc_null_path(job["estimator_seed"]))) for job in pending]
                    for future in as_completed(futures):
                        result = future.result()
                        completed[result.case_id.iloc[0]] = result
                        if len(completed) % 6 == 0 or len(completed) == len(jobs):
                            write_status("running", len(completed), len(jobs), started)
                            print(f"Completed {len(completed)}/{len(jobs)} crossed cases.", flush=True)
        results = collect_results(jobs, completed)
        save_outputs(results, jobs, process_diagnostics, clock_diagnostics)
        state = "complete" if len(completed) == len(jobs) else "partial"
        write_status(state, len(completed), len(jobs), started)
        print(f"Process/chronology cross {state}: {len(completed)}/{len(jobs)} cases.", flush=True)
    except BaseException as error:
        if jobs:
            write_table(collect_results(jobs, completed), DATA_DIR / "member_results.csv")
        write_status("interrupted_or_failed", len(completed), len(jobs), started, error=f"{type(error).__name__}: {error}")
        raise


def _plot_03_1_cell_data(states, jobs):
    """Count exact XY decisions against all requested inputs, including pending ones."""

    cases = pd.DataFrame(jobs)
    figure_rows = []
    structures = ["h0_common_driver", "h1_direct"]
    for structure in structures:
        for method in ["PCMCI+", "LKIF", "GC"]:
            for realization in [1, 2, 3]:
                for age_scheme in [0, 1, 2, 3]:
                    condition = "irregular_accurate" if age_scheme == 0 else "z_only"
                    requested = cases[(cases.structure == structure)
                                      & (cases.realization_id == realization)
                                      & (cases.age_scheme_id == age_scheme)
                                      & (cases.condition == condition)]
                    if len(requested) == 0:
                        raise ValueError(f"No requested inputs: {structure}/{realization}/{age_scheme}.")
                    selected = states[(states.analysis == "conditional") & (states.method == method)
                                      & states.case_id.isin(requested.case_id)]
                    if selected.case_id.duplicated().any():
                        raise ValueError("Duplicate conditional XY states in figure inputs.")
                    selected = selected.set_index("case_id").reindex(requested.case_id)
                    missing = selected.state.isna()
                    correct = selected.correct_xy.astype("boolean").fillna(False)
                    unresolved = selected.unresolved.astype("boolean").fillna(False)
                    failed = selected.failed.astype("boolean").fillna(False)
                    pending = selected.pending.astype("boolean").fillna(False) | missing
                    figure_rows.append({"structure": structure, "method": method,
                        "realization_id": realization, "age_scheme_id": age_scheme,
                        "condition": condition, "requested": len(requested),
                        "correct_xy": int(correct.sum()), "unresolved": int(unresolved.sum()),
                        "failed": int(failed.sum()), "pending": int(pending.sum()),
                        "recovery_percent": 100 * int(correct.sum()) / len(requested)})
    summary = pd.DataFrame(figure_rows)
    return {"cell_recovery": summary[summary.age_scheme_id > 0].reset_index(drop=True),
            "accurate_recovery": summary[summary.age_scheme_id == 0].reset_index(drop=True)}


def _plot_03_1_design(jobs, data_dir, output_dir):
    """Illustrate different trajectories and independently assigned Z-age schemes."""

    cases = pd.DataFrame(jobs)
    data_dir = Path(data_dir)
    fig, axes = plt.subplots(2, 2, figsize=(19, 12))
    ax = axes[0, 0]
    boxes = [(0.60, "Realizations 1, 2, 3", "Different innovations generate different X, Y and Z trajectories.\n"
              "The linear model and its memory settings stay fixed.", "#E9F3FA"),
             (0.21, "Age schemes 1, 2, 3", "Each has one observed Z-age curve and 10 ordered members.\n"
              "The whole scheme is reused across the three realizations.", "#F3EDF8")]
    for bottom, title, description, color in boxes:
        ax.add_patch(Rectangle((0.03, bottom), 0.94, 0.29, transform=ax.transAxes,
                               facecolor=color, edgecolor="0.65", linewidth=1))
        ax.text(0.5, bottom + 0.22, title, transform=ax.transAxes, ha="center", va="center",
                fontsize=14, fontweight="bold")
        ax.text(0.5, bottom + 0.105, description, transform=ax.transAxes, ha="center", va="center",
                fontsize=11, linespacing=1.5)
    ax.text(0.5, 0.55, "Cross every realization with every age scheme", transform=ax.transAxes,
            ha="center", va="center", fontsize=12)
    ax.text(0.5, 0.065, "Within each row: X and Y stay identical; only Z timing changes.",
            transform=ax.transAxes, ha="center", fontsize=11)
    ax.set_axis_off()
    ax.set_title("Two separate sources of variation", fontsize=14, pad=14)

    ax = axes[0, 1]
    realization_colors = ["#0072B2", "#D55E00", "#009E73"]
    for realization, color in zip([1, 2, 3], realization_colors):
        reference = cases[(cases.structure == "h0_common_driver")
                          & (cases.realization_id == realization)
                          & (cases.condition == "irregular_accurate")]
        if len(reference) != 1:
            raise ValueError("Expected one irregular accurate reference per realization.")
        series = pd.read_csv(data_dir / "inputs" / f"{reference.case_id.iloc[0]}.csv",
                             float_precision="round_trip")
        window = series.age.between(13, 17)
        ax.plot(series.loc[window, "age"], series.loc[window, "z"], color=color,
                linewidth=1.3, label=f"Realization {realization}")
    ax.set_xlim(17, 13)
    ax.set_xlabel("Age (kyr BP)")
    ax.set_ylabel("Z (unforced background SD units)")
    ax.set_title("Different Z trajectories under the same model\nAccurate ages; fixed 17–13 kyr BP illustration window", fontsize=13, pad=12)
    ax.legend(frameon=False, ncol=3, fontsize=10, loc="upper right")
    ax.grid(color="0.92", linewidth=0.6)

    ax = axes[1, 0]
    ax.imshow(np.eye(3), cmap=ListedColormap(["#F4F4F4", "#D8EAF6"]), vmin=0, vmax=1)
    for realization in range(3):
        for age_scheme in range(3):
            reused = realization == age_scheme
            label = "10 members\nReused\nfrom 03" if reused else "10 members\nNew\ncombination"
            ax.text(age_scheme, realization, label, ha="center", va="center", fontsize=11,
                    linespacing=1.6)
            if reused:
                ax.add_patch(Rectangle((age_scheme - 0.46, realization - 0.46), 0.92, 0.92,
                                       fill=False, edgecolor="0.15", linewidth=1.8, zorder=4))
    ax.set_xticks(range(3), [f"Age scheme\n{number}" for number in [1, 2, 3]], fontsize=11)
    ax.set_yticks(range(3), [f"Realization {number}" for number in [1, 2, 3]], fontsize=11)
    ax.set_xticks(np.arange(-0.5, 3, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 3, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="both", length=0)
    ax.set_title("The same 3 × 3 crossing in each system", fontsize=14, pad=14)

    ax = axes[1, 1]
    scheme_colors = ["#6A51A3", "#A6761D", "#1B9E77"]
    for age_scheme, color in zip([1, 2, 3], scheme_colors):
        scheme_rows = cases[(cases.age_scheme_id == age_scheme) & (cases.condition == "z_only")]
        seeds = scheme_rows.chronology_seed.unique()
        if len(seeds) != 1:
            raise ValueError("Each displayed age scheme must identify exactly one saved Z clock.")
        with np.load(data_dir / "clocks" / f"s{int(seeds[0])}_z_C.npz") as clock:
            true_age = clock["age_true"]
            observed_error = 1000 * (clock["age_obs"] - true_age)
            member_errors = 1000 * (clock["ages"] - true_age[None, :])
        window = (true_age >= 1) & (true_age <= 29)
        low, high = np.min(member_errors, axis=0), np.max(member_errors, axis=0)
        ax.fill_between(true_age[window], low[window], high[window], color=color, alpha=0.13,
                        linewidth=0)
        ax.plot(true_age[window], observed_error[window], color=color, linewidth=1.15,
                label=f"Age scheme {age_scheme}")
    ax.axhline(0, color="0.4", linewidth=0.8, zorder=0)
    ax.set_xlim(29, 1)
    ax.set_xlabel("True age (kyr BP)")
    ax.set_ylabel("Assigned age − true age (yr)")
    ax.set_title("Three complete Z-age schemes\nLine: observed centre; band: min–max of 10 members", fontsize=13, pad=12)
    ax.legend(frameon=False, ncol=3, fontsize=10, loc="upper right")
    ax.grid(color="0.92", linewidth=0.6)
    for ax in [axes[0, 1], axes[1, 1]]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("0. Separate simulated-history variation from Z-age uncertainty", fontsize=18, y=0.99)
    fig.text(0.5, 0.047, "Both systems have Z → X and Z → Y; the second also has X → Y. "
             "Only the complete X–Y decision is assessed here.", ha="center", fontsize=11)
    fig.text(0.5, 0.025, "Three histories × three age schemes are not nine independent histories. "
             "Pairwise methods, including CCM, do not change when only Z's dates change.", ha="center", fontsize=10.5)
    add_panel_labels(axes, x=-0.11, y=1.04)
    fig.tight_layout(rect=(0.025, 0.075, 0.99, 0.945), h_pad=3.4, w_pad=3.8)
    save_figure(fig, output_dir, "00_cross_design")
    plt.close(fig)


def _plot_03_1_recovery(cell_data, accurate_data, output_dir):
    """Show exact XY recovery in each realization-by-age-scheme cell."""

    systems = [("h0_common_driver", "Common driver only", "Correct: neither X–Y direction"),
               ("h1_direct", "Common driver + direct X → Y", "Correct: X → Y only")]
    fig, axes = plt.subplots(2, 3, figsize=(20, 11.8))
    for row, (structure, system_name, correct_label) in enumerate(systems):
        for column, method in enumerate(["PCMCI+", "LKIF", "GC"]):
            ax = axes[row, column]
            selected = cell_data[(cell_data.structure == structure) & (cell_data.method == method)]
            values = selected.pivot(index="realization_id", columns="age_scheme_id",
                                    values="recovery_percent").reindex(index=[1, 2, 3], columns=[1, 2, 3])
            if values.isna().any().any():
                raise ValueError("Incomplete requested 3 × 3 crossing in figure summary.")
            image = ax.imshow(values.to_numpy(), cmap="Blues", vmin=0, vmax=100)
            row_labels = []
            for realization in [1, 2, 3]:
                reference = accurate_data[(accurate_data.structure == structure)
                                          & (accurate_data.method == method)
                                          & (accurate_data.realization_id == realization)]
                if len(reference) != 1:
                    raise ValueError("Expected one accurate reference summary per method and realization.")
                reference = reference.iloc[0]
                baseline_label = f"Acc: {int(reference.correct_xy)}/{int(reference.requested)}"
                exceptions = [f"{label}:{int(reference[field])}" for field, label in
                              [("unresolved", "U"), ("failed", "F"), ("pending", "P")]
                              if reference[field] > 0]
                if exceptions:
                    baseline_label = baseline_label + "; " + ",".join(exceptions)
                row_labels.append(f"Realization {realization}\n{baseline_label}")
                for age_scheme in [1, 2, 3]:
                    cell = selected[(selected.realization_id == realization)
                                    & (selected.age_scheme_id == age_scheme)].iloc[0]
                    label = f"{int(cell.correct_xy)}/{int(cell.requested)}\nU: {int(cell.unresolved)}"
                    if cell.failed or cell.pending:
                        label = label + f"\nF: {int(cell.failed)}  P: {int(cell.pending)}"
                    ax.text(age_scheme - 1, realization - 1, label, ha="center", va="center",
                            fontsize=12, linespacing=1.5,
                            color="white" if cell.recovery_percent >= 60 else "#111111")
                    if realization == age_scheme:
                        ax.add_patch(Rectangle((age_scheme - 1.46, realization - 1.46), 0.92, 0.92,
                                               fill=False, edgecolor="0.2", linewidth=1.7))
            ax.set_xticks(range(3), [f"Age scheme {number}" for number in [1, 2, 3]], fontsize=10.5)
            ax.set_yticks(range(3), row_labels, fontsize=10.5)
            ax.set_xticks(np.arange(-0.5, 3, 1), minor=True)
            ax.set_yticks(np.arange(-0.5, 3, 1), minor=True)
            ax.grid(which="minor", color="white", linewidth=2)
            ax.tick_params(which="both", length=0)
            ax.set_title(f"{method}\n{system_name}\n{correct_label}", fontsize=12, pad=14)
    fig.suptitle("1. Does the complete X–Y decision depend on the realization or the Z-age scheme?", fontsize=17, y=0.99)
    fig.text(0.5, 0.945, "Conditional X,Y,Z analyses | Cells: correct / requested, followed by unresolved count U | "
             "Acc: own-realization irregular accurate-age reference", ha="center", fontsize=10.5)
    fig.subplots_adjust(left=0.105, right=0.90, bottom=0.14, top=0.855, wspace=0.55, hspace=0.53)
    colorbar_axis = fig.add_axes([0.935, 0.28, 0.015, 0.43])
    colorbar = fig.colorbar(image, cax=colorbar_axis)
    colorbar.set_ticks(np.arange(0, 101, 20))
    colorbar.set_ticklabels([f"{number}%" for number in range(0, 101, 20)])
    colorbar.set_label("Exact X–Y recovery / requested members", fontsize=11, labelpad=12)
    add_panel_labels(axes, x=-0.16, y=1.04)
    fig.text(0.5, 0.086, "Outlined diagonal cells reuse the corresponding 03 results. "
             "All other cells apply an existing age scheme to a different realization.", ha="center", fontsize=10.5)
    fig.text(0.5, 0.059, "Unresolved decisions are not correct; failures (F) and pending cases (P), where present, remain in the denominator. "
             "Saved states distinguish missed and extra arrows.", ha="center", fontsize=10.5)
    fig.text(0.5, 0.032, "Row and column patterns are descriptive: there are three simulated histories, three complete age schemes, "
             "and fixed estimator seeds within each realization.", ha="center", fontsize=10.5)
    save_figure(fig, output_dir, "01_xy_recovery_matrix")
    plt.close(fig)


def plot_results(results, jobs, process_diagnostics, clock_diagnostics, data_dir, output_dir):
    """Render two focused figures without changing or writing scientific results."""
    states = summarize_results(results)["xy_case_states"]
    figure_data = _plot_03_1_cell_data(states, jobs)
    _plot_03_1_design(jobs, data_dir, output_dir)
    _plot_03_1_recovery(figure_data["cell_recovery"], figure_data["accurate_recovery"], output_dir)
    return figure_data


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
