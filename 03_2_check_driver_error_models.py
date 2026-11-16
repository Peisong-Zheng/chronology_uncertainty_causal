"""Compare C/L/E/J uncertainty in the observed common driver's chronology.

X/Y retain accurate ages. Reuse all C-type and accurate results from 03_1;
generate only L/E/J Z chronologies, then fit the three conditional methods.
Dynamics, sampling, memories and estimator settings remain those of 03.
See agent/03_2_driver_error_models_plan.md.
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
from matplotlib.patches import Rectangle
from toolbox.plotting import add_panel_labels, save_figure
from toolbox.observation import make_age_uncertainty, simulate_observations, draw_age_ensemble

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
source_spec = importlib.util.spec_from_file_location("_driver_error_cross_source", ROOT / "03_1_check_process_chronology_cross.py")
cross = importlib.util.module_from_spec(source_spec)
sys.modules[source_spec.name] = cross
source_spec.loader.exec_module(cross)
base = cross.base

RUN_CAUSAL = True             # False reads committed outputs without new fits or imports.
PREPARE_ONLY = False          # True prepares/validates inputs and stops before case work.
RESUME = True                 # Retain committed results, including inherited fit failures.
N_WORKERS = 32                # Outer Windows processes; one numerical thread per worker.
PROCESS_SEEDS = [20261009, 20261301, 20261302]  # Existing process realization identifiers.
CHRONOLOGY_SEEDS = PROCESS_SEEDS.copy()       # Existing Z schemes, crossed with each history.
STRUCTURES = ["h0_common_driver", "h1_direct"]  # Without/with the existing direct X-to-Y link.
CONDITIONS = ["original_uniform", "irregular_accurate", "z_only"]
N_AGE_MEMBERS = 10            # Ordered members around each uncertain observed Z chronology.
AGE_TYPES = ["C", "L", "E", "J"]  # Existing uncertainty-envelope families.
TYPE_SCALES_YEARS = {"C": 40.0, "L": 40.0, "E": 40.0, "J": 10.0}  # Peak proposal SDs, not RMS matched.
AGE_EXP_SHAPE = 2.0           # Fixed E-envelope curvature, unchanged from 02.
JUMP_COUNT = 5               # Level-3 local intervals; 50% of the 28-kyr analysis window.
JUMP_INTERVAL_WIDTH = 2.8    # Full width of each compact J interval, kyr.
JUMP_CENTER_FRACTIONS = [0.5, 0.2, 0.8, 0.35, 0.65]  # Fractions of the 29-1 kyr window.
AGE_CORRELATION_LENGTH = 1 / 60  # Squared-exponential length in normalized sample order.
AGE_MAX_DRAWS = 20000        # Whole-curve proposal budget; no sorting or repair.
METHODS = ["PCMCI+", "LKIF", "GC", "CCM"]
ANALYSES = {"pairwise": METHODS, "conditional": ["PCMCI+", "LKIF", "GC"]}
SUPPORT_THRESHOLD = 0.90      # Descriptive within-cell directional support threshold.
MIN_VALID_FRACTION = 0.90     # Required fraction of jointly usable requested X/Y decisions.
IMPLEMENTATION_REVISION = "driver_error_types_saved_xy_v1"
STEM = Path(__file__).stem
SOURCE_DIR = ROOT / "data/processed/03_1_check_process_chronology_cross"
PROCESS_DIR = ROOT / "data/processed/03_check_common_driver"
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
    """Require the completed, unchanged 03_1 comparison and original 03 source."""
    source_dir = Path(source_dir)
    cross.check_parameters(cross.execution_parameters(PROCESS_DIR), source_dir)
    saved = json.loads((source_dir / "parameters.json").read_text(encoding="utf-8"))
    if (saved["PROCESS_SEEDS"] != PROCESS_SEEDS or saved["CHRONOLOGY_SEEDS"] != CHRONOLOGY_SEEDS
            or saved["STRUCTURES"] != STRUCTURES or saved["N_AGE_MEMBERS"] != N_AGE_MEMBERS
            or saved["ANALYSES"] != ANALYSES):
        raise ValueError("The type comparison does not match completed experiment 03_1.")
    status = json.loads((source_dir / "run_status.json").read_text(encoding="utf-8"))
    if status["state"] != "complete" or status["completed_cases"] != status["requested_cases"]:
        raise ValueError("Source experiment 03_1 is not complete.")
    cases = pd.read_csv(source_dir / "cases.csv")
    required = set()
    for seed in PROCESS_SEEDS:
        for structure in STRUCTURES:
            required.update(cross.case_name(seed, structure, condition) for condition in CONDITIONS[:2])
            required.update(cross.case_name(seed, structure, "z_only", chronology_seed, member)
                            for chronology_seed in CHRONOLOGY_SEEDS for member in range(N_AGE_MEMBERS))
    if cases.case_id.duplicated().any() or set(cases.case_id) != required or not cases.input_status.eq("ok").all():
        raise ValueError("Source 03_1 references or C cells are missing, duplicate or invalid.")
    for job in cases.to_dict("records"):
        if (job["estimator_seed"] != job["process_seed"] + base.METHOD_SEED_OFFSET
                or cross.read_completed_case(job, source_dir) is None):
            raise ValueError(f"Incomplete source case or changed seed: {job['case_id']}")
    return saved, cases


def source_file_digests(source_dir=SOURCE_DIR):
    """Fingerprint used source evidence without changing either source directory."""
    source_dir = Path(source_dir)
    _, cases = validate_source(source_dir)
    paths = [source_dir / name for name in ["parameters.json", "cases.csv", "sampling.csv",
             "process_diagnostics.csv", "clock_diagnostics.csv"]]
    paths.extend(source_dir / "clocks" / f"s{seed}_z_C.npz" for seed in CHRONOLOGY_SEEDS)
    for case_id in cases.case_id:
        paths.extend(source_case_paths(case_id, source_dir))
    hashes = {"cross/" + str(path.relative_to(source_dir)).replace("\\", "/"): file_digest(path)
              for path in sorted(set(paths))}
    for seed in PROCESS_SEEDS:
        for structure in STRUCTURES:
            for suffix in [".csv", "_truth.json"]:
                path = PROCESS_DIR / "generated" / f"s{seed}_{structure}{suffix}"
                hashes["process/" + str(path.relative_to(PROCESS_DIR)).replace("\\", "/")] = file_digest(path)
    return hashes


def execution_parameters(source_dir=SOURCE_DIR):
    source_parameters, _ = validate_source(source_dir)
    controls = {"RUN_CAUSAL", "PREPARE_ONLY", "RESUME", "N_WORKERS"}
    return {"RUN_CAUSAL": RUN_CAUSAL, "PREPARE_ONLY": PREPARE_ONLY, "RESUME": RESUME,
        "N_WORKERS": N_WORKERS, "PROCESS_SEEDS": PROCESS_SEEDS, "CHRONOLOGY_SEEDS": CHRONOLOGY_SEEDS,
        "STRUCTURES": STRUCTURES, "CONDITIONS": CONDITIONS, "N_AGE_MEMBERS": N_AGE_MEMBERS,
        "AGE_TYPES": AGE_TYPES, "TYPE_SCALES_YEARS": TYPE_SCALES_YEARS,
        "AGE_EXP_SHAPE": AGE_EXP_SHAPE, "JUMP_COUNT": JUMP_COUNT,
        "JUMP_INTERVAL_WIDTH": JUMP_INTERVAL_WIDTH, "JUMP_CENTER_FRACTIONS": JUMP_CENTER_FRACTIONS,
        "AGE_CORRELATION_LENGTH": AGE_CORRELATION_LENGTH, "AGE_MAX_DRAWS": AGE_MAX_DRAWS,
        "ANALYSES": ANALYSES, "SUPPORT_THRESHOLD": SUPPORT_THRESHOLD,
        "MIN_VALID_FRACTION": MIN_VALID_FRACTION, "IMPLEMENTATION_REVISION": IMPLEMENTATION_REVISION,
        "source_parameters": {k: v for k, v in source_parameters.items() if k not in controls},
        "source_sha256": source_file_digests(source_dir), "calculation_sha256": calculation_digest(),
        "age_range_kyr_bp": [0.0, base.DURATION],
        "scale_definition": "02 level-3 envelopes on the retained 0-30 kyr range; not RMS matched; J intervals in 29-1 kyr.",
        "pairing": "Same observed/ensemble RNG seeds across types; whole-curve rejection can break accepted-member pairing.",
        "interpretation": "Three saved realizations crossed with three schemes per error type; X/Y accurate; X-Y outcomes only."}


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


def case_name(process_seed, structure, condition, chronology_seed=0, member=-1, age_type="accurate"):
    return f"{cross.case_name(process_seed, structure, condition, chronology_seed, member)}_{age_type}"


def age_envelope(age, age_type):
    """Use the existing level-3 envelope on the full retained age range, in kyr."""
    centers = base.ANALYSIS_YOUNG + (base.ANALYSIS_OLD - base.ANALYSIS_YOUNG) * np.array(JUMP_CENTER_FRACTIONS[:JUMP_COUNT])
    intervals = [(float(center - JUMP_INTERVAL_WIDTH / 2), float(center + JUMP_INTERVAL_WIDTH / 2))
                 for center in centers]
    return make_age_uncertainty(age, kind=age_type, scale=TYPE_SCALES_YEARS[age_type] / 1000,
        age_range=(0.0, base.DURATION), exp_shape=AGE_EXP_SHAPE,
        jump_intervals=intervals if age_type == "J" else None)


def prepare_clocks(output_dir, source_dir, sampling, grid):
    """Generate nine new schemes around observed centers and verify the three C schemes."""
    indices = sampling.native_index.to_numpy(dtype=int)
    full_age = base.DURATION - np.arange(int(round(base.DURATION / base.DT)) + 1) * base.DT
    age_true = full_age[indices]
    if not np.array_equal(age_true, sampling.age_true.to_numpy()):
        raise ValueError("Saved sampling ages do not match the original native grid.")
    clocks, clock_rows, acceptance_rows = {}, [], []
    for scheme, seed in enumerate(CHRONOLOGY_SEEDS, start=1):
        age_seed = int(np.random.SeedSequence([seed, 3103, 2]).generate_state(1)[0])
        ensemble_seed = int(np.random.SeedSequence([seed, 3104, 2]).generate_state(1)[0])
        for age_type in AGE_TYPES:
            observation = simulate_observations(full_age, np.zeros(len(full_age)), indices,
                age_envelope(age_true, age_type), age_correlation_length=AGE_CORRELATION_LENGTH,
                measurement_noise_ratio=0, rng_age=age_seed, rng_measurement=0)
            observed = observation["observed"]
            ensemble = draw_age_ensemble(observed, n_members=N_AGE_MEMBERS,
                age_correlation_length=AGE_CORRELATION_LENGTH, max_draws=AGE_MAX_DRAWS, rng=ensemble_seed)
            if not ensemble["diagnostics"]["complete"]:
                raise RuntimeError(f"Incomplete chronology generation: {age_type}/scheme {scheme}")
            clock = {"age_true": age_true, "age_obs": observed.age_obs.to_numpy(),
                "age_std": observed.age_std.to_numpy(), "ages": ensemble["ages"],
                "sample_id": ensemble["sample_id"], "native_indices": indices}
            filename = f"s{seed}_z_{age_type}.npz"
            destination = Path(output_dir) / "clocks" / filename
            if age_type == "C":
                with np.load(Path(source_dir) / "clocks" / filename, allow_pickle=False) as saved:
                    if set(saved.files) != set(clock) or any(not np.array_equal(saved[key], clock[key]) for key in clock):
                        raise ValueError("Native C regeneration differs from saved 03_1 chronology.")
                if destination.exists():
                    base.save_clock(clock, destination)
                copy_exact(Path(source_dir) / "clocks" / filename, destination)
            else:
                base.save_clock(clock, destination)
            clocks[(seed, age_type)] = clock
            member_errors = []
            for ages in clock["ages"]:
                if (not np.isfinite(ages).all() or np.any(np.diff(ages) >= 0)
                        or ages[0] < grid[0] or ages[-1] > grid[-1]):
                    raise ValueError(f"Chronology order/coverage failed for {age_type}/scheme {scheme}; no extrapolation.")
                member_errors.append(np.interp(-grid, -age_true, ages - age_true) * 1000)
            member_errors = np.asarray(member_errors)
            observed_error = np.interp(-grid, -age_true, clock["age_obs"] - age_true) * 1000
            proposal_rms = float(np.sqrt(np.mean((1000 * age_envelope(grid, age_type))**2)))
            observed_rms = float(np.sqrt(np.mean(observed_error**2)))
            effective_sd = float(np.mean(np.std(member_errors, axis=0, ddof=0)))
            common = {"age_type": age_type, "age_scheme_id": scheme, "chronology_seed": seed,
                "record": "z", "condition": "z_only", "age_seed": age_seed, "ensemble_seed": ensemble_seed,
                "proposal_rms_years": proposal_rms, "observed_rms_years": observed_rms,
                "effective_member_sd_years": effective_sd, "analysis_old_kyr": base.ANALYSIS_OLD,
                "analysis_young_kyr": base.ANALYSIS_YOUNG,
                "effective_sd_definition": "Mean pointwise population SD of 10 age-error curves on true 29-1 kyr grid; ddof=0"}
            acceptance_rows.append({**common, **ensemble["diagnostics"],
                "source_c_verified": age_type == "C", "proposal_peak_sd_years": TYPE_SCALES_YEARS[age_type]})
            for member, (ages, error) in enumerate(zip(clock["ages"], member_errors)):
                rms = float(np.sqrt(np.mean(error**2)))
                clock_rows.append({**common, "member": member, "rms_error_years": rms,
                    "mean_error_years": float(np.mean(error)), "relative_x_y_rms_years": 0.0,
                    "relative_x_z_rms_years": rms, "relative_y_z_rms_years": rms,
                    "minimum_gap_years": float(np.min(-np.diff(ages))) * 1000, "covers_analysis_window": True,
                    "accepted_members": ensemble["diagnostics"]["accepted_members"],
                    "attempted_draws": ensemble["diagnostics"]["attempted_draws"],
                    "order_rejected": ensemble["diagnostics"]["order_rejected"],
                    "age_order_violations": ensemble["diagnostics"]["age_order_violations"],
                    "acceptance_rate": ensemble["diagnostics"]["acceptance_rate"]})
    write_table(pd.DataFrame(acceptance_rows), Path(output_dir) / "chronology_acceptance.csv")
    return clocks, pd.DataFrame(clock_rows)


def prepare_inputs(output_dir=DATA_DIR, source_dir=SOURCE_DIR):
    """Prepare 12 references and 720 Z-only inputs; only L/E/J clocks are new."""
    output_dir, source_dir = Path(output_dir), Path(source_dir)
    _, source_cases = validate_source(source_dir)
    source_jobs = source_cases.set_index("case_id").to_dict("index")
    source_hashes = {case_id: source_case_digest(case_id, source_dir) for case_id in source_cases.case_id}
    sampling = pd.read_csv(source_dir / "sampling.csv", float_precision="round_trip")
    indices = sampling.native_index.to_numpy(dtype=int)
    age_true = sampling.age_true.to_numpy()
    if not np.array_equal(indices, sampling.native_index.to_numpy()) or np.any(np.diff(indices) <= 0):
        raise ValueError("Source sampling indices are not ordered integers.")
    copy_exact(source_dir / "sampling.csv", output_dir / "sampling.csv")
    grid = np.linspace(base.ANALYSIS_OLD, base.ANALYSIS_YOUNG,
                       int(round((base.ANALYSIS_OLD - base.ANALYSIS_YOUNG) / base.DT)) + 1)
    clocks, clock_diagnostics = prepare_clocks(output_dir, source_dir, sampling, grid)
    jobs = []
    for realization, seed in enumerate(PROCESS_SEEDS, start=1):
        for structure in STRUCTURES:
            full = pd.read_csv(PROCESS_DIR / "generated" / f"s{seed}_{structure}.csv", float_precision="round_trip")
            if not np.array_equal(full.age.to_numpy()[indices], age_true):
                raise ValueError("Saved process and sampling ages disagree.")
            accurate_id = cross.case_name(seed, structure, "irregular_accurate")
            accurate = pd.read_csv(source_dir / "inputs" / f"{accurate_id}.csv", float_precision="round_trip")
            accurate_native = pd.read_csv(source_dir / "inputs" / f"{accurate_id}.csv")
            if not np.array_equal(accurate.age.to_numpy(), grid):
                raise ValueError("Source reference uses a different analysis grid.")
            assignments = [(condition, 0, -1, "accurate") for condition in CONDITIONS[:2]]
            assignments.extend(("z_only", chronology_seed, member, age_type)
                for age_type in AGE_TYPES for chronology_seed in CHRONOLOGY_SEEDS for member in range(N_AGE_MEMBERS))
            for condition, chronology_seed, member, age_type in assignments:
                scheme = CHRONOLOGY_SEEDS.index(chronology_seed) + 1 if chronology_seed else 0
                identity = case_name(seed, structure, condition, chronology_seed, member, age_type)
                complete_reuse = condition != "z_only" or age_type == "C"
                source_id = cross.case_name(seed, structure, condition, chronology_seed, member) if complete_reuse else accurate_id
                if condition == "z_only":
                    series = accurate.copy()
                    series["z"] = np.interp(-grid, -clocks[(chronology_seed, age_type)]["ages"][member], full.z.to_numpy()[indices])
                    if complete_reuse:
                        original = pd.read_csv(source_dir / "inputs" / f"{source_id}.csv", float_precision="round_trip")
                        if not np.array_equal(series.to_numpy(), original.to_numpy()):
                            raise ValueError(f"Reconstructed C input differs from original 03_1: {identity}")
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
                    raise ValueError("Default-parser X/Y estimator inputs changed during type preparation.")
                jobs.append({"case_id": identity, "process_seed": seed, "chronology_seed": chronology_seed,
                    "realization_id": realization, "age_scheme_id": scheme, "age_type": age_type,
                    "structure": structure, "condition": condition, "member": member,
                    "estimator_seed": int(source_jobs[source_id]["estimator_seed"]),
                    "source_case_id": source_id, "pairwise_source_case": source_id,
                    "input_origin": "source03_1_reference" if condition != "z_only" else "source03_1_C" if complete_reuse else "new_Z_error_type",
                    "full_result_reused": complete_reuse, "input_status": "ok",
                    "source_case_sha256": source_hashes[source_id], "input_sha256": file_digest(destination)})
    write_table(pd.DataFrame(jobs), output_dir / "cases.csv")
    process_diagnostics = pd.read_csv(source_dir / "process_diagnostics.csv", float_precision="round_trip")
    write_table(process_diagnostics, output_dir / "process_diagnostics.csv")
    write_table(clock_diagnostics, output_dir / "clock_diagnostics.csv")
    return jobs, process_diagnostics, clock_diagnostics


def validate_case_source(job, output_dir, source_dir=SOURCE_DIR):
    """Check actual estimator inputs and complete source identity before reuse."""
    output_dir, source_dir = Path(output_dir), Path(source_dir)
    if job["condition"] not in CONDITIONS or job["age_type"] not in AGE_TYPES + ["accurate"]:
        raise ValueError("Unknown chronology condition or error type.")
    uncertain = job["condition"] == "z_only"
    if uncertain != (job["age_type"] in AGE_TYPES):
        raise ValueError("Accurate and uncertain case identities disagree.")
    expected_id = case_name(job["process_seed"], job["structure"], job["condition"],
                           job["chronology_seed"], job["member"], job["age_type"])
    expected_reuse = not uncertain or job["age_type"] == "C"
    expected_source = cross.case_name(job["process_seed"], job["structure"], job["condition"],
        job["chronology_seed"], job["member"]) if expected_reuse else cross.case_name(
        job["process_seed"], job["structure"], "irregular_accurate")
    if (job["case_id"] != expected_id or job["full_result_reused"] != expected_reuse
            or job["source_case_id"] != expected_source or job["pairwise_source_case"] != expected_source
            or job["realization_id"] != PROCESS_SEEDS.index(job["process_seed"]) + 1
            or job["age_scheme_id"] != (CHRONOLOGY_SEEDS.index(job["chronology_seed"]) + 1 if uncertain else 0)
            or (uncertain and job["member"] not in range(N_AGE_MEMBERS))
            or (not uncertain and (job["chronology_seed"] != 0 or job["member"] != -1))):
        raise ValueError("Prepared case identity or source assignment changed.")
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
                edges["result_origin"] = "source03_1_complete" if job["full_result_reused"] else "source03_1_pairwise" if reused else "new_conditional"
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
            "age_type", "condition", "member", "analysis", "method"]
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
            print(f"Prepared {len(jobs)} driver-error-type inputs.", flush=True)
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
                print(f"Starting {len(pending)} L/E/J cases, three conditional methods each; workers={N_WORKERS}.", flush=True)
                context = multiprocessing.get_context("spawn")
                with ProcessPoolExecutor(max_workers=min(N_WORKERS, len(pending)), mp_context=context) as pool:
                    futures = [pool.submit(run_case, job, str(DATA_DIR), str(gpdc_null_path(job["estimator_seed"]))) for job in pending]
                    for future in as_completed(futures):
                        result = future.result()
                        completed[result.case_id.iloc[0]] = result
                        if len(completed) % 6 == 0 or len(completed) == len(jobs):
                            write_status("running", len(completed), len(jobs), started)
                            print(f"Completed {len(completed)}/{len(jobs)} driver-error cases.", flush=True)
        results = collect_results(jobs, completed)
        save_outputs(results, jobs, process_diagnostics, clock_diagnostics)
        state = "complete" if len(completed) == len(jobs) else "partial"
        write_status(state, len(completed), len(jobs), started)
        print(f"Driver-error comparison {state}: {len(completed)}/{len(jobs)} cases.", flush=True)
    except BaseException as error:
        if jobs:
            write_table(collect_results(jobs, completed), DATA_DIR / "member_results.csv")
        write_status("interrupted_or_failed", len(completed), len(jobs), started, error=f"{type(error).__name__}: {error}")
        raise



def _plot_03_2_count(states, requested, structure, method):
    """Count XY states against requested inputs without dropping unavailable cases."""
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
    count = len(requested)
    if count == 0:
        raise ValueError("No requested inputs for a figure cell.")
    return {"structure": structure, "method": method, "requested": count,
            "correct_xy": int(correct.sum()), "unresolved": int(unresolved.sum()),
            "failed": int(failed.sum()), "pending": int(pending.sum()),
            "recovery_percent": 100 * int(correct.sum()) / count}


def _plot_03_2_cell_data(states, jobs):
    """Keep the three age schemes visible beneath each type-level count."""
    cases = pd.DataFrame(jobs)
    scheme_rows, accurate_rows = [], []
    for structure in ["h0_common_driver", "h1_direct"]:
        for method in ["PCMCI+", "LKIF", "GC"]:
            for realization in [1, 2, 3]:
                references = cases[(cases.structure == structure)
                                   & (cases.realization_id == realization)
                                   & (cases.condition == "irregular_accurate")]
                if len(references) != 1:
                    raise ValueError("Expected one accurate reference per system and realization.")
                reference = _plot_03_2_count(states, references, structure, method)
                reference["realization_id"] = realization
                accurate_rows.append(reference)
                for kind in ["C", "L", "E", "J"]:
                    for scheme in [1, 2, 3]:
                        requested = cases[(cases.structure == structure)
                                          & (cases.realization_id == realization)
                                          & (cases.condition == "z_only")
                                          & (cases.age_type == kind)
                                          & (cases.age_scheme_id == scheme)]
                        cell = _plot_03_2_count(states, requested, structure, method)
                        cell.update(realization_id=realization, age_type=kind,
                                    age_scheme_id=scheme)
                        scheme_rows.append(cell)
    scheme_data = pd.DataFrame(scheme_rows)
    group_keys = ["structure", "method", "realization_id", "age_type"]
    count_fields = ["requested", "correct_xy", "unresolved", "failed", "pending"]
    type_data = scheme_data.groupby(group_keys, as_index=False, sort=False)[count_fields].sum()
    type_data["recovery_percent"] = 100 * type_data.correct_xy / type_data.requested
    return {"scheme_recovery": scheme_data, "type_recovery": type_data,
            "accurate_recovery": pd.DataFrame(accurate_rows)}


def _plot_03_2_reference_labels(accurate_data, structure, method):
    """Show each realization's own accurate-age decision alongside perturbations."""
    labels = []
    for realization in [1, 2, 3]:
        selected = accurate_data[(accurate_data.structure == structure)
                                 & (accurate_data.method == method)
                                 & (accurate_data.realization_id == realization)]
        if len(selected) != 1:
            raise ValueError("Missing or repeated accurate-reference figure entry.")
        reference = selected.iloc[0]
        label = f"Realization {realization}\nAcc: {int(reference.correct_xy)}/{int(reference.requested)}"
        exceptions = [f"{symbol}:{int(reference[field])}" for field, symbol in
                      [("unresolved", "U"), ("failed", "F"), ("pending", "P")]
                      if reference[field] > 0]
        if exceptions:
            label = label + "; " + ",".join(exceptions)
        labels.append(label)
    return labels


def _plot_03_2_clock_data(clock_diagnostics, jobs, data_dir):
    """Extract the saved age envelopes and one diagnostic row per complete scheme."""
    diagnostics = clock_diagnostics.copy()
    fields = ["proposal_rms_years", "observed_rms_years", "effective_member_sd_years",
              "accepted_members", "attempted_draws", "acceptance_rate"]
    rows, envelopes = [], []
    cases = pd.DataFrame(jobs)
    for kind in ["C", "L", "E", "J"]:
        for scheme in [1, 2, 3]:
            selected = diagnostics[(diagnostics.age_type == kind)
                                   & (diagnostics.age_scheme_id == scheme)]
            if len(selected) != 10 or selected.member.duplicated().any():
                raise ValueError("Each plotted age scheme requires its ten member diagnostics.")
            if not np.isfinite(selected.rms_error_years).all():
                raise ValueError("Age-error RMS diagnostics must be finite.")
            for field in fields:
                if selected[field].nunique(dropna=False) != 1:
                    raise ValueError(f"Scheme diagnostic is inconsistent across members: {field}.")
            row = {"age_type": kind, "age_scheme_id": scheme,
                   "median_rms_years": float(selected.rms_error_years.median()),
                   "minimum_rms_years": float(selected.rms_error_years.min()),
                   "maximum_rms_years": float(selected.rms_error_years.max())}
            row.update({field: selected[field].iloc[0] for field in fields})
            rows.append(row)
            seeds = cases[(cases.age_type == kind) & (cases.age_scheme_id == scheme)
                          & (cases.condition == "z_only")].chronology_seed.unique()
            if len(seeds) != 1:
                raise ValueError("Figure age scheme has an ambiguous source seed.")
            with np.load(Path(data_dir) / "clocks" / f"s{int(seeds[0])}_z_{kind}.npz") as clock:
                if scheme == 1:
                    envelopes.append(pd.DataFrame({"age_type": kind,
                        "age": clock["age_true"], "proposal_sd_years": 1000 * clock["age_std"]}))
    return pd.DataFrame(rows), pd.concat(envelopes, ignore_index=True)


def _plot_03_2_error_models(clock_diagnostics, clock_data, envelopes, output_dir):
    """Compare the actual error scales as well as the prescribed envelope shapes."""
    kinds = ["C", "L", "E", "J"]
    type_colors = {"C": "#0072B2", "L": "#D55E00", "E": "#009E73", "J": "#AA4499"}
    scheme_colors = ["#6A51A3", "#A6761D", "#1B9E77"]
    offsets = [-0.23, 0, 0.23]
    fig, axes = plt.subplots(1, 3, figsize=(21, 8))
    ax = axes[0]
    for kind in kinds:
        selected = envelopes[envelopes.age_type == kind]
        ax.plot(selected.age, selected.proposal_sd_years, color=type_colors[kind],
                linewidth=2, label=kind)
    ax.axvspan(30, 29, color="0.92", zorder=0)
    ax.axvspan(1, 0, color="0.92", zorder=0)
    ax.set_xlim(30, 0)
    ax.set_ylim(-1, 45)
    ax.set_xlabel("Age (kyr BP)", fontsize=12)
    ax.set_ylabel("Proposal SD (yr)", fontsize=12)
    ax.set_title("Prescribed Z-age uncertainty\nC/L/E: max 40 yr; J: five peaks of 10 yr\nDefined over 0–30 kyr; analysis uses 29–1 kyr", fontsize=12.5, pad=14)
    ax.legend(frameon=False, ncol=4, loc="upper center", fontsize=11)

    ax = axes[1]
    for kind_number, kind in enumerate(kinds):
        for scheme, color, offset in zip([1, 2, 3], scheme_colors, offsets):
            row = clock_data[(clock_data.age_type == kind)
                             & (clock_data.age_scheme_id == scheme)].iloc[0]
            x = kind_number + offset
            ax.vlines(x, row.minimum_rms_years, row.maximum_rms_years, color=color, linewidth=1.8)
            ax.scatter(x, row.median_rms_years, color=color, s=55, zorder=4)
            ax.scatter(x, row.observed_rms_years, facecolors="white", edgecolors=color,
                       marker="D", linewidths=1.4, s=48, zorder=5)
        proposal = clock_data[clock_data.age_type == kind].proposal_rms_years.unique()
        if len(proposal) != 1:
            raise ValueError("The same type has different proposal RMS across age schemes.")
        ax.plot([kind_number - 0.39, kind_number + 0.39], [proposal[0], proposal[0]],
                color="0.25", linewidth=1.2, linestyle="--", zorder=2)
    ax.set_xticks(range(4), kinds, fontsize=12)
    ax.set_xlim(-0.6, 3.6)
    ax.set_ylim(0, max(clock_data.maximum_rms_years.max(), clock_data.observed_rms_years.max()) * 1.16)
    ax.set_ylabel("RMS relative to true ages (yr)", fontsize=12)
    ax.set_title("Actual Z-age errors\nDots/ranges: 10 members; diamonds: observed centre", fontsize=13, pad=14)
    ax.text(0.5, -0.16, "Dashed: proposal-envelope RMS", transform=ax.transAxes,
            ha="center", fontsize=11)

    ax = axes[2]
    for kind_number, kind in enumerate(kinds):
        for scheme, color, offset in zip([1, 2, 3], scheme_colors, offsets):
            row = clock_data[(clock_data.age_type == kind)
                             & (clock_data.age_scheme_id == scheme)].iloc[0]
            x = kind_number + offset
            ax.bar(x, row.effective_member_sd_years, width=0.19, color=color, alpha=0.85)
            label = f"{int(row.accepted_members)}/{int(row.attempted_draws)}"
            ax.annotate(label, (x, row.effective_member_sd_years), xytext=(0, 6),
                        textcoords="offset points", ha="center", fontsize=9.5, rotation=90)
    ax.set_xticks(range(4), kinds, fontsize=12)
    ax.set_xlim(-0.6, 3.6)
    ax.set_ylim(0, clock_data.effective_member_sd_years.max() * 1.36)
    ax.set_ylabel("Mean pointwise member SD (yr)", fontsize=12)
    ax.set_title("Effective ensemble spread\nLabels: accepted / attempted whole curves", fontsize=13, pad=14)
    ax.text(0.5, -0.16, "SD across the 10 accepted members (ddof = 0)", transform=ax.transAxes,
            ha="center", fontsize=11)
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", color="0.92", linewidth=0.7)
        ax.set_axisbelow(True)
    handles = [Rectangle((0, 0), 1, 1, facecolor=color, label=f"Age scheme {scheme}")
               for scheme, color in zip([1, 2, 3], scheme_colors)]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.66, 0.91),
               ncol=3, frameon=False, fontsize=11)
    fig.suptitle("0. Z-age models differ in both error pattern and overall error size", fontsize=18, y=0.985)
    fig.text(0.5, 0.105, "X and Y keep accurate ages. Each type has three complete age schemes; "
             "each scheme is reused across three simulated histories and both systems.", ha="center", fontsize=11)
    fig.text(0.5, 0.068, "Types are not RMS-matched. J's 50% coverage refers to five defined intervals, "
             "not the fraction of samples exceeding an error threshold.", ha="center", fontsize=11)
    fig.text(0.5, 0.031, "Member error ranges are descriptive min–max ranges, not confidence intervals. "
             "Members are drawn around uncertain observed ages; all diagnostics use 29–1 kyr BP.", ha="center", fontsize=10.5)
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.30, top=0.78, wspace=0.31)
    add_panel_labels(axes, x=-0.10, y=1.04)
    save_figure(fig, output_dir, "00_error_models")
    plt.close(fig)


def _plot_03_2_recovery(type_data, accurate_data, output_dir):
    """Show all three realizations separately, with thirty requested members per type."""
    kinds = ["C", "L", "E", "J"]
    systems = [("h0_common_driver", "Common driver only", "Correct: neither X–Y direction"),
               ("h1_direct", "Common driver + direct X → Y", "Correct: X → Y only")]
    fig, axes = plt.subplots(2, 3, figsize=(21, 12))
    for row, (structure, system_name, correct_label) in enumerate(systems):
        for column, method in enumerate(["PCMCI+", "LKIF", "GC"]):
            ax = axes[row, column]
            selected = type_data[(type_data.structure == structure) & (type_data.method == method)]
            values = selected.pivot(index="realization_id", columns="age_type",
                                    values="recovery_percent").reindex(index=[1, 2, 3], columns=kinds)
            if values.isna().any().any():
                raise ValueError("Incomplete requested realization-by-type figure summary.")
            image = ax.imshow(values.to_numpy(), cmap="Blues", vmin=0, vmax=100, aspect="auto")
            for realization in [1, 2, 3]:
                for kind_number, kind in enumerate(kinds):
                    cell = selected[(selected.realization_id == realization)
                                    & (selected.age_type == kind)].iloc[0]
                    label = f"{int(cell.correct_xy)}/{int(cell.requested)}\nU: {int(cell.unresolved)}"
                    if cell.failed or cell.pending:
                        label = label + f"\nF:{int(cell.failed)} P:{int(cell.pending)}"
                    ax.text(kind_number, realization - 1, label, ha="center", va="center",
                            fontsize=11.5, linespacing=1.5,
                            color="white" if cell.recovery_percent >= 60 else "#111111")
            ax.add_patch(Rectangle((-0.47, -0.47), 0.94, 2.94, fill=False,
                                   edgecolor="0.2", linewidth=1.6))
            ax.set_xticks(range(4), kinds, fontsize=12)
            ax.set_yticks(range(3), _plot_03_2_reference_labels(accurate_data, structure, method), fontsize=10.5)
            ax.set_xticks(np.arange(-0.5, 4, 1), minor=True)
            ax.set_yticks(np.arange(-0.5, 3, 1), minor=True)
            ax.grid(which="minor", color="white", linewidth=2)
            ax.tick_params(which="both", length=0)
            ax.set_title(f"{method}\n{system_name}\n{correct_label}", fontsize=12, pad=14)
    fig.suptitle("1. Does the Z-age error model change the complete X–Y decision?", fontsize=18, y=0.99)
    fig.text(0.5, 0.945, "Conditional X,Y,Z analyses | Cells: correct / requested, with unresolved count U | "
             "Acc: own-realization irregular accurate-age reference", ha="center", fontsize=11)
    fig.subplots_adjust(left=0.105, right=0.90, bottom=0.15, top=0.855, wspace=0.55, hspace=0.53)
    colorbar_axis = fig.add_axes([0.935, 0.29, 0.015, 0.42])
    colorbar = fig.colorbar(image, cax=colorbar_axis)
    colorbar.set_ticks(np.arange(0, 101, 20))
    colorbar.set_ticklabels([f"{number}%" for number in range(0, 101, 20)])
    colorbar.set_label("Exact X–Y recovery / requested members", fontsize=11, labelpad=12)
    add_panel_labels(axes, x=-0.16, y=1.04)
    fig.text(0.5, 0.096, "Each cell contains 3 age schemes × 10 members for one simulated history. "
             "Outlined C columns reuse 03_1 exactly; Figure 02 separates the three schemes.", ha="center", fontsize=11)
    fig.text(0.5, 0.066, "Unresolved decisions are not correct. Failed (F) and pending (P) cases, where present, "
             "remain in the requested denominator.", ha="center", fontsize=10.5)
    fig.text(0.5, 0.036, "The same three simulated histories are used for every type. Types are not RMS-matched; "
             "these fractions are descriptive, not independent-history error rates.", ha="center", fontsize=10.5)
    save_figure(fig, output_dir, "01_xy_recovery")
    plt.close(fig)


def _plot_03_2_scheme_detail(scheme_data, accurate_data, output_dir):
    """Expose all three scheme results in every type-by-realization cell."""
    kinds = ["C", "L", "E", "J"]
    scheme_colors = ["#6A51A3", "#A6761D", "#1B9E77"]
    systems = [("h0_common_driver", "Common driver only", "Correct: neither X–Y direction"),
               ("h1_direct", "Common driver + direct X → Y", "Correct: X → Y only")]
    fig, axes = plt.subplots(2, 3, figsize=(24, 13.5))
    for row, (structure, system_name, correct_label) in enumerate(systems):
        for column, method in enumerate(["PCMCI+", "LKIF", "GC"]):
            ax = axes[row, column]
            selected = scheme_data[(scheme_data.structure == structure) & (scheme_data.method == method)]
            for realization in [1, 2, 3]:
                bottom = realization - 1
                if realization % 2 == 1:
                    ax.add_patch(Rectangle((0, bottom), 4, 1, facecolor="#F7F7F7", edgecolor="none", zorder=0))
                for kind_number, kind in enumerate(kinds):
                    for scheme, color in zip([1, 2, 3], scheme_colors):
                        cell = selected[(selected.realization_id == realization)
                                        & (selected.age_type == kind)
                                        & (selected.age_scheme_id == scheme)].iloc[0]
                        text_y = bottom + 0.14 + (scheme - 1) * 0.30
                        bar_y = text_y + 0.09
                        left, width = kind_number + 0.08, 0.84
                        ax.barh(bar_y, width, left=left, height=0.055, color="0.86", edgecolor="none", zorder=1)
                        ax.barh(bar_y, width * cell.correct_xy / cell.requested, left=left,
                                height=0.055, color=color, edgecolor="none", zorder=2)
                        ax.text(kind_number + 0.08, text_y, f"A{scheme}", color=color, fontsize=9.5,
                                ha="left", va="center", fontweight="bold")
                        ax.text(kind_number + 0.47, text_y, f"{int(cell.correct_xy)}/{int(cell.requested)}",
                                fontsize=9.5, ha="center", va="center")
                        ax.text(kind_number + 0.92, text_y, f"U:{int(cell.unresolved)}", fontsize=8.8,
                                ha="right", va="center")
                        if cell.failed or cell.pending:
                            ax.text(kind_number + 0.5, bar_y + 0.075,
                                    f"F:{int(cell.failed)}  P:{int(cell.pending)}", fontsize=7.5,
                                    ha="center", va="center", color="0.35")
            for edge in range(5):
                ax.axvline(edge, color="0.80", linewidth=0.8, zorder=0)
            for edge in range(4):
                ax.axhline(edge, color="0.65", linewidth=0.9, zorder=0)
            ax.add_patch(Rectangle((0.025, 0.025), 0.95, 2.95, fill=False, edgecolor="0.25", linewidth=1.4))
            ax.set_xlim(0, 4)
            ax.set_ylim(3, 0)
            ax.set_xticks(np.arange(4) + 0.5, kinds, fontsize=13)
            ax.set_yticks(np.arange(3) + 0.5,
                          _plot_03_2_reference_labels(accurate_data, structure, method), fontsize=10.5)
            ax.tick_params(which="both", length=0, pad=7)
            ax.set_title(f"{method}\n{system_name}\n{correct_label}", fontsize=12.5, pad=15)
            for spine in ax.spines.values():
                spine.set_visible(False)
    handles = [Rectangle((0, 0), 1, 1, facecolor=color, label=f"A{scheme}: age scheme {scheme}")
               for scheme, color in zip([1, 2, 3], scheme_colors)]
    fig.suptitle("2. Do the three age schemes agree within each error type?", fontsize=19, y=0.99)
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=3,
               frameon=False, fontsize=11)
    fig.text(0.5, 0.91, "Each short bar spans 0–10 correct members; a full bar means 10/10. "
             "Labels give correct/requested and unresolved count U.", ha="center", fontsize=11)
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.13, top=0.81, wspace=0.37, hspace=0.53)
    add_panel_labels(axes, x=-0.10, y=1.04)
    fig.text(0.5, 0.084, "Read across C/L/E/J within the same realization. Within each cell, "
             "A1/A2/A3 retain their fixed top-to-bottom order and colours.", ha="center", fontsize=11)
    fig.text(0.5, 0.055, "Age-scheme labels pair random-generation settings across types; whole-curve rejection "
             "can make accepted member indices differ. C results are exact copies from 03_1.", ha="center", fontsize=10.5)
    fig.text(0.5, 0.026, "Failed (F) and pending (P) cases remain in denominators. "
             "These are repeated chronology assignments to three histories, not 90 independent histories per type.", ha="center", fontsize=10.5)
    save_figure(fig, output_dir, "02_scheme_detail")
    plt.close(fig)


def plot_results(results, jobs, process_diagnostics, clock_diagnostics, data_dir, output_dir):
    """Draw three explanatory figures without editing saved scientific evidence."""
    states = summarize_results(results)["xy_case_states"]
    figure_data = _plot_03_2_cell_data(states, jobs)
    clock_data, envelopes = _plot_03_2_clock_data(clock_diagnostics, jobs, data_dir)
    figure_data["clock_summary"] = clock_data
    figure_data["proposal_envelopes"] = envelopes
    _plot_03_2_error_models(clock_diagnostics, clock_data, envelopes, output_dir)
    _plot_03_2_recovery(figure_data["type_recovery"], figure_data["accurate_recovery"], output_dir)
    _plot_03_2_scheme_detail(figure_data["scheme_recovery"], figure_data["accurate_recovery"], output_dir)
    return figure_data


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
