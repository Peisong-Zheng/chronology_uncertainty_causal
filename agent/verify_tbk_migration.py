"""Check the extracted TBK model directly against the archived notebook core.

Run from the project root with .venv/Scripts/python.exe agent/verify_tbk_migration.py.
This check does not generate observations, fit causal methods or modify old results.
"""

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import textwrap
import time as timer

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from toolbox.tbk_model import simulate_tbk

SOURCE = ROOT / "data/archives/tbk_original_20261110T214805311811Z/One_D_model_TBK.ipynb"
SOURCE_SHA256 = "7f8e2a5cf055ea56fc5d7f90a8f76d88c6bcc20cdea04a3a182331e253fd51ee"
SEED = 20261009
EXTENDED_YEARS = 39000
STARTUP_YEARS = 9000
ARRAY_NAMES = {
    "time": "time",
    "temperature": "Temp",
    "transport_target": "arb_bpss_target_time",
    "transport_actual": "arb_bpss_actual_time",
    "diffusion": "F",
    "transport": "Fadv",
    "noise": "wnoiseint",
    "restoring_global": "Frectifier_global",
    "restoring_area": "Frectifier_area",
    "gain_flux": "Fgain",
    "global_temperature": "global_temp_arr",
    "area_temperature": "area_temp_arr",
    "event_intervals": "timetrip",
}


def archived_core(runlength, south_diffusion, noise_sd, extend_events=False):
    """Execute the actual archived loop body, changing only requested knobs.

    The notebook's nine noisy realizations and one noiseless realization share
    the same body. Selecting i=1 or i=10 extracts one such realization. Core
    equations and their arithmetic order are read from the archive, not copied
    into this verifier. The original global RNG state is restored afterwards.
    """
    notebook = json.loads(SOURCE.read_text(encoding="utf-8"))
    source = "".join(notebook["cells"][12]["source"])
    prefix = source[:source.index("num_simulations =")]
    prefix = prefix.replace("import matplotlib.pyplot as plt\n", "")
    loop_start = source.index("for i in range(1, num_simulations + 1):")
    body_start = source.index("\n", loop_start) + 1
    body_end = source.index("    # Save Temp, F, Fadv", body_start)
    body = textwrap.dedent(source[body_start:body_end])
    code = prefix + f"\ni = {10 if noise_sd == 0 else 1}\n" + body
    assert code.count("runlength = 10000") == 1
    assert code.count("Dy_var = [0.005, 0.005, 0.1, 0.1, 0.1]") == 1
    code = code.replace("runlength = 10000", f"runlength = {runlength}")
    code = code.replace("Dy_var = [0.005, 0.005, 0.1, 0.1, 0.1]",
                        f"Dy_var = [{south_diffusion}, {south_diffusion}, 0.1, 0.1, 0.1]")
    assert noise_sd in (0, 0.005)
    if extend_events:
        event_count = max(20, int(np.ceil(runlength / 1500)))
        assert code.count("timetrip_pre = np.arange(1, 21)") == 1
        assert code.count("for k in range(20)") == 1
        code = code.replace("timetrip_pre = np.arange(1, 21)",
                            f"timetrip_pre = np.arange(1, {event_count + 1})")
        code = code.replace("for k in range(20)", f"for k in range({event_count})")
    namespace = {}
    random_state = np.random.get_state()
    try:
        np.random.seed(SEED)
        exec(compile(code, str(SOURCE) + "::cell12_core", "exec"), namespace)
    finally:
        np.random.set_state(random_state)
    return {name: namespace[original_name] for name, original_name in ARRAY_NAMES.items()}


def compare_arrays(actual, expected):
    comparison = {}
    for name in ARRAY_NAMES:
        left = actual[name]
        right = expected[name]
        assert left.shape == right.shape, (name, left.shape, right.shape)
        assert left.dtype == right.dtype, (name, left.dtype, right.dtype)
        assert np.array_equal(left, right), name
        comparison[name] = {
            "shape": list(left.shape), "dtype": str(left.dtype),
            "exact": True, "max_absolute_difference": 0.0,
            "sha256": hashlib.sha256(left.tobytes()).hexdigest(),
        }
    return comparison


def run_checks():
    start = timer.perf_counter()
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256
    cases = []
    for south_diffusion in (0.005, 0.0005):
        for noise_sd in (0.005, 0):
            before = np.random.get_state()
            result = simulate_tbk(south_diffusion=south_diffusion, noise_sd=noise_sd, seed=SEED)
            after = np.random.get_state()
            assert before[0] == after[0]
            assert np.array_equal(before[1], after[1])
            assert before[2:] == after[2:]
            reference = archived_core(10000, south_diffusion, noise_sd)
            comparison = compare_arrays(result, reference)
            cases.append({
                "runlength": 10000, "south_diffusion": south_diffusion,
                "noise_sd": noise_sd, "arrays": comparison,
                "global_random_state_unchanged": True,
            })
            print(f"Exact archived-core agreement: D_south={south_diffusion}, noise={noise_sd}", flush=True)

    # All 39 kyr comparisons use the same full draw shape. A shorter run with
    # the same seed has different draw offsets for boxes 2..6; it is not the
    # appropriate same-input reference for testing the extended integration.
    result = simulate_tbk(runlength=EXTENDED_YEARS, seed=SEED)
    reference = archived_core(EXTENDED_YEARS, 0.005, 0.005, extend_events=True)
    extended_comparison = compare_arrays(result, reference)
    legacy = archived_core(EXTENDED_YEARS, 0.005, 0.005, extend_events=False)
    prefix_matches = {}
    for name in ARRAY_NAMES:
        if name == "event_intervals":
            assert np.array_equal(result[name][:20], legacy[name])
        else:
            assert np.array_equal(result[name][..., :30000], legacy[name][..., :30000]), name
        prefix_matches[name] = True

    target = result["transport_target"]
    time = result["time"]
    intervals = result["event_intervals"]
    expected_target = np.zeros(EXTENDED_YEARS)
    inside_event = np.any(
        (time[:-2, None] > intervals[None, :, 0])
        & (time[:-2, None] < intervals[None, :, 1]), axis=1,
    )
    expected_target[1:-1] = np.where(inside_event, 1, -1)
    assert np.array_equal(target, expected_target)
    assert intervals.shape == (26, 2)
    assert np.any(target[30000:-1] == 1)
    assert not np.any(legacy["transport_target"][30000:-1] == 1)

    temperature = result["temperature"]
    for name in ARRAY_NAMES:
        assert np.all(np.isfinite(result[name])), name
    assert temperature.shape == (1, 6, EXTENDED_YEARS)
    assert np.array_equal(time, np.arange(1, EXTENDED_YEARS + 1))
    assert np.array_equal(temperature[0, 0], temperature[0, 1])
    assert np.array_equal(temperature[0, 5], temperature[0, 4])
    assert np.all(result["gain_flux"] == 0)
    assert np.all(result["transport"][0, 1:5].sum(axis=0) == 0)
    diffusion_balance = float(np.max(np.abs(result["diffusion"][0, 1:5].sum(axis=0))))
    assert diffusion_balance < 1e-12
    for name in ("transport_target", "global_temperature", "area_temperature"):
        assert result[name][-1] == 0  # The original final diagnostics were not filled.

    retained = time >= STARTUP_YEARS
    ages = (EXTENDED_YEARS - time[retained]) / 1000
    assert len(ages) == 30001
    assert ages[0] == 30 and ages[-1] == 0
    analysis_years = np.arange(10000, 38001, 20)
    assert len(analysis_years) == 1401
    assert np.all(np.isin(analysis_years, time[retained]))
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256

    return {
        "passed": True,
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": str(SOURCE.relative_to(ROOT)), "source_sha256": SOURCE_SHA256,
        "model_sha256": hashlib.sha256((ROOT / "toolbox/tbk_model.py").read_bytes()).hexdigest(),
        "python": sys.version, "numpy": np.__version__, "platform": platform.platform(),
        "seed": SEED, "default_length_cases": cases,
        "extended_case": {
            "runlength": EXTENDED_YEARS, "arrays": extended_comparison,
            "first_30000_values_match_original_20_event_model": prefix_matches,
            "event_count": len(intervals), "first_twenty_intervals_unchanged": True,
            "target_matches_original_strict_intervals_and_update_order": True,
            "continues_events_after_30000_years": True,
            "all_arrays_finite": True, "neumann_boundaries_exact": True,
            "advective_increment_sums_to_zero": True,
            "maximum_absolute_diffusive_increment_sum": diffusion_balance,
            "maximum_absolute_temperature": float(np.max(np.abs(temperature))),
            "retained_years": 30000, "retained_annual_points": len(ages),
            "analysis_grid_points": len(analysis_years),
            "legacy_final_unfilled_diagnostics_preserved": True,
        },
        "elapsed_seconds": timer.perf_counter() - start,
        "interpretation": "Exact numerical migration for identical local inputs; no causal estimation tested.",
    }


if __name__ == "__main__":
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output = ROOT / "data/benchmarks/tbk_migration_checks" / stamp
    output.mkdir(parents=True, exist_ok=False)
    report = run_checks()
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Passed: {output / 'report.json'}", flush=True)
    print(f"Elapsed: {report['elapsed_seconds']:.3f} seconds", flush=True)
