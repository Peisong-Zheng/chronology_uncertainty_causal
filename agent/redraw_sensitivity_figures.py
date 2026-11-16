"""Redraw uncertainty and sampling figures from committed tables, without fitting."""

from pathlib import Path
import importlib.util
import json
import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pandas as pd


def load_analysis(filename):
    specification = importlib.util.spec_from_file_location(Path(filename).stem, ROOT / filename)
    analysis = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(analysis)
    return analysis


def main():
    uncertainty = load_analysis("02_3_check_uncertainty_escalation.py")
    parameters = json.loads((uncertainty.DATA_DIR / "parameters.json").read_text(encoding="utf-8"))
    # These constants affect envelope display, not estimator execution.
    for name in ["DURATION", "ANALYSIS_OLD", "ANALYSIS_YOUNG", "AGE_EXP_SHAPE",
                 "J_INTERVAL_WIDTH", "J_CENTER_FRACTIONS", "J_PEAK_COUNTS", "LEVELS",
                 "N_CAUSAL_MEMBERS", "SUPPORT_THRESHOLD", "MIN_VALID_FRACTION", "CCM_LIBRARY_SIZE"]:
        if name in parameters:
            setattr(uncertainty, name, parameters[name])
    scenarios = pd.read_csv(uncertainty.DATA_DIR / "scenarios.csv", float_precision="round_trip")
    summary = pd.read_csv(uncertainty.DATA_DIR / "ensemble_summary.csv", float_precision="round_trip")
    uncertainty.plot_uncertainty_envelopes(scenarios)
    uncertainty.plot_results(summary, scenarios)

    sampling = load_analysis("02_4_check_sampling_sensitivity.py")
    saved = json.loads((sampling.DATA_DIR / "parameters.json").read_text(encoding="utf-8"))
    for name in ["ANALYSIS_OLD", "ANALYSIS_YOUNG", "IRREGULAR_GAPS_YEARS"]:
        if name in saved:
            setattr(sampling, name, saved[name])
    # Use the existing results regardless of the current production CCM default.
    # A visual revision neither invokes main() nor bypasses a resume guard.
    originals = {process: pd.read_csv(sampling.DATA_DIR / f"{process}_original.csv", float_precision="round_trip")
                 for process in ["linear", "nonlinear", "independent"]}
    clocks = pd.read_csv(sampling.DATA_DIR / "sampling_clocks.csv", float_precision="round_trip")
    points = pd.read_csv(sampling.DATA_DIR / "sampling_points.csv", float_precision="round_trip")
    results = pd.read_csv(sampling.DATA_DIR / "member_results.csv", float_precision="round_trip")
    summaries = pd.read_csv(sampling.DATA_DIR / "scenario_summary.csv", float_precision="round_trip")
    contrasts = pd.read_csv(sampling.DATA_DIR / "paired_score_contrasts.csv", float_precision="round_trip")
    sampling.plot_sampling_designs(originals, clocks, points)
    sampling.plot_sampling_results(results, summaries, contrasts, sampling.FIGURE_DIR)
    sampling.plot_sampling_individual_scores(results, sampling.FIGURE_DIR)
    print("Redrew 4 uncertainty and 7 sampling figures from saved tables.")


if __name__ == "__main__":
    main()
