"""Redraw the existing main/strength figures from saved results only."""
from pathlib import Path
import importlib.util
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from toolbox.plotting import save_figure


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Any accidental estimator/generator call during redraw is a hard error.
    def forbidden(*args, **kwargs):
        raise RuntimeError("Scientific calculation is disabled in this plot-only redraw.")
    for symbol in ["generate_process", "draw_age_ensemble", "simulate_observations", "run_pcmciplus", "run_lkif", "run_ccm", "run_granger"]:
        if hasattr(module, symbol):
            setattr(module, symbol, forbidden)
    return module


def read(path):
    return pd.read_csv(path, float_precision="round_trip")


def main():
    pipeline = load_script("02_check_causal_pipeline")
    data = pipeline.DATA_DIR
    parameters = json.loads((data / "baseline_parameters.json").read_text(encoding="utf-8"))
    pipeline.CCM_LIBRARY_SIZE = parameters["CCM_LIBRARY_SIZE"]
    uniform = read(data / "baseline_original_uniform_results.csv")
    irregular = read(data / "baseline_irregular_accurate_results.csv")
    truths = {name: {"truth": truth} for name, truth in parameters["process_truth"].items()}
    fig, _ = pipeline.plot_process_truth(truths)
    save_figure(fig, pipeline.FIGURE_DIR, "00_model_truth")
    plt.close(fig)
    sampling = {}
    for process in ["linear", "nonlinear", "independent"]:
        sampling[process] = {
            "uniform": read(data / f"r0_{process}_original_uniform_series.csv"),
            "irregular": read(data / f"r0_{process}_irregular_accurate_series.csv"),
            "observed": {record: read(data / f"r0_{process}_{record}_irregular_observed.csv") for record in ["x", "y"]},
        }
    pipeline.plot_checks(uniform, irregular, sampling)
    diagnostics = read(data / "age_diagnostics.csv")
    chronology = {}
    for kind in ["C", "L", "E", "J"]:
        for record in ["x", "y"]:
            key = f"{record}_{kind}"
            with np.load(data / f"r0_{key}_chronology.npz", allow_pickle=False) as saved:
                arrays = {name: saved[name].copy() for name in saved.files}
            chronology[key] = {
                "observation": {"truth": pd.DataFrame({"age_true": arrays["age_true"]}), "observed": pd.DataFrame({"age_obs": arrays["age_obs"]})},
                "ensemble": {"ages": arrays["ages"], "diagnostics": diagnostics.loc[diagnostics.record == key].iloc[0].to_dict()},
            }
    examples = {
        "baseline": irregular,
        "accurate": sampling["linear"]["irregular"],
        "observed": {record: read(data / f"r0_linear_L-J_{record}_observed.csv") for record in ["x", "y"]},
        "resampled": [read(data / f"r0_linear_L-J_m{member}_series.csv") for member in range(parameters["N_CAUSAL_MEMBERS"])],
    }
    pipeline.plot_age_checks(read(data / "ensemble_summary.csv"), diagnostics, chronology, examples,
        scenario_names=[item[0] for item in parameters["AGE_SCENARIOS"]], methods=parameters["METHODS"],
        processes_order=["linear", "nonlinear", "independent"], figure_dir=pipeline.FIGURE_DIR,
        n_causal_members=parameters["N_CAUSAL_MEMBERS"], support_threshold=parameters["SUPPORT_THRESHOLD"])
    pipeline.plot_ccm_filter_comparison(read(data / "ccm_filter_comparison.csv"))
    strength = load_script("02_1_check_coupling_strength")
    saved = json.loads((strength.DATA_DIR / "parameters.json").read_text(encoding="utf-8"))
    # Plotting uses the saved categorical gains, never the current estimator settings.
    strength.GAINS = saved["GAINS"]
    diagnostics = json.loads((strength.DATA_DIR / "method_diagnostics.json").read_text(encoding="utf-8"))
    library_sizes = {item["lib_size"] for item in diagnostics if item["method"] == "CCM"}
    if len(library_sizes) != 1:
        raise ValueError("Strength figure requires one saved CCM library size.")
    fig = strength.plot_strength_scores(read(strength.DATA_DIR / "strength_results.csv"),
                                       ccm_library_size=library_sizes.pop())
    save_figure(fig, strength.FIGURE_DIR, "01_strength_ranking")
    plt.close(fig)
    print("Redrew 12 figures (24 PNG/PDF exports) from saved arrays and tables; no estimates recomputed.")


if __name__ == "__main__":
    main()
