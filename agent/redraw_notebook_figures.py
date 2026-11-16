"""Redraw the seven cached notebook figures without running any estimator.

Only the explicitly selected plotting cells from 02_5, 02_6 and 02_7 execute.
All numerical inputs come from saved CSV/NPZ artifacts. The unexecuted 01 and
02_2 tutorials are not run by this entry point.
"""
from pathlib import Path
import json
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from toolbox.plotting import add_panel_labels


def plotting_cell(notebook, index, namespace):
    """Execute one reviewed plotting cell, never the notebook's fit cells."""
    document = json.loads((ROOT / notebook).read_text(encoding="utf-8"))
    cell = document["cells"][index]
    assert cell["cell_type"] == "code"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="FigureCanvasAgg is non-interactive.*", category=UserWarning)
        exec(compile("".join(cell["source"]), f"{notebook}:cell{index}", "exec"), namespace)


def base_namespace():
    return {"Path": Path, "np": np, "pd": pd, "plt": plt,
            "display": lambda value: None, "add_panel_labels": add_panel_labels}


def main():
    plt.rcdefaults()
    processes = ["linear", "nonlinear", "independent"]

    stem = "02_5_compare_gpdc_backends"
    data = ROOT / "data/processed" / stem
    ns = base_namespace()
    ns.update(PROCESSES=processes, DT=0.020,
              DATA_DIR=ROOT / "data/processed/02_check_causal_pipeline",
              FIG=ROOT / "figures" / stem,
              summary=pd.read_csv(data / "backend_summary.csv", float_precision="round_trip"),
              comparison=pd.read_csv(data / "regression_comparison.csv", float_precision="round_trip"))
    info = json.loads((data / "run_info.json").read_text(encoding="utf-8"))
    ns["ALIGN_ATOL"] = info["atol"]
    ns["backends"] = ns["summary"].backend.drop_duplicates().tolist()
    ns["cases"], ns["results"] = {}, {}
    with np.load(data / "regression_arrays.npz", allow_pickle=False) as arrays:
        for case in ns["summary"].case.drop_duplicates():
            ns["cases"][case] = {"age": arrays[f"{case}__age"]}
            for backend in ns["backends"]:
                key = f"{case}__{backend}__residuals"
                if key in arrays:
                    ns["results"][(case, backend)] = {"residuals": arrays[key]}
    plotting_cell(stem + ".ipynb", 4, ns)
    plotting_cell(stem + ".ipynb", 16, ns)
    plt.close("all")
    print(f"Redrew three figures: {stem}", flush=True)

    stem = "02_6_check_ccm_library_length"
    data = ROOT / "data/processed" / stem
    settings = json.loads((data / "settings.json").read_text(encoding="utf-8"))
    ns = base_namespace()
    with np.load(data / "lag_curves.npz", allow_pickle=False) as arrays:
        ns.update(PROCESSES=arrays["processes"].tolist(),
                  DIRECTIONS=arrays["directions"].tolist(),
                  LIBRARY_SIZES=arrays["library_sizes"].tolist(),
                  lag_years=arrays["lag_years"], lag_rho=arrays["rho"])
    ns.update(COLORS=["#4682B4", "#E69F00", "#009E73"],
              N_LIBRARIES=settings["n_libraries"], DT_KYR=settings["dt_kyr"],
              FIGURE_DIR=ROOT / "figures" / stem,
              peaks=pd.read_csv(data / "peak_lags.csv", float_precision="round_trip"))
    plotting_cell(stem + ".ipynb", 10, ns)
    plotting_cell(stem + ".ipynb", 12, ns)
    plt.close("all")
    print(f"Redrew two figures: {stem}", flush=True)

    stem = "02_7_check_ccm_neighbor_reuse"
    data = ROOT / "data/processed" / stem
    settings = json.loads((data / "settings.json").read_text(encoding="utf-8"))
    ns = base_namespace()
    ns.update(PROCESSES=processes, LAG_STEPS=np.array(settings["lag_steps"]),
              DT_KYR=settings["dt_kyr"], N_TIMING_PASSES=settings["timing_passes"],
              DIRECTIONS=["X -> Y", "Y -> X"], FIGURE_DIR=ROOT / "figures" / stem,
              timings=pd.read_csv(data / "timings.csv", float_precision="round_trip"),
              cache_stats=pd.read_csv(data / "cache_stats.csv", float_precision="round_trip"),
              results={})
    for process in processes:
        for backend in ["native", "reuse"]:
            with np.load(data / f"{process}_{backend}_pass1.npz", allow_pickle=False) as arrays:
                ns["results"][(0, process, backend)] = {"rho": arrays["rho"]}
    plotting_cell(stem + ".ipynb", 14, ns)
    plt.close("all")
    print(f"Redrew two figures: {stem}", flush=True)


if __name__ == "__main__":
    main()
