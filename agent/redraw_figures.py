"""Redraw the 30 existing research figures from saved results, without fitting.

Run from the project root with the production Python environment. The 01 and
02_2 notebooks have no saved results here; this launcher does not execute them.
"""

from pathlib import Path
import os
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ["redraw_main_figures.py", "redraw_sensitivity_figures.py",
           "redraw_notebook_figures.py"]


def main():
    environment = os.environ.copy()
    for name in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"]:
        environment[name] = "1"
    environment["MPLBACKEND"] = "Agg"
    environment["PYTHONUTF8"] = "1"
    environment["PYTHONIOENCODING"] = "utf-8"
    for name in SCRIPTS:
        print(f"Redrawing saved figures: {name}", flush=True)
        subprocess.run([sys.executable, "-u", str(ROOT / "agent" / name)],
                       cwd=ROOT, env=environment, check=True)
    print("Redrew 30 saved figures. The 01 and 02_2 plot cells are styled but were not executed.",
          flush=True)


if __name__ == "__main__":
    main()
