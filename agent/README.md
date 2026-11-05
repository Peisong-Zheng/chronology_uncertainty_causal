# Project support tools

This directory contains maintained, versioned support code and pinned dependencies. Temporary agent work belongs in `tmp/`. Scientific scripts and notebooks remain at the project root, and the existing `toolbox/`, `data/` and `figures/` locations are unchanged.

- `run_experiment.py`: launch the root scientific scripts sequentially with one numerical thread per worker and a separate log for each script.
- `verify_optimizations.py`: compare complete native numerical outputs before and after the existing performance optimizations.
- `benchmark_workers.py`: compare worker counts on the same prepared experiment-02 inputs.
- `explore_cuda.py`: run an optional CUDA matrix-operation probe in a separate environment; this does not move the causal estimators to the GPU.
- `requirements.txt`: pinned scientific dependencies, including the corrected LKIF commit.
- `requirements-notebook.txt`: optional notebook tools plus the scientific dependencies.
- `requirements-gpu.txt`: optional CUDA probe dependencies plus the scientific dependencies.

The notebook and GPU requirements include their sibling `requirements.txt` using `-r requirements.txt`; keep the three requirements files together.

Run these Windows commands from the project root:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r agent\requirements.txt
.\.venv\Scripts\python.exe agent\verify_optimizations.py
```

The current launcher selects `02_4_check_sampling_sensitivity.py` alone, with 32 workers in the scientific script. Its [experiment plan](02_4_sampling_sensitivity_plan.md) defines the 492 paired sampling cases. Run it from the project root:

```powershell
.\.venv\Scripts\python.exe agent\run_experiment.py
```

For a new chronology-study setup on different hardware, set `EXPERIMENT="02_check_causal_pipeline.py"`, `FOLLOWUP_EXPERIMENTS=[]` in the launcher and `PREPARE_ONLY=True` in 02. Prepare the baseline and chronology inputs, then benchmark:

```powershell
.\.venv\Scripts\python.exe agent\run_experiment.py
.\.venv\Scripts\python.exe agent\benchmark_workers.py
```

For the chronology production sequence, set `PREPARE_ONLY=False`, select the measured worker count in both 02 and 02_3, and set `FOLLOWUP_EXPERIMENTS=["02_3_check_uncertainty_escalation.py"]` in the launcher. Then run:

```powershell
.\.venv\Scripts\python.exe agent\run_experiment.py
```

Keep `EXPERIMENT="02_check_causal_pipeline.py"` for this sequence. The launcher stops on a nonzero exit. Read the project-root [AGENTS.md](../AGENTS.md) for scientific constraints, optional environments and resume rules.
