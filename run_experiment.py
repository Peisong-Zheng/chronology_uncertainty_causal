r"""Launch experiments sequentially with single-thread libraries and saved logs.

Windows: .venv\Scripts\python.exe run_experiment.py
Edit EXPERIMENT below and the scientific settings in that experiment first.
Optional FOLLOWUP_EXPERIMENTS run in order after the preceding script succeeds.
Each script has its own log; a nonzero exit stops the remaining sequence.
The launcher stays in the foreground; keep its terminal open until it finishes.
Run only one instance of an experiment at a time because its output folder is shared.
On Windows it temporarily prevents automatic system sleep while it runs.
Display sleep remains allowed, and no persistent power settings are changed.
"""

import ctypes
from datetime import datetime, timezone
from pathlib import Path
import os
import signal
import subprocess
import sys


EXPERIMENT = "02_check_causal_pipeline.py"  # Main study; use 02_3_check_uncertainty_escalation.py for four levels.
FOLLOWUP_EXPERIMENTS = ["02_3_check_uncertainty_escalation.py"]  # Complete the other levels after 02 succeeds.
NUMERICAL_THREADS = 1                     # Per-process BLAS/OpenMP threads; workers are set in the experiment.

ROOT = Path(__file__).resolve().parent
LOG_DIR = ROOT / "data" / "logs"


def set_windows_system_awake(required):
    """Request or release a temporary system-sleep assertion on this thread."""
    if os.name != "nt":
        return False
    # ES_CONTINUOUS keeps the request active; ES_SYSTEM_REQUIRED blocks system sleep.
    # Deliberately omit ES_DISPLAY_REQUIRED, so the screen may still turn off.
    flags = 0x80000001 if required else 0x80000000
    try:
        set_state = ctypes.WinDLL("kernel32", use_last_error=True).SetThreadExecutionState
        set_state.argtypes = [ctypes.c_uint]
        set_state.restype = ctypes.c_uint
        if set_state(flags) == 0:
            print("Warning: Windows could not update the temporary system-sleep request.", file=sys.stderr, flush=True)
            return False
    except (AttributeError, OSError) as error:
        print(f"Warning: Windows system-sleep control is unavailable: {error}", file=sys.stderr, flush=True)
        return False
    return True


def stop_child(process):
    """Forward Ctrl+C to the experiment group, then stop it if it cannot exit."""
    if process.poll() is not None:
        return
    if os.name == "nt":
        try:
            process.send_signal(signal.CTRL_BREAK_EVENT)
        except OSError:
            # A launcher attached only to pipes may have no console for CTRL_BREAK.
            # Recheck before taskkill: the child can exit while signal delivery fails.
            if process.poll() is None:
                subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"], check=False)
            process.wait()
            return
    else:
        os.killpg(process.pid, signal.SIGTERM)
    try:
        process.wait(timeout=15)
    except subprocess.TimeoutExpired:
        if os.name == "nt":
            if process.poll() is None:
                subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"], check=False)
        else:
            os.killpg(process.pid, signal.SIGKILL)
        process.wait()


def run_experiment(experiment=None):
    """Launch one script; calls without an argument still use EXPERIMENT."""
    script = ROOT / (EXPERIMENT if experiment is None else experiment)
    if not script.is_file() or script.suffix != ".py" or script.parent != ROOT:
        raise ValueError("Each experiment must name a Python script in the project root.")
    if NUMERICAL_THREADS != 1:
        raise ValueError("Keep one numerical thread per worker for the current process-pool design.")
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    log_path = LOG_DIR / f"{script.stem}_{stamp}.log"

    # These settings reach the child before it imports NumPy, SciPy or the estimators.
    environment = os.environ.copy()
    for name in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"]:
        environment[name] = str(NUMERICAL_THREADS)
    environment["PYTHONUTF8"] = "1"
    environment["PYTHONIOENCODING"] = "utf-8"
    environment["PYTHONUNBUFFERED"] = "1"
    environment["MPLBACKEND"] = "Agg"
    environment["MPLCONFIGDIR"] = str(ROOT / "data" / "cache" / "matplotlib")

    group_options = {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP} if os.name == "nt" else {"start_new_session": True}
    with log_path.open("w", encoding="utf-8") as log:
        header = f"Started UTC: {datetime.now(timezone.utc).isoformat()}\nInterpreter: {sys.executable}\nExperiment: {script.name}\n"
        log.write(header)
        log.flush()
        print(header, end="")
        print(f"Log: {log_path}", flush=True)
        process = subprocess.Popen(
            [sys.executable, "-u", str(script)], cwd=ROOT, env=environment,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace", bufsize=1, **group_options,
        )
        pid_message = f"Experiment PID: {process.pid}\n"
        log.write(pid_message)
        log.flush()
        print(pid_message, end="", flush=True)
        try:
            for line in process.stdout:
                log.write(line)
                log.flush()
                print(line, end="", flush=True)
            return_code = process.wait()
        except KeyboardInterrupt:
            print("\nStopping the experiment and its workers; completed checkpoints are retained.", flush=True)
            stop_child(process)
            log.write("\nInterrupted from the launcher.\n")
            return_code = 130
        finally:
            process.stdout.close()
        log.write(f"\nExit code: {return_code}\n")
    return return_code


def main():
    experiments = [EXPERIMENT, *FOLLOWUP_EXPERIMENTS]
    keep_awake = False
    try:
        keep_awake = set_windows_system_awake(True)
        if keep_awake:
            print("Windows automatic system sleep is temporarily prevented; display sleep remains allowed.", flush=True)
        # Keep the sleep request active across script boundaries and stop on failure.
        for position, experiment in enumerate(experiments, start=1):
            if len(experiments) > 1:
                print(f"Experiment {position}/{len(experiments)}: {experiment}", flush=True)
            return_code = run_experiment(experiment)
            if return_code != 0:
                print(f"Launcher stopped after {experiment} with exit code {return_code}.", flush=True)
                return return_code
        return 0
    finally:
        if keep_awake:
            set_windows_system_awake(False)


if __name__ == "__main__":
    raise SystemExit(main())
