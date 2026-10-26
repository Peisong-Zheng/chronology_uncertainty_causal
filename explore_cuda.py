"""Probe CUDA and time one float64 GP-related matrix operation.

This is NOT a GPU implementation of GPDC, PCMCI+, CCM, LKIF or Granger.
Run in the optional .venv-gpu environment on the NVIDIA Windows PC.
The production causal methods continue using their existing CPU packages.
"""

from pathlib import Path
import json
import os
import time

# Set numerical threads before importing the CPU numerical libraries.
for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[name] = "1"

N_SAMPLES = 1391       # Current effective sample count in GPDC after lag trimming.
N_FEATURES = 3         # Illustrative GP conditioning dimension; not a complete PCMCI test.
REPEATS = 5            # Timed repetitions after one warm-up for each backend.
SEED = 20261009        # Fixed input for comparing CPU and GPU arithmetic on this machine.
RTOL = 1e-9           # Required relative agreement of the solved coefficient vector.
ATOL = 1e-10          # Absolute agreement tolerance near zero.
ROOT = Path(__file__).resolve().parent
RESULT_DIR = ROOT / "data" / "benchmarks" / "cuda_exploration"


def main():
    try:
        import cupy as cp
        from cupyx.scipy.linalg import solve_triangular as gpu_triangular
        count = cp.cuda.runtime.getDeviceCount()
        if count == 0:
            raise RuntimeError("No CUDA device was found.")
        device = cp.cuda.runtime.getDeviceProperties(0)
    except Exception as error:
        print(f"CUDA exploration unavailable: {type(error).__name__}: {error}")
        print("Use a separate Windows NVIDIA environment; the CPU experiment is unaffected.")
        return 2

    import numpy as np
    from scipy.linalg import solve_triangular

    rng = np.random.default_rng(SEED)
    x = rng.normal(size=(N_SAMPLES, N_FEATURES)).astype(np.float64)
    y = rng.normal(size=N_SAMPLES).astype(np.float64)

    def cpu_solve():
        norms = np.sum(x * x, axis=1)
        distance = np.maximum(norms[:, None] + norms[None, :] - 2 * (x @ x.T), 0)
        covariance = np.exp(-0.5 * distance)
        covariance.flat[::N_SAMPLES + 1] = covariance.diagonal() + 1.0
        factor = np.linalg.cholesky(covariance)
        intermediate = solve_triangular(factor, y, lower=True)
        return solve_triangular(factor.T, intermediate, lower=False)

    def gpu_solve():
        # Include both input and output transfers in the timing.
        gx = cp.asarray(x)
        gy = cp.asarray(y)
        norms = cp.sum(gx * gx, axis=1)
        distance = cp.maximum(norms[:, None] + norms[None, :] - 2 * (gx @ gx.T), 0)
        covariance = cp.exp(-0.5 * distance)
        diagonal = cp.arange(N_SAMPLES)
        covariance[diagonal, diagonal] = covariance[diagonal, diagonal] + 1.0
        factor = cp.linalg.cholesky(covariance)
        intermediate = gpu_triangular(factor, gy, lower=True)
        result = gpu_triangular(factor.T, intermediate, lower=False)
        result = cp.asnumpy(result)
        cp.cuda.get_current_stream().synchronize()
        return result

    cpu_result = cpu_solve()
    gpu_result = gpu_solve()
    matches = bool(np.allclose(cpu_result, gpu_result, rtol=RTOL, atol=ATOL))
    timings = {}
    for backend, solve in [("cpu", cpu_solve), ("cuda_with_transfers", gpu_solve)]:
        elapsed = []
        for _ in range(REPEATS):
            start = time.perf_counter()
            solve()
            elapsed.append(time.perf_counter() - start)
        timings[backend] = {"seconds": elapsed, "median_seconds": float(np.median(elapsed))}
    device_name = device["name"]
    if isinstance(device_name, bytes):
        device_name = device_name.decode("utf-8")
    report = {
        "scope": "Float64 RBF covariance construction and Cholesky solve only; not a causal-method backend.",
        "device": device_name, "device_count": count, "cupy_version": cp.__version__,
        "cuda_runtime_version": cp.cuda.runtime.runtimeGetVersion(),
        "numpy_version": np.__version__, "n_samples": N_SAMPLES, "n_features": N_FEATURES,
        "seed": SEED, "repeats": REPEATS, "dtype": "float64", "rtol": RTOL, "atol": ATOL,
        "matches_cpu": matches, "max_absolute_difference": float(np.max(np.abs(cpu_result - gpu_result))),
        "timings": timings,
        "primitive_speed_ratio": timings["cpu"]["median_seconds"] / timings["cuda_with_transfers"]["median_seconds"],
        "interpretation": "Do not extrapolate this ratio to PCMCI+ or CCM. Existing production packages remain CPU-based.",
    }
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    (RESULT_DIR / "matrix_probe.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if matches else 1


if __name__ == "__main__":
    raise SystemExit(main())
