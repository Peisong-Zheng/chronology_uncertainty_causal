"""Regularly sampled, stationary processes with a known causal direction."""

import numpy as np
import pandas as pd


def generate_process(
    *,
    n_steps: int,
    dt: float,
    age_start: float,
    lag: float,
    ar_x: float = 0.8,
    ar_y: float = 0.6,
    coupling: float = 0.4,
    coupling_form: str = "linear",
    noise_std=(1.0, 1.0),
    burn_in: int = 1000,
    rng=None,
) -> dict:
    """Generate X -> Y, or independent X and Y when coupling is zero.

    The discrete equations are X[t] = ar_x * X[t-1] + u[t] and
    Y[t] = ar_y * Y[t-1] + coupling * g(X[t-lag_steps]) + v[t].
    Independent Gaussian innovations u and v have the two ``noise_std``
    values. X must have positive innovation SD; Y may have zero SD.

    ``g`` is X / driver_std for linear coupling, or
    tanh(2 * X / driver_std) / coupling_transform_std for nonlinear coupling.
    Both transforms have unit SD under X's theoretical stationary normal
    distribution. No sample-dependent normalization is applied.

    All time arguments use kyr. ``lag / dt`` must be a positive integer.
    ``age_start`` is the oldest retained age in kyr BP, referenced to 1950 CE.
    Rows advance forward in time: elapsed ``time`` increases, ``age`` decreases.
    Negative ages, if requested, represent time after 1950 CE. AR coefficients
    refer to one discrete step; changing dt does not rescale the equations.

    For positive AR memory times T_x and T_y, callers can set
    ar_x = exp(-dt / T_x) and ar_y = exp(-dt / T_y). To keep the stationary
    background SDs S_x and S_y fixed, use noise_std =
    (S_x * sqrt(1-ar_x**2), S_y * sqrt(1-ar_y**2)). S_y describes Y without
    the driver, not its total SD. A sustained unit g(X) has long-run response
    K when coupling = (1-ar_y) * K. These are discrete-model parameters;
    the delayed coupled system is not an exact joint SDE discretization.

    ``rng`` accepts a seed, a numpy Generator, or None. Return a dictionary
    containing the DataFrame ``series`` and the model metadata ``truth``.
    """
    if isinstance(n_steps, (bool, np.bool_)) or not isinstance(n_steps, (int, np.integer)):
        raise ValueError("n_steps must be an integer of at least 2.")
    if n_steps < 2:
        raise ValueError("n_steps must be an integer of at least 2.")
    if isinstance(burn_in, (bool, np.bool_)) or not isinstance(burn_in, (int, np.integer)):
        raise ValueError("burn_in must be an integer.")

    parameters = np.asarray([dt, age_start, lag, ar_x, ar_y, coupling], dtype=float)
    if not np.all(np.isfinite(parameters)):
        raise ValueError("Time and equation parameters must be finite.")
    if dt <= 0 or lag <= 0:
        raise ValueError("dt and lag must be positive, in kyr.")
    if abs(ar_x) >= 1 or abs(ar_y) >= 1:
        raise ValueError("Stationary processes require abs(ar_x) and abs(ar_y) below 1.")
    if coupling_form not in ("linear", "nonlinear"):
        raise ValueError("coupling_form must be 'linear' or 'nonlinear'.")

    lag_ratio = lag / dt
    lag_steps = int(round(lag_ratio))
    if lag_steps < 1 or not np.isclose(lag_ratio, lag_steps, rtol=0, atol=1e-10):
        raise ValueError("lag / dt must be a positive integer.")
    if burn_in < lag_steps:
        raise ValueError("burn_in must be at least lag / dt steps.")

    noise_std = np.asarray(noise_std, dtype=float)
    if noise_std.shape != (2,) or not np.all(np.isfinite(noise_std)):
        raise ValueError("noise_std must contain two finite standard deviations.")
    if noise_std[0] <= 0 or noise_std[1] < 0:
        raise ValueError("X noise SD must be positive and Y noise SD nonnegative.")

    rng = np.random.default_rng(rng)
    driver_std = noise_std[0] / np.sqrt(1 - ar_x**2)
    coupling_transform_std = 1.0
    if coupling_form == "nonlinear":
        # Fixed normal-distribution integral, independent of the generated data.
        nodes, weights = np.polynomial.hermite.hermgauss(256)
        transformed_nodes = np.tanh(2 * np.sqrt(2) * nodes)
        transformed_variance = np.sum(weights * transformed_nodes**2) / np.sqrt(np.pi)
        coupling_transform_std = np.sqrt(transformed_variance)

    # The extra lag history is discarded separately from the requested burn-in.
    start = lag_steps + burn_in
    total_steps = start + n_steps
    noise_x = rng.normal(0, noise_std[0], total_steps)
    noise_y = rng.normal(0, noise_std[1], total_steps)
    x = np.empty(total_steps)
    y = np.zeros(total_steps)

    # Initialize X from its stationary distribution; Y warms up from zero.
    x[0] = noise_x[0] / np.sqrt(1 - ar_x**2)
    for step in range(1, total_steps):
        x[step] = ar_x * x[step - 1] + noise_x[step]

    driver = x / driver_std
    if coupling_form == "nonlinear":
        driver = np.tanh(2 * driver) / coupling_transform_std

    for step in range(lag_steps, total_steps):
        y[step] = (
            ar_y * y[step - 1]
            + coupling * driver[step - lag_steps]
            + noise_y[step]
        )

    time = np.arange(n_steps) * dt
    series = pd.DataFrame(
        {
            "time": time,
            "age": age_start - time,
            "x": x[start:],
            "y": y[start:],
        }
    )
    if not np.all(np.isfinite(series.to_numpy())):
        raise ValueError("The chosen parameters produced nonfinite times or process values.")

    truth = {
        "edges": [] if coupling == 0 else [["x", "y"]],
        "coupling_form": coupling_form,
        "coupling": coupling,
        "lag": lag,
        "lag_steps": lag_steps,
        "ar_x": ar_x,
        "ar_y": ar_y,
        "noise_std": tuple(noise_std),
        "driver_std": driver_std,
        "y_background_std": noise_std[1] / np.sqrt(1 - ar_y**2),
        "coupling_gain": coupling / (1 - ar_y),
        # For negative AR coefficients these are decay envelopes, not monotone memory.
        "x_ar_decay_time": -dt / np.log(abs(ar_x)) if ar_x != 0 else 0.0,
        "y_intrinsic_ar_decay_time": -dt / np.log(abs(ar_y)) if ar_y != 0 else 0.0,
        "coupling_transform_std": coupling_transform_std,
        "dt": dt,
        "n_steps": n_steps,
        "burn_in": burn_in,
        "burn_in_duration": burn_in * dt,
        "duration": (n_steps - 1) * dt,
        "age_start": age_start,
        "time_unit": "kyr",
        "age_unit": "kyr BP",
        "value_unit": "arbitrary",
        "reference_year": 1950,
    }
    return {"series": series, "truth": truth}
