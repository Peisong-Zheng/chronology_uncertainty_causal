"""Four-box temperature model adapted from Thomas Bauska's 14 October 2024 draft.

The numerical core is extracted from cell 12 of One_D_model_TBK.ipynb.
The original notebook is archived in data/archives/tbk_original_20261110T214805311811Z/.
Only the prescribed event table is extended for integrations beyond 30 kyr.
Dates, observations, causal analysis and plotting belong in calling code.
"""

import numpy as np


def simulate_tbk(*, runlength=10000, south_diffusion=0.005, noise_sd=0.005,
                 seed=20261009):
    """Return one annual four-box realization with the original update ordering.

    Temperature-like arrays have shape (1, 6, runlength), including the two
    ghost boxes. Time is the legacy coordinate 1..runlength in nominal years.
    Noise is a per-step temperature increment in model units, not a flux or
    a continuous-time diffusion coefficient. A local RandomState preserves
    the legacy np.random.randn sequence for a specified seed and run length.

    The last entries of target/global/area diagnostic arrays remain unfilled
    as in the original loop. The complete temperature and actual-transport
    arrays include their final updates. Callers can compute separately named
    complete mean-temperature diagnostics from the returned temperature.
    """
    if isinstance(runlength, (bool, np.bool_)) or not isinstance(runlength, (int, np.integer)) or runlength < 2:
        raise ValueError("runlength must be an integer number of years of at least two.")
    if not np.isfinite(south_diffusion) or not 0 < south_diffusion <= 0.1:
        raise ValueError("south_diffusion must be positive and no larger than the original maximum 0.1.")
    if not np.isfinite(noise_sd) or noise_sd < 0:
        raise ValueError("noise_sd must be finite and nonnegative.")
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or not 0 <= seed < 2**32:
        raise ValueError("seed must be an integer in the legacy RandomState range.")
    runlength = int(runlength)
    rng = np.random.RandomState(int(seed))

    # Basic setup and parameter values from the original model.
    nx = 1
    ny = 4
    lengthx = 1
    lengthy = 4
    dt = 1
    nsteps = int(runlength / dt)
    bpss_reoccurence = 1000
    bpss_length = 500

    Dy_var = [south_diffusion, south_diffusion, 0.1, 0.1, 0.1]
    maxD = 0.1
    maxw = 0.0

    wnoise_level = noise_sd
    wnoise_sigma = np.array([1, 1, 1, 1, 1, 1]) * wnoise_level
    wnoise_dt = 1

    arb_bpss_target = 0
    arb_bpss_tau = 10
    rectifier_global_temp = 0
    rectifier_global_tau = 100
    rectifier_area_temp = 0
    rectifier_area_tau = 100000
    gain = [0, 0, 0, 0, 0, 0]

    # Preserve the first twenty events; continue the same cycle for longer runs.
    event_count = max(20, int(np.ceil(runlength / (bpss_reoccurence + bpss_length))))
    timetrip_pre = np.arange(1, event_count + 1)
    timetrip = np.column_stack((
        bpss_reoccurence * timetrip_pre + bpss_length * (timetrip_pre - 1),
        (bpss_reoccurence + bpss_length) * timetrip_pre
    ))

    hx = lengthx / nx
    hy = lengthy / ny
    sigma = 0.7
    dt_diff = (1 / 4) * min(hx ** 2, hy ** 2) / maxD
    if maxw == 0:
        dt_conv = np.inf
    else:
        dt_conv = 4 * maxD / (maxw ** 2)
    dtstabilitycheck = sigma * min(dt_diff, dt_conv)
    if dt > dtstabilitycheck:
        raise ValueError("The annual model step exceeds the original stability check.")

    time = np.arange(1, nsteps + 1) * dt
    Temp = np.zeros((nx, ny + 2, nsteps))
    arb_bpss_target_time = np.zeros(nsteps)
    arb_bpss_actual_time = np.zeros(nsteps)
    Fadv = np.zeros((nx, ny + 2, nsteps))
    F = np.zeros((nx, ny + 2, nsteps))
    Fgain = np.zeros((nx, ny + 2, nsteps))
    Frectifier_global = np.zeros((nx, ny + 2, nsteps))
    Frectifier_area = np.zeros((nx, ny + 2, nsteps))
    wnoiseint = np.zeros((nx, ny + 2, nsteps))

    # The shape and draw order include ghost boxes, as in the notebook.
    wnoiserawtime = np.arange(0, time[-1] + wnoise_dt, wnoise_dt)
    wnoiseraw_shape = (nx, ny + 2, len(wnoiserawtime))
    wnoiseraw = wnoise_sigma.reshape(1, -1, 1) * rng.randn(*wnoiseraw_shape)
    for y in range(ny + 2):
        for x in range(nx):
            wnoiseint[x, y, :] = np.interp(time, wnoiserawtime, wnoiseraw[x, y, :])

    global_temp_arr = np.zeros(nsteps)
    area_temp_arr = np.zeros(nsteps)

    # Core equations: keep the original arithmetic and update order.
    for t in range(nsteps - 1):
        for y in range(1, 5):
            x = 0
            F[x, y, t + 1] = +(Dy_var[y] * dt / hy ** 2) * (Temp[x, y + 1, t] - Temp[x, y, t]) \
                            + (Dy_var[y - 1] * dt / hy ** 2) * (Temp[x, y - 1, t] - Temp[x, y, t])

        x = 0
        y = 4
        arb_bpss_target_time[t] = arb_bpss_target
        arb_bpss_actual_time[t + 1] = arb_bpss_actual_time[t] \
                                      + (arb_bpss_target_time[t] - arb_bpss_actual_time[t]) * dt * (1 / arb_bpss_tau)

        Fadv[x, y - 2, t + 1] = -arb_bpss_actual_time[t] * dt
        Fadv[x, y - 1, t + 1] = +0.5 * arb_bpss_actual_time[t] * dt
        Fadv[x, y, t + 1] = +0.5 * arb_bpss_actual_time[t] * dt

        global_temp_arr[t] = Temp[0, 1:5, t].mean()
        area_temp_arr[t] = Temp[0, 2:4, t].mean()
        Frectifier_global[:, :, t + 1] = (rectifier_global_temp - global_temp_arr[t]) * (1 / rectifier_global_tau) * dt
        Frectifier_area[x, 2:4, t + 1] = (rectifier_area_temp - area_temp_arr[t]) * (1 / rectifier_area_tau) * dt
        Temp[:, :, t + 1] = Temp[:, :, t] \
                           + F[:, :, t + 1] \
                           + Fadv[:, :, t + 1] \
                           + wnoiseint[:, :, t + 1] \
                           + Frectifier_global[:, :, t + 1] \
                           + Frectifier_area[:, :, t + 1] \
                           + Fgain[:, :, t]
        Fgain[:, :, t + 1] = np.array(gain).reshape(1, -1) * (Temp[:, :, t + 1] - Temp[:, :, t])

        # The switch is intentionally updated after the temperature step.
        if any((time[t] > timetrip[k, 0]) and (time[t] < timetrip[k, 1]) for k in range(len(timetrip))):
            arb_bpss_target = 1
        else:
            arb_bpss_target = -1

        Temp[0, 0, t + 1] = Temp[0, 1, t + 1]
        Temp[0, 5, t + 1] = Temp[0, 4, t + 1]

    return {
        "time": time,
        "temperature": Temp,
        "transport_target": arb_bpss_target_time,
        "transport_actual": arb_bpss_actual_time,
        "diffusion": F,
        "transport": Fadv,
        "noise": wnoiseint,
        "restoring_global": Frectifier_global,
        "restoring_area": Frectifier_area,
        "gain_flux": Fgain,
        "global_temperature": global_temp_arr,
        "area_temperature": area_temp_arr,
        "event_intervals": timetrip,
        "parameters": {
            "model": "TBK_four_box_cell12",
            "source_notebook": "One_D_model_TBK.ipynb",
            "source_cell_zero_based": 12,
            "native_step_years": dt,
            "runlength_years": runlength,
            "south_diffusion": float(south_diffusion),
            "diffusion_coefficients": list(Dy_var),
            "noise_sd": float(noise_sd),
            "noise_kind": "independent_gaussian_temperature_increment",
            "noise_grid_years": wnoise_dt,
            "seed": int(seed),
            "random_generator": "numpy.random.RandomState(MT19937)",
            "forcing_off_years": bpss_reoccurence,
            "forcing_on_years": bpss_length,
            "forcing_response_years": arb_bpss_tau,
            "forcing_event_count": event_count,
            "global_restoring_years": rectifier_global_tau,
            "area_restoring_years": rectifier_area_tau,
            "gain": gain,
            "legacy_final_diagnostics_unfilled": [
                "transport_target", "global_temperature", "area_temperature"
            ],
        },
    }
