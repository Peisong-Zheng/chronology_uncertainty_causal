"""Fixed TBK records and independent, ordered southern/northern chronologies.

The existing observation module supplies the Gaussian covariance factor. The
explicit proposal loop adds candidate IDs and separate order/coverage counts;
it never changes a temperature, sorts an age curve or repairs a rejected draw.
"""

from functools import lru_cache

import numpy as np
import pandas as pd

from toolbox.causal import resample_pair
from toolbox.observation import _age_error_factor
from toolbox.tbk_model import simulate_tbk


@lru_cache(maxsize=2)
def age_error_factor(n_samples, correlation_length):
    """Reuse the unchanged native factor across records and amplitude levels."""
    return _age_error_factor(n_samples, correlation_length)


def draw_ordered_chronologies(center, age_std, factor, rng, n_members,
                             analysis_old, analysis_young, max_draws):
    """Draw around the supplied center, with the native Gaussian arithmetic."""
    eigenvectors, eigenvalue_roots = factor
    ages, candidate_ids = [], []
    order_rejected = coverage_rejected = 0
    attempted = 0
    while len(ages) < n_members and attempted < max_draws:
        standard_error = eigenvectors @ (eigenvalue_roots * rng.standard_normal(len(center)))
        candidate = center + age_std * standard_error
        attempted = attempted + 1
        if not np.all(np.isfinite(candidate)):
            raise ValueError("Nonfinite TBK chronology proposal.")
        if not np.all(np.diff(candidate) < 0):
            order_rejected = order_rejected + 1
        elif candidate[0] < analysis_old or candidate[-1] > analysis_young:
            coverage_rejected = coverage_rejected + 1
        else:
            ages.append(candidate)
            candidate_ids.append(attempted)  # one-based, including rejected draws
    counts = dict(requested=n_members, accepted=len(ages), attempted=attempted,
                  order_rejected=order_rejected, coverage_rejected=coverage_rejected,
                  acceptance_rate=len(ages) / attempted,
                  complete=len(ages) == n_members)
    if len(ages) != n_members:
        raise RuntimeError(f"TBK chronology draw budget exhausted: {counts}")
    return np.asarray(ages), np.asarray(candidate_ids, dtype=int), counts


def error_diagnostics(age, center, age_true, age_grid):
    """Errors on the common true-time analysis window, in years."""
    error = np.interp(-age_grid, -age_true, (age - age_true) * 1000)
    around_center = np.interp(-age_grid, -age_true, (age - center) * 1000)
    spacing = -np.diff(age) * 1000
    stretch = np.diff(age) / np.diff(age_true)
    return dict(mean_error_years=float(error.mean()),
                rms_error_years=float(np.sqrt(np.mean(error ** 2))),
                max_abs_error_years=float(np.max(np.abs(error))),
                rms_about_center_years=float(np.sqrt(np.mean(around_center ** 2))),
                min_spacing_years=float(spacing.min()),
                min_stretch=float(stretch.min()), max_stretch=float(stretch.max()),
                covers_grid=bool(age[0] >= age_grid[0] and age[-1] <= age_grid[-1]))


def prepare_tbk_records(process_seed, sigma_years=(20, 50, 100, 200),
                        n_members=10, runlength=39000, spinup_years=9000,
                        south_diffusion=0.005, noise_sd=0.005,
                        sample_step_years=20, analysis_old=29., analysis_young=1.,
                        correlation_years=1000., max_draws=10000):
    """Prepare one physical history; X is T1 (South), Y is T4 (North).

    Observed centers are drawn around truth; members are drawn around those
    observed centers. Independent SeedSequences identify the two records and
    two uncertainty layers. Recreating each stream at each level pairs candidate
    Gaussian draws, not necessarily the accepted member indices.
    """
    if not 0 < spinup_years < runlength:
        raise ValueError("Require a positive spinup shorter than the integration.")
    if not isinstance(sample_step_years, int) or sample_step_years < 1:
        raise ValueError("sample_step_years must be a positive integer.")
    if (runlength - spinup_years) % sample_step_years:
        raise ValueError("Sampling must include both retained endpoints.")
    if not 0 <= analysis_young < analysis_old <= (runlength - spinup_years) / 1000:
        raise ValueError("Analysis window must lie inside the retained support.")
    if n_members < 1 or max_draws < n_members or correlation_years <= 0:
        raise ValueError("Invalid chronology draw count or correlation length.")
    if any(not np.isfinite(s) or s <= 0 for s in sigma_years):
        raise ValueError("Every uncertain level requires a finite positive SD.")
    model = simulate_tbk(runlength=runlength, south_diffusion=south_diffusion,
                         noise_sd=noise_sd, seed=process_seed)
    sampled_time = np.arange(spinup_years, runlength + 1, sample_step_years)
    indices = sampled_time - 1  # original model time is 1..runlength
    age_true = (runlength - sampled_time) / 1000.
    temperatures = model['temperature'][0]
    sampled = pd.DataFrame(dict(sample_id=np.arange(len(indices)),
                                model_time_years=sampled_time, age_true=age_true,
                                x=temperatures[1, indices], y=temperatures[4, indices],
                                t2=temperatures[2, indices], t3=temperatures[3, indices],
                                q=model['transport_actual'][indices]))
    selected = (age_true <= analysis_old) & (age_true >= analysis_young)
    age_grid = age_true[selected].copy()
    if age_grid[0] != analysis_old or age_grid[-1] != analysis_young:
        raise ValueError("Analysis endpoints must coincide with the sampling grid.")
    x, y = sampled['x'].to_numpy(), sampled['y'].to_numpy()
    accurate = resample_pair(age_true, x, age_true, y, age_grid=age_grid)['series']
    np.testing.assert_array_equal(accurate['x'].to_numpy(), x[selected])
    np.testing.assert_array_equal(accurate['y'].to_numpy(), y[selected])
    case_base = dict(process_seed=int(process_seed), estimator_seed=int(process_seed) + 400,
                     input_status='ok')
    cases = [dict(case_base, case_id=f'tbk_s{process_seed}_accurate', level=0,
                  sigma_years=0., member=-1, condition='accurate', series=accurate)]
    correlation_length = correlation_years / (runlength - spinup_years)
    factor = age_error_factor(len(age_true), correlation_length)
    chronologies, chronology_rows, relative_rows, acceptance_rows = {}, [], [], []
    for level, sigma in enumerate(sigma_years, start=1):
        age_std = np.full(len(age_true), sigma / 1000.)
        for record_index, (record, value) in enumerate([('x', x), ('y', y)]):
            center_seed = [int(process_seed), 4101, record_index]
            member_seed = [int(process_seed), 4102, record_index]
            center, center_ids, center_counts = draw_ordered_chronologies(
                age_true, age_std, factor, np.random.default_rng(center_seed),
                1, analysis_old, analysis_young, max_draws)
            center = center[0]
            ages, member_ids, member_counts = draw_ordered_chronologies(
                center, age_std, factor, np.random.default_rng(member_seed),
                n_members, analysis_old, analysis_young, max_draws)
            base = dict(process_seed=int(process_seed), level=level,
                        sigma_years=float(sigma), record=record)
            chronologies[f'level{level}_{record}'] = dict(
                base, age_true=age_true.copy(), age_obs=center, ages=ages,
                age_std=age_std.copy(), value=value.copy(),
                member_candidate_ids=member_ids, center_candidate_id=int(center_ids[0]),
                center_seed=center_seed, member_seed=member_seed,
                center_counts=center_counts, member_counts=member_counts,
                correlation_length=correlation_length, correlation_years=correlation_years)
            for phase, counts in [('center', center_counts), ('member', member_counts)]:
                acceptance_rows.append(dict(base, phase=phase, **counts))
            chronology_rows.append(dict(base, phase='center', member=-1,
                                        candidate_id=int(center_ids[0]),
                                        **error_diagnostics(center, center, age_true, age_grid)))
            for member, ages_m in enumerate(ages):
                chronology_rows.append(dict(base, phase='member', member=member,
                                            candidate_id=int(member_ids[member]),
                                            **error_diagnostics(ages_m, center, age_true, age_grid)))
        southern = chronologies[f'level{level}_x']
        northern = chronologies[f'level{level}_y']
        for member in range(-1, n_members):
            ages_x = southern['age_obs'] if member < 0 else southern['ages'][member]
            ages_y = northern['age_obs'] if member < 0 else northern['ages'][member]
            relative = np.interp(-age_grid, -age_true, (ages_y - ages_x) * 1000)
            relative_rows.append(dict(process_seed=int(process_seed), level=level,
                                      sigma_years=float(sigma), member=member,
                                      phase='center' if member < 0 else 'member',
                                      relative_rms_years=float(np.sqrt(np.mean(relative ** 2))),
                                      relative_mean_years=float(relative.mean()),
                                      relative_max_abs_years=float(np.max(np.abs(relative)))))
            if member >= 0:
                series = resample_pair(ages_x, x, ages_y, y, age_grid=age_grid)['series']
                cases.append(dict(case_base, case_id=f'tbk_s{process_seed}_level{level}_m{member:02d}',
                                  level=level, sigma_years=float(sigma), member=member,
                                  condition='chronology', series=series))
    process_diagnostics = dict(process_seed=int(process_seed), native_points=len(model['time']),
                               retained_annual_points=runlength - spinup_years + 1,
                               sampled_points=len(sampled), analysis_points=len(age_grid),
                               sample_step_years=sample_step_years, spinup_years=spinup_years,
                               runlength=runlength, south_diffusion=south_diffusion,
                               noise_sd=noise_sd, correlation_years=correlation_years,
                               correlation_length=correlation_length,
                               all_model_values_finite=bool(np.isfinite(temperatures).all()),
                               south_sd=float(np.std(x[selected])), north_sd=float(np.std(y[selected])),
                               zero_error_exact=True)
    return dict(model=model, sampled=sampled, age_grid=age_grid, cases=cases,
                chronologies=chronologies,
                diagnostics=dict(process=process_diagnostics,
                                 chronology=pd.DataFrame(chronology_rows),
                                 relative=pd.DataFrame(relative_rows),
                                 acceptance=pd.DataFrame(acceptance_rows)))
