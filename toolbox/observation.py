"""Irregular proxy observations with dating errors, in kyr before 1950 CE."""

import numpy as np
import pandas as pd


def _check_age_grid(age):
    age = np.asarray(age, dtype=float)
    if age.ndim != 1 or len(age) < 2 or not np.all(np.isfinite(age)):
        raise ValueError("age must be a finite one-dimensional grid with at least two points.")
    spacing = -np.diff(age)
    if np.any(spacing <= 0) or not np.allclose(spacing, spacing[0], rtol=1e-8, atol=1e-12):
        raise ValueError("age must be regularly spaced and strictly decreasing (old to young).")
    return age, spacing[0]


def sample_irregular_indices(
    age,
    *,
    young_spacing,
    old_spacing,
    jitter=0.2,
    rng=None,
):
    """Select grid indices with larger mean gaps toward older ages.

    Ages and spacings are in kyr. The target spacing grows linearly from
    young_spacing at the youngest endpoint to old_spacing at the oldest.
    Each gap is multiplied by Uniform(1-jitter, 1+jitter), then rounded to
    grid steps. This quantization makes the target mean spacing approximate.
    The oldest point is included; the youngest is not forcibly appended.
    rng accepts a NumPy Generator, an integer seed, or None.
    """
    age, dt = _check_age_grid(age)
    settings = [young_spacing, old_spacing, jitter]
    if not np.all(np.isfinite(settings)):
        raise ValueError("Sampling settings must be finite.")
    if not 0 <= jitter < 1 or not 0 < young_spacing <= old_spacing:
        raise ValueError("Require 0 <= jitter < 1 and 0 < young_spacing <= old_spacing.")
    if young_spacing * (1 - jitter) < dt * (1 - 1e-10):
        raise ValueError("The shortest possible sampling gap must be at least one grid step.")

    rng = np.random.default_rng(rng)
    indices = [0]
    index = 0
    while True:
        relative_age = (age[index] - age[-1]) / (age[0] - age[-1])
        target_gap = young_spacing + (old_spacing - young_spacing) * relative_age
        gap = target_gap * rng.uniform(1 - jitter, 1 + jitter)
        index = index + int(np.rint(gap / dt))
        if index >= len(age):
            break
        indices.append(index)
    return np.asarray(indices, dtype=int)


def make_age_uncertainty(
    age,
    *,
    kind,
    scale,
    age_range,
    exp_shape=2.0,
    jump_intervals=None,
):
    """Evaluate a C/L/E/J envelope of dating-error standard deviations.

    All age quantities are in kyr BP; scale is in kyr (one standard deviation,
    not a hard bound). age_range is (young, old), shared across records.
    C is constant, L linear, and E exponential toward older ages.
    For J, supply non-overlapping (young, old) intervals. Each contains a
    symmetric smooth bell with height scale at its center, smoothly reaching
    zero at both boundaries. Outside these intervals it is exactly zero.
    With z = (age - center) / half_width, the bell is exp(-z**2/(1-z**2))
    for abs(z) < 1. It is Gaussian-like near the center, not a normal PDF.
    Coverage is the total interval length, not a probability or a threshold
    exceedance fraction. Samples need not hit the exact peak center.
    """
    age = np.asarray(age, dtype=float)
    bounds = np.asarray(age_range, dtype=float)
    if age.ndim != 1 or age.size == 0 or not np.all(np.isfinite(age)):
        raise ValueError("age must be a nonempty finite one-dimensional array.")
    if bounds.shape != (2,) or not np.all(np.isfinite(bounds)) or bounds[0] >= bounds[1]:
        raise ValueError("age_range must be finite (young, old), with young < old.")
    if np.any(age < bounds[0]) or np.any(age > bounds[1]):
        raise ValueError("All ages must fall within the shared age_range.")
    if not np.isfinite(scale) or scale < 0:
        raise ValueError("scale must be a finite nonnegative standard deviation.")

    relative_age = (age - bounds[0]) / (bounds[1] - bounds[0])
    if kind == "C":
        return np.full(age.shape, scale, dtype=float)
    if kind == "L":
        return scale * relative_age
    if kind == "E":
        if not np.isfinite(exp_shape) or exp_shape <= 0:
            raise ValueError("exp_shape must be finite and positive.")
        # Equivalent to expm1(k*u)/expm1(k), without overflow for large k.
        growth = np.exp(exp_shape * (relative_age - 1))
        growth = growth * (-np.expm1(-exp_shape * relative_age)) / (-np.expm1(-exp_shape))
        return scale * growth
    if kind != "J":
        raise ValueError("kind must be 'C', 'L', 'E', or 'J'.")

    intervals = np.asarray(jump_intervals, dtype=float)
    if intervals.ndim != 2 or intervals.shape[1] != 2 or len(intervals) == 0:
        raise ValueError("J requires a nonempty list of (young, old) intervals.")
    if not np.all(np.isfinite(intervals)) or np.any(intervals[:, 0] >= intervals[:, 1]):
        raise ValueError("Each J interval must have finite young < old boundaries.")
    ordered = intervals[np.argsort(intervals[:, 0])]
    if ordered[0, 0] < bounds[0] or ordered[-1, 1] > bounds[1]:
        raise ValueError("J intervals must lie within age_range.")
    if np.any(ordered[1:, 0] < ordered[:-1, 1]):
        raise ValueError("J intervals must not overlap.")

    age_std = np.zeros_like(age)
    for young, old in intervals:
        center = (young + old) / 2
        half_width = (old - young) / 2
        distance = (age - center) / half_width
        inside = np.abs(distance) < 1
        squared_distance = distance[inside] ** 2
        # The bell and its slopes approach zero at the interval boundaries.
        age_std[inside] = scale * np.exp(-squared_distance / (1 - squared_distance))
    return age_std


def _age_error_factor(n_samples, correlation_length):
    """Factor the sample-order kernel without adding independent dating noise."""
    sample_position = np.linspace(0.0, 1.0, n_samples)
    distance = sample_position[:, None] - sample_position[None, :]
    correlation = np.exp(-0.5 * (distance / correlation_length) ** 2)
    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
    tolerance = 100 * np.finfo(float).eps * n_samples
    if np.min(eigenvalues) < -tolerance:
        raise ValueError("Dating-error correlation matrix is not positive semidefinite.")
    # Smooth kernels can be numerically rank deficient; remove roundoff only.
    eigenvalues = np.maximum(eigenvalues, 0.0)
    return eigenvectors, np.sqrt(eigenvalues)


def _draw_age_error(age_std, correlation_length, rng):
    """Draw a smooth Gaussian error without adding diagonal noise to the kernel."""
    eigenvectors, eigenvalue_roots = _age_error_factor(len(age_std), correlation_length)
    standard_error = eigenvectors @ (eigenvalue_roots * rng.standard_normal(len(age_std)))
    return age_std * standard_error


def draw_age_ensemble(
    observed, *, n_members, age_correlation_length, max_draws, rng=None,
):
    """Draw strictly ordered working chronologies around observed ages, in kyr BP.

    Only sample_id, age_obs and age_std are read from the observed DataFrame.
    Rows must retain increasing integer sample_id, from old to young; age_obs
    itself may contain inversions. Proposals use the same squared-exponential
    kernel on normalized sample order as simulate_observations, centered at
    age_obs with unconstrained marginal SD age_std. Reject an entire proposal
    unless every adjacent age decreases strictly. Zero-SD ages remain fixed.

    Return ages (accepted member by sample), copied sample_id and diagnostics.
    max_draws counts candidate curves, including deterministic zero-error
    candidates. A depleted budget returns the accepted subset with complete=False;
    contradictory fixed ages return an empty, explicitly infeasible result.
    No sorting, clipping, support-window conditioning or value adjustment occurs.

    This is a truncated Gaussian working uncertainty distribution, not an
    automatically calibrated age posterior. Conditioning on order can change
    its mean and marginal SD. rng accepts a Generator, integer seed, or None.
    """
    required = {"sample_id", "age_obs", "age_std"}
    if not isinstance(observed, pd.DataFrame) or not required.issubset(observed.columns):
        raise ValueError("observed must be a DataFrame with sample_id, age_obs and age_std.")
    if not observed.columns.is_unique:
        raise ValueError("observed must have unique column names.")
    sample_id = observed["sample_id"].to_numpy(copy=True)
    age_obs = observed["age_obs"].to_numpy(dtype=float, copy=True)
    age_std = observed["age_std"].to_numpy(dtype=float, copy=True)
    if len(sample_id) < 2 or not np.issubdtype(sample_id.dtype, np.integer):
        raise ValueError("Require at least two integer sample IDs in known stratigraphic order.")
    if sample_id[0] < 0 or np.any(sample_id[1:] <= sample_id[:-1]):
        raise ValueError("sample_id must be nonnegative and strictly increasing; rows are not sorted.")
    if not np.all(np.isfinite(age_obs)):
        raise ValueError("age_obs must contain finite ages.")
    if not np.all(np.isfinite(age_std)) or np.any(age_std < 0):
        raise ValueError("age_std must contain finite nonnegative standard deviations.")
    for name, count in [("n_members", n_members), ("max_draws", max_draws)]:
        if isinstance(count, (bool, np.bool_)) or not isinstance(count, (int, np.integer)) or count < 1:
            raise ValueError(f"{name} must be a positive integer.")
    n_members = int(n_members)
    max_draws = int(max_draws)
    if not np.isfinite(age_correlation_length) or age_correlation_length <= 0:
        raise ValueError("age_correlation_length must be finite and positive.")

    diagnostics = {
        "requested_members": int(n_members),
        "accepted_members": 0,
        "attempted_draws": 0,
        "order_rejected": 0,
        "acceptance_rate": np.nan,
        "complete": False,
        "status": "infeasible_fixed_ages",
        "age_order_violations": int(np.count_nonzero(np.diff(age_obs) >= 0)),
        "min_age_spacing_kyr": np.nan,
    }
    fixed_age = age_obs[age_std == 0]
    if np.any(np.diff(fixed_age) >= 0):
        return {"ages": np.empty((0, len(sample_id))), "sample_id": sample_id,
                "diagnostics": diagnostics}

    rng = np.random.default_rng(rng)
    has_uncertainty = np.any(age_std > 0)
    if has_uncertainty:
        # Reuse this factor for every candidate, including rejected curves.
        eigenvectors, eigenvalue_roots = _age_error_factor(len(sample_id), age_correlation_length)
    ages = np.empty((min(n_members, max_draws), len(sample_id)))
    accepted = 0
    attempted = 0
    while accepted < n_members and attempted < max_draws:
        candidate = age_obs.copy()
        if has_uncertainty:
            standard_error = eigenvectors @ (eigenvalue_roots * rng.standard_normal(len(sample_id)))
            candidate = candidate + age_std * standard_error
        attempted = attempted + 1
        if not np.all(np.isfinite(candidate)):
            raise ValueError("Dating-error settings produced nonfinite candidate ages.")
        if np.all(np.diff(candidate) < 0):
            ages[accepted] = candidate
            accepted = accepted + 1

    ages = ages[:accepted].copy()
    diagnostics.update({
        "accepted_members": accepted,
        "attempted_draws": attempted,
        "order_rejected": attempted - accepted,
        "acceptance_rate": accepted / attempted,
        "complete": accepted == n_members,
        "status": "complete" if accepted == n_members else "max_draws_reached",
        "min_age_spacing_kyr": float(np.min(-np.diff(ages, axis=1))) if accepted else np.nan,
    })
    return {"ages": ages, "sample_id": sample_id, "diagnostics": diagnostics}


def simulate_observations(
    age,
    value,
    sample_indices,
    age_std,
    *,
    age_correlation_length=0.05,
    measurement_noise_ratio=0.0,
    rng_age=None,
    rng_measurement=None,
):
    """Generate one irregular record, keeping observed information separate from truth.

    age/value describe the full regular process, ordered old to young.
    age_std gives one dating-error SD in kyr for EACH selected sample.
    Dating errors have a squared-exponential correlation kernel on normalized
    sample order, with age_correlation_length measured as a fraction of that
    order (not kyr). Errors are centered at zero across repeated draws, not
    forced to have zero mean within one realization.

    Measurement-noise SD is measurement_noise_ratio times the population SD
    (ddof=0) of the FULL supplied value series. Reuse that series and recreate
    the same measurement RNG to pair observations across dating scenarios.
    The two RNG arguments accept independent Generators, seeds, or None.

    Returns observed/truth DataFrames and a metadata dict. age_obs is a raw
    estimate, not a reconstructed chronology: inversions and negative BP ages
    are retained. No sorting, clipping, interpolation, or monotonic repair.
    """
    age, _ = _check_age_grid(age)
    value = np.asarray(value, dtype=float)
    indices = np.asarray(sample_indices)
    age_std = np.asarray(age_std, dtype=float)
    if value.shape != age.shape or not np.all(np.isfinite(value)):
        raise ValueError("value must be finite and have the same shape as age.")
    if indices.ndim != 1 or indices.size == 0 or not np.issubdtype(indices.dtype, np.integer):
        raise ValueError("sample_indices must be a nonempty one-dimensional integer array.")
    if indices[0] < 0 or indices[-1] >= len(age) or np.any(indices[1:] <= indices[:-1]):
        raise ValueError("sample_indices must be increasing, unique, and within the process grid.")
    if age_std.shape != indices.shape or not np.all(np.isfinite(age_std)) or np.any(age_std < 0):
        raise ValueError("age_std must contain one finite nonnegative SD per sample.")
    if not np.isfinite(age_correlation_length) or age_correlation_length <= 0:
        raise ValueError("age_correlation_length must be finite and positive.")
    if not np.isfinite(measurement_noise_ratio) or measurement_noise_ratio < 0:
        raise ValueError("measurement_noise_ratio must be finite and nonnegative.")

    age_true = age[indices]
    value_true = value[indices]
    rng_age = np.random.default_rng(rng_age)
    rng_measurement = np.random.default_rng(rng_measurement)
    age_error = np.zeros(len(indices))
    if np.any(age_std > 0):
        age_error = _draw_age_error(age_std, age_correlation_length, rng_age)
    measurement_reference_std = float(np.std(value, ddof=0))
    measurement_std = measurement_noise_ratio * measurement_reference_std
    measurement_error = measurement_std * rng_measurement.standard_normal(len(indices))
    age_obs = age_true + age_error
    sample_id = np.arange(len(indices))

    observed = pd.DataFrame({
        "sample_id": sample_id,
        "age_obs": age_obs,
        "value_obs": value_true + measurement_error,
        "age_std": age_std,
        "measurement_std": measurement_std,
    })
    truth = pd.DataFrame({
        "sample_id": sample_id,
        "source_index": indices,
        "time_true": age[0] - age_true,
        "age_true": age_true,
        "value_true": value_true,
        "age_error": age_error,
        "measurement_error": measurement_error,
    })
    metadata = {
        "age_unit": "kyr BP",
        "time_unit": "kyr",
        "value_unit": "arbitrary",
        "reference_year": 1950,
        "age_start": float(age[0]),
        "age_std_definition": "unconstrained Gaussian standard deviation",
        "age_correlation_kernel": "squared_exponential",
        "age_correlation_coordinate": "normalized_sample_order",
        "age_correlation_length": float(age_correlation_length),
        "measurement_noise_ratio": float(measurement_noise_ratio),
        "measurement_reference_std": measurement_reference_std,
        "age_order_violations": int(np.count_nonzero(np.diff(age_obs) >= 0)),
    }
    return {"observed": observed, "truth": truth, "metadata": metadata}
