"""Keep chronology support separate from accuracy across independent records."""

import numpy as np
import pandas as pd
from statsmodels.stats.proportion import proportion_confint


def summarize_ensemble(
    member_results,
    *,
    n_members,
    baseline=None,
    support_threshold=0.90,
    min_valid_fraction=0.90,
):
    """Summarize two directed edges without reading the simulation truth.

    Input has two rows per requested member and method, including failed rows.
    Required columns are member, method, source, target, score, edge_detected,
    and status. Optional unoriented marks unresolved contemporaneous evidence.
    Missing members remain in the requested denominator for validity checks.
    Scores use every finite score with status='ok', not significant members
    alone; a successful direction score survives a failure in its partner.
    A baseline table has method/source/target/score columns. If it also has
    status, only status='ok' supplies a reference score.

    q uses members with valid decisions in BOTH directions. Joint fractions
    require no unresolved edge; their sum plus fraction_unoriented is one.
    q >= support_threshold is a reporting rule, not a calibrated probability.
    The same threshold marks high unoriented support. A thresholded summary
    does not imply that every individual member has the reported graph.
    """
    required = {"member", "method", "source", "target", "score", "edge_detected", "status"}
    if not isinstance(member_results, pd.DataFrame) or not required.issubset(member_results.columns):
        raise ValueError(f"member_results requires columns {sorted(required)}.")
    if not isinstance(n_members, (int, np.integer)) or isinstance(n_members, (bool, np.bool_)) or n_members < 1:
        raise ValueError("n_members must be a positive integer.")
    n_members = int(n_members)
    if not 0.5 < support_threshold <= 1 or not 0 < min_valid_fraction <= 1:
        raise ValueError("Require 0.5 < support_threshold <= 1 and 0 < min_valid_fraction <= 1.")
    if member_results.duplicated(["method", "member", "source", "target"]).any():
        raise ValueError("Duplicate member/direction results cannot be counted twice.")
    if member_results[["method", "member", "source", "target", "status"]].isna().any().any():
        raise ValueError("Every result needs method, member ID, direction and explicit status.")
    directions = list(zip(member_results["source"], member_results["target"]))
    if any(direction not in [("x", "y"), ("y", "x")] for direction in directions):
        raise ValueError("Only x->y and y->x result rows are supported.")
    for column in ["edge_detected", "unoriented"]:
        if column in member_results:
            values = member_results[column].dropna()
            if any(not isinstance(value, (bool, np.bool_)) for value in values):
                raise ValueError(f"{column} must be boolean or missing, not a numeric or text vote.")
    if baseline is not None:
        baseline_required = {"method", "source", "target", "score"}
        if not isinstance(baseline, pd.DataFrame) or not baseline_required.issubset(baseline.columns):
            raise ValueError("baseline requires method, source, target and score columns.")
        if baseline.duplicated(["method", "source", "target"]).any():
            raise ValueError("baseline must contain one score per method and direction.")

    summaries = []
    for method, rows in member_results.groupby("method", sort=False):
        if rows["member"].nunique() > n_members:
            raise ValueError("More member IDs were supplied than requested.")
        forward_votes = []
        reverse_votes = []
        unresolved_votes = []
        for member, pair in rows.groupby("member", sort=False):
            forward = pair[(pair["source"] == "x") & (pair["target"] == "y")]
            reverse = pair[(pair["source"] == "y") & (pair["target"] == "x")]
            if len(forward) != 1 or len(reverse) != 1:
                continue
            forward = forward.iloc[0]
            reverse = reverse.iloc[0]
            if forward["status"] != "ok" or reverse["status"] != "ok":
                continue
            if pd.isna(forward["edge_detected"]) or pd.isna(reverse["edge_detected"]):
                continue
            forward_votes.append(bool(forward["edge_detected"]))
            reverse_votes.append(bool(reverse["edge_detected"]))
            unresolved_votes.append(bool(pair["unoriented"].fillna(False).any()) if "unoriented" in pair else False)

        forward_votes = np.asarray(forward_votes, dtype=bool)
        reverse_votes = np.asarray(reverse_votes, dtype=bool)
        unresolved_votes = np.asarray(unresolved_votes, dtype=bool)
        n_valid = len(forward_votes)
        row = {
            "method": method,
            "n_requested": n_members,
            "n_reported": rows["member"].nunique(),
            "n_valid": n_valid,
            "valid_fraction": n_valid / n_members,
            "support_threshold": support_threshold,
            "q_x_to_y": np.nan,
            "q_y_to_x": np.nan,
            "fraction_x_only": np.nan,
            "fraction_y_only": np.nan,
            "fraction_both": np.nan,
            "fraction_neither": np.nan,
            "fraction_unoriented": np.nan,
        }
        if n_valid:
            resolved = ~unresolved_votes
            row["q_x_to_y"] = float(np.mean(forward_votes))
            row["q_y_to_x"] = float(np.mean(reverse_votes))
            row["fraction_x_only"] = float(np.mean(forward_votes & ~reverse_votes & resolved))
            row["fraction_y_only"] = float(np.mean(reverse_votes & ~forward_votes & resolved))
            row["fraction_both"] = float(np.mean(forward_votes & reverse_votes & resolved))
            row["fraction_neither"] = float(np.mean(~forward_votes & ~reverse_votes & resolved))
            row["fraction_unoriented"] = float(np.mean(unresolved_votes))

        row["decision"] = "insufficient_valid_members"
        row["high_unoriented"] = bool(n_valid and row["fraction_unoriented"] >= support_threshold)
        if n_valid and row["valid_fraction"] >= min_valid_fraction:
            high_forward = row["q_x_to_y"] >= support_threshold
            high_reverse = row["q_y_to_x"] >= support_threshold
            if high_forward and high_reverse:
                row["decision"] = "both"
            elif high_forward:
                row["decision"] = "x_to_y"
            elif high_reverse:
                row["decision"] = "y_to_x"
            elif row["high_unoriented"]:
                row["decision"] = "unoriented"
            else:
                row["decision"] = "no_high_support"

        for source, target in [("x", "y"), ("y", "x")]:
            direction = f"{source}_to_{target}"
            edge_rows = rows[(rows["source"] == source) & (rows["target"] == target) & (rows["status"] == "ok")]
            scores = pd.to_numeric(edge_rows["score"], errors="coerce").to_numpy()
            scores = scores[np.isfinite(scores)]
            row[f"n_scores_{direction}"] = len(scores)
            row[f"median_{direction}"] = np.nan
            row[f"q25_{direction}"] = np.nan
            row[f"q75_{direction}"] = np.nan
            row[f"iqr_{direction}"] = np.nan
            row[f"baseline_{direction}"] = np.nan
            row[f"change_{direction}"] = np.nan
            if len(scores):
                q25, median, q75 = np.percentile(scores, [25, 50, 75])
                row[f"median_{direction}"] = float(median)
                row[f"q25_{direction}"] = float(q25)
                row[f"q75_{direction}"] = float(q75)
                row[f"iqr_{direction}"] = float(q75 - q25)
            if baseline is not None:
                reference = baseline[(baseline["method"] == method) & (baseline["source"] == source) & (baseline["target"] == target)]
                if "status" in reference:
                    reference = reference[reference["status"] == "ok"]
                if len(reference) == 1:
                    reference_score = pd.to_numeric(reference["score"], errors="coerce").iloc[0]
                    if np.isfinite(reference_score):
                        row[f"baseline_{direction}"] = float(reference_score)
                        row[f"change_{direction}"] = row[f"median_{direction}"] - row[f"baseline_{direction}"]
        summaries.append(row)
    return pd.DataFrame(summaries)


def evaluate_replicates(summaries, *, truth_edges):
    """Score independent records for ONE true process, grouped by method/scenario.

    summaries has one row per replicate/scenario/method, produced above.
    truth_edges is [], [('x', 'y')], or [('y', 'x')]. No feedback experiment is
    implemented here. Wilson intervals use independent records, never members.
    Failed ensembles are counted separately, not scored as absent causal links.
    Exactness applies to the thresholded ensemble graph: unoriented support
    must remain below that ensemble's support_threshold. Individual correct
    member fractions exclude every member with any unresolved evidence.
    For independent truth, any high-support directed or unoriented connection
    contributes to null_false_edge_rate.
    missed_true_edge_rate is the complement of true_edge_rate among valid
    independent records; unavailable fits are failures, not missed edges.
    """
    truth_edges = {tuple(edge) for edge in truth_edges}
    if truth_edges not in (set(), {("x", "y")}, {("y", "x")}):
        raise ValueError("This evaluation supports an independent or single-direction bivariate truth.")
    required = {"replicate", "scenario", "method", "decision", "q_x_to_y", "q_y_to_x",
                "support_threshold", "fraction_unoriented", "fraction_x_only",
                "fraction_y_only", "fraction_neither"}
    if not isinstance(summaries, pd.DataFrame) or not required.issubset(summaries.columns):
        raise ValueError(f"summaries requires columns {sorted(required)}.")
    if summaries[["replicate", "scenario", "method"]].isna().any().any():
        raise ValueError("Every summary needs replicate, scenario and method labels.")
    decisions = {"insufficient_valid_members", "x_to_y", "y_to_x", "both", "unoriented", "no_high_support"}
    if not summaries["decision"].isin(decisions).all():
        raise ValueError("Unknown or missing ensemble decisions cannot be treated as empty graphs.")
    if summaries.duplicated(["replicate", "scenario", "method"]).any():
        raise ValueError("Each replicate/scenario/method must have exactly one summary.")
    valid_rows = summaries[summaries["decision"] != "insufficient_valid_members"]
    probabilities = valid_rows[["q_x_to_y", "q_y_to_x", "fraction_unoriented", "fraction_x_only",
                                "fraction_y_only", "fraction_neither"]].to_numpy(dtype=float)
    if not np.all(np.isfinite(probabilities)) or np.any((probabilities < 0) | (probabilities > 1)):
        raise ValueError("Valid ensemble summaries require finite support and member fractions in [0, 1].")
    thresholds = valid_rows["support_threshold"].to_numpy(dtype=float)
    if not np.all(np.isfinite(thresholds)) or np.any((thresholds <= 0.5) | (thresholds > 1)):
        raise ValueError("Valid summaries require 0.5 < support_threshold <= 1.")

    results = []
    for (method, scenario), group in summaries.groupby(["method", "scenario"], sort=False):
        valid = group[group["decision"] != "insufficient_valid_members"]
        n_valid = len(valid)
        row = {"method": method, "scenario": scenario, "n_replicates": len(group), "n_valid": n_valid}
        row["failure_fraction"] = 1 - n_valid / len(group)
        row["mean_unoriented_fraction"] = valid["fraction_unoriented"].mean()
        forward = valid["q_x_to_y"] >= valid["support_threshold"]
        reverse = valid["q_y_to_x"] >= valid["support_threshold"]
        resolved = valid["fraction_unoriented"] < valid["support_threshold"]
        if truth_edges == {("x", "y")}:
            exact = forward & ~reverse & resolved
            true_positive = forward
            false_reverse = reverse
            row["mean_correct_member_fraction"] = valid["fraction_x_only"].mean()
        elif truth_edges == {("y", "x")}:
            exact = reverse & ~forward & resolved
            true_positive = reverse
            false_reverse = forward
            row["mean_correct_member_fraction"] = valid["fraction_y_only"].mean()
        else:
            exact = ~forward & ~reverse & resolved
            true_positive = None
            false_reverse = None
            row["mean_correct_member_fraction"] = valid["fraction_neither"].mean()
        events = {"exact_direction_rate": exact, "true_edge_rate": true_positive,
                  "missed_true_edge_rate": ~true_positive if true_positive is not None else None,
                  "reverse_error_rate": false_reverse,
                  "null_false_edge_rate": forward | reverse | ~resolved if not truth_edges else None}
        for name, event in events.items():
            row[name] = np.nan
            row[f"{name}_low"] = np.nan
            row[f"{name}_high"] = np.nan
            if event is not None and n_valid:
                successes = int(np.count_nonzero(event))
                low, high = proportion_confint(successes, n_valid, alpha=0.05, method="wilson")
                row[name] = successes / n_valid
                row[f"{name}_low"] = float(low)
                row[f"{name}_high"] = float(high)
        results.append(row)
    return pd.DataFrame(results)
