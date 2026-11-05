"""Per-target coverage -> pool rebalancing recommendations (the feedback loop's input).

rel_mean = target mean coverage / panel mean coverage; adjustment = 1/rel_mean clamped to 0.1-10.
"""


def calculate_relative_coverage(targets, panel_mean):
    for t in targets:
        t["rel_mean"] = t["mean_coverage"] / panel_mean if panel_mean > 0 else 0.0
    return targets


def recommend(rel_mean):
    """Return (adjustment_factor, action) for a relative coverage value."""
    if rel_mean == 0:
        return 10.0, "FLAG FOR REDESIGN"
    adj = max(0.1, min(10.0, 1.0 / rel_mean))
    if rel_mean < 0.5:
        return adj, f"INCREASE concentration to {adj:.1f}x"
    if rel_mean > 2.0:
        return adj, f"DECREASE concentration to {adj:.2f}x"
    return adj, "Maintain current concentration"


def build_feedback(panel, targets):
    """JSON consumed by `design.py --feedback`. `panel` comes from metrics.panel_metrics."""
    mean_cov = panel["mean_target_coverage"]
    targets = calculate_relative_coverage(targets, mean_cov)
    recs = []
    for t in sorted(targets, key=lambda t: t["mean_coverage"]):
        adj, action = recommend(t["rel_mean"])
        recs.append({"target_id": t["name"], "mean_coverage": t["mean_coverage"], "rel_mean": t["rel_mean"],
                     "recommended_adjustment": adj, "action": action})
    uniform = sum(t["rel_mean"] > 0.2 for t in targets) / len(targets) * 100 if targets else 0.0
    return {
        "overall_metrics": {"mean_coverage": mean_cov, "on_target_pct": panel["pct_on_target_bases"],
                            "fold_80": panel["fold_80"], "uniformity_0.2x": uniform},
        "target_recommendations": recs,
    }
