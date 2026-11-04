"""Per-target coverage -> pool rebalancing recommendations.

Outputs a JSON consumed by `design.py --feedback` and a TSV table for humans.
rel_mean = target mean coverage / panel mean coverage; adjustment = 1/rel_mean clamped to 0.1-10.
"""
import argparse
import csv
import json

from analysis_utils import (calculate_relative_coverage, parse_per_target_coverage,
                            parse_picard_pcr_metrics, to_float)


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


def build_feedback(metrics, targets):
    mean_cov = to_float(metrics.get('MEAN_TARGET_COVERAGE'))
    targets = calculate_relative_coverage(targets, mean_cov)
    recs = []
    for t in sorted(targets, key=lambda t: t['mean_coverage']):
        adj, action = recommend(t['rel_mean'])
        recs.append({'target_id': t['name'], 'mean_coverage': t['mean_coverage'],
                     'rel_mean': t['rel_mean'], 'recommended_adjustment': adj, 'action': action})
    uniform = sum(t['rel_mean'] > 0.2 for t in targets) / len(targets) * 100 if targets else 0.0
    return {
        'overall_metrics': {
            'mean_coverage': mean_cov,
            'on_target_pct': to_float(metrics.get('PCT_AMPLIFIED_BASES')) * 100,
            'fold_80': to_float(metrics.get('FOLD_80_BASE_PENALTY')),
            'uniformity_0.2x': uniform,
        },
        'target_recommendations': recs,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('-m', '--metrics', required=True, help="Picard TargetedPcrMetrics file")
    ap.add_argument('-c', '--coverage', required=True, help="Picard per-target coverage file")
    ap.add_argument('-j', '--json', required=True, help="Output feedback JSON")
    ap.add_argument('-t', '--tsv', required=True, help="Output summary TSV")
    args = ap.parse_args()

    fb = build_feedback(parse_picard_pcr_metrics(args.metrics), parse_per_target_coverage(args.coverage))
    with open(args.json, 'w') as f:
        json.dump(fb, f, indent=2)
    with open(args.tsv, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['target', 'mean_coverage', 'rel_mean', 'adjustment', 'action'])
        for r in fb['target_recommendations']:
            w.writerow([r['target_id'], f"{r['mean_coverage']:.1f}", f"{r['rel_mean']:.3f}",
                        f"{r['recommended_adjustment']:.2f}", r['action']])
    m = fb['overall_metrics']
    print(f"mean coverage {m['mean_coverage']:.1f}x | fold-80 {m['fold_80']:.2f} | "
          f"{m['uniformity_0.2x']:.1f}% of targets > 0.2x mean | {len(fb['target_recommendations'])} targets")


if __name__ == '__main__':
    main()
