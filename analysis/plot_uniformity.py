"""Log-scale scatter of per-amplicon coverage relative to the panel mean."""
import argparse
import logging

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from analysis_utils import (calculate_relative_coverage, parse_per_target_coverage,  # noqa: E402
                            parse_picard_pcr_metrics, to_float)

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


def plot_uniformity(metrics_path, coverage_path, output_path, sample_name):
    mean = to_float(parse_picard_pcr_metrics(metrics_path).get('MEAN_TARGET_COVERAGE'))
    targets = parse_per_target_coverage(coverage_path)
    if not targets or mean == 0:
        logging.error("No coverage data; plot not written.")
        return
    rel = [t['rel_mean'] for t in calculate_relative_coverage(targets, mean)]
    pct = sum(r > 0.2 for r in rel) / len(rel) * 100

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.scatter(range(len(rel)), [max(r, 1e-3) for r in rel], s=18,
               c=['#2563eb' if r > 0.2 else '#dc2626' for r in rel])
    ax.set_yscale('log')
    for y, style in ((1.0, '--'), (0.2, '-'), (5.0, '--')):
        ax.axhline(y, linestyle=style, linewidth=1, color='gray')
    ax.set_title(f"{sample_name}: {pct:.1f}% of amplicons > 0.2x mean")
    ax.set_xlabel("Amplicon index")
    ax.set_ylabel("Relative coverage (log)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    print(f"Uniformity plot saved to {output_path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-m", "--metrics", required=True)
    ap.add_argument("-c", "--coverage", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-s", "--sample", default="Sample")
    a = ap.parse_args()
    plot_uniformity(a.metrics, a.coverage, a.output, a.sample)
