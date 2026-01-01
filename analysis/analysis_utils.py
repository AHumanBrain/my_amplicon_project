import csv
import logging

def parse_picard_pcr_metrics(metrics_file):
    """Parse Picard TargetedPcrMetrics output into a dictionary."""
    metrics = {}
    try:
        with open(metrics_file) as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if line.strip().startswith('CUSTOM_AMPLICON_SET'):
                    headers = lines[i].strip().split('\t')
                    values = lines[i+1].strip().split('\t')
                    metrics = dict(zip(headers, values))
                    break
    except Exception as e:
        logging.warning(f"Could not parse Picard metrics {metrics_file}: {e}")
    return metrics

def parse_per_target_coverage(coverage_file):
    """Parse Picard per-target coverage output into a list of dictionaries."""
    targets = []
    try:
        with open(coverage_file) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                targets.append({
                    'name': row['name'],
                    'chrom': row['chrom'],
                    'start': int(row['start']),
                    'end': int(row['end']),
                    'mean_coverage': float(row['mean_coverage'])
                })
    except Exception as e:
        logging.warning(f"Could not parse per-target coverage {coverage_file}: {e}")
    return targets

def calculate_relative_coverage(targets, global_mean):
    """Normalizes target coverage against the panel mean."""
    for t in targets:
        t['rel_mean'] = t['mean_coverage'] / global_mean if global_mean > 0 else 0
    return targets
