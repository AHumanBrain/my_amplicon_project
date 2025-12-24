import argparse
import subprocess
import sys
import json
import csv
from pathlib import Path


def parse_picard_pcr_metrics(metrics_file):
    """Parse Picard TargetedPcrMetrics output."""
    metrics = {}
    try:
        with open(metrics_file) as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if line.strip().startswith('CUSTOM_AMPLICON_SET'):
                    # This line is headers, next line is values
                    headers = lines[i].strip().split('\t')
                    values = lines[i+1].strip().split('\t')
                    metrics = dict(zip(headers, values))
                    break
    except Exception as e:
        print(f"Warning: Could not parse Picard metrics: {e}")
    return metrics


def parse_per_target_coverage(coverage_file):
    """Parse Picard per-target coverage output."""
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
        print(f"Warning: Could not parse per-target coverage: {e}")
    return targets


def generate_html_report(picard_metrics, per_target_stats, output_file, feedback_file):
    """Generate an enhanced HTML report with balancing recommendations."""
    
    total_amplicons = len(per_target_stats)
    mean_overall = float(picard_metrics.get('MEAN_TARGET_COVERAGE', 0))
    on_target_pct = float(picard_metrics.get('PCT_AMPLIFIED_BASES', 0)) * 100
    fold_80 = float(picard_metrics.get('FOLD_80_BASE_PENALTY', 0))
    
    # Calculate Uniformity (% > 0.2x mean)
    threshold_20 = 0.2 * mean_overall
    pct_20 = (len([t for t in per_target_stats if t['mean_coverage'] > threshold_20]) / total_amplicons * 100) if total_amplicons > 0 else 0

    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Amplicon Performance Report</title>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 40px; background: #0f172a; color: #f8fafc; line-height: 1.6; }}
        .header {{ border-bottom: 2px solid #334155; padding-bottom: 20px; margin-bottom: 30px; }}
        h1 {{ color: #38bdf8; margin: 0; }}
        .summary-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin-bottom: 40px; }}
        .stat-card {{ background: #1e293b; padding: 20px; border-radius: 12px; border: 1px solid #334155; }}
        .stat-value {{ font-size: 28px; font-weight: 800; color: #38bdf8; }}
        .stat-label {{ color: #94a3b8; font-size: 14px; text-transform: uppercase; letter-spacing: 0.05em; }}
        
        h2 {{ color: #f1f5f9; border-left: 4px solid #38bdf8; padding-left: 15px; margin-top: 40px; }}
        
        table {{ width: 100%; border-collapse: collapse; margin-top: 20px; background: #1e293b; border-radius: 8px; overflow: hidden; }}
        th, td {{ padding: 12px 15px; text-align: left; border-bottom: 1px solid #334155; }}
        th {{ background: #334155; color: #38bdf8; text-transform: uppercase; font-size: 13px; }}
        tr:hover {{ background: #2d3748; }}
        
        .badge {{ padding: 4px 8px; border-radius: 4px; font-size: 12px; font-weight: 600; }}
        .badge-up {{ background: #065f46; color: #34d399; }}
        .badge-down {{ background: #7f1d1d; color: #f87171; }}
        .badge-ok {{ background: #1e3a8a; color: #60a5fa; }}
        .badge-redesign {{ background: #d946ef; color: #ffffff; }}

        .chart-container {{ margin-top: 20px; }}
        .bar-row {{ display: flex; align-items: center; margin-bottom: 8px; }}
        .bar-label {{ width: 250px; font-size: 13px; color: #cbd5e1; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }}
        .bar-outer {{ flex: 1; background: #334155; height: 12px; border-radius: 6px; position: relative; }}
        .bar-inner {{ height: 100%; border-radius: 6px; transition: width 0.5s ease; }}
        .low {{ background: #ef4444; }}
        .med {{ background: #f59e0b; }}
        .high {{ background: #10b981; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 Amplicon Performance & Balancing Report</h1>
        <p>Advanced metrics powered by Picard & custom balancing logic.</p>
    </div>

    <div class="summary-grid">
        <div class="stat-card">
            <div class="stat-value">{total_amplicons}</div>
            <div class="stat-label">Total Amplicons</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{mean_overall:.1f}x</div>
            <div class="stat-label">Mean Coverage</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{on_target_pct:.1f}%</div>
            <div class="stat-label">On-Target Bases</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{fold_80:.2f}</div>
            <div class="stat-label">Fold-80 Penalty</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{pct_20:.1f}%</div>
            <div class="stat-label">Uniformity (>0.2x)</div>
        </div>
    </div>

    <h2>🛠️ Primer Balancing Recommendations</h2>
    <p>Adjust these concentrations in your next master pool to improve uniformity.</p>
    <table>
        <thead>
            <tr>
                <th>Target ID</th>
                <th>Coverage</th>
                <th>Relative to Mean</th>
                <th>Adjustment Factor</th>
                <th>Recommendation (v8.0 Ready)</th>
            </tr>
        </thead>
        <tbody>
"""
    
    # Sort by coverage ascending to show problems first
    sorted_stats = sorted(per_target_stats, key=lambda x: x['mean_coverage'])
    
    max_cov = max([s['mean_coverage'] for s in sorted_stats]) if sorted_stats else 1

    feedback_data = []

    for stat in sorted_stats:
        cov = stat['mean_coverage']
        rel_mean = cov / mean_overall if mean_overall > 0 else 0
        
        # Calculation for balancing
        # If rel_mean is 0.5, we need 2x concentration to bring it to mean
        # We cap adjustment between 0.1x and 10x for practical reasons
        adj_factor = 1.0 / rel_mean if rel_mean > 0 else 10.0
        adj_factor = max(0.1, min(10.0, adj_factor))
        
        if rel_mean < 0.2:
            rec = f"INCREASE concentration to {adj_factor:.1f}x"
            badge = "badge-up"
        elif rel_mean < 0.5:
            rec = f"INCREASE concentration to {adj_factor:.1f}x"
            badge = "badge-up"
        elif rel_mean > 2.0:
            rec = f"DECREASE concentration to {adj_factor:.2f}x"
            badge = "badge-down"
        else:
            rec = "Maintain current concentration"
            badge = "badge-ok"
            
        # Redesign Check
        if rel_mean == 0:
            rec = "🔥 FLAG FOR REDESIGN"
            badge = "badge-redesign"

        html += f"""
            <tr>
                <td>{stat['name']}</td>
                <td>{cov:.1f}x</td>
                <td>{rel_mean:.2f}x</td>
                <td>{adj_factor:.2f}x</td>
                <td><span class="badge {badge}">{rec}</span></td>
            </tr>
"""
        feedback_data.append({
            'target_id': stat['name'],
            'mean_coverage': cov,
            'rel_mean': rel_mean,
            'recommended_adjustment': adj_factor,
            'action': rec
        })

    html += """
        </tbody>
    </table>

    <h2>📊 Coverage Visualization</h2>
    <div class="chart-container">
"""

    for stat in sorted_stats:
        width = (stat['mean_coverage'] / max_cov) * 100 if max_cov > 0 else 0
        color = "low" if width < 20 else ("med" if width < 50 else "high")
        html += f"""
        <div class="bar-row">
            <div class="bar-label">{stat['name']}</div>
            <div class="bar-outer">
                <div class="bar-inner {color}" style="width: {width}%"></div>
            </div>
        </div>
"""

    html += """
    </div>
    <div style="margin-top: 50px; font-size: 12px; color: #64748b;">
        Report generated by <code>plot_coverage.py</code>. 
        Feedback data exported to <code>primer_balancing_feedback.json</code> for Design v8.0.
    </div>
</body>
</html>
"""
    
    with open(output_file, 'w') as f:
        f.write(html)
        
    # Also write JSON feedback for v8.0 design script
    with open(feedback_file, 'w') as f:
        json.dump({
            'overall_metrics': {
                'mean_coverage': mean_overall,
                'on_target_pct': on_target_pct,
                'fold_80': fold_80,
                'uniformity_0.2x': pct_20
            },
            'target_recommendations': feedback_data
        }, f, indent=4)


def main():
    parser = argparse.ArgumentParser(description='Enhanced Amplicon Performance Report')
    parser.add_argument('-m', '--metrics', required=True, help='Picard TargetedPcrMetrics file')
    parser.add_argument('-c', '--coverage', required=True, help='Picard per-target coverage file')
    parser.add_argument('-o', '--output', default='coverage_report.html', help='Output HTML report')
    parser.add_argument('-j', '--json', default='primer_balancing_feedback.json', help='Output JSON feedback file')
    
    args = parser.parse_args()
    
    print(f"Parsing Picard metrics: {args.metrics}")
    picard_metrics = parse_picard_pcr_metrics(args.metrics)
    
    print(f"Parsing per-target coverage: {args.coverage}")
    per_target_stats = parse_per_target_coverage(args.coverage)
    
    print("Generating enhanced report...")
    generate_html_report(picard_metrics, per_target_stats, args.output, args.json)
    print(f"Completed! Reports at {args.output} and {args.json}")


if __name__ == '__main__':
    main()
