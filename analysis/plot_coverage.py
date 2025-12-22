#!/usr/bin/env python3
"""
plot_coverage.py: Generate per-amplicon coverage plots from BAM and BED files.

Usage:
    python plot_coverage.py -b sample.bam -t targets.bed -o coverage_report.html
"""

import argparse
import subprocess
import sys
from pathlib import Path


def parse_bed(bed_file):
    """Parse BED file and return list of amplicon regions."""
    amplicons = []
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
                amplicons.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': name
                })
    return amplicons


def get_coverage_for_region(bam_file, chrom, start, end):
    """Use samtools depth to get average coverage for a region."""
    try:
        cmd = ['samtools', 'depth', '-r', f'{chrom}:{start}-{end}', bam_file]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        depths = []
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                if len(parts) >= 3:
                    depths.append(int(parts[2]))
        
        return sum(depths) / len(depths) if depths else 0
    except subprocess.CalledProcessError:
        return 0


def generate_html_report(amplicons, coverages, output_file):
    """Generate an HTML report with a bar chart of coverage per amplicon."""
    
    # Determine color thresholds
    max_cov = max(coverages) if coverages else 100
    
    html = """<!DOCTYPE html>
<html>
<head>
    <title>Amplicon Coverage Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; background: #1a1a2e; color: #eee; }
        h1 { color: #00d4ff; }
        .container { max-width: 1200px; margin: 0 auto; }
        .bar-container { margin: 10px 0; display: flex; align-items: center; }
        .bar-label { width: 200px; font-size: 12px; overflow: hidden; text-overflow: ellipsis; }
        .bar-wrapper { flex: 1; background: #2a2a4a; border-radius: 4px; height: 24px; }
        .bar { height: 100%; border-radius: 4px; transition: width 0.3s; }
        .bar-value { margin-left: 10px; font-size: 12px; min-width: 60px; }
        .low { background: linear-gradient(90deg, #ff4757, #ff6b81); }
        .medium { background: linear-gradient(90deg, #ffa502, #ffcd00); }
        .high { background: linear-gradient(90deg, #2ed573, #7bed9f); }
        .summary { background: #2a2a4a; padding: 15px; border-radius: 8px; margin-bottom: 20px; }
        .stat { display: inline-block; margin-right: 30px; }
        .stat-value { font-size: 24px; font-weight: bold; color: #00d4ff; }
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Amplicon Coverage Report</h1>
        <div class="summary">
            <div class="stat">
                <div class="stat-value">""" + str(len(amplicons)) + """</div>
                <div>Total Amplicons</div>
            </div>
            <div class="stat">
                <div class="stat-value">""" + f"{sum(coverages)/len(coverages):.1f}x" + """</div>
                <div>Mean Coverage</div>
            </div>
            <div class="stat">
                <div class="stat-value">""" + f"{min(coverages):.1f}x" + """</div>
                <div>Min Coverage</div>
            </div>
            <div class="stat">
                <div class="stat-value">""" + f"{max(coverages):.1f}x" + """</div>
                <div>Max Coverage</div>
            </div>
        </div>
        <h2>Per-Amplicon Coverage</h2>
"""
    
    for amp, cov in zip(amplicons, coverages):
        # Determine color class
        if cov < 30:
            color_class = 'low'
        elif cov < 100:
            color_class = 'medium'
        else:
            color_class = 'high'
        
        # Calculate bar width (max 100%)
        width = min(100, (cov / max_cov) * 100) if max_cov > 0 else 0
        
        html += f"""
        <div class="bar-container">
            <div class="bar-label" title="{amp['name']}">{amp['name']}</div>
            <div class="bar-wrapper">
                <div class="bar {color_class}" style="width: {width}%;"></div>
            </div>
            <div class="bar-value">{cov:.1f}x</div>
        </div>
"""
    
    html += """
    </div>
</body>
</html>
"""
    
    with open(output_file, 'w') as f:
        f.write(html)
    
    print(f"Coverage report written to: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Generate per-amplicon coverage report')
    parser.add_argument('-b', '--bam', required=True, help='Input BAM file')
    parser.add_argument('-t', '--targets', required=True, help='BED file with amplicon targets')
    parser.add_argument('-o', '--output', default='coverage_report.html', help='Output HTML file')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not Path(args.bam).exists():
        print(f"Error: BAM file not found: {args.bam}")
        sys.exit(1)
    if not Path(args.targets).exists():
        print(f"Error: BED file not found: {args.targets}")
        sys.exit(1)
    
    # Parse amplicons
    print(f"Parsing amplicons from: {args.targets}")
    amplicons = parse_bed(args.targets)
    print(f"Found {len(amplicons)} amplicons")
    
    # Calculate coverage for each amplicon
    print("Calculating coverage...")
    coverages = []
    for amp in amplicons:
        cov = get_coverage_for_region(args.bam, amp['chrom'], amp['start'], amp['end'])
        coverages.append(cov)
        print(f"  {amp['name']}: {cov:.1f}x")
    
    # Generate report
    generate_html_report(amplicons, coverages, args.output)


if __name__ == '__main__':
    main()
