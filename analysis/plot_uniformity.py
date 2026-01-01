import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import logging
from pathlib import Path

# --- Modular Imports ---
from analysis_utils import parse_picard_pcr_metrics, parse_per_target_coverage, calculate_relative_coverage

# --- Logging Configuration ---
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def plot_uniformity(pcr_metrics_path, target_cov_path, output_path, sample_name):
    """Replicate Paragon Genomics uniformity plot style."""
    
    # Use modular utilities
    metrics_dict = parse_picard_pcr_metrics(pcr_metrics_path)
    global_mean = float(metrics_dict.get('MEAN_TARGET_COVERAGE', 0))
    targets_list = parse_per_target_coverage(target_cov_path)
    
    if not targets_list or global_mean == 0:
        logging.error("Required data missing. Aborting plot.")
        return

    # Convert to DataFrame for easier plotting
    df = pd.DataFrame(targets_list)
    df = pd.DataFrame(calculate_relative_coverage(targets_list, global_mean))
    
    # Calculate Uniformity Metric: % of amplicons > 0.2x mean
    uniformity_02 = (df['rel_mean'] > 0.2).sum() / len(df) * 100
    
    # Visual Setup (Dark Theme for Premium Feel)
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Plotting
    x = range(len(df))
    y = df['rel_mean']
    
    # Color-coding
    colors = ['#38bdf8' if val > 0.2 else '#ef4444' for val in y] # Blue for good, Red for low
    
    ax.scatter(x, y, c=colors, s=15, alpha=0.7, edgecolors='none')
    
    # Log Scale
    ax.set_yscale('log')
    
    # Reference Lines
    ax.axhline(y=1.0, color='#94a3b8', linestyle='--', linewidth=1, alpha=0.5, label='Mean (100%)')
    ax.axhline(y=0.2, color='#ef4444', linestyle='-', linewidth=1.5, alpha=0.8, label='0.2x Threshold (20%)')
    ax.axhline(y=5.0, color='#f59e0b', linestyle='--', linewidth=1, alpha=0.5, label='5.0x Threshold (500%)')
    
    # Annotate the Thresholds (Paragon Style)
    ax.annotate('500%', xy=(len(df)*0.02, 5.0), xytext=(len(df)*0.02, 6.5),
                arrowprops=dict(arrowstyle='->', color='#f59e0b'), color='#f59e0b', fontsize=10)
    ax.annotate('20%', xy=(len(df)*0.02, 0.2), xytext=(len(df)*0.02, 0.1),
                arrowprops=dict(arrowstyle='->', color='#ef4444'), color='#ef4444', fontsize=10)

    # Formatting
    ax.set_title(f"{sample_name}: Coverage Uniformity {uniformity_02:.1f}%", 
                 fontsize=18, fontweight='bold', pad=20, color='#f8fafc')
    ax.set_xlabel("Amplicon Index", fontsize=12, labelpad=10, color='#94a3b8')
    ax.set_ylabel("Relative Coverage (Log10)", fontsize=12, labelpad=10, color='#94a3b8')
    
    # Clean up grid and spines
    ax.grid(True, which="both", ls="-", alpha=0.1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('#334155')
    ax.spines['left'].set_color('#334155')
    
    # Y-axis ticks
    ax.set_yticks([0.01, 0.1, 0.2, 1.0, 5.0, 10.0])
    ax.get_yaxis().set_major_formatter(plt.ScalarFormatter())
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Uniformity plot saved to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Paragon-style Uniformity Plots")
    parser.add_argument("-m", "--metrics", required=True, help="Picard PCR metrics file")
    parser.add_argument("-c", "--coverage", required=True, help="Picard per-target coverage file")
    parser.add_argument("-o", "--output", required=True, help="Output image file (PNG/PDF)")
    parser.add_argument("-s", "--sample", default="Sample", help="Sample Name")
    
    args = parser.parse_args()
    plot_uniformity(args.metrics, args.coverage, args.output, args.sample)
