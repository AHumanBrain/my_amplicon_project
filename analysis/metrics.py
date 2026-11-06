"""Panel metrics from samtools output (replaces Picard TargetedPcrMetrics).

On the same BAM these reproduce Picard exactly (checked on a 21-amplicon panel): mean target coverage,
per-target mean coverage, fold-80 base penalty and aligned/on-target bases.
"""
import subprocess

import numpy as np


def samtools(*args):
    """Run samtools and return stdout; raise a readable RuntimeError on failure."""
    try:
        return subprocess.run(["samtools", *args], capture_output=True, text=True, check=True).stdout
    except FileNotFoundError as e:
        raise RuntimeError("samtools not found on PATH (activate the conda env from environment.yml)") from e
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"samtools {args[0]} failed: {e.stderr.strip()}") from e


def parse_bedcov(text):
    """`samtools bedcov` output (the BED columns + summed depth) -> per-amplicon mean coverage."""
    targets = []
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) < 4:
            continue
        start, end = int(f[1]), int(f[2])
        name = f[3] if len(f) > 4 else f"{f[0]}:{start}-{end}"
        targets.append({"name": name, "chrom": f[0], "start": start, "end": end,
                        "mean_coverage": int(f[-1]) / (end - start)})
    return targets


def parse_stats(text):
    """Selected `samtools stats` summary numbers."""
    wanted = {"raw total sequences": "reads_total", "reads mapped": "reads_mapped",
              "bases mapped (cigar)": "bases_mapped"}
    out = {}
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) >= 3 and f[0] == "SN" and f[1].rstrip(":") in wanted:
            out[wanted[f[1].rstrip(":")]] = int(f[2].split()[0])
    return out


def fold_80(depths):
    """Fold over-coverage needed to raise 80% of target bases to the mean: mean / 20th percentile.

    None when the 20th percentile is 0 (undefined; Picard prints '?').
    """
    depths = np.asarray(depths, dtype=float)
    if depths.size == 0:
        return None
    p20 = np.percentile(depths, 20)
    return float(depths.mean() / p20) if p20 > 0 else None


def panel_metrics(targets, depths, stats):
    """Whole-panel numbers from per-target coverage, per-base target depths and alignment stats."""
    depths = np.asarray(depths, dtype=float)
    mean_cov = float(depths.mean()) if depths.size else 0.0
    aligned = stats.get("bases_mapped", 0)
    return {
        "n_targets": len(targets), "target_territory_bp": int(depths.size),
        "mean_target_coverage": mean_cov, "fold_80": fold_80(depths),
        "pct_on_target_bases": 100.0 * float(depths.sum()) / aligned if aligned else 0.0,
        "reads_total": stats.get("reads_total", 0), "reads_mapped": stats.get("reads_mapped", 0),
        "pct_reads_mapped": 100.0 * stats.get("reads_mapped", 0) / stats["reads_total"] if stats.get("reads_total") else 0.0,
    }


def depth_profile(text):
    """`samtools depth -a -b bed` output -> (contigs, positions (1-based), depths) as numpy arrays."""
    rows = [line.split("\t") for line in text.splitlines() if line]
    return (np.array([r[0] for r in rows]), np.array([int(r[1]) for r in rows]), np.array([int(r[2]) for r in rows]))


def measure(bam, bed, profile=False):
    """Run samtools on `bam` for the amplicons in `bed`; returns (targets, depths, stats),
    plus the per-base (contigs, positions, depths) profile when profile=True."""
    targets = parse_bedcov(samtools("bedcov", bed, bam))
    prof = depth_profile(samtools("depth", "-a", "-b", bed, bam))
    result = (targets, prof[2], parse_stats(samtools("stats", bam)))
    return result + (prof,) if profile else result
