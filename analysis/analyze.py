"""Amplicon sequencing analysis: primer trimming -> alignment -> per-amplicon coverage -> rebalancing feedback.

    python analysis/analyze.py --reads "raw_data/S_R{1,2}.fastq.gz" --genome ref.fa \\
        --bed design/output/panel_pool_1.bed --primers design/output/panel_pool_1_trimming.fasta --outdir results

Tools: cutadapt, bwa, samtools (plus numpy, matplotlib). One run per pool (each pool has its own BED and
primer FASTA). `--reads` uses `{1,2}` for the read number (paired-end) or is a plain glob (single-end);
a `*` in the head matches several samples. Per sample it writes, in <outdir>/<sample>/:
  trimmed FASTQ + cutadapt.log, sorted BAM (+ .bai), target_coverage.tsv, balancing_feedback.json,
  target_summary.tsv, uniformity.png, summary.json
"""
import argparse
import csv
import glob
import json
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import balance  # noqa: E402
import metrics  # noqa: E402
import plot_uniformity  # noqa: E402

TOOLS = ("cutadapt", "bwa", "samtools")


def check_tools():
    missing = [t for t in TOOLS if not shutil.which(t)]
    if missing:
        raise RuntimeError(f"missing on PATH: {', '.join(missing)}. Activate the conda env from environment.yml.")


def find_samples(pattern):
    """[(sample, [read files])]. '{1,2}' marks the read number (paired-end); otherwise single-end glob.

    The text after '{1,2}' must be literal (e.g. '.fastq.gz'); '*' is allowed before it.
    """
    if "{1,2}" not in pattern:
        files = sorted(glob.glob(pattern))
        return [(re.sub(r"(\.fastq|\.fq)(\.gz)?$", "", os.path.basename(f)), [f]) for f in files]
    head, tail = pattern.split("{1,2}")
    if any(c in tail for c in "*?[]"):
        raise ValueError("the text after {1,2} must not contain wildcards")
    samples = []
    for r1 in sorted(glob.glob(head + "1" + tail)):
        r2 = r1[:len(r1) - len(tail) - 1] + "2" + tail
        if not os.path.exists(r2):
            raise FileNotFoundError(f"mate not found for {r1}: expected {r2}")
        stem = os.path.basename(r1)[:len(os.path.basename(r1)) - len(tail) - 1]
        samples.append((re.sub(r"[_.-]?[Rr]?$", "", stem) or stem, [r1, r2]))
    return samples


def split_primers(fasta, outdir):
    """Write fwd/rev primer FASTAs from a trimming FASTA whose names end in _F / _R."""
    paths = {}
    for tag in ("F", "R"):
        paths[tag] = os.path.join(outdir, f"primers_{tag}.fa")
        keep, lines = False, []
        for line in open(fasta):
            if line.startswith(">"):
                keep = line.strip().endswith(f"_{tag}")
            if keep:
                lines.append(line if line.endswith("\n") else line + "\n")
        with open(paths[tag], "w") as f:
            f.writelines(lines)
    if os.path.getsize(paths["F"]) == 0:
        raise ValueError(f"{fasta}: no forward primers (names must end in _F / _R)")
    return paths["F"], paths["R"]


def run(cmd, log=None):
    try:
        return subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"{cmd[0]} failed: {e.stderr.strip()[-500:]}") from e


def trim(sample, reads, fwd, rev, outdir, threads):
    """Anchored 5' primer trimming (R1: forward primers, R2: reverse primers). Returns (files, stats)."""
    base = ["cutadapt", "--discard-untrimmed", "-j", str(threads), "-g", f"^file:{fwd}"]
    if len(reads) == 2:
        out = [os.path.join(outdir, f"{sample}.trimmed_R{i}.fastq.gz") for i in (1, 2)]
        cmd = base + ["-G", f"^file:{rev}", "-o", out[0], "-p", out[1], *reads]
    else:
        out = [os.path.join(outdir, f"{sample}.trimmed.fastq.gz")]
        cmd = base + ["-o", out[0], reads[0]]
    log = run(cmd).stdout
    with open(os.path.join(outdir, "cutadapt.log"), "w") as f:
        f.write(log)
    total = re.search(r"Total (?:read pairs|reads) processed:\s+([\d,]+)", log)
    kept = re.search(r"(?:Pairs|Reads) written \(passing filters\):\s+([\d,]+)", log)
    n_in = int(total.group(1).replace(",", "")) if total else 0
    n_out = int(kept.group(1).replace(",", "")) if kept else 0
    return out, {"reads_in": n_in, "reads_kept": n_out, "pct_retained": 100.0 * n_out / n_in if n_in else 0.0}


def index_genome(genome, outdir):
    """bwa-index a symlink inside <outdir>/reference so the reference folder stays untouched."""
    ref_dir = os.path.join(outdir, "reference")
    os.makedirs(ref_dir, exist_ok=True)
    link = os.path.join(ref_dir, "genome.fa")
    if not os.path.exists(link):
        os.symlink(os.path.abspath(genome), link)
    if not all(os.path.exists(link + ext) for ext in (".bwt", ".pac", ".ann", ".amb", ".sa")):
        run(["bwa", "index", link])
    return link


def align(sample, reads, genome, outdir, threads):
    """bwa mem piped into samtools sort; returns the indexed BAM path."""
    bam = os.path.join(outdir, f"{sample}.sorted.bam")
    rg = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"
    bwa = subprocess.Popen(["bwa", "mem", "-t", str(threads), "-R", rg, genome, *reads],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    srt = subprocess.Popen(["samtools", "sort", "-@", str(threads), "-o", bam, "-"], stdin=bwa.stdout,
                           stderr=subprocess.PIPE)
    bwa.stdout.close()
    sort_err = srt.communicate()[1]
    bwa_err = bwa.stderr.read()
    bwa.wait()
    if bwa.returncode or srt.returncode:
        raise RuntimeError(f"alignment failed: {(bwa_err or sort_err).decode()[-500:]}")
    run(["samtools", "index", bam])
    return bam


def analyze_sample(sample, reads, genome_index, bed, fwd, rev, outdir, threads):
    sdir = os.path.join(outdir, sample)
    os.makedirs(sdir, exist_ok=True)
    trimmed, trim_stats = trim(sample, reads, fwd, rev, sdir, threads)
    if trim_stats["reads_kept"] == 0:
        raise RuntimeError(f"{sample}: no reads started with a primer from the primer FASTA (all discarded). "
                           "Check that --primers belongs to this pool and that reads begin at the primer.")
    bam = align(sample, trimmed, genome_index, sdir, threads)
    targets, depths, stats = metrics.measure(bam, bed)
    panel = metrics.panel_metrics(targets, depths, stats)
    feedback = balance.build_feedback(panel, targets)

    with open(os.path.join(sdir, "target_coverage.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["name", "chrom", "start", "end", "mean_coverage"])
        for t in targets:
            w.writerow([t["name"], t["chrom"], t["start"], t["end"], f"{t['mean_coverage']:.3f}"])
    with open(os.path.join(sdir, f"{sample}_balancing_feedback.json"), "w") as f:
        json.dump(feedback, f, indent=2)
    with open(os.path.join(sdir, f"{sample}_target_summary.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["target", "mean_coverage", "rel_mean", "adjustment", "action"])
        for r in feedback["target_recommendations"]:
            w.writerow([r["target_id"], f"{r['mean_coverage']:.1f}", f"{r['rel_mean']:.3f}",
                        f"{r['recommended_adjustment']:.2f}", r["action"]])
    plotted = plot_uniformity.plot_uniformity(targets, os.path.join(sdir, f"{sample}_uniformity.png"), sample)
    summary = {"sample": sample, **trim_stats, **panel, "uniformity_0.2x": feedback["overall_metrics"]["uniformity_0.2x"],
               "uniformity_plot": plotted}
    with open(os.path.join(sdir, f"{sample}_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    return summary


def analyze(reads, genome, bed, primers, outdir, threads=None, sample=None):
    """Analyse every sample matching `reads`. Returns a list of summary dicts."""
    check_tools()
    for path, what in ((genome, "--genome"), (bed, "--bed"), (primers, "--primers")):
        if not os.path.exists(path):
            raise FileNotFoundError(f"{what} not found: {path}")
    samples = find_samples(reads)
    if not samples:
        raise FileNotFoundError(f"no reads match {reads}")
    if sample and len(samples) != 1:
        raise ValueError("--sample can only rename a single sample")
    threads = threads or os.cpu_count() or 1
    os.makedirs(outdir, exist_ok=True)
    fwd, rev = split_primers(primers, outdir)
    genome_index = index_genome(genome, outdir)
    return [analyze_sample(sample or name, files, genome_index, bed, fwd, rev, outdir, threads)
            for name, files in samples]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reads", required=True, help="FASTQ glob; use {1,2} for the read number")
    ap.add_argument("--genome", required=True, help="reference FASTA the amplicons were designed on")
    ap.add_argument("--bed", required=True, help="amplicon BED from design.py (<prefix>_pool_N.bed)")
    ap.add_argument("--primers", required=True, help="trimming FASTA from design.py (<prefix>_pool_N_trimming.fasta)")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--threads", type=int)
    ap.add_argument("--sample", help="override the sample name (single sample only)")
    a = ap.parse_args()
    try:
        for s in analyze(a.reads, a.genome, a.bed, a.primers, a.outdir, a.threads, a.sample):
            f80 = f"{s['fold_80']:.2f}" if s["fold_80"] is not None else "n/a"
            print(f"{s['sample']}: {s['pct_retained']:.1f}% reads kept after trimming, {s['pct_reads_mapped']:.1f}% mapped, "
                  f"mean target coverage {s['mean_target_coverage']:.0f}x, fold-80 {f80}, "
                  f"{s['uniformity_0.2x']:.1f}% of {s['n_targets']} amplicons > 0.2x mean -> {os.path.join(a.outdir, s['sample'])}/")
    except (RuntimeError, ValueError, FileNotFoundError) as e:
        sys.exit(f"error: {e}")


if __name__ == "__main__":
    main()
