"""Coverage along the genome (log scale) with the amplicon layout underneath."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

DEEP, BLUE, AQUA, INK, MUTED, GRID, SURFACE = "#0d4f87", "#2a78d6", "#1baf7a", "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def plot_coverage_track(profile, targets, output_path, sample_name, max_contigs=4):
    """profile = (contigs, positions, depths) per target base; targets = [{'name','chrom','start','end'}, ...].
    One panel per contig (up to max_contigs, longest coverage first). Returns False if nothing to plot."""
    contigs, pos, depth = profile
    if not len(depth) or depth.max() == 0:
        return False
    names = sorted(set(contigs), key=lambda c: -int((contigs == c).sum()))[:max_contigs]
    fig, axes = plt.subplots(2 * len(names), 1, figsize=(11, 2.9 * len(names) + 0.6), facecolor=SURFACE,
                             gridspec_kw={"height_ratios": [3, 0.55] * len(names)})
    axes = np.atleast_1d(axes)
    for i, contig in enumerate(names):
        top, band = axes[2 * i], axes[2 * i + 1]
        m = contigs == contig
        x, y = pos[m], np.maximum(depth[m], 1)
        lo, hi = int(x.min()), int(x.max())
        top.set_facecolor(SURFACE)
        top.fill_between(x, 1, y, step="mid", color=DEEP, linewidth=0)
        top.set_yscale("log")
        top.set_ylim(1, max(10, y.max() * 1.6))
        top.set_xlim(lo - 0.01 * (hi - lo), hi + 0.01 * (hi - lo))
        top.set_ylabel("Log$_{10}$ coverage", color=MUTED, fontsize=11)
        top.tick_params(colors=MUTED, labelbottom=False)
        top.grid(axis="y", which="major", color=GRID, linewidth=0.6)
        top.set_axisbelow(True)
        top.set_title(f"{sample_name}: {contig}" if len(names) > 1 or i == 0 else contig, loc="left", fontsize=11, color=INK)
        band.set_facecolor(SURFACE)
        for k, t in enumerate((t for t in targets if t["chrom"] == contig)):
            band.add_patch(plt.Rectangle((t["start"], 0.55 if k % 2 == 0 else 0.1), t["end"] - t["start"], 0.32,
                                         color=BLUE if k % 2 == 0 else AQUA, linewidth=0))
        band.set_xlim(top.get_xlim()); band.set_ylim(0, 1); band.set_yticks([])
        band.set_xlabel("Position", color=MUTED, fontsize=11)
        band.tick_params(colors=MUTED)
        for ax in (top, band):
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            for sp in ("left", "bottom"):
                ax.spines[sp].set_color(MUTED)
        band.spines["left"].set_visible(False)
    fig.tight_layout(h_pad=0.4)
    fig.savefig(output_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return True
