"""Normalized coverage per amplicon (log scale) against reference lines at multiples of the panel mean."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BLUE, ORANGE, INK, MUTED, GRID, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
# (multiple of the mean, colour, linestyle, linewidth)
LINES = ((5.0, "#2a78d6", "--", 2.2), (1.0, "#8a4b2a", "--", 1.6), (0.2, "#d6342a", "--", 2.2),
         (0.1, "#1baf7a", ":", 1.8), (0.05, "#e8a0b4", "--", 1.8))


def plot_uniformity(targets, output_path, sample_name):
    """`targets` need a 'rel_mean' (see balance.calculate_relative_coverage). Returns False if nothing to plot."""
    rel = [t["rel_mean"] for t in targets]
    if not rel or max(rel) == 0:
        return False
    pct = sum(r > 0.2 for r in rel) / len(rel) * 100
    fig, ax = plt.subplots(figsize=(10, 5.2), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.scatter(range(len(rel)), [max(r, 0.03) for r in rel], s=22, alpha=0.85, linewidths=0,
               c=[BLUE if r > 0.2 else ORANGE for r in rel], zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(0.035, 8)
    for y, color, style, lw in LINES:
        ax.axhline(y, linestyle=style, linewidth=lw, color=color, zorder=2)
        ax.text(len(rel) * 0.012, y * 1.07, f"{y:g} * mean cov", fontsize=9, color=INK, va="bottom")
    ax.grid(axis="y", which="both", color=GRID, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.set_xticks([])
    ax.set_xlabel(f"Amplicons 1-{len(rel)}", color=MUTED, fontsize=12)
    ax.set_ylabel("Normalized coverage per amplicon", color=MUTED, fontsize=12)
    ax.set_title(f"{sample_name}: {pct:.1f}% of amplicons > 0.2x mean", loc="left", fontsize=11, color=INK)
    ax.tick_params(colors=MUTED)
    for sp in ax.spines.values():
        sp.set_color(MUTED)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return True
