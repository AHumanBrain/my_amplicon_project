"""Log-scale scatter of per-amplicon coverage relative to the panel mean."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def plot_uniformity(targets, output_path, sample_name):
    """`targets` need a 'rel_mean' (see balance.calculate_relative_coverage). Returns False if nothing to plot."""
    rel = [t["rel_mean"] for t in targets]
    if not rel or max(rel) == 0:
        return False
    pct = sum(r > 0.2 for r in rel) / len(rel) * 100
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.scatter(range(len(rel)), [max(r, 1e-3) for r in rel], s=18, c=["#2563eb" if r > 0.2 else "#dc2626" for r in rel])
    ax.set_yscale("log")
    for y, style in ((1.0, "--"), (0.2, "-"), (5.0, "--")):
        ax.axhline(y, linestyle=style, linewidth=1, color="gray")
    ax.set_title(f"{sample_name}: {pct:.1f}% of amplicons > 0.2x mean")
    ax.set_xlabel("Amplicon index")
    ax.set_ylabel("Relative coverage (log)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return True
