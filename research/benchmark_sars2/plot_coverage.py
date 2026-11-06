"""Figure: where each whole-genome SARS-CoV-2 design places its amplicons and what it leaves uncovered.

  python plot_coverage.py <out.png> <genome.fa> artic <amplicons.bed> <primers.fasta> \
      design "<label>" <prefix> <targets.tsv>  [design "<label>" <prefix> <targets.tsv> ...]
Pools alternate blue/aqua; orange bands are bases no amplicon insert covers after primer trimming.
"""
import csv, json, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coverage_report as cr

BLUE, AQUA, ORANGE, INK, MUTED, GRID, SURFACE = "#2a78d6", "#1baf7a", "#eb6834", "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

out, fa = sys.argv[1], sys.argv[2]
L = cr.genome_length(fa)
rows = []                                    # (label, amplicons [(s, e, fl, rl, pool)], analysis)
args = sys.argv[3:]
if args[0] == "artic":
    amps = cr.read_artic(args[1], args[2])
    pools = [i % 2 for i, _ in enumerate(sorted(amps))]
    by = sorted(amps)
    rows.append(("ARTIC v3 scheme", [(*a[:4], i % 2) for i, a in enumerate(by)], cr.analyse(amps, L)))
    args = args[3:]
while args:
    assert args[0] == "design"
    label, prefix, targets = args[1:4]
    amps, pool_of = [], {}
    import glob
    for pi, f in enumerate(sorted(glob.glob(f"{prefix}_pool_*.csv"))):
        for r in csv.DictReader(open(f)):
            amps.append((int(r["start"]), int(r["end"]), len(r["fwd_primer"]), len(r["rev_primer"]), r["target_id"]))
            pool_of[r["target_id"]] = pi
    m = json.load(open(f"{prefix}_manifest.json"))
    wins = {r["name"]: (int(r["start"]) - 1, int(r["end"])) for r in csv.DictReader(open(targets), delimiter="\t")}
    ok = {t for p in m["results"]["pools"] for t in p["targets"]}
    failed = [wins[n] for n in wins if n not in ok]
    res = cr.analyse(amps, L, failed)
    rows.append((label, [(s, e, fl, rl, pool_of[t]) for s, e, fl, rl, t in amps], res))
    args = args[4:]

fig, axes = plt.subplots(len(rows), 1, figsize=(12, 1.45 * len(rows) + 0.9), sharex=True, facecolor=SURFACE)
for ax, (label, amps, res) in zip(axes, rows):
    ax.set_facecolor(SURFACE)
    for a, b in [x for d in res["detail"].values() for x in d]:
        ax.axvspan(a, b, color=ORANGE, alpha=0.9, lw=0)
    for s, e, fl, rl, pool in amps:
        y = 0.36 if pool == 0 else 0.03
        ax.add_patch(plt.Rectangle((s + fl, y), (e - rl) - (s + fl), 0.28, color=BLUE if pool == 0 else AQUA, lw=0))
        ax.add_patch(plt.Rectangle((s, y + 0.08), fl, 0.12, color=MUTED, lw=0))
        ax.add_patch(plt.Rectangle((e - rl, y + 0.08), rl, 0.12, color=MUTED, lw=0))
    ax.set_xlim(0, L); ax.set_ylim(0, 1.0)
    ax.set_yticks([])
    ax.text(0.003, 1.03, f"{label}", transform=ax.transAxes, va="bottom", ha="left", fontsize=10.5, color=INK, weight="bold")
    ax.text(0.997, 1.03, f"{res['pct_covered']:.1f}% covered · {res['n_amplicons']} amplicons · {res['uncovered']:,} bp uncovered",
            transform=ax.transAxes, va="bottom", ha="right", fontsize=9.5, color=MUTED)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9)
axes[-1].set_xlabel("SARS-CoV-2 position (bp)", color=MUTED, fontsize=10)
fig.suptitle("Whole-genome SARS-CoV-2 panels: amplicon inserts by pool (blue / aqua), primers (grey), uncovered bases (orange)",
             x=0.012, ha="left", fontsize=11, color=INK)
fig.tight_layout(h_pad=1.4, rect=[0, 0, 1, 0.945])
fig.savefig(out, dpi=200, facecolor=SURFACE)
print("wrote", out)
