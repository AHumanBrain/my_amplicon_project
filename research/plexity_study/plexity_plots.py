"""Draw the plexity-study figures from results.json (and offtarget_*.json when present).

    python research/plexity_study/plexity_plots.py --workdir <study dir> --outdir <figure dir>

Colours follow the dataviz palette (validated): blue = naive primers, aqua = brute-force screened (repair),
orange = brute-force first-fit; ordinal blue ramp = chemistry ladder steps (lighter = less chemistry).
"""
import argparse
import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
LADDER_BLUE = {"L0": "#86b6ef", "L1": "#2a78d6", "L2": "#1c5cab", "L3": "#184f95", "L4": "#0d366b"}
LADDER_NAME = {"L0": "plain: no hot-start, 55 °C, 250 nM", "L1": "+ hot-start", "L2": "+ 60 °C annealing",
               "L3": "+ 50 nM primers", "L4": "+ rebalanced conc."}
MARK = {"L0": "o", "L1": "s", "L2": "D", "L3": "^", "L4": "v"}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": GRID, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 10, "axes.titlesize": 11,
    "axes.titleweight": "semibold", "axes.titlelocation": "left", "legend.frameon": False,
})


def series(res, scenario, level, key):
    rows = res["scenarios"][scenario]
    return [r["N"] for r in rows], [r["levels"][level][key] for r in rows]


def line(ax, x, y, color, marker, label, ls="-", lw=2.0, ms=6, ring=True):
    ax.plot(x, y, ls=ls, lw=lw, color=color, marker=marker, ms=ms, label=label,
            markeredgecolor=SURFACE if ring else color, markeredgewidth=1.2, clip_on=False)


def end_label(ax, x, y, text, color, dy=1.0, dx=1.06, fs=8.5):
    ax.annotate(text, xy=(x, y), xytext=(x * dx, y * dy), color=INK2, fontsize=fs, va="center",
                annotation_clip=False)


def plain_xticks(ax):
    """Readable plexity ticks on the log x axis (10, 30, 100, 300, 1000)."""
    ax.set_xticks([10, 30, 100, 300, 1000])
    ax.set_xticklabels(["10", "30", "100", "300", "1,000"])
    ax.minorticks_off()


def finish(fig, name, outdir, caption):
    fig.text(0.01, 0.005, caption, fontsize=7.6, color=INK2, va="bottom", ha="left", wrap=True)
    path = os.path.join(outdir, name)
    fig.savefig(path, dpi=190)
    plt.close(fig)
    return path


def fig_dimers(res, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.6))
    fig.subplots_adjust(left=0.06, right=0.84, top=0.86, bottom=0.2, wspace=0.52)
    panels = [("risky_pairs", "Extendable-dimer pairs in the pools", "pairs (count)", "symlog"),
              ("mean_seq", "Mean fraction of primer 3′ ends tied up in dimers", "% of primers", "log")]
    for ax, (key, title, ylab, scale) in zip(axes, panels):
        mult = 100.0 if key == "mean_seq" else 1.0
        for lv, color, label in (("L0", LADDER_BLUE["L0"], "naive, plain"), ("L1", LADDER_BLUE["L1"], "naive + hot-start"),
                                 ("L3", LADDER_BLUE["L4"], "naive + full chemistry")):
            x, y = series(res, "naive", lv, key)
            y = np.array(y) * mult
            if key == "mean_seq":
                y = np.maximum(y, 0.0008)
            line(ax, x, y, color, MARK[lv], label)
            end_label(ax, x[-1], y[-1], label, color)
        x, y = series(res, "repair_plain", "L0", key)
        y = np.array(y) * mult
        line(ax, x, np.maximum(y, 0.0008) if key == "mean_seq" else y, AQUA, "s", "brute-force screened, plain", ls="--")
        end_label(ax, x[-1], max(y[-1], 0.0008), "brute-force screened, plain", AQUA, dy=1.0)
        ax.set_xscale("log")
        ax.set_yscale(scale, **({"linthresh": 1.0} if scale == "symlog" else {}))
        ax.set_xlim(10, 1000)
        plain_xticks(ax)
        ax.set_xlabel("attempted plexity (targets, tiled design, 2 pools)")
        ax.set_ylabel(ylab)
        ax.set_title(title)
        if key == "risky_pairs":
            n2 = np.array([10, 1000])
            ax.plot(n2, 155 * (n2 / 100.0) ** 2, color=INK2, lw=1, ls=":", zorder=1)
            ax.annotate("∝ N²", xy=(300, 155 * 9), color=INK2, fontsize=8.5)
        else:
            ax.set_ylim(0.0008, 120)
            ax.set_yticks([0.001, 0.01, 0.1, 1, 10, 100])
            ax.set_yticklabels(["0.001", "0.01", "0.1", "1", "10", "100"])
    fig.suptitle("Primer-dimer potential grows ~N²: screening flattens it by discarding targets, chemistry flattens it without losing any",
                 x=0.06, ha="left", y=0.985, fontsize=13, fontweight="semibold")
    return finish(fig, "fig1_dimer_potential_vs_plexity.png", outdir,
                  "Model: nearest-neighbour 3′-end ΔG (primer3, KAPA salts) → pair occupancy c/(c+Kd). 'Risky' = occupancy ≥ 0.42% (the tool's −6 kcal/mol at 37 °C, 250 nM). "
                  "Plain = dimers evaluated at 37 °C (no hot-start), 250 nM per primer at EVERY plexity, a worst case. Zero values drawn at the axis floor.\n"
                  "Brute-force screened = best primer pair per tile, chosen to keep every pair under the plain limit (only 290/935 tiles can be placed, see Fig 2). In silico only.")


def fig_capacity(res, outdir):
    fig, ax = plt.subplots(figsize=(8.8, 5.4))
    fig.subplots_adjust(left=0.09, right=0.70, top=0.86, bottom=0.25)
    resc = [100.0 * res["scenarios"][sc][-1]["placed"] / res["scenarios"][sc][-1]["N"]
            for sc in ("firstfit_rescue", "repair_rescue") if sc in res["scenarios"]]
    ax.axhline(100, color=INK2, lw=1.2, ls="--", zorder=1)
    ax.text(11.2, 103.2, f"hot-start / 60 °C / 50 nM limit: {min(resc):.1f}–{max(resc):.0f}% placed, both searches" if resc else "",
            color=INK2, fontsize=8.8, va="bottom")
    specs = [("firstfit_plain", ORANGE, "o", "first-fit"), ("repair_plain", AQUA, "s", "first-fit + repair")]
    for sc, color, mk, label in [sp for sp in specs if sp[0] in res["scenarios"]]:
        rows = res["scenarios"][sc]
        x = [r["N"] for r in rows]
        y = [100.0 * r["placed"] / r["N"] for r in rows]
        line(ax, x, y, color, mk, label + ", plain limit (−6 kcal/mol, 37 °C, 250 nM)")
        dy = 5.5 if sc == "repair_plain" else -5.5          # keep the two end labels apart
        ax.annotate(f"{label}: {y[-1]:.0f}% placed", xy=(x[-1], y[-1]), xytext=(x[-1] * 1.07, y[-1] + dy), color=INK2,
                    fontsize=8.8, va="center", annotation_clip=False)
    ax.set_xscale("log")
    ax.set_xlim(10, 1000)
    plain_xticks(ax)
    ax.set_ylim(0, 112)
    ax.set_xlabel("attempted plexity (tiled tiles)")
    ax.set_ylabel("% of targets that can be placed")
    ax.set_title("Under the plain dimer limit a brute-force search places only ~30% of a 1,000-plex tiling", fontsize=10.5)
    ax.legend(loc="lower left", fontsize=8.5, bbox_to_anchor=(0.0, -0.02))
    return finish(fig, "fig2_capacity_vs_plexity.png", outdir,
                  "Greedy per-pool selection from 25 candidates/tile (genome order); 'repair' moves one blocking member to an alternative. Tiled windows leave few distinct primers per tile\n"
                  "(tile 11: 10 left / 3 right primers sharing one motif), so one early pick can block a neighbour. A joint search could place more. In silico only.")


def fig_uniformity(res, outdir):
    fig, (ax, zx) = plt.subplots(1, 2, figsize=(13.4, 5.6), gridspec_kw={"width_ratios": [1.0, 1.0]})
    fig.subplots_adjust(left=0.065, right=0.86, top=0.84, bottom=0.22, wspace=0.42)
    sets = [("naive", "L0", LADDER_BLUE["L0"], "o", "-", "naive, plain"), ("naive", "L1", LADDER_BLUE["L1"], "s", "-", "naive + hot-start"),
            ("naive", "L3", LADDER_BLUE["L3"], "^", "-", "naive + 50 nM primers"), ("naive", "L4", LADDER_BLUE["L4"], "v", "-", "naive + rebalanced conc."),
            ("repair_plain", "L0", AQUA, "s", "--", "brute-force screened, plain")]
    zoom_dy = {"naive + hot-start": -1.1, "naive + 50 nM primers": 0.0, "naive + rebalanced conc.": 0.9, "brute-force screened, plain": 2.2}
    for sc, lv, color, mk, ls, label in sets:
        x, y = series(res, sc, lv, "uniform_frac_ge_0p2")
        y = np.array(y) * 100
        for a_ in (ax, zx):
            line(a_, x, y, color, mk, label, ls=ls)
        if label == "naive, plain":
            ax.annotate(label, xy=(x[-1], y[-1]), xytext=(x[-1] * 1.07, y[-1]), color=INK2, fontsize=8.6, va="center", annotation_clip=False)
        else:
            zx.annotate(label, xy=(x[-1], y[-1]), xytext=(x[-1] * 1.07, y[-1] + zoom_dy[label]), color=INK2, fontsize=8.6, va="center", annotation_clip=False)
    ax.set_ylim(0, 105)
    zx.set_ylim(85, 104)
    for a_ in (ax, zx):
        a_.set_xscale("log")
        a_.set_xlim(10, 1000)
        plain_xticks(a_)
        a_.set_xlabel("attempted plexity (tiled tiles)")
        a_.set_ylabel("% of amplicons ≥ 0.2× mean coverage (model)")
    ax.set_title("Full range")
    zx.set_title("Zoom on the 85–100% band")
    ax.legend(loc="lower left", fontsize=8, bbox_to_anchor=(0.0, -0.02))
    fig.suptitle("Uniformity collapses with dimers; stringency costs some back until primers are rebalanced", x=0.065, ha="left", y=0.985,
                 fontsize=13, fontweight="semibold")
    return finish(fig, "fig3_uniformity_vs_plexity.png", outdir,
                  "MODEL, uncalibrated: per-cycle efficiency e = 0.95 × target occupancy(F) × occupancy(R) × (1−dimer sequestration F)(1−dimer sequestration R); yield ∝ (1+e)^30; no plateau compression. "
                  "Only relative comparisons are meaningful.\n"
                  "Rebalancing sets each primer's concentration (10–250 nM) to equalise target occupancy, as the coverage-feedback loop does empirically. Brute-force screened places only 290/935 tiles (Fig 2).")


def fig_ladder(res, outdir):
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.8))
    fig.subplots_adjust(left=0.2, right=0.97, top=0.84, bottom=0.24, wspace=0.22)
    N = res["plex"][-1]
    entry = [e for e in res["scenarios"]["naive"] if e["N"] == N][0]
    keys = ["L0", "L1", "L2", "L3", "L4"]
    metrics = [("risky_pairs", "Extendable-dimer pairs", "symlog", 1.0), ("mean_seq", "Mean primer sequestered (%)", "log", 100.0),
               ("uniform_frac_ge_0p2", "Amplicons ≥ 0.2× mean (%)", "linear", 100.0)]
    for i, (ax, (key, title, scale, mult)) in enumerate(zip(axes, metrics)):
        vals = [max(entry["levels"][k][key] * mult, 0.002 if key == "mean_seq" else 0) for k in keys]
        ypos = np.arange(len(keys))[::-1]
        for yp, k, v in zip(ypos, keys, vals):
            ax.plot([ax.get_xlim()[0] if scale == "linear" else 1e-4 if scale == "log" else 0, v], [yp, yp], color=GRID, lw=3, solid_capstyle="round", zorder=1)
            ax.plot(v, yp, "o", ms=11, color=LADDER_BLUE[k], markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=3)
            txt = f"{v:,.0f}" if key == "risky_pairs" else (f"{v:.2g}" if key == "mean_seq" else f"{v:.0f}")
            ax.annotate(txt, xy=(v, yp), xytext=(8, 0), textcoords="offset points", fontsize=8.7, color=INK, va="center")
        ax.set_xscale(scale, **({"linthresh": 1.0} if scale == "symlog" else {}))
        if scale == "linear":
            ax.set_xlim(0, 115)
        elif scale == "log":
            ax.set_xlim(1e-3, 300)
        else:
            ax.set_xlim(0, 4e4)
        ax.set_ylim(-0.6, len(keys) - 0.4)
        ax.set_yticks(ypos)                                   # after the scales, or matplotlib resets the ticks
        ax.set_yticklabels([LADDER_NAME[k] for k in keys] if i == 0 else [""] * len(keys))
        ax.set_title(title)
        ax.grid(axis="y", visible=False)
    fig.suptitle(f"The chemistry ladder at {N}-plex (naive primers, 2 tiled pools): hot-start does most of the work",
                 x=0.02, ha="left", y=0.985, fontsize=13, fontweight="semibold")
    return finish(fig, "fig4_rescue_ladder.png", outdir,
                  "Each step adds one lever to the previous. 'plain' evaluates dimers at 37 °C (reaction setup, no hot-start). Most commercial master mixes are already hot-start, so the realistic baseline is the 2nd row.\n"
                  "These levers are commodity reagents/cycling choices; I cannot vouch for the IP status of any of them. Predictions of potential, not measurements.")


def own_product_counts(workdir, scenario, plex, max_product=1000):
    """Products formed purely by facing on-target (designed) primers of different targets in one pool.

    In a tiled 2-pool layout these are geometric (A's forward primer faces C's reverse primer ~600 bp away)
    and do not depend on chemistry, so they are subtracted from the off-target product counts.
    """
    import pickle
    import plexity_study as ps
    cands = pickle.load(open(os.path.join(workdir, "candidates.pkl"), "rb"))["cands"]
    if scenario == "naive":
        order = sorted(cands)
        chosen = {t: 0 for t in order}
    else:
        tag, rule = scenario.split("_")
        sel = pickle.load(open(os.path.join(workdir, "selection_repair.pkl" if tag == "repair" else "selection.pkl"), "rb"))
        order = sel["order"]
        chosen = {t: r for t, r in sel[f"screened_{rule}"].items() if r is not None}
    out = {}
    for n in plex:
        cutoff = order[n - 1]
        total = 0
        for pool in (0, 1):
            sites = []
            for t in chosen:
                if t <= cutoff and (int(t[1:]) - 1) % 2 == pool:
                    f, r = ps.design.priming_footprints(cands[t][chosen[t]])
                    for fp in (f, r):
                        sites.append(dict(contig="NC_000913.3", lo=fp[0], hi=fp[1], strand=fp[2], target=t, own=True))
            total += ps._count_products(sites, max_product)
        out[n] = total
    return out


def fig_offtarget(res, outdir, workdir):
    data = {}
    for sc in ("naive", "repair_plain"):
        p = os.path.join(workdir, f"offtarget_{sc}.json")
        if os.path.exists(p):
            data[sc] = json.load(open(p))
    if "naive" not in data:
        return None
    plex = [r["N"] for r in data["naive"]]
    own = {sc: own_product_counts(workdir, sc, plex) for sc in data}
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.4))
    fig.subplots_adjust(left=0.06, right=0.84, top=0.86, bottom=0.2, wspace=0.52)
    for ax, (key, title, ylab) in zip(axes, [("competitive_sites", "Competitive off-target priming sites", "sites (count)"),
                                               ("unintended_products", "Off-target-driven products (≤ 1 kb)", "products (count)")]):
        def val(sc, r, lv):
            v = r["levels"][lv][key]
            return v - own[sc][r["N"]] if key == "unintended_products" else v
        for lv in ("L0", "L1", "L3"):
            x = [r["N"] for r in data["naive"]]
            y = [val("naive", r, lv) for r in data["naive"]]
            color = LADDER_BLUE["L0"] if lv == "L0" else LADDER_BLUE["L1"] if lv == "L1" else LADDER_BLUE["L4"]
            label = {"L0": "naive, plain", "L1": "naive + hot-start", "L3": "naive + full chemistry"}[lv]
            line(ax, x, y, color, MARK[lv], label)
            end_label(ax, x[-1], max(y[-1], 0.5), label, color)
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=1.0)
        ax.set_xlim(10, 1000)
        ax.set_ylim(bottom=0)
        plain_xticks(ax)
        ax.set_xlabel("attempted plexity (tiled tiles)")
        ax.set_ylabel(ylab)
        ax.set_title(title)
    fig.suptitle("Off-target priming: hot-start cuts competitive sites ~3×, the full ladder ~28× (935-plex)", x=0.06, ha="left", y=0.985, fontsize=13, fontweight="semibold")
    return finish(fig, "fig5_offtarget_vs_plexity.png", outdir,
                  "Sites: BLAST near-matches (≥14 aligned bp incl. exact 3′ end, ≤2 mismatches) scored by primer–site duplex ΔG; 'competitive' = occupancy ≥ 5% of the primer's real target at that temperature and concentration.\n"
                  "Products = facing sites in one pool within 1 kb that involve at least one off-target site; the geometric products between neighbouring designed primers (A faces C, ~600 bp) are subtracted. Dimer screening does not select against off-target sites, so it is not shown. Genome: E. coli K-12 only (4.6 Mb); a human-scale genome has ~1000× more candidate sites. In silico only.")


def write_table(res, outdir):
    path = os.path.join(outdir, "plexity_study_table.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["scenario", "plexity", "placed", "ladder", "risky_pairs", "mean_sequestration_pct", "p95_sequestration_pct",
                    "pct_amplicons_ge_0.2x_mean"])
        for sc, rows in res["scenarios"].items():
            for r in rows:
                for lv, v in r["levels"].items():
                    w.writerow([sc, r["N"], r["placed"], lv, v["risky_pairs"], f"{v['mean_seq'] * 100:.4f}", f"{v['p95_seq'] * 100:.4f}",
                                f"{v['uniform_frac_ge_0p2'] * 100:.1f}"])
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    res = json.load(open(os.path.join(a.workdir, "results.json")))
    for fn in os.listdir(a.workdir):                       # newest per-scenario caches win over results.json
        if fn.startswith("scenario_") and fn.endswith(".json"):
            res["scenarios"][fn[len("scenario_"):-5]] = json.load(open(os.path.join(a.workdir, fn)))
    paths = [fig_dimers(res, a.outdir), fig_capacity(res, a.outdir), fig_uniformity(res, a.outdir), fig_ladder(res, a.outdir)]
    off = fig_offtarget(res, a.outdir, a.workdir)
    if off:
        paths.append(off)
    paths.append(write_table(res, a.outdir))
    print("\n".join(paths))


if __name__ == "__main__":
    main()
