"""Whole-genome coverage of a designed panel once primers are trimmed, and why bases are uncovered.

  python coverage_report.py design  <prefix> <targets.tsv> <genome.fa>   a design.py run (--output-prefix <prefix>)
  python coverage_report.py artic   <amplicons.bed> <primers.fasta> <genome.fa>   the ARTIC scheme, scored the same way

A base counts as covered if some amplicon's INSERT (amplicon minus its two primers) contains it, which is what is
left after primer trimming. Uncovered stretches are labelled: genome end, failed tile (a target that was not
designed), or overlap gap (neighbouring amplicons overlap by less than two primer lengths).
"""
import csv, glob, json, sys


def genome_length(fa):
    return sum(len(l.strip()) for l in open(fa) if not l.startswith(">"))


def read_design(prefix):
    amps = []
    for f in sorted(glob.glob(f"{prefix}_pool_*.csv")):
        for r in csv.DictReader(open(f)):
            amps.append((int(r["start"]), int(r["end"]), len(r["fwd_primer"]), len(r["rev_primer"]), r["target_id"]))
    return amps


def read_artic(bed, fasta):
    plen, k = {}, None
    for l in open(fasta):
        l = l.strip()
        if l.startswith(">"):
            k = l[1:]; plen[k] = 0
        else:
            plen[k] += len(l)
    out = []
    for l in open(bed):
        c, s, e, n = l.split()[:4]
        out.append((int(s), int(e), plen[f"{n}_F"], plen[f"{n}_R"], n))
    return out


def analyse(amps, L, failed_windows=()):
    cov = [0] * L
    for s, e, fl, rl, _ in amps:
        for i in range(max(s + fl, 0), min(e - rl, L)):
            cov[i] += 1
    stretches, i = [], 0
    while i < L:
        if cov[i] == 0:
            j = i
            while j < L and cov[j] == 0:
                j += 1
            stretches.append((i, j)); i = j
        else:
            i += 1
    causes = {"genome end": 0, "failed tile": 0, "overlap gap": 0}
    detail = {"genome end": [], "failed tile": [], "overlap gap": []}
    for a, b in stretches:
        if a == 0 or b == L:
            c = "genome end"
        elif any(max(a, ws) < min(b, we) for ws, we in failed_windows):
            c = "failed tile"
        else:
            c = "overlap gap"
        causes[c] += b - a
        detail[c].append((a, b))
    unc = sum(b - a for a, b in stretches)
    return dict(n_amplicons=len(amps), genome=L, uncovered=unc, pct_covered=100 * (L - unc) / L,
                stretches=len(stretches), causes=causes, detail=detail,
                redundant=sum(1 for c in cov if c >= 2), lengths=sorted(e - s for s, e, *_ in amps))


def show(name, r):
    c = r["causes"]
    print(f"{name}: {r['n_amplicons']} amplicons cover {r['pct_covered']:.1f}% of {r['genome']} bp "
          f"({r['uncovered']} bp uncovered in {r['stretches']} stretches); uncovered bp by cause: "
          f"failed tile {c['failed tile']}, overlap gap {c['overlap gap']}, genome end {c['genome end']}; "
          f"{100 * r['redundant'] / r['genome']:.0f}% of bases covered by 2+ amplicons")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "design":
        prefix, targets, fa = sys.argv[2:5]
        m = json.load(open(f"{prefix}_manifest.json"))
        wins = {r["name"]: (int(r["start"]) - 1, int(r["end"])) for r in csv.DictReader(open(targets), delimiter="\t")}
        ok = {t for p in m["results"]["pools"] for t in p["targets"]}
        failed = [wins[n] for n in wins if n not in ok]
        show(prefix.split("/")[-1], analyse(read_design(prefix), genome_length(fa), failed))
        print(f"  tiles: {len(ok)} designed, {len(failed)} not designed")
    else:
        bed, fasta, fa = sys.argv[2:5]
        show("ARTIC v3", analyse(read_artic(bed, fasta), genome_length(fa)))
