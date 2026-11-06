"""Do our thermodynamic metrics predict which real ARTIC v3 amplicons drop out? usage: predict_dropouts.py <workdir of scripts/real_data_check.sh>"""
import csv, sys
import numpy as np, primer3
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from core import biophysics as bp
W = sys.argv[1]   # workdir of scripts/real_data_check.sh (holds out/sample1/target_coverage.tsv)
EX = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'examples', 'sars2')
ref = "".join(l.strip() for l in open(EX + "/NC_045512.2.fa") if l[0] != ">").upper()
prim = {}
k = None
for l in open(EX + "/artic_v3_primers.fasta"):
    l = l.strip()
    if l[0] == ">": k = l[1:]; prim[k] = ""
    else: prim[k] += l
cov = {r["name"]: float(r["mean_coverage"]) for r in csv.DictReader(open(f"{W}/out/sample1/target_coverage.tsv"), delimiter="\t")}
bed = {l.split()[3]: (int(l.split()[1]), int(l.split()[2])) for l in open(EX + "/artic_v3_amplicons.bed")}
tm = lambda s: primer3.calc_tm(s, mv_conc=50, dv_conc=2.5, dntp_conc=1.2, dna_conc=50)
gc = lambda s: (s.count("G") + s.count("C")) / len(s)
pool = lambda n: int(n[3:].split("_")[0]) % 2        # ARTIC convention: odd amplicons pool 1, even pool 2
main = {n: s for n, s in prim.items() if "_alt" not in n}
rows = []
for a, (s, e) in bed.items():
    f, r = main[f"{a}_F"], main[f"{a}_R"]
    mine = [x for n, x in main.items() if pool(n) == pool(a) and not n.startswith(a + "_")]
    others = [min(bp.get_heterodimer_dg(p, o) for o in mine) for p in (f, r)]
    rows.append(dict(amp=a, cov=cov[a], length=e - s, amp_gc=gc(ref[s:e]), tm_f=tm(f), tm_r=tm(r), tm_min=min(tm(f), tm(r)),
                     tm_diff=abs(tm(f) - tm(r)), gc_f=gc(f), gc_r=gc(r),
                     self_dg=min(bp.get_homodimer_dg(f), bp.get_homodimer_dg(r)), pair_dg=bp.get_heterodimer_dg(f, r),
                     end3=min(bp.get_end_stability_dg(f, f), bp.get_end_stability_dg(r, r)),
                     worst_cross_dg=min(others), len_f=len(f), len_r=len(r), n_alt=sum(1 for n in prim if n.startswith(a + "_alt") or n.startswith(a + "_alt"))))
y = np.log10([r["cov"] for r in rows])
def rank(x):
    o = np.argsort(x, kind="stable"); r = np.empty(len(x)); r[o] = np.arange(len(x))
    for v in np.unique(x):
        m = x == v
        if m.sum() > 1: r[m] = r[m].mean()
    return r
def spearman(x, y): return float(np.corrcoef(rank(np.asarray(x)), rank(y))[0, 1])
rnd = np.random.default_rng(0)
print(f"{len(rows)} amplicons; coverage {min(r['cov'] for r in rows):.0f}-{max(r['cov'] for r in rows):.0f}x")
print(f"{'metric':16}{'rho':>7}{'perm p':>9}")
out = []
for m in ("length", "amp_gc", "tm_min", "tm_diff", "gc_f", "gc_r", "self_dg", "pair_dg", "end3", "worst_cross_dg", "len_f", "len_r"):
    x = np.array([r[m] for r in rows]); rho = spearman(x, y)
    perm = np.mean([abs(spearman(x, rnd.permutation(y))) >= abs(rho) for _ in range(2000)])
    out.append((m, rho, perm)); print(f"{m:16}{rho:7.2f}{perm:9.3f}")
# how many ARTIC primers would pass OUR design limits?
ok = lambda s: 19 <= len(s) <= 21 and 59 <= tm(s) <= 61 and bp.get_homodimer_dg(s) >= -9
print("ARTIC primers passing our size/Tm/self-dimer limits:", sum(ok(s) for s in main.values()), "of", len(main))
print("ARTIC primer Tm (KAPA conditions): median", round(float(np.median([tm(s) for s in main.values()])), 1), "length median", int(np.median([len(s) for s in main.values()])))
low = sorted(rows, key=lambda r: r["cov"])[:10]
print("10 lowest amplicons:", [(r["amp"], round(r["cov"]), round(r["amp_gc"], 2), round(r["tm_min"], 1)) for r in low])
print("amplicon GC overall median", round(float(np.median([r["amp_gc"] for r in rows])), 2))
