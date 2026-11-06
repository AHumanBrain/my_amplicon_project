# SARS-CoV-2 benchmark (unmaintained)

A benchmark with ground truth: the public ARTIC v3 scheme (`examples/sars2/`) and a real MiSeq run of it (the data behind
`scripts/real_data_check.sh`). Run that script first and pass its workdir below. Not linted or tested.
The findings are now in the main README (Verified) and the fix is in `design.py` (`--primer-*` options).

| Script | What it does |
|---|---|
| `make_tiles.py` | (output is committed as `design/targets/sars2_tiles.tsv`) 99 overlapping 400 bp windows (step 300) over NC_045512.2 as a tiled targets table |
| `predict_dropouts.py <workdir>` | correlates our thermodynamic metrics for the ARTIC primers with the real per-amplicon coverage |
| `relaxed_design.py` | `design.py` with wider primer size/GC/Tm limits, to test why the default design fails |

## Results (2026-11-06; everything in silico except the real read counts)

**1. Our tiled design on SARS-CoV-2, default settings** (`design.py`, tiled layout, amplicons 330-400 bp, reference NC_045512.2):
66 of 99 tiles designed in 2 pools [36, 30]; 33 failed (8 had no primer pair at all, 25 no compatible pair in the tiled pool);
**11,680 bp (39% of the genome) uncovered after primer trimming**, in 51 stretches. ARTIC covers the genome with 98 amplicons.

**2. Why:** `design.py` fixes primers at 19-21 nt and 40-60% GC. The 196 ARTIC primers are 22-30 nt (median 24): none is 19-21 nt.
SARS-CoV-2 is 38% GC. With 20-30 nt, GC >= 30%, Tm 57-64 (`relaxed_design.py 20,30,30,57,64`): **90 of 99 designed, 3,467 bp uncovered** (48 stretches).
So primer size/GC limits caused most failures on an AT-rich genome; the rest of the gap is that tiling does not optimise neighbour overlap (a known open item).

**3. Do our metrics predict which real amplicons drop out?** No. Real per-amplicon coverage spans 44x-1522x (fold-80 5.6).
Spearman correlation of log coverage with 12 metrics (amplicon length and GC, primer Tm, Tm difference, primer GC and length,
self-dimer, F/R heterodimer, 3' end stability, worst cross-dimer in the pool): every |rho| <= 0.19, none significant at p < 0.05
(permutation test, 98 amplicons, no multiple-comparison correction; amplicon length rho = -0.19, p = 0.051 is the closest).
The 10 weakest amplicons (44-85x) look unremarkable on these metrics. One real dataset, one scheme, one chemistry (ARTIC uses
a different polymerase and annealing temperature than our KAPA model).
