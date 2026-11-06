# SARS-CoV-2 benchmark (unmaintained)

A benchmark with ground truth: the public ARTIC v3 scheme (`examples/sars2/`) and a real MiSeq run of it (the data behind
`scripts/real_data_check.sh`). Run that script first and pass its workdir below. Not linted or tested.
The findings are now in the main README (Verified) and the fix is in `design.py` (`--primer-*` options).

| Script | What it does |
|---|---|
| `make_tiles.py` | (output is committed as `design/targets/sars2_tiles.tsv`) 99 overlapping 400 bp windows (step 300) over NC_045512.2 as a tiled targets table |
| `predict_dropouts.py <workdir>` | correlates our thermodynamic metrics for the ARTIC primers with the real per-amplicon coverage |
| `coverage_report.py` | whole-genome coverage of a design once primers are trimmed, with each uncovered stretch labelled (failed tile / overlap gap / genome end); scores the ARTIC scheme the same way |

## Results (2026-11-06; everything in silico except the real read counts)

**1. Our tiled design on SARS-CoV-2, default settings** (`design.py`, tiled layout, amplicons 330-400 bp, reference NC_045512.2):
66 of 99 tiles designed in 2 pools [36, 30]; 33 failed (8 had no primer pair at all, 25 no compatible pair in the tiled pool);
**11,680 bp (39% of the genome) uncovered after primer trimming**, in 51 stretches. ARTIC covers the genome with 98 amplicons.

**2. Why:** `design.py` fixes primers at 19-21 nt and 40-60% GC. The 196 ARTIC primers are 22-30 nt (median 24): none is 19-21 nt.
SARS-CoV-2 is 38% GC. With 20-30 nt, GC >= 30%, Tm 57-64 (`--primer-min-size 20 --primer-max-size 30 --primer-min-gc 30 --primer-min-tm 57 --primer-max-tm 64`, now `panel.example.toml`): **90 of 99 designed, 3,467 bp uncovered** (48 stretches).
So primer size/GC limits caused most failures on an AT-rich genome; the rest of the gap is that tiling does not optimise neighbour overlap (a known open item).

**3. Do our metrics predict which real amplicons drop out?** No. Real per-amplicon coverage spans 44x-1522x (fold-80 5.6).
Spearman correlation of log coverage with 12 metrics (amplicon length and GC, primer Tm, Tm difference, primer GC and length,
self-dimer, F/R heterodimer, 3' end stability, worst cross-dimer in the pool): every |rho| <= 0.19, none significant at p < 0.05
(permutation test, 98 amplicons, no multiple-comparison correction; amplicon length rho = -0.19, p = 0.051 is the closest).
The 10 weakest amplicons (44-85x) look unremarkable on these metrics. One real dataset, one scheme, one chemistry (ARTIC uses
a different polymerase and annealing temperature than our KAPA model).

## Whole-genome panel: the complete problem (2026-11-06, in silico)

Goal: every base of SARS-CoV-2 (29,903 bp) inside some amplicon insert (amplicon minus primers) in at most 2 pools.
`coverage_report.py` scores a design that way. Targets: `design/targets/sars2_tiles.tsv` (99 windows of 400 bp, step 300).

| Design | Tiles designed | Genome covered | Uncovered bp: failed tiles / overlap gaps / ends | Bases covered 2+ times |
|---|---|---|---|---|
| ARTIC v3 scheme (reference point) | 98 of 98 | **99.6%** | 0 / 0 / 121 | 13% |
| This tool, default primer limits | 66 of 99 | 60.0% | 10,343 / 1,337 / 283 | 2% |
| This tool, `panel.example.toml` (wider primer limits) | 89 of 99 | 86.1% | 3,082 / 873 / 212 | 3% |
| Same, `--num-candidates 100` | 94 of 99 | 90.6% | 1,463 / 1,141 / 212 | 3% |

What stands between 90.6% and ARTIC's 99.6%, in order of how it showed up:

1. **Primer limits** (fixed by the `--primer-*` options): 60% to 86%.
2. **Pool compatibility.** With the wider limits every failed tile still has primer pairs, but none fits with the other
   amplicons in its pool. The rejections are all primer-dimer checks against primers already in the pool (heterodimer
   dG < -8 or 3' end dimer dG < -6; never cross-product or amplicon overlap). The dimer limits are absolute while the primers
   are now 20-30 nt, and each tile has only 25 candidates (the packer reached rank 24). More candidates help
   (89 to 94 tiles) but do not finish it.
3. **Tile placement.** Amplicons sit anywhere inside fixed windows, so neighbours often overlap by less than two primer
   lengths: the tool's own gap report shows about 1,100 bp of overlap gaps even when most tiles design, and only 3% of bases
   are covered twice. ARTIC overlaps neighbours on purpose (13%), so one lost amplicon does not open a hole.
4. **Genome ends:** about 200 bp uncovered (ARTIC 121 bp).

All numbers are predicted from the design alone. They say nothing about wet-lab dropout: in the real ARTIC run no amplicon
fell below 44x, even though our metrics could not predict which would be weakest.
