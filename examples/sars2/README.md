# SARS-CoV-2 example (the development genome)

A 30 kb genome with a real, published tiled amplicon scheme and real reads: small enough that design, BLAST and the
full check run in about a minute, and it has ground truth for tiled designs and for coverage.

| File | What it is |
|---|---|
| `NC_045512.2.fa` | SARS-CoV-2 Wuhan-Hu-1 reference (NCBI RefSeq NC_045512.2, 29,903 bp, 38% GC), as shipped in the nf-core `test-datasets` repository |
| `artic_v3_primers.fasta` | the 196 ARTIC v3 primers (98 amplicons; alternates named `ampN_altX_F/R`), renamed `ampN_F` / `ampN_R` so `analysis/analyze.py` can use them |
| `artic_v3_amplicons.bed` | the 98 amplicons as intervals from the outer primer edges (`ampN`), derived from the ARTIC v3 primer BED |

Derived from the ARTIC v3 files in nf-core `test-datasets` (branch `viralrecon`, `genome/NC_045512.2/amplicon/`); these are
public data, check upstream terms before redistributing outside this repository.

Used by: `panel.example.toml` and `design/targets/sars2_{tiles,genes}.tsv` (designs), `scripts/check_all.sh`,
`scripts/real_data_check.sh` (a real MiSeq ARTIC v3 run, downloaded), and `research/benchmark_sars2/`.
