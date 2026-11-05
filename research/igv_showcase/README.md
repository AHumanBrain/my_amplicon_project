# IGV showcase (unmaintained)

Helper scripts used to make `docs/igv/*.png`; not tested or linted. Needs `npm i igv` (3.8.10 used), Playwright with Chromium, and a design plus its simulated and analysed BAMs.

1. `simulate_and_analyze.sh <design-prefix> <name>`: simulates one tube per pool and runs `analysis/analyze.py` (paths inside are session-specific; edit them).
2. Copy the genome (+ `.fai`), `genes.gff`, amplicons BED (from `scripts/design_viewer.py`), targets BED and BAMs or bedGraph coverage into `site/data/`, plus `igv.esm.min.js`, and write the track configs.
3. `python serve_range.py 8765` (a static server with HTTP Range support, which igv.js needs), then `node shot.js <outdir>`.

For IGV desktop, load `<prefix>_view/session.xml` from `design_viewer.py` (not tested here).
